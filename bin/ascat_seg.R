{
    library(optparse)
    library(devtools)
    library(khtools)
    
    ## DO NOT FAIL SILENTLY
    options(error = function() {traceback(2); quit("no", 1)})

    ## prepared following these comments:
    ## https://github.com/VanLoo-lab/ascat/issues/19#issuecomment-371866049

    if (!exists('opt')) {
        option_list = list(
            make_option(c("-i", "--id"), type = "character", help = "Sample / pair ID (required)"),
            make_option("--libdir", type = "character", help = "libdir", default=""),
            make_option("--outdir", type = "character", help = "outdir", default="./"),
            make_option("--variants", type = "character", help = "path to variants file",
                        default="./"),
            make_option("--gc", type = "logical", help = "gc correct?", default=TRUE),
            make_option("--from_maf", type = "logical", help = "starting from MAF?", default=FALSE),
            make_option("--rebin_width", type = "numeric",
                        help = "rebin coverage width?", default=5e4),
            make_option("--snp_file", type = "character",
                        help = "path to file with SNVs",
                        default = ""),
            make_option("--coverage", type = "character", help = "path to coverage file",
                        default="./"),
            make_option("--seg", type = "character", help = "path to seg file",
                        default="/dev/null"),
            make_option("--field", type = "character", help = "field for ratio", default="ratio"),
            make_option("--penalty", type = "numeric", help = "penalty for aspcf", default=70),
            make_option("--hets_thresh", type = "numeric",
                        help = "only include sites with ALT allele fraction above this value",
                        default=0)
        )
        parseobj = OptionParser(option_list=option_list)
        opt = parse_args(parseobj)

        print(opt)
        
        print(.libPaths())
        options(error=function()traceback(2))

        ## keep record of run
        writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
        saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
    }

    library(ASCAT)
    library(data.table)
    library(gUtils)
    library(GenomicRanges)
    library(skitools)

    source("~/modules/ascat_seg/utils.R")

    if (grepl(pattern = "txt$", x = opt$variants)) {
        variants.dt = fread(opt$variants)
        variants.dt[, ":="(alt.count.n = as.numeric(as.character(alt.count.n)),
                           ref.count.n = as.numeric(as.character(ref.count.n)),
                           alt.count.t = as.numeric(as.character(alt.count.t)),
                           ref.count.t = as.numeric(as.character(ref.count.t)))]

        hets.dt = variants.dt[(alt.count.n > 0 | ref.count.n > 0) & (!is.na(alt.count.t)) & (!is.na(ref.count.t))]

        ## recompute totals and allele fractions
        hets.dt[, ":="(Chr = seqnames, Position = start,
                       normal.total = alt.count.n + ref.count.n,
                       tumor.total = alt.count.t + ref.count.t)]

        hets.dt[, ":="(alt.frac.t = alt.count.t / tumor.total,
                       ref.frac.t = ref.count.t / tumor.total,
                       alt.frac.n  = alt.count.n / normal.total,
                       ref.frac.n = ref.count.n / normal.total)]

        hets.dt[, ":="(normal.logR = 0)]
        hets.dt = hets.dt[(!is.infinite(tumor.total)) & (!is.na(tumor.total)) & (!is.na(normal.total))  & (!is.infinite(normal.total))]
        hets.dt[, tumor.adj := tumor.total / mean(tumor.total, na.rm = TRUE)]
        hets.dt[, normal.adj := normal.total / mean(normal.total, na.rm = TRUE)]
        hets.dt[, tumor.R := tumor.adj / normal.adj]
        hets.dt[, tumor.logR := log2(tumor.R / median(tumor.R, na.rm = TRUE))]
        hets.dt[, adjusted.ratio := tumor.total / normal.total]

    } else if (grepl(pattern = "rds$", x = opt$variants) || opt$from_maf) {
        variants.dt = as.data.table(readRDS(opt$variants))
        variants.dt[, alt.frac.n := 0.5]
        variants.dt[, ref.frac.n := 0.5]
        variants.dt[, alt.frac.t := alt.count.t / (ref.count.t + alt.count.t)]
        variants.dt[, ref.frac.t := ref.count.t / (ref.count.t + alt.count.t)]
        variants.dt[, Chr := seqnames]
        variants.dt[, Position := start]
        all.sites.dt = fread(opt$snp_file)
        all.sites.dt[, alt.frac.n := 0.99]
        all.sites.dt[, ref.frac.n := 0.01]
        all.sites.dt[, seqnames := Chr]
        all.sites.dt[, start := Position]
        all.sites.dt[, end := Position]
        hets.dt = rbind(variants.dt, all.sites.dt, fill = TRUE)
        hets.dt = unique(hets.dt, by = c("Chr", "Position"))
    } else {
        stop("either txt or rds must be supplied.")
    }

    ## read coverage file (expected to be .rds)
    cov.gr = readRDS(opt$coverage)


    ## ## transfer ratio
    message("Transferring ratio")
    ## gc correct tumor and normal
    if (opt$gc) {
        if ("tumor" %in% names(values(cov.gr)) & "normal" %in% names(values(cov.gr))) {
            message("Applying GC correction")
            tum.gr = khtools::.gc(cov.gr, "tumor")
            norm.gr = khtools::.gc(cov.gr, "normal")
            ratio.gr = khtools::.gc(cov.gr, "ratio")
            values(cov.gr)[, opt$field] = values(ratio.gr)[, "ratio"]
        } else {
            message("Skipping GC correction!!")
        }
    }
    if (!(any(grepl("^chr", hets.dt[, seqnames]), na.rm = TRUE))) {
        ov = gr.findoverlaps(GRanges(seqnames = hets.dt[, seqnames],
                                     ranges = IRanges(start = hets.dt[, start], width = 1)),
                             gr.nochr(cov.gr[, opt$field]),
                             return.type = "data.table")
    } else {
        ov = gr.findoverlaps(GRanges(seqnames = hets.dt[, seqnames],
                                     ranges = IRanges(start = hets.dt[, start], width = 1)),
                             cov.gr[, opt$field],
                             return.type = "data.table")
    }
    ov[, val := values(cov.gr)[, opt$field][subject.id]]
    hets.dt[, ratio := ov$val[match(1:.N, ov$query.id)]]


    message("Filtering SNP sites")
    good = hets.dt[, which((!is.infinite(ratio)) & (!is.na(ratio)) & ratio > 0)]
    hets.dt = hets.dt[good,]

    message("Normalizing log ratio...")
    ## try rerunning with mean instead of median?
    hets.dt[, ratio := ratio / mean(ratio, na.rm = TRUE)] ## should it be mean or median??

    normal.logr.dt = hets.dt[, .(Chr, Position, N = 0)]
    tumour.logr.dt = hets.dt[, .(Chr, Position, S = log2(ratio))]
    normal.BAF.dt = hets.dt[, .(Chr, Position, N = alt.frac.n)]
    tumour.BAF.dt = hets.dt[, .(Chr, Position, S = alt.frac.t)]

    ## get data frames
    normal.logr.df = as.data.frame(normal.logr.dt)
    tumour.logr.df = as.data.frame(tumour.logr.dt)
    normal.BAF.df = as.data.frame(normal.BAF.dt)
    tumour.BAF.df = as.data.frame(tumour.BAF.dt)

    ## set row names
    rownames(normal.logr.df) = paste0("SNP", 1:normal.logr.dt[, .N])
    rownames(tumour.logr.df) = paste0("SNP", 1:normal.logr.dt[, .N])
    rownames(normal.BAF.df) = paste0("SNP", 1:normal.logr.dt[, .N])
    rownames(tumour.BAF.df) = paste0("SNP", 1:normal.logr.dt[, .N])

    ## make temporary output directory
    message("Dumping temporary files...")
    temp.dir = file.path(opt$outdir, "tmp")
    dir.create(temp.dir, recursive = TRUE)

    ## dup files
    fwrite(normal.logr.df, file.path(temp.dir, "normal.logR.txt"), row.names = T, sep = "\t")
    fwrite(tumour.logr.df, file.path(temp.dir, "tumour.logR.txt"), row.names = T, sep = "\t")
    fwrite(normal.BAF.df, file.path(temp.dir, "normal.BAF.txt"), row.names = T, sep = "\t")
    fwrite(tumour.BAF.df, file.path(temp.dir, "tumour.BAF.txt"), row.names = T, sep = "\t")

    ##browser()

    ## grab gender
    gender = ifelse(hets.dt[(seqnames == "X") | (seqnames == "chrX"),
                            median(ratio, na.rm = T)] > 0.8,
                    "XX",
                    "XY")
    message("The gender of this sample: ", gender)

    message("Starting ASCAT!!")
    message("Reading LogR and BAF files")

    ## check if we need chromosome prefix
    if (!(any(grepl("^chr", hets.dt[, seqnames]), na.rm = TRUE))) {
        ascat.dat = ASCAT::ascat.loadData(Tumor_LogR_file = file.path(temp.dir, "tumour.logR.txt"),
                                          Tumor_BAF_file = file.path(temp.dir, "tumour.BAF.txt"),
                                          Germline_LogR_file = file.path(temp.dir, "normal.logR.txt"),
                                          Germline_BAF_file = file.path(temp.dir, "normal.BAF.txt"),
                                          gender = gender)
    } else {
        ascat.dat = ASCAT::ascat.loadData(Tumor_LogR_file = file.path(temp.dir, "tumour.logR.txt"),
                                      Tumor_BAF_file = file.path(temp.dir, "tumour.BAF.txt"),
                                      Germline_LogR_file = file.path(temp.dir, "normal.logR.txt"),
                                      Germline_BAF_file = file.path(temp.dir, "normal.BAF.txt"),
                                      gender = gender,
                                      chrs = paste0("chr", c(as.character(1:22), "X", "Y")))
    }

    ascat.plotRawData(ascat.dat, opt$outdir)

    ascat.dat = ASCAT::ascat.aspcf(ascat.dat, penalty = opt$penalty) ## higher penalty??

    ascat.plotSegmentedData(ascat.dat, opt$outdir)##"~/public_html")
    ascat.res = ASCAT::ascat.runAscat(ASCATobj = ascat.dat, gamma = 1)

    if (is.null(ascat.res$ploidy)) {
        stop("ASCAT failed to find a good purity/ploidy!")
    } else {
        message("ASCAT ploidy: ", ascat.res$ploidy)
        message("ASCAT cellularity: ", ascat.res$aberrantcellfraction)
    }

    message("Preparing output!")

    segments.gr = GRanges(seqnames = ascat.res$segments[, "chr"],
                      ranges = IRanges(start = ascat.res$segments[, "startpos"],
                                       end = ascat.res$segments[, "endpos"]),
                      nMajor = ascat.res$segments[, "nMajor"],
                      nMinor = ascat.res$segments[, "nMinor"],
                      nTotal = ascat.res$segments[, "nMajor"] + ascat.res$segments[, "nMinor"])

    pp.dt = data.table(id = opt$id,
                       purity = ascat.res$aberrantcellfraction,
                       ploidy = ascat.res$ploidy)

    saveRDS(pp.dt, "ascat_pp.rds")
    saveRDS(segments.gr, "ascat_seg.rds")

    message("Cleaning up temporary files")
    file.remove(file.path(temp.dir, "tumour.logR.txt"))
    file.remove(file.path(temp.dir, "normal.logR.txt"))
    file.remove(file.path(temp.dir, "tumour.BAF.txt"))
    file.remove(file.path(temp.dir, "normal.BAF.txt"))

    message("Done")

    quit("no")
}
