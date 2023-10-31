withAutoprint(
{
    library(optparse)
    ## RLIBDIR = '/cga/meyerson/home/marcin/Software/R/x86_64-unknown-linux-gnu-library/3.1/'
    ## GIT.HOME = '/cga/meyerson/home/marcin/DB/git'

    option_list = list(
        make_option(c("-t", "--tbam"), type = "character", help = "Path to tumor .bam file"),
        make_option(c("-n", "--nbam"), type = "character", default = 'NULL', help = "Path to normal .bam file"),
        make_option(c("-s", "--sites"), type = "character", default = 1, help = "Path to .vcf file containing sites to probe eg hapmap, dbSNP"),
        make_option(c("-c", "--ncores"), type = "integer", default = 1, help = "Number of cores"),
        make_option(c("-l", "--libdir"), type = "character", default = paste(Sys.getenv('GIT_HOME'), 'isva', sep = '/'), help = "Directory containing karyoMIP.R file (eg default GIT.HOME/isva)"),
        make_option(c("-m", "--maxdepth"), type = "integer", default = 500, help = "maximum read depth for each SNP site (default 500)"),
        make_option(c("-o", "--outdir"), type = "character", default = './', help = "Directory to dump output into"),
        make_option(c("-k", "--filter"), type = "character", default = 'TRUE', help = "TRUE set to true if only keep candidate hets in output, otherwise will keep everything")
    )

    parseobj = OptionParser(option_list=option_list)
    opt = parse_args(parseobj)

    options(error = function() { traceback(2); quit("no", 1) })

    if (is.null(opt$tbam) | is.null(opt$nbam) | is.null(opt$sites))
        stop(print_help(parseobj))

    print(opt)

    clock = function(expr)
    {
        now = Sys.time()
        eval(expr)
        return(Sys.time()-now)
    }

                                        #Sys.setenv(GIT_HOME = GIT.HOME)
                                        #print(Sys.getenv('GIT_HOME'))

    print('LIBPATHS!');
    print(.libPaths())
    library(parallel)
    library(methods)
    library(abind)
    library(skitools)
    library(gUtils)
    library(bamUtils)
    require(Rsamtools)
    require(GenomicRanges)
    require(RColorBrewer)

    #source(paste0(opt$libdir, "/", "overwrite.R")) ## contains varcount which adjusts seqlevels in hapmap query if there is "chr" in bamfiles


    path = BiocGenerics::path

                                        #source(paste(CODE.HOME, 'functions.R', sep = "/"))
                                        #source(paste(opt$libdir, 'grUtils.R', sep = "/")) ## source local version of grUtils which contains the function hets (called below)

    ## mafcount = function (tum.bam, norm.bam = NULL, maf, chunk.size = 100, verbose = TRUE,
    ##           mc.cores = 1, ...)
    ## {
    ##     if (is.character(tum.bam)) {
    ##         if (!file.exists(tum.bam))
    ##             return("tumor bam file doesn't exist")
    ##         tum.bam = BamFile(tum.bam)
    ##     }
    ##     bams = BamFileList(tum.bam)
    ##     if (!is.null(norm.bam)) {
    ##         if (is.character(norm.bam)) {
    ##             if (!file.exists(norm.bam))
    ##                 return("normal bam file doesn't exist")
    ##             norm.bam = BamFile(norm.bam)
    ##         }
    ##         if (identical(seqlengths(bams), seqlengths(norm.bam))) {
    ##             bams = c(bams, BamFileList(norm.bam))
    ##         }
    ##         else {
    ##             bams2 = BamFileList(norm.bam)
    ##         }
    ##     }
    ##     if (is.null(dim(maf))) {
    ##         chunks = chunk(1, length(maf), chunk.size)
    ##     }
    ##     else {
    ##         chunks = chunk(1, nrow(maf), chunk.size)
    ##     }
    ##     if (is.null(maf$Tumor_Seq_Allele1)) {
    ##         maf$Tumor_Seq_Allele1 = maf$alt
    ##     }
    ##     if (is.null(maf$Tumor_Seq_Allele1)) {
    ##         maf$Tumor_Seq_Allele1 = maf$ALT
    ##     }
    ##     if (is.null(maf$Reference_Allele)) {
    ##         maf$Reference_Allele = maf$ref
    ##     }
    ##     if (is.null(maf$Reference_Allele)) {
    ##         maf$Reference_Allele = maf$REF
    ##     }
    ##     if (is.null(maf$Reference_Allele) | is.null(maf$Tumor_Seq_Allele1)) {
    ##         stop("Error: Cannot locate variant columns in input GRanges, please check input to make sure it either has standard VCF ALT / REF columns or MAF file columns specifying alt and ref allele")
    ##     }
    ##     if (!all(is.character(maf$Tumor_Seq_Allele1))) {
    ##         maf$Tumor_Seq_Allele1 = sapply(maf$Tumor_Seq_Allele1,
    ##                                        function(x) as.character(x)[1])
    ##     }
    ##     if (!all(is.character(maf$Reference_Allele))) {
    ##         maf$Reference_Allele = as.character(maf$Reference_Allele)
    ##     }
    ##     maf$alt.count.t = maf$ref.count.t = NA
    ##     if (!is.null(norm.bam)) {
    ##         maf$alt.count.n = maf$ref.count.n = NA
    ##     }
    ##     if (verbose) {
    ##         cat("Initialized\n")
    ##     }
    ##     if (is.data.frame(maf)) {
    ##         maf = seg2gr(maf)
    ##     }
    ##     tmp = do.call("rbind", mclapply(1:nrow(chunks), function(i) {
    ##         if (verbose) {
    ##             cat("Starting chunk ", chunks[i, 1], " to ", chunks[i,
    ##                                                                 2], "\n")
    ##         }
    ##         ix = chunks[i, 1]:chunks[i, 2]
    ##         if (verbose) {
    ##             now = Sys.time()
    ##         }
    ##         vc = varcount(bams, maf[ix], ...)
    ##         if (exists("bams2")) {
    ##             vc2 = varcount(bams2, maf[ix], ...)
    ##         }
    ##         if (verbose) {
    ##             print(Sys.time() - now)
    ##         }
    ##         tum.count = vc$counts[, , 1]
    ##         if (exists("bams2")) {
    ##             norm.count = vc2$counts[, , 1]
    ##         }
    ##         if (is.null(dim(tum.count))) {
    ##             tum.count = cbind(tum.count)
    ##         }
    ##         out = cbind(tum.count[cbind(match(maf$Tumor_Seq_Allele1[ix],
    ##                                           rownames(tum.count)), 1:length(ix))], tum.count[cbind(match(maf$Reference_Allele[ix],
    ##                                                                                                       rownames(tum.count)), 1:length(ix))])
    ##         if (verbose) {
    ##             cat("Num rows:", nrow(out), "\n")
    ##         }
    ##         if (!is.null(norm.bam)) {
    ##             if (identical(seqlengths(bams), seqlengths(norm.bam))) {
    ##                 norm.count = vc$counts[, , 2]
    ##             }
    ##             else {
    ##                 norm.count = vc2$counts[, , 1]
    ##             }
    ##             if (is.null(dim(norm.count))) {
    ##                 norm.count = cbind(norm.count)
    ##             }
    ##             out = cbind(out, norm.count[cbind(match(maf$Tumor_Seq_Allele1[ix],
    ##                                                     rownames(norm.count)), 1:length(ix))], norm.count[cbind(match(maf$Reference_Allele[ix],
    ##                                                                                                                   rownames(norm.count)), 1:length(ix))])
    ##         }
    ##         return(out)
    ##     }, mc.cores = mc.cores))
    ##     if (all(is.na(tmp))) {
    ##         return(GRanges())
    ##     }
    ##     else {
    ##         maf$alt.count.t = tmp[, 1]
    ##         maf$ref.count.t = tmp[, 2]
    ##         maf$alt.frac.t = maf$alt.count.t/(maf$alt.count.t + maf$ref.count.t)
    ##         maf$ref.frac.t = 1 - maf$alt.frac.t
    ##         if (!is.null(norm.bam)) {
    ##             maf$alt.count.n = tmp[, 3]
    ##             maf$ref.count.n = tmp[, 4]
    ##             maf$alt.frac.n = maf$alt.count.n/(maf$alt.count.n +
    ##                                               maf$ref.count.n)
    ##             maf$ref.frac.n = 1 - maf$alt.frac.n
    ##         }
    ##         return(maf)
    ##     }
    ## }


#############################
                                        # hets
                                        #
                                        # generates allele fraction at all possible hets at sites specified by vcf (eg hapmap) input
                                        # for tumor and normal
                                        #
#############################
    hets = function(tum.bam, norm.bam = NULL, out.file, vcf.file = '/cga/meyerson/home/marcin/DB/dbSNP/hapmap_3.3.b37.vcf', chunk.size1 = 1e3, chunk.size2 = 1e2, mc.cores = 14, verbose = T, na.rm = TRUE, max.depth = 500,
                    filt.norm = T ## if TRUE will remove any sites that have allele fraction 0 or 1 or NA in MAF
                    )
    {
        f = file(vcf.file, 'r')
        if (grepl('VCF', readLines(f, 1)))
            vcf = TRUE
        else
            vcf = FALSE

        sl = hg_seqlengths()

        if (verbose)
            st = Sys.time()

        nprocessed = 0
        nhets = 0
        first = T
        ## get past headers
        while (grepl('^#', last.line <<- readLines(f, n=1))){}

        if (verbose)
            cat('Opened vcf, writing hets to text file', out.file, '\n')

        out.cols = c('seqnames', 'start', 'end', 'Tumor_Seq_Allele1', 'Reference_Allele', 'ref.count.t', 'alt.count.t', 'ref.count.n', 'alt.count.n', 'alt.frac.t', 'ref.frac.t', 'alt.frac.n', 'ref.frac.n')


        if (vcf)
            col.ix = 1:5
        else
        {
            col.ix = match(c("Chromosome", "Start_position", "End_position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"), strsplit(last.line, '\t')[[1]])
            if (any(is.na(col.ix)))
                stop('Error processing variant file: must be valid VCF or MAF')
        }

        while (!is.null(tmp <- tryCatch(read.delim(file = f, as.is = T, header = F, nrows = chunk.size1)[, col.ix], error = function(x) NULL)))
        {
            if (vcf)
                names(tmp) = c('chr', 'start', 'name', 'ref', 'alt')
            else
            {
                names(tmp) = c('chr', 'start', 'name', 'ref', 'alt', 'alt2')
                ## just in case the first tumor seq allele is equal to reference .. which happens in mafs
                tmp$alt = ifelse(tmp$alt==tmp$ref, tmp$alt2, tmp$alt)
            }

            loc = seg2gr(tmp, seqlengths = sl)
            clock({loc.count = bamUtils::mafcount(tum.bam, norm.bam, loc, indel = T, chunk.size = chunk.size2, mc.cores = mc.cores, max.depth = max.depth)})
            nprocessed = nprocessed + length(loc.count)
            if (filt.norm & !is.null(loc.count$alt.frac.n))
                loc.count = loc.count[which(loc.count$alt.frac.n != 1 & loc.count$alt.frac.n != 0)]

            nhets = nhets + length(loc.count)
            if (length(loc.count)>0)
            {
                df = as.data.frame(loc.count)
                if (na.rm) ## remove any entries with 0 ref or alt reads in tumor or normal
                {
                    if (!is.null(norm.bam))
                        naix = apply(df[, c('alt.count.t', 'ref.count.t', 'alt.count.n', 'ref.count.n')], 1, function(x) all(is.na(x)))
                    else
                        naix = apply(df[, c('alt.count.t', 'ref.count.t')], 1, function(x) all(is.na(x)))
                    df = df[which(!naix), ]
                }
                out.cols = intersect(out.cols, names(df))
                if (first)
                {

                    write.tab(df[, out.cols], out.file, append = F, col.names = T)
                    first = F
                }
                else
                    write.tab(df[, out.cols], out.file, append = T, col.names = F)
            }

            if (verbose)
                cat(sprintf('Processed %s sites, wrote %s candidate hets\n', nprocessed, nhets))

            if (verbose)
            {
                cat('Time elapsed:\n')
                print(Sys.time() - st)
            }
        }

        close(f)

        if (verbose)
            cat('Finished het processing wrote to file', out.file, '\n')
    }


    ## keep record of run w opts
    writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'args.txt', sep = '/'))
    saveRDS(opt, paste(opt$outdir, 'opts.rds', sep = '/'))

    if (!is.null(opt$nbam))
    {
        if (!file.exists(opt$nbam))
            opt$nbam = NULL
        else
        {
            nbami = paste(opt$nbam, '.bai', sep = '')
            if (!file.exists(nbami))
                nbami = gsub('\\.bam$', '.bai', opt$nbam)
            opt$nbam = BamFile(opt$nbam, index = nbami)
        }
    }


    tbami = paste(opt$tbam, '.bai', sep = '')
    if (!file.exists(tbami))
        tbami = gsub('\\.bam$', '.bai', opt$tbam)
    opt$tbam = BamFile(opt$tbam, index = tbami)

    set.seed(42);

    hets(tum.bam = opt$tbam, norm.bam = opt$nbam, out.file = paste(opt$outdir, 'sites.txt', sep = '/'), vcf.file = opt$sites, chunk.size1 = 1e4, chunk.size2 = 2.5e3, max.depth = opt$maxdepth, mc.cores = opt$ncores, filt.norm = grepl('true', opt$filter, ignore.case = T))



    cat('Done\n')


    quit("no", 0)
}, echo = FALSE)
