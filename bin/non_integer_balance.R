{
    library(optparse)
    options(bitmapType='cairo')
    ## options(error = function() {traceback(2); quit("no", 1)})

    if (!exists('opt'))
    {
      option_list = list(
        make_option("--id", type = "character", help = "sample id"),
        make_option("--jab", type = "character", help = "contains jabba_rds"),
        make_option("--cov", type = "character", help = "path to coverage file"),
        make_option("--field", type = "character", help = "signal field in cov", default = "foreground"),
        make_option("--hets", type = "character", help = "het pileups", default = "/dev/null"),
        make_option("--hets_thresh", type = "numeric", help = "threshold for calling a site a het",
                    default = 0.2),
        make_option("--mask", type = "character", help = "path to mask file"),
        make_option("--overwrite", type = "logical", help = "overwrite?", default = FALSE),
        make_option("--lambda", type = "numeric", default = 20, help = "slack penalty"),
        make_option("--allin", type = "logical", default = TRUE, help = "force all ALT?"),
        make_option("--fix_thresh", type = "numeric", default = 10,
                    help = "penalty threshold for fixing nodes"),
        make_option("--nodebounds", type = "logical", default = TRUE,
                    help = "set boundaries on heavy nodes?"),
        make_option("--ism", type = "logical", default = FALSE, help = "ism?"),
        make_option("--build", type = "character", default = "hg19", help = "either hg19 or hg38"),
        make_option("--epgap", type = "numeric", default = 1e-6, help = "epgap"),
        make_option("--tilim", type = "numeric", default = 100, help = "time limit for optimization"),
        make_option("--gurobi", type = "logical", default = TRUE, help = "use gurobi?"), ## because why the hell not
        make_option("--fasta", type = "character", help = "path to fasta file",
                    default = "/gpfs/commons/groups/imielinski_lab/DB/GATK/human_g1k_v37_decoy.fasta"),
        make_option("--pad", type = "numeric", help = "how much padding around junction for map", default = 101),
        make_option(c("-o", "--outdir"), type = "character", default = './', help = "Directory to dump output into")
        )
        parseobj = OptionParser(option_list=option_list)
        opt = parse_args(parseobj)
        saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
    }

    ## updated balance? yikes haha
    library(gUtils)
    library(skitools)
    library(JaBbA)
    library(gGnome)
    library(zitools)
    ## devtools::load_all("~/git/gGnome") ## noninteger capacity now added to master branch

    #' @name grab.hets
    #' @title grab.hets
    #'
    #' @description
    #'
    #' returns allele gtrack given sites.txt from het pileup
    #'
    #' @param agt.fname (character) path to sites.txt
    #' @param min.frac (numeric) between 0 and 1, min frequency in normal to count as het site
    #' @param max.frac (numeric) between 0 and 1, max frequency in normal to count as het site
    #'
    #' @return allele gTrack
    grab.hets = function(agt.fname = NULL,
                        min.frac = 0.2,
                        max.frac = 0.8)
    {
        if (is.null(agt.fname) || !file.exists(agt.fname)) {
            stop("agt.fname does not exist")
        }

        ## prepare and filter
        agt.dt = fread(agt.fname)[alt.frac.n > min.frac & alt.frac.n < max.frac,]
        ## add major and minor
        agt.dt[, which.major := ifelse(alt.count.t > ref.count.t, "alt", "ref")]
        agt.dt[, major.count := ifelse(which.major == "alt", alt.count.t, ref.count.t)]
        agt.dt[, minor.count := ifelse(which.major == "alt", ref.count.t, alt.count.t)]

        ## melt the data frame
        agt.melted = rbind(agt.dt[, .(seqnames, start, end, count = major.count, allele = "major")],
                        agt.dt[, .(seqnames, start, end, count = minor.count, allele = "minor")]
                        )

        ## make GRanges
        agt.gr = dt2gr(agt.melted[, .(seqnames, start, end, count, allele)])

        return (agt.gr)
    }

    #' @name grab.hets.from.maf
    #' @title grab.hets.from.maf
    #'
    #' @description
    #'
    #' get hets into format needed by phased.binstats from HMF maf_approx
    #'
    #' @param agt.fname (character)
    #' @param min.frac (numeric)
    #'
    #' @param GRanges
    grab.hets.from.maf = function(agt.fname, min.frac = 0.2) {

        if (!file.exists(agt.fname)) {
            stop("invalid file")
        }

        gr = readRDS(agt.fname)

        if (!inherits(gr, "GRanges")) {
            stop("rds file must contain GRanges object.")
        }

        dt = as.data.table(gr)
        dt[, major.count := ifelse(alt.count.t > ref.count.t, alt.count.t, ref.count.t)]
        dt[, minor.count := ifelse(alt.count.t > ref.count.t, ref.count.t, alt.count.t)]

        ## out.dt = rbind(dt[, .(seqnames, start, end, count = major.count, allele = "major")],
        ##                dt[, .(seqnames, start, end, count = minor.count, allele = "minor")])

        return(dt)
    }

    ## define output file names
    balanced.gg.fn = paste0(opt$outdir, '/', 'balanced.gg.rds')
    binstats.gg.fn = paste0(opt$outdir, '/', 'binstats.gg.rds')
    hets.gg.fn = paste0(opt$outdir, '/', 'hets.gg.rds')

    message("Reading gGraph!")
    jab = gG(jabba = opt$jab)

    message("Reading coverage!")
    cov = readRDS(opt$cov)

    message("Loading reference!")
    ref = RSeqLib::BWA(fasta = opt$fasta)

    if (file.exists(opt$mask) & file.info(opt$mask)$size) {
        message("Checking coverage mask")
        mask = readRDS(opt$mask)
        ## remove masked bins
        cov = cov %Q% (!(cov %^% mask))
    }

    ## add ncn field to bins
    ## first get gender of the sample from the JaBbA graph
    ncn.x = jab$nodes$dt[(seqnames == "X" | seqnames == "chrX"),
                         weighted.mean(cn,
                                       w = end - start + 1,
                                       na.rm = TRUE)]
    message("mean CN of X: ", ncn.x)
    ncn.vec = rep(2, length(cov))
    if (ncn.x < 1.4) {
        message("Adjusting ncn for XY")
        ncn.vec[which(as.character(seqnames(cov)) %in% c("chrX", "chrY", "X", "Y"))] = 1
    }

    values(cov)[, "ncn"] = ncn.vec

    if ((opt$overwrite) | (!file.exists(binstats.gg.fn))) {

        message("Starting binstats")
        binstats.gg = gGnome::binstats(jab, bins = cov, field = opt$field, lp = TRUE)

        ## save binstats
        saveRDS(binstats.gg, binstats.gg.fn)

    } else {
        binstats.gg = readRDS(binstats.gg.fn)
    }

    if (opt$overwrite | (!file.exists(balanced.gg.fn))) {

        if (opt$allin) {
            message("Setting edge lower bounds")

            ## identify only edges between "standard chromosomes"
            if (binstats.gg$edges$dt[type == "ALT", .N] > 0) {
                edt = binstats.gg$edges$dt[type == "ALT"]
                edt[, seqnames.1 := binstats.gg$nodes$dt[n1, seqnames]]
                edt[, seqnames.2 := binstats.gg$nodes$dt[n2, seqnames]]

                ## define standard chromosomes
                all.chrs = unique(seqnames(binstats.gg$nodes$gr))
                std.chrs = as.character(all.chrs)[grepl("^(chr)*[0-9XY]+$",
                                                        as.character(all.chrs))]
                message("chrs: ", std.chrs)

                ## get edges to fix
                if (edt[seqnames.1 %in% std.chrs &
                        seqnames.2 %in% std.chrs, .N] > 0) {
                    efix = edt[seqnames.1 %in% std.chrs &
                               seqnames.2 %in% std.chrs, edge.id]
                } else {
                    efix = numeric()
                }

                ## get edges to remove
                if (edt[(!seqnames.1 %in% std.chrs) |
                        (!seqnames.2 %in% std.chrs), .N] > 0) {
                    eremove = edt[(!seqnames.1 %in% std.chrs) |
                                  (!seqnames.2 %in% std.chrs), edge.id]
                } else {
                    eremove = numeric()
                }

                ## also get some more edges to remove
                ## only check standard chromosomes for mappability
                jj = binstats.gg$junctions[type == "ALT" & cn > 0 & edge.id %in% efix]

                if (length(jj)) {
                    message("Getting mappability...")
                    jj = zitools:::junction_mappability(junctions = jj,
                                              ref = ref,
                                              pad = opt$pad,
                                              build = opt$build,
                                              verbose = TRUE)
                    if (is.null(jj$dt$FILTER)) {
                        jj$set(FILTER = "PASS") ## no choice i guess
                    }
                    if (jj$dt[((bp1.mapq < 60) | (bp2.mapq < 60)) & (!FILTER == "PASS"), .N]) {
                        message("Number of additional junctions to remove: ",
                                jj$dt[((bp1.mapq < 60) | (bp2.mapq < 60)) & (!FILTER == "PASS"), .N])
                        eremove = c(eremove, jj$dt[((bp1.mapq < 60) | (bp2.mapq < 60)) & (!FILTER == "PASS"), edge.id])
                    }
                }

                ## don't fix any edges that should be removed
                efix = setdiff(efix, eremove)

                if (length(efix)) {
                    binstats.gg$edges[efix]$mark(lb = 1)
                    binstats.gg$junctions[efix]$set(lb = 1)
                } else {
                    message("no edges to fix!")
                }

                if (length(eremove)) {
                    binstats.gg$edges[eremove]$mark(ub = 0)
                    binstats.gg$junctions[eremove]$set(ub = 0)
                } else {
                    message("no edges to remove!")
                }
            } else {
                efix = numeric(); eremove = numeric()
            }

            message("Number of fixed edges: ", length(efix))
            message("Number of removed edges: ", length(eremove))
            binstats.gg$edges$mark(cn = NA)
            binstats.gg$junctions$set(cn = NA)
        } else {
            efix = NULL
            eremove = NULL
        }

        binstats.gg$nodes[cn > 999]$mark(cn = NA)


        message("Starting balance")
        res = balance(binstats.gg,
                      tilim = opt$tilim,
                      lambda = opt$lambda,
                      epgap = opt$epgap,
                      ##efix = ##efix,
                      lp = TRUE,
                      verbose = 2,
                      use.gurobi = opt$gurobi,
                      ism = opt$ism,
                      nonintegral = TRUE,
                      debug = TRUE,
                      trelim = 64000,
                      nodefileind = 3)

        message("Finished with balance!")

        balanced.gg = res$gg
        sol = res$sol
        balanced.gg$set(epgap = sol$epgap)
        balanced.gg$set(status = sol$status)

        ## save purity and ploidy
        message("Transfering purity and ploidy metadata")
        balanced.gg$set(purity = jab$meta$purity, ploidy = jab$meta$ploidy)

        balanced.gg = gGnome:::loosefix(balanced.gg)
        ## save results
        saveRDS(balanced.gg, balanced.gg.fn)


    } else {
        balanced.gg = readRDS(balanced.gg.fn)
    }

    ## if (opt$overwrite | (!file.exists(hets.gg.fn))) {
    ## fuck lol regenerate graphs with hets
    if (TRUE) {
        if (file.exists(opt$hets) & file.info(opt$hets)$size > 0) {

            message("Found hets!")
            if (grepl(pattern = "txt$", x = opt$hets)) {
                hets = grab.hets(opt$hets)

                hets.dt = fread(opt$hets)

                ## remove sites with zero counts
                hets.dt = hets.dt[(alt.count.n > 0 | ref.count.n > 0) & (alt.count.t > 0 | ref.count.t > 0),]

                ## filter hets for heterozygous sites
                ## mark with fractions
                message("Filtering sites for hets")
                hets.dt[, ":="(alt.frac.t = alt.count.t / (alt.count.t + ref.count.t),
                               ref.frac.t = ref.count.t / (alt.count.t + ref.count.t),
                               alt.frac.n = alt.count.n / (alt.count.n + ref.count.n),
                               ref.frac.n = ref.count.n / (alt.count.n + ref.count.n))]

                hets.dt = hets.dt[(alt.frac.n > opt$hets_thresh & ref.frac.n > opt$hets_thresh),]
                hets.dt[, major.frac := pmax(alt.frac.t, ref.frac.t)]
                hets.dt[, minor.frac := pmin(alt.frac.t, ref.frac.t)]
                hets.dt[, major.count := pmax(alt.count.t, ref.count.t)]
                hets.dt[, minor.count := pmin(alt.count.t, ref.count.t)]

            } else if (grepl(pattern = "rds$", x = opt$hets)) {
                hets.dt = grab.hets.from.maf(opt$hets)
            } else {
                stop("either txt or rds must be supplied.")
            }

            hets.gr = dt2gr(hets.dt[, .(seqnames, start, end, strand = "*", alt = alt.count.t, ref = ref.count.t)])

            jab = zitools:::gg2jab(balanced.gg, purity = balanced.gg$meta$purity, ploidy = balanced.gg$meta$ploidy)
            jab = c(jab,
                    JaBbA:::jabba.alleles(jab, hets.gr, verbose = TRUE, uncoupled = TRUE)[c("asegstats", "aadj", "agtrack")])

        } else {
            message("Hets not supplied! Skipping!")
            jab = zitools:::gg2jab(balanced.gg, purity = balanced.gg$meta$purity, ploidy = balanced.gg$meta$ploidy)
        }

        jab = gGnome:::loosefix(gG(jabba = jab))
        saveRDS(jab, hets.gg.fn)

    } else {
        message("DONE!")
    }

    quit("no", status = 0)
}
