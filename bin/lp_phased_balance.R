{
    library(optparse)
    options(bitmapType='cairo')
    ## options(error = function() {traceback(2); quit("no", 1)})

    if (!exists('opt'))
    {
      option_list = list(
        make_option("--id", type = "character", help = "sample id"),
        make_option("--jab", type = "character", help = "contains jabba_rds"),
        make_option("--hets", type = "character", help = "hets_pileup_wgs"),
        make_option("--lambda", type = "numeric", default = 100, help = "slack penalty"),
        make_option("--cnloh", type = "logical", default = FALSE, help = "allow cnloh?"),
        make_option("--major", type = "logical", default = TRUE, help = "force major?"),
        make_option("--allin", type = "logical", default = TRUE, help = "force all ALT?"),
        make_option("--marginal", type = "logical", default = TRUE, help = "force marginal?"),
        make_option("--from_maf", type = "logical", default = FALSE, help = "maf format? (HMF)"),
        make_option("--mask", type = "character", default = "/dev/null", help = "path to coverage mask"),
        make_option("--ism", type = "logical", default = TRUE, help = "ism?"),
        make_option("--epgap", type = "numeric", default = 1e-3, help = "epgap"),
        make_option("--hets_thresh", type = "numeric", default = 0.2, help = "threshold for c~alling het"),
        make_option("--min_bins", type = "numeric", default = 5, help = "min bins for marking node as unphased"),
        make_option("--min_width", type = "numeric", default = 5e3, help = "min width for marking node as unphased"),
        make_option("--trelim", type = "numeric", default = 64000, help = "max uncompressed tree size (MB)"),
        make_option("--reward", type = "numeric", default = 10, help = "penalty on CNLOH"),
        make_option("--nodefileind", type = "numeric", default = 1, help = "node file indicator (1 to keep in memory, 3 to write to disk)"),
        make_option("--tilim", type = "numeric", default = TRUE, help = "time limit for optimization"),
        make_option(c("-o", "--outdir"), type = "character", default = './', help = "Directory to dump output into")
        )
      parseobj = OptionParser(option_list=option_list)
        opt = parse_args(parseobj)
        saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
    }


    ## updated balance? yikes haha
    library(JaBbA)
    devtools::load_all("~/git/gGnome_ZC")
    devtools::load_all("~/git/zitools")
    library(gUtils)
    library(skitools)
    library(DNAcopy)


    ## the following are neede for grab.hets and grab.hets.from.maf but not doing that for now..
    ## source("~/projects/phasing/phasing_utils.R")
    ## source("~/projects/phasing/tmp.R")

    source("~/utils.R")

    message("Reading input!")

    ## read input
    jab = readRDS(opt$jab)
    if (inherits(jab, 'gGraph')) {
        message("JaBbA is gGraph")
        gg = readRDS(opt$jab)
    } else {
        gg = gG(jabba = opt$jab)
    }

    ## ## check whether hets or MAF approximation was provided
    if (grepl(pattern = "txt$", x = opt$hets)) {
        hets.dt = fread(opt$hets)
    } else if (grepl(pattern = "rds$", x = opt$hets)) {
        hets.dt = grab.hets.from.maf(opt$hets)
    } else {
        stop("either txt or rds must be supplied.")
    }

    ## prepare hets for jabba.alleles format
    ## hets.dt = fread(opt$hets)



    ## filter hets for heterozygous sites
    ## mark with fractions
    message("Filtering sites for hets")
    if ((!is.null(hets.dt$alt.count.n)) & (!is.null(hets.dt$ref.count.n))) {
        ## remove sites with zero counts
        hets.dt[, alt.count.n := as.numeric(as.character(alt.count.n))]
        hets.dt[, ref.count.n := as.numeric(as.character(ref.count.n))]
        hets.dt[, alt.count.t := as.numeric(as.character(alt.count.t))]
        hets.dt[, ref.count.t := as.numeric(as.character(ref.count.t))]
        hets.dt = hets.dt[(alt.count.n > 0 | ref.count.n > 0) & (alt.count.t > 0 | ref.count.t > 0),]
        hets.dt[, ":="(alt.frac.t = alt.count.t / (alt.count.t + ref.count.t),
                       ref.frac.t = ref.count.t / (alt.count.t + ref.count.t),
                       alt.frac.n = alt.count.n / (alt.count.n + ref.count.n),
                       ref.frac.n = ref.count.n / (alt.count.n + ref.count.n))]

        hets.dt = hets.dt[(alt.frac.n > opt$hets_thresh & ref.frac.n > opt$hets_thresh),]
    } else {
        hets.dt[, alt.count.t := as.numeric(as.character(alt.count.t))]
        hets.dt[, ref.count.t := as.numeric(as.character(ref.count.t))]
        hets.dt = hets.dt[(alt.count.t > 0 | ref.count.t > 0),]
        hets.dt [, ":="(alt.frac.t = alt.count.t / (alt.count.t + ref.count.t),
                        ref.frac.t = ref.count.t / (alt.count.t + ref.count.t))]
    }
    hets.dt[, major.frac := pmax(alt.frac.t, ref.frac.t)]
    hets.dt[, minor.frac := pmin(alt.frac.t, ref.frac.t)]
    hets.dt[, major.count := pmax(alt.count.t, ref.count.t)]
    hets.dt[, minor.count := pmin(alt.count.t, ref.count.t)]

    hets.gr = dt2gr(hets.dt[, .(seqnames, start, end, strand = "*", alt = alt.count.t, ref = ref.count.t)])

    ## get rid of masked SNPs
    if (file.exists(opt$mask) && file.info(opt$mask)$size) {
        message("reading mask!!")
        mask = readRDS(opt$mask)
        sel = which(!(hets.gr %^% mask))
        hets.gr = hets.gr[sel]
        hets.dt = hets.dt[sel,]
    } else {
        message("not using mask!!")
    }

    ## compute the BAF
    hets.dt[, BAF := minor.count / (minor.count + major.count)]

    ## start segmentation
    ## we want to only segment nodes with at least 10 SNPs
    ## and with CN > 0
    ## get stretches of constant total CN
    if (is.null(gg$meta$ploidy)) {
        max.cn = 3
    } else {
        max.cn = ceiling(gg$meta$ploidy + 1)
    }

    nds = gg$nodes$gr[, "cn"] %Q% (width(gg$nodes$gr) > 1e7) %Q% (cn > 0) %Q% (cn <= max.cn)
    spl = unlist(reduce(split(nds, ~ cn)))
    names(spl) = NULL

    nodes.ov.hets = gr.findoverlaps(spl,
                                    dt2gr(hets.dt[, .(seqnames, start, end, strand = "*")]),
                                    return.type = "data.table")

    nodes.ov.hets[, count := .N, by = query.id]

    ## we should only segment relatively wide nodes with low copy number
    nodes.ov.hets = nodes.ov.hets[count >= 100,]

    nodes.to.segment = nodes.ov.hets[, unique(query.id)]

    new.sl = seqlengths(gg)

    segs = lapply(nodes.to.segment,
                  function(qid) {
                      message("Starting segmentation for : ", qid)
                      hets.subset.dt = hets.dt[nodes.ov.hets[query.id == qid, subject.id],]
                      cna = CNA(hets.subset.dt[, BAF],
                                hets.subset.dt[, as.character(seqnames)],
                                hets.subset.dt[, start],
                                data.type = "logratio")
                      seg = segment(smooth.CNA(cna),
                                    alpha = 1e-5,
                                    verbose = TRUE)
                      utils::capture.output({seg_dt = print(seg); setDT(seg_dt)},
                                            type = "output",
                                            file = "/dev/null")
                      out = seg2gr(seg_dt[!(is.na(seg.mean) | is.na(loc.start) | is.na(loc.end))], new.sl)
                      out = gr.fix(out, new.sl, drop = T)
                      ## get the number of hets per segment
                      values(out)[, "nhets"] = out %N% hets.gr
                      ## make sure the segment is on the order of high kbp
                      out = out %Q% (nhets > 50)
                      message("Number of segments: ", length(out))
                      names(out) = NULL
                      if (length(out) > 1) {
                          return(gr.start(out[2:length(out)]))
                      }
                      return(GRanges())
                      message("Finished!")
                  })

    segs = do.call(grbind, segs)

    if (is.null(segs)) {
        message("No extra breakends")
        cnloh.bnds = GRanges()
    } else {
        message("Number of extra breakends: ", length(segs))
        cnloh.bnds = segs
    }


    ## add CNLOH breakends ## it doesn't help if there's too many of these
    if (length(cnloh.bnds) > 0) {
        ## create junctions corresponding to CNLOH breakpoints
        message("Creating CNLOH junctions")
        cnloh.bp1 = GRanges(seqnames = seqnames(cnloh.bnds),
                            ranges = IRanges(start = GenomicRanges::start(cnloh.bnds),
                                             width = 1),
                            strand = "-")
        cnloh.bp2 = GRanges(seqnames = seqnames(cnloh.bnds),
                            ranges = IRanges(start = GenomicRanges::start(cnloh.bnds) + 1,
                                             width = 1),
                            strand = "+")
        cnloh.grl = grl.pivot(GRangesList(cnloh.bp1, cnloh.bp2))
        cnloh.jj = jJ(cnloh.grl)
        cnloh.jj$set(type = "ALT", cnloh = TRUE)
        gg$edges$mark(cnloh = FALSE)
        gg = gg$copy$add(junctions = cnloh.jj)
    } else {
        message("No CNLOH junctions!")
        gg$edges$mark(cnloh = FALSE)
    }

    if(!opt$from_maf) {
        jab = zitools:::gg2jab(gg, purity = gg$meta$purity, ploidy = gg$meta$ploidy)
        jab = jabba.alleles2(jab, hets.gr, verbose = TRUE, uncoupled = TRUE, marginal = opt$marginal)

        ## transfer back the allelic annotations
        new.nodes = gg$nodes$gr[, "node.id"] %$% (jab$segstats %Q% (snode.id > 0))
        new.gg = gG(nodes = new.nodes, edges = gg$edges$dt)
        binstats.gg = diploid2haploid(new.gg)

        ## transfer back edge rewards
        ##binstats.gg$edges$mark(reward = new.gg$edges$dt[match(binstats.gg$edges$dt$og.edge.id, edge.id), reward])
        binstats.gg$edges$mark(cnloh = new.gg$edges$dt[match(binstats.gg$edges$dt$og.edge.id, edge.id), cnloh])
        ## binstats.gg$edges$mark(cnloh = gg$edges$dt[match(binstats.gg$edges$dt$og.edge.id, edge.id), cnloh])
    } else {
        tmp = as.data.table(readRDS(opt$hets))
        tmp = rbind(tmp[, .(seqnames, start, end, strand = "*",
                            allele = "major",
                            count = ifelse(baf >= 0.5, baf, 1 - baf))],
                    tmp[, .(seqnames, start, end, strand = "*",
                            allele = "minor",
                            count = ifelse(baf >= 0.5, 1 - baf, baf))])
        tmp.gr = dt2gr(tmp)
        if (file.exists(opt$mask) && file.info(opt$mask)$size) {
            message("reading mask!!")
            mask = readRDS(opt$mask)
            tmp.gr = tmp.gr %Q% (!tmp.gr %^% mask)
        }
        binstats.gg = phased.binstats(gg, bins = tmp.gr,
                                      purity = gg$meta$purity,
                                      ploidy = gg$meta$ploidy,
                                      count.field = "count",
                                      allele.field = "allele")
    }

    binstats.gg$edges[(cnloh)]$mark(reward = -opt$lambda * 1.5, type = "ALT")

    message("Generating marginals")

    marginal.gr = gg$nodes$gr[, "cn"]
    if (opt$marginal) {
        ##marginal.gr$fix = as.numeric(width(marginal.gr) > 1e7)
        values(marginal.gr)[, "nhets"] = marginal.gr %N% hets.gr
        values(marginal.gr)[, "fix"] = FALSE##as.numeric(marginal.gr$nhets > 100) ## more than 100 snps
        ## values(marginal.gr)[, "fix"] = TRUE
        values(marginal.gr)[, "weight"] = ifelse(values(marginal.gr)[, "nhets"] > 10, marginal.gr$nhets, 1)
        ## for non-fixed marginals, make the weight equal to number of hets
    } else {
        marginal.gr$fix = 0
    }

    ## M needs to be at least as big as the biggest node...
    binstats.gg$nodes[cn >= 999]$mark(cn = NA) ## NA anything bigger than M
    binstats.gg$edges$mark(cn = NA) ## no edge CN

    ## stash binstats gg
    message("Stashing binstats.gg")
    saveRDS(binstats.gg, paste0(opt$outdir, "/", "binstats.gg.rds"))

    message("Starting balance")


    res = balance(binstats.gg,
                  lambda = opt$lambda,
                  marginal = marginal.gr,
                  ism = opt$ism, ## false for now? because should be TRUE just by virtue of parent graph TRU
                  lp = TRUE,
                  M = 1000, ## importantly needs to be big enough or else this will be infeasible
                  verbose = 2,
                  tilim = opt$tilim,
                  epgap = opt$epgap,
                  cnloh = opt$cnloh,
                  force.major = opt$major,
                  force.alt = opt$allin,
                  phased = TRUE,
                  trelim = 128000,
                  nodefileind = opt$nodefileind,
                  debug = TRUE)

    balanced.gg = res$gg
    sol = res$sol
    balanced.gg$nodes$mark(epgap = sol$epgap, status = sol$status)
    message("Number of cnloh edges: ", length(balanced.gg$junctions[(cnloh) & cn > 0]))
    saveRDS(balanced.gg, paste0(opt$outdir, '/', 'balanced.gg.rds'))

    message("Marking nodes in unphased graph")
    allele.annotated.gg = copy(gg)
    major.nodes.gr = res$gg$gr %Q% (allele == "major")
    minor.nodes.gr = res$gg$gr %Q% (allele == "minor")
    major.pmt = match(values(allele.annotated.gg$nodes$gr)[, "node.id"], values(major.nodes.gr)[, "og.node.id"])
    minor.pmt = match(values(allele.annotated.gg$nodes$gr)[, "node.id"], values(minor.nodes.gr)[, "og.node.id"])
    allele.annotated.gg$nodes$mark(cn.high = values(major.nodes.gr)[major.pmt, "cn"],
                                   cn.low = values(minor.nodes.gr)[minor.pmt, "cn"],
                                   cn.old.high = values(major.nodes.gr)[major.pmt, "cn.old"],
                                   cn.old.low = values(minor.nodes.gr)[minor.pmt, "cn.old"])

    saveRDS(allele.annotated.gg, paste0(opt$outdir, "/", "unphased.gg.rds"))
    quit("no", status = 0)
}
