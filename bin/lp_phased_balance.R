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
    library(zitools)
    library(gGnome)
    library(gUtils)
    library(skitools)
    library(DNAcopy)


    # utils.R has been pasted below
    # source("~/utils.R")

    ## start of utils.R

    #' @name diploid2haploid
    #' @title diploid2haploid
    #'
    #' @param gg (gGraph) diploid gGraph (expect cn.low and cn.high as node metadata)
    #' @param verbose (logical) default FALSE
    #'
    #' @return gGraph with field allele and cn
    diploid2haploid = function(gg, verbose = FALSE) {
        og.nodes.gr = gg$nodes$gr[, c("cn.low", "cn.high", "var.low", "var.high", "cn", "node.id", "nhets")]
        values(og.nodes.gr)[, "og.node.id"] = values(gg$nodes$gr)[, "node.id"]
        names(values(og.nodes.gr))[names(values(og.nodes.gr)) == "cn"] = "cn.total"

        ## prepare nodes for melted graph
        phased.gg.nodes = c(og.nodes.gr, og.nodes.gr)
        values(phased.gg.nodes)[, "cn"] = c(values(og.nodes.gr)[, "cn.high"], values(og.nodes.gr)[, "cn.low"])
        values(phased.gg.nodes)[, "allele"] = c(rep("major", length(og.nodes.gr)), rep("minor", length(og.nodes.gr)))
        values(phased.gg.nodes)[, "variance"] = c(values(og.nodes.gr)[, "var.high"], values(og.nodes.gr)[, "var.low"])
        values(phased.gg.nodes)[, "variance"] = c(values(og.nodes.gr)[, "var.high"], values(og.nodes.gr)[, "var.low"])
        ## bound the variance
        values(phased.gg.nodes)[, "var.adj"] = pmax(pmax(values(phased.gg.nodes)[, "variance"], values(phased.gg.nodes)[, "cn"]), 1)

        ## compute the weight
        values(phased.gg.nodes)[, "weight"] = values(phased.gg.nodes)[, "nhets"] / values(phased.gg.nodes)[, "var.adj"]

        ## prepare edges for melted graph
        phased.gg.edges = rbind(
            gg$edges$dt[, .(n1, n2, n1.side, n2.side, type,
                og.edge.id = edge.id,
                n1.allele = "major",
                n2.allele = "major")],
            gg$edges$dt[, .(n1 = n1 + length(og.nodes.gr), n2 = n2 + length(og.nodes.gr), type,
                n1.side, n2.side,
                og.edge.id = edge.id,
                n1.allele = "minor",
                n2.allele = "minor")],
            gg$edges$dt[, .(n1, n2 = n2 + length(og.nodes.gr), type,
                n1.side, n2.side,
                og.edge.id = edge.id,
                n1.allele = "major",
                n2.allele = "minor")],
            gg$edges$dt[, .(n1 = n1 + length(og.nodes.gr), n2, type,
                n1.side, n2.side,
                og.edge.id = edge.id,
                n1.allele = "minor",
                n2.allele = "major")]
        )

        ## add n1/n2 chromosome information
        phased.gg.edges[, ":="(n1.chr = seqnames(phased.gg.nodes)[n1] %>% as.character,
            n2.chr = seqnames(phased.gg.nodes)[n2] %>% as.character)]

        ## add edge connection type (straight/cross)
        phased.gg.edges[n1.chr == n2.chr & n1.allele == n2.allele, connection := "straight"]
        phased.gg.edges[n1.chr == n2.chr & n1.allele != n2.allele, connection := "cross"]

        phased.gg = gG(nodes = phased.gg.nodes, edges = phased.gg.edges)

        ref.edge.col = alpha("blue", 0.3)
        alt.edge.col = alpha("red", 0.3)
        ref.edge.lwd = 0.5
        alt.edge.lwd = 1.0
        phased.gg$edges$mark(col = ifelse(phased.gg$edges$dt$type == "REF", ref.edge.col, alt.edge.col),
            lwd = ifelse(phased.gg$edges$dt$type == "REF", ref.edge.lwd, alt.edge.lwd))

        major.node.col = alpha("red", 0.5)
        minor.node.col = alpha("blue", 0.5)
        phased.gg$nodes$mark(col = ifelse(phased.gg$nodes$dt$allele == "major", major.node.col, minor.node.col),
            ywid = 0.8)

        return(phased.gg)
    }

    ####################
    #' @name jabba.alleles2
    #' @title jabba.alleles2
    #' @rdname internal
    #' jabba.alleles
    #'
    #' @description
    #' Populates allelic value s for JaBbA object.  This does not explicitly impose junction balance constraints on alleles, but rather just computes
    #' the maximum likelihood estimate given allelic counts and the inferred total copy number on a given segment according to JaBbA
    #'
    #' @param jab JaBbA object
    #' @param het.sites GRanges with meta data fields (see below) for alt and rref count
    #' @param alt.count.field character specifying alt.count meta data field in input het.sites (default $alt)
    #' @param ref.count.field character specifying ref.count meta data field in input het.sites (default $ref)
    #' @param split.ab logical flag whether to split aberrant segmetns (segmentss with ab edge entering or leaving prior to computing allelic states (default FALSE)
    #' @param uncoupled logical flag whether to not collapse segments after inferring MLE estimate (default FALSE), if FALSE will try to merge adjacent segments and populate allele-specific junctions with copy numbers on the basis of the MLE fit on individual allelic segments
    #' @param conservative if TRUE then will leave certain allelic segments "unphased" if one cannot sync the high / low interval state with the incoming and / or outgoing junction state
    #' @param marginal fix marginal? default TRUE
    #' @param verbose logical flag
    #' @return
    #' list with following fields:
    #' $segstats = GRanges of input segments with $cn.high and $cn.low segments populated
    #' $asegstats = GRanges of allelic segments (length is 2*length(segstats)) with high and low segments each having $cn, this is a "melted" segstats GRAnges
    #' $agtrack = gTrack of allelic segments and supporting input het.sites
    #' $aadj = allelic adjacency matrix of allele specific junctions
    #' $ab.ix = indices of aberrant edges in $aadj
    #' $ref.ix = indices of reference edges in $aadj
    ############################################
    jabba.alleles2 = function(jab,
                              het.sites, ## granges with meta data fields for alt.count and
                              alt.count.field = 'alt',
                              ref.count.field = 'ref',
                              baf.field = 'baf.t',
                              split.ab = F, ## if split.ab == T, then will split across any "aberrant" segment (i.e. segment with ab edge entering or leaving prior to computing allelic states (note: this might create gaps)
                              uncoupled = FALSE, ## if uncoupled, we just assign each high low allele the MLE conditioning on the total copy number
                              conservative = FALSE, ## if TRUE then will leave certain allelic segments "unphased" if one cannot sync the high / low interval state with the incoming and / or outgoing junction state
                              marginal = TRUE,
                              verbose = F
                              )
    {
        if (!all(c(alt.count.field, ref.count.field) %in% names(values(het.sites)))){
            jwarning('count fields not found in meta data of het.sites input, trying BAF...')
            if (!(baf.field %in% names(values(het.sites))))
            jerror('BAF field not found in meta data of het.sites input either!')
            else{
                ## outputs are re.seg$low and re.seg$high
                ## test deviations of observed BAF from expected by beta distribution
                if (verbose)
                message('Processing', length(het.sites),
                    'het sites using fields', baf.field, '\n')

            }
        } else {
            ## jerror('count fields not found in meta data of het.sites input')

            if (verbose)
            {
                message('Processing ', length(het.sites), ' het sites using fields ', alt.count.field, ' and ', ref.count.field)
            }

            het.sites$low.count = pmin(values(het.sites)[, alt.count.field], values(het.sites)[, ref.count.field])
            het.sites$high.count = pmax(values(het.sites)[, alt.count.field], values(het.sites)[, ref.count.field])

            het.sites = het.sites[!is.na(het.sites$low.count) & !is.na(het.sites$high.count)]

            ## stretch out het sites
            bin.gaps = gaps(het.sites)
            bin.gaps = bin.gaps %Q% (strand(bin.gaps) == "*")
            bin.gaps = resize(bin.gaps, width = width(bin.gaps) + 1, fix = "start")
            values(bin.gaps)[, "low.count"] = gr.val(query = bin.gaps[, c()],
                target = het.sites,
                val = "low.count",
                mean = TRUE,
                na.rm = TRUE)$low.count
            values(bin.gaps)[, "high.count"] = gr.val(query = bin.gaps[, c()],
                target = het.sites,
                val = "high.count",
                mean = TRUE,
                na.rm = TRUE)$high.count

            het.sites = bin.gaps ## use these stretched out sites

            het.sites = het.sites[!is.na(het.sites$low.count) & !is.na(het.sites$high.count)]

            ss.p = jab$segstats[ as.logical( strand(jab$segstats)=='+' ) ]

            ## find the reference junctions
            ord.ix = order(jab$segstats)
            rev.ix = as.logical(strand(jab$segstats[ord.ix]) == '-')
            ord.ix = c(ord.ix[!rev.ix], rev(ord.ix[rev.ix]))

            ref.jun = cbind(ord.ix[-length(ord.ix)], ord.ix[-1])
            ref.jun = ref.jun[which(jab$adj[ref.jun]>0), ]

            has.ab.rand = 0
            if (split.ab)
            {
                ab.adj = jab$adj
                ab.adj[ref.jun] = 0
                has.ab = as.numeric(Matrix::rowSums(ab.adj!=0)!=0 | Matrix::colSums(ab.adj!=0)!=0)[which( as.logical( strand(jab$segstats)=='+')) ]
                has.ab.rand = runif(length(ss.p)) * 1e-6 * has.ab
            }

            ss.p = ss.p[!is.na(ss.p$cn)]
            ## browser()
            #' zchoo Wednesday, Apr 27, 2022 08:34:21 AM
            ## don't do this since this merges ranges with the same score!
            ## re.seg = as(coverage(ss.p, weight = ss.p$cn + has.ab.rand), 'GRanges')
            re.seg = ss.p[, "cn"]
            ## re.seg$cn = round(re.seg$score)

            het.sites$ix = gr.match(het.sites, re.seg)

            if (verbose)
            {
                message('Computed high / low counts and matched to segs')
            }


            highs = split(het.sites$high.count, het.sites$ix)[as.character(seq_along(re.seg))]
            lows = split(het.sites$low.count, het.sites$ix)[as.character(seq_along(re.seg))]

            het.sites$cn = re.seg$cn[het.sites$ix]
            purity = jab$purity
            ploidy = mean(het.sites$cn, na.rm = T) ## ploidy may be slightly different from "global ploidy" depending on the distribution of sites

            sw = length(het.sites)
            total = sum(as.numeric(c(het.sites$high.count, het.sites$low.count)))

            cn = re.seg$cn
            ## gamma = 2*(1-purity)/purity  ## gammas and betas need to be recomputed for
            ## beta = (2*(1-purity)*sw + purity*ploidy*sw) / (purity * total)
            gamma = 1*(1-purity)/purity  ## gammas and betas need to be recomputed for  (1 since we are looking at het alleles)
            beta = (1*(1-purity)*sw + purity*ploidy*sw) / (purity * total)
            centers = (0:(max(cn)+1) + gamma)/beta

            if (verbose)
            {
                message('Computed SNP ploidy and allelic copy centers')
            }

            ## now test deviation from each absolute copy combo using poisson model
            ## i.e. counts ~ poisson(expected mean)
            ##
            re.seg.tmp = lapply(seq_along(re.seg), function(i)
                {
                    ##        if (verbose)
                    ##          cat('.')
                    x = lows[[i]]
                    if (length(x)==0)
                    return(list(low = NA, high = NA))
                    y = highs[[i]]
                    tot.cn = cn[i]
                    ## allow for +/- errors
                    ## use negative binomial??
                    ## ll = sapply(0:(floor(tot.cn/2)), function(j) sum(pnbinom(x, mu = centers[j+1],
                    ##                                                          size = 0,##centers[j+1] / 2,
                    ##                                                          log.p = T) +
                    ##                                                  pnbinom(y, mu = centers[tot.cn - j + 1],
                    ##                                                          size = centers[tot.cn - j + 1],
                    ##                                                          log.p = T)))
                    ll = sapply(0:(floor(tot.cn/2)), function(j) sum(ppois(x,centers[j+1], log.p = T) +
                        ppois(y,centers[tot.cn-j+1],log.p = T)))
                    ll = ll - min(ll)
                    curr.best = max(ll)
                    curr.cn = which.max(ll) - 1
                    curr.hcn = tot.cn - curr.cn
                    if (!marginal) {
                        if (tot.cn > 1) {
                            ## ll = sapply(0:(floor(tot.cn/2)), function(j) sum(pnbinom(x, mu = centers[j+1],
                            ##                                                          size = centers[j+1] / 2,
                            ##                                                          log.p = T) +
                            ##                                                  pnbinom(y, mu = centers[tot.cn - j],
                            ##                                                          size = centers[tot.cn - j] / 2,
                            ##                                                          log.p = T)))
                            ll = sapply(0:(floor((tot.cn - 1)/2)), function(j) sum(ppois(x,centers[j+1], log.p = T) +
                                ppois(y,centers[tot.cn-j],log.p = T)))
                            ll = ll - min(ll)
                            if (max(ll) > curr.best) {
                                curr.best = max(ll)
                                curr.cn = which.max(ll) - 1
                                curr.hcn = tot.cn - curr.cn - 1
                            }
                        }
                        ## ll = sapply(0:(floor(tot.cn/2)), function(j) sum(pnbinom(x, mu = centers[j+1],
                        ##                                                          size = centers[j+1],
                        ##                                                          log.p = T) +
                        ##                                                  pnbinom(y, mu = centers[tot.cn - j + 2],
                        ##                                                          size = centers[tot.cn - j + 2],
                        ##                                                          log.p = T)))
                        ll = sapply(0:(floor((tot.cn + 1)/2)), function(j) sum(ppois(x,centers[j+1], log.p = T) +
                            ppois(y,centers[tot.cn-j+2],log.p = T)))
                        ll = ll - min(ll)
                        if (max(ll) > curr.best) {
                            curr.best = max(ll)
                            curr.cn = which.max(ll) - 1
                            curr.hcn = tot.cn - curr.cn + 1
                        }
                    }
                    ##return(list(low = curr.cn, high = tot.cn - curr.cn))
                    return(list(low = curr.cn, high = curr.hcn))
                })

            re.seg$low = lapply(re.seg.tmp, function(x) {x$low}) %>% unlist
            re.seg$high = lapply(re.seg.tmp, function(x) {x$high}) %>% unlist
        }
        ## #########################################################################
        ## borderline, below are common to both methods
        ## no need to round and NA segments are fine
        jab$segstats$cn.low = gr.val(jab$segstats, re.seg, 'low', na.rm = TRUE)$low
        jab$segstats$cn.high = gr.val(jab$segstats, re.seg, 'high', na.rm = TRUE)$high

        ## also get number of hets per segment
        jab$segstats$nhets = jab$segstats %N% het.sites

        ## and the variance
        jab$segstats$var.low = gr.val(query = jab$segstats,
            target = het.sites,
            val = 'low.count',
            na.rm = TRUE,
            FUN = function(x, w, na.rm) {var(x, na.rm = na.rm)})$low.count

        jab$segstats$var.high = gr.val(query = jab$segstats,
            target = het.sites,
            val = 'high.count',
            na.rm = TRUE,
            FUN = function(x, w, na.rm) {var(x, na.rm = na.rm)})$high.count

        ## NA-out some nodes
        na.ix = (!gr.val(jab$segstats, re.seg, 'low', FUN = function(x,w,na.rm) any(!is.na(x)))$low) |
        (!gr.val(jab$segstats, re.seg, 'high', FUN = function(x,w,na.rm) any(!is.na(x)))$high)
        jab$segstats$cn.low[na.ix] = jab$segstats$cn.high[na.ix] = NA
        jab$segstats$var.low[na.ix] = jab$segstats$var.high[na.ix] = NA
        jab$segstats$nhets[na.ix] = 0

        ## ## ###########
        ## ## phasing
        ## ## ###########

        ## ## iterate through all reference junctions and apply (wishful thinking) heuristic
        ## ##
        ## ## populate n x n x 2 adjacency matrix, which we will later expand to a bigger matrix
        ## adj.ab = jab$adj
        ## adj.ab[ref.jun] = 0
        ## adj.ref = jab$adj*0
        ## adj.ref[ref.jun] = jab$adj[ref.jun]
        ## high = low = jab$segstats[, c()]
        ## high$cn = jab$segstats$cn.high
        ## low$cn = jab$segstats$cn.low
        ## high$parent = low$parent = seq_along(jab$segstats)
        ## high$type = 'high'
        ## low$type = 'low'
        ## high$id = seq_along(jab$segstats)
        ## low$id = length(jab$segstats) + seq_along(jab$segstats)
        ## asegstats = c(high, low)
        ## amap = cbind(high$id, low$id) ## maps segstats id x allele combos to asegstats id

        ## aadj = sparseMatrix(1, 1, x = 0, dims = c(length(asegstats), length(asegstats)))

        ## .flip = function(x) x %% 2+1

        ## asegstats = c(high, low)
        ## acn = cbind(high$cn, low$cn)

        ## phased.out = phased.in = rep(TRUE, length(asegstats))

        ## str = strand(asegstats)

        ## if (verbose)
        ##     message('Starting phasing ')

        ## for (k in 1:nrow(ref.jun))
        ## {
        ##     i = ref.jun[k, 1]
        ##     j = ref.jun[k, 2]
        ##     a = acn[ref.jun[k,1],]
        ##     b = acn[ref.jun[k,2],]

        ##     phased.out[amap[i, ]] = FALSE
        ##     phased.in[amap[j, ]] = FALSE

        ##     pairs.ij = cbind(rep(c(1:2), 2), rep(c(1:2), each = 2)) ## 4 possible matches
        ##     m = setdiff(which(a[pairs.ij[,1]] == b[pairs.ij[,2]]), NA)

        ##     if (!(length(m) %in% c(0, 4))) ## 1,2, and 3 matches are fine (3 matches occur if one interval is in allelic balance, and the other not
        ##     {
        ##         if (length(m)==2) ## pick the phase that the alleles can handle
        ##             m = rev(m[order(as.numeric(sum(adj.ab[i, ])<=a[pairs.ij[m,1]]) + as.numeric(sum(adj.ab[, j])<=b[pairs.ij[m,2]]))])

        ##         m.ij = pairs.ij[m[1], ]
        ##         fm.ij = .flip(m.ij)
        ##         aadj[amap[i, m.ij[1]], amap[j, m.ij[2]]] = min(a[m.ij[1]], jab$adj[i, j])
        ##         aadj[amap[i, fm.ij[1]], amap[j, fm.ij[2]]] = jab$adj[i, j] - aadj[amap[i, m.ij[1]], amap[j, m.ij[2]]]

        ##         phased.out[amap[i, ]] = TRUE
        ##         phased.in[amap[j, ]] = TRUE

        ##         if (length(a.ab <- Matrix::which(adj.ab[i,]!=0))>0)
        ##         {
        ##             ## if a.ab (partner) is already phased then unpopulate the non-ab allelic junction, otherwise populate both alleles of partner
        ##             ## BUG: a.ab is length 2????
        ##             ## hack: replace a.ab with a.ab[1]
        ##             if (any(ph <- aadj[amap[i, fm.ij[1]], amap[a.ab[1], ]] !=0))
        ##             {
        ##                 aadj[amap[i, fm.ij[1]], amap[a.ab[1], ph]] = adj.ab[i, a.ab[1]]
        ##                 aadj[amap[i, m.ij[1]], amap[a.ab[1], ph]] = 0
        ##             }
        ##             else
        ##                 ## otherwise diffuse copy into both alleles of the partner (will be resolved when we resolve phase for the partner interval)
        ##                 ## or collapse unphased nodes back
        ##                 aadj[amap[i, fm.ij[1]], amap[a.ab[1], ]] = adj.ab[i, a.ab[1]]/2

        ##             if (!conservative)
        ##                 if (a[fm.ij[1]] < adj.ab[i, a.ab]) # if the allelic node can't handle the outgoing allelic edge flux, so unphase
        ##                     phased.out[amap[i, ]] = FALSE
        ##         }

        ##         if (length(b.ab <- Matrix::which(adj.ab[,j]!=0))>0)
        ##         {
        ##             ## if b.ab (partner) is already phased then concentrate all of the junction copy into the aberrant allele of this interval
        ##             ## BUG: why b.ab is length 2???? I thought we resolved this long ago
        ##             ## hack: replace a.ab with a.ab[1]
        ##             if (any(ph <- aadj[amap[b.ab[1], ], amap[j, fm.ij[2]]] !=0))
        ##             {
        ##                 aadj[amap[b.ab[1], ph], amap[j, fm.ij[2]]] = adj.ab[b.ab[1], j]
        ##                 aadj[amap[b.ab[1], ph], amap[j, m.ij[2]]] = 0
        ##             }
        ##             else
        ##                 ## otherwise diffuse copy into both alleles of the partner (will be resolved when we resolve phase for the partner interval)
        ##                 ## or collapse unphased nodes back
        ##                 aadj[amap[b.ab[1],], amap[j, fm.ij[2]]] = adj.ab[b.ab[1], j]/2

        ##             if (!conservative)
        ##                 if (b[fm.ij[2]] < adj.ab[b.ab, j]) # the allelic node cn can't handle the incoming allelic edge flux, so unphase
        ##                     phased.in[amap[j, ]] = FALSE
        ##         }
        ##     }
        ## }

        ## if (verbose)
        ##     message('Finished phasing, finalizing.')

        ## asegstats$phased.in = phased.in
        ## asegstats$phased.out = phased.out

        ## if (uncoupled)
        ##     unphased = rep(FALSE, length(asegstats))
        ## else
        ##     unphased = !asegstats$phased.out | !asegstats$phased.in

        ## unphased.parents = unique(asegstats$parent[unphased])
        ## aadj.unphunph = jab$adj[unphased.parents, unphased.parents]
        ## aadj.phph = aadj[!unphased, !unphased]

        ## asegstats$new.ind = NA
        ## asegstats$new.ind[!unphased] = 1:sum(!unphased)
        ## asegstats$new.ind[unphased] = as.integer(factor(asegstats$parent[unphased], unphased.parents))
        ## mat.collapse = sparseMatrix(which(unphased), asegstats$new.ind[unphased], x = 1, dims = c(nrow(aadj), length(unphased.parents)))

        ## aadj.phunph = aadj[!unphased, ] %*% mat.collapse
        ## aadj.unphph = t(mat.collapse) %*% aadj[, which(!unphased)]

        ## aadj.final = rbind(
        ##     cbind(aadj.phph, aadj.phunph),
        ##     cbind(aadj.unphph, aadj.unphunph)
        ## )

        ## asegstats.unphased = asegstats[match(unphased.parents, asegstats$parent)]
        ## asegstats.unphased$cn = jab$segstats$cn[asegstats.unphased$parent]
        ## asegstats.final = c(asegstats[!unphased], asegstats.unphased)
        ## asegstats.final$phased = c(rep(T, sum(!unphased)), rep(F, length(asegstats.unphased)))
        ## asegstats.final$type[!asegstats.final$phased] = 'total'

        ## tmp.str = gr.string(gr.stripstrand(asegstats), mb = F, other.cols = 'type');
        ## asegstats$tile.id = as.integer(factor(tmp.str, unique(tmp.str)))

        ## ix = order(asegstats.final)
        ## asegstats.final = asegstats.final[ix]
        ## aadj.final = aadj.final[ix, ix]

        ## if (verbose)
        ##     message('Annotating allelic vertices')

        ## tmp.string = gr.string(asegstats, mb = F, other.cols = 'type'); tmp.string2 = gr.string(gr.flipstrand(asegstats), mb = F, other.cols = 'type')
        ## asegstats$flip.ix = match(tmp.string, tmp.string2)
        ## asegstats$phased = !unphased

        ## asegstats.final$edges.in = sapply(seq_along(asegstats.final),
        ##                                   function(x) {ix = Matrix::which(aadj.final[,x]!=0); paste(ix, '(', aadj.final[ix,x], ')', '->', sep = '', collapse = ',')})
        ## asegstats.final$edges.out = sapply(seq_along(asegstats.final),
        ##                                    function(x) {ix = Matrix::which(aadj.final[x, ]!=0); paste('->', ix, '(', aadj.final[x,ix], ')', sep = '', collapse = ',')})

        ## asegstats.final$slack.in = asegstats.final$cn - Matrix::colSums(aadj.final)
        ## asegstats.final$slack.out = asegstats.final$cn - Matrix::rowSums(aadj.final)

        ## asegstats.final$new.ind = asegstats.final$phased.out = asegstats.final$phased.in = asegstats.final$id = NULL
        ## asegstats.final$tile.id = as.integer(factor(gr.string(gr.stripstrand(asegstats.final), mb = F, other.cols = 'type')))

        ## m = sparseMatrix(seq_along(asegstats.final), asegstats.final$parent, x = 1);

        ## hh = rep(het.sites[, c()], 2)
        ## hh$count = c(het.sites$high.count, het.sites$low.count)
        ## hh$type = rep(c('high', 'low'), each = length(het.sites))

        ## hh$ywid = 0.5
        ## atd = c(
        ##     gTrack(hh, angle = 0, y.field = 'count', y0 = 0,
        ##            colormaps = list(type = c('high' = alpha('red', 0.3), 'low' = alpha('blue', 0.3))), name = 'hets', y.quantile = 0.001, lwd.border = 2),
        ##     gTrack(asegstats.final, angle = 0, y.field = 'cn', y0 = 0,
        ##            colormaps = list(type = c('high' = alpha('red', 0.3), 'low' = alpha('blue', 0.3), 'total' = alpha('purple', 0.3))), name = 'alleles')
        ## )

        ## out = list(
        ##     segstats = jab$segstats,
        ##     asegstats = asegstats.final,
        ##     atd = atd,
        ##     agtrack = atd,
        ##     aadj = aadj.final,
        ##     ab.ix = Matrix::which((m %*% adj.ab %*% t(m))!=0, arr.ind = T),
        ##     ref.ix = Matrix::which((m %*% adj.ref %*% t(m))!=0, arr.ind = T)
        ## )

        return(jab)
    }

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

        ## return(dt2gr(out.dt))
        return(dt)
    }

    ## end of utils.R

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

    # handle case where nodes.ov.hets is empty
    if (nrow(nodes.ov.hets) == 0) {
        message("No nodes overlap with hets")
        segs = NULL
    } else {
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
    }

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
        gg$meta$ploidy = ifelse(is.null(gg$meta$ploidy), NA, gg$meta$ploidy)
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
        print('tracer2')
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

    # override gGnome balance with zc_dev gGnome balance function)

    #' @name balance
    #' @title balance gGnome graphs
    #' @description
    #'
    #' Here we analyze gGraphs with "cn" (copy number) field to enforce integer
    #' cn and junction balance, ie sum of incoming (or outgoing) edge
    #' cn should be equal to node copy cn.
    #'
    #' The goal is to find a balaned assignment of "cn" to the nodes and edges of the gGraph
    #' that maximally resemble the input weights while minimizing the loose end penalty.
    #' The similarity / distance function can be weighted by optional node / edge
    #' metadata field $weight (when weighted = TRUE).
    #'
    #' To output this gGraph, we design a MIP with
    #' (1) objective function that minimizes (weighted sum of square) distance of fit node and junction copy number to provided values in
    #'     $cn field
    #' (2) and lambda* the sum of copy number at non terminal loose ends subject to
    #' (3) junction balance constraint
    #' (4) fixing copy number of a subset of nodes and junctions
    #'
    #' Objective weight can be modulated at nodes and edges with $weight metadata
    #' field (default node weight is node width, and edge weight is 1).
    #' These fields will then set the penalty incurred to a fit of x to that node / edge
    #' with copy number c and weight w as (x-c)^2/w.
    #'
    #' Lambda can be modulated at nodes with $lambda node metadata field (default 1)
    #'
    #' For "haplographs" ie graphs that have more than one node overlapping a given location, it may
    #' be important to constrain total copy number using a haploid read depth signal.
    #' The marginal parameter enables this through a GRanges annotated with $cn and optionally $weight
    #' field that provides a target total copy number that the optimization will attempt to satisfy.
    #' This provided copy number c and weight w (default 1) will be evaluated against the
    #' sum s of the fit copy numbers of all nodes overlapping that location by adding a penalty
    #' of (c-s)^2/w to the corresponding solution. marginal can also have an optional logical field
    #' $fix that will actually constrain the marginal copy number to be equal to the provided value
    #' (note: that the optimization may be infeasible, and function will error out)
    #'
    #' Additional controls can be inputted by changing the graph metadata - e.g. adding fields
    #' $lb and $ub to nodes and edges will constrain their fit copy number to those bounds.
    #' Adding $reward field to edges will add a reward for each copy of that edge in the solution.
    #'
    #'
    #' @param gg gGraph with field $cn, can be NA for some nodes and edges, optional field $weight which will adjust the quadratic penalty on the fit to x as (x-$cn)^2/weight
    #' @param lambda positive number specifying loose end penalty, note if gg$node metadata contain $lambda field then this lambda will be multiplied by the node level lambda (default 10)
    #' @param marginal GRanges with field $cn and optional $weight field will be used to fit the summed values at each base of the genome to optimally fit the marginal value, optional field $fix will actually constrain the marginal to be the provided value
    #' @param emarginal Junctions object with marginal CN in the $cn field (and optionally $weight in the weight field). optional field $fix will actually constrain the marginal to be the provided value.
    #' @param tight indices or epxression on node metadata specifying at which nodes to disallow loose ensd
    #' @param nfix indices or expression on node metadata specifying which node cn to fix
    #' @param efix indices or expression on edge metadata specifying which edge cn to fix
    #' @param nrelax indices or expression on node metadata specifying which nodes cn to relax
    #' @param erelax  indices or expression on edge metadata specifying which edges cn to relax
    #' @param L0  flag whether to apply loose end penalty as L1 (TRUE)
    #' @param loose.collapse (parameter only relevant if L0 = TRUE) will count all unique (by coordinate) instances of loose ends in the graph as the loose end penalty, rather than each instance alone ... useful for fitting a metagenome graph   (FALSE)
    #' @param phased (logical) indicates whether to run phased/unphased. default = FALSE
    #' @param ism  (logical) additional ISM constraints (FALSE)
    #' @param force.major (logical) force major allele CN to be >= minor allele CN (default FALSE)
    #' @param force.alt (logical) force incorporation of ALT edges, only applicable for phasing (default TRUE)
    #' @param cnloh (logical) allow CN LOH? only relevant if phasing = TRUE. default FALSE.
    #' @param lp (logical) solve as linear program using abs value (default TRUE)
    #' @param M  (numeric) big M constraint for L0 norm loose end penalty (default 1e3)
    #' @param verbose (integer)scalar specifying whether to do verbose output, value 2 will spit out MIP (1)
    #' @param tilim (numeric) time limit on MIP in seconds (10)
    #' @param epgap (numeric) relative optimality gap threshhold between 0 and 1 (default 1e-3)

    #' @param trelim (numeric) max size of uncompressed tree in MB (default 32e3)
    #' @param nodefileind (numeric) one of 0 (no node file) 1 (in memory compressed) 2 (on disk uncompressed) 3 (on disk compressed) default 1
    #' @param debug (logical) returns list with names gg and sol. sol contains full RCPLEX solution. (default FALSE)
    #' @param force.haplotypes (logical) default TRUE
    #' @param max.span (numeric) the maximum span of an edge below which both endpoints must be on the same parental haplotype default 1e9
    #' @param use.gurobi (logical) use gurobi optimizer? if TRUE uses gurobi instead of cplex. default FALSE.
    #' @param nonintegral (logical) run without integer constraints on REF edges and nodes? default FALSE.
    #'
    #' @return balanced gGraph maximally resembling input gg in CN while minimizing loose end penalty lambda.
    #' @author Marcin Imielinski
    #'
    #' @export

    balance = function(gg,
                       lambda = 10,
                       marginal = NULL,
                       emarginal = NULL,
                       tight = NULL,
                       nfix = NULL, efix = NULL, nrelax = NULL, erelax = NULL,
                       L0 = TRUE,
                       loose.collapse = FALSE,
                       M = 1e3,
                       phased = FALSE,
                       ism = TRUE,
                       force.major = TRUE,
                       force.alt = FALSE,
                       cnloh = FALSE,
                       lp = TRUE,
                       verbose = 1,
                       tilim = 10,
                       trelim = 32e3,
                       nodefileind = 1,
                       epgap = 1e-3,
                       max.span = 1e6, ## max span in bp
                       debug = FALSE,
                       use.gurobi = FALSE,
                       nonintegral = FALSE)
    {
        if (verbose) {
            message("creating copy of input gGraph")
        }

        if (use.gurobi) {
            if (!requireNamespace("gurobi", quietly = TRUE)) {
                stop("use.gurobi is TRUE but gurobi is not installed")
            }
        }

        gg = gg$copy

        if (verbose) {
            message("Checking inputs")
        }

        if (nodefileind) {
            if (!(nodefileind %in% c(0,1,2,3))) {
                warning("Invalid choice for nodefileind, resetting to default 1")
                nodefileind = 1
            }
        }
        nodefileind = as.integer(nodefileind)

        if (ism) {
            if (!L0) {
                stop("ISM can only be set to true if using L0 penalty")
            }
        }

        if (!('cn' %in% names(gg$nodes$dt)))
        {
            warning('cn field not defined on nodes, setting to NA')
            gg$nodes$mark(cn = NA_real_)
        }

        if (!('cn' %in% names(gg$edges$dt)))
        {
            warning('cn not defined on edges, providing NA')
            gg$edges$mark(cn = NA_real_)
        }

        if (phased) {
            if (!("allele" %in% names(gg$nodes$dt))) {
                stop("cannot run phased balance without $allele field in nodes")
            }
        }

        if (!is.null(marginal)) {
            if (!inherits(marginal, 'GRanges') || is.null(marginal$cn)) {
                stop('marginal must be a GRanges with field $cn')
            }
            if (is.null(marginal$fix)) {
                if (verbose) {
                    message("$fix not supplied. marginals not fixed by default.")
                }
                marginal$fix = 0
            }
            if (is.null(marginal$weight)) {
                if (verbose) {
                    message("$weight not supplied. set to range width in Mbp by default.")
                }
                marginal$weight = width(marginal)
            }
        }

        if (!is.null(emarginal)) {
            if (!inherits(emarginal, 'Junction') || is.null(emarginal$dt$cn)) {
                stop('emarginal must be Junction with field $cn')
            }
            ## don't mutate?
            ## emarginal = emarginal$copy
            if (is.null(emarginal$dt$fix)) {
                if (verbose) {
                    message('$fix not supplied in emarginal. not fixed by default')
                }
                emarginal$set(fix = 0)
            }
            if (is.null(emarginal$dt$weight)) {
                if (verbose) {
                    message("$weight not supplied in emarginal. set to 1 by default")
                }
                emarginal$set(weight = 1)
            }
        }

        ## default local lambda: default local lambda is 1 for consistency with JaBbA
        if (!('lambda' %in% names(gg$nodes$dt)))
            gg$nodes$mark(lambda = 1)


        ## default node weight is its width
        if (!('weight' %in% names(gg$nodes$dt)))
        {
            gg$nodes$mark(weight = width(gg$nodes$gr))
        }

        ## default edge weight is its width
        if (!('weight' %in% names(gg$edges$dt)))
        {
            gg$edges$mark(weight = 1)
        }

        ## default reward is 0
        if (!('reward' %in% names(gg$edges$dt)))
        {
            gg$edges$mark(reward = 0)
        }

        ## handle parsing of efix, nfix, nrelax, erelax
        if (!any(deparse(substitute(nfix)) == "NULL")) ## R voodo to allow "with" style evaluation
            nfix = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(nfix)))), parent.frame()), gg$nodes$dt, parent.frame(2)), error = function(e) NULL)

        if (!any(deparse(substitute(nrelax)) == "NULL")) ## R voodo to allow "with" style evaluation
            nrelax = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(nrelax)))), parent.frame()), gg$nodes$dt, parent.frame(2)), error = function(e) NULL)

        if (!any(deparse(substitute(efix)) == "NULL")) ## R voodo to allow "with" style evaluation
            efix = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(efix)))), parent.frame()), gg$edges$dt, parent.frame(2)), error = function(e) NULL)

        if (!any(deparse(substitute(erelax)) == "NULL")) ## R voodo to allow "with" style evaluation
            erelax = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(erelax)))), parent.frame()), gg$edges$dt, parent.frame(2)), error = function(e) NULL)

        if (!any(deparse(substitute(tight)) == "NULL")) ## R voodo to allow "with" style evaluation
            tight = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(tight)))), parent.frame()), gg$nodes$dt, parent.frame(2)), error = function(e) NULL)


        if (is.logical(nfix))
            nfix = which(nfix)

        if (is.logical(efix))
            efix = which(efix)

        if (is.logical(nrelax))
            nrelax = which(nrelax)

        if (is.logical(erelax))
            erelax = which(erelax)

        if (length(nfix) & verbose)
            message('Fixing ', length(nfix), ' nodes')

        if (length(efix) & verbose)
            message('Fixing ', length(efix), ' edges')

        if (length(nrelax) & verbose)
            message('Relaxing ', length(nrelax), ' nodes')

        gg$nodes[nrelax]$mark(weight = 0)

        if (length(erelax) & verbose)
            message('Relaxing ', length(erelax), ' edges')

        gg$nodes[erelax]$mark(weight = 0)

        if (!is.logical(tight))
            tight = 1:length(gg$nodes) %in% tight

        if (any(tight) & verbose)
            message('Leaving ', sum(tight), ' nodes tight')

        gg$nodes$mark(tight = tight)

        if (is.null(gg$nodes$dt$lb))
            gg$nodes$mark(lb = 0)

        if (is.null(gg$nodes$dt$ub))
            gg$nodes$mark(ub = Inf)

        if (is.null(gg$edges$dt$lb))
            gg$edges$mark(lb = 0)

        if (is.null(gg$edges$dt$ub))
            gg$edges$mark(ub = Inf)

        if (loose.collapse)
        {
            if (verbose)
                message('Collapsing loose ends')

            uleft = unique(gr.start(gg$nodes$gr))
            uright = unique(gr.end(gg$nodes$gr))

            gg$nodes$mark(loose.left.id = paste0(gr.match(gr.start(gg$nodes$gr), uleft), 'l'))
            gg$nodes$mark(loose.right.id = paste0(gr.match(gr.end(gg$nodes$gr), uright), 'r'))
        }
        else
        {
            gg$nodes$mark(loose.left.id = paste0(1:length(gg$nodes), 'l'))
            gg$nodes$mark(loose.right.id = paste0(1:length(gg$nodes), 'r'))
        }

    ########
        ## VARIABLES
    ########

        ## create state space, keeping track of graph ids
        vars = rbind(
            gg$dt[, .(cn, snode.id, lb, ub, weight, gid = index, type = 'node', vtype = 'I')],
            gg$sedgesdt[, .(from, to, lb, ub, sedge.id,  cn, reward,
                            gid = sedge.id, type = 'edge', vtype = 'I',
                            span = gg$junctions$span[match(edge.id, gg$junctions$dt[, edge.id])])],
            ## for loose ends lid marks all "unique" loose ends (which if loose.collapse = TRUE
            ## will be defined on the basis of coordinate overlap)
            gg$dt[tight == FALSE, .(cn = NA, snode.id, lambda, gid = index,
                                    ulid = paste0(index, 'i'),
                                    lid = ifelse(strand == '+', loose.left.id, paste0('-', loose.right.id)),
                                    type = 'loose.in', vtype = 'I')], ## incoming loose ends
            gg$dt[tight == FALSE, .(cn = NA, snode.id, lambda, gid = index,
                                    ulid = paste0(index, 'o'),
                                    lid = ifelse(strand == '+', loose.right.id, paste0('-', loose.left.id)),
                                    type = 'loose.out', vtype = 'I')], ## outgoing loose ends
            gg$dt[tight == FALSE, .(gid = index, cn,
                                    weight, type = 'nresidual',
                                    vtype = 'C')], ## node residual
            gg$sedgesdt[, .(gid = sedge.id, cn,
                            weight, type = 'eresidual',
                            vtype = 'C')], ## edge residual
            fill = TRUE)


        ## add "slush" variables - there will be one per chromosome
        if (nonintegral) {

            if (verbose) { message("Adding slush variables") }
            ## grab standard chromosomes from gGraph
            chr.names = grep("^(chr)*[0-9XY]+$", as.character(seqlevels(gg)), value = TRUE)
            slush.dt = data.table(chr = chr.names,
                                  lb = -0.49999,
                                  ub = 0.49999,
                                  type = "slush",
                                  vtype = "C",
                                  gid = 1:length(chr.names))

            vars = rbind(vars, slush.dt, fill = TRUE)

            vars[type == "node", chr := as.character(gg$dt$seqnames)[match(snode.id, gg$dt$snode.id)]]
            vars[!(chr %in% slush.dt[, chr]), chr := NA_character_]

        }

        if (L0)
        {
            ## loose ends are labeled with lid and ulid, lid is only relevant if loose.collapse is true
            ## (i.e. we need indicator.sum and indicator.sum.indicator
            if (verbose) {
                message("adding l0 penalty indicator")
            }

            vars = rbind(vars,
                         rbind(
                             vars[type == 'loose.in', ][ , type := 'loose.in.indicator'][, vtype := 'B'][, gid := lid],
                             vars[type == 'loose.out', ][ , type := 'loose.out.indicator'][, vtype := 'B'][, gid := lid]
                         ))

            if (loose.collapse)
            {
                ## sum will sum all the loose ends assocaited with the same lid
                vars = rbind(vars,
                             unique(rbind(
                                 vars[type == 'loose.in', ][ , type := 'loose.in.indicator.sum'][, vtype := 'I'][, gid := lid],
                                 vars[type == 'loose.out', ][ , type := 'loose.out.indicator.sum'][, vtype := 'I'][, gid := lid]
                             ), by = 'gid'))

                ## sum.indicator is an binary indicator on the sum
                vars = rbind(vars,
                             rbind(
                                 vars[type == 'loose.in.indicator.sum', ][ , type := 'loose.in.indicator.sum.indicator'][, vtype := 'B'][, gid := lid],
                                 vars[type == 'loose.out.indicator.sum', ][ , type := 'loose.out.indicator.sum.indicator'][, vtype := 'B'][, gid := lid]
                             ))
            }
        }

        if (!is.null(marginal)) {
            ## first disjoin marginal against the nodes
            ## ie wee ned to create a separate residual variable for every unique
            ## disjoint overlap of marginal with the nodes
            dmarginal = gg$nodes$gr %>% gr.stripstrand %*% grbind(marginal %>% gr.stripstrand) %>%
                disjoin %$% marginal[, c('cn', 'weight', 'fix')] %Q%
                (!is.na(cn)) %Q% (!is.na(weight)) %Q% (!is.infinite(weight))

            vars = rbind(vars,
                         gr2dt(dmarginal)[, .(cn, weight, mfix = fix>0,
                                              rid = 1:.N, type = 'mresidual', vtype = 'C')][, gid := rid],
                         fill = TRUE
                         )
        }


        if (!is.null(emarginal)) {
            ## we need to identify which junction in the marginal each junction in the phased graph corresponds to
            junction.map = merge.Junction(
                phased = gg$junctions[, c()],
                emarginal = emarginal[, c("cn", "weight", "fix")],
                cartesian = TRUE,
                all.x = TRUE)$dt
            ## match this back with edge id and add this to vars
            ## vars[type == "edge", emarginal.id := junction.map[abs(sedge.id), seen.by.emarginal]]
            vars[type == "edge", emarginal.id := junction.map[abs(sedge.id), subject.id]]
            ## add weight and target total CN
            emarginal = merge.data.table(unique(
                vars[type == "edge" & !is.na(emarginal.id),][, type := "emresidual"][, .(emarginal.id, sedge.id, lb = -M, ub = M, gid, type, vtype = "C", from, to)],
                by = "emarginal.id"),
                junction.map[, .(subject.id, weight, cn, fix)],
                by.x = "emarginal.id",
                by.y = "subject.id")
            vars = rbind(vars, emarginal, fill = TRUE)
        }

        if (lp) {
            ## need delta plus and delta minus for nodes and edges
            delta.node = gg$dt[tight == FALSE, .(gid = index, cn, weight, vtype = 'C')]
            delta.edge = gg$sedgesdt[, .(gid = sedge.id, cn, weight, reward, vtype = 'C')]

            deltas = rbind(
                delta.node[, .(gid, weight, vtype, type = "ndelta.plus")],
                delta.node[, .(gid, weight, vtype, type = "ndelta.minus")],
                delta.edge[, .(gid, weight, vtype, type = "edelta.plus")],
                delta.edge[, .(gid, weight, vtype, type = "edelta.minus")]
            )

            vars = rbind(
                vars,
                deltas,
                fill = TRUE
            )

            ## add deltas for marginals if marginals are supplied
            if (!is.null(marginal)) {
                mdeltas = rbind(
                    vars[type == "mresidual", .(rid, weight, vtype, type = "mdelta.plus")][, gid := rid],
                    vars[type == "mresidual", .(rid, weight, vtype, type = "mdelta.minus")][, gid := rid]
                )
                vars = rbind(vars, mdeltas, fill = TRUE)
            }

            ## add deltas for emresiduals if emarginals are supplied
            if (!is.null(emarginal)) {
                emdeltas = rbind(
                    vars[type == "emresidual", .(emarginal.id, weight, vtype, type = "emdelta.plus")][, gid := emarginal.id],
                    vars[type == "emresidual", .(emarginal.id, weight, vtype, type = "emdelta.minus")][, gid := emarginal.id]
                )
                vars = rbind(vars, emdeltas, fill = TRUE)
            }
        }

        ## moved from being ISM-specific annotation as this is more generally useful information
        ## add telomeric annotation
        qtips = gr.end(si2gr(seqlengths(gg$nodes))) ## location of q arm tips
        term.in = c(which(start(gg$nodes$gr) == 1), ## beginning of chromosome
                    -which(gg$nodes$gr %^% qtips)) ## flip side of chromosome end
        term.out = -term.in ## out is reciprocal of in

        ## annotate loose indicators with this
        vars[!is.na(snode.id), telomeric := ifelse(snode.id %in% term.in | snode.id %in% term.out,
                                                   TRUE,
                                                   FALSE)]



        if (phased) {
            ## add allele information and og.node.id
            node.match = match(vars[, snode.id], gg$dt$snode.id)
            vars[, ":="(allele = gg$dt$allele[node.match],
                        og.node.id = gg$dt$og.node.id[node.match])]




            ## add ref/alt information and og.edge.id
            edge.match = match(vars[, sedge.id], gg$sedgesdt$sedge.id)
            gg$edges$mark(span = gg$junctions$span)

            vars[, ":="(ref.or.alt = gg$sedgesdt$type[edge.match],
                        connection = gg$sedgesdt$connection[edge.match],
                        og.edge.id = gg$sedgesdt$og.edge.id[edge.match],
                        span = gg$sedgesdt$span[edge.match], ## NEED SPAN
                        n1 = gg$dt$snode.id[gg$sedgesdt$from[edge.match]],
                        n2 = gg$dt$snode.id[gg$sedgesdt$to[edge.match]])]

            vars[, n1.side := ifelse(n1 > 0, "right", "left")]
            vars[, n2.side := ifelse(n2 > 0, "left", "right")]
            vars[, n1 := abs(n1)]
            vars[, n2 := abs(n2)]

            edge.indicator.vars = vars[type == "edge"][, type := "edge.indicator"][, vtype := "B"][, gid := sedge.id]
            vars = rbind(vars, edge.indicator.vars, fill = TRUE)

            #' zchoo Thursday, Sep 02, 2021 05:30:53 PM
            #' moved to earlier
            ## REF edge configuration constraint (added by default basically)
            ## only add this if there are no unphased nodes
            if (cnloh) {

                ## if allow CNLOH, the sum of edge indicators corresponding with og edge id is LEQ 2
                ## this is only allowed in constant CN regions and if breakpoint is not shared with any ALT edges

                ## penalize CNLOH edges

                if (!is.null(gg$edges$dt$cnloh)) {
                    cnloh.edges = gg$edges$dt[cnloh == TRUE & type == "ALT", edge.id] %>% unique
                    if (verbose) {
                        message("Number of marked CNLOH edges: ", length(cnloh.edges))
                    }

                    ## add CNLOH annotation to variables
                    ## browser()
                    vars[, cnloh := FALSE]
                    vars[(type == "edge.indicator" | type == "edge" | type == "eresidual") &
                         ref.or.alt == "ALT" & (abs(sedge.id) %in% cnloh.edges) & (sedge.id > 0),
                         ":="(cnloh = TRUE)]

                } else {
                    warning("CNLOH not specified on edges. Disallowing!")
                    cnloh.og.edges = c()
                    vars[, cnloh := FALSE]
                }
            } else {

                cnloh.og.edges = c()
                vars[, cnloh := FALSE]

            }

            ## add node haplotype indicators
            ## these are binary indicators that determine whether a node belongs to H1
            ## only constrain positive-stranded nodes due to skew symmetry
            haplotype.indicators = gg$dt[(allele == "major" | allele == "minor") & snode.id > 0,
                                         .(cn, snode.id, lb, ub, weight, og.node.id,
                                           allele, gid = index, type = 'haplotype', vtype = 'B')]
            vars = rbind(vars, haplotype.indicators, fill = TRUE)

            ## add H1 and H2 'AND' indicators which should have n1/n1.side/n2/n2.side metadata
            ## only add these for low-span edges
            ## and where CNLOH is FALSE
            h1.and.indicators = vars[sedge.id > 0 & type == "edge" & (connection == "straight" | connection == "cross") & span < max.span & cnloh == FALSE,][, vtype := "B"][, type := "h1.and.indicator"][, gid := sedge.id]
            h2.and.indicators = vars[sedge.id > 0 & type == "edge" & (connection == "straight" | connection == "cross") & span < max.span & cnloh == FALSE,][, vtype := "B"][, type := "h2.and.indicator"][, gid := sedge.id]

            vars = rbind(vars, h1.and.indicators, h2.and.indicators, fill = TRUE)

        }


        ## if we want to implement edge reward or ISM, we need to add edge indicators
        ## these are binary variables that allow us to reward/penalize L0 norms
        ##if (TRUE) { ##(ism | any(gg$edges$dt$reward != 0, na.rm = TRUE)) {
        ## if not phased, must add edge indicators (for just the ALT edges)
        if (!phased) {
            edge.match = match(vars[, sedge.id], gg$sedgesdt$sedge.id)
            vars[, ":="(ref.or.alt = gg$sedgesdt$type[edge.match])] ## need ref.or.alt information
            edge.indicator.vars = vars[type == "edge" & ref.or.alt == "ALT"][, type := "edge.indicator"][, vtype := "B"][, gid := sedge.id]
            vars = rbind(vars, edge.indicator.vars, fill = TRUE)
        }

        ## loose indicators are only required if running with ISM
        if (ism) {
            vars[type == "loose.in.indicator" & sign(snode.id) == 1, ee.id := paste(snode.id, "left")]
            vars[type == "loose.out.indicator" & sign(snode.id) == 1, ee.id := paste(snode.id, "right")]
        }

        ## but even without ISM, if we want to add reward, we must add edge indicators
        vars[type == "edge.indicator" & sign(sedge.id) == 1 & ref.or.alt == "ALT",
             ":="(ee.id.n1 = paste(gg$edges$dt$n1[match(sedge.id, gg$edges$dt$sedge.id)],
                                   gg$edges$dt$n1.side[match(sedge.id, gg$edges$dt$sedge.id)]),
                  ee.id.n2 = paste(gg$edges$dt$n2[match(sedge.id, gg$edges$dt$sedge.id)],
                                   gg$edges$dt$n2.side[match(sedge.id, gg$edges$dt$sedge.id)]))]

        if (phased) {
            ## homologous extremity exclusivity (only for phased graphs)
            ## get stranded breakpoint ID's associated with the start and end of each node

            ## number of unique starts should be equal to number of snodes in the original unphased graph
            ## aka 2 * number of og edge ids

            vars[type == "loose.in.indicator", hee.id := paste(og.node.id, "in")]
            vars[type == "loose.out.indicator", hee.id := paste(og.node.id, "out")]

            vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                 ":="(og.n1 = gg$dt$og.node.id[from],
                      og.n1.side = gg$edges$dt$n1.side[match(abs(sedge.id), gg$edges$dt$edge.id)],
                      og.n2 = gg$dt$og.node.id[to],
                      og.n2.side = gg$edges$dt$n2.side[match(abs(sedge.id), gg$edges$dt$edge.id)])]

            vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                 ":="(hee.id.n1 = ifelse(og.n1.side == "left",
                                         paste(og.n1, "in"),
                                         paste(og.n1, "out")),
                      hee.id.n2 = ifelse(og.n2.side == "left",
                                         paste(og.n2, "in"),
                                         paste(og.n2, "out")))]

            ## reciprocal homologous extremity exclusivity
            ## implement config indicators. there is one per og.edge.id per configuration
            straight.config = unique(vars[type == "edge.indicator" & ref.or.alt == "REF" & sedge.id > 0, ][, type := "straight.config"][, config.id := paste("straight", og.edge.id)], by = "og.edge.id")
            cross.config = unique(vars[type == "edge.indicator" & ref.or.alt == "REF" & sedge.id > 0, ][, type := "cross.config"][, config.id := paste("cross", og.edge.id)], by = "og.edge.id")

            ## add these to vars
            vars = rbind(vars, straight.config, cross.config, fill = TRUE)

            ## add straight/cross to REF edges
            vars[type == "edge.indicator" & ref.or.alt == "REF",
                 connection := gg$sedgesdt$connection[match(sedge.id, gg$sedgesdt$sedge.id)]]

            ## add config ID's to corresponding edge indicators
            vars[type == "edge.indicator" & ref.or.alt == "REF" & sedge.id > 0,
                 config.id := paste(connection, og.edge.id)]


            ## add config ID's to corresponding edge indicators
            vars[type == "edge.indicator" & ref.or.alt == "REF" & sedge.id > 0,
                 config.id := paste(connection, og.edge.id)]

            ## add straight edge id e.g. for each n1 and n2, add the sedge.id of the corresponding straight edge
            straight.sedges = gg$edges$dt[type == "REF" & connection == "straight" & sedge.id > 0,
                                          .(n1.full = paste(n1, n1.side), n2.full = paste(n2, n2.side), sedge.id)]
            cross.sedges = gg$edges$dt[type == "REF" & connection == "cross" & sedge.id > 0,
                                       .(n1.full = paste(n1, n1.side), n2.full = paste(n2, n2.side), sedge.id)]


            ## pull alt edges from sedgesdt
            alt.sedges = gg$edges$dt[type == "ALT" & sedge.id > 0,
                                     .(n1.full = paste(n1, n1.side), n2.full = paste(n2, n2.side), sedge.id)]

            alt.sedges[, ":="(s1 = straight.sedges$sedge.id[match(n1.full, straight.sedges$n1.full)],
                              s2 = straight.sedges$sedge.id[match(n2.full, straight.sedges$n1.full)],
                              s3 = straight.sedges$sedge.id[match(n1.full, straight.sedges$n2.full)],
                              s4 = straight.sedges$sedge.id[match(n2.full, straight.sedges$n2.full)])]

            alt.sedges[, ":="(c1 = cross.sedges$sedge.id[match(n1.full, cross.sedges$n1.full)],
                              c2 = cross.sedges$sedge.id[match(n2.full, cross.sedges$n1.full)],
                              c3 = cross.sedges$sedge.id[match(n1.full, cross.sedges$n2.full)],
                              c4 = cross.sedges$sedge.id[match(n2.full, cross.sedges$n2.full)])]

            ## pull loose ends
            vars[type == "loose.in.indicator" & snode.id > 0, n2.full := paste(snode.id, "left")]
            vars[type == "loose.out.indicator" & snode.id > 0, n1.full := paste(snode.id, "right")]

            ## merge sedge.id
            vars[type == "loose.in.indicator" & snode.id > 0, ":="(s = straight.sedges$sedge.id[match(n2.full, straight.sedges$n2.full)])]
            vars[type == "loose.out.indicator" & snode.id > 0, ":="(s = straight.sedges$sedge.id[match(n1.full, straight.sedges$n1.full)])]

            vars[type == "loose.in.indicator" & snode.id > 0, ":="(c = cross.sedges$sedge.id[match(n2.full, cross.sedges$n2.full)])]
            vars[type == "loose.out.indicator" & snode.id > 0, ":="(c = cross.sedges$sedge.id[match(n1.full, cross.sedges$n1.full)])]


            ## merge this info into vars
            vars = merge(vars,
                         alt.sedges[, .(sedge.id, s1, s2, s3, s4, c1, c2, c3, c4)],
                         by = "sedge.id",
                         all.x = TRUE,
                         all.y = FALSE)

        }

        vars[, id := 1:.N] ## set id in the optimization
        vars[is.na(lb), lb := -Inf]
        vars[is.na(ub), ub := Inf]
        vars[, relax := FALSE][, fix := FALSE]
        if ("mresidual" %in% vars$type) {
            vars[type == 'mresidual' & mfix == TRUE, ":="(lb = 0, ub = 0)]
            message("Number of fixed marginals: ", nrow(vars[type == 'mresidual' & mfix == TRUE,]))
        }
        if ("emresidual" %in% vars$type) {
            vars[type == "emresidual" & fix == TRUE, ":="(lb = 0, ub = 0)]
        }

        ## redo setting lb and ub
        vars[type %in% c('node', 'edge') & is.na(lb), lb := 0]
        vars[type %in% c('node', 'edge') & is.na(ub), ub := M]
        vars[type %in% c('node', 'edge') & lb < 0, lb := 0]
        vars[type %in% c('node', 'edge') & ub > M, ub := M]
        ## vars[type %in% c('node', 'edge'), lb := ifelse(is.na(lb), 0, pmax(lb, 0, na.rm = TRUE)]
        ## vars[type %in% c('node', 'edge'), ub := ifelse(is.na(ub), M, pmin(ub, M, na.rm = TRUE))]
        vars[type %in% c('loose.in', 'loose.out'), ":="(lb = 0, ub = Inf)]

        ## reward shouldn't have to be positive
        ## vars[type %in% c('edge'), reward := pmax(reward, 0, na.rm = TRUE)]


        ## figure out junctions and nodes to fix

        vars[!is.na(cn) & type == 'node' & abs(snode.id) %in% nfix, ":="(lb = cn, ub = cn, fix = TRUE)]
        vars[!is.na(cn) & type == 'edge' & abs(sedge.id) %in% efix, ":="(lb = cn, ub = cn, fix = TRUE)]

        ## figure out terminal node sides for in and out loose ends
        ## these will not have loose ends penalized
        qtips = gr.end(si2gr(seqlengths(gg$nodes))) ## location of q arm tips
        term.in = c(which(start(gg$nodes$gr) == 1), ## beginning of chromosome
                    -which(gg$nodes$gr %^% qtips)) ## flip side of chromosome end
        term.out = -term.in
        vars$terminal = FALSE
        vars[(type %in% c('loose.in', 'loose.in.indicator')) & (snode.id %in% term.in), terminal := TRUE]
        vars[(type %in% c('loose.out', 'loose.out.indicator')) & (snode.id %in% term.out), terminal := TRUE]

        ## if not using integral constraints,
        ## change the vtype of terminal loose ends, nodes, and REF edges
        ## additionally relax their lower bound to -0.5
        if (nonintegral) {
            vars[type == "loose.in" & (snode.id %in% term.in), vtype := "C"]
            vars[type == "loose.out" & (snode.id %in% term.out), vtype := "C"]
            vars[type == "loose.in" & (snode.id %in% term.in), lb := -0.4999]
            vars[type == "loose.out" & (snode.id %in% term.out), lb := -0.4999]
            vars[type == "edge" & ref.or.alt == "REF", vtype := "C"]
            vars[type == "edge" & ref.or.alt == "REF", lb := -0.4999]
            ## vars[type == "node", lb := -0.4999]
            ## vars[type == "node", vtype := "C"]
        }

        ## fix the heaviest node
        ## browser()
        ## maxn = vars[type == "node"][which.max(weight), snode.id]
        ## vars[snode.id == maxn & type == "node", ":="(ub = (cn), lb = (cn))]

        ## browser()
        ## table(vars[, .(type, vtype)])

    ########
        ## CONSTRAINTS
        ## the key principle behind this "melted" form of constraint building is the cid
        ## (constraint id) which is the key that will group coefficients into constraints
        ## when we finally build the matrices.  So all we need to do is make sure that
        ## that value / cid pairs make sense and that every cid has an entry in b
    ########

        ## we need one junction balance constraint per loose end

        ## constraints indexed by cid
        if (nonintegral) {

            ## if running without integer constraints we have to include the slush variables
            ## for each node
            slush.sub = vars[type == "slush"]
            node.slush.in = vars[type == "node", .(value = -1,
                                                   id = slush.sub$id[match(chr, slush.sub$chr)],
                                                   cid = paste("in", gid))]
            node.slush.out = vars[type == "node", .(value = -1,
                                                   id = slush.sub$id[match(chr, slush.sub$chr)],
                                                   cid = paste("out", gid))]

            node.slush = rbind(node.slush.in, node.slush.out)[!is.na(id)]

            constraints = rbind(
                node.slush,
                vars[type == 'loose.in', .(value = 1, id, cid = paste('in', gid))],
                vars[type == 'edge', .(value = 1, id, cid = paste('in', to))],
                vars[type == 'node', .(value = -1, id, cid = paste('in', gid))],
                vars[type == 'loose.out', .(value = 1, id, cid = paste('out', gid))],
                vars[type == 'edge', .(value = 1, id, cid = paste('out', from))],
                vars[type == 'node', .(value = -1, id, cid = paste('out', gid))],
                fill = TRUE)


        } else {
            constraints = rbind(
                vars[type == 'loose.in', .(value = 1, id, cid = paste('in', gid))],
                vars[type == 'edge', .(value = 1, id, cid = paste('in', to))],
                vars[type == 'node', .(value = -1, id, cid = paste('in', gid))],
                vars[type == 'loose.out', .(value = 1, id, cid = paste('out', gid))],
                vars[type == 'edge', .(value = 1, id, cid = paste('out', from))],
                vars[type == 'node', .(value = -1, id, cid = paste('out', gid))],
                fill = TRUE)

        }

        b = rbind(
            vars[type == 'node', .(value = 0, sense = 'E', cid = paste('in', gid))],
            vars[type == 'node', .(value = 0, sense = 'E', cid = paste('out', gid))],
            fill = TRUE)

        ## add to the constraints the definitions of the node and edge
        if (nonintegral) {

            ## if running without integer constraints we have to include the slush variables
            ## for each node
            slush.sub = vars[type == "slush"]
            node.slush = vars[type == "node", .(value = 1,
                                                id = slush.sub$id[match(chr, slush.sub$chr)],
                                                cid = paste("nresidual", gid))]
            constraints = rbind(
                constraints,
                rbind(
                    node.slush,
                    vars[type == 'node', .(value = 1, id, cid = paste('nresidual', gid))],
                    vars[type == 'nresidual', .(value = -1, id, cid = paste('nresidual', gid))],
                    vars[type == 'edge', .(value = 1, id, cid = paste('eresidual', gid))],
                    vars[type == 'eresidual', .(value = -1, id, cid = paste('eresidual', gid))],
                    fill = TRUE)
            )
        } else {
            constraints = rbind(
                constraints,
                rbind(
                    vars[type == 'node', .(value = 1, id, cid = paste('nresidual', gid))],
                    vars[type == 'nresidual', .(value = -1, id, cid = paste('nresidual', gid))],
                    vars[type == 'edge', .(value = 1, id, cid = paste('eresidual', gid))],
                    vars[type == 'eresidual', .(value = -1, id, cid = paste('eresidual', gid))],
                    fill = TRUE)
            )
        }

        b = rbind(b,
                  vars[type == 'node', .(value = cn, sense = 'E', cid = paste('nresidual', gid))],
                  vars[type == 'edge', .(value = cn, sense = 'E', cid = paste('eresidual', gid))],
                  fill = TRUE)

        ## add the reverse complement equality constraints on nodes and edges
        constraints = rbind(
            constraints,
            rbind( ## +1 coefficient for positive nodes, -1 for negative nodes, matched by abs (snode.id)
                vars[type == 'node', .(value = sign(snode.id), id, cid = paste('nrc', abs(snode.id)))],
                vars[type == 'edge', .(value = sign(sedge.id), id, cid = paste('erc', abs(sedge.id)))],
                fill = TRUE)
        )

        b = rbind(b,
                  vars[type == 'node' & snode.id>0, .(value = 0, sense = 'E', cid = paste('nrc', abs(snode.id)))],
                  vars[type == 'edge' & sedge.id>0, .(value = 0, sense = 'E', cid = paste('erc', abs(sedge.id)))],
                  fill = TRUE)


        ## if solving as LP, add deltas constraints (absolute value trick)

        if (lp) {
            if (verbose) {
                message("adding delta constraints for LP")
            }

        vars[type %like% "delta.plus" | type %like% "delta.minus", ":="(ub = M, lb = 0)]

            ## add the residual constraints
            ndelta.slack = rbind(
                vars[type == "nresidual", .(value = -1, id, cid = paste("ndelta.minus.slack", gid))],
                vars[type == "ndelta.minus", .(value = -1, id, cid = paste("ndelta.minus.slack", gid))],
                vars[type == "nresidual", .(value = 1, id, cid = paste("ndelta.plus.slack", gid))],
                vars[type == "ndelta.plus", .(value = -1, id, cid = paste("ndelta.plus.slack", gid))]
            )

            ndelta.slack.rhs = rbind(
                vars[type == "ndelta.minus", .(value = 0, sense = "L", cid = paste("ndelta.minus.slack", gid))],
                vars[type == "ndelta.plus", .(value = 0, sense = "L", cid = paste("ndelta.plus.slack", gid))]
            )

            edelta.slack = rbind(
                vars[type == "eresidual", .(value = -1, id, cid = paste("edelta.minus.slack", gid))],
                vars[type == "edelta.minus", .(value = -1, id, cid = paste("edelta.minus.slack", gid))],
                vars[type == "eresidual", .(value = 1, id, cid = paste("edelta.plus.slack", gid))],
                vars[type == "edelta.plus", .(value = -1, id, cid = paste("edelta.plus.slack", gid))]
            )

            edelta.slack.rhs = rbind(
                vars[type == "edelta.minus", .(value = 0, sense = "L", cid = paste("edelta.minus.slack", gid))],
                vars[type == "edelta.plus", .(value = 0, sense = "L", cid = paste("edelta.plus.slack", gid))]
            )

            mdelta.slack = rbind(
                vars[type == "mresidual", .(value = -1, id, cid = paste("mdelta.minus.slack", gid))],
                vars[type == "mdelta.minus", .(value = -1, id, cid = paste("mdelta.minus.slack", gid))],
                vars[type == "mresidual", .(value = 1, id, cid = paste("mdelta.plus.slack", gid))],
                vars[type == "mdelta.plus", .(value = -1, id, cid = paste("mdelta.plus.slack", gid))]
            )

            mdelta.slack.rhs = rbind(
                vars[type == "mdelta.minus", .(value = 0, sense = "L", cid = paste("mdelta.minus.slack", gid))],
                vars[type == "mdelta.plus", .(value = 0, sense = "L", cid = paste("mdelta.plus.slack", gid))]
            )

            constraints = rbind(constraints, ndelta.slack, edelta.slack, mdelta.slack, fill = TRUE)
            b = rbind(b, ndelta.slack.rhs, edelta.slack.rhs, mdelta.slack.rhs, fill = TRUE)

        }

        if (phased) {

            ## add haplotype indicator constraints
            ## e.g. the haplotype indicators corresponding to the same og node must add up to 1
            iconstraints = vars[type == "haplotype" & snode.id > 0,
                                .(value = 1, id, cid = paste("haplotype.node", og.node.id))]
            rhs = unique(vars[type == "haplotype",
                              .(value = 1, sense = "E", cid = paste("haplotype.node", og.node.id))],
                         by = "cid")

            constraints = rbind(constraints, iconstraints, fill = TRUE)
            b = rbind(b, rhs, fill = TRUE)

            ## check that there are two per og edge id
            ## browser()
            ## tmp = iconstraints[, .(count = .N), by = cid]
            ## all(tmp[, count] == 2, na.rm = TRUE)
            ## length(unique(rhs[, cid])) == length(unique(iconstraints[, cid]))
            ## length(unique(rhs[, cid])) == length(unique(gg$nodes$dt[allele == "major", og.node.id]))

            ## add H1 AND constraint
            h1.and.ids = merge.data.table(vars[type == "h1.and.indicator", .(n1, n2, edge.id = id, sedge.id)],
                                          vars[type == "haplotype", .(n1.snode.id = snode.id, n1.id = id)],
                                          by.x = "n1",
                                          by.y = "n1.snode.id") %>%
                merge.data.table(vars[type == "haplotype", .(n2.snode.id = snode.id, n2.id = id)],
                                 by.x = "n2",
                                 by.y = "n2.snode.id")

            h2.and.ids = merge.data.table(vars[type == "h2.and.indicator", .(n1, n2, edge.id = id, sedge.id)],
                                          vars[type == "haplotype", .(n1.snode.id = snode.id, n1.id = id)],
                                          by.x = "n1",
                                          by.y = "n1.snode.id") %>%
                merge.data.table(vars[type == "haplotype", .(n2.snode.id = snode.id, n2.id = id)],
                                 by.x = "n2",
                                 by.y = "n2.snode.id")

            ## verify only + sedge id
            ## browser()

            ## there are four constraints that are needed to implement this first edge constraint (c1-3)
            iconstraints = rbind(h1.and.ids[, .(value = 1, id = edge.id, cid = paste("h1.and.c1", sedge.id))],
                                 h1.and.ids[, .(value = -1, id = n1.id, cid = paste("h1.and.c1", sedge.id))],
                                 h1.and.ids[, .(value = 1, id = edge.id, cid = paste("h1.and.c2", sedge.id))],
                                 h1.and.ids[, .(value = -1, id = n2.id, cid = paste("h1.and.c2", sedge.id))],
                                 h1.and.ids[, .(value = 1, id = edge.id, cid = paste("h1.and.c3", sedge.id))],
                                 h1.and.ids[, .(value = -1, id = n1.id, cid = paste("h1.and.c3", sedge.id))],
                                 h1.and.ids[, .(value = -1, id = n2.id, cid = paste("h1.and.c3", sedge.id))])

            rhs = rbind(h1.and.ids[, .(value = 0, sense = "L", cid = paste("h1.and.c1", sedge.id))],
                        h1.and.ids[, .(value = 0, sense = "L", cid = paste("h1.and.c2", sedge.id))],
                        h1.and.ids[, .(value = -1, sense = "G", cid = paste("h1.and.c3", sedge.id))])

            constraints = rbind(constraints, iconstraints, fill = TRUE)
            b = rbind(b, rhs, fill = TRUE)

            ## tmp = iconstraints[, .(count = .N), by = cid]
            ## all(tmp[cid %like% "c1", count] == 2)
            ## all(tmp[cid %like% "c2", count] == 2)
            ## all(tmp[cid %like% "c3", count] == 3)


            iconstraints = rbind(h2.and.ids[, .(value = 1, id = edge.id, cid = paste("h2.and.c1", sedge.id))],
                                 h2.and.ids[, .(value = 1, id = n1.id, cid = paste("h2.and.c1", sedge.id))],
                                 h2.and.ids[, .(value = 1, id = edge.id, cid = paste("h2.and.c2", sedge.id))],
                                 h2.and.ids[, .(value = 1, id = n2.id, cid = paste("h2.and.c2", sedge.id))],
                                 h2.and.ids[, .(value = 1, id = edge.id, cid = paste("h2.and.c3", sedge.id))],
                                 h2.and.ids[, .(value = 1, id = n1.id, cid = paste("h2.and.c3", sedge.id))],
                                 h2.and.ids[, .(value = 1, id = n2.id, cid = paste("h2.and.c3", sedge.id))])

            rhs = unique(rbind(h2.and.ids[, .(value = 1, sense = "L", cid = paste("h2.and.c1", sedge.id))],
                               h2.and.ids[, .(value = 1, sense = "L", cid = paste("h2.and.c2", sedge.id))],
                               h2.and.ids[, .(value = 1, sense = "G", cid = paste("h2.and.c3", sedge.id))]),
                         by = "cid")

            constraints = rbind(constraints, iconstraints, fill = TRUE)
            b = rbind(b, rhs, fill = TRUE)

            ## verify that there are no weird NA's and that there is only one set of constraints per sedge.id
            ## tmp = iconstraints[, .(count = .N), by = cid]
            ## all(tmp[cid %like% "c1", count] == 2)
            ## all(tmp[cid %like% "c2", count] == 2)
            ## all(tmp[cid %like% "c3", count] == 3)

            ## connect edge indicators to the haplotype configuration of connected edges
            iconstraints = rbind(vars[type == "h1.and.indicator",
                                      .(value = -1, id, cid = paste("haplotype.indicator", sedge.id))],
                                 vars[type == "h2.and.indicator",
                                      .(value = -1, id, cid = paste("haplotype.indicator", sedge.id))],
                                 vars[type == "edge.indicator" &
                                      (sedge.id %in% vars[type == "h1.and.indicator",]$sedge.id),
                                      .(value = 1, id, cid = paste("haplotype.indicator", sedge.id))])
            rhs = unique(iconstraints[, .(value = 0, sense = "L", cid)], by = "cid")

            constraints = rbind(constraints, iconstraints, fill = TRUE)
            b = rbind(b, rhs, fill = TRUE)

            ## verify that there are three of these per sedge.id!
            ## tmp = iconstraints[, .(count = .N), by = cid]
            ## all(tmp[, count] == 3)

            ## add constraints that force indicators to be 1 if edge CN > 0

            ## add constraints for upper bound (same setup as L0 penalty) - one per edge
            iconstraints = vars[type == "edge", .(value = 1, id,
                                                  sedge.id,
                                                  cid = paste("edge.indicator.ub", sedge.id))]

            ## add matching indicator variables, matching by cid
            iconstraints = rbind(
                iconstraints,
                vars[type == "edge.indicator", ][
                    sedge.id %in% iconstraints$sedge.id, .(value = -M, id, cid = iconstraints$cid, sedge.id)],
                fill = TRUE)

            ## upper bound is M if indicator is positive, and zero otherwise
            constraints = rbind(
                constraints,
                iconstraints,
                fill = TRUE)

            ## add the RHS of this constraint (upper bound)
            b = rbind(
                b,
                vars[type == "edge", .(value = 0, sense = "L", cid = paste("edge.indicator.ub", sedge.id))],
                fill = TRUE
            )

            ## add constraints for the lower bound
            iconstraints = vars[type == "edge",
                                .(value = 1, id, sedge.id, cid = paste("edge.indicator.lb", sedge.id))]

            ## add matching indicator variables for LB
            iconstraints = rbind(
                iconstraints,
                vars[type == "edge.indicator", ][sedge.id %in% iconstraints$sedge.id,
                                                 .(value = -0.1, id, cid = iconstraints$cid, sedge.id)],
                fill = TRUE)

            constraints = rbind(
                constraints,
                iconstraints,
                fill = TRUE)

            ## add the RHS of this constraint (upper bound)
            b = rbind(
                b,
                vars[type == "edge", .(value = 0, sense = "G", cid = paste("edge.indicator.lb", sedge.id))],
                fill = TRUE
            )
        }

        ## implement edge indicators for ISM and edge reward
        if (ism | any(gg$edges$dt$reward != 0, na.rm = TRUE)) {

            ## implement edge edge indicators if not already (e.g. if not doing phasing)
            if (!phased) {

                ## importantly, we only want to add these for ALT edges
                iconstraints = vars[type == "edge" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                                    .(value = 1, id,
                                      sedge.id,
                                      cid = paste("edge.indicator.ub", sedge.id))]

                ## add matching indicator variables, matching by cid
                iconstraints = rbind(
                    iconstraints,
                    vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id) == 1, ][
                        sedge.id %in% iconstraints$sedge.id,
                        .(value = -M, id, cid = iconstraints$cid, sedge.id)],
                    fill = TRUE)

                ## upper bound is M if indicator is positive, and zero otherwise
                constraints = rbind(
                    constraints,
                    iconstraints,
                    fill = TRUE)

                ## add the RHS of this constraint (upper bound)
                b = rbind(
                    b,
                    vars[type == "edge" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                         .(value = 0, sense = "L", cid = paste("edge.indicator.ub", sedge.id))],
                    fill = TRUE
                )

                ## add constraints for the lower bound
                iconstraints = vars[type == "edge" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                                    .(value = 1, id, sedge.id, cid = paste("edge.indicator.lb", sedge.id))]

                ## add matching indicator variables for LB
                iconstraints = rbind(
                    iconstraints,
                    vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id) == 1, ][
                        sedge.id %in% iconstraints$sedge.id,
                        .(value = -0.1, id, cid = iconstraints$cid, sedge.id)],
                    fill = TRUE)

                constraints = rbind(
                    constraints,
                    iconstraints,
                    fill = TRUE)

                ## add the RHS of this constraint (upper bound)
                b = rbind(
                    b,
                    vars[type == "edge" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                         .(value = 0, sense = "G", cid = paste("edge.indicator.lb", sedge.id))],
                    fill = TRUE
                )
            }

            ## fix loose ends at zero if there's a junction there (only valid if not phasing)
            #' zchoo Tuesday, Jun 15, 2021 11:53:15 AM
            #' this constraint appears to be valid even if running phasing.
            ## if (!phased) {
            ## extremity exclusivity (relevant for ALL graphs)

            if (ism) {
                loose.constraints = rbind(
                    vars[type == "loose.in.indicator" & sign(snode.id) == 1 & telomeric == FALSE,
                         .(value = 1, id, cid = paste("extremity.exclusivity", ee.id))],
                    vars[type == "loose.out.indicator" & sign(snode.id) == 1 & telomeric == FALSE,
                         .(value = 1, id, cid = paste("extremity.exclusivity", ee.id))]
                )

                edge.constraints = rbind(
                    vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                         .(value = 1, id, cid = paste("extremity.exclusivity", ee.id.n1))],
                    vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                         .(value = 1, id, cid = paste("extremity.exclusivity", ee.id.n2))]
                )

                constraints = rbind(constraints, loose.constraints, edge.constraints, fill = TRUE)

                loose.b = unique(loose.constraints[, .(cid, value = 1, sense = "L")], by = "cid")
                edge.b = unique(edge.constraints[, .(cid, value = 1, sense = "L")], by = "cid")

                b = rbind(b, edge.b, loose.b, fill = TRUE)

                ## fix loose ends at zero if they coincide with a called junction
                edge.ee.ids = unique(c(vars[type == "edge.indicator", ee.id.n1], vars[type == "edge.indicator", ee.id.n2]))
                edge.ee.ids = edge.ee.ids[!is.na(edge.ee.ids)]

                loose.zeros = rbind(
                    vars[type == "loose.in.indicator" & sign(snode.id) == 1 & ee.id %in% edge.ee.ids,
                         .(value = 1, id, cid = paste("extremity.exclusivity", ee.id))],
                    vars[type == "loose.out.indicator" & sign(snode.id) == 1 & ee.id %in% edge.ee.ids,
                         .(value = 1, id, cid = paste("extremity.exclusivity", ee.id))]
                )

                loose.zeros.rhs = unique(loose.zeros[, .(cid, value = 0, sense = "E")], by = "cid")

                constraints = rbind(constraints, loose.zeros, fill = TRUE)
                b = rbind(b, loose.zeros.rhs, fill = TRUE)
            }
        }

        if (phased) {
            ## homologous extremity exclusivity
            ## this is actually redundant with previous constraints

            if (ism) {
                loose.constraints = rbind(
                    vars[type == "loose.in.indicator" & sign(snode.id)==1 & telomeric == FALSE,
                         .(value = 1, id, cid = paste("homol.extremity.exclusivity", hee.id))],
                    vars[type == "loose.out.indicator" & sign(snode.id)==1 & telomeric == FALSE,
                         .(value = 1, id, cid = paste("homol.extremity.exclusivity", hee.id))]
                )

                edge.constraints = rbind(
                    vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id)==1,
                         .(value = 1, id, cid = paste("homol.extremity.exclusivity", hee.id.n1))],
                    vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id)==1,
                         .(value = 1, id, cid = paste("homol.extremity.exclusivity", hee.id.n2))]
                )

                ## we should allow loose ends to violate this as some loose ends are germline
                ## therefore, loose ends with id's not shared with edge constraints
                loose.constraints = loose.constraints[(cid %in% edge.constraints[, cid])]

                ## add these constraints to the existing table
                constraints = rbind(constraints, loose.constraints, edge.constraints, fill = TRUE)

                rhs = unique(rbind(
                    vars[type == "loose.in.indicator" & sign(snode.id)==1 & telomeric == FALSE,
                         .(value = 1, sense = "L", cid = paste("homol.extremity.exclusivity", hee.id))],
                    vars[type == "loose.out.indicator" & sign(snode.id)==1 & telomeric == FALSE,
                         .(value = 1, sense = "L", cid = paste("homol.extremity.exclusivity", hee.id))]
                ), by = "cid")

                rhs = rhs[cid %in% edge.constraints[, cid]]

                b = rbind(b, rhs, fill = TRUE)

                if (verbose) {
                    message("Number of homologous extremity exclusivity constraints: ",
                            nrow(rhs))
                }

                config.dt = vars[type == "straight.config" | type == "cross.config",]

                config.constraints.lt = rbind(
                    vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "REF",
                         .(value = -1, id, cid = paste("config lt", sedge.id))],
                    vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "REF",
                         .(value = 1, id = config.dt$id[match(config.id, config.dt$config.id)],
                           cid = paste("config lt", sedge.id))])

                config.constraints.gt = rbind(
                    vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "REF",
                         .(value = -1, id, cid = config.id)],
                    vars[type == "straight.config" & sedge.id > 0 & ref.or.alt == "REF",
                         .(value = 1, id, cid = config.id)],
                    vars[type == "cross.config" & sedge.id > 0 & ref.or.alt == "REF",
                         .(value = 1, id, cid = config.id)])

                rhs = unique(rbind(
                    config.constraints.lt[, .(cid, value = 0, sense = "G")],
                    config.constraints.gt[, .(cid, value = 0, sense = "L")]),
                    by = "cid")

                constraints = rbind(constraints, config.constraints.lt, config.constraints.gt, fill = TRUE)
                b = rbind(b, rhs, fill = TRUE)

                ## implement reciprocal homologous extremity exclusivity
                straight.config.dt = vars[type == "straight.config",]
                cross.config.dt = vars[type == "cross.config",]

                rhomol.constraints = rbind(
                    ## corresponding cross indicator
                    vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "REF" & connection == "straight",
                         .(value = 1, id = cross.config.dt$id[match(og.edge.id, cross.config.dt$og.edge.id)],
                           cid = paste("rhee", sedge.id))],

                    ## corresponding cross indicator
                    vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "REF" & connection == "cross",
                         .(value = 1, id = straight.config.dt$id[match(og.edge.id, straight.config.dt$og.edge.id)],
                           cid = paste("rhee", sedge.id))],

                    ## actual ALT edges
                    vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(s1),
                         .(value = 1, id, cid = paste("rhee", s1))],
                    vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(s2),
                         .(value = 1, id, cid = paste("rhee", s2))],
                    vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(s3),
                         .(value = 1, id, cid = paste("rhee", s3))],
                    vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(s4),
                         .(value = 1, id, cid = paste("rhee", s4))],
                    vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(c1),
                         .(value = 1, id, cid = paste("rhee", c1))],
                    vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(c2),
                         .(value = 1, id, cid = paste("rhee", c2))],
                    vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(c3),
                         .(value = 1, id, cid = paste("rhee", c3))],
                    vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(c4),
                         .(value = 1, id, cid = paste("rhee", c4))],

                    ## loose indicators
                    vars[type == "loose.in.indicator" & snode.id > 0 & !is.na(s) & telomeric == FALSE,
                         .(value = 1, id, cid = paste("rhee", s))],
                    vars[type == "loose.in.indicator" & snode.id > 0 & !is.na(c) & telomeric == FALSE,
                         .(value = 1, id, cid = paste("rhee", c))],
                    vars[type == "loose.out.indicator" & snode.id > 0 & !is.na(s) & telomeric == FALSE,
                         .(value = 1, id, cid = paste("rhee", s))],
                    vars[type == "loose.out.indicator" & snode.id > 0 & !is.na(c) & telomeric == FALSE,
                         .(value = 1, id, cid = paste("rhee", c))]
                )

                rhs = unique(rhomol.constraints[, .(value = 2, sense = "L", cid)], by = "cid")

                if (verbose) {
                    message("Number of reciprocal homologous constraints: ", nrow(rhs))
                }

                constraints = rbind(constraints, rhomol.constraints, fill = TRUE)
                b = rbind(b, rhs, fill = TRUE)

            }

            ## add the edge indicator sum constraints (ISM consistency)
            iconstraints = unique(
                vars[type == "edge.indicator" & ref.or.alt == "ALT",
                     .(value = 1, id, og.edge.id,
                       edge.id = abs(sedge.id),
                       cid = paste("edge.indicator.sum.ub", og.edge.id))],
                by = "edge.id"
            )

            constraints = rbind(
                constraints,
                iconstraints[, .(value, id, cid)],
                fill = TRUE)

            edge.indicator.b = unique(
                vars[type == "edge.indicator" & ref.or.alt == "ALT",
                     .(value = 1, sense = "L", cid = paste("edge.indicator.sum.ub", og.edge.id))],
                by = "cid"
            )

            b = rbind(b, edge.indicator.b, fill = TRUE)

            ## force major allele to have higher CN than minor allele
            ## may not work for phased blocks
            if (force.major) {

                iconstraints = rbind(
                    vars[type == "node" & allele == "major" & snode.id > 0,
                         .(value = 1, id, cid = paste("force.major", og.node.id))],
                    vars[type == "node" & allele == "minor" & snode.id > 0,
                         .(value = -1, id, cid = paste("force.major", og.node.id))])

                rhs = unique(vars[type == "node" & snode.id > 0 & allele == "major",
                                  .(value = 0, sense = "G", cid = paste("force.major", og.node.id))],
                             by = "cid")

                constraints = rbind(constraints, iconstraints, fill = TRUE)
                b = rbind(b, rhs, fill = TRUE)

            }


            ## force nonzero CN for ALT edges (because these have nonzero CN in original JaBbA output)
            ## can become infeasible if original graph is not compatible with ISM
            if (force.alt) {

                if (ism) {
                    warning("Forcing ALT edges while running ISM can make some problems infeasible!")
                }

                iconstraints = unique(
                    vars[type == "edge.indicator" & ref.or.alt == "ALT" & cnloh != TRUE & sedge.id > 0,
                         .(value = 1, id, og.edge.id,
                           edge.id = abs(sedge.id),
                           cid = paste("edge.indicator.sum.lb", og.edge.id))],
                    by = "edge.id"
                )

                constraints = rbind(
                    constraints,
                    iconstraints[, .(value, id, cid)],
                    fill = TRUE)

                edge.indicator.b = unique(
                    vars[type == "edge.indicator" & ref.or.alt == "ALT" & cnloh != TRUE & sedge.id > 0,
                         .(value = 1, sense = "G", cid = paste("edge.indicator.sum.lb", og.edge.id))],
                    by = "cid"
                )

                b = rbind(b, edge.indicator.b, fill = TRUE)
            }
        }


        if (L0) ## add "big M" constraints
        {
            ## indicator constraints ie on ulids
            iconstraints = rbind(
                vars[type == 'loose.out', .(value = 1, id, ulid, cid = paste('loose.out.indicator.ub', ulid))],
                vars[type == 'loose.in', .(value = 1, id, ulid, cid = paste('loose.in.indicator.ub', ulid))],
                fill = TRUE)

            ## add the matching indicator variables, matching to the cid from above
            iconstraints = rbind(
                iconstraints,
                vars[type %in% c('loose.out.indicator', 'loose.in.indicator'), ][
                    match(iconstraints$ulid, ulid), .(value = -M, id, cid = iconstraints$cid)],
                fill = TRUE)

            ## upper bounds "infinity" ie M if indicator positive, 0 otherwise
            constraints = rbind(
                constraints,
                iconstraints,
                fill = TRUE)

        }

        if (L0) ## add "big M" constraints
        {
            ## indicator constraints ie on ulids
            iconstraints = rbind(
                vars[type == 'loose.out', .(value = 1, id, ulid, cid = paste('loose.out.indicator.ub', ulid))],
                vars[type == 'loose.in', .(value = 1, id, ulid, cid = paste('loose.in.indicator.ub', ulid))],
                fill = TRUE)

            ## add the matching indicator variables, matching to the cid from above
            iconstraints = rbind(
                iconstraints,
                vars[type %in% c('loose.out.indicator', 'loose.in.indicator'), ][
                    match(iconstraints$ulid, ulid), .(value = -M, id, cid = iconstraints$cid)],
                fill = TRUE)

            ## upper bounds "infinity" ie M if indicator positive, 0 otherwise
            constraints = rbind(
                constraints,
                iconstraints,
                fill = TRUE)

            ## upper bound sense is 'L' i.e. less than because -M on left hand side
            b = rbind(b,
                      vars[type == 'loose.in', .(value = 0, sense = 'L', cid = paste('loose.in.indicator.ub', ulid))],
                      vars[type == 'loose.out', .(value = 0, sense = 'L', cid = paste('loose.out.indicator.ub', ulid))],
                      fill = TRUE)

            ## lower bound 0.1 if indicator positive, 0 otherwise
            iconstraints = rbind(
                vars[type == 'loose.out', .(value = 1, id, ulid, cid = paste('loose.out.indicator.lb', ulid))],
                vars[type == 'loose.in', .(value = 1, id, ulid, cid = paste('loose.in.indicator.lb', ulid))],
                fill = TRUE)

            ## add the matching indicator variables, matching to the cid from above
            iconstraints = rbind(
                iconstraints,
                vars[type %in% c('loose.out.indicator', 'loose.in.indicator'), ][
                    match(iconstraints$ulid, ulid), .(value = -.1, id, cid = iconstraints$cid)],
                fill = TRUE)

            ## upper bounds "infinity" ie M if indicator positive, 0 otherwise
            constraints = rbind(
                constraints,
                iconstraints,
                fill = TRUE)

            ## lower bound sense is 'G' i.e. greater than because -M on left hand side
            b = rbind(b,
                      vars[type == 'loose.in', .(value = 0, sense = 'G', cid = paste('loose.in.indicator.lb', ulid))],
                      vars[type == 'loose.out', .(value = 0, sense = 'G', cid = paste('loose.out.indicator.lb', ulid))],
                      fill = TRUE)

            if (loose.collapse)
            {
    ##################
                ## loose indicator sum  = sum of indicators
    ##################
                iconstraints = rbind(
                    vars[type == 'loose.out.indicator', .(value = 1, id, lid, cid = paste('loose.out.indicator.sum', lid))],
                    vars[type == 'loose.in.indicator', .(value = 1, id, lid, cid = paste('loose.in.indicator.sum', lid))],
                    fill = TRUE)

                ## indicator sum is the sum of all indicators mapping to that loose end
                iconstraints = rbind(
                    iconstraints,
                    unique(vars[type %in% c('loose.out.indicator.sum', 'loose.in.indicator.sum'), ][
                        match(iconstraints$lid, lid), .(value = -1, id, lid, cid = iconstraints$cid)], by = 'lid'),
                    fill = TRUE)

                constraints = rbind(
                    constraints,
                    iconstraints,
                    fill = TRUE)

                b = rbind(b,
                          vars[type == 'loose.in.indicator.sum', .(value = 0, sense = 'E', cid = paste('loose.in.indicator.sum', lid))],
                          vars[type == 'loose.out.indicator.sum', .(value = 0, sense = 'E', cid = paste('loose.out.indicator.sum', lid))],
                          fill = TRUE)

    ##################
                ## now we make new indicator variables on the sum of the individual loose end indicators
                ## upper bound bound 0.1 if indicator positive, 0 otherwise
    ##################

                iconstraints = rbind(
                    vars[type == 'loose.out.indicator.sum', .(value = 1, id, lid, cid = paste('loose.out.indicator.sum.indicator.ub', lid))],
                    vars[type == 'loose.in.indicator.sum', .(value = 1, id, lid, cid = paste('loose.in.indicator.sum.indicator.ub', lid))],
                    fill = TRUE)

                ## add the matching indicator variables, matching to the cid from above
                iconstraints = rbind(
                    iconstraints,
                    vars[type %in% c('loose.out.indicator.sum.indicator', 'loose.in.indicator.sum.indicator'), ][
                        match(iconstraints$lid, lid), .(value = -M, id, lid, cid = iconstraints$cid)],
                    fill = TRUE)

                ## upper bounds "infinity" ie M if indicator positive, 0 otherwise
                constraints = rbind(
                    constraints,
                    iconstraints,
                    fill = TRUE)

                ## upper bound sense is 'L' i.e. less than because -M on left hand side
                b = rbind(b,
                          vars[type == 'loose.in.indicator.sum', .(value = 0, sense = 'L', cid = paste('loose.in.indicator.sum.indicator.ub', lid))],
                          vars[type == 'loose.out.indicator.sum', .(value = 0, sense = 'L', cid = paste('loose.out.indicator.sum.indicator.ub', lid))],
                          fill = TRUE)

                ## lower bound 0.1 if indicator positive, 0 otherwise
                iconstraints = rbind(
                    vars[type == 'loose.out.indicator.sum', .(value = 1, id, lid, cid = paste('loose.out.indicator.sum.indicator.lb', lid))],
                    vars[type == 'loose.in.indicator.sum', .(value = 1, id, lid, cid = paste('loose.in.indicator.sum.indicator.lb', lid))],
                    fill = TRUE)

                ## add the matching indicator variables, matching to the cid from above
                iconstraints = rbind(
                    iconstraints,
                    vars[type %in% c('loose.out.indicator.sum', 'loose.in.indicator.sum'), ][
                        match(iconstraints$lid, lid), .(value = -.1, id, lid, cid = iconstraints$cid)],
                    fill = TRUE)

                ## upper bounds "infinity" ie M if indicator positive, 0 otherwise
                constraints = rbind(
                    constraints,
                    iconstraints,
                    fill = TRUE)

                ## lower bound sense is 'G' i.e. greater than because -M on left hand side
                b = rbind(b,
                          vars[type == 'loose.in.indicator.sum', .(value = 0, sense = 'G', cid = paste('loose.in.indicator.sum.indicator.lb', lid))],
                          vars[type == 'loose.out.indicator.sum', .(value = 0, sense = 'G', cid = paste('loose.out.indicator.sum.indicator.lb', lid))],
                          fill = TRUE)

            }
        }


        if (!is.null(marginal) && length(dmarginal))
        {
            ## match against nodes and store query.id as rid
            ## this will be the constraint id that will allow us
            ## to sum the appropriate nodes to constrain to the residual
            ov = dmarginal[, c('cn', 'weight')] %*% gg$nodes$gr %>% gr2dt

            ov[, rid := query.id]

            constraints = rbind(
                constraints,
                rbind(
                    ## match up vars and marginal by snode.id and populate coefficients
                    merge.data.table(vars[type == 'node', !"rid"], ov, by = 'snode.id')[, .(value = 1, id , cid = paste('mresidual', rid))],
                    ## the residual is the difference between the sum and marginal cn
                    vars[type == 'mresidual' & rid %in% ov$rid, .(value = -1, id, cid = paste('mresidual', rid))],
                    fill = TRUE),
                fill = TRUE
            )

            b = rbind(b,
                      vars[type == 'mresidual' & rid %in% ov$rid, .(value = cn, sense = 'E', cid = paste('mresidual', rid))],
                      fill = TRUE)
        }

        if (!is.null(emarginal)) {

            emconstraints = rbind(
                vars[type == "edge", .(value = 1, id, cid = paste("emresidual", emarginal.id))],
                vars[type == "emresidual", .(value = -1, id, cid = paste("emresidual", emarginal.id))]
            )

            constraints = rbind(constraints, emconstraints, fill = TRUE)

            emb = vars[type == "emresidual", .(value = cn, sense = "E", cid = paste("emresidual", emarginal.id))]

            b = rbind(emb, b, fill = TRUE)
        }

    ########
        ## MAKE MATRICES
    ########

        ## now Rcplex time
        ## remove any rows with b = NA
        ## get rid of any constraints with NA values
        keep.constraints = intersect(b[!is.na(value), cid], constraints[!is.na(value), cid])
        b = b[cid %in% keep.constraints,]
        constraints = constraints[cid %in% keep.constraints,]

        ## convert constraints to integers
        ucid = unique(b$cid)
        b[, cid.char := cid]
        b[, cid := cid %>% factor(ucid) %>% as.integer]
        constraints[, cid.char := cid]
        constraints[, cid := cid %>% factor(ucid) %>% as.integer]

        pmt = match(ucid, b$cid.char) ## get right permutation
        bvec = b[pmt, value]
        sense = b[pmt, sense]
        if (verbose) {
            message("Unique cids (A): ", length(unique(constraints$cid)))
            message("Unique cids (b): ", length(unique(b$cid)))
            message("Number of variables: ", length(unique(constraints$id)))
        }

        ## create constraint matrix, Qmat, and cobj, lb, ub from vars and constraints  lambda = 10
        Amat = sparseMatrix(constraints$cid, constraints$id, x = constraints$value, dims = c(length(ucid), nrow(vars)))
        vars[is.na(weight), weight := 0]

        if (verbose) {

            message("bvec length: ", length(bvec))
            message("Amat nrow: ", nrow(Amat))

        }
        if (any(ix <- is.infinite(vars$weight)))
        {
            warning('nodes with infinite weight, setting to 0, please check inputs')
            vars[ix, weight := 0]
        }
        Qmat = vars[, weight * (type %in% c('nresidual', 'eresidual', 'mresidual'))] %>% as.numeric %>% Diagonal(x = .) %>% as('CsparseMatrix')

        ## set lambda to 0 at terminal or other non NA nodes
        vars[is.na(lambda), lambda := 0]


        ## set cvec by multiplying global lambda by local lambda for non-terminal loose end
        ## vars (or their indicators if L0 is TRUE)
        if (L0)
        {
            if (loose.collapse)
            {
                cvec = lambda*(vars[, lambda*(type %in% c('loose.in.indicator.sum.indicator', 'loose.out.indicator.sum.indicator') & !terminal)] %>% as.numeric)
                ## cvec = lambda*(vars[, lambda*(type %in% c('loose.in.indicator.sum.indicator', 'loose.out.indicator.sum.indicator', 'loose.in.indicator', 'loose.out.indicator') & !terminal)] %>% as.numeric)
            }
            else
            {
                cvec = lambda*(vars[, lambda * (type %in% c('loose.in.indicator', 'loose.out.indicator') & !terminal)] %>% as.numeric)
            }
        } else {
            cvec = lambda*(vars[, lambda*(type %in% c('loose.in', 'loose.out') & !terminal)] %>% as.numeric)
        }

        ## message("CVEC: ", length(cvec))

        if (length(indices <- which(vars[, type == "edge.indicator" & reward != 0])))
        {
            if (verbose) {
            }
            message('Applying reward')
            ## grab edge indicator variables with reward
            cvec[indices] = -vars$reward[indices]
        }

        if (lp) {
            ## add weights of stuff
            indices = which(vars$type %in% c("mdelta.plus", "mdelta.minus",
                                             "ndelta.plus", "ndelta.minus",
                                             "edelta.plus", "edelta.minus"))
            wts = vars$weight[indices]
            cvec[indices] = wts
            Qmat = NULL ## no Q if solving LP
        }

        ## browser()
        if (cnloh) {

            if ("cnloh" %in% colnames(vars)) {
                indices = which(vars$type == "edge.indicator" & !is.na(vars$cnloh) & vars$cnloh == TRUE)
                cvec[indices] = lambda

                message("Number of penalized CNLOH edges: ", length(indices))
            }
        }

        ## check constraints of CNLOH
        ## browser()
        ## vars[type == "edge.indicator" & cnloh == TRUE]
        ## vars[type == "edge.indicator" & cnloh == TRUE, .N, by = og.edge.id]


        lb = vars$lb
        ub = vars$ub

        control = list(trace = ifelse(verbose>=2, 1, 0), tilim = tilim, epgap = epgap, round = 1, trelim = trelim, nodefileind = nodefileind, method = 4)

        ## call our wrapper for CPLEX
        if (use.gurobi) {

            if (verbose) { message("Starting optimization with gurobi!") }

            sol = run_gurobi(cvec = cvec,
                             Amat = Amat,
                             bvec = bvec,
                             Qmat = Qmat,
                             lb = lb,
                             ub = ub,
                             sense = sense,
                             vtype = vars$vtype,
                             objsense = "min",
                             control = control)
        } else {

            if (verbose) { message("Starting optimization with CPLEX!") }

            sol =  gGnome:::Rcplex2(cvec,
                                    Amat,
                                    bvec,
                                    Qmat = Qmat,
                                    lb = lb,
                                    ub = ub,
                                    sense = sense,
                                    vtype = vars$vtype,
                                    objsense = "min",
                                    control = control,
                                    tuning = FALSE)
        }

        vars$cvec = cvec
        vars$x = sol$x

        ## for debugging
        ppc = function(x) (x %>% merge(vars, by = 'id') %>% merge(b, by = 'cid.char'))[, paste(paste(round(value.x, 1), '*', paste(type, gid, sep=  '_'), '(', signif(x, 2), ')', collapse = ' + '), ifelse(sense[1] == 'E', '=', ifelse(sense[1] == 'G', '>=', '<=')), round(value.y[1],2)), by = cid.char]

        ppv = function(x) {tmp = x %>% merge(constraints, by = 'id'); constraints[cid %in% tmp$cid, ] %>% ppc}

        .check = function(x) data.table(obs = sign(as.numeric(round(Amat %*% x - bvec))),
                                        sense)
        chk = .check(sol$x)

        if (any(is.na(sol$x)))
            stop('Rcplex did not converge or failed to find a solution, please run with verbose = 2 to get more detailed output')

        if (chk[sense == 'E', any(obs != 0, na.rm = TRUE)] |
            chk[sense == 'G', any(obs < 0, na.rm = TRUE)] |
            chk[sense == 'L', any(obs > 0, na.rm = TRUE)])
            stop('Constraint violation likely due to M parameter being too large for problem causing CPLEX numerical instability, consider lowering M parameter')

        ##.obj = function(x) 0.5 * rbind(x) %*% Qmat %*% cbind(x) + cvec %*% x

        ## mark haplotypes if phasing
        if (phased) {
            haplotypes.dt = vars[type == "haplotype" & snode.id > 0,
                                 .(node.id = snode.id,
                                   haplotype = ifelse(x == 1, "h1", "h2"),
                                   col = ifelse(x == 1, alpha("red", 0.5), alpha("blue", 0.5)))]
            gg$nodes[haplotypes.dt$node.id]$mark(haplotype = haplotypes.dt$haplotype)
            gg$nodes[haplotypes.dt$node.id]$mark(col = haplotypes.dt$col)
        }

        ## update graph
        nmark = vars[type == 'node', .(nid = abs(snode.id), cn = round(x))]
        emark = vars[type == 'edge', .(eid = abs(sedge.id), cn = round(x))]

        loosei = vars[type == 'loose.in' & snode.id>0, .(cn = round(x)), keyby = snode.id]
        looseo = vars[type == 'loose.out' & snode.id>0, .(cn = round(x)), keyby = snode.id]

        nodes = gg$nodes[loosei$snode.id] ## need to do this to use nodes active binding settings
        nodes$loose.left = loosei$cn>0

        nodes = gg$nodes[looseo$snode.id] ## need to do this to use nodes active binding settings
        nodes$loose.right = looseo$cn>0

        gg$nodes$mark(loose.cn.left = 0, loose.cn.right = 0)
        gg$nodes[loosei$snode.id]$mark(loose.cn.left = loosei$cn)
        gg$nodes[looseo$snode.id]$mark(loose.cn.right = looseo$cn)

        ## cache old cn values
        gg$nodes$mark(cn.old = gg$nodes$dt$cn)
        gg$edges$mark(cn.old = gg$edges$dt$cn)
        gg$nodes$mark(cn = NULL) ## reset to avoid weird type casting issue
        gg$edges$mark(cn = NULL) ## reset to avoid weird type casting issue
        gg$nodes[nmark$nid]$mark(cn = nmark$cn)
        gg$edges[emark$eid]$mark(cn = emark$cn)
        gg$set(y.field = 'cn')

        gg$set(obj = sol$obj)
        gg$set(status = sol$status)
        gg$set(epgap = sol$epgap)
        if (!use.gurobi) {
            gg$set(code = readRDS(system.file('extdata', 'cplex_codes.rds', package="gGnome"))[.(sol$status), code])
        }

        if (verbose) {
          message("CPLEX epgap ", sol$epgap, " with solution status ", gg$meta$code)
        }

        ##  fix loose ends
        nodes = gg$nodes
        nodes$loose.left = nodes$dt$loose.cn.left>0
        nodes$loose.right = nodes$dt$loose.cn.right>0

        ## if phased, mark edges with different colors to make it easier to visualize
        if (phased) {
            if (verbose) {
                message("formatting phased graph...")
            }
            ## edge formatting
            ref.edge.col = alpha("blue", 0.5)
            alt.edge.col = alpha("red", 0.5)
            ref.edge.lwd = 1.0
            alt.edge.lwd = 1.0
            edge.col = ifelse(gg$edges$dt$type == "REF", ref.edge.col, alt.edge.col)
            edge.lwd = ifelse(gg$edges$dt$type == "REF", ref.edge.lwd, alt.edge.lwd)
            gg$edges$mark(col = edge.col, lwd = edge.lwd)

            ## mark zero cn edges
            zero.cn.col = alpha("gray", 0.1)
            zero.cn.lwd = 0.5
            zero.cn.edges = which(gg$edges$dt$cn == 0)
            gg$edges[zero.cn.edges]$mark(col = zero.cn.col, lwd = zero.cn.lwd)
        } else {

            ## edge formatting
            ref.edge.col = alpha("blue", 0.2)
            alt.edge.col = alpha("red", 0.4)
            ref.edge.lwd = 0.5
            alt.edge.lwd = 1.0
            edge.col = ifelse(gg$edges$dt$type == "REF", ref.edge.col, alt.edge.col)
            edge.lwd = ifelse(gg$edges$dt$type == "REF", ref.edge.lwd, alt.edge.lwd)
            gg$edges$mark(col = edge.col, lwd = edge.lwd)

            ## mark zero cn edges
            zero.cn.col = alpha("gray", 0)
            zero.cn.lwd = 0.5
            zero.cn.edges = which(gg$edges$dt$cn == 0)
            gg$edges[zero.cn.edges]$mark(col = zero.cn.col, lwd = zero.cn.lwd)
        }

        ## if nonintegral also return the offsets as graph metadata
        ## maybe it will be useful
        if (nonintegral) {
            gg$set(meta = vars[type == "slush", .(chr, offset = x)])
        }

        if (debug) {
            return(list(gg = gg, sol = sol))
        }
        return(gg)
    }

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
