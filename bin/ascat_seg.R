{
    library(optparse)
    
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

    library(devtools)
    library(khtools)
    library(ASCAT)
    library(data.table)
    library(gUtils)
    library(GenomicRanges)
    library(skitools)

    #source("~/modules/ascat_seg/utils.R")

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

    if (grepl(pattern = "txt$", x = opt$variants)) {
        variants.dt = fread(opt$variants)
    	variants.dt[[1]] <- gsub("chr","",variants.dt[[1]])
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
    ## Edit by Tanubrata: Adds a fix to the column names to do GC correction when passing CBS coverge, else
    ## ASCAT breaks when passing raw drycleaned coverage without GC correction
    if (opt$gc) {
        if ("tum.counts" %in% names(values(cov.gr)) & "norm.counts" %in% names(values(cov.gr))) {
            message("Applying GC correction")
            tum.gr = khtools::.gc(cov.gr, "tum.counts")
            norm.gr = khtools::.gc(cov.gr, "norm.counts")
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
