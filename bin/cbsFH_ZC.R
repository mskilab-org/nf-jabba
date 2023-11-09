library(optparse)
options(bitmapType='cairo')
## RLIBDIR = '/cga/meyerson/home/marcin/Software/R/x86_64-unknown-linux-gnu-library/3.1/'
#GIT.HOME = '/cga/meyerson/home/marcin/DB/git'

# run cbs

if (!exists('opt')) ## if opt already exists allow over-ride of command line arg processing for debugging purposes
{
    option_list = list(
        make_option(c("-l", "--libdir"), type = "character", help = "Libdir", default = './'),
        make_option(c("-t", "--tcov"), type = "character", help = "Tumor coverage (GRanges or tab delimited file) OR coverage file"),
        make_option(c("-n", "--ncov"), type = "character", help = "Normal coverage (GRanges or tab delimited file) OR NULL"),
        make_option("--hets", type = "character", help = "text file containing hets", default = "/dev/null"),
        make_option("--hets_thresh", type = "numeric", help = "threshold for counting hets", default = 0.2),
        make_option("--distance_thresh", type = "numeric", help = "distance threshold", default = 1e5),
        make_option(c("--cnsignif"), type = "numeric", default = 1e-5, help = "alpha value for CBS"),
        make_option(c("--mask"), type = "character", default = "/dev/null", help = "path to coverage mask"),
        make_option(c("-m", "--name"), type = "character", default = 'tumor', help = "Sample / Individual name"),
        make_option("--undo_splits", type = "character", default = "none", help = "CBS undo split method"),
        make_option("--undo_prune", type = "numeric", default = 0.1, help = "undo prune SSE thres"),
        make_option("--undo_SD", type = "numeric", default = 3, help = "undo splits SD thres"),
        make_option(c("-o", "--outdir"), type = "character", default = './', help = "Directory to dump output into"),
        make_option(c("-f", "--field"), type = "character", default = 'records', help = "Column of tumor and normal coverage which to use for downstream processing")
    )
    parseobj = OptionParser(option_list=option_list)
    opt = parse_args(parseobj)
}

print(opt)

if (is.null(opt$tcov))
  stop(print_help(parseobj))

writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
saveRDS(opt, paste(opt$outdir, 'cmd.opts', sep = '/'))

#.libPaths(c(.libPaths(), opt$rlibdir))

## print(GIT.HOME)
## CODE.HOME = paste(GIT.HOME,  '/marcin/R/', sep = "");
## ISVA.HOME = paste(GIT.HOME,  '/isva/', sep = "");
## library(methods)
## source(paste(CODE.HOME, 'functions.R', sep = "/"))
## source(paste(GIT.HOME, 'grUtils/grUtils.R', sep = "/"))
## source(paste(GIT.HOME, 'trackData/trackData.R', sep = "/"))
source(paste(opt$libdir, 'JaBbA.R', sep = "/")) ## source local version of karyomip

library(gUtils)
library(gTrack)
##library(JaBbA)
library(skitools)
library(skidb)
library(DNAcopy)
require(Rsamtools)
require(GenomicRanges)
require(parallel)
require(sequenza)
library(data.table)

set.seed(42)

seqinfo2gr = si2gr

if (nchar(Sys.getenv('HG.FFT'))==0)
    Sys.setenv(HG.FFT = '/nethome/mimielinski/DB/ffTracks//hg19.rds')

# tcov$ratio = values(tcov)[, opt$field] / values(ncov)[, opt$field]
out.file.cov = paste(opt$outdir, 'cov.rds', sep = '/')
out.file.rds = paste(opt$outdir, 'seg.rds', sep = '/')
out.file.nseg.rds = paste(opt$outdir, 'nseg.rds', sep = '/')
out.file.txt = paste(opt$outdir, 'seg.txt', sep = '/')
out.file.nozzle = paste(opt$outdir, 'nozzle', sep = '/')

.getcov = function(fn, field, mask)
    {
        if (grepl('\\.rds$', mask) && file.exists(mask) && file.info(mask)$size) {
            mask.gr = readRDS(mask)
        } else {
            mask.gr = GRanges()
        }
        if (grepl('\\.rds$', fn)) ## GRanges rds of coverage
            cov = readRDS(fn)
        else
            {
                if (grepl('\\.csv(\\.gz)?$', fn))  ## assume bedgraph-like format w fourth column having counts
                    {        
                        tmp = read.delim(fn, sep = ',', header = FALSE)
                        names(tmp)[1:4] = c('chr', 'start', 'end', field)
                    }
                else if (grepl('(bedgraph$)', fn) | grepl('\\.tsv(\\.gz)?$', fn))
                    {
                        tmp = read.delim(fn, sep = '\t', header = FALSE)
                        names(tmp)[1:4] = c('chr', 'start', 'end', field)
                    }
                else if (grepl('txt$', fn))
                    tmp = read.delim(fn, header = TRUE)

                cov = seg2gr(tmp)
            }

        if (length(mask.gr)) {
            mask.ix = which(cov %^% mask.gr)
            values(cov)[mask.ix, field] = NA
        }
        
        return(cov)
    }

tcov = .getcov(opt$tcov, opt$field, opt$mask)

if (is.null(opt$ncov))
    opt$ncov = ''
    
if (file.exists(opt$ncov) & file.size(opt$ncov) > 0)
    {
        ncov = .getcov(opt$ncov, opt$field, opt$mask)

        cat('reading', opt$tcov, 'and', opt$ncov, '\n')

#        if (!identical(width(tcov), width(ncov)))
#            stop('Tumor and normal Coverage GRanges objects must be of equal length and width')

        if (!(opt$field %in% names(values(tcov)))| !(opt$field %in% names(values(ncov))))
            stop(sprintf('Opt$Field %s not included in either tumor or normal coverage GRanges object', opt$field))
        
        tmp.s = as.character(seqnames(ncov))
        tmp.e = end(ncov)
        sl = sapply(split(tmp.e, tmp.s), max)
        sl = sl[order(names(sl))] ## order characters
        sl = sl[order(suppressWarnings(as.numeric(names(sl))))] ## order numeric components
        sl = sl[which(!grepl("_", names(sl)))] ### current fix to remove junk chromosomes while remaining relatively reference agnostic
        ## sl = sl[nchar(names(sl))<=3] ### current fix to remove junk chromosomes while remaining relatively reference agnostic
        rm(list=ls(pattern="tmp"))
        gc()

        ## limit chroms to those with read support in at least a "chunk" in the normal
        tmp.chunktile = gr.tile(sl, 5e6) ## broad chromosome pieces to assess whether a chromosome has any read support in any "chunk"
        map = gr.tile.map(tmp.chunktile, ncov, verbose = T); gc() ## Removed mc.cores = 1
        ## Valo = values(ncov)[, opt$field]; gc()
        val = values(ncov)[, opt$field]; gc()
        val[val==0] = NA
        tmp.chunktile$mean = sapply(map, function(x,y) if (length(x)==0) NA else if ((sum(is.na(y[x]))/length(x))<.10) mean(y[x], na.rm = T) else NA, val)
        new.sl = sl[names(sl) %in% as.character(seqnames(tmp.chunktile))[!is.na(tmp.chunktile$mean)]] ## limit genome to only chromosomes with at least a "piece" seen in normal
        rm(list=ls(pattern="tmp"))
        gc()

        ## now compute coverage mean x variance in a (for now) super coarse segmentation of normal (chromosome only) 
        gr.nseg = seqinfo2gr(new.sl)
        tmp.ncov = ncov
        val = values(tmp.ncov)[, opt$field]
        val[val==0] = NA
        values(tmp.ncov)[, opt$field] = val
        gc()
        ss.n = segstats(gr.nseg, tmp.ncov, opt$field, mc.cores = 1, na.thresh = 0.9)
        gc()
        meds = sapply(split(val, as.character(seqnames(tmp.ncov))), median, na.rm = T)
        ss.n$mean = meds[as.character(seqnames(ss.n))] ## replace mean with median, i.e. ghetto robust mean to get rid of crapola, keep the segstats variance
        rm(list=ls(pattern="tmp"))
        gc()
        
        ## pp.n = ppgrid(ss.n, verbose = T, purity.min = 0.5, purity.max = 1, ploidy.max = 6)
        ## nudging it to take the ploidy where purity is closest to 1
        ## pp.n = ppgrid(ss.n, verbose = T, purity.min = 0.5, purity.max = 1, ploidy.max = 3) 
        pp.n = ppgrid(ss.n, verbose = T, purity.min = 0.5,
                      purity.max = 1.5, ploidy.max = 2.5) 
        ## pp.n = ppgrid2(cov = ncov,
        ##                segstats = ss.n,
        ##                verbose = T,
        ##                purity.min = 0.5,
        ##                purity.max = 1,
        ##                ploidy.max = 3,
        ##                mc.cores=5)
        tmp.i = which.min(abs(pp.n$purity-1))
        ss.n$cn = round(rel2abs(ss.n, field = 'mean', gamma = pp.n$gamma[tmp.i], beta = pp.n$beta[tmp.i]))
##        td.back = gTrack(gr.nseg, y.field = 'mean', col = 'black', y0 = 0, y1 = 1600)

        map = gr.tile.map(ss.n, ncov, verbose = T)
        gc()
        ix = rep(NA, length(ncov))
        ix[unlist(map)] = rep(1:length(map), sapply(map, length))
                                        # ncounts = values(ncov)[, opt$field]

        gc()
                                        #normFactor = ncounts / (1+ss.n$cn[ix]); ## normalize local copy number from normFactor (ie allow 'ratio' to maintain copy number)


        normFactor = values(ncov)[, opt$field] / (ss.n$cn[ix]); ## Inf when ncn is 0
        inf.ix = which(is.infinite(normFactor))
        ## normalize local copy number from normFactor (ie allow 'ratio' to maintain copy number)
        gc()
        tcounts = values(tcov)[, opt$field]
        gc()

        six = sample(setdiff(seq_along(tcov), inf.ix), 1e4)
        reScale = mean(normFactor[six], na.rm = T)/mean(tcounts[six], na.rm = T) ## scalar to get values around 1 if total read depth for tumor != normal
        values(tcov) = data.frame(ratio = tcounts / normFactor * reScale,
                                  tum.counts = values(tcov)[, opt$field],
                                  norm.counts = values(ncov)[, opt$field])



        
        ## normFactor = values(ncov)[, opt$field] / (1+ss.n$cn[ix]); ## normalize local copy number from normFactor (ie allow 'ratio' to maintain copy number)
        ## gc()
        ## tcounts = values(tcov)[, opt$field]
        ## gc()

        ## six = sample(length(tcov), 1e4)
        ## reScale = mean(normFactor[six], na.rm = T)/mean(tcounts[six], na.rm = T) ## scalar to get values around 1 if total read depth for tumor != normal
        ## values(tcov) = data.frame(ratio = tcounts / normFactor * reScale,
        ## tum.counts = values(tcov)[, opt$field], norm.counts = values(ncov)[, opt$field])
        rm(ncov)
        gc()

        cat('removing normal coverage from memory\n')

        ## remove any pieces of tcov not in new.sl (ie chromosomes that have (virtually) 0 read support in normal)

        cat('trimming tcov to final seqlengths\n')

        tcov = gr.fix(tcov, new.sl, drop = T)

        ## save "normal segs" for use downstream
        names(ss.n) = NULL
        saveRDS(ss.n, out.file.nseg.rds)
        gc()
    } else { 
        tcov$ratio = values(tcov)[, opt$field]
        ss.n = NULL
        new.sl = seqlengths(tcov)
        saveRDS(GRanges(), out.file.nseg.rds)
        if (any(is.na(new.sl)))
            {
                library(data.table)
                tmp.sl = data.table(sn = as.character(seqnames(tcov)), end = end(tcov))[, max(end, na.rm = T), by = sn][,  structure(V1, names = sn)]
                new.sl[is.na(new.sl)] = tmp.sl[is.na(new.sl)]
                new.sl = new.sl[!is.na(new.sl)]
            }
    }

cat('saving RDS files\n')
saveRDS(tcov, out.file.cov)

## read hets if supplied
if (file.exists(opt$hets) && (file.size(opt$hets) > 0)) {
    hets.dt = data.table::fread(opt$hets)
    hets.dt = hets.dt[(alt.count.n > 0 | ref.count.n > 0) & (alt.count.t > 0 | ref.count.t > 0),]
    hets.dt[, ":="(alt.frac.t = alt.count.t / (alt.count.t + ref.count.t),
                   ref.frac.t = ref.count.t / (alt.count.t + ref.count.t),
                   alt.frac.n = alt.count.n / (alt.count.n + ref.count.n),
                   ref.frac.n = ref.count.n / (alt.count.n + ref.count.n))]

    hets.dt = hets.dt[(alt.frac.n > opt$hets_thresh & ref.frac.n > opt$hets_thresh),]
    hets.dt[, major.frac := pmax(alt.frac.t, ref.frac.t)]
    hets.dt[, minor.frac := pmin(alt.frac.t, ref.frac.t)]
    hets.dt[, major.count := pmax(alt.count.t, ref.count.t)]
    hets.dt[, minor.count := pmin(alt.count.t, ref.count.t)]
    hets.dt[, BAF := minor.count / (minor.count + major.count)]

    ## create GRanges
    hets.gr = GRanges(seqnames = hets.dt[, seqnames],
                      ranges = IRanges(start = hets.dt[, start],
                                       width = 1),
                      strand = "*",
                      baf = hets.dt[, BAF])
} else {
    hets.gr = GRanges()
}


if (!file.exists(out.file.rds))
    {
        ix = which(!is.na(tcov$ratio))
        cat('sending ', length(ix), ' segments\n')
        cna = CNA(log(tcov$ratio[ix]), as.character(seqnames(tcov))[ix], start(tcov)[ix], data.type = 'logratio')
        gc()
        cat('finished making cna\n')
        seg = segment(smooth.CNA(cna), alpha = opt$cnsignif, verbose = T,
                      undo.splits = opt$undo_splits,
                      undo.prune = opt$undo_prune,
                      undo.SD = opt$undo_SD) ## also 1e-5!!! TODO URGENT
        utils::capture.output({seg_dt = print(seg); setDT(seg_dt)}, type = "output", file = "/dev/null") #### KH
        out = seg2gr(seg_dt[!(is.na(seg.mean) | is.na(loc.start) | is.na(loc.end))], new.sl) ## remove seqlengths that have not been segmented #### KH
        cat('finished segmenting\n')
        if (length(hets.gr)) {
            cat("Segmenting BAF\n")
            hets.cna = CNA(hets.dt[, BAF],
                           hets.dt[, as.character(seqnames)],
                           hets.dt[, start],
                           data.type = "logratio")
            hets.seg = segment(smooth.CNA(cna),
                               alpha = opt$cnsignif,
                               verbose = TRUE)
            utils::capture.output({hets_seg_dt = print(hets.seg); setDT(hets_seg_dt)},
                                  type = "output",
                                  file = "/dev/null")
            hets.out = seg2gr(hets_seg_dt[!(is.na(seg.mean) | is.na(loc.start) | is.na(loc.end))], new.sl)
            hets.out = gr.fix(hets.out, new.sl, drop = T)
            ## keep any breakends that are more than 1e5 from exising start sites
            distances = GenomicRanges::distance(gr.start(hets.out),
                                                gr.start(out),
                                                ignore.strand = TRUE)
            hets.out = hets.out[which(distances > opt$distance_thresh)]
            out = GenomicRanges::disjoin(gUtils::grbind(gUtils::gr.stripstrand(out),
                                                        gUtils::gr.stripstrand(hets.out)))
        }
        ## out = seg2gr(print(seg), new.sl) ## remove seqlengths that have not been segmented
        out = gr.fix(out, new.sl, drop = T)
        cat(length(out), ' segments produced\n')
        names(out) = NULL
        saveRDS(out, out.file.rds)
        write.tab(as.data.frame(out), out.file.txt)
        cat('finished dumping to file, now making plots\n')
        gc()
    } else
    out = readRDS(out.file.rds)

## ix = sample(1:length(tcov), 1e5)
## if (!is.null(ss.n))
##     {
## #        y1.n = mean(sample(tcov$norm.counts, 1e5), na.rm = T) + 3*sd(sample(tcov$norm.counts, 1e5), na.rm = T)
## #        td.tcounts = gTrack(tcov[ix], y.field = 'tum.counts', y1 = mean(sample(tcov$tum.counts, 1e5), na.rm = T) + 3*sd(sample(tcov$tum.counts, 1e5), na.rm = T), name = 'tCov');
## #        td.ncounts = gTrack(tcov[ix], y.field = 'norm.counts', y1 = y1.n, name = 'nCov')
## #        td.nseg = gTrack(ss.n, y.field = 'mean', y1 = y1.n, name = 'n')
## #        td.ncn = gTrack(ss.n, y.field = 'cn', y1 = 5, name = 'nSeg')
##     } else {
##         td.ncounts = td.nseg = td.ncn = td.tcounts = gTrack()
##     }

## td.ratio = gTrack(tcov[ix], y.field = 'ratio', y1 = 10, name = 'covR')
## out$mean = exp(out$seg.mean)
## td.seg = gTrack(out, y.field = 'mean', y1 = 10, name = 'tSeg')

#library("Nozzle.R1")
#cat('making nozzle report\n')

## e = tryCatch({
##     out.file.png = paste(opt$outdir, 'seg.png', sep = '/')
##     png(out.file.png, height = 1500, width = 1500)
##         plot(c(td.ncounts, td.tcounts, td.ratio, td.nseg, td.ncn, td.seg), y0 = 0, border = alpha('black', 0.5))
##     dev.off()
##     report = newReport(paste('CovCBS:', opt$name))
##     fig = newFigure(file.name(out.file.png), fileHighRes = file.name(out.file.png), exportId = "FIGURE_1", "Tracks displaying raw coverage counts for tumor (tCov) and normal (nCov), normalized tumor coverage (CovR), normal tissue integer segmentation (nSeg), and segmented tumor data (tSeg)")
##     report = addToResults(report, fig)
##     writeReport(report, out.file.nozzle)}, error = function(e) 'error')

## if (!is.null(e))
##     {
##         print('failed making with png now trying pdf')
##         tryCatch({
##             out.file.pdf = paste(opt$outdir, 'seg.pdf', sep = '/')
##             pdf(out.file.pdf, height = 1500, width = 1500)
##             plot(c(td.ncounts, td.tcounts, td.ratio, td.nseg, td.ncn, td.seg), y0 = 0, border = alpha('black', 0.5))
##             dev.off()
##             report = newReport(paste('CovCBS:', opt$name))
##             fig = newFigure(file.name(out.file.pdf), fileHighRes = file.name(out.file.pdf), exportId = "FIGURE_1", "Tracks displaying raw coverage counts for tumor (tCov) and normal (nCov), normalized tumor coverage (CovR), normal tissue integer segmentation (nSeg), and segmented tumor data (tSeg)")
##             report = addToResults(report, fig)
##             writeReport(report, out.file.nozzle)}, error = function(e) 'error')
##     }
    

    
cat('done\n')
