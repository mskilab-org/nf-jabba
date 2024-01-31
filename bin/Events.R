withAutoprint(
{
    library(optparse)
    options(bitmapType='cairo')
    options(error = function() {traceback(2); quit("no", 1)})

    if (!exists('opt'))
    {
        option_list = list(
            make_option(c("-l", "--libdir"), type = "character", help = "libdir"),
            make_option(c("-i", "--id"), type = "character", help = "sample id"),
            make_option(c("-g", "--gGraph"), type = "character", help = "an RDS file contains a gGraph or JaBbA graph with cn annotation on nodes and edges"),
            make_option(c("-r", "--ref"), type = "character", help = "path to reference file"),
            make_option(c("-o", "--outdir"), type = "character", default = './', help = "Directory to dump output into")
        )
        parseobj = OptionParser(option_list=option_list)
        opt = parse_args(parseobj)
        saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
    }


    ## source("~/lab/home/khadi/git/khscripts/.Rprofile")
    library(gGnome)
    library(gUtils)
    library(skitools)
    library(khtools)
    reduce = GenomicRanges::reduce

    print(.libPaths())
    print(events)

    setDTthreads(1)

    registerS3method("lengths", "gWalk", gGnome::lengths.gWalk, envir = environment())
    registerS3method("merge", "data.table", merge.data.table, envir = environment())


    formals(tic)[c('max.insert', 'min.cushion', 'min.span', 'max.small')] =
        alist(max.insert = 5e5,
              min.cushion = 1e6,
              min.span = 1e6,
              max.small = 1e4)
    newtic = tic
    overwritefun("newtic", "tic", "gGnome")


    ## call complex events
    ## source(paste0(opt$libdir, "/fix.R")); relib3(gGnome); source(paste0(opt$libdir, "/fix.R")) # gUtils gr.val bug?
    gg = gG(jab = opt$gGraph) %>% events(QRP = TRUE)
    
    ## gg = recip_caller(gg)
    ## ev = rbind(gg$meta$events, gg$meta$qrdup, gg$meta$qrdel, fill = T)[
    ##    ,`:=`(ev.id, seq_len(.N))]
    ## gg$set(events = ev)                                  

    ## tryCatch since sometimes get goofy reference mismatch
    gg = tryCatch(microhomology(gg, hg = opt$ref),
                  error = function(e) gg)

    ## update events with sample id
    ev = gg$meta$events[, id := opt$id]
    gg$set(events = ev)

    saveRDS(gg, paste0(opt$outdir, '/', 'complex.rds'))
    
    quit("no", status = 0)
}, echo = FALSE)
