withAutoprint(
    {
        library(optparse)
        options(bitmapType = "cairo")

        options(error = function() {
            traceback(2)
            quit("no", 1)
        })
        if (!exists("opt")) {
            option_list <- list(
                make_option(c("-i", "--id"), type = "character", help = "sample id"),
                make_option(c("-g", "--gGraph"), type = "character", help = "an RDS file contains a gGraph or JaBbA graph with cn annotation on nodes and edges"),
                make_option(c("-r", "--gencode"), type = "character", help = "an RDS or GTF file of GENCODE"),
                make_option(c("-o", "--outdir"), type = "character", default = "./", help = "Directory to dump output into"),
                make_option(c("--cores"), type = "integer", default = 1L, help = "Number of cores")
            )
            parseobj <- OptionParser(option_list = option_list)
            opt <- parse_args(parseobj)
            saveRDS(opt, paste(opt$outdir, "cmd.args.rds", sep = "/"))
        }

        library(gGnome)
        library(gUtils)
        library(parallel) ## needed for mc.cores

        ## setDTthreads(10)
        if (grepl(".rds$", opt$gencode)) {
            gencode <- readRDS(as.character(opt$gencode))
        } else {
            gencode <- rtracklayer::import(opt$gencode)
        }


        ## call complex events
        ## fus = fusions(gG(jab = opt$gGraph), gencode, verbose = TRUE, opt$cores)
        fus <- fusions(gG(jab = opt$gGraph), gencode, verbose = TRUE, mc.cores = opt$cores)

        ## update events with sample id
        if (length(fus)) {
            fus$set(id = opt$id)
            fus$set(mincn = fus$eval(edge = min(cn, na.rm = TRUE)))
        }

        saveRDS(fus, paste0(opt$outdir, "/", "fusions.rds"))

        quit("no", 0)
    },
    echo = FALSE
)
