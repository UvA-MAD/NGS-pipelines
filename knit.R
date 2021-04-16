design.file= snakemake@input[["design_file"]]
figpath = snakemake@params[["knit_figdir"]]
fq.path = snakemake@params[["fastq_store"]]
result.dir=snakemake@params[["result_dir"]]
rmd.file=snakemake@params[["rmd_file"]]
outfile=snakemake@output[["md"]]

#options(useHTTPS=FALSE)
library(knitr)

opts_knit$set(base.dir = result.dir)
knit(rmd.file, output = outfile,quiet=TRUE)

