## testing WGS analysis
library(ngsflows)
source('~/projects/papers_flow/analysis/fastq_mutect.R')


## ------------------   WGS samples ------ ------ ------ ------  ##
set_opts(fastq_format = "{{samplename}}_{{num}}.fastq.gz",
         .parse = FALSE,
         paired_end = TRUE,
         split_by_chr = "TRUE",
         bwa_method = "mem",
         verbose = 1,
         flow_run_path = "/rsrch2/iacs/iacs_dep/sseth/flowr/wgs/runs")
tumor = "/rsrch2/iacs/iacs_dep/sseth/flowr/wgs/coad/TCGA-A6-6141-01A-11D-1771-10"
normal = "/rsrch2/iacs/iacs_dep/sseth/flowr/wgs/coad/TCGA-A6-6141-10A-01D-1771-10"
deffile = "~/Dropbox/projects/papers_flow/analysis/fastq_mutect_mem.def"

## -------   running a single node mode:
flowmat = fastq_mutect(tumor, normal, 
                       sample1 = "TCGA-A6-6141-01A",
                       sample2 = "TCGA-A6-6141-10A")


flowdef = as.flowdef(deffile)
plot_flow(deffile, pdffile = gsub("def$", "pdf", deffile))
## submission
fobj = to_flow(flowmat, flowdef, 
               flowname = "fastq_mutect", 
               flow_run_path = get_opts("flow_run_path"))


write_sheet(flowmat, file.path(get_opts('flow_run_path'), "flowmat.tsv"))
write_sheet(flowdef, file.path(get_opts('flow_run_path'), "flowdef.tsv"))
#fobj_submit = submit_flow(fobj, execute = FALSE)
#fobj_submit = submit_flow(fobj, execute = TRUE)
