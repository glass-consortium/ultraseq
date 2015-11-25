
## Rscript ~/projects/papers_flow/analysis/run_fastq_haplotyper.R

#library(params)
library(ngsflows)

## for haplotyper
submit <- function(){

  flowmat = rbind(f_merge$flowmat, f_preproc$flowmat, f_hapl$flowmat)
  opts_flow$set(verbose = 2)

  ## submission
  fobj = to_flow(flowmat, flowdef,
                 flowname = "fastq_haplotyper",
                 flow_run_path = "/rsrch2/iacs/iacs_dep/sseth/flowr/runs/illumina_plat")
  opts_flow$set(verbose = 1)
  submit_flow(fobj, execute = TRUE)

}

# this create a flowmat which includes bwa alignment from fastq to bam
fastq_preproc <- function(wd,
                          samplename,
                          fastq_format = get_opts("fastq_format"),
                          jobname_suffix = "_a",
                          split_by_chr = get_opts("split_by_chr")){

  out = create_fq_sheet(x=wd, format = fastq_format)
  out = subset(out, !is.na(num))

  set_opts(samplename = samplename)

  ##   opts_flow$load("~/Dropbox/public/github_ngsflows/inst/conf/ngsflows_mda.conf",
  ##                  check = FALSE)
  #source("~/Dropbox/public/github_ngsflows/inst/pipelines/fastq_bam_bwa.R")
  source(fetch_pipes('fastq_bam_bwa')$pip)

  #fqs1 = head(subset(out, read == 1)$file, 2)
  #fqs2 = head(subset(out, read == 2)$file, 2)

  if('read' %in% colnames(out)){
    fqs1 = subset(out, read == 1)$file
    fqs2 = subset(out, read == 2)$file
    f_merge = fastq_bam_bwa(fqs1, fqs2)
  }else{
    fqs1 = out$file
    f_merge = fastq_bam_bwa(fqs1 = fqs1)
  }

  ## the following steps, are chained
  ## we supply the output of one as
  ## input to the other
  f_preproc = preprocess(f_merge$outfile, split_by_chr = split_by_chr)

  flowmat = rbind(f_merge$flowmat, f_preproc$flowmat)
  flowmat$jobname = paste0(flowmat$jobname, jobname_suffix)

  return(list(flowmat = flowmat, outfile = f_preproc$outfile))
}


fastq_mutect <- function(wd1, wd2, sample1, sample2,
                         split_by_chr = get_opts("split_by_chr")){


  ## bam file names depend on the sample names we supply
  out_a = fastq_preproc(wd1, sample1, jobname_suffix = "_a")
  out_b = fastq_preproc(wd2, sample2, jobname_suffix = "_b")

  paired_name = paste0(sample1, "___", sample2)
  ## -- run mutect on each of the smaller chrs
  mut = mutect(x = out_a$outfile,
               y = out_b$outfile,
               outfile = paired_name,
               samplename = paired_name,
               is_merged = FALSE,
               split_by_chr = split_by_chr)

  ## generate a bigger flowdef with all the elements
  flowmat = rbind(out_a$flowmat, out_b$flowmat, mut$flowmat)

  ## since the three modules are to be submitted as a single pipeline
  ## we will put the same sample name throughout
  flowmat$samplename = paired_name

  #def = to_flowdef(flowmat)
  #write_sheet(def, "analysis/fastq_mutect.def")

  return(flowmat)
}




if(FALSE){

  ## format for NA files
  set_opts(fastq_format = "{{samplename}}_{{read}}.{{num}}.fastq")

  #mypack:::reload('flowr')
  source('~/Dropbox/projects/papers_flow/analysis/fastq_mutect.R')
  ## WEX example: C835.HCC1143_BL.4_2.009.fastq
  ## notice samplename has a _; should confirm that this parses properly
  ## set_opts tries to complete, values in {{}}, which is undesired here.
  ## so turning OFF parsing.
  set_opts(fastq_format = "{{samplename}}_{{read}}.{{num}}.fastq",
           verbose = 1, .parse = FALSE)
  tumor = "/rsrch2/iacs/iacs_dep/sseth/flowr/wex/89324e86-5b7a-4f69-92c1-3b67293f8748"
  normal = "/rsrch2/iacs/iacs_dep/sseth/flowr/wex/ebdb53ae-6386-4bc4-90b1-4f249ff9fcdf"
  set_opts(flow_run_path = "/rsrch2/iacs/iacs_dep/sseth/flowr/wex/runs")


  ## -------   running a single node mode:
  set_opts(split_by_chr = "TRUE") ## off for all
  flowmat = fastq_mutect(tumor, normal,
                      sample1 = "C835.HCC1143.2",
                      sample2 = "C835.HCC1143_BL.4")
  flowdef = as.flowdef("~/Dropbox/projects/papers_flow/analysis/fastq_mutect.def")
  
  ## submission
  fobj = to_flow(flowmat, flowdef,
                 flowname = "fastq_mutect",
                 flow_run_path = get_opts("flow_run_path"))
  fobj_submit = submit_flow(fobj, execute = TRUE)

  ## -------   running a single node mode:
  ## - split_by_chr: off
  ## - previous step finishes completely before starting the next
  ## - steps supporting multicore use 24 cores
  set_opts(split_by_chr = "FALSE") ## off for all
  flowmat = fastq_mutect(tumor, normal,
                         sample1 = "C835.HCC1143.2",
                         sample2 = "C835.HCC1143_BL.4")
  flowdef = as.flowdef("~/Dropbox/projects/papers_flow/analysis/fastq_mutect_node.def")
  ## submission
  fobj = to_flow(flowmat, flowdef,
                 flowname = "fastq_mutect",
                 flow_run_path = get_opts("flow_run_path"))
  fobj_submit = submit_flow(fobj, execute = TRUE)


  #plot_flow(fobj)

  ## rerun
  wd = "/rsrch2/iacs/iacs_dep/sseth/flowr/wex/runs/fastq_mutect-C835.HCC1143.2___C835.HCC1143_BL.4-20150906-03-21-49-s43B8hA1"
  #debug(rerun)
  rerun(x = wd, mat = flowmat, def = flowdef, start_from = "mutect", kill = FALSE, execute = FALSE)

  ## check commands:
  fobj@jobs$mutect@cmds
  fobj_submit@jobs$mutect@cmds
  tail(flowmat$cmd, 1)



  ## in house illumina examples
  base_wd = "/rsrch2/iacs/ngs_runs/150819_SN1222_0289_BC790RACXX/fastqs/Y76I8n1Y76/Project_BreastPDX-MS10batch2"
  tumor = "Sample_MBernstam-BreastPDX-MS10batch2-BCX080P0-T"
  normal = "Sample_MBernstam-BreastPDX-MS10batch2-BCX080Normal"
  set_opts(fastq_format = "{{samplename}}_{{index}}_L00{{lane}}_R{{read}}_{{num}}.fastq.gz")



}
