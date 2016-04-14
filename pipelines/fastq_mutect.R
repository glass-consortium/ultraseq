
## Rscript ~/projects/papers_flow/analysis/run_fastq_haplotyper.R

#library(params)
library(funr)
library(ultraseq)

# this workflow imports a few sub-workflows
# fetch_pipes would search
# current wd,
# ~/flowr/pipes and
# {{{ultraseq_home}}}/pipes
# if you want you replace this with absolute paths
script_home = detect_home()
source(script_home, "fastq_preprocess.R")
source(script_home, "fastq_bam_bwa.R")

# optionally, find files in flowr pipelines space
#source(fetch_pipes("fastq_preprocess.R")$pip)
#source(fetch_pipes("fastq_bam_bwa.R")$pip)

# in the future this workflow would be renamed and changed
# to include a cnv and trsnslocation caller

#' fastq to mutect
#' either specify
#' @param wd1
#' @param wd2
#' @param fastq_path_tum
#' @param fastq_path_ref
fastq_mutect <- function(
  fqs_tum_1, fqs_tum_2,
  fqs_ref_1, fqs_ref_2,
  fastq_path_tum,
  fastq_path_ref,
  samplename_tum, samplename_ref
  
  ){

  # bam file names depend on the sample names we supply
  #fq_pre = fetch_pipes("fastq_preprocess")
  paired_name = paste0(samplename_tum, "___", samplename_ref)

  preproc_tum = fastq_preprocess(
    fqs1 = fqs_tum_1,
    fqs2 = fqs_tum_2,
    fastq_path = fastq_path_ref,
    samplename =  paired_name)

  preproc_ref = fastq_preprocess(
    fqs1 = fqs_ref_1,
    fqs2 = fqs_ref_1,
    fastq_path = fastq_path_ref,
    samplename = paired_name)

  # add suffix to jobnme
  preproc_tum$flowmat$jobname = paste0(preproc_tum$flowmat$jobname, "_tum")
  preproc_ref$flowmat$jobname = paste0(preproc_ref$flowmat$jobname, "_ref")

  # -- run mutect on each of the smaller chrs
  print(preproc_tum$outfiles$recalibed_bam)

  mut = mutect(tumor_bam = preproc_tum$outfiles$recalibed_bam,
               normal_bam = preproc_ref$outfiles$recalibed_bam,
               outfile = paired_name,
               samplename = paired_name)

  ## generate a bigger flowdef with all the elements
  flowmat = rbind(preproc_tum$flowmat, preproc_ref$flowmat, mut$flowmat)

  ## since the three modules are to be submitted as a single pipeline
  ## we will put the same sample name throughout
  flowmat$samplename = paired_name

  #def = to_flowdef(flowmat)
  #write_sheet(def, "analysis/fastq_mutect.def")

  return(list(flowmat = flowmat, outfiles = mut$outfiles))

}

# ------------------ example with single fastq files ---------------------------

if(FALSE){

 # test example with split_chr ON:
 opts_flow$set(
   split_by_chr = TRUE # false is experimental
   )



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
                 flow_run_path = opts_flow$get("flow_run_path"))
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
                 flow_run_path = opts_flow$get("flow_run_path"))
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
