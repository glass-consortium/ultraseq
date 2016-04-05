#

ultraseq.fq <- function(

){
  # align
  # preprocess
  # ultraseq.bam


  return(flowmat)
}



ultraseq.bam <- function(

){




  return(flowmat)
}

ultraseq_help_text =
"
do this and this etc...

flowr run ultraseq tbam=$tbam nbam=$nbam oprefix=$oprefix

"

ultraseq <- function(

  # options for starting with fastq
  start_from = c("fastq", "bam"),
  samplesheet, # df from create_sample_sheet
  pairs,       # two column table with sample ids

  # options for starting with bam
  tbam,
  nbam,
  out_prefix,

  # options for starting with
  aligner = c("bwa", "bwa.mem"),
  sv_caller = c("mutect"),
  sv_caller = c("lumpy"),
  cnv_caller = c("cnvkit"),

  # subset and check def. Should be compatible with this
  def,


){
  if(start_from == "fastq")
    ultraseq.fq(samplesheet, pairs)

  # incorrect inputs
  if(missing(tbam) | missing(nbam)) & missing(samplesheet))
  message(ultraseq_help_text)

  return(flowmat)
}
