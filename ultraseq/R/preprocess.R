
#' Flow following Broad's best practices for variant calling, starting from sorted bam
#'
#' @title Pre-process bam files following Broad's best practices for variant calling, starting from aligned BAM file
#' @description This function provides a wrapper around the best practices described on \href{https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1}{GATK's website}.
#' If the link is broken google 'GATK best practices'
#'
#' This aims to perform the following steps ( for DNA ):
#'
#' \itemize{
#' \item mark duplicates
#' \item realign indels
#' \item recalibrate bases
#' \item current version: \emph{3.4-46}
#' }
#'
#' For RNA GATK recommends a additional step of split n trim, which is not currently supported (contributions welcome !).
#'
#' \strong{NOTE}:
#'
#' Some GATK tools use \href{https://www.broadinstitute.org/gatk/guide/article?id=1975}{CPU threads while others use data threads},
#' flowr tries to use efficiently make the best use of both/either depending on tool's compatibility.
#'
#' @param java_exe path to java
#' @param java_tmp path to java tmp, can leave blank
#' @param cpu_markdup not used.
#' @param cpu_target number of threads used for GATK target creation step
#' @param bam 
#' @param samplename 
#' @param split_by_chr 
#' @param gatk_jar 
#' @param picard_dir 
#' @param samtools_exe 
#' @param mem_markdup 
#' @param mem_target 
#' @param cpu_realign 
#' @param ref_fasta 
#' @param gatk_target_opts 
#' @param gatk_realign_opts 
#' @param gatk_baserecalib_opts 
#' @param gatk_printreads_opts 
#' @param outfile 
#' @param bam 
#' @param outfile 
#' @param samplename 
#' @param split_by_chr 
#' @param java_exe 
#' @param java_tmp 
#' @param gatk_jar 
#' @param picard_dir 
#' @param samtools_exe 
#' @param cpu_markdup 
#' @param mem_markdup 
#' @param cpu_target 
#' @param mem_target 
#' @param cpu_realign 
#' @param mem_realign 
#' @param cpu_baserecalib 
#' @param mem_baserecalib 
#' @param cpu_printreads 
#' @param mem_printreads 
#' @param ref_fasta 
#' @param gatk_target_opts 
#' @param gatk_realign_opts 
#' @param gatk_baserecalib_opts 
#' @param gatk_printreads_opts 
#' @param execute_cmds run commands, after creation. Useful for testing/debugging and running on local platforms.
#' @export
#'
#' @examples \dontrun{
#' ## load options, including paths to tools and other parameters
#' load_opts(fetch_conf("ngsflows.conf"), check = FALSE)
#' out = preprocess("my_wex.bam", samplename = "samp", split_by_chr = TRUE)
#'
#' }
preprocess <- function(bam,
                       outfile,
                       samplename = opts_flow$get("samplename"),
                       split_by_chr = opts_flow$get("split_by_chr"),

                       java_exe = opts_flow$get("java_exe"),
                       java_tmp = opts_flow$get("java_tmp"),

                       gatk_jar = opts_flow$get('gatk_jar'),
                       picard_jar = opts_flow$get('picard_jar'),
                       samtools_exe = opts_flow$get('samtools_exe'),

                       cpu_markdup = 1,
                       mem_markdup = opts_flow$get("mem_markdup"),
                       
                       cpu_target = opts_flow$get("cpu_target"),  ## not used
                       mem_target = opts_flow$get("mem_target"),

                       cpu_realign = opts_flow$get("cpu_realign"),
                       mem_realign= opts_flow$get("mem_realign"),

                       ## scatter 8 per node nct=8
                       cpu_baserecalib = opts_flow$get("cpu_baserecalib"),
                       mem_baserecalib = opts_flow$get("mem_baserecalib"),

                       ## scatter 8 per node nct=8
                       cpu_printreads = opts_flow$get("cpu_printreads"),
                       mem_printreads = opts_flow$get("mem_printreads"),

                       ref_fasta = opts_flow$get('ref_fasta'),

                       picard_markdup_opts = opts_flow$get('picard_markdup_opts'),
                       gatk_target_opts = opts_flow$get('gatk_target_opts'),
                       gatk_realign_opts = opts_flow$get('gatk_realign_opts'),
                       gatk_baserecalib_opts = opts_flow$get('gatk_baserecalib_opts'),
                       gatk_printreads_opts = opts_flow$get('gatk_printreads_opts'), 

                       execute_cmds = FALSE
                       ){

  check_args(ignore = "outfile")
  
  bamset = bam_set(bam = bam, outfile = outfile, ref_fasta = ref_fasta, split_by_chr = split_by_chr)

  # get the name of the function
  pipename = match.call()[[1]]
  message("Generating a ", pipename, " flowmat for sample: ", samplename)

  # ------------ dedup; SINGLE FILE
  dedupbam <- paste0(bamset$out_prefix, ".marked.bam")
  metricsfile <- paste0(bamset$out_prefix, ".marked.metrics")
  cmd_markdup <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s %s",
                         java_exe, mem_markdup, java_tmp, picard_jar, bamset$bam,
                         dedupbam, metricsfile, 
                         picard_markdup_opts)
  cmd_markdup

  ## ------------ realign; SINGLE FILE
  intervalsfiles <- paste0(bamset$out_prefix, ".realign.intervals")
  ## ------------ do this for all chrs
  cmd_target <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T RealignerTargetCreator -R %s -I %s -o %s -nt %s %s",
                        java_exe, mem_target, java_tmp, gatk_jar, ref_fasta, dedupbam,
                        intervalsfiles, cpu_target, gatk_target_opts)

  realignedbams <- paste0(bamset$out_prefix_chr ,".realigned.bam")
  cmd_realign <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s %s %s",
                         java_exe, mem_realign, java_tmp, gatk_jar, ref_fasta, dedupbam,
                         intervalsfiles, realignedbams, gatk_realign_opts, bamset$gatk_intervals)

  ## ------------ base recalibration
  ## explicity define intervals here as well, though not required.
  recalibbams <- paste0(bamset$out_prefix_chr, ".recalibed.bam")
  recalibtabfile <- paste0(bamset$out_prefix_chr, ".recalib.tab")
  cmd_baserecalib <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T BaseRecalibrator -R %s -I %s -o %s -nct %s %s %s",
                             java_exe, mem_baserecalib, java_tmp, gatk_jar, ref_fasta,
                             realignedbams, recalibtabfile, cpu_baserecalib,
                             gatk_baserecalib_opts, bamset$gatk_intervals)
  
  cmd_printreads1 <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T PrintReads -R %s -I %s -BQSR %s -o %s -nct %s %s %s",
                            java_exe, mem_printreads, java_tmp, gatk_jar, ref_fasta, realignedbams,
                            recalibtabfile, recalibbams, cpu_printreads,
                            gatk_printreads_opts, bamset$gatk_intervals)
  cmd_printreads2 <- sprintf("%s index %s", samtools_exe, recalibbams)
  cmd_printreads = sprintf("%s;%s", cmd_printreads1, cmd_printreads2)

  cmds <- list(markdup = cmd_markdup,
               target = cmd_target, realign = cmd_realign,
               baserecalib = cmd_baserecalib, printreads = cmd_printreads)
  sapply(cmds, length)
  

  if(execute_cmds)
    sapply(cmds, system)
  
  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(flowmat=flowmat, outfiles = recalibbams))

}

#' @export
get_bam_chrs <- function(x){
  if(file.exists(x)){
    out = Rsamtools::scanBamHeader(x)
    chrs = names(out[[1]]$targets)
  }else{
    message("bam does not exists, returning hg19 chrs")
    chrs = c(1:22,"X","Y","MT")
  }
  return(chrs)
}


#' Read the associated dictionary file and return a list of chromosome names
#'
#' @param fa a reference genome fasta file
#' @param length provide length of each chromosome as well [FALSE]
#'
#' @export
#'
#' @importFrom params read_sheet
#' @importFrom tools file_path_sans_ext
#'
get_fasta_chrs <- function(fa = opts_flow$get("ref_fasta"),
                           length = FALSE){
  check_args()
  
  #dict = gsub("fasta$", "dict", x)
  dict = paste0(file_path_sans_ext(fa), ".dict")
  
  # use default chrs if they are already set
  if(!is.null(opst_flow$get("ref_fasta_chrs")))
    return(list(chrs = opts_flow$get("ref_fasta_chrs"), lens = NA))
    
  if(!file.exists(dict)){
    message(c("We need a dictionary (extension: .dict) for the reference fasta file to proceed. ",
              "Follow this link to learn more: http://lmgtfy.com/?q=create+dict+fasta"))
    warning("By default using hg19 chrs.")
    chrs = c(1:22,"X","Y","MT")
    lens = rep(NA, length(chrs))

  }else{
    seqs = read_sheet(dict, ext = "tsv", skip = 1, head = FALSE)
    chrs = gsub("SN:", "", seqs[, 2], fixed = TRUE)
    lens = gsub("LN:", "", seqs[, 3], fixed = TRUE)
    # 
    # if(length)
    #   return(list(chrs = chrs, lens = lens))
  }
  
  return(list(chrs = chrs, lens = lens))
}

.get_bam_chrs <- function(bam, samtools_exe = "samtools") {
  cmd <- sprintf("%s view -H %s | grep  '\\@SQ' | cut -f 2,3",
                 samtools_exe, bam)
  chrs_info <- system(cmd, intern = TRUE)
  chrs_info <- do.call(rbind, lapply(chrs_info, function(x){
    y = strsplit(x, "\t|:")[[1]]
    y[c(2,4)]
  }))
  return (chrs_info)
}
