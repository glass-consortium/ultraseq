

#' @rdname samtools
#' @title set of wrapper functions around samtools
#' 
#' @description several wrappers around samtools
#'
#' @param x input sam/bam file to be sorted and indexed
#' @param samplename sample name
#' @param outbam output bam file
#' @param samtools_exe path to samtools
#' 
#' @export
samtools_sort_index <- function(x, samplename, outbam, samtools_exe){

  cmd_sort = sprintf("%s sort %s -o %s",
    samtools_exe, x, gsub(outbam, ".bam", ""))

  cmd_index = sprintf("%s index %s",
    samtools_exe, outbam)

  cmds = list(sort = cmd_sort, index = cmd_index)
  flowmat = to_flowmat(cmds, samplename)
  return(list(flowmat = flowmat, outfiles = outbam) )

}
