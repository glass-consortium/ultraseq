bam_set <- function(bam, outfile, ref_fasta, split_by_chr){
  
  obj = list()
  
  bam = tools::file_path_as_absolute(bam)
  assert_that(has_extension(bam, "bam"))
  obj$bam = bam
  
  # determine output file name
  if(missing(outfile))
    obj$out_prefix <- gsub(".bam", "", basename(bam))
  else
    obj$out_prefix <- gsub(".bam", "", basename(outfile))
  
  ## if file is available determine whether to split for faster processing
  tmp <- get_fasta_chrs(ref_fasta)
  obj$chr_names = tmp$chrs
  obj$chr_lengths = tmp$lens
  
  if(split_by_chr){
    obj$out_prefix_chr <- paste(obj$out_prefix , obj$chr_names, sep = "_") ## bam names
    obj$gatk_intervals = paste0(" -L ", obj$chr_names)             ## interval files
    
  }else{
    
    obj$out_prefix_chr = obj$out_prefix 
    obj$gatk_intervals = ""
    
  }
  
  
  obj
  
  
}

# single: tumor/normal
#   |__: can create a out-prefix
#
# multiple: tumor/normal
#   |__: assumed one for each chr.
#   |__: length   -  of both files should be same
#   |__: sequence -  of files and those of chrs, is assumed to be the same
#   |__: length of chr and tumor normal files should be the same
#   | 
#   |__: DICT is available. 
#   |__: length and sequence of dict, is assumed (lexicographically sorted)
#   |__: 
#   |__: 
#   |__: 


# create a definitive out prefix
# samplename, may not always be unique, using bams is better
# determine output file name, if out_prefix is not provided


# this is shared, by ALL callers
#' Title
#'
#' @param tumor_bam 
#' @param normal_bam 
#' @param out_prefix 
#' @param is_merged 
#' @param split_by_chr 
#' 
#' @details
#' 
#' \code{
#' # single: tumor/normal
#'   |__: can create a out-prefix
#'
#' multiple: tumor/normal
#'   |__: assumed one for each chr.
#'   |__: length   -  of both files should be same
#'   |__: sequence -  of files and those of chrs, is assumed to be the same
#'   |__: length of chr and tumor normal files should be the same
#'   | 
#'   |__: DICT is available. 
#'   |__: length and sequence of dict, is assumed (lexicographically sorted)
#'   |__: 
#'   |__: 
#'   |__: 
#' create a definitive out prefix
#' samplename, may not always be unique, using bams is better
#' determine output file name, if out_prefix is not provided
#' }
#' 
#' @import assertthat
#'
paired_bam_set <- function(tumor_bam, 
                           normal_bam,
                           out_prefix,
                           is_merged,
                           split_by_chr,
                           ref_fasta = opts_flow$get('ref_fasta')
                           
){
  
  check_args(ignore = "out_prefix")
  
  obj = list()
  
  
  obj$split_by_chr = split_by_chr
  
  obj$out_prefix_chr = ""
  
  # check number of tumor and normal: should be same
  # if number of tumor/normal more than 1, is_merged should be FALSE
  tumor_bam = sapply(tumor_bam, tools::file_path_as_absolute)
  normal_bam = sapply(normal_bam, tools::file_path_as_absolute)

  sapply(sapply(tumor_bam, has_extension, "bam"), assert_that)
  sapply(sapply(normal_bam, has_extension, "bam"), assert_that)
  obj$bam$tumor = tumor_bam
  obj$bam$normal = normal_bam
  
  assert_that(length(tumor_bam) == length(tumor_bam))
  obj$bam$len = length(tumor_bam)
  
  
  if(split_by_chr & is_merged & length(tumor_bam) > 1){
    stop("multiple bams supplied, expected 1; perhaps is_merged should be FALSE?")
    
  }else if(split_by_chr & !is_merged & length(tumor_bam) == 1){
    stop("single bam supplied, expected multiple (one for each chromosome); perhaps is_merged should be TRUE")
  }
  
  
  if(obj$bam$len > 1 & missing(out_prefix))
    stop(">> multiple bams supplied, and out_prefix is missing")
  
  if(obj$bam$len == 1 & missing(out_prefix))
    obj$out_prefix <- paste0(gsub(".bam", "", basename(tumor_bam)), "__", gsub(".bam", "", basename(normal_bam)))
  else
    obj$out_prefix = out_prefix
  
  tmp <- get_fasta_chrs(ref_fasta)
  obj$chrs_names = tmp$chrs
  obj$chrs_lengths = tmp$lens
  
  if(split_by_chr){
    obj$out_prefix_chr <- paste0(obj$out_prefix, "_", obj$chrs_names) # bam names
    obj$gatk_intervals = paste0(" -L ", obj$chr_names)             ## interval files
    
  }else{
    obj$out_prefix_chr = obj$out_prefix 
    obj$gatk_intervals = ""
  }
  
  obj
  
}
