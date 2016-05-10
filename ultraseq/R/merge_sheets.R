#' merge tables by row, and filter rows
#'
#' @param x a character vector of file to be merged
#' @param outfile path to the merged output file
#' @param .filter a filter string
#' @param ... other arguments supplied to params::read_sheet
#' 
#' @importFrom dplyr bind_rows
#' @importFrom params write_sheet
#'
#' @export
#' 
#' @examples \dontrun{
#' merge_sheets(c("mutect.chr1.txt", "mutect.chr2.txt"), 
#'    outfile = "mutect.merged.txt", .filter = "judgement==KEEP")
#' }
#' 
merge_sheets <- function(x, outfile, .filter = NA, ...){
  tmp <- lapply(x, function(fl){
    message(".", appendLF = FALSE)
    tab = read_sheet(fl, ...)
    
    # convert all columns to character, fool proofing
    tab[] <- lapply(tab, as.character)
    dim(tab)

    # filtering each file, before merging
    if(!is.na(.filter)){
      tab2 = dplyr::filter_(tab, .filter)
      return(tab2)
    }else{
      return(tab)
    }
  })
  
  # mrgd = try(do.call(rbind, tmp))
  # if fails try using dplyr
  mrgd = try(bind_rows(tmp))
  
  if(!missing(outfile))
    write_sheet(mrgd, outfile)
  
  invisible(mrgd)
}
