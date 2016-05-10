.onAttach <- function(lib, pkg){
  packageStartupMessage("ultraseq: genomic flows made faster")

  
  conf_dir = c(system.file(package = "ultraseq", "conf"),
          system.file(package = "ultraseq", "inst/conf"))
  conf_dir = conf_dir[file.exists(conf_dir)]
  
  fls = flowr::fetch_conf("*", places = conf_dir)
  
  suppressMessages(flowr::opts_flow$load(fls, check = FALSE))

  if(flowr::opts_flow$get('verbose') > 1)
    packageStartupMessage("\nverbose level: 2 (debug mode)\nfollowing files are being loaded:\n", paste(fls, "\n"))

}
