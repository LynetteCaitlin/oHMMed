.onAttach <- function(libname, pkgname) {
  
  v <- getNamespaceVersion("oHMMed")
  packageStartupMessage(paste0("oHMMed ", v, " loaded"))
  packageStartupMessage("For more information please visit: ")
  packageStartupMessage("https://lynettecaitlin.github.io/oHMMed/")
}


.onDetach <- function(libpath) { 
  
  rule <- paste0(rep("=", getOption("width")), collapse = "")
  packageStartupMessage(rule)
  packageStartupMessage("Thank you for using the oHMMed package!")
  packageStartupMessage(rule)
}