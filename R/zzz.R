#THISPKG <- "GenometriCorr"

.gnmCorr <- new.env(parent=emptyenv())

.onAttach <- function(libname, pkgname) {
       version <- packageDescription("GenometriCorr", fields="Version")
       packageStartupMessage(paste("Welcome to GenometriCorr version ", version))
}

#.onUnload <- function(libpath){
#       library.dynam.unload(THISPKG, libpath)
#}

