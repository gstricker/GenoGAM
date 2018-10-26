## clean up
.onUnload <- function (libpath) {
  library.dynam.unload("GenoGAM", libpath)
}
