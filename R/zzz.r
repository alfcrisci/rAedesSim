
.onLoad = function(libname, pkgname) {

library.dynam(pkgname, package=pkgname, lib.loc=.libPaths())
options(digits=10)
}
.onUnload = function(libpath) {
  pkgname <- "rAedesSim";
  library.dynam.unload(chname=pkgname, libpath=libpath)
}
