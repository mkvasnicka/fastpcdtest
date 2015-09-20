.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to fastpcdtest package.",
                        "Jak citovat.")
}

.onUnload <- function (libpath) {
  library.dynam.unload("fastpcdtest", libpath)
}
