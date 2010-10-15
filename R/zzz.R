".onLoad" <- function (lib, pack) {
    library.dynam(pack, pack, lib)
    packageStartupMessage("RcppDE package")
    packageStartupMessage("  C++ Implementation of Differntial Evolution Optimisation")
    packageStartupMessage("  Author: Dirk Eddelbuettel")
    packageStartupMessage("Extending:")
    packageStartupMessage("  DEoptim package")
    packageStartupMessage("  Differential Evolution algorithm in R")
    packageStartupMessage("  Authors: David Ardia and Katharine Mullen")
}
