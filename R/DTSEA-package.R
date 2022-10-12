#' DTSEA
#' @title The Drug target set enrichment analysis (DTSEA)
#' @name DTSEA-package
#' @aliases DTsEA
#' @description The DTSEA implements a novel application to GSEA and extends the
#' adoption of GSEA.
#
NULL

utils::globalVariables(c("."))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "This package SHOULD NOT BE USED UNDER INTEL MATH KERNEL LIBRARY ON ANY OCCASION.
There is an avoidable but critical bug with Intel Math Kernel Library (MKL) on various operating systems.

======================================
For better performance, we recommend not using RStudio on Windows because RStudio cannot take advantage of the multi-core capabilities available on modern computers. ")
}
