#!/usr/bin/Rscript

##
## Refresh and re-compile an R package
##
## Run with sudo to install the package globally (in Linux, into 
## ‘/usr/local/lib/R/site-library’)

## Author(s)     : Arho Virkki
## Copyright     : VTT Technical Reseach Centre of Finland
## Original date : 2013-08-27

PACKAGE_NAME <- "Locida"
PACKAGE_VERSION <- "1.0"

## Generate the documentation with Roxygen2
library(utils)
library(methods)
library(roxygen2)
roxygenize(PACKAGE_NAME, overwrite=TRUE)

## The actual build
system(paste("R CMD build", PACKAGE_NAME))

PACKAGE_FILE <- paste(PACKAGE_NAME, "_", PACKAGE_VERSION, ".tar.gz", sep="")
system(paste("R CMD INSTALL", PACKAGE_FILE))

