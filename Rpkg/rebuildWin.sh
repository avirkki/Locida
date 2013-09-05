#!/bin/bash

## Refresh and re-compile the Windows version of an R package using wine
##
## Author(s) : Arho Virkki 
## Copyright : VTT Technical Reseach Centre of Finland
## Original date : 2013-08-27

## Instructions to build for Windows under R can be found at:
## http://thomas.zumbrunn.name/blog/2011/12/17/buildung-r-binary-packages-
##        for-windows-with-wine-under-gnulinux/
##
## In summary, all one needs to do under Ubuntu 12.04 is to
##
## 1. Install wine with "sudo apt-get wine"
## 2. Download the Windows version of R and install it with e.g.
##    "wine ~/Downloads/R-3.0.1-win.exe"
## 3. Download zip from ftp://ftp.info-zip.org/pub/infozip/win32/zip300xn.zip
##    and install it.
## 4. Issue "wine regedit" to edit HKEY_CURRENT_USER -> Environment -> PATH 
##    to read "C:\Program Files\R\R-3.0.1\bin\i386;C:\Program Files\infozip"
##    (or similar)
## 5. Build the R packages with "wine R CMD INSTALL --build <package>.tar.gz"



PACKAGE_NAME="Locida"
PACKAGE_VERSION="1.0"

## Delete the old Build
rm $PACKAGE_NAME\_$PACKAGE_VERSION.zip

## Build the Windows version. (Using less shield the Terminal from the 
## harmfull escape characters from Wine).
wine R --vanilla CMD INSTALL --build $PACKAGE_NAME\_$PACKAGE_VERSION.tar.gz | less



