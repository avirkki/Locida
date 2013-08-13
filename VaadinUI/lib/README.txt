The bundled Eclipse .classpath points to this directory to find
the required jar packages.

Rserve libraries are needed for communication between Java and R. The libraries 
are developed by Simon Urbanek and the sources are licensed under GPL.

This script will download the latest libraries from 
http://www.rforge.net/Rserve/files/

We also need to download the RVaadin package by Arho Virkki et al to make the 
Vaadin user interface to work with R. The library is distributed under Apache 2
license, and is available through GitHub (similar to this package) at:
https://github.com/avirkki/RVaadin

To download the files (under Linux), issue:
./getlibs.sh

