#!/bin/bash

echo
echo "The bundled Eclipse .classpath points to this directory to find"
echo "the required Rserve jar packages REngine.jar and RserveEngine.jar"
echo
echo "Rserve is needed for the communication between Java and R. The library is"
echo "developed by Simon Urbanek and the sources are licensed under GPL."
echo
echo "This script will download the latest libraries from" 
echo "http://www.rforge.net/Rserve/files/"
echo 
echo "We also need to download the RVaadin package by Arho Virkki et al"
echo "to make the Vaadin user interface to work with R. The library is"
echo "distributed under Apache 2 license, and is available through GitHub"
echo "(similar to this package)."
echo
echo "Press Ctrl+C to cancel, or any other key to proceed"
read

wget -N http://www.rforge.net/Rserve/files/RserveEngine.jar
wget -N http://www.rforge.net/Rserve/files/REngine.jar

wget -N https://github.com/avirkki/RVaadin/blob/master/jar/RVaadin.jar?raw=true
mv "RVaadin.jar?raw=true" RVaadin.jar
