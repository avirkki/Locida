## Author(s)     : Arho Virkki
## Copyright     : VTT Technical Reseach Centre of Finland
##
## An example Locida run with a pre-packaged VTT SR1 Nicotiana Tabacum dataset.
## This example is identical to example 1, except that now we read the data
## from an Excel file.

library(Locida)

## Read the data set from Excel file (using the gdata package).
library(gdata)
fname <- "SR1_July2011_calibr_plain.xlsx"
W <- read.xls(fname, check.names=FALSE,
    stringsAsFactors=FALSE)

## This is a convention: Suppose that non-numeric column names correspond to 
## factors, sample names and other descriptive data, and the rest is the 
## spectrum data

suppressWarnings( numericColNames <- as.numeric(colnames(W)) )
factorColumns <- is.na(numericColNames)

## Now, let us separate the data frame into sample description and the
## statistical data matrix X, and we are done.
samples <- W[,factorColumns]
X <- W[,!factorColumns]


## Next, we choose to factor the samples along the GM group descriptions
sampleCategory <- samples[["gm_group_description"]]
factors <- unique(sampleCategory)
print(factors)

## We choose the factors by hand and construct a boolean 
## nrow(samples) x 2 matrix designFrame. The first colums indicates the treated
## and the second row the control samples. Of course, the labels only shown in 
## the plot can be any character strings.
## 
treated <- factors[1]
controls <- factors[5]

designFrame <- data.frame( 
    samples[["gm_group_description"]] %in% treated,
    samples[["gm_group_description"]] %in% controls)

## Set the labels to show in the plot
colnames(designFrame) <- c(treated, paste(controls, collapse="\n"))

## Initialize the locida object
loc <- newLocida(X)
loc$setDesignFrame( designFrame )

## Now, choose the range and plot the image
loc$setXlim( c(10,0) )
loc$setYlim(loc$getMaxYlim())

loc$getPlot()

##
## res <- loc$getDifferenceMatrices( sampleCategory, factors )
## loc$getHeatmapDisplay( res, factors )













