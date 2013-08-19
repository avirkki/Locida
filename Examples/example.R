## Author(s)     : Arho Virkki
## Copyright     : VTT Technical Reseach Centre of Finland
##
## An example Locida run with a pre-packaged VTT SR1 Nicotiana Tabacum dataset

library(Locida)

## The SR1 data set consists of the sample description data frame 'samples' and
## the data matrix 'X':
data(SR1)
attach(SR1)

## Explore the dimensions of the bucketed data and see the sample names
dim(X)
colnames(samples)

## Next, we choose to factor the samples along the GM group descriptions
sampleCategory <- samples[["gm_group_description"]]
factors <- unique(sampleCategory)
print(factors)

## We choose the factors by hand and construct a boolean 
## nrow(samples) x 2 matrix designFrame. The first and second colum indicate the
## treated and the control samples, respectively. The labels are only used to
## label the plot, and can be any character strings.
## 
treated <- factors[1]
controls <- factors[5:6]

designFrame <- data.frame( 
    samples[["gm_group_description"]] %in% treated,
    samples[["gm_group_description"]] %in% controls)

## Set the labels in the plot
colnames(designFrame) <- c(treated, paste(controls, collapse="\n"))

## Initialize the Locida object
loc <- newLocida(X)
loc$setDesignFrame( designFrame )

## For Windows, we suggest using the Cairo package for faster graphics
## library(Cairo)
## CairoWin()

## Now, choose the range and plot the image
loc$setXlim( c(10,0) )
loc$setYlim(loc$getMaxYlim())
loc$getPlot()

## Copy the resulting image into a pdf file
## dev.copy2pdf(file="LocidaTest.pdf")

## res <- loc$getDifferenceMatrices( sampleCategory, factors )
## loc$getHeatmapDisplay( res, factors )













