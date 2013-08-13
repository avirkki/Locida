## Author(s)     : Arho Virkki
## Copyright     : VTT Technical Reseach Centre of Finland
##
## An example Locida run with a pre-packaged VTT SR1 Nicotiana Tabacum dataset

library(Locida)

## The SR1 data set consists of the sample description data frame 'samples' and
## the data matrix 'X':
data(SR1)
attach(SR1)

dim(X)
colnames(samples)



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
controls <- factors[5:6]

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













