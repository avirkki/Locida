
## Local intensity difference analysis (Locida) for NMR data.
## 
## A method originally developed for the EU SmartCell plant science project
## and with the from the European Union Seventh Framework Programme 
## FP7/2007‐2013 under grant agreement number 222716 – SMARTCELL.

## Author(s)     : Arho Virkki
## Copyright     : VTT Technical Reseach Centre of Finland
## Original date : 2012-04-24


## For analysis, we need
##
## samples - The data frame describing the plant samples.
## X - The standard statistical data matric
## designFrame - A boolean N x g matrix with dim(designFrame) ==  dim(samples),
## where N is the number of samples and g the number of groups (here g = 2).

## The following structure is called closure in the R nomeclature, which 
## closely resembles the Java class structure.


##' Generate a new Locida object as R closure 
##'
##' Blah...
##'
##' @details Details...
##' 
##' \itemize{
##' \item Item 1
##' \item Item 2
##' \item Item 3}
##'
##' @references The format of 
##' \href{http://www.moleculardevices.com/Documents/general-documents/mkt-appnotes/microarray-appnotes/GenePix_Pro_AppNote_Making_GAL_Files_rev_B.pdf}{
##' \emph{Jack Ye Zhai: Making GenePix Array List (GAL) Files}}.
##'
##' @param X A standard statistical data matrix describing the NMR data
##' @param spectLim Optional bucket limits (if not deduceable from X)
##' @return R closure with methods similar to Java/C++ objects 
##' @author Arho Virkki <arho.virkki@@vtt.fi>
##' @export
newLocida <- function( X, spectLim=NULL ) {
  
  
  ## Inspect and set the limits of the spectral data. By convention, the bucket
  ## (column) name equals to its middle point. This computation should yield 
  ## something like 'spectLim = c(9.98, 0.02)' 
  ##
  if( is.null(spectLim) ) {
    spectLim <- as.numeric(colnames(X)[c(1,ncol(X))]) 
    if( any(is.na(spectLim)) ) {
      stop("Cannot deduce the spectral limits from the data matrix")
    }
  }
  
  ## Saturation point for the nonlinear detection method
  I1 <- 0.75
 
  ## Initially empty designframe to compare two groups
  designFrame <- data.frame("Treatment" = rep(FALSE, nrow(X)), 
      "Control" = rep(FALSE, nrow(X)))

  ## Default measure and method:
  measure <- "absolute"
  method <- "saturating"
  
  ## Public methods for measure and computational method
  getMeasures <- function() {
    return( c("absolute", "relative") )
  } 
  getMeasure <- function() {
    return( measure )
  } 
  setMeasure <- function( measure ) {
    measure <<- measure
  }
 
  getMethods <- function() {
    return( c("linear", "cut-off", "saturating") )
  }
  getMethod <- function() {
    return( method )
  }
  setMethod <- function( method ) {
    method <<- method
  }
  
  ## Detection model parameters for the default selection
  d0 <- 0.001
  d1 <- 0.003
  
  ## Detection model parameters for the relative case
  d0rel <- 30
  
  
  ## Statistical detection parameters
  p0 <- 0.05
  p1 <- signif(p0/ncol(X),3)

  
  ## The color scheme
  alpha <- "55"
  blue <- "#0000BB"
  blendedBlue <- "#6391C4"
  red <- "#FF0000"
  cols <- c(red, blue)
  alpha_cols <- paste(c(red, blendedBlue), alpha, sep="") 
  
  ## Unfortunately, R cannot handle reversed coordinates like in the typical
  ## NMR spectrum (10...0). For this reason, we need to set up own our
  ## coordinate system based on the lenght of the data vectors, where k=1...N
  ##
  ## 1                 k                                 N
  ## |-----------------|---------------------------------|
  ## a                 x                                 b
  ##
  ## Where
  ##
  ## (1) x = a + (k-1)h
  ## (2) b = a + (N-1)h   =>   h = (b-a)/(N-1)   "step size"
  ## (3) k = (x-a)/h + 1
  
  N <- ncol(X)
  a <- spectLim[1]
  b <- spectLim[2]
  h <- (b-a)/(N-1)
  
  ## The default xlim (something like c(10,0))
  xlim <- c(spectLim[1]-h/2, spectLim[2]+h/2)
  
  ## designFrame-related member variables
  Ngroups <- NULL
  inDesign <- NULL
  groupAverages <- NULL 
  meanDiff <- NULL
  relSymDiff <- NULL
  relDiff <- NULL
  
  ## By setting the designFrame, we also set several related measures including
  ## Ngroups, inDesign, groupAverages, meanDiff
  ## 
  setDesignFrame <- function( designFrame ) {
    
    ## Change the member variable value
    designFrame <<- designFrame
    
    ## Inspect the number of groups and construct a booleand vector being
    ## true if the sample appears in one of the groups
    Ngroups <<- ncol(designFrame)
    inDesign <<- rowSums(designFrame) > 0
    
    ## Compute the group averages and standard deviations
    groupAverages <<- matrix(nrow=Ngroups, ncol=N)
    
    for( grp in 1:Ngroups ) {
      groupAverages[grp,] <<- colMeans(X[designFrame[,grp],])
    }
    
    ## Choose all the samples that appear in the designFrame and compute the
    ## standard deviation for each column (spectrum bucket) individually.
    ## (not needed at the moment)
    ## 
    ## groupSd <- apply(X[inDesign,], 2, sd) 

    ## In the following, we assume N x 2 designFrame. The extension for
    ## N x 3+ groups needs to be considered separately.
    
    treatedMean <- groupAverages[1,]
    controlMean <- groupAverages[2,]
    
    ## Absolute difference
    meanDiff <<- treatedMean - controlMean

    ## Relative % difference computed in the symmetric way: 
    ## -100% decrease (conventionally 50% decrease) can be compesated  with
    ## +100% increase.
    ##
    relDiff2 <- (treatedMean - controlMean)/controlMean * 100
    relDiff1 <- (controlMean - treatedMean)/treatedMean * 100
  
    relSymDiff <<- relDiff2 * (relDiff2 > 0) - relDiff1 * (relDiff1 > 0)
   
    ## Relative % difference computed in the conventional way:
    ## 50 % decrease can be compensated only with 100% increase.
    relDiff <<- relDiff2
  }
  
  ## Initialize these member variables
  setDesignFrame(designFrame) 
  
  
  getDesignFrame <- function() {
    return( designFrame )
  }
   
  getMeanDiff <- function( onlyForXlim=FALSE ) {
    if( onlyForXlim ) {
      return( meanDiff[getK(xlim[1]):getK(xlim[2])] )
    } else {
      return( meanDiff )
    }
  }
  
  ## Does the treatment increase the signal intensity?
  isIncreased <- function() {
    return( meanDiff > 0 )
  }
  
  ## Compute the index k = 1...N for given frequency x
  getK <- function( x ) {
    return( pmax(1, pmin(N, round( (x - a)/h ) + 1 )) )
  }
  
  ## Compute the frequency x = [a,b] for given index k.
  getFreq <- function( k ) {
    return( a + (k-1)*h )
  }
  
  getMaxYlim <- function() {
    ## Choose all the involved values from X and compute their range
    if(all(inDesign == FALSE)) {
      return(c(0,0.1))
    } else {
      return( range(X[inDesign,
                  getK(xlim[1]):getK(xlim[2])],na.rm=TRUE) )
    }
  }
  ylim <- getMaxYlim()
  
  
  ##
  ## Introduce the member functions
  ##
  
  ## Setters for the graph limits
  setYlim <- function( ylim ) {
    ylim <<- ylim
  }
  setXlim <- function( xlim ) {
    
    ## Sanity check that we are not setting the spectrum wider than
    ## the data range
    xlim[1] <- min(spectLim[1]-h/2, xlim[1])
    xlim[2] <- max(spectLim[2]+h/2, xlim[2])
    
    xlim <<- xlim 
  }

  ##
  ## The Mathematical considerations to compute the color intensity
  ##
  
  ## Approximate the noise from empty signal range 
  ##
  getNoiseSd <- function( xlim=c(10,9.5) ) {
        
    ## Choose those samples which belong to the design, and compute the
    ## combined standard deviation for a range of the spectrum which 
    ## should contain only noise
    combinedData <- as.vector(as.matrix(X[inDesign,getK(xlim[1]):getK(xlim[2])]))
    
    noiseSd <- sd(na.omit(combinedData))
    if( is.null(noiseSd) ) {
      # Return something feasible
      return( 0.001 )
    } else {
      return( noiseSd )
    }
  }
  
  ## Get a feasible estimate for d0 and d1 (using criteria which differ between
  ## different methods)
  ##
  getEstimateD <- function( xlim=c(10,9.5) ) {
    
    if( measure == "absolute" ) {
      
      ## Approximate the noise levels
      noiseLevel <- 100*getNoiseSd( xlim=xlim )
      
      ## Three digits is enought (and looks better in the UI) for 
      ## these ballpark estimates 
      d0 <- signif(noiseLevel,3)
      
      if( method == "cut-off" ) {
        ## Both limits are the same
        d1 <- d0
      }
      
      if( method == "linear" ) {
        ## The largest difference gets 100% color intensity
        d1 <- max( abs(meanDiff), na.rm=TRUE)
      }
      
      if( method == "saturating") {
        
        ## Choose the 10% of greatest difference, compute its mean, and do
        ## not believe too small estimated values.
        NmaxDiff <- ceiling(length(meanDiff)/10)
        dEst <- mean(sort(abs(meanDiff), decreasing=TRUE)[1:NmaxDiff])
        d1 <- signif(max(dEst, 3*d0),3)
      } 
      
    }
    
    if( measure == "relative" ) {
      
      if( method == "cut-off" ) {
        ## 95% quantile as cut-off d0 == d1
        ## qVal <- quantile(abs(relSymDiff), probs=0.95, na.rm=TRUE)
        ## I'd rather keep this fixed instead of d0 <- qVal[1]
        d0 <- d0rel
        d1 <- d0rel
      }      
      
      if( method == "linear" ) {
        ## 95% and 100% quantiles as d0 and d1
        qVal <- quantile(abs(relSymDiff), probs=c(0.95, 1), na.rm=TRUE)
        
        ## I'd rather keep this fixed instead of d0 <- qVal[1]
        d0 <- d0rel
        d1 <- qVal[2]
      }
      
      if( method == "saturating" ) {
        ## These percentual values can be hard-coded, as we are speaking about 
        ## relative saturating method (which adapts to different ranges of 
        ## intesity values as such).
        d0 <- d0rel
        d1 <- 100
      }
    }
    
    return( c(d0,d1) )
  }
  
  ## Setters and getters for the model parameters
  setD <- function( d ) {
    d0 <<- d[1]
    d1 <<- d[2]
    
    ## Fine tuning the UI
    if( measure == "relative") {
      d0rel <<- d[1]
    }
  }
  getD <- function() {
    return( c(d0,d1) ) 
  }
  setP <- function( p ) {
    p0 <<- p[1]
    p1 <<- p[2]
  }
  getP <- function() {
    return( c(p0,p1) ) 
  }
    
  ## The function to compute the relevance of hit in the range [0,1] that
  ## will also be used as the color intensity. If onlyForXlim == TRUE, the 
  ## vector consists only values for the current xlim
  ##
  getColorIntensity <- function( onlyForXlim=FALSE ) {
    
    if( onlyForXlim ) {
      klim <- getK(xlim)
    } else {
      klim <- c(1,N)
    }
    
    ## In case of 2 groups, do Wilcoxon group comparison for them
    if( Ngroups == 2 ) {
      
      ## Suppose that 1st column is treated and 2nd is for the controls
      Xt <- X[designFrame[,1],]
      Xc <- X[designFrame[,2],]
      
      statistic <- rep(as.numeric(NA),N)
      p.value <- rep(as.numeric(NA),N)
      
      for( k in klim[1]:klim[2] ) {
        
        ## Both vectors must contain non-NA values
        if( !all(is.na(Xc[,k])) && !all(is.na(Xt[,k])) ) {
          
          res <- wilcox.test(Xc[,k],Xt[,k], exact=FALSE)
          statistic[k] <- res$statistic
          p.value[k] <- res$p.value 
        }     
      }
      
      ## Construct a rank based on the Wilcox statistic and two threshold
      ## p-values
      
      ## Compute the logical threshold vector for the lower limit p0
      pTresh <- p.value < p0
      
      ## Use logarithmic model scaled by the upper threshold p1 limited by 1
      pLevel <- pmin(-log(p.value)/-log(p1),1)*pTresh
     
      ##
      ## Detection measure: "absolute" or "relative"
      ##
      
      if( measure == "absolute" ) {
        d <- abs(meanDiff)  
      
      } else if( measure == "relative" ) {
        d <- abs(relSymDiff)
      
      } else {
        stop("Invalid measure")
      }
      
      ##
      ## Detection methods: "linear", "cut-off", "saturating"
      ##
      
      if( method == "cut-off" ) {

        ## Just a binary {0,1} response
        Id <- as.numeric( d > d0 )
        
        
      } else if( method == "linear" ) {
        
        ## A ramp function 
        ##
        ## I(d) := 0 for d <= d0, 1 for d > d1, and linear increase
        ##         (d-d0)/(d1-d0) between d0 and d1.
        ##   
        Id <- (d > d0) * pmin( (d-d0)/(d1-d0), 1 )
        
      } else if( method == "saturating" ) {
          
        ## Use the following non-linear index for the intensity change
        ##
        ## I(d) := k(d-d0) / ( k(d-d0) + 1 ) for  d > d0, otherwise I(d) := 0
        ##
        ## where k = 1/(d1-d0) * I1 / (1-I1)  can be computed from the upper
        ## limit d1, where I1 is a matter of definition.
     
        k <- 1/(d1-d0)*I1/(1-I1)
        Id <- k*(d-d0)/(k*(d-d0) + 1) * (d>d0)
      
      } else {
        stop("Invalid method")
      }
      
      return( list("intensity"=(pLevel*Id)[klim[1]:klim[2]],
              "p.value"=p.value[klim[1]:klim[2]]) )
    }
  }
  
  getDetectionCurve <- function() {
    
    ## Store the graphics parameters
    old.par <- par(c("mar","xpd","cex"))
    par("mar"=c(5.1,4.1,3.1,1.1))
    par("xpd"=NA)
    par("cex"=1.0)
  
    ## Illustrate all the 6 different detection modes here
    
    if( method == "linear" ) {
      ## Ramp function
      plot( x=c(0,d0,d1,1.1*d1), y=c(0,0,100,100), type="l", bty="L",
            xlim=c(0,1.1*d1), ylim=c(0,100), xlab="", ylab="")
        
      lines(x=c(d0,d0), y=c(-5,105), lty="dashed")
      lines(x=c(d1,d1), y=c(-5,105 ), lty="dashed")
      text(x=d0, y=112, labels="d0")
      text(x=d1, y=112, labels="d1")
    }
    
    if( method == "cut-off" ) {
      ## Step function
      plot( x=c(0,d0,d0,2*d0), y=c(0,0,100,100), type="l", bty="L",
          xlim=c(0,1.5*d1), ylim=c(0,100),
          xlab="", ylab="")
      
      lines(x=c(d0,d0), y=c(100,105), lty="dashed")
      text(x=d0, y=112, labels="d0")
    }
    
    if( method == "saturating" ) {

      ## Saturation curve
      k <- 1/(d1-d0)*I1/(1-I1)

      x <- c(0, seq(d0, 1.5*d1, len=100))
      y <- k*(x-d0)/(k*(x-d0) + 1) * (x>d0) * 100
      
      plot( x, y, type="l", bty="L", xlim=c(0,1.5*d1), ylim=c(0,100),
          xlab="", ylab="")
      
      lines(x=c(d0,d0), y=c(-5,105), lty="dashed")
      lines(x=c(d1,d1), y=c(-5,105 ), lty="dashed")
      lines(x=c(0, 1.5*d1), y=c(100,100), lty="dashed" )
      lines(x=c(0.85,1.5)*d1,y=c(I1*100, I1*100), lty="dashed")
      
      text(x=d0, y=112, labels="d0")
      text(x=d1, y=112, labels="d1")
      text(x=0.72*d1, y=I1*100, labels=paste(as.character(I1*100),"%"))
    }
    
     
    title( ylab="Color intensity (%)" )
    if( measure == "absolute" ) {
      title(xlab="Difference of Group Means (Intensity)")
    }
    if( measure == "relative" ) {
      title(xlab="Difference of Group Means (%)")
    }  
    
    
    
    ## Restore the default graphics parameters
    par(old.par)
  }
  

  ## Summarize the hits (i.e. the detected regions of the spectrum) into a 
  ## data frame to be scrutinized later
  ##
  getResultFrame <- function( onlyForXlim = TRUE ) {

    ## Index range for the current xlim
    kRange <- seq(getK(xlim[1]), getK(xlim[2]))

    ## Axis for the chemical shift
    csRange <- seq( from=spectLim[1], to=spectLim[2], length.out=ncol(X))
    cs <- csRange[kRange]
    
    ## Intensity and P values
    cInt <- getColorIntensity( onlyForXlim )
    intensityVec <- cInt$intensity
    intensityVec[is.na(intensityVec)] <- 0
    p.value <- cInt$p.value

    ## This detection is based on the concept of increasing and decreasing
    ## signal slopes arising e.g. in analog electronics
    hit <- intensityVec != 0
    hitStartIdx <- which(diff(c(FALSE, hit)) == 1)
    hitEndIdx <- which(diff(c(hit,FALSE)) == -1)
    
    nHits <- length(hitStartIdx)
    
    resFrame <- data.frame( 
        "Start"=numeric(nHits),
        "End"=numeric(nHits),
        "RelChange (%)"=numeric(nHits),
        "RelSymChange"=numeric(nHits),
        "AbsChange"=numeric(nHits), 
        "P"=numeric(nHits), check.names=FALSE)
    
    
    ## Now all the indices are known, and we can start slicing the results. 
    ## For simplicity, if the hit consists of multiple buckets, only the maximal
    ## difference and minimal p-value are recorded (usually co-inciding to the 
    ## same bucket)
    ##
    for( i in seq_along(hitStartIdx) ) {
      
      ## Record the hit statistics into the resFrame
      resFrame[i,"Start"] = cs[hitStartIdx[i]]
      resFrame[i,"End"] = cs[hitEndIdx[i]]
  
      resFrame[i,"RelChange (%)"] = max(relDiff[kRange]
              [hitStartIdx[i]:hitEndIdx[i]])
      
      resFrame[i,"RelSymChange"] = max(relSymDiff[kRange]
              [hitStartIdx[i]:hitEndIdx[i]])

      resFrame[i,"AbsChange"] = max(meanDiff[kRange]
              [hitStartIdx[i]:hitEndIdx[i]])
      
      resFrame[i,"P"] = min(p.value[hitStartIdx[i]:hitEndIdx[i]])
    }
    
    return( resFrame )
  }

  
  ## Computes the difference matrices for the curren xlim, and 
  ## given factorVec and factor classes 'factorClass'.
  ##
  ## For example, viable choices are
  ## factorVec <- samples$gm_group_description
  ## factorClass <- unique(samples$gm_group_description)
  
  getDifferenceMatrices <- function(factorVec, factorClass) {
    
    Nfac <- length(factorClass)
    
    ## Two matrices, other for the increased signal because of the treatment 
    ## (redM) and the other vice versa (blueM)
    ##
    redM <- matrix(rep(0, Nfac^2), nrow=Nfac)
    blueM <- matrix(rep(0, Nfac^2), nrow=Nfac)
    
    eps <- sqrt(.Machine$double.eps)
    
    ## Save the old design
    dfOld <- getDesignFrame()
    
    
    ## Of course, we need Nfac >= 2:
    
    for( i in 1:Nfac ) {
      samples_i <- factorVec %in% factorClass[i]
      for( j in setdiff(1:Nfac,i) ) {
        samples_j <- factorVec %in% factorClass[j]
        ## Set the corresponding designFrame
        designFrame <- data.frame(
            "T" = samples_i,
            "C" = samples_j )
        
        colnames(designFrame) <- c(as.character(i),as.character(j))
        
        setDesignFrame( designFrame )
        
        intensity <- getColorIntensity( onlyForXlim=TRUE )$intensity
        intensity[is.na(intensity)] <- 0
        
        d <- getMeanDiff( onlyForXlim=TRUE )
        
        ## d > 0  ==> Treatment increases signal (which should be show as *red* 
        ## according to our convention)
        
        hits <- (intensity > eps )
        
        ## NA values do not decrease nor increase:
        doesIncrease <- (d > 0)
        doesIncrease[is.na(doesIncrease)] <- FALSE
        
        increasingHits <- intensity[doesIncrease & hits]
        if( length(increasingHits) > 0 ) {
          redM[i,j] <- sum(increasingHits)
        }
      }
    }
    
    ## Put the old frame back
    setDesignFrame( dfOld )
    
    maxVal <- max(redM,blueM)
    TM <- redM/maxVal
    
    return( list("TM"=TM, "CM"=t(TM)) )
  } 
  
  
  ## The plot routine which produces the graph to the default device
  ##
  getPlot <- function() {
    
    ## Ask for the intensities
    colorIntensity <- getColorIntensity()$intensity
    
    ## Store the graphics parameters
    old.par <- par(c("mar","xpd"))
    par("mar"=c(6.1,5.1,2.1,0.6))
    par("xpd"=NA)
    
    ## Open the plot window   
    plot.new()
    klim <- getK(xlim)
    plot.window( xlim=c(klim[1]-1,klim[2]+1), ylim=ylim )
    
    ## Draw the Axis and labels
    
    ## A simple method to determine the tic distance between xlim[1] > xlim[2]
    ## on the x-axis (abscissa)
    ticStep <- 1
    scaleUp <- 5
    tmp <- xlim
    while( (tmp[1] - tmp[2]) < 5 ) {
      tmp <- tmp*scaleUp
      ticStep <- ticStep/scaleUp
    }
    
    labels <- seq(round(tmp[1])*ticStep, xlim[2], by=-ticStep)
    ##labels <- seq(round(tmp[1])*ticStep, by=-ticStep)
    axis(1, at=getK(labels), labels=labels )
    

    
    
    ## The default routine is fine for the y-axis (ordinate)
    axis(2)
    title(xlab="Chemical shift (ppm)", ylab="Signal intesity")
    
    ## The x range for all signals
    xIdx <- klim[1]:klim[2]
         
    notNA <- !is.na(colorIntensity)
    shadeCol <- rep(as.numeric(NA),N)
    
    increased <- isIncreased()
    
    ## If treatment increases signal, it should appear red
    shadeCol[notNA & increased] <- rgb(1,0.5,0.5, 
        alpha=colorIntensity[notNA & increased])
    
    ## Decreasement appears as blue
    shadeCol[notNA & !increased] <- rgb(0.5,0.5,1, 
        alpha=colorIntensity[notNA & !increased])
    
    
    for( k in xIdx[!is.na(xIdx)] ) {
      bucketSpan <- c(k-0.5,k+0.5)
      
      polygon(border=NA,
          x=c(bucketSpan[1],bucketSpan[1],bucketSpan[2], bucketSpan[2]),
          y=c(ylim[1],ylim[2],ylim[2],ylim[1]), col=shadeCol[k])
    }
    
    ## The samples to be included are given by the columns of the design frame
    for( grp in 1:Ngroups ) {
      
      ## Draw each signal individually:   
      groupIdx <- which(designFrame[grp] == TRUE)
      for( sgIdx in groupIdx ) {
        
        ## This nested loop *is* heavy!
        #for( i in xIdx ) {
        #  bucketSpan <- c(i-0.4,i+0.4)
        #  lines(bucketSpan, rep(X[sgIdx,i],2), col=alpha_cols[grp], lwd=2)
        #}
        
        ## To avoid the loop, introduce a heuristic measure to choose correct 
        ## character expansion 'cex'!
        A <- (xlim[1]-xlim[2])
        cex <- 11/A
        
        lines(xIdx, X[sgIdx,xIdx], col=alpha_cols[grp], lwd=1, type="p", pch="-", 
            cex=cex)  
      }
    }  
    
    for( grp in 1:Ngroups ) {
      ## Draw the average spectrum for the group (on top of all graphics)
      lines(xIdx, groupAverages[grp,xIdx], col=cols[grp], lwd=1, type="l", pch=1)
      points(xIdx, groupAverages[grp,xIdx], col="black", lwd=1, pch="-", cex=1)
    }
    
    ## Construct labels for the image
    ##
    yStep <- -diff(ylim)/10
    for( grp in 1:Ngroups ) {
      text(x=klim[1], y=ylim[2]+(grp-1)*yStep, col=cols[grp],
          labels=colnames(designFrame)[grp], pos=4, cex=0.9)
    }
    
    ## Restore the default graphics parameters
    par(old.par)
  }
  
  ## Helper function to produce an individual difference matrix from 
  ## difference matrices T and C (computed by 'getDifferenceMatrices')
  ## 
  getHeatmap <- function( T, C=NULL, rgbCode=c(1,0,0), rgbCodeC=c(0,0,1), ... ) {
    
    par.old <- par("mgp")
    N <- ncol(T)
    ## Open the plot window   
    plot.new()
    plot.window( xlim=c(0,N), ylim=c(0,N) )
    
    ## Draw the Axis and labels
    
    labels <- as.character(1:N)
    for(k in c(1,3)) {
      axis(k, at=seq(0.5,N-0.5), labels=labels, lwd=0, line=-1 )
    }
    for(k in c(2,4)) {
      axis(k, at=seq(0.5,N-0.5), labels=rev(labels), lwd=0, line=-1 )
    }
    
    par("mgp"=c(1.5,0.8,0))
    title(ylab="Treated", xlab="Control", ... )
    
    for( i in 1:N ) {
      for( j in 1:N ) {
        if( is.null(C) ) {
          polygon(border=NA,
              y=N-c(i-1,i-1,i,i), x=c(j-1,j,j,j-1),
              col=rgb(rgbCode[1],rgbCode[2],rgbCode[3],T[i,j]))
          
        } else {
          polygon(border=NA,
              y=N-c(i-1,i-1,i,i), x=c(j-1,j,j,j-1),
              col=rgb(rgbCode[1],rgbCode[2],rgbCode[3],max(T[i,j]-C[i,j],0)))
          polygon(border=NA,
              y=N-c(i-1,i-1,i,i), x=c(j-1,j,j,j-1),
              col=rgb(rgbCodeC[1],rgbCodeC[2],rgbCodeC[3],max(C[i,j]-T[i,j],0)))
        }
      }
    }
    
    par("mgp"=par.old)
  }
  
  ## Produces a three-image heat map display from the output of 
  ## 'getDifferenceMatrices', res$TM and res$CM. 
  ## 
  getHeatmapDisplay <- function( res, factors ) {
    layout(matrix(c(1,2,4,3), 2, 2, byrow = TRUE), c(1,1), c(1,1))
    getHeatmap( res$TM,rgbCode=c(1,0,0),main="Intensity increase")
    getHeatmap( res$CM,rgbCode=c(0,0,1),main="Intensity decrease")
    getHeatmap( res$TM, res$CM,
        rgbCode=c(1,0,0), rgbCodeC=c(0,0,1), main="Intesity change")
    
    plot.new()
    par(xpd=NA)
    plot.window( xlim=c(0,1), ylim=c(0,1) )
    
    legendText <- 
        paste(paste( seq_along(factors),"=", factors ), collapse="\n")
    text(x=-0.2, y=0.3, labels=legendText, pos=4)
    
    ## Restore default layout
    layout(1)
  }
  
  ## A file name for the resultFrame
  ## 
  getFileName <- function() {
    return( paste(measure, '-', method, '_d0=', d0, '_d1=', d1,'_p0=', 
            p0, '_p1=',p1, '.csv', sep='') )
  }
  
  ## Return the public methods
  return( list( "setXlim"= setXlim, "setYlim"= setYlim, 
          "getMaxYlim"=getMaxYlim, 
          "setDesignFrame"=setDesignFrame, "getDesignFrame"=getDesignFrame,
          "getMeasures"=getMeasures, "getMeasure"=getMeasure, 
          "setMeasure"=setMeasure,
          "getMethods"=getMethods, "getMethod"=getMethod, "setMethod"=setMethod,
          "getPlot"=getPlot, 
          "getColorIntensity"=getColorIntensity, "getMeanDiff"=getMeanDiff,
          "isIncreased"=isIncreased, "getNoiseSd"=getNoiseSd,
          "getEstimateD"=getEstimateD, 
          "setD"=setD, "getD"=getD, "setP"=setP, "getP"=getP,
          "getDifferenceMatrices"=getDifferenceMatrices,
          "getHeatmap"=getHeatmap, "getHeatmapDisplay"=getHeatmapDisplay,
          "getResultFrame"=getResultFrame, 
          "getDetectionCurve"=getDetectionCurve,
          "getFileName"=getFileName) )
}
