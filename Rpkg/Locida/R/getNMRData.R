## Analysis for Hannu's NMR data for SmartCell 2012

## Author(s)     : Arho Virkki
## Copyright     : VTT Technical Reseach Centre of Finland
## Original date : 2012-04-24


## Open the connection and read the SR1 data The function returns a list of 
## objects, where
##
## 'samples' describe the samples as the corresponding entries from the
##           plant table from the database
##
## 'sampleGroups' is a ready-made grouping of the samples (which is redundant
##                information, and can be re-constructed from  the sample data)
##
## 'X' The standard data matrix with possibly NA entries
## 'U' The polished data matrix with no NA entries
## 
##'@export                 
getData <- function( measurement.label = "SR1_July2011_calibr") {

  cat("Reading data from the SmartCell database...\n")
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, host="localhost", user="sc-forge",
      dbname="smartcell", password="periwinkle098")
  
  ## Anneli's SR1 dataset
  sampleFilter <- paste(
      "measurement.label =", paste("'",measurement.label,"'",sep=""), "AND",
      "plant.label != 'VTT-2011-05-13-SR1-hr-506' AND",
      "plant.label != 'VTT-2011-05-31-SR1-hr-553'")
  
  samples <- dbGetQuery(con, paste(
          "SELECT DISTINCT plant.id, plant.label, species, experiment,",
          "cell_culture, genotype,", 
          "gm_group_description, gm_group_code, gm_experiment_code, gm_comment",
          "FROM plant, result, measurement", 
          "WHERE plant.id = result.plant_id AND",
          "result.measurement_id = measurement.id AND",
          sampleFilter,
          "ORDER BY plant.label;"))
  
  ## Query also the groups 
  #  sampleGroups <- dbGetQuery(con, paste(
  #          "SELECT DISTINCT gm_group_code, gm_group_description, gm_comment",
  #          "FROM plant, result, measurement",
  #          "WHERE plant.id = result.plant_id AND",
  #          "result.measurement_id = measurement.id AND",
  #          sampleFilter,
  #          "ORDER BY gm_group_code;"))
  
  ## Query the NMR spectrum data for this dataset. In this query, we already
  ## know that samples 
  ## 
  ## VTT-2011-05-31-SR1-hr-553 and
  ## VTT-2011-05-13-SR1-hr-506
  
  nmr_spectrum <- dbGetQuery(con, paste(
          "SELECT plant.label, chemical_shift_ppm, value",
          "FROM plant, result, measurement, metadata, nmr_spectrum",
          "WHERE plant.id = result.plant_id AND",
          "result.metadata_id = metadata.id AND",
          "result.measurement_id = measurement.id AND",
          "metadata.id = nmr_spectrum.metadata_id AND",
          sampleFilter,
          "ORDER BY label, chemical_shift_ppm DESC;"))
  
  ## Align the data as an ordinary statistical data matrix 'X'
  ##
  spectrumPoints <- length(nmr_spectrum$label)/length(samples$label)
  X <- matrix(nmr_spectrum$value, ncol=spectrumPoints, byrow=TRUE)
  rownames(X) <- samples$label
  colnames(X) <- nmr_spectrum$chemical_shift_ppm[1:spectrumPoints]
  
  ## Remove the coluns with NA content
  #validColumns <- apply(X, 2, function(x) all(!is.na(x)) )
  #U <- X[,validColumns]
  
  dbDisconnect(con)
  
  #return( list("samples"=samples, "sampleGroups" = sampleGroups, "X"=X, "U"=U))
  return( list("samples"=samples, "X"=X, "factorName"="gm_group_description" ) )
}

