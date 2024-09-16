run_occtest <- function(occ_species, output_env, species) {
  names(occ_species)[4] <- 'eventDate'
  default_settings <- defaultSettings()
  my_analysis_settings <- default_settings$analysisSettings
  
  # Human and land use detection
  my_analysis_settings$humanDetection$doHumanDetection <- FALSE
  my_analysis_settings$landUse$doLandUse <- FALSE
  my_analysis_settings$institutionLocality$doInstitutionLocality <- T # for natural points
  
  # Country and range analysis
  my_analysis_settings$centroidDetection$doCentroidDetection <- FALSE
  my_analysis_settings$countryStatusRange$doCountryRecordAnalysis <- FALSE
  my_analysis_settings$countryStatusRange$doRangeAnalysis <- FALSE
  
  # Time accuracy test
  my_analysis_settings$doTimeAccuracy <- FALSE
  
  # Environmental outliers test
  my_analysis_settings$envOutliers$doEnvOutliers <- FALSE
  
  # geographical outliers
  myGeoOutliers <- list(
    doGeoOutliers = TRUE,
    methodGeoOutliers = c('alphaHull', 'distance', 'median', 'grubbs', 'dixon', 'rosner', 'mcp'),
    alpha.parameter = 2,
    mcp_percSample = 95
  )
  
  my_analysis_settings$geoOutliers <- myGeoOutliers
  my_analysis_settings$geoOutliers$doGeoOutliers <- FALSE  # whether to use geographical outliers
  
  occ_species_test <- occTest(
    sp.name = species,
    sp.table = occ_species,
    r.env = output_env,
    analysisSettings = my_analysis_settings,
    verbose = TRUE
  )
  clean_occ <- occFilter(df = occ_species_test, errorAcceptance = 'strict')
  # if (errorAcceptance == 'strict')    {colMetaData$errorThreshold =0.2}
  # if (errorAcceptance == 'majority')  {colMetaData$errorThreshold =0.5}
  # if (errorAcceptance == 'relaxed')   {colMetaData$errorThreshold =0.7}
  list(
    occ_species_test = occ_species_test,
    clean_occ = clean_occ
  )
}