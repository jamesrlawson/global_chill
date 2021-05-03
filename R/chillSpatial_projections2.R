require(chillR)
require(sp)
require(tidyverse)
require(raster)
require(lubridate)
require(terra)
require(pbapply)
require(rasterVis)
require(future)
require(future.apply)
require(terra)

source('../../R/MCA_utilities.R')
source('./horticulture/MCA_cherry/R/dev/chillSpatial_functions.R')
# source('chillSpatial_functions.R')

plan(multiprocess, workers=6, gc=TRUE)

models <- tolower(c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"))
scenarios <- c("historical", "ssp126", "ssp585")

historical_GFDL_ESM4 <- getChillWorld(scenario=scenarios[1], model=models[1], year_range=1991:2010)
historical_IPSL_CM6A_LR <- getChillWorld(scenario=scenarios[1], model=models[2], year_range=1991:2010)
historical_MPI_ESM1_2_HR <- getChillWorld(scenario=scenarios[1], model=models[3], year_range=1991:2010)
historical_MRI_ESM2_0 <- getChillWorld(scenario=scenarios[1], model=models[4], year_range=1991:2010)
historical_UKESM1_0_LL <- getChillWorld(scenario=scenarios[1], model=models[5], year_range=1991:2010)

ssp126_GFCDL_ESM4 <- getChillWorld(scenario=scenarios[2], model=models[1], year_range=2041:2060)
ssp126_IPSL_CM6A_LR <- getChillWorld(scenario=scenarios[2], model=models[2], year_range=2041:2060)
ssp126_MPI_ESM1_2_HR <- getChillWorld(scenario=scenarios[2], model=models[3], year_range=2041:2060)
  ssp126_MRI_ESM2_0 <- getChillWorld(scenario=scenarios[2], model=models[4], year_range=2041:2060)
  ssp126_UKESM1_0_LL <- getChillWorld(scenario=scenarios[2], model=models[5], year_range=2041:2060)
  
  ssp585_GFDL_ESM4 <- getChillWorld(scenario=scenarios[3], model=models[1], year_range=2041:2060)
  ssp585_IPSL_CM6A_LR <- getChillWorld(scenario=scenarios[3], model=models[2], year_range=2041:2060)
  ssp585_MPI_ESM1_2_HR <- getChillWorld(scenario=scenarios[3], model=models[3], year_range=2041:2060)
  ssp585_MRI_ESM2_0 <- getChillWorld(scenario=scenarios[3], model=models[4], year_range=2041:2060)
  ssp585_UKESM1_0_LL <- getChillWorld(scenario=scenarios[3], model=models[5], year_range=2041:2060)
  
  ssp126_GFCDL_ESM4 <- getChillWorld(scenario=scenarios[2], model=models[1], year_range=2081:2100)
  ssp126_IPSL_CM6A_LR <- getChillWorld(scenario=scenarios[2], model=models[2], year_range=2081:2100)
  ssp126_MPI_ESM1_2_HR <- getChillWorld(scenario=scenarios[2], model=models[3], year_range=2081:2100)
  ssp126_MRI_ESM2_0 <- getChillWorld(scenario=scenarios[2], model=models[4], year_range=2081:2100)
  ssp126_UKESM1_0_LL <- getChillWorld(scenario=scenarios[2], model=models[5], year_range=2081:2100)
  
  ssp585_GFDL_ESM4 <- getChillWorld(scenario=scenarios[3], model=models[1], year_range=2081:2100) # files not found
  ssp585_IPSL_CM6A_LR <- getChillWorld(scenario=scenarios[3], model=models[2], year_range=2081:2100) # error in firstday:lastday
  ssp585_MPI_ESM1_2_HR <- getChillWorld(scenario=scenarios[3], model=models[3], year_range=2081:2100)
  ssp585_MRI_ESM2_0 <- getChillWorld(scenario=scenarios[3], model=models[4], year_range=2081:2100)
  ssp585_UKESM1_0_LL <- getChillWorld(scenario=scenarios[3], model=models[5], year_range=2081:2100)