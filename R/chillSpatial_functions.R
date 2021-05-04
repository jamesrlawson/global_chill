# MHT1 interpolates hourly temperatures from daily tmin and tmax

MHT1 <- function (latitude, tmin=tmin, tmax=tmax, dates, keep_sunrise_sunset = FALSE, ...) 
{  
  
  this.year = unique(dates$Year)[1]
  
  message(paste0('Calculating hourly temperatures for year: ', this.year))
  
  datcells <- which(!is.na(tmin[[1]][[1]][]))
  
  message('  Loading temperature data into memory..')
  
  lats <- coordinates(raster(tmin))[,2]
  
  tmax <- as.matrix(tmax)
  tmin <- as.matrix(tmin)
  
  JDay <- dates$JDay
  
  if(max(JDay)>365) {
    JDay[JDay>365] <- JDay[JDay>365] -365
  }
  

  message('  Processing daylengths..')
  
  hourly_temps.list <- vector('list', length(datcells))
  pb <- txtProgressBar(max=length(datcells), style=3)
  for(n in datcells) {
    setTxtProgressBar(pb, n)

    dates.cell <- data.frame(Year = dates$Year,
                             JDay = JDay,
                             Tmax = tmax[n,],
                             Tmin = tmin[n,])

    hourly_temps.list[[n]] <- round(chillR::make_hourly_temps(lats[n], dates.cell),2)

  }

  message("  Reticulating splines..")
  # remove null elements of list
  hourly_temps.list[sapply(hourly_temps.list, is.null)] <- NULL
  return(hourly_temps.list)
}

#### calculate hourly temps across a series of years ####
# wrapper for MHT()
# years is a vector of years corresponding to the list objects of tmin and tmax data
# daylight is calculated by DH()
# JDay is a vector of julian days describing the dormancy period of the plant
# writeToDisk is a string describing an output directory where hourly_temps files can be written to disk
  # this will prevent these outputs using too much system memory

# MHT.years <- function(years, lat, JDay, tmin, tmax, writeToDisk=FALSE,...) {
#   
#   message('Annualizing temperature data..')
#   
#   if(class(tmin)=='list') {
#     
#     MHThourly_temps <- future_lapply(seq_along(years), function(x) MHT(lat,
#                                                                     tmin=tmin[[x]],
#                                                                     tmax=tmax[[x]],
#                                                                     dates = data.frame(Year = years[x],
#                                                                                        JDay = JDay)))
#     
#   } else {
#     
#     # should check that tmin and tmax represent the same date range (or are the same dim's)
#     
#     # subset input data to dormancy JDay range
#     # change this to convert to SpatRaster if input data is in raster:: format
#     if(class(tmin[[1]]) %in% c("RasterBrick", "RasterStack")) {
#       dates <- as.Date(1:nlayers(tmin), origin=as.Date(paste0(years[1], "-01-01"))-1)
#     } else {
#       if(class(tmin[[1]]) %in% c("SpatRaster")) {
#         dates <- as.Date(1:nlyr(tmin), origin=as.Date(paste0(years[1], "-01-01"))-1)
#       }
#     }
#     
#     dates.year <- lubridate::year(dates)
#     
#     tmin.out <- vector("list", length(unique(dates.year)))
#     tmax.out <- vector("list", length(unique(dates.year)))
#     period.dates <- vector("list", length(unique(dates.year)))
#     
#     # if(max(JDay)>365) {
#     yearRange <- unique(dates.year)[1:(length(unique(dates.year))-1)]
#     # } else {
#     # yearRange <- unique(dates.year)
#     # }
#     
#     for(i in seq_along(yearRange)) {
#       firstDay <- which(dates.year %in% years[i])[JDay][1]
#       lastDay <- firstDay + length(JDay)-1
#       period <- firstDay:lastDay
#       period.dates[[i]] <- dates[period]
#       tmin.out[[i]] <- tmin[[period]]
#       tmax.out[[i]] <- tmax[[period]]
#     }
# 
#     if(!writeToDisk) {
#       MHThourly_temps <- lapply(seq_along(yearRange), function(x) MHT(lat,
#                                                                       tmin=tmin.out[[x]],
#                                                                       tmax=tmax.out[[x]],
#                                                                       dates = data.frame(Year = year(period.dates[[x]]),
#                                                                                          JDay = yday(period.dates[[x]]))))
#       return(MHThourly_temps)
#       
#     } else {
#       
#       lapply(seq_along(yearRange), function(x) MHT(lat,
#                                                    tmin=tmin.out[[x]],
#                                                    tmax=tmax.out[[x]],
#                                                    dates = data.frame(Year = year(period.dates[[x]]),
#                                                                       JDay = yday(period.dates[[x]]))))
#       return(NULL)
#     }
#     
#   }
#   
#   
#   
# }


#### calculate chill portions  ####

# takes output from MHT.years as input and calculates chill portions values for a specified raster cell
# sits inside wrapper chill_spatial.years which runs chill_spatial for each year of data

chill_spatial <- function(cell, hourly_temps, JDay) {
  out <- chillR::chilling(chillR::stack_hourly_temps(cell), Start_JDay=JDay[1], End_JDay=JDay[length(JDay)])$Chill_portions
  return(out)
}

# wrapper for chill_spatial to run across multuple years
chill_spatial.years <- function(hourly_temps, template, readFromDisk=FALSE, JDay=NULL) {

  # run chill spatial on hourly_temps list object
    out <- future_lapply(hourly_temps, chill_spatial, JDay=JDay)
    
    # take output values and assign to a template raster
    out.unlist <- unlist(out)
    chill.ras <- template
    cells <- which(!is.na(chill.ras[]))
    chill.ras[cells] <- out.unlist
    return(chill.ras)
}


# parent function for calculating hourly temperatures and then chill portions
getChillSpatial <- function(years, lat, JDay, tmin, tmax, template, writeToDisk=FALSE,...) {
  message('Annualizing temperature data..')
  # should check that tmin and tmax represent the same date range (or are the same dim's)
  
  # subset input data to dormancy JDay range
  # change this to convert to SpatRaster if input data is in raster:: format
  
  if(class(tmin)=='list') {
    
    # interpolate hourly temperatures  

    MHThourly_temps <- future_lapply(seq_along(years), function(x) MHT1(lat,
                                                                tmin=tmin[[x]],
                                                                tmax=tmax[[x]],
                                                                dates = data.frame(Year = years[x],
                                                                                   JDay = JDay)))
    
  } else {
      if(class(tmin[[1]]) %in% c("RasterBrick", "RasterStack", "RasterLayer")) {
        dates <- as.Date(1:nlayers(tmin), origin=as.Date(paste0(years[1], "-01-01"))-1)
      } else {
        if(class(tmin[[1]]) %in% c("SpatRaster")) {
          dates <- as.Date(1:nlyr(tmin), origin=as.Date(paste0(years[1], "-01-01"))-1)
        }
      }
      
      dates.year <- lubridate::year(dates)
      
      tmin.out <- vector("list", length(unique(dates.year)))
      tmax.out <- vector("list", length(unique(dates.year)))
      period.dates <- vector("list", length(unique(dates.year)))
      
      # if(max(JDay)>365) {
      # yearRange <- unique(dates.year)[1:(length(unique(dates.year))-1)]
      yearRange <- years
      
      # } else {
      # yearRange <- unique(dates.year)
      # }
      for(i in seq_along(yearRange)) {
        firstDay <- which(dates.year %in% years[i])[JDay][1]
        lastDay <- firstDay + length(JDay)-1
        period <- firstDay:lastDay
        period.dates[[i]] <- dates[period]
        tmin.out[[i]] <- tmin[[period]]
        tmax.out[[i]] <- tmax[[period]]
      }
      
      # run MHT function across years range
      
      if(max(JDay)>365) {
        JDay[JDay>365] <- JDay[JDay>365]-365
      }
      
      daytimes <-  DL(lat, JDay)
      browser()
      # MHThourly_temps <- future_lapply(seq_along(years), function(x) MHT1(lat,
      #                                                                     tmin=tmin.out[[x]],
      #                                                                     tmax=tmax.out[[x]],
      #                                                                     dates = data.frame(Year = year(period.dates[[x]]),
      #                                                                                        JDay = yday(period.dates[[x]]))))
      
      MHThourly_temps <- future_lapply(seq_along(years), function(x) MHT(lat,
                                                                         Day_times=daytimes,
                                                                          tmin=tmin.out[[x]],
                                                                          tmax=tmax.out[[x]],
                                                                          dates = data.frame(Year = year(period.dates[[x]]),
                                                                                             JDay = yday(period.dates[[x]]))))
      
    #   for(i in seq_along(yearRange)) {
    #     
    #     MHThourly_temps <- MHT1(lat,
    #                             tmin=tmin.out[[i]],
    #                             tmax=tmax.out[[i]],
    #                             dates = data.frame(Year = year(period.dates[[i]]),
    #                                                JDay = yday(period.dates[[i]])))
    #   
    # }
  
  }

  # calculate chill portions from hourly temperatures object
    
    plan(multiprocess,  gc=TRUE)

    out <- lapply(MHThourly_temps, chill_spatial.years, JDay=JDay, template=template[[1]][[1]])
    
    future:::ClusterRegistry("stop")

  return(out)

}

getChillWorld <- function(scenario, model, year_range) {
  # convert year_range to vector
  years <- year_range[1]:year_range[length(year_range)-1] # added in the -1 to allow for northern hemisphere dormancy running over 2 calendar years
  # message(years)
  
  # get data files from scenario and model arguments
  dat.dir <- paste0('../../chillPortions/data/', scenario, '/')
  dat.files <- list.files(dat.dir, pattern=model, full.names = TRUE)
  dat.files <- grep('xml', dat.files, invert = TRUE, value = TRUE)
  # dat.files <- grep(year_range, dat.files, invert = TRUE, value = TRUE)
  tmin.files <- grep("tasmin", dat.files, value=TRUE)
  tmax.files <- grep("tasmax", dat.files, value=TRUE)
  
  
  if(year_range==1991:2010) {
    tmin.in <- c(rast(tmin.files[1]), rast(tmin.files[2]))# %>% aggregate(fact=2)
    tmax.in <- c(rast(tmax.files[1]), rast(tmax.files[2]))# %>% aggregate(fact=2)
  } 
  if(year_range==2041:2060) {
    tmin.in <- c(rast(tmin.files[1]), rast(tmin.files[2]))# %>% aggregate(fact=2)
    tmax.in <- c(rast(tmax.files[1]), rast(tmax.files[2]))# %>% aggregate(fact=2)
  } 
  if(year_range==2081:2100) {
    tmin.in <- c(rast(tmin.files[3]), rast(tmin.files[4]))# %>% aggregate(fact=2)
    tmax.in <- c(rast(tmax.files[3]), rast(tmax.files[4]))# %>% aggregate(fact=2)
  }
  
  # set dormancy period as a vector of julian days (for northern hemisphere dormancy, days in second calendar year should be represented as 365+n)
  JDay.south <- 92:306
  JDay.north <- 92:306 + 180
  
  # chill portions spatial calculations
  
  # produce spatial layer of latitudes and use to calculate day lengths
  lat <- tmin.in[[JDay.south]]
  xy <- coordinates(raster(tmin.in[[1]][[1]]))
  values(lat) <- xy[,2]
  
  #### run calculations on northern hemisphere ####
  
  tmin.north <- crop(tmin.in, ext(c(-180,180,0,90)))
  tmax.north <- crop(tmax.in, ext(c(-180,180,0,90)))
  
  rm(tmin.in, tmax.in)
  gc()
  
  # set a template raster to define the attributes of the final output raster
  # this can be a single layer of temperature data, for example
  template.ras <- tmin.north[[1]]
  
  # process daily temperatures
  t <- proc.time()
  files <- list.files('../../chillPortions/worldchill', full.names = TRUE, pattern='hourly')
  file.remove(files)
  chill_portions.north <- getChillSpatial(years, lat, JDay.north, tmin=tmin.north, tmax=tmax.north, template=template.ras)
  t1 <- t-proc.time()
  t1/60
  
  dir.create(paste0('../../chillPortions/worldchill/chill_portions/', scenario, '/', model, '/', year_range[10]), recursive = TRUE)
  
  chill_portions.north.raster <- rast(stack(lapply(chill_portions.north, raster)))
  
  writeRaster(chill_portions.north.raster,
              paste0('../../chillPortions/worldchill/chill_portions/',
                     scenario, '/', model, '/', scenario, '_', model, '_', year_range[10], '_', 'chill_portions_north.tif'),
              overwrite=TRUE)
  
  rm(hourly_temps.north, tmin.north, tmax.north)
  gc()
  
  #### run calculations on southern hemisphere ####

  if(year_range==1991:2010) {
    tmin.in <- c(rast(tmin.files[1]), rast(tmin.files[2]))# %>% aggregate(fact=2)
    tmax.in <- c(rast(tmax.files[1]), rast(tmax.files[2]))# %>% aggregate(fact=2)
  } 
  
  if(year_range==2041:2060) {
    tmin.in <- c(rast(tmin.files[1]), rast(tmin.files[2]))# %>% aggregate(fact=2)
    tmax.in <- c(rast(tmax.files[1]), rast(tmax.files[2]))# %>% aggregate(fact=2)
  } 
  if(year_range==2081:2100) {
    tmin.in <- c(rast(tmin.files[3]), rast(tmin.files[4]))# %>% aggregate(fact=2)
    tmax.in <- c(rast(tmax.files[3]), rast(tmax.files[4]))# %>% aggregate(fact=2)
  }
  
  tmin.south <- crop(tmin.in, ext(c(-180,180,-90,0)))
  tmax.south <- crop(tmax.in, ext(c(-180,180,-90,0)))
  
  rm(tmin.in, tmax.in)
  gc()
  
  # set a template raster to define the attributes of the final output raster
  # this can be a single layer of temperature data, for example
  template.ras <- tmin.south[[1]]#[[1]]
  
  # process daily temperatures
  t <- proc.time()
  files <- list.files('../../chillPortions/worldchill', full.names = TRUE, pattern='hourly')
  file.remove(files)
  chill_portions.south <- getChillSpatial(years, lat, JDay.south, tmin=tmin.south, tmax=tmax.south, template=template.ras)
  t1 <- t-proc.time()
  t1/60
  
  chill_portions.south.raster <- rast(stack(lapply(chill_portions.south, raster)))
  writeRaster(chill_portions.south.raster, 
              paste0('../../chillPortions/worldchill/chill_portions/', 
                     scenario, '/', model, '/', scenario, '_', model, '_', year_range[10], '_', 'chill_portions_south.tif'),
              overwrite=TRUE)
  rm(hourly_temps.south, tmin.south, tmax.south)
  gc()
  
  # stitch northern and southern hemisphere rasters together
  
  world <- lapply(seq_along(1:length(chill_portions.north)), function(x) {
    north.x <- raster(chill_portions.north[[x]])
    south.x <- raster(chill_portions.south[[x]])
    world.x <- mosaic(north.x, south.x, fun=min)
    return(world.x)
  })
  
  writeRaster(rast(stack(world)), 
              paste0('../../chillPortions/worldchill/chill_portions/', 
                     scenario, '/', model, '/', scenario, '_', model, '_', year_range[10], '_', 'chill_portions_world.tif'),
              overwrite=TRUE)
  
  return(world)
  
}



