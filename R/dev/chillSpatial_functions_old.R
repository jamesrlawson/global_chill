# ## split continuous daily data stack into a list of annual data stacks
# # years is the range of years that the continuous stack spans
# # dat is a 
# annualise <- function(years, dat, startDay=0, useTerra=TRUE, writeDir=NULL, writeName=NULL) {
#   
#   # get date range using first year in year range provided
#   
#   if(class(dat) %in% c("RasterBrick", "RasterStack")) {
#     dates <- as.Date(1:nlayers(dat), origin=as.Date(paste0(years[1], "-01-01"))-1)
#   } else {
#     if(class(dat) %in% c("SpatRaster")) {
#       dates <- as.Date(1:nlyr(dat), origin=as.Date(paste0(years[1], "-01-01"))-1)
#     }
#   }
#   
#   dates.year <- lubridate::year(dates)
#   
#   if(!all(unique(dates.year)==years)) {
#     stop("Year range provided does not match input data.")
#   }
#   
#   dat.out <- vector("list", length(unique(dates.year)))
#   for(i in seq_along(unique(dates.year))) {
#     
#     # extract a year's worth of daily data and store as list object
#     message(years[i])
#     subset <- which(dates.year %in% years[i]) + startDay
#     
#     # if startDay is not == 0, trim date range for the final year
#     if(years[i] == max(years) & startDay != 0) {
#       subset <- subset[subset <= max(subset)-startDay]
#     }
#     
#     out <- dat[[subset]]
#     
#     if(useTerra & class(out) != 'SpatRaster') out <- terra::rast(out)
#     
#     dat.out[[i]] <- out
#     
#     # write to file (speeds up processing in later steps)
#     if(!is.null(writeDir)){
#       if(!dir.exists(writeDir)) dir.create(writeDir, recursive=TRUE, showWarnings=FALSE)
#       if(useTerra) terra::writeRaster(dat.out[[i]], paste0(writeDir, "/", writeName, years[i], '.tif'), overwrite=TRUE)
#       if(!useTerra) raster::writeRaster(dat.out[[i]], paste0(writeDir, "/", writeName, years[i], '.tif'), format="GTiff", overwrite=TRUE)
#       dat.out[[i]] <- paste0(writeDir, "/", writeName, years[i], '.tif')
#     }
#     
#   }
#   
#   return(dat.out)
#   
# }
# 
# 
# #### Compute sunrise and sunset times and daylength spatially ####
# # modified from chillR::daylength()
# 
# DL <- function (latitude, JDay, notimes.as.na = FALSE) 
# {  
#   
#   # if (missing(latitude)) 
#   #   stop("'latitude' not specified")
#   # if (missing(JDay)) 
#   #   stop("'JDay' not specified")
#   # if (!isTRUE(all(is.numeric(JDay)))) 
#   #   stop("'JDay' contains non-numeric values")
#   # if (length(latitude) > 1) 
#   #   stop("'latitude' has more than one element")
#   # if (!is.numeric(latitude)) 
#   #   stop("'latitude' is not numeric")
#   # if (latitude > 90 | latitude < (-90)) 
#   #   warning("'latitude' is usually between -90 and 90")
#   Gamma <- 2 * pi/365 * ((JDay) - 1)
#   Delta <- 180/pi * (0.006918 - 0.399912 * cos(Gamma) + 0.070257 * 
#                        sin(Gamma) - 0.006758 * cos(Gamma) + 0.000907 * sin(Gamma) - 
#                        0.002697 * cos(3 * (Gamma)) + 0.00148 * sin(3 * (Gamma)))
#   
#   CosWo.1 <- latitude
#   values(CosWo.1) <- sin(-0.8333/360 * 2 * pi)
#   
#   CosWo <- (CosWo.1 - sin(latitude/360 * 2 * pi) * sin(Delta/360 * 2 * pi))/(cos(latitude/360 * 
#                                                                                    2 * pi) * cos(Delta/360 * 2 * pi))
#   normal_days <- as.vector(CosWo[] >= -1 & CosWo[] <= 1)
#   
#   Sunrise <- rep(-99, length(CosWo))
#   
#   Sunrise[normal_days[]] <- 12 - acos(CosWo[][normal_days])/(15/360 *
#                                                                2 * pi)
#   
#   Sunset <- rep(-99, length(CosWo))
#   Sunset[normal_days[]] <- 12 + acos(CosWo[][normal_days])/(15/360 * 
#                                                               2 * pi)
#   Daylength <- Sunset - Sunrise
#   Daylength[which(CosWo[] > 1)] <- 0
#   Daylength[which(CosWo[] < (-1))] <- 24
#   Sunrise[which(Daylength == 24)] <- 99
#   Sunset[which(Daylength == 24)] <- 99
#   if (notimes.as.na) {
#     Sunrise[which(Sunrise %in% c(-99, 99))] <- NA
#     Sunset[which(Sunset %in% c(-99, 99))] <- NA
#   }
#   Sunset[which(is.na(JDay))] <- NA
#   Sunrise[which(is.na(JDay))] <- NA
#   Daylength[which(is.na(JDay))] <- NA
#   
#   Sunrise[which(Sunrise == 99)] <- 0
#   Sunrise[which(Sunrise == -99)] <- 12
#   Sunset[which(Sunset == 99)] <- 24
#   Sunset[which(Sunset == -99)] <- 12
#   
#   times=rep(JDay, each=length(Sunrise)/length(JDay))
#   cell <- rep(1:ncell(lat), times=length(unique(JDay)))
#   
#   Sunrise <- split(Sunrise, f=cell, drop=TRUE)
#   Sunset <- split(Sunset, f=cell, drop=TRUE)
#   Daylength <- split(Daylength, f=cell, drop=TRUE)
#   
#   return(list(Sunrise = Sunrise, Sunset = Sunset, Daylength = Daylength, JDay=rep(JDay, each=length(Sunrise))))#/length(JDay))))
# }

#### Make hourly temperature record from daily data ####
# modified from chillR::make_hourly_temps()

MHT <- function (latitude, tmin=tmin, tmax=tmax, dates, keep_sunrise_sunset = FALSE, ...) 
{
  # browser()
  
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
    # message(n)
    # if(n==2236) browser()
    setTxtProgressBar(pb, n)
    
    dates.cell <- data.frame(Year = dates$Year,
                             JDay = JDay,
                             Tmax = tmax[n,],
                             Tmin = tmin[n,])
    
    hourly_temps.list[[n]] <- chillR::make_hourly_temps(lats[n], dates.cell)

  }
  
  message("  Reticulating splines..")
  # remove null elements of list
  hourly_temps.list[sapply(hourly_temps.list, is.null)] <- NULL
  
  close(pb)
  
  # if(writeToDisk) {
    dir.create('../../chillPortions/worldchill', recursive = TRUE, showWarnings = FALSE)
    message('  Writing data to disk..')
    write_rds(hourly_temps.list, paste0('../../chillPortions/worldchill', '/hourly_temps', this.year, '.rds'))
  # } else {
    # return(hourly_temps.list)
  # }
    return(NULL)
  
}

MHT1 <- function (latitude, tmin=tmin, tmax=tmax, dates, keep_sunrise_sunset = FALSE, ...) 
{
  # browser()
  
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
  # options(future.globals.maxSize= 991289600)
  # plan(multiprocess, multiprocess, gc=TRUE)
  # hourly_temps.list <- future_lapply(datcells, function(n, dates.year=dates$Year, JDay1=JDay, tmax1=tmax, tmin1=tmin) {
  #   # message(n)
  #     dates.cell <- data.frame(Year = dates.year,
  #                              JDay = JDay1,
  #                              Tmax = tmax1[n,],
  #                              Tmin = tmin1[n,])
  #     chillR::make_hourly_temps(lats[n], dates.cell)
  # })
  # 
  # future:::ClusterRegistry("stop")
  
  hourly_temps.list <- vector('list', length(datcells))
  pb <- txtProgressBar(max=length(datcells), style=3)
  for(n in datcells) {
    # message(n)
    # if(n==2236) browser()
    setTxtProgressBar(pb, n)

    dates.cell <- data.frame(Year = dates$Year,
                             JDay = JDay,
                             Tmax = tmax[n,],
                             Tmin = tmin[n,])

    hourly_temps.list[[n]] <- chillR::make_hourly_temps(lats[n], dates.cell)

  }

  message("  Reticulating splines..")
  # remove null elements of list
  hourly_temps.list[sapply(hourly_temps.list, is.null)] <- NULL
  return(hourly_temps.list)
  # close(pb)
  # 
  # # if(writeToDisk) {
  # dir.create('../../chillPortions/worldchill', recursive = TRUE, showWarnings = FALSE)
  # message('  Writing data to disk..')
  # write_rds(hourly_temps.list, paste0('../../chillPortions/worldchill', '/hourly_temps', this.year, '.rds'))
  # # } else {
  # # return(hourly_temps.list)
  # # }
  # return(NULL)
  
}


#### calculate hourly temps across a series of years ####
# wrapper for MHT()
# years is a vector of years corresponding to the list objects of tmin and tmax data
# daylight is calculated by DH()
# JDay is a vector of julian days describing the dormancy period of the plant
# writeToDisk is a string describing an output directory where hourly_temps files can be written to disk
  # this will prevent these outputs using too much system memory

MHT.years <- function(years, lat, JDay, tmin, tmax, writeToDisk=FALSE,...) {
  
  message('Annualizing temperature data..')
  
  if(class(tmin)=='list') {
    
    MHChourly_temps <- future_lapply(seq_along(years), function(x) MHT(lat,
                                                                    tmin=tmin[[x]],
                                                                    tmax=tmax[[x]],
                                                                    dates = data.frame(Year = years[x],
                                                                                       JDay = JDay)))
    
  } else {
    
    # should check that tmin and tmax represent the same date range (or are the same dim's)
    
    # subset input data to dormancy JDay range
    # change this to convert to SpatRaster if input data is in raster:: format
    if(class(tmin[[1]]) %in% c("RasterBrick", "RasterStack")) {
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
    yearRange <- unique(dates.year)[1:(length(unique(dates.year))-1)]
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
    
    # MHChourly_temps <- lapply(seq_along(yearRange), function(x) MHT(lat,
    #                                                      tmin=tmin.out[[x]],
    #                                                      tmax=tmax.out[[x]],
    #                                                      dates = data.frame(Year = year(period.dates[[x]]), 
    #                                                                         JDay = yday(period.dates[[x]]))))
    
    
    
    
    if(!writeToDisk) {
      MHChourly_temps <- lapply(seq_along(yearRange), function(x) MHT(lat,
                                                                      tmin=tmin.out[[x]],
                                                                      tmax=tmax.out[[x]],
                                                                      dates = data.frame(Year = year(period.dates[[x]]),
                                                                                         JDay = yday(period.dates[[x]]))))
      return(MHChourly_temps)
      
    } else {
      
      lapply(seq_along(yearRange), function(x) MHT(lat,
                                                   tmin=tmin.out[[x]],
                                                   tmax=tmax.out[[x]],
                                                   dates = data.frame(Year = year(period.dates[[x]]),
                                                                      JDay = yday(period.dates[[x]]))))
      return(NULL)
    }
    
  }
  
  
  
}


#### calculate chill portions  ####

# takes output from MHT.years as input and calculates chill portions values for a specified raster cell
# sits inside wrapper chill_spatial.years which runs chill_spatial for each year of data

chill_spatial <- function(cell, hourly_temps, JDay) {
  out <- chillR::chilling(chillR::stack_hourly_temps(cell), Start_JDay=JDay[1], End_JDay=JDay[length(JDay)])$Chill_portions
  return(out)
}

# # runs chill_spatial() across an entire raster, across a sequence of years
# # use tmin or tmax as a source raster
# chill_spatial.years <- function(years, hourly_temps=NULL, template, readFromDisk=FALSE, JDay=NULL) {
#   
#   browser()
#   
#   if(max(JDay)>365) {
#     JDay[JDay>365] <- JDay[JDay>365]-365
#   }
#   # years.seq <- 1:length(hourly_temps)
#   if(readFromDisk) {
#     message('Reading hourly temperature data from disk..')
#     files <- list.files('../../chillPortions/worldchill', full.names = TRUE, pattern='hourly')
#     # files <- grep('chill_portions', files, invert = TRUE, value = TRUE)
#     # hourly_temps <- lapply(files, read_rds)
#     years.seq <- 1:length(files)
#   } else {
#     years.seq <- seq_along(years)
#   }
#   
#   
#   out.list <- vector('list', length(years.seq))
#   for(i in years.seq) {
#     message(paste0("Calculating chill portions for year: ", years[i]))
#     hourly_temps.list <- read_rds(files[i])
#     # run chill_spatial calculations over every dataframe in the above list
#     out <- future_lapply(hourly_temps.list, chill_spatial, JDay=JDay)
#     # take output values and assign to a template raster
#     out.unlist <- unlist(out)
#     chill.ras <- template
#     cells <- which(!is.na(chill.ras[]))
#     chill.ras[cells] <- out.unlist
#     # output chill portions raster to a list
#     out.list[[i]] <- chill.ras
#   }
#   return(out.list)
# }

chill_spatial.years <- function(hourly_temps, template, readFromDisk=FALSE, JDay=NULL) {
    # out.list <- vector('list', length(hourly_temps))
  # for(i in seq_along(hourly_temps)) {
    # message(paste0("Calculating chill portions for year: ", years[i]))
    # run chill_spatial calculations over every dataframe in the above list
    out <- future_lapply(hourly_temps, chill_spatial, JDay=JDay)
    # take output values and assign to a template raster
    out.unlist <- unlist(out)
    chill.ras <- template
    cells <- which(!is.na(chill.ras[]))
    chill.ras[cells] <- out.unlist
    # output chill portions raster to a list
    # out.list[[i]] <- chill.ras
  # }
  # return(out.list)
    return(chill.ras)
}



getChillSpatial <- function(years, lat, JDay, tmin, tmax, template, writeToDisk=FALSE,...) {
  message('Annualizing temperature data..')
  
  # should check that tmin and tmax represent the same date range (or are the same dim's)
  
  # subset input data to dormancy JDay range
  # change this to convert to SpatRaster if input data is in raster:: format
  
  if(class(tmin)=='list') {
    
    MHChourly_temps <- future_lapply(seq_along(years), function(x) MHT1(lat,
                                                                tmin=tmin[[x]],
                                                                tmax=tmax[[x]],
                                                                dates = data.frame(Year = years[x],
                                                                                   JDay = JDay)))
    
  } else {

      if(class(tmin[[1]]) %in% c("RasterBrick", "RasterStack")) {
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
      
      chillPortions.list <- vector('list', length(yearRange))
      for(i in seq_along(yearRange)) {
        
        MHChourly_temps <- MHT1(lat,
                                tmin=tmin.out[[i]],
                                tmax=tmax.out[[i]],
                                dates = data.frame(Year = year(period.dates[[i]]),
                                                   JDay = yday(period.dates[[i]])))
      
    }
  
  }
    
    plan(multiprocess,  gc=TRUE)
    # plan(multiprocess)
    # out <- future_lapply(MHChourly_temps, chill_spatial, JDay=JDay)
    
    
    # out <- vector("list", length=length(MHChourly_temps[[1]]))
    # 
    # for(i in seq_along(MHChourly_temps[[1]])) {
    #   message(i)
    #   cell <- MHChourly_temps[[1]][[i]]
    #   out[[i]] <- chillR::chilling(chillR::stack_hourly_temps(cell), Start_JDay=JDay[1], End_JDay=JDay[length(JDay)])$Chill_portions
    # }

    # out <- lapply(MHChourly_temps, chill_spatial, JDay=JDay)
    out <- lapply(MHChourly_temps, chill_spatial.years, JDay=JDay, template=template[[1]][[1]])
    
    future:::ClusterRegistry("stop")
    # # take output values and assign to a template raster
    # out.unlist <- unlist(out)
    # chill.ras <- template
    # cells <- which(!is.na(chill.ras[]))
    # chill.ras[cells] <- out.unlist
    # 
    # chillPortions.list[[i]] <- chill.ras
    # 
 
  
  # chillPortions <- lapply(seq_along(yearRange), function(x) {
  #   browser()
  #   MHChourly_temps <- MHT(lat,
  #           tmin=tmin.out[[x]],
  #           tmax=tmax.out[[x]],
  #           dates = data.frame(Year = year(period.dates[[x]]),
  #                              JDay = yday(period.dates[[x]])))
  #   
  #   # run chill_spatial calculations over every dataframe in the above list
  #   out <- future_lapply(MHChourly_temps, chill_spatial, JDay=JDay)
  #   # take output values and assign to a template raster
  #   out.unlist <- unlist(out)
  #   chill.ras <- template
  #   cells <- which(!is.na(chill.ras[]))
  #   chill.ras[cells] <- out.unlist
  #   
  #   chill.ras
  #   
  #   })
  
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
  # JDay.north[JDay.north > 365] <- JDay.north[JDay.north > 365] - 365
  # JDay.north <- JDay.north[1:365]
  # JDay.north <- JDay.north[!is.na(JDay.north)]
  
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
  template.ras <- tmin.north[[1]]#[[1]]
  
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
  
  # tmin.in <- c(rast(tmin.files[1]), rast(tmin.files[2])) %>% aggregate(fact=2)
  # tmax.in <- c(rast(tmax.files[1]), rast(tmax.files[2])) %>% aggregate(fact=2)
  
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
  
  # write_rds(chill_portions.south, '../../chillPortions/worldchill/chill_portions/chill_portions.south.rds')
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



