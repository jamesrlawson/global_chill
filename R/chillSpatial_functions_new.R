# calculate day lengths from latitudes

DL <- function (latitude, JDay, notimes.as.na = FALSE) 
{  
  
  # if (missing(latitude)) 
  #   stop("'latitude' not specified")
  # if (missing(JDay)) 
  #   stop("'JDay' not specified")
  # if (!isTRUE(all(is.numeric(JDay)))) 
  #   stop("'JDay' contains non-numeric values")
  # if (length(latitude) > 1) 
  #   stop("'latitude' has more than one element")
  # if (!is.numeric(latitude)) 
  #   stop("'latitude' is not numeric")
  # if (latitude > 90 | latitude < (-90)) 
  #   warning("'latitude' is usually between -90 and 90")
  Gamma <- 2 * pi/365 * ((JDay) - 1)
  Delta <- 180/pi * (0.006918 - 0.399912 * cos(Gamma) + 0.070257 * 
                       sin(Gamma) - 0.006758 * cos(Gamma) + 0.000907 * sin(Gamma) - 
                       0.002697 * cos(3 * (Gamma)) + 0.00148 * sin(3 * (Gamma)))
  
  CosWo.1 <- latitude
  values(CosWo.1) <- sin(-0.8333/360 * 2 * pi)
  
  Delta.sin <- sin(Delta/360 * 2 * pi)
  Delta.cos <- cos(Delta/360 * 2 * pi)
  
  CosWo.a <- (CosWo.1 - sin(latitude/360 * 2 * pi))
  CosWo.b <- (cos(latitude/360 * 2 * pi))
  
  termA <- rast(stack(lapply(seq_along(Delta.sin), function(x) Delta.sin[x] * CosWo.a[[x]])))
  termB <- rast(stack(lapply(seq_along(Delta.cos), function(x) Delta.cos[x] * CosWo.b[[x]])))
  
  CosWo <- termA/termB
  
  # CosWo <- (CosWo.1 - sin(latitude/360 * 2 * pi) * sin(Delta/360 * 2 * pi))/(cos(latitude/360 * 
  # 2 * pi) * cos(Delta/360 * 2 * pi))
  normal_days <- as.vector(CosWo[] >= -1 & CosWo[] <= 1)
  
  Sunrise <- rep(-99, length(CosWo[]))
  Sunrise[normal_days[]] <- 12 - acos(CosWo[][normal_days])/(15/360 *
                                                               2 * pi)
  Sunset <- rep(-99, length(CosWo[]))
  Sunset[normal_days[]] <- 12 + acos(CosWo[][normal_days])/(15/360 * 
                                                              2 * pi)
  Daylength <- Sunset - Sunrise
  Daylength[which(CosWo[] > 1)] <- 0
  Daylength[which(CosWo[] < (-1))] <- 24
  Sunrise[which(Daylength == 24)] <- 99
  Sunset[which(Daylength == 24)] <- 99
  if (notimes.as.na) {
    Sunrise[which(Sunrise %in% c(-99, 99))] <- NA
    Sunset[which(Sunset %in% c(-99, 99))] <- NA
  }
  Sunset[which(is.na(JDay))] <- NA
  Sunrise[which(is.na(JDay))] <- NA
  Daylength[which(is.na(JDay))] <- NA
  
  Sunrise[which(Sunrise == 99)] <- 0
  Sunrise[which(Sunrise == -99)] <- 12
  Sunset[which(Sunset == 99)] <- 24
  Sunset[which(Sunset == -99)] <- 12
  
  times=rep(JDay, each=length(Sunrise)/length(JDay))
  cell <- rep(1:ncell(lat), times=length(unique(JDay)))
  
  Sunrise <- round(Sunrise,2)
  Sunset <- round(Sunset,2)
  Daylength <- round(Daylength,2)
  
  Sunrise <- split(Sunrise, f=cell, drop=TRUE)
  Sunset <- split(Sunset, f=cell, drop=TRUE)
  Daylength <- split(Daylength, f=cell, drop=TRUE)
  
  return(list(Sunrise = Sunrise, Sunset = Sunset, Daylength = Daylength, JDay=rep(JDay, each=length(Sunrise))))#/length(JDay))))
}

# calculate chill portions from daily tmin and tmax
getCP <- function (latitude, tmin=tmin, tmax=tmax, dates, Day_times=daytimes, keep_sunrise_sunset = FALSE, template=NULL, ...) 
{
  # if (missing(latitude)) 
  #   stop("'latitude' not specified")
  # if (length(latitude) > 1) 
  #   stop("'latitude' has more than one element")
  # if (!is.numeric(latitude)) 
  #   stop("'latitude' is not numeric")
  # if (latitude > 90 | latitude < (-90)) 
  #   warning("'latitude' is usually between -90 and 90")
  
  datcells <- which(!is.na(tmin[[1]][[1]][]))
  
  message('  Loading temperature data into memory..')
  
  tmax <- as.matrix(tmax)
  tmin <- as.matrix(tmin)
  
  message('  Processing daylengths..')
  
  dates.list <- vector('list', ncell(lat))
  pb <- txtProgressBar(max=length(datcells), style=3)
  for(n in datcells) {
    # message(n)
    setTxtProgressBar(pb, n)
    dates.cell <- data.frame(Year = dates$Year,
                             JDay = dates$JDay,
                             Tmax = tmax[n,][JDay],
                             Tmin = tmin[n,][JDay],
                             cell = n,
                             Sunrise = NA,
                             Sunset = NA,
                             Daylength = NA,
                             prev_Sunset = NA,
                             next_Sunrise = NA,
                             prev_max = NA,
                             next_min = NA,
                             prev_min = NA)
    
    Day_times.cell <- list(Sunrise=Day_times$Sunrise[[n]], Sunset=Day_times$Sunset[[n]], Daylength=Day_times$Daylength[[n]])
    dates.cell$Sunrise[2:(length(Day_times.cell$Sunrise) - 1)] <- Day_times.cell$Sunrise[2:(length(Day_times.cell$Sunrise) - 1)]
    dates.cell$Sunset[2:(length(Day_times.cell$Sunset) - 1)] <- Day_times.cell$Sunset[2:(length(Day_times.cell$Sunset) - 1)]
    dates.cell$Daylength[2:(length(Day_times.cell$Daylength) - 1)] <- Day_times.cell$Daylength[2:(length(Day_times.cell$Daylength) - 1)]
    dates.cell$prev_Sunset[1:(length(Day_times.cell$Sunset) - 2)] <- Day_times.cell$Sunset[1:(length(Day_times.cell$Sunset) - 2)]
    dates.cell$next_Sunrise[3:length(Day_times.cell$Sunrise)] <- Day_times.cell$Sunrise[3:length(Day_times.cell$Sunrise)]
    dates.cell$prev_max[c(1:(nrow(dates.cell) - 1))] <- dates.cell$Tmax[c(1:(nrow(dates.cell) - 1))]
    dates.cell$next_min[c(2:nrow(dates.cell))] <- dates.cell$Tmin[c(2:nrow(dates.cell))]
    dates.cell$prev_min[c(1:(nrow(dates.cell) - 1))] <- dates.cell$Tmin[c(1:(nrow(dates.cell) - 1))]
    dates.cell$Tsunset <- dates.cell$Tmin + (dates.cell$Tmax - dates.cell$Tmin) *
      sin((pi * (dates.cell$Sunset - dates.cell$Sunrise)/(dates.cell$Daylength + 4)))
    dates.cell$prev_Tsunset <- dates.cell$prev_min + (dates.cell$prev_max -  dates.cell$prev_min) * sin((pi * (dates.cell$Daylength)/(dates.cell$Daylength + 4)))
    
    dates.list[[n]] <- dates.cell
  }
  
  close(pb)
  
  dates <- do.call(rbind, dates.list)
  
  message("  Interpolating hourly temperatures..")
  
  pb <- txtProgressBar(max=23, style=3)
  
  colnum <- ncol(dates) + 1
  hourcol <- c(colnum:(colnum + 23))
  for (hour in 0:23) {
    setTxtProgressBar(pb, hour)
    hourcount <- hour + 1
    no_riseset <- which(dates$Daylength %in% c(0, 24, 
                                               -99))
    dates[no_riseset, colnum + hour] <- ((dates$Tmax + 
                                            dates$Tmin)/2)[no_riseset]
    c_morn <- which(hour <= dates$Sunrise)
    if (1 %in% c_morn)
      if (!length(c_morn) == 1)
        c_morn <- c_morn[2:length(c_morn)]
    else c_morn <- c()
    c_day <- which(hour > dates$Sunrise & hour <= dates$Sunset)
    c_eve <- which(hour >= dates$Sunset)
    if (nrow(dates) %in% c_eve) 
      c_eve <- c_eve[1:(length(c_eve) - 1)]
    dates[c_morn, colnum + hour] <- dates$prev_Tsunset[c_morn] - 
      ((dates$prev_Tsunset[c_morn] - dates$Tmin[c_morn])/log(max(1, 24 - (dates$prev_Sunset[c_morn] - dates$Sunrise[c_morn]),na.rm=TRUE)) * 
         log(hour + 24 - dates$prev_Sunset[c_morn] + 
               1))
    dates[c_day, colnum + hour] <- dates$Tmin[c_day] + 
      (dates$Tmax[c_day] - dates$Tmin[c_day]) * 
      sin((pi * (hour - dates$Sunrise[c_day])/(dates$Daylength[c_day] + 
                                                 4)))
    dates[c_eve, colnum + hour] <- dates$Tsunset[c_eve] - 
      ((dates$Tsunset[c_eve] - dates$next_min[c_eve])/log(24 - 1) * log(hour - dates$Sunset[c_eve] + 1))
  }
  
  preserve_columns <- c('Year', 'JDay', 'Tmin', 'Tmax', 'cell')
  
  colnames(dates)[(ncol(dates) - 23):(ncol(dates))] <- c(paste("Hour_", 
                                                               0:23, sep = ""))
  if (!keep_sunrise_sunset) 
    dates <- dates[, c(preserve_columns, paste("Hour_", 
                                               0:23, sep = ""))]
  if (keep_sunrise_sunset) 
    dates <- dates[, c(preserve_columns, "Sunrise", 
                       "Sunset", "Daylength", paste("Hour_", 
                                                    0:23, sep = ""))]
  dates[1, (ncol(dates) - 23):(ncol(dates))][which(is.na(dates[1, 
                                                               (ncol(dates) - 23):(ncol(dates))]))] <- dates[1,"Tmin"]
  dates[nrow(dates), (ncol(dates) - 23):(ncol(dates))][which(is.na(dates[nrow(dates), 
                                                                         (ncol(dates) - 23):(ncol(dates))]))] <- dates[nrow(dates), 
                                                                                                                       "Tmin"]
  # dates <- dates[dates$JDay %in% JDay.provided,]
  close(pb)
  
  
  dates <- round(dates,2)
  dates <- split(dates, f=dates$cell)
  
  message('  Calculating chill portions..')
  sourceCpp('../../chillPortions/global_chill/Cpp/chill_func.cpp')
  CP <- lapply(dates, function(x) chill_func(na.omit(as.vector(t(x[,6:29]))))) # 6:29 refers to columns the hourly temperature data
  # na.omit(as.vector(t())) converts the dataframe of values to a vector by concatenating rows, giving a vector of hourly temperatures
  template[datcells] <- NA
  template[datcells] <- unlist(CP)
  return(template)
  # }
  
}

# from https://github.com/RPertille/ChillModels
dynamic_model <- function(x,total=TRUE){
  e0 <- 4153.5
  e1 <- 12888.8
  a0 <- 139500
  a1 <- 2.567e+18
  slp <- 1.6
  tetmlt <- 277
  aa <- a0/a1
  ee <- e1 - e0
  TK <- x + 273
  ftmprt <- slp * tetmlt * (TK - tetmlt)/TK
  sr <- exp(ftmprt)
  xi <- sr/(1 + sr)
  xs <- aa * exp(ee/TK)
  ak1 <- a1 * exp(-e1/TK)
  interE <- 0
  memo <- new.env(hash = TRUE)
  posi <- 1
  assign(x = paste(1), value = 0, envir = memo)
  E = 0
  S <- ak1
  S[1] <- 0
  E <- S
  options(scipen = 30)
  for (l in 2:length(x)) {
    if (E[l - 1] < 1) {
      S[l] <- E[l - 1]
      E[l] <- xs[l] - (xs[l] - S[l]) * exp(-ak1[l])
    }
    else {
      S[l] <- E[l - 1] - E[l - 1] * xi[l - 1]
      E[l] <- xs[l] - (xs[l] - S[l]) * exp(-ak1[l])
    }
  }
  interE <- E
  y <- rep(0, length(x))
  y[which(interE >= 1)] <- interE[which(interE >= 1)] * 
    xi[which(interE >= 1)]
  if (total == TRUE) 
    return(tail(cumsum(y),n=1))
  else return(y)
}


# parent function for calculating chill portions
getChillSpatial <- function(years, lat, JDay, tmin, tmax, template, writeToDisk=FALSE,...) {
  message('Annualizing temperature data..')
  # should check that tmin and tmax represent the same date range (or are the same dim's)
  
  # subset input data to dormancy JDay range
  # change this to convert to SpatRaster if input data is in raster:: format
  if(class(tmin)=='list') {

    # interpolate hourly temperatures  

    daytimes <-  DL(lat, JDay)
    
    CP <- lapply(seq_along(years), function(x) getCP(lat,
                                                   tmin=tmin[[x]],
                                                   tmax=tmax[[x]],
                                                   Day_times=daytimes,
                                                   template=tmin[[1]][[1]],
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
      
      # run getCP function across years range
      
      if(max(JDay)>365) {
        JDay[JDay>365] <- JDay[JDay>365]-365
      }
      
      daytimes <-  DL(lat, JDay)

      CP <- future_lapply(seq_along(years), function(x) getCP(lat,
                                                            Day_times=daytimes,
                                                            tmin=tmin.out[[x]],
                                                            tmax=tmax.out[[x]],
                                                            dates = data.frame(Year = year(period.dates[[x]]),
                                                                               JDay = yday(period.dates[[x]]))))
      
  }

  return(CP)

}

getChillWorld <- function(scenario, model, year_range) {
  # convert year_range to vector
  # years <- year_range[1]:year_range[length(year_range)-1] # added in the -1 to allow for northern hemisphere dormancy running over 2 calendar years
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



