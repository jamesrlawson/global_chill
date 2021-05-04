## split continuous daily data stack into a list of annual data stacks
# years is the range of years that the continuous stack spans
# dat is a 
annualise <- function(years, dat, startDay=0, useTerra=TRUE, writeDir=NULL, writeName=NULL) {
  
  # get date range using first year in year range provided
  
  if(class(dat) %in% c("RasterBrick", "RasterStack")) {
    dates <- as.Date(1:nlayers(dat), origin=as.Date(paste0(years[1], "-01-01"))-1)
  } else {
    if(class(dat) %in% c("SpatRaster")) {
      dates <- as.Date(1:nlyr(dat), origin=as.Date(paste0(years[1], "-01-01"))-1)
    }
  }
  
  dates.year <- lubridate::year(dates)
  
  if(!all(unique(dates.year)==years)) {
    stop("Year range provided does not match input data.")
  }
  
  dat.out <- vector("list", length(unique(dates.year)))
  for(i in seq_along(unique(dates.year))) {
    
    # extract a year's worth of daily data and store as list object
    message(years[i])
    subset <- which(dates.year %in% years[i]) + startDay
    
    # if startDay is not == 0, trim date range for the final year
    if(years[i] == max(years) & startDay != 0) {
      subset <- subset[subset <= max(subset)-startDay]
    }
    
    out <- dat[[subset]]
    
    if(useTerra & class(out) != 'SpatRaster') out <- terra::rast(out)
    
    dat.out[[i]] <- out
    
    # write to file (speeds up processing in later steps)
    if(!is.null(writeDir)){
      if(!dir.exists(writeDir)) dir.create(writeDir, recursive=TRUE, showWarnings=FALSE)
      if(useTerra) terra::writeRaster(dat.out[[i]], paste0(writeDir, "/", writeName, years[i], '.tif'), overwrite=TRUE)
      if(!useTerra) raster::writeRaster(dat.out[[i]], paste0(writeDir, "/", writeName, years[i], '.tif'), format="GTiff", overwrite=TRUE)
      dat.out[[i]] <- paste0(writeDir, "/", writeName, years[i], '.tif')
    }
    
  }
  
  return(dat.out)
  
}


#### Compute sunrise and sunset times and daylength spatially ####
# modified from chillR::daylength()

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
  
  CosWo <- (CosWo.1 - sin(latitude/360 * 2 * pi) * sin(Delta/360 * 2 * pi))/(cos(latitude/360 * 
                                                                                   2 * pi) * cos(Delta/360 * 2 * pi))
  normal_days <- as.vector(CosWo[] >= -1 & CosWo[] <= 1)
  
  Sunrise <- rep(-99, length(CosWo))
  
  Sunrise[normal_days[]] <- 12 - acos(CosWo[][normal_days])/(15/360 *
                                                               2 * pi)
  
  Sunset <- rep(-99, length(CosWo))
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
  
  Sunrise <- split(Sunrise, f=cell, drop=TRUE)
  Sunset <- split(Sunset, f=cell, drop=TRUE)
  Daylength <- split(Daylength, f=cell, drop=TRUE)
  
  return(list(Sunrise = Sunrise, Sunset = Sunset, Daylength = Daylength, JDay=rep(JDay, each=length(Sunrise))))#/length(JDay))))
}

#### Make hourly temperature record from daily data ####
# modified from chillR::make_hourly_temps()

MHT <- function (latitude, tmin=tmin, tmax=tmax, dates, daytimes=Day_times, lowMem=writeToDisk, keep_sunrise_sunset = FALSE, ...) 
{
  
  # if (missing(latitude)) 
  #   stop("'latitude' not specified")
  # if (length(latitude) > 1) 
  #   stop("'latitude' has more than one element")
  # if (!is.numeric(latitude)) 
  #   stop("'latitude' is not numeric")
  # if (latitude > 90 | latitude < (-90)) 
  #   warning("'latitude' is usually between -90 and 90")

  year = unique(dates$Year)
  JDay.provided = unique(dates$JDay)
  
  JDay <- c(JDay.provided[1]-1, JDay.provided, JDay.provided[length(JDay.provided)]+1, JDay.provided[length(JDay.provided)]+2)
  
  # JDay <- (min(JDay.provided)-1):(max(JDay.provided)+2) # expand JDay to account for truncation of 1 day at beginning and 2 days at end of provided dormancy period
  
  message(paste0('Calculating hourly temperatures for year: ', year))
  
  # lat <- tmin[[dates$JDay]]
  # xy <- coordinates(raster(tmin[[1]]))
  # lat[] <- xy[, 2]
  
  # message('  Calculating daylengths for period...')
  
  # Day_times <-  DL(lat, dates$JDay)
  
  # allcells <- 1:ncell(lat)
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
    dates.cell <- data.frame(Year = year,
                             JDay = JDay,
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
  
  message("  Reticulating splines..")
  
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
  dates <- dates[dates$JDay %in% JDay.provided,]
  close(pb)
  
  if(!dir.exists(lowMem)) {
    dir.create(lowMem, recursive=TRUE)
  }
  
  if(!is.null(lowMem)) {
    message('  Writing data to disk..')
    write_rds(dates, paste0(lowMem, '/hourly_temps', year, '.rds'))
  } else {
    return(dates)
  }
  
}

#### calculate hourly temps across a series of years ####
# wrapper for MHT()
# years is a vector of years corresponding to the list objects of tmin and tmax data
# daylight is calculated by DH()
# JDay is a vector of julian days describing the dormancy period of the plant
# writeToDisk is a string describing an output directory where hourly_temps files can be written to disk
  # this will prevent these outputs using too much system memory

MHT.years <- function(years, lat, JDay, tmin, tmax, writeToDisk=NULL) {
  
  years.seq <- seq_along(years)
  
  MHChourly_temps <- lapply(years.seq, function(x) MHT(lat,
                                                       tmin=tmin[[x]],
                                                       tmax=tmax[[x]],
                                                       dates = data.frame(Year = years[x], JDay = JDay)))
  if(is.null(writeToDisk)) {
    return(MHChourly_temps)
  }
  
}


#### calculate chill portions  ####

# takes output from MHT.years as input and calculates chill portions values for a specified raster cell
# sits inside wrapper chill_spatial.years which runs chill_spatial for each year of data

chill_spatial <- function(cell, hourly_temps, JDay) {

    out <- chillR::chilling(chillR::stack_hourly_temps(cell), Start_JDay=JDay[1], End_JDay=JDay[length(JDay)])$Chill_portions

  return(out)
}

# runs chill_spatial() across an entire raster, across a sequence of years
# use tmin or tmax as a source raster
chill_spatial.years <- function(years, hourly_temps=NULL, template, readFromDisk=NULL) {
  
  years.seq <- seq_along(years)
  
  if(!is.null(readFromDisk)) {
    message('Reading hourly temperature data from disk..')
    files <- list.files(writeToDisk, full.names = TRUE)
    hourly_temps <- lapply(files, read_rds)
  }
  
  out.list <- vector('list', length(years.seq))
  for(i in years.seq) {
    message(paste0("Calculating chill portions for year: ", years[i]))
    
    #split hourly_temps object into a list containing a dataframe describing hourly temperature values for each cell
    hourly_temps.list <- split(hourly_temps[[i]], f=hourly_temps[[i]]$cell,drop=TRUE)

    # run chill_spatial calculations over every dataframe in the above list
    out <- future_lapply(hourly_temps.list, chill_spatial, JDay=JDay)

    # take output values and assign to a template raster
    out.unlist <- unlist(out)
    chill.ras <- template
    cells <- which(!is.na(chill.ras[]))
    chill.ras[cells] <- out.unlist
    # output chill portions raster to a list
    out.list[[i]] <- chill.ras
  }
  return(out.list)
}

