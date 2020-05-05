
#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Calculate GDD for each wheat growing season

#====================================================================================== -

.libPaths("T:/R3UserLibs")

dirfrom <- "O:/Projects/KP0011/2/Data/data_raw/"
dirto <- "O:/Projects/KP0011/2/Data/data_prep/"

setwd(dirfrom)

#====================================================================================== -

#2018 ----

data_climate <- read.csv("data_climate_2018.csv", sep = ";", na.strings = "?")
data_climate$time_string <- strptime(data_climate$time_string, "%d.%m.%Y")

# calcualte mean daily temperature, using hourly means 
day_temp <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), mean, na.rm = TRUE)
names(day_temp)[2] <- "mean"

# get min and max hourly temperatures and calculate meanMM temperature  
day_temp$min <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), min, na.rm = TRUE)$x
day_temp$max <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), max, na.rm = TRUE)$x
day_temp$meanMM <- (day_temp$max + day_temp$min)/2

day_temp$day <- as.Date(day_temp$day)

# set mean values < 0 to 0 (base temperature)
day_temp$mean_calc <- ifelse(day_temp$mean > 0, day_temp$mean, 0)
day_temp$meanMM_calc <- ifelse(day_temp$meanMM > 0, day_temp$meanMM, 0)

# calcualte GDD and GDDMM
day_temp$GDD <- cumsum(day_temp$mean_calc)
day_temp$GDDMM <- cumsum(day_temp$meanMM_calc)

day_temp18 <- day_temp

# save GDD data
write.csv(day_temp[,-c(5,7)], paste0(dirto, "gdd_2018.csv"), row.names = FALSE)

#====================================================================================== -
#2019 ----

# data_climate <- read.csv("data_climate_2019.csv", sep = ";", na.strings = "?")
data_climate <- read.csv("data_climate_2019_db.csv")
data_climate$timestamp <- gsub("T", " ", data_climate$timestamp)
colnames(data_climate) <- c("time_string", "mean_temp")
# data_climate$time_string <- strptime(data_climate$time_string, "%d.%m.%Y")
data_climate$time_string <- strptime(data_climate$time_string, "%Y-%m-%d")

# calcualte mean daily temperature, using hourly means 
day_temp <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), mean, na.rm = TRUE)
names(day_temp)[2] <- "mean"

# get min and max hourly temperatures and calculate meanMM temperature  
day_temp$min <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), min, na.rm = TRUE)$x
day_temp$max <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), max, na.rm = TRUE)$x
day_temp$meanMM <- (day_temp$max + day_temp$min)/2

day_temp$day <- as.Date(day_temp$day)

# set mean values < 0 to 0 (base temperature)
day_temp$mean_calc <- ifelse(day_temp$mean > 0, day_temp$mean, 0)
day_temp$meanMM_calc <- ifelse(day_temp$meanMM > 0, day_temp$meanMM, 0)

# calcualte GDD and GDDMM
day_temp$GDD <- cumsum(day_temp$mean_calc)
day_temp$GDDMM <- cumsum(day_temp$meanMM_calc)

day_temp19 <- day_temp

# save GDD data
write.csv(day_temp[,-c(5,7)], paste0(dirto, "gdd_2019.csv"), row.names = FALSE)

#====================================================================================== -

day_temp <- bind_rows(day_temp18, day_temp19) %>% 
  dplyr::select(-5, -7)

data.table::fwrite(day_temp, paste0(dirto, "gdd.csv"))

#====================================================================================== -