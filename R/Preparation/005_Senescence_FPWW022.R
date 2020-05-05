
#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Conversion of senescence scoring data to GDDAH and GDDAS 

#====================================================================================== -

rm(list = ls())
.libPaths("T:R3UserLibs")
library(dplyr)

base_path <- "O:/Projects/KP0011/2/Data/"

#senescence data
data2018 <- read.csv(paste0(base_path,"data_raw/Senescence_FPWW022.csv"), sep = ";", check.names = FALSE)
#heading data
d.head2018 <- read.csv(paste0(base_path, "data_prep/Heading_FPWW022_DAS.csv"))
#gdd data
day_temp <- read.csv(paste0(base_path, "data_prep/gdd_2018.csv"))

#====================================================================================== -

data <- right_join(d.head2018, data2018)
   
data_Fl0 <- dplyr::select(data, c(1:8, starts_with("Fl0"))) %>% 
  gather(object_rating, SnsFl0, `Fl0 14.06.2018`:`Fl0 11.07.2018`, factor_key = TRUE)
data_Cnp <- dplyr::select(data, c(1:8, starts_with("Can"))) %>% 
  gather(object_rating, SnsCnp, `Can 14.06.2018`:`Can 11.07.2018`, factor_key = TRUE)
data <- bind_cols(data_Fl0, data_Cnp)

data$grading_date <- sapply(lapply(lapply(as.character(data$object_rating), strsplit, " ", fixed = TRUE), unlist), function(x) x=as.character(x[2])) %>% 
  as.Date(., "%d.%m.%Y")
data$grading_DAS <- data$grading_date - as.Date("17.10.2017", "%d.%m.%Y")
data$heading_DAS <- as.difftime(data$HeadingDAS, units = "days")
data$grading_DAH <- data$grading_DAS - data$heading_DAS

#calculate GDDAH
grading_GDDAH <- NULL
for (i in 1:nrow(data)) {
  if (!is.na(data$heading_date)[i]) {
    grading_GDDAH[i] <- day_temp[day_temp$day == paste(data$grading_date)[i], 6] -
      day_temp[day_temp$day == paste(data$heading_date)[i], 6]
  }
  else paste("NA")
  print(grading_GDDAH[i])
}
data$grading_GDDAH <- as.numeric(round(grading_GDDAH, 0))

#calculate GDDAS
grading_GDDAS <- NULL
for(i in 1:nrow(data)){
  if(!is.na(data$heading_date)[i]){
    grading_GDDAS[i] <- as.numeric(data$heading_GDDAS[i]) + as.numeric(data$grading_GDDAH[i])
  }
  else {
    grading_GDDAS[i] <- day_temp[day_temp$day == paste(data$grading_date)[i], 6]
  }
  print(grading_GDDAS[i])
}
data$grading_GDDAS <- as.numeric(round(grading_GDDAS, 0))

# tidy up data
data <- data %>% dplyr::select(Plot_ID, Lot, RangeLot, RowLot, Gen_Name, grading_date, 
                               heading_date, heading_DAS, heading_GDDAS, grading_DAS, 
                               grading_DAH, grading_GDDAH, grading_GDDAS, SnsFl0, SnsCnp)

write.csv(as.data.frame(data), paste0(base_path, "data_prep/Senescence_FPWW022_DAS.csv"), row.names = FALSE)

#====================================================================================== -
