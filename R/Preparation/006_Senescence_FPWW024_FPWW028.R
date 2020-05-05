
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
data_FPWW024 <- read.csv(paste0(base_path, "data_raw/Senescence_FPWW024.csv"), sep = ";", check.names = FALSE) %>% dplyr::select(-plot_number, -row, -range)
data_FPWW028 <- read.csv(paste0(base_path, "data_raw/Senescence_FPWW028.csv"), sep = ";", check.names = FALSE) %>% dplyr::select(-row, -range)
data2019 <- rbind(data_FPWW024, data_FPWW028) %>% 
  rename("Plot_ID" = plot_label,
         "RangeLot" = range_in_lot,
         "RowLot" = row_in_lot) %>% 
  dplyr::select(Plot_ID, Lot, RangeLot, RowLot, Gen_Name, everything())
#heading data
d.head2019 <- read.csv(paste0(base_path, "data_prep/Heading_FPWW024_FPWW028_DAS.csv"))
#gdd data
day_temp <- read.csv(paste0(base_path, "data_prep/gdd_2019.csv"))

#====================================================================================== -

data <- right_join(d.head2019, data2019)

data_Fl0 <- dplyr::select(data, c(1:8, starts_with("Fl0"))) %>% 
  gather(object_rating, SnsFl0, `Fl0 25.06.2019`:`Fl0 17.07.2019`, factor_key = TRUE)
data_Cnp <- dplyr::select(data, c(1:8, starts_with("Cnp"))) %>% 
  gather(object_rating, SnsCnp, `Cnp 25.06.2019`:`Cnp 17.07.2019`, factor_key = TRUE)
data <- bind_cols(data_Fl0, data_Cnp)

data$grading_date <- sapply(lapply(lapply(as.character(data$object_rating), strsplit, " ", fixed = TRUE), unlist), function(x) x=as.character(x[2])) %>% 
  as.Date(., "%d.%m.%Y")
data$grading_DAS <- data$grading_date - as.Date("17.10.2018", "%d.%m.%Y")
data$heading_DAS <- as.difftime(data$HeadingDAS, units = "days")
data$grading_DAH <- data$grading_DAS - data$heading_DAS

#calculate GDDAH
grading_GDDAH <- NULL
for (i in 1:nrow(data)) {
  if (!is.na(data$heading_date)[i]) {
    grading_GDDAH[i] <- day_temp[day_temp$day == paste(data$grading_date)[i], 6] -
      day_temp[day_temp$day == paste(data$heading_date)[i], 6]
  }
  else {
    grading_GDDAH[i] <- paste("NA")
  }
  print(grading_GDDAH[i])
}
data$grading_GDDAH <- as.numeric(round(as.numeric(grading_GDDAH), 0)) 

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
data$grading_GDDAS <- as.numeric(round(as.numeric(grading_GDDAS), 0))

# tidy up data
data <- data %>% dplyr::select(Plot_ID, Lot, RangeLot, RowLot, Gen_Name, grading_date, 
                               heading_date, heading_DAS, heading_GDDAS, grading_DAS, 
                               grading_DAH, grading_GDDAH, grading_GDDAS, SnsFl0, SnsCnp)

write.csv(as.data.frame(data), paste0(base_path, "data_prep/Senescence_FPWW024_FPWW028_DAS.csv"), row.names = FALSE)

#====================================================================================== -
