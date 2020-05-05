
#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Estimation of heading date, based on interpolation of individual assessments
# convertion of heading date in days after sowing (heading_DAS) and growing degree days after sowing (heading_GDDAS)

#====================================================================================== -

rm(list = ls())
.libPaths("T:/R3UserLibs")
library(tidyverse)

dirfrom <- "O:/Projects/KP0011/2/Data/data_raw/"
dirto <- "O:/Projects/KP0011/2/Data/data_prep/"

setwd(dirfrom)

filenames <- as.list(list.files(pattern = "Heading_FPWW024_FPWW028.csv"))

data <- lapply(filenames, read_csv) %>% bind_rows()

# heading data
d.head <-data %>% 
  dplyr::select(plot.UID, year_site.UID, plot.range_in_lot, plot.row_in_lot, genotype.name, trait_value.timestamp, trait_value.value) %>% 
  mutate(Lot = strsplit(as.character(year_site.UID), "_") %>% lapply("[[", 2) %>% unlist() %>% gsub("lot", "", .)) %>% 
  mutate(trait_value.timestamp = as.Date(str_sub(trait_value.timestamp, 1,10), "%Y-%m-%d")) %>% 
  mutate(HeadingDAS = as.factor(trait_value.timestamp - as.Date("2018-10-17", "%Y-%m-%d"))) %>% 
  dplyr::select(plot.UID, Lot, plot.range_in_lot, plot.row_in_lot, genotype.name, trait_value.timestamp, HeadingDAS, trait_value.value)

names <- c("Plot_ID", "Lot", "RangeLot", "RowLot", "Gen_Name", "heading_date", "HeadingDAS", "heading_GDDAS")
names(d.head) <- names

#====================================================================================== -

write.csv(d.head, paste0(dirto, "Heading_FPWW024_FPWW028_DAS.csv"), row.names = FALSE)

#====================================================================================== -