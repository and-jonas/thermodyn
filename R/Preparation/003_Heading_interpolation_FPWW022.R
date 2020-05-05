
#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Estimation of heading date, based on interpolation of individual assessments
# convertion of heading date in days after sowing (heading_DAS) and growing degree days after sowing (heading_GDDAS)

#====================================================================================== -

rm(list = ls())
.libPaths("T:/R3UserLibs")
library(tidyverse)

base_path <- "O:/Projects/KP0011/2/Data/"

#heading data
d.head <- read.csv(paste0(base_path, "data_raw/Heading_FPWW022.csv"), na.strings = c("NA",""), sep = ";") %>% 
  mutate_at(vars(matches("X")), as.numeric)

#GDD data
day_temp <- read.csv(paste0(base_path, "data_prep/gdd_2018.csv"), sep = ",")

#====================================================================================== -

#heading_DAS ----

plots <- d.head[,1]

#estimate heading date in DAS
HeadingDAS <- NULL
for (i in 1:length(plots)){
  x<- d.head[i,-c(1:6)]
  if(is.na(x[3])){
    print("NA")
    HeadingDAS<-c(HeadingDAS, NA)
  } else if(x[1] == 60){
    print("213")
    HeadingDAS<-c(HeadingDAS, 213)
  } else if (x[2] == 60){
    print("216")
    HeadingDAS<-c(HeadingDAS, 216)
  } else if (x[3] == 60){
    print("218")
    HeadingDAS<-c(HeadingDAS, 218)    
  } else if (x[4] == 60){
    print("220")
    HeadingDAS<-c(HeadingDAS, 220)    
  } else if (x[5] == 60){
    print("223")
    HeadingDAS<-c(HeadingDAS, 223)    
  } else if (x[6] == 60){
    print("225")
    HeadingDAS<-c(HeadingDAS, 225)    
  } else if (x[1] > 56 & x[1] < 60 & x[2] >= 60){
    print("214")
    HeadingDAS<-c(HeadingDAS, 214)
  } else if (x[1] > 53 & x[1] < 60 & x[2] >= 60){
    print("215")
    HeadingDAS<-c(HeadingDAS, 215)
  } else if (x[2] > 56 & x[2] < 60 & x[3] >= 60){
    print("217")
    HeadingDAS<-c(HeadingDAS, 217)
  } else if (x[2] > 50 & x[2] < 60 & x[3] >= 60){
    print("218")
    HeadingDAS<-c(HeadingDAS, 218)
  } else if (x[3] > 56 & x[3] < 60 & x[4] >= 60){
    print("219")
    HeadingDAS<-c(HeadingDAS, 219)
  } else if (x[3] > 50 & x[3] < 60 & x[4] >= 60){
    print("220")
    HeadingDAS<-c(HeadingDAS, 220)
  } else if (x[4] > 56 & x[4] < 60 & x[5] >= 60){
    print("221")
    HeadingDAS<-c(HeadingDAS, 221)
  } else if (x[4] > 53 & x[4] < 60 & x[5] >= 60){
    print("222")
    HeadingDAS<-c(HeadingDAS, 222)
  } else if (x[4] > 50 & x[4] < 60 & x[5] >= 60){
    print("223")
    HeadingDAS<-c(HeadingDAS, 223)
  } else if (x[5] > 56 & x[5] < 60 & x[6] >= 60){
    print("224")
    HeadingDAS<-c(HeadingDAS, 224) 
  } else if (x[5] > 50 & x[5] < 60 & x[6] >= 60){
    print("225")
    HeadingDAS<-c(HeadingDAS, 225) 
  } else if (x[6] > 56 & x[6] < 60 & x[7] >= 60){
    print("226")
    HeadingDAS<-c(HeadingDAS, 226) 
  } else{
    print("NA")
    HeadingDAS<-c(HeadingDAS, NA)
  }
}

HeadingDAS <- cbind(d.head[,1:5], HeadingDAS)

#====================================================================================== -

#heading_date ----

HeadingDAS$heading_date <- as.Date("17.10.2017", "%d.%m.%Y") + HeadingDAS$HeadingDAS

#====================================================================================== -

#Quality Check
HeadingDAS[HeadingDAS$Gen_Name=="CH CLARO", ] #OK
HeadingDAS[HeadingDAS$Gen_Name=="SURETTA", ] #OK
HeadingDAS[HeadingDAS$Gen_Name=="CH NARA", ] #OK

#====================================================================================== -

#heading_GDDAS ----

#convert to dates
day_temp$day <- as.Date(day_temp$day, "%Y-%m-%d")

#calculate headingGDDAS for each plot
heading_GDDAS <- NULL
for (i in 1:nrow(HeadingDAS)) {
  if (!is.na(HeadingDAS$heading_date)[i]) {
    heading_GDDAS[i] <- day_temp[day_temp$day == paste(HeadingDAS$heading_date)[i], 6]
  }
  else {heading_GDDAS[i] <- NA}
  print(heading_GDDAS[i])
}

HeadingDAS$heading_GDDAS <- round(heading_GDDAS, 2)

# #====================================================================================== -
# 
# write.csv(HeadingDAS, paste0(dirto, "Heading_FPWW022_DAS_ANOVA.csv"), row.names = FALSE)
# 
# #====================================================================================== -

#Only one Lot was assessed:
#add heading dates to Lot 6,
#for checks, use the average heading date recorded in Lot 2
checks <- c("SURETTA", "CH CLARO", "CH NARA") #define which cultivars represent checks
d <- HeadingDAS[c(1:378),c(5:8)] %>% filter(!Gen_Name %in% checks) %>% 
  add_row(., Gen_Name = "SURETTA", HeadingDAS = 219, heading_date = as.Date("2018-05-24", "%Y-%m-%d"), heading_GDDAS = 1279.45) %>% 
  add_row(., Gen_Name = "CH NARA", HeadingDAS = 219, heading_date = as.Date("2018-05-24", "%Y-%m-%d"), heading_GDDAS = 1279.45) %>% 
  add_row(., Gen_Name = "CH CLARO", HeadingDAS = 218, heading_date = as.Date("2018-05-23", "%Y-%m-%d"), heading_GDDAS = 1260.91)

##add derived mean values to the df again
data <- right_join(d, HeadingDAS[1:5]) %>% 
  select(Plot_ID, Lot, RangeLot, RowLot, Gen_Name, heading_date, HeadingDAS, heading_GDDAS)

#====================================================================================== -

write.csv(data, paste0(base_path, "data_prep/Heading_FPWW022_DAS.csv"), row.names = FALSE)

#====================================================================================== -