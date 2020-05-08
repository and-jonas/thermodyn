
#HEADER ----

rm(list = ls())
.libPaths("T:/R3UserLibs")
library(data.table)
library(tidyverse)
base_path <- "O:/Projects/KP0011/2/"
setwd(base_path)

source("thermodyn/R/Utils/001_spectra_utils.R")

############################################################################################# -

#design ----
exp_design <- read_csv("Data/data_prep/exp_design.csv") %>% 
  as_tibble() %>% 
  mutate_at(vars(Lot, RangeLot, RowLot, Range, Row, RowBL, RangeBL), funs(as.numeric)) %>% 
  mutate_at(vars(Gen_ID, Rep, Xf, Yf), funs(as.factor))

#heading ----
head <- read_csv("Data/data_prep/heading.csv") %>% # heading date
  dplyr::select(Plot_ID, heading_date, HeadingDAS, heading_GDDAS) %>% 
  mutate(heading_date = as.Date(heading_date),
         heading_DAS = as.numeric(HeadingDAS)) %>% as_tibble()

#onsen ----
onsen <- read_csv("Data/data_prep/dynpars_scr.csv") %>% 
  dplyr::select(Plot_ID, level, contains("onsen_gom"), contains("onsen_lin")) %>% 
  gather(par, value, contains("onsen")) %>% 
  mutate(Trait = paste(level, par, sep = "_")) %>% 
  dplyr::select(-par, -level) %>% 
  spread(Trait, value)

#gy ----
grnyld22 <- read.csv("Data/data_from_db/FPWW022_gy.csv") %>% as_tibble() %>% 
  dplyr::select(plot.UID, trait_value.value) %>% 
  rename("Plot_ID" = plot.UID,
         "GY" = trait_value.value)
grnyld24 <- read.csv("Data/data_from_db/FPWW024_lot2_lot4_gy.csv") %>% as_tibble() %>% 
  dplyr::select(plot.UID, trait_value.value) %>% 
  rename("Plot_ID" = plot.UID,
         "GY" = trait_value.value)
grnyld28 <- read.csv("Data/data_from_db/FPWW028_lot2_lot4_gy.csv") %>% as_tibble() %>% 
  dplyr::select(plot.UID, trait_value.value) %>% 
  rename("Plot_ID" = plot.UID,
         "GY" = trait_value.value)
grnyld <- bind_rows(grnyld22, grnyld24, grnyld28)

#gpc ----
grnprt <- read.csv("O:/FIP/2018/WW022/RefTraits/FPWW022_protein_NIRS.csv") %>% as_tibble() %>% 
  dplyr::select(plot_label, value) %>% 
  dplyr::rename("Plot_ID" =  plot_label, 
                "GPC" = value)

#tkw ----
tkw <- read_csv("O:/FIP/2018/WW022/RefTraits/FPWW022_marvin_TKW.csv") %>% 
  mutate(TKW = TKM.g./100) %>% 
  dplyr::select(Kennung, TKW) %>% 
  rename("Plot_ID" = Kennung)#gy, gpc ----

#glc ---- 
glc1 <- readxl::read_excel("O:/FIP/2017/WW018/RefTraits/Fl0Glc_20170614.xlsx", skip = 3) %>% 
  dplyr::select(Plot_ID, Fl0Glc)
glc2 <- read.csv("O:/FIP/2018/WW022/RefTraits/GlcEye/glaucousness_2018-06-20.csv") %>% 
  dplyr::select(plot_label, value) %>% 
  dplyr::rename(Plot_ID = 1,
                Fl0Glc = 2)
glc <- rbind(glc1, glc2)

#angle ----
angle <- read.csv("O:/FIP/2018/WW022/RefTraits/Leaf_Angle/2018-06-18-02-08-20_171002_Design_FPWW022_FieldBook_database.csv", skip = 3) %>% 
  dplyr::select(Plot_ID, value) %>% 
  dplyr::rename(angle = 2) %>% 
  mutate(angle = as.numeric(angle)) %>% 
  filter(!grepl("Daten wurden gespeichert", angle))

#Fl0Wid, Fl0Len ----
widlen <- read.csv("O:/FIP/2018/WW022/RefTraits/Wid_Len/2018-06-29-09-07-23_FPWW022_FieldBook_Plt_Rul_Fl0WidLen_database.csv") %>% 
  mutate(value = as.numeric(as.character(value))) %>% as_tibble() %>% 
  group_by(Plot_ID, Trait_ID) %>% dplyr::summarise(mean = mean(value, na.rm = TRUE)) %>% 
  spread(Trait_ID, mean)

#  data ----
spc <- data.table::fread("O:/Projects/KP0011/Spectral_modelling/Data/spc_raw.csv") %>% 
  mutate(Plot_ID = as.factor(Plot_ID),
         meas_date = as.Date(meas_date),
         replicate = as.factor(replicate)) %>% 
  as_tibble()
spc <- spc %>% filter(grepl("^FPWW022", Plot_ID))

#GDDAH for measurements and scorings (plot-specific)
gddah <- read.csv("O:/Projects/KP0011/1/Data_Andereggetal2019/Data/helper_files/gddah_data.csv") %>% 
  mutate(meas_date = meas_date %>% as.Date()) %>% 
  as_tibble()

#lodging plots, to be excluded from the analysis
lodg_plots <- read.csv("O:/Projects/KP0011/1/Data_Andereggetal2019/Data/other_data/lodging.csv", sep = ";") %>% 
  filter(Lodging %in% c("l", "ll")) %>% pull(Plot_ID)

#remove lodging plots
spc <- spc %>% 
  filter(!Plot_ID %in% lodg_plots)

#create dataset to extract dynamics parameters
##This dataset comprises all measurement dates available,
##except for those carried out prior to heading
SI <- spc %>%
  #smooth spectra using the Savitzky-Golay smoothing filter
  f_spc_smooth(3, 11, 0) %>%
  #average spectra for each plot
  f_spc_avg() %>% 
  #calculate spectral vegetation indices
  f_calc_si() %>%
  #add gddah data
  left_join(., gddah, by = c("Plot_ID", "meas_date")) %>% 
  #rearrange
  dplyr::filter(Exp == "FPWW022") %>% 
  filter(meas_date == as.Date("2018-06-23")) %>% 
  dplyr::select(Plot_ID, SI_1200)

#join by Plot_ID, where plot-level data is available
## MUST BE 2,232 x 30
d <- list(exp_design, head, onsen, grnyld, grnprt, tkw, glc, angle, widlen, SI) %>% 
  Reduce(function(df1, df2) left_join(df1, df2, by = "Plot_ID"),.)

des_names <- names(exp_design)
names <- c(des_names) %>% unique()
dd <- d %>% dplyr::select(one_of(names), everything())

############################################################################################# -

#JOIN BY GEN_NAME OR GEN_ID

#load awn data
awn <- readxl::read_excel("O:/Projects/KP0011/1/Data preparation/Data/raw/Additional_PhenInfo_ww012_2.xls") %>% 
  dplyr::select(Gen, Awns) %>% 
  mutate(Awns = ifelse(is.na(Awns), 0, 1) %>% as.factor()) %>% unique()

#load final height data
load("O:/Projects/KP0011/Spectral_modelling/AllDFList.RDA")
hgt <- AllDFList$BLUEs3Y[,c("Gen_ID","FH")] %>% as_tibble() %>% 
  mutate(Gen_ID = as.factor(Gen_ID))

#join by Gen_ID or Gen_Name, where genotype-level data is available
Data <- left_join(dd, awn, by = c("Gen_Name" = "Gen")) %>% 
  left_join(., hgt, by = "Gen_ID")

saveRDS(Data, "Data/data_ready/1_design_phenodata.rds")

############################################################################################# -

Data <- readRDS("Data/data_ready/1_design_phenodata.rds")

# THERM data ----
#>2018 ----
##experimental plots
dir <- "O:/Projects/KP0011/2/Data/Thermal_data/"
files <- as.list(list.files(dir, pattern = "^FPWW022", full.names = TRUE))
d0 <- files %>% lapply(., readr::read_csv)
data <- d0 %>% bind_rows(.) %>% 
  # dplyr::select(UID, timestamp, percentile_050) %>% #UID does not exist in this file anymore ?!
  dplyr::select(plot_label, timestamp, percentile_050) %>% rename("UID" = plot_label) %>% 
  mutate(temp = percentile_050/100) %>% dplyr::select(-percentile_050) %>% 
  mutate(temp = ifelse(temp == 0, NA, temp)) %>% 
  mutate(meas_date = as.Date(timestamp, "%Y-%m-%d"))

#standards
files <- as.list(list.files(dir, pattern = "^REF", full.names = TRUE))
d1 <- files %>% lapply(., readr::read_csv)
std <- d1 %>% bind_rows(.) %>% 
  # dplyr::select(UID, timestamp, percentile_050) %>% #UID does not exist in this file anymore ?!
  dplyr::select(plot_label, timestamp, percentile_050) %>% rename("UID" = plot_label) %>% 
  mutate(temp = percentile_050/100) %>% dplyr::select(-percentile_050) %>% 
  mutate(temp = ifelse(temp == 0, NA, temp)) %>% 
  mutate(meas_date = as.Date(timestamp, "%Y-%m-%d"))

#get bush reference average temperature
ref_bush <- std %>% filter(grepl("bush", UID)) %>% 
  group_by(timestamp) %>% 
  summarize(mean_temp = mean(temp, na.rm = TRUE),
            stddev = sd(temp)) #should be OK (?)                   <==== -
ref_bush <- ref_bush %>% transmute(timestamp, mean_bsh_temp = mean_temp)

#get street reference average temperature
ref_street <- std %>% filter(grepl("street", UID)) %>% 
  group_by(timestamp) %>% 
  summarize(mean_temp = mean(temp, na.rm = TRUE),
            stddev = sd(temp)) #should be OK (?)                   <==== -
ref_street <- ref_street %>% transmute(timestamp, mean_str_temp = mean_temp)

normalized <- data %>% 
  #add reference termperature
  full_join(., ref_bush, by = "timestamp") %>% 
  full_join(., ref_street, by = "timestamp") %>% 
  #calculate deltaT
  mutate(deltaT_bsh = temp - mean_bsh_temp,
         deltaT_str = temp - mean_str_temp) %>% 
  rename("Plot_ID" = UID)

############################################################################################# -

## >2019 ----
trait_th_path <- paste0("O:/Projects/KP0011/2/Data/Masterarbeit/", "01_Data/CroPyDB/Thermal/")

df_traits_th <- data.frame()
files <- list.files(path = trait_th_path, pattern="*.csv")
for (file in files) {
  df_trait_th <- read_csv(paste0(trait_th_path, file))
  df_traits_th <- bind_rows(df_traits_th, df_trait_th)
  remove(df_trait_th)
}

df_traits_th <- df_traits_th %>% mutate(date = date(as.POSIXlt(timestamp, tz="Etc/GMT-1")))

#extract percentiles values from json
df_traits_th <- df_traits_th %>% mutate(value_json = gsub("'", '"', value_json))
df_traits_th <- df_traits_th %>%  
  mutate(json = purrr::map(value_json, fromJSON)) %>% 
  mutate(json = purrr::map(json, as.data.frame)) %>%
  unnest(json)

unique(df_traits_th$date)

data <- df_traits_th %>%  
  dplyr::select(plot.UID, timestamp, percentile_050) %>% 
  mutate(temp = percentile_050/100) %>% dplyr::select(-percentile_050) %>% 
  mutate(temp = ifelse(temp == 0, NA, temp)) %>% 
  mutate(meas_date = as.Date(timestamp, "%Y-%m-%d", tz="Etc/GMT-1")) %>% 
  rename("Plot_ID" = plot.UID)

#combine years 
data_therm <- bind_rows(normalized, data) %>% arrange(Plot_ID, timestamp)

############################################################################################# -

#Merge with covariate data
dat <- left_join(Data, data_therm, by = "Plot_ID")

############################################################################################# -

#>meas_gddah, meas_gddas ----

#Temperature data
temp <- read.csv("Data/data_prep/gdd.csv")

#unique heading_date x meas_date combinations
sub <- dat %>% dplyr::select(heading_date, meas_date) %>% unique()

#Add meas_GDDAH 
meas_GDDAH <- NULL
for (i in 1:nrow(sub)) {
  if (!is.na(sub$heading_date)[i]) {
    meas_GDDAH[i] <- temp[temp$day == paste(sub$meas_date)[i], 6] -
      temp[temp$day == paste(sub$heading_date)[i], 6]
  }
  else meas_GDDAH[i] <- NA
  print(meas_GDDAH[i])
}
sub$meas_GDDAH <- meas_GDDAH
dd <- dat %>% full_join(., sub, by = c("heading_date", "meas_date")) %>% 
  arrange(Plot_ID)

#Add meas_GDDAS
sub <- dd %>% dplyr::select(meas_date) %>% unique()
meas_GDDAS <- NULL
for (i in 1:nrow(sub)) {
  meas_GDDAS[i] <- temp[temp$day == paste(sub$meas_date)[i], 6]
}
sub$meas_GDDAS <- meas_GDDAS
ddd <- dd %>% full_join(., sub, by = c("meas_date"))

############################################################################################# -

#create nested tibble

design <- ddd %>% dplyr::select(one_of(names(exp_design))) %>% unique() %>% 
  group_by(Plot_ID) %>% group_nest(.key = "design")

covars <- grep("Exp|Plot_ID|^heading|^Cnp|^Fl0|angle|Awns|FH|GY|GPC|TKW|SI_1200", names(ddd), value = TRUE)

covars <- ddd %>% dplyr::select(one_of(covars)) %>% unique() %>%
  group_by(Plot_ID) %>% group_nest(.key = "covariates")

meas <- ddd %>% dplyr::select(one_of(c(names(normalized), "meas_GDDAH", "meas_GDDAS", "Plot_ID"))) %>% 
  group_by(Plot_ID, timestamp) %>% group_nest(.key = "thermodata")
  
FULL <- full_join(design, covars, by = "Plot_ID") %>% 
  full_join(., meas, by = "Plot_ID")
saveRDS(FULL, "Data/data_ready/2_thermodata.rds")

############################################################################################# -                        