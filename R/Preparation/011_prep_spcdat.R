

rm(list = ls())
.libPaths("T:/R3UserLibs")
library(SpATS)
library(tidyverse)

workdir <- "O:/Projects/KP0011/2/"
setwd(workdir)

#functions
source("thermodyn/R/Utils/001_spectra_utils.R")
source("thermodyn/R/Utils/002_thermo_utils.R")
source("thermodyn/R/Utils/003_h2_BLUEs_utils.R")

#============================================================================================================== -

#extract BLUE from SpATS
spat_corr_spats <- function(
  object, # a fitted SpATS object
  response, #the modelled trait
  element_return = "minimal") #either "minimal" (only corrected values per plot) or "full" (including design) 
{
  fitted <- object$fitted
  intercept <- object$coeff['Intercept']
  gen_mod_mat <- construct_genotype_prediction_matrix(object, object$data)
  gen_coeff <- object$coeff[1:ncol(gen_mod_mat)]
  geno_pred <- as.vector(gen_mod_mat %*% gen_coeff)
  residual <- object$residuals
  if(element_return == "full"){
    plot_corr <- as.data.frame(intercept + geno_pred + residual) %>% dplyr::rename(spat_corr = 1) %>% 
      bind_cols(object$data, .) %>% as_tibble() %>% dplyr::select(-one_of(response), -weights)
  } else if(element_return == "minimal"){
    plot_corr <- as.data.frame(intercept + geno_pred + residual) %>% dplyr::rename(spat_corr = 1) %>% 
      as_tibble()
  }
}

#============================================================================================================== -

#load data ----

#load and filter spectral data
spc <- data.table::fread("O:/Projects/KP0011/Spectral_modelling/Data/spc_raw.csv") %>% 
  mutate(Plot_ID = as.factor(Plot_ID),
         meas_date = as.Date(meas_date),
         replicate = as.factor(replicate)) %>% 
  as_tibble() %>% 
  #only FPWW022
  filter(grepl("^FPWW022", Plot_ID))

#load temperature data and heading dates
day_temp <- read.csv("Data/data_prep/gdd_2018.csv") %>% 
  mutate(day = as.Date(day, "%Y-%m-%d"))
heading <- readRDS("Data/data_ready/1_design_phenodata.rds") %>% 
  dplyr::select(Plot_ID, heading_date, heading_GDDAS)

#load design
#design ----
design <- read_csv("Data/data_prep/exp_design.csv") %>% 
  as_tibble() %>% 
  mutate_at(vars(Lot, RangeLot, RowLot, Range, Row, RowBL, RangeBL), funs(as.numeric)) %>% 
  mutate_at(vars(Gen_ID, Rep, Xf, Yf), funs(as.factor))

#============================================================================================================== -

#Prep SVI dataset ----
#Calculate SVI
SI <- spc %>%
  #smooth spectra using the Savitzky-Golay smoothing filter
  f_spc_smooth(3, 11, 0) %>%
  #average spectra for each plot
  f_spc_avg() %>% 
  #calculate spectral vegetation indices
  f_calc_si()

#add heading date
data <- left_join(SI, heading, by = "Plot_ID") %>% dplyr::select(-starts_with("SI_"), everything())

#Convert measurement dates to meas_GDDAH and meas_GDDAS
##GDDAH
meas_GDDAH <- NULL
for (i in 1:nrow(data)) {
  if (!is.na(data$heading_date)[i]) {
    meas_GDDAH[i] <- day_temp[day_temp$day == paste(data$meas_date)[i], 6] -
      day_temp[day_temp$day == paste(data$heading_date)[i], 6]
  }
  else {
    meas_GDDAH[i] <- paste("NA")
  }
  print(meas_GDDAH[i])
}
data$meas_GDDAH <- as.numeric(round(as.numeric(meas_GDDAH), 0)) 

##GDDAS
meas_GDDAS <- NULL
for(i in 1:nrow(data)){
  if(!is.na(data$heading_date)[i]){
    meas_GDDAS[i] <- as.numeric(data$heading_GDDAS[i]) + as.numeric(data$meas_GDDAH[i])
  }
  else {
    meas_GDDAS[i] <- day_temp[day_temp$day == paste(data$meas_date)[i], 6]
  }
  print(meas_GDDAS[i])
}
data$meas_GDDAS <- as.numeric(round(as.numeric(meas_GDDAS), 0))

#rearrange output
SI <- data %>% dplyr::select(-starts_with("SI_"), everything())

SIready <- SI %>% 
  #select a subset of SI
  dplyr::select(Plot_ID, meas_date, meas_GDDAH, meas_GDDAS, matches("SI_PSRI|SI_NDVI_nb_ASD|MCARI2")) %>%
  #INVERT PSRI VALUES TO HELP INTERPRETATION
  mutate(SI_PSRI = -SI_PSRI) %>% 
  #select only first five measurements
  group_by(Plot_ID) %>% slice(1:5) %>% group_by(Plot_ID, meas_date) %>% 
  group_nest(.key = "Spc")

#============================================================================================================== -

#spat corr ----

corr <- SI %>% 
  dplyr::select(Plot_ID, meas_date, meas_GDDAH, meas_GDDAS, starts_with("SI_")) %>% 
  right_join(design,. ) %>% 
  gather(., SI, value, starts_with("SI_")) %>% 
  filter(meas_date < as.Date("2018-06-26")) %>%
  filter(grepl("SI_PSRI|SI_NDVI_nb_ASD|MCARI2", SI)) %>% 
  group_by(SI, meas_date) %>% group_nest() %>% 
  # calculate within-year repeatablity
  mutate(w2 = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", 
                         fixed = "~ NULL", genotype.as.random = TRUE, genotype = "Gen_Name",
                         .f = possibly(f_spats, otherwise = NA_real_)) %>%
           purrr::map_dbl(.x = .,
                          .f = possibly(get_h2, otherwise = NA_real_)))  %>%
  # extract BLUEs and spatially corrected plot values
  mutate(obj = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", 
                          fixed = "~ NULL", genotype.as.random = FALSE, genotype = "Gen_Name",
                          .f = possibly(f_spats, otherwise = NA_real_))) %>%
  mutate(BLUE =  purrr::map(.x = obj,
                            .f = possibly(get_BLUE_spats, otherwise = NA_real_))) %>%
  mutate(spat_corr = purrr::map(.x = obj, response = "value", element_return = "minimal",
                                .f = possibly(spat_corr_spats, otherwise = NA_real_))) %>% 
  #get spatial component
  mutate(sp = purrr::map(.x = obj,
                         .f = possibly(get_spatial, otherwise = NA_real_))) %>% 
  #plot the spatial trend 
  mutate(plot = purrr::map(.x = sp,  
                           form = formula(spatial ~ RangeLot + RowLot | Lot),
                           .f = possibly(desplot::desplot, otherwise = NA_real_)))

corr$plot

post <- corr %>% dplyr::select(meas_date, SI, data, spat_corr) %>% 
  unnest(c(data, spat_corr))
#============================================================================================================== -

#create output ----

pp <- post %>% dplyr::select(Plot_ID, SI, meas_date, meas_GDDAH, meas_GDDAS, spat_corr) %>% 
  spread(., SI, spat_corr) %>% 
  group_by(Plot_ID, meas_date) %>% group_nest(.key = "Spc_corr")

SPC <- full_join(SIready, pp)

saveRDS(SPC, "Data/data_ready/A_spc.rds")

#============================================================================================================== -

#cor raw ~ corr ----
spc <- readRDS("Data/A_spc.rds") 
Spc <- spc %>% dplyr::select(Plot_ID, meas_date, Spc) %>% unnest(Spc) %>% 
  gather(SI, val, starts_with("SI_"))
Spcc <- spc %>% dplyr::select(Plot_ID, meas_date, Spc_corr) %>% unnest(Spc_corr) %>% 
  gather(., SI, val_c, starts_with("SI"))

all <- full_join(Spc, Spcc) %>% 
  filter(!meas_date == "2018-06-27") #wrf???
ggplot(all) +
  geom_point(aes(x = val, y = val_c)) +
  facet_wrap(~interaction(SI, meas_date), scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank())

#============================================================================================================== -
