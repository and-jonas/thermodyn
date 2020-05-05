
rm(list = ls())
.libPaths("T:R3UserLibs")
library(SpATS)
detach("package:plyr")
library(tidyverse)

workdir <- "O:/Projects/KP0011/2/"
setwd(workdir)

#functions
source("thermodyn/R/Utils/003_h2_BLUEs_utils.R")
source("thermodyn/R/Utils/002_thermo_utils.R")

#============================================================================================================== -

#read data
dat <- readRDS("Data/data_ready/2_thermodata.rds")
covariates <- dplyr::select(dat, Plot_ID, covariates)
design <- dplyr::select(dat, Plot_ID, design) %>% unique()

#define covariates to use 
vars <- names(dat$covariates[[1]])[c(3, 5, 11, 7, 13, 9, 17:26)]

#select best flights
cs <- c("2018-06-04 12:18:34", "2018-06-07 13:47:57", 
        "2018-06-16 13:27:18", "2018-06-20 11:59:41", 
        "2018-06-23 13:08:47", "2018-06-26 11:48:37", 
        "2018-06-30 09:01:21", "2018-07-04 11:35:08", 
        "2018-07-09 11:15:07", "2018-07-12 12:29:44",
        grep("^2019", unique(dat$timestamp), value = TRUE))
cs <- paste(cs, collapse = "|")
data0 <- dat %>% filter(grepl(cs, timestamp)) 

#============================================================================================================== -

#reshape data
d1.1 <- data0 %>%   
  extract_covars_from_nested("thermodata", "temp") %>% 
  dplyr::select(Plot_ID, design, timestamp, temp) %>% 
  #modifications to allow processing in one step
  mutate(timestamp = as.factor(timestamp)) %>%
  rename(value = temp) %>% 
  unnest(design) %>% 
  group_by(timestamp) %>% group_nest() %>% 
  rename(trait = timestamp)
d1.2 <- data0 %>% 
  extract_covars_from_nested("covariates", vars) %>%
  dplyr::select(Plot_ID, design, vars) %>% 
  unique() %>% 
  unnest(design) %>% 
  gather(., trait, value, vars) %>% 
  group_by(trait) %>% group_nest()
newd <- rbind(d1.1, d1.2) 

#============================================================================================================== -

#perform spatial correction
corrected_all <- newd %>% 
  # calculate within-year repeatablity
  mutate(w2 = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", 
                         fixed = "~ check", genotype.as.random = TRUE, genotype = "Gen_Name",
                         .f = possibly(f_spats, otherwise = NA_real_)) %>%
           purrr::map_dbl(.x = .,
                          .f = possibly(get_h2, otherwise = NA_real_)))  %>%
  # extract BLUEs and spatially corrected plot values
  mutate(obj = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", 
                          fixed = "~ NULL", genotype.as.random = FALSE, genotype = "Gen_Name",
                          .f = possibly(f_spats, otherwise = NA_real_))) %>%
  mutate(BLUE =  purrr::map(.x = obj,
                            .f = possibly(get_BLUE_spats, otherwise = NA_real_))) %>%
  mutate(spat_corr = purrr::map(.x = obj, response = "value", element_return = "full",
                                .f = possibly(get_spat_corr_spats, otherwise = NA_real_)))

W2 <- corrected_all %>% dplyr::select(trait, w2)

#============================================================================================================== -

#extract corrected covars
covariates_corr <- corrected_all %>% 
  #exclude thermal measurements
  dplyr::filter(!grepl("^20", trait)) %>% 
  dplyr::select(trait, spat_corr) %>% unnest(spat_corr) %>% 
  spread(., trait, spat_corr) %>% 
  dplyr::select(Plot_ID, vars) %>% 
  group_by(Plot_ID) %>% group_nest(.key = "covariates_corr")

#extract corrected temps
temp_corr <- corrected_all %>% 
  #select thermal measurements
  dplyr::filter(grepl("^20", trait)) %>% 
  dplyr::select(trait, spat_corr) %>% unnest(spat_corr) %>% 
  ungroup() %>%
  mutate(trait = as.POSIXct(trait, "%Y-%m-%d %H:%M:%OS")) %>%
  group_by(trait) %>%
  dplyr::select(trait, Plot_ID, spat_corr) %>% 
  rename("temp" = spat_corr) %>% 
  arrange(Plot_ID)
meas_GDDAH <- data0 %>% extract_covars_from_nested("thermodata", c("meas_GDDAH", "meas_GDDAS")) %>%
  dplyr::select(meas_GDDAH, meas_GDDAS)
temp_corr <- bind_cols(temp_corr, meas_GDDAH) %>% 
  group_by(trait, Plot_ID) %>% group_nest(.key = "thermodata_corr") %>% 
  rename("timestamp" = trait)

#add corrected traits to data
DD <- data0 %>% 
  full_join(., covariates_corr) %>% 
  full_join(., temp_corr)
  
#rearrange
out <- DD %>% dplyr::select(1:3, 6, 4, 5, 7) 

#save
saveRDS(out, "Data/data_ready/3_data.rds")

#============================================================================================================== -

# corr vs. raw ----

d <- out %>% dplyr::select(Plot_ID, design, thermodata, thermodata_corr) %>% 
  unnest(c(design, thermodata, thermodata_corr), names_repair = "unique")

plot <- ggplot(d) + 
  geom_point(aes(x = temp...18, y = temp...26, color = as.factor(Lot))) +
  xlab("CT_raw") + ylab("CT_corr")+
  facet_wrap(~meas_date, scales = "free") #OK

png("Figures/2year/02_CT_raw_corr.png", units = "in", width = 10, height = 8, res = 300)
plot(plot)
dev.off()

#============================================================================================================== -