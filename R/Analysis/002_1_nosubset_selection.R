
rm(list = ls())
.libPaths("T:/R3UserLibs")
library(tidyverse)
library(plotly)
library(ggplot2)
library(ggpmisc)
library(scales)
library(data.table)

workdir <- "O:/Projects/KP0011/2/"
setwd(workdir)

source("thermodyn/R/Utils/001_spectra_Utils.R")
source("thermodyn/R/Utils/002_thermo_utils.R")

#==================================================================================== -

dat <- readRDS("Data/data_ready/3_data.rds")
design <- dplyr::select(dat, Plot_ID, design) %>% unique()

#Selection on Heading/Onsen ----

#>Compare LP01 and GABI ----
expanded <- dat %>% 
  extract_covars_from_nested("design", vars = c("Exp", "Set", "harvest_year")) %>% 
  extract_covars_from_nested("covariates", vars = c("heading_GDDAS",
                                                    "Cnp_onsen_gom_GDDAS_fitted", "Fl0_onsen_gom_GDDAS_fitted",
                                                    "Cnp_onsen_gom_GDDAH", "Cnp_onsen_gom_GDDAS"))
FPWW022 <- expanded[expanded$Exp == "FPWW022", c("Plot_ID", "heading_GDDAS", 
                                                "Cnp_onsen_gom_GDDAS_fitted", "Cnp_onsen_gom_GDDAH",
                                                "Exp", "Set", "harvest_year")]
FPWW024028 <- expanded[expanded$Exp == "FPWW024"|expanded$Exp == "FPWW028",c("Plot_ID", "heading_GDDAS", 
                                                                             "Cnp_onsen_gom_GDDAS_fitted", "Cnp_onsen_gom_GDDAH", 
                                                                             "Exp", "Set", "harvest_year")]
FPWW <- bind_rows(FPWW022, FPWW024028)

plotdat_hd <- FPWW[c("Plot_ID", "heading_GDDAS", "Set", "harvest_year")] %>% unique()
plothd <- ggplot(plotdat_hd) +
  geom_histogram(aes(x = heading_GDDAS, fill = Set)) +
  facet_wrap(~harvest_year)

plotdat_sen <- FPWW[c("Plot_ID", "Cnp_onsen_gom_GDDAS_fitted", "Cnp_onsen_gom_GDDAH", "Set", "harvest_year")] %>% unique()
plotsen <- ggplot(plotdat_sen) +
  geom_histogram(aes(x = Cnp_onsen_gom_GDDAS_fitted, fill = Set)) +
  facet_wrap(~harvest_year)

plotdat_sen_gddah <- FPWW[c("Plot_ID", "Cnp_onsen_gom_GDDAH", "Set", "harvest_year")] %>% unique()
plotsen_gddah <- ggplot(plotdat_sen) +
  geom_histogram(aes(x = Cnp_onsen_gom_GDDAH, fill = Set)) +
  facet_wrap(~harvest_year)

png("Figures/2year/hd_lp01_gabi.png", units = "in", width = 10, height = 8, res = 300)
plot(plothd)
dev.off()

png("Figures/2year/sen_lp01_gabi.png", units = "in", width = 10, height = 8, res = 300)
plot(plotsen)
dev.off()

png("Figures/2year/stg_lp01_gabi.png", units = "in", width = 10, height = 8, res = 300)
plot(plotsen_gddah)
dev.off()

#==================================================================================== -
#==================================================================================== -

#>Select subset for each year ----

expanded <- expanded %>% 
  extract_covars_from_nested("thermodata", vars = c("meas_GDDAS", "meas_GDDAH"))

#select flights
dates19 <- c("2019-06-13 12:14:00",
             "2019-06-18 11:38:00",
             "2019-06-25 11:28:00",
             "2019-06-27 11:26:00",
             "2019-06-29 11:31:00",
             "2019-07-01 11:47:00") %>% 
  paste(., collapse = "|")
selected <- expanded %>% filter(timestamp < "2018-06-26 11:48:37 UTC" | grepl(dates19, timestamp))

data_selected <- selected %>% ungroup() %>% dplyr::select(Plot_ID:thermodata_corr)
saveRDS(data_selected, "Data/data_all_plots/4_data_selected.rds")

#==================================================================================== -
#==================================================================================== -

#>Select subset phen + stg ----

year_list <- split(expanded, expanded$harvest_year)

selected <- list()
for (i in 1:length(year_list)){
  #get data
  data <- year_list[[i]]
  #get LP01 range in heading
  range_hd <- data %>% dplyr::filter(Set == "LP01") %>% 
    dplyr::select(heading_GDDAS) %>% range(., na.rm = TRUE)
  #select subset for phenology
  subset_phen <- data %>% dplyr::filter(Set == "LP01" | heading_GDDAS %inrange% range_hd)
  #select subset for "sufficient" stay-green (i.e. at least 5 flights prior to onset)
  subset_phen_stg <- subset_phen %>% filter(meas_GDDAS < Cnp_onsen_gom_GDDAS_fitted) %>%
    arrange(Plot_ID) %>% group_by(Plot_ID) %>% group_nest() %>% 
    filter(map_int(data, nrow) >= 5) %>% unnest(c(data)) %>% 
    group_by(timestamp)
  #get number of flights meeting the criterium at each flight
  nplots <- subset_phen_stg %>% nest() %>% mutate(plots = map_dbl(data, nrow)) 
  #get flights with at least 90% of the plots pre-senescence
  flights <- nplots[which(nplots$plots > 0.9*max(nplots$plots)), ] %>% pull(timestamp)
  #exclude measurements carried out prior to end of heading
  prec_flights <- subset_phen_stg %>% slice(which.min(meas_GDDAH)) %>% dplyr::filter(meas_GDDAH < 0) %>% pull(timestamp)
  #FILTER FLIGHTS
  selected[[i]] <- subset_phen_stg %>% filter(timestamp %in% flights) %>% 
    filter(!timestamp %in% prec_flights)
} #five flights in 2018, 6 flights in 2019

data_selected <- bind_rows(selected) %>% ungroup() %>% dplyr::select(Plot_ID:thermodata_corr)
saveRDS(data_selected, "Data/data_ready/4_data_selected.rds")

#==================================================================================== -

#>Select subset phen ----

year_list <- split(selected, selected$harvest_year)

selected <- list()
for (i in 1:length(year_list)){
  #get data
  data <- year_list[[i]]
  #get LP01 range in heading
  range_hd <- data %>% dplyr::filter(Set == "LP01") %>% 
    dplyr::select(heading_GDDAS) %>% range(., na.rm = TRUE)
  #select subset for phenology
  subset_phen <- data %>% dplyr::filter(Set == "LP01" | heading_GDDAS %inrange% range_hd) %>% 
    group_by(timestamp)
  #get number of flights meeting the criterium at each flight
  nplots <- subset_phen %>% nest() %>% mutate(plots = map_dbl(data, nrow)) 
  #get flights with at least 90% of the plots pre-senescence
  flights <- nplots[which(nplots$plots > 0.9*max(nplots$plots)), ] %>% pull(timestamp)
  #exclude measurements carried out prior to end of heading
  prec_flights <- subset_phen %>% slice(which.min(meas_GDDAH)) %>% dplyr::filter(meas_GDDAH < 0) %>% pull(timestamp)
  #FILTER FLIGHTS
  selected[[i]] <- subset_phen %>% filter(timestamp %in% flights) %>% 
    filter(!timestamp %in% prec_flights)
} #five flights in 2018, 6 flights in 2019

data_selected <- bind_rows(selected) %>% ungroup() %>% dplyr::select(Plot_ID:thermodata_corr)
saveRDS(data_selected, "Data/data_ready/4_data_selected_pheno.rds")

#==================================================================================== -
#==================================================================================== -
