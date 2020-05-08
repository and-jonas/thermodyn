
rm(list = ls())
.libPaths("T:R3UserLibs")
library(SpATS)
detach("package:plyr")
library(tidyverse)
library(asreml)
library(desplot)
workdir <- "O:/Projects/KP0011/2/"
setwd(workdir)

#functions
source("thermodyn/R/Utils/002_thermo_utils.R")
source("thermodyn/R/Utils/003_h2_BLUEs_utils.R")

#============================================================================================================== -

# dat <- readRDS("Data/data_ready/5_data_dyntraits.rds") ## use only selected subset
dat <- readRDS("Data/data_all_plots/5_data_dyntraits.rds") ## use all available data
design <- dat %>% dplyr::select(Plot_ID, design) %>% unnest(c(design)) %>% unique()

# reshape
newnames <- c(names(dat$dyntraits[[1]]), paste0(names(dat$dyntraits[[1]]), "_fc"))
data0 <- dat %>% 
  extract_covars_from_nested(from = "covariates", vars = c("heading_GDDAS", "Cnp_onsen_gom_GDDAS_fitted")) %>% 
  dplyr::select(Plot_ID, design, heading_GDDAS, Cnp_onsen_gom_GDDAS_fitted, dyntraits, dyntraits_fc) %>% 
  unnest(c(design, dyntraits, dyntraits_fc), names_repair = "unique") 
names(data0)[20:25] <- newnames

data0 <- data0%>% 
  gather(., dyntrait, value, slp_ct:sc_midstg_fc) %>% 
  unique()

#============================================================================================================== -

#perform spatial correction
#for dyntraits and dyntraits_fc
corrected_all <- data0 %>% 
  group_by(dyntrait, harvest_year) %>% group_nest() %>% 
  # filter(grepl("_fc", dyntrait)) %>% 
  # calculate within-year repeatablity
  mutate(w2 = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", 
                         fixed = "~ check + heading_GDDAS + Cnp_onsen_gom_GDDAS_fitted", genotype.as.random = TRUE, genotype = "Gen_Name",
                         .f = possibly(f_spats, otherwise = NA_real_)) %>%
           purrr::map_dbl(.x = .,
                          .f = possibly(get_h2, otherwise = NA_real_)))  %>%
  # extract BLUEs and spatially corrected plot values
  mutate(obj = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", 
                          fixed = "~ heading_GDDAS", genotype.as.random = FALSE, genotype = "Gen_Name",
                          .f = possibly(f_spats, otherwise = NA_real_))) %>%
  mutate(plot_obj = purrr::map(.x = obj, .f = plot)) %>% 
  mutate(BLUE =  purrr::map(.x = obj,
                            .f = possibly(get_BLUE_spats, otherwise = NA_real_))) %>%
  mutate(spat_corr = purrr::map(.x = obj, response = "value", element_return = "full",
                                .f = possibly(get_spat_corr_spats, otherwise = NA_real_))) %>% 
  #get spatial component
  mutate(sp = purrr::map(.x = obj,
                         .f = possibly(get_spatial, otherwise = NA_real_))) %>% 
  #plot the spatial trend 
  mutate(plot = purrr::map(.x = sp,  
                           form = formula(spatial ~ RangeLot + RowLot | Lot),
                           .f = possibly(desplot::desplot, otherwise = NA_real_)))

corrected_all$dyntrait
corrected_all$plot

#============================================================================================================== -

#extract corrected dyntraits (calculated from raw NRCT)
dyntraits_corr <- corrected_all %>% 
  dplyr::select(dyntrait, spat_corr) %>% unnest(spat_corr) %>% 
  spread(., dyntrait, spat_corr) %>% 
  dplyr::select(Plot_ID,  sc_midstg, sc_onsen, slp_ct) %>% 
  group_by(Plot_ID) %>% group_nest(.key = "dyntraits_corr")

#extract corrected dyntraits (calculated from spatially corrected NRCT)
dyntraits_fc_corr <- corrected_all %>% 
  dplyr::select(dyntrait, spat_corr) %>% unnest(spat_corr) %>% 
  spread(., dyntrait, spat_corr) %>% 
  dplyr::select(Plot_ID,  sc_midstg_fc, sc_onsen_fc, slp_ct_fc) %>% 
  #REVERT to standard names
  rename_at(vars(contains("_fc")), .funs = list(~gsub("_fc", "", .))) %>% 
  group_by(Plot_ID) %>% group_nest(.key = "dyntraits_fc_corr")

out <- dyntraits_corr %>% full_join(dyntraits_fc_corr) %>% 
  full_join(dat,.)

cor <- out %>% 
  extract_covars_from_nested("covariates", "heading_GDDAS") %>% 
  extract_covars_from_nested("dyntraits_fc_corr", "slp_ct") %>% 
  extract_covars_from_nested("design", "harvest_year") %>% 
  dplyr::select(Plot_ID, harvest_year, heading_GDDAS, slp_ct) %>% unique() %>% 
  dplyr::select(-Plot_ID) %>% group_by(harvest_year) %>% group_nest() %>% 
  mutate(cor = purrr::map_dbl(data, do_cor_test, x = "heading_GDDAS", y = "slp_ct"))


saveRDS(out, "Data/data_all_plots/6_data_dyntraits_corr.rds")

#============================================================================================================== -

# compare dyntraits extracted from raw NRCT and from spatially corrected NRCT
ding <- out %>% dplyr::select(Plot_ID, design, dyntraits, dyntraits_fc) %>% unique() %>% 
  unnest(c(design, dyntraits, dyntraits_fc), names_repair = "unique")
p1 <- ggplot(ding) + 
  geom_point(aes(x = slp_ct...18, y =slp_ct...21, color = as.factor(Lot))) +
  xlab("slp_ct_fr") + ylab("slp_ct_fc")

# compare dyntraits extracted from spatially corrected NRCT and dyntraits extracted from raw NRCT with subsequent spatial correction
ding <- out %>% dplyr::select(Plot_ID, design, dyntraits_fc, dyntraits_corr) %>% unique() %>% 
  unnest(c(design, dyntraits_fc, dyntraits_corr), names_repair = "unique")
p2 <- ggplot(ding) + 
  geom_point(aes(x = slp_ct...18, y =slp_ct...23, color = as.factor(Lot))) +
  xlab("slp_ct_fc") + ylab("slp_ct_fr_corr")

# compare dyntraits extracted from spatially corrected NRCT and dyntraits extracted from spatially corrected NRCT with subsequent repeated spatial correction
ding <- out %>% dplyr::select(Plot_ID, design, dyntraits_fc, dyntraits_fc_corr) %>% unique() %>% 
  unnest(c(design, dyntraits_fc, dyntraits_fc_corr), names_repair = "unique")
p3 <- ggplot(ding) + 
  geom_point(aes(x = slp_ct...18, y =slp_ct...23, color = as.factor(Lot))) +
  xlab("slp_ct_fc") + ylab("slp_ct_fc_corr")

#arrange plots and save
png("Figures/2year/dynpars_cor_corrmethods.png", width = 7, height = 5, units = 'in', res = 300)
gridExtra::grid.arrange(grobs=list(p1, p2, p3), ncol=2)
dev.off()

#====================================================================================== -

#TWO STAGE ----
#heritability from best linear unbiased estimators
h2_BLUE <- corrected_all %>%
  dplyr::select(dyntrait, harvest_year, BLUE) %>% unnest(BLUE) %>% 
  # add design
  # full_join(design, ., by = c("Gen_Name", "harvest_year")) %>% 
  mutate(Gen_Name = as.factor(Gen_Name),
         harvest_year = as.factor(harvest_year)) %>% 
  group_by(dyntrait) %>% group_nest() %>% 
  mutate(h2_2stage = purrr::map_dbl(data, possibly(get_h2_asreml2, otherwise = NA_real_), #after bug-fix, replace by map_dbl
                                  fixed = "BLUE ~ harvest_year",
                                  random = "~Gen_Name",
                                  residual = "~NULL",
                                  cullis = FALSE)) %>% 
  mutate(h2_2stage_c = purrr::map_dbl(data, possibly(get_h2_asreml2, otherwise = NA_real_), #after bug-fix, replace by map_dbl
                                    fixed = "BLUE ~ harvest_year",
                                    random = "~Gen_Name",
                                    residual = "~NULL",
                                    cullis = TRUE))

#====================================================================================== -

#ONE STAGE ----
#heritability from spatially corrected plot values
#scorings and agronomic traits
h2_spatcorr <- corrected_all %>%
  dplyr::select(dyntrait, harvest_year, spat_corr) %>% unnest(spat_corr) %>% 
  mutate(Gen_Name = as.factor(Gen_Name),
         harvest_year = as.factor(harvest_year)) %>% 
  group_by(dyntrait) %>% group_nest() %>% 
  mutate(h2_1stage_c = purrr::map_dbl(data, possibly(get_h2_asreml1, otherwise = NA_real_),
                                      fixed = "spat_corr ~ 1",
                                      random  = "~ harvest_year + Gen_Name + Gen_Name:harvest_year + Rep:at(harvest_year)",
                                      residual = "~dsum(~id(units) | harvest_year)",
                                      cullis = TRUE))


data <- h2_spatcorr$data[[6]]

data[is.na(data$spat_corr),] %>% nrow()

ggplot(data) +
  geom_histogram(aes(x = spat_corr)) +
  facet_wrap(~harvest_year, scales = "free")


#====================================================================================== -

## THESE ARE FINAL RESULTS 

#====================================================================================== -


#perform spatial correction
#for dyntraits and dyntraits_fc
corrected_all <- data0 %>% 
  group_by(dyntrait, harvest_year) %>% group_nest() %>% 
  filter(grepl("_fc", dyntrait)) %>% 
  # calculate within-year repeatablity
  mutate(w2_0 = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", 
                         fixed = "~ NULL", genotype.as.random = TRUE, genotype = "Gen_Name",
                         .f = possibly(f_spats, otherwise = NA_real_)) %>%
           purrr::map_dbl(.x = .,
                          .f = possibly(get_h2, otherwise = NA_real_))) %>% 
  mutate(w2_1 = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", 
                         fixed = "~ check", genotype.as.random = TRUE, genotype = "Gen_Name",
                         .f = possibly(f_spats, otherwise = NA_real_)) %>%
           purrr::map_dbl(.x = .,
                          .f = possibly(get_h2, otherwise = NA_real_))) %>% 
  mutate(w2_21 = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", 
                         fixed = "~ check + heading_GDDAS", genotype.as.random = TRUE, genotype = "Gen_Name",
                         .f = possibly(f_spats, otherwise = NA_real_)) %>%
           purrr::map_dbl(.x = .,
                          .f = possibly(get_h2, otherwise = NA_real_))) %>% 
  mutate(w2_22 = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", 
                         fixed = "~ check + Cnp_onsen_gom_GDDAS_fitted", genotype.as.random = TRUE, genotype = "Gen_Name",
                         .f = possibly(f_spats, otherwise = NA_real_)) %>%
           purrr::map_dbl(.x = .,
                          .f = possibly(get_h2, otherwise = NA_real_))) %>% 
  mutate(w2_3 = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", 
                           fixed = "~ check + heading_GDDAS + Cnp_onsen_gom_GDDAS_fitted", genotype.as.random = TRUE, genotype = "Gen_Name",
                           .f = possibly(f_spats, otherwise = NA_real_)) %>%
           purrr::map_dbl(.x = .,
                          .f = possibly(get_h2, otherwise = NA_real_)))

W2 <- corrected_all %>% dplyr::select(-data)


