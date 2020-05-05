
rm(list = ls())
.libPaths("T:/R3UserLibs")
library(tidyverse)
workdir <- "O:/Projects/KP0011/2/"
setwd(workdir)

source("thermodyn/R/Utils/002_thermo_utils.R")

#==================================================================================== -

dat <- readRDS("Data/data_ready/3_data.rds")
design <- dplyr::select(dat, Plot_ID, design) %>% unique()

#analyzed covariates
covars <- c("FH", "Heading_GDDAS", "OnSen_GDDAH","Onsen_GDDAS",
            "Fl0Ang", "Fl0Len", "Fl0Wid", "Fl0Glc", "Ratio1200")

#reshape data
d1 <- dat %>%   
  dplyr::filter(grepl("^FPWW022", Plot_ID)) %>% 
  #get temp data
  extract_covars_from_nested("thermodata", "temp") %>%
  dplyr::select(Plot_ID, covariates, timestamp, temp) %>%
  mutate(timestamp = as.factor(timestamp)) %>%
  dplyr::rename(value = temp) %>%
  #get design
  unnest(covariates) %>%
  #rename covars
  rename(Fl0Ang = "angle",
         Heading_GDDAS = "heading_GDDAS",
         Ratio1200 = "SI_1200",
         OnSen_GDDAH = "Cnp_onsen_gom_GDDAH",
         Onsen_GDDAS = "Cnp_onsen_gom_GDDAS_fitted")

#Simple cor ----
d1.1 <- d1 %>% 
  #reshape
  dplyr::select(timestamp, value, covars) %>% 
  gather(., Covariate, covarval, FH:Ratio1200) %>% 
  group_by(timestamp, Covariate) %>% nest() %>% 
  #calculate correlations
  mutate(cor = purrr::map_dbl(data, do_cor_test, x = "covarval", y = "value", return = "estimate")) %>%
  mutate(p.val = purrr::map_dbl(data, do_cor_test, x = "covarval", y = "value", return = "p.value")) %>% 
  # mutate(cor = abs(cor)) %>% 
  dplyr::select(-data) %>% ungroup() %>% 
  mutate(trait = as.Date(timestamp))

p <- ggplot(d1.1) +
  # geom_point(aes(y = cor, x = trait), shape = 2, size = 3) + 
  # scale_shape_discrete(solid = F) +
  geom_line(aes(y = cor, x = trait, group = Covariate, color = Covariate, lty = Covariate), size = 1) + 
  xlab("Date") + ylab("Pearson correlation coefficient") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

png("Figures/2year/cors_CT_covars_atTPs.png", width = 6, height = 4, units = 'in', res = 300)
plot(p)
dev.off()

#==================================================================================== -

#MuLinRegr @TPs ----

moddata <- d1 %>% 
  dplyr::select(timestamp, value, one_of(covars)) %>% 
  group_by(timestamp) %>% group_nest() %>% 
  mutate(fit = purrr::map(data, lm, formula = as.formula("value~."))) %>% 
  #extract rsq and pval
  mutate(glance = purrr::map(fit,possibly(broom::glance, NA_real_)),
         rsq = purrr::map_dbl(glance,  possibly(`$`, NA_real_), "adj.r.squared"),
         pval = purrr::map_dbl(glance, possibly(`$`, NA_real_), "p.value"))

#==================================================================================== -