
rm(list = ls())
.libPaths("T:/R3UserLibs")
library(tidyverse)
library(GGally)

workdir <- "O:/Projects/KP0011/2/"
setwd(workdir)

source("thermodyn/R/Utils/002_thermo_utils.R")

#==================================================================================== -

#load data for selected plots and measurement time points
dat <- readRDS("Data/data_ready/4_data_selected.rds")
# dat <- readRDS("Data/data_all_plots/4_data_selected.rds")
design <- dplyr::select(dat, Plot_ID, design) %>% unique()

y <- c("thermodata", "thermodata_corr")
dat_dyn <- thermodata_sc <- list()

for(i in y){
  
  #reshape data
  ddd <- dat %>%   
    #get temp data
    extract_covars_from_nested(from = i, "temp") %>%
    dplyr::select(Plot_ID, covariates, i, timestamp, temp) %>%
    arrange(timestamp) %>% group_by(timestamp)
  
  #Scale (and revert?) CT values
  scaled <- ddd %>% 
    mutate(temp = f_scale(temp)) %>% 
    # mutate(temp = 1-temp) %>% 
    ungroup()
  
  #check distributions
  plotdat <- scaled %>%
    dplyr::select(timestamp, temp) %>%
    mutate(timestamp = as.factor(strptime(timestamp, "%Y-%m-%d"))) %>%
    mutate(year = str_sub(timestamp, 1, 4))
  library(plyr)
  plotdat.cor <- plyr::ddply(.data = plotdat, .(timestamp),summarize, n=paste("n =", length(temp)))
  p <- ggplot(plotdat) +
    geom_histogram(aes(x = temp, fill = year)) +
    xlab("scaled CT") +
    facet_wrap(~timestamp) +
    geom_text(data=plotdat.cor, aes(x=0.8, y=40, label=n), colour="black") +
    theme_bw() +
    theme(panel.grid = element_blank())
  detach("package:plyr")

  # png(paste0("Figures/2year/", i, "_scaled_distr.png"), width = 7, height = 7, units = 'in', res = 300)
  # plot(p)
  # dev.off()
    
  #reshape data
  scaled_ <- scaled %>% 
    extract_covars_from_nested(from = i, vars = "meas_GDDAH") %>%
    dplyr::select(Plot_ID, covariates, timestamp, meas_GDDAH, temp) %>%
    extract_covars_from_nested("covariates", vars = c("heading_GDDAS","Cnp_onsen_gom_GDDAH", "Cnp_onsen_gom_GDDAS_fitted"))
  
  # assemble output
  suffix <- ifelse(i == "thermodata_corr", "_fc", "")
  thermodata_sc[[i]] <- scaled_ %>% dplyr::select(Plot_ID, timestamp, meas_GDDAH, temp) %>% 
    group_by(Plot_ID, timestamp) %>% group_nest(.key = paste0("thermodata_sc", suffix))
  
  # get intercepts
  ## fit linear trends
  fit_scaled <- scaled_ %>% 
    #reshape
    ungroup() %>% dplyr::select(Plot_ID, temp, Cnp_onsen_gom_GDDAH, meas_GDDAH) %>% 
    arrange(Plot_ID) %>% 
    group_by(Plot_ID) %>% nest() %>% 
    #fit linear trend
    mutate(fit = purrr::map(data, possibly(lm, NA_real_), formula = as.formula(temp~meas_GDDAH))) %>% 
    unnest(data) %>% dplyr::select(Plot_ID, Cnp_onsen_gom_GDDAH, fit) %>% unique()
  ##predict scaled CT from linear trends
  predicted_ct <- fit_scaled %>%
    #exclude plots without heading date (and thus, no Cnp_onsen_gom_GDDAH)
    filter(!is.na(Cnp_onsen_gom_GDDAH)) %>% 
    dplyr::mutate(sc_onsen = map2_dbl(fit, 3, ~predict(.x, newdata = list(meas_GDDAH = Cnp_onsen_gom_GDDAH)))) %>% 
    dplyr::mutate(sc_midstg =  map2_dbl(fit, 3, ~predict(.x, newdata = list(meas_GDDAH = Cnp_onsen_gom_GDDAH/2))))
  pred_sc_uncorr <- predicted_ct %>% dplyr::select(Plot_ID, starts_with("sc_"))
  
  # get slopes
  fit_scaled <- scaled_ %>% 
    #reshape
    ungroup() %>% dplyr::select(Plot_ID, temp, meas_GDDAH) %>% 
    arrange(Plot_ID) %>% 
    group_by(Plot_ID) %>% group_nest() %>% 
    #fit linear trend
    mutate(fit = purrr::map(data, possibly(lm, NA_real_), formula = as.formula(temp~meas_GDDAH))) %>% 
    #extract slope
    mutate(slope = purrr::map_dbl(fit, possibly(get_lm_pars, NA_real_), pred = "meas_GDDAH")) %>% 
    #extract rsq and pval and rmse
    mutate(glance = purrr::map(fit,possibly(broom::glance, NA_real_)),
           rsq = purrr::map_dbl(glance,  possibly(`$`, NA_real_), "r.squared"),
           pval = purrr::map_dbl(glance, possibly(`$`, NA_real_), "p.value"),
           RMSE = purrr::map_dbl(fit, possibly(get_RMSE, NA_real_)))
  slopes <- fit_scaled %>% dplyr::select(Plot_ID, slope) %>% dplyr::rename("slp_ct" = slope)
           
  # assemble output
  dat_dyn[[i]] <- full_join(slopes, pred_sc_uncorr) %>% 
    group_by(Plot_ID) %>% group_nest(.key = paste0("dyntraits", suffix))
  
  # Model evaluation summary
  plotdata <- fit_scaled %>% 
    left_join(., design, by = "Plot_ID") %>% extract_covars_from_nested("design", "harvest_year") %>% 
    mutate(harvest_year = as.factor(harvest_year)) %>% 
    dplyr::select(Plot_ID, harvest_year, slope, rsq, pval, RMSE) %>% 
    gather(., par, value, slope:RMSE) %>% 
    arrange(match(par, c("slope", "rsq", "pval", "rmse")))
  ### create separate name vectors
  am_names <- c(
    pval = "p-value",
    rsq = "italic(R^{2})",
    slope = "Slope[NRCT%~%GDD]",
    RMSE = "italic(RMSE)"
  )
  p <- ggplot(plotdata) +
    geom_histogram(aes(x = value, fill = harvest_year), alpha = 0.5) +
    facet_wrap(~par, scales = "free_x",
               labeller = labeller(par  = as_labeller(am_names,  label_parsed))) +
    labs(y = "Count", x = "Value") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white")) +
    ggtitle("A")
  
  # png(paste0("Figures/2year/mod_eval", suffix,".png"), height = 4, width = 6, units = 'in', res = 300)
  # plot(p)
  # dev.off()
  
}  

#assemble output
out_dat1 <- thermodata_sc %>% Reduce(function(df1, df2) left_join(df1, df2, by = c("Plot_ID", "timestamp")),.)
out_dat2 <- dat_dyn %>% Reduce(function(df1, df2) left_join(df1, df2, by = "Plot_ID"),.)
OUT <- list(dat, out_dat1, out_dat2) %>% Reduce(function(df1, df2) left_join(df1, df2),.)

saveRDS(OUT, "Data/data_ready/5_data_dyntraits.rds")
# saveRDS(OUT, "Data/data_all_plots/5_data_dyntraits.rds")

#==================================================================================== -

# ##correlation between predicted ranks at different timepoints
# cordat <- predicted_ct %>% ungroup() %>% 
#   mutate(Exp = str_sub(Plot_ID, 1, 7) %>% as.factor()) %>% 
#   dplyr::select(Plot_ID, Exp, starts_with("sc_"))
# 
# pdf("Figures/pheno_cor_pred_sc.pdf")
# pairs(cordat[3:4], upper.panel = panel.cor, col = cordat$Exp)
# dev.off()

#==================================================================================== -

#plot fits ---- 

for (i in y){
  
  suffix <- ifelse(i == "thermodata_corr", "_fc", "")
  
  plotdat <- OUT %>% 
    dplyr::select(Plot_ID, design, paste0("thermodata_sc", suffix)) %>% 
    extract_covars_from_nested("design", c("Gen_Name", "Lot")) %>% 
    unnest(c(paste0("thermodata_sc", suffix))) %>% 
    mutate(Plot_ID = as.factor(Plot_ID),
           Gen_Name = as.factor(Gen_Name))
  
  formula <- y ~ x 
  p <- ggplot(plotdat, aes(x = meas_GDDAH, y = temp, color = Lot)) +
    geom_point(alpha = 0.5) +
    facet_wrap_paginate(~Gen_Name, ncol = 3, nrow = 3) +
    geom_smooth(method = "lm", formula = formula, se = F) +
    stat_poly_eq(aes(label = paste(..rr.label..)), 
                 label.x.npc = "right", label.y.npc = c(0.1, 0.2),
                 formula = formula, parse = TRUE, size = 3) +
    labs(x = "Thermal time after heading (GDD)", y = "Canopy Temperature") +
    # scale_color_manual(values = c("darkgreen", "red")) + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white"))
  required_pages <- n_pages(p)
  p <- list()  
  for(i in 1:required_pages){
    p[[i]] <- ggplot(plotdat, aes(x = meas_GDDAH, y = temp, color = Lot)) +
      geom_point(alpha = 0.5) +
      facet_wrap_paginate(~Gen_Name, ncol = 3, nrow = 3, page = i) +
      geom_smooth(method = "lm", formula = formula, se = F) +
      stat_poly_eq(aes(label = paste(..rr.label..)), 
                   label.x.npc = "right", label.y.npc = c(0.1, 0.2),
                   formula = formula, parse = TRUE, size = 3) +
      labs(x = "Thermal time after heading (GDD)", y = "Canopy Temperature") +
      # scale_color_manual(values = c("darkgreen", "red")) + 
      theme_bw() +
      theme(panel.grid = element_blank(),
            strip.background = element_rect(fill = "white"))
  }  
  pdf(paste0("Figures/2year/linear_fits", suffix, ".pdf"))
  print(p)
  dev.off()

}

#==================================================================================== -
