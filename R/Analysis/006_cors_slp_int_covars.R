
rm(list = ls())
.libPaths("T:R3UserLibs")
detach("package:plyr")
library(tidyverse)
library(GGally)
library(gridExtra)
workdir <- "O:/Projects/KP0011/2/"
setwd(workdir)

#functions
source("thermodyn/R/Utils/002_thermo_utils.R")

# dat <- readRDS("Data/data_ready/6_data_dyntraits_corr.rds")
dat <- readRDS("Data/data_all_plots/6_data_dyntraits_corr.rds") ## use all available data


#============================================================================================================== -

# cor slope~rank ----
dyntraits <- grep("^dyntraits", names(dat), value = TRUE)
p <- list()
for (j in dyntraits){
  newdat <- dat %>% 
    extract_covars_from_nested(., from = "design", vars = c("harvest_year")) %>% 
    dplyr::select(Plot_ID, harvest_year, j) %>% unnest(j)
  cordat <- split(newdat, newdat$harvest_year)
  corplot <- list()
  for(i in names(cordat)){
    corplot[[i]] <- ggpairs(cordat[[i]][-c(1:2)],
                            lower = list(continuous = wrap("smooth", alpha = 0.3, size=1))) +
      ggtitle(paste(j, i, sep = "_")) +
      theme_bw() +
      theme(panel.grid = element_blank())
  }
  p[[j]] <- corplot
}

plist <- unlist(p, recursive = FALSE)

pdf("Figures/2year/pairs_int_slp.pdf", onefile = TRUE)
plist
dev.off()

#============================================================================================================== -

# cor slope~covars ----

#analyzed covariateS
covars <- c("FH", "heading_GDDAS", "Cnp_onsen_gom_GDDAH", "Cnp_onsen_gom_GDDAS_fitted",
            "angle", "Fl0Len", "Fl0Wid", "Fl0Glc", "SI_1200")

dyntraits <- grep("^dyntraits", names(dat), value = TRUE)

cors <- list()
for(i in dyntraits){
  print(i)
  #if slp derived from uncorrected temperature values, use raw covariate values
  if(i == "dyntraits"){
    d2 <- dat %>% 
      extract_covars_from_nested(from = "design", vars = "harvest_year") %>% 
      dplyr::select(Plot_ID, harvest_year, covariates, i) %>% 
      unnest(c(covariates, i)) %>% unique() %>% 
      dplyr::select(slp_ct, harvest_year, one_of(covars))
  # if slp derived from corrected temperature values or similar, use corrected covariate values (??)
  } else {
    d2 <- dat %>%
      extract_covars_from_nested(from = "design", vars = "harvest_year") %>%
      dplyr::select(Plot_ID, harvest_year, covariates_corr, i) %>%
      unnest(c(covariates_corr, i)) %>% unique() %>%
      dplyr::select(slp_ct, harvest_year, one_of(covars))
  }
  cors[[i]] <- d2 %>% 
    gather(., Covariate, covarval, FH:SI_1200) %>% 
    group_by(harvest_year, Covariate) %>% group_nest() %>% 
    #calculate correlations
    mutate(cor = purrr::map_dbl(data, possibly(do_cor_test, NA_real_), x = "covarval", y = "slp_ct", return = "estimate")) %>%
    mutate(p.val = purrr::map_dbl(data, possibly(do_cor_test, NA_real_), x = "covarval", y = "slp_ct", return = "p.value")) %>% 
    # mutate(cor = abs(cor)) %>% 
    dplyr::select(-data) %>% ungroup()
}

result <- cors %>% lapply(., dplyr::select, -p.val) %>%  
  Reduce(function(df1, df2) full_join(df1, df2, by = c("harvest_year", "Covariate")), .)

#============================================================================================================== -



