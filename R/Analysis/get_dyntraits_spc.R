
rm(list = ls())

# Helper functions ----

# couloured  diagonal panel
diag.col <- function(x, ...){
  ll <- par("usr") 
  rect(ll[1], ll[3], ll[2], ll[4], 
       col="grey")
}

#============================================================================================================== -

spc <- readRDS("Data/data_ready/A_spc.rds")
data <- readRDS("Data/006_data.rds")
design <- dplyr::select(data, Plot_ID, design) %>% unique()

# bm <- data %>% select(Plot_ID, covariates) %>% unnest(covariates) %>% unique() %>% 
#   dplyr::select(Plot_ID, SI_1200, FH, SnsCnp_onsen_gom_GDDAS) %>% 
#   rename("Ratio1200" = 2,
#          "OnSen" = 4)
# write.csv(bm, "Data/covariates.csv", row.names = F)

# extract ranks based on spatially corrected index value
scaled <- spc %>% 
  dplyr::select(Plot_ID, meas_date, Spc_corr) %>% unnest(Spc_corr) %>% 
  gather(., SI, value, starts_with("SI_")) %>% 
  group_by(meas_date, SI) %>% arrange(desc(value), .by_group = TRUE) %>% 
  mutate(value = f_scale(value))

sub <- scaled %>% 
  arrange(Plot_ID) %>% 
  slice(1:20)

ggplot(sub) +
  geom_point(aes(x = meas_GDDAH, y = value, color = SI)) +
  facet_wrap(~ Plot_ID)

#============================================================================================================== -

## Fit linear temporal trend ----
# to ranks
fit_scaled <- scaled %>% ungroup() %>% dplyr::select(Plot_ID, SI, value, meas_GDDAH) %>% 
  arrange(Plot_ID) %>% 
  group_by(Plot_ID, SI) %>% nest() %>% 
  mutate(fit = purrr::map(data, possibly(lm, NA_real_), formula = as.formula(value~meas_GDDAH))) %>% 
  mutate(slope = purrr::map_dbl(fit, possibly(get_lm_pars, NA_real_), pred = "meas_GDDAH")) %>% 
  mutate(glance = purrr::map(fit,possibly(broom::glance, NA_real_)),
         rsq = purrr::map_dbl(glance,  possibly(`$`, NA_real_), "r.squared"),
         pval = purrr::map_dbl(glance, possibly(`$`, NA_real_), "p.value"))

#PSRI results in highest Rsq
ord <- fit_scaled %>% ungroup() %>% group_by(SI) %>% 
  summarize(avg_rsq = mean(rsq)) %>% 
  arrange()

#============================================================================================================== -

#slopes correlation ----
dat_thm <- data %>%
  dplyr::select(Plot_ID, dyntraits_fc) %>% 
  unnest(dyntraits_fc) %>% unique() %>% 
  dplyr::select(Plot_ID, slope_scale)
  
dat_spc <- fit_scaled %>% 
  dplyr::select(Plot_ID, SI, slope)

ALL <- full_join(dat_thm, dat_spc)

wide <- ALL %>% 
  slice(1:1133) %>% 
  group_by(Plot_ID) %>% 
  spread(., SI, slope)

names(wide)[-1] <- c("slope_therm", "slope_MCARI2", "slope_NDVI", "slope_PSRI")

tiff("Figures/forpub/cors_slopes.tiff", width = 6, height = 5, units = 'in', res = 400)
pairs(wide[-1], upper.panel = panel.cor, lower.panel = panel.smooth, adj = 0, diag.panel = NULL)
dev.off()

#============================================================================================================== -

#Relative slopes ----
ddd <- dat_spc %>% 
  spread(., SI, slope) %>% 
  full_join(dat_thm,. , by = "Plot_ID") %>% 
  #calculate relative slopes
  mutate_at(vars(starts_with("SI_")), funs(.-slope_scale)) %>% 
  dplyr::select(-slope_scale) %>% 
  full_join(design, .) %>% 
  unnest(design) %>% 
  filter(check != "_NoCheck") %>% 
  #Remove outliers
  split(.$Gen_ID) %>%
  # lapply(., findOutlier, var = "SI_PSRI", cutoff = 3) %>%
  bind_rows() # No genotypic differences among the check varieties, NOISY!

dall <- dat_spc %>% 
  spread(., SI, slope) %>% 
  full_join(dat_thm,. , by = "Plot_ID") %>% 
  #calculate relative slopes
  mutate_at(vars(starts_with("SI_")), funs(.-slope_scale)) %>% 
  dplyr::select(-slope_scale) %>% 
  full_join(design, .) %>% 
  unnest(design) %>% 
  # #Remove outliers
  findOutlier(var = "SI_PSRI", cutoff = 3) %>%
  findOutlier(var = "SI_MCARI2", cutoff = 3) %>%
  findOutlier(var = "SI_NDVI_nb_ASD", cutoff = 3) %>%
  gather(., SI, value, starts_with("SI_")) %>% 
  mutate(SVI = strsplit(SI, "_") %>% lapply(., "[[", 2) %>% unlist() %>% as.factor()) %>% 
  droplevels()

meas_dates <- ddd %>% dplyr::select(Gen_Name, starts_with("SI")) %>% 
  gather(., SI, value, starts_with("SI_")) %>% 
  mutate(SVI = strsplit(SI, "_") %>% lapply(., "[[", 2) %>% unlist() %>% as.factor()) %>% 
  group_by(Gen_Name, SVI) %>% 
  summarise(x1 = min(value, na.rm = T),
            x2 = max(value, na.rm = T)) %>% 
  ungroup() %>% arrange(SVI) %>% 
  mutate(y1 = c(1.50, 1.47, 1.44,
                1.41, 1.38, 1.35,
                1.32, 1.29, 1.26)) %>% 
  mutate(y2 = c(1.53, 1.495, 1.465, 
                1.435, 1.405, 1.375,
                1.345, 1.315, 1.285)) %>% 
  mutate(Gen_Name = as.factor(Gen_Name))

tiff("Figures/forpub/distr_rel_slopes.tiff", width = 7, height = 6, units = 'in', res = 400)
ggplot(dall) +
  geom_density(aes(x = value, group = SVI, fill = SVI), alpha = 0.2) +
  geom_rect(meas_dates, mapping = aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = SVI), alpha=0.6) +
  geom_text(meas_dates, mapping = aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=Gen_Name), size=3, fontface = "bold") +
  labs(y = "Density", x = expression(Slope["SVI"] - Slope["CT"]))+
  scale_fill_manual(values = c("red", "seagreen", "yellow")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.9, 0.5))
dev.off()


outs1 <- dall %>% 
  arrange(value) %>% 
  slice(1:3) %>% pull(Gen_Name)

outs2 <- dall %>% unique() %>% 
  arrange(desc(value)) %>% 
  slice(1:3) %>% pull(Gen_Name)
  
##################################################################################### -

dat_thm <- readRDS("modeldata.rds") %>%
  #for checks only
  filter(!grepl("_NoCheck", check)) %>% 
  dplyr::select(Plot_ID, Gen_ID, slope_ranks) %>% 
  rename("slope" = slope_ranks)

dd <- dat_spc %>% 
  spread(., Index, slope) %>% 
  left_join(dat_thm, .) %>% 
  dplyr::select(-Plot_ID)

list <- split(dd, dd$Gen_ID)

lapply(list, "[", -1) %>% lapply(., 
                                 pairs, upper.panel = panel.cor)

##################################################################################### -

#Gen Eff
##Thermo trend
ggplot(dat_thm) +
  geom_boxplot(aes(x = Gen_ID, y = slope)) +
  geom_point(aes(x = Gen_ID, y = slope))

ggplot(dat_thm) +
  geom_histogram(aes(x = slope), bins = 8) + 
  facet_wrap(~Gen_ID, scales = "free")

mod <- lm(slope ~ Gen_ID, data = dat_thm)
anova(mod)

##Spectral trend
dd <- dat_spc %>% 
  filter(grepl("MCARI2|PSRI", Index)) %>% 
  left_join(dat_thm[1:2], ., by = "Plot_ID") %>% 
  dplyr::select(-Plot_ID)

ggplot(dd) +
  geom_boxplot(aes(x = Gen_ID, y = slope, group = interaction(Index, Gen_ID), color = Index))

ggplot(dd) +
  geom_histogram(aes(x = slope), bins = 12) + 
  facet_wrap(~Gen_ID*Index, scales = "free")

mod <- lm(slope ~ Gen_ID*Index, data = dd)
anova(mod)

##################################################################################### -
