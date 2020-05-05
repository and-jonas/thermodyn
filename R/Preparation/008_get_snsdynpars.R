#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Extract senescence dynamics parameters from scorings
# using linear interpolation and gompertz models.


# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#====================================================================================== -

rm(list = ls())
.libPaths("T:/R3UserLibs") #Specify path to R Libraries

library(tidyverse)
library(nls.multstart)
library(plotly)

path_to_data <- "O:/Projects/KP0011/2/Data/data_prep/" #Specify path to data
dirto <- "O:/Projects/KP0011/2/Data/data_prep/"
path_to_utils <- "O:/Projects/KP0011/2/thermodyn/R/Utils/" #Specify path to functions

source(paste0(path_to_utils, "001_spectra_utils.R"))

#====================================================================================== -

#FIT USING X = GDDAH ----

#load data
data <- read.csv(paste0(path_to_data, "scr_sca.csv"))

# extract DynPars: Scorings ----
data <- data %>% 
  # slice(1:150) %>% 
  tidyr::gather(Trait, Score, SnsFl0:SnsCnp, factor_key = TRUE) %>% 
  filter(complete.cases(.)) %>% 
  data.frame()

# do linear interpolation of scorings
# and fit nls using multstart
# Performs 250 times repeated NLS fitting (Levenberg-Marquardt algorithm)
# with random-search start parameter sets randomly sampled from a uniform
# distribution between upper and lower starting parameter bounds
data_fits <- data %>%
  # slice(1:150) %>% 
  group_by(Plot_ID, Trait) %>%
  nest() %>% 
  mutate(fit_lin = purrr::map(data, lin_approx, x = "grading_GDDAH")) %>% 
  mutate(fit_cgom = purrr::map(data, ~ nls_multstart(Score ~ Gompertz_constrained(b, M, tempsum = grading_GDDAH),
                                                     data = .x,
                                                     iter = 250,
                                                     start_lower = c(b = -0.1, M = 1600),
                                                     start_upper = c(b = 0, M = 2200),
                                                     convergence_count = 100)))

# new data frame of predictions
new_preds <- data %>%
  do(., data.frame(grading_GDDAH = seq(min(.$grading_GDDAH), max(.$grading_GDDAH), length.out = 1000), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(data, Plot_ID) %>%
  summarise(., min_gGDDAH = min(grading_GDDAH), max_gGDDAH = max(grading_GDDAH)) %>%
  ungroup()

# create new predictions
preds2 <- data_fits %>%
  mutate(x = purrr::map(fit_cgom, broom::augment, newdata = new_preds)) %>% unnest(x) %>% 
  dplyr::select(-data, -fit_cgom) %>% 
  merge(., max_min, by = "Plot_ID") %>%
  group_by(., Plot_ID) %>%
  filter(., grading_GDDAH > unique(min_gGDDAH) & grading_GDDAH < unique(max_gGDDAH)) %>%
  arrange(., Plot_ID, Trait, grading_GDDAH) %>%
  rename(., Score = .fitted) %>%
  ungroup()

# check whether model converged
convInfo <- data_fits %>%
  transmute(convInfo = purrr::map(purrr::map(fit_cgom, "convInfo"), "isConv"))

#extract parameters from nls fits
preds_gom <- preds2 %>% 
  dplyr::select(-fit_lin) %>% 
  group_by(Plot_ID, Trait) %>%
  group_nest() %>%
  mutate(pars_gom = purrr::map(data, extract_pars)) %>%
  select(-data)

#extract parameters from linear interpolation
preds_lin <- data_fits %>%
  unnest(fit_lin) %>% 
  mutate(data_lin = purrr::map(fit_lin, bind_cols)) %>%
  transmute(pars_lin = purrr::map(data_lin, extract_pars))

#combine predictions and convInfo
pred <- list(convInfo, preds_gom, preds_lin) %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("Plot_ID", "Trait")), .)

#combine predictions and convInfo
df <- pred %>% unnest(pars_gom, pars_lin, convInfo)
names(df)[4:length(df)] <- c("onsen_gom", "midsen_gom", "endsen_gom", "tsen_gom",
                             "onsen_lin", "midsen_lin", "endsen_lin", "tsen_lin")

#replace parameter value with NA if model did not converge
df[df$convInfo == FALSE, c(4:7)] <- NA

fwrite(df, file = paste0(dirto, "dynpars_scr_gddah.csv"))

#====================================================================================== -

#FIT USING X = GDDAS ----

rm(data)

#load data
data <- read.csv(paste0(path_to_data, "scr_sca.csv"))

# extract DynPars: Scorings ----
data <- data %>% 
  tidyr::gather(Trait, Score, SnsFl0:SnsCnp, factor_key = TRUE) %>% 
  data.frame()

data <- data[-which(is.na(data$Score)),]

# do linear interpolation of scorings
# and fit nls using multstart
# Performs 250 times repeated NLS fitting (Levenberg-Marquardt algorithm)
# with random-search start parameter sets randomly sampled from a uniform
# distribution between upper and lower starting parameter bounds
data_fits <- data %>%
  # slice(1:150) %>% 
  group_by(Plot_ID, Trait) %>%
  nest() %>% 
  mutate(fit_lin = purrr::map(data, lin_approx, x = "grading_GDDAS")) %>% 
  mutate(fit_cgom = purrr::map(data, ~ nls_multstart(Score ~ Gompertz_constrained(b, M, tempsum = grading_GDDAS),
                                                     data = .x,
                                                     iter = 250,
                                                     start_lower = c(b = -0.1, M = 1600),
                                                     start_upper = c(b = 0, M = 2200),
                                                     convergence_count = 100)))

# new data frame of predictions
new_preds <- data %>%
  do(., data.frame(grading_GDDAS = seq(min(.$grading_GDDAS), max(.$grading_GDDAS), length.out = 1000), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(data, Plot_ID) %>%
  summarise(., min_gGDDAS = min(grading_GDDAS), max_gGDDAS = max(grading_GDDAS)) %>%
  ungroup()

# create new predictions
preds2 <- data_fits %>%
  mutate(x = purrr::map(fit_cgom, broom::augment, newdata = new_preds)) %>% unnest(x) %>% 
  dplyr::select(-data, -fit_cgom) %>% 
  merge(., max_min, by = "Plot_ID") %>%
  group_by(., Plot_ID) %>%
  filter(., grading_GDDAS > unique(min_gGDDAS) & grading_GDDAS < unique(max_gGDDAS)) %>%
  arrange(., Plot_ID, Trait, grading_GDDAS) %>%
  rename(., Score = .fitted) %>%
  ungroup()

# check whether model converged
convInfo <- data_fits %>%
  transmute(convInfo = purrr::map(purrr::map(fit_cgom, "convInfo"), "isConv"))

#extract parameters from nls fits
preds_gom <- preds2 %>% 
  dplyr::select(-fit_lin) %>% 
  group_by(Plot_ID, Trait) %>%
  group_nest() %>%
  mutate(pars_gom = purrr::map(data, extract_pars)) %>%
  select(-data)

#extract parameters from linear interpolation
preds_lin <- data_fits %>%
  unnest(fit_lin) %>% 
  mutate(data_lin = purrr::map(fit_lin, bind_cols)) %>%
  transmute(pars_lin = purrr::map(data_lin, extract_pars))

#combine predictions and convInfo
pred <- list(convInfo, preds_gom, preds_lin) %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("Plot_ID", "Trait")), .)

#combine predictions and convInfo
df <- pred %>% unnest(pars_gom, pars_lin, convInfo)
names(df)[4:length(df)] <- c("onsen_gom_GDDAS", "midsen_gom_GDDAS", "endsen_gom_GDDAS", "tsen_gom_GDDAS",
                             "onsen_lin_GDDAS", "midsen_lin_GDDAS", "endsen_lin_GDDAS", "tsen_lin_GDDAS")

#replace parameter value with NA if model did not converge
df[df$convInfo == FALSE, c(4:7)] <- NA

fwrite(df, file = paste0(dirto, "dynpars_scr_gddas.csv"))

#====================================================================================== -

#MERGE ----

#convert to DATE and to GDDAS

df <- read_csv(paste0(dirto, "dynpars_scr_gddah.csv"))
sub <- data %>% dplyr::select(Plot_ID, heading_date, heading_DAS, heading_GDDAS) %>% unique()

all <- full_join(df, sub, by = "Plot_ID")

all <- all %>% 
  rename_at(vars(contains( "sen") ), list(~paste(., "GDDAH", sep = "_"))) %>% 
  mutate_at(vars(contains("sen")), funs(GDDAS = . + heading_GDDAS)) %>% 
  rename_at(vars(contains("GDDAH_GDDAS")), list(~gsub("GDDAH_", "", .))) %>% 
  dplyr::select(-contains("heading"))

SnsFl0 <- all %>% arrange(Plot_ID) %>% 
  filter(Trait == "SnsFl0")

SnsCnp <- all %>% arrange(Plot_ID) %>% 
  filter(Trait == "SnsCnp")

design1 <- read_csv(paste0(path_to_data, "exp_design.csv")) %>% filter(Exp %in% c("FPWW022", "FPWW024", "FPWW028")) %>% mutate(Trait = "SnsFl0")
design2 <- read_csv(paste0(path_to_data, "exp_design.csv")) %>% filter(Exp %in% c("FPWW022", "FPWW024", "FPWW028")) %>% mutate(Trait = "SnsCnp")

final <- left_join(design1, SnsFl0) %>% arrange(Plot_ID)
final2 <- left_join(design2, SnsCnp) %>% arrange(Plot_ID)

final <- bind_rows(final, final2) %>% dplyr::select(Plot_ID, Trait, contains("sen")) %>% 
  rename("level" = Trait) %>% 
  mutate(level = gsub("Sns", "", level))

dd <- read_csv(paste0(dirto, "dynpars_scr_gddas.csv")) %>% 
  ungroup() %>% 
  rename_at(vars(contains("sen")), list(~paste0(., "_fitted"))) %>% 
  rename("level" = Trait) %>% 
  mutate(level = gsub("Sns", "", level)) %>% 
  dplyr::select(Plot_ID, level, contains("sen"))

ff <- full_join(final, dd)

data.table::fwrite(ff, paste0(dirto, "dynpars_scr.csv"))

#====================================================================================== -

#CHECKS----

ff$dings <- ff$onsen_gom_GDDAS - ff$onsen_gom_GDDAS_fitted

dd <- ff %>% dplyr::select(onsen_gom_GDDAS, onsen_gom_GDDAS_fitted) %>% 
  mutate(diff = onsen_gom_GDDAS-onsen_gom_GDDAS_fitted) %>% 
  arrange(desc(abs(diff)))

p <- ggplot(ff) +
  geom_point(aes(x = onsen_gom_GDDAS,  y = onsen_gom_GDDAS_fitted))

ggplotly(p)

#====================================================================================== -

  