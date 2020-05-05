

rm(list = ls())
library(tidyverse)

source("O:/Projects/KP0011/2/thermodyn/R/Utils/001_spectra_utils.R")

scr <- f_sen_read(dir = "O:/Projects/KP0011/2/Data/data_prep/",
                  file_names = c("Senescence_FPWW022_DAS.csv", "Senescence_FPWW024_FPWW028_DAS.csv")) %>%
  f_invert_sen() %>%
  f_scale_sen()

data.table::fwrite(scr, "O:/Projects/KP0011/2/Data/data_prep/scr_sca.csv", row.names = FALSE)