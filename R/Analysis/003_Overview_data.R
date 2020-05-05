
#==================================================================================== -

rm(list = ls())
.libPaths("T:/R3UserLibs")
library(tidyverse)
library(data.table)
workdir <- "O:/Projects/KP0011/2/"
setwd(workdir)

source("thermodyn/R/Utils/002_thermo_utils.R")

#==================================================================================== -

#prepare data (full and subsets)
dat_all <- readRDS("Data/data_ready/3_data.rds")
dat_all <- dat_all %>% 
  mutate(subset = "None")
dat_sub <- readRDS("Data/data_ready/4_data_selected.rds") %>% 
  mutate(subset = "Sel")
data <- bind_rows(dat_all, dat_sub)
year_list <- data %>% extract_covars_from_nested(., from = "design", vars = c("harvest_year")) %>% 
  split(., .$harvest_year)
subsets <- list()
for(i in 1:length(year_list)){
  subsets[[i]] <- split(year_list[[i]], year_list[[i]]$subset)
}
list <- purrr::flatten(subsets)

#Create overview plots for different subsets and years ----
plot <- list()
for(i in 1:length(list)){
  
  #get data subset
  dd <- list[[i]] %>% 
    # filter(timestamp < "2018-06-26 11:48:37 UTC") %>% 
    arrange(Plot_ID) %>%
    extract_covars_from_nested(., from = "thermodata", c("meas_GDDAH", "meas_GDDAS", "temp")) %>%
    extract_covars_from_nested(., from = "covariates",
                               vars = c("Cnp_onsen_gom_GDDAS", "Fl0_onsen_gom_GDDAS",
                                        "heading_GDDAS"))
 
 #get measurement dates
 meas_dates <- dd %>% dplyr::select(timestamp, meas_GDDAS, starts_with("Cnp")) %>% 
    group_by(., timestamp) %>% 
    summarise(meas_GDDAS = mean(meas_GDDAS)) %>% 
    mutate(Date = as.factor(as.Date(timestamp, "%Y-%m-%d"))) %>% 
    rename("r" = timestamp)
  
  #get phenology ranges
  onsen_cnp <- dd %>% dplyr::select(Plot_ID, starts_with("Cnp")) %>% 
    summarise(x1 = min(Cnp_onsen_gom_GDDAS, na.rm = TRUE),
              x2 = max(Cnp_onsen_gom_GDDAS, na.rm = TRUE),
              y1 = 44,
              y2 = 44.5)
  onsen_fl0 <- dd %>% dplyr::select(Plot_ID, starts_with("Fl0")) %>% 
    summarise(x1 = min(Fl0_onsen_gom_GDDAS, na.rm = TRUE),
              x2 = max(Fl0_onsen_gom_GDDAS, na.rm = TRUE),
              y1 = 47.5,
              y2 = 48)
  heading <- dd %>% dplyr::select(Plot_ID, heading_GDDAS) %>% 
      summarise(x1 = min(heading_GDDAS, na.rm = TRUE),
                x2 = max(heading_GDDAS, na.rm = TRUE),
                y1 = 45,
                y2 = 46)
  
  #Create inset histograms to show the distribution in phenologcy
  inset_1 <- ggplot(dd) +
    geom_histogram(aes(x= heading_GDDAS), fill = "darkgreen", alpha = 0.5) +
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.margin = margin(0,0,0,0))
  inset_2 <- ggplot(dd) +
    geom_histogram(aes(x= Cnp_onsen_gom_GDDAS), fill = "tan4", alpha = 0.5) +
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.margin = margin(0,0,0,0))
  inset_3 <- ggplot(dd) +
    geom_histogram(aes(x= Fl0_onsen_gom_GDDAS), fill = "goldenrod4", alpha = 0.5) +
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.margin = margin(0,0,0,0))

  #create plot
  plot[[i]] <- ggplot() +
    geom_line(dd, mapping = aes(x = meas_GDDAS, y = temp, group = Plot_ID), alpha = 0.1) +
    geom_point(meas_dates, mapping = aes(x = meas_GDDAS, y = 25), shape = 17, size = 2, color="darkgrey", fill = "darkgrey") +
    geom_rect(onsen_cnp, mapping = aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha = 0.75, fill = "tan4") +
    geom_rect(onsen_fl0, mapping = aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha = 0.75, fill = "goldenrod4") +
    geom_rect(heading, mapping = aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha = 0.75, fill = "darkgreen") +
    xlim(1250, 2250) + ylim(25, 50) +
    annotation_custom(grob=ggplotGrob(inset_1), 
                      ymin = 46, ymax=49, xmin=min(dd$heading_GDDAS, na.rm = T), xmax=max(dd$heading_GDDAS, na.rm = T)) +
    annotation_custom(grob=ggplotGrob(inset_2), 
                      ymin = 44.5, ymax=47, xmin=min(dd$Cnp_onsen_gom_GDDAS, na.rm = T), xmax=max(dd$Cnp_onsen_gom_GDDAS, na.rm = T)) +
    annotation_custom(grob=ggplotGrob(inset_3), 
                      ymin = 48, ymax=50.5, xmin=min(dd$Fl0_onsen_gom_GDDAS, na.rm = T), xmax=max(dd$Fl0_onsen_gom_GDDAS, na.rm = T)) +
    xlab("Thermal time after sowing (GDD)") + ylab ("Canopy temperature (°C)") +
    # scale_y_continuous(limits = c(-4, 12.5)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank()) 
  
}

#arrange plots and save
png("Figures/2year/CT_trends.png", width = 7, height = 7, units = 'in', res = 300)
gridExtra::grid.arrange(arrangeGrob(label,
                                     grobs=plot, ncol=2))
dev.off()

#==================================================================================== -

#Create meteoplot ----

airT <- read.csv("Data/data_raw/airT_db_17_19.csv")
Precip <- read.csv("Data/data_raw/Precip_db_17_19.csv")
data_climate <- full_join(airT, Precip, by = "timestamp")
names(data_climate) <- c("timestamp", "mean_temp", "precipitation")

nrow(airT); nrow(Precip); nrow(data_climate) #interessant

# calcualte mean daily temperature, using hourly means 
data_climate$timestamp <- gsub("T", " ", data_climate$timestamp)
colnames(data_climate) <- c("time_string", "mean_temp","precipitation")
# data_climate$time_string <- strptime(data_climate$time_string, "%d.%m.%Y")
data_climate$time_string <- strptime(data_climate$time_string, "%Y-%m-%d")

# calcualte mean daily temperature, using hourly means 
day_temp <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), mean, na.rm = TRUE)
names(day_temp)[2] <- "mean"

# get min and max hourly temperatures and calculate meanMM temperature  
day_temp$min <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), min, na.rm = TRUE)$x
day_temp$max <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), max, na.rm = TRUE)$x
day_temp$meanMM <- (day_temp$max + day_temp$min)/2

# calcualte cumulative daily rainfall, using hourly values 
day_rainfall <- aggregate(data_climate$precipitation, list(day = cut(data_climate$time_string, "days")), sum, na.rm = TRUE)
names(day_rainfall)[2] <- "precipitation"

meteo_daily <- full_join(day_temp, day_rainfall) %>% as_tibble() %>% 
  mutate(day = as.Date(strptime(day, "%Y-%m-%d"))) %>% 
  mutate(year = str_sub(day, 1,4))

#get period of interest
subset <- meteo_daily %>% filter(day %between% c("2018-05-01", "2018-07-13") | day %between% c("2019-05-01", "2019-07-23"))

#get heading date and add to plot
phenodat <- dat_all %>% 
  extract_covars_from_nested(., from = "design", vars = c("harvest_year")) %>% 
  extract_covars_from_nested(., from = "covariates", vars = c("heading_date"))
heading <- phenodat %>% group_by(harvest_year) %>% 
  dplyr::select(Plot_ID, harvest_year, heading_date) %>% 
  summarise(x1 = min(heading_date, na.rm = TRUE),
            x2 = max(heading_date, na.rm = TRUE),
            y1 = 29,
            y2 = 30) %>% 
  dplyr::rename("year" = harvest_year) %>% 
  mutate(year = as.character(year))

#create figure
p <- ggplot(subset) +
  geom_line(aes(day, mean), col = "red") +
  facet_wrap(~year, scales = "free") +
  geom_bar(aes(day, precipitation/2), stat = "identity", col = "dodgerblue4", fill = "dodgerblue4") +  
  geom_rect(heading, mapping = aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha = 0.75, fill = "darkgreen") +
  scale_y_continuous(sec.axis = sec_axis(~.*2, name = "Rainfall [mm/m2]\n"), limits = c(0,30)) + # trick to move axis labels add \n to introduce a line break
  labs(x = "", y = "Temperature [°C]") +  #re-insert x = "\nDate" if showing "Date" is needed
  theme_bw() +
  theme(panel.grid = element_blank())

png("Figures/2year/meteo_heading.png", width = 8, height = 4, units = 'in', res = 300)
plot(p)
dev.off()

#==================================================================================== -
