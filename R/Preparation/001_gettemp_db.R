
.libPaths("T:R3UserLibs")
#install.packages("devtools")
# devtools::install_git('https://gitlab.ethz.ch/crop_phenotyping/RCroPyDB/')

library(RCroPyDB)
library(openssl)
library(jsonlite)
library(dplyr)

# A json file has been downloaded from this link: 
# https://crop-phenotyping-db.ethz.ch/api/covariate_value?include=covariate&filter[objects]=[{"name":"frequency","op":"eq","val":"1:00:00"},{"name":"measurement","op":"has","val":{"name":"campaign","op":"has","val":{"name":"position_id","op":"in","val":[1]}}},{"name":"covariate_id","op":"eq","val":1},{"name":"timestamp","op":"gt","val":"2018-10-17"},{"name":"timestamp","op":"lt","val":"2019-07-25"}]&fields[covariate_value]=value,frequency,timestamp,position_z,covariate&fields[covariate]=name,label,si_unit
# "NaN" has been manually replaced by "null" in the downloaded json file
data <- jsonlite::fromJSON("O:/Projects/KP0011/2/from_db_2019/weather")
attributes <- data$data$attributes

dd <- attributes %>% filter(position_z == 2.0)

dd <- dd[1:6720,]

write.csv(dd[3:4], "O:/Projects/KP0011/1/Data preparation/Data/raw/data_climate_2019_db.csv", row.names = F)

#####

# A json file has been downloaded from this link: 
# https://crop-phenotyping-db.ethz.ch/api/covariate_value?include=covariate&filter[objects]=[{"name":"frequency","op":"eq","val":"1:00:00"},{"name":"measurement","op":"has","val":{"name":"campaign","op":"has","val":{"name":"position_id","op":"in","val":[1]}}},{"name":"covariate_id","op":"eq","val":1},{"name":"timestamp","op":"gt","val":"2017-10-16"},{"name":"timestamp","op":"lt","val":"2018-07-14"}]&fields[covariate_value]=value,frequency,timestamp,position_z,covariate&fields[covariate]=name,label,si_unit
# "NaN" has been manually replaced by "null" in the downloaded json file
data <- jsonlite::fromJSON("O:/Projects/KP0011/2/from_db_2019/weather_2018.txt")
attributes <- data$data$attributes

dd <- attributes %>% filter(position_z == 2.0)

dd <- dd[25:6500,]

write.csv(dd[3:4], "O:/Projects/KP0011/1/Data preparation/Data/raw/data_climate_2018_db.csv", row.names = F)

#####

AirT <- jsonlite::fromJSON("O:/Projects/KP0011/2/Data/data_raw/covars_from_db/airT_17_19")
attributes <- AirT$data$attributes

dd <- attributes %>% filter(position_z == 2.0)
write.csv(dd[3:4], "O:/Projects/KP0011/2/Data/data_raw/airT_db_17_19.csv", row.names = F)

#####

Precip <- jsonlite::fromJSON("O:/Projects/KP0011/2/Data/data_raw/covars_from_db/Precip_17_19")
attributes <- Precip$data$attributes

write.csv(attributes[3:4], "O:/Projects/KP0011/2/Data/data_raw/Precip_db_17_19.csv", row.names = F)

#####