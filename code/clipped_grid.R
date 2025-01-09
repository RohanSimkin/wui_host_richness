# This script creates a version of the Equi7Grid that contains no overlaps across regions
# This can be used to summarize data to ensure that areas are not counted more than once
# as per Simkin et al 2024 - Global exposure to zoonotic diseases in the wildland-urban interface


# The following must be saved locally:
# - Natural Earth country boundaries can be downloaded here: https://www.naturalearthdata.com/downloads/50m-cultural-vectors/
# - World Bank income groups can be downloaded here: https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups

# The script contains paths to these files that need to be updated
# Ouptput location must also be updated manually



library(sf)
sf_use_s2(FALSE)
library(fs)
library(tidyverse)

#Provide path to Natural Earth country data:
countries_path = "<insert path>" #e.g. "countries/ne_50m_admin_0_countries.shp"
countries = st_read(countries_path) %>%
  filter(NAME != "Antarctica")
#Correct the names of French territories so that they are not all grouped in the analysis
french_territories = countries %>% filter(ADMIN == "France") %>%
  st_cast("POLYGON") %>%
  mutate(area = st_area(geometry)) %>%
  arrange(desc(area)) %>%
  mutate(NAME_LONG = c("France", "French Guiana", "France", "Réunion", "Martinique", "Guadeloupe", "Guadeloupe", "Mayotte", "France", "Guadeloupe"),
         ADM0_A3 = c("FRA", NA, "FRA", NA, NA, NA, NA, NA, "FRA", NA)) %>%
  select(-area) %>%
  distinct(NAME_LONG, ADM0_A3, .keep_all = TRUE)
countries = bind_rows(countries %>% filter(ADMIN != "France"), french_territories)



#Income data - See link above to download

#provide location of the World Bank Income Group data:
gni_path = "<insert_path"> # e.g. "gni/"

#The CLASS.csv file contains income group categorization
gni = read_csv(paste0(gni_path, "CLASS.csv")) %>%
  select("Economy", "Code", "Income group")
names(gni) <- c("Country_Name", "Country_Code", "Income_Group")

#edit these country names so they can be matched with the Natural Earth Countries data
gni = gni %>% mutate(Country_Code = case_when(
  Country_Code == "SSD" ~ "SDS", #South Sudan
  Country_Code == "XKX" ~ "KOS", #Kosovo
  TRUE ~ Country_Code
))

#Join to countries data
gni = countries %>% tibble() %>% select(NAME_LONG, ADM0_A3) %>%
  left_join(gni %>% select(Country_Code, Income_Group), by = c("ADM0_A3" = "Country_Code")) %>%
  mutate(lmic = ifelse(Income_Group == "High income", "High", "Low-Middle"))

#Note: There are several countries in the Natural Earth Data that don't have a direct match in the World Bank Income Group data
#E.g. Taiwan, Somaliland Palestine, French Guiana, Western Sahara
  


#Get the data that demarcates the Equi7Grid zones

#insert path to data stored locally
#This should be a folder that contains all of the Equi7Grid ZONE files e.g. "EQUI7_V14_SA_GEOG_ZONE.shp"
geog_zone_path = "<insert path>" #e.g. "equi7/"


#Create new zone polygons with overlaps removed
zone_fl = dir_ls(geog_zone_path, recurse = TRUE, regexp = "GEOG_ZONE.shp") %>%
  #Remove Antarctica
  str_subset("_AN_", negate = TRUE)
geog_zones = lapply(zone_fl, function(x){
  st_read(x)
}) %>% bind_rows() %>%
  st_difference() #Remove overlapping portions

#Get the full set of Equi7Grid 100km tile grid files

#Update file path to folder that contains the Equi7Grid 100km tile grids in WGS84 e.g. "EQUI7_V14_SA_GEOG_TILE_T1.shp"
gg_file_path = "<insert path>" #e.g. "equi7/"

gg_fl = dir_ls(gg_file_path, regexp = "_GEOG_TILE_T1.shp")

#Identify folder to write results to
out_path = "<insert path>" #e.g. "equi7/"

#Loop over the tile datasets
#Note that this script relies on the Equi7Grid file names not being changed 
#i.e. each filename must resemble "EQUI7_V14_SA_GEOG_TILE_T1.shp"
lapply(gg_fl, function(x){
  print(x)
  reg_2ch = x %>% str_split("_") %>% list_c() %>% .[[length(.) - 3]]
  reg = case_when(
    reg_2ch == "AF" ~ "Africa",
    reg_2ch == "AS" ~ "Asia",
    reg_2ch == "NA" ~ "Northamerica",
    reg_2ch == "SA" ~ "Southamerica",
    reg_2ch == "EU" ~ "Europe",
    reg_2ch == "OC" ~"Oceania"
  )
  gg = st_read(x) #Read the tiled grid file
  geog_zone = geog_zones %>% filter(ZONE == reg) #Filter to the relevant zone polygon
  #Identify the tiles that are not covered by the zone file and remove them
  containers = st_covered_by(gg, geog_zone)
  gg$ol <- lengths(containers) > 0
  gg$ol <- gg$ol == FALSE
  #Clip to the zone boundaries
  clipped = st_intersection(gg, geog_zone)
  #filter to the country boundaries
  clipped = st_filter(clipped, countries)
  #Identify just those tiles that are clipped to the zone boundaries and calculate thier new area
  clipped_ol = clipped %>% filter(ol == TRUE) 
  clipped_ol$a = as.numeric(st_area(clipped_ol))
  #Assign the area of tiles that do NOT overlap the zone boundaries
  clipped_nol = clipped %>% filter(ol == FALSE)
  clipped_nol$a = 100000^2
  #Merge the tiles back into a single dataset
  clipped = bind_rows(clipped_ol, clipped_nol)
  #Add country names - based on country with largest overlap
  clipped = st_join(clipped, countries %>% select(NAME_LONG, REGION_UN, SUBREGION), largest = TRUE)
  #Add teh income group based on country with largest overlap
  clipped = left_join(clipped, gni %>% select(NAME_LONG, Income_Group, lmic), by = "NAME_LONG")
  #Write to file
  write_sf(clipped, paste0(out_path, "geog_grid_", reg_2ch, "_clipped_100k.gpkg"))
})

