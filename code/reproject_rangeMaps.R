# This script reprojects mammal rangemaps into the Equi7Grid projection
# as per Simkin et al 2024 - Global exposure to zoonotic diseases in the wildland-urban interface

# It operates on a single mammal order in a single region defined by the Equi7Grid
# It may be edited to iterate over mammal orders and/or regions
# Each mammal range map dataset is first filtered to only include those that are recognized as hosts for zoonotic diseases


# Mammal range maps separated by order are availbale here: https://mol.org/datasets/ec694c34-bddd-4111-ba99-926a5f7866e8

# Host species names are stored in "SI_Table_3_host_pathogen_data.csv" 

# A new .gpkg file containing transformed range maps is produced for the mammal order and region




library(fs)
# library(stars)
library(sf)
sf_use_s2(FALSE)
library(tidyverse)

#Path to folder containing the Equi7Grid projection files for the region being assessed
#This sets the region for the analysis
#Note the last characters in the string must denote the region e.g. "AF/"
#Edit as appropriate
grid_path = "<insert folder path>" #e.g. "TUW-GEO-Equi7Grid-9ed362e/src/equi7grid/grids/AF/"

#get region as a character string
reg_chr = grid_path %>% str_sub(-3, -2) #either AF, AS, NA, SA, EU or OC


#Read in the land shapefiles
proj_land = dir_ls(grid_path, recurse = TRUE, regexp = "_PROJ_LAND.shp")[[1]] %>%
  st_read() %>%
  #Only keep the geometry column
  select(geometry)
# Get Equi7Grid projection for the region  
e7_crs = st_crs(proj_land)

#Get the land polygon in wgs84
geog_land = dir_ls(grid_path, recurse = TRUE, regexp = "_GEOG_LAND.shp")[[1]] %>%
  st_read() %>%
  #Only keep the geometry column
  select(geometry)

# Get a list of host species names
host_data_path = "<insert file path here" # e.g. "SI_Table_3_host_pathogen_data.csv"
host_data = read_csv(host_data_path)
host_list = host_data$mdd_sciName %>% unique()


#Get mammal range map polygon
#This sets the mammal order to be assessed
#Range map file path
rm_path = "<insert file path>"
p = st_read(rm_path, quiet = TRUE) %>%
  mutate(sciname = str_replace(sciname, " ", "_")) %>%
  filter(sciname %in% host_list) # filter to the host species

#Get the species that occur in the region by clipping to the land polygons
p_int = st_intersection(p, geog_land)
if(nrow(p_int) < 1) {
  message("No species in this region")
}else{
  #reproject
  p_e7 = st_transform(p_int, e7_crs)
  # Save in the same place as the mammal range maps
  w_fn = paste0(rm_path %>% str_sub(end = -6), "_", reg_chr, "_e7.gpkg")
  st_write(p_e7, w_fn, driver = "GPKG", append = FALSE, quiet = TRUE)
}
  