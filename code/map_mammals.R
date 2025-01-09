# This script creates habitat maps for each mammal species 
#  as per Simkin et al 2024 - Host richness in the global wildland-urban interface
# This script is designed to be run for a single global region as defined by the Equi7Grid projection (i.e. Africa, Asia, Northamerica, Southamerica, Europe or Oceania)


# Host species names are stored in the supplementary information file "host_pathogen_data.csv" 

# Habitat preferences for host species are available from Map of Life 
# Leaving habitat and elevation preferences blank will assign all "natural" land cover categories to a species and will not apply elevation limits.



#Note that this script creates individual species habitat maps clipped to each tile
#This means many thousands of files if this script is run for the whole dataset.
#The script will create a set of folders:
# - a folder is created for each species
# - subfolders contain habitat rasters for each global region as defined by the Equi7Grid projection
# We recommend running this on a subset of the full dataset and scaling up as necessary
# Ensure that you have adequate storage space


# This script requires that you have the required raster datasets tiled and 
#   reprojected into the Equi7Grid projection at the required (~300m) resolution stored locally on your machine
# The following are required:
# - CCI Landcover data: see script "reproject_cci.R"
# - Digital elevation model: see script "reproject_dem.R")
# - WUI data aggregated to ~300m resolution: see script "aggregate_wui.R"
# - Mammal rangemaps reprojected into the Equi7Grid projection: see script "reproject_rangeMaps.R"


# Note that the script below must be updated with file paths to local input and output locations




#Use install.packages("<packagename>") if the below are not installed
library(sf)
sf_use_s2(FALSE)
library(stars)
library(fs)
library(tidyverse)

#Read in the region 100k tile grid for the region being assessed

proj_100k_path = "<insert file path here>" #e.g. "EQUI7_V14_AF_PROJ_TILE_T1.shp"
proj_100k = st_read(proj_100k_path)

#get region as a character string
reg_chr = proj_100k %>% slice(1) %>% pull(ZONE)


#A list of "natural" cci landcover classes - assign these to any species that do not have prefs
#Mosaic cropland (> 50%) is excluded
#bare ground and permanent snow/ice is excluded
cci_nat = c(40, 50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, 100, 110, 120, 121, 122, 130, 140, 150, 151, 152, 153, 160, 170, 180)

#Read in the host pathogen data from supplementary information data "SI_Table_3_host_pathogen_data.csv" - update path as necessary
host_pathogen_data_path = "<insert file path here>" # e.g. "host_pathogen_data.csv"
host_pathogen_data = read_csv(host_pathogen_data_path)

#get a list of host names
host_names = host_pathogen_data$mdd_sciName %>% unique() 

#Read in the host habitat preferences from Map of Life - update path as necessary
prefs_path = "<insert file path here>" # e.g. "habitat_prefs.csv"
prefs = read_csv(prefs_path) %>% distinct() %>%
  mutate(prefs_esa_landcover = str_split(prefs_esa_landcover, ", "))
prefs$prefs_esa_landcover <- lapply(prefs$prefs_esa_landcover, function(x){as.numeric(x)})

#Set path to the range map data set
rmaps_path = "<insert path>"

#Read in range maps and join the habitat preferences
#Note that rangemaps must be reprojected into the Equi7Grid format prior to this stage (see reproject_rangeMaps.R)
rmaps = st_read(rmaps_path) %>%
  select(sciname, order) %>%
  rename(mdd_sciName = sciname, mdd_order = order) %>%
  mutate(mdd_sciName = str_replace(mdd_sciName, " ", "_")) %>%
  mutate(mdd_order = str_to_upper(mdd_order)) %>%
  filter(mdd_sciName %in% host_names) %>%
  st_collection_extract("POLYGON") %>%
  left_join(prefs, by = c("mdd_sciName", "mdd_order"))

# Workflow - filter the rangemaps by habitat and elevation
# Loop over species
# For each species:
# - get a list of the TILE variable for each 100k grid that intersects the habitat polygon
# - loop over the TILE list
#  - read the cci and dem into memory
#  - apply habitat filter
#  - write as tif

for(i in 1:nrow(rmaps)){
  p_e7 = rmaps %>% slice(i)
  mdd_sp = p_e7 %>% pull(mdd_sciName)
  print(mdd_sp)
  #create a directory for each species
  dir_create(mdd_sp)
  #Create a subdirectory for each region
  folderName = paste0(mdd_sp, "/", reg_chr)
  dir_create(folderName)
  #get habitat classes as a list
  lc_list = p_e7$prefs_esa_landcover %>% 
    unlist()
  #Get character vector of 100km grids names that intersect the species polygon
  gr_l = st_filter(proj_100k, p_e7) %>%
    pull(TILE)
  #loop over the grids
      for(g in 1:length(gr_l)) {
        fn = paste0(folderName, "/", mdd_sp, "_", reg_chr, "_", gr_l[[g]], "_hab.tiff")
        #cci raster file name - edit path as necessary
        g_cci_fn = paste0("cci/cci_100k_e7/cci_", reg_chr, "_", gr_l[[g]], "_e7.tif")
        #cci raster
        g_cci = read_stars(g_cci_fn, proxy = FALSE) %>% setNames("lccs_class")
        #dem filename - update path as nessesary
        g_dem_fn = paste0("dem/dem_100k_e7/dem_", reg_chr, "_", gr_l[[g]], "_e7.tif")
        #dem raster
        g_dem = read_stars(g_dem_fn, proxy = FALSE)
        #reclassify the cci raster so that any pixels that match a habitat preference become 1 and all other pixels set to 0
        c = g_cci %>% st_crop(p_e7, crop = FALSE) %>%
          mutate(lccs_class = case_when(
            lccs_class %in% lc_list ~ 1,
            TRUE ~ 0
          ))
        # #mask cells out of elevation
        e_min = p_e7 %>% pull(elev_min) %>% as.numeric()
        e_max = p_e7 %>% pull(elev_max) %>% as.numeric()
        d = g_dem %>% st_crop(p_e7) %>%
          setNames("x") %>%
          mutate(x = case_when(
            x < e_min ~ 0,
            x > e_max ~ 0,
            TRUE ~ 1
          ))
          #The final raster - cci cropped to elevation
          c[d == 0] <- 0
          
          #Write the clipped species raster - note that this is where the lazy operations are carried out.
          write_stars(c, fn, overwrite = TRUE)
      }
}

