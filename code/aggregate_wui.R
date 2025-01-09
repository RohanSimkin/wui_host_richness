# This script aggregates raster maps of the wildland urban interface into the 
#   resolution used in Simkin et al 2024 - Global exposure to zoonotic diseases in the wildland-urban interface

# WUI data can be downloaded from: https://zenodo.org/records/7941460 as described in:
#   Schug, F., Bar-Massada, A., Carlson, A.R. et al. The global wildland–urban interface. Nature 621, 94–99 (2023). https://doi.org/10.1038/s41586-023-06320-0

# WUI data is stored in a file structure that relates to the Equi7Grid
# We provide a file that matches the WUI file structure to the ZONE and TILE identifiers in the Equi7Grid (see "xx_nameMatch.csv")


#This script iterates over WUI tiles from a single region
#To reproject multiple regions, the script must be run separately for each region, or could be edited to iterate over regions

# Note that this script could run on a single computer - However it will be slow and produces a lot of data.
# We recommend that this code is run in a cluster environment where possible, so that the task of aggregating the full WUI dataset 
#   can be split across multiple nodes

#A copy of each wui raster file (.tif) is produced in a lower resolution



library(fs)
library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)


#The path to the folder containing WUI files - must only contain the WUI files from the region of interest
wui_dir = "<insert path here>" #e.g "wui/AF/"


#path to the match list - matching between the 100k grid and the WUI folder structure
matches_path = "<insert path here>" #e.g. "af_nameMatch.csv"
matches = read_csv(matches_path)

#Get a list off ALL the wui rasters for that region
#Note the $ selects only files that end in WUI.tif and excludes aux.xml files
files_l = dir_ls(wui_dir, recurse = TRUE, regexp = "WUI.tif$")

for(i in 1:length(files_l)) {
  f = files_l[[i]]
  print(f)
  #get the WUI folder name e.g. "X0050_Y0072"
  # This assumes that the folder structure of the WUI data has not been changed
  #i.e. the path looks like e.g. "wui/AF/X0050_Y0072/WUI.tiff"
  folder = f %>% str_split("/") %>% list_c() %>% .[[length(.) -1]]
  #Then match this to the "TILE" identifier in the Equi7Grid
  grid_no = matches %>% filter(wui_name == folder) %>% pull(TILE)
  
  #name to write aggregated raster - Note results are saved in the same folder structure as the WUI dataset
  fn = paste0(wui_dir, grid_no, "_WUI_300m.tif")
  #read wui raster
  r = read_stars(f, proxy = FALSE)
  #get wui raster bounding box
  r_bbox = r %>% st_bbox()
  # create a template for reprojecting
  t = r_bbox %>% st_as_stars(nx = 333, ny = 333, pretty = FALSE, inside = FALSE, values = 0)
  #Warp to the aggregated resolution  
  r_300 = st_warp(r, dest = t, method = "mode", use_gdal = TRUE, no_data_value = 0)
  #write result
  write_stars(r_300, dsn = fn, update = FALSE) 
}

