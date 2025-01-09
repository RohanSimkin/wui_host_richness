# This script produces a set of host count rasters at 300m scale in equi7grid projection
#   by summing the number of hosts with suitable habitat at the pixel level
#  as per Simkin et al 2024 - Global exposure to zoonotic diseases in the wildland-urban interface

# It requires that you have produced a set of tiled habitat rasters for each species of interest as 
#   detailed in the script "map_mammals.R"

#Equi7Grid files can be downloaded from here: https://zenodo.org/records/8252376.






library(fs)
library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)




#enter region name
reg_ch = "<insert region name>" # one of "Africa", "Asia", "Northamerica", "Southamerica", "Europe" or "Oceania"


# Update file paths with local file locations
# Get the projected land outline
proj_land_path = "<insert path>" # e.g. "equi7/EQUI7_V14_AF_PROJ_LAND.shp"
proj_land = st_read(proj_land_path)
#get the projected tile grid 
proj_grid_path = "<insert path>" # e.g. "equi7/EQUI7_V14_AF_PROJ_TILE_T1.shp"
proj_grid = st_read(proj_grid_path) %>%
  st_filter(proj_land)
e7_crs = st_crs(proj_grid)


#Specify the folder to write results into
write_path = "<insert path>" #e.g. "results/hostCount_rasters"

#Specify the folder where host habitat rasters are stored
hab_path = "<insert path>" #e.g. "results/hab_rasters"


#loop over the tiles
for(i in 1:nrow(proj_grid)){
  x = proj_grid %>% slice(i) %>% pull(TILE)
  print(x)
  #Update the destination for file to be written to local path
  write_fn = paste0(write_path,  "/host_count_", reg_ch, "_", x, ".tif")
  #get a list of all species habitat tiles in the region
  #Note the regexp command may need to be changed to match the file naming convention
  hab_fl = dir_info(hab_path, recurse = TRUE, regexp = paste0(reg_ch, "_", x), type = "file") %>%
      filter(size > 0) %>%
      pull(path)
    #read each raster and join as bands
    if(length(hab_fl) > 0){
      rList = lapply(hab_fl, function(f) {
        tryCatch({
          read_stars(f, proxy = FALSE)
        }, error=function(e){
          message("Error reading habitat file")
          print(f)
          NA
        })
      })
      r = do.call("c", rList)
      r = st_redimension(r)
      r_sum = st_apply(r, MARGIN = 1:2, FUN = sum, na.rm = TRUE)
      write_stars(r_sum, write_fn)  
    }else{
      message("No host species present - nothing written")
    }
}
