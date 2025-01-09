#This script reprojects and tiles the DEM data into the Equi7Grid projection.
#As per Simkin et al 2024 - Global exposure to zoonotic diseases in the wildland-urban interface



#The Equi7Grid can be downloaded from here: https://zenodo.org/records/8252376.

#DEM can be downloaded here: https://www.eorc.jaxa.jp/ALOS/en/dataset/aw3d30/aw3d30_e.htm

#The code below requires that you have a copy of the Equi7Grid projection files and DEM raster stored locally

#Note that this script produces a large amount of data - ensure you have adequate storage


#Use install.packages(<packagename>) if these packages are not installed
library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)


#Read the global DEM as a proxy object
#Update path to DEM file
dem_path = "<insert path here>"
#Read the DEM
dem = read_stars(dem_path, proxy = TRUE)

#The grid in wgs84 contains multipolygons that are split at the dateline.
#This function splits them into polygons so that the raster can be cropped for each part independently 
#then merged back together in the new crs

#A function for splitting multipolygons and merging the resulting raster
r_fun = function(x) {
  if(st_is(x, "MULTIPOLYGON")){
    split_g = x %>% st_cast("POLYGON")
    rList = c()
    for(i in 1:nrow(split_g)) {
      g_i = split_g %>% slice(i) %>% st_bbox() %>% st_as_sfc()
      r_g_i = dem[g_i] %>% st_as_stars()
      rList[[i]] <- r_g_i
    }
    if(length(rList) == 1){
      rList[[1]]
    }else{
      st_mosaic(rList[[1]], rList[[2]])  
    }
  }else{
    dem[x] %>% st_as_stars()
  }
}

#A global function for reading the dem raster that catches the error caused when 
#the dem raster does not overlap with the tile
err_fun = function(x) {
  tryCatch({
    r_fun(x)
  },
  error = function(e) {
    message("Error on read/crop")
    print(e)
    message("returning a grid of zeros")
    g %>% st_bbox() %>% st_as_stars(dx = 0.00277778, dy = -0.00277778, pretty = FALSE, inside = FALSE, values = 0)
  })
}


#Loop over global regions
regions = c("Africa", "Asia", "Northamerica", "Southamerica", "Oceania", "Europe")

lapply(regions, function(x){
  reg = x #e.g. "Africa"
  print(reg)
  reg_2ch = case_when(
    reg == "Africa" ~ "AF",
    reg == "Asia" ~ "AS",
    reg == "Northamerica" ~ "NA",
    reg == "Southamerica" ~ "SA",
    reg == "Europe" ~ "EU",
    reg == "Oceania" ~ "OC"
  )

  #Path to the Equi7grid files - update as appropriate
  grid_path = "TUW-GEO-Equi7Grid-9ed362e/src/equi7grid/grids/"
  
  # read the 100k grid in e7 and get the crs
  proj_100k = st_read(paste0(grid_path, reg_2ch, "/PROJ/EQUI7_V14_", reg_2ch, "_PROJ_TILE_T1.shp"))
  e7_crs =  proj_100k %>% st_crs()
  
  #read land polygons in wgs84
  geog_land = st_read(paste0(grid_path, reg_2ch, "/GEOG/EQUI7_V14_", reg_2ch, "_GEOG_LAND.shp"))
  
  #read the 100km tile polygons in wgs_84
  geog_100k = st_read(paste0(grid_path, reg_2ch, "/GEOG/EQUI7_V14_", reg_2ch, "_GEOG_TILE_T1.shp")) %>%
    st_filter(geog_land)
  
  #Loop over the 100k polygons and write to reprojected rasters
  for(i in 1:nrow(geog_100k)) {
    g = geog_100k %>% slice(i)
    z = g %>% pull(ZONE)
    t = g %>% pull(TILE)
    #filename to write projected dem to
    dem_fn = paste0("dem_", z, "_", t, "_e7.tif")
    #crop the dem proxy to the 100k grid - return grid of zeros and warn if there is an error
    r = err_fun(g)
    #get the corresponding projected grid cell
    g_proj = proj_100k %>% filter(TILE == t)
    #create a template for reprojecting
    proj_temp = g_proj %>% st_bbox() %>% st_as_stars(nx = 333, ny = 333, pretty = FALSE, inside = FALSE)
    #Warp
    w = st_warp(r, dest = proj_temp, method = "near", no_data_value = NA_real_)
    #write
    write_stars(w, dem_fn, overwrite = TRUE)
  }
})




