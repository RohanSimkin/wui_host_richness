# This script reprojects and tiles the ESA CCI Landcover data into the Equi7Grid projection.
# As per Simkin et al 2024 - Global exposure to zoonotic diseases in the wildland-urban interface



#The Equi7Grid can be downloaded from here: https://zenodo.org/records/8252376.

#CCI Landcover can be downloaded here: https://www.esa-landcover-cci.org/?q=node/164
#Follow the instructions to download and transform the 2020 Land cover layer from NetCDF to .tif file using gdal

#The code below requires that you have a copy of the Equi7Grid projection files and 2020 CCI landcover raster stored locally

#The output is a reprojected raster for each 100km tile in the Equi7Grid
#Note that this script produces a large amount of data - ensure you have adequate storage


#Use install.packages(<packagename>) if these packages are not installed
library(stars)
library(sf)
sf_use_s2(FALSE)
library(tidyverse)


#Loop over global regions
regions = c("Africa", "Asia", "Northamerica", "Southamerica", "Oceania", "Europe")

#Path to the directory containing Equi7grid files - update as appropriate
grid_path = "<insert directory path>" # e.g.  "TUW-GEO-Equi7Grid-9ed362e/src/equi7grid/grids/"


lapply(regions, function(x){
  reg = x #e.g. "Africa"
  reg_2ch = case_when(
    reg == "Africa" ~ "AF",
    reg == "Asia" ~ "AS",
    reg == "Northamerica" ~ "NA",
    reg == "Southamerica" ~ "SA",
    reg == "Europe" ~ "EU",
    reg == "Oceania" ~ "OC"
  )
  #Read the Equi7Grid
  proj_grid = st_read(paste0(grid_path, reg_2ch, "/PROJ/EQUI7_V14_", reg_2ch, "_PROJ_TILE_T1.shp"))
  #get the crs
  e7_crs = st_crs(proj_grid)
  #read land polygons in wgs84
  geog_land = st_read(paste0(grid_path, reg_2ch, "/GEOG/EQUI7_V14_", reg_2ch, "_GEOG_LAND.shp"))
  
  #read the 100km tile polygons in wgs_84 and filter to the land
  geog_grid = st_read(paste0(grid_path, reg_2ch, "/GEOG/EQUI7_V14_", reg_2ch, "_GEOG_TILE_T1.shp")) %>%
    st_filter(geog_land)
  
  # read CCI raster data 
  cci_path = "<insert path to the 2020 CCI raster layer here>"
  cci = read_stars(cci_path, proxy = TRUE)
  
  
  #The grid in wgs84 contains multipolygons that are split at the dateline.
  #This function splits them into polygons so that the raster can be cropped for each part independently 
  #then merged back together in the new crs
  
  #A function for splitting multipolygons, cropping the cci raster, and merging the resulting raster
  r_fun = function(x) {
    if(st_is(x, "MULTIPOLYGON")){
      split_g = x %>% st_cast("POLYGON")
      rList = c()
      for(i in 1:nrow(split_g)) {
        g_i = split_g %>% slice(i) %>% st_bbox() %>% st_as_sfc()
        r_g_i = cci[g_i] %>% st_as_stars()
        rList[[i]] <- r_g_i
      }
      if(length(rList) == 1){
        rList[[1]]
      }else{
        st_mosaic(rList[[1]], rList[[2]])  
      }
    }else{
      cci[x] %>% st_as_stars()
    }
  }
  
  #Loop over the tile polygons and write to reprojected rasters
  #Specify output path here
  out_path = "<insert output path>"
  for(i in 1:nrow(geog_grid)) {
    g = geog_grid %>% slice(i)
    z = g %>% pull(ZONE)
    t = g %>% pull(TILE)
    fn = paste0(out_path, "cci_", z, "_", t, "_e7.tif")
    print(fn)
    #crop the cci raster to the extent of the tile
    r = r_fun(g)
    #get the corresponding projected tile
    g_proj = proj_grid %>% filter(TILE == t)
    #A template for reprojecting
    proj_temp = g_proj %>% st_bbox() %>% st_as_stars(nx = 333, ny = 333, pretty = FALSE, inside = FALSE)
    #Reproject
    w = st_warp(st_normalize(r), proj_temp, method = "near", no_data_value = NA_real_)
    write_stars(w, dsn = fn, overwrite = TRUE)
  }
})

