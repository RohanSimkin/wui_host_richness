#This script reprojects and tiles the GHSL population raster data into the Equi7Grid projection.
#As per Simkin et al 2024 - Global exposure to zoonotic diseases in the wildland-urban interface



#The Equi7Grid can be downloaded from here: https://zenodo.org/records/8252376.

#GHSL population data can be downloaded here: https://ghsl.jrc.ec.europa.eu/datasets.php

#The code below requires that you have a copy of the Equi7Grid projection files and population raster stored locally

#Note that this script produces a large amount of data - ensure you have adequate storage


#Use install.packages("<packagename>") if these packages are not installed
library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)

#Read population raster
pop_path = "<insert file path here>"
pop = read_stars(pop_path, proxy = TRUE)


#Path to the Equi7grid files - update as appropriate
grid_path = "TUW-GEO-Equi7Grid-9ed362e/src/equi7grid/grids/"


r_fun = function(x) {
  split_g = x %>% st_cast("POLYGON")
  rList = c()
  for(i in 1:nrow(split_g)) {
    g_i = split_g %>% slice(i) %>% st_bbox() %>% st_as_sfc()
    pop_g_i = pop[g_i]
    rList[[i]] <- pop_g_i
  }
  return(rList)
}



#Loop over global regions
regions = c("Africa", "Asia", "Northamerica", "Southamerica", "Oceania", "Europe")

lapply(regions, function(x){

  #get region name
  reg = x #e.g. "Africa"
  reg_2ch = case_when(
    reg == "Africa" ~ "AF",
    reg == "Asia" ~ "AS",
    reg == "Northamerica" ~ "NA",
    reg == "Southamerica" ~ "SA",
    reg == "Europe" ~ "EU",
    reg == "Oceania" ~ "OC"
  )

  # read the 100k grid in e7 and get the crs
  proj_100k = st_read(paste0(grid_path, reg_2ch, "/PROJ/EQUI7_V14_", reg_2ch, "_PROJ_TILE_T1.shp"))
  e7_crs =  proj_100k %>% st_crs()

  #read land polygons in wgs84
  geog_land = st_read(paste0(grid_path, reg_2ch, "/GEOG/EQUI7_V14_", reg_2ch, "_GEOG_LAND.shp"))
  
  #read the 100km tile polygons in wgs_84
  geog_100k = st_read(paste0("equi7/EQUI7_V14_", reg_2ch, "_GEOG_TILE_T1.shp")) %>%
    st_filter(geog_land)
  
  #Loop over the 100k polygons and write to reprojected rasters
  for(i in 1:nrow(geog_100k)) {
    g = geog_100k %>% slice(i)
    z = g %>% pull(ZONE)
    t = g %>% pull(TILE)
    #filename to write projected population raster to
    pop_fn = paste0("ghslPop_", z, "_", t, "_e7.tif")
    print(pop_fn)
    #crop the pop proxy to the 100k grid
    #r is a list of 1 unless the grid crosses the dateline, then it will be length 2
    r = r_fun(g)
    #get the corresponding projected grid cell
    g_proj = proj_100k %>% filter(TILE == t)
    message("warping-------------------")
    w_list = lapply(r, function(x){
      w = st_warp(x %>% st_as_stars(), crs = e7_crs)
      return(w[g_proj])})
    if(length(w_list) == 1){
      write_stars(w_list[[1]], dsn = paste0("ghsl/ghslPop_e7/", pop_fn), overwrite = TRUE)
      }else{
        w1 = w_list[[1]] %>% st_as_stars()
        w2 = w_list[[2]] %>% st_as_stars()
        w = st_mosaic(w1, w2)
        write_stars(w, pop_fn, overwrite = TRUE)
      }
    } 
})