#This script produces a summary of the number of people living in WUI areas across binned host richness levels
#  as per Simkin et al 2024 - Global exposure to zoonotic diseases in the wildland-urban interface

#It runs for an individual region as defined by the Equi7Grid projection
#The following are required:
# - The global GHSL population data reprojected and tiled into the equi7grid projection as per "reproject_pop.R"
# - a set of host count rasters produced by script "host_counts.R" for the area of interest
# - a version of the Equi7Grid with all overlaps removed.  See script "clipped_grid.R"
# - aggregated WUI raster files. See script "aggregate_wui.R"
# - data tables provided in this download that identify matches between the WUI file structure and ZONE & TILE identifiers in the Equi7Grid data. File name is in the format "xx_namesMatch.csv"


library(fs)
library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)



#Specify region name
reg_ch = "<insert name>" # One of "Africa", "Asia", "Northamerica", "Southamerica", "Europe", "Oceania"
reg_2ch = case_when(
  reg_ch == "Africa" ~ "AF",
  reg_ch == "Asia" ~ "AS",
  reg_ch == "Northamerica" ~ "NA",
  reg_ch == "Southamerica" ~ "SA",
  reg_ch == "Europe" ~ "EU",
  reg_ch == "Oceania" ~ "OC"
)

# read the clipped geog tile grid

#Path to clipped raster
c_path = "<insert path>" #e.g. "equi7/"
#Read the clipped raster - update filename if nessesary
c = st_read(paste0(c_path, "geog_grid_", reg_2ch, "_clipped_1k.gpkg")) %>%
  #join multipolygons that are split across the dateline
  group_by(ZONE, TILE) %>%
  summarise()#This joins any polygons split across the dateline


#Proj grid - just to get the projected crs, then removed
proj_grid_path = "<insert path>" #e.g. "equi7/"
proj_grid = st_read(paste0(proj_grid_path, "EQUI7_V14_", reg_2ch, "_PROJ_TILE_T1.shp"))
e7 =  st_crs(proj_grid)
rm(proj_grid)

#Transform the clipped grid into the Equi7Grid projection
c_e7 = st_transform(c, crs = e7)

#Get the csv that shows matches between the wui dataset folder names and the 100k grid TILE ids
#Specify the folder containing the file
match_path = "<insert_path>" #e.g "wui/"
wuiMatch_df = read_csv(paste0(match_path, str_to_lower(reg_2ch), "_nameMatch.csv"))


#Specify location of host count rasters
hc_path = "<insert path>" #e.g. "mammals/host_counts/" 
#Specify location of the WUI files - aggregated as per "aggregate_wui.R"
wui_file_path = "<insert_path>" # e.g. "wui/"
#Specify location of GHSL population raster files
pop_file_path = "<insert path>" # e.g. "ghsl/ghslPop_e7/"

#Loop over the tiles in the clipped tile grid
hPops_l = lapply(1:nrow(c_e7), function(i){
  #Filter to a single tile
  g = c_e7 %>% slice(i)
  t = g$TILE
  print(t)
  #get host_counts
  #If no hosts are present return a blank raster
  hc_fn = paste0(hc_path, "host_count_", reg_ch, "_", t, ".tif")
  if(file_exists(hc_fn)){
    hc = read_stars(hc_fn) %>%
      setNames("h")}
  else{hc = st_as_stars(st_bbox(g), nx = 333, ny = 333, pretty = FALSE, inside = FALSE, values = 0) %>% setNames("h")}
  #create a raster mask to make cropping to the tile faster
  gTemplate = hc %>% setNames("g") %>% mutate(g = NA_real_)
  g_rast = st_rasterize(g %>% select() %>% mutate(g = 1), template = gTemplate)
  #Function that rounds to the highest multiple of x
  mround <- function(x,base){
      base*ceiling(x/base)
    }
    #create a raster with binned host counts
    h_bins = seq(0, mround(max(hc$h, na.rm = TRUE), 5), by = 5)
    if(length(h_bins) > 1){
      h_labs = tail(h_bins, -1)
      hc_binned = hc %>%
        mutate(h = cut(h, breaks = h_bins, right = TRUE, include.lowest = FALSE, labels = h_labs)) %>%
        mutate(h = as.numeric(as.character(h)))
      hc_binned[hc == 0] <- 0}
    else{hc_binned = hc}
    
    #mask habitat map to WUI - set non wui areas to NA
    #Update wui raster file naming convention as necessary
    wui_name = wuiMatch_df %>% filter(TILE == t) %>% pull(wui_name)
    if(length(wui_name > 0)){
      wuiRast = read_stars(paste0(wui_file_path, reg_2ch, "/", wui_name, "/", t, "_WUI_300m.tif"), proxy = FALSE) %>%
        setNames("w") %>%
        mutate(w = case_when(
          w == 0 ~ NA_real_, #water
          w >= 1 & w <= 4 ~ 1, #wui
          w >= 5 ~ NA_real_, #non-wui land
          TRUE ~ NA_real_ #anything else
        ))
    }else{
      wuiRast = st_as_stars(st_bbox(g), nx = 333, ny = 333, pretty = FALSE, inside = FALSE, values = NA_real_)
    }
    hc_binned[is.na(wuiRast)] <- NA_real_
    
    #Mask to tile extent
    hc_binned[is.na(g_rast)] <- NA_real_
    
    #get population raster
    #mask population raster to just those pixels that have >= 1 host
    #First need to disaggregate the habitat raster to 100m so that it will perfectly overlap with the population raster
    popFile = paste0(pop_file_path, "ghslPop_", reg_ch, "_", t, "_e7.tif") #Update file naming convention if nessesary
    # Read the raster or provide a blank raster if none exists
    if(file_exists(popFile)){
      tryCatch({
        popRast = read_stars(popFile, proxy = FALSE) %>%
          setNames("p")
      }, error=function(e){message("POPULATION ERROR - returning blank population raster")
        popRast = st_as_stars(st_bbox(g), nx = 333, ny = 333, pretty = FALSE, inside = FALSE, values = 0)})
    }else{
      message("NO POPULATION FILE - returning blank raster")
      popRast = st_as_stars(st_bbox(g), nx = 333, ny = 333, pretty = FALSE, inside = FALSE, values = 0)
    }
    #Mask to tile extent
    popRast[is.na(g_rast)] <- NA_real_
    #Create a template for warping
    t100 = popRast %>% mutate(p = NA_real_)
    #Warp to match resolution of the population raster
    hc_binned = hc_binned %>% st_warp(dest = t100, method = "near")
    #Get a list of hte unique host bins present
    h_unq = hc_binned$h %>% as.vector() %>% unique()
    h_unq = h_unq[!is.na(h_unq)]
    #Sum the population for each host bin
    hPops = lapply(h_unq, function(l){
      temp_popRast = popRast
      temp_popRast[hc_binned != l | is.na(hc_binned)] <- 0
      popSum = sum(temp_popRast$p %>% as.vector(), na.rm = TRUE)
      tibble(ZONE = reg_ch, TILE = t, h_bin = l, pop = popSum)
    }) %>% bind_rows()
    return(hPops)
})
  
reg_hPops = bind_rows(hPops_l)

#specify location to write files
write_path = "<insert path>" #e.g. "mammals/binnedPopCounts/"
#Write to csv file
write_fn = paste0(write_path, "binnedPop_summary_", reg_ch, "_", s, "_", e, ".csv")
write_csv(reg_hPops, write_fn)


