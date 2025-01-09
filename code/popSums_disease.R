#Get a csv file that shows the number of people exposed to >=1 host for a given disease for each tile in a given region
#  as per Simkin et al 2024 - Global exposure to zoonotic diseases in the wildland-urban interface

# Requires the following:
# - The global GHSL population data reprojected and tiled into the equi7grid projection as per "reproject_pop.R"
# - Tile grid in WGS84 with overlapping tiles removed - See "clipped_grid.R"
# - Projected Equi7Grid tile grid available: https://zenodo.org/records/8190018
# - File with host pathogen associations "host_pathogen_data.csv" available with this download
# - File that shows matches between the WUI file structure and the Equi7Grid ZONE and TILE identifiers available with this download
# - Individual species habitat rasters as per "map_mammals.R"


library(fs)
library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)


#Set region name
reg_ch = "<insert region>" #One of "Northamerica", "Southamerica", "Asia", "Africa", "Europe", "Oceania"
reg_2ch = case_when(
  reg_ch == "Africa" ~ "AF",
  reg_ch == "Asia" ~ "AS",
  reg_ch == "Northamerica" ~ "NA",
  reg_ch == "Southamerica" ~ "SA",
  reg_ch == "Europe" ~ "EU",
  reg_ch == "Oceania" ~ "OC"
)


# read the clipped geog tile grid
c_path = "<insert path>" #e.g. "equi7/"
c = st_read(paste0(c_path, "geog_grid_", reg_2ch, "_clipped_1k.gpkg")) %>%
  #join multipolygons that are split across the dateline
  group_by(ZONE, TILE) %>%
  summarise()#join any polygons that are split across the dateline
  
  
  
#Proj grid - just to get the CRS then remove
proj_grid_path = "<insert path>" #e.g. "equi7/"
proj_grid = st_read(paste0("equi7/EQUI7_V14_", reg_2ch, "_PROJ_TILE_T1.shp"))
e7 =  st_crs(proj_grid)
rm(proj_grid)

#Transform the clipped grid into the Equi7Grid projection  
c_e7 = st_transform(c, crs = e7)
  

#Get host pathogen associations
#Update with local file path
dis_data_path = "<insert path>" #e.g. "mammals/"
dis_data = read_csv(paste0(dis_data_path, "host_pathogen_data.csv")) %>%
  select(mdd_sciName, gideon_disease) %>%
  distinct()
  
#Get the csv that shows matches between the wui dataset folder names and the 100k grid TILE ids
#Update file naming convention as necessary
match_path = "<insert path>" #e.g. "wui/"
wuiMatch_df = read_csv(paste0(match_path, reg_2ch, "/", str_to_lower(reg_2ch), "_nameMatch.csv"))
  

#Identify folder containing host habitat files
host_file_path = "<insert path>" #e.g. "mammals/hab_rasters/"
#Get full list of host habitat files for the region
habRast_fl_reg = list.files(host_file_path, recursive = TRUE, pattern = reg_ch, full.names = TRUE) %>%
  lapply(function(x){
    file.info(x) %>% select(size) %>%
      rownames_to_column(var = "path")
  }) %>% bind_rows() %>%
  filter(size > 0) %>%
  pull(path)
  
#Provide path to folder containing population rasters
pop_file_path = "<insert path>" #e.g. "pop/ghslPop_e7/"

#Provide path to folder containing aggregated wui rasters
wui_file_path = "<insert path>" #e.g. "wui/"

#Loop over tiles in the clipped, projected tile grid  
d_pop_df_l = lapply(1:nrow(c_e7), function(i){
  g = c_e7 %>% slice(i)
  t = g$TILE
  print(t)
    
  #get species list for the tile
  habRast_fl = habRast_fl_reg %>% str_subset(pattern = t)
  #get folder for WUI raster that corresponds to the tile  
  wui_name = wuiMatch_df %>% filter(TILE == t) %>% pull(wui_name)
  # wui_name = wuiMatch_df %>% filter(TILE == t) %>% pull(wui_name)
  #Get teh population raster - update file naming convention as nessesary 
  popFile = paste0(pop_file_path, "ghslPop_", reg_ch, "_", t, "_e7.tif")
    
    #First check that all nessesary files are present
    if(length(habRast_fl > 0) & length(wui_name > 0) & file_exists(popFile)){
      #make a dataframe with species name and disease - update file name convention as necessary
      spp_names = lapply(habRast_fl, function(f){
        spp = f %>% str_extract(pattern = "([^/]+$)") %>%
          str_split("_") %>%
          list_c() %>%
          .[1:2] %>%
          paste(collapse = '_')
        tibble(mdd_sciName = spp, fn = f)
      }) %>% bind_rows() %>%
        left_join(dis_data, by = "mdd_sciName")
      
      #group by disease and loop over disease
      dis = spp_names$gideon_disease %>% unique()
      d_pop_df = lapply(dis, function(d){
        print(d)
        #read and join hab rasters for all species
        d_fl = spp_names %>% filter(gideon_disease == d) %>%
          pull(fn)
        rList = lapply(d_fl, function(d_f) {
          tryCatch({
            read_stars(d_f) %>% setNames("h")
          }, error=function(e){
            message("Error reading habitat file")
            print(d_f)
            NA
          })
        })
        rList = rList[!is.na(rList)]
        
        #Combine rasters
        r = do.call("c", rList)
        r = st_redimension(r)
        #get sum of hosts at each pixel
        r = st_apply(r, MARGIN = 1:2, FUN = sum)
        r[r == 0] <- NA_real_
        
        #mask habitat map to WUI - update file name convention as appropriate
        wuiRast = read_stars(paste0(wui_file_path, reg_2ch, "/", wui_name, "/", t, "_WUI_300m.tif"), proxy = FALSE) %>%
          setNames("w") %>%
          mutate(w = case_when(
            w == 0 ~ NA_real_, #water
            w >= 1 & w <= 4 ~ 1, #wui
            w >= 5 ~ NA_real_, #non-wui land
            TRUE ~ NA_real_ #anything else
          ))
        # Set all non wui pixels to NA
        r[is.na(wuiRast)] <- NA_real_
        
        #crop to the tile boundary
        #Rasterize the tile to make this much faster than cropping to a polygon
        gTemplate = r %>% setNames("g") %>% mutate(g = NA_real_)
        g_rast = st_rasterize(g %>% select() %>% mutate(g = 1), template = gTemplate)
        r[is.na(g_rast)] <- NA_real_
        
        
        #mask population raster to just those pixels that have >= 1 host and sum population
        #First need to disaggregate the habitat raster to 100m so that it will perfectly overlap with the population raster
        tryCatch({
          popRast = read_stars(popFile, proxy = FALSE) %>%
            setNames("p")
          t100 = popRast %>% mutate(p = 0)
          r = r %>% st_warp(dest = t100, method = "near")
          popRast[is.na(r)] <- NA_real_
          #sum population
          d_pop = popRast %>% pull(p) %>% sum(na.rm = TRUE)
          tibble(ZONE = reg_ch, TILE = t, gideon_disease = d, pop = d_pop)
        }, error=function(e){message("POPULATION ERROR")
          tibble(ZONE = reg_ch, TILE = t, gideon_disease = NA, pop = 0)})
      })
    }else{
      message("missing data - skipping")
      tibble(ZONE = reg_ch, TILE = t, gideon_disease = NA, pop = 0)
    }
  }) 
  
  d_pop_df = bind_rows(d_pop_df_l)
  #Write file
  write_path = "<insert path>" #e.g. "pop_disSums/"
  write_fn = paste0(write_path, "pop_disSums_", reg_ch, ".csv")
  write_csv(d_pop_df, write_fn)

