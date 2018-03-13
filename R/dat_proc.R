library(tidyverse)
library(readxl)
library(raster)
library(maptools)
library(rgdal)
library(sf)
library(rmapshaper)
library(proj4shortcut)

prj <- geo_wgs84

##
# csci data, statewide
csci_raw <- read_csv('data/raw/csci_raw.txt')
save(csci_raw, file = 'data/csci_raw.RData', compress = 'xz')

##
# psa 6
calipsa <- readOGR('S:/Spatial_Data/RCMP_needs editting/Inputs/PSA6_090111/PSA6_2011.shp') %>% 
  spTransform(prj) %>% 
  st_as_sf
save(calipsa, file = 'data/calipsa.RData', compress = 'xz')

##
# psa labels by centroid, deserts modoc has two
psalab <- calipsa %>% 
  as('Spatial') %>% 
  rmapshaper::ms_explode() %>% 
  st_as_sf %>% 
  mutate(
    AREA = st_area(.), 
    AREA = gsub('\\sm\\^2$', '', AREA), 
    AREA = as.numeric(AREA)
  ) %>% 
  group_by(PSA6) %>% 
  top_n(2, AREA) %>% 
  ungroup %>% 
  filter(!(PSA6 == 'South Coast' & AREA < 1e10)) %>% 
  filter(!(PSA6 == 'Chaparral' & AREA < 1e10)) %>% 
  st_centroid

psalab <- psalab %>% 
  st_coordinates %>% 
  data.frame %>% 
  rename(
    long = X, 
    lat = Y
  ) %>% 
  mutate(
    Region = psalab$PSA6,
    Region = factor(Region, 
      levels = c('Central Valley', 'Chaparral', 'Deserts Modoc', 'North Coast', 'Sierra Nevada', 'South Coast'),
      labels = c('CV', 'C', 'DM', 'NC', 'SN', 'SC')
    )
  )

save(psalab, file = 'data/psalab.RData', compress = 'xz')

##
# land use data from gap
# https://gapanalysis.usgs.gov/gaplandcover/data/download/

# master raster
gapland <- raster('Z:/MarcusBeck/GIS/gaplf2011lc_v30_CA/gaplf2011lc_v30_ca.tif')

# # gap key
# # ag is 555 - 557, urban is 580 - 584
# gapkey <- read_delim(
#   'C:/Users/Marcus.SCCWRP2K/Desktop/gaplf2011lc_v30_CA/GAP_LANDFIRE_National_Terrestrial_Ecosystems_2011_Attributes.txt',
#   delim = '\t')

# reclassify matrix, right closed
rcmat <- c(-1, 554, 1, 554, 557, 2, 557, 579, 1, 579, 584, 3) %>%
  matrix(., ncol = 3, byrow = T)

# get cali shapefile in format for clip
data(calishp)
toprj <- crs(gapland) %>% 
  as.character
calishp <- calishp %>% 
  st_transform(toprj)

# reclassify, aggregate, vectorize, clip, simplify, dissolve by feature, transform
# takes a minute
ludat <- gapland %>% 
  reclassify(rcmat) %>% 
  aggregate(fact = 10, fun = modal) %>% 
  rasterToPolygons(dissolve = T) %>% 
  st_as_sf %>% 
  st_intersection(calishp) %>% 
  ms_simplify(as(., 'Spatial')) %>% 
  st_as_af %>% 
  st_buffer(0) %>% 
  group_by(layer) %>% 
  summarize %>% 
  st_cast %>% 
  st_transform(prj)

save(ludat, file = 'data/ludat.RData', compress = 'xz')