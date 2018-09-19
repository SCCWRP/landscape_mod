# get high res, unsimplified map of stream classes for all Cali NHD

library(tidyverse)
library(sf)
library(leaflet)
library(RColorBrewer)
library(rgdal)

# data(caliclsplo)
load(file = '../data/calipsa.RData')
load(file = '../data/comid_prd.RData')

source('funcs.R')

prj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# comid predictions to join with calinhd
comid_prd <- comid_prd %>% 
  dplyr::select(COMID, matches('^core|^COMID$'))

# all cali hydro, really simplified
calinhd <- readOGR('S:/Spatial_Data/NHDPlus/All_Ca_NHDPlus.shp') %>% 
  spTransform(prj) %>% 
  st_as_sf %>% 
  # st_simplify(dTolerance = 0.5, preserveTopology = T) %>% 
  filter(COMID %in% unique(comid_prd$COMID)) %>% 
  dplyr::select(COMID) %>%
  inner_join(comid_prd, by = 'COMID')

# get biological condition expectations
cls <- getcls2(calinhd, thrsh = 0.79, tails = 0.1, modls = 'core')

calicls <- calinhd %>% 
  left_join(cls, by = 'COMID')

caliclsid <- calicls %>% 
  as('Spatial') %>% 
  .@data %>% 
  dplyr::select(COMID, strcls) %>% 
  rownames_to_column('id') 
caliclsplo <- calicls %>% 
  as('Spatial') %>% 
  fortify %>% 
  left_join(caliclsid, by = 'id') 

# color palette for stream expectations
pal_exp <- colorFactor(
  palette = brewer.pal(9, 'Paired')[c(2, 1, 5, 6)],
  na.color = 'yellow',
  levels = c('likely unconstrained', 'possibly unconstrained', 'possibly constrained', 'likely constrained'))

# state base
calista <- map_data('state', region = 'california')

# psa data fortified
psadat <- calipsa %>%
  as('Spatial') %>%
  as.data.frame %>%
  rownames_to_column('id') %>%
  dplyr::select(id, PSA6) %>%
  dplyr::rename(PSA = PSA6) %>%
  mutate(Region = factor(PSA,
                         levels = c('Central Valley', 'Chaparral', 'Deserts Modoc', 'North Coast', 'Sierra Nevada', 'South Coast'),
                         labels = c('CV', 'Ch', 'DM', 'NC', 'SN', 'SC')
  ))
psafrt <- calipsa %>%
  dplyr::select(-PSA_REGION, -AREA, -PERIMETER, -PSA8) %>%
  as('Spatial') %>%
  fortify %>%
  left_join(psadat, by = 'id')

# all comid class
toplo <- caliclsplo %>% 
  na.omit

# comid class
p1 <- ggplot(toplo) +
  geom_path(aes(x = long, y = lat, group = id, colour = strcls), size = 0.30) +
  geom_polygon(data = psafrt, aes(x = long, y = lat, group = group), fill = NA, size = 0.5, colour = scales::alpha('black', 0.9)) +
  coord_map() +
  scale_colour_manual(values = pal_exp(levels(toplo$strcls))) +
  guides(colour = guide_legend('Segment classification', title.position = 'top', ncol = 2)) +
  theme_void(base_family = 'serif', base_size = 44) +
  theme(
    legend.position = 'top',
    legend.direction = 'horizontal',
    plot.margin = grid::unit(c(0, 0, 0, 0), units = 'pt')
  )

jpeg('C:/Users/Marcus.SCCWRP2K/Desktop/calires.jpeg', height = 36, width = 24, units = 'in', res = 600, family = 'serif')
p1
dev.off()
