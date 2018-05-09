library(tidyverse)
library(readxl)
library(raster)
library(maptools)
library(rgdal)
library(sf)
library(rmapshaper)
library(proj4shortcut)
library(leaflet)
library(RColorBrewer)
library(doParallel)
library(foreach)
library(gridExtra)

source('R/funcs.R')
prj <- geo_wgs84

######
# setup

# setwd('C:/proj/manuscripts/landscape_mod/')

data(calipsa)
data(comid_prd)

# state base
calista <- map_data('state', region = 'california')

# psa data fortified
psadat <- calipsa %>%
  as('Spatial') %>%
  as.data.frame %>%
  rownames_to_column('id') %>%
  dplyr::select(id, PSA6) %>%
  rename(PSA = PSA6) %>%
  mutate(Region = factor(PSA,
                         levels = c('Central Valley', 'Chaparral', 'Deserts Modoc', 'North Coast', 'Sierra Nevada', 'South Coast'),
                         labels = c('CV', 'Ch', 'DM', 'NC', 'SN', 'SC')
  ))
psafrt <- calipsa %>%
  dplyr::select(-PSA_REGION, -AREA, -PERIMETER, -PSA8) %>%
  as('Spatial') %>%
  fortify %>%
  left_join(psadat, by = 'id')

# all cali hydro, really simplified
calinhd <- readOGR('S:/Spatial_Data/NHDPlus/All_Ca_NHDPlus.shp') %>% 
  spTransform(prj) %>% 
  st_as_sf %>% 
  st_simplify(dTolerance = 0.5, preserveTopology = T) %>% 
  filter(COMID %in% unique(comid_prd$COMID)) %>% 
  dplyr::select(COMID) %>%
  inner_join(comid_prd, by = 'COMID')


# comid predictions to join with calinhd
comid_prd <- comid_prd %>% 
  dplyr::select(COMID, matches('^full|^core|^COMID$'))

######
# full vs core model

# full 
clsfl <- getcls2(calinhd, thrsh = 0.79, tails = 0.05, modls = 'full')

caliclsfl <- calinhd %>% 
  left_join(clsfl, by = 'COMID')

caliclsidfl <- caliclsfl %>% 
  as('Spatial') %>% 
  .@data %>% 
  dplyr::select(COMID, strcls) %>% 
  rownames_to_column('id') 
caliclsplofl <- caliclsfl %>% 
  as('Spatial') %>% 
  fortify %>% 
  left_join(caliclsidfl, by = 'id') 

toplofl <- caliclsplofl %>% 
  na.omit

# core

clscr <- getcls2(calinhd, thrsh = 0.79, tails = 0.05, modls = 'core')

caliclscr <- calinhd %>% 
  left_join(clscr, by = 'COMID')

caliclsidcr <- caliclscr %>% 
  as('Spatial') %>% 
  .@data %>% 
  dplyr::select(COMID, strcls) %>% 
  rownames_to_column('id') 
caliclsplocr <- caliclscr %>% 
  as('Spatial') %>% 
  fortify %>% 
  left_join(caliclsidcr, by = 'id') 

toplocr <- caliclsplocr %>% 
  na.omit

## 
# plots

# full plot
pfl <- ggplot(toplofl) +
  geom_path(aes(x = long, y = lat, group = id, colour = strcls), size = 0.30) +
  geom_polygon(data = psafrt, aes(x = long, y = lat, group = group), fill = NA, size = 0.5, colour = scales::alpha('black', 0.9)) +
  coord_map() +
  scale_colour_manual(values = pal_exp(levels(toplofl$strcls))) +
  guides(colour = guide_legend('Reach classification', title.position = 'top', ncol = 2)) +
  theme_void(base_family = 'serif') +
  theme(
    legend.position = 'none',
    plot.margin = grid::unit(c(0, 0, 0, 0), units = 'pt')
  ) + 
  ggtitle("Full")


# core plot
pcr <- ggplot(toplocr) +
  geom_path(aes(x = long, y = lat, group = id, colour = strcls), size = 0.30) +
  geom_polygon(data = psafrt, aes(x = long, y = lat, group = group), fill = NA, size = 0.5, colour = scales::alpha('black', 0.9)) +
  coord_map() +
  scale_colour_manual(values = pal_exp(levels(toplocr$strcls))) +
  guides(colour = guide_legend('Reach classification', title.position = 'top', ncol = 2)) +
  theme_void(base_family = 'serif') +
  theme(
    legend.position = 'top',
    legend.direction = 'horizontal',
    plot.margin = grid::unit(c(0, 0, 0, 0), units = 'pt')
  ) + 
  ggtitle("Core")
pleg <- g_legend(pcr)
pcr <- pcr + theme(legend.position = 'none')

setEPS(height = 6, width = 8, family = 'serif')
postscript('C:/Users/Marcus.SCCWRP2K/Desktop/p1.eps')
# png('C:/Users/Marcus.SCCWRP2K/Desktop/p1.png', height = 6, width = 8, units = 'in', res = 400, family = 'serif')
grid.arrange(
  pleg,
  arrangeGrob(pfl, pcr, ncol = 2), 
  ncol = 1, heights = c(0.2, 1)
)
dev.off()
######
# different tails

# 5th/95th
cls1 <- getcls2(calinhd, thrsh = 0.79, tails = 0.05, modls = 'full')

calicls1 <- calinhd %>% 
  left_join(cls1, by = 'COMID')

caliclsid1 <- calicls1 %>% 
  as('Spatial') %>% 
  .@data %>% 
  dplyr::select(COMID, strcls) %>% 
  rownames_to_column('id') 
caliclsplo1 <- calicls1 %>% 
  as('Spatial') %>% 
  fortify %>% 
  left_join(caliclsid1, by = 'id') 

toplo1 <- caliclsplo1 %>% 
  na.omit

# 25th/75th
cls2 <- getcls2(calinhd, thrsh = 0.79, tails = 0.25, modls = 'full')

calicls2 <- calinhd %>% 
  left_join(cls2, by = 'COMID')

caliclsid2 <- calicls2 %>% 
  as('Spatial') %>% 
  .@data %>% 
  dplyr::select(COMID, strcls) %>% 
  rownames_to_column('id') 
caliclsplo2 <- calicls2 %>% 
  as('Spatial') %>% 
  fortify %>% 
  left_join(caliclsid2, by = 'id') 

toplo2 <- caliclsplo2 %>% 
  na.omit

# 45th/55th
cls3 <- getcls2(calinhd, thrsh = 0.79, tails = 0.45, modls = 'full')

calicls3 <- calinhd %>% 
  left_join(cls3, by = 'COMID')

caliclsid3 <- calicls3 %>% 
  as('Spatial') %>% 
  .@data %>% 
  dplyr::select(COMID, strcls) %>% 
  rownames_to_column('id') 
caliclsplo3 <- calicls3 %>% 
  as('Spatial') %>% 
  fortify %>% 
  left_join(caliclsid3, by = 'id') 

toplo3 <- caliclsplo3 %>% 
  na.omit

##
# plots

# 5th/95th
p1 <- ggplot(toplo1) +
  geom_path(aes(x = long, y = lat, group = id, colour = strcls), size = 0.30) +
  geom_polygon(data = psafrt, aes(x = long, y = lat, group = group), fill = NA, size = 0.5, colour = scales::alpha('black', 0.9)) +
  coord_map() +
  scale_colour_manual(values = pal_exp(levels(toplo1$strcls))) +
  guides(colour = guide_legend('Reach classification', title.position = 'top', ncol = 2)) +
  theme_void(base_family = 'serif') +
  theme(
    legend.position = 'none',
    plot.margin = grid::unit(c(0, 0, 0, 0), units = 'pt')
  ) + 
  ggtitle("5th/95th")

# 25th/75th
p2 <- ggplot(toplo2) +
  geom_path(aes(x = long, y = lat, group = id, colour = strcls), size = 0.30) +
  geom_polygon(data = psafrt, aes(x = long, y = lat, group = group), fill = NA, size = 0.5, colour = scales::alpha('black', 0.9)) +
  coord_map() +
  scale_colour_manual(values = pal_exp(levels(toplo2$strcls))) +
  guides(colour = guide_legend('Reach classification', title.position = 'top', ncol = 2)) +
  theme_void(base_family = 'serif') +
  theme(
    legend.position = 'none',
    plot.margin = grid::unit(c(0, 0, 0, 0), units = 'pt')
  ) + 
  ggtitle("25th/75th")

# 45th/55th
p3 <- ggplot(toplo3) +
  geom_path(aes(x = long, y = lat, group = id, colour = strcls), size = 0.30) +
  geom_polygon(data = psafrt, aes(x = long, y = lat, group = group), fill = NA, size = 0.5, colour = scales::alpha('black', 0.9)) +
  coord_map() +
  scale_colour_manual(values = pal_exp(levels(toplo3$strcls))) +
  guides(colour = guide_legend('Reach classification', title.position = 'top', ncol = 2)) +
  theme_void(base_family = 'serif') +
  theme(
    legend.position = 'top',
    legend.direction = 'horizontal',
    plot.margin = grid::unit(c(0, 0, 0, 0), units = 'pt')
  ) + 
  ggtitle("45th/55th")

pleg <- g_legend(p3)
p3 <- p3 + theme(legend.position = 'none')

setEPS(height = 6, width = 10, family = 'serif')
postscript('C:/Users/Marcus.SCCWRP2K/Desktop/p2.eps')
# png('C:/Users/Marcus.SCCWRP2K/Desktop/p2.png', height = 6, width = 10, units = 'in', res = 400, family = 'serif')
grid.arrange(
  pleg,
  arrangeGrob(p1, p2, p3, ncol = 3), 
  ncol = 1, heights = c(0.2, 1)
)
dev.off()

######
# different thresholds

# 0.63
cls1 <- getcls2(calinhd, thrsh = 0.63, tails = 0.05, modls = 'full')

calicls1 <- calinhd %>% 
  left_join(cls1, by = 'COMID')

caliclsid1 <- calicls1 %>% 
  as('Spatial') %>% 
  .@data %>% 
  dplyr::select(COMID, strcls) %>% 
  rownames_to_column('id') 
caliclsplo1 <- calicls1 %>% 
  as('Spatial') %>% 
  fortify %>% 
  left_join(caliclsid1, by = 'id') 

toplo1 <- caliclsplo1 %>% 
  na.omit

# 0.79
cls2 <- getcls2(calinhd, thrsh = 0.79, tails = 0.05, modls = 'full')

calicls2 <- calinhd %>% 
  left_join(cls2, by = 'COMID')

caliclsid2 <- calicls2 %>% 
  as('Spatial') %>% 
  .@data %>% 
  dplyr::select(COMID, strcls) %>% 
  rownames_to_column('id') 
caliclsplo2 <- calicls2 %>% 
  as('Spatial') %>% 
  fortify %>% 
  left_join(caliclsid2, by = 'id') 

toplo2 <- caliclsplo2 %>% 
  na.omit

# 0.92
cls3 <- getcls2(calinhd, thrsh = 0.92, tails = 0.05, modls = 'full')

calicls3 <- calinhd %>% 
  left_join(cls3, by = 'COMID')

caliclsid3 <- calicls3 %>% 
  as('Spatial') %>% 
  .@data %>% 
  dplyr::select(COMID, strcls) %>% 
  rownames_to_column('id') 
caliclsplo3 <- calicls3 %>% 
  as('Spatial') %>% 
  fortify %>% 
  left_join(caliclsid3, by = 'id') 

toplo3 <- caliclsplo3 %>% 
  na.omit

##
# plots

# 0.63
p1 <- ggplot(toplo1) +
  geom_path(aes(x = long, y = lat, group = id, colour = strcls), size = 0.30) +
  geom_polygon(data = psafrt, aes(x = long, y = lat, group = group), fill = NA, size = 0.5, colour = scales::alpha('black', 0.9)) +
  coord_map() +
  scale_colour_manual(values = pal_exp(levels(toplo1$strcls))) +
  guides(colour = guide_legend('Reach classification', title.position = 'top', ncol = 2)) +
  theme_void(base_family = 'serif') +
  theme(
    legend.position = 'none',
    plot.margin = grid::unit(c(0, 0, 0, 0), units = 'pt')
  ) + 
  ggtitle("0.63")

# 0.79
p2 <- ggplot(toplo2) +
  geom_path(aes(x = long, y = lat, group = id, colour = strcls), size = 0.30) +
  geom_polygon(data = psafrt, aes(x = long, y = lat, group = group), fill = NA, size = 0.5, colour = scales::alpha('black', 0.9)) +
  coord_map() +
  scale_colour_manual(values = pal_exp(levels(toplo2$strcls))) +
  guides(colour = guide_legend('Reach classification', title.position = 'top', ncol = 2)) +
  theme_void(base_family = 'serif') +
  theme(
    legend.position = 'none',
    plot.margin = grid::unit(c(0, 0, 0, 0), units = 'pt')
  ) + 
  ggtitle("0.79")

# 0.92
p3 <- ggplot(toplo3) +
  geom_path(aes(x = long, y = lat, group = id, colour = strcls), size = 0.30) +
  geom_polygon(data = psafrt, aes(x = long, y = lat, group = group), fill = NA, size = 0.5, colour = scales::alpha('black', 0.9)) +
  coord_map() +
  scale_colour_manual(values = pal_exp(levels(toplo3$strcls))) +
  guides(colour = guide_legend('Reach classification', title.position = 'top', ncol = 2)) +
  theme_void(base_family = 'serif') +
  theme(
    legend.position = 'top',
    legend.direction = 'horizontal',
    plot.margin = grid::unit(c(0, 0, 0, 0), units = 'pt')
  ) + 
  ggtitle("0.92")

pleg <- g_legend(p3)
p3 <- p3 + theme(legend.position = 'none')

setEPS(height = 6, width = 10, family = 'serif')
postscript('C:/Users/Marcus.SCCWRP2K/Desktop/p3.eps')
# png('C:/Users/Marcus.SCCWRP2K/Desktop/p3.png', height = 6, width = 10, units = 'in', res = 400, family = 'serif')
grid.arrange(
  pleg,
  arrangeGrob(p1, p2, p3, ncol = 3), 
  ncol = 1, heights = c(0.2, 1)
)
dev.off()

