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

data(calipsa)
data(comid_prd)

# south coast poly
polysc <- calipsa %>% 
  dplyr::filter(PSA_REGION %in% 'South_Coast')

# state base
calista <- map_data('state', region = 'california')

# all cali hydro, really simplified
calinhd <- readOGR('S:/Spatial_Data/NHDPlus/All_Ca_NHDPlus.shp') %>% 
  spTransform(prj) %>% 
  st_as_sf %>% 
  filter(COMID %in% unique(comid_prd$COMID)) %>% 
  dplyr::select(COMID) %>%
  inner_join(comid_prd, by = 'COMID')

# south coast clips
scnhd <- calinhd[polysc, ] %>% 
  dplyr::select(COMID, core0.50) %>% 
  mutate(bcgcat = base::cut(core0.50, 
                            breaks = c(-Inf, 0.325, 0.625, 0.825, 1.025, Inf), 
                            labels = c('six', 'five', 'four', 'three', 'one/two'), 
                            right = F), 
         bcgcat = factor(bcgcat, levels = rev(levels(bcgcat)))
         )

# color scale
cols <- RColorBrewer::brewer.pal(6, 'Paired')[c(2, 1, 4, 3, 6)]

p1 <- ggplot() +
  geom_sf(data = scnhd, aes(colour = bcgcat, fill = bcgcat), size = 0.3) +
  geom_sf(data = polysc, fill = NA, size = 0.5, colour = scales::alpha('black', 0.9)) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  guides(
    colour = guide_legend('BCG bin', title.position = 'top'),
    fill = guide_legend('BCG bin', title.position = 'top')
    ) +
  theme_void(base_family = 'serif') +
  theme(
    plot.margin = grid::unit(c(0, 0, 0, 0), units = 'pt'), 
    panel.grid.major = element_line(colour = NA)
  )


png('C:/Users/Marcus.SCCWRP2K/Desktop/bcg_sc.png', height = 4, width = 6, units = 'in', res = 600, family = 'serif')
p1
dev.off()
