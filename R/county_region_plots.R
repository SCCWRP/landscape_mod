library(tidyverse)
library(sf)
library(leaflet)
library(RColorBrewer)
library(rgdal)

source('R/funcs.R')

prj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

strcls_all <- st_read('L:/Channels in developed landscapes_RM/Data/strm_constraints/strm_constraints.shp') 
library(tidyverse)
library(sf)
library(leaflet)
library(RColorBrewer)
library(rgdal)

# color palette for stream expectations
pal_exp <- colorFactor(
  palette = brewer.pal(9, 'Paired')[c(2, 1, 5, 6)],
  na.color = 'yellow',
  levels = c('likely unconstrained', 'possibly unconstrained', 'possibly constrained', 'likely constrained'))

source('R/funcs.R')

prj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

strcls_all <- st_read('L:/Channels in developed landscapes_RM/Data/strm_constraints/strm_constraints.shp') 

# rb 3 boundary
rbnds <- st_read('S:/Spatial_Data/RWQCBdistricts/rwqcbnda.shp') %>% 
  st_transform(crs = prj) %>% 
  filter(RB == 3) %>% 
  st_simplify(dTolerance = 0.001, preserveTopology = T)

# SB county
sbco <- st_read('S:/Spatial_Data/CA_Counties/cnty24k97.shp') %>% 
  st_transform(crs = prj) %>% 
  filter(NAME %in% 'Santa Barbara') %>% 
  st_simplify(dTolerance = 0.001, preserveTopology = T)

# sb county, by class -----------------------------------------------------

strcls <- strcls_all[sbco, ] %>% 
  dplyr::select(COMID, Ref10) 
toplo <- strcls %>% 
  # st_simplify(dTolerance = 0.01, preserveTopology = T) %>% 
  rename(strcls = Ref10) %>% 
  mutate(
    strcls = factor(strcls, 
                    levels = c('likely unconstrained', 'possibly unconstrained', 'possibly constrained', 'likely constrained')
    )
  ) %>% 
  filter(!is.na(strcls)) %>% 
  .[sbco, ]

# comid class
p1 <- ggplot() + 
  geom_sf(data = sbco) + 
  geom_sf(data = toplo, aes(colour = strcls, fill = strcls)) +
  scale_colour_manual(values = pal_exp(levels(toplo$strcls))) +
  scale_fill_manual(values = pal_exp(levels(toplo$strcls))) +
  guides(
    colour = guide_legend('Segment classification', title.position = 'top', ncol = 2),
    fill = guide_legend('Segment classification', title.position = 'top', ncol = 2)
  ) +
  theme_void(base_family = 'serif', base_size = 12) +
  theme(
    legend.position = 'top',
    legend.direction = 'horizontal',
    plot.margin = grid::unit(c(0, 0, 0, 0), units = 'pt'), 
    panel.grid.major = element_line(colour = 'white')
  )


jpeg('C:/Users/Marcus.SCCWRP2K/Desktop/sb_streamclass.jpeg', height = 8, width = 6, units = 'in', res = 400, family = 'serif')
# pdf('C:/Users/Marcus.SCCWRP2K/Desktop/calires.pdf', height = 36, width = 24, family = 'serif')
p1
dev.off()

# rb3 by class ------------------------------------------------------------

strcls <- strcls_all[rbnds, ] %>% 
  dplyr::select(COMID, Ref10) 
toplo <- strcls %>% 
  st_simplify(dTolerance = 0.1, preserveTopology = T) %>%
  rename(strcls = Ref10) %>% 
  mutate(
    strcls = factor(strcls, 
                    levels = c('likely unconstrained', 'possibly unconstrained', 'possibly constrained', 'likely constrained')
    )
  ) %>% 
  filter(!is.na(strcls)) %>% 
  .[rbnds, ]

# comid class
p2 <- ggplot() + 
  geom_sf(data = rbnds) + 
  geom_sf(data = toplo, aes(colour = strcls, fill = strcls)) +
  scale_colour_manual(values = pal_exp(levels(toplo$strcls))) +
  scale_fill_manual(values = pal_exp(levels(toplo$strcls))) +
  guides(
    colour = guide_legend('Segment classification', title.position = 'top', ncol = 2),
    fill = guide_legend('Segment classification', title.position = 'top', ncol = 2)
  ) +
  theme_void(base_family = 'serif', base_size = 12) +
  theme(
    legend.position = 'top',
    legend.direction = 'horizontal',
    plot.margin = grid::unit(c(0, 0, 0, 0), units = 'pt'), 
    panel.grid.major = element_line(colour = 'white')
  )


jpeg('C:/Users/Marcus.SCCWRP2K/Desktop/rb3_streamclass.jpeg', height = 8, width = 6, units = 'in', res = 400, family = 'serif')
# pdf('C:/Users/Marcus.SCCWRP2K/Desktop/calires.pdf', height = 36, width = 24, family = 'serif')
p2
dev.off()

# sb county, by median prediction -----------------------------------------

strcls <- strcls_all[sbco, ] %>% 
  dplyr::select(COMID, qt50) 
toplo <- strcls %>% 
  # st_simplify(dTolerance = 0.01, preserveTopology = T) %>% 
  filter(!is.na(qt50)) %>% 
  filter(qt50 >= 0) %>% 
  .[sbco, ]

p1 <- ggplot() + 
  geom_sf(data = sbco) + 
  geom_sf(data = toplo, aes(colour = qt50)) +
  scale_colour_gradientn('Median CSCI prediction', colours = pal(seq(0, 1.4, by = 0.1))) +
  theme_void(base_family = 'serif', base_size = 12) +
  theme(
    legend.position = 'top',
    legend.direction = 'horizontal',
    plot.margin = grid::unit(c(10, 0, 0, 0), units = 'pt'), 
    panel.grid.major = element_line(colour = 'white')
  )


jpeg('C:/Users/Marcus.SCCWRP2K/Desktop/sb_medianpred.jpeg', height = 8, width = 6, units = 'in', res = 400, family = 'serif')
# pdf('C:/Users/Marcus.SCCWRP2K/Desktop/calires.pdf', height = 36, width = 24, family = 'serif')
p1
dev.off()

# rb 3 boundary
rbnds <- st_read('S:/Spatial_Data/RWQCBdistricts/rwqcbnda.shp') %>% 
  st_transform(crs = prj) %>% 
  filter(RB == 3) %>% 
  st_simplify(dTolerance = 0.001, preserveTopology = T)

# rb3 median prediction ---------------------------------------------------

strcls <- strcls_all[rbnds, ] %>% 
  dplyr::select(COMID, qt50) 
toplo <- strcls %>% 
  st_simplify(dTolerance = 0.01, preserveTopology = T) %>%
  filter(!is.na(qt50)) %>% 
  filter(qt50 >= 0) %>% 
  .[rbnds, ]

# comid class
p2 <- ggplot() + 
  geom_sf(data = rbnds) + 
  geom_sf(data = toplo, aes(colour = qt50)) +
  scale_colour_gradientn('Median CSCI prediction', colours = pal(seq(0, 1.4, by = 0.1))) +
  theme_void(base_family = 'serif', base_size = 12) +
  theme(
    legend.position = 'top',
    legend.direction = 'horizontal',
    plot.margin = grid::unit(c(10, 0, 0, 0), units = 'pt'), 
    panel.grid.major = element_line(colour = 'white')
  )


jpeg('C:/Users/Marcus.SCCWRP2K/Desktop/rb3_medianpred.jpeg', height = 8, width = 6, units = 'in', res = 400, family = 'serif')
# pdf('C:/Users/Marcus.SCCWRP2K/Desktop/calires.pdf', height = 36, width = 24, family = 'serif')
p2
dev.off()
