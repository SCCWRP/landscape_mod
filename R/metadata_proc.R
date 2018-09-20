# create the data

library(tidyverse)
library(sf)
library(leaflet)
library(RColorBrewer)
library(rgdal)

# data(caliclsplo)
load(file = 'data/calipsa.RData')
load(file = 'data/comid_prd.RData')

source('R/funcs.R')

prj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# comid predictions to join with calinhd
comid_prd <- comid_prd %>% 
  dplyr::select(COMID, matches('^core|^COMID$'))

# all cali hydro, really simplified
calinhd <- readOGR('S:/Spatial_Data/NHDPlus/All_Ca_NHDPlus.shp') %>% 
  spTransform(prj) %>% 
  st_as_sf %>% 
  dplyr::select(COMID) %>%
  left_join(comid_prd, by = 'COMID')

# get biological condition expectations
cls1 <- getcls2(calinhd, thrsh = 0.63, tails = 0.1, modls = 'core')
cls2 <- getcls2(calinhd, thrsh = 0.79, tails = 0.1, modls = 'core')
cls3 <- getcls2(calinhd, thrsh = 0.92, tails = 0.1, modls = 'core')

calicls <- calinhd %>%   
  replace(is.na(.), -999999) %>% 
  left_join(cls1, by = 'COMID') %>%
  rename(
    Ref01 = strcls
  ) %>% 
  left_join(cls2, by = 'COMID') %>% 
  rename(
    Ref10 = strcls
  ) %>% 
  left_join(cls3, by = 'COMID') %>% 
  rename(
    Ref30 = strcls
  ) %>% 
  select(-strcls_int.x, - strcls_int.y, -strcls_int) %>% 
  rename_at(vars(contains('core')), funs(gsub('^core0\\.', 'qt', .)))

st_write(calicls, 'data/metadata/strm_constraints.shp', delete_layer = TRUE)

# ##
# # create metadata
# 
# library(EML)
# f <- system.file("examples/hf205-methods.md", package = "EML")
# set_methods(methods_file = f)
# 
# title <- 'Constrained streams for biological integrity in California'
# abstract <- as(set_TextType("data/metadata/abstract.md"), "abstract")
# intellectualRights <- "This work is licensed under a Creative Commons Attribution 4.0 International License."
# pubDate <- '2018'
# keywordSet <- c(
#   new("keywordSet", 
#       keyword = c('bioassessment', 'biotic integrity', 'streams', 'urbanization', 'modified channels', 'landscape stressors'))
# )
# marcus <- as.person("Marcus W. Beck <marcusb@sccwrp.org> [cre]")
# creator <- as(marcus, "creator")
# contact <- as(marcus, "contact")
# setwd('L:/Channels in developed landscapes_RM/Marcus/landscape_mod/')
# methods <- set_methods(methods_file = "data/metadata/methods.md")


