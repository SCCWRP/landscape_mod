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
library(randomForest)

source('R/funcs.R')
prj <- geo_wgs84

# ##
# # csci data, statewide
# csci_raw <- read_csv('data/raw/csci_raw.txt')
# save(csci_raw, file = 'data/csci_raw.RData', compress = 'xz')
# 
# ##
# # psa 6
# calipsa <- readOGR('S:/Spatial_Data/RCMP_needs editting/Inputs/PSA6_090111/PSA6_2011.shp') %>% 
#   spTransform(prj) %>% 
#   st_as_sf
# save(calipsa, file = 'data/calipsa.RData', compress = 'xz')
# 
# ##
# # psa labels by centroid, deserts modoc has two
# psalab <- calipsa %>% 
#   as('Spatial') %>% 
#   rmapshaper::ms_explode() %>% 
#   st_as_sf %>% 
#   mutate(
#     AREA = st_area(.), 
#     AREA = gsub('\\sm\\^2$', '', AREA), 
#     AREA = as.numeric(AREA)
#   ) %>% 
#   group_by(PSA6) %>% 
#   top_n(2, AREA) %>% 
#   ungroup %>% 
#   filter(!(PSA6 == 'South Coast' & AREA < 1e10)) %>% 
#   filter(!(PSA6 == 'Chaparral' & AREA < 1e10)) %>% 
#   st_centroid
# 
# psalab <- psalab %>% 
#   st_coordinates %>% 
#   data.frame %>% 
#   rename(
#     long = X, 
#     lat = Y
#   ) %>% 
#   mutate(
#     Region = psalab$PSA6,
#     Region = factor(Region, 
#       levels = c('Central Valley', 'Chaparral', 'Deserts Modoc', 'North Coast', 'Sierra Nevada', 'South Coast'),
#       labels = c('CV', 'CH', 'DM', 'NC', 'SN', 'SC')
#     )
#   )
# 
# save(psalab, file = 'data/psalab.RData', compress = 'xz')
# 
# ##
# # cali NHD simplify, fortify, and save
# 
# load(file = 'data/comid_prd.RData')
# 
# # all cali hydro, really simplified
# calinhd <- readOGR('S:/Spatial_Data/NHDPlus/All_Ca_NHDPlus.shp') %>% 
#   spTransform(prj) %>% 
#   st_as_sf %>% 
#   st_simplify(dTolerance = 0.5, preserveTopology = T)
# 
# # comid predictions to join with calinhd
# comid_prd <- comid_prd %>% 
#   dplyr::select(COMID, core0.50)
# 
# # fortified calinhd, joind with comid pred
# nhdplo <- calinhd %>% 
#   filter(COMID %in% unique(comid_prd$COMID)) %>% 
#   dplyr::select(COMID) %>% 
#   as('Spatial')
# comidid <- nhdplo@data %>% 
#   dplyr::select(COMID) %>% 
#   rownames_to_column('id')
# nhdplo <- nhdplo %>% 
#   fortify %>% 
#   left_join(comidid, by = 'id') %>%  
#   inner_join(comid_prd, by = 'COMID')
# 
# save(nhdplo, file = 'data/nhdplo.RData', compress = 'xz')
# 
######
# cali nhd stream classes, simplified, fortified, and saved

data(comid_prd)

# comid predictions to join with calinhd
comid_prd <- comid_prd %>% 
  dplyr::select(COMID, matches('^core|^COMID$'))

# all cali hydro, really simplified
calinhd <- readOGR('S:/Spatial_Data/NHDPlus/All_Ca_NHDPlus.shp') %>% 
  spTransform(prj) %>% 
  st_as_sf %>% 
  st_simplify(dTolerance = 0.5, preserveTopology = T) %>% 
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

save(calicls, file = 'data/calicls.RData', compress = 'xz')
save(caliclsplo, file = 'data/caliclsplo.RData', compress = 'xz')

######
# all cali CSCI site expectations

# all comid performance
data(csci_comid)

# all cali hydro, really simplified
calinhd <- readOGR('S:/Spatial_Data/NHDPlus/All_Ca_NHDPlus.shp') %>% 
  spTransform(prj) %>% 
  st_as_sf %>% 
  st_simplify(dTolerance = 0.5, preserveTopology = T) %>% 
  filter(COMID %in% unique(comid_prd$COMID)) %>% 
  dplyr::select(COMID) %>%
  inner_join(comid_prd, by = 'COMID')

# format csci data for site_exp function, takes averages of repeats at each site
csci_comid <- csci_comid %>% 
  dplyr::select(COMID, StationCode, CSCI, SampleDate, FieldReplicate, New_Lat, New_Long) %>% 
  group_by(COMID, StationCode) %>% 
  summarise(
    csci = mean(CSCI, na.rm = T), 
    lat = mean(New_Lat, na.rm = T),
    long = mean(New_Long, na.rm = T)
  ) %>% 
  ungroup

# get site expectations
caliexp <- site_exp(calinhd, csci_comid, thrsh = 0.79, tails = 0.1, modls = 'core')

save(caliexp, file = 'data/caliexp.RData', compress = 'xz')

######
# get stream length in each stream class by PSA

load(file = 'data/calicls.RData')
load(file = 'data/caliexp.RData')
load(file = 'data/csci_comid.RData')
load(file = 'data/calipsa.RData')

strclslen <- calicls %>% 
  dplyr::select(COMID, strcls) %>% 
  st_intersection(calipsa)

lens <- strclslen %>% 
  as('Spatial') %>% 
  sp::SpatialLinesLengths(longlat = T)

strclslen <- strclslen %>% mutate(lens = lens) 

save(strclslen, file = 'data/strclslen.RData', compress = 'xz')

# ######
# # expected range of scores for urban, ag, other by region
# data(calicls)
# data(strclslen)
# 
# # streamcat data
# strmcat <- rbind(read.csv("Z:/MarcusBeck/Landscape models from rafi/Streamcat_v2_AllCOMID_030117/exp_1.csv", stringsAsFactors = F),
#                  read.csv("Z:/MarcusBeck/Landscape models from rafi/Streamcat_v2_AllCOMID_030117/exp_2.csv", stringsAsFactors = F),
#                  read.csv("Z:/MarcusBeck/Landscape models from rafi/Streamcat_v2_AllCOMID_030117/exp_3.csv", stringsAsFactors = F))
# comid.nats <- read.csv("Z:/MarcusBeck/Landscape models from rafi/ALL_COMID_Nats.csv", stringsAsFactors = F)
# strmcat <- plyr::join(strmcat, comid.nats[,setdiff(names(comid.nats), "WsAreaSqKm")])
# 
# strmcat$TotUrb2011Ws<-  rowSums(strmcat[,c("PctUrbOp2011Ws","PctUrbLo2011Ws","PctUrbMd2011Ws","PctUrbHi2011Ws")])
# strmcat$TotUrb2011Cat<-  rowSums(strmcat[,c("PctUrbOp2011Cat","PctUrbLo2011Cat","PctUrbMd2011Cat","PctUrbHi2011Cat")])
# strmcat$TotUrb2011WsRp100<-  rowSums(strmcat[,c("PctUrbOp2011WsRp100","PctUrbLo2011WsRp100","PctUrbMd2011WsRp100","PctUrbHi2011WsRp100")])
# strmcat$TotUrb2011CatRp100<-  rowSums(strmcat[,c("PctUrbOp2011CatRp100","PctUrbLo2011CatRp100","PctUrbMd2011CatRp100","PctUrbHi2011CatRp100")])
# 
# strmcat$TotAg2011Ws<-  rowSums(strmcat[,c("PctHay2011Ws","PctCrop2011Ws")])
# strmcat$TotAg2011Cat<-  rowSums(strmcat[,c("PctHay2011Cat","PctCrop2011Cat")])
# strmcat$TotAg2011WsRp100<-  rowSums(strmcat[,c("PctHay2011WsRp100","PctCrop2011WsRp100")])
# strmcat$TotAg2011CatRp100<-  rowSums(strmcat[,c("PctHay2011CatRp100","PctCrop2011CatRp100")])
# 
# strmcat <- strmcat %>% 
#   dplyr::select(COMID, TotUrb2011Ws, TotAg2011Ws)
# 
# # PSA by all comid
# psaall <- strclslen %>% 
#   dplyr::select(COMID, PSA6)
# st_geometry(psaall) <- NULL
# 
# # data to summarize
# scrdist <- calicls
# st_geometry(scrdist) <- NULL
# scrdist <- scrdist %>% 
#   dplyr::select(-strcls, -strcls_int) %>% 
#   left_join(strmcat, by = 'COMID') %>% 
#   left_join(psaall, by = 'COMID') %>% 
#   filter(!is.na(PSA6)) %>% 
#   rename(Region = PSA6)
# 
# # get range of scores for average urban, ag, open locations, by region
# scrdistreg <- scrdist %>% 
#   group_by(Region) %>% 
#   nest %>% 
#   mutate(grps = purrr::map(data, function(x){
#     
#     # get kmeans, centers
#     tomod <- x %>% 
#       mutate(
#         Urb = log10(1 + TotUrb2011Ws), 
#         Ag = log10(1 + TotAg2011Ws)
#       ) %>% 
#       dplyr::select(Urb, Ag) 
#     ngrps <- 8
#     kmod <- kmeans(tomod, centers = ngrps, nstart = 10, iter.max = 100, algorithm = 'MacQueen')
#     grps <- kmod$cluster
#     cent <- kmod$centers %>% 
#       data.frame
#     
#     # find the typical ag, urb group
#     urb <- which.max(cent[, 'Urb'])
#     ag <- which.max(cent[, 'Ag'])
#     oth <- which.min(rowSums(cent))
#     
#     grplabs <- list(
#       urb = urb,
#       ag = ag,
#       oth = oth
#     )
#     
#     ests <- x %>% 
#       mutate(
#         grps = grps
#       ) %>% 
#       filter(grps %in% unlist(grplabs)) %>% 
#       mutate(
#         grps = factor(grps, levels = unlist(grplabs), labels = names(grplabs))
#       ) %>% 
#       group_by(grps) %>% 
#       summarise(
#         lo05 = mean(full0.05, na.rm = T), 
#         hi95 = mean(full0.95, na.rm = T),
#         lo25 = mean(full0.25, na.rm = T), 
#         hi75 = mean(full0.75, na.rm = T), 
#         lo45 = mean(full0.45, na.rm = T),
#         hi55 = mean(full0.55, na.rm = T)
#       )
#     
#     return(ests)
#     
#   })) %>% 
#   dplyr::select(-data) %>% 
#   unnest %>% 
#   mutate(
#     Region = factor(Region, 
#          levels = c('Central Valley', 'Chaparral', 'Deserts Modoc', 'North Coast', 'Sierra Nevada', 'South Coast'),
#          labels = c('CV', 'CH', 'DM', 'NC', 'SN', 'SC')),
#     Region = as.character(Region), 
#     grps = as.character(grps)
#   )
# 
# # typical statewide
# scrdistall <- scrdist %>% 
#   mutate(Region = 'Statewide') %>% 
#   group_by(Region) %>% 
#   nest %>% 
#   mutate(grps = purrr::map(data, function(x){
#     
#     # get kmeans, centers
#     tomod <- x %>% 
#       mutate(
#         Urb = log10(1 + TotUrb2011Ws), 
#         Ag = log10(1 + TotAg2011Ws)
#       ) %>% 
#       dplyr::select(Urb, Ag) 
#     ngrps <- 8
#     kmod <- kmeans(tomod, centers = ngrps, nstart = 10, iter.max = 100, algorithm="MacQueen")
#     grps <- kmod$cluster
#     cent <- kmod$centers %>% 
#       data.frame
#     
#     # find the typical ag, urb group
#     urb <- which.max(cent[, 'Urb'])
#     ag <- which.max(cent[, 'Ag'])
#     oth <- which.min(rowSums(cent))
#     
#     grplabs <- list(
#       urb = urb,
#       ag = ag,
#       oth = oth
#     )
#     
#     ests <- x %>% 
#       mutate(
#         grps = grps
#       ) %>% 
#       filter(grps %in% unlist(grplabs)) %>% 
#       mutate(
#         grps = factor(grps, levels = unlist(grplabs), labels = names(grplabs))
#       ) %>% 
#       group_by(grps) %>% 
#       summarise(
#         lo05 = mean(full0.05, na.rm = T), 
#         hi95 = mean(full0.95, na.rm = T),
#         lo25 = mean(full0.25, na.rm = T), 
#         hi75 = mean(full0.75, na.rm = T), 
#         lo45 = mean(full0.45, na.rm = T),
#         hi55 = mean(full0.55, na.rm = T)
#       )
#     
#     return(ests)
#     
#   })) %>% 
#   dplyr::select(-data) %>% 
#   unnest %>% 
#   mutate(
#     Region = as.character(Region),
#     grps = as.character(grps)
#     )
# 
# 
# typscrs <- bind_rows(scrdistall, scrdistreg)
# 
# save(typscrs, file = 'data/typscrs.RData', compress = 'xz')
# 
##
# rf model importance for constraints in each region

load(file = 'data/caliclsplo.RData')
load(file = 'data/strclslen.RData')

# streamcat data
strmcat <- rbind(read.csv("Z:/MarcusBeck/Landscape models from rafi/Streamcat_v2_AllCOMID_030117/exp_1.csv", stringsAsFactors = F),
                 read.csv("Z:/MarcusBeck/Landscape models from rafi/Streamcat_v2_AllCOMID_030117/exp_2.csv", stringsAsFactors = F),
                 read.csv("Z:/MarcusBeck/Landscape models from rafi/Streamcat_v2_AllCOMID_030117/exp_3.csv", stringsAsFactors = F))
comid.nats <- read.csv("Z:/MarcusBeck/Landscape models from rafi/ALL_COMID_Nats.csv", stringsAsFactors = F)
strmcat <- plyr::join(strmcat, comid.nats[,setdiff(names(comid.nats), "WsAreaSqKm")])

# PSA by all comid
psaall <- strclslen %>% 
  filter(!is.na(strcls)) %>% 
  dplyr::select(COMID, PSA6)
st_geometry(psaall) <- NULL

allenv <- caliclsplo %>% 
  dplyr::select(COMID, strcls) %>% 
  unique %>% 
  left_join(strmcat, by = 'COMID') %>% 
  left_join(psaall, by = 'COMID') %>% 
  filter(!is.na(PSA6)) %>% 
  mutate(strcls = gsub('^possibly\\s|^likely\\s', '', strcls))


# remove wshed vars, focus on watershed only, remove COMID from model
cnstrfrst <- allenv[, !grepl('Cat|^COMID$', names(allenv))] %>% 
  group_by(PSA6) %>% 
  nest %>% 
  mutate(mods = purrr::map(data, function(x){
    
    tmp <- randomForest(as.factor(strcls) ~ .,
                        data = x, 
                        importance = TRUE, 
                        ntree = 1000, na.action = na.omit)
    
    tmp$importance
    
  })) %>% 
  dplyr::select(-data)

save(cnstrfrst, file = 'data/cnstrfrst.RData', compress = 'xz')

# ##
# # SMC watershed 
# 
# shd_pth <- 'S:/Spatial_Data/SMCBasefiles/Boundaries/SMCSheds/SMCSheds2009.shp'
# 
# # watersheds
# shed <- readOGR(shd_pth) %>%
#   spTransform(CRS(prj)) %>%
#   st_as_sf %>% 
#   filter(SMC_Name %in% 'San Gabriel') %>% 
#   dplyr::select(SMC_Name)
# 
# save(shed, file = '/data/shed.RData', compress = 'xz')
# 
# ##
# # land use data from gap
# # https://gapanalysis.usgs.gov/gaplandcover/data/download/
# 
# # master raster
# gapland <- raster('Z:/MarcusBeck/GIS/gaplf2011lc_v30_CA/gaplf2011lc_v30_ca.tif')
# 
# # # gap key
# # # ag is 555 - 557, urban is 580 - 584
# # gapkey <- read_delim(
# #   'C:/Users/Marcus.SCCWRP2K/Desktop/gaplf2011lc_v30_CA/GAP_LANDFIRE_National_Terrestrial_Ecosystems_2011_Attributes.txt',
# #   delim = '\t')
# 
# # reclassify matrix, right closed
# rcmat <- c(-1, 554, 1, 554, 557, 2, 557, 579, 1, 579, 584, 3) %>%
#   matrix(., ncol = 3, byrow = T)
# 
# # get cali shapefile in format for clip
# data(calishp)
# toprj <- crs(gapland) %>% 
#   as.character
# calishp <- calishp %>% 
#   st_transform(toprj)
# 
# # reclassify, aggregate, vectorize, clip, simplify, dissolve by feature, transform
# # takes a minute
# ludat <- gapland %>% 
#   reclassify(rcmat) %>% 
#   aggregate(fact = 10, fun = modal) %>% 
#   rasterToPolygons(dissolve = T) %>% 
#   st_as_sf %>% 
#   st_intersection(calishp) %>% 
#   ms_simplify(as(., 'Spatial')) %>% 
#   st_as_af %>% 
#   st_buffer(0) %>% 
#   group_by(layer) %>% 
#   summarize %>% 
#   st_cast %>% 
#   st_transform(prj)
# 
# save(ludat, file = 'data/ludat.RData', compress = 'xz')
# 
# ######
# # SGR lu clip
# 
# # master raster
# gapland <- raster('Z:/MarcusBeck/GIS/gaplf2011lc_v30_CA/gaplf2011lc_v30_ca.tif')
# 
# # # gap key
# # # ag is 555 - 557, urban is 580 - 584
# # gapkey <- read_delim(
# #   'C:/Users/Marcus.SCCWRP2K/Desktop/gaplf2011lc_v30_CA/GAP_LANDFIRE_National_Terrestrial_Ecosystems_2011_Attributes.txt',
# #   delim = '\t')
# 
# # reclassify matrix, right closed
# rcmat <- c(-1, 161, 1, 161, 162, 2, 162, 182, 1, 182, 183, 1, 183, 580, 1, 580, 582, 3, 582, 583, 4, 583, 584, 5) %>%
#   matrix(., ncol = 3, byrow = T)
# 
# # master raster
# gapland <- raster('Z:/MarcusBeck/GIS/gaplf2011lc_v30_CA/gaplf2011lc_v30_ca.tif')
# 
# 
# # get shed shapefile in format for clip
# data(shed)
# toprj <- crs(gapland) %>% 
#   as.character
# tocrp <- shed %>% 
#   st_transform(toprj)# %>% 
#   # st_buffer(dist = 50000)
# 
# ludat <- tocrp %>% 
#   as('Spatial') %>% 
#   raster::crop(gapland, .) %>% 
#   aggregate(fact = 4, fun = modal)
# ludat <- tocrp %>% 
#   as('Spatial') %>% 
#   raster::mask(ludat, .) %>% 
#   reclassify(rcmat) %>% 
#   projectRaster(crs = prj)
# 
# ludat@data@values <- ludat@data@values %>% 
#   pmax(., 1) %>% 
#   pmin(., 5) %>% 
#   round(0)
# 
# sgrlu <- ludat
# save(sgrlu, file = 'data/sgrlu.RData', compress = 'xz')
# 
# ######
# # sensitivity analysis, statewide
# 
# data(comid_prd)
# data(csci_comid)
# data(strclslen)
# 
# source('R/funcs.R')
# prj <- geo_wgs84
# 
# # calipsa
# calipsa <- readOGR('S:/Spatial_Data/RCMP_needs editting/Inputs/PSA6_090111/PSA6_2011.shp') %>% 
#   spTransform(prj) %>% 
#   st_as_sf
# 
# # comid predictions to join with calinhd
# comid_prd <- comid_prd %>% 
#   dplyr::select(COMID, matches('^full|^COMID$'))
# 
# # all cali hydro, really simplified
# calinhd <- readOGR('S:/Spatial_Data/NHDPlus/NHDPLusCalifornia/NHDPlusCalifornia.shp') %>% 
#   spTransform(prj) %>% 
#   st_as_sf %>% 
#   st_simplify(dTolerance = 0.5, preserveTopology = T) %>% 
#   filter(COMID %in% unique(comid_prd$COMID)) %>% 
#   dplyr::select(COMID) %>%
#   inner_join(comid_prd, by = 'COMID')
# 
# # format csci data for site_exp function, takes averages of repeats at each site
# csci_agg <- csci_comid %>% 
#   dplyr::select(COMID, StationCode, CSCI, SampleDate, FieldReplicate, New_Lat, New_Long) %>% 
#   group_by(COMID, StationCode) %>% 
#   summarise(
#     csci = mean(CSCI, na.rm = T), 
#     lat = mean(New_Lat, na.rm = T),
#     long = mean(New_Long, na.rm = T)
#   ) %>% 
#   ungroup
# 
# # psa by csci comid
# cscipsa <- csci_comid %>% 
#   dplyr::select(COMID, PSA6c) %>% 
#   unique %>% 
#   mutate(COMID = as.character(COMID))
# 
# # search grid
# grd_chk <- list(
#   thrsh = c(0.63, 0.79, 0.92), 
#   tails = seq(0.05, 0.45, by = 0.05)
#   ) %>% 
#   cross_df
# 
# # setup parallel backend
# ncores <- detectCores() - 3  
# cl<-makeCluster(ncores)
# registerDoParallel(cl)
# strt<-Sys.time()
# 
# sensres <- foreach(i = 1:nrow(grd_chk), 
#                .export = c('calipsa', 'comid_prd', 'calinhd', 'csci_agg', 'cscipsa'), 
#                .packages = c('sp', 'mapview', 'leaflet', 'tidyverse', 'raster', 'maptools', 'rgdal', 'sf')) %dopar% {
#                  
#   sink('log.txt')
#   cat(i, 'of', nrow(grd_chk), '\n')
#   print(Sys.time()-strt)
#   sink()
#   
#   source("R/funcs.R")
#   
#   thrshi <- grd_chk[[i, 'thrsh']]
#   tailsi <- grd_chk[[i, 'tails']]
#   
#   # get biological condition expectations
#   cls <- getcls2(calinhd, thrsh = thrshi, tails = tailsi, modls = 'full')
#   
#   calicls <- calinhd %>% 
#     left_join(cls, by = 'COMID')
#   
#   ######
#   # all cali CSCI site expectations
#   
#   # get site expectations
#   caliexp <- site_exp(calinhd, csci_agg, thrsh = thrshi, tails = tailsi, modls = 'full')
#   
#   ######
#   # get stream length in each stream class by PSA
#   
#   strclsleni <- calicls %>% 
#     dplyr::select(COMID, strcls) %>% 
#     st_intersection(calipsa)
#     
#   lens <- strclsleni %>% 
#     as('Spatial') %>% 
#     sp::SpatialLinesLengths(longlat = T)
#     
#   strclsleni <- strclsleni %>% mutate(lens = lens) 
#   
#   st_geometry(strclsleni) <- NULL
#   lentot <- strclsleni %>% 
#     filter(!is.na(strcls)) %>% 
#     mutate(
#     Region = factor(PSA6, 
#                    levels = c('Central Valley', 'Chaparral', 'Deserts Modoc', 'North Coast', 'Sierra Nevada', 'South Coast'),
#                    labels = c('CV', 'CH', 'DM', 'NC', 'SN', 'SC')
#       )
#     ) %>% 
#     group_by(Region, strcls) %>% 
#     summarise(len = sum(lens)) %>% 
#     spread(strcls, len, fill = 0) %>% 
#     gather('strcls', 'len', -Region) %>% 
#     group_by(Region) %>% 
#     mutate(
#     perc = 100 * len / sum(len), 
#     strcls = factor(strcls, levels = c('likely constrained', 'possibly constrained', 'possibly unconstrained', 'likely unconstrained'))
#     )
#   
#   lentotall <- strclsleni %>% 
#     filter(!is.na(strcls)) %>% 
#     group_by(strcls) %>% 
#     summarise(len = sum(lens)) %>%
#     ungroup %>% 
#     mutate(
#     Region = 'statewide',
#     perc = 100 * len / sum(len), 
#     strcls = factor(strcls, levels = c('likely constrained', 'possibly constrained', 'possibly unconstrained', 'likely unconstrained'))
#     ) 
#   
#   totab1 <- bind_rows(lentotall, lentot)
#   
#   # format table of site counts by region
#   reltot <- caliexp %>% 
#     dplyr::select(COMID, csci, StationCode, strcls, perf) %>% 
#     filter(!is.na(strcls)) %>% 
#     left_join(cscipsa, by = 'COMID') %>% 
#     rename(Region = PSA6c) %>% 
#     group_by(Region, perf) %>% 
#     summarise(
#     n = n(), 
#     csciave = mean(csci, na.rm = T), 
#     cscisd = sd(csci, na.rm = T)
#     ) %>% 
#     group_by(Region) %>% 
#     mutate(
#     perc = n / sum(n)
#     ) %>% 
#     complete(Region, perf)
#   
#   # whole state
#   reltotall <- caliexp %>% 
#     dplyr::select(COMID, csci, StationCode, strcls, perf) %>% 
#     filter(!is.na(strcls)) %>% 
#     left_join(cscipsa, by = 'COMID') %>% 
#     rename(Region = PSA6c) %>% 
#     group_by(perf) %>% 
#     summarise(
#     n = n(), 
#     csciave = mean(csci, na.rm = T), 
#     cscisd = sd(csci, na.rm = T)
#     ) %>% 
#     ungroup %>% 
#     mutate(
#     perc = n / sum(n)
#     ) %>% 
#     mutate(Region = 'Statewide') 
#   
#   totab2 <- bind_rows(reltotall, reltot)
#   
#   out <- list(totab1, totab2)
#   return(out)
#   
# }
# 
# save(sensres, file = 'data/sensres.RData', compress = 'xz')
# 
# ######
# # sensitivity analysis, SGR only
# 
# data(scrs)
# data(spat)
# 
# source('R/funcs.R')
# 
# # search grid
# grd_chk <- list(
#   thrsh = c(0.63, 0.79, 0.92),
#   tails = seq(0.05, 0.45, by = 0.05)
# ) %>%
#   cross_df
# 
# # setup parallel backend
# ncores <- detectCores() - 3
# cl<-makeCluster(ncores)
# registerDoParallel(cl)
# strt<-Sys.time()
# 
# senssgr <- foreach(i = 1:nrow(grd_chk),
#                    .export = c('scrs', 'spat'),
#                    .packages = c('sp', 'mapview', 'leaflet', 'tidyverse', 'raster', 'maptools', 'rgdal', 'sf')) %dopar% {
# 
#    sink('log.txt')
#    cat(i, 'of', nrow(grd_chk), '\n')
#    print(Sys.time()-strt)
#    sink()
# 
#    source("R/funcs.R")
# 
#    thrshi <- grd_chk[[i, 'thrsh']]
#    tailsi <- grd_chk[[i, 'tails']]
# 
#    # get biological condition expectations
#    cls <- getcls2(spat, thrsh = thrshi, tails = tailsi, modls = 'full')
# 
#    spat2 <- spat %>%
#      left_join(cls, by = 'COMID')
# 
#    ######
#    # get stream length total and percent
# 
#    strclsleni <- spat2 %>%
#      dplyr::select(COMID, strcls)
# 
#    lens <- strclsleni %>%
#      as('Spatial') %>%
#      sp::SpatialLinesLengths(longlat = T)
# 
#    strclsleni <- strclsleni %>% mutate(lens = lens)
# 
#    st_geometry(strclsleni) <- NULL
#    lentot <- strclsleni %>%
#      filter(!is.na(strcls)) %>%
#      group_by(strcls) %>%
#      summarise(len = sum(lens)) %>%
#      spread(strcls, len, fill = 0) %>%
#      gather('strcls', 'len') %>%
#      mutate(
#        perc = 100 * len / sum(len),
#        strcls = factor(strcls, levels = c('likely constrained', 'possibly constrained', 'possibly unconstrained', 'likely unconstrained'))
#      )
# 
#    out <- lentot
#    return(out)
# 
#   }
# 
# save(senssgr, file = 'data/senssgr.RData', compress = 'xz')
# 
# 
# ######
# # create data from Rafi
# # create csci comid data, rf core and full models for csci quantiles, quantile preds for all comid in CA
# # original file from Z:/MarcusBeck/Landscape models from rafi/modlu_120117.R
# 
# setwd("Z:/MarcusBeck/Landscape models from rafi/")
# library(plyr)
# library(dplyr)
# library(reshape2)
# 
# # csci<-read.csv("CSCI_LU_temp_022817.csv", stringsAsFactors = F) #Load data
# csci<-join(read.csv("CSCI_LU_temp_040617.csv", stringsAsFactors = F), #Load data
#            read.csv("CSCI_Nat_temp_061217.csv", stringsAsFactors = F))
# nat.vars.full<-setdiff(names(read.csv("CSCI_Nat_temp_061217.csv")), names(read.csv("CSCI_LU_temp_040617.csv")))
# nat.vars.full<-unique(c("WsAreaSqKm",nat.vars.full))
# 
# csci$SampleID<-paste(csci$StationCode, csci$SampleDate, csci$CollectionMethodCode, csci$FieldReplicate, sep="_")
# csci<-csci[which(!is.na(csci$PctFrstLossWs)),] #Drop rows with missing STREAMCAT data
# 
# 
# csci$RdDensCatRp100[is.na(csci$RdDensCatRp100)]<-0
# 
# csci$TotUrb2011Ws<-  rowSums(csci[,c("PctUrbOp2011Ws","PctUrbLo2011Ws","PctUrbMd2011Ws","PctUrbHi2011Ws")])
# csci$TotUrb2011Cat<-  rowSums(csci[,c("PctUrbOp2011Cat","PctUrbLo2011Cat","PctUrbMd2011Cat","PctUrbHi2011Cat")])
# csci$TotUrb2011WsRp100<-  rowSums(csci[,c("PctUrbOp2011WsRp100","PctUrbLo2011WsRp100","PctUrbMd2011WsRp100","PctUrbHi2011WsRp100")])
# csci$TotUrb2011CatRp100<-  rowSums(csci[,c("PctUrbOp2011CatRp100","PctUrbLo2011CatRp100","PctUrbMd2011CatRp100","PctUrbHi2011CatRp100")])
# 
# csci$TotAg2011Ws<-  rowSums(csci[,c("PctHay2011Ws","PctCrop2011Ws")])
# csci$TotAg2011Cat<-  rowSums(csci[,c("PctHay2011Cat","PctCrop2011Cat")])
# csci$TotAg2011WsRp100<-  rowSums(csci[,c("PctHay2011WsRp100","PctCrop2011WsRp100")])
# csci$TotAg2011CatRp100<-  rowSums(csci[,c("PctHay2011CatRp100","PctCrop2011CatRp100")])
# 
# gis.area<-read.csv("csci_AREA.csv", stringsAsFactors = F)
# 
# setdiff(csci$StationCode, gis.area$StationCode)
# #Let's look at PCA-space:
# 
# 
# 
# #COMIDs
# all.comid<-rbind(read.csv("Streamcat_v2_AllCOMID_030117/exp_1.csv", stringsAsFactors = F),
#                  read.csv("Streamcat_v2_AllCOMID_030117/exp_2.csv", stringsAsFactors = F),
#                  read.csv("Streamcat_v2_AllCOMID_030117/exp_3.csv", stringsAsFactors = F))
# comid.nats<-read.csv("ALL_COMID_Nats.csv", stringsAsFactors = F)
# all.comid<-join(all.comid, comid.nats[,setdiff(names(comid.nats), "WsAreaSqKm")])
# 
# all.comid$TotUrb2011Ws<-  rowSums(all.comid[,c("PctUrbOp2011Ws","PctUrbLo2011Ws","PctUrbMd2011Ws","PctUrbHi2011Ws")])
# all.comid$TotUrb2011Cat<-  rowSums(all.comid[,c("PctUrbOp2011Cat","PctUrbLo2011Cat","PctUrbMd2011Cat","PctUrbHi2011Cat")])
# all.comid$TotUrb2011WsRp100<-  rowSums(all.comid[,c("PctUrbOp2011WsRp100","PctUrbLo2011WsRp100","PctUrbMd2011WsRp100","PctUrbHi2011WsRp100")])
# all.comid$TotUrb2011CatRp100<-  rowSums(all.comid[,c("PctUrbOp2011CatRp100","PctUrbLo2011CatRp100","PctUrbMd2011CatRp100","PctUrbHi2011CatRp100")])
# 
# all.comid$TotAg2011Ws<-  rowSums(all.comid[,c("PctHay2011Ws","PctCrop2011Ws")])
# all.comid$TotAg2011Cat<-  rowSums(all.comid[,c("PctHay2011Cat","PctCrop2011Cat")])
# all.comid$TotAg2011WsRp100<-  rowSums(all.comid[,c("PctHay2011WsRp100","PctCrop2011WsRp100")])
# all.comid$TotAg2011CatRp100<-  rowSums(all.comid[,c("PctHay2011CatRp100","PctCrop2011CatRp100")])
# 
# all.comid$RdDensWsRp100<-all.comid$RoadDensityRipBuf100_CA.RdDensWsRp100
# all.comid$RdDensCatRp100<-all.comid$RoadDensityRipBuf100_CA.RdDensCatRp100
# 
# library(plyr)
# csci<-join(csci, gis.area)
# 
# 
# 
# library(dplyr)
# # csci2<-csci %>% group_by(COMID) %>% slice(which.max(AREA_SQKM))
# 
# csci$Largest<-
#   sapply(1:nrow(csci), function(i){
#     comid<-csci$COMID[i]
#     area<-csci$AREA_SQKM[i]
#     mydf<-csci[which(csci$COMID==comid),]
#     area==max(mydf$AREA_SQKM)
#   })
# 
# # csci[which(csci$COMID==8264738),c("StationCode","COMID","AREA_SQKM","Largest")]
# csci.full<-csci
# csci<-csci.full[csci.full$Largest,]
# 
# library(quantregForest)
# 
# 
# 
# #Create cal and val data
# sites<-unique(csci[,c("StationCode","PSA6c","PctImp2006Ws")]) #Create a sites DF
# library(plyr)
# sites.d<-ddply(sites, .(PSA6c), summarize, ImpT1=quantile(PctImp2006Ws, probs=c(0.25)),ImpT2=quantile(PctImp2006Ws, probs=c(0.5)),ImpT3=quantile(PctImp2006Ws, probs=c(0.75)))
# write.table(sites.d, "clipboard", sep="\t")
# 
# sites$RegImpQ<-   #Divide regions into thirds based on imperviousness
#   sapply(1:nrow(sites), function(i){
#     region<-sites$PSA6c[i]
#     imp<-sites$PctImp2006Ws[i]
#     imp.t1<-sites.d$ImpT1[which(sites.d$PSA6c==region)]
#     imp.t2<-sites.d$ImpT2[which(sites.d$PSA6c==region)]
#     imp.t3<-sites.d$ImpT3[which(sites.d$PSA6c==region)]
#     ifelse(imp<imp.t1,"T1", ifelse(imp<imp.t2,"T2", ifelse(imp<imp.t3,"T3","T4")))
#   })
# 
# table(sites$PSA6c, sites$RegImpQ)
# sites$Stratum<-paste(sites$PSA6c, sites$RegImpQ, sep="_") #Create strata by dividing each region into thirds
# table(sites$Stratum)
# 
# #Random assignment into cal (80%) and val (20%) data sets by Stratum field
# library(dplyr)
# sites.t<- sites %>% group_by(Stratum)
# set.seed(500)
# sites.cal<-sample_frac(sites.t, 0.75, replace=F)
# sites$SiteSet<-ifelse(sites$StationCode %in% sites.cal$StationCode, "Cal","Val")
# csci<-join(csci, sites[,c("StationCode","RegImpQ", "Stratum", "SiteSet")])
# 
# #
# #Select a single sample for assessment
# csci.t<- csci %>% group_by(StationCode)
# set.seed(501)
# samps.sel<-sample_n(csci.t, 1)
# csci$SelectedSample<-ifelse(csci$SampleID %in% samps.sel$SampleID, "Selected","NotSelected")
# 
# #PotentialVars
# metrics<-read.csv("metrics.csv", stringsAsFactors = F)
# 
# 
# core.candidates<-c("CanalDensCat","CanalDensWs", 
#                    "PctImp2006Cat","PctImp2006Ws","PctImp2006CatRp100","PctImp2006WsRp100",
#                    "TotUrb2011Ws","TotUrb2011Cat","TotUrb2011WsRp100","TotUrb2011CatRp100",
#                    "TotAg2011Ws","TotAg2011Cat","TotAg2011WsRp100","TotAg2011CatRp100",
#                    "RdDensCat", "RdDensWs", "RdDensCatRp100", "RdDensWsRp100", "RdCrsCat","RdCrsWs")
# core.candidates.disag<-c("CanalDensCat","CanalDensWs",
#                          "PctImp2006Cat","PctImp2006Ws","PctImp2006CatRp100","PctImp2006WsRp100",
#                          "PctUrbOp2011Cat", "PctUrbLo2011Cat",
#                          "PctUrbMd2011Cat", "PctUrbHi2011Cat", "PctHay2011Cat", "PctCrop2011Cat",
#                          "PctUrbOp2011Ws", "PctUrbLo2011Ws", "PctUrbMd2011Ws", "PctUrbHi2011Ws",
#                          "PctHay2011Ws", "PctCrop2011Ws", "PctUrbOp2011CatRp100", "PctUrbLo2011CatRp100",
#                          "PctUrbMd2011CatRp100", "PctUrbHi2011CatRp100", "PctHay2011CatRp100",
#                          "PctCrop2011CatRp100", "PctUrbOp2011WsRp100", "PctUrbLo2011WsRp100",
#                          "PctUrbMd2011WsRp100", "PctUrbHi2011WsRp100", "PctHay2011WsRp100",
#                          "PctCrop2011WsRp100",
#                          "RdDensCat", "RdDensWs", "RdDensCatRp100", "RdDensWsRp100", "RdCrsCat", "RdCrsSlpWtdCat","RdCrsWs","RdCrsSlpWtdWs")
# 
# core.candidates.mines.dams<-c(core.candidates,
#                               "MineDensCat", "MineDensWs", "MineDensCatRp100", "MineDensWsRp100",
#                               "DamDensCat", "DamDensWs",  "DamNrmStorM3Cat",  "DamNrmStorM3Ws")
# 
# core.candidates.disag.mines.dams<-c(core.candidates.disag,
#                                     "MineDensCat", "MineDensWs", "MineDensCatRp100", "MineDensWsRp100",
#                                     "DamDensCat", "DamDensWs",  "DamNrmStorM3Cat",  "DamNrmStorM3Ws")
# 
# stressors.atm<-c(#"NH4_2008Cat","NO3_2008Cat","InorgNWetDep_2008Cat","SN_2008Cat", #Cats have a small number of NAs
#   "NH4_2008Ws","NO3_2008Ws","InorgNWetDep_2008Ws","SN_2008Ws")
# 
# stressors.other<-c("PctNonAgIntrodManagVegCat","PctNonAgIntrodManagVegWs","PctNonAgIntrodManagVegCatRp100","PctNonAgIntrodManagVegWsRp100")
# 
# core.candidates.mines.dams.atm<-c(core.candidates.mines.dams, stressors.atm)
# core.candidates.disag.mines.dams.atm<-c(core.candidates.disag.mines.dams, stressors.atm)
# 
# core.candidates.mines.dams.atm.veg<-c(core.candidates.mines.dams.atm,stressors.other)
# nat.cands<-setdiff(nat.vars.full, c("MAST_2008","MAST_2009","MAST_2013","MAST_2014",
#                                     "MSST_2008","MSST_2009","MSST_2013","MSST_2014",
#                                     "MWST_2008","MWST_2009","MWST_2013","MWST_2014"))
# const.vars<-c("PctGlacTilClayCat","PctHydricCat","PctGlacTilLoamCat","PctGlacLakeFineCat","PctGlacLakeCrsCat","PctGlacTilLoamWs","PctGlacLakeCrsWs","PctCoastCrsCat","PctGlacTilClayWs","PctHydricWs","PctCoastCrsWs",
#               "MineDensCat","MineDensCatRp100", "DamDensCat","DamNrmStorM3Cat",
#               "PctCarbResidCat","PctAlkIntruVolCat","PctColluvSedCat","PctEolCrsCat","PctEolCrsWs", "PctEolFineCat","PctSalLakeCat","PctWaterCat","PctAlkIntruVolWs",
#               "PctColluvSedWs","PctEolFineWs","PctGlacLakeFineWs") #Includes nearly-constant vars, other rejects
# full.vars<-setdiff(c(core.candidates.mines.dams.atm.veg, nat.cands), const.vars)
# setdiff(full.vars, names(csci))
# 
# csci.rf.dat<-csci[which(csci$SiteSet=="Cal" & csci$SelectedSample=="Selected"),]
# 
# # junk<-melt(csci.rf.dat[,core.candidates])
# # nrow(ddply(junk, .(variable), summarize, n_unique=length(unique(value))))
# # 
# # ggplot(data=junk[which(junk$variable %in% unique(junk$variable)[1:9]),],aes(x=value))+
# #   geom_histogram()+
# #   facet_wrap(~variable, scales="free")+scale_y_sqrt()
# 
# # summary(csci.rf.dat[,full.vars])
# 
# 
# 
# 
# library(quantregForest)
# 
# rf_full.dat<-na.omit(csci.rf.dat[,c("SampleID", "COMID", "CSCI",full.vars)])
# 
# 
# 
# 
# set.seed(10101)
# rf_full<-quantregForest(y=rf_full.dat$CSCI,
#                         x=as.matrix(rf_full.dat[,c(full.vars)]),
#                         keep.inbag=T, importance=T,proximity=T)
# 
# set.seed(10012)
# rf_core<-quantregForest(y=csci.rf.dat$CSCI,
#                         x=as.matrix(csci.rf.dat[,core.candidates]),
#                         keep.inbag=T, importance=T,proximity=T)
# 
# 
# revisited.sites<-  unique(csci$StationCode[which(csci$SelectedSample=="NotSelected")])
# csci$Revisited<-csci$StationCode %in% revisited.sites
# 
# ######
# # get model predictions, have to separate calibration oob from statewide
#   
# ##
# # core model
# 
# # prediction data w/o calibration dataset
# newdatcr <- all.comid %>% 
#   filter(!COMID %in% csci.rf.dat$COMID) %>% 
#   select_(.dots = c('COMID', core.candidates)) %>% 
#   na.omit
# 
# # out of bag predictions for calibration dataset
# # estimates for same comid must be averaged with oob estimates
# predcore_oob <- predict(rf_core, what=seq(from=0.05, to=.95, by=.05), na.rm=T) %>% 
#   as.data.frame %>% 
#   mutate(COMID = csci.rf.dat$COMID) %>% 
#   gather('var', 'val', -COMID) %>% 
#   group_by(COMID, var) %>% 
#   summarize(val = mean(val)) %>% 
#   spread(var, val) %>% 
#   .[, c(2:20, 1)]
# predcore_all <- predict(rf_core, newdata = newdatcr[, -1], what=seq(from=0.05, to=.95, by=.05), na.rm=T) %>% 
#   as.data.frame %>% 
#   mutate(COMID = newdatcr$COMID)
# 
# # join calibration oob with statewide
# predcore <- bind_rows(predcore_oob, predcore_all)
# names(predcore) <- c(paste0("core",formatC(as.numeric(seq(from=0.05, to=.95, by=.05)), format = 'f', flag='0', digits = 2)), 'COMID')
# 
# ##
# # full model
# 
# # prediction data w/o calibration dataset
# newdatfl <- all.comid %>% 
#   select_(.dots = c('COMID', core.candidates.mines.dams.atm.veg, nat.cands)) %>% 
#   na.omit %>% 
#   filter(!COMID %in% rf_full.dat$COMID) 
# 
# # out of bag predictions for calibration dataset
# predfull_oob <- predict(rf_full, what=seq(from=0.05, to=.95, by=.05), na.rm=T) %>% 
#   as.data.frame %>% 
#   mutate(COMID = rf_full.dat$COMID) %>% 
#   gather('var', 'val', -COMID) %>% 
#   group_by(COMID, var) %>% 
#   summarize(val = mean(val)) %>% 
#   spread(var, val) %>% 
#   .[, c(2:20, 1)]
# predfull_all <- predict(rf_core, newdata = newdatfl[, -1], what=seq(from=0.05, to=.95, by=.05), na.rm=T) %>% 
#   as.data.frame %>% 
#   mutate(COMID = newdatfl$COMID)
# 
# # join calibration oob with statewide
# predfull <- bind_rows(predfull_oob, predfull_all)
# names(predfull) <- c(paste0("full",formatC(as.numeric(seq(from=0.05, to=.95, by=.05)), format = 'f', flag='0', digits = 2)), 'COMID')
# 
# pred_all <- predcore %>% 
#   left_join(predfull, by = 'COMID') %>% 
#   left_join(all.comid[,c("COMID", core.candidates, setdiff(full.vars,core.candidates))], by = 'COMID') %>% 
#   as.data.frame
#               
# pred_all$DevData<-
#   ifelse(pred_all$COMID %in% csci$COMID[which(csci$SiteSet=="Cal")],"Cal",
#          ifelse(pred_all$COMID %in% csci$COMID[which(csci$SiteSet=="Val")],"Val","No"))
# comid_prd <- pred_all
# 
# # csci data for comparison with stream comid
# csci_comid <- csci
# 
# ## 
# # save all
# 
# save(csci_comid, file = 'C:/proj/manuscripts/landscape_mod/data/csci_comid.RData', compress = 'xz')
# save(comid_prd, file = 'C:/proj/manuscripts/landscape_mod/data/comid_prd.RData', compress = 'xz')
# save(rf_core,file = "C:/proj/manuscripts/landscape_mod/data/rf_core.Rdata", compress = 'xz')
# save(rf_full,file = "C:/proj/manuscripts/landscape_mod/data/rf_full.Rdata", compress = 'xz')
