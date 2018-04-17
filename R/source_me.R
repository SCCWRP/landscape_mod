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

######
# sensitity analysis, statewide

data(comid_prd)
data(csci_comid)
data(strclslen)

source('R/funcs.R')
prj <- geo_wgs84

# calipsa
calipsa <- readOGR('S:/Spatial_Data/RCMP_needs editting/Inputs/PSA6_090111/PSA6_2011.shp') %>% 
  spTransform(prj) %>% 
  st_as_sf

# comid predictions to join with calinhd
comid_prd <- comid_prd %>% 
  dplyr::select(COMID, matches('^full|^COMID$'))

# all cali hydro, really simplified
calinhd <- readOGR('S:/Spatial_Data/NHDPlus/NHDPLusCalifornia/NHDPlusCalifornia.shp') %>% 
  spTransform(prj) %>% 
  st_as_sf %>% 
  st_simplify(dTolerance = 0.5, preserveTopology = T) %>% 
  filter(COMID %in% unique(comid_prd$COMID)) %>% 
  dplyr::select(COMID) %>%
  inner_join(comid_prd, by = 'COMID')

# format csci data for site_exp function, takes averages of repeats at each site
csci_agg <- csci_comid %>% 
  dplyr::select(COMID, StationCode, CSCI, SampleDate, FieldReplicate, New_Lat, New_Long) %>% 
  group_by(COMID, StationCode) %>% 
  summarise(
    csci = mean(CSCI, na.rm = T), 
    lat = mean(New_Lat, na.rm = T),
    long = mean(New_Long, na.rm = T)
  ) %>% 
  ungroup

# psa by csci comid
cscipsa <- csci_comid %>% 
  dplyr::select(COMID, PSA6c) %>% 
  unique %>% 
  mutate(COMID = as.character(COMID))

# search grid
grd_chk <- list(
  thrsh = c(0.63, 0.79, 0.92), 
  tails = c(0.05, 0.25, 0.45)
) %>% 
  cross_df

# setup parallel backend
ncores <- detectCores() - 2  
cl<-makeCluster(ncores)
registerDoParallel(cl)
strt<-Sys.time()

sensres <- foreach(i = 1:nrow(grd_chk), 
                   .export = c('calipsa', 'comid_prd', 'calinhd', 'csci_agg', 'cscipsa'), 
                   .packages = c('sp', 'mapview', 'leaflet', 'tidyverse', 'raster', 'maptools', 'rgdal', 'sf')) %dopar% {
                     
                     sink('log.txt')
                     cat(i, 'of', nrow(grd_chk), '\n')
                     print(Sys.time()-strt)
                     sink()
                     
                     source("R/funcs.R")
                     
                     thrshi <- grd_chk[[i, 'thrsh']]
                     tailsi <- grd_chk[[i, 'tails']]
                     
                     # get biological condition expectations
                     cls <- getcls2(calinhd, thrsh = thrshi, tails = tailsi, modls = 'full')
                     
                     calicls <- calinhd %>% 
                       left_join(cls, by = 'COMID')
                     
                     ######
                     # all cali CSCI site expectations
                     
                     # get site expectations
                     caliexp <- site_exp(calinhd, csci_agg, thrsh = thrshi, tails = tailsi, modls = 'full')
                     
                     ######
                     # get stream length in each stream class by PSA
                     
                     strclsleni <- calicls %>% 
                       dplyr::select(COMID, strcls) %>% 
                       st_intersection(calipsa)
                     
                     lens <- strclsleni %>% 
                       as('Spatial') %>% 
                       sp::SpatialLinesLengths(longlat = T)
                     
                     strclsleni <- strclsleni %>% mutate(lens = lens) 
                     
                     st_geometry(strclsleni) <- NULL
                     lentot <- strclsleni %>% 
                       filter(!is.na(strcls)) %>% 
                       mutate(
                         Region = factor(PSA6, 
                                         levels = c('Central Valley', 'Chaparral', 'Deserts Modoc', 'North Coast', 'Sierra Nevada', 'South Coast'),
                                         labels = c('CV', 'CH', 'DM', 'NC', 'SN', 'SC')
                         )
                       ) %>% 
                       group_by(Region, strcls) %>% 
                       summarise(len = sum(lens)) %>% 
                       spread(strcls, len, fill = 0) %>% 
                       gather('strcls', 'len', -Region) %>% 
                       group_by(Region) %>% 
                       mutate(
                         perc = 100 * len / sum(len), 
                         strcls = factor(strcls, levels = c('likely constrained', 'possibly constrained', 'possibly unconstrained', 'likely unconstrained'))
                       )
                     
                     lentotall <- strclsleni %>% 
                       filter(!is.na(strcls)) %>% 
                       group_by(strcls) %>% 
                       summarise(len = sum(lens)) %>%
                       ungroup %>% 
                       mutate(
                         Region = 'statewide',
                         perc = 100 * len / sum(len), 
                         strcls = factor(strcls, levels = c('likely constrained', 'possibly constrained', 'possibly unconstrained', 'likely unconstrained'))
                       ) 
                     
                     totab1 <- bind_rows(lentotall, lentot)
                     
                     # format table of site counts by region
                     reltot <- caliexp %>% 
                       dplyr::select(COMID, csci, StationCode, strcls, perf) %>% 
                       filter(!is.na(strcls)) %>% 
                       left_join(cscipsa, by = 'COMID') %>% 
                       rename(Region = PSA6c) %>% 
                       group_by(Region, perf) %>% 
                       summarise(
                         n = n(), 
                         csciave = mean(csci, na.rm = T), 
                         cscisd = sd(csci, na.rm = T)
                       ) %>% 
                       group_by(Region) %>% 
                       mutate(
                         perc = n / sum(n)
                       ) %>% 
                       complete(Region, perf)
                     
                     # whole state
                     reltotall <- caliexp %>% 
                       dplyr::select(COMID, csci, StationCode, strcls, perf) %>% 
                       filter(!is.na(strcls)) %>% 
                       left_join(cscipsa, by = 'COMID') %>% 
                       rename(Region = PSA6c) %>% 
                       group_by(perf) %>% 
                       summarise(
                         n = n(), 
                         csciave = mean(csci, na.rm = T), 
                         cscisd = sd(csci, na.rm = T)
                       ) %>% 
                       ungroup %>% 
                       mutate(
                         perc = n / sum(n)
                       ) %>% 
                       mutate(Region = 'Statewide') 
                     
                     totab2 <- bind_rows(reltotall, reltot)
                     
                     out <- list(totab1, totab2)
                     return(out)
                     
                   }

save(sensres, file = 'data/sensres.RData', compress = 'xz')
