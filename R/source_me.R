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

source('R/funcs.R')
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
                    labels = c('CV', 'CH', 'DM', 'NC', 'SN', 'SC')
    )
  )

save(psalab, file = 'data/psalab.RData', compress = 'xz')

##
# cali NHD simplify, fortify, and save

load(file = 'data/comid_prd.RData')

# all cali hydro, really simplified
calinhd <- readOGR('S:/Spatial_Data/NHDPlus/All_Ca_NHDPlus.shp') %>% 
  spTransform(prj) %>% 
  st_as_sf %>% 
  st_simplify(dTolerance = 0.5, preserveTopology = T)

# comid predictions to join with calinhd
comid_prd <- comid_prd %>% 
  dplyr::select(COMID, full0.50)

# fortified calinhd, joind with comid pred
nhdplo <- calinhd %>% 
  filter(COMID %in% unique(comid_prd$COMID)) %>% 
  dplyr::select(COMID) %>% 
  as('Spatial')
comidid <- nhdplo@data %>% 
  dplyr::select(COMID) %>% 
  rownames_to_column('id')
nhdplo <- nhdplo %>% 
  fortify %>% 
  left_join(comidid, by = 'id') %>%  
  inner_join(comid_prd, by = 'COMID')

save(nhdplo, file = 'data/nhdplo.RData', compress = 'xz')

######
# cali nhd stream classes, simplified, fortified, and saved

data(comid_prd)

# comid predictions to join with calinhd
comid_prd <- comid_prd %>% 
  dplyr::select(COMID, matches('^full|^COMID$'))

# all cali hydro, really simplified
calinhd <- readOGR('S:/Spatial_Data/NHDPlus/All_Ca_NHDPlus.shp') %>% 
  spTransform(prj) %>% 
  st_as_sf %>% 
  st_simplify(dTolerance = 0.5, preserveTopology = T) %>% 
  filter(COMID %in% unique(comid_prd$COMID)) %>% 
  dplyr::select(COMID) %>%
  inner_join(comid_prd, by = 'COMID')

# get biological condition expectations
cls <- getcls2(calinhd, thrsh = 0.79, tails = 0.05, modls = 'full')

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
caliexp <- site_exp(calinhd, csci_comid, thrsh = 0.79, tails = 0.05, modls = 'full')

save(caliexp, file = 'data/caliexp.RData', compress = 'xz')

######
# get stream length in each stream class by PSA

load(file = 'data/calicls.RData')
load(file = 'data/caliexp.RData')
load(file = 'data/csci_comid.RData')

strclslen <- calicls %>% 
  dplyr::select(COMID, strcls) %>% 
  st_intersection(calipsa)

lens <- strclslen %>% 
  as('Spatial') %>% 
  sp::SpatialLinesLengths(longlat = T)

strclslen <- strclslen %>% mutate(lens = lens) 

save(strclslen, file = 'data/strclslen.RData', compress = 'xz')

######
# expected range of scores for urban, ag, other by region
data(calicls)
data(strclslen)

# streamcat data
strmcat <- rbind(read.csv("Z:/MarcusBeck/Landscape models from rafi/Streamcat_v2_AllCOMID_030117/exp_1.csv", stringsAsFactors = F),
                 read.csv("Z:/MarcusBeck/Landscape models from rafi/Streamcat_v2_AllCOMID_030117/exp_2.csv", stringsAsFactors = F),
                 read.csv("Z:/MarcusBeck/Landscape models from rafi/Streamcat_v2_AllCOMID_030117/exp_3.csv", stringsAsFactors = F))
comid.nats <- read.csv("Z:/MarcusBeck/Landscape models from rafi/ALL_COMID_Nats.csv", stringsAsFactors = F)
strmcat <- plyr::join(strmcat, comid.nats[,setdiff(names(comid.nats), "WsAreaSqKm")])

strmcat$TotUrb2011Ws<-  rowSums(strmcat[,c("PctUrbOp2011Ws","PctUrbLo2011Ws","PctUrbMd2011Ws","PctUrbHi2011Ws")])
strmcat$TotUrb2011Cat<-  rowSums(strmcat[,c("PctUrbOp2011Cat","PctUrbLo2011Cat","PctUrbMd2011Cat","PctUrbHi2011Cat")])
strmcat$TotUrb2011WsRp100<-  rowSums(strmcat[,c("PctUrbOp2011WsRp100","PctUrbLo2011WsRp100","PctUrbMd2011WsRp100","PctUrbHi2011WsRp100")])
strmcat$TotUrb2011CatRp100<-  rowSums(strmcat[,c("PctUrbOp2011CatRp100","PctUrbLo2011CatRp100","PctUrbMd2011CatRp100","PctUrbHi2011CatRp100")])

strmcat$TotAg2011Ws<-  rowSums(strmcat[,c("PctHay2011Ws","PctCrop2011Ws")])
strmcat$TotAg2011Cat<-  rowSums(strmcat[,c("PctHay2011Cat","PctCrop2011Cat")])
strmcat$TotAg2011WsRp100<-  rowSums(strmcat[,c("PctHay2011WsRp100","PctCrop2011WsRp100")])
strmcat$TotAg2011CatRp100<-  rowSums(strmcat[,c("PctHay2011CatRp100","PctCrop2011CatRp100")])

strmcat <- strmcat %>% 
  dplyr::select(COMID, TotUrb2011Ws, TotAg2011Ws)

# PSA by all comid
psaall <- strclslen %>% 
  dplyr::select(COMID, PSA6)
st_geometry(psaall) <- NULL

# data to summarize
scrdist <- calicls
st_geometry(scrdist) <- NULL
scrdist <- scrdist %>% 
  dplyr::select(-strcls, -strcls_int) %>% 
  left_join(strmcat, by = 'COMID') %>% 
  left_join(psaall, by = 'COMID') %>% 
  filter(!is.na(PSA6)) %>% 
  rename(Region = PSA6)

# get range of scores for average urban, ag, open locations, by region
scrdistreg <- scrdist %>% 
  group_by(Region) %>% 
  nest %>% 
  mutate(grps = purrr::map(data, function(x){
    
    # get kmeans, centers
    tomod <- x %>% 
      mutate(
        Urb = log10(1 + TotUrb2011Ws), 
        Ag = log10(1 + TotAg2011Ws)
      ) %>% 
      dplyr::select(Urb, Ag) 
    ngrps <- 8
    kmod <- kmeans(tomod, centers = ngrps, nstart = 10, iter.max = 100, algorithm = 'MacQueen')
    grps <- kmod$cluster
    cent <- kmod$centers %>% 
      data.frame
    
    # find the typical ag, urb group
    urb <- which.max(cent[, 'Urb'])
    ag <- which.max(cent[, 'Ag'])
    oth <- which.min(rowSums(cent))
    
    grplabs <- list(
      urb = urb,
      ag = ag,
      oth = oth
    )
    
    ests <- x %>% 
      mutate(
        grps = grps
      ) %>% 
      filter(grps %in% unlist(grplabs)) %>% 
      mutate(
        grps = factor(grps, levels = unlist(grplabs), labels = names(grplabs))
      ) %>% 
      group_by(grps) %>% 
      summarise(
        lo05 = mean(full0.05, na.rm = T), 
        hi95 = mean(full0.95, na.rm = T),
        lo25 = mean(full0.25, na.rm = T), 
        hi75 = mean(full0.75, na.rm = T), 
        lo45 = mean(full0.45, na.rm = T),
        hi55 = mean(full0.55, na.rm = T)
      )
    
    return(ests)
    
  })) %>% 
  dplyr::select(-data) %>% 
  unnest %>% 
  mutate(
    Region = factor(Region, 
                    levels = c('Central Valley', 'Chaparral', 'Deserts Modoc', 'North Coast', 'Sierra Nevada', 'South Coast'),
                    labels = c('CV', 'CH', 'DM', 'NC', 'SN', 'SC')),
    Region = as.character(Region), 
    grps = as.character(grps)
  )

# typical statewide
scrdistall <- scrdist %>% 
  mutate(Region = 'Statewide') %>% 
  group_by(Region) %>% 
  nest %>% 
  mutate(grps = purrr::map(data, function(x){
    
    # get kmeans, centers
    tomod <- x %>% 
      mutate(
        Urb = log10(1 + TotUrb2011Ws), 
        Ag = log10(1 + TotAg2011Ws)
      ) %>% 
      dplyr::select(Urb, Ag) 
    ngrps <- 8
    kmod <- kmeans(tomod, centers = ngrps, nstart = 10, iter.max = 100, algorithm="MacQueen")
    grps <- kmod$cluster
    cent <- kmod$centers %>% 
      data.frame
    
    # find the typical ag, urb group
    urb <- which.max(cent[, 'Urb'])
    ag <- which.max(cent[, 'Ag'])
    oth <- which.min(rowSums(cent))
    
    grplabs <- list(
      urb = urb,
      ag = ag,
      oth = oth
    )
    
    ests <- x %>% 
      mutate(
        grps = grps
      ) %>% 
      filter(grps %in% unlist(grplabs)) %>% 
      mutate(
        grps = factor(grps, levels = unlist(grplabs), labels = names(grplabs))
      ) %>% 
      group_by(grps) %>% 
      summarise(
        lo05 = mean(full0.05, na.rm = T), 
        hi95 = mean(full0.95, na.rm = T),
        lo25 = mean(full0.25, na.rm = T), 
        hi75 = mean(full0.75, na.rm = T), 
        lo45 = mean(full0.45, na.rm = T),
        hi55 = mean(full0.55, na.rm = T)
      )
    
    return(ests)
    
  })) %>% 
  dplyr::select(-data) %>% 
  unnest %>% 
  mutate(
    Region = as.character(Region),
    grps = as.character(grps)
  )


typscrs <- bind_rows(scrdistall, scrdistreg)

save(typscrs, file = 'data/typscrs.RData', compress = 'xz')

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
ncores <- detectCores() - 1  
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
