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
      labels = c('CV', 'Ch', 'DM', 'NC', 'SN', 'SC')
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

######
# create data from Rafi
# create csci comid data, rf core and full models for csci quantiles, quantile preds for all comid in CA
# original file from Z:/MarcusBeck/Landscape models from rafi/modlu_120117.R

setwd("Z:/MarcusBeck/Landscape models from rafi/")
library(plyr)
library(dplyr)
library(reshape2)

# csci<-read.csv("CSCI_LU_temp_022817.csv", stringsAsFactors = F) #Load data
csci<-join(read.csv("CSCI_LU_temp_040617.csv", stringsAsFactors = F), #Load data
           read.csv("CSCI_Nat_temp_061217.csv", stringsAsFactors = F))
nat.vars.full<-setdiff(names(read.csv("CSCI_Nat_temp_061217.csv")), names(read.csv("CSCI_LU_temp_040617.csv")))
nat.vars.full<-unique(c("WsAreaSqKm",nat.vars.full))

csci$SampleID<-paste(csci$StationCode, csci$SampleDate, csci$CollectionMethodCode, csci$FieldReplicate, sep="_")
csci<-csci[which(!is.na(csci$PctFrstLossWs)),] #Drop rows with missing STREAMCAT data


csci$RdDensCatRp100[is.na(csci$RdDensCatRp100)]<-0

csci$TotUrb2011Ws<-  rowSums(csci[,c("PctUrbOp2011Ws","PctUrbLo2011Ws","PctUrbMd2011Ws","PctUrbHi2011Ws")])
csci$TotUrb2011Cat<-  rowSums(csci[,c("PctUrbOp2011Cat","PctUrbLo2011Cat","PctUrbMd2011Cat","PctUrbHi2011Cat")])
csci$TotUrb2011WsRp100<-  rowSums(csci[,c("PctUrbOp2011WsRp100","PctUrbLo2011WsRp100","PctUrbMd2011WsRp100","PctUrbHi2011WsRp100")])
csci$TotUrb2011CatRp100<-  rowSums(csci[,c("PctUrbOp2011CatRp100","PctUrbLo2011CatRp100","PctUrbMd2011CatRp100","PctUrbHi2011CatRp100")])

csci$TotAg2011Ws<-  rowSums(csci[,c("PctHay2011Ws","PctCrop2011Ws")])
csci$TotAg2011Cat<-  rowSums(csci[,c("PctHay2011Cat","PctCrop2011Cat")])
csci$TotAg2011WsRp100<-  rowSums(csci[,c("PctHay2011WsRp100","PctCrop2011WsRp100")])
csci$TotAg2011CatRp100<-  rowSums(csci[,c("PctHay2011CatRp100","PctCrop2011CatRp100")])

gis.area<-read.csv("csci_AREA.csv", stringsAsFactors = F)

setdiff(csci$StationCode, gis.area$StationCode)
#Let's look at PCA-space:



#COMIDs
all.comid<-rbind(read.csv("Streamcat_v2_AllCOMID_030117/exp_1.csv", stringsAsFactors = F),
                 read.csv("Streamcat_v2_AllCOMID_030117/exp_2.csv", stringsAsFactors = F),
                 read.csv("Streamcat_v2_AllCOMID_030117/exp_3.csv", stringsAsFactors = F))
comid.nats<-read.csv("ALL_COMID_Nats.csv", stringsAsFactors = F)
all.comid<-join(all.comid, comid.nats[,setdiff(names(comid.nats), "WsAreaSqKm")])

all.comid$TotUrb2011Ws<-  rowSums(all.comid[,c("PctUrbOp2011Ws","PctUrbLo2011Ws","PctUrbMd2011Ws","PctUrbHi2011Ws")])
all.comid$TotUrb2011Cat<-  rowSums(all.comid[,c("PctUrbOp2011Cat","PctUrbLo2011Cat","PctUrbMd2011Cat","PctUrbHi2011Cat")])
all.comid$TotUrb2011WsRp100<-  rowSums(all.comid[,c("PctUrbOp2011WsRp100","PctUrbLo2011WsRp100","PctUrbMd2011WsRp100","PctUrbHi2011WsRp100")])
all.comid$TotUrb2011CatRp100<-  rowSums(all.comid[,c("PctUrbOp2011CatRp100","PctUrbLo2011CatRp100","PctUrbMd2011CatRp100","PctUrbHi2011CatRp100")])

all.comid$TotAg2011Ws<-  rowSums(all.comid[,c("PctHay2011Ws","PctCrop2011Ws")])
all.comid$TotAg2011Cat<-  rowSums(all.comid[,c("PctHay2011Cat","PctCrop2011Cat")])
all.comid$TotAg2011WsRp100<-  rowSums(all.comid[,c("PctHay2011WsRp100","PctCrop2011WsRp100")])
all.comid$TotAg2011CatRp100<-  rowSums(all.comid[,c("PctHay2011CatRp100","PctCrop2011CatRp100")])

all.comid$RdDensWsRp100<-all.comid$RoadDensityRipBuf100_CA.RdDensWsRp100
all.comid$RdDensCatRp100<-all.comid$RoadDensityRipBuf100_CA.RdDensCatRp100

library(plyr)
csci<-join(csci, gis.area)



library(dplyr)
# csci2<-csci %>% group_by(COMID) %>% slice(which.max(AREA_SQKM))

csci$Largest<-
  sapply(1:nrow(csci), function(i){
    comid<-csci$COMID[i]
    area<-csci$AREA_SQKM[i]
    mydf<-csci[which(csci$COMID==comid),]
    area==max(mydf$AREA_SQKM)
  })

# csci[which(csci$COMID==8264738),c("StationCode","COMID","AREA_SQKM","Largest")]
csci.full<-csci
csci<-csci.full[csci.full$Largest,]

library(quantregForest)



#Create cal and val data
sites<-unique(csci[,c("StationCode","PSA6c","PctImp2006Ws")]) #Create a sites DF
library(plyr)
sites.d<-ddply(sites, .(PSA6c), summarize, ImpT1=quantile(PctImp2006Ws, probs=c(0.25)),ImpT2=quantile(PctImp2006Ws, probs=c(0.5)),ImpT3=quantile(PctImp2006Ws, probs=c(0.75)))
write.table(sites.d, "clipboard", sep="\t")

sites$RegImpQ<-   #Divide regions into thirds based on imperviousness
  sapply(1:nrow(sites), function(i){
    region<-sites$PSA6c[i]
    imp<-sites$PctImp2006Ws[i]
    imp.t1<-sites.d$ImpT1[which(sites.d$PSA6c==region)]
    imp.t2<-sites.d$ImpT2[which(sites.d$PSA6c==region)]
    imp.t3<-sites.d$ImpT3[which(sites.d$PSA6c==region)]
    ifelse(imp<imp.t1,"T1", ifelse(imp<imp.t2,"T2", ifelse(imp<imp.t3,"T3","T4")))
  })

table(sites$PSA6c, sites$RegImpQ)
sites$Stratum<-paste(sites$PSA6c, sites$RegImpQ, sep="_") #Create strata by dividing each region into thirds
table(sites$Stratum)

#Random assignment into cal (80%) and val (20%) data sets by Stratum field
library(dplyr)
sites.t<- sites %>% group_by(Stratum)
set.seed(500)
sites.cal<-sample_frac(sites.t, 0.75, replace=F)
sites$SiteSet<-ifelse(sites$StationCode %in% sites.cal$StationCode, "Cal","Val")
csci<-join(csci, sites[,c("StationCode","RegImpQ", "Stratum", "SiteSet")])

#
#Select a single sample for assessment
csci.t<- csci %>% group_by(StationCode)
set.seed(501)
samps.sel<-sample_n(csci.t, 1)
csci$SelectedSample<-ifelse(csci$SampleID %in% samps.sel$SampleID, "Selected","NotSelected")

#PotentialVars
metrics<-read.csv("metrics.csv", stringsAsFactors = F)


core.candidates<-c("CanalDensCat","CanalDensWs", 
                   "PctImp2006Cat","PctImp2006Ws","PctImp2006CatRp100","PctImp2006WsRp100",
                   "TotUrb2011Ws","TotUrb2011Cat","TotUrb2011WsRp100","TotUrb2011CatRp100",
                   "TotAg2011Ws","TotAg2011Cat","TotAg2011WsRp100","TotAg2011CatRp100",
                   "RdDensCat", "RdDensWs", "RdDensCatRp100", "RdDensWsRp100", "RdCrsCat","RdCrsWs")
core.candidates.disag<-c("CanalDensCat","CanalDensWs",
                         "PctImp2006Cat","PctImp2006Ws","PctImp2006CatRp100","PctImp2006WsRp100",
                         "PctUrbOp2011Cat", "PctUrbLo2011Cat",
                         "PctUrbMd2011Cat", "PctUrbHi2011Cat", "PctHay2011Cat", "PctCrop2011Cat",
                         "PctUrbOp2011Ws", "PctUrbLo2011Ws", "PctUrbMd2011Ws", "PctUrbHi2011Ws",
                         "PctHay2011Ws", "PctCrop2011Ws", "PctUrbOp2011CatRp100", "PctUrbLo2011CatRp100",
                         "PctUrbMd2011CatRp100", "PctUrbHi2011CatRp100", "PctHay2011CatRp100",
                         "PctCrop2011CatRp100", "PctUrbOp2011WsRp100", "PctUrbLo2011WsRp100",
                         "PctUrbMd2011WsRp100", "PctUrbHi2011WsRp100", "PctHay2011WsRp100",
                         "PctCrop2011WsRp100",
                         "RdDensCat", "RdDensWs", "RdDensCatRp100", "RdDensWsRp100", "RdCrsCat", "RdCrsSlpWtdCat","RdCrsWs","RdCrsSlpWtdWs")

core.candidates.mines.dams<-c(core.candidates,
                              "MineDensCat", "MineDensWs", "MineDensCatRp100", "MineDensWsRp100",
                              "DamDensCat", "DamDensWs",  "DamNrmStorM3Cat",  "DamNrmStorM3Ws")

core.candidates.disag.mines.dams<-c(core.candidates.disag,
                                    "MineDensCat", "MineDensWs", "MineDensCatRp100", "MineDensWsRp100",
                                    "DamDensCat", "DamDensWs",  "DamNrmStorM3Cat",  "DamNrmStorM3Ws")

stressors.atm<-c(#"NH4_2008Cat","NO3_2008Cat","InorgNWetDep_2008Cat","SN_2008Cat", #Cats have a small number of NAs
  "NH4_2008Ws","NO3_2008Ws","InorgNWetDep_2008Ws","SN_2008Ws")

stressors.other<-c("PctNonAgIntrodManagVegCat","PctNonAgIntrodManagVegWs","PctNonAgIntrodManagVegCatRp100","PctNonAgIntrodManagVegWsRp100")

core.candidates.mines.dams.atm<-c(core.candidates.mines.dams, stressors.atm)
core.candidates.disag.mines.dams.atm<-c(core.candidates.disag.mines.dams, stressors.atm)

core.candidates.mines.dams.atm.veg<-c(core.candidates.mines.dams.atm,stressors.other)
nat.cands<-setdiff(nat.vars.full, c("MAST_2008","MAST_2009","MAST_2013","MAST_2014",
                                    "MSST_2008","MSST_2009","MSST_2013","MSST_2014",
                                    "MWST_2008","MWST_2009","MWST_2013","MWST_2014"))
const.vars<-c("PctGlacTilClayCat","PctHydricCat","PctGlacTilLoamCat","PctGlacLakeFineCat","PctGlacLakeCrsCat","PctGlacTilLoamWs","PctGlacLakeCrsWs","PctCoastCrsCat","PctGlacTilClayWs","PctHydricWs","PctCoastCrsWs",
              "MineDensCat","MineDensCatRp100", "DamDensCat","DamNrmStorM3Cat",
              "PctCarbResidCat","PctAlkIntruVolCat","PctColluvSedCat","PctEolCrsCat","PctEolCrsWs", "PctEolFineCat","PctSalLakeCat","PctWaterCat","PctAlkIntruVolWs",
              "PctColluvSedWs","PctEolFineWs","PctGlacLakeFineWs") #Includes nearly-constant vars, other rejects
full.vars<-setdiff(c(core.candidates.mines.dams.atm.veg, nat.cands), const.vars)
setdiff(full.vars, names(csci))

csci.rf.dat<-csci[which(csci$SiteSet=="Cal" & csci$SelectedSample=="Selected"),]

# junk<-melt(csci.rf.dat[,core.candidates])
# nrow(ddply(junk, .(variable), summarize, n_unique=length(unique(value))))
# 
# ggplot(data=junk[which(junk$variable %in% unique(junk$variable)[1:9]),],aes(x=value))+
#   geom_histogram()+
#   facet_wrap(~variable, scales="free")+scale_y_sqrt()

# summary(csci.rf.dat[,full.vars])




library(quantregForest)

rf_full.dat<-na.omit(csci.rf.dat[,c("SampleID", "CSCI",full.vars)])




set.seed(10101)
rf_full<-quantregForest(y=rf_full.dat$CSCI,
                        x=as.matrix(rf_full.dat[,c(full.vars)]),
                        keep.inbag=T, importance=T,proximity=T)

set.seed(10012)
rf_core<-quantregForest(y=csci.rf.dat$CSCI,
                        x=as.matrix(csci.rf.dat[,core.candidates]),
                        keep.inbag=T, importance=T,proximity=T)


revisited.sites<-  unique(csci$StationCode[which(csci$SelectedSample=="NotSelected")])

csci$Revisited<-csci$StationCode %in% revisited.sites



# all.comid2<-na.omit(all.comid[,stressors.hiq])

start.time <- Sys.time()
comid.pred.core<-predict(rf_core, newdata=as.matrix(na.omit(all.comid[, core.candidates])), what=seq(from=0.05, to=.95, by=.05), na.rm=T)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

comid.pred.core.df<-as.data.frame(comid.pred.core)
names(comid.pred.core.df)<-paste0("core",formatC(as.numeric(seq(from=0.05, to=.95, by=.05)), format = 'f', flag='0', digits = 2))
comid.pred.core.df$COMID<-  na.omit(all.comid[,c("COMID",core.candidates)])$COMID
head(comid.pred.core.df)


# junk<-((all.comid[, c(core.candidates.mines.dams.atm.veg, nat.cands)]))
# head(junk)

comid.pred.full<-predict(rf_full, newdata=as.matrix(na.omit(all.comid[, c(core.candidates.mines.dams.atm.veg, nat.cands)])), what=seq(from=0.05, to=.95, by=.05), na.rm=T)
comid.pred.full.df<-as.data.frame(comid.pred.full)
names(comid.pred.full.df)<-paste0("full",formatC(as.numeric(seq(from=0.05, to=.95, by=.05)), format = 'f', flag='0', digits = 2))
comid.pred.full.df$COMID<-  na.omit(all.comid[,c("COMID",core.candidates.mines.dams.atm.veg, nat.cands)])$COMID
head(comid.pred.full.df)

comid.pred.modeled<-join(comid.pred.core.df[,c(20,1:19)], comid.pred.full.df)
summary(comid.pred.modeled)

comid.pred.modeled<-join(comid.pred.modeled, all.comid[,c("COMID", core.candidates, setdiff(full.vars,core.candidates))])

comid.pred.modeled$DevData<-
  ifelse(comid.pred.modeled$COMID %in% csci$COMID[which(csci$SiteSet=="Cal")],"Cal",
         ifelse(comid.pred.modeled$COMID %in% csci$COMID[which(csci$SiteSet=="Val")],"Val","No"))
comid_prd <- comid.pred.modeled

# csci data for comparison with stream comid
csci_comid <- csci

## 
# save all

save(csci_comid, file = 'C:/proj/manuscripts/landscape_mod/ata/csci_comid.RData', compress = 'xz')
save(comid_prd, file = 'C:/proj/manuscripts/landscape_mod/data/comid_prd.RData', compress = 'xz')
save(rf_core,file = "C:/proj/manuscripts/landscape_mod/data/rf_core.Rdata", compress = 'xz')
save(rf_full,file = "C:/proj/manuscripts/landscape_mod/data/rf_full.Rdata", compress = 'xz')
