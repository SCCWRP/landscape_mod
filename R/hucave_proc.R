# get HUC12 median CSCI predictions from landscape model

# create the data

library(tidyverse)
library(sf)

# from https://knb.ecoinformatics.org/view/urn:uuid:75411f50-32ed-42a5-bbfd-26833c7a441f
# created in R/metadata_proc.R
prd <- st_read('Z:/MarcusBeck/GIS/strm_constraints.shp')

# HUC polygons, https://ftp.wildlife.ca.gov/?ShareToken=E1ED000FE07BE1FD24B40CDAEA8C454569499F8C
# see email for password
huc <- st_read(dsn = 'Z:/MarcusBeck/GIS/HUCs/baseLayers.gdb') %>% 
  st_transform(st_crs(prd))

# intersect comids with predictions  
int <- st_intersection(prd, huc)

# get COMID lengths
lens <- st_length(int)

# remove geometry
hucave <- int
st_geometry(hucave) <- NULL

# add lengths
hucave <- hucave %>% 
  select(COMID, qt50, Ref10, HUC12) %>% 
  mutate(
    lens = lens,
    qt50 = ifelse(qt50 == -999999, NA, qt50)
    )

# summarize by HUC12 as weighted mean, reach count, and reach count as likely unconstrained
hucave <- hucave %>% 
  group_by(HUC12) %>%
  summarise(
    qt50 = weighted.mean(qt50, w = lens, na.rm = T),
    n_reach = n(),
    n_lu = sum(Ref10 %in% 'likely unconstrained')
  )

# combine aggregates with polygon shapefile
huc <- huc %>% 
  select(HUC12) %>% 
  left_join(hucave, by = 'HUC12')

st_write(huc, 'data/metadata/hucave.shp', delete_layer = TRUE)

