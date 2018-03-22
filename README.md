# README

Materials for landcape modelling manuscript

## Data 

All files created using `R/dat_proc.R` (sometimes not from local source)

* `calipsa.RData` sf object of california PSA regions
* `comid_prd.RData` predicted csci quantiles for all comid in california
* `csci_raw.RData` data.frame of current csci scores by lat/lon
* `csci_comid.RData` csci scores by comid
* `ludat.RData` sf object of california land use by ag, urban, other
* `psalab.RData` data.frame and psa centroids, labels
* `rf_core.RData` core quantile random forest model of csci scores
* `rf_full.RData` full quantile random forest model of csci scores
* `scrs.RData` data.frame of CSCI scores for sites in SGR, from SGRRMP .rproj file
* `shed.RData` sf object of San Gabriel watershed
* `spat.RData` sf object of nhd hydrolines with CSCI predictions for SGR, from SGRRMP .rproj file
* `typetab.RData` data.frame of table for site types