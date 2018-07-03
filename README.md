# README

Materials for landcape modelling manuscript

## Data 

All files created using `R/dat_proc.R` (sometimes not from local source)

* `calicls.RData` statewide stream classifcations by nhd comid, sf object
* `caliclsplo.RData` statewide stream classifications by nhd comid, simplified and fortified for figures
* `caliexp.RData` all CSCI site expectations and comid stream classes for cali
* `calipsa.RData` sf object of california PSA regions
* `cnstrfrst.RData` results from random forest models of constraints in each region, importance only
* `comid_prd.RData` predicted csci quantiles for all comid in california
* `csci_raw.RData` data.frame of current csci scores by lat/lon
* `csci_comid.RData` csci scores by comid
* `ludat.RData` sf object of california land use by ag, urban, other
* `nhdplo.RData` statewide median csci prediction by nhd comid, simplified and fortified for figures
* `psalab.RData` data.frame and psa centroids, labels
* `rf_core.RData` core quantile random forest model of csci scores
* `rf_full.RData` full quantile random forest model of csci scores
* `scrs.RData` data.frame of CSCI scores for sites in SGR, from SGRRMP .rproj file
* `sensres.RData` results of sensitivity analysis statewide, list of lists for each region and statewide, classification and relative score changes
* `senssgr.RData` results of sensitivity analysis, SGR watershed only
* `sgrclslen.RData` length of stream classifications for all reaches in SGR watershed
* `shed.RData` sf object of San Gabriel watershed
* `spat.RData` sf object of nhd hydrolines with CSCI predictions for SGR, from SGRRMP .rproj file
* `strclslen.RData` sf object of all cali stream classes with length and PSA
* `typescrs.RData` summary of typical score expectations by region, landuse, certainty
* `typetab.RData` data.frame of table for site types