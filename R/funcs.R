# function for formatting p-values in tables
p_ast <- function(x){
  
  sig_cats <- c('**', '*', 'ns')
  sig_vals <- c(-Inf, 0.005, 0.05, Inf)
  
  out <- cut(x, breaks = sig_vals, labels = sig_cats, right = FALSE)
  out <- as.character(out)
  
  return(out)
  
}

######
# size palettes
pal_siz <- function(levin){
  
  out <- tibble(
    levels = c('High', 'Medium', 'Low', 'baseline'), 
    sizes = c(4, 3, 2, 1)
  ) %>% 
    filter(levels %in% levin) %>% 
    .$sizes
  
  return(out)
  
}

######
# color palettes

# color palette for priorities
pal_pri <- colorFactor(
  palette = RColorBrewer::brewer.pal(9, 'Greys')[c(9, 6, 4, 1)],
  na.color = 'yellow',
  levels = c('High', 'Medium', 'Low', 'baseline'))

# color palette for csci scores
pal <- colorNumeric(
  palette = c('#d7191c', '#abd9e9', '#2c7bb6'),
  na.color = 'yellow',
  domain = c(0, 1.4))

# color palette for priorities
pal_pri <- colorFactor(
  palette = RColorBrewer::brewer.pal(9, 'Greys')[c(9, 6, 4, 1)],
  na.color = 'yellow',
  levels = c('High', 'Medium', 'Low', 'baseline'))

# color palette for stream expectations
pal_exp <- colorFactor(
  palette = brewer.pal(9, 'Paired')[c(2, 1, 5, 6)],
  na.color = 'yellow',
  levels = c('likely unconstrained', 'possibly unconstrained', 'possibly constrained', 'likely constrained'))

# color palette for CSCI scoring performance
pal_prf <- colorFactor(
  palette = c(
    colorRampPalette(brewer.pal(9, 'Blues'))(100)[c(90, 74, 58)],
    colorRampPalette(brewer.pal(9, 'Blues'))(100)[c(42, 26, 10)],
    colorRampPalette(brewer.pal(9, 'Reds'))(100)[c(42, 26, 10)],
    colorRampPalette(brewer.pal(9, 'Reds'))(100)[c(90, 74, 58)]
  ),
  na.color = 'yellow',
  levels = c(
    'over scoring (lu)', 'expected (lu)', 'under scoring (lu)',
    'over scoring (pu)', 'expected (pu)', 'under scoring (pu)',
    'over scoring (pc)', 'expected (pc)', 'under scoring (pc)',
    'over scoring (lc)', 'expected (lc)', 'under scoring (lc)')
)

# color palette for CSCI type
pal_typ <- colorFactor(
  palette = RColorBrewer::brewer.pal(11, 'Spectral'),#hue_pal()(100),
  na.color = 'yellow',
  domain = paste0('Type', sprintf('%02d', seq(1:16)))
)

# color palette for csci score differences
pal_difr <- colorNumeric(
  palette = c('black', 'purple', 'white', 'darkgreen', 'black'),
  na.color = 'yellow',
  domain = c(-0.6, 0.6))

# pal lu
pal_lu <- colorFactor(
  palette = c('grey20', 'grey40', 'grey60', 'khaki3', 'khaki2'),
  domain = c('Urban: hi', 'Urban: md', 'Urban: lo', 'Open: forest', 'Open: scrub')
)

######
# get legend from an existing ggplot object
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#' Get biological expectation with changing threshold and likelihoods, same as getcls but only returns class designation
plot_strcls <- function(datin, thrsh = 0.79, tails = 0.05, lbs = list('likely unconstrained' = 0, 'possibly unconstrained' = 1, 'possibly constrained' = 2, 'likely constrained' = 3)){
  
  # sanity check
  if(tails > 0.5)
    stop('tails must be less than 0.5')
  
  dat <- datin %>% 
    gather('var', 'val', -Site) %>% 
    group_by(Site) %>% 
    nest %>% 
    mutate(
      datcut = purrr::map(data, function(x){
        
        tofilt <- data.frame(
          val = seq(x$val[1], x$val[2], length = 19), 
          qts = round(seq(0.05, 0.95, by = 0.05), 2)
        )
        
        # get quantile labels to filter
        lovl <- round(tails, 2)
        hivl <- round(1 - tails, 2)
        
        # filter by quantile labels and median
        rngs <- tofilt %>% 
          filter(qts %in% c(lovl, 0.5, hivl))
        
        return(rngs)
        
      }),
      strcls_int = purrr::map(datcut, function(x){
        
        # get threshold range
        rngs <- x$val
        
        # find if above/below or covering thrsh
        cls <- findInterval(thrsh, rngs)
        
        return(cls)
        
      })
      
    ) 
  
  # subset lbs by those in interval
  lbs <- unique(dat$strcls_int) %>% 
    na.omit %>% 
    as.numeric %>% 
    match(unlist(lbs)) %>% 
    lbs[.] %>% 
    .[names(sort(unlist(.)))]
  
  # strcls as correct factor levels
  dat <- dat %>%
    mutate(strcls = factor(strcls_int, levels = unlist(lbs), labels = names(lbs)))
  
  return(dat)
  
}

#' Get CSCI StationCode expectations and performance classication
plot_exp <- function(datin, scrs, thrsh = 0.79, tails = 0.05, lbs = list('over scoring' = 2, 'expected' = 1, 'under scoring' = 0),
                     ...
){
  
  # join plot strcls data and scrs
  incl <- datin %>% 
    plot_strcls(thrsh = thrsh, tails = tails, ...) %>% 
    left_join(scrs, by = 'Site') 
  
  # get CSCI performance (over/under)
  incl <- incl %>% 
    mutate(
      perf = pmap(list(datcut, csci), function(datcut, csci){
        
        # within datcut interval
        prf <- findInterval(csci, range(datcut$val))
        
        return(prf)
        
      })
    ) %>% 
    unnest(perf) %>%
    mutate(bythrsh = ifelse(csci < thrsh, 0, 1)) %>%
    unite('typeprf', strcls_int, perf, bythrsh, remove = FALSE) %>% 
    mutate(
      typelv = typ_lbs(typeprf, thrsh = thrsh, tails = tails),
      typeoc = typ_lbs(typeprf, thrsh = thrsh, tails = tails, obs_sc = T)
    ) %>% 
    dplyr::select(-typeprf, -bythrsh)
  
  # subset lbs by those in interval
  lbs <- unique(incl$perf) %>% 
    na.omit %>% 
    as.numeric %>% 
    match(unlist(lbs)) %>% 
    lbs[.] %>% 
    .[names(sort(unlist(.)))]
  
  # perf as correct factor levels
  incl <- incl %>%
    mutate(perf = factor(perf, levels = unlist(lbs), labels = names(lbs)))
  
  return(incl)
  
}

#' Get type labels from three level code 
#'
#' @param vec chr string vector of codes to typify 
#' @param thrsh numeric for CSCI scoring thresholds
#' @param tails numeric for tails to truncate expectations for overlap with thrsh
#' @param obs_sc logical if observed score text qualifiers are returned
#' @param get_cds logical indicating if three level codes as list is returned
#' @details  The three level codes for type are (0, 1, 2, 3), (0, 1, 2), and (0, 1).  The first level is likely unconstrained (0), possibly unconstrained (1), possibly constrained (2), likely constrianed (3); the second level is under-scoring (0), as expected (1), and over-scoring (2); the third level is below threshold (0) and above threshold (1)
typ_lbs <- function(vec = NULL, thrsh = 0.79, tails = 0.05, obs_sc = FALSE, get_cds = FALSE){
  
  # type labels from codes in vec
  lbs <- list(
    Type01 = '0_2_1', 
    Type02 = '0_1_1',
    Type03 = '0_0_1',
    Type04 = '0_0_0',
    Type05 = '1_2_1',
    Type06 = '1_1_1',
    Type07 = '1_1_0',  
    Type08 = '1_0_0',
    Type09 = '2_2_1',
    Type10 = '2_1_1',
    Type11 = '2_1_0',
    Type12 = '2_0_0',
    Type13 = '3_2_1',
    Type14 = '3_2_0', 
    Type15 = '3_1_0', 
    Type16 = '3_0_0'
  )
  
  if(get_cds) return(lbs)
  
  # kill NA entries
  vec <- ifelse(grepl('NA', vec), NA, vec)
  
  # get tails in chr format
  lovl <- 100 * tails
  hivl <- 100 * (1 -  tails) 
  vls <- c(lovl, hivl) %>% 
    round(., 0) %>% 
    paste0(., 'th')
  
  # subset lbs by those in vec
  lbs <- unique(vec) %>% 
    na.omit %>%
    match(unlist(lbs)) %>% 
    lbs[.] %>% 
    .[sort(names(.))]
  
  # assign factor levels to vec
  vec <- factor(vec, levels = unlist(lbs), labels = names(lbs))
  
  if(!obs_sc) return(vec)
  
  # observed score chr types
  obs_sc <- list(
    Type01 = paste('>=', vls[2]), 
    Type02 = paste(vls[1], 'to', vls[2]),
    Type03 = paste(thrsh, 'to', vls[1]),
    Type04 = paste('<', thrsh),
    Type05 = paste('>=', vls[2]),
    Type06 = paste(thrsh, 'to', vls[2]),
    Type07 = paste(vls[1], 'to', thrsh), 
    Type08 = paste('<', vls[1]),
    Type09 = paste('>=', vls[2]),
    Type10 = paste(thrsh, 'to', vls[2]),
    Type11 = paste(vls[1], 'to', thrsh), 
    Type12 = paste('<', vls[1]),
    Type13 = paste('>=', thrsh),
    Type14 = paste(vls[2], 'to', thrsh),
    Type15 = paste(vls[1], 'to', vls[2]),
    Type16 = paste('<', vls[1])
  ) %>% 
    enframe('vec', 'obs_sc') %>% 
    unnest
  
  # assign factor levels to vec
  vec <- vec %>% 
    data.frame(vec = .) %>% 
    mutate(vec = as.character(vec)) %>% 
    left_join(., obs_sc, by = 'vec') %>% 
    .$obs_sc
  
  return(vec)
  
}

#' Get performance multi classifications
#' Input is scr_exp() reactive from app
get_perf_mlt <- function(scr_exp, lbs = c('over scoring (lu)', 'expected (lu)', 'under scoring (lu)',
                                          'over scoring (pu)', 'expected (pu)', 'under scoring (pu)',
                                          'over scoring (pc)', 'expected (pc)', 'under scoring (pc)',
                                          'over scoring (lc)', 'expected (lc)', 'under scoring (lc)')
){
  
  # format perf_mlt as column combos of perf and strcls
  out <- scr_exp %>% 
    mutate(
      torm = purrr::map(strcls, function(chr){
        
        if(!is.na(chr)){
          as.character(chr) %>% 
            strsplit(., ' ') %>% 
            lapply(function(x) substr(x, 1, 1) %>% paste(collapse = '') %>% paste0('(', ., ')')) %>% 
            .[[1]]
        } else {
          NA
        }
        
      }),
      torm = unlist(torm)
    ) %>% 
    unite('perf_mlt', perf, torm, sep = ' ', remove = F) %>% 
    mutate(perf_mlt = ifelse(grepl('^NA|NA$', perf_mlt), NA, perf_mlt))
  
  # subset lbs by those in data
  lbs2 <- unique(out$perf_mlt) %>% 
    na.omit %>% 
    match(lbs) %>% 
    lbs[.] 
  
  # get correct order
  lbs <- lbs[lbs %in% lbs2]
  lbs2 <- lbs2[match(lbs, lbs2)]
  
  # reassign factor levels
  out <- out %>% 
    mutate(perf_mlt = factor(perf_mlt, levels = lbs2))
  
  return(out)
  
}

# master funcion that classifies everything using only thrsh and tails
proc_all <- function(datin, scrs, thrsh, tails){
  
  out <- plot_exp(datin, scrs, thrsh = thrsh, tails = tails) %>% 
    get_perf_mlt %>% 
    mutate(datcut = purrr::map(datcut, function(x){
      
      out <- x %>% 
        filter(qts != 0.5) %>% 
        mutate(lev = c('minv_qt', 'maxv_qt')) %>% 
        dplyr::select(-qts) %>% 
        spread(lev, val)
      
      return(out)
      
    })
    ) %>% 
    mutate(data = purrr::map(data, spread, var, val)) %>% 
    unnest %>% 
    rename(
      `Relative\nperformance` = perf_mlt,
      `Stream class` = strcls, 
      `CSCI score` = csci
    )
  
  return(out)
  
}

#' Summarize data from site_exp by counts and all types, used for table in app
#'
#' @param datin output data.frame from site_exp
#' @param thrsh numeric for CSCI scoring thresholds
#' @param tails numeric for tails to truncate expectations for overlap with thrsh
#' @param obs_sc logical if observed score text qualifiers are returned
#' @param lbs_str chr string labels for stream comid expectation as likely constrained, possibly constrained, possibly unconstrained, or likely unconstrained
#' @param lbs_sta chr string labels for site/station performance as over, expected, or under scoring
get_tab <- function(datin, thrsh = 0.79, tails = 0.05, lbs_str = list('likely unconstrained' = 0, 'possibly unconstrained' = 1, 'possibly constrained' = 2, 'likely constrained' = 3), lbs_sta = list('over scoring' = 2, 'expected' = 1, 'under scoring' = 0)){ 
  
  # typeoc labels
  typeoc <- typ_lbs(get_cds = T) %>% 
    unlist %>% 
    typ_lbs(., thrsh = thrsh, tails = tails, obs_sc = T)
  
  # type labels from codes
  lbs <- typ_lbs(get_cds = T) %>%
    enframe('typelv', 'codes') %>%
    unnest %>%
    separate(codes, c('strcls', 'perf', 'thrsh'), sep = '_', remove = FALSE) %>% 
    mutate(
      typelv = factor(typelv, levels = typelv),
      strcls = factor(strcls, levels = unlist(lbs_str), labels = names(lbs_str)),
      perf = factor(perf, levels = unlist(lbs_sta), labels = names(lbs_sta))
    ) %>% 
    mutate(typeoc = typeoc) %>% 
    dplyr::select(-codes, -thrsh)
  
  # get type levels and convert all to character
  lvs <- levels(lbs$typelv)
  lbs <- lbs %>% 
    mutate_all(as.character)
  
  # types from observed data, join with complete table
  totab <- datin %>%
    dplyr::select(strcls, perf, typeoc, typelv) %>%
    group_by(strcls, perf, typeoc, typelv) %>%
    summarise(Sites = n()) %>%
    na.omit %>%
    ungroup %>% 
    mutate_all(as.character) %>% 
    left_join(lbs, ., by = c('strcls', 'perf', 'typeoc', 'typelv')) %>% 
    mutate(
      Sites = ifelse(is.na(Sites), 0, Sites), 
      Sites = as.numeric(Sites),
      typelv = factor(typelv, levels = lvs)
    ) %>% 
    arrange(`typelv`) %>%
    rename(
      `Reach expectation` = strcls,
      `Site performance` = perf,
      `Observed score` = typeoc,
      Type = typelv
    ) %>%
    mutate(
      Type = gsub('^Type|^Type0', '', Type),
      Type = as.numeric(Type)
    ) %>% 
    dplyr::select(`Reach expectation`, `Site performance`, `Observed score`, Type, Sites)
  
  return(totab)
  
}

#' Get biological expectation with changing threshold and likelihoods
#' 
#' @param datin sf object with stream COMIDS and quantile expectations
#' @param thrsh numeric for CSCI scoring thresholds
#' @param tails numeric for tails to truncate expectations for overlap with thrsh
#' @param modls chr string for selecting core (simple) or full model for expectations
#' @param lbs chr string labels for interval classifications
#' 
#' @return a nested data frame sorted by increaesing median value of expected score of COMID and nested columns as original data, cut data by tails, and sream classification (strcls).  The strcls column indicates if the ranges in datcut are within, above, or below those defind by thrsh.
getcls <- function(datin, thrsh = 0.79, tails = 0.05,  modls = c('core', 'full'), lbs = list('likely unconstrained' = 0, 'possibly unconstrained' = 1, 'possibly constrained' = 2, 'likely constrained' = 3)){
  
  # sanity check
  if(tails >= 0.5)
    stop('tails must be less 0.5')
  
  # models argument
  modls <- match.arg(modls)
  
  dat <- datin
  st_geometry(dat) <- NULL
  dat <- dat %>% 
    dplyr::select(matches(paste0('^COMID$|^', modls, '0'))) %>% 
    gather('var', 'val', -COMID) %>% 
    arrange(COMID, var) %>% 
    group_by(COMID) %>% 
    nest %>% 
    mutate(
      datcut = purrr::map(data, function(x){
        
        # get quantile labels to filter
        lovl <- 100 * tails
        hivl <- 100 * (1 -  tails) 
        vls <- c(lovl, 50, hivl) %>% 
          str_pad(2, pad = '0') %>% 
          paste0(modls, '0.', .)
        
        # filter by quantile labels and median
        x <- x %>% 
          filter(var %in% vls)
        
        return(x)
        
      }),
      
      medv = purrr::map(data, ~ filter(.x, var %in% paste0(modls, '0.50')) %>% .$val), 
      strcls_int = purrr::map(datcut, function(x){
        
        # return NA if any zero values in predictions
        if(any(is.na(x$val))){
          
          cls <- NA
          return(cls)
          
        } 
        
        # find if above/below or covering thrsh
        cls <- findInterval(thrsh, x$val)
        
        return(cls)
        
      })
      
    ) %>% 
    unnest(medv) %>% 
    arrange(medv) %>% 
    mutate(COMID = factor(COMID, levels = COMID)) %>% 
    unnest(strcls_int) 
  
  # subset lbs by those in interval
  lbs <- unique(dat$strcls_int) %>% 
    na.omit %>% 
    as.numeric %>% 
    match(unlist(lbs)) %>% 
    lbs[.]
  
  # strcls as correct factor levels
  dat <- dat %>%
    mutate(strcls = factor(strcls_int, levels = unlist(lbs), labels = names(lbs)))
  
  return(dat)
  
}

#' Get biological expectation with changing threshold and likelihoods, same as getcls but only returns class designation
#' 
#' @param datin sf object with stream COMIDS and quantile expectations
#' @param thrsh numeric for CSCI scoring thresholds
#' @param tails numeric for tails to truncate expectations for overlap with thrsh
#' @param modls chr string for selecting core (simple) or full model for expectations
#' @param lbs chr string labels for interval classifications
#' 
#' @return a nested data frame sorted by increaesing median value of expected score of COMID and nested columns as original data, cut data by tails, and sream classification (strcls).  The strcls column indicates if the ranges in datcut are within, above, or below those defind by thrsh.
getcls2 <- function(datin, thrsh = 0.79, tails = 0.05, modls = c('core', 'full'), lbs = list('likely unconstrained' = 0, 'possibly unconstrained' = 1, 'possibly constrained' = 2, 'likely constrained' = 3)){
  
  # sanity check
  if(tails >= 0.5)
    stop('tails must be less than 0.5')
  
  # models argument
  modls <- match.arg(modls)
  
  dat <- datin
  st_geometry(dat) <- NULL
  dat <- dat %>% 
    dplyr::select(matches(paste0('^COMID$|^', modls, '0'))) %>% 
    gather('var', 'val', -COMID) %>% 
    arrange(COMID, var) %>% 
    group_by(COMID) %>% 
    nest %>% 
    mutate(
      strcls_int = purrr::map(data, function(x){
        
        # return NA if any zero values in predictions
        if(any(is.na(x$val))){
          
          cls <- NA
          return(cls)
          
        } 
        
        # get quantile labels to filter
        lovl <- 100 * tails
        hivl <- 100 * (1 -  tails) 
        vls <- c(lovl, 50, hivl) %>% 
          str_pad(2, pad = '0') %>% 
          paste0(modls, '0.', .)
        
        # filter by quantile labels and median
        ints <- x %>% 
          filter(var %in% vls) %>% 
          .$val 
        
        cls <- findInterval(thrsh, ints)
        
        return(cls)
        
      })
      
    ) %>% 
    dplyr::select(-data) %>% 
    unnest 
  
  # subset lbs by those in interval
  lbs <- unique(dat$strcls_int) %>% 
    na.omit %>% 
    as.numeric %>% 
    match(unlist(lbs)) %>% 
    lbs[.] %>% 
    .[names(sort(unlist(.)))]
  
  # strcls as correct factor levels
  dat <- dat %>%
    mutate(strcls = factor(strcls_int, levels = unlist(lbs), labels = names(lbs)))
  
  return(dat)
  
}

#' Get CSCI StationCode expectations and performance classication
#'
#' @param datin sf object with stream COMIDS and quantile expectations
#' @param scrs CSCI scores by COMID and StationCode
#' @param thrsh numeric for CSCI scoring thresholds
#' @param tails numeric for tails to truncate expectations for overlap with thrsh
#' @param lbs chr string labels for site performance as over, expected, or under scoring
#' @param ... additional arguments passed to getcls
#' 
site_exp <- function(datin, scrs, thrsh = 0.79, tails = 0.05, lbs = list('over scoring' = 2, 'expected' = 1, 'under scoring' = 0),
                     ...
){
  
  # site csci scores
  scrs <- scrs %>% 
    mutate(COMID = as.character(COMID)) %>% 
    dplyr::select(COMID, StationCode, csci, lat, long)
  
  # filter comids with csci scores, classify, join with scores
  incl <- datin %>% 
    filter(COMID %in% scrs$COMID) %>% 
    getcls(thrsh = thrsh, tails = tails, ...) %>% 
    mutate(COMID = as.character(COMID)) %>% 
    left_join(scrs, by = 'COMID') %>% 
    arrange(medv) %>% 
    mutate(StationCode = factor(StationCode, levels = unique(StationCode))) %>% 
    dplyr::select(-medv)
  
  # get CSCI performance (over/under)
  incl <- incl %>% 
    mutate(
      perf = pmap(list(datcut, csci), function(datcut, csci){
        
        # return NA if any zero values in predictions
        if(any(is.na(datcut$val))){
          
          prf <- NA
          return(prf)
          
        } 
        
        # within datcut interval
        prf <- findInterval(csci, range(datcut$val))
        
        return(prf)
        
      })
    ) %>% 
    unnest(perf) %>%
    mutate(bythrsh = ifelse(csci < thrsh, 0, 1)) %>%
    unite('typeprf', strcls_int, perf, bythrsh, remove = FALSE) %>% 
    mutate(
      typelv = typ_lbs(typeprf, thrsh = thrsh, tails = tails),
      typeoc = typ_lbs(typeprf, thrsh = thrsh, tails = tails, obs_sc = T)
    ) %>% 
    dplyr::select(-typeprf, -bythrsh)
  
  # subset lbs by those in interval
  lbs <- unique(incl$perf) %>% 
    na.omit %>% 
    as.numeric %>% 
    match(unlist(lbs)) %>% 
    lbs[.] %>% 
    .[names(sort(unlist(.)))]
  
  # perf as correct factor levels
  incl <- incl %>%
    mutate(perf = factor(perf, levels = unlist(lbs), labels = names(lbs)))
  
  return(incl)
  
}

#' Get priority action from inputs mapped to data
#' 
#' @param plot_ex output reactive for plot selector
#' @param scr_exp_map output reactive for score expectations on map
#' 
get_pri_inp <- function(input, plot_ex, scr_exp_map){
  
  # get plot example site, types
  ex_jn <- plot_ex %>% 
    dplyr::select(Site, typelv)
  
  # site classications
  scr_exp_map <- scr_exp_map %>%
    mutate(typelv = as.character(typelv))
  
  # format site priorities from input
  scr_pri <- list(
    `Site 1` = input$`Site 1`,
    `Site 2` = input$`Site 2`,
    `Site 3` = input$`Site 3`,
    `Site 4` = input$`Site 4`,
    `Site 5` = input$`Site 5`,
    `Site 6` = input$`Site 6`,
    `Site 7` = input$`Site 7`,
    `Site 8` = input$`Site 8`,
    `Site 9` = input$`Site 9`,
    `Site 10` = input$`Site 10`,
    `Site 11` = input$`Site 11`,
    `Site 12` = input$`Site 12`,
    `Site 13` = input$`Site 13`,
    `Site 14` = input$`Site 14`,
    `Site 15` = input$`Site 15`,
    `Site 16` = input$`Site 16`
  ) %>% 
    enframe('Site', 'Priority') %>%
    mutate(Priority = purrr::map(Priority, ~ ifelse(is.null(.x), 'Do nothing', .x))) %>%
    unnest %>%
    mutate(
      Site = factor(Site, levels = levels(ex_jn$Site))
    ) %>%
    left_join(ex_jn, by = 'Site') %>%
    mutate(typelv = as.character(typelv)) %>%
    dplyr::select(-Site) %>%
    left_join(scr_exp_map, ., by = 'typelv') %>%
    dplyr::select(StationCode, COMID, csci, perf_mlt, typelv, lat, long, Priority) %>% 
    split(.$Priority) %>% 
    enframe('Priority')
  
  return(scr_pri)
  
}