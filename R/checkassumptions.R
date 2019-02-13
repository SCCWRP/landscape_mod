# compare constraint classes to physical channel alteration
library(tidyverse)
library(readxl)
library(leaflet)
library(RColorBrewer)

load(file = '../data/caliexp.RData')

# phab data
phab <- read_excel('../data/raw/phab_all.xlsx') %>% 
  mutate(
    VariableResult = ifelse(VariableResult %in% 'Not Recorded', NA, VariableResult),
    val = ifelse(is.na(Result), VariableResult, Result), 
    val = as.numeric(val)
  ) %>% 
  dplyr::select(-Result, -VariableResult, -VariableName) %>% 
  group_by(StationCode, Variable) %>% 
  summarise(val = mean(val, na.rm = T)) %>% 
  filter(!is.na(val)) %>% 
  filter(val > -87) 

# phab, channel alteration
phab_ch <- phab %>% 
  filter(Variable %in% 'RBP_CHN') %>% 
  dplyr::select(-Variable) %>% 
  rename(ch_alt = val)

# all csci scores and classes
allscrs <- caliexp %>% 
  dplyr::select(COMID, csci, StationCode, strcls, perf) %>% 
  filter(!is.na(strcls)) %>%  
  mutate(StationCode = as.character(StationCode)) %>% 
  left_join(phab_ch, by = 'StationCode') %>% 
  filter(!is.na(ch_alt))

ggplot(allscrs, aes(x = strcls, y = ch_alt, fill = strcls, colour = strcls)) + 
  geom_boxplot(colour = 'black', alpha = 0.7, outlier.color = NA) +
  geom_point(position = position_jitter(w = 0.15, h = 0.15), alpha = 0.6) +
  scale_y_continuous('PHAB: Channel alteration') +
  theme_bw(base_family = 'serif', base_size = 16) +
  theme(
    axis.title.y = element_blank(), 
    legend.position = 'none'
  ) +
  scale_fill_manual(values = pal_exp(levels(allscrs$strcls))) +
  scale_colour_manual(values = pal_exp(levels(allscrs$strcls))) +
  coord_flip()

# tecolote creek 906M23302
allscrs %>% 
  filter(strcls %in% 'likely constrained' & ch_alt > 19)

##
# same analysis as above but using SMC channel mod data
# evaluate channel mod with constraints via contingency table

load(file = '../data/caliexp.RData')

chdat <- read.csv('C:/proj/SMC_report/ignore/csci.alg.df.csv', stringsAsFactors = FALSE) %>%
  dplyr::select(StationCode, ChannelClass_3ef) %>%
  rename(
    chcls = ChannelClass_3ef
  ) %>% 
  mutate(
    chcls = case_when(
      chcls %in% 'All or partially earthen' ~ 'Natural', 
      T ~ chcls
    )
  )

cmp <- caliexp %>% 
  dplyr::select(StationCode, strcls) %>% 
  inner_join(chdat, by = 'StationCode')

table(cmp[, -1])

###
# other strmcat variables that may explain deviation from predictions
library(tidyverse)

load(file = '../data/cnstrfrst.RData')
load(file = '../data/comid_prd.RData')
load(file = '../data/csci_comid.RData')

# observed and predicted data
obsdat <- csci_comid %>%
  dplyr::select(COMID, StationCode, SiteSet, SelectedSample, CSCI, PSA6c, RegImpQ) %>% 
  filter(SelectedSample %in% 'Selected')
prddat <- comid_prd %>% 
  dplyr::select(COMID, core0.50)

# combined data for residuals and streamcat vars
dat <- obsdat %>% 
  inner_join(prddat, by = 'COMID') %>% 
  dplyr::rename(
    observed = CSCI,
    predicted = core0.50
  ) %>% 
  mutate(resid = observed - predicted) %>% 
  dplyr::select(COMID, StationCode, SiteSet, PSA6c, resid) %>% 
  left_join(strmcat, by = 'COMID')

# corests
corests <- dat %>%
  gather('var', 'val', -COMID, -StationCode, -SiteSet, -PSA6c, -resid) %>% 
  group_by(PSA6c) %>% 
  nest %>% 
  mutate(
    corv = map(data, function(x){
      
      x %>% 
        group_by(var) %>% 
        summarise(
          cor = cor.test(resid, val, method = 'spearman')$estimate 
        ) %>% 
        arrange(cor) %>% 
        .[1:5, ]
      
    })
  ) %>% 
  dplyr::select(-data) %>% 
  unnest %>% 
  arrange(PSA6c, cor)

p1 <- ggplot(corests, aes(x = var, y = cor)) + 
  geom_bar(stat = 'identity') +
  facet_wrap(~PSA6c, ncol = 1, scales = 'free_y') + 
  coord_flip()


toplo2 <- dat %>% 
  dplyr::select(StationCode, PSA6c, resid, MWST_2008, DamNrmStorM3Ws, BFIWs, PctFire2006Ws) %>% 
  gather('var', 'val', -StationCode, -PSA6c, -resid) %>% 
  filter(
    PSA6c %in% 'SN' & var %in% 'DamNrmStorM3Ws' |
      PSA6c %in% 'NC' & var %in% 'MWST_2008' |
      PSA6c %in% 'SC' & var %in% 'BFIWs' |
      PSA6c %in% 'SC' & var %in% 'PctFire2006Ws'
  )

p2 <- ggplot(dat, aes(x = DamNrmStorM3Ws, y = resid)) + 
  geom_point() + 
  stat_smooth() #+
# scale_x_log10()

##
# in stream measures

library(tidyverse)
library(readxl)

load(file = '../data/cnstrfrst.RData')
load(file = '../data/comid_prd.RData')
load(file = '../data/csci_comid.RData')

# phab data
phab <- read_excel('../data/raw/phab_all.xlsx') %>% 
  mutate(
    VariableResult = ifelse(VariableResult %in% 'Not Recorded', NA, VariableResult),
    val = ifelse(is.na(Result), VariableResult, Result), 
    val = as.numeric(val)
  ) %>% 
  dplyr::select(-Result, -VariableResult, -VariableName) %>% 
  group_by(StationCode, Variable) %>% 
  summarise(val = mean(val, na.rm = T)) %>% 
  filter(!is.na(val)) %>% 
  filter(val > -87)

# nitrogen
nitr <- read.csv('../data/raw/smc.swamp_mastertable_SMR_CRAM.CSV', stringsAsFactors = F) %>% 
  rename(
    val = Nitrogen_Total_mgPerL
  ) %>% 
  filter(!is.na(val)) %>% 
  group_by(StationCode) %>% 
  summarise(val = mean(val, na.rm = T)) %>% 
  mutate(Variable = 'tn') %>% 
  dplyr::select(StationCode, Variable, val)

# combined extra data for phab and nitrogen, long format
exdat <- phab %>% 
  bind_rows(nitr)

# observed and predicted data
obsdat <- csci_comid %>%
  dplyr::select(COMID, StationCode, SiteSet, SelectedSample, CSCI, PSA6c, RegImpQ) %>% 
  filter(SelectedSample %in% 'Selected')
prddat <- comid_prd %>% 
  dplyr::select(COMID, core0.50)

# combined data for residuals and streamcat vars
dat <- obsdat %>% 
  inner_join(prddat, by = 'COMID') %>% 
  dplyr::rename(
    observed = CSCI,
    predicted = core0.50
  ) %>% 
  mutate(resid = observed - predicted) %>% 
  dplyr::select(COMID, StationCode, SiteSet, PSA6c, resid)

# joinnnn  
toeval <- dat %>% 
  inner_join(exdat, by = 'StationCode')

ggplot(toeval, aes(x = val, y = resid)) + 
  geom_point() +
  facet_wrap(PSA6c ~ Variable, scales = 'free_x') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_smooth(method = 'lm') +
  # scale_x_log10() +
  theme_bw()

######
# look at score differences between perennial/intermittent
# flow indicators from Chad, RB9

source('R/funcs.R')

data(caliexp)

# format flow indication data (filter those where sop status not met)
flos <- read_excel('raw/FromChadTable 1b.xlsx', sheet = 'Sheet2') %>% 
  select(StationCode, Flow_Simple, SOP_Status) %>% 
  rename(flow = Flow_Simple) %>% 
  filter(grepl('^SOP met', SOP_Status)) %>% 
  select(-SOP_Status) %>% 
  unique %>% 
  mutate(flow = factor(flow, levels = c('P', 'I'), labels = c('Perennial', 'Intermittent')))

# combine flow indication with csci, landscape mod ests, get pvals from t.tests
dat <- caliexp %>% 
  select(StationCode, datcut, strcls, csci) %>% 
  unnest %>% 
  mutate(StationCode = as.character(StationCode)) %>% 
  inner_join(flos, by = 'StationCode') %>% 
  spread(var, val) %>% 
  gather('var', 'val', csci, core0.10, core0.50, core0.90) %>% 
  mutate(
    StationCode = as.character(StationCode),
    var = factor(var, 
                 levels = c('csci', 'core0.10', 'core0.50', 'core0.90'),
                 labels = c('CSCI', 'tenth', 'median', 'ninetieth')
    )
  ) %>% 
  group_by(var) %>% 
  nest %>% 
  mutate(
    pvl = purrr::map(data, function(x){
      int <- x %>% filter(flow %in% 'Intermittent')
      per <- x %>% filter(flow %in% 'Perennial')
      
      tst <- t.test(int$val, per$val, var.equal = T)
      pvl <- tst$p.value %>% 
        round(2) %>% 
        paste('p =', .)
      
      return(pvl)
      
    })
  ) %>% 
  unnest(pvl) %>% 
  unnest(data) %>% 
  unite('varp', var, pvl, remove = F, sep = ', ') 

# boxplot comparisons  
ggplot(dat, aes(x = flow, y = val)) + 
  # geom_violin(draw_quantiles = 0.5, alpha = 0.2, fill = 'lightgrey') + 
  geom_boxplot(alpha = 0.5, fill = 'lightgrey', outlier.colour = NA) + 
  geom_point(aes(fill = strcls), position = 'jitter', col = 'lightgrey', alpha = 0.9, pch = 21, size = 3) + 
  geom_hline(yintercept = 0.79, linetype = 'dotted', colour = 'tomato1', size = 2) + 
  theme_bw(base_family = 'serif', base_size = 16) +
  facet_wrap(~varp, ncol = 4) + 
  scale_fill_manual(values = pal_exp(levels(dat$strcls)), drop = F) +
  ylab('CSCI') +
  theme(
    axis.title.x = element_blank(), 
    strip.background = element_blank(), 
    axis.text.x = element_text(angle = 20, hjust = 1), 
    legend.title = element_blank()
  )

# this looks at relationship between observed/predicted csci between perennial/intermittent
# combine flow indication with csci, format for plots
dat <- caliexp %>% 
  select(StationCode, datcut, strcls, csci) %>% 
  unnest %>% 
  mutate(StationCode = as.character(StationCode)) %>% 
  inner_join(flos, by = 'StationCode') %>% 
  spread(var, val) %>% 
  rename(
    CSCI = csci, 
    median = core0.50, 
    tenth = core0.10, 
    ninetieth = core0.90
  )

# plot  
ggplot(dat, aes(x = median, y = CSCI)) + 
  geom_errorbar(aes(ymin = tenth, ymax = ninetieth, colour = strcls), width = 0, size = 1) + 
  geom_point(aes(shape = flow), size = 2.5) +
  geom_smooth(method = 'lm', aes(linetype = flow), se = F, colour = 'black') + 
  theme_bw(base_family = 'serif', base_size = 16) +
  scale_colour_manual('', values = pal_exp(levels(dat$strcls)), drop = F) +
  xlab('Median predicted CSCI') +
  ylab('Observed CSCI') +
  theme(
    legend.title = element_blank()
  )


