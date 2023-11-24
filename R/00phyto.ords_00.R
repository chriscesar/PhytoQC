# 00phyto.ords_00.R
### Side project on phyto data ordination

## load packages
# load packages #
ld_pkgs <- c("tidyverse","purrr","readxl", "vegan")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

# set themes/global presets
cbPalette <- c("#0072B2","#e79f00","#009E73", "#9ad0f3", "#000000",
               "#D55E00", "#CC79A7", "#F0E442")
ppi <- 300
perms <- 9999 #number of permutations of analyses
theme_set(ggthemes::theme_few())###set theme for all ggplot objects

# load data 
df0 <- as_tibble(readxl::read_xlsx("data/Copy of PhytoChl_2000-2020 (WB+RegionalSeas).xlsx",
                                   sheet = "Abundance2000_2020USE"))

# filter, convert to wide to prep for ordination ####
df0 %>% 
  dplyr::select(.,-aphiaID, -Time) %>% #remove unneccessary cols
  filter(Year >= 2012) %>%  ##keep only samples gathered from 2012 onwards
  #drop NA values in TypologyCode to keep only WBs selected for analysis
  filter(!is.na(TypologyCode)) %>% 
  # names(.)
  # group_by(across(c(-Abundance))) %>% count(.) %>% View(.)
  pivot_wider(names_from = TaxonUSE,
              values_from = Abundance,
              values_fill=0) -> dfw

saveRDS(dfw, file = "data/out/phyto_WBs_sideProj.Rdat")

## ordinate data ###
## remove metadata
dfw %>% 
  dplyr::select(-c(`Regional Sea`:WBArea)) -> dftmp

### RUN ORDINATION ####
set.seed(pi);ord <- metaMDS(dftmp)
