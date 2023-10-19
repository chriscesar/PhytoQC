# QAWBsummary.R ####
### Summaries of workbook QA tests

#### load packages ####
ld_pkgs <- c("tidyverse")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

### set metadata
theme_set(ggthemes::theme_few())###set theme for all ggplot objects
ppi <- 300 #image resolution
# colourblind friendly colour palette (RGB values also commented)
cbPalette <- c("#999999", #153/153/153
               "#E69F00",#230/159/000
               "#56B4E9",#086/180/233
               "#CC79A7", #204/121/167
               "#009E73",#000/158/115
               "#F0E442",#240/228/066
               "#0072B2",#000/114/178
               "#D55E00",#213/094/000
               
               "#444444", 
               "#C34D55",
               "#33A2C4",
               "#554C31",
               "#C5C221",
               "#5531A1",
               "#B32C55",
               "#BB3593" 
               
)

### load data ####
### taxon data
df0 <- as_tibble(read.csv("data/out/extracted_data_ALL.csv"))

### to do ###
# assign values to correct dates


#remove taxon data
df <- df0 %>% 
  dplyr::select(.,-c(Tax_Qual,dens,prop,
                     SD14_WaterSampleVol_ml,
                     SD15_SubSampleVol_ml,
                     AnalysisType)) %>% 
  filter(.,!SD02_LabSwap=="Yes") %>% 
  distinct()

# summarise Passes/fails ###
## overall
df %>% 
  group_by(SD01_AnaylsisLab,QA04_Final) %>% 
  count()

## Total abundance
df %>% 
  group_by(SD01_AnaylsisLab,QA01_TotAbund) %>% 
  count()

## Dominant taxa
df %>% 
  group_by(SD01_AnaylsisLab,QA02_DomTax) %>% 
  count()

## Shared taxa
df %>% 
  group_by(SD01_AnaylsisLab,QA03_SharedTax) %>% 
  count()


