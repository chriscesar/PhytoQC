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

### keep only non labswaps
df0NoSwp <- df0 %>%
  filter(., SD02_LabSwap == "No")

# assign values to correct dates
df0NoSwp %>% 
  filter(nchar(SD06_SampleDate)==5) -> df_tmp_5

df_tmp_5$SD06_SampleDate <- as.Date(as.numeric(df_tmp_5$SD06_SampleDate),
                                    origin = "1899-12-30")

df0 %>% 
  filter(nchar(SD06_SampleDate)!=5) -> df_tmp_n5

df_tmp_n5$SD06_SampleDate <- as.Date(df_tmp_n5$SD06_SampleDate,
                                     format = "%d/%m/%Y")

### join the data together
df0NoSwp <- rbind(df_tmp_n5,df_tmp_5); rm(df_tmp_n5,df_tmp_5)

#remove taxon data
df <- df0NoSwp %>% 
  dplyr::select(.,-c(Tax_Qual,dens,prop,
                     SD14_WaterSampleVol_ml,
                     SD15_SubSampleVol_ml,
                     AnalysisType)) %>% 
  filter(.,!SD02_LabSwap=="Yes") %>% 
  distinct()

# summarise Passes/fails ###
## how many samples by lab?
df %>% 
  group_by(SD01_AnaylsisLab) %>% 
  count()

# by year
df %>% 
  group_by(year = format(SD06_SampleDate, "%Y"),
           SD01_AnaylsisLab) %>% 
  count()

## overall
df %>% 
  group_by(SD01_AnaylsisLab,QA04_Final) %>% 
  count()

# by year
df %>% 
  group_by(year = format(SD06_SampleDate, "%Y"),
           SD01_AnaylsisLab,
           QA04_Final) %>% 
  count()

###

## Total abundance
df %>% 
  group_by(SD01_AnaylsisLab,QA01_TotAbund) %>% 
  count()

# by year
df %>% 
  group_by(year = format(SD06_SampleDate, "%Y"),
           SD01_AnaylsisLab,
           QA01_TotAbund) %>% 
  count()

###

## Dominant taxa
df %>% 
  group_by(SD01_AnaylsisLab,
           QA02_DomTax) %>% 
  count()

# by year
df %>% 
  group_by(year = format(SD06_SampleDate, "%Y"),
           SD01_AnaylsisLab,
           QA02_DomTax) %>% 
  count()

###

## Shared taxa
df %>% 
  group_by(SD01_AnaylsisLab,QA03_SharedTax) %>% 
  count()

#by year
df %>% 
  group_by(year = format(SD06_SampleDate, "%Y"),
           SD01_AnaylsisLab,
           QA03_SharedTax) %>% 
  count()

###

