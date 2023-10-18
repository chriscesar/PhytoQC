# 101_init_data_expl.R ####
# initial look at non-labswap data

# admin ####
# load packages #
ld_pkgs <- c("tidyverse","ggthemes","vegan")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

# set themes/global presets
cbPalette <- c("#0072B2","#e79f00","#009E73", "#9ad0f3", "#000000",
               "#D55E00", "#CC79A7", "#F0E442")
ppi <- 300
theme_set(ggthemes::theme_few())###set theme for all ggplot objects

# load data 
df0 <- as_tibble(read.csv("data/out/extracted_data_NonLabSwap.csv"))

# Widen data and run MDS ####

### filter by AnalysisType = Original & widen
df0_orig_w <- df0 %>% 
  filter(.,!is.na(SD00FileName)) %>% 
  filter(AnalysisType == "Original") %>% 
  dplyr::select(.,
                !c(prop,
                   QA01_TotAbund,
                   QA02_DomTax,
                   QA03_SharedTax,
                   QA04_Final)) %>% 
  dplyr::group_by(
    SD00FileName,
    SD01_AnaylsisLab,
    SD02_LabSwap,
    SD03_OriginalAnalyst,
    SD04_AnalysisDateOrig,
    SD05_NameOfSurvey_WFD,
    SD06_SampleDate,
    SD07_EAOldSiteCode,
    SD08_EAWIMSCode,
    SD09_InternalSampleID,
    SD10_AuditAnalyst,
    SD11_AuditDateBaseRec,
    SD12_AuditDateRepSub,
    SD13_Comments,
    SD14_WaterSampleVol_ml,
    SD15_SubSampleVol_ml,
    AnalysisType,
    Tax_Qual
  ) %>% 
  pivot_wider(names_from = Tax_Qual, values_from = dens,
              values_fill = 0) %>% 
  ungroup()

### create temporary file for ordination
# remove metadata
dftmp <- df0_orig_w %>% 
  select(!c(SD00FileName:AnalysisType))
set.seed(pi);ord <- vegan::metaMDS(dftmp,
                      try = 100,
                      trymax = 1000,
                      k = 3)
plot(ord)
