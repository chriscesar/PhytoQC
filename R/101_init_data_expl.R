# 101_init_data_expl.R ####
# initial look at non-labswap data

# admin ####
# load packages #
ld_pkgs <- c("tidyverse","vegan", "rgl", "mvabund","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

tictoc::tic.clearlog();tic("TOTAL TIME")
tic("ADMIN: Set metadata; load data")
source("R/00_setMeta.R")
# load data 
df0 <- as_tibble(read.csv("data_processed/extracted_data_NonLabSwap.csv"))

# Widen data and run MDS ####

### filter by AnalysisType = "Original" & widen
df0_orig_w <- df0 %>% 
  filter(.,!is.na(SD00_FileName)) %>% 
  filter(AnalysisType == "Original") %>% 
  dplyr::select(., # remove unnecessary columns
                !c(prop,
                   QA01_TotAbund,
                   QA02_DomTax,
                   QA03_SharedTax,
                   QA04_Final)) %>% 
  dplyr::group_by(
    SD00_FileName,
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
              values_fill = 0) %>% # widen data
  ungroup()

## save wide data ###
write.csv(df0_orig_w, file="data_processed/extracted_data_NonLabSwap_wide.csv")
toc(log=TRUE)

tic("Run ordination")
### create temporary file for ordination
# remove metadata
dftmp <- df0_orig_w %>% 
  dplyr::select(!c(SD00_FileName:AnalysisType))
set.seed(pi);ord <- vegan::metaMDS(dftmp,
                      try = 20,
                      trymax = 1000,
                      k = 3)
#plot(ord)
toc(log=TRUE)

tic("Plot ordination using ggplot")
## extract data for ggplot2 ####
#### extract ordination axes ####
scores_site <- df0_orig_w %>% 
  dplyr::select(c(1:17))
tmp_sites <- as_tibble(as.data.frame(scores(ord,display = "site")))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
scores_site$NMDS1 <- tmp_sites$NMDS1 #location of individual samples in NMDS space
scores_site$NMDS2 <- tmp_sites$NMDS2 #location of individual samples in NMDS space
scores_site$NMDS3 <- tmp_sites$NMDS3 #location of individual samples in NMDS space
#assign colours to variable 'SD01_AnalysisLab'
unq <- levels(as.factor(scores_site$SD01_AnaylsisLab)) ## extract unique Lab values
scores_site$LabCol <- case_when(scores_site$SD01_AnaylsisLab == unq[1] ~ cbPalette[1],
                                scores_site$SD01_AnaylsisLab == unq[2] ~ cbPalette[2],
                                   TRUE ~ NA_character_)
saveRDS(scores_site, file = "data_processed/scores_site3d.Rdata")
rm(tmp_sites)

# Using the scores function from vegan to extract the species scores and
# convert to a data.frame
scores_species <- as.data.frame(scores(ord,display = "species"))
scores_species$lbfull <-  row.names(scores_species)
scores_species$lb <-  make.cepnames(row.names(scores_species))#shorten names
scores_species$lb <- gsub("µ","u",scores_species$lb)

saveRDS(scores_species, file = "data_processed/scores_species3d.Rdata")

plot3d(x = scores_site$NMDS1,
       y = scores_site$NMDS2,
       z = scores_site$NMDS3,
       col=scores_site$LabCol,
       type="s", size = 1,
       xlab = "NMDS1",ylab = "NMDS2",zlab = "NMDS3")

text3d(x = scores_species$NMDS1,
       y = scores_species$NMDS2,
       z = scores_species$NMDS3,
       scores_species$lb, size=1,
       cex = .75, col="grey")

axes3d();title3d(xlab="NMDS1",
                 ylab="NMDS2",
                 zlab="NMDS3",
                 font=2)
toc(log=TRUE)

# statistical comparisons ####
### do assemblages differ by lab?
## vegan version
tic("Run 1 way PERMANOVA")
(fit01 <- vegan::adonis2(dftmp ~ df0_orig_w$SD01_AnaylsisLab,
                         permutations = perms))
toc(log=TRUE)

## mvabund
tic("Run mvabund model: raw counts")
mvdftmp <- mvabund(dftmp)
fit_raw <- manyglm(mvdftmp ~ df0_orig_w$SD01_AnaylsisLab,family="negative.binomial")
anova_fit_raw <- mvabund::anova.manyglm(fit_raw,p.uni = "adjusted")
toc(log=TRUE)

tic("Run mvabund model: sqrt transformed")
fit_sqrt <- manyglm(sqrt(mvdftmp) ~ df0_orig_w$SD01_AnaylsisLab,family="negative.binomial")
anova_fit_sqrt <- mvabund::anova.manyglm(fit_sqrt,p.uni = "adjusted")
toc(log=TRUE)

### mvabund based on presence/absence data
tic("Run mvabund model: binary")
dftmp_bin <- dftmp
dftmp_bin[dftmp_bin>0] = 1
mvdftmp_bin <- mvabund(dftmp_bin)
fit_bin <- manyglm(mvdftmp_bin ~ df0_orig_w$SD01_AnaylsisLab, family=binomial())
anova_fit_bin <- mvabund::anova.manyglm(fit_bin,p.uni = "adjusted")
toc(log=TRUE)

## save outputs
saveRDS(fit_raw, file = "model_outputs/fit_raw_model.Rdat")
saveRDS(anova_fit_raw, file = "model_outputs/fit_raw_model_anova.Rdat")
saveRDS(fit_sqrt, file = "model_outputs/fit_sqrt_model.Rdat")
saveRDS(anova_fit_sqrt, file = "model_outputs/fit_sqrt_model_anova.Rdat")
saveRDS(fit_bin, file = "model_outputs/fit_bin_model.Rdat")
saveRDS(anova_fit_bin, file = "model_outputs/fit_bin_model_anova.Rdat")

(x <- unlist(tictoc::tic.log()))
saveRDS(x,file = "model_outputs/model_times.Rdat")
########################################
########################################
########################################

# summary(fit03 <- manyglm(mvdftmp_bin~df0_orig_w$SD01_AnaylsisLab,
#                  family="binomial"))
# anova_fit03 <- mvabund::anova.manyglm(fit03,
#                                       p.uni = "adjusted")
# saveRDS(anova_fit03, file = "data/out/anova_fit03_binomial.Rdat")
anova_fit03 <- readRDS("data/out/anova_fit03_binomial.Rdat")

# which taxa differ?
tx <- as_tibble(t(anova_fit03$uni.p))
tx$nm <- colnames(anova_fit03$uni.p)
tx <- tx[,-1]
txsig <- tx[tx$`df0_orig_w$SD01_AnaylsisLab`<0.06,]

# what proportion significantly differ?
nrow(txsig)/nrow(tx)*100

### prevalence of these by lab
dftmp_bin$Lab <- df0_orig_w$SD01_AnaylsisLab
dftmp_bin %>% 
  group_by(Lab) %>% 
  summarise(across(`Cerataulina pelagica`:`Indet. naked dinoflagellate_>50`,
                   mean, na.rm=TRUE)) %>% 
  t(.) -> mean_bin

mean_bin <- as_tibble(mean_bin)
mean_bin <- mean_bin[-1,]
colnames(mean_bin) <- c(unique(df0_orig_w$SD01_AnaylsisLab)[2],
                        unique(df0_orig_w$SD01_AnaylsisLab)[1])
mean_bin$tx <- tx$nm

### retain only 'significant' taxa
(mean_bin <- mean_bin[mean_bin$tx %in% txsig$nm,])

mean_bin$lab_more <- apply(mean_bin, 1, function(row) {
  colnames(mean_bin)[which.max(row)]
  })
View(mean_bin)

### quick plot
mean_bin %>% 
  pivot_longer(cols=1:2,
               names_to = "Lab",values_to = "Mean") %>% 
  dplyr::select(-lab_more) %>% 
  ggplot(.,aes(x=tx, y= Mean, group=Lab, fill=Lab,shape=Lab))+
  scale_shape_manual(values = c(22,24))+
  scale_fill_manual(values = c("#0072B2","#e79f00"))+
  geom_point(size=4)+
  geom_vline(xintercept = seq(from = .5, to = 20.5, by = 1),
             col= "grey", lty=2)+
  theme(axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  coord_flip()

# to do: ####
# convert to binary data.  Which taxa are more likely to be found/missed by diff labs?
