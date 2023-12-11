# 101_init_data_expl.R ####
# initial look at non-labswap data #
## looking at presence-absence data #

# admin ####
# load packages #
ld_pkgs <- c("tidyverse","ggthemes","vegan", "rgl", "mvabund")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

# set themes/global presets
cbPalette <- c("#0072B2","#e79f00","#009E73", "#9ad0f3", "#000000",
                        "#D55E00", "#CC79A7", "#F0E442"
                        )
ppi <- 300
perms <- 9999 #number of permutations of analyses
theme_set(ggthemes::theme_few())###set theme for all ggplot objects

# load data 
df0 <- as_tibble(read.csv("data/out/extracted_data_NonLabSwap.csv"))

### add variable to ID presence-absence of taxa
## NB only 'present' taxa are recorded in the long format data
df0$dens <- 1

# Widen data and run MDS ####
### filter by AnalysisType = "Original" & widen
df0_orig_w <- df0 %>% 
  filter(.,!is.na(SD00FileName)) %>% 
  filter(AnalysisType == "Original") %>% 
  dplyr::select(., # remove unnecessary columns
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
              values_fill = 0) %>% # widen data
  ungroup()
                        
### create temporary file for ordination
# remove metadata
dftmp <- df0_orig_w %>%
  dplyr::select(!c(SD00FileName:AnalysisType))
set.seed(pi)
ord <- vegan::metaMDS(dftmp,
                      try = 20,
                      trymax = 1000,
                      k = 3)
#plot(ord)

# statistical comparisons ####
### do assemblages differ by lab?
## vegan version
(fit01 <-
   vegan::adonis2(dftmp ~ df0_orig_w$SD01_AnaylsisLab,
                  permutations = perms))

## mvabund
mvdftmp <- mvabund(dftmp)
fit02 <-
  manyglm(mvdftmp ~ df0_orig_w$SD01_AnaylsisLab)
# anova_fit02 <- mvabund::anova.manyglm(fit02,p.uni = "adjusted")
# saveRDS(anova_fit02, file = "data/out/anova_fit02_PA.Rdat")

### mvabund based on presence/absence data
dftmp_bin <- dftmp
summary(fit03 <- manyglm(mvdftmp~df0_orig_w$SD01_AnaylsisLab,
                 family="binomial"))
# anova_fit03 <- mvabund::anova.manyglm(fit03,
#                                       p.uni = "adjusted")
# saveRDS(anova_fit03, file = "data/out/anova_fit03_binomial_PA.Rdat")
anova_fit03 <- readRDS("data/out/anova_fit03_binomial_PA.Rdat")

### extract scores
scores_species <-
  as.data.frame(scores(ord, display = "species"))

scores_site <- df0_orig_w %>%
  dplyr::select(c(1:17))

# which taxa differ?
tx <- as_tibble(t(anova_fit03$uni.p))
tx$nm <- rownames(scores_species)
tx <- tx[, -1]
txsig <- tx[tx$`df0_orig_w$SD01_AnaylsisLab` < 0.06, ]

# what proportion significantly differ?
nrow(txsig) / nrow(tx) * 100

### prevalence of these by lab
dftmp_bin$Lab <- df0_orig_w$SD01_AnaylsisLab
dftmp_bin %>%
  group_by(Lab) %>%
  summarise(across(
    `Cerataulina pelagica`:`Indet. naked dinoflagellate_>50`,
    mean,
    na.rm = TRUE
  )) %>%
  t(.) -> mean_bin

mean_bin <- as_tibble(mean_bin)
mean_bin <- mean_bin[-1, ]
colnames(mean_bin) <- c(unique(df0_orig_w$SD01_AnaylsisLab)[2],
                        unique(df0_orig_w$SD01_AnaylsisLab)[1])
mean_bin$tx <- tx$nm

### retain only 'significant' taxa
(mean_bin <- mean_bin[mean_bin$tx %in% txsig$nm, ])

mean_bin$lab_more <- apply(mean_bin, 1, function(row) {
  colnames(mean_bin)[which.max(row)]
})
View(mean_bin)

# extract data for ggplot2 ####
#### extract ordination axes ####
tmp_sites <- as_tibble(as.data.frame(scores(ord, display = "site")))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
scores_site$NMDS1 <- tmp_sites$NMDS1 #location of individual samples in NMDS space
scores_site$NMDS2 <- tmp_sites$NMDS2 #location of individual samples in NMDS space
scores_site$NMDS3 <- tmp_sites$NMDS3 #location of individual samples in NMDS space

#assign colours to variable 'SD01_AnalysisLab'
unq <- levels(as.factor(scores_site$SD01_AnaylsisLab)) ## extract unique Lab values
scores_site$LabCol <- case_when(
  scores_site$SD01_AnaylsisLab == unq[1] ~ cbPalette[1],
  scores_site$SD01_AnaylsisLab == unq[2] ~ cbPalette[2],
  TRUE ~ NA_character_
  )
saveRDS(scores_site, file = "data/out/scores_site3d_PA.Rdata")
rm(tmp_sites)

# Using the scores function from vegan to extract the species scores and
# convert to a data.frame

scores_species$lbfull <-
  row.names(scores_species)
scores_species$lb <-
  make.cepnames(row.names(scores_species))#shorten names
scores_species$lb <-
  gsub("µ", "u", scores_species$lb)

### include only tx in 'significant' list (run later code first!)
scores_species$lb2 <- ifelse(scores_species$lbfull %in% txsig$nm,
                             scores_species$lbfull, " ")
scores_species$lb2 <-
  gsub("µ", "u", scores_species$lb2)
scores_species$lb2 <-
  gsub("≤", "<", scores_species$lb2)
scores_species$lb2 <-
  gsub("≥", ">", scores_species$lb2)

saveRDS(scores_species, file = "data/out/scores_species3d_PA.Rdata")

plot3d(
  x = scores_site$NMDS1,
  y = scores_site$NMDS2,
  z = scores_site$NMDS3,
  col = scores_site$LabCol,
  type = "s",
  size = 1,
  xlab = "NMDS1",
  ylab = "NMDS2",
  zlab = "NMDS3"
)

text3d(
  x = scores_species$NMDS1,
  y = scores_species$NMDS2,
  z = scores_species$NMDS3,
  scores_species$lb,#ALL taxa
  # scores_species$lb2,#'SIG' taxa only
  size = 1,
  cex = .75,
  col = "grey"
)

text3d(
  x = scores_species$NMDS1,
  y = scores_species$NMDS2,
  z = scores_species$NMDS3,
  scores_species$lb2,#'SIG' taxa only
  size = 12,
  cex = .9,
  col = "red"
)

axes3d()
title3d(
  xlab = "NMDS1",
  ylab = "NMDS2",
  zlab = "NMDS3",
  font = 2
)

### quick plot
mean_bin %>%
  pivot_longer(cols = 1:2,
               names_to = "Lab",
               values_to = "Mean") %>%
  group_by(tx) %>% 
  mutate(Max = max(Mean),
         Min = min(Mean)) %>% 
  ungroup() %>% 
  mutate(Mean = as.numeric(Mean)) %>% 
  dplyr::select(-lab_more) %>%
  arrange(desc(Mean),Lab) %>% 
  ggplot(., aes(
    # x = tx,
    x = reorder(tx,Mean),
    y = Mean,
    group = Lab,
    fill = Lab,
    shape = Lab
  )) +
  ylab("Mean prevalence")+
  ylim(0,1)+
  scale_shape_manual(values = c(22, 24)) +
  scale_fill_manual(values = c("#0072B2", "#e79f00")) +
  geom_point(size = 4) +
  geom_vline(
    xintercept = seq(from = .5, to = 20.5, by = 1),
    col = "grey",
    lty = 2
  ) +
  theme(
    # axis.text.x = element_blank(),
    #axis.title.x = element_blank(),
    axis.title.y = element_blank()
        ) +
  coord_flip()
