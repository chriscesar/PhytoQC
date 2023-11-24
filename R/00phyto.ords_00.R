# 00phyto.ords_00.R
### Side project on phyto data ordination

## load packages
# load packages #
ld_pkgs <- c("tidyverse","purrr","readxl", "vegan","seas")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

# set themes/global presets
cbPalette <- c("#0072B2","#e79f00","#009E73", "#9ad0f3", "#000000",
               "#D55E00", "#CC79A7", "#F0E442")
ppi <- 300
perms <- 9999 #number of permutations of analyses
theme_set(ggthemes::theme_few())###set theme for all ggplot objects

#load data
df0 <- as_tibble(readxl::read_xlsx("data/Copy of PhytoChl_2000-2020 (WB+RegionalSeas).xlsx",
                                   sheet = "Abundance2000_2020USE"))

## highlight problematic samples which skew outputs
probsiteIDs <- c(808147)

# filter, convert to wide to prep for ordination ####
df0 %>%
  dplyr::select(.,-TorN,-aphiaID, -Time) %>% #remove unneccessary cols
  filter(!(`Sample Id` %in% probsiteIDs)) %>% 
  filter(Year >= 2012) %>%  ##keep only samples gathered from 2012 onwards
  #drop NA values in TypologyCode to keep only WBs selected for analysis
  filter(!is.na(TypologyCode)) %>%
  # mutate(Abundance = 1) %>% ## convert all abundances to Presence/Absence
  pivot_wider(names_from = TaxonUSE,
              values_from = Abundance,
              values_fill=0) -> dfw

saveRDS(dfw, file = "data/out/phyto_WBs_sideProj_2012_Abund.Rdat")
# write.csv(dfw, file = "data/out/phyto_WBs_sideProj_2012_Abund.csv", row.names = FALSE)
# dfw <- readRDS("data/out/phyto_WBs_sideProj_2019_Bin.Rdat")

## remove singletons ####
# removes samples which only contain a single taxon
dfwtrim <- dfw %>% 
  filter(rowSums(dplyr::select(., -c(`Regional Sea`:WBArea)) != 0, na.rm = TRUE) > 1)

## ordinate data ###
## remove metadata
dfwtrim %>% 
  dplyr::select(-c(`Regional Sea`:WBArea)) -> dftmp

### pre-processing data ####
dftmptrns <- sqrt(dftmp) ## square root trans
dftmptrns <- wisconsin(dftmptrns)

### RUN ORDINATION ####
start.time <- Sys.time() # start timer
set.seed(pi);ord <- metaMDS(dftmptrns, autotransform = FALSE)
saveRDS(ord,file = "data/out/phyto_WBs_sideProj_2012_Abund_ORD.Rdat")
plot(ord)

time.taken <- round(end.time - start.time,2)
time.taken

#### plot in ggplot2 ####
scores_site <- as_tibble(as.data.frame(scores(ord,display = "site")))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
scores_site$SampID <- dfwtrim$`Sample Id`
scores_site$year <- dfwtrim$Year
scores_site$date <- dfwtrim$Date
scores_site$DJF <- as.factor(mkseas(dfwtrim$Date, width="DJF"))#convert dates to 3month seasonal block
scores_site$Region <- dfwtrim$`Regional Sea`
scores_site$WB <- dfwtrim$Waterbody
scores_site$CW_TW <- dfwtrim$`CW/TW`
scores_site$Typology <- dfwtrim$Typology
scores_site$TypologyCode <- dfwtrim$TypologyCode
scores_site$WBType <- dfwtrim$WBType
scores_site$WBArea <- dfwtrim$WBArea
scores_site$Class_DIN <- dfwtrim$`Nut status`
scores_site$Class_Phyto <- dfwtrim$`Phyto Status`
scores_site$Class_OppMac <- dfwtrim$`Opp Mac`

# Using the scores function from vegan to extract the species scores and
# convert to a data.frame
scores_species <- as.data.frame(scores(ord,display = "species"))
scores_species$lb <-  make.cepnames(row.names(scores_species))#shorten names

#### generate mean centroids ####
## by Typology
scores_site %>% 
  group_by(Typology) %>%
  summarise(mn_ax1_Typology=mean(NMDS1),mn_ax2_Typology=mean(NMDS2)) %>%
  ungroup() -> centr

scores_site <- left_join(scores_site,centr,by="Typology");rm(centr)

## by TypologyCode
scores_site %>% 
  group_by(TypologyCode) %>%
  summarise(mn_ax1_TypologyCode=mean(NMDS1),mn_ax2_TypologyCode=mean(NMDS2)) %>%
  ungroup() -> centr

scores_site <- left_join(scores_site,centr,by="TypologyCode");rm(centr)

## by WB
scores_site %>% 
  group_by(WB) %>%
  summarise(mn_ax1_WB=mean(NMDS1),mn_ax2_WB=mean(NMDS2)) %>%
  ungroup() -> centr

scores_site <- left_join(scores_site,centr,by="WB");rm(centr)

## by Region
scores_site %>% 
  group_by(Region) %>%
  summarise(mn_ax1_Region=mean(NMDS1),mn_ax2_Region=mean(NMDS2)) %>%
  ungroup() -> centr

scores_site <- left_join(scores_site,centr,by="Region");rm(centr)

## by Classification: DIN
scores_site %>% 
  group_by(Class_DIN) %>%
  summarise(mn_ax1_Class_DIN=mean(NMDS1),mn_ax2_Class_DIN=mean(NMDS2)) %>%
  ungroup() -> centr

scores_site <- left_join(scores_site,centr,by="Class_DIN");rm(centr)

## by Classification: Phyto
scores_site %>% 
  group_by(Class_Phyto) %>%
  summarise(mn_ax1_Class_Phyto=mean(NMDS1),mn_ax2_Class_Phyto=mean(NMDS2)) %>%
  ungroup() -> centr

scores_site <- left_join(scores_site,centr,by="Class_Phyto");rm(centr)

## by Classification: OppMac
scores_site %>% 
  group_by(Class_OppMac) %>%
  summarise(mn_ax1_Class_OppMac=mean(NMDS1),mn_ax2_Class_OppMac=mean(NMDS2)) %>%
  ungroup() -> centr

scores_site <- left_join(scores_site,centr,by="Class_OppMac");rm(centr)

#### generate plots! ####
ggplot()+
  geom_hline(colour="grey",yintercept = 0, lty=2)+
  geom_vline(colour="grey",xintercept = 0, lty=2)+
  geom_text(data=scores_species, aes(x = NMDS1, y=NMDS2, label=lb),
            size=3,
            alpha=0.2)+
  geom_segment(data=scores_site,aes(x=NMDS1,y=NMDS2,
                                  colour=Region,
                                  xend=mn_ax1_Region,yend=mn_ax2_Region),
             show.legend = FALSE)+
  geom_point(data=scores_site, show.legend=TRUE,
             aes(x=NMDS1, y=NMDS2,
                 fill = Region,
                 shape = Region),
             size=3)+
  scale_fill_manual(values=c(cbPalette))+
  scale_colour_manual(values=c(cbPalette))+
  # scale_shape_manual(values = rep(c(25:21),each=2))+
  scale_shape_manual(values = rep(c(24,25,23),each=2))+
  # coord_equal()+
  # labs(title="Non-metric Multidimensional Scaling of zooplankton taxon abundances",
  #      subtitle="Colours & shapes indicate region",
  #      caption=paste0("Stress = ",round(ord$stress,3),"\nSamples gathered between ",
  #                     min(dfw$sample.date)," & ",max(dfw$sample.date)))+
  theme(legend.title = element_blank(),
        axis.title = element_text(face="bold"))
