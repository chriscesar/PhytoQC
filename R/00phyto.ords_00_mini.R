# 00phyto.ords_00_mini.R
### Side project on phyto data ordination
### this version trims 'big' data set to a single year and a single region

## load packages
# load packages #
ld_pkgs <- c("tidyverse","purrr","readxl", "vegan","seas","ggtext")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

# set themes/global presets
cbPalette <- c("#0072B2","#e79f00","#009E73", "#9ad0f3", "#000000",
               "#D55E00", "#CC79A7", "#F0E442","#0021a0","#305841")
ppi <- 300
perms <- 9999 #number of permutations of analyses
theme_set(ggthemes::theme_few())###set theme for all ggplot objects

# #load data
# df0 <- as_tibble(readxl::read_xlsx("data/Copy of PhytoChl_2000-2020 (WB+RegionalSeas).xlsx",
#                                    sheet = "Abundance2000_2020USE"))
# 
# ## highlight problematic samples which skew outputs
# probsiteIDs <- c(808147)
# 
# # filter, convert to wide to prep for ordination ####
# df0 %>%
#   dplyr::select(.,-TorN,-aphiaID, -Time) %>% #remove unneccessary cols
#   filter(!(`Sample Id` %in% probsiteIDs)) %>% 
#   filter(Year >= 2012) %>%  ##keep only samples gathered from 2012 onwards
#   #drop NA values in TypologyCode to keep only WBs selected for analysis
#   filter(!is.na(TypologyCode)) %>%
#   # mutate(Abundance = 1) %>% ## convert all abundances to Presence/Absence
#   pivot_wider(names_from = TaxonUSE,
#               values_from = Abundance,
#               values_fill=0) -> dfw
# 
# saveRDS(dfw, file = "data/out/phyto_WBs_sideProj_2012_Abund.Rdat")
# write.csv(dfw, file = "data/out/phyto_WBs_sideProj_2012_Abund.csv", row.names = FALSE)
dfw <- readRDS("data/out/phyto_WBs_sideProj_2012_Abund.Rdat")

## assign colours to Typology
col.pal <- c("1C" = cbPalette[1],
             "1T" = cbPalette[2],
             "2C" = cbPalette[3],
             "2T" = cbPalette[4],
             "3T" = cbPalette[5],
             "4C" = cbPalette[6],
             "4T" = cbPalette[7],
             "5C" = cbPalette[8],
             "7C" = cbPalette[9],
             "8C" = cbPalette[10]
             )
sh.pal <- c( "1C" = 21,
             "1T" = 24,
             "2C" = 21,
             "2T" = 24,
             "3T" = 24,
             "4C" = 21,
             "4T" = 24,
             "5C" = 21,
             "7C" = 21,
             "8C" = 21
             )
unique(dfw$`Regional Sea`)

dfw %>% 
  filter(Year == 2019) %>% 
  filter(`Regional Sea` %in% c("Northern North Sea","Southern North Sea")) -> dfwNS2019
  # filter(`Regional Sea` %in% c("Irish Sea")) -> dfwNS2019
  # filter(`Regional Sea` %in% c("Southern North Sea")) -> dfwNS2019
  # filter(`Regional Sea` %in% c("Northern North Sea")) -> dfwNS2019
  
dfwNS2019 <- dfwNS2019 %>% #convert to long format, remove zero values
  pivot_longer(cols = Bacillariaceae:"Ceratium longipes",
               names_to = "taxon", values_to="Abundance") %>% 
  filter(Abundance>0) %>% #remove zeroes
  pivot_wider(names_from = taxon, values_from = Abundance,
              values_fill = 0)

## ordinate data ###
## remove metadata
dfwNS2019 %>%
  dplyr::select(-c(`Regional Sea`:WBArea)) -> dftmp

### pre-processing data ####
dftmptrns <- sqrt(dftmp) ## square root trans
dftmptrns <- wisconsin(dftmptrns)

### RUN ORDINATION ####
start.time <- Sys.time() # start timer
set.seed(pi);ord <- metaMDS(dftmptrns, autotransform = FALSE, trymax = 500)
saveRDS(ord,file = "data/out/phyto_WBs_sideProj_2019_NS_Abund_ORD.Rdat")
end.time <- Sys.time() # stop timer
time.taken <- round(end.time - start.time,2)
time.taken

# ord <- readRDS("data/out/phyto_WBs_sideProj_2019_NS_Abund_ORD.Rdat")
plot(ord)

#### plot in ggplot2 ####
# remove problematics
dfwtrim <- dfwNS2019

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

#### calc centroids ####
## by Typology
scores_site %>% 
  group_by(Typology) %>%
  summarise(mn_ax1_Typology=mean(NMDS1),mn_ax2_Typology=mean(NMDS2)) %>%
  ungroup() -> centr
centr_Typol <- centr
names(centr_Typol)[2] <- "NMDS1";names(centr_Typol)[3] <- "NMDS2"

scores_site <- left_join(scores_site,centr,by="Typology");rm(centr)

## by TypologyCode
scores_site %>% 
  group_by(TypologyCode) %>%
  summarise(mn_ax1_TypologyCode=mean(NMDS1),mn_ax2_TypologyCode=mean(NMDS2)) %>%
  ungroup() -> centr
centr_TypolCode <- centr
names(centr_TypolCode)[2] <- "NMDS1";names(centr_TypolCode)[3] <- "NMDS2"

scores_site <- left_join(scores_site,centr,by="TypologyCode");rm(centr)

## by WB
scores_site %>% 
  group_by(WB) %>%
  summarise(mn_ax1_WB=mean(NMDS1),mn_ax2_WB=mean(NMDS2)) %>%
  ungroup() -> centr
centr_WB <- centr
names(centr_WB)[2] <- "NMDS1";names(centr_WB)[3] <- "NMDS2"

scores_site <- left_join(scores_site,centr,by="WB");rm(centr)

## by Region
scores_site %>% 
  group_by(Region) %>%
  summarise(mn_ax1_Region=mean(NMDS1),mn_ax2_Region=mean(NMDS2)) %>%
  ungroup() -> centr
centr_Region <- centr
names(centr_Region)[2] <- "NMDS1";names(centr_Region)[3] <- "NMDS2"

scores_site <- left_join(scores_site,centr,by="Region");rm(centr)

## by Classification: DIN
scores_site %>% 
  group_by(Class_DIN) %>%
  summarise(mn_ax1_Class_DIN=mean(NMDS1),mn_ax2_Class_DIN=mean(NMDS2)) %>%
  ungroup() -> centr
centr_Class_DIN <- centr
names(centr_Class_DIN)[2] <- "NMDS1";names(centr_Class_DIN)[3] <- "NMDS2"

scores_site <- left_join(scores_site,centr,by="Class_DIN");rm(centr)

## by Classification: Phyto
scores_site %>% 
  group_by(Class_Phyto) %>%
  summarise(mn_ax1_Class_Phyto=mean(NMDS1),mn_ax2_Class_Phyto=mean(NMDS2)) %>%
  ungroup() -> centr
centr_Class_Phyto <- centr
names(centr_Class_Phyto)[2] <- "NMDS1";names(centr_Class_Phyto)[3] <- "NMDS2"

scores_site <- left_join(scores_site,centr,by="Class_Phyto");rm(centr)

## by Classification: OppMac
scores_site %>% 
  group_by(Class_OppMac) %>%
  summarise(mn_ax1_Class_OppMac=mean(NMDS1),mn_ax2_Class_OppMac=mean(NMDS2)) %>%
  ungroup() -> centr
centr_Class_oppMac <- centr
names(centr_Class_oppMac)[2] <- "NMDS1";names(centr_Class_oppMac)[3] <- "NMDS2"

scores_site <- left_join(scores_site,centr,by="Class_OppMac");rm(centr)

#### generate plots! ####
ggplot()+
  geom_hline(colour="grey",yintercept = 0, lty=2)+
  geom_vline(colour="grey",xintercept = 0, lty=2)+
  geom_text(data=scores_species, aes(x = NMDS1, y=NMDS2, label=lb),
            size=3,
            alpha=0.2)+
  geom_segment(data=scores_site,aes(x=NMDS1,y=NMDS2,
                                  colour=Typology,
                                  xend=mn_ax1_Typology,yend=mn_ax2_Typology),
             show.legend = FALSE)+
  geom_point(data=scores_site, show.legend=TRUE,
             aes(x=NMDS1, y=NMDS2,
                 fill = Typology,
                 shape = Typology),
             size=3)+
  scale_fill_manual(values=c(col.pal))+
  scale_colour_manual(values=c(col.pal))+
  # scale_shape_manual(values = rep(c(25:21),each=2))+
  # scale_shape_manual(values = rep(c(24,25,23),each=3))+
  scale_shape_manual(values = sh.pal)+
  geom_textbox(data=centr_Typol,aes(x=NMDS1,y=NMDS2,label=Typology),
               hjust=0.5,width = unit(0.045, "npc"),
               fontface="bold")+
  coord_equal()+
  labs(title="Non-metric Multidimensional Scaling of phyplankton abundances in the North Sea",
       subtitle="Based on data gathered in 2019. Colours & shapes indicate waterbody typology",
       caption=paste0("Stress = ",round(ord$stress,3)))+
  theme(legend.title = element_blank(),
        axis.title = element_text(face="bold"))
