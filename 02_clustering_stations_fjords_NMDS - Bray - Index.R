
##########################################################
# Load working envt and packages
##########################################################
setwd("C:/Users/loicj/OneDrive - UBC/01_Fjord_eDNA/05_Scratch/NEW")
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(vegan)
library(gclus)
library(ape)
library(ggdendro)
library(dendextend)
library(ggpubr)
library(ggrepel)

sessionInfo()
packageVersion("vegan")

# k value for NMDS
k <- 4

########################
#import data 
########################

# load table
station_df <- read.csv("01_Metadata_site.csv", sep = ",")
species_df <- read.csv("01_Tables/12S_species_table_eDNAIndex_BC-db_2023.csv", sep = ",", row.names = "scientificName")

# remove columns
colnames(species_df)
species_df <- subset(species_df, select = -c(IA03_Summer_100m_2023, IA03_Summer_bottom_2023))# remove this column because no metadata

###############################
# apply filter
###############################
# colnames(species_df)

# # keep only data from the fjord (without historical data)
species_df <- dplyr::select(species_df, contains(c("DB","Phylum","Class","Order","Family","Genus","Species", "2023")))

# then filter by depth
species_df_deep <- dplyr::select(species_df, contains(c("DB","Phylum","Class","Order","Family","Genus","Species",
                                                         "100m", "bottom")))

species_df_surf <- dplyr::select(species_df, contains(c("DB","Phylum","Class","Order","Family","Genus","Species",
                                                        "Surface", "10m")))



################################################
# Subset environmental variables to fit on NMDS
################################################
colnames(station_df)
env <- station_df %>% dplyr::select(Site_ID_Depth_Season, Depth, Temperature, Dissolved_O2, Salinity, Fluorometry)
row.names(env) <- env$Site_ID_Depth_Season

# filter to keep similar number of samples between species and env table
list <- colnames(species_df[7:ncol(species_df)])
list_surf <- colnames(species_df_surf[7:ncol(species_df_surf)])
list_deep <- colnames(species_df_deep[7:ncol(species_df_deep)])

env <- filter(env, Site_ID_Depth_Season %in% list)
env_surf <- filter(env, Site_ID_Depth_Season %in% list_surf)
env_deep <- filter(env, Site_ID_Depth_Season %in% list_deep)

env <- env %>% dplyr::select(Depth, Temperature, Fluorometry, Dissolved_O2, Salinity) # select numerical variables
env_surf <- env_surf %>% dplyr::select(Depth, Temperature, Fluorometry, Dissolved_O2, Salinity) # select numerical variables
env_deep  <- env_deep  %>% dplyr::select(Depth, Temperature, Fluorometry, Dissolved_O2, Salinity) # select numerical variables

# Standardize env values
# suppress differences between scales and units between variables
env.z<-decostand(env, method="standardize",na.rm = T)
env.z_surf<-decostand(env_surf, method="standardize",na.rm = T)
env.z_deep<-decostand(env_deep, method="standardize",na.rm = T)

apply(env.z, 2, mean, na.rm = T) # center (mean~0)...
apply(env.z_surf, 2, mean, na.rm = T) # center (mean~0)...
apply(env.z_deep, 2, mean, na.rm = T) # center (mean~0)...

apply(env.z, 2, sd, na.rm = T)   # reduce (sd=1)
apply(env.z_surf, 2, sd, na.rm = T)   # reduce (sd=1)
apply(env.z_deep, 2, sd, na.rm = T)   # reduce (sd=1)

#########################################
# NMDS with all data
########################################

######################### Prepare and standardize data:
# see: http://www.statsoft.fr/concepts-statistiques/glossaire/c/centrer.html
colnames(species_df)

# subset table to only have numeric columns
df_sp <- species_df[,7:ncol(species_df)]# for index

# turn table
df_sp <- as.data.frame(t(df_sp))

# order table, for vegan so row match between tables
df_sp <- df_sp[order(as.character(row.names(df_sp))), ]
env <- env[order(as.character(row.names(env))), ]

# Hellinger transformation: correct influence of rare species
# gives low weights to variables with low counts (rare species) and many zeros
# ! USE ONLY IF WORKING WITH RELATIVE ABUNDANCE MATRIX (index is already standardized)
# df_hel_sp <-decostand(df_sp, method="hellinger")

#df_hel <- na.omit(df_hel)

# calcul de la matrice de Bray-Curtis
df_BC_sp<-vegdist(df_sp, method="bray", na.rm = F)
#df_BC_sp<-vegdist(df_sp, method="jaccard", na.rm = F, binary = T) # use for presence/absence matrix


######################## Calculate NMDS##################################
# order df_sp and env.z so row.names order match for envfit
df_sp <- df_sp[order(row.names(df_sp)),] # reorder table by order
env.z <- env.z[order(row.names(env.z)),] # reorder table by order

# calcule nMDS
# best goodness fitting is lower values as possible
otus_MDS_sp <- metaMDS(df_sp, distance = "bray", k = k, trymax = 100, na.rm = T)
otus_MDS_sp # give information on MDS

# plot Shepard diagram
# x ==> original distance
# y ==> distance output by a dimension reduction algorithm (here NMDS)
# a perfect dimension reduction will produce a straight line. Impossible, but the closest the better
stressplot(otus_MDS_sp)

# establish relationship between MDS axes and envt variables (dataframe format):
vec.envfit <- envfit(otus_MDS_sp, env.z, na.rm = T, permutations = 999) 
vec.envfit #donne les axes significatifs

#basic plot
plot(otus_MDS_sp)
plot(vec.envfit)

# nb: adjust the multiplier (ordiArrowMul) based on the open plot
vec.envfit.df <- as.data.frame(scores(vec.envfit, "vectors")*ordiArrowMul(vec.envfit))

# Extract NMDS results for ggplot
#build a data frame with NMDS coordinates and metadata
#NMDS_sp = as.data.frame(scores(otus_MDS_sp)) score for other axis
MDS1_sp = otus_MDS_sp$points[,1]
MDS2_sp = otus_MDS_sp$points[,2]
NMDS_sp = data.frame(MDS1 = MDS1_sp, MDS2 = MDS2_sp)
NMDS_sp$Site_ID_Depth_Season <- row.names(NMDS_sp)
row.names(NMDS_sp) <- NMDS_sp$Site_ID_Depth_Season

# Join metadata
NMDS_sp <- left_join(NMDS_sp, station_df, by = "Site_ID_Depth_Season")

# create a column season-depth for NMDS centroids links
NMDS_sp$Season_Depth <- paste(NMDS_sp$Season, NMDS_sp$Depth_feature, sep = "_")
  
#####################################################
# extract centroids for NMDS
#####################################################
# extract NMDS score
scrs <- scores(otus_MDS_sp, display = 'sites')

# compute centroids
cent_Depth <- aggregate(cbind(MDS1, MDS2) ~ Depth, data = NMDS_sp, FUN = mean)
cent_Season <- aggregate(cbind(MDS1, MDS2) ~ Season, data = NMDS_sp, FUN = mean)
cent_Fjord <- aggregate(cbind(MDS1, MDS2) ~ Location.ID, data = NMDS_sp, FUN = mean)
cent_Feature <- aggregate(cbind(MDS1, MDS2) ~ Feature, data = NMDS_sp, FUN = mean)

cent_Season_Depth <- aggregate(cbind(MDS1, MDS2) ~ Season+Depth_feature, data = NMDS_sp, FUN = mean)
cent_Season_Depth$Season_Depth <- paste(cent_Season_Depth$Season, cent_Season_Depth$Depth_feature, sep = "_")
cent_Season_Depth <- cent_Season_Depth[, c("Season_Depth","MDS1", "MDS2")]

# extract segment
segs_Depth <- merge(NMDS_sp, setNames(cent_Depth, c('Depth','oMDS1','oMDS2')),
              by = 'Depth', sort = FALSE)

segs_Season <- merge(NMDS_sp, setNames(cent_Season, c('Season','oMDS1','oMDS2')),
                            by = 'Season', sort = FALSE)

segs_Fjord <- merge(NMDS_sp, setNames(cent_Fjord, c('Location.ID','oMDS1','oMDS2')),
                     by = 'Location.ID', sort = FALSE)

segs_Feature <- merge(NMDS_sp, setNames(cent_Feature, c('Feature','oMDS1','oMDS2')),
                    by = 'Feature', sort = FALSE)

segs_Season_Depth <- merge(NMDS_sp, setNames(cent_Season_Depth, c('Season_Depth','oMDS1','oMDS2')),
                      by = 'Season_Depth', sort = FALSE)

####################### ANOSIM test #############################
# check if NMDS clustering is significant
#https://jkzorz.github.io/2019/06/11/ANOSIM-test.html
# An R value close to "1.0" suggests dissimilarity between groups
# while an R value close to "0" suggests an even distribution of
# high and low ranks within and between groups
#################################################################

ano_depth = anosim(df_sp, NMDS_sp$Feature, distance = "bray", permutations = 999)
ano_fjord = anosim(df_sp, NMDS_sp$Location.ID, distance = "bray", permutations = 999)
ano_season = anosim(df_sp, NMDS_sp$Season, distance = "bray", permutations = 999)

ano_depth
ano_fjord
ano_season

# PERMANOVA with multiple features
ado = adonis2(df_sp~., data=env.z, distance = "bray", permutations = 999, na.action=na.omit)
ado

###################################
# Plot NMDS clustering by categories
###################################

# color by depth
NMDS_plot_depth <- ggplot() + 
   geom_point(data=NMDS_sp,aes(x = MDS1, y = MDS2, col = as.character(Feature)
                               #shape = Location.ID
                               ), 
              size = 3, alpha = 0.9, stroke = 2) +
   labs(title = paste("NMDS"), 
        subtitle = paste("stress:", otus_MDS_sp$stress))+
  scale_color_manual(values=c("#8c6bb1","orange","black","red"))+
   scale_fill_manual(values=c("#8c6bb1","orange","black","red"))+
  # plot the centroids
   geom_segment(data = segs_Feature,
                linetype = "dashed",
                mapping = aes(x = MDS1, y = MDS2, xend = oMDS1, yend = oMDS2,
                              col = as.character(Feature)))+ #link samples to centroids
   geom_label(data = cent_Feature, aes(x = MDS1, y = MDS2, label = Feature, fill = as.character(Feature)),
              size = 3,
              color = "white") + # centroids
  geom_segment(data = vec.envfit.df, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = NULL, colour = "blue")+ # envfit variable
  geom_text_repel(data = vec.envfit.df, aes(x = NMDS1, y = NMDS2, label = row.names(vec.envfit.df)), size = 5, colour = "blue")+#envfit label
  ggtitle(label = "Depth")+
  xlim(-0.9, 1.2)+
  ylim(-1, 1)+
   theme_bw()
NMDS_plot_depth


# create color vectors
season_color <- c("Fall"="#F8766D", 
                  "Spring"="#7CAE00",
                  "Summer"="#007643",
                  "Winter"="#00BFC4",
                  "Others"="grey8") # season color vector

season_color2 <- c("Fall_top_10m"="#F8766D", 
                  "Spring_top_10m"="#7CAE00",
                  "Summer_top_10m"="#007643",
                  "Winter_top_10m"="#00BFC4",
                  "Fall_bottom"="#F8766D", 
                  "Spring_bottom"="#7CAE00",
                  "Summer_bottom"="#007643",
                  "Winter_bottom"="#00BFC4",
                  "Fall"="#F8766D", 
                  "Spring"="#7CAE00",
                  "Summer"="#007643",
                  "Winter"="#00BFC4") # season color vector

# filter centroids and links table to remove 100m later
cent_Season_Depth <- filter(cent_Season_Depth,)
cent_Season_Depth_filter <- cent_Season_Depth %>% filter(grepl("bottom|top", Season_Depth))
segs_Season_Depth_filter <- segs_Season_Depth %>% filter(grepl("bottom|top", Season_Depth))

# plot
NMDS_plot_season <- ggplot() + 
  geom_point(data=NMDS_sp,aes(x = MDS1, y = MDS2, col = Season, shape = Location.ID), 
             size = 3, alpha = 0.9, stroke = 2) +
  labs(title = paste("NMDS"), 
       subtitle = paste("stress:", otus_MDS_sp$stress))+
  scale_shape_manual(values = c(15, 1, 2, 8, 17))+
  # plot the centroids
  geom_segment(data = segs_Season,
               linetype = "dashed",
               mapping = aes(x = MDS1, y = MDS2, xend = oMDS1, yend = oMDS2,
                             col = as.character(Season)))+
  geom_label(data = cent_Season, aes(x = MDS1, y = MDS2, label = Season, fill = as.character(Season)),
             size = 3,
             color = "white") + 
  scale_colour_manual(values=season_color2)+
  scale_fill_manual(values=season_color2)+
  ggtitle(label = "Season")+
  xlim(-1, 1.2)+
  ylim(-1, 1)+
  theme_bw()
NMDS_plot_season

NMDS_plot_season_depth <- ggplot() + 
   geom_point(data=NMDS_sp,aes(x = MDS1, y = MDS2, col = Season, shape = Location.ID), 
              size = 3, alpha = 0.9, stroke = 2) +
   labs(title = paste("NMDS"), 
        subtitle = paste("stress:", otus_MDS_sp$stress))+
  scale_shape_manual(values = c(15, 1, 2, 8, 17))+
   # plot the centroids
   geom_segment(data = segs_Season_Depth_filter,
               linetype = "dashed",
               mapping = aes(x = MDS1, y = MDS2, xend = oMDS1, yend = oMDS2,
                              col = as.character(Season_Depth)))+
   geom_label(data = cent_Season_Depth_filter, aes(x = MDS1, y = MDS2, label = Season_Depth, fill = as.character(Season_Depth)),
              size = 3,
              color = "white") + 
   scale_colour_manual(values=season_color2)+
   scale_fill_manual(values=season_color2)+
   ggtitle(label = "Season")+
  xlim(-1, 1.2)+
  ylim(-1, 1)+
   theme_bw()
NMDS_plot_season_depth

# color by fjord
NMDS_plot_fjord <- ggplot() + 
   geom_point(data=NMDS_sp,aes(x = MDS1, y = MDS2, col = as.character(Location.ID), shape = Feature), 
              size = 2, alpha = 0.9, stroke = 2) +
   #geom_text(data=NMDS_sp,aes(x=MDS1,y=MDS2,label=Site_ID_Depth_Season), size=3, vjust=0, check_overlap = T)+
   labs(title = paste("NMDS"), 
        subtitle = paste("stress:", otus_MDS_sp$stress))+
   scale_shape_manual(values = c(16, 16, 16, 16))+
   scale_color_manual(values=c("#7fc97f", "#beaed4", "#fdc086","#386cb0", "#f0027f"))+
   scale_fill_manual(values=c("#7fc97f", "#beaed4", "#fdc086","#386cb0", "#f0027f"))+
   # plot the centroids
   geom_segment(data = segs_Fjord,
                linetype = "dashed",
                mapping = aes(x = MDS1, y = MDS2, xend = oMDS1, yend = oMDS2,
                              col = as.character(Location.ID)))+
   geom_label(data = cent_Fjord, aes(x = MDS1, y = MDS2, label = Location.ID, fill = as.character(Location.ID)),
              size = 3,
              color = "white") + # centroids
   
   ggtitle(label = "Fjord")+
   theme_bw()
NMDS_plot_fjord


#########################################
# NMDS with surface/10m data
########################################

######################### Prepare and standardize data:

# see: http://www.statsoft.fr/concepts-statistiques/glossaire/c/centrer.html
colnames(species_df_surf)

# subset table to only have numeric columns
df_sp_surf <- species_df_surf[,7:ncol(species_df_surf)]# for index

# turn table
df_sp_surf <- as.data.frame(t(df_sp_surf))

# transformation de Hellinger pour corriger l'influence des especes rares
# gives low weights to variables with low counts (rare species) and many zeros
# ! USE ONLY IF WORKING WITH RELATIVE ABUNDANCE MATRIX (index already standardized)
#df_hel_sp_surf <-decostand(df_sp_surf, method="hellinger")

#df_hel <- na.omit(df_hel)

# calcul de la matrice de distance (choose one)
df_BC_sp_surf<-vegdist(df_sp_surf, method="bray", na.rm = F)#Bray-Curtis
#df_BC_sp_surf<-vegdist(df_sp_surf, method="jaccard", na.rm = F, binary = T)#Jaccard


# order table, for vegan so row match between tables
df_sp_surf <- df_sp_surf[order(as.character(row.names(df_sp_surf))), ]
env_surf <- env[order(as.character(row.names(env_surf))), ]


######################## Calculate NMDS:

# calcul MDS (=PCoA) (unconstrained): only bio
#MDS_sp <- cmdscale(vegdist(df_sp_surf, method = "bray"), k = 2, eig = T, add = T)

# calcule nMDS
# best goodness fitting is lower values as possible
otus_MDS_sp_surf <- metaMDS(df_sp_surf, distance = "bray", k = k, trymax = 100, na.rm = T)

otus_MDS_sp_surf # give information on MDS
stressplot(otus_MDS_sp_surf)

#build a data frame with NMDS coordinates and metadata
MDS1_sp_surf = otus_MDS_sp_surf$points[,1]
MDS2_sp_surf = otus_MDS_sp_surf$points[,2]
NMDS_sp_surf = data.frame(MDS1 = MDS1_sp_surf, MDS2 = MDS2_sp_surf)
NMDS_sp_surf$Site_ID_Depth_Season <- row.names(NMDS_sp_surf)

# Join metadata
NMDS_sp_surf <- left_join(NMDS_sp_surf, station_df, by = "Site_ID_Depth_Season")

row.names(otus_MDS_sp_surf)
# establish relationship between MDS axes and envt variables (dataframe format):
vec.envfit_surf <- envfit(otus_MDS_sp_surf, env.z_surf, na.rm = T, permutations = 999) 
vec.envfit_surf #donne les axes significatifs

plot(otus_MDS_sp_surf)
plot(vec.envfit_surf)

# the multiplier 'ordiArrowMul' based on the open plot, dont forget to plot the NMDS before running that command
vec.envfit.df_surf <- as.data.frame(scores(vec.envfit_surf, "vectors")*ordiArrowMul(vec.envfit_surf))

############### extract centroids for NMDS

# extract NMDS score
scrs <- scores(otus_MDS_sp_surf, display = 'sites')

# compute centroids
cent_Depth_surf <- aggregate(cbind(MDS1, MDS2) ~ Depth, data = NMDS_sp_surf, FUN = mean)
cent_Season_surf <- aggregate(cbind(MDS1, MDS2) ~ Season, data = NMDS_sp_surf, FUN = mean)
cent_Fjord_surf <- aggregate(cbind(MDS1, MDS2) ~ Location.ID, data = NMDS_sp_surf, FUN = mean)
cent_Feature_surf <- aggregate(cbind(MDS1, MDS2) ~ Feature, data = NMDS_sp_surf, FUN = mean)

# extract segment
segs_Depth_surf <- merge(NMDS_sp_surf, setNames(cent_Depth_surf, c('Depth','oMDS1','oMDS2')),
                    by = 'Depth', sort = FALSE)

segs_Season_surf <- merge(NMDS_sp_surf, setNames(cent_Season_surf, c('Season','oMDS1','oMDS2')),
                     by = 'Season', sort = FALSE)

segs_Fjord_surf <- merge(NMDS_sp_surf, setNames(cent_Fjord_surf, c('Location.ID','oMDS1','oMDS2')),
                    by = 'Location.ID', sort = FALSE)

segs_Feature_surf <- merge(NMDS_sp_surf, setNames(cent_Feature_surf, c('Feature','oMDS1','oMDS2')),
                      by = 'Feature', sort = FALSE)


#######################: ANOSIM test
# see if NMDS clustering is significant
#https://jkzorz.github.io/2019/06/11/ANOSIM-test.html

#An R value close to "1.0" suggests dissimilarity between groups
# while an R value close to "0" suggests an even distribution of
# high and low ranks within and between groups
ano_surf_fjord = anosim(df_sp_surf, NMDS_sp_surf$Location.ID, distance = "bray", permutations = 999)
ano_surf_season = anosim(df_sp_surf, NMDS_sp_surf$Season, distance = "bray", permutations = 999)

ano_surf_fjord
ano_surf_season

# PERMANOVA with multiple features
ado = adonis2(df_sp_surf~., data=env_surf, distance = "bray", permutations = 999, na.action=na.omit)
ado


########################## Plot it:
# color by depth
NMDS_plot_depth_surf <- ggplot() + 
  geom_point(data=NMDS_sp_surf,aes(x = MDS1, y = MDS2, col = as.character(Feature), shape = Season), 
             size = 2, alpha = 0.9, stroke = 2) +
  #geom_text(data=NMDS_sp,aes(x=MDS1,y=MDS2,label=Site_ID_Depth_Season), size=3, vjust=0, check_overlap = T)+
  labs(title = paste("NMDS"), 
       subtitle = paste("stress:", otus_MDS_sp_surf$stress))+
  scale_shape_manual(values = c(16, 16, 16, 16))+
  scale_color_manual(values=c("orange","red"))+
  scale_fill_manual(values=c("orange","red"))+
  # plot the centroids
  geom_segment(data = segs_Feature_surf,
               linetype = "dashed",
               mapping = aes(x = MDS1, y = MDS2, xend = oMDS1, yend = oMDS2,
                             col = as.character(Feature)))+
  geom_label(data = cent_Feature_surf, aes(x = MDS1, y = MDS2, label = Feature, fill = as.character(Feature)),
             size = 3,
             color = "white") + # centroids
  ggtitle(label = "Depth / top")+
  theme_bw()
NMDS_plot_depth_surf

# color by season
season_color <- c("Fall"="#F8766D", 
                  "Spring"="#7CAE00",
                  "Summer"="#007643",
                  "Winter"="#00BFC4",
                  "Others"="grey8") # season color vector

NMDS_plot_season_surf <- ggplot() + 
  geom_point(data=NMDS_sp_surf, aes(x = MDS1, y = MDS2,
                                    #shape = Location.ID,
                                    col = Season), 
             size = 4, alpha = 0.7, stroke = 2) +
  labs(title = paste("NMDS"), 
       subtitle = paste("stress:", otus_MDS_sp_surf$stress))+
  # plot the centroids
  geom_segment(data = segs_Season_surf,
               linetype = "dashed",
               mapping = aes(x = MDS1, y = MDS2, xend = oMDS1, yend = oMDS2,
                             col = as.character(Season)))+
  geom_label(data = cent_Season_surf, aes(x = MDS1, y = MDS2, label = Season, fill = as.character(Season)),
             size = 3,
             color = "white") + 
  scale_colour_manual(values=season_color)+
  scale_fill_manual(values=season_color)+
  ggtitle(label = "Top 10m")+
  xlim(-1, 1.6)+
  ylim(-1, 1.1)+
  theme_bw()
NMDS_plot_season_surf

# color by fjord
NMDS_plot_fjord_surf <- ggplot() + 
  geom_point(data=NMDS_sp_surf,aes(x = MDS1, y = MDS2, col = as.character(Location.ID)), 
             size = 4, alpha = 0.7, stroke = 2) +
  #geom_text(data=NMDS_sp,aes(x=MDS1,y=MDS2,label=Site_ID_Depth_Season), size=3, vjust=0, check_overlap = T)+
  labs(title = paste("NMDS"), 
       subtitle = paste("stress:", otus_MDS_sp_surf$stress))+
  #scale_shape_manual(values = c(16, 17, 18, 15))+
  scale_color_manual(values=c("#7fc97f", "#beaed4", "#fdc086","#386cb0", "#f0027f"))+
  scale_fill_manual(values=c("#7fc97f", "#beaed4", "#fdc086","#386cb0", "#f0027f"))+
  # plot the centroids
  geom_segment(data = segs_Fjord_surf,
               linetype = "dashed",
               mapping = aes(x = MDS1, y = MDS2, xend = oMDS1, yend = oMDS2,
                             col = as.character(Location.ID)))+
  geom_label(data = cent_Fjord_surf, aes(x = MDS1, y = MDS2, label = Location.ID, fill = as.character(Location.ID)),
             size = 3,
             color = "white") + # centroids
  ggtitle(label = "Surf")+
  #xlim(-1, 1.6)+
  #ylim(-1, 1.1)+
  theme_bw()
NMDS_plot_fjord_surf

#########################################
# NMDS with bottom
########################################

######################### Prepare and standardize data:

# see: http://www.statsoft.fr/concepts-statistiques/glossaire/c/centrer.html
colnames(species_df_deep)

# subset table to only have numeric columns
df_sp_bt <- species_df_deep[,7:ncol(species_df_deep)]# for index

# turn table
df_sp_bt <- as.data.frame(t(df_sp_bt))

# transformation de Hellinger pour corriger l'influence des especes rares
# gives low weights to variables with low counts (rare species) and many zeros
# ! USE ONLY IF WORKING WITH RELATIVE ABUNDANCE MATRIX (index already standardized)
# df_hel_sp_bt <-decostand(df_sp_bt, method="hellinger")

#df_hel <- na.omit(df_hel)

# calcul de la matrice de Bray-Curtis
df_BC_sp_bt<-vegdist(df_sp_bt, method="bray", na.rm = F)
#df_BC_sp_bt<-vegdist(df_sp_surf, method="jaccard", na.rm = F, binary = T)#Jaccard


######################## Calculate NMDS:

# calcul MDS (=PCoA) (unconstrained): only bio
#MDS_sp <- cmdscale(vegdist(df_sp_bt, method = "bray"), k = 2, eig = T, add = T)

# calcule nMDS
# best goodness fitting is lower values as possible
otus_MDS_sp_bt <- metaMDS(df_sp_bt, distance = "bray", k = k, trymax = 100, na.rm = T)

otus_MDS_sp_bt # give information on MDS
stressplot(otus_MDS_sp_bt)

#build a data frame with NMDS coordinates and metadata
MDS1_sp_bt = otus_MDS_sp_bt$points[,1]
MDS2_sp_bt = otus_MDS_sp_bt$points[,2]
NMDS_sp_bt = data.frame(MDS1 = MDS1_sp_bt, MDS2 = MDS2_sp_bt)
NMDS_sp_bt$Site_ID_Depth_Season <- row.names(NMDS_sp_bt)

# Join metadata
NMDS_sp_bt <- left_join(NMDS_sp_bt, station_df, by = "Site_ID_Depth_Season")


plot(otus_MDS_sp_bt)

############### extract centroids for NMDS

# extract NMDS score
scrs <- scores(otus_MDS_sp_bt, display = 'sites')

# compute centroids
cent_Depth_bt <- aggregate(cbind(MDS1, MDS2) ~ Depth, data = NMDS_sp_bt, FUN = mean)
cent_Season_bt <- aggregate(cbind(MDS1, MDS2) ~ Season, data = NMDS_sp_bt, FUN = mean)
cent_Fjord_bt <- aggregate(cbind(MDS1, MDS2) ~ Location.ID, data = NMDS_sp_bt, FUN = mean)
cent_Feature_bt <- aggregate(cbind(MDS1, MDS2) ~ Feature, data = NMDS_sp_bt, FUN = mean)

# extract segment
segs_Depth_bt <- merge(NMDS_sp_bt, setNames(cent_Depth_bt, c('Depth','oMDS1','oMDS2')),
                         by = 'Depth', sort = FALSE)

segs_Season_bt <- merge(NMDS_sp_bt, setNames(cent_Season_bt, c('Season','oMDS1','oMDS2')),
                          by = 'Season', sort = FALSE)

segs_Fjord_bt <- merge(NMDS_sp_bt, setNames(cent_Fjord_bt, c('Location.ID','oMDS1','oMDS2')),
                         by = 'Location.ID', sort = FALSE)

segs_Feature_bt <- merge(NMDS_sp_bt, setNames(cent_Feature_bt, c('Feature','oMDS1','oMDS2')),
                           by = 'Feature', sort = FALSE)

#######################: ANOSIM test
# see if NMDS clustering is significant
#https://jkzorz.github.io/2019/06/11/ANOSIM-test.html

#An R value close to "1.0" suggests dissimilarity between groups
# while an R value close to "0" suggests an even distribution of
# high and low ranks within and between groups
ano_bt_fjord = anosim(df_sp_bt, NMDS_sp_bt$Location.ID, distance = "bray", permutations = 999)
ano_bt_season = anosim(df_sp_bt, NMDS_sp_bt$Season, distance = "bray", permutations = 999)

ano_bt_fjord
ano_bt_season

########################## Plot it:

# color by depth
NMDS_plot_depth_bt <- ggplot() + 
  geom_point(data=NMDS_sp_bt,aes(x = MDS1, y = MDS2,
                                 #shape = Season,
                                 col = as.character(Feature)), 
             size = 2, alpha = 0.9, stroke = 2) +
  #geom_text(data=NMDS_sp,aes(x=MDS1,y=MDS2,label=Site_ID_Depth_Season), size=3, vjust=0, check_overlap = T)+
  labs(title = paste("NMDS"), 
       subtitle = paste("stress:", otus_MDS_sp_bt$stress))+
  #scale_shape_manual(values = c(16, 16, 16, 16))+
  scale_color_manual(values=c("grey","black"))+
  scale_fill_manual(values=c("grey","black"))+
  # plot the centroids
  geom_segment(data = segs_Feature_bt,
               linetype = "dashed",
               mapping = aes(x = MDS1, y = MDS2, xend = oMDS1, yend = oMDS2,
                             col = as.character(Feature)))+
  geom_label(data = cent_Feature_bt, aes(x = MDS1, y = MDS2, label = Feature, fill = as.character(Feature)),
             size = 3,
             color = "white") + # centroids
  ggtitle(label = "Depth")+
  xlim(-1, 1.3)+
  ylim(-1, 1.3)+
  theme_bw()
NMDS_plot_depth_bt

# color by season
season_color <- c("Fall"="#F8766D", 
                  "Spring"="#7CAE00",
                  "Summer"="#007643",
                  "Winter"="#00BFC4",
                  "Others"="grey8") # season color vector


NMDS_plot_season_bt <- ggplot() + 
  geom_point(data=NMDS_sp_bt,aes(x = MDS1, y = MDS2, 
                                 #shape = Location.ID,
                                 col = Season), 
             size = 4, alpha = 0.7, stroke = 2) +
  labs(title = paste("NMDS"), 
       subtitle = paste("stress:", otus_MDS_sp_bt$stress))+
  scale_color_manual(values=season_color)+
  scale_fill_manual(values=season_color)+
  # plot the centroids
  geom_segment(data = segs_Season_bt,
               linetype = "dashed",
               mapping = aes(x = MDS1, y = MDS2, xend = oMDS1, yend = oMDS2,
                             col = as.character(Season)))+
  geom_label(data = cent_Season_bt, aes(x = MDS1, y = MDS2, label = Season, fill = as.character(Season)),
             size = 3,
             color = "white") + 
  ggtitle(label = "100m/bottom")+
  xlim(-1.2, 1.3)+
  ylim(-1.5, 1.3)+
  theme_bw()
NMDS_plot_season_bt

# color by fjord
fjords_color <- c("Rivers"="#386cb0", 
                  "Bute"="#7fc97f",
                  "Toba"="#f0027f",
                  "Howe_Sound"="#beaed4",
                  "Indian"="#fdc086") # fjords color vector

NMDS_plot_fjord_bt <- ggplot() + 
  geom_point(data=NMDS_sp_bt,aes(x = MDS1, y = MDS2, col = as.character(Location.ID)), 
             size = 4, alpha = 0.7, stroke = 2) +
  #geom_text(data=NMDS_sp,aes(x=MDS1,y=MDS2,label=Site_ID_Depth_Season), size=3, vjust=0, check_overlap = T)+
  labs(title = paste("NMDS"), 
       subtitle = paste("stress:", otus_MDS_sp_bt$stress))+
  #scale_shape_manual(values = c(16, 17, 18, 15))+
  scale_color_manual(values=c("#7fc97f", "#beaed4", "#fdc086","#386cb0", "#f0027f"))+
  scale_fill_manual(values=c("#7fc97f", "#beaed4", "#fdc086","#386cb0", "#f0027f"))+
  # plot the centroids
  geom_segment(data = segs_Fjord_bt,
               linetype = "dashed",
               mapping = aes(x = MDS1, y = MDS2, xend = oMDS1, yend = oMDS2,
                             col = as.character(Location.ID)))+
  geom_label(data = cent_Fjord_bt, aes(x = MDS1, y = MDS2, label = Location.ID, fill = as.character(Location.ID)),
             size = 3,
             color = "white") 
  ggtitle(label = "Deep fjords")+
  xlim(-1.2, 1.3)+
  ylim(-1.5, 1.3)+
  theme_bw()
NMDS_plot_fjord_bt

# Plot them all, season
season <- ggarrange(NMDS_plot_season_surf,
          NMDS_plot_season_bt,
          nrow = 2, ncol = 1,  legend = F)

ggarrange(NMDS_plot_depth, season,
          nrow = 1, ncol = 2,
          align = "v" , common.legend = T)
