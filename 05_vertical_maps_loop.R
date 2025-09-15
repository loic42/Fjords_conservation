setwd("C:/Users/loicj/OneDrive - UBC/01_Fjord_eDNA/05_Scratch/NEW")
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(vegan)
library(ggbeeswarm)


########################
#import data 
########################

# load table
station_df <- read.csv("01_Metadata_site.csv", sep = ",")
species_df <- read.csv("01_Tables/12S_species_table_eDNAIndex_BC-db_2023.csv")#12S or COISal
#species_df <- read.csv("01_Tables/COISal_species_table_eDNAIndex_fjords.csv")

# make a tidy table
species_df_long <- species_df %>% pivot_longer(-c(Phylum,Class,Order,Family,Genus,Species, scientificName), names_to = "Site_ID_Depth_Season", values_to = "count_reads") # turn reads df to tidy

# merge metadata
species_df_merge <- left_join(species_df_long, station_df, by = "Site_ID_Depth_Season")

# reorder table by depth and Season
species_df_merge$Feature <- factor(species_df_merge$Feature, levels=c("bottom", "100m", "10m", "Surface"))
species_df_merge$Season <- factor(species_df_merge$Season, levels=c("Winter","Spring","Summer", "Fall"))

# filter by Fjord
species_df_merge_IA <- filter(species_df_merge, Location.ID == "Indian")
species_df_merge_HS <- filter(species_df_merge, Location.ID == "Howe_Sound")
species_df_merge_TO <- filter(species_df_merge, Location.ID == "Toba")
species_df_merge_BU <- filter(species_df_merge, Location.ID == "Bute")
species_df_merge_RI <- filter(species_df_merge, Location.ID == "Rivers")

# reorder IA and HS so it goes same "flowing" order than the other fjords
species_df_merge_IA$Site_ID <- factor(species_df_merge_IA$Site_ID, levels=c("IA03", "IA02", "IA01"))
species_df_merge_HS$Site_ID <- factor(species_df_merge_HS$Site_ID, levels=c("HS03","HSGLASSR1", "HS02","HS01","HSPLANT"))

#########################
# Select variables
########################

# color by season
season_color <- c("Fall"="#F8766D", 
                  "Spring"="#7CAE00",
                  "Summer"="#007643",
                  "Winter"="#00BFC4",
                  "Others"="grey8") # season color vector

############################################
# Make a loop function to generate multiple plots
############################################
library(purrr)

# Define a function to generate the plot
make_species_plot <- function(data, region_name, species_name) {
  ggplot(data = filter(data, scientificName == species_name),
         aes(x = Site_ID, y = Feature, shape = as.character(Year))) +
    geom_point(aes(size = count_reads), col = "white", alpha = 1) +# to keep the y axis even when empty
    geom_beeswarm(data = ~ dplyr::filter(.x, count_reads > 0.001),
                  aes(size = count_reads, fill = Season),
                  alpha = 0.8, cex = 3) +
    scale_fill_manual(values = season_color) +
    scale_size_area(max_size = 10, limits = c(0, 1), name = "eDNA Index Value") +
    scale_shape_manual(values = c("2019" = 24, "2022" = 22,"2023" =  21)) +
    labs(x = "Latitude", y = "depth (m)") +
    ggtitle(region_name, subtitle = species_name) +
    theme_bw()
}

# Loop over species and regions
species_list <- unique(c(
  species_df_merge_IA$scientificName,
  species_df_merge_HS$scientificName,
  species_df_merge_BU$scientificName,
  species_df_merge_TO$scientificName,
  species_df_merge_RI$scientificName
))

region_dfs <- list(
  "Indian" = species_df_merge_IA,
  "Howe Sound" = species_df_merge_HS,
  "Bute" = species_df_merge_BU,
  "Toba" = species_df_merge_TO,
  "Rivers" = species_df_merge_RI
)

# Create a nested list of plots: region -> species -> plot
all_plots <- imap(region_dfs, function(region_df, region_name) {
  map(species_list, function(sp) {
    make_species_plot(region_df, region_name, sp)
  }) %>% set_names(species_list)
})

# here is how to call a plot
all_plots$`Howe Sound`$`Salmo salar`

####################################################################
# plot them all
####################################################################
library(ggpubr)

# List CPs species
# forage fishes
sp1 <- "Ammodytes personatus"
sp2 <- "Clupea pallasii"
sp3 <- "Engraulis mordax"
sp4 <- "Cymatogaster aggregata"
sp5 <- "Mallotus villosus"
sp7 <- "Thaleichthys pacificus"

# demersal/"deep-sea" fish
sp9 <- "Anarrhichthys ocellatus"
sp10 <- "Anoplopoma fimbria"
sp11 <- "Cryptacanthodes giganteus"
sp12 <- "Microstomus pacificus"
sp13 <- "Ophiodon elongatus"
sp14 <- "Stenobrachius leucopsarus"
sp15 <- "Leuroglossus schmidti"
sp6 <- "Merluccius productus"
sp8 <- "Gadus chalcogrammus"

# flatfish
sp16 <- "Myoxocephalus polyacanthocephalus"
sp17 <- "Platichthys stellatus"
sp18 <- "Pleuronectes quadrituberculatus"
sp19 <- "Atheresthes evermanni"

# coastal/tidal species
sp20 <- "Leptoclinus maculatus"
sp21 <- "Hemilepidotus jordani"
sp22 <- "Sebastolobus macrochir"

# shark and Marine mammals
sp23 <- "Hexanchus griseus"
sp24 <- "Caliraja rhina"
sp25 <- "Squalus suckleyi"
sp26 <- "Lagenorhynchus obliquidens"
sp27 <- "Megaptera novaeangliae"
sp28 <- "Balaenoptera musculus"
sp29 <- "Phoca vitulina"
sp30 <- "Phocoena phocoena"

# salmonids
sp31 <- "Oncorhynchus gorbuscha"
sp32 <- "Oncorhynchus keta"
sp33 <- "Oncorhynchus kisutch"
sp34 <- "Oncorhynchus mykiss"
sp35 <- "Oncorhynchus nerka"
sp36 <- "Oncorhynchus tshawytscha"
sp37 <- "Salmo salar"
sp38 <- "Salvelinus confluentus"
sp39 <- "Salvelinus malma"

all <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9, sp10, sp11, sp12,
         sp13, sp14, sp15, sp16, sp17, sp18, sp19, sp20, sp21, sp22, sp23,
         sp24, sp25, sp26, sp27, sp28, sp29, sp30, sp31, sp32, sp33, sp34,
         sp35, sp36, sp37, sp38, sp39)

################################
# single plot
################################
library(ggpubr)

ggarrange(all_plots$`Rivers`$`Oncorhynchus kisutch`, all_plots$`Rivers`$`Ammodytes personatus`,
          nrow = 2, ncol = 1,
          common.legend = T)

ggarrange(all_plots$`Bute`$`Gadus chalcogrammus`, all_plots$`Bute`$`Merluccius productus`,
          all_plots$`Bute`$`Squalus suckleyi`, all_plots$`Bute`$`Hexanchus griseus`,
          nrow = 2, ncol = 2,
          common.legend = T)


ggarrange(all_plots$`Howe Sound`$`Gadus chalcogrammus`, all_plots$`Howe Sound`$`Merluccius productus`,
          all_plots$`Howe Sound`$`Leuroglossus schmidti`, all_plots$`Howe Sound`$`Anoplopoma fimbria`,
          all_plots$`Howe Sound`$`Hexanchus griseus`, all_plots$`Howe Sound`$`Squalus suckleyi`,
          nrow = 3, ncol = 2,
          common.legend = T)

############################
# Plot and export 15x18
############################
library(ggpubr)

# vertical salmons
ggarrange(all_plots$`Rivers`$`Oncorhynchus gorbuscha`,all_plots$`Rivers`$`Oncorhynchus keta`,
          all_plots$`Rivers`$`Oncorhynchus kisutch`,all_plots$`Rivers`$`Oncorhynchus tshawytscha`,
          all_plots$`Rivers`$`Oncorhynchus nerka`,all_plots$`Rivers`$`Oncorhynchus mykiss`,
          all_plots$`Bute`$`Oncorhynchus gorbuscha`,all_plots$`Bute`$`Oncorhynchus keta`,
          all_plots$`Bute`$`Oncorhynchus kisutch`,all_plots$`Bute`$`Oncorhynchus tshawytscha`,
          all_plots$`Bute`$`Oncorhynchus nerka`,all_plots$`Bute`$`Oncorhynchus mykiss`,
          all_plots$`Toba`$`Oncorhynchus gorbuscha`,all_plots$`Toba`$`Oncorhynchus keta`,
          all_plots$`Toba`$`Oncorhynchus kisutch`,all_plots$`Toba`$`Oncorhynchus tshawytscha`,
          all_plots$`Toba`$`Oncorhynchus nerka`,all_plots$`Toba`$`Oncorhynchus mykiss`,
          all_plots$`Howe Sound`$`Oncorhynchus gorbuscha`,all_plots$`Howe Sound`$`Oncorhynchus keta`,
          all_plots$`Howe Sound`$`Oncorhynchus kisutch`,all_plots$`Howe Sound`$`Oncorhynchus tshawytscha`,
          all_plots$`Howe Sound`$`Oncorhynchus nerka`,all_plots$`Howe Sound`$`Oncorhynchus mykiss`,
          all_plots$`Indian`$`Oncorhynchus gorbuscha`,all_plots$`Indian`$`Oncorhynchus keta`,
          all_plots$`Indian`$`Oncorhynchus kisutch`,all_plots$`Indian`$`Oncorhynchus tshawytscha`,
          all_plots$`Indian`$`Oncorhynchus nerka`,all_plots$`Indian`$`Oncorhynchus mykiss`,
          nrow = 5, ncol = 6,
          common.legend = T
)

# Figure: vertical forage fish export 15x15
ggarrange(all_plots$`Rivers`$`Ammodytes personatus`,all_plots$`Rivers`$`Clupea pallasii`,
          all_plots$`Rivers`$`Engraulis mordax`, all_plots$`Rivers`$`Thaleichthys pacificus`,
          all_plots$`Bute`$`Ammodytes personatus`,all_plots$`Bute`$`Clupea pallasii`,
          all_plots$`Bute`$`Engraulis mordax`,all_plots$`Bute`$`Thaleichthys pacificus`,
          all_plots$`Toba`$`Ammodytes personatus`,all_plots$`Toba`$`Clupea pallasii`,
          all_plots$`Toba`$`Engraulis mordax`,all_plots$`Toba`$`Thaleichthys pacificus`,
          all_plots$`Howe Sound`$`Ammodytes personatus`,all_plots$`Howe Sound`$`Clupea pallasii`,
          all_plots$`Howe Sound`$`Engraulis mordax`,all_plots$`Howe Sound`$`Thaleichthys pacificus`,
          all_plots$`Indian`$`Ammodytes personatus`,all_plots$`Indian`$`Clupea pallasii`,
          all_plots$`Indian`$`Engraulis mordax`,all_plots$`Indian`$`Thaleichthys pacificus`,
          nrow = 5, ncol = 4,
          common.legend = T
)

# Figure: vertical mesopelagic fish export 15x15
ggarrange(all_plots$`Rivers`$`Stenobrachius leucopsarus`,all_plots$`Rivers`$`Gadus chalcogrammus`,
          all_plots$`Rivers`$`Merluccius productus`, all_plots$`Rivers`$`Squalus suckleyi`,
          all_plots$`Bute`$`Stenobrachius leucopsarus`,all_plots$`Bute`$`Gadus chalcogrammus`,
          all_plots$`Bute`$`Merluccius productus`,all_plots$`Bute`$`Squalus suckleyi`,
          all_plots$`Toba`$`Stenobrachius leucopsarus`,all_plots$`Toba`$`Gadus chalcogrammus`,
          all_plots$`Toba`$`Merluccius productus`,all_plots$`Toba`$`Squalus suckleyi`,
          all_plots$`Howe Sound`$`Stenobrachius leucopsarus`,all_plots$`Howe Sound`$`Gadus chalcogrammus`,
          all_plots$`Howe Sound`$`Merluccius productus`,all_plots$`Howe Sound`$`Squalus suckleyi`,
          all_plots$`Indian`$`Stenobrachius leucopsarus`,all_plots$`Indian`$`Gadus chalcogrammus`,
          all_plots$`Indian`$`Merluccius productus`,all_plots$`Indian`$`Squalus suckleyi`,
          nrow = 5, ncol = 4,
          common.legend = T
)

# vertical  fish
ggarrange(all_plots$`Rivers`$`Gadus chalcogrammus`,all_plots$`Rivers`$`Merluccius productus`,
          all_plots$`Rivers`$`Leuroglossus schmidti`,all_plots$`Rivers`$`Anoplopoma fimbria`,
          all_plots$`Rivers`$`Ophiodon elongatus`,all_plots$`Rivers`$`Cryptacanthodes giganteus`,
          all_plots$`Bute`$`Gadus chalcogrammus`,all_plots$`Bute`$`Merluccius productus`,
          all_plots$`Bute`$`Leuroglossus schmidti`,all_plots$`Bute`$`Anoplopoma fimbria`,
          all_plots$`Bute`$`Ophiodon elongatus`,all_plots$`Bute`$`Cryptacanthodes giganteus`,
          all_plots$`Toba`$`Gadus chalcogrammus`,all_plots$`Toba`$`Merluccius productus`,
          all_plots$`Toba`$`Leuroglossus schmidti`,all_plots$`Toba`$`Anoplopoma fimbria`,
          all_plots$`Toba`$`Ophiodon elongatus`,all_plots$`Toba`$`Cryptacanthodes giganteus`,
          all_plots$`Howe Sound`$`Gadus chalcogrammus`,all_plots$`Howe Sound`$`Merluccius productus`,
          all_plots$`Howe Sound`$`Leuroglossus schmidti`,all_plots$`Howe Sound`$`Anoplopoma fimbria`,
          all_plots$`Howe Sound`$`Ophiodon elongatus`,all_plots$`Howe Sound`$`Cryptacanthodes giganteus`,
          all_plots$`Indian`$`Gadus chalcogrammus`,all_plots$`Indian`$`Merluccius productus`,
          all_plots$`Indian`$`Leuroglossus schmidti`,all_plots$`Indian`$`Anoplopoma fimbria`,
          all_plots$`Indian`$`Ophiodon elongatus`,all_plots$`Indian`$`Cryptacanthodes giganteus`,
          nrow = 5, ncol = 6,
          common.legend = T
)


ggarrange(all_plots$`Howe Sound`$`Oncorhynchus gorbuscha`, all_plots$`Indian`$`Oncorhynchus gorbuscha`,
          all_plots$`Howe Sound`$`Oncorhynchus keta`, all_plots$`Indian`$`Oncorhynchus keta`,
          all_plots$`Howe Sound`$`Oncorhynchus kisutch`, all_plots$`Indian`$`Oncorhynchus kisutch`,
           all_plots$`Howe Sound`$`Oncorhynchus tshawytscha`, all_plots$`Indian`$`Oncorhynchus tshawytscha`,
          all_plots$`Howe Sound`$`Oncorhynchus nerka`, all_plots$`Indian`$`Oncorhynchus nerka`,
          nrow = 5, ncol = 2,
          common.legend = T
)


#######################################################
# Cumulated index calculation
#######################################################
forage_fish_list <- c("Ammodytes personatus","Clupea pallasii",
                      "Engraulis mordax","Cymatogaster aggregata",
                      "Mallotus villosus","Thaleichthys pacificus")

# filter table
df <- filter(species_df_merge, scientificName %in% all)
unique(df$scientificName)

# subset table
colnames(df)
df <- df %>% select(Site_ID,Feature, count_reads, Location.ID, Season)

# aggregate edna index by location and depth
df_sum <- aggregate(df$count_reads, by = list(df$Site_ID, df$Feature, df$Location.ID, df$Season), FUN = sum)
colnames(df_sum) <- c("Site_ID", "Feature", "Location.ID", "Season", "sum_index")

# Define a function to generate the plot
cum_RI_plot <- ggplot(data = filter(df_sum, Location.ID == "Rivers"),
       aes(x = Site_ID, y = Feature)) +
  geom_point(aes(size = sum_index), col = "white", alpha = 1) +# to keep the y axis even when empty
  geom_beeswarm(data = ~ dplyr::filter(.x, sum_index > 0.001),
                aes(size = sum_index, col = Season),
                alpha = 0.8, cex = 3) +
  scale_color_manual(values = season_color) +
  scale_size_area(max_size = 10,
                  limits = c(0, 3),
                  name = "Cumulated eDNA Index Value") +
  labs(x = "Latitude", y = "depth (m)") +
  ggtitle("Rivers") +
  theme_bw()

cum_BU_plot <- ggplot(data = filter(df_sum, Location.ID == "Bute"),
         aes(x = Site_ID, y = Feature)) +
    geom_point(aes(size = sum_index), col = "white", alpha = 1) +# to keep the y axis even when empty
    geom_beeswarm(data = ~ dplyr::filter(.x, sum_index > 0.001),
                  aes(size = sum_index, col = Season),
                  alpha = 0.8, cex = 3) +
  scale_color_manual(values = season_color) +
    scale_size_area(max_size = 10,
                    limits = c(0, 3),
                    name = "Cumulated eDNA Index Value") +
    labs(x = "Latitude", y = "depth (m)") +
    ggtitle("Bute") +
    theme_bw()

cum_TO_plot <- ggplot(data = filter(df_sum, Location.ID == "Toba"),
       aes(x = Site_ID, y = Feature)) +
  geom_point(aes(size = sum_index), col = "white", alpha = 1) +# to keep the y axis even when empty
  geom_beeswarm(data = ~ dplyr::filter(.x, sum_index > 0.001),
                aes(size = sum_index, col = Season),
                alpha = 0.8, cex = 3) +
  scale_color_manual(values = season_color) +
  scale_size_area(max_size = 10,
                  limits = c(0, 3),
                  name = "Cumulated eDNA Index Value") +
  labs(x = "Latitude", y = "depth (m)") +
  ggtitle("Toba") +
  theme_bw()

cum_HS_plot <- ggplot(data = filter(df_sum, Location.ID == "Howe_Sound"),
       aes(x = Site_ID, y = Feature)) +
  geom_point(aes(size = sum_index), col = "white", alpha = 1) +# to keep the y axis even when empty
  geom_beeswarm(data = ~ dplyr::filter(.x, sum_index > 0.001),
                aes(size = sum_index, col = Season),
                alpha = 0.8, cex = 3) +
  scale_color_manual(values = season_color) +
  scale_size_area(max_size = 10,
                  limits = c(0, 3),
                  name = "Cumulated eDNA Index Value") +
  labs(x = "Latitude", y = "depth (m)") +
  ggtitle("Howe Sound") +
  theme_bw()


cum_IA_plot <- ggplot(data = filter(df_sum, Location.ID == "Indian"),
       aes(x = Site_ID, y = Feature)) +
  geom_point(aes(size = sum_index), col = "white", alpha = 1) +# to keep the y axis even when empty
  geom_beeswarm(data = ~ dplyr::filter(.x, sum_index > 0.001),
                aes(size = sum_index, col = Season),
                alpha = 0.8, cex = 3) +
  scale_color_manual(values = season_color) +
  scale_size_area(max_size = 10,
                  limits = c(0, 3),
                  name = "Cumulated eDNA Index Value") +
  labs(x = "Latitude", y = "depth (m)") +
  ggtitle("Indian") +
  theme_bw()

library(ggpubr)
ggarrange(cum_RI_plot, cum_BU_plot, cum_TO_plot, cum_HS_plot, cum_IA_plot,
          ncol=2, nrow = 3, common.legend = T)
