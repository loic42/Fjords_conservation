library(dplyr)
library(tidyverse)
library(ggplot2)

setwd("C:/Users/loicj/OneDrive - UBC/01_Fjord_eDNA/05_Scratch/NEW")

# load table
df <- read.csv("01_Tables/12S_species_table_eDNAIndex_BC-db_2023.csv", sep = ",")
df <- select(df, contains(c("scientificName", "2023")))
metadata <- read.csv("01_Metadata_site.csv", sep = ",")
taxo <- read.csv("01_taxo_all_genus.csv", sep = ",")
df_CPs <- read.csv("02_Conservation_priorities/01_conservation_priorities_final_table.csv")

# filter table to contains only CPs species
df_CPs <- filter(df_CPs, Final_score == 2)
list_CPs <- as.data.frame(df_CPs$scientificName)
colnames(list_CPs) <- c("scientificName")

# Select CPs species only
df <- filter(df, scientificName %in% list_CPs$scientificName)

# pivot table
df_long <- df %>% pivot_longer(cols = 2:ncol(df) , names_to = "Site_ID_Depth_Season", values_to = "count") # turn APC df to tidy

# add a column presence-absence
df_long$count_pres <- ifelse(df_long$count > 0, 1, 0)

#merge metadata
df_long <- left_join(df_long, metadata, by = "Site_ID_Depth_Season")

# aggregate table (mean relative abundance) by fjords and species
df_long_merged_fjord <- aggregate(df_long$count, by = list(df_long$scientificName,
                                                     df_long$Location.ID), FUN = mean)

df_long_merged_season <- aggregate(df_long$count, by = list(df_long$scientificName,
                                                           df_long$Season), FUN = mean)

df_long_merged_depth <- aggregate(df_long$count, by = list(df_long$scientificName,
                                                           df_long$Feature), FUN = mean)

# For the depth table, keep only depth at surf, 10, 100 and bottom to keep consistency in the mean calculation
df_long_merged_depth <- filter(df_long_merged_depth, Group.2 %in% c("Surface", "10m", "100m", "bottom"))

# aggregate table (presence_absence) by fjords and species
df_long_merged_fjord_pres <- aggregate(df_long$count_pres, by = list(df_long$scientificName,
                                                           df_long$Location.ID), FUN = max)

df_long_merged_season_pres <- aggregate(df_long$count_pres, by = list(df_long$scientificName,
                                                            df_long$Season), FUN = max)

df_long_merged_depth_pres <- aggregate(df_long$count_pres, by = list(df_long$scientificName,
                                                           df_long$Feature), FUN = max)

# calculate maximum relative abundance value
colnames(df_long_merged_fjord) <- c("scientificName", "Location.ID", "count")
colnames(df_long_merged_season) <- c("scientificName", "Season", "count")
colnames(df_long_merged_depth) <- c("scientificName", "Feature", "count")

colnames(df_long_merged_fjord_pres) <- c("scientificName", "Location.ID", "count")
colnames(df_long_merged_season_pres) <- c("scientificName", "Season", "count")
colnames(df_long_merged_depth_pres) <- c("scientificName", "Feature", "count")

# calculate maximum value for mean eDNA Index calculation per category
max_relative_ab_df_fjord <- df_long_merged_fjord %>%
  group_by(scientificName) %>%
  summarise(max = max(count))

max_relative_ab_df_season <- df_long_merged_season %>%
  group_by(scientificName) %>%
  summarise(max = max(count))

max_relative_ab_df_depth <- df_long_merged_depth %>%
  group_by(scientificName) %>%
  summarise(max = max(count))

# merge maximum value
df_long_merged_fjord <- left_join(df_long_merged_fjord, max_relative_ab_df_fjord, by = "scientificName")
df_long_merged_season <- left_join(df_long_merged_season, max_relative_ab_df_season, by = "scientificName")
df_long_merged_depth <- left_join(df_long_merged_depth, max_relative_ab_df_depth, by = "scientificName")

# calculate a mean DNA Index value
df_long_merged_fjord$mean_index <- df_long_merged_fjord$count/df_long_merged_fjord$max
df_long_merged_season$mean_index <- df_long_merged_season$count/df_long_merged_season$max
df_long_merged_depth$mean_index <- df_long_merged_depth$count/df_long_merged_depth$max


# reorder x-axis according to North-South gradient for the heatmap
df_long_merged_fjord$Location.ID <- factor(df_long_merged_fjord$Location.ID,
                                     levels = c("Rivers", "Bute", "Toba", "Howe_Sound", "Indian")) #lock in the defined order

df_long_merged_season$Season <- factor(df_long_merged_season$Season,
                                           levels = c("Winter", "Spring", "Summer", "Fall")) #lock in the defined order

df_long_merged_depth$Feature <- factor(df_long_merged_depth$Feature,
                                           levels = c("Surface", "10m", "100m", "bottom")) #lock in the defined order

df_long_merged_fjord_pres$Location.ID <- factor(df_long_merged_fjord$Location.ID,
                                           levels = c("Rivers", "Bute", "Toba", "Howe_Sound", "Indian")) #lock in the defined order

df_long_merged_season_pres$Season <- factor(df_long_merged_season$Season,
                                       levels = c("Winter", "Spring", "Summer", "Fall")) #lock in the defined order

df_long_merged_depth_pres$Feature <- factor(df_long_merged_depth$Feature,
                                       levels = c("Surface", "10m", "100m", "bottom")) #lock in the defined order


# reorder y-axis accodring to taxonomy for the heatmap
df_long_merged_fjord <- df_long_merged_fjord %>% separate(scientificName, c("Genus", "Species"), remove = F)
df_long_merged_fjord <- left_join(df_long_merged_fjord, taxo, by = "Genus")
df_long_merged_fjord$Order <- factor(df_long_merged_fjord$Order, levels = unique(df_long_merged_fjord$Order)) #lock in the defined order

df_long_merged_season <- df_long_merged_season %>% separate(scientificName, c("Genus", "Species"), remove = F)
df_long_merged_season <- left_join(df_long_merged_season, taxo, by = "Genus")
df_long_merged_season$Order <- factor(df_long_merged_season$Order, levels = unique(df_long_merged_season$Order)) #lock in the defined order

df_long_merged_depth <- df_long_merged_depth %>% separate(scientificName, c("Genus", "Species"), remove = F)
df_long_merged_depth <- left_join(df_long_merged_depth, taxo, by = "Genus")
df_long_merged_depth$Order <- factor(df_long_merged_depth$Order, levels = unique(df_long_merged_depth$Order)) #lock in the defined order


# define the fjord color palette
fjord_pal <- c("Indian" = "#fdc086",
               "Howe_Sound" = "#beaed4",
               "Toba" = "#f0027f",
               "Bute" = "#7fc97f",
               "Rivers" = "#386cb0")

# define an Order (taxonmy) color vector
order_pal <- c("Clupeiformes" = "#e31a1c",
               "Salmoniformes" = "#fb9a99",
               "Scorpaeniformes" = "#b2df8a",
               "Gadiformes" = "#33a02c",
               "Argentiniformes" = "#1f78b4",
               "Pleuronectiformes" = "#a6cee3",
               "Myctophiformes" = "#8c2d04",
               "Squaliformes" = "#ff7f00",
               "Ovalentaria" = "#00868b",
               "Perciformes" = "#cab2d6",
               "Osmeriformes" = "#6a3d9a",
               "Trachiniformes" = "#ffff99",
               "Labriformes" = "grey70",
               "Chimaeriformes" = "grey40",
               "Scombriformes"="grey25",
               "Carnivora" = "black",
               "Artiodactyla" = "Grey",
               "Hexanchiformes" = "Grey",
               "Syngnathiformes" = "#005824",
               "Rajiformes" = "#fec44f",
               "Gobiesociformes" = "grey90",
               "Gobiiformes" = "Grey",
               "Cypriniformes" = "yellow",
               "Blenniiformes" = "Grey",
               "Ophidiiformes" = "Grey")

# order y name by the column order
df_long_merged_fjord <- df_long_merged_fjord %>% mutate(scientificName = fct_reorder(scientificName, as.integer(Order)))
df_long_merged_depth <- df_long_merged_depth %>% mutate(scientificName = fct_reorder(scientificName, as.integer(Order)))
df_long_merged_season <- df_long_merged_season %>% mutate(scientificName = fct_reorder(scientificName, as.integer(Order)))

# plot the heatmap
heatmap_fjord <- ggplot(df_long_merged_fjord, aes(x = Location.ID, y = scientificName))+
  geom_point(aes(size = mean_index))+
  scale_colour_manual(values=order_pal)+
  labs(title = "Fjord", 
       x = "Fjord",
       y = "Species")+
  scale_size_area()+
  scale_alpha_continuous(range = c(0,1))+
  coord_fixed(ratio = 0.4)+ # space between points
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust = 1))
heatmap_fjord

heatmap_season <- ggplot(df_long_merged_season, aes(x = Season, y = scientificName))+
  geom_point(aes(size = mean_index))+
  scale_colour_manual(values=order_pal)+
  labs(title = "Season", 
       x = "season",
       y = "Species")+
  scale_size_area()+
  scale_alpha_continuous(range = c(0,1))+
  coord_fixed(ratio = 0.4)+ # space between points
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust = 1))
heatmap_season

heatmap_depth <- ggplot(df_long_merged_depth, aes(x = Feature, y = scientificName))+
  geom_point(aes(size = mean_index))+
  scale_colour_manual(values=order_pal)+
  labs(title = "Depth", 
       x = "depth",
       y = "Species")+
  scale_size_area()+
  scale_alpha_continuous(range = c(0,1))+
  coord_fixed(ratio = 0.4)+ # space between points
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust = 1))
heatmap_depth

################################################
# Then calculate Indval for depth and Season
##################################################

# format df to Indval
df2 <- df
row.names(df2) <- df2$scientificName
df2 <- select(df2, contains(c("2023")))
df2 <- as.data.frame(t(df2))

# load the taxonomy
taxo <- read.csv("01_taxo_all_genus.csv", sep = ",")

# filter metadata for similar samples than edna index top50 samples
metadata <- filter(metadata, Site_ID_Depth_Season %in% row.names(df2))
row.names(metadata) <- metadata$Site_ID_Depth_Season

# reorder table to make order match for Indval
df2 <- df2[order(row.names(df2)),] # reorder table by order
metadata <- metadata[order(row.names(metadata)),] # reorder table by order

# calculate indicator species of depth, season and fjord
library(indicspecies)

indval_sp_depth <- multipatt(df2, metadata$Feature, func = "IndVal.g", control = how(nperm = 999))
indval_sp_fjord <- multipatt(df2, metadata$Location.ID, func = "IndVal.g", control = how(nperm = 999))
indval_sp_season <- multipatt(df2, metadata$Season, func = "IndVal.g", control = how(nperm = 999))

indval_sp_depth_df <- indval_sp_depth$sign
#indval_sp_depth_df <- filter(indval_sp_depth_df, p.value <= 0.05)
indval_sp_depth_list <- row.names(indval_sp_depth_df)
indval_sp_depth_list
indval_sp_depth_df$scientificName <- row.names(indval_sp_depth_df)
indval_sp_depth_df <- indval_sp_depth_df %>% separate(scientificName, c("Genus", "Species"))
indval_sp_depth_df2 <- left_join(indval_sp_depth_df, taxo, by = "Genus")
indval_sp_depth_df3 <- filter(indval_sp_depth_df2, stat>0.5)# filter only significant IndVal species
indval_sp_depth_df3$scientificName <- paste(indval_sp_depth_df3$Genus,indval_sp_depth_df3$Species) # recreate colunm scientificName
indval_sp_depth_df4 <- select(indval_sp_depth_df3, contains(c("scientificName", "Surface", "10m", "100m", "bottom")))


indval_sp_season_df <- indval_sp_season$sign
#indval_sp_season_df <- filter(indval_sp_season_df, p.value <= 0.05)
indval_sp_season_list <- row.names(indval_sp_season_df)
indval_sp_season_list
indval_sp_season_df$scientificName <- row.names(indval_sp_season_df)
indval_sp_season_df <- indval_sp_season_df %>% separate(scientificName, c("Genus", "Species"))
indval_sp_season_df2 <- left_join(indval_sp_season_df, taxo, by = "Genus")
indval_sp_season_df3 <- filter(indval_sp_season_df2, stat>0.5)# filter only significant IndVal species
indval_sp_season_df3$scientificName <- paste(indval_sp_season_df3$Genus,indval_sp_season_df3$Species) # recreate colunm scientificName
indval_sp_season_df4 <- select(indval_sp_season_df3, contains(c("scientificName", "Spring", "Summer", "Fall", "Winter")))


indval_sp_fjord_df <- indval_sp_fjord$sign
#indval_sp_fjord_df <- filter(indval_sp_fjord_df, p.value <= 0.05)
indval_sp_fjord_list <- row.names(indval_sp_fjord_df)
indval_sp_fjord_list
indval_sp_fjord_df$scientificName <- row.names(indval_sp_fjord_df)
indval_sp_fjord_df <- indval_sp_fjord_df %>% separate(scientificName, c("Genus", "Species"))
indval_sp_fjord_df2 <- left_join(indval_sp_fjord_df, taxo, by = "Genus")
indval_sp_fjord_df3 <- filter(indval_sp_fjord_df2, stat>0.5)# filter only significant IndVal species
indval_sp_fjord_df3$scientificName <- paste(indval_sp_fjord_df3$Genus,indval_sp_fjord_df3$Species) # recreate colunm scientificName
indval_sp_fjord_df4 <- select(indval_sp_fjord_df3, contains(c("scientificName", "Rivers", "Bute", "Toba","Howe", "Indian")))

# join table
indval_df <- left_join(df_long_merged_fjord, indval_sp_fjord_df4, by = "scientificName")
indval_df <- left_join(indval_df, indval_sp_season_df4, by = "scientificName")
indval_df <- left_join(indval_df, indval_sp_depth_df4, by = "scientificName")

# pivot column depth, season and fjord
colnames(indval_df)
indval_df <- indval_df %>% pivot_longer(cols = c("s.Rivers", "s.Bute", "s.Toba","s.Howe_Sound", "s.Indian"),  names_to = "Indval_fjord", values_to = "Indval_fjord_value")
indval_df <- indval_df %>% pivot_longer(cols = c("s.Spring", "s.Summer", "s.Fall","s.Winter"),  names_to = "Indval_season", values_to = "Indval_season_value")
indval_df <- indval_df %>% pivot_longer(cols = c("s.Surface", "s.10m", "s.100m","s.bottom"),  names_to = "Indval_depth", values_to = "Indval_depth_value")

# replace NA value by zero
indval_df <- indval_df %>% mutate(Indval_fjord_value = if_else(is.na(Indval_fjord_value), 0, Indval_fjord_value))
indval_df <- indval_df %>% mutate(Indval_season_value = if_else(is.na(Indval_season_value), 0, Indval_season_value))
indval_df <- indval_df %>% mutate(Indval_depth_value = if_else(is.na(Indval_depth_value), 0, Indval_depth_value))

# change s. name
indval_df$Indval_fjord <- gsub('s.','',indval_df$Indval_fjord)
indval_df$Indval_season <- gsub('s.','',indval_df$Indval_season)
indval_df$Indval_depth <- gsub('s.','',indval_df$Indval_depth)

# order y name by the column order
indval_df <- indval_df %>% mutate(scientificName = fct_reorder(scientificName, as.integer(Order)))

heatmap_indval_fjord <- ggplot(indval_df, aes(x = Indval_fjord, scientificName))+
  geom_point(aes(alpha = Indval_fjord_value, color = Order), size = 4)+
  scale_colour_manual(values=order_pal)+
  labs(title = "B", 
       x = "Fjord",
       y = "Species",
       color = "Order")+
  scale_size_area()+
  scale_alpha_continuous(range = c(0,1))+
  coord_fixed(ratio = 0.4)+ # space between points
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
heatmap_indval_fjord


# for season indval
season_pal <- c("Fall"="#F8766D", 
                  "Spring"="#7CAE00",
                  "Summer"="#007643",
                  "Winter"="#00BFC4",
                  "Others"="grey8") # season color vector

heatmap_indval_season <- ggplot(indval_df, aes(x = Indval_season, scientificName))+
  geom_point(aes(alpha = Indval_season_value, color = Order), size = 4)+
  scale_colour_manual(values=order_pal)+
  labs(title = "C", 
       x = "season",
       y = "Species",
       color = "Order")+
  scale_size_area()+
  scale_alpha_continuous(range = c(0,1))+
  coord_fixed(ratio = 0.4)+ # space between points
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
heatmap_indval_season

# for depth indval
depth_pal <- c("Surface"="red", 
                "10m"="orange",
                "100m"="grey",
                "bottom"="black") # depth color vector

# reorder x-axis for the heatmap
indval_df$Indval_depth <- factor(indval_df$Indval_depth, levels = c("Surface", "10m", "100m", "bottom")) #lock in the defined order

heatmap_indval_depth <- ggplot(indval_df, aes(x = Indval_depth, scientificName))+
  geom_point(aes(alpha = Indval_depth_value, color = Order), size = 4)+
  scale_colour_manual(values=order_pal)+
  labs(title = "D", 
       x = "depth",
       y = "Species",
       color = "Order")+
  scale_size_area()+
  scale_alpha_continuous(range = c(0,1))+
  coord_fixed(ratio = 0.4)+ # space between points
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
heatmap_indval_depth


# plot them all!
library(ggpubr)

ggarrange(heatmap_fjord,heatmap_indval_fjord, heatmap_indval_season, heatmap_indval_depth,
          nrow = 1, ncol = 4, common.legend = T, align = "hv")

# plot the heatmap

# filter data for points with value higher than 0
df_long_merged_fjord <- df_long_merged_fjord %>%
  filter(mean_index > 0.01)

# import and merge common names
common_name_df <- read.csv("01_updated_common_names.csv")
df_long_merged_fjord <- left_join(df_long_merged_fjord, common_name_df, by = "scientificName")
colnames(df_long_merged_fjord)

# order y name by the column order
df_long_merged_fjord <- df_long_merged_fjord %>% mutate(scientificName = fct_reorder(scientificName, as.integer(Order)))
df_long_merged_depth <- df_long_merged_depth %>% mutate(scientificName = fct_reorder(scientificName, as.integer(Order)))
df_long_merged_season <- df_long_merged_season %>% mutate(scientificName = fct_reorder(scientificName, as.integer(Order)))

# plot it
heatmap_fjord <- ggplot(df_long_merged_fjord, aes(x = Location.ID, y = scientificName))+
  geom_point(aes(size = mean_index))+
  geom_point(data = filter(indval_df, Indval_fjord_value > 0),
             aes(x = Indval_fjord, y = scientificName, alpha = Indval_fjord_value),
             shape = 18, color = "red")+
  scale_y_discrete(labels = setNames(df_long_merged_fjord$FullNameAbbreviated, df_long_merged_fjord$scientificName))+
  scale_colour_manual(values=order_pal)+
  labs(title = "Fjord", 
       x = "Fjord",
       y = "Species")+
  scale_size_area()+
  scale_alpha_continuous(range = c(0,1))+
  coord_fixed(ratio = 0.4)+ # space between points
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust = 1))
heatmap_fjord

# filter data for points with value higher than 0
df_long_merged_season <- df_long_merged_season %>%
  filter(mean_index > 0.01)

heatmap_season <- ggplot(df_long_merged_season, aes(x = Season, y = scientificName))+
  geom_point(aes(size = mean_index))+
  #geom_point(data = df_long_merged_season_with_border, aes(x = Season, y = scientificName), 
  #           shape = 21, color = "black", size = 4, stroke = 0.5) +
  geom_point(data = filter(indval_df, Indval_season_value > 0) , aes(x = Indval_season, y = scientificName, alpha = Indval_season_value), shape = 18, color = "red")+
  scale_colour_manual(values=order_pal)+
  labs(title = "Season", 
       x = "season",
       y = "Species")+
  scale_size_area()+
  scale_alpha_continuous(range = c(0,1))+
  coord_fixed(ratio = 0.4)+ # space between points
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust = 1), axis.text.y = element_blank())
heatmap_season

# filter data for points with value higher than 0
df_long_merged_depth <- df_long_merged_depth %>%
  filter(mean_index > 0.01)

heatmap_depth <- ggplot(df_long_merged_depth, aes(x = Feature, y = scientificName))+
  geom_point(aes(size = mean_index))+
  #geom_point(data = df_long_merged_depth_with_border, aes(x = Feature, y = scientificName), 
  #           shape = 21, color = "black", size = 4, stroke = 0.5) +
  geom_point(data = filter(indval_df, Indval_depth_value > 0), aes(x = Indval_depth, y = scientificName, alpha = Indval_depth_value), shape = 18, color = "red")+
  scale_colour_manual(values=order_pal)+
  labs(title = "Depth", 
       x = "depth",
       y = "Species")+
  scale_size_area()+
  scale_alpha_continuous(range = c(0,1))+
  coord_fixed(ratio = 0.4)+ # space between points
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust = 1), axis.text.y = element_blank())
heatmap_depth

# Plot them all!
# Save as 13*13
ggarrange(heatmap_fjord,heatmap_season, heatmap_depth,
          nrow = 1, ncol = 3, common.legend = T, align = "hv")

