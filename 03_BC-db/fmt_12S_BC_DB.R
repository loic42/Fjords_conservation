# Import MIDORI database and filter it for Pacific taxa only using OBIS as reference

##########################################################
# Load working envt and packages
##########################################################
setwd("C:/Users/loicj/OneDrive - UBC/01_Fjord_eDNA/12_BC_Databases")

library(tidyr)
library(tidyverse)
library(robis)
library(worrms)
library(dplyr)

########################
#import data
########################

# import raw species table from MIDORI (with only genus, species and total as header)
# this list was extracted from the MIDORI fasta file
midori_12S_base <- read.csv("files/midori.12S.up.list.csv", header = T, sep = ";")

###################################################
# find ambiguous taxa
###################################################
# see https://api.obis.org/#/Occurrence/get_occurrence__id_
# see the map: https://mapper.obis.org/?areaid=31906,31912

# keep only row with non multiple hits
list_midori_base <- as.data.frame(midori_12S_base$ScientificName)

colnames(list_midori_base) <- "scientificName"

# get a list of area ID
areas <- area(verbose = F)
    
# extract species information for regions of interest (run just once)
# see here: https://mapper.obis.org/?areaid=35,265
# for BC: canada + alaska identifier
species_OBIS_df <- checklist(areaid = c("35", "265"))

# Filter the df to extract group of interest (fish, birds, mammals and sharks)
species_OBIS_df_BC <- filter(species_OBIS_df, class %in% (c("Aves", "Elasmobranchii", "Actinopteri", "Mammalia")))

write.csv(species_OBIS_df_BC, "OBIS_BC_Alaska_table.csv")

