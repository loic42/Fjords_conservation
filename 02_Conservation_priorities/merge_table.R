setwd("C:/Users/loicj/OneDrive - UBC/01_Fjord_eDNA/05_Scratch/NEW")

library(dplyr)

# load list Fishes from the Salish Sea
df_FishBase <- read.csv("06_Conservation_priorities/Tables/BCFjords-FBSLB-20240704.csv", sep = ",", header = T)
df_Gale <- read.csv("06_Conservation_priorities/Tables/Gale et al. - Table.csv", sep = ",", header = T)

# filter table
colnames(df_Gale)
df_FishBase_filter <- df_FishBase[c("Genus", "scientificName", "Troph","Max_length", "Vulnerability_fisheries", 
                     "Vulnerability_fisheries_category", "Vulnerability_climate", "Vulnerability_climate_category")]

df_Gale_filter <- df_Gale[c("Common_Name", "scientificName", "Upper_level_Predators","Forage_Species","Any_CC",
                            "Global_CC","National_CC", "Regional_CC", "Recommended_CP", "Nutrient_Transporting_Species")]

# load list of species from the fjords
df_fjord <- read.csv("01_Tables/12S_species_table_BC-db.csv", sep = ",", header = T)
df_fjord <- filter(df_fjord, Species != "sp") # filter out species annotated as sp.
fjord_sp <- as.data.frame(df_fjord$scientificName) # extract list of detected species
colnames(fjord_sp) <- "scientificName"

# find intersection and unique species of both list
# for these species, criteria will be checked and added manually
inter_sp <- intersect(fjord_sp$scientificName, df_FishBase_filter$scientificName)
unique_fjord_sp_Gale <- as.data.frame(setdiff(fjord_sp$scientificName, df_Gale_filter$scientificName))
unique_fjord_sp_Fishbase <- as.data.frame(setdiff(fjord_sp$scientificName, df_FishBase_filter$scientificName))

#write.csv(list_combined, "test.csv")

# merge fjords and CPs table
df_CPs_merged <- left_join(fjord_sp, df_FishBase_filter, by = "scientificName")
df_CPs_merged2 <- left_join(df_CPs_merged, df_Gale_filter, by = "scientificName")

# add column to know where to manually check
unique_fjord_sp_Gale$check_gale <- "Yes"
unique_fjord_sp_Fishbase$check_fishbase <- "Yes"

colnames(unique_fjord_sp_Gale) <- c("scientificName", "check_gale")
colnames(unique_fjord_sp_Fishbase) <- c("scientificName", "check_fishbase")


df_CPs_merged3 <- left_join(df_CPs_merged2,unique_fjord_sp_Gale, by = "scientificName")
df_CPs_merged4 <- left_join(df_CPs_merged3,unique_fjord_sp_Fishbase, by = "scientificName")


# write csv
write.csv(df_CPs_merged4, "06_Conservation_priorities/01_conservation_priorities_final_table.csv")


############################################################
# Assess the status of species not list in Gale et al.(2018)
############################################################
sp_list <- read.csv("06_Conservation_priorities/Tables/remaining_taxa.csv")

BC_list <- read.csv("06_Conservation_priorities/Tables/Conservation_lists/BC_list.csv")
IUCN_list <- read.csv("06_Conservation_priorities/Tables/Conservation_lists/IUCN_list.csv")
SARA_list <- read.csv("06_Conservation_priorities/Tables/Conservation_lists/SARA_list.csv")
WILD_list <- read.csv("06_Conservation_priorities/Tables/Conservation_lists/wild_species_list.csv")

intersect(sp_list$scientificName, BC_list$Scientific.Name)
intersect(sp_list$scientificName,IUCN_list$scientificName)
intersect(sp_list$scientificName,SARA_list$Scientific.name)
intersect(sp_list$scientificName,WILD_list$scientific_name)

################################################################
# Add common name to the table
################################################################
sp_df <- read.csv("06_Conservation_priorities/01_conservation_priorities_final_table.csv")
common_name_df <- read.csv("01_updated_common_names.csv")

x <- left_join(sp_df, common_name_df, by = "scientificName")

write.csv(x, "test.csv")
