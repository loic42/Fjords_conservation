# Objective: Turn raw species table into exploitable table for subsequent analysis
# Steps:
# 1) Import eDNA species raw table
# 2) Aggregate multiple hits
# 3) Filter by taxonomy
# 4) Field negative control correction. Turn to zero value that are below the maximum reads count in negative controls
# 5) Extraction negative control correction. Turn to zero value that are below the maximum reads count in negative controls (library specific).
# 6) Turn to zero occurrence of 1 (singletons)
# 7) Merge replicates
# 8) Remove samples with low counts (<1000)
# 9) Remove species with total reads < 10
# 10) Calculate relative abundance and eDNA Index table
# 11) Export tables as csv
# 12) Format table for shiny applications

# NB: Date of Negative control in species table should be in the following format: "Summer_EC-16-01-2024"

##########################################################
# Loading working environment and packages 
##########################################################
setwd("C:/Users/loicj/OneDrive - UBC/01_Fjord_eDNA/05_Scratch/NEW")

library(readr)
library(devtools)
library(ggplot2)
library(ape)
library(vegan)
library(picante)
library(permute)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(robis)
library(worrms)
library(dplyr)
library(matrixStats)

sessionInfo()

########################
# Import data
# original fjords data, raw reads processed with BARQUE and annotated with BC database at 97% identity
# NB: Sebastes species as been turned manually to Sebastes sp.
########################

asv_table <- read.csv("00_Data/12s_MiFish_species_table_multiple_added_2023.csv", header = T, sep = ",")
asv_taxonomy <- read.csv("01_taxo_all_genus.csv", header = T, sep = ",")# import taxonomy
asv_df <- left_join(asv_table, asv_taxonomy, by = "Genus") # merge ASV table and taxonomy
metadata <- read.csv("01_Metadata_all.csv", header = T, sep = ",")
  

############################  
# aggregate species table
############################
colnames(asv_df)

# species level
sp_df <- aggregate(asv_df[5:(ncol(asv_df)-5)], by = list(asv_df$Genus, asv_df$scientificName), FUN = sum)
sp_df <- sp_df %>% rename(Genus = Group.1, scientificName = Group.2)

# verify if duplicated name
sp_df[duplicated(sp_df$scientificName)]
colnames(sp_df)

########################
# merge table with taxo
########################

# merge species table with taxonomy
asv_taxonomy_unique_sp <- asv_taxonomy[!duplicated(asv_taxonomy$scientificName), ]# remove duplicate from taxo table
sp_merged <- left_join(sp_df, asv_taxonomy, by = "Genus")


#############################
# filter table by taxonomy
#############################
colnames(sp_merged)
colnames(sp_df)

# Creating a not in operator:
`%notin%` <- Negate(`%in%`)

# keep only phylum chordata (include fish and marine mammals)
unique(sp_merged$Phylum)
phylum_list <- c("Chordata")
sp_merged_phylum <- filter(sp_merged, Phylum %in% phylum_list)

# filter at order level
# remove human dna and birds (to not overinterpret vertical effect)
unique(sp_merged_phylum$Order)
order_list <- c("Primates","Anseriformes", "Gaviiformes", "Charadriiformes", "Passeriformes")
sp_merged_order <- filter(sp_merged_phylum, Order %notin% order_list)

# filter at class level 
# remove insecta, reptiles and human food consumption
# NB: Not needed if using curated fish-specific database
class_list <- c("Insecta", "Amphibia", "Chilopoda", "Arachnida", "Collembola", "NA", 
                "Lepidosauria", "Clitellata", "Reptilia", "Rhinolophidae", "Chlamyphoridae",
               "Cercopithecidae", "Elephantidae", "Gomphotheriidae", "Capromyidae")

family_list <- c("Bovidae", "Felidae","Equidae", "Mustelidae", "Procyonidae", "Vespertilionidae","Erinaceidae", "Leporidae",
                 "Soricidae", "Hominidae","Phasianidae", "Castoridae", "Cricetidae", "Muridae", "Sciuridae", "Tupaiidae", "Dasyuridae",
                 "Suidae", "Cervidae", "NA", "Scincidae", "Australasian", "Mammutidae", "Hominidae", "Vespertilionidae", "Canidae")

# remove non-marine species at Order level
sp_merged_class <- filter(sp_merged_order, Class %notin% class_list)
sp_merged_family <- filter(sp_merged_class, Family %notin% family_list)
sp_df <- sp_merged_family

# keep a list of removed taxa
removed_df1 <- filter(sp_merged, Phylum %notin% phylum_list)
removed_df2 <- filter(sp_merged, Order %in% order_list)
removed_df3 <- filter(sp_merged, Class %in% class_list)
removed_df4 <- filter(sp_merged, Family %in% family_list)

# extract a list of removed taxa
list_removed <- unique(c(removed_df1$scientificName, removed_df2$scientificName,
                  removed_df3$scientificName, removed_df4$scientificName))

list_removed_taxo <- filter(sp_merged, scientificName %in% list_removed)
list_removed_taxo <- dplyr::select(list_removed_taxo, contains(c("scientificName", "Phylum", "Class", "Order",
                                                                 "Genus", "Species", "Total")))

write.csv(list_removed_taxo, "01_Tables/list_removed_taxa_2023.csv")


###################################################
# Field Negative control: Library - specific
###################################################
colnames(sp_df)

# subset negative control columns (one table per library)
# NB: fjedna one are field blanks
neg_12S_winter <- dplyr::select(sp_df, contains(c("scientificName", "winter")))
neg_12S_spring <- dplyr::select(sp_df, contains(c("scientificName", "spring")))
neg_12S_summer <- dplyr::select(sp_df, contains(c("scientificName", "summer","FjeDNA441", "FjeDNA484", "FjeDNA521", "FjeDNA578", "FjeDNA481")))
neg_12S_fall <- dplyr::select(sp_df, contains(c("scientificName", "fall","FjeDNA670", "FjeDNA671")))

# keep only field negative controls
neg_12S_winter <- dplyr::select(neg_12S_winter, contains(c("scientificName", "field", "FjeDNA")))
neg_12S_spring <- dplyr::select(neg_12S_spring, contains(c("scientificName","field", "FjeDNA")))
neg_12S_summer <- dplyr::select(neg_12S_summer, contains(c("scientificName","field", "FjeDNA")))
neg_12S_fall <- dplyr::select(neg_12S_fall, contains(c("scientificName", "field", "FjeDNA")))

# keep the maximum value among each  negative control
neg_12S_winter$max_negative <- rowMaxs(as.matrix(neg_12S_winter[,2:ncol(neg_12S_winter)]))
neg_12S_spring$max_negative <- rowMaxs(as.matrix(neg_12S_spring[,2:ncol(neg_12S_spring)]))
neg_12S_summer$max_negative <- rowMaxs(as.matrix(neg_12S_summer[,2:ncol(neg_12S_summer)]))
neg_12S_fall$max_negative <- rowMaxs(as.matrix(neg_12S_fall[,2:ncol(neg_12S_fall)]))

# keep only species-specific maximum negative column and scientificName
neg_12S_winter <- dplyr::select(neg_12S_winter, contains(c("scientificName","max_negative")))
neg_12S_spring <- dplyr::select(neg_12S_spring, contains(c("scientificName","max_negative")))
neg_12S_summer <- dplyr::select(neg_12S_summer, contains(c("scientificName","max_negative")))
neg_12S_fall <- dplyr::select(neg_12S_fall, contains(c("scientificName","max_negative")))

# subset metadata by library
unique(metadata$Survey)
metadata_winter <- filter(metadata, Survey == "FJORDS_BIODIVERSITY" & Season == "Winter")
metadata_spring <- filter(metadata, Survey == "FJORDS_BIODIVERSITY" & Season == "Spring")
metadata_summer <- filter(metadata, Survey == "FJORDS_BIODIVERSITY" & Season == "Summer")
metadata_fall <- filter(metadata, Survey == "FJORDS_BIODIVERSITY" & Season == "Fall")

# subset table by library
sp_df_winter_library <- dplyr::select(sp_df, contains(c("scientificName", metadata_winter$Hakai.ID)))
sp_df_spring_library <- dplyr::select(sp_df, contains(c("scientificName", metadata_spring$Hakai.ID)))
sp_df_summer_library <- dplyr::select(sp_df, contains(c("scientificName", metadata_summer$Hakai.ID)))
sp_df_fall_library <- dplyr::select(sp_df, contains(c("scientificName", metadata_fall$Hakai.ID)))

# merge tables
sp_df_winter_long <- left_join(sp_df_winter_library, neg_12S_winter, by = "scientificName")
sp_df_spring_long <- left_join(sp_df_spring_library, neg_12S_spring, by = "scientificName")
sp_df_summer_long <- left_join(sp_df_summer_library, neg_12S_summer, by = "scientificName")
sp_df_fall_long <- left_join(sp_df_fall_library, neg_12S_fall, by = "scientificName")

# turn table by keeping max_negative as a column
sp_df_winter_long <- sp_df_winter_long %>% pivot_longer(cols = 2:(ncol(sp_df_winter_long)-1) , names_to = "station", values_to = "count") # turn APC df to tidy
sp_df_spring_long <- sp_df_spring_long %>% pivot_longer(cols = 2:(ncol(sp_df_spring_long)-1) , names_to = "station", values_to = "count") # turn APC df to tidy
sp_df_summer_long <- sp_df_summer_long %>% pivot_longer(cols = 2:(ncol(sp_df_summer_long)-1) , names_to = "station", values_to = "count") # turn APC df to tidy
sp_df_fall_long <- sp_df_fall_long %>% pivot_longer(cols = 2:(ncol(sp_df_fall_long)-1) , names_to = "station", values_to = "count") # turn APC df to tidy

# apply the filter: turn to zero value that are below the maximum reads count threshold in negative controls
sp_df_winter_long$count_reads_corrected <- ifelse(sp_df_winter_long$count > sp_df_winter_long$max_negative, sp_df_winter_long$count, 0) # loop to correct
sp_df_spring_long$count_reads_corrected <- ifelse(sp_df_spring_long$count > sp_df_spring_long$max_negative, sp_df_spring_long$count, 0) # loop to correct
sp_df_summer_long$count_reads_corrected <- ifelse(sp_df_summer_long$count > sp_df_summer_long$max_negative, sp_df_summer_long$count, 0) # loop to correct
sp_df_fall_long$count_reads_corrected <- ifelse(sp_df_fall_long$count > sp_df_fall_long$max_negative, sp_df_fall_long$count, 0) # loop to correct

# Keep only corrected reads columns
sp_df_winter_long <- subset(sp_df_winter_long, select=c(scientificName, station, count_reads_corrected)) # remove column used for correction
sp_df_spring_long <- subset(sp_df_spring_long, select=c(scientificName, station, count_reads_corrected)) # remove column used for correction
sp_df_summer_long <- subset(sp_df_summer_long, select=c(scientificName, station, count_reads_corrected)) # remove column used for correction
sp_df_fall_long <- subset(sp_df_fall_long, select=c(scientificName, station, count_reads_corrected)) # remove column used for correction

# re-order table according to station
sp_df_winter_long <- sp_df_winter_long[order(sp_df_winter_long$station),]
sp_df_spring_long <- sp_df_spring_long[order(sp_df_spring_long$station),]
sp_df_summer_long <- sp_df_summer_long[order(sp_df_summer_long$station),]
sp_df_fall_long <- sp_df_fall_long[order(sp_df_fall_long$station),]

###################################################################
# Extraction Negative Controls correction
###################################################################

# subset metadata to keep only extraction dates
neg_extr_df <- subset(metadata,select = c('Hakai.ID','Extraction_date'))
neg_extr_df <- neg_extr_df[!is.na(neg_extr_df$Extraction_date), ]

# subset negative controls columns
neg_12S_winter <- dplyr::select(sp_df, contains(c("scientificName", "winter")))
neg_12S_spring <- dplyr::select(sp_df, contains(c("scientificName", "spring")))
neg_12S_summer <- dplyr::select(sp_df, contains(c("scientificName", "summer","FjeDNA441", "FjeDNA484", "FjeDNA521", "FjeDNA578")))
neg_12S_fall <- dplyr::select(sp_df, contains(c("scientificName", "fall", "FjeDNA670", "FjeDNA671")))

# keep only extraction negative controls
neg_12S_winter <- dplyr::select(neg_12S_winter, contains(c("scientificName", "EC")))
neg_12S_spring <- dplyr::select(neg_12S_spring, contains(c("scientificName","EC")))
neg_12S_summer <- dplyr::select(neg_12S_summer, contains(c("scientificName","EC")))
neg_12S_fall <- dplyr::select(neg_12S_fall, contains(c("scientificName", "EC")))

# pivot negative extraction table
neg_12S_winter_long <- neg_12S_winter %>% pivot_longer(cols = 2:(ncol(neg_12S_winter)) , names_to = "Date_EC", values_to = "count_EC") # turn APC df to tidy
neg_12S_spring_long <- neg_12S_spring %>% pivot_longer(cols = 2:(ncol(neg_12S_spring)) , names_to = "Date_EC", values_to = "count_EC") # turn APC df to tidy
neg_12S_summer_long <- neg_12S_summer %>% pivot_longer(cols = 2:(ncol(neg_12S_summer)) , names_to = "Date_EC", values_to = "count_EC") # turn APC df to tidy
neg_12S_fall_long <- neg_12S_fall %>% pivot_longer(cols = 2:(ncol(neg_12S_fall)) , names_to = "Date_EC", values_to = "count_EC") # turn APC df to tidy

# separate the date column
neg_12S_winter_long <- neg_12S_winter_long %>% separate_wider_delim(Date_EC, "_EC.", names = c("season", "Extraction_date"))
neg_12S_spring_long <- neg_12S_spring_long %>% separate_wider_delim(Date_EC, "_EC.", names = c("season", "Extraction_date"))
neg_12S_summer_long <- neg_12S_summer_long %>% separate_wider_delim(Date_EC, "_EC.", names = c("season", "Extraction_date"))
neg_12S_fall_long <- neg_12S_fall_long %>% separate_wider_delim(Date_EC, "_EC.", names = c("season", "Extraction_date"))

# replace "." by "-"
neg_12S_winter_long$Extraction_date <- gsub('\\.', '-', neg_12S_winter_long$Extraction_date)
neg_12S_spring_long$Extraction_date <- gsub('\\.', '-', neg_12S_spring_long$Extraction_date)
neg_12S_summer_long$Extraction_date <- gsub('\\.', '-', neg_12S_summer_long$Extraction_date)
neg_12S_fall_long$Extraction_date <- gsub('\\.', '-', neg_12S_fall_long$Extraction_date)

# merge metadata with negative extraction table
neg_12S_winter_long <- left_join(neg_12S_winter_long,neg_extr_df, by = "Extraction_date")
neg_12S_spring_long <- left_join(neg_12S_spring_long,neg_extr_df, by = "Extraction_date")
neg_12S_summer_long <- left_join(neg_12S_summer_long,neg_extr_df, by = "Extraction_date")
neg_12S_fall_long <- left_join(neg_12S_fall_long,neg_extr_df, by = "Extraction_date")

# remove row with "NA" (means no extraction negative controls for this sample)
neg_12S_winter_long <- neg_12S_winter_long[!is.na(neg_12S_winter_long$Hakai.ID), ]
neg_12S_spring_long <- neg_12S_spring_long[!is.na(neg_12S_spring_long$Hakai.ID), ]
neg_12S_summer_long <- neg_12S_summer_long[!is.na(neg_12S_summer_long$Hakai.ID), ]
neg_12S_fall_long <- neg_12S_fall_long[!is.na(neg_12S_fall_long$Hakai.ID), ]

# resubset negative control table
neg_12S_winter_long <- subset(neg_12S_winter_long,select = c("scientificName", "Extraction_date", "count_EC", "Hakai.ID"))
neg_12S_spring_long <- subset(neg_12S_spring_long,select = c("scientificName", "Extraction_date", "count_EC", "Hakai.ID"))
neg_12S_summer_long <- subset(neg_12S_summer_long,select = c("scientificName", "Extraction_date", "count_EC", "Hakai.ID"))
neg_12S_fall_long <- subset(neg_12S_fall_long,select = c("scientificName", "Extraction_date", "count_EC", "Hakai.ID"))

# merge species table with negative control table
sp_df_winter_long <- left_join(sp_df_winter_long, neg_12S_winter_long, by = c("station" = "Hakai.ID", "scientificName"))
sp_df_spring_long <- left_join(sp_df_spring_long, neg_12S_spring_long, by = c("station" = "Hakai.ID", "scientificName"))
sp_df_summer_long <- left_join(sp_df_summer_long, neg_12S_summer_long, by = c("station" = "Hakai.ID", "scientificName"))
sp_df_fall_long <- left_join(sp_df_fall_long, neg_12S_fall_long, by = c("station" = "Hakai.ID", "scientificName"))

# Turn count_EC with "NA" values to zero
sp_df_winter_long$count_EC[is.na(sp_df_winter_long$count_EC)] <- 0
sp_df_spring_long$count_EC[is.na(sp_df_spring_long$count_EC)] <- 0
sp_df_summer_long$count_EC[is.na(sp_df_summer_long$count_EC)] <- 0
sp_df_fall_long$count_EC[is.na(sp_df_fall_long$count_EC)] <- 0

# apply the filter: turn to zero value that are below the maximum reads count in negative controls
sp_df_winter_long$count_reads_corrected <- ifelse(sp_df_winter_long$count_reads_corrected > sp_df_winter_long$count_EC, sp_df_winter_long$count_reads_corrected, 0) # loop to correct
sp_df_spring_long$count_reads_corrected <- ifelse(sp_df_spring_long$count_reads_corrected > sp_df_spring_long$count_EC, sp_df_spring_long$count_reads_corrected, 0) # loop to correct
sp_df_summer_long$count_reads_corrected <- ifelse(sp_df_summer_long$count_reads_corrected > sp_df_summer_long$count_EC, sp_df_summer_long$count_reads_corrected, 0) # loop to correct
sp_df_fall_long$count_reads_corrected <- ifelse(sp_df_fall_long$count_reads_corrected > sp_df_fall_long$count_EC, sp_df_fall_long$count_reads_corrected, 0) # loop to correct

# resubset long table
sp_df_winter_long <- subset(sp_df_winter_long,select = c("scientificName", "station", "count_reads_corrected"))
sp_df_spring_long <- subset(sp_df_spring_long,select = c("scientificName", "station", "count_reads_corrected"))
sp_df_summer_long <- subset(sp_df_summer_long,select = c("scientificName", "station", "count_reads_corrected"))
sp_df_fall_long <- subset(sp_df_fall_long,select = c("scientificName", "station", "count_reads_corrected"))

###################################################################
# Filter: Turn to zero occurrence of 1 and re-assemble final table
###################################################################

sp_df_winter_long$count_reads_corrected <- ifelse(sp_df_winter_long$count_reads_corrected < 2, 0, sp_df_winter_long$count_reads_corrected) # loop to correct
sp_df_spring_long$count_reads_corrected <- ifelse(sp_df_spring_long$count_reads_corrected < 2, 0, sp_df_spring_long$count_reads_corrected) # loop to correct
sp_df_summer_long$count_reads_corrected <- ifelse(sp_df_summer_long$count_reads_corrected < 2, 0, sp_df_summer_long$count_reads_corrected) # loop to correct
sp_df_fall_long$count_reads_corrected <- ifelse(sp_df_fall_long$count_reads_corrected < 2, 0, sp_df_fall_long$count_reads_corrected) # loop to correct

# re-pivot table
sp_df_winter_library_filtered <- sp_df_winter_long %>% pivot_wider(names_from = station, values_from = count_reads_corrected) # turn df to tidy
sp_df_spring_library_filtered <- sp_df_spring_long %>% pivot_wider(names_from = station, values_from = count_reads_corrected) # turn df to tidy
sp_df_summer_library_filtered <- sp_df_summer_long %>% pivot_wider(names_from = station, values_from = count_reads_corrected) # turn df to tidy
sp_df_fall_library_filtered <- sp_df_fall_long %>% pivot_wider(names_from = station, values_from = count_reads_corrected) # turn df to tidy

# merge distinct library table into one
sp_df <- Reduce(function (...) { merge(..., all = F) },   # full join
                        list(sp_df_winter_library_filtered,
                             sp_df_spring_library_filtered,
                             sp_df_summer_library_filtered,
                             sp_df_fall_library_filtered))

########################################################
# Merge replicates
########################################################
colnames(sp_df)
colnames(metadata)

# turn table by keeping max_negative as a column
sp_df_long <- sp_df %>% pivot_longer(cols = 2:(ncol(sp_df)) , names_to = "Hakai.ID", values_to = "count") # turn APC df to tidy

# filter and merge metadata
metadata_filter <- dplyr::select(metadata, contains(c("Hakai.ID", "Site_ID_Depth_Season")))

sp_df_long_merged <- left_join(sp_df_long, metadata_filter, by = "Hakai.ID")# if error, check duplicated Hakai.ID in metadata file

# SUM replicates together
sp_df_long_merged <- aggregate(sp_df_long_merged$count, by = list(sp_df_long_merged$scientificName,
                                                                   sp_df_long_merged$Site_ID_Depth_Season),
                               FUN = sum, na.rm = T)

# change colnames
colnames(sp_df_long_merged) <- c("scientificName", "ID", "count")

# re-pivot table
sp_df_merged <- sp_df_long_merged %>% pivot_wider(names_from = ID, values_from = count) # turn df to tidy
colnames(sp_df_merged)

#########################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# USE ONLY IF WANT NO REPLICATES MERGING
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#sp_df_merged <- sp_df
#########################################

###################################################
# Visualize and Remove low count stations (<1000)
###################################################

sp_df_sum <- as.data.frame(rbind(colnames(sp_df_merged[,2:ncol(sp_df_merged)]), colSums(sp_df_merged[,2:ncol(sp_df_merged)])))
colnames(sp_df_sum) <- sp_df_sum[1,]
sp_df_sum <- sp_df_sum[-1,]
sp_df_sum

# pivot table 
sp_df_sum <- sp_df_sum %>% pivot_longer(cols = 1:ncol(sp_df_sum) , names_to = "station", values_to = "count") # turn APC df to tidy

# filter to find low abundant samples
sp_df_sum_low <- filter(sp_df_sum, as.numeric(count) < 1000)

# remove samples with less than 1000 reads
sp_df_merged <- dplyr::select(sp_df_merged, -matches(sp_df_sum_low$station))

####################################################
# remove species that have less than 2 total reads
####################################################

# make a column sum
sp_df_merged$sum <- rowSums(sp_df_merged[,2:ncol(sp_df_merged)])
  
# filter species with total reads < 10
sp_df_merged <- filter(sp_df_merged, sum >= 10)

# remove the sum column
sp_df_merged <- select(sp_df_merged, -c(sum))

######################
# turn into relative
######################
colnames(sp_df_merged)
sp_df_relative <- as.data.frame(sapply(sp_df_merged[,2:ncol(sp_df_merged)], prop.table))

# change row names
rownames(sp_df_relative) <- sp_df_merged$scientificName

################################################
# filter empty species lines
################################################
colnames(sp_df_relative)

# add a column total reads
sp_df_relative$Total <- rowSums(sp_df_relative)

# add a column scientificName
sp_df_relative$scientificName <- row.names(sp_df_relative)

# filter  empty species lines
sp_df_relative_rm <- filter(sp_df_relative, Total == 0)
sp_df_relative <- filter(sp_df_relative, Total > 0)

# Remove column "Total" using subset
sp_df_relative <- subset(sp_df_relative, select = -Total)

######################
# Calculate edna index
######################
colnames(sp_df_relative)

# vector of maximum relative abundance for each row
Maximum_vec <- apply(sp_df_relative[,1:(ncol(sp_df_relative)-1)],1,max)

# divide by the vector to create an abundance index ranging from 0 to 1
df_index <- as.data.frame(as.matrix(sp_df_relative[,1:(ncol(sp_df_relative)-1)]) / Maximum_vec)

###############################
# add full taxonomy to tables
###############################
taxonomy <- read.csv("01_taxo_all_genus.csv", sep = ",")

# recreate a column scientificName
df_index$scientificName <- sp_df_relative$scientificName

# split into Genus and Species
df_index <- df_index %>% separate(scientificName, c("Genus", "Species"))
sp_df_relative <- sp_df_relative %>% separate(scientificName, c("Genus", "Species"))
sp_df_merged <- sp_df_merged  %>% separate(scientificName, c("Genus", "Species"))

# Join taxonomy
sp_df_relative <- left_join(sp_df_relative, taxonomy, by = "Genus")
sp_df_merged <- left_join(sp_df_merged, taxonomy, by = "Genus")
df_index <- left_join(df_index, taxonomy, by = "Genus")

# recreate column scientificname
df_index$scientificName <- paste(df_index$Genus,df_index$Species)
sp_df_relative$scientificName <- paste(sp_df_relative$Genus,sp_df_relative$Species)
sp_df_merged$scientificName <- paste(sp_df_merged$Genus,sp_df_merged$Species)

###############################
# write species and genus table
###############################
# NB: reiterate the script with MIDORI annotated to compare with BC database annotation

# export table before turn into relative and filter at 0.0001% of cumulative relative abundance
write.csv(sp_df_merged, "01_Tables/12S_species_table_BC-db_2023.csv")

# export table after turn into relative
write.csv(sp_df_relative, "01_Tables/HS_12S_species_table_relative_BC-db.csv")

# export edna index table
#colnames(df_index) <- gsub(x = colnames(df_index), pattern = ".1", replacement = "", fixed = T)  
write.csv(df_index, "01_Tables/12S_species_table_eDNAIndex_MIDORI_2023.csv")

