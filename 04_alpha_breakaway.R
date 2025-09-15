
##########################################################
# load envt and packages
##########################################################
setwd("C:/Users/loicj/OneDrive - UBC/01_Fjord_eDNA/05_Scratch/NEW")
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(vegan)
library(phyloseq)
library(ggpubr)
library(ggpmisc)
library(breakaway)

########################
#import data 
########################
# https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
# https://adw96.github.io/breakaway/articles/intro-diversity-estimation.html#species-richness-estimation-with-breakaway-1

df <- read.csv("01_Tables/12S_species_table_BC-db_2023.csv", sep = ",")
rownames(df) <- df$scientificName
df <- dplyr::select(df, contains(c("2023")))
metadata <- read.csv("01_Metadata_site.csv", sep = ",")

# build a frequency table
frequencytablelist <- build_frequency_count_tables(df)

# calculate richness estimate using breakaway
index_div <- breakaway(df)

# basic breakaway plot
plot(index_div)

# extract data
index_div_df <- summary(index_div)

# remove row with NAs
index_div_df <- na.omit(index_div_df)

# rename table
colnames(index_div_df) <- c("estimate", "error", "lower", "upper", "Site_ID_Depth_Season", "name", "model")

# join metadata
index_div_df <- left_join(index_div_df, metadata, by = "Site_ID_Depth_Season")

##################################
# Plot global
##################################
colnames(index_div_df)
unique(index_div_df$Season)

# reorder x axis on boxplot
index_div_df$Season <- factor(index_div_df$Season , levels=c("Winter", "Spring", "Summer", "Fall"))
index_div_df$Location.ID <- factor(index_div_df$Location.ID , levels=c("Indian", "Howe_Sound", "Toba", "Bute", "Rivers"))
index_div_df$Feature <- factor(index_div_df$Feature , levels=c("Surface", "10m", "100m", "bottom"))

#########################
# Boxplot per Depth
#########################
depth_pal <- c("orange", "red","#8c6bb1","black")

# boxplot
richness_boxplot_depth <- ggplot(index_div_df, aes(x=Feature, y=estimate , fill = Feature)) + 
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  theme_bw()+
  labs(x = "Depth", 
       y = "Richness estimate")+
  #ylim(0, 100)+
  scale_fill_manual(values=c("Surface"="orange", 
                             "10m"="red",
                             "100m"="#8c6bb1",
                             "bottom"="black"))+
  theme(axis.text.x=element_blank())
richness_boxplot_depth

# Comparison of each group against each others
# http://www.sthda.com/english/wiki/comparing-means-in-r
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

# one way anova test: in this test, indicate that some of the group means are different,
# but we don't know which pairs of groups are different
res.aov_depth <- aov(estimate ~ Location.ID, data = index_div_df)
summary(res.aov_depth)

# if the One-Way ANOVA is significant, we can compute a Tukey test to perform multiple pairwise-comparisons
# between the means of groups:
# diff show the difference between mean of the two groups, with padj is the the p-value (sign if < 0.05)
TukeyHSD(res.aov_depth)

# 1. test homogeneity of variances: if there is relatinship between residuals and fitted-value (mean of each groups),
# data do not fit the homogeneity of variance condition
plot(res.aov_depth, 1)

# 2. Check the normarlity assumption: if the normality is respected, the plot should show
# approximately a straight line
plot(res.aov_depth, 2)

# A shapiro test can support normality test
# Extract the residuals
aov_depth_residuals <- residuals(object = res.aov_depth)

# Run Shapiro-Wilk test
# if p-value < 0.05, reject the null hypothesis that the population is normally distributed
# if pvalue > 0.05, ok
shapiro.test(x = aov_depth_residuals)

# When assumption 1 (variance) or 2 (normality) are not met, we can use a Kruskal-Wallis test (non parametric)
# p-value < 0.05, ==> significant differences between the treatment groups.
kruskal.test(estimate ~ Feature, data = index_div_df)

# to know which pairs of group are different, use pairwise.wilcox.test()
# only treatment with pvalue < 0.05 are significantly differents
pairwise.wilcox.test(index_div_df$estimate, index_div_df$Feature,
                     p.adjust.method = "BH")

# select the significant comparaison from pairwise.wilcox.test and add them to the plot
sign_comparisons_depth <- list( c("Surface", "10m"), c("bottom", "10m"))

# boxplot_depth
boxplot_depth <- richness_boxplot_depth + stat_compare_means(comparisons = sign_comparisons_depth)
boxplot_depth


# OR Comparison of each group against the base-mean
# default is wilcoxon test
compare_means(estimate ~ Feature, index_div_df, ref.group = ".all.", method = "wilcox.test")

boxplot_depth <- ggboxplot(index_div_df, x = "Feature", y = "estimate", add = "jitter",
                           color = "Feature", group="Feature",
                           palette = depth_pal, 
                           legend.title="Depth")+
  geom_hline(yintercept = mean(index_div_df$estimate), linetype = 2) 
boxplot_depth

boxplot_depth <- boxplot_depth + stat_compare_means(label = "p.signif",
                                                    ref.group = ".all.", # compare each boxplot with means
                                                    method = "wilcox.test", 
                                                    label.y = 50) # select label position on y axis
boxplot_depth

# and multiple pairwise comparison between group
pairwise_wilcox_depth <- pairwise.wilcox.test(index_div_df$estimate, index_div_df$Feature,
                     p.adjust.method = "BH")
#write.csv(pairwise_wilcox_depth$p.value, "08_Alpha_div/pairwise_wilcox_depth.csv")

# Boxplot per season
richness_boxplot_season <- ggplot(index_div_df, aes(x=Season, y=estimate , fill = Season)) + 
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  theme_bw()+
  labs(x = "Season", 
       y = "Richness estimate")+
   #ylim(0, 100)+
   scale_fill_manual(values=c("Fall"="#F8766D", 
                              "Spring"="#7CAE00",
                              "Summer"="#007643",
                              "Winter"="#00BFC4"))+
  theme(axis.text.x=element_blank())
richness_boxplot_season

season_pal <- c("#00BFC4", "#7CAE00","#007643","#F8766D")

# compare against base-mean
compare_means(estimate ~ Season, index_div_df, ref.group = ".all.", method = "wilcox.test")

boxplot_season <- ggboxplot(index_div_df, x = "Season", y = "estimate", add = "jitter",
                           color = "Season", group="Season",
                           palette = season_pal, 
                           legend.title="season")+
  geom_hline(yintercept = mean(index_div_df$estimate), linetype = 2) 
boxplot_season

boxplot_season <- boxplot_season + stat_compare_means(label = "p.signif",
                                                    ref.group = ".all.", # compare each boxplot with means
                                                    method = "wilcox.test", 
                                                    label.y = 50) # select label position on y axis
boxplot_season

# and multiple pairwise comparison between group
pairwise_wilcox_season <- pairwise.wilcox.test(index_div_df$estimate, index_div_df$Season,
                                              p.adjust.method = "BH")
pairwise_wilcox_season$p.value
#write.csv(pairwise_wilcox_season$p.value, "08_Alpha_div/pairwise_wilcox_season.csv")


# boxplot per fjord
richness_boxplot_fjord <- ggplot(index_div_df, aes(x=Location.ID, y=estimate , fill = Location.ID)) + 
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  theme_bw()+
  labs(x = "Fjord", 
       y = "Richness estimate")+
  #ylim(0, 100)+
  scale_fill_manual(values=c(c("Bute" = "#7fc97f",
                               "Howe_Sound"="#beaed4", 
                               "Indian"="#fdc086",
                               "Rivers"="#386cb0", 
                               "Toba"="#f0027f")))+
  theme(axis.text.x=element_blank())
richness_boxplot_fjord

fjord_pal <- c("#fdc086","#beaed4","#f0027f","#7fc97f","#386cb0")

# Comparison of each group against base-mean
compare_means(estimate ~ Location.ID, index_div_df, ref.group = ".all.", method = "wilcox.test")

boxplot_fjord <- ggboxplot(index_div_df, x = "Location.ID", y = "estimate", add = "jitter",
                           color = "Location.ID", group="Location.ID",
                           palette = fjord_pal, 
                           legend.title="fjord")+
  geom_hline(yintercept = mean(index_div_df$estimate), linetype = 2) 
boxplot_fjord

boxplot_fjord <- boxplot_fjord + stat_compare_means(label = "p.signif",
                                                    ref.group = ".all.", # compare each boxplot with means
                                                    method = "wilcox.test", 
                                                    label.y = 50) # select label position on y axis
boxplot_fjord

# and multiple pairwise comparison between group
pairwise_wilcox_fjord <- pairwise.wilcox.test(index_div_df$estimate, index_div_df$Location.ID,
                                               p.adjust.method = "BH")
pairwise_wilcox_fjord$p.value
#write.csv(pairwise_wilcox_fjord$p.value, "08_Alpha_div/pairwise_wilcox_fjord.csv")


# plot them all
ggarrange(boxplot_depth, boxplot_season, boxplot_fjord,
          ncol = 1, nrow = 3,
          common.legend = ,
          legend = "left",
          align = "hv")



##################################
# Plot diversity boxplot per fjords
##################################
colnames(index_div_df)
unique(index_div_df$Season)
# reorder x axis on boxplot

# Fjords filter
index_div_df_BU <- filter(index_div_df, Location.ID == "Bute")
index_div_df_TO <- filter(index_div_df, Location.ID == "Toba")
index_div_df_RI <- filter(index_div_df, Location.ID == "Rivers")
index_div_df_HS <- filter(index_div_df, Location.ID == "Howe_Sound")
index_div_df_IA <- filter(index_div_df, Location.ID == "Indian")

#########################
# Boxplot per Depth
#########################

# Boxplot per season
richness_boxplot_RI <- ggplot(index_div_df_RI, aes(x=Season, y=estimate , fill = Feature)) + 
  geom_boxplot()+
  #geom_jitter(width = 0.1)+
  theme_bw()+
  labs(y = "Richness estimate")+
  #ylim(0, 100)+
  scale_fill_manual(values=c("Surface"="red", 
                             "10m"="orange",
                             "100m"="#8c6bb1",
                             "bottom"="black"))+
  ylim(0, 60)+
  ggtitle("Rivers")
richness_boxplot_RI

richness_boxplot_BU <- ggplot(index_div_df_BU, aes(x=Season, y=estimate , fill = Feature)) + 
  geom_boxplot()+
  #geom_jitter(width = 0.1)+
  theme_bw()+
  labs(y = "Richness estimate")+
  #ylim(0, 100)+
  scale_fill_manual(values=c("Surface"="red", 
                             "10m"="orange",
                             "100m"="#8c6bb1",
                             "bottom"="black"))+
  ylim(0, 60)+
  ggtitle("Bute")
richness_boxplot_BU

richness_boxplot_TO <- ggplot(index_div_df_TO, aes(x=Season, y=estimate , fill = Feature)) + 
  geom_boxplot()+
  #geom_jitter(width = 0.1)+
  theme_bw()+
  labs(y = "Richness estimate")+
  #ylim(0, 100)+
  scale_fill_manual(values=c("Surface"="red", 
                             "10m"="orange",
                             "100m"="#8c6bb1",
                             "bottom"="black"))+
  ylim(0, 60)+
  ggtitle("Toba")
richness_boxplot_TO

richness_boxplot_HS <- ggplot(index_div_df_HS, aes(x=Season, y=estimate , fill = Feature)) + 
  geom_boxplot()+
  #geom_jitter(width = 0.1)+
  theme_bw()+
  labs(y = "Richness estimate")+
  #ylim(0, 100)+
  scale_fill_manual(values=c("Surface"="red", 
                             "10m"="orange",
                             "100m"="#8c6bb1",
                             "bottom"="black"))+
  ylim(0, 60)+
  ggtitle("Howe Sound")
richness_boxplot_HS

richness_boxplot_IA <- ggplot(index_div_df_IA, aes(x=Season, y=estimate , fill = Feature)) + 
  geom_boxplot()+
  theme_bw()+
  labs(y = "Richness estimate")+
  #ylim(0, 100)+
  scale_fill_manual(values=c("Surface"="red", 
                             "10m"="orange",
                             "100m"="#8c6bb1",
                             "bottom"="black"))+
  ylim(0, 60)+
  ggtitle("Indian Arm")
richness_boxplot_IA


# plot them all
library(ggpubr)

ggarrange(richness_boxplot_RI, richness_boxplot_BU,
          richness_boxplot_TO, richness_boxplot_HS,
          richness_boxplot_IA,
          nrow = 3, ncol = 2,
          common.legend = T)
