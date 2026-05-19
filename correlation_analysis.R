#### this script will perform preliminary analysis of the antibody titre dataset ####
## the excel sheet used in the script was shared by Eva on 24/06/2024 ##
## the paper linked to the dataset is stored in the supplementary_material folder ##

#### importing relevant libraries ####
library(readxl)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(factoextra)
library(cluster)
library(RColorBrewer)
library(ggalluvial)
library(corrplot)
library(smplot2)

#### importing the dataset ####

#reading the dataset
#antibody_titre_info = read_xlsx("../data_Martin/R updated overview table lifelines cohort.xlsx")
ab_info_mn_neg = read_xlsx("input_files/!!correct Overview table Lifelines cohort MN- MN+ (1).xlsx", sheet = 2)
ab_info_mn_pos = read_xlsx("input_files/!!correct Overview table Lifelines cohort MN- MN+ (1).xlsx", sheet = 3)

#subsetting both
ab_info_mn_neg = ab_info_mn_neg[,c("participant number", "assessment", "age cohort", "Seks", "WIV-specific IgG (mean)", "IgG1", "IgG2", 
                                   "IgG4", "fcyrlllam(mean)", "% CD107+")]

ab_info_mn_pos = ab_info_mn_pos[,c("participant number", "assessment", "age cohort", "Seks", "WIV-specific IgG (mean)", "IgG1", "IgG2", 
                                   "IgG4", "fcyrlllam(mean)", "% CD107+")]

#changing column names for both
colnames(ab_info_mn_neg) = c("ID", "timepoint", "age", "gender", "total_IgG", "IgG1", "IgG2", "IgG4",  "FCG", "NK")
colnames(ab_info_mn_pos) = c("ID", "timepoint", "age", "gender", "total_IgG", "IgG1", "IgG2", "IgG4",  "FCG", "NK")

#### Question 1 - WIV vs Day 1 and Day 2 MN titres per strain ####
#### writing a function that will take both the MN negative and positive info and make correlation plots
#### for both timepoints, all age groups and also overall
ab_info_mn_correlate = function(ab_subset, extension){

  #### normalizing the antibody titres to log 10 ####
  
  ab_subset$IgG1 = log10(ab_subset$IgG1)
  ab_subset$IgG2 = log10(ab_subset$IgG2)
  ab_subset$IgG4 = log10(ab_subset$IgG4)
  ab_subset$FCG = log10(ab_subset$FCG)
  ab_subset$total_IgG = log10(ab_subset$total_IgG)
  
  ## IgG titres - all D1##
  #subsetting the D1 titers
  ab_subset_d1 = ab_subset[ab_subset$timepoint == 1,c("total_IgG", "IgG1", "IgG2", "IgG4", "FCG", "NK", "age")]
  
  #removing NA subset
  ab_subset_d1 = ab_subset_d1 %>% drop_na()
  
  #changing all Inf values to zero
  ab_subset_d1$IgG1[ab_subset_d1$IgG1 == -Inf] = 0
  ab_subset_d1$IgG2[ab_subset_d1$IgG2 == -Inf] = 0
  ab_subset_d1$IgG4[ab_subset_d1$IgG4 == -Inf] = 0
  ab_subset_d1$FCG[ab_subset_d1$FCG == -Inf] = 0
  ab_subset_d1$total_IgG[ab_subset_d1$total_IgG == -Inf] = 0 
  
  ## IgG titres - all D2##
  #subsetting the D1 titers
  ab_subset_d2 = ab_subset[ab_subset$timepoint == 2,c("total_IgG", "IgG1", "IgG2", "IgG4", "FCG", "NK", "age")]
  
  #removing NA subset
  ab_subset_d2 = ab_subset_d2 %>% drop_na()
  
  #changing all Inf values to zero
  ab_subset_d2$IgG1[ab_subset_d2$IgG1 == -Inf] = 0
  ab_subset_d2$IgG2[ab_subset_d2$IgG2 == -Inf] = 0
  ab_subset_d2$IgG4[ab_subset_d2$IgG4 == -Inf] = 0
  ab_subset_d2$FCG[ab_subset_d2$FCG == -Inf] = 0
  ab_subset_d2$total_IgG[ab_subset_d2$total_IgG == -Inf] = 0
  
  
  #### writing a function that will take in the ab day-specific data and making age specific correlation ####
  ab_correlate_plot = function(ab_info, day, age_group){
    
    #subsetting the ab_info
    ab_info = ab_info[ab_info$age == age_group,]
    
    #subsetting to remove age
    ab_info = subset(ab_info, select = -c(age))
    
    #making p value matrix
    ab_info_cor = cor.mtest(ab_info)
    
    #filename
    filename = paste("correlation_plots/", extension, "/", age_group, "_TP", day, ".jpeg", sep = "")
    
    #title
    title = paste(extension, " ; Correlation for : ", age_group, " & Day : ", day, sep = "")
    
    #jpeg device
    jpeg(filename = filename, width = 10, height = 10, res = 300, units = "in")
    
    #making correlation plot
    corrplot(cor(ab_info), method = "color", type = "upper", insig = "label_sig", tl.cex = 2, 
             p.mat = ab_info_cor$p, sig.level = c(.001, .01, .05), title = title, mar=c(0,0,2,0))
    
    #printing the plot
    #save_plot(plot = bp, filename = filename)
    
    #dev off
    dev.off()
    
    
  }
  
  ### making the correlaiton plots ###
  
  #tp1
  lapply(unique(ab_subset_d1$age), function(x){ab_correlate_plot(ab_info = ab_subset_d1, day = 1, age_group = x)})
  
  #tp2
  lapply(unique(ab_subset_d2$age), function(x){ab_correlate_plot(ab_info = ab_subset_d2, day = 2, age_group = x)})
  
  
  ### making the correlation plots - overall ###
  
  #tp1
  
  #making p value matrix
  ab_info_cor_d1 = cor.mtest(ab_subset_d1[,c(1,2,3,4,5,6)])
  
  #filename
  filename = paste("correlation_plots/", extension, "/overall_tp1.jpeg", sep = "")
  
  #title
  title = paste(extension, " ; Correlation TP1 Overall", sep = "")
  
  #jpeg device
  jpeg(filename = filename, width = 10, height = 10, res = 300, units = "in")
  
  #making correlation plot
  corrplot(cor(ab_subset_d1[,c(1,2,3,4,5,6)]), method = "color", type = "upper", insig = "label_sig", tl.cex = 2, 
           p.mat = ab_info_cor_d1$p, sig.level = c(.001, .01, .05), title = title, mar=c(0,0,2,0))
  
  #printing the plot
  #save_plot(plot = bp, filename = filename)
  
  #dev off
  dev.off()

  
  #tp2
  
  #making p value matrix
  ab_info_cor_d2 = cor.mtest(ab_subset_d2[,c(1,2,3,4,5,6)])
  
  #filename
  filename = paste("correlation_plots/", extension, "/overall_tp2.jpeg", sep = "")
  
  #title
  title = paste(extension, " ; Correlation TP2 Overall", sep = "")
  
  #jpeg device
  jpeg(filename = filename, width = 10, height = 10, res = 300, units = "in")
  
  #making correlation plot
  corrplot(cor(ab_subset_d2[,c(1,2,3,4,5,6)]), method = "color", type = "upper", insig = "label_sig", tl.cex = 2, 
           p.mat = ab_info_cor_d2$p, sig.level = c(.001, .01, .05), title = title, mar=c(0,0,2,0))
  
  #printing the plot
  #save_plot(plot = bp, filename = filename)
  
  #dev off
  dev.off()
  
  
}

## giving call to the function that will make the plots for me directly for both MN neg and MN pos
ab_info_mn_correlate(ab_subset = ab_info_mn_neg, extension = "MN_negative")
ab_info_mn_correlate(ab_subset = ab_info_mn_pos, extension = "MN_positive")


