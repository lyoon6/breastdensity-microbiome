#################################################
# The association between breast density and    #
#   gut microbiota composition at 2 years post  #
#    menarche: A cross-sectional study of       #
#    adolescents in Santiago, Chile             #
#                                               #
# MaAsLin Analyses                              #
#                                               #
# Lara Yoon                                     #
#    PhD Student, UCLA Dept of Epidemiology     #
#    lyoon6@ucla.edu                            #
#                                               #
# Last updated 03/16/2021                       #
# Created with R version 4.0.4                  #
#################################################

## Note: This script uses data file generated in the 'BD-Micro_DiversityAnalyses_Script.R' file


#############################################
##                                         ##
##       Package loading and options       ##
##                                         ##
#############################################

# Refer to individual package documentation for installation if necessary
pacman::p_load("devtools","here", "BiocManager", # tools 
               "assertr", "tidyverse", "reshape2", "janitor", "data.table", # data management
               "dada2", "phyloseq", "qiime2R", "ggplot2", "DESeq2", "GUniFrac", "vegan", "microbiome", "Maaslin2", "boot", "ade4", # analysis
               "scales", "grid", "ggsci", "ggpubr", "viridis", "picante", "knitr", "gridExtra" # data viz
)  

pacman::p_loaded()

# Use standard notation
options(scipen=999)


#############################################
##                                         ##
##       Loading the data sets             ##
##                                         ##
#############################################

# ASV Table (samples are rows; taxa are columns)
otutab <- read_delim("Data/FinalData/OTUtable_filt.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
otutab <- column_to_rownames(otutab, var="SampleID")

# Tax Table 
taxtab <- read_delim("Data/FinalData/TAXtable_filt.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
taxtab <- column_to_rownames(taxtab, var="Taxonomy")
taxtab <- t(taxtab)

# Sample data 
samptab <- read_csv("Data/FinalData/Meta_after_filter_2.csv")
samptab <- column_to_rownames(samptab, var="SampleID")

# Convert imported files to phyloseq objects 
otu <- otu_table(otutab, taxa_are_rows=FALSE)
tax <- tax_table(taxtab)
met <- sample_data(samptab)
phyobj <- merge_phyloseq(otu, tax, met)



#############################################
##                                         ##
##      PREP FOR Maaslin                   ##
##                                         ##
#############################################

# Agglomerate to the genus level
phyobjsubglom <- tax_glom(phyobj, taxrank="Genus", NArm=TRUE)
ntaxa(phyobj); ntaxa(phyobjsubglom) #1600; 269
# Note that if we use NArm=TRUE, we get 223 taxa as opposed to 269 WITH naRM=false

# Define individual objects as data frames
taxtable<-as.data.frame(tax_table(phyobjsubglom))
asvtable<-as.data.frame(otu_table(phyobjsubglom))
metatable<-meta(phyobjsubglom)

# Merge tax table and otu table 
names(taxtable)
taxtable$taxa <- paste0(taxtable$Phylum, "; ", taxtable$Class, "; ", taxtable$Order, "; ", taxtable$Family, "; ", taxtable$Genus)
taxtabsub <- taxtable %>%  select(c("taxa"))
asvtable <- t(asvtable)
asvnames <-  merge(taxtabsub, asvtable, by="row.names")
rownames(asvnames) <- asvnames$taxa
asvnames <- select(asvnames, -c("Row.names", "taxa"))
asvnamest <- as.data.frame(t(asvnames)) #samples are rows


###########################################  
# Sensitivity analysis: subset to top 100 taxa 
# Do not run this first time around
###########################################  

# Prune to top 100 genera & redo 
poglomtop <- prune_taxa(names(sort(taxa_sums(phyobjsubglom), TRUE))[1:100], phyobjsubglom)
ntaxa(poglomtop) #100 taxa 

# define individual objects as data frames
taxtable<-as.data.frame(tax_table(poglomtop)) #get taxonomy
asvtable<-as.data.frame(otu_table(poglomtop))
metatable<-meta(poglomtop) #get asvs

# merge tax table and otu table 
names(taxtable)
taxtable$taxa <- paste0(taxtable$Phylum, "; ", taxtable$Class, "; ", taxtable$Order, "; ", taxtable$Family, "; ", taxtable$Genus)
taxtabsub <- taxtable %>%  select(c("taxa"))
asvtable <- t(asvtable)
asvnames <-  merge(taxtabsub, asvtable, by="row.names")
rownames(asvnames) <- asvnames$taxa
asvnames <- select(asvnames, -c("Row.names", "taxa"))
asvnamest <- as.data.frame(t(asvnames)) #samples are rows



#############################################
##                                         ##
##       MAASLIN                           ##
##                                         ##
#############################################


###################################
## Treating outcomes as categories
###################################

# pFGV 
############################################
# MODEL 1
maas1_pfgv <- Maaslin2(
  input_data=asvnamest, 
  input_metadata=metatable, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("pFGV_tert"), 
  reference=c("pFGV_tert,T1"), 
  min_prevalence=0)

# MODEL 2
maas2_pfgv_fat2 <- Maaslin2(
  input_data=asvnamest, 
  input_metadata=metatable, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("pFGV_tert", "fat2Y", "age_2pm"), 
  reference=c("pFGV_tert,T1"), 
  min_prevalence=0)

# MODEL 2, no fat
maas2_pfgv_fat2 <- Maaslin2(
  input_data=asvnamest, 
  input_metadata=metatable, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("pFGV_tert", "age_2pm"), 
  reference=c("pFGV_tert,T1"), 
  min_prevalence=0)

# MODEL 3
maas3_pfgv <- Maaslin2(
  input_data=asvnamest, 
  input_metadata=metatable, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("pFGV_tert","fat2Y","abx6mo","meducatbin","AvgCalQ", "tv", "ethnicity", "birth_mode", "age_2pm", "bfeedcat"), 
  reference=c("pFGV_tert,T1","meducatbin,1","AvgCalQ,1","agemencat,1", "tv,1", "ethnicity,1", "birth_mode,1", "bfeedcat,1"), 
  min_prevalence=0)

# MODEL 3, no fat
maas3_pfgv <- Maaslin2(
  input_data=asvnamest, 
  input_metadata=metatable, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("pFGV_tert", "abx6mo","meducatbin","AvgCalQ", "tv", "ethnicity", "birth_mode", "age_2pm", "bfeedcat"), 
  reference=c("pFGV_tert,T1","meducatbin,1","AvgCalQ,1","agemencat,1", "tv,1", "ethnicity,1", "birth_mode,1", "bfeedcat,1"), 
  min_prevalence=0)


# aFGV 
############################################
maas1_afgv <- Maaslin2(
  input_data=asvnamest, 
  input_metadata=metatable, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("aFGV_tert"), 
  reference=c("aFGV_tert,T1"), 
  min_prevalence=0)

# MODEL 2
maas2_afgv <- Maaslin2(
  input_data=asvnamest, 
  input_metadata=metatable, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("aFGV_tert","fat2Y", "age_2pm"), 
  reference=c("aFGV_tert,T1"), 
  min_prevalence=0)

# MODEL 3
maas3_afgv <- Maaslin2(
  input_data=asvnamest, 
  input_metadata=metatable, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("aFGV_tert","fat2Y","abx6mo","meducatbin","AvgCalQ", "tv", "ethnicity", "birth_mode", "age_2pm", "bfeedcat"), 
  reference=c("aFGV_tert,T1","meducatbin,1","AvgCalQ,1","agemencat,1", "tv,1", "ethnicity,1", "birth_mode,1", "bfeedcat,1"), 
  min_prevalence=0)



###################################
## Treating outcomes as continuous
###################################

# MODEL 3
maas3_pfgv <- Maaslin2(
  input_data=asvnamest, 
  input_metadata=metatable, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("DXA_FGVpercent","fat2Y","abx6mo","meducatbin","AvgCalQ", "tv", "ethnicity", "birth_mode", "age_2pm", "bfeedcat"), 
  reference=c("meducatbin,1","AvgCalQ,1","agemencat,1", "tv,1", "ethnicity,1", "birth_mode,1", "bfeedcat,1"), 
  min_prevalence=0)

# MODEL 3
maas3_afgv <- Maaslin2(
  input_data=asvnamest, 
  input_metadata=metatable, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("DXA_FGVabsolute","fat2Y","abx6mo","meducatbin","AvgCalQ", "tv", "ethnicity", "birth_mode", "age_2pm", "bfeedcat"), 
  reference=c("meducatbin,1","AvgCalQ,1","agemencat,1", "tv,1", "ethnicity,1", "birth_mode,1", "bfeedcat,1"), 
  min_prevalence=0)




