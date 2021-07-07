#################################################
# The association between breast density and    #
#   gut microbiota composition at 2 years post  #
#    menarche: A cross-sectional study of       #
#    adolescents in Santiago, Chile             #
#                                               #
# PICRUSt Analyses                              #
#                                               #
# Lara Yoon                                     #
#    PhD Student, UCLA Dept of Epidemiology     #
#    lyoon6@ucla.edu                            #
#                                               #
# Last updated 03/16/2021                       #
# Created with R version 4.0.4                  #
#################################################


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
##         THIS STEP FOR EXPORTING         ##
##     PROCESSED PHYLOSEQ FILE TO QIIME    ##
##                                         ##  
#############################################

##  Using phyloseq object generated in 'BD-Micro_DiversityAnalyses_Script.R'
subpo2

## Getting .BIOM file    
otu <- as(otu_table(subpo2), "matrix") # taxa are rows 
otu_biom <- make_biom(data=otu)
write_biom(otu_biom, "<REPLACE WITH FILEPATH>") 



#############################################
##                                         ##
##         THIS STEP FOR IMPORTING         ##
##     PICRUST FILES FOR DOWNSTREAMS       ##
##                                         ##  
#############################################

# Note: Use options skip=1 and comment.char= to read in the .tsv file 

# Import Pathway Abundance file from PICRUSt output
metacyc <- read.table(file="<REPLACE WITH FILEPATH>", sep = '\t', header = TRUE, check.names = FALSE)
# Samples are columns; transpose 
rownames(metacyc) <- metacyc$'OTU ID'
metacyc <- select(metacyc, -c('OTU ID'))
metacyct <- as.data.frame(t(metacyc)) #samples are rows


# Import EC metagnome file from PICRUSt output (note physical deletion of # on OTU row)
ec <- read.table(file="<REPLACE WITH FILEPATH>", sep = '\t', skip=1, header = TRUE, comment.char = "#",  check.names = FALSE)
# Samples are columns; transpose 
rownames(ec) <- ec$'OTU ID'
ec <- select(ec, -c('OTU ID'))
ect <- as.data.frame(t(ec)) #samples are rows

# Import KO metagnome file from PICRUSt output
ko <- read.table(file="<REPLACE WITH FILEPATH>", sep = '\t', skip=1, header = TRUE, comment.char = "#",  check.names = FALSE)
# Samples are columns; transpose 
rownames(ko) <- ko$'OTU ID'
ko <- select(ko, -c('OTU ID'))
kot <- as.data.frame(t(ko)) #samples are rows

# Import metadata 
metapi <- read_csv("<REPLACE WITH FILEPATH>")  # Samples are rows 
# Make 'Sample ID' the rowname 
rownames(metapi) <- metapi$SampleID 


# Import quick description of each functional category using PICRUSt mapping files 
ec_desc <- read.table(file="<REPLACE WITH FILEPATH>", sep = '\t', check.names = FALSE)
ko_desc <- read.table(file="<REPLACE WITH FILEPATH>", sep = '\t', check.names = FALSE)
metacyc_desc <- read.table(file="<REPLACE WITH FILEPATH>", sep = '\t', check.names = FALSE)

#############################################
##                                         ##
##         CREATING PHYLOSEQ OBJECTS       ##
##                                         ##  
#############################################

FEATmetacyc <- otu_table(metacyct, taxa_are_rows = FALSE)
view(FEATmetacyc)
FEATec <- otu_table(ect, taxa_are_rows = FALSE)
FEATko <- otu_table(kot, taxa_are_rows = FALSE)
META <- sample_data(metapi)
rownames(META) <- META$SampleID 

POmetacyc <- phyloseq(FEATmetacyc,META)
summarize_phyloseq(POmetacyc)

POec <- phyloseq(FEATec, META)
summarize_phyloseq(POec)

POko <- phyloseq(FEATko, META)
summarize_phyloseq(POec)


#############################################
##                                         ##
##        PREPROCESSING                    ##
##                                         ##  
#############################################

logt  = transform_sample_counts(POmetacyc, function(x) log(1 + x) )

# PCoA with bray to look for outliers (for both pfgv and afgv)
out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues
plot_ordination(logt, out.pcoa.logt, type = "samples", 
                title='PCoA (bray) to look for outliers') + 
  coord_fixed(sqrt(evals[2] / evals[1])) 

rel_abund <- t(apply(otu_table(POmetacyc), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram") +
  xlab("Relative abundance")


#############################################
##                                         ##
##         BETA DIVERSITY                  ##
##                                         ##  
#############################################

########### 
# MetaCyc #
###########

# log transforming the abundance data for the analyses
qplot(log10(rowSums(otu_table(POmetacyc)))) + xlab("Logged counts-per-sample")
POmetacyc.log <- transform_sample_counts(POmetacyc, function(x) log(1 + x))

#plot PCoA using Bray-Curtis as distance metric 
pcoa_bray_ord <- ordinate(POmetacyc.log, method="PCoA", distance="bray")

# TOtal distance captured 
plot_scree(pcoa_bray_ord,'Scree plot, Bray/PCoA')
# How much variation do the first two axes (ones we will plot) explain?
(100*sum(pcoa_bray_ord$values$Relative_eig[1:2])) #53%

#plotting 
pcoa_bray_pfgv <- plot_ordination(POmetacyc.log, pcoa_bray_ord, color="pFGV_tert") + theme(aspect.ratio=1) + 
  theme_bw()+
  scale_color_viridis(discrete=TRUE)+ 
  labs(color="Percent FGV (Tercile)")+
  stat_ellipse(size=1)+
  theme(legend.position="bottom")
pcoa_bray_pfgv

pcoa_bray_afgv <- plot_ordination(POmetacyc.log, pcoa_bray_ord, color="aFGV_tert") + theme(aspect.ratio=1) + 
  theme_bw()+
  scale_color_viridis(discrete=TRUE)+ 
  labs(color="Absolute FGV (Tercile)")+
  stat_ellipse(size=1)+
  theme(legend.position="bottom")
pcoa_bray_afgv


# group together [FIGURE 4]
ggarrange(pcoa_bray_pfgv, pcoa_bray_afgv,
          labels=c("A", "B"),
          ncol=2, nrow=1, common.legend = FALSE, legend = "bottom")


##########################
## PERMANOVA
##########################

# calculate distances 
bray_dist <- phyloseq::distance(POmetacyc.log, method="bray")

#PERMANOVA to test whether BD tercile differ significantly from each other
# model 1
permbray12 <- adonis2(bray_dist ~ sample_data(POmetacyc.log)$pFGV_tert)

# model 2
permbray22 <- adonis2(bray_dist ~ sample_data(POmetacyc.log)$pFGV_tert 
                      + sample_data(POmetacyc.log)$fatcat)

# model 3 [TABLE 3]
permbray32 <- adonis2(bray_dist ~ sample_data(POmetacyc.log)$pFGV_tert 
                      + sample_data(POmetacyc.log)$fatcat
                      + sample_data(POmetacyc.log)$age_2pm
                      + sample_data(POmetacyc.log)$abx6mo
                      + sample_data(POmetacyc.log)$birth_mode
                      + sample_data(POmetacyc.log)$bfeedcat
                      + sample_data(POmetacyc.log)$AvgCalQ
                      + sample_data(POmetacyc.log)$meducatbin 
                      + sample_data(POmetacyc.log)$ethnicity )

# model 4
permbray42 <- adonis2(bray_dist ~ sample_data(POmetacyc.log)$aFGV_tert)

# model 5
permbray52 <- adonis2(bray_dist ~ sample_data(POmetacyc.log)$aFGV_tert 
                      + sample_data(POmetacyc.log)$fatcat)

# model 6 [TABLE 3]
permbray62 <- adonis2(bray_dist ~ sample_data(POmetacyc.log)$aFGV_tert 
                      + sample_data(POmetacyc.log)$fatcat
                      + sample_data(POmetacyc.log)$age_2pm
                      + sample_data(POmetacyc.log)$abx6mo
                      + sample_data(POmetacyc.log)$birth_mode
                      + sample_data(POmetacyc.log)$bfeedcat
                      + sample_data(POmetacyc.log)$AvgCalQ
                      + sample_data(POmetacyc.log)$meducatbin 
                      + sample_data(POmetacyc.log)$ethnicity )


# homogeneity of dispersion 
beta_p <- betadisper(bray_dist, sample_data(POmetacyc.log)$pFGV_tert)
permutest(beta_p)
beta_a <- betadisper(bray_dist, sample_data(POmetacyc.log)$aFGV_tert)
permutest(beta_a)
# not significant, meaning we cannot reject the null hypothesis that our groups have the same dispersions.



#############################################
##                                         ##
##      MaAsLin2                           ##
##                                         ##
#############################################

# Inputs: Feature table and metadata 
featuretab <- as.data.frame(metacyct)
metadata <- as.data.frame(sample_data(POmetacyc))

# pFGV 
############################################
# MODEL 1
maas1_pfgv <- Maaslin2(
  input_data=featuretab, 
  input_metadata=metadata, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("pFGV_tert"), 
  reference=c("pFGV_tert,T1"), 
  min_prevalence=0)


# MODEL 2
maas2_pfgv_fat2 <- Maaslin2(
  input_data=featuretab, 
  input_metadata=metadata, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("pFGV_tert", "fat2Y", "age_2pm"), 
  reference=c("pFGV_tert,T1"), 
  min_prevalence=0)

# MODEL 3
maas3_pfgv <- Maaslin2(
  input_data=featuretab, 
  input_metadata=metadata, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("pFGV_tert","fat2Y","abx6mo","meducatbin","AvgCalQ", "tv", "ethnicity", "birth_mode", "age_2pm", "bfeedcat"), 
  reference=c("pFGV_tert,T1","meducatbin,1","AvgCalQ,1","agemencat,1", "tv,1", "ethnicity,1", "birth_mode,1", "bfeedcat,1"), 
  min_prevalence=0)

# aFGV 
############################################
# MODEL 1
maas1_afgv <- Maaslin2(
  input_data=featuretab, 
  input_metadata=metadata, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("aFGV_tert"), 
  reference=c("aFGV_tert,T1"), 
  min_prevalence=0)


# MODEL 2
maas2_afgv_fat2 <- Maaslin2(
  input_data=featuretab, 
  input_metadata=metadata, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("aFGV_tert", "fat2Y", "age_2pm"), 
  reference=c("aFGV_tert,T1"), 
  min_prevalence=0)

# MODEL 3
maas3_afgv <- Maaslin2(
  input_data=featuretab, 
  input_metadata=metadata, 
  output= "<REPLACE WITH FILEPATH>", 
  fixed_effects=c("aFGV_tert","fat2Y","abx6mo","meducatbin","AvgCalQ", "tv", "ethnicity", "birth_mode", "age_2pm", "bfeedcat"), 
  reference=c("aFGV_tert,T1","meducatbin,1","AvgCalQ,1","agemencat,1", "tv,1", "ethnicity,1", "birth_mode,1", "bfeedcat,1"), 
  min_prevalence=0)
