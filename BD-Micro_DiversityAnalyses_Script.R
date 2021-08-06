#################################################
# The association between breast density and    #
#   gut microbiota composition at 2 years post  #
#    menarche: A cross-sectional study of       #
#    adolescents in Santiago, Chile             #
#                                               #
# Diversity Analyses                            #
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
library(pacman)
pacman::p_load("devtools","here", "BiocManager", "ape","biomformat", "Biostrings", # tools 
               "assertr", "tidyverse", "reshape2", "janitor", "data.table",  # data management
               "dada2", "phyloseq", "qiime2R", "ggplot2", "DESeq2", "GUniFrac", "vegan", "microbiome", "Maaslin2", "boot", "ade4", "car", # analysis
               "scales", "grid", "ggsci", "ggpubr", "viridis", "picante", "knitr", "gridExtra") # data viz


# Use standard notation
options(scipen=999)


#############################################
##                                         ##
##       Loading the data sets             ##
##                                         ##
#############################################

# Files for import 
feat <- here("Microbiome/Input", "tablesilva.qza")
taxo <- here("Microbiome/Input", "silva_97_taxonomy.qza" )
treeroot <- here("Microbiome/Input", "rooted_tree_masked_alignment.qza")
meta <- here("Data/FinalData/MetadataImputed.txt")   # Use metadata after dropping controls and samples without breast density

# Convert imported  files to phyloseq object 
phyobj<-qza_to_phyloseq(features=feat,
                        taxonomy=taxo,
                        tree=treeroot,
                        metadata=meta)

# Define individual objects
taxtable<-as.data.frame(tax_table(phyobj)) 
phytree<-phy_tree(phyobj) 
metadata<-sample_data(phyobj) 
asvs<-otu_table(phyobj) 


# OTU Table
phyobjasv <- as(otu_table(phyobj), "matrix")
if(taxa_are_rows(phyobj)){phyobjasv <- t(phyobjasv)}
phyobjasvdf <- as.data.frame(phyobjasv)
phyobjasvdf <- rownames_to_column(phyobjasvdf, "SampleID")
write.csv(phyobjasvdf, file="Data/FinalData/ASVtable.csv", 
            sep="\t", row.names=FALSE)

# Tax Table 
phyobjtax <- as(tax_table(phyobj), "matrix")
if(taxa_are_rows(phyobj)){phyobjtax <- t(phyobjtax)}
phyobjtaxdf <- as.data.frame(phyobjtax)
phyobjtaxdf <- rownames_to_column(phyobjtaxdf, "Taxonomy")
write.table(phyobjtaxdf, file="Data/FinalData/TAXtable.csv", 
            sep="\t", row.names=FALSE)




#############################################
##                                         ##
##       Summarize data                    ##
##                                         ##
#############################################

# Summarize phyloseq object
summarize_phyloseq(phyobj)

# Sum number of distinct ASVs observed by sample (different than 0)
no_asvs <- apply( otu_table(phyobj), 2, function(x) sum(x != 0) )
no_asvs_df <- data.frame("Sample" = names(no_asvs),
                         "Type_alpha" = rep("Observed", length(no_asvs)),
                         "Alpha_div" = no_asvs, 
                         row.names = names(no_asvs))
colnames(no_asvs_df) <- c("Samples", "Type of alpha-diversity", "Alpha-diversity measure") 



## Creating a function to summarize taxa (Function author:  Paul J. McMurdie )
fast_melt = function(phyobj){
  # supports "naked" otu_table as `phyobj` input.
  otutab = as(otu_table(phyobj), "matrix")
  if(!taxa_are_rows(phyobj)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(phyobj, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(phyobj, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

summarize_taxa = function(phyobj, Rank, GroupBy = NULL){
  Rank <- Rank[1]
  if(!Rank %in% rank_names(phyobj)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(phyobj), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){
    GroupBy <- GroupBy[1]
    if(!GroupBy %in% sample_variables(phyobj)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(phyobj), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt
  mdt = fast_melt(phyobj)
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    sdt = data.table(SampleID = sample_names(phyobj),
                     var1 = get_variable(phyobj, GroupBy))
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  summarydt = mdt[, list(meanRA = mean(RelativeAbundance),
                         sdRA = sd(RelativeAbundance),
                         minRA = min(RelativeAbundance),
                         maxRA = max(RelativeAbundance)),
                  by = c(Rank, GroupBy)]
  return(summarydt)
}

plot_taxa_summary = function(phyobj, Rank, GroupBy = NULL){
  # Get taxa summary table 
  dt1 = summarize_taxa(phyobj, Rank = Rank, GroupBy = GroupBy)
  # Set factor appropriately for plotting
  RankCol = which(colnames(dt1) == Rank)
  setorder(dt1, -meanRA)
  dt1[, RankFac := factor(dt1[[Rank]], 
                          levels = rev(dt1[[Rank]]))]
  dt1[, ebarMax := max(c(0, min(meanRA + sdRA))), by = eval(Rank)]
  dt1[, ebarMin := max(c(0, min(meanRA - sdRA))), by = eval(Rank)]
  # Set zeroes to one-tenth the smallest value
  ebarMinFloor = dt1[(ebarMin > 0), min(ebarMin)]
  ebarMinFloor <- ebarMinFloor / 10
  dt1[(ebarMin == 0), ebarMin := ebarMinFloor]
  
  pRank = ggplot(dt1, aes(x = meanRA, y = RankFac)) +
    scale_x_log10() +
    xlab("Mean Relative Abundance") +
    ylab(Rank) +
    theme_bw()
  if(!is.null(GroupBy)){
    # pRank <- pRank + facet_wrap(facets = as.formula(paste("~", GroupBy)))
    pRank <- pRank + geom_point(mapping = aes_string(colour = GroupBy),
                                size = 5)
  } else {
    # Don't include error bars for faceted version
    pRank <- pRank + geom_errorbarh(aes(xmax = ebarMax,
                                        xmin = ebarMin))
  }
  return(pRank)
}

plot_taxa_summary(phyobj, "Phylum")
plot_taxa_summary(phyobj, "Class")
summarize_taxa(phyobj, "Phylum")
summarize_taxa(phyobj, "Class")



#############################################
##                                         ##
##      PRE-PROCESSING                     ##
##                                         ##
#############################################

##########################
## TAXANOMIC FILTERING
##########################

# What phyla are represented? 
get_taxa_unique(phyobj, "Phylum")

# Create table, number of features for each phyla 
table.features_phyla <- as.data.frame(table(tax_table(phyobj)[,"Phylum"], exclude=NULL))

# Remove NA and ambiguous phylum annotation 
phyobj <- subset_taxa(phyobj, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
table.features_phyla2 <- as.data.frame(table(tax_table(phyobj)[,"Phylum"], exclude=NULL))

# Examine prevalence of each feature
prevdf <- apply(X=otu_table(phyobj), 
               MARGIN = ifelse(taxa_are_rows(phyobj), yes=1, no=2), 
               FUN=function(x) {sum(x>0)})

# Add taxonomy and total read counts
prevdf <- data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phyobj),
                    tax_table(phyobj))

# Compute the total and average prevalences of the features in each phylum
prevdfSUMS <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
prevdfSUMS
names(prevdfSUMS)[2] <- "meanprev"
names(prevdfSUMS)[3] <- "sumprev"

# Filter out low prevalence phyla, where prevalence <=1
PhyForFilt <- c("Chlamydiae", "Kiritimatiellaeota", "Deinococcus-Thermus", "Deferribacteres")
subpo <- subset_taxa(phyobj, !Phylum %in% PhyForFilt)


##########################
## SAMPLE PREV. FILTERING
##########################

# Subset to remaining phyla 
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(subpo, "Phylum"))

# Plot taxa prevalence versus total count
prevfilterplot <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(subpo),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none") + 
  ggtitle("Taxa prevalence versus total count (filtering step) ")

# Define prevalence threshold as 2% of total samples (.02*218=4.36)
prevthresh <- 0.02 * nsamples(subpo) 

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevthresh)]
subpo2 <- prune_taxa(keepTaxa, subpo)


##########################
## SAMPLE  FILTERING
##########################

# Variance stabilize with log transformation 
logt  = transform_sample_counts(subpo2, function(x) log(1 + x) )

# PCoA with bray to look for outliers 
out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues
plot_ordination(logt, out.pcoa.logt, type = "samples", title='PCoA (bray) to look for outliers') + 
  coord_fixed(sqrt(evals[2] / evals[1])) 

# MDS with bray to look for outliers 
out.mds.logt <- ordinate(logt, method = "MDS", distance = "bray")
plot_ordination(logt, out.mds.logt, type = "samples",  title='MDS (bray) to look for outliers') + 
  coord_fixed(sqrt(evals[2] / evals[1])) 

# PCoA with bray for phylum
plot_ordination(logt, out.pcoa.logt, type = "species", color = "Phylum", 
                title='PCoA (bray) by phylum') +
  coord_fixed(sqrt(evals[2] / evals[1]))

# No clear outliers --> no samples dropped


##########################
## CHECK READS
##########################

# Get total number of reads 
sum(colSums(otu_table(subpo2)))
sample_sum_df <- data.frame(sum=sample_sums(subpo2))

# Plot read count
ggplot(sample_sum_df, aes(x=sum))+ 
  geom_histogram(color="black", fill="indianred", binwidth = 2500) + 
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") + 
  theme(axis.title.y = element_blank()) 

# Get read stats
min(sample_sums(subpo2)) 
mean(sample_sums(subpo2)) 
max(sample_sums(subpo2))
median(sample_sums(subpo2)) 

# Look for low performing samples
qplot(colSums(otu_table(logt)),bins=10) + xlab("Logged counts-per-sample")

#############################################
##                                         ##
##      NORMALIZATION AND TRANSFORMATION   ##
##                                         ##
#############################################

# Is log transformation sufficient for normalizing the data? 
qplot(sample_data(subpo2)$age_2pm, geom = "histogram", binwidth=1) + xlab("age")
qplot(log10(rowSums(otu_table(subpo2))), binwidth=1) +  xlab("Logged counts-per-sample")

# Yes. 


#############################################
##                                         ##
##     ABUNDANCE                           ##
##                                         ##
#############################################

# Setting theme for plot
theme_set(theme_bw())
cols <- c("Firmicutes"= "dodgerblue2", 
          "Bacteroidetes"="#E31A1C", 
          "Actinobacteria"= "green4", 
          "Proteobacteria"="#6A3D9A",
          "Epsilonbacteraeota"="#FF7F00",
          "Euryarchaeota"="skyblue2", 
          "Verrucomicrobia"="#FB9A99", 
          "Synergistetes"="palegreen2", 
          "Fusobacteria"= "#CAB2D6", 
          "Cyanobacteria"="#FDBF6F", 
          "Lentisphaerae"="maroon", 
          "Spirochaetes"="orchid1", 
          "Elusimicrobia"="gold1", 
          "Tenericutes"="blue1", 
          "Acidobacteria"= "darkorange4", 
          "Patescibacteria" = "blueviolet")


##########################
## BARPLOTS
##########################

# Barplot of relative abundance (phylum) by pFGV tercile 
subpo2.phy.pfgv <- subpo2 %>% 
  tax_glom(taxrank="Phylum") %>%  #agglomerate to level of phylum
  transform_sample_counts(function(x) x / sum(x)) %>%  #convert to proportions 
  merge_samples("pFGV_tert") %>%  # sum abundances within categories
  transform_sample_counts(function(x) x / sum(x)) %>%  #convert to proportions within category
  psmelt() %>% # prep for ggplot 
  arrange(Phylum) # sort by phylum 
ra_phy_pfgv <- ggplot(subpo2.phy.pfgv, aes(x=pFGV_tert, y=Abundance, fill=Phylum)) + 
  scale_fill_manual(values=cols)+ 
  theme_classic() + 
  geom_bar(stat = "identity") +
  ylab("Relative Abundance \n") + 
  xlab("Percent FGV (Tercile)")

# Barplot of relative abundance (phylum) by aFGV tercile 
subpo2.phy.afgv <- subpo2 %>% 
  tax_glom(taxrank="Phylum") %>% 
  transform_sample_counts(function(x) x / sum(x)) %>% 
  merge_samples("aFGV_tert") %>%  
  transform_sample_counts(function(x) x / sum(x)) %>%  
  # filter(Abundance > 0.01) %>% 
  psmelt() %>% 
  arrange(Phylum)
ra_phy_afgv <- ggplot(subpo2.phy.afgv, aes(x=aFGV_tert, y=Abundance, fill=Phylum)) + 
  scale_fill_manual(values=cols)+ 
  theme_classic() + 
  geom_bar(stat = "identity") +
  ylab("Relative Abundance \n") + 
  xlab("Absolute FGV (Tercile)")

# Group together  [FIGURE 1]
ggarrange(ra_phy_pfgv, ra_phy_afgv,
          labels=c("A", "B"),
          ncol=2, nrow=1, common.legend = TRUE, legend = "right")

# Summarize by phylum and tercile 
raphyavg_pfgv <- subpo2.phy.pfgv %>% 
  group_by(Phylum, pFGV_tert) %>% 
  summarise(mean=mean(Abundance))

raphyavg_afgv <- subpo2.phy.afgv %>% 
  group_by(Phylum, aFGV_tert) %>% 
  summarise(mean=mean(Abundance))


##########################
## BOXPLOTS AND STATS
##########################

# Transformation 
subpo2ra <- transform_sample_counts(subpo2, function(x) x/sum(x))
subpo2_phyglom <- tax_glom(subpo2ra, taxrank = 'Phylum')
dat <- data.table(psmelt(subpo2_phyglom))
dat$Phylum <- as.character(dat$Phylum)
dat[, median := median(Abundance, na.rm = TRUE), 
    by = "Phylum"]
dat[(median <= 0.001), Phylum := "Other"]

# Order Phylum by abundance
dat$Phylum <- factor(dat$Phylum, levels=unique(as.character(dat$Phylum)))
dat <- transform(dat, Phylum=reorder(Phylum, Median)) 

# Set pairwise comparisons to make
comparisonsra <- list(c("T1", "T2"), c("T2","T3"), c("T1","T3"))

# Boxplot 
relabun_pfgv <- ggboxplot(dat[Abundance > 0], x="pFGV_tert", y="Abundance", 
                          fill="pFGV_tert", facet.by="Phylum", nrow=1) + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_y_log10() + 
  scale_fill_viridis(discrete=TRUE, name="Percent FGV (Tercile)") +
  ylab("Relative Abundance") + 
  theme(legend.position="bottom") +
  stat_compare_means(aes(group=pFGV_tert),method="kruskal.test", label="p.format",  label.y.npc="bottom", label.x.npc="left", hide.ns = TRUE) + 
  stat_compare_means( comparisons=comparisonsra, label="p.signif", tip.length = .01) 



relabun_afgv <- ggboxplot(dat[Abundance > 0], x="aFGV_tert", y="Abundance", 
                          fill="aFGV_tert", facet.by="Phylum", nrow=1) + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_y_log10() + 
  scale_fill_viridis(discrete=TRUE, name="Aboslute FGV (Tercile)") +
  ylab("Relative Abundance") + 
  theme(legend.position="bottom") +
  stat_compare_means(aes(group=aFGV_tert),method="kruskal.test", label="p.format", label.y.npc="bottom", label.x.npc="left", hide.ns = TRUE) + 
  stat_compare_means( comparisons=comparisonsra, label="p.signif", tip.length = .01) 

#[ADDITIONAL FIGURE 1]
relabun_pfgv 
#[ADDITIONAL FIGURE 2]
relabun_afgv


#############################################
##                                         ##
##      ALPHA DIVERSITY ANALYSES           ##
##                                         ##
#############################################

##########################
## RAREFACTION
##########################

# Plot of the number of species as a fn of the number of samples
rarecurve(t(otu_table(subpo2)), step=50, cex=0.5, label = FALSE)

# Rarefy without replacement to minimum 
subpo2rf <- rarefy_even_depth(subpo2, sample.size=min(sample_sums(subpo2)), rngseed=801, replace=F)

# Check sample sums before and after rarefying 
sample_sums(subpo2)
sample_sums(subpo2rf)


##########################
## PLOT ALPHA DIVERSITY
##########################

# Define alpha diversity metrics 
alpha_meas = c("Observed","Shannon")

# Plot alpha diversity for percent FGV (pFGV) and absolute FGV (aFGV)
adiv.rich.pfgv <- plot_richness(subpo2rf, "pFGV_tert", measures="Observed")
adiv.shan.pfgv <- plot_richness(subpo2rf, "pFGV_tert", measures="Shannon")
adiv.rich.afgv <- plot_richness(subpo2rf, "aFGV_tert", measures='Observed')
adiv.shan.afgv <- plot_richness(subpo2rf, "aFGV_tert", measures='Shannon')

# Set pairwise comparisons
comparisons <- list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3"))

# Observed Richness, percent FGV
adiv.rich.pfgv.box <- adiv.rich.pfgv +
  geom_boxplot(data=adiv.rich.pfgv$data,
               aes(pFGV_tert, value, fill=pFGV_tert), 
               show.legend=FALSE,
               size=.5, alpha=.5,
               outlier.shape=NA)+
  theme_bw() + 
  scale_fill_viridis(discrete=TRUE)+ 
  xlab("Percent FGV (Tercile)")+ 
  ylab("Observed Richness")+
  stat_compare_means(method="kruskal.test", label.y=2.3 , label.x.npc="left", label="p.format") + 
  stat_compare_means(comparisons=comparisons, label="p.format", hide.ns = FALSE, method = "wilcox.test", tip.length = .01)+ 
  theme(strip.background = element_blank(), strip.text = element_blank()) + 
  coord_cartesian(ylim=c(0,420))

# Shannon index, percent FGV
adiv.shan.pfgv.box <- adiv.shan.pfgv +
  geom_boxplot(data=adiv.shan.pfgv$data,
               aes(pFGV_tert, value, fill=pFGV_tert), 
               show.legend=FALSE,
               size=.5, alpha=.5,
               outlier.shape=NA)+
  theme_bw() + 
  scale_fill_viridis(discrete=TRUE)+ 
  xlab("Percent FGV (Tercile)")+ 
  ylab("Shannon Index")+
  stat_compare_means(method="kruskal.test", label.y=2.3 , label.x.npc="left", label="p.format")+
  theme(strip.background = element_blank(), strip.text = element_blank())

# Observed Richness, absolute FGV
adiv.rich.afgv.box <- adiv.rich.afgv +
  geom_boxplot(data=adiv.rich.afgv$data,
               aes(aFGV_tert, value, fill=aFGV_tert), 
               show.legend=FALSE,
               size=.5, alpha=.5,
               outlier.shape=NA)+
  theme_bw() + 
  scale_fill_viridis(discrete=TRUE)+ 
  xlab("Absolute FGV (Tercile)")+ 
  ylab("Observed Richness")+
  stat_compare_means(method="kruskal.test", label.y=2.3 , label.x.npc="left", label="p.format") +
  theme(strip.background = element_blank(), strip.text = element_blank())+ 
  coord_cartesian(ylim=c(0,420))

# Shannon index, absolute FGV
adiv.shan.afgv.box <- adiv.shan.afgv +
  geom_boxplot(data=adiv.shan.afgv$data,
               aes(aFGV_tert, value, fill=aFGV_tert), 
               show.legend=FALSE,
               size=.5, alpha=.5,
               outlier.shape=NA)+
  theme_bw() + 
  scale_fill_viridis(discrete=TRUE)+ 
  xlab("Absolute FGV (Tercile)")+ 
  ylab("Shannon Index")+
  stat_compare_means(method="kruskal.test", label.y=2.3 , label.x.npc="left", label="p.format")+
  theme(strip.background = element_blank(), strip.text = element_blank())

# Combine to single image  [FIGURE 2]
ggarrange(adiv.rich.pfgv.box, adiv.shan.pfgv.box, adiv.rich.afgv.box, adiv.shan.afgv.box,
          labels=c("A", "B", "C", "D"),
          ncol=2, nrow=2, common.legend = TRUE, legend = "bottom")


##########################
## ALPHA DIVERSITY STATS
##########################

# Get datasets
adiv.rich.df <- as.data.frame(adiv.rich.pfgv.box$data)
adiv.shan.df <- as.data.frame(adiv.shan.pfgv.box$data)

# Get mean by group (replace df as necessary)
adivsumall <- adiv.rich.df %>% 
  group_by(variable) %>% 
  summarise(mean=mean(value))

adivsumpfgv <- adiv.rich.df %>% 
  group_by(variable, pFGV_tert) %>% 
  summarise(mean=mean(value))

adivsumafgv <- adiv.rich.df %>% 
  group_by(variable, aFGV_tert) %>% 
  summarise(mean=mean(value))

# Get median by group 
adivsumpfgv <- adivdf %>% 
  group_by(variable) %>% 
  summarise(median=median(value))

adivsumafgv <- adivdf %>% 
  group_by(variable) %>% 
  summarise(median=median(value))


#############################################
##                                         ##
##      BETA DIVERSITY ANALYSES            ##
##                                         ##
#############################################

##########################
## PCoA BRAY CURTIS 
##########################

# Ordinate, PCoA using Bray-Curtis as metric 
pcoa_bray_ord <- ordinate(subpo2rf, method="PCoA", distance="bray")

# TOtal distance captured 
plot_scree(pcoa_bray_ord,'Scree plot, Bray/PCoA')

# How much variation do the first two axes explain? 15.02%
(100*sum(pcoa_bray_ord$values$Relative_eig[1:2])) 

# PCoA plots
pcoa_bray_pfgv <- plot_ordination(subpo2rf, pcoa_bray_ord, color="pFGV_tert") + theme(aspect.ratio=1) + 
  theme_bw()+
  scale_color_viridis(discrete=TRUE)+ 
  labs(color="Percent FGV (Tercile)")+
  stat_ellipse(size=1)+
  theme(legend.position="bottom")

pcoa_bray_afgv <- plot_ordination(subpo2rf, pcoa_bray_ord, color="aFGV_tert") + theme(aspect.ratio=1) + 
  theme_bw()+
  scale_color_viridis(discrete=TRUE)+ 
  labs(color="Absolute FGV (Tercile)")+
  stat_ellipse(size=1)+
  theme(legend.position="bottom")

# Group together [FIGURE 3]
ggarrange(pcoa_bray_pfgv, pcoa_bray_afgv,
          labels=c("A", "B"),
          ncol=2, nrow=1, common.legend = FALSE, legend = "bottom")

##########################
## PERMANOVA
##########################

# Calculate distances 
bray_dist <- phyloseq::distance(subpo2rf, method="bray")

# Model 1: percent FGV
permbray1 <- adonis2(bray_dist ~ sample_data(subpo2rf)$pFGV_tert)

# Model 2: percent FGV ~ Fat
permbray2 <- adonis2(bray_dist ~ sample_data(subpo2rf)$pFGV_tert 
                     + sample_data(subpo2rf)$fatcat_2pm)

# Model 3: percent FGV ~ fat, age, antibiotic use, birth mode, breast feeding, caloric intake, medu, ethnicity  [TABLE 2]
permbray3 <- adonis2(bray_dist ~ sample_data(subpo2rf)$pFGV_tert 
                     + sample_data(subpo2rf)$fatcat_2pm
                     + sample_data(subpo2rf)$age_2pm
                     + sample_data(subpo2rf)$abx_6mo
                     + sample_data(subpo2rf)$birthmode
                     + sample_data(subpo2rf)$bfeedcat
                     + sample_data(subpo2rf)$AvgCalQ
                     + sample_data(subpo2rf)$medu 
                     + sample_data(subpo2rf)$ethnicity
                     + sample_data(subpo2rf)$tvcat)

# Model 4: absolute FGV
permbray4 <- adonis2(bray_dist ~ sample_data(subpo2rf)$aFGV_tert)

# Model 5: absolute FGV ~ Fat
permbray5 <- adonis2(bray_dist ~ sample_data(subpo2rf)$aFGV_tert 
                     + sample_data(subpo2rf)$fatcat_2pm)

# Model 6: absolute FGV ~ fat, age, antibiotic use, birth mode, breast feeding, caloric intake, medu, ethnicity [TABLE 2]
permbray6 <- adonis2(bray_dist ~ sample_data(subpo2rf)$aFGV_tert 
                     + sample_data(subpo2rf)$fatcat_2pm
                     + sample_data(subpo2rf)$age_2pm
                     + sample_data(subpo2rf)$abx_6mo
                     + sample_data(subpo2rf)$birthmode
                     + sample_data(subpo2rf)$bfeedcat
                     + sample_data(subpo2rf)$AvgCalQ
                     + sample_data(subpo2rf)$medu 
                     + sample_data(subpo2rf)$ethnicity
                     + sample_data(subpo2rf)$tvcat)

# Evaluate homogeneity of dispersion 
beta_p <- betadisper(bray_dist, sample_data(subpo2rf)$pFGV_tert)
permutest(beta_p)
beta_a <- betadisper(bray_dist, sample_data(subpo2rf)$aFGV_tert)
permutest(beta_a)
# Not significant, meaning we cannot reject the null hypothesis that our groups have the same dispersion



#############################################
##                                         ##
##   SAVE PHYLOSEQ OBJECTS                 ##
##                                         ##
#############################################

# OTU Table
subpo2OTU <- as(otu_table(subpo2), "matrix")
if(taxa_are_rows(subpo2)){subpo2OTU <- t(subpo2OTU)}
subpo2OTUdf <- as.data.frame(subpo2OTU)
subpo2OTUdf <- rownames_to_column(subpo2OTUdf, "SampleID")
write.table(subpo2OTUdf, file="Data/FinalData/OTUtable_filt.csv", sep="\t", row.names=FALSE)

# Tax Table 
subpo2TAX <- as(tax_table(subpo2), "matrix")
if(taxa_are_rows(subpo2)){subpo2TAX <- t(subpo2TAX)}
subpo2TAXdf <- as.data.frame(subpo2TAX)
subpo2TAXdf <- rownames_to_column(subpo2TAXdf, "Taxonomy")
write.table(subpo2TAXdf, file="Data/FinalData/TAXtable_filt.csv", sep="\t", row.names=FALSE)

# Sample data 
view(sample_data(subpo2))
subpo2samp <- subpo2 %>%  
  sample_data() %>%  
  data.frame() %>%  
  rownames_to_column() %>%  
  as_tibble()
names(subpo2samp)
subpo2samp$SampleID <- subpo2samp$rowname
subpo2sampdf <- subpo2samp %>%
  select(SampleID, everything()) %>%  
  select(-c("rowname"))
write.csv(subpo2sampdf, "Data/FinalData/Meta_after_filter_2.csv",  row.names = FALSE )
