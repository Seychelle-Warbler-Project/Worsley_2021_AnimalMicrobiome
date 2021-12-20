#Sarah Worsley

# R code and analysis for: Gut microbiome composition, not alpha diversity, is associated with survival in a natural vertebrate population

#load packages
  library(ggplot2)
  library(phyloseq)
  library(ape)
  library(philr)
  library(microbiome)
  library(pairwiseAdonis)
  library(vegan)
  library(ggordiplots)
  library(BiodiversityR)

#Upload ASV, taxonomy, tree and metadata files to phyloseq:

#ASV table
  asv_table <- read.csv ("ASV_table.csv", row.names=1) #read in asv table with feature names as rownames
  str (asv_table) #should be 56690 features, 636 samples (note 2 of these samples were filtered out during QIIME steps as low read number)
  head (asv_table)
  asv_table <- as.matrix (asv_table) #make into a matrix

#taxonomy file
  taxonomy <- read.csv ("taxonomy.csv", row.names=1)
  str (taxonomy) #should be 56690 observations, 7 taxonomic groupings, feature names as row names.
  taxonomy <- as.matrix (taxonomy)

#read in tree as a phyloseq object
  phy_tree <- read_tree ("tree.nwk")

#load metadata file
  metadatafull <- read.csv("Metadata_Full.csv")
  str(metadatafull)  # contains 638 samples- need to eliminate the 2 samples filtered by QIIME
  row.names(metadatafull)<-metadatafull$Sample.ID
  
#remove filtered samples
  samples_remove <- c("Sample.102.115", "Sample.221.427")
  metadata_filter <- metadatafull[!metadatafull$Sample.ID %in% samples_remove,]  
  str(metadata_filter) # 636 samples
  write.csv(metadata_filter, "metadata_filtered.csv", row.names=TRUE)
  
#load filtered metadata with sample names as rownames
  metadata<-read.csv("metadata_filtered.csv", row.names = 1)
  str(metadata)  
  
  head(metadata)

#import all as phyloseq objects
  ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
  TAX <- tax_table(taxonomy)
  META <- sample_data(metadata)
  head(META)

#check that the ASV and sample names are consistent across objects
  str(taxa_names(TAX))
  str(taxa_names(ASV))
  str(taxa_names(phy_tree))
  
  str(sample_names(ASV))
  str(sample_names(META))

#### MERGE INTO PHYLOSEQ OBJECT ####
  physeq <- phyloseq(ASV, TAX, META, phy_tree)
  physeq #56690 taxa, 636 samples 

#check the number of reads per sample
  sample_reads<-data.frame(reads=sample_sums(physeq))
  head(sample_reads)

##################
### FILTERING ####
##################
  
####filter to remove instances where features are not assigned as bacteria/ unassigned/ not assigned at phylum level/ #####

#summarise number of features at kingdom level- need to remove cases of archaea and taxa unassigned.  
  table(tax_table(physeq)[,"Kingdom"]) 
  # removes all labelled as unassigned at kingdom level
  physeq<-subset_taxa(physeq, !is.na(Kingdom) & !Kingdom %in% c("", "Unassigned")) 
  table(tax_table(physeq)[,"Kingdom"])

#remove those features labelled as Archaea  
  physeq<-subset_taxa(physeq, !Kingdom %in% "Archaea") 
  table(tax_table(physeq)[,"Kingdom"]) 

#print features per phylum- can see that 506 are phylum unassigned (have a blank name)
  table(tax_table(physeq)[,"Phylum"])
# removes all labelled as unassigned at phylum level
  physeq<-subset_taxa(physeq, !is.na(Phylum) & !Phylum %in% c("", "Unassigned"))
  table(tax_table(physeq)[,"Phylum"]) #check removed- gives the number of features per phylum

#check adjusted numbers - now 55664 taxa in 636 samples  
  physeq


### removal of samples from controls ####
# remove control samples (positive/negative controls, collection controls)  
  physeq2 <- subset_samples(physeq, SampleType == "F") 
  physeq2 #610 samples

### remove samples with less than 10000 reads 
  physeq2
  physeq3<-prune_samples(sample_sums(physeq2)>=10000, physeq2)
  physeq3 # removes 23 samples (some of these are sequencing duplicates)- leaves 587 samples



#### REMOVE SEQUENCING DUPLICATES ####
  
  #See SequencingExtractionRepeatability.R for an assessment of reproducibility
  
  #### extract samples sequenced/extracted twice (keep the one with the highest read number)###
  
  physeqNoDup<- subset_samples(physeq3, SequencingDuplicate=="No") #528 remaining
  physeqNoDup
  
  #### Remove catch duplicates ####
  
  # Keeping 1 sample per catch in following order of priority: 1) from tray 2)  from bag, or if neither keep highest read count
  
  physeqFiltered<- subset_samples(physeqNoDup, CatchRepeat=="No") #471 samples, 55664 taxa 
  physeqFiltered

###SAMPLE FILTERING###

#remove outlier sample (Sample.71.81)- very tiny sample with extremely low alpha diversity

  physeqFiltered<- subset_samples( physeqFiltered, Sample.ID != "Sample.71.81")
  physeqFiltered #470 samples, 55664 taxa  

#remove nestlings
  physeqFiltered<- subset_samples( physeqFiltered,Ageclass != "N")
  physeqFiltered #458 samples

#remove floaters with no TQ
  physeqFiltered<- subset_samples( physeqFiltered,Status!= "FLOAT")#leaves 450 samples
  physeqFiltered #450 samples, from 309 individuals
  datFiltered<- sample_data(physeqFiltered)
  unique(datFiltered$BirdID)
  

##########################
### PRUNE RARE TAXA ######
##########################
  
# Compute prevalence of each feature (total number of samples in which a taxon appears at least once), store as data.frame
  prevdf = apply(X = otu_table(physeqFiltered),
                 MARGIN = ifelse(taxa_are_rows(physeqFiltered), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts for each phylum to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(physeqFiltered),
                      tax_table(physeqFiltered))
  head(prevdf)
  str(prevdf)
  
  plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) #average prevelance of features within each phylum and the sum of feature prevalence within each phylum
  
# Plot the unique phyla: each dot will be a feature- total abundance of that feature across samples on x axis and the prevalance (the fraction of all samples it occurs in on the y axis).
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeqFiltered, "Phylum"))
  ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(physeqFiltered),color=Phylum)) +
  # Include filtering threshold here; 2%, total abundance across all samples set to 50
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_vline(xintercept = 50, alpha = 0.5, linetype = 2)+  geom_point(size = 1, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# Define prevalence threshold as 2% samples and abundance as 50 total reads across samples

  prevalenceThreshold<-2*(450/100)

  abundanceThreshold<-50

# Execute the prevalence filter, using `prune_taxa()` function

  head(prevdf1)

  KeepTaxa1<-prevdf1[prevdf1$Prevalence>=prevalenceThreshold,]
  str(KeepTaxa1)
  head(KeepTaxa1)

  KeepTaxa2<- rownames(KeepTaxa1)[(KeepTaxa1$TotalAbundance >= abundanceThreshold)]
  str(KeepTaxa2) #3247 taxa retained after filtering at 2% sample level

  physeqBeta<- prune_taxa(KeepTaxa2,physeqFiltered)
  physeqBeta # 3247 taxa (down from 56029)


  
######################################
#### CLR TRANSFORM ASV ABUNDANCES#####
######################################

  library(microbiome)
  physeq_clr <- microbiome::transform(physeqBeta, "clr")
  head(otu_table(physeq_clr))


#function to extract an ASV matrix
  vegan_otu <- function(physeq) {
    OTU <- otu_table(physeq)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }

#Extract ASV Matrix and Sample Data  
  clr_Matrix<-vegan_otu(physeq_clr)
  clr_SampleData<-as(sample_data(physeq_clr),"data.frame")


#Principal Components Analysis 
  library(vegan)
  clr_pca<-rda(clr_Matrix)
  
  library(BiodiversityR)
  sigCLR <- PCAsignificance (clr_pca, axes = 8) #how much does each principal component explain
  sigCLR #PC1= 6.91%, PC2= 4.14%

###############  
###PERMANOVA###
###############
  
#test for influence of age, sex, sampling season and territory quality on GM composition
  
  #correct data formats
  clr_SampleData$BirdID<- as.factor(clr_SampleData$BirdID)
  clr_SampleData$SexEstimate<- as.factor(clr_SampleData$SexEstimate)
  clr_SampleData$Ageclass<- as.factor(clr_SampleData$Ageclass)
  clr_SampleData$FieldPeriodIDFactor<- as.factor(clr_SampleData$FieldPeriodIDFactor)
  clr_SampleData$Ageclass <- factor(clr_SampleData$Ageclass, levels= c("FL", "OFL", "SA", "A"))
  clr_SampleData$FieldPeriodIDFactor <- factor(clr_SampleData$FieldPeriodIDFactor, levels= c("Major17", "Minor18", "Major18", "Minor19", "Major19", "Minor20"))
  str(clr_SampleData)
  
  #set strata blocks- stratify by BirdID to account for repeat sampling of individuals.
  perm <- how(nperm = 9999)
  set.seed(549888)
  setBlocks(perm) <- with(clr_SampleData, BirdID)
  permanovaAllData<- adonis2(clr_Matrix ~ Ageclass + SexEstimate + FieldPeriodIDFactor + TerritoryQuality, data=clr_SampleData, permutations = perm, method = "euclidean", by= "margin")
  permanovaAllData
  
  str(permanovaAllData)
  resPermAll<- data.frame(permanovaAllData[c(1,2,3,4,5)])
  resPermAll
  write.csv(resPermAll, "AllSamplesPermanova.csv") # CLR results for Table 1

############################  
### Pairwise permanovas ####
############################ 

### Pairwise permanova for age class - results reported in Table S1
  
  library(pairwiseAdonis)
  
  perm <- how(nperm = 9999)
  set.seed(549838)
  setBlocks(perm) <- with(clr_SampleData, BirdID)
  AgePairwise<- pairwise.adonis2(clr_Matrix ~ Ageclass + SexEstimate + FieldPeriodIDFactor + TerritoryQuality,
                   data = clr_SampleData, method="euclidean", by="margin", nperm=perm) 
  AgePairwise
  p.adjustedAge <- p.adjust(c(0.083,0.01,0.105,0.045,0.110, 0.001),method="BH")
  p.adjustedAge
  

## Pairwise permanova for field period - results reported in Table S1

  perm <- how(nperm = 9999)
  set.seed(54844)
  setBlocks(perm) <- with(clr_SampleData, BirdID)
  PairwiseFP<- pairwise.adonis2(clr_Matrix ~ FieldPeriodIDFactor + Ageclass + SexEstimate + TerritoryQuality,
                   data = clr_SampleData, method="euclidean", by="margin", nperm=perm) 
  PairwiseFP
  
  p.adjustedFP <- p.adjust(c(0.001,0.003,0.001,0.001,0.001, 0.007,0.001,0.001, 0.001, 0.001, 0.001, 0.001,
                             0.001,0.001, 0.020), method="BH")
 
  
### Checking group dispersions with betadisper analysis for age and field period
  #i.e are sig p values due to differences in dispersion?
  
  #age - reported in Fig S3
  distMatrixBeta <- vegdist(clr_Matrix, method="euclidean")
  BDAge<-betadisper(distMatrixBeta, clr_SampleData$Ageclass)
  perm <- how(nperm = 9999)
  set.seed(54)
  permutest(BDAge, permutations = perm)
  par(mar=c(6,6,6,6))
  plotDistancesAge<- boxplot(BDAge,ylab="Distance to centroid\n", xlab="\nAge class", ylim=c(20,130), cex.axis=1.5, cex.lab=1.5)
  plotDistancesAge
  set.seed(567)
  TukeyHSD(BDAge, which = "group", ordered = FALSE)
  
  #sex
  distMatrixBeta <- vegdist(clr_Matrix, method="euclidean")
  BDSex<-betadisper(distMatrixBeta, clr_SampleData$SexEstimate)
  perm <- how(nperm = 9999)
  set.seed(54)
  permutest(BDSex, permutations = perm)
  par(mar=c(6,6,6,6))
  plotDistancesSex<- boxplot(BDSex,ylab="Distance to centroid\n", xlab="\nSex Estimate", ylim=c(20,140), cex.axis=1.5, cex.lab=1.5)
  plotDistancesSex

  
  #Field period - reported in Fig S5
  distMatrixBeta <- vegdist(clr_Matrix, method="euclidean")
  BDFP<-betadisper(distMatrixBeta, clr_SampleData$FieldPeriodIDFactor)
  perm <- how(nperm = 9999)
  set.seed(55)
  permutest(BDFP, permutations = perm)
  par(mar=c(6,6,6,6))
  boxplot(BDFP,ylab="Distance to centroid\n", xlab="\nSampling period", ylim=c(20,140), cex.axis=1.5, cex.lab=1.5)
  TukeyHSD(BDFP, which = "group", ordered = FALSE)


###################  
#### PCA plots ####
###################
  
  library(ggordiplots)
    
#by sex - Fig S4
  plotPCASex<-gg_ordiplot(clr_pca,groups=clr_SampleData$SexEstimate,plot=FALSE, conf=0.95, pt.size = 2) 
  SexPlot<-plotPCASex$plot
  SexPlot + theme_bw()  + guides(color=guide_legend("Sex")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 6.91%") + ylab("PC2 4.14%\n") + theme(axis.text = element_text(size=20), axis.title = element_text(size=20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
    scale_color_manual(values= c("#E69F00","#0072B2"))
  
  table(clr_SampleData$Ageclass)
  
#by age class - Fig S3
  plotPCAAge<-gg_ordiplot(clr_pca,groups=clr_SampleData$Ageclass,plot=FALSE, conf = 0.95, pt.size = 2) 
  AgeDat<- plotPCAAge$df_ord
  AgeDat$Ageclass<- as.factor(clr_SampleData$Ageclass)
  str(AgeDat)
  AgePlot<- ggplot(data= AgeDat, aes(x=x, y=y, fill=Ageclass, shape=Ageclass, color=Ageclass)) + theme_bw() +
    geom_point(size=4)+
    guides(color=guide_legend("Ageclass")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 6.91%") + ylab("PC2 4.14%\n") + scale_shape_manual(values=c(21,22,23,24)) + 
    theme(axis.text = element_text(size=20), axis.title = element_text(size=20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
    scale_fill_manual(values= c("purple4","goldenrod","aquamarine4", "grey65")) +
    scale_color_manual(values= c("purple4","goldenrod","aquamarine4", "grey65"))
  
  AgePlot
  
  table(clr_SampleData$Ageclass)
  
#plot by season(major/minor) and year (2017,2018,2019, 2020) - Figure 2
  
  plotPCASeason<-gg_ordiplot(clr_pca,groups=clr_SampleData$FieldPeriodIDFactor,plot=FALSE, conf = 0.95, pt.size = 2)  
  ord.data <- plotPCASeason$df_ord
  ord.data$SampleYear<-as.factor(clr_SampleData$SampleYear)
  ord.data$Season<-as.factor(clr_SampleData$FieldPeriodIDBinary)
  head(ord.data)
  
  FPAndYear<- ggplot(data = ord.data, aes(x = x, y = y, color = SampleYear, fill=SampleYear, shape = Season)) +
    theme_bw() + geom_point(size = 3, stroke=1.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 6.91%") + ylab("PC2 4.14%\n") + scale_shape_manual(values=c(2,1)) + 
    scale_colour_manual(values=c("goldenrod","aquamarine4", "purple4", "grey65")) +
    scale_fill_manual(values=c("goldenrod","aquamarine4", "purple4", "grey65")) + 
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22), legend.title = element_text(size=18), legend.text = element_text(size=18))
  
  FPAndYear
  

########################################
### PhILR transformed asv abundances####
########################################

  library(philr)

# Add a pseudocount of 1 to ASVs in the filtered phyloseq object- avoids calculating log-ratios involving zeros

  physeqBetaPseudo <- transform_sample_counts(physeqBeta, function(x) x+1)
  otu_table(physeqBetaPseudo)

# Next check that the phylogenetic tree is rooted and binary (all multichotomies have been resolved). 

  is.rooted(phy_tree(physeqBetaPseudo)) # Is the tree Rooted? Yes
  is.binary(phy_tree(physeqBetaPseudo)) # All multichotomies resolved? - false
  physeqBetaTree <- multi2di(phy_tree(physeqBetaPseudo)) # use multi2di in ape to replace multichotomies with a series of dichotomies with one (or several) branch(es) of zero length
  
#number the internal nodes of the tree so it's easier to work with- prefix the node number with "n"
  physeqBetaTreeN <- makeNodeLabel(physeqBetaTree, method="number", prefix='n')
  name.balance(physeqBetaTreeN, tax_table(physeqBetaPseudo), 'n86') #find a consensus name for the two clades that descend from a given balance (out of x/y, x= numerator, y= denominator in balance)

#transpose data: ensures taxa are columns, samples are rows
  otu.table <- t(otu_table(physeqBetaPseudo))
  tree <- physeqBetaTreeN
  metadata <- sample_data(physeqBetaPseudo)
  tax <- tax_table(physeqBetaPseudo)
  
  otu.table[1:2,1:2] # OTU Table
  tree # Phylogenetic Tree
  head(metadata,2) # Metadata
  head(tax,2) # taxonomy table

#transform the data using philr (wrapper function)
#uses taxon weighting:down weight the influence of taxa with many zeros or non-zero counts
#branch length weighting: scale balances using phylogenetic distance (communities differing in abundance of closely related species are more similar than those differing in distantly related species)
  phy.philr <- philr(otu.table, tree, 
                   part.weights='enorm.x.gm.counts', 
                   ilr.weights='blw.sqrt')
  phy.philr[1:5,1:5] #expressed as balances: these represent the log-ratio of the geometric mean abundance of the two groups of taxa that descend from a given internal node


######################  
###PERMANOVA PhILR ###
######################
  
  #test for influence of age, sex, sampling season and territory quality on GM composition
  
  #correct data formats
  philr_SampleData<- data.frame(metadata)
  philr_SampleData$BirdID<- as.factor(philr_SampleData$BirdID)
  philr_SampleData$SexEstimate<- as.factor(philr_SampleData$SexEstimate)
  philr_SampleData$Ageclass<- as.factor(philr_SampleData$Ageclass)
  philr_SampleData$FieldPeriodIDFactor<- as.factor(philr_SampleData$FieldPeriodIDFactor)
  philr_SampleData$Ageclass <- factor(philr_SampleData$Ageclass, levels= c("FL", "OFL", "SA", "A"))
  philr_SampleData$FieldPeriodIDFactor <- factor(philr_SampleData$FieldPeriodIDFactor, levels= c("Major17", "Minor18", "Major18", "Minor19", "Major19", "Minor20"))
  str(philr_SampleData)
  
  #set strata blocks- stratify by BirdID to account for repeat sampling of individuals.
  perm <- how(nperm = 9999)
  set.seed(54998)
  setBlocks(perm) <- with(philr_SampleData, BirdID)
  permanovaAllDataPhilr<- adonis2(phy.philr ~ Ageclass + SexEstimate + FieldPeriodIDFactor +
                                    TerritoryQuality, data=philr_SampleData, permutations = perm,
                                  method = "euclidean", by= "margin")
  permanovaAllDataPhilr
  
  resPermAllphilr<- data.frame(permanovaAllDataPhilr[c(1,2,3,4,5)])
  resPermAllphilr
  write.csv(resPermAllphilr, "AllSamplesPermanovaPhilr.csv") #reported in Table 1

  
############################  
### Pairwise permanovas ####
############################ 
  
  ### Pairwise permanova for field period - reported in Table S1
  
  library(pairwiseAdonis)

  perm <- how(nperm = 9999)
  set.seed(634)
  setBlocks(perm) <- with(philr_SampleData, BirdID)
  PairwiseFPPhilr<- pairwise.adonis2(phy.philr ~ FieldPeriodIDFactor + Ageclass + SexEstimate + TerritoryQuality,
                                data = philr_SampleData, method="euclidean", by="margin", nperm = perm) 
  PairwiseFPPhilr
  
  p.adjustedFPphilr <- p.adjust(c(0.034,0.163,0.015,0.001,0.002, 0.001,0.001,0.001, 0.001, 0.180, 0.006, 0.077,
                             0.001,0.001, 0.040), method="BH")
  
  
### Checking group dispersions with betadisper analysis for sex and field period###
  
  #Sex
  distMatrixBetaPD <- vegdist(phy.philr, method="euclidean")
  BDSexPD<-betadisper(distMatrixBetaPD, philr_SampleData$SexEstimate)
  set.seed(2699)
  permutest(BDSexPD, permutations = 9999) 
  boxplot(BDSexPD,ylab="Distance to centroid\n", xlab="\nSex", cex.axis=1.5, cex.lab=1.5)

  str(philr_SampleData)
  
  #Field period
  distMatrixBeta <- vegdist(phy.philr, method="euclidean")
  BDFPPD<-betadisper(distMatrixBeta, philr_SampleData$FieldPeriodIDFactor)
  set.seed(269)
  permutest(BDFPPD, permutations = 9999)
  boxplot(BDFPPD,ylab="Distance to centroid\n", xlab="\nSampling period", cex.axis=1.5, cex.lab=1.5)
  TukeyHSD(BDFPPD, which = "group", ordered = FALSE)

##################
#### PCA plots####
##################  
  
  ## Ordination in PhILR Space
  
  library(vegan)
  ilr_pca<-rda(phy.philr)
  
  library(ggordiplots)
  
  sigPhILR <- PCAsignificance (ilr_pca, axes = 8) #how much does each principal component explain
  sigPhILR #PC1= 20.58, PC2= 10.44

  library(ggordiplots)
  
### by sex - Fig S4
  plotPCASex<-gg_ordiplot(ilr_pca,groups=philr_SampleData$SexEstimate,plot=FALSE, conf=0.95, pt.size = 2) 
  SexPlot<-plotPCASex$plot
  SexPlot + theme_bw()  + guides(color=guide_legend("Sex")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 20.58%") + ylab("PC2 10.44%\n") + theme(axis.text = element_text(size=20), axis.title = element_text(size=20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
    scale_color_manual(values= c("#E69F00","#0072B2"))
  
### by field period - Figure 2
  #plot with minor/major season and year
  plotPCAphilFP<-gg_ordiplot(ilr_pca,groups=metadata$FieldPeriodIDFactor,plot=FALSE, conf = 0.95, pt.size=2) 
  ord.dataphilr <- plotPCAphilFP$df_ord
  ord.dataphilr$SampleYear<-as.factor(philr_SampleData$SampleYear)
  ord.dataphilr$Season<-as.factor(philr_SampleData$FieldPeriodIDBinary)
  head(ord.dataphilr)
  
  FPplotPhilr<- ggplot(data = ord.dataphilr, aes(x = x, y = y, color = SampleYear, fill=SampleYear, shape = Season)) +
    theme_bw() + geom_point(size = 3, stroke=1.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 20.58%") + ylab("PC2 10.44%\n") + scale_shape_manual(values=c(2,1)) + 
    scale_colour_manual(values=c("goldenrod","aquamarine4", "purple4", "grey65")) +
    scale_fill_manual(values=c("goldenrod","aquamarine4", "purple4", "grey65")) + 
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22), legend.title = element_text(size=18), legend.text = element_text(size=18))
  FPplotPhilr
  

  

  

