#Sarah Worsley

#R code and analysis for: Gut microbiome composition, not alpha diversity, is associated with survival in a natural vertebrate population

#load packages
  library(ggplot2)
  library(phyloseq)
  library(ape)
  library(microbiome)
  library(vegan)
  library(ggordiplots)
  library(data.table)


#Upload files to phyloseq:

#ASV table
  asv_table <- read.csv ("ASV_table.csv", row.names=1) #read in asv table with feature names as rownames
  str (asv_table) #should be 56690 features, 636 samples
  head (asv_table)
  asv_table <- as.matrix (asv_table) #make into a matrix

#taxonomy
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

#import all as phyloseq objects
  ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
  TAX <- tax_table(taxonomy)
  META <- sample_data(metadata)
  head(META)

#check that the ASV and sample names are consistent across objects (e.g have dashes been replaced with dots?)
  str(taxa_names(TAX))
  str(taxa_names(ASV))
  str(taxa_names(phy_tree))

  str(sample_names(ASV))
  str(sample_names(META))

#### MERGE INTO PHYLOSEQ OBJECT ####
  physeq <- phyloseq(ASV, TAX, META, phy_tree)
  physeq #check there are the right numbers of samples/features/variables etc...  

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

#check adjusted numbers - now 55664 taxa and 636 samples 
  physeq


#### removal of samples from controls####
# remove control samples (positive/negative controls, collection controls)  
  physeq2 <- subset_samples(physeq, SampleType == "F") 
  physeq2 #610 samples

### remove samples with less than 10000 reads (didn't reach completion: see SequencingExtractionRepeatability.R for iNext code) 
### BUT DO NOT rarefy for beta diversity analysis
  physeq2
  physeq3<-prune_samples(sample_sums(physeq2)>=10000, physeq2)
  physeq3 # removes 23 samples (some of these are sequencing duplicates)- leaves 587 samples


######################################
#### REMOVE SEQUENCING DUPLICATES ####
######################################
  
#See SequencingExtractionRepeatability.R for assessment of repeatability

#### extract samples sequenced/extracted twice (keep the one with the highest read number)###
  
  physeqNoDup<- subset_samples(physeq3, SequencingDuplicate=="No") #528 remaining
  physeqNoDup
  
#### Remove catch duplicates ####
  
  # Keeping 1 sample per catch in following order of priority: 1) from tray 2)  from bag, or if neither keep highest read count
  
  physeqFiltered<- subset_samples(physeqNoDup, CatchRepeat=="No") #471 samples, 55664 taxa 
  physeqFiltered

##########################
#### SAMPLE FILTERING ####
##########################
  
  
#remove outlier sample (Sample.71.81)- very tiny sample with extremely low alpha diversity
  
  physeqFiltered<- subset_samples( physeqFiltered, Sample.ID != "Sample.71.81")
  physeqFiltered #470 samples, 55664 taxa  
  
  range(sample_sums(physeqFiltered)) #range of read numbers across samples
  observationThreshold <- 1
  
  #check ASV number per sample in unrarefied samples
  ASVno470SamplesUnrarefied<-data.frame(estimate_richness(physeqFiltered, split=TRUE, measures= c("Observed")))
  head(ASVno470SamplesUnrarefied)
  str(ASVno470SamplesUnrarefied) #470 samples
  summary(ASVno470SamplesUnrarefied) #mean ASV no= 399.2
  sd(ASVno470SamplesUnrarefied$Observed) #sd 285.24
  
  #how many birds?
  dataFiltered<- sample_data(physeqFiltered)
  str(dataFiltered)
  str(unique(dataFiltered$BirdID)) #how many unique birds are there = 320 individuals


### For survival analysis- remove 2020 samples as no follow-up census to assess survival

  filteredMeta<-sample_data(physeqFiltered)
  str(filteredMeta) #470 samples
  Data_no2020 <- filteredMeta[!is.na(filteredMeta$SurvivedNextSeason),]
  str(Data_no2020) #398 samples
  SurvivalID<-Data_no2020$Sample.ID

  physeqSurvival<- prune_samples(SurvivalID, physeqFiltered)
  physeqSurvival 
  

###Only keep the latest sample per bird as not enough with multiple samples (keep the same ones as for alpha analysis)
  
  SurvivalData<- sample_data(physeqSurvival)
  latestBird<- SurvivalData[SurvivalData$LatestSample=="yes"] 
  LatestBirdID<- latestBird$Sample.ID
  
  physeqSurvivalBeta<- prune_samples(LatestBirdID, physeqSurvival)
  physeqSurvivalBeta #278 birds
  
  
### remove nestlings (too few samples)
  
  phyBetaSurvivalnoCH<- subset_samples(physeqSurvivalBeta, Ageclass != "N")
  phyBetaSurvivalnoCH #269 samples
  

### remove floaters to be consistent with survival alpha diversity analysis (no assigned territory)
  phyBetaSurvivalnoFloat<- subset_samples(phyBetaSurvivalnoCH, Status != "FLOAT")
  phyBetaSurvivalnoFloat #264 samples
  

#how many survived?
  datSurvival<-data.frame(sample_data(phyBetaSurvivalnoFloat))
  str(datSurvival[datSurvival$SurvivedNextSeasonFactor=="no",]) #38 samples from birds that died
  

  
##########################  
### PRUNE RARE TAXA ######
##########################
  
# Compute prevalence of each feature (total number of samples in which a taxon appears at least once), store as data.frame
  prevdf = apply(X = otu_table(phyBetaSurvivalnoFloat),
                 MARGIN = ifelse(taxa_are_rows(phyBetaSurvivalnoFloat), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts for each phylum to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(phyBetaSurvivalnoFloat),
                      tax_table(phyBetaSurvivalnoFloat))
  head(prevdf)
  str(prevdf)
  
  plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) #average prevelance of features within each phylum and the sum of feature prevalence within each phylum
  
# Plot the unique phyla: each dot will be a feature- total abundance of that feature across samples on x axis and the prevalance (the fraction of all samples it occurs in on the y axis).
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(phyBetaSurvivalnoFloat, "Phylum"))
  ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(phyBetaSurvivalnoFloat),color=Phylum)) +
    # Include filtering threshold here; 2% = approx 5 samples, total abundance across all samples set to 50
    geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_vline(xintercept = 50, alpha = 0.5, linetype = 2)+  geom_point(size = 1, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + theme(legend.position="none")
  
# Define prevalence threshold as in at least 2% of samples and abundance threshold as 50 total reads across samples
 
  prevalenceThreshold<-2*(264/100)
  
  abundanceThreshold<-50

# Execute the prevalence filter, using `prune_taxa()` function
  
  head(prevdf1)
  
  KeepTaxa1<-prevdf1[prevdf1$Prevalence>=prevalenceThreshold,]
  str(KeepTaxa1)
  head(KeepTaxa1)
  
  KeepTaxa2<- rownames(KeepTaxa1)[(KeepTaxa1$TotalAbundance >= abundanceThreshold)]
  str(KeepTaxa2) #2900 taxa retained after filtering at 2% sample level
  
  
  physeqBeta<- prune_taxa(KeepTaxa2, phyBetaSurvivalnoFloat)
  physeqBeta # 2900 taxa (down from 55664)

##################################
##Split up juveniles and adults###
##################################
  
#need to do this to avoid confounding age and survival (as mortality is higher in juvenile age classes- SA, OFL, FL)

  physeqBeta   

#Create age group category 
  
  juv <- c("FL", "OFL", "SA")
    
  sample_data(physeqBeta)$AgeJuvAdult<- ifelse(sample_data(physeqBeta)$Ageclass %in% juv, "juvenile", "adult")
  head(sample_data(physeqBeta))

######################
##### Juveniles ######
######################
  
  JuvenilePhyseqBeta<- subset_samples(physeqBeta, AgeJuvAdult == "juvenile")
  JuvenilePhyseqBeta # 116 samples
  JuvDat<- data.frame(sample_data(JuvenilePhyseqBeta))
  table(JuvDat$SurvivedNextSeasonFactor) #17 died, 99 survived


##CLR transform
  Juvphyseq_clr <- microbiome::transform(JuvenilePhyseqBeta, "clr")
  head(otu_table(Juvphyseq_clr))

#Extract ASV Matrix and Sample Data  
  Juvclr_Matrix<-vegan_otu(Juvphyseq_clr)
  Juvclr_SampleData<-as(sample_data(Juvphyseq_clr),"data.frame")
  str(Juvclr_SampleData)
  Juvclr_SampleData$FieldPeriodIDFactor <- factor(Juvclr_SampleData$FieldPeriodIDFactor, levels= c("Major17", "Minor18", "Major18", "Minor19", "Major19"))
  Juvclr_SampleData$Ageclass <- factor(Juvclr_SampleData$Ageclass, levels= c("FL", "OFL", "SA"))


###Principal Components Analysis

  Juvclr_pca<-rda(Juvclr_Matrix)


### Plot the PCA

#by survival to next season
  plotPCASurvivalJuveniles<-gg_ordiplot(Juvclr_pca,groups=Juvclr_SampleData$SurvivedNextSeasonFactor,ellipse=FALSE, plot=FALSE, pt.size=3, choices=c(1,2)) 
  JuvSurvivalPlot<-plotPCASurvivalJuveniles$plot
  JuvSurvivalPlot + geom_line(alpha=0.1)+ theme_bw()  + guides(color=guide_legend("Survived")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 6.62%") + ylab("PC2 5.51%\n") + theme(axis.text = element_text(size=20), axis.title = element_text(size=20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
    scale_color_manual(values= c("#E69F00","#0072B2")) 
  

####Permanova analysis - Table 2

  perm <- how(nperm = 9999)
  set.seed(5334)
  JuvenilePermanovaSurvival<- adonis2(Juvclr_Matrix ~ SurvivedNextSeasonFactor + Ageclass + 
                                        SexEstimate + TerritoryQuality + FieldPeriodIDFactor, data=Juvclr_SampleData, 
                                      permutations = perm, method = "euclidean", by= "margin")
  JuvenilePermanovaSurvival


######################
######## Adults ######
######################

  AdultPhyseqBeta<- subset_samples(physeqBeta, AgeJuvAdult == "adult")
  AdultPhyseqBeta # 148 samples
  AdultDat<- data.frame(sample_data(AdultPhyseqBeta))
  table(AdultDat$SurvivedNextSeasonFactor) #21 died, 127 survived
  
##CLR transform
  Adultphyseq_clr <- microbiome::transform(AdultPhyseqBeta, "clr")
  head(otu_table(Adultphyseq_clr))
  
#Extract ASV Matrix and Sample Data  
  Adultclr_Matrix<-vegan_otu(Adultphyseq_clr)
  Adultclr_SampleData<-as(sample_data(Adultphyseq_clr),"data.frame")
  str(Adultclr_SampleData)
  Adultclr_SampleData$FieldPeriodIDFactor <- factor(Adultclr_SampleData$FieldPeriodIDFactor, levels= c("Major17", "Minor18", "Major18", "Minor19", "Major19"))
  Adultclr_SampleData$Ageclass <- factor(Adultclr_SampleData$Ageclass, levels= c("FL", "OFL", "SA"))
  
  
###Principal Components Analysis###
  
  Adultclr_pca<-rda(Adultclr_Matrix)
  
  
  ### Plot the PCA- Figure 3
  
  #by survival to next season
  plotPCASurvivalAdults<-gg_ordiplot(Adultclr_pca,groups=Adultclr_SampleData$SurvivedNextSeasonFactor,ellipse =FALSE, plot=FALSE, pt.size=3, choices=c(1,2)) 
  AdultSurvivalPlot<-plotPCASurvivalAdults$plot
  AdultSurvivalPlot + geom_line(alpha=0.1)+ theme_bw()  + guides(color=guide_legend("Survived")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 7.51%") + ylab("PC2 4.01%\n") + theme(axis.text = element_text(size=20), axis.title = element_text(size=20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
    scale_color_manual(values= c("#E69F00","#0072B2")) 
  
  
####Permanova analysis - Table 2
  
  perm <- how(nperm = 9999)
  set.seed(1458)
  AdultPermanovaSurvival<- adonis2(Adultclr_Matrix ~ SurvivedNextSeasonFactor +
                                        SexEstimate + FieldPeriodIDFactor +TerritoryQuality, data=Adultclr_SampleData, 
                                      permutations = perm, method = "euclidean", by= "margin")
  AdultPermanovaSurvival

  
  
### Checking group dispersions with betadisper analysis
  
  #survival
  distMatrixBetaA <- vegdist(Adultclr_Matrix, method="euclidean")
  BDA<-betadisper(distMatrixBetaA, Adultclr_SampleData$SurvivedNextSeason)
  set.seed(25)
  permutest(BDA, permutations = 9999) # not sig
  
  
  
#############################
### PhILR transformation ####
#############################   
    
#includes phylogenetic information: PHYLOGENETIC ISOMETRIC LOG RATIO TRANSFORMATION ####
  
## INSTALLATION

    #if(!requireNamespace("BiocManager")){
    # install.packages("BiocManager")
    #}
    
    #library(BiocManager)
    
    #BiocManager::install("philr")
    
    
    library(philr)
    
    
##################################
##Split up juveniles and adults###
##################################
    
    
######################
##### Juveniles ######
######################
  
  JuvenilePhyseqBeta
    
# Add a pseudocount of 1 to ASVs in the filtered phyloseq object- avoids calculating log-ratios involving zeros
    
    physeqBetaPseudoJuv <- transform_sample_counts(JuvenilePhyseqBeta, function(x) x+1)
    otu_table(physeqBetaPseudoJuv)
    
# Next we check that the phylogenetic tree is rooted and binary (all multichotomies have been resolved). 
    
    is.rooted(phy_tree(physeqBetaPseudoJuv)) # Is the tree Rooted? Yes
    is.binary(phy_tree(physeqBetaPseudoJuv)) # All multichotomies resolved? - false
    physeqBetaTreeJuv <- multi2di(phy_tree(physeqBetaPseudoJuv)) # use multi2di in ape to replace multichotomies with a series of dichotomies with one (or several) branch(es) of zero length
    
#number the internal nodes of the tree so it's easier to work with- prefix the node number with "n"
    physeqBetaTreeNJuv <- makeNodeLabel(physeqBetaTreeJuv, method="number", prefix='n')
    name.balance(physeqBetaTreeNJuv, tax_table(physeqBetaPseudoJuv), 'n86') #find a consensus name for the two clades that descend from a given balance (out of x/y, x= numerator, y= denominator in balance)
    
#transpose data: ensures taxa are columns, samples are rows
    otu.table <- t(otu_table(physeqBetaPseudoJuv))
    tree <- physeqBetaTreeNJuv
    metadata <- sample_data(physeqBetaPseudoJuv)
    tax <- tax_table(physeqBetaPseudoJuv)
    
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
    
    
## Ordination in PhILR Space
    
    Juvilr_pca<-rda(phy.philr)
    
    library(ggordiplots)
    metadata$FieldPeriodIDFactor <- factor(metadata$FieldPeriodIDFactor, levels= c("Major17", "Minor18", "Major18", "Minor19", "Major19"))
    
    #survival
    JuvplotPCAphilSurv<-gg_ordiplot(Juvilr_pca,groups=metadata$SurvivedNextSeasonFactor,plot=FALSE, ellipse=FALSE, pt.size=3) 
    JuvplotPCAphilSurv<-JuvplotPCAphilSurv$plot
    JuvplotPCAphilSurv + theme_bw()  + guides(color=guide_legend("Survived")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      xlab("\nPC1 21.91%") + ylab("PC2 11.34%\n") + theme(axis.text = element_text(size=20), axis.title = element_text(size=20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
      scale_color_manual(values= c("#E69F00","#0072B2"))
    
####Permanova analysis - Table 2
    
    JuvPhilrDat<- data.frame(metadata)
    JuvPhilrDat$SurvivedNextSeasonFactor<- as.factor(JuvPhilrDat$SurvivedNextSeasonFactor)
    JuvPhilrDat$SexEstimate<- as.factor(JuvPhilrDat$SexEstimate)
    JuvPhilrDat$Ageclass<- as.factor(JuvPhilrDat$Ageclass)
    str(JuvPhilrDat)
    perm <- how(nperm = 9999)
    set.seed(9934)
    JuvenilePermanovaSurvival_philr<- adonis2(phy.philr ~ SurvivedNextSeasonFactor + Ageclass + 
                                          SexEstimate + TerritoryQuality + FieldPeriodIDFactor, data=JuvPhilrDat, 
                                        permutations = perm, method = "euclidean", by= "margin")
    JuvenilePermanovaSurvival_philr
    
    
######################
##### Adults #########
######################
    
    AdultPhyseqBeta # 148 samples
    
# Add a pseudocount of 1 to ASVs in the filtered phyloseq object- avoids calculating log-ratios involving zeros
    
    physeqBetaPseudoAdult <- transform_sample_counts(AdultPhyseqBeta, function(x) x+1)
    otu_table(physeqBetaPseudoAdult)
    
# Next we check that the phylogenetic tree is rooted and binary (all multichotomies have been resolved). 
    
    is.rooted(phy_tree(physeqBetaPseudoAdult)) # Is the tree Rooted? Yes
    is.binary(phy_tree(physeqBetaPseudoAdult)) # All multichotomies resolved? - false
    physeqBetaTreeAdult <- multi2di(phy_tree(physeqBetaPseudoAdult)) # use multi2di in ape to replace multichotomies with a series of dichotomies with one (or several) branch(es) of zero length
    
#number the internal nodes of the tree so it's easier to work with- prefix the node number with "n"
    physeqBetaTreeNAdult <- makeNodeLabel(physeqBetaTreeAdult, method="number", prefix='n')
    name.balance(physeqBetaTreeNAdult, tax_table(physeqBetaPseudoAdult), 'n86') #find a consensus name for the two clades that descend from a given balance (out of x/y, x= numerator, y= denominator in balance)
    
#transpose data: ensures taxa are columns, samples are rows
    Adultotu.table <- t(otu_table(physeqBetaPseudoAdult))
    Adulttree <- physeqBetaTreeNAdult
    Adultmetadata <- sample_data(physeqBetaPseudoAdult)
    Adulttax <- tax_table(physeqBetaPseudoAdult)
    
    Adultotu.table[1:2,1:2] # OTU Table
    Adulttree # Phylogenetic Tree
    head(Adultmetadata,2) # Metadata
    head(Adulttax,2) # taxonomy table
    
#transform the data using philr (wrapper function)
#uses taxon weighting:down weight the influence of taxa with many zeros or non-zero counts
#branch length weighting: scale balances using phylogenetic distance (communities differing in abundance of closely related species are more similar than those differing in distantly related species)
    Adultphy.philr <- philr(Adultotu.table, Adulttree, 
                       part.weights='enorm.x.gm.counts', 
                       ilr.weights='blw.sqrt')
    Adultphy.philr[1:5,1:5] #expressed as balances: these represent the log-ratio of the geometric mean abundance of the two groups of taxa that descend from a given internal node
    
    
## Ordination in PhILR Space
    
    Adultilr_pca<-rda(Adultphy.philr)
    
    library(ggordiplots)
    Adultmetadata$FieldPeriodIDFactor <- factor(Adultmetadata$FieldPeriodIDFactor, levels= c("Major17", "Minor18", "Major18", "Minor19", "Major19"))
    
    #survival
    AplotPCAphilSurv<-gg_ordiplot(Adultilr_pca,groups=Adultmetadata$SurvivedNextSeasonFactor,plot=FALSE, ellipse = FALSE, pt.size=2) 
    AplotPCAphilSurv<-AplotPCAphilSurv$plot
    AplotPCAphilSurv + theme_bw()  + guides(color=guide_legend("Survived")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      xlab("\nPC1 21.86%") + ylab("PC2 10.83%\n") + theme(axis.text = element_text(size=20), axis.title = element_text(size=20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
      scale_color_manual(values= c("#E69F00","#0072B2"))
    
####Permanova analysis - Table 2
    AdultPhilrDat<- data.frame(Adultmetadata)
    AdultPhilrDat$SurvivedNextSeasonFactor<- as.factor(AdultPhilrDat$SurvivedNextSeasonFactor)
    AdultPhilrDat$SexEstimate<- as.factor(AdultPhilrDat$SexEstimate)
    str(AdultPhilrDat)
    perm <- how(nperm = 9999)
    set.seed(1359)
    AdultPermanovaSurvival_philr<- adonis2(Adultphy.philr ~ SurvivedNextSeasonFactor +
                                                SexEstimate + FieldPeriodIDFactor +TerritoryQuality, data=AdultPhilrDat, 
                                              permutations = perm, method = "euclidean", by= "margin")
    AdultPermanovaSurvival_philr
    
    
    
    
#######################################    
#### Differential abundance testing ###
####################################### 

##### Adults only - juveniles did not show compositional differences  
##### Use the non-rarefied (but prevalence filtered) reads.    
##### Using Ancom BC https://www.nature.com/articles/s41467-020-17041-7
    
    library(dplyr)
    library(nloptr)
    
#install latest release of AncomBC
    #install.packages("remotes")
    #remotes::install_github("FrederickHuangLin/ANCOMBC")
    
    library(ANCOMBC)
    library(corncob)
    
    AdultPhyseqBeta # unrarefied reads, but filtered to remove rare taxa: 2900 taxa, 270 samples
    
    str(sample_data(AdultPhyseqBeta))
    sample_data(AdultPhyseqBeta)$SurvivedNextSeasonFactor <- as.factor(sample_data(AdultPhyseqBeta)$SurvivedNextSeasonFactor)
    sample_data(AdultPhyseqBeta)$Ageclass <- as.factor(sample_data(AdultPhyseqBeta)$Ageclass)
    sample_data(AdultPhyseqBeta)$SexEstimate <- as.factor(sample_data(AdultPhyseqBeta)$SexEstimate)
    sample_data(AdultPhyseqBeta)$FieldPeriodID <- as.factor(sample_data(AdultPhyseqBeta)$FieldPeriodID)
    
# running Ancom BC
    
    SurvivalAncom<-ancombc(AdultPhyseqBeta, formula= "SurvivedNextSeasonFactor + SexEstimate + FieldPeriodID", p_adj_method = "BH", lib_cut = 10000, alpha = 0.05, neg_lb = TRUE, group = "FieldPeriodID", struc_zero = TRUE,global = TRUE)
    str(SurvivalAncom)
    
    results<-SurvivalAncom$res
    str(results)
    head(results)
    
    
    diffAbund<-as.data.frame(results$diff_abn)
    head(diffAbund)
    str(diffAbund)
    DiffSurvival<-diffAbund[diffAbund$SurvivedNextSeasonFactoryes=="TRUE",]
    str(DiffSurvival)
    rows<-row.names(DiffSurvival)
    
    results2<-as.data.frame(results)
    head(results2)
    significantSurvival<-subset(results2, rownames(results2) %in% rows)
    head(significantSurvival)
    str(significantSurvival) # 28 observations
    
    taxon<-otu_to_taxonomy(OTU = row.names(DiffSurvival), data = physeqBeta)
    taxon2<-as.data.frame(taxon)
    
    significantSurvival$taxonomy<-taxon2$taxon
    #write.csv(significantSurvival,"AncomSurvivalAdult.csv")
    significantSurvival<- read.csv("AncomSurvivalAdult.csv")
    str(significantSurvival)

  
    significantSurvival<-data.frame(significantSurvival$beta.SurvivedNextSeasonFactoryes,significantSurvival$Order, significantSurvival$Phylum,  significantSurvival$Family, significantSurvival$se.SurvivedNextSeasonFactoryes)
    str(significantSurvival)
    colnames(significantSurvival)<-c("Beta", "Order", "Phylum", "Family", "se")
    significantSurvival$Order<-as.factor(significantSurvival$Order)
    significantSurvival$Phylum<-as.factor(significantSurvival$Phylum)
    significantSurvival$Family<-as.factor(significantSurvival$Family)
    significantSurvival$Species<-as.factor(seq(from=1, to=28, by=1))
    
## plot of all differentially abundant taxa- Figure 4
    
    ggplot(significantSurvival, aes(x= Beta, y=Order, colour=Phylum, group=Species)) + 
      geom_point(position=position_dodge(width = 0.8), size=2) +
      geom_errorbar(aes(xmin=Beta -se, xmax=Beta + se), width=0,lwd=1, position= position_dodge(width=0.8)) +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_x_continuous(name="\nLog fold change", breaks=c(-3,-2, -1,0, 1, 2)) + 
      ylab("Bacterial order\n")  + geom_vline(xintercept = 0, linetype="dotted", size=1.5) +
      theme(axis.text = element_text(size=20), axis.title = element_text(size=22), legend.text = element_text(size=20), legend.title = element_text(size=22)) +
      scale_colour_manual(values=c("goldenrod2","#CC79A7", "#56B4E9",  "aquamarine4", "darkgrey"))
    
    
# plotting differentially abundant ASVs individually  
    my_subset <- subset(otu_table(physeqBeta), rownames(otu_table(physeqBeta)) %in% c(rows))
    new_physeq <- merge_phyloseq(my_subset, tax_table(physeqBeta), sample_data(physeqBeta))
    mphyseq = psmelt(new_physeq)
    ggplot(data = mphyseq, mapping = aes_string(x = "SurvivedNextSeasonFactor",y = "Abundance")) +
      geom_boxplot() +
      geom_point(size = 1, alpha = 0.3,
                 position = position_jitter(width = 0.3)) +
      scale_y_log10()+ facet_wrap(facets = vars(Genus))+
      theme(legend.position="none")
    