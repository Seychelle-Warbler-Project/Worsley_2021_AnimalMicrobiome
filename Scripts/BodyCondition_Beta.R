#Sarah Worsley

#R code and analysis for: Gut microbiome composition, not alpha diversity, is associated with survival in a natural vertebrate population

#load packages
  library(ggplot2)
  library(phyloseq)
  library(ape)
  library(vegan)
  library(dplyr)
  library(BiodiversityR)
  library(data.table)
  library(microbiome)
  library(ggordiplots)

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

### remove samples with less than 10000 reads - did not reach completion- see SequencingExtractionRepeatability.R for iNext Code
  physeq2
  physeq3<-prune_samples(sample_sums(physeq2)>=10000, physeq2)
  physeq3 # removes 23 samples (some of these are sequencing duplicates)- leaves 587 samples

  
###############################
#### REMOVE SEQ DUPLICATES ####
###############################

#See SequencingExtractionRepeatability.R for assessment of repeatability
  
#### extract samples sequenced/extracted twice (keep the one with the highest read number)###
  
  physeqNoDup<- subset_samples(physeq3, SequencingDuplicate=="No") #528 remaining
  physeqNoDup
  
#### Remove catch duplicates ####
  
  # Keeping 1 sample per catch in following order of priority: 1) from tray 2)  from bag, or if neither keep highest read count
  
  physeqFiltered<- subset_samples(physeqNoDup, CatchRepeat=="No") #471 samples, 55664 taxa 
  physeqFiltered
  
  
##############################
### Sample filtering #########
##############################
  
#remove outlier sample (Sample.71.81)- very tiny sample (see extraction notes) with extremely low alpha diversity

  physeqFiltered<- subset_samples( physeqFiltered, Sample.ID != "Sample.71.81")
  physeqFiltered #470 samples, 55664 taxa  

#remove nestlings
  physeqFiltered<- subset_samples( physeqFiltered,Ageclass != "N")
  physeqFiltered #458 samples


#remove samples where no body mass is recorded
  physeqFiltered<- subset_samples( physeqFiltered,BodyMass != "NA")#leaves 436 samples
  physeqFiltered
  
#remove females with eggs as outliers for mass
  physeqFiltered<- subset_samples( physeqFiltered, Sample.ID != "Sample.40.46")  #Notes say possible egg present
  physeqFiltered<- subset_samples( physeqFiltered, Sample.ID != "Sample.48.54") #egg present
  physeqFiltered # 434 samples

#remove floaters with no TQ
  physeqFiltered<- subset_samples( physeqFiltered,Status!= "FLOAT")#leaves 426 samples
  physeqFiltered #426 samples, from 296 individuals

#remove bird with no catch time 
  physeqFiltered<- subset_samples(physeqFiltered,MinutesSinceSunrise != "NA")#leaves 425 samples
  physeqFiltered #425 samples, from 296 individuals
  

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
  # Include filtering threshold here; 2% = approx 5 samples, total abundance across all samples set to 50
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_vline(xintercept = 50, alpha = 0.5, linetype = 2)+  geom_point(size = 1, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# Define prevalence threshold as 2% samples and abundance as 50 total reads across samples

  prevalenceThreshold<-2*(425/100)

  abundanceThreshold<-50

# Execute the prevalence filter, using `prune_taxa()` function

  head(prevdf1)

  KeepTaxa1<-prevdf1[prevdf1$Prevalence>=prevalenceThreshold,]
  str(KeepTaxa1)
  head(KeepTaxa1)

  KeepTaxa2<- rownames(KeepTaxa1)[(KeepTaxa1$TotalAbundance >= abundanceThreshold)]
  str(KeepTaxa2) #3080 taxa retained after filtering at 2% sample level


  physeqBeta<- prune_taxa(KeepTaxa2,physeqFiltered)
  physeqBeta # 3080 taxa (down from 56029)

  
#########################################################################
#########  Generate residuals/categories for Body Condition #############
#########################################################################
    
#Make category for juveniles (less than one year old, i.e FL, OFL and SA) and adults
#This is so that you can test for compositional differences without confounding condition and age (as body mass is greater in adults)
  
  Juveniles<- c("FL", "OFL", "SA")
  sample_data(physeqBeta)$Ageclass2Way<- ifelse(sample_data(physeqBeta)$Ageclass %in% Juveniles, "Juveniles", "Adults")
  
  PhySampleData<- data.frame(sample_data(physeqBeta))
  str(PhySampleData)
  ggplot(PhySampleData, aes(Ageclass, BodyMass)) + geom_boxplot() + geom_point(position=position_jitter(width = 0.2))
  ggplot(PhySampleData, aes(Ageclass2Way, BodyMass)) + geom_boxplot() + geom_point(position=position_jitter(width = 0.2))
  
  
# generate body condition measure (residuals from model relating body mass and tarsus length)
  #do separately for the 2 sexes

  ggplot(PhySampleData, aes(RightTarsus, BodyMass, col=SexEstimate)) + geom_point()
  
  malesMass<- PhySampleData[PhySampleData$SexEstimate=="M",]
  femalesMass<- PhySampleData[PhySampleData$SexEstimate=="F",]
  
  #model for males- extract residuals
  modelResidualsMale<- lm(BodyMass~ RightTarsus + MinutesSinceSunrise, data= malesMass)
  summary(modelResidualsMale)
  residualsMales<- data.frame(residuals(modelResidualsMale))
  residualsMales$Sample.ID<- row.names(residualsMales)
  colnames(residualsMales)<- c("ConditionResiduals", "Sample.ID")
  str(residualsMales)
  
  #model for females- extract residuals
  modelResidualsFemale<- lm(BodyMass~ RightTarsus + MinutesSinceSunrise, data= femalesMass)
  summary(modelResidualsFemale)
  residualsFemales<- data.frame(residuals(modelResidualsFemale))
  residualsFemales$Sample.ID<- row.names(residualsFemales)
  colnames(residualsFemales)<- c("ConditionResiduals", "Sample.ID")
  str(residualsFemales)
  
  #merge and add to the metadata
  residualsAll<- rbind(residualsMales, residualsFemales)
  str(residualsAll)  
  
  PhySampleData2<- merge(PhySampleData, residualsAll, by="Sample.ID")
  str(PhySampleData2)

  ggplot(PhySampleData2, aes(Ageclass2Way, ConditionResiduals)) + geom_boxplot() + geom_point(position = position_jitter(width=0.2))

  
##############################################################
#######  PERMANOVA ANALYSIS - CLR Transformed Abundances #####
##############################################################  
  
### split up analyses for juveniles and adults to avoid influence of conflating condition and age on GM composition
  
  
#############################  
#### Juveniles permanova ####
#############################
  
  JuvenilesOnlyPhyseq<- subset_samples(physeqBeta, Ageclass2Way=="Juveniles")
  JuvenilesOnlyPhyseq #205 samples
    
  JuvenilesConditionDat<- data.frame(sample_data(JuvenilesOnlyPhyseq))
  unique(JuvenilesConditionDat$BirdID) #175 individuals
  freq<-data.frame(table(JuvenilesConditionDat$BirdID))
  str(freq[freq$Freq>1,]) #25 individuals have repeat samples

  JuvenilesConditionDat<- merge(JuvenilesConditionDat, residualsAll, by="Sample.ID") #add residuals
  str(JuvenilesConditionDat)
  

#### CLR TRANSFORM ASV ABUNDANCES#####

  juveniles_clr <- microbiome::transform(JuvenilesOnlyPhyseq, "clr")
  head(otu_table(juveniles_clr))
  
# function to extract an ASV matrix
  vegan_otu <- function(physeq) {
    OTU <- otu_table(physeq)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  
# Extract ASV Matrix and Sample Data  
  Juvenilesclr_Matrix<-vegan_otu(juveniles_clr)
  
#### PERMANOVA ####

#Include Bird ID as a blocking factor to control for repeat samples
  
  perm <- how(nperm = 9999)
  set.seed(5498)
  setBlocks(perm) <- with(JuvenilesConditionDat, BirdID)
  permanovaJuvenilesCondition<- adonis2(Juvenilesclr_Matrix ~ Ageclass + SexEstimate +  TerritoryQuality +FieldPeriodIDFactor +
                                    ConditionResiduals, data=JuvenilesConditionDat, 
                                  permutations = perm, method = "euclidean", by= "margin")
  permanovaJuvenilesCondition
  

 
###########################
##### Adults permanova ####
###########################

  AdultsOnlyPhyseq<- subset_samples(physeqBeta, Ageclass2Way=="Adults")
  AdultsOnlyPhyseq #220 samples
  
  AdultsConditionDat<- data.frame(sample_data(AdultsOnlyPhyseq))
  str(unique(AdultsConditionDat$BirdID)) #165 individuals
  
  AdultsConditionDat<- merge(AdultsConditionDat, residualsAll, by="Sample.ID") #add residuals
  str(AdultsConditionDat)
  
  
#### CLR TRANSFORM ASV ABUNDANCES#####
  
  Adults_clr <- microbiome::transform(AdultsOnlyPhyseq, "clr")
  head(otu_table(Adults_clr))
  
#Extract ASV Matrix and Sample Data  
  Adultsclr_Matrix<-vegan_otu(Adults_clr)
  
#### PERMANOVA ####
  
  perm <- how(nperm = 9999)
  set.seed(549899)
  setBlocks(perm) <- with(AdultsConditionDat, BirdID)
  permanovaAdultsCondition<- adonis2(Adultsclr_Matrix ~ SexEstimate + TerritoryQuality + FieldPeriodIDFactor +
                                          ConditionResiduals, data=AdultsConditionDat, 
                                        permutations = perm, method = "euclidean", by= "margin")
  permanovaAdultsCondition
  
  
################################################################
###  PERMANOVA ANALYSIS - Philr Transformed Abundances #########
################################################################  

############################# 
##### Juveniles permanova ###
#############################
  
#### PhILR TRANSFORM ASV ABUNDANCES #####
  
  library(philr)
  
  # Add a pseudocount of 1 to ASVs in the filtered phyloseq object- avoids calculating log-ratios involving zeros
  
  JuvenilesPseudo <- transform_sample_counts(JuvenilesOnlyPhyseq, function(x) x+1)
  otu_table(JuvenilesPseudo)
  
  # Next we check that the phylogenetic tree is rooted and binary (all multichotomies have been resolved). 
  
  is.rooted(phy_tree(JuvenilesPseudo)) # Is the tree Rooted? Yes
  is.binary(phy_tree(JuvenilesPseudo)) # All multichotomies resolved? - false
  JuvenilesTree <- multi2di(phy_tree(JuvenilesPseudo)) # use multi2di in ape to replace multichotomies with a series of dichotomies with one (or several) branch(es) of zero length
  
  #number the internal nodes of the tree so it's easier to work with- prefix the node number with "n"
  JuvenilesTreeN <- makeNodeLabel(JuvenilesTree, method="number", prefix='n')
  name.balance(JuvenilesTreeN, tax_table(JuvenilesPseudo), 'n856') #find a consensus name for the two clades that descend from a given balance (out of x/y, x= numerator, y= denominator in balance)
  
  #transpose data: ensures taxa are columns, samples are rows
  otu.table <- t(otu_table(JuvenilesPseudo))
  tree <- JuvenilesTreeN
  metadata <- sample_data(JuvenilesConditionDat)
  tax <- tax_table(JuvenilesPseudo)
  
  otu.table[1:2,1:2] # OTU Table
  tree # Phylogenetic Tree
  head(metadata,2) # Metadata
  head(tax,2) # taxonomy table
  
  #transform the data using philr (wrapper function)
  #uses taxon weighting:down weight the influence of taxa with many zeros or non-zero counts
  #branch length weighting: scale balances using phylogenetic distance (communities differing in abundance of closely related species are more similar than those differing in distantly related species)
  Juvenilephy.philr <- philr(otu.table, tree, 
                     part.weights='enorm.x.gm.counts', 
                     ilr.weights='blw.sqrt')
  Juvenilephy.philr[1:5,1:5] #expressed as balances: these represent the log-ratio of the geometric mean abundance of the two groups of taxa that descend from a given internal node
  
  
### permanova ###
  
  perm <- how(nperm = 9999)
  set.seed(7898)
  setBlocks(perm) <- with(JuvenilesConditionDat, BirdID)
  permanovaJuvenilesPhilr<- adonis2(Juvenilephy.philr ~ Ageclass + SexEstimate + + TerritoryQuality+ FieldPeriodIDFactor +
                                        ConditionResiduals, data=JuvenilesConditionDat, 
                                        permutations = perm, method = "euclidean", by= "margin")
  permanovaJuvenilesPhilr
  

############################ 
##### Adults permanova #####
############################  
  
#### Philr TRANSFORM ASV ABUNDANCES#####
  
  
  # Add a pseudocount of 1 to ASVs in the filtered phyloseq object- avoids calculating log-ratios involving zeros
  
  AdultsPseudo <- transform_sample_counts(AdultsOnlyPhyseq, function(x) x+1)
  otu_table(AdultsPseudo)
  
  # Next we check that the phylogenetic tree is rooted and binary (all multichotomies have been resolved). 
  
  is.rooted(phy_tree(AdultsPseudo)) # Is the tree Rooted? Yes
  is.binary(phy_tree(AdultsPseudo)) # All multichotomies resolved? - false
  AdultsTree <- multi2di(phy_tree(AdultsPseudo)) # use multi2di in ape to replace multichotomies with a series of dichotomies with one (or several) branch(es) of zero length
  
  #number the internal nodes of the tree so it's easier to work with- prefix the node number with "n"
  AdultsTreeN <- makeNodeLabel(AdultsTree, method="number", prefix='n')
  name.balance(AdultsTreeN, tax_table(AdultsPseudo), 'n856') #find a consensus name for the two clades that descend from a given balance (out of x/y, x= numerator, y= denominator in balance)
  
  #transpose data: ensures taxa are columns, samples are rows
  otu.table <- t(otu_table(AdultsPseudo))
  tree <- AdultsTreeN
  metadata <- sample_data(AdultsConditionDat)
  tax <- tax_table(AdultsPseudo)
  
  otu.table[1:2,1:2] # OTU Table
  tree # Phylogenetic Tree
  head(metadata,2) # Metadata
  head(tax,2) # taxonomy table
  
  #transform the data using philr (wrapper function)
  #uses taxon weighting:down weight the influence of taxa with many zeros or non-zero counts
  #branch length weighting: scale balances using phylogenetic distance (communities differing in abundance of closely related species are more similar than those differing in distantly related species)
  Adultphy.philr <- philr(otu.table, tree, 
                             part.weights='enorm.x.gm.counts', 
                             ilr.weights='blw.sqrt')
  Adultphy.philr[1:5,1:5] #expressed as balances: these represent the log-ratio of the geometric mean abundance of the two groups of taxa that descend from a given internal node
  
  
### PERMANOVA adults ###
  
  perm <- how(nperm = 9999)
  set.seed(7898)
  setBlocks(perm) <- with(AdultsConditionDat, BirdID)
  permanovaAdultPhilr<- adonis2(Adultphy.philr ~ SexEstimate +  TerritoryQuality+ FieldPeriodIDFactor +
                                      ConditionResiduals, data=AdultsConditionDat, 
                                    permutations = perm, method = "euclidean", by= "margin")
  permanovaAdultPhilr
  
  
