# Sarah Worsley

# R code and analysis for: Gut microbiome composition, not alpha diversity, is associated with survival in a natural vertebrate population

#load packages
  library(ggplot2)
  library(phyloseq)
  library(ape)
  library(dplyr)
  library(FSA)
  library(iNEXT)
  library(boot)
  library(microbiome)
  library(vegan)

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
  row.names(metadatafull)<- metadatafull$Sample.ID
  
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
  physeq 

#check the number of reads per sample
  sample_reads<-data.frame(reads=sample_sums(physeq))
  head(sample_reads)

##################
### FILTERING ####
##################

#filter to remove instances where features are not assigned as bacteria/ unassigned/ not assigned at phylum level/ #####

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

#############################
#### SAMPLE COMPLETENESS#####
#############################
  
#Completeness alpha rarefaction plotting (species accumulation) using iNEXT

  head(otu_table(physeq2))

#make the otu abundance table into a matrix and check by printing first 2 rows/columns
  abund <- as(otu_table(physeq2), "matrix") 
  abund[1:2,1:2]

#convert to a dataframe
  abund2 <- as.data.frame(abund)
  str(abund2)

#iNEXT only takes numeric values, so change all values in the dataframe to numeric values instead of integers.
  df2 <- mutate_all(abund2, function(x) as.numeric(x)) 
  str(df2)

# Install and load inext
#iNEXT= https://cran.r-project.org/web/packages/iNEXT/iNEXT.pdf
#install.packages("iNEXT")

  library(iNEXT)

#Run inext
#q=0 specifies that the function should use species richness (rather than shannon (1) or simpson (2) indices) for the rarefaction. 
#Datatype= abundance because you have raw abundances. 
#endpoint=20000 specifies the number of reads (sample size) that should be used as an endpoint for the rarefaction/extrapolation. 

  inext<-iNEXT(df2, q=0, datatype="abundance", endpoint=20000)

#plot rarefaction curve
  rarefaction<- ggiNEXT(inext, type=1, se=TRUE, facet.var="none", grey= TRUE) + theme(legend.position = "none")+ xlab("Sequencing depth") + ylab("Observed ASVs")
  rarefaction

#Plot sample completeness curve - Figure S1
  completeness<-ggiNEXT(inext, type=2, se=TRUE, facet.var="none", grey=TRUE)+scale_shape_manual(values=rep(20,45))+ theme(legend.position = "none") +xlab("Read count") +ylab("Sample completeness")
  completeness + geom_vline(xintercept=9000, alpha=0.5, linetype=2)


### remove samples with less than 10000 reads ###
  physeq2
  physeq3<-prune_samples(sample_sums(physeq2)>=10000, physeq2)
  physeq3 # removes 23 samples (note some of these are sequencing duplicates)- leaves 587 samples

  
#Find mean number of ASVs across samples prior to rarefaction (don't include outlier in this calculation)
  physeq3noOutlier <-subset_samples(physeq3, Sample.ID != "Sample.71.81")
  physeq3noOutlier #586 samples once outlier removed (unusual diversity measure as very small sample)

  range(sample_sums(physeq3noOutlier)) #range of read numbers across samples
  observationThreshold <- 1
  sumASVs<- (apply(otu_table(physeq3noOutlier) > observationThreshold, 2, sum))
  str(sumASVs)
  sumASVs <- data.frame(sumASVs)
  str(sumASVs)
  mean(sumASVs[,1]) #mean ASVs per sample
  sd(sumASVs[,1]) #Standard deviation of ASVs per sample
 
  
#################################################  
###Repeatability of alpha diversity measures ####
#################################################  
  

##### RAREFY READS TO MIN SAMPLING DEPTH ######
  
#rarefy to 10000 and set seed before rarefying (28367), so that results are reproducible  
  physeqRare<-rarefy_even_depth(physeq3, 10000, rngseed = 28367)

# 6548 ASVs removed after subsampling- leaves 49116 across 587 samples
  physeqRare
  sample_sums(physeqRare)

# Remove one outlier sample "really tiny" in extraction notes= 586 samples

  physeqRare <-subset_samples(physeqRare, Sample.ID != "Sample.71.81")
  physeqRare


##### CALCULATE ALPHA DIVERSITY METRICS ####

# calculate shannon diversity using estimate_richness()
  richnessEstRare<-estimate_richness(physeqRare, split=TRUE, measures= c("Chao1", "Shannon", "observed"))
  head(richnessEstRare)
  str(richnessEstRare)


#add alpha diversity metrics to metadata
  physeqRareMeta <- as.data.frame(sample_data(physeqRare))
  head(physeqRareMeta)
  physeqRareMeta$Chao1 <- richnessEstRare$Chao1
  physeqRareMeta$Shannon <- richnessEstRare$Shannon
  #write.csv(physeqRareMeta, "MetaForMatrix.csv")

###Extract samples sequenced/extracted twice ######

  MetaForMatrix<- read.csv("MetaForMatrix.csv", row.names = 1)
  head(MetaForMatrix)
  str(MetaForMatrix)

#extract those sequenced twice - make a column in the excel sheet entitled comparison and add these as sequence duplicates
  row_namesSeqDup <- as.vector(grep("b", rownames(physeqRareMeta), value=TRUE)) 
  str(row_namesSeqDup) #49 seq reps- 6 were filtered out as low read no (<10000).

#extract samples with DNA extracted twice - add these as extraction duplicates to excel
  row_namesExDup <- as.vector(grep("R", rownames(physeqRareMeta), value=TRUE))
  str(row_namesExDup) #10 extraction reps

#Filter metadata to only include duplicated tube numbers  
  Dups<- c(row_namesExDup,row_namesSeqDup)
  str(Dups) #59
  
  DupsMeta<- MetaForMatrix[MetaForMatrix$Sample.ID %in% Dups,]
  TubeNosDups<-as.vector(DupsMeta$TubeNumber) # get tube numbers of all duplicates
  
  TubesForMatrix<- MetaForMatrix[MetaForMatrix$TubeNumber %in% TubeNosDups,] #subset to just those samples that are duplicated
  str(TubesForMatrix)
  Tubenos<- TubesForMatrix[,c(2,5)]
  str(Tubenos)

###### calculate the reproducibility of alpha diversity metrics #####

#Make a distance matrix with euclidean distances of shannon diversity
  shannon <- TubesForMatrix$Shannon
  samples <- TubesForMatrix$OriginalSampleIDUnique
  names(shannon) <- samples
  library(vegan)
  distMatrix <- vegdist(shannon, method="euclidean")
  distMatrix <- as.matrix(distMatrix)
  distMatrix[1:3,1:3]


#extract as a dataframe (just upper triangle- one way comparisons)
  pairwiseDist <- t(combn(colnames(distMatrix), 2))
  pairwiseDist<- data.frame(pairwiseDist, dist=distMatrix[pairwiseDist])
  head(pairwiseDist)
  colnames(pairwiseDist)<- c("OriginalSampleIDUnique","ID2", "dist")
  str(pairwiseDist)
  pairwiseDist<- merge(pairwiseDist,Tubenos,by="OriginalSampleIDUnique")
  colnames(pairwiseDist)<- c("ID1","OriginalSampleIDUnique", "dist", "Tube1")
  pairwiseDist<- merge(pairwiseDist,Tubenos,by="OriginalSampleIDUnique")
  colnames(pairwiseDist)<- c("ID2","ID1", "dist", "Tube1", "Tube2")
  
  
#Copy the two Id columns twice so that we can check for matches: ID1seq, ID2seq, ID1extr, ID2extr
  pairwiseDist$ID1Seq<- pairwiseDist$ID1
  pairwiseDist$ID2Seq<- pairwiseDist$ID2
  pairwiseDist$ID1extr<- pairwiseDist$ID1
  pairwiseDist$ID2extr<- pairwiseDist$ID2
  str(pairwiseDist)
  
#Remove "b" from two of them (ID1seq,ID2seq) to check for tube number matches (the second duplicate is the tube number+b) between seq dup
  pairwiseDist$ID1Seq<-gsub("b","",as.character(pairwiseDist$ID1Seq))
  pairwiseDist$ID2Seq<-gsub("b","",as.character(pairwiseDist$ID2Seq))
  
#remove "R" from ID1extr and ID2extr to check for extraction replicates (the second extraction has R after the tube number)
  pairwiseDist$ID1extr<-gsub("R","",as.character(pairwiseDist$ID1extr))
  pairwiseDist$ID2extr<-gsub("R","",as.character(pairwiseDist$ID2extr))
  
  str(pairwiseDist)

  
#name type of comparison

  pairwiseDist$SeqDup<-ifelse(pairwiseDist$ID1Seq==pairwiseDist$ID2Seq, "SeqDuplicate", "Between")
  pairwiseDist$extrRep<-ifelse(pairwiseDist$ID1extr==pairwiseDist$ID2extr, "ExDuplicate", "Between")
  
  pairwiseDist$Comparison<- pairwiseDist$SeqDup
  which(pairwiseDist$extrRep == "ExDuplicate")
  pairwiseDist$Comparison[which(pairwiseDist$extrRep == "ExDuplicate")] <- "ExtrDuplicate"

#plot the distances- Figure S2
  pairs<- pairwiseDist
  head(pairs)
  str(pairs)
  pairs$Comparison<-as.factor(pairs$Comparison)
  samplePairs<- c("Different samples", "Extraction duplicates", "Sequencing duplicates")
  ggplot(pairs, aes(Comparison, dist)) +geom_jitter(col="darkgrey", width=0.2, shape=21, size=2)+ 
    geom_boxplot(alpha=0.1, outlier.shape=NA) + stat_summary(fun=mean, geom="point", shape=23, size=4)+
    xlab("\nPairwise comparison") +ylab("Euclidean distance\n") +theme_minimal() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +theme (axis.line = element_line(colour="black"))+
    scale_x_discrete(labels=samplePairs) + theme(element_blank())


####statistical test- bootstrapped Kruskal-wallis statified by tube number (i.e sample ID)####
  
  str(pairs)
  head(pairs)

  
  library(boot)
  
  Kruskal_BootstrapAlpha <- function(d,i) {
    d2 <- d[i,]
    objectKW<- (kruskal.test(d2$dist~d2$Comparison))
    return(c(p= objectKW$p.value, stat =objectKW$statistic,
             df = objectKW$parameter))
  }

  
  set.seed(678)
  bootKW<- boot(data = pairs,statistic = Kruskal_BootstrapAlpha, R = 99999, 
                parallel = "multicore", strata = pairs$Tube1, ncpus=2)
  print(bootKW) #t1 is P value, t2 is test stat, t3 is df
  str(bootKW)
  mean(bootKW$t[,1]) #bootstrap estimate of P value = 2.560999e-13
  mean(bootKW$t[,2]) #bootstrap estimate of test statistic (chisq)= 99.32
  mean(bootKW$t[,3]) #df = 2
  
  
  # output: original is the same as $t0 (original value of statistic in full dataset)
  # bias is difference between the mean of bootstrap realizations (mean from $t- called a bootstrap estimate of T) and value in original dataset (the one from $t0).
  # bootMed is the median of the bootstrapped values
  
  #calculate confidence intervals
  bootKW_CI<-boot.ci(bootKW,type = c("norm", "basic", "perc"))
  
  
  
### posthoc test- bootstrapped Dunn's test ###
  

  Dunn_Bootstrap <- function(d,i) {
    d2 <- d[i,]
    objectDunn<- (dunnTest(d2$dist~d2$Comparison, method= "bh"))
    return(c(p= objectDunn$res$P.unadj, padj =objectDunn$res$P.adj,
             teststat = objectDunn$res$Z))
  }
  
  set.seed(678)
  bootDunn<- boot(data = pairs,statistic = Dunn_Bootstrap, R = 99999, 
                  parallel = "multicore", strata = pairs$Tube1, ncpus=2)
  print(bootDunn) #t1 is P value, t2 is test stat, t3 is df
  str(bootDunn)
  head(bootDunn)
  mean(bootDunn$t[,4]) #bootstrap p.adj for between-extrDuplicate = 0.009
  mean(bootDunn$t[,5]) #bootstrap p.adj for between-SeqDuplicate = 3.046078e-12
  mean(bootDunn$t[,6]) #bootstrap p.adj for ExtrDuplicate-SeqDuplicate = 0.583

  
  #calculate confidence intervals
  bootDunn_CI<-boot.ci(bootDunn,type = c("norm", "basic", "perc"))
  bootDunn_CI$normal
  

  table(pairs$Comparison)

  
#######################################
##### Beta diversity repeatability ####
#######################################
  
  physeq3 #use the unrarefied reads

# Remove one outlier sample "really tiny" in extraction notes= 586 samples

  physeq4 <-subset_samples(physeq3, Sample.ID != "Sample.71.81")
  physeq4

# filter rare taxa from phyloseq object
  
# Compute prevalence of each feature (total number of samples in which a taxon appears at least once), store as data.frame
  prevdf = apply(X = otu_table(physeq4),
               MARGIN = ifelse(taxa_are_rows(physeq4), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts for each phylum to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(physeq4),
                    tax_table(physeq4))
  head(prevdf)
  str(prevdf)

  plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) #average prevelance of features within each phylum and the sum of feature prevalence within each phylum

# Plot the unique phyla: each dot will be a feature- total abundance of that feature across samples on x axis and the prevalance (the fraction of all samples it occurs in on the y axis).
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq4, "Phylum"))
  ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(physeq4),color=Phylum)) +
    # Include a guess for the threshold- here 2% = approx 5 samples, total abundance of 50 reads across all samples
    geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_vline(xintercept = 50, alpha = 0.5, linetype = 2)+  geom_point(size = 1, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + theme(legend.position="none")

# Define prevalence threshold as 2% samples and abundance as 50 total reads across samples
prevalenceThreshold<-2*(586/100)

abundanceThreshold<-50

# Execute the prevalence filter, using `prune_taxa()` function
  head(prevdf1)

  KeepTaxa1<-prevdf1[prevdf1$Prevalence>=prevalenceThreshold,]
  str(KeepTaxa1)
  head(KeepTaxa1)

  KeepTaxa2<- rownames(KeepTaxa1)[(KeepTaxa1$TotalAbundance >= abundanceThreshold)]
  str(KeepTaxa2) #3057 taxa retained after filtering at 2% sample level

  physeq4
  physeqfiltered<- prune_taxa(KeepTaxa2, physeq4)
  physeqfiltered # 3057 taxa, 586 samples

#Filter the phyloseq object to contain just the sequencing and repeat extractions

  TubeNosDups #tube nos of seq dups and repeat extractions
  datBeta<- sample_data(physeqfiltered)
  samplesduplicated <- datBeta[datBeta$TubeNumber %in% TubeNosDups,]
  sampleIDsBetaDup<- samplesduplicated$Sample.ID

  physeqFilteredDups <- prune_samples(sampleIDsBetaDup, physeqfiltered)
  physeqFilteredDups #59 samples duplicated so 118 in total

#CLR transformation of ASV abundances
  library(microbiome)
  physeq_clr <- microbiome::transform(physeqFilteredDups, "clr")

#function to extract an ASV matrix
  vegan_otu <- function(physeq) {
    OTU <- otu_table(physeq)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }

#Extract OTU Matrix and Sample Data  
  clr_v<-vegan_otu(physeq_clr)
  clr_s<-as(sample_data(physeq_clr),"data.frame")  

  library(vegan)
  distMatrixBeta <- vegdist(clr_v, method="euclidean")
  str(distMatrixBeta)
  distMatrixBeta <- as.matrix(distMatrixBeta)
  distMatrixBeta[1:3,1:3]

#extract as a dataframe (just upper triangle- one way comparisons)
  pairwiseDistBeta <- t(combn(colnames(distMatrixBeta), 2))
  head(pairwiseDistBeta)
  pairwiseDistBeta<- data.frame(pairwiseDistBeta, dist=distMatrixBeta[pairwiseDistBeta])
  head(pairwiseDistBeta)
  str(pairwiseDistBeta)
  
#match sample ids to unique tube IDs
  TubeMeta<- data.frame(sample_data(physeq_clr))
  str(TubeMeta)
  TubeMetaNos<- TubeMeta[,c(1,2,5)]
  str(TubeMetaNos)
  
  colnames(pairwiseDistBeta)<- c("Sample.ID","ID2", "dist")
  pairwiseDistBeta2<- merge(pairwiseDistBeta, TubeMetaNos, by="Sample.ID")
  str(pairwiseDistBeta2)
  
  colnames(pairwiseDistBeta2)<- c("ID1","Sample.ID", "dist", "TubeNoID1", "Tube1")
  pairwiseDistBeta3<- merge(pairwiseDistBeta2, TubeMetaNos, by="Sample.ID")
  str(pairwiseDistBeta3)
  colnames(pairwiseDistBeta3)<- c("ID2","ID1", "dist", "TubeNoID1", "Tube1","TubeNoID2","Tube2")
  head(pairwiseDistBeta3)
  
#Copy the two Id columns twice so that we can check for matches between seq or extraction dups:
  pairwiseDistBeta3$ID1Seq<- pairwiseDistBeta3$TubeNoID1
  pairwiseDistBeta3$ID2Seq<- pairwiseDistBeta3$TubeNoID2
  pairwiseDistBeta3$ID1extr<- pairwiseDistBeta3$TubeNoID1
  pairwiseDistBeta3$ID2extr<- pairwiseDistBeta3$TubeNoID2
  str(pairwiseDistBeta3)
  
  #Remove "b" from two of them (ID1seq,ID2seq) to check for tube number matches (the second duplicate is the tube number+b) between seq dup
  pairwiseDistBeta3$ID1Seq<-gsub("b","",as.character(pairwiseDistBeta3$ID1Seq))
  pairwiseDistBeta3$ID2Seq<-gsub("b","",as.character(pairwiseDistBeta3$ID2Seq))
  
  #remove "R" from ID1extr and ID2extr to check for extraction replicates (the second extraction has R after the tube number)
  pairwiseDistBeta3$ID1extr<-gsub("R","",as.character(pairwiseDistBeta3$ID1extr))
  pairwiseDistBeta3$ID2extr<-gsub("R","",as.character(pairwiseDistBeta3$ID2extr))
  
  str(pairwiseDistBeta3)
  
#name type of comparison
  
  pairwiseDistBeta3$SeqDup<-ifelse(pairwiseDistBeta3$ID1Seq==pairwiseDistBeta3$ID2Seq, "SeqDuplicate", "Between")
  pairwiseDistBeta3$extrRep<-ifelse(pairwiseDistBeta3$ID1extr==pairwiseDistBeta3$ID2extr, "ExDuplicate", "Between")
  
  pairwiseDistBeta3$Comparison<- pairwiseDistBeta3$SeqDup
  which(pairwiseDistBeta3$extrRep == "ExDuplicate")
  pairwiseDistBeta3$Comparison[which(pairwiseDistBeta3$extrRep == "ExDuplicate")] <- "ExtrDuplicate"
  

#plot of within and between samples- Figure S2
  pairsbeta<- pairwiseDistBeta3
  head(pairsbeta)
  str(pairsbeta)
  pairsbeta$Comparison<-as.factor(pairsbeta$Comparison)
  ggplot(pairsbeta, aes(Comparison, dist)) +geom_jitter(col="darkgrey", width=0.2, shape=21, size=2)+ 
    geom_boxplot(alpha=0.1, outlier.shape=NA) + stat_summary(fun=mean, geom="point", shape=23, size=4) +
    xlab("\nPairwise comparison") +ylab("Euclidean distance\n") +theme_minimal() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +theme (axis.line = element_line(colour="black"))+
    scale_x_discrete(labels=samplePairs) + theme(element_blank())

####statistical test - bootstrapped Kruskal-Wallis#####
  
  library(boot)
  
  Kruskal_BootstrapBeta <- function(d,i) {
    d2 <- d[i,]
    objectKW<- (kruskal.test(d2$dist~d2$Comparison))
    return(c(p= objectKW$p.value, stat =objectKW$statistic,
             df = objectKW$parameter))
  }
  
  
  set.seed(678)
  bootKW<- boot(data = pairsbeta,statistic = Kruskal_BootstrapBeta, R = 99999, 
                parallel = "multicore", strata = pairsbeta$Tube1, ncpus=2)
  print(bootKW) #t1 is P value, t2 is test stat, t3 is df
  str(bootKW)
  mean(bootKW$t[,1]) #bootstrap estimate of P value = 1.103661e-14
  mean(bootKW$t[,2]) #bootstrap estimate of test statistic (chisq)= 111.2322
  mean(bootKW$t[,3]) #df = 2
  
  
  # output: original is the same as $t0 (original value of statistic in full dataset)
  # bias is difference between the mean of bootstrap realizations (mean from $t- called a bootstrap estimate of T) and value in original dataset (the one from $t0).
  # bootMed is the median of the bootstrapped values
  
  #calculate confidence intervals
  bootKW_CI<-boot.ci(bootKW,type = c("norm", "basic", "perc"))
  
  
### posthoc test- bootstrapped Dunn's test ###
  
  
  Dunn_BootstrapBeta <- function(d,i) {
    d2 <- d[i,]
    objectDunn<- (dunnTest(d2$dist~d2$Comparison, method= "bh"))
    return(c(p= objectDunn$res$P.unadj, padj =objectDunn$res$P.adj,
             teststat = objectDunn$res$Z))
  }
  
  set.seed(678)
  bootDunn<- boot(data = pairsbeta,statistic = Dunn_Bootstrap, R = 99999,
                  parallel = "multicore", strata = pairsbeta$Tube1, ncpus=2)
  print(bootDunn) #t1 is P value, t2 is test stat, t3 is df
  str(bootDunn)
  head(bootDunn)
  mean(bootDunn$t[,4]) #bootstrap p.adj for between-extrDuplicate = 0.00136
  mean(bootDunn$t[,5]) #bootstrap p.adj for between-SeqDuplicate = <0.001
  mean(bootDunn$t[,6]) #bootstrap p.adj for ExtrDuplicate-SeqDuplicate = 0.710

  
