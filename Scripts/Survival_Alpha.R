# Sarah Worsley

# R code and analysis for: Gut microbiome composition, not alpha diversity, is associated with survival in a natural vertebrate population

#load packages
  library(ggplot2)
  library(phyloseq)
  library(ape)
  library(dplyr)
  library(picante)
  library(plyr)
  library(RColorBrewer)
  library(forcats)
  library(patchwork)
  library(GGally)
  library(car)
  library(arm)
  library(multcomp)
  library(DHARMa)
  library(picante)


#Upload files to phyloseq:

#ASV table
  asv_table <- read.csv ("ASV_table.csv", row.names=1) #read in asv table with feature names as rownames
  str (asv_table) #should be 56690 features, 636 samples outputed from QIIME analysis
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
  physeq #check there are the right numbers of samples/features/variables (636/56690/91)

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
  

### removal of samples from controls####
# remove control samples (positive/negative controls, collection controls)  
  physeq2 <- subset_samples(physeq, SampleType == "F") 
  physeq2 #610 samples


################################################  
### remove samples with less than 10000 reads###
################################################  
  
  #Samples with fewer reads did not reach completion- see SequencingExtractionRepeatability.R for iNext code.
  
  physeq2
  physeq3<-prune_samples(sample_sums(physeq2)>=10000, physeq2)
  physeq3 # removes 23 samples (some of these are sequencing duplicates)- leaves 587 samples
  
  #check the number of reads per sample
  sample_reads<-data.frame(reads=sample_sums(physeq3))
  head(sample_reads)
  str(sample_reads)
  outlier<- c("Sample.71.81") #outlier sample removed from downstream analysis
  sample_reads2<-sample_reads[!(row.names(sample_reads) %in% outlier),]
  str(sample_reads2) #586 samples
  summary(sample_reads2) #reads range from 10,979 to 744,600 per sample (excluding the outlier sample)
  
  #check ASV number per sample- calc mean and SD
  ASVno<-data.frame(estimate_richness(physeq3, split=TRUE, measures= c("Observed"))) #observed ASVs per sample
  head(ASVno)
  ASVno2<-ASVno[!(row.names(ASVno) %in% outlier),]
  str(ASVno2) #586
  summary(ASVno2) #mean ASV no= 383
  sd(ASVno2) #sd 277
  
###############################################  
##### RAREFY READS TO MIN SAMPLING DEPTH ######
###############################################
  
#rarefy to 10000 and set seed before rarefying (28367), so that results are reproducible  
  physeqRare<-rarefy_even_depth(physeq3, 10000, rngseed = 28367)

# 6548 ASVs removed after subsampling- leaves 49116 across 587 samples
  physeqRare
  sample_sums(physeqRare)
  
  
######################################
#### REMOVE SEQUENCING DUPLICATES ####
######################################
  
  #See SequencingExtractionRepeatability.R for assessment of repeatability
  
  #### extract samples sequenced/extracted twice (keep the one with the highest read number)###
  
  physeqRareNoDup<- subset_samples(physeqRare, SequencingDuplicate=="No") #528 remaining
  physeqRareNoDup

#################################    
#### Remove catch duplicates ####
################################# 
  
  # Keeping 1 sample per catch in following order of priority: 1) from tray 2)  from bag, or if neither keep highest read count
  
  physeqRareFiltered<- subset_samples(physeqRareNoDup, CatchRepeat=="No")  
  physeqRareFiltered #471 samples, 49116 taxa  

############################################
##### CALCULATE ALPHA DIVERSITY METRICS ####
############################################ 
  
# calculate chao1, shannon using estimate_richness()
  richnessEstRare<-estimate_richness(physeqRareFiltered, split=TRUE, measures= c("Chao1", "Shannon", "Observed"))
  head(richnessEstRare)
  str(richnessEstRare)
  
#install.packages("picante") for calculating Faith's PD
  library(picante)
  faiths_asv <- as.data.frame(physeqRareFiltered@otu_table)
  faiths_tree <- physeqRareFiltered@phy_tree
  faithPD <- pd(t(faiths_asv),faiths_tree, include.root = T) # t transposes the otu table for use in picante
  head(faithPD)
  richnessEstRare$FaithPD <- faithPD$PD 
  head(richnessEstRare)

#add alpha diversity metrics to metadata
  physeqRareMeta <- as.data.frame(sample_data(physeqRareFiltered))
  head(physeqRareMeta)
  physeqRareMeta$Chao1 <- richnessEstRare$Chao1
  physeqRareMeta$Shannon <- richnessEstRare$Shannon
  physeqRareMeta$FaithPD <- richnessEstRare$FaithPD
  physeqRareMeta$Observed <- richnessEstRare$Observed
  head(physeqRareMeta)
  str(physeqRareMeta) #471 samples, 321 individuals (106 birds have multiple samples)
  str(unique(physeqRareMeta$BirdID))
  
  #write.csv(physeqRareMeta, "FilteredRarefiedSamples_alphaDiversity.csv")

#########################################
##### PLOT RELATIVE ABUNDANCES ##########
#########################################

  
#remove outlier sample (Sample.71.81)- very tiny sample with extremely low alpha diversity
  
  physeqFilteredoutlier<- subset_samples( physeqRareFiltered, Sample.ID != "Sample.71.81")
  physeqFilteredoutlier #470 samples, 49116 taxa  
  
  dataFiltered<- sample_data(physeqFilteredoutlier)
  str(dataFiltered)
  str(unique(dataFiltered$BirdID)) #how many unique birds are there = 320 individuals (not including outlier)
  
#check ASV number per sample in rarefied samples
  ASVnofiltered<-data.frame(estimate_richness(physeqFilteredoutlier, split=TRUE, measures= c("Observed")))
  head(ASVnofiltered)
  str(ASVnofiltered)
  summary(ASVnofiltered) #mean ASV no= 368
  sd(ASVnofiltered$Observed) #sd 253.112
  
#transform counts to relative abundances
  physeqFilteredoutlier
  relabund <- transform_sample_counts(physeqFilteredoutlier, function(x){(x / sum(x))*100})
  relabundTable<-otu_table(relabund)
  head(relabundTable)
  colSums(relabundTable) 
  
#agglomerate to phylum level
  ps.phylum<- tax_glom(relabund, taxrank="Phylum", NArm=FALSE)
  ps.phylum #detected 42 taxa
  plot_bar(ps.phylum, fill="Phylum")
  
  dat <- psmelt(ps.phylum) #create a dataframe of agglomerated data at phylum level.
  str(dat)
  medians <- ddply(dat, ~Phylum, function(x) c(median=median(x$Abundance))) #find the median count per phylum
  
#Find mean abundance of proteobacteria
  prot<-dat[dat$Phylum=="Proteobacteria",]
  mean(prot$Abundance)
  sd(prot$Abundance)
  
#Find mean abundance of firmicutes
  Firm<-dat[dat$Phylum=="Firmicutes",]
  mean(Firm$Abundance)
  sd(Firm$Abundance)
  
#Find mean abundance of actinobacteria
  Actino<-dat[dat$Phylum=="Actinobacteria",]
  mean(Actino$Abundance)
  sd(Actino$Abundance)

# find Phyla whose median relabundance is less than 1%
  other <- medians[medians$median <= 0.01,]$Phylum
  str(other)
# change their name to "Other" to make plot less complex
  dat[dat$Phylum %in% other,]$Phylum <- 'Other'
  
  
#barplot:relative abundance of each phylum per sample.
  #order the bars by the abundance of proteobacteria
  
  unique(dat$Phylum)
  dat$Phylum<-factor(dat$Phylum, c("Proteobacteria", "Firmicutes", "Actinobacteria", "Chloroflexi", "Planctomycetes", "Bacteroidetes", "Cyanobacteria", "Acidobacteria", "Patescibacteria",  "Verrucomicrobia","Other"))
  dfprot<-dat[dat$Phylum=="Proteobacteria",]
  str(dfprot)
  dfprot<-dfprot[order(-dfprot$Abundance),]
  sampleOrder<-dfprot$Sample
  head(sampleOrder)
  sampleOrder[1:37]
  
  dat$Sample<-factor(dat$Sample, sampleOrder)
  str(dat)
  colourCount <- length(unique(dat$Phylum))
  pal<- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FFFF99", "#FB9A99","#B15928", "#E31A1C", "#FDBF6F","#FF7F00", "#CAB2D6", "#6A3D9A")
  getPalette <- colorRampPalette(pal)
  
  p <- ggplot(data=dat, aes(x=Sample, y=Abundance, fill=Phylum))
  p + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount))  + guides(fill=guide_legend(nrow=5)) + ylab("Relative abundance (%)")
  
  
#barplot: phyla rel abundance per samples, per ageclass- Figure 1
  dat$Ageclass<-as.factor(dat$Ageclass)
  dat$Ageclass<-factor(dat$Ageclass, c("N", "FL", "OFL", "SA", "A"))
  levels(dat$Ageclass)
  levels(dat$Ageclass) <- c("Nestlings", "Fledglings", "Old Fledlings", "Sub-Adults", "Adults")
  p <- ggplot(data=dat, aes(x=Sample, y=Abundance, fill=Phylum))
  p + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount))  + guides(fill=guide_legend(nrow=11)) + facet_wrap(~Ageclass, scales="free_x", nrow=1)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))+ ylab("Relative Abundance (%)")+
    theme(axis.title = element_text(size=18), axis.text = element_text(size=16), legend.text = element_text(size=14), legend.title = element_text(size=16), strip.text = element_text(size=16))
  
  
#number of samples per ageclass
  dataAge<- sample_data(physeqFilteredoutlier)
  table(dataAge$Ageclass)
  
  
###################################################
####### SURVIVAL ANALYSIS- alpha diversity ########
###################################################
  
  alphaData <- read.csv("FilteredRarefiedSamples_alphaDiversity.csv") #read in alpha diversity metadata
  head(alphaData)
  str(alphaData) #471 samples
  
#remove the 72 samples from 2020 as no survival data for these birds, leaves 399 samples
  
  alphaData_no2020 <- alphaData[!is.na(alphaData$SurvivedNextSeason),]
  str(alphaData_no2020)
  str(unique(alphaData_no2020$BirdID)) #399 samples, 279 birds when 2020 samples removed.
  head(alphaData_no2020$SampleDate) 
  
#Take the latest sample per bird as not enough with multiple samples to control for this
  
  alphaDataLatestSample<- alphaData_no2020[alphaData_no2020$LatestSample=="yes",]
  str(alphaDataLatestSample) # 279 samples/birds
  
#check survival numbers
  str(alphaDataLatestSample[alphaDataLatestSample$SurvivedNextSeason==1,]) #235 survived to next season
  str(alphaDataLatestSample[alphaDataLatestSample$SurvivedNextSeason==0,]) #44 did not survive to next season
  write.csv(alphaDataLatestSample,"SurvivalAlphaData_latestBird.csv")
  
  
#extract variable from metadata that will be used in the model
  
  SurvivalAlphaDiversity<- read.csv("SurvivalAlphaData_latestBird.csv")
  str(SurvivalAlphaDiversity)
  SurvivalData <- SurvivalAlphaDiversity[, c("Sample.ID", "SampleYear", "BirdID","Ageclass" ,"SexEstimate", "SurvivedNextSeason", "SurvivedNextSeasonFactor", "TerritoryQualityRaw","TerritoryQuality", "FieldPeriodID","SequenceRunNo", "Shannon", "Chao1", "FaithPD", "Status")]
  
  
#correct the data format
  
  SurvivalData$SampleYear <- as.factor(SurvivalData$SampleYear)
  SurvivalData$Ageclass <- as.factor(SurvivalData$Ageclass)
  SurvivalData$SexEstimate <- as.factor(SurvivalData$SexEstimate)
  SurvivalData$FieldPeriodID <- as.factor(SurvivalData$FieldPeriodID)
  SurvivalData$SequenceRunNo <- as.factor(SurvivalData$SequenceRunNo)
  SurvivalData$SurvivedNextSeasonFactor <- as.factor(SurvivalData$SurvivedNextSeasonFactor)
  SurvivalData$logChao <- log(SurvivalData$Chao1)
  SurvivalData$logFaith <- log(SurvivalData$FaithPD)  
  str(SurvivalData)
  
#remove nestlings
  SurvivalData2<-SurvivalData[SurvivalData$Ageclass != "N",]
  str(SurvivalData2) #270 samples
  str(SurvivalData2[SurvivalData2$SurvivedNextSeason==0,]) #41 did not survive to next season (3 removed)
  
#remove floaters as no TQ values assigned
  SurvivalData2<-SurvivalData2[SurvivalData2$Status != "FLOAT",]
  str(SurvivalData2) #265 samples
  str(SurvivalData2[SurvivalData2$SurvivedNextSeason==0,]) #38 did not survive to next season (3 removed)
  
#remove outlier with low alpha diversity
  SurvivalData2<-SurvivalData2[SurvivalData2$Sample.ID != "Sample.71.81",] #outlier: Extraction notes- "really tiny sample". After removal= 264 birds, 38 didn't survive
  
#final numbers died vs survived
  table(SurvivalData2$SurvivedNextSeasonFactor)
  
#plot survival to the next season 
  
  S_Survival<- ggplot(SurvivalData2, aes(SurvivedNextSeasonFactor, Shannon, fill=SurvivedNextSeasonFactor))+ geom_boxplot(outlier.shape = NA, fill=c("white","white")) +  geom_point(pch = 21, position = position_jitterdodge(), size=2) +theme(axis.text = element_text(size=14), axis.title = element_text(size=17), legend.position = "none") +
    xlab("\nSurvived to the following season") +ylab("Shannon diversity\n")
  
  C_Survival<- ggplot(SurvivalData2, aes(SurvivedNextSeasonFactor, logChao, fill=SurvivedNextSeasonFactor))+ geom_boxplot(outlier.shape = NA, fill=c("white","white")) +  geom_point(pch = 21, position = position_jitterdodge(), size=2) +theme(axis.text = element_text(size=14), axis.title = element_text(size=17), legend.position = "none") +
    xlab("\nSurvived to the following season") +ylab("Log Chao1 richness\n")
  
  F_Survival<- ggplot(SurvivalData2, aes(SurvivedNextSeasonFactor, logFaith, fill=SurvivedNextSeasonFactor))+ geom_boxplot(outlier.shape = NA, fill=c("white","white")) +  geom_point(pch = 21, position = position_jitterdodge(), size=2) +theme(axis.text = element_text(size=14), axis.title = element_text(size=17), legend.position = "none") + 
    xlab("\nSurvived to the following season") +ylab("Log Faith's PD\n")
  
  library(patchwork)
  (S_Survival| C_Survival| F_Survival)
  
### Plot proportions surviving per age class - Figure S6
  
  SurvivalByAge<- SurvivalData2[,c(4,7)]
  str(SurvivalByAge)
  
  SurvPercent<- data.frame(prop.table(table(SurvivalByAge$Ageclass, SurvivalByAge$SurvivedNextSeasonFactor),
                                      margin = 1))
  colnames(SurvPercent)<- c("Ageclass", "Survived","Percent")
  SurvPercent<-SurvPercent[!SurvPercent$Ageclass == "N",]
  SurvPercent<-SurvPercent[order(SurvPercent$Ageclass),]

  count<- data.frame(SurvivalByAge %>%
    group_by(Ageclass, SurvivedNextSeasonFactor) %>%
    tally())
  count
  
  SurvPercent$Counts<- count$n
  
  SurvPercent$Ageclass <- factor(SurvPercent$Ageclass,c("FL", "OFL", "SA", "A"))
  
  SurvivalProportions<- ggplot(SurvPercent,aes(Ageclass, Percent, fill=Survived,label=Counts)) +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=c("#E69F00","#0072B2"))+
    geom_text(size = 5, position = position_stack(vjust = 0.5))+
    xlab("\nAge class")+ ylab("Proportion survived\n")+ theme_bw()+
    theme(axis.text = element_text(size=18), axis.title = element_text(size=18), 
          legend.text = element_text(size=18), legend.title=element_text(size=18))+
    theme(panel.grid = element_blank())
    
  SurvivalProportions + guides(fill=guide_legend(title="Survived"))
  
#######################  
#####SHANNON MODEL##### 
#######################
  
#check for correlations
  pairsplot<-SurvivalData2[, c("Ageclass", "TerritoryQuality", "Shannon", "SexEstimate", "SurvivedNextSeason", "SampleYear", "FieldPeriodID",  "SequenceRunNo")]
  ggpairs(pairsplot)

#Check the VIFs
  
  ShannonVIF<- glm(SurvivedNextSeason~ Shannon +  Ageclass + SexEstimate + TerritoryQuality + SampleYear, data= SurvivalData2, family = binomial("logit"))
  plot(ShannonVIF)
  vif(ShannonVIF)

# All VIFs are below 2
# Note that you can't have sequence run number and the FieldPeriodID in the same model as they are highly correlated.

#Run the full model

  m1Shannon<- glm(SurvivedNextSeason~ Shannon + Ageclass + SexEstimate + TerritoryQuality + SampleYear +
                  I(Shannon^2) + Shannon*Ageclass + Shannon*SexEstimate + Shannon*TerritoryQuality + Shannon*SampleYear, data= SurvivalData2, family = binomial("logit"))
  summary(m1Shannon) 

  M1shannon.std <- standardize(m1Shannon, standardize.y=FALSE) #standardise continuous variables
  summary(M1shannon.std)
  car::Anova(M1shannon.std)


#drop Ageclass interaction as not significant (enables interpretation of main effect)

  m2Shannon<-  glm(SurvivedNextSeason~ Shannon + Ageclass + SexEstimate + TerritoryQuality + SampleYear +
                     I(Shannon^2) +  Shannon*SexEstimate + Shannon*TerritoryQuality + Shannon*SampleYear, data= SurvivalData2, family = binomial("logit"))
  
  summary(m2Shannon)  

  M2shannon.std <- standardize(m2Shannon, standardize.y=FALSE)
  summary(M2shannon.std)
  car::Anova(M2shannon.std)


#drop Sex interaction term

  m3Shannon<-glm(SurvivedNextSeason~ Shannon + Ageclass + SexEstimate + TerritoryQuality + SampleYear +
                   I(Shannon^2) +  Shannon*TerritoryQuality + Shannon*SampleYear, data= SurvivalData2, family = binomial("logit"))
  
  summary(m3Shannon)  

  M3shannon.std <- standardize(m3Shannon, standardize.y=FALSE)
  summary(M3shannon.std)
  car::Anova(M3shannon.std)


#drop Sample Year interaction term

  m4Shannon<- glm(SurvivedNextSeason~ Shannon + Ageclass + SexEstimate + TerritoryQuality + SampleYear +
                    I(Shannon^2) +  Shannon*TerritoryQuality, data= SurvivalData2, family = binomial("logit"))
  
  summary(m4Shannon)  

  M4shannon.std <- standardize(m4Shannon, standardize.y=FALSE)
  summary(M4shannon.std)

  car::Anova(M4shannon.std)

#drop TQ interaction

  m5Shannon<-  glm(SurvivedNextSeason~ Shannon + Ageclass + SexEstimate + TerritoryQuality + SampleYear +
                     I(Shannon^2) , data= SurvivalData2, family = binomial("logit"))
  
  summary(m5Shannon)  

  M5shannon.std <- standardize(m5Shannon, standardize.y=FALSE)
  summary(M5shannon.std)
  car::Anova(M5shannon.std)
  

#drop shannon squared term to enable interpretation of main effects

  m6Shannon<- glm(SurvivedNextSeason~ Shannon + Ageclass + SexEstimate + TerritoryQuality + SampleYear, data= SurvivalData2, family = binomial("logit"))
  
  summary(m6Shannon)  
  
  M6shannon.std <- standardize(m6Shannon, standardize.y=FALSE)
  summary(M6shannon.std)
  car::Anova(M6shannon.std)

  summary(glht(M6shannon.std, mcp(Ageclass="Tukey"))) ## SA and OFL significantly different


#check model diagnostics using DHARMa
  
  simulationOutputShannon <- simulateResiduals(fittedModel = M6shannon.std, plot = F)
  plot(simulationOutputShannon)
  plotResiduals(simulationOutputShannon, form = SurvivalData2$Shannon)
  testDispersion(simulationOutputShannon)
  
  
#extract coefficients
  outShannonSurv <- summary(M6shannon.std)
  coeffSSurv<- data.frame(outShannonSurv$coefficients[,c(1:4)])
  coeffSSurv$term <- rownames(coeffSSurv)
  PValsShannonSurv<- coeffSSurv[,c(4,5)]
  #write.csv(coeffSSurv, "Shannon_coefficients.csv")
  
  
#####################
#####Chao1 MODEL#####  
#####################
  
#Check the VIFs
  ChaoVIF<- glm(SurvivedNextSeason~ logChao + Ageclass + SexEstimate + TerritoryQuality +SampleYear, data= SurvivalData2, family = binomial("logit"))
  summary(ChaoVIF)  

  #library(car)
  vif(ChaoVIF) #all below 2
  

#Run the full model

  m1Chao<- glm(SurvivedNextSeason~ logChao + Ageclass + SexEstimate + TerritoryQuality + SampleYear+
                  I(logChao^2) + logChao*Ageclass + logChao*SexEstimate +logChao*TerritoryQuality  + logChao*SampleYear, data= SurvivalData2, family = binomial("logit"))
  summary(m1Chao)  

  M1chao.std <- standardize(m1Chao, standardize.y=FALSE)
  summary(M1chao.std)  
  car::Anova(M1chao.std)


#drop TQ interaction- not significant so remove to interpret main effects

  m2Chao<- glm(SurvivedNextSeason~ logChao + Ageclass + SexEstimate + TerritoryQuality + SampleYear+
                 I(logChao^2) + logChao*Ageclass + logChao*SexEstimate  + logChao*SampleYear, data= SurvivalData2, family = binomial("logit"))
  
  summary(m2Chao)  

  M2chao.std <- standardize(m2Chao, standardize.y=FALSE)
  summary(M2chao.std)  

  car::Anova(M2chao.std)

#drop age interaction

  m3Chao<- glm(SurvivedNextSeason~ logChao + Ageclass + SexEstimate + TerritoryQuality + SampleYear+
                 I(logChao^2) + logChao*SexEstimate  + logChao*SampleYear, data= SurvivalData2, family = binomial("logit"))
  
  summary(m3Chao)  

  M3chao.std <- standardize(m3Chao, standardize.y=FALSE)
  summary(M3chao.std)  

  car::Anova(M3chao.std)


#drop year interaction

  m4Chao<-glm(SurvivedNextSeason~ logChao + Ageclass + SexEstimate + TerritoryQuality + SampleYear+
              I(logChao^2) + logChao*SexEstimate, data= SurvivalData2, family = binomial("logit"))

  summary(m4Chao)  

  M4chao.std <- standardize(m4Chao, standardize.y=FALSE)
  summary(M4chao.std)  

  car::Anova(M4chao.std)

#drop sex interaction
  m5Chao<- glm(SurvivedNextSeason~ logChao + Ageclass + SexEstimate + TerritoryQuality + SampleYear+
                 I(logChao^2), data= SurvivalData2, family = binomial("logit"))
  summary(m5Chao)  

  M5chao.std <- standardize(m5Chao, standardize.y=FALSE)
  summary(M5chao.std)  

  car::Anova(M5chao.std)


#drop quadratic to interpret linear effect

  m6Chao<- glm(SurvivedNextSeason~ logChao + Ageclass + SexEstimate + TerritoryQuality + SampleYear, data= SurvivalData2, family = binomial("logit"))
  summary(m6Chao)  

  M6chao.std <- standardize(m6Chao, standardize.y=FALSE)
  summary(M6chao.std)  

  car::Anova(M6chao.std)

  
  library(multcomp)
  summary(glht(M6chao.std, mcp(Ageclass="Tukey"))) ## SA and OFL significantly different
  
#check model diagnostics using dharma
  library(DHARMa)
  
  simulationOutputChao <- simulateResiduals(fittedModel = M6chao.std, plot = F)
  plot(simulationOutputChao)
  plotResiduals(simulationOutputChao, form = SurvivalData2$logChao)
  testDispersion(simulationOutputChao)


#extract coefficients
  outChaoSurv <- summary(M6chao.std)
  coeffCSurv<- data.frame(outChaoSurv$coefficients[,c(1:4)])
  coeffCSurv$term <- rownames(coeffCSurv)
  PValsChaoSurv<- coeffCSurv[,c(4,5)]
  
  #write.csv(coeffCSurv, "Chao_coefficients.csv")

######################
#####Faiths MODEL#####
######################

#Check the VIFs
  faithVIF<- glm(SurvivedNextSeason~ logFaith + Ageclass + SexEstimate + TerritoryQuality + SampleYear, data= SurvivalData2, family = binomial("logit"))
  summary(faithVIF)  

#library(car)
  vif(faithVIF) #all below 3

#Run the full model

  m1fpd<- glm(SurvivedNextSeason~ logFaith + Ageclass + SexEstimate + TerritoryQuality + SampleYear +
               I(logFaith^2) + logFaith*Ageclass + logFaith*SexEstimate +logFaith*TerritoryQuality + logFaith*SampleYear, data= SurvivalData2, family = binomial("logit"))
  summary(m1fpd)  
  
  std.FPD <-standardize(m1fpd, standardize.y=FALSE)
  summary(std.FPD)
  car::Anova(std.FPD)
  
#drop TQ interaction term

  m2fpd<-  glm(SurvivedNextSeason~ logFaith + Ageclass + SexEstimate + TerritoryQuality + SampleYear +
               I(logFaith^2) + logFaith*Ageclass + logFaith*SexEstimate + logFaith*SampleYear, data= SurvivalData2, family = binomial("logit"))

  summary(m2fpd)  

  std.FPD2 <-standardize(m2fpd, standardize.y=FALSE)
  summary(std.FPD2)

  car::Anova(std.FPD2)

#drop age interaction

  m3fpd<- glm(SurvivedNextSeason~ logFaith + Ageclass + SexEstimate + TerritoryQuality + SampleYear +
                I(logFaith^2) + logFaith*SexEstimate + logFaith*SampleYear, data= SurvivalData2, family = binomial("logit"))
  
  summary(m3fpd)  

  std.FPD3 <-standardize(m3fpd, standardize.y=FALSE)
  summary(std.FPD3)

  car::Anova(std.FPD3)


#drop year interaction term

  m4fpd<- glm(SurvivedNextSeason~ logFaith + Ageclass + SexEstimate + TerritoryQuality + SampleYear +
                 I(logFaith^2) + logFaith*SexEstimate, data= SurvivalData2, family = binomial("logit"))
  summary(m4fpd)  

  std.FPD4 <-standardize(m4fpd, standardize.y=FALSE)
  summary(std.FPD4)

  car::Anova(std.FPD4)

#drop Sex interaction
  
  m5fpd<-glm(SurvivedNextSeason~ logFaith + Ageclass + SexEstimate + TerritoryQuality + SampleYear +
               I(logFaith^2), data= SurvivalData2, family = binomial("logit"))
  summary(m5fpd)  

  std.FPD5 <-standardize(m5fpd, standardize.y=FALSE)
  summary(std.FPD5)

  car::Anova(std.FPD5)


#drop quadratic

  m6fpd<-glm(SurvivedNextSeason~ logFaith + Ageclass + SexEstimate + TerritoryQuality + SampleYear, data= SurvivalData2, family = binomial("logit"))
  summary(m6fpd)  

  std.FPD6 <-standardize(m6fpd, standardize.y=FALSE)
  summary(std.FPD6)
  car::Anova(std.FPD6)

  library(multcomp)
  summary(glht(std.FPD6, mcp(Ageclass="Tukey"))) ## SA and OFL significantly different

#check model diagnostics using dharma
  library(DHARMa)
  
  simulationOutputfpd <- simulateResiduals(fittedModel = std.FPD6, plot = F)
  plot(simulationOutputfpd)
  plotResiduals(simulationOutputfpd, form = SurvivalData2$logFaith)
  testDispersion(simulationOutputfpd)
  
#extract coefficients
  outFaithsSurv <- summary(std.FPD6)
  coeffFSurv<- data.frame(outFaithsSurv$coefficients[,c(1:4)])
  coeffFSurv$term <- rownames(coeffFSurv)
  PValsFaithSurv<- coeffFSurv[,c(4,5)]
  #write.csv(coeffFSurv, "Faiths_coefficients.csv")
  
  
#### adjust P values for multiple testing- control for test of same hypothesis using three different alpha diversity metrics ####
  
  PValsShannonSurv
  PValsChaoSurv
  PValsFaithSurv
  
  PvaluesAllSurvival<- data.frame(PValsShannonSurv$Pr...z..,PValsChaoSurv$Pr...z.., PValsFaithSurv$Pr...z..)
  row.names(PvaluesAllSurvival)<- c("Intercept","alphaDiversity", "FL", "OFL", "SA", "Sex",
                            "TerrQuality","yr2018", "yr2019")
  
  PvaluesAllSurvival<- data.frame(t(PvaluesAllSurvival))
  str(PvaluesAllSurvival)
  
  lapply(PvaluesAllSurvival, \(x) p.adjust (x, method = "BH", n = length(x)))
  



