#Sarah Worsley

#R code and analysis for: Gut microbiome composition, not alpha diversity, is associated with survival in a natural vertebrate population

#load packages
  library(ggplot2)
  library(phyloseq)
  library(ape)
  library(picante)
  library(MuMIn)
  library(lme4)
  library(lmerTest)
  library(car)
  library(arm)
  library(DHARMa)

#Upload ASV, taxonomy, tree and metadata files to phyloseq:

#ASV table
  asv_table <- read.csv ("ASV_table.csv", row.names=1) #read in asv table with feature names as rownames
  str (asv_table) #should be 56690 features, 636 samples (note 2 samples werefiltered out during QIIME steps as low read number)
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
  physeq2 #610 samples, 55664 taxa


### remove samples with less than 10000 reads (did not reach completion- see file SequencingExtractionRepeatability.R for iNEXT code)
  physeq2
  physeq3<-prune_samples(sample_sums(physeq2)>=10000, physeq2)
  physeq3 # removes 23 samples (some of these are sequencing duplicates)- leaves 587 samples,55664 taxa

  
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

# extract samples sequenced/extracted twice (keep the one with the highest read number)

    physeqRarenoDup<- subset_samples(physeqRare, SequencingDuplicate=="No") #528 remaining
    physeqRarenoDup

#### REMOVE CATCH DUPLICATES ####

# Keeping 1 sample per catch in following order of priority: 1) from tray 2)  from bag, or if neither keep highest read count

    physeqRareFiltered<- subset_samples(physeqRarenoDup, CatchRepeat=="No") #471 samples, 49116 taxa 
    physeqRareFiltered
  
############################################
##### CALCULATE ALPHA DIVERSITY METRICS ####
############################################

# calculate chao1, shannon using estimate_richness()
  set.seed(787)
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

  #write.csv(physeqRareMeta, "FilteredRarefiedSamples_alphaDiversity.csv")

########################################
####### BODY CONDITION ANALYSIS ########
########################################
  
  alphaData <- read.csv("FilteredRarefiedSamples_alphaDiversity.csv") #read in alpha diversity metadata-
  head(alphaData)
  str(alphaData) #471 samples

#remove alpha diversity outlier 
  alphaData <- alphaData[alphaData$Sample.ID != "Sample.71.81",] #outlier: Extraction notes- "really tiny sample". After removal= 470 samples
  
#remove samples where no body mass is recorded

  alphaData_BMass <- alphaData[!is.na(alphaData$BodyMass),] #leaves 448 samples
  str(alphaData_BMass)
  str(unique(alphaData_BMass$BirdID)) #448 samples, 308 birds
  
#remove nestlings
  alphaData_BMNoNestling<-alphaData_BMass[alphaData_BMass$Ageclass != "N",]
  str(alphaData_BMNoNestling) # 436 samples
  str(unique(alphaData_BMNoNestling$BirdID)) #299 individuals

#remove females with eggs as outliers
  alphaData_BMNoNestling<- alphaData_BMNoNestling[alphaData_BMNoNestling$Sample.ID != "Sample.40.46",] #Notes say possible egg present
  alphaData_BMNoNestling<- alphaData_BMNoNestling[alphaData_BMNoNestling$Sample.ID != "Sample.48.54",] #egg present
  str(alphaData_BMNoNestling) # 434 samples
  
#remove floaters with no TQ
  alphaData_BMfiltered<-alphaData_BMNoNestling[alphaData_BMNoNestling$Status != "FLOAT",]
  str(alphaData_BMfiltered) #426 samples
  
#check all have a catch time

  alphaData_BMfiltered<- alphaData_BMfiltered[!is.na(alphaData_BMfiltered$MinutesSinceSunrise),] #leaves 425 samples (one sample does not have a recorded catch time)
  str(alphaData_BMfiltered)
  birdsBM<- as.vector(alphaData_BMfiltered$BirdID)
  str(unique(birdsBM)) #296 individuals
  birds<-data.frame(table(alphaData_BMfiltered$BirdID))
  str(birds[birds$Freq>1,]) #100 birds have > 1 sample = 34%
  
#look at tarsus measurements per observer- check no systematic bias
  
  ggplot(alphaData_BMfiltered, aes(Observer, RightTarsus)) + geom_point()
  ggplot(alphaData_BMfiltered, aes(Observer, BodyMass)) + geom_point()
  summary(alphaData_BMfiltered$RightTarsus)
  
#extract variables from metadata that will be used in the model
  
  ConditionAlphaDiversity<- alphaData_BMfiltered[, c("TerritoryID", "SampleYear" ,"Sample.ID", "BirdID","Ageclass" ,"SexEstimate","TerritoryQuality","FieldPeriodID","SequenceRunNo", "Shannon", "Chao1", "FaithPD", "BodyMass", "RightTarsus", "MinutesSinceSunrise", "Status", "Observer")]
  str(ConditionAlphaDiversity)
  ConditionAlphaDiversity$SexEstimate <- as.factor(ConditionAlphaDiversity$SexEstimate)
  ConditionAlphaDiversity$Ageclass <- as.factor(ConditionAlphaDiversity$Ageclass)
  ConditionAlphaDiversity$BirdID <- as.factor(ConditionAlphaDiversity$BirdID)
  ConditionAlphaDiversity$FieldPeriodID<- as.factor(ConditionAlphaDiversity$FieldPeriodID)
  ConditionAlphaDiversity$SampleYear<- as.factor(ConditionAlphaDiversity$SampleYear)
  ConditionAlphaDiversity$TerritoryID<- as.factor(ConditionAlphaDiversity$TerritoryID)
  ConditionAlphaDiversity$Observer<-as.factor(ConditionAlphaDiversity$Observer)
  ConditionAlphaDiversity$logchao <- log(ConditionAlphaDiversity$Chao1)
  ConditionAlphaDiversity$logfaith <- log(ConditionAlphaDiversity$FaithPD)

  #write.csv(ConditionAlphaDiversity, "BodyConditionAlpha.csv")
  
  


###############################
## LMM for Shannon diversity###
###############################
  
  Vifmodelcheck<- lmer(BodyMass ~ RightTarsus + Shannon + Ageclass + SexEstimate + 
                         TerritoryQuality + FieldPeriodID + MinutesSinceSunrise + 
                         (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  vif(Vifmodelcheck) #all VIFs are under 3, also age and diversity are not correlated (vifs approx 1)
  
  
#Model1- all interactions and squared terms included but removed sequentialy if not significant (enables interpretation of main effects)
  
  model1Shannon<-lmer(BodyMass ~ RightTarsus + Shannon + I(Shannon^2) + Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise +
                        Shannon*SexEstimate + Shannon*Ageclass + (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  
  m1Shannon.std <- standardize(model1Shannon, standardize.y=FALSE) #scale the continuous variables
  summary(m1Shannon.std)
  vcov(m1Shannon.std)
  car::Anova(m1Shannon.std)

  
#drop Ageclass interaction as non-significant interaction term (least significant)
  
  model2Shannon<-lmer(BodyMass ~ RightTarsus + Shannon + I(Shannon^2) + Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise +
                        Shannon*SexEstimate + (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  
  m2Shannon.std <- standardize(model2Shannon, standardize.y=FALSE)
  summary(m2Shannon.std)
  car::Anova(m2Shannon.std)


#drop Sex interaction term as also non-significant interaction term
  
  model3Shannon<- lmer(BodyMass ~ RightTarsus + Shannon + I(Shannon^2) + Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise +
                          (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  
  m3Shannon.std <- standardize(model3Shannon, standardize.y=FALSE)
  summary(m3Shannon.std)
  car::Anova(m3Shannon.std)

  
#drop shannon^2 (not significant) to aid in interpretation of main effect
  
  model4Shannon<-lmer(BodyMass ~ RightTarsus + Shannon + Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise +
                        (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  
  m4Shannon.std <- standardize(model4Shannon, standardize.y=FALSE)
  summary(m4Shannon.std)
  vcov(m4Shannon.std)
  car::Anova(m4Shannon.std)


  library(MuMIn)
  r.squaredGLMM(m4Shannon.std) #R2 calculations
  
#model diagnostics
  
  library(DHARMa)
  
  simulationOutputShannon <- simulateResiduals(fittedModel = m4Shannon.std, plot = F)
  plot(simulationOutputShannon)
  plotResiduals(simulationOutputShannon, form = ConditionAlphaDiversity$Shannon)
  testDispersion(simulationOutputShannon)
  
  
#extract coefficients
  outShannonBM <- summary(m4Shannon.std)
  coeffSBM <- data.frame(outShannonBM$coefficients[,c(1:5)])
  coeffSBM$term<- rownames(coeffSBM)
  str(coeffSBM)
  PValsShannonBM<- coeffSBM[,c(5,6)]
  
  #write.csv(coeffSBM, "Shannon_coefficients.csv")
  
  
#####################  
#### Chao1 model ####
#####################
  

  Vifmodelcheck<- lmer(BodyMass ~ RightTarsus + logchao + Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise + (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  vif(Vifmodelcheck) #all <3


# Model1 chao1 richness (note log transform Chao1 richness as skewed)
  model1Chao<-lmer(BodyMass ~ RightTarsus + logchao + I(logchao^2) + Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise +
                     logchao*SexEstimate + logchao*Ageclass  + (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  summary(model1Chao)
  
  m1chao.std <- standardize(model1Chao, standardize.y=FALSE)
  summary(m1chao.std)
  car::Anova(m1chao.std)
  
  
#drop age*chao interaction as most non-significant term- do this to enable interpretation of main effects
  
  model2Chao<-lmer(BodyMass ~ RightTarsus + logchao + I(logchao^2) + Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise +
                     logchao*SexEstimate + (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  summary(model2Chao)
  
  m2chao.std <- standardize(model2Chao, standardize.y=FALSE)
  summary(m2chao.std)
  car::Anova(m2chao.std)
  

#drop chao*sex interaction as most non-significant term- do this to enable interpretation of main effects
  
  model3Chao<-lmer(BodyMass ~ RightTarsus + logchao + I(logchao^2)+ Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise +
                      (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  
  m3chao.std <- standardize(model3Chao, standardize.y=FALSE)
  summary(m3chao.std)
  car::Anova(m3chao.std)
  
  
#drop chao squared to interpret main effect
  
  model4Chao<- lmer(BodyMass ~ RightTarsus + logchao +  Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise +
                      (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  
  m4chao.std <- standardize(model4Chao, standardize.y=FALSE)
  summary(m4chao.std)
  car::Anova(m4chao.std)
  
  library(MuMIn)
  r.squaredGLMM(m4chao.std) # R2
  
  
#check model diagnostics
  library(DHARMa)
  
  simulationOutputchao <- simulateResiduals(fittedModel = m4chao.std, plot = F)
  plot(simulationOutputchao)
  plotResiduals(simulationOutputchao, form = ConditionAlphaDiversity$logchao)
  testDispersion(simulationOutputchao)


#extract coefficients
  outChaoBM <- summary(m4chao.std)
  coeffCBM<- data.frame(outChaoBM$coefficients[,c(1:5)])
  coeffCBM$term<- rownames(coeffCBM)
  PValsChaoBM<- coeffCBM[,c(5,6)]
  
  
  #write.csv(coeffCBM, "Chao_coefficients.csv")
 

##################### 
#### Faiths model ###
#####################
  
  Vifmodelcheck<- lmer(BodyMass ~ RightTarsus + logfaith + Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise + (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  vif(Vifmodelcheck) #all <3
  
# model1 faithsPD
  
  model1faith<-lmer(BodyMass ~ RightTarsus + logfaith + I(logfaith^2) + Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise +
                      logfaith*SexEstimate + logfaith*Ageclass  + (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  summary(model1faith)
  m1faiths.std <- standardize(model1faith, standardize.y=FALSE)
  summary(m1faiths.std)
  car::Anova(m1faiths.std)

  
#drop age*faith interaction to enable interpretation of main effects
  
  model2faith<-lmer(BodyMass ~ RightTarsus + logfaith + I(logfaith^2) + Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise +
                      logfaith*SexEstimate  + (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  summary(model2faith)
  
  m2faiths.std <- standardize(model2faith, standardize.y=FALSE)
  summary(m2faiths.std)
  car::Anova(m2faiths.std)
  
  
#drop sex*faith interaction
  
  model3faith<-lmer(BodyMass ~ RightTarsus + logfaith + I(logfaith^2) + Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise +
                       (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  summary(model3faith)
 
  m3faiths.std <- standardize(model3faith, standardize.y=FALSE)
  summary(m3faiths.std)
  car::Anova(m3faiths.std)
  

#drop faiths2
  
  model4faith<-lmer(BodyMass ~ RightTarsus + logfaith +  Ageclass + SexEstimate + TerritoryQuality + FieldPeriodID + MinutesSinceSunrise +
                      (1|BirdID) + (1|TerritoryID) + (1|Observer), data= ConditionAlphaDiversity)
  summary(model4faith)
  
  m4faiths.std <- standardize(model4faith, standardize.y=FALSE)
  summary(m4faiths.std)
  car::Anova(m4faiths.std)
  
  library(MuMIn)
  r.squaredGLMM(m4faiths.std)
  
  
#check model diagnostics
  library(DHARMa)
  
  simulationOutputfpd <- simulateResiduals(fittedModel = m4faiths.std, plot = F)
  plot(simulationOutputfpd)
  plotResiduals(simulationOutputfpd, form = ConditionAlphaDiversity$logfaith)
  testDispersion(simulationOutputfpd)
  
  plot(m4faiths.std)
  resid<-residuals(m4chao.std, type= "deviance")
  hist(resid)

#extract coefficients
  
  outFPDBM <- summary(m4faiths.std)
  coeffFBM<- data.frame(outFPDBM$coefficients[,c(1:5)])
  coeffFBM$term <- rownames(coeffFBM)
  PValsFaithBM<- coeffFBM[,c(5,6)]
  
  #write.csv(coeffFBM, "faiths_coefficients.csv")
  summary(m4faiths.std)
  

##############################################  
#### adjust P values for multiple testing ####
##############################################  
  
#Controls for testing the same hypothesis with different diversity metrics 
  
  PValsShannonBM
  PValsChaoBM
  PValsFaithBM
  
  PvaluesAll<- data.frame(PValsShannonBM$Pr...t..,PValsChaoBM$Pr...t.., PValsFaithBM$Pr...t..)
  row.names(PvaluesAll)<- c("Intercept", "RightTarsus", "alphaDiversity", "FL", "OFL", "SA", "Sex",
                            "TerrQuality","FP166", "FP167", "FP171", "FP173", "FP174", "MinutesSinceSunrise")
  
  PvaluesAll<- data.frame(t(PvaluesAll))
  str(PvaluesAll)
  
  lapply(PvaluesAll, \(x) p.adjust (x, method = "BH", n = length(x)))
  

  
  