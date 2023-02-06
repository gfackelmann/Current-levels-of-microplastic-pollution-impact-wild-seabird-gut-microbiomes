#Copyright (c) 2023 Insitute of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany. See the LICENSE file at the top-level directory of this distribution (https://github.com/gfackelmann/Current-levels-of-microplastic-pollution-impact-wild-seabird-gut-microbiomes/blob/main/LICENSE).
#Author: Gloria Fackelmann (gloria.fackelmann@uni-ulm.de; https://www.researchgate.net/profile/Gloria-Fackelmann)

####################################################################
# Load libraries
####################################################################
library("phyloseq") #1.28.0
library("btools") #0.0.1
library("dplyr") #1.0.5
library("nlme") #3.1.141
library("HTSSIP") #1.4.1
library("geiger") #2.0.6.2
library("ape") #5.3
library("deal") #1.2.39
library("data.tree") #1.0.0

##############################################################################################################
#Calculate alpha diversity metrics: Faith's PD, Observed number of ASVs, Shannon index
##############################################################################################################
PD_df <- c()
Richness <- c()
PD_df <- estimate_pd(phylo_obj)
Richness <- estimate_richness(phylo_obj, measures=c("Observed","Shannon"))

sample_data(phylo_obj)$NotRareObserved <- Richness$Observed
sample_data(phylo_obj)$NotRareShannon <- Richness$Shannon
sample_data(phylo_obj)$NotRareChao1 <- Richness$Chao1
sample_data(phylo_obj)$NotRareInvSimpson <- Richness$InvSimpson
sample_data(phylo_obj)$NotRarePD <- PD_df$PD

remove(Richness)
remove(PD_df)

##############################################################################################################
#Calculate alpha diversity metrics: Allen's H
##############################################################################################################
#Code adapted from: Alberdi & Gilbert (2019), "A guide to the application of Hill numbers to DNA-based diversity analyses"
#https://github.com/cylove1112/hilldiv/blob/master/R/hill.div.r

abund <- otu_table(phylo_obj)
tree <- phy_tree(phylo_obj)
samples <- colnames(abund)
sample.vector <- c()
for (s in samples){
  vector <- abund[,s]
  names(vector) <- rownames(abund)
  Li <- tree$edge.length #Get branch lengths
  ltips <- sapply(tree$edge[, 2], function(node) tips(tree, node)) #Sum relative abundances per lineage
  ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector]))) #Sum relative abundances per lineage
  T <- sum(Li * ai) #Get total tree depth
  Li <- Li[ai != 0] #Remove zeros
  ai <- ai[ai != 0] #Remove zeros
  allen <- -sum(Li*ai*log(ai))
  names(allen) <- s
  sample.vector <- c(sample.vector,allen)
}      
return(sample.vector)

sample_data(phylo_obj)$SampleID==names(sample.vector)
sample_data(phylo_obj)$allen <- sample.vector

##############################################################################################################
#Linear mixed models
##############################################################################################################
#Faith's PD
model_PD <- nlme::lme(NotRarePD  ~ 
                  species +
                  scale(propMPcount)*GITlocation +
                  scale(propMPmass)*GITlocation +
                  scale(SequencingDepth),
                random = ~1 | birdID, 
                method = "REML",
                weights = varComb(varIdent(form = ~1 | GITlocation), varExp(form = ~ (scale(SequencingDepth)))),
                na.action=na.fail, data=phylo_obj_df)

#Observed number of ASVs
model_observed <- nlme::lme(sqrt(NotRareObserved)  ~ 
                  species +
                  scale(propMPcount)*GITlocation +
                  scale(propMPmass)*GITlocation +
                  scale(SequencingDepth),
                random = ~1 | birdID, 
                method = "REML",
                weights = varComb(varIdent(form = ~1 | GITlocation), varExp(form = ~ (scale(SequencingDepth)))),
                na.action=na.fail, data=phylo_obj_df)

#Shannon index
model_shannon <- nlme::lme(NotRareShannon  ~ 
                  species +
                  scale(propMPcount)*GITlocation +
                  scale(propMPmass)*GITlocation +
                  scale(SequencingDepth),
                random = ~1 | birdID, 
                method = "REML",
                weights = varComb(varIdent(form = ~1 | GITlocation), varExp(form = ~ (scale(SequencingDepth)))),
                na.action=na.fail, data=phylo_obj_df)

#Allen's H
model_allen <- nlme::lme(sqrt(allen)  ~ 
                         species +
                         scale(propMPcount)*GITlocation +
                         scale(propMPmass)*GITlocation +
                         scale(SequencingDepth),
                       random = ~1 | birdID, 
                       method = "REML",
                       weights = varComb(varIdent(form = ~1 | GITlocation), varExp(form = ~ (scale(SequencingDepth)))),
                       na.action=na.fail, data=phylo_obj_df)