#Copyright (c) 2023 Insitute of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany. See the LICENSE file at the top-level directory of this distribution (https://github.com/gfackelmann/Current-levels-of-microplastic-pollution-impact-wild-seabird-gut-microbiomes/blob/main/LICENSE).
#Author: Gloria Fackelmann (gloria.fackelmann@uni-ulm.de; https://www.researchgate.net/profile/Gloria-Fackelmann)

####################################################################
# Load libraries
####################################################################
library("phyloseq") #1.28.0
library("vegan") #2.5.5

####################################################################
# Calculate the various distance matrices
####################################################################
#weighted unifrac
DistW <- phyloseq::distance(phylo_obj, method="wUniFrac")
#unweighted unifrac
DistUW <- phyloseq::distance(phylo_obj, method="uunifrac")

####################################################################
# Permanova
####################################################################
propMPcount <- as(sample_data(phylo_obj), "data.frame")$propMPcount
propMPmass <- as(sample_data(phylo_obj), "data.frame")$propMPmass
species <- as(sample_data(phylo_obj), "data.frame")$species
GITlocation <- as(sample_data(phylo_obj), "data.frame")$GITlocation
sex <- as(sample_data(phylo_obj), "data.frame")$sex
SequencingDepth <- as(sample_data(phylo_obj), "data.frame")$SequencingDepth
birdID <- as(sample_data(phylo_obj), "data.frame")$birdID


perm_W <- adonis(DistW ~
                         species +
                         GITlocation +
                         scale(SequencingDepth) + 
                         sex +
                         scale(propMPmass)*GITlocation +
                         scale(propMPcount)*GITlocation +
                         scale(propMPmass)*species +
                         scale(propMPcount)*species,
                       strata = birdID,
                       data = as(sample_data(phylo_obj), "data.frame"),
                       permutations=9999)

perm_UW <- adonis(DistUW ~
                         species +
                         GITlocation +
                         scale(SequencingDepth) + 
                         sex +
                         scale(propMPmass)*GITlocation +
                         scale(propMPcount)*GITlocation +
                         scale(propMPmass)*species +
                         scale(propMPcount)*species,
                       strata = birdID,
                       data = as(sample_data(phylo_obj), "data.frame"),
                       permutations=9999)