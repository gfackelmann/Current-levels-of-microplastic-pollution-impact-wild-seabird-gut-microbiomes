#Copyright (c) 2023 Insitute of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany. See the LICENSE file at the top-level directory of this distribution (https://github.com/gfackelmann/Current-levels-of-microplastic-pollution-impact-wild-seabird-gut-microbiomes/blob/main/LICENSE).
#Author: Gloria Fackelmann (gloria.fackelmann@uni-ulm.de; https://www.researchgate.net/profile/Gloria-Fackelmann)

####################################################################
# Load libraries
####################################################################
library("microbiome") #1.9.16
library("nlme") #3.1.141
source("utils/ancom.R")

####################################################################
# prep for ANCOm function
####################################################################
1formula <- "scale(propMPcount) + scale(propMPmass) + scale(propMPcount)*GITlocation + scale(propMPmass)*GITlocation + scale(propMPcount)*species + scale(propMPmass)*species + scale(SequencingDepth) + sex"

random_formula <- "~1 | birdID"

var_list <- c("scale(propMPcount)", "scale(propMPmass)", "scale(propMPcount):GITlocation", "scale(propMPmass):GITlocation", "scale(propMPcount):species", "scale(propMPmass):species")

var_list_betas <- c("scale(propMPcount)", "scale(propMPmass)", "scale(propMPcount):GITlocationCloacalSwab", "scale(propMPmass):GITlocationCloacalSwab", "scale(propMPcount):speciesNorthernFulmar", "scale(propMPmass):speciesNorthernFulmar")

####################################################################
# run ANCOM
####################################################################

phylo_obj_df_t <- as.data.frame(t(otu_table(phylo_obj)))
colnames(phylo_obj_df_t) <- paste("X", colnames(phylo_obj_df_t), sep="")
colnames(phylo_obj_df_t) <- gsub('-','.', colnames(phylo_obj_df_t))
colnames(phylo_obj_df_t) <- gsub('[(]','', colnames(phylo_obj_df_t))
colnames(phylo_obj_df_t) <- gsub('[)]','', colnames(phylo_obj_df_t))
phylo_obj_df_t$Sample.ID <- rownames(phylo_obj_df_t)

system.time({ancom_res <- ANCOM.main(
  OTUdat=phylo_obj_df_t, 
  Vardat=meta_df,
  fformula=1formula,
  random_formula=random_formula,
  var_list = var_list,
  var_list_betas = var_list_betas,
  multcorr=2, 
  sig=0.05,
  prev.cut=1)})


