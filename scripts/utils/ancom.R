#Code adapted from: https://github.com/sidhujyatha/ANCOM and modified by Gloria Fackelmann (gloria.fackelmann@uni-ulm.de).

ancom.W = function(	otu_data,
					var_data,
					fformula,
					random_formula,
					var_list,
					var_list_betas,
					multcorr,
					sig){

	n_otu=dim(otu_data)[2]-1
	otu_ids=head(colnames(otu_data), -1)
	data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
	fformula = paste0("lr ~ ", fformula)

	num_errors <- 0
	logratio_mats <- lapply(1:length(var_list), matrix, data= NA, nrow=n_otu, ncol=n_otu)
	betas_stats <- lapply(1:length(var_list), matrix, data=0.0, nrow=n_otu, ncol=8)
	for(ii in 1:n_otu){
	  betas_vals <- lapply(1:length(var_list), matrix, data=NA, nrow=2, ncol=n_otu)
	  for(jj in 1:n_otu){
	    if(ii == jj){
	      for(cur_var in 1:length(var_list)){
	        logratio_mats[[cur_var]][[ii,jj]] <- 1.
	      }
	    }else{
    		data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
    		lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
    		if(sum(abs(lr)) == 0){
    			for(cur_var in 1:length(var_list)){
    				logratio_mats[[cur_var]][ii,jj] <- 1.
    			}
    		}else{
    			lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
    			tryCatch(
            	{
    				model = nlme::lme(
    							formula(fformula),
    							data = lr_dat,
    							random = formula(random_formula),
    							na.action=na.omit,
    							control = lmeControl(sing.tol=1e-20)) 
    				anova_res <- anova(model)
    				betas_model <- summary(model)[[4]][[1]]
    				
    				for(cur_var in 1:length(var_list)){
    				  picker=which(gsub(" ","",row.names(anova_res))== var_list[cur_var])
    				  logratio_mats[[cur_var]][ii,jj] <- anova_res[["p-value"]][picker]
    				  if(anova_res[["p-value"]][picker]<sig){
    				    betas_vals[[cur_var]][1,jj] <- betas_model[[var_list_betas[[cur_var]]]]
    				    betas_vals[[cur_var]][2,jj] <- betas_model[[var_list_betas[[cur_var]]]] + betas_model[["(Intercept)"]]
    				  }
    				}
    
    			},
    			error=function(cond){
    				for(cur_var in 1:length(var_list)){
    					logratio_mats[[cur_var]][ii,jj] <- 1.
    				}
    				num_errors <- num_errors + 1
    			})
    		}
	    }
	  }
	  if(ii %% 1 == 0){
  	  print("Current")
  	  print(ii)
	  }
	  for(cur_var in 1:length(var_list_betas)){
	    betas_stats[[cur_var]][ii,1] <- mean(betas_vals[[cur_var]][1,], na.rm = TRUE)
	    betas_stats[[cur_var]][ii,2] <- sd(betas_vals[[cur_var]][1,], na.rm = TRUE)
	    betas_stats[[cur_var]][ii,3] <- quantile(betas_vals[[cur_var]][1,], na.rm=TRUE)[['25%']]
	    betas_stats[[cur_var]][ii,4] <- quantile(betas_vals[[cur_var]][1,], na.rm=TRUE)[['75%']]
	    
	    betas_stats[[cur_var]][ii,5] <- mean(betas_vals[[cur_var]][2,], na.rm = TRUE)
	    betas_stats[[cur_var]][ii,6] <- sd(betas_vals[[cur_var]][2,], na.rm = TRUE)
	    betas_stats[[cur_var]][ii,7] <- quantile(betas_vals[[cur_var]][2,], na.rm=TRUE)[['25%']]
	    betas_stats[[cur_var]][ii,8] <- quantile(betas_vals[[cur_var]][2,], na.rm=TRUE)[['75%']]
	  }
	}

	print("Number of errors:")
	print(num_errors)
	
	W_total <- list()
	for(cur_var in 1:length(var_list)){
		cur_mat <- logratio_mats[[cur_var]]
		
		ind <- lower.tri(cur_mat)
		cur_mat[ind] <- t(cur_mat)[ind] 

		cur_mat[which(is.finite(cur_mat)==FALSE)] <- 1

		mc.pval <- t(apply(cur_mat,1,function(x){
			s <- p.adjust(x, method = "BH")
			return(s)
		}))

		a <- cur_mat[upper.tri(cur_mat,diag=FALSE)==TRUE]

		b <- matrix(0,ncol=n_otu,nrow=n_otu)

		b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
		diag(b)  <- NA
		ind.1    <- lower.tri(b)
		b[ind.1] <- t(b)[ind.1]

		#########################################
		### Code to extract surrogate p-value
		surr.pval <- apply(mc.pval,1,function(x){
			s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
			return(s0)
		})
		#########################################
		### Conservative
		if(multcorr==1){
			W <- apply(b,1,function(x){
				  subp <- length(which(x<sig))
				})
		### Moderate
		} else if(multcorr==2){
			W <- apply(mc.pval,1,function(x){
				  subp <- length(which(x<sig))
				})
		### No correction
		} else if(multcorr==3){
			W <- apply(cur_mat,1,function(x){
				  subp <- length(which(x<sig))
				})
		}
		
		W_total[[length(W_total)+1]] <- W
		
	}
	
	beta_frame <- data.frame(otu_ids)
	for(cur_var in 1:length(var_list_betas)){
	  column_name <- paste(var_list_betas[[cur_var]], "mean", sep="_")
	  beta_frame[,column_name] <- betas_stats[[cur_var]][,1]
	  column_name <- paste(var_list_betas[[cur_var]], "sd", sep="_")
	  beta_frame[,column_name] <- betas_stats[[cur_var]][,2]
	  column_name <- paste(var_list_betas[[cur_var]], "q25", sep="_")
	  beta_frame[,column_name] <- betas_stats[[cur_var]][,3]
	  column_name <- paste(var_list_betas[[cur_var]], "q75", sep="_")
	  beta_frame[,column_name] <- betas_stats[[cur_var]][,4]
	  
	  column_name <- paste(var_list_betas[[cur_var]], "inter", "mean", sep="_")
	  beta_frame[,column_name] <- betas_stats[[cur_var]][,5]
	  column_name <- paste(var_list_betas[[cur_var]], "inter", "sd", sep="_")
	  beta_frame[,column_name] <- betas_stats[[cur_var]][,6]
	  column_name <- paste(var_list_betas[[cur_var]], "inter", "q25", sep="_")
	  beta_frame[,column_name] <- betas_stats[[cur_var]][,7]
	  column_name <- paste(var_list_betas[[cur_var]], "inter", "q75", sep="_")
	  beta_frame[,column_name] <- betas_stats[[cur_var]][,8]
	  
	}
	return(list(W_total, beta_frame))
  }



ANCOM.main = function(	OTUdat,
						Vardat,
						fformula,
						random_formula,
						var_list,
						var_list_betas,
						multcorr,
						sig,
						prev.cut){
  library("microbiome")
  library("nlme")
	p.zeroes=apply(OTUdat[,-1],2,function(x){
		s=length(which(x==0))/length(x)
		})
	zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
	colnames(zeroes.dist)=c("Taxon","Proportion_zero")
	zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
			  geom_histogram(binwidth=0.1,colour="black",fill="white") + 
			  xlab("Proportion of zeroes") + ylab("Number of taxa") +
			  theme_bw()
	OTUdat.thinned=OTUdat
	OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<=prev.cut))]
	otu.names=head(colnames(OTUdat.thinned), -1)
	res_ancom   <- ancom.W(	OTUdat.thinned,
							Vardat,
							fformula,
							random_formula,
							var_list,
							var_list_betas,
							multcorr,
							sig)
  W.detected <- res_ancom[[1]]
  final_betas <- res_ancom[[2]]
	W_stat_list       <- W.detected
	

	### Bubble plot

	W_frame_list <- list()
	for(cur_var in 1:length(var_list)){
		W_stat <- W_stat_list[[cur_var]]
		W_frame = data.frame(otu.names,W_stat,row.names=NULL)
		W_frame = W_frame[order(-W_frame$W_stat),]
		W_frame$detected_0.95=rep(FALSE,dim(W_frame)[1])
		W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
		W_frame$detected_0.85=rep(FALSE,dim(W_frame)[1])
		W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
		W_frame$detected_0.75=rep(FALSE,dim(W_frame)[1])
		W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
		W_frame$detected_0.65=rep(FALSE,dim(W_frame)[1])
		W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])

		W_frame$detected_0.95[which(W_frame$W_stat>0.95*dim(OTUdat.thinned[,-1])[2])]=TRUE
		W_frame$detected_0.9[which(W_frame$W_stat>0.9*dim(OTUdat.thinned[,-1])[2])]=TRUE
		W_frame$detected_0.85[which(W_frame$W_stat>0.85*dim(OTUdat.thinned[,-1])[2])]=TRUE
		W_frame$detected_0.8[which(W_frame$W_stat>0.8*dim(OTUdat.thinned[,-1])[2])]=TRUE
		W_frame$detected_0.75[which(W_frame$W_stat>0.75*dim(OTUdat.thinned[,-1])[2])]=TRUE
		W_frame$detected_0.7[which(W_frame$W_stat>0.7*dim(OTUdat.thinned[,-1])[2])]=TRUE
		W_frame$detected_0.65[which(W_frame$W_stat>0.65*dim(OTUdat.thinned[,-1])[2])]=TRUE
		W_frame$detected_0.6[which(W_frame$W_stat>0.6*dim(OTUdat.thinned[,-1])[2])]=TRUE

		W_frame_list[[length(W_frame_list)+1]] <- W_frame
	}
	
	final_results=list(W_frame_list, zero.plot, final_betas)
	names(final_results)=c("W.taxa.list","PLot.zeroes", "Betas")
	return(final_results)
}
