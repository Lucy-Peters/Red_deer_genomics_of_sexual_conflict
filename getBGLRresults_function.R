#function to extract individual and SNP based results from BGLR model

getBGLRresults<-function(model, posterior, geno_mat){
  
  #get credible intervals 95%
  SNP_CIs_data<-getCIS(posterior)
  pos<-length(model$ETA)
  
  #test if credible interval contains 0
  SNP_CIs_data$significant<-ifelse(SNP_CIs_data$lower.bound< 0 & SNP_CIs_data$upper.bound > 0, "no", "yes")
  SNP_CIs_data$bhat<-model$ETA[[pos]]$b
  SNP_CIs_data$SD.bhat<-model$ETA[[pos]]$SD.b
  SNP_CIs_data$SNP.name<-model$ETA[[pos]]$colNames
  SNP_CIs_data$var.b<-model$ETA[[pos]]$varB
  SNP_CIs_data$SD.var.b<-model$ETA[[pos]]$SD.varB
  
  SNP_CIs_data<-SNP_CIs_data[, c("SNP.name", "bhat", "SD.bhat","lower.bound", "upper.bound","var.b", "SD.var.b",  "significant")]
  
  #response variable estimates
  ID_PhenoGeno_data<-as.data.frame(model$probs)
  ID_PhenoGeno_data$yhat<-model$yHat
  ID_PhenoGeno_data$SD.yhat<-model$SD.yHat
  ID_PhenoGeno_data$probs.yhat<-pnorm(model$yHat)
  ID_PhenoGeno_data$y<-model$y
  
  #get estimated genomic values (for phenotype)
  gHat<-as.vector(geno_mat%*%model$ETA[[pos]]$b)
  
  #add genomic value to phenotypic estimates
  ID_PhenoGeno_data$ghat<-gHat
  ID_PhenoGeno_data$probs.ghat<-pnorm(ID_PhenoGeno_data$ghat)
  colnames(ID_PhenoGeno_data)[c(1:2)]<-c("probability.failure", "probability.survival")
  
  ID_PhenoGeno_data<-ID_PhenoGeno_data[, c("y", "yhat", "SD.yhat", "probs.yhat", "probability.survival", "ghat", "probs.ghat")]
  
  BGLR_results_list<-list(SNP_CIs_data, ID_PhenoGeno_data)
  
  BGLR_results_list
  
}