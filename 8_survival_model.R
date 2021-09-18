#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Red deer: yearling survival model - sex-specific allele dropout; Bayesian GWAS  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# library(RODBC)
# library(plyr)
# library(dplyr)
# library(ggplot2)
# library(BGLR)
# library(data.table)
# library(pROC)
# 
# db<-"C:\\Users\\s1767711\\Documents\\Red_deer\\database\\current\\RedDeer1.93.accdb" 
# con<-odbcConnectAccess2007(db)
# queries<-sqlTables(con, tableType = "VIEW")
# by_dy_data<-sqlFetch(con, "Birth_death_year_data")
# odbcClose(con)
# 
# head(by_dy_data)
# write.table(by_dy_data, file="data/Deer_survival_data-raw.txt", sep="\t", col.names = T, row.names = F, quote=F)
# 
# #create variable for lifespan (if dead)
# by_dy_data<-mutate(by_dy_data, lifespan=DeathYear-BirthYear)
# colnames(by_dy_data)[c(1:4)]<-c("ID", "sex", "birthyear", "deathyear")
# by_dy_data<-subset(by_dy_data, !is.na(birthyear))
# 
# #exclude calves born in 2019 - data doesn't include 2020 data, don't know if calves alive in 2019
# #made it to 2020
# by_dy_data<-subset(by_dy_data, !birthyear==2019)
# 
# #exclude deer with unknown sex
# table(by_dy_data$sex)
# by_dy_data<-subset(by_dy_data, !sex==3)
# 
# by_dy_data<-by_dy_data %>%
#   mutate(sex = ifelse(sex==1, "female", "male"))
# 
# #by_dy_data<-subset(by_dy_data, birthyear>=1980)
# 
# #data set only including animals with high quality data (non-fringe, study area)
# by_dy_data<-subset(by_dy_data, DataQuality == 1)
# 
# #create variable that represents lifespan or current age if still alive
# max(by_dy_data$birthyear, na.rm = T)
# by_dy_data$lifespan_comp<-ifelse(is.na(by_dy_data$lifespan)==TRUE, (2019-by_dy_data$birthyear), by_dy_data$lifespan)
# 
# #separate date last seen into month and year and add as variables to data frame
# last.seen<-stringr::str_split(by_dy_data$LastSeen, "-")
# 
# last.year<-character()
# last.month<-character()
# 
# for (n in c(1:nrow(by_dy_data))) {
#   year<-last.seen[[n]][1]
#   last.year<-append(last.year, year)
#   
#   month<-last.seen[[n]][2]
#   last.month<-append(last.month, month)
# }
# 
# by_dy_data<-by_dy_data %>%
#   mutate(last.seen.year=last.year, last.seen.month=last.month)
# 
# #if death year or month is unobserved use year or month animal was last seen
# by_dy_data$deathyear.comp<-ifelse(is.na(by_dy_data$deathyear)==TRUE & by_dy_data$Status=="D", 
#                                   by_dy_data$last.seen.year, by_dy_data$deathyear)
# 
# by_dy_data$deathmonth.comp<-ifelse(is.na(by_dy_data$DeathMonth)==TRUE & by_dy_data$Status=="D", 
#                                   by_dy_data$last.seen.month, by_dy_data$DeathMonth)
# 
# #add yearling survival variable
# by_dy_data<-by_dy_data %>%
#   mutate(yearling_survival= ifelse(lifespan_comp == 0, 0, 1))
# 
# by_dy_data<-by_dy_data %>%
#   mutate(yearling_survival= ifelse(yearling_survival == 1 & lifespan_comp < 2 & deathmonth.comp <= 5, 0, yearling_survival))
# 
# #exclude calves that died before first winter 
# by_dy_data<-subset(by_dy_data, !(yearling_survival == 0 & deathmonth.comp >= 5 & deathmonth.comp < 10 ))
# 
# #make variable indicating if calf has survived to 2 years
# by_dy_data<-by_dy_data %>%
#   mutate(secondyear_survival= ifelse(lifespan_comp <= 1, 0, 1))
# 
# by_dy_data<-by_dy_data %>%
#   mutate(secondyear_survival= ifelse(secondyear_survival == 1 & lifespan_comp < 3 & Status == "D" & deathmonth.comp <= 5, 
#                                      0, secondyear_survival))
# 
# #make separate data frames for yearling and second year survival 
# by_dy_data_yearling<-by_dy_data
# by_dy_data_secondyear<-by_dy_data
# 
# #exclude calves younger than or 2 years that were shot during shooting season from secondyear table
# by_dy_data_secondyear<-subset(by_dy_data_secondyear, !(DeathType == "S" & lifespan_comp < 3 & deathmonth.comp < 3))
# 
# #exclude calves that died in first year due to being shot
# by_dy_data_yearling<-subset(by_dy_data_yearling, !(DeathType=="S" & yearling_survival == 0))
# 
# #how many males and females are in data set?
# ggplot(by_dy_data_yearling, aes(sex))+geom_histogram(fill="grey", stat = "count")
# ggplot(by_dy_data_secondyear, aes(sex))+geom_histogram(fill="grey", stat = "count")
# 
# #how many males and female calves died in first winter?
# survival_count1<-by_dy_data_yearling %>%
#   group_by(sex, yearling_survival) %>%
#   summarise(count = n())
# 
# #how many male and female died after 2 years?
# survival_count2<-by_dy_data_secondyear %>%
#   group_by(sex, secondyear_survival) %>%
#   summarise(count = n())
# 
# ggplot(by_dy_data_yearling, aes(as.factor(yearling_survival), fill=sex))+
#   geom_histogram(position = "dodge", bins = 2, stat= "count")
# 
# ggplot(by_dy_data_secondyear, aes(as.factor(secondyear_survival), fill=sex))+
#   geom_histogram(position = "dodge", bins = 2, stat= "count")
# 
# 
# #read in and collate imputed deer genotype files (AlphaImpute)

# #make large tables with genotype and marker information of all chromosomes
# 
# geno_all_chr<-data.frame()
# impute_qc_all<-data.frame()
# 
# for (chr in c(1:2,4:13, 15:22, 24:32)) {
#   
#   #AlphaImpute imputed genotype files
#   assign("geno_chr_in", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaImpute_run/Results/ImputeGenotypes_chr",
#                                            chr, ".txt"), header=F))
#   
#   assign("impute_qc_chr", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaImpute_run/Results/ImputationQuality_chr",
#                                           chr, ".txt"), header=F))
#   
#   
#   if(chr==1){
#     
#     geno_all_chr<-geno_chr_in
#     
#   } else {
#     
#     geno_all_chr<-cbind(geno_all_chr, geno_chr_in[, -1])
#     
#   }
#   
#   impute_qc_all<-rbind(impute_qc_all, impute_qc_chr)
#   
# }
# 
# 
# 
# #count how many chromosomes per individual have imputation quality of at least 0.95
# impute_qc_summary<-impute_qc_all %>%
#   group_by(ID) %>%
#   summarise(count = sum(V7>=0.95))
# 
# table(impute_qc_summary$count)
# 
# #only keep individuals where all chromosomes (29) have imputation quality of at least 0.95
# impute_qc_ID_ok<-impute_qc_summary %>%
#   filter(count==29)
# 
# linkmap<-read.table("data/Cervus_elaphus_linkage_map_data.ordered.txt", header=T)
# linkmap_sub<-subset(linkmap, !(CEL.LG == 3 | CEL.LG == 14 | CEL.LG == 23| CEL.LG == 33 | CEL.LG == 34))

# 
# #read in imputed genotype, imputation quality and marker map file
# 
# #AlphaImpute imputed genotype files
# assign("geno_all_chr", read.table("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaImpute_run/Results/ImputeGenotypes.txt", header=F))
# 
# assign("impute_qc_all", read.table("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaImpute_run/Results/ImputationQuality.txt", header=F))
# 
# assign("marker_map_all", read.table("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaImpute_run/Input_files/Deer_AlphaImpute_marker_file_all.txt", header=T))
# 
# #exclude IDs not in phenotypic data set
# colnames(geno_all_chr)[1]<-"ID"
# geno_all_chr<-subset(geno_all_chr, ID %in% by_dy_data_secondyear$ID)
# 
# #only include individuals with high imputation quality overall
# colnames(impute_qc_all)[1]<-"ID"
# impute_qc_all<-subset(impute_qc_all, ID %in% geno_all_chr$ID)
# impute_qc_ID_ok<-subset(impute_qc_all, V7>=0.95)
# 
# geno_all_chr<-subset(geno_all_chr, ID %in% impute_qc_ID_ok$ID)
# 
# #make sure IDs in phenotype data are the same and is ordered same as genotype data
# by_dy_in<-subset(by_dy_data_secondyear, ID %in% geno_all_chr$ID)
# geno_all_chr$ID.order<-c(1:nrow(geno_all_chr))
# by_dy_in<-join(by_dy_in, geno_all_chr[, c("ID", "ID.order")], by="ID")
# by_dy_in<-arrange(by_dy_in, ID.order)
# 
# #make ID names row names in genotype matrix and add SNP names as column names
# geno_all_chr$ID.order<-NULL
# rownames(geno_all_chr)<-geno_all_chr$ID
# geno_all_chr$ID<-NULL
# geno_all_chr<-as.matrix(geno_all_chr)
# 
# colnames(geno_all_chr)<-marker_map_all$SNP.Name
# head(geno_all_chr)[, c(1:10)]
# 
# #subset of genotype matrix
# geno_all_chr<-geno_all_chr[, c(1:1000)]
# geno_all_chr<-as.matrix(geno_all_chr)
# 
# #add second winter as variable to phenotype df
# by_dy_in<-by_dy_in %>%
#   mutate(second_winter=birthyear+1)
# 
# 
# write.table(by_dy_in[, -8], file="results/Deer_secondyear_survival_data.txt", sep = "\t", col.names = T, row.names = F, quote = F)
# save(geno_all_chr, file = "results/Deer_imputed_genotype_matrix.RData")

#model testing run on ashworth server

library(plyr)
library(dplyr)
library(BGLR)
library(pROC)


args <- commandArgs(trailingOnly = TRUE)
resultpath <- as.character(args[1])

#read in phenotype data
by_dy_in<-read.table("results/Deer_secondyear_survival_data.txt", header=T, stringsAsFactors = F)
yBin<-by_dy_in$secondyear_survival

#load genotype data
load("results/Deer_imputed_genotype_matrix.RData")

#basic BayesB model

#model and prior specifications
ETA_BB<-list(list(~factor(sex), data = by_dy_in, model = "FIXED"),
          list(X = geno_all_chr, model = "BayesB", saveEffects=TRUE))

model_BB<-BGLR(y=yBin, response_type='ordinal', ETA=ETA_BB,
          nIter=120000, burnIn=20000, thin = 5, saveAt = paste0(resultpath, "BGLR_BayesB_"))

save(model_BB, file=paste0(resultpath, "BGLR_model_BayesB.RData"))

#get estimated genomic values (for phenotype)
gHat_BB<-as.vector(pnorm(geno_all_chr%*%model_BB$ETA[[2]]$b))

roc_BB<-roc(yBin, gHat_BB, auc.polygon=TRUE, grid=TRUE, plot = TRUE )

# #basic BayesA model
# 
# #model and prior specifications
# ETA_BA<-list(list(~factor(sex), data = by_dy_in, model = "FIXED"),
#              list(X = geno_all_chr, model = "BayesA", saveEffects=TRUE))
# 
# model_BA<-BGLR(y=yBin, response_type='ordinal', ETA=ETA_BA,
#                nIter=120000, burnIn=20000, thin = 5, saveAt = paste0(resultpath, "BGLR_BayesA_"))
# 
# save(model_BA, file=paste0(resultpath, "BGLR_model_BayesA.RData"))
#  
# #get estimated genomic values (for phenotype)
# gHat_BA<-as.vector(pnorm(geno_all_chr%*%model_BA$ETA[[2]]$b))
# 
# roc_BA<-roc(yBin, gHat_BA, auc.polygon=TRUE, grid=TRUE, plot = TRUE )
# 
# #test which model is better at predicting phenotype from genotype
# roc.test.geno<-roc.test(roc_BB, roc_BA, alternative = "greater")
# 
# save(roc.test.geno, file = paste0(resultpath, "ROC_test_BayesB_vs_BayesA.RData"))
# 
# #decide on genetic model based on roc.test result
# if(roc.test.geno$p.value<0.05){
#   model.geno<-"BayesB"
# } else {
#   model.geno<-"BayesA"
# }
# 

model.geno<-"BayesB"

#run model with selected genetic model and add environmental effect (birthyear)

ETA_by<-list(list(~factor(sex), data = by_dy_in, model = "FIXED"),
             list(~factor(birthyear), data = by_dy_in, model = "BRR"),
             list(X = geno_all_chr, model = model.geno, saveEffects=TRUE))

model_by<- BGLR(y=yBin, response_type='ordinal', ETA=ETA_by,
                nIter=120000, burnIn=20000, thin = 5, saveAt = paste0(resultpath, "BGLR_by_"))

save(model_by, file = paste0(resultpath, "BGLR_model_birthyear.RData"))

#get estimated genomic values (for phenotype)
gHat_by<-as.vector(pnorm(geno_all_chr%*%model_by$ETA[[3]]$b))

roc_by<-roc(yBin, gHat_by, auc.polygon=TRUE, grid=TRUE, plot = TRUE )

#test which model is better at predicting phenotype from genotype
if(model.geno == "BayesB"){
  roc.test.by<-roc.test(roc_by, roc_BB)
} else {
  roc.test.by<-roc.test(roc_by, roc_BA, alternative = "greater")
}

save(roc.test.by, file = paste0(resultpath, "ROC_test_birthyear_vs_genetic.RData"))
save(roc.test.by, file =  "results/ROC_test_birthyear_vs_genetic.RData")

#make new ETA based on roc.test result
if(roc.test.by$p.value<0.05){

  roc_plus<- roc_by
  model_plus<-model_by

  ETA_sy<- list(list(~factor(sex), data = by_dy_in, model = "FIXED"),
                list(~factor(birthyear), data = by_dy_in, model = "BRR"),
                list(~factor(second_winter), data = by_dy_in, model = "BRR"),
                list(X = geno_all_chr, model = model.geno, saveEffects=TRUE))

} else if ( (roc.test.by$p.value>0.05) && (model.geno == "BayesB") ) {

  roc_plus<- roc_BB
  model_plus<-model_BB

  ETA_sy<- list(list(~factor(sex), data = by_dy_in, model = "FIXED"),
                list(~factor(second_winter), data = by_dy_in, model = "BRR"),
                list(X = geno_all_chr, model = model.geno, saveEffects=TRUE))

} else if ( (roc.test.by$p.value>0.05) && (model.geno == "BayesA") ) {

  roc_plus<- roc_BA
  model_plus<-model_BA

  ETA_sy<- list(list(~factor(sex), data = by_dy_in, model = "FIXED"),
                list(~factor(second_winter), data = by_dy_in, model = "BRR"),
                list(X = geno_all_chr, model = model.geno, saveEffects=TRUE))

}


model_sy<- BGLR(y=yBin, response_type='ordinal', ETA=ETA_sy,
                nIter=120000, burnIn=20000, thin = 5, saveAt = paste0(resultpath, "BGLR_sy_"))

save(model_sy, file = paste0(resultpath, "BGLR_model_secondyear.RData"))

pos<-length(model_sy$ETA)

#get estimated genomic values (for phenotype)
gHat_sy<-as.vector(pnorm(geno_all_chr%*%model_sy$ETA[[pos]]$b))

roc_sy<-roc(yBin, gHat_sy, auc.polygon=TRUE, grid=TRUE, plot = TRUE )

roc.test.sy<-roc.test(roc_sy, roc_plus, alternative = "greater")

save(roc.test.sy, file = paste0(resultpath, "ROC_test_secondyear_addition.RData"))

print(paste0("Final model is ..."))

if(roc.test.sy$p.value<0.05){

  print(paste0(str(model_sy)))

} else {

  print(paste0(str(model_plus)))

}


#test if there is a difference in genomic prediction accuracy between birth year and second year model
load("results/ROC_test_birthyear_vs_genetic.RData")
load("results/ROC_test_secondyear_addition.RData")

roc_by<-roc.test.by$roc1
roc_sy<-roc.test.sy$roc1

#test maternal (Mum ID) effect
#from previous models we know BayesB with no additional random effects is best model

#read in pedigree
ped<-read.table("results/Pedigree_deer_2017_19_merged.txt", header =T)
head(ped)

#join Mum ID and phenotype data
by_dy_in<-join(by_dy_in, ped[, c("ID", "Mother")], by="ID")

#exclude entries with unknown mothers
by_dy_in<-subset(by_dy_in, !is.na(Mother))
yBin<-by_dy_in$secondyear_survival

#subset genotype matrix to IDs in phenotype data
geno_all_chr<-subset(geno_all_chr, rownames(geno_all_chr) %in% by_dy_in$ID)

#Mum effect model

ETA_mum<-list(list(~factor(sex), data = by_dy_in, model = "FIXED"),
             list(~factor(Mother), data = by_dy_in, model = "BRR"),
             list(X = geno_all_chr, model = "BayesA", saveEffects=TRUE))

model_mum<- BGLR(y=yBin, response_type='ordinal', ETA=ETA_mum,
                nIter=120000, burnIn=20000, thin = 5, saveAt = paste0(resultpath, "BGLR_mum_"))

save(model_mum, file = paste0(resultpath, "BGLR_model_Mum.RData"))

#get Mum model roc object
gHat_mum<-as.vector(pnorm(geno_all_chr%*%model_by$ETA[[3]]$b))
roc_mum<-roc(yBin, gHat_mum, auc.polygon=TRUE, grid=TRUE, plot = TRUE )

#get roc from previously established best model (BayesA, sex and gene effects only)
load("results/ROC_test_BayesB_vs_BayesA.RData")
roc_plus<-roc.test.geno[["roc2"]]

roc.test.mum<-roc.test(roc_mum, roc_plus, alternative = "greater")
save(roc.test.mum, file = paste0(resultpath, "ROC_test_Mum_addition.RData"))

print(paste0("Final model is ..."))

if(roc.test.mum$p.value<0.05){

  print(paste0("BayesA model, sex, Mum ID and gene effects "))

} else {

  print(paste0("BayesA model, sex and gene effects only"))

}


#
#full BayesB model
#
##model and prior specifications
# ETA_BB<-list(list(~factor(sex), data = by_dy_in, model = "FIXED"),
#              list(~factor(birthyear), data = by_dy_in, model = "BRR"),
#              list(~factor(second_winter), data = by_dy_in, model = "BRR"),
#              list(~factor(Mother), data = by_dy_in, model = "BRR"),
#              list(X = geno_all_chr, model = "BayesB", saveEffects=TRUE))
#
# model_BB<-BGLR(y=yBin, response_type='ordinal', ETA=ETA_BB,
#          nIter=120000, burnIn=20000, thin = 5, saveAt = paste0(resultpath, "BGLR_BayesB_MumEnv_"))
#
#save(model_BB, file=paste0(resultpath, "BGLR_model_BayesB_MumEnv.RData"))
#
#
##run sex specific models
# 
#phenotype data
by_dy_fem<-subset(by_dy_in, !sex == "male")
yBin_fem<-by_dy_fem$secondyear_survival

by_dy_male<-subset(by_dy_in, !sex == "female")
yBin_male<-by_dy_male$secondyear_survival

#genotype data
geno_all_chr_fem<-subset(geno_all_chr, rownames(geno_all_chr) %in% by_dy_fem$ID)
geno_all_chr_male<-subset(geno_all_chr, rownames(geno_all_chr) %in% by_dy_male$ID)

# #model structure
# ETA_fem<-list(list(~factor(birthyear), data = by_dy_fem, model = "BRR"),
#               list(~factor(second_winter), data = by_dy_fem, model = "BRR"),
#               list(~factor(Mother), data = by_dy_fem, model = "BRR"),
#               list(X = geno_all_chr_fem, model = "BayesB", saveEffects=TRUE))
#
#
# ETA_male<-list(list(~factor(birthyear), data = by_dy_male, model = "BRR"),
#               list(~factor(second_winter), data = by_dy_male, model = "BRR"),
#               list(~factor(Mother), data = by_dy_male, model = "BRR"),
#               list(X = geno_all_chr_male, model = "BayesB", saveEffects=TRUE))
# 
##sex specific models
#model_fem<-BGLR(y=yBin_fem, response_type='ordinal', ETA=ETA_fem,
#                 nIter=120000, burnIn=20000, thin = 5, saveAt = paste0(resultpath, "BGLR_fem_MumEnv_"))
# 
#save(model_fem, file = paste0(resultpath, "BGLR_model_female_MumEnv.RData"))
#
# model_male<-BGLR(y=yBin_male, response_type='ordinal', ETA=ETA_male,
#                  nIter=120000, burnIn=20000, thin = 5, saveAt = paste0(resultpath, "BGLR_male_MumEnv_"))
# 
# save(model_male, file = paste0(resultpath, "BGLR_model_male_MumEnv.RData"))
# 

# 
# #heritability model
# 
# #Computing the genomic relationship matrix
# geno_all_chr_scaled<-scale(geno_all_chr,center=TRUE,scale=TRUE)
# G<-tcrossprod(geno_all_chr_scaled)/ncol(geno_all_chr_scaled)
# #Computing the eigen-value decomposition of G
# EVD <-eigen(G)
# #3# Setting the linear predictor
# ETA_h2<-list(list(~factor(sex), data = by_dy_in, model = "FIXED"),
#              list(V=EVD$vectors,d=EVD$values, model='RKHS'))
# 
# #h2 model
# model_h2<-BGLR(y=yBin, response_type='ordinal', ETA=ETA_h2,
#                nIter=120000, burnIn=20000, thin = 5, saveAt = paste0(resultpath, "BGLR_h2_"))
# 
# save(model_h2, file = paste0(resultpath, "BGLR_model_h2.RData"))
# 
# 
#Evaluate best model (general, not sex specific)
load("results/BGLR_model_BayesB_MumEnv.RData")

plot(model_BA$ETA[[2]]$SD.b ~ model_BA$ETA[[2]]$b,col=2,
     main='Vulcano Plot',
     xlab='Estimated Effect',ylab='Est. Posterior SD')

bHat<- model_BA$ETA[[2]]$b
SD.bHat<- model_BA$ETA[[2]]$SD.b

plot(bHat^2, ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')

plot(abs(0.5-pnorm(bHat)), ylab='Estimated binomial Effect',
     type='o',cex=.5,col=4,main='Marker Effects')


#get marker effects and SD on data scale (probit link) - evaluate at standard normal cumulative distribution function
bHat_binom<-pnorm(bHat)
SD.bHat_binom<-pnorm(SD.bHat)

plot(SD.bHat_binom~bHat_binom)

#filter out 'significant SNPs
#calculate probability that effect is different from intercept 0

load("results/BGLR_model_BayesB_MumEnv.RData")

#get posterior distributions
BayesB.post<-readBinMat("results/BGLR_BayesB_MumEnv_ETA_5_b.bin")

#hist(BayesB.post[, 1], breaks = 100)
#mean(BayesB.post[, 1])

#test
mat<-BayesB.post[c(1:nrow(BayesB.post)), c(1:4)]

getCIS<-function(mat){

  mat.ordered<-as.data.frame(apply(mat, 2, sort))

  lower.conf<-t(mat.ordered %>%
    filter(row_number()==(nrow(mat)*0.025)))

  upper.conf<-t(mat.ordered %>%
    filter(row_number()==(nrow(mat)*0.975)))

  CI_df<-data.frame(lower.bound = lower.conf[,1], upper.bound = upper.conf[,1])

  CI_df

}


source("scripts/getBGLRresults.R")


#sex effect
bhat.sex<-model_BB$ETA[[1]]$b
bhat.sex.post<-read.table("results/BGLR_BayesB_MumEnv_ETA_1_b.dat", header=T)
bhat.sexCIS<-getCIS(bhat.sex.post)

bhat.sexCIS$bhat<-bhat.sex
bhat.sexCIS.data<-as.data.frame(sapply(bhat.sexCIS, pnorm))
bhat.sexCIS.effect<-as.data.frame(bhat.sexCIS.data$`sapply(bhat.sexCIS, pnorm)` - 0.5)
  
#female BGLR model
#load model and posterior SNP distributions
load("results/BGLR_model_female_MumEnv.RData")
female.post<-readBinMat("results/BGLR_fem_MumEnv_ETA_4_b.bin")

#phenotype data
by_dy_fem<-subset(by_dy_in, !sex == "male")

#genotype data
geno_all_chr_fem<-subset(geno_all_chr, rownames(geno_all_chr) %in% by_dy_fem$ID)

BGLR_fem_results<-getBGLRresults(model_fem, female.post, geno_all_chr_fem)


#male BGLR model
#load model and posterior SNP distributions
load("results/BGLR_model_male_MumEnv.RData")
male.post<-readBinMat("results/BGLR_male_MumEnv_ETA_4_b.bin")

#phenotype data
by_dy_male<-subset(by_dy_in, !sex == "female")

#genotype data
geno_all_chr_male<-subset(geno_all_chr, rownames(geno_all_chr) %in% by_dy_male$ID)

BGLR_male_results<-getBGLRresults(model_male, male.post, geno_all_chr_male)

#load result compiled on server
load("results/BGLR_male_results.RData")
load("results/BGLR_fem_results.RData")
load("results/BGLR_BayesB_results.RData")

BGLR_male_SNPs<-BGLR_male_results[[1]]
table(BGLR_male_SNPs$significant)
BGLR_male_SNPs$model<-"male_specific"

BGLR_fem_SNPs<-BGLR_fem_results[[1]]
table(BGLR_fem_SNPs$significant)
BGLR_fem_SNPs$model<-"female_specific"

BGLR_BayesB_SNPs<-BGLR_BayesB_results[[1]]
table(BGLR_BayesB_SNPs$significant)
BGLR_BayesB_SNPs$model<-"both_sexes"

BGLR_SNP_effects<-rbind(BGLR_BayesB_SNPs, BGLR_fem_SNPs)
BGLR_SNP_effects<-rbind(BGLR_SNP_probs, BGLR_male_SNPs)
BGLR_SNP_effects<-BGLR_SNP_effects %>%
  select(-"significant")

write.table(BGLR_SNP_effects, file="results/BGLR_model_SNP_effects_all_models.txt", sep="\t", col.names = T,
            row.names = F, quote = F)

#phenotypic probabilities
BGLR_fem_phen<-BGLR_fem_results[[2]]
BGLR_fem_phen$ID<-model_fem$ETA[[1]]$data$ID
BGLR_fem_phen$model<-"female_specific"

BGLR_BayesB_phen<-BGLR_BayesB_results[[2]]
BGLR_BayesB_phen$ID<-model_BB$ETA[[1]]$data$ID
BGLR_BayesB_phen$model<-"both_sexes"

BGLR_male_phen<-BGLR_male_results[[2]]
BGLR_male_phen$ID<-model_male$ETA[[1]]$data$ID
BGLR_male_phen$model<-"male_specific"

BGLR_phen_probs<-rbind(BGLR_BayesB_phen, BGLR_fem_phen)
BGLR_phen_probs<-rbind(BGLR_phen_probs, BGLR_male_phen)

BGLR_phen_probs<-BGLR_phen_probs[, c("ID", "y","yhat", "probs.yhat", "ghat", "probs.ghat", "model")]

write.table(BGLR_phen_probs, file="results/BGLR_model_probabilities_phenotype_all_models.txt", sep="\t", col.names = T,
            row.names = F, quote = F)


#plot SNP effects from male and female BGLR models
linkmap_cow_pos<-read.table("data/Cervus_elaphus_linkage_map_additional_BT_pos.txt", header=T)
head(linkmap_cow_pos)
names(linkmap_cow_pos)[names(linkmap_cow_pos)=="SNP.Name"]<-"SNP.name"

#df for SNP effects
SNP_effects<-data.frame(SNP.name = c(model_fem$ETA[[4]]$colNames, model_male$ETA[[4]]$colNames), 
                        effect_size = c(model_fem$ETA[[4]]$b, model_male$ETA[[4]]$b),
                        sex = c( rep(c("female", "male"), each = length(model_fem$ETA[[4]]$colNames))))

SNP_effects<-join(SNP_effects, linkmap_cow_pos[, c("SNP.name","CEL.LG")])
names(SNP_effects)[names(SNP_effects)=="CEL.LG"]<-"Chromosome"

SNP_effects$effect_size<-abs(0.5-pnorm(SNP_effects$effect_size))
linkmap_sub<-linkmap_cow_pos[, c("SNP.name", "CEL.order")]

source("scripts/ggplot_manhattan.R")

scales_y <- list(female = scale_y_continuous(limits = c(0, 0.004)),
                 male = scale_y_continuous(limits = c(0, 0.0015)))

library(RColorBrewer)
palette<-brewer.pal(2, "Set2")
palette<-palette[c(1:2)]

p_effect_size<-manhattan_plot(df_plot = SNP_effects, link_order = linkmap_sub, y.var = "effect_size",
                              chr_no = 33, chr_colours= palette, xlab_text = "CEL linkage group")

ggsave("plots/BGLR_SNP_effect_sizes_split_sex.png", p_effect_size, width = 25, height = 25, units = "cm")

#heritability of survival
load("results/BGLR_model_h2_MumEnv.RData")
varA<-model_h2$ETA[[5]]$varU

#get posterior distributions for all random effects
Bayesh2.post<-read.table("results/BGLR_h2_MumEnv_ETA_5_varU.dat", header=F)
Bayesh2Mum.post<-read.table("results/BGLR_h2_MumEnv_ETA_4_varB.dat", header=F)
Bayesh2sy.post<-read.table("results/BGLR_h2_MumEnv_ETA_3_varB.dat", header=F)
Bayesh2by.post<-read.table("results/BGLR_h2_MumEnv_ETA_2_varB.dat", header=F)


#CIs for h2
h2.CIs<-getCIS(Bayesh2.post)
h2.post<-(Bayesh2.post)/(Bayesh2.post+Bayesh2Mum.post+Bayesh2sy.post+Bayesh2by.post+1)
h2.postCIs<-getCIS(h2.post)
h2_genomic<-mean(h2.post$V1)
h2.CIs$lower.h2<-h2.postCIs$lower.bound
h2.CIs$upper.h2<-h2.postCIs$upper.bound
h2.CIs$variance<-varA
h2.CIs$h2<-h2_genomic
h2.CIs$variable<-"genetic"

#CIs for proportion of variance explained of additional random effects
Bayesh2MumCIs<-getCIS(Bayesh2Mum.post)
Mum.h2.post<-(Bayesh2Mum.post)/(Bayesh2.post+Bayesh2Mum.post+Bayesh2sy.post+Bayesh2by.post+1)
Mum.h2.postCIs<-getCIS(Mum.h2.post)
Mum_h2<-mean(Mum.h2.post$V1)
Bayesh2MumCIs$lower.h2<-Mum.h2.postCIs$lower.bound
Bayesh2MumCIs$upper.h2<-Mum.h2.postCIs$upper.bound
Bayesh2MumCIs$variance<-model_h2$ETA[[4]]$varB
Bayesh2MumCIs$h2<-Mum_h2
Bayesh2MumCIs$variable<-"mother"

Bayesh2syCIs<-getCIS(Bayesh2sy.post)
Sy.h2.post<-(Bayesh2sy.post)/(Bayesh2.post+Bayesh2Mum.post+Bayesh2sy.post+Bayesh2by.post+1)
Sy.h2.postCIs<-getCIS(Sy.h2.post)
Sy_h2<-mean(Sy.h2.post$V1)
Bayesh2syCIs$lower.h2<-Sy.h2.postCIs$lower.bound
Bayesh2syCIs$upper.h2<-Sy.h2.postCIs$upper.bound
Bayesh2syCIs$variance<-model_h2$ETA[[3]]$varB
Bayesh2syCIs$h2<-Sy_h2
Bayesh2syCIs$variable<-"second_year"

Bayesh2byCIs<-getCIS(Bayesh2by.post)
By.h2.post<-(Bayesh2by.post)/(Bayesh2.post+Bayesh2Mum.post+Bayesh2sy.post+Bayesh2by.post+1)
By.h2.postCIs<-getCIS(By.h2.post)
By_h2<-mean(By.h2.post$V1)
Bayesh2byCIs$lower.h2<-By.h2.postCIs$lower.bound
Bayesh2byCIs$upper.h2<-By.h2.postCIs$upper.bound
Bayesh2byCIs$variance<-model_h2$ETA[[2]]$varB
Bayesh2byCIs$h2<-By_h2
Bayesh2byCIs$variable<-"birthyear"

h2.CIs<-rbind(h2.CIs, Bayesh2MumCIs)
h2.CIs<-rbind(h2.CIs, Bayesh2byCIs)
h2.CIs<-rbind(h2.CIs, Bayesh2syCIs)
h2.CIs<-h2.CIs[, c("variable","variance", "lower.bound", "upper.bound", "h2", "lower.h2", "upper.h2")]

write.table(h2.CIs, file ="results/Survival_h2.txt", sep = "\t", row.names = F, col.names = T, quote = F)

library(xtable)
print(xtable(h2.CIs, type = "latex", 
             digits= c(1,1,4,4,4,4,4,4) ),
      file = "results/BGLR_h2_variance_effect_tbl.tex")


#sex effect
bhat.sex<-model_h2$ETA[[1]]$b
bhat.sex.post<-read.table("results/BGLR_h2_MumEnv_ETA_1_b.dat", header=T)
bhat.sexCIS<-getCIS(bhat.sex.post)

bhat.sexCIS$bhat<-bhat.sex
bhat.sexCIS.data<-as.data.frame(sapply(bhat.sexCIS, pnorm))
bhat.sexCIS.effect<-bhat.sexCIS.data$`sapply(bhat.sexCIS, pnorm)` - 0.5


#trace plot for Va and h2
varU<-scan('results/BGLR_h2_MumEnv_ETA_5_varU.dat')
plot(varU,type='o',col=2,cex=.5,ylab=expression(var[A]))
abline(h=model_h2$ETA[[5]]$varU,col=4,lwd=2);
abline(v=model_h2$burnIn/model_h2$thin,col=4)

plot((varU/(varU+Bayesh2by.post$V1+Bayesh2Mum.post$V1+Bayesh2sy.post$V1+1)),type='o',col=2,cex=.5,ylab=expression(h^2))
abline(h=h2_genomic,col=4,lwd=2);
abline(v=model_h2$burnIn/model_h2$thin,col=4)








# #get variances
# B<-readBinMat("results/BGLR_BayesB_ETA_4_b.bin")
# plot(B[,1], type = "o", col=4)
# X<-scale(geno_sub)
# VAR<-getVariances(B=B, X=X, sets = sample(1:20, size=1000, replace=TRUE))
# head(VAR)
# 
# 
# demo(ordinal)
# demo(BL)
# 
# #demo getVariances
# data(wheat)
# y=wheat.Y[,1] ; X=scale(wheat.X)
# dir.create('test_saveEffects')
# setwd('test_saveEffects')
# fm=BGLR(y=y,ETA=list(list(X=X,model='BayesB',saveEffects=TRUE)),nIter=12000,thin=2,burnIn=2000)
# B=readBinMat('ETA_1_b.bin')
# plot(B[,1],type='o',col=4)
# VAR=getVariances(B=B,X=X,sets=c(1:1279))
# head(VAR)
# plot(VAR[,"total"])
# 
# VAR.G<-sum(fm$ETA[[1]]$varB)
