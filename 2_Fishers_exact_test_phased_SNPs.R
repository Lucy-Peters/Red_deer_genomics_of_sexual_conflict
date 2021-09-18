## Apply Fisher's exact test to test for transmission distortion on phased SNP data (linkphase and AlphaPeel)##

library(tidyr)
library(plyr)
library(dplyr)
library(reshape2)
library(RODBC)
library(ggplot2)
library(xtable)

#create tables to use as input for Fisher's exact test 

#get sex and birthyear data (by maybe needed in later models)
db<-"C:\\Users\\s1767711\\Documents\\Red_deer\\database\\current\\RedDeer1.93.accdb"
con<-odbcConnectAccess2007(db)
queries<-sqlTables(con, tableType = "VIEW")
by_sex_data<-sqlFetch(con, "Birthyear_and_sex_data")
odbcClose(con)

head(by_sex_data)
colnames(by_sex_data)<-c("ID", "sex", "birthyear") #sex in database file (2=male, 1=female)
#noticed ID change between database birthyear/sex data and new pedigree: MINNI in pedigree is DIL01 in birthyear/sex data
#make ID same as in pedigree file
by_sex_data$ID<-gsub("^DIL01$", "MINNI", by_sex_data$ID)

write.table(by_sex_data, file="Birthyear_sex_data_2019.txt", sep = "\t", col.names = T, row.names = F, quote = F)

by_sex_data<-read.table("data/Birthyear_sex_data_2019.txt", header=TRUE)

#ID_index_link<-read.table("ID_conversion_linkphase.txt", header=T)
ID_index_link<-read.table("ID_conversion_linkphase_new_ped.txt", header=T)

#list of genotyped IDs
IDs_genotyped<-read.table("results/IDs_genotyped_reclustered_SNPs.txt", header=T)

#get information about individuals in simulations
load("results/Deer_genedrop_objects_sim_part1.RData")
genedrop_obj<-genedrop_list[[1]]
rm(genedrop_list)
genedrop_ped<-genedrop_obj@pedigree %>%data.frame()
colnames(genedrop_ped)<-c("ID", "Father", "Mother", "sex", "Cohort")
rm(genedrop_obj)

#ped<-read.csv("Pedigree_2019-07-02_DB.csv", header = T)
ped<-read.table("results/Pedigree_deer_2017_19_merged.txt", header=T)
head(ped)
ped_dup<-subset(ped, duplicated(ID))
ped<-subset(ped, !duplicated(ID))

#MAF_compare<-read.table("Compare_MAF_pre_post_phasing_reclustered.txt", header=T)
#MAF_compare_sub<-subset(MAF_compare, unphased_fraction<0.16)

#loop for making TD tables

#lists with compiled data tables

master_tab_list<-list()
Mumtrans_list<-list()
Dadtrans_list<-list()
offspringtrans_list<-list()
Mum_off_sextrans_list<-list()
Dad_off_sextrans_list<-list()

#subset df lists
# Mumtrans_sub_list<-list()
# Dadtrans_sub_list<-list()
# offspringtrans_sub_list<-list()
# Mum_off_sextrans_sub_list<-list()
# Dad_off_sextrans_sub_list<-list()

SNP_error_all<-data.frame()

for (chr in c(1:33)) {
  
  #LINKPHASE phase files
  #assign("phase_chr", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/Linkphase_run/Results/phasesLINK_chr", 
                                         #chr), header=F))
  
  #LINKPHASE marker files
  #assign("marker_map_chr", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/Linkphase_run/Input_files/Deer_LINKPHASE_marker_file_new_chr", 
                                             #chr,".txt"), header=F))
  
  #AlphaPeel phase files
  assign("phase_chr", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaPeel_run/Results/AlphaPeel_out_chr", 
                                        chr, ".called_phase.0.95"), header=F))
  
  #AlphaPeel marker files not required for AlphaPeel but nedded to assign SNP names
  assign("marker_map_chr", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaPeel_run/Input_files/Deer_AlphaPeel_marker_file_chr", 
                                        chr,".txt"), header=T))
  
  colnames(marker_map_chr)<-c("marker.order", "SNP.name", "Est.Mb.pos")
  
  #LINKPHASE
  # colnames(phase_chr)<-c("ID_new", "hap.origin", as.character(marker_map_chr$SNP.name))
  # SNP.names<-c(as.character(marker_map_chr$SNP.name))
  # phase_chr<-join(phase_chr, ID_index_link, by="ID_new")
  # phase_chr<-phase_chr[, c("ID", "hap.origin", SNP.names)]
  # phase_chr$hap.origin<-ifelse(phase_chr$hap.origin == 1, "paternal", "maternal")
  
  #AlphaPeel
  colnames(phase_chr)<-c("ID", as.character(marker_map_chr$SNP.name))
  SNP.names<-c(as.character(marker_map_chr$SNP.name))
  phase_chr<-subset(phase_chr, ID %in% genedrop_ped$ID)
  phase_chr<-subset(phase_chr, ID %in% IDs_genotyped$ID)
  hap.origin.names<-as.vector(c("paternal", "maternal"))
  phase_chr$hap.origin<-rep(hap.origin.names, (nrow(phase_chr)/2))
  
  phase_chr_melt<- reshape2::melt(phase_chr, id.vars = c("ID", "hap.origin"))
  colnames(phase_chr_melt)[3]<-"SNP.name"
  phase_chr_cast<- reshape2::dcast(phase_chr_melt, SNP.name+ID~hap.origin, value.var = "value" )
  
  phase_chr_cast<-join(phase_chr_cast, ped, by="ID")
  phase_chr_cast<-join(phase_chr_cast, by_sex_data, by="ID")
  
  phase_chr_cast<-phase_chr_cast%>%
    mutate(Chromosome=paste0(chr))
  
  phase_chr_cast$sex<-ifelse(phase_chr_cast$sex==1, "female", 
                             ifelse(phase_chr_cast$sex==2,"male", "unknown"))
  
  #allocate random alle labels 
  
  phase_chr_random<-data.frame()
  
  for (snp in SNP.names) {
    phase_chr_cast_sub<-subset(phase_chr_cast, SNP.name==snp)
    
    draw<-rbinom(1, 1, 0.5)
    if(draw==1){
      
      #LINKPHASE
      #allele 1 = A, allele 2 = B
      # phase_chr_cast_sub$maternal<-ifelse(phase_chr_cast_sub$maternal==1, "A",
      #                                 ifelse(phase_chr_cast_sub$maternal==2, "B", NA))
      # 
      # phase_chr_cast_sub$paternal<-ifelse(phase_chr_cast_sub$paternal==1, "A",
      #                                 ifelse(phase_chr_cast_sub$paternal==2, "B", NA))
      # 
      # phase_chr_cast_sub$minor.allele<-"A"
       
      #AlphaPeel
      #allele 0 = A, allele 1 = B
      phase_chr_cast_sub$maternal<-ifelse(phase_chr_cast_sub$maternal==0, "A",
                                          ifelse(phase_chr_cast_sub$maternal==1, "B", NA))
      
      phase_chr_cast_sub$paternal<-ifelse(phase_chr_cast_sub$paternal==0, "A",
                                          ifelse(phase_chr_cast_sub$paternal==1, "B", NA))
      
      phase_chr_cast_sub$minor.allele<-"A"
      
      
      
    } else {
      
      #LINKPHASE
      #allele 1 = B, allele 2 = A
      # phase_chr_cast_sub$maternal<-ifelse(phase_chr_cast_sub$maternal==1, "B",
      #                                 ifelse(phase_chr_cast_sub$maternal==2, "A", NA))
      # 
      # phase_chr_cast_sub$paternal<-ifelse(phase_chr_cast_sub$paternal==1, "B",
      #                                 ifelse(phase_chr_cast_sub$paternal==2, "A", NA))
      # 
      # phase_chr_cast_sub$minor.allele<-"B"
      
      #AlphaPeel
      #allele 0 = B, allele 1 = A
      phase_chr_cast_sub$maternal<-ifelse(phase_chr_cast_sub$maternal==0, "B",
                                          ifelse(phase_chr_cast_sub$maternal==1, "A", NA))
      
      phase_chr_cast_sub$paternal<-ifelse(phase_chr_cast_sub$paternal==0, "B",
                                          ifelse(phase_chr_cast_sub$paternal==1, "A", NA))
      
      phase_chr_cast_sub$minor.allele<-"B"
      
    }
    
    phase_chr_random<-rbind(phase_chr_random, phase_chr_cast_sub)
  }
  
  
  phase_chr_random_mum<-phase_chr_random[, c("SNP.name", "ID", "maternal", "paternal")]
  colnames(phase_chr_random_mum)<-c("SNP.name", "Mother", "maternal_mother", "paternal_mother")
  phase_chr_random<-join(phase_chr_random, phase_chr_random_mum)
  
  phase_chr_random_dad<-phase_chr_random[, c("SNP.name", "ID", "maternal", "paternal")]
  colnames(phase_chr_random_dad)<-c("SNP.name", "Father", "maternal_father", "paternal_father")
  phase_chr_random<-join(phase_chr_random, phase_chr_random_dad)
  
  phase_chr_random<-subset(phase_chr_random, !(is.na(maternal)&is.na(paternal)))
  phase_chr_random[is.na(phase_chr_random)]<-"Z"
  
  phase_chr_random$informative<-ifelse(phase_chr_random$maternal_mother!=phase_chr_random$paternal_mother & 
                                         phase_chr_random$maternal_mother!="Z" & phase_chr_random$paternal_mother!="Z" &
                                         phase_chr_random$maternal_father==phase_chr_random$paternal_father, "yes_mat", 
                                       ifelse(phase_chr_random$maternal_father!=phase_chr_random$paternal_father & 
                                                phase_chr_random$maternal_father!="Z" & phase_chr_random$paternal_father!="Z" &
                                                phase_chr_random$maternal_mother==phase_chr_random$paternal_mother, "yes_pat",
                                              ifelse(phase_chr_random$maternal_mother!=phase_chr_random$paternal_mother &
                                                       phase_chr_random$maternal_mother!="Z" & phase_chr_random$paternal_mother!="Z" &
                                                       phase_chr_random$maternal_father!=phase_chr_random$paternal_father &
                                                       phase_chr_random$maternal_father!= "Z" & phase_chr_random$paternal_father!="Z", "yes_both", "no")))
  
  phase_chr_random_clean<-subset(phase_chr_random, informative=="yes_mat"|informative=="yes_pat"|informative=="yes_both")
  
  phase_chr_random_clean$mend_error<-ifelse(phase_chr_random_clean$maternal != phase_chr_random_clean$maternal_mother &
                                              phase_chr_random_clean$maternal != phase_chr_random_clean$paternal_mother &
                                              phase_chr_random_clean$maternal !="Z"& phase_chr_random_clean$maternal_mother !="Z" &
                                              phase_chr_random_clean$paternal_mother !="Z","mat_error",
                                            ifelse(phase_chr_random_clean$paternal != phase_chr_random_clean$maternal_father &
                                                     phase_chr_random_clean$paternal != phase_chr_random_clean$paternal_father &
                                                     phase_chr_random_clean$paternal !="Z" & phase_chr_random_clean$maternal_father !="Z" &
                                                     phase_chr_random_clean$paternal_father !="Z","pat_error", "no_error"))
  
  
  phase_chr_random_clean$maternal[which(phase_chr_random_clean$maternal=="Z")]<-NA
  phase_chr_random_clean$paternal[which(phase_chr_random_clean$paternal=="Z")]<-NA
  
  SNP_error<-phase_chr_random_clean %>%
    filter(mend_error!="no_error") %>%
    group_by(SNP.name)%>%
    summarize(error.count=n()) %>%
    na.omit
  
  SNP_error<-as.data.frame(SNP_error)
  SNP_error_all<-rbind(SNP_error_all, SNP_error)
  
  # #exclude any transmissions with Mendelian error
  phase_chr_random_clean<-subset(phase_chr_random_clean, !mend_error=="mat_error" & !mend_error=="pat_error")


  #do all the following transmission frequency calculations both on full data set and on subset filtered for SNP
  #phasing rate

  #summarise transmission frequencies by transmitting parent and convert into format that can be used for
  #Fisher's exact test (dcast - 1 row per SNP, columns are alleles A and B, values are frequencies )

  Mumtrans<-subset(phase_chr_random_clean, informative=="yes_mat"|informative=="yes_both")
  # #Mumtrans_sub<-subset(Mumtrans, SNP.name %in% MAF_compare_sub$SNP.name)

  Mumtrans<-Mumtrans %>%
    group_by(SNP.name, maternal)%>%
    summarise(freq=n())%>%
    na.omit

  Mumtrans <- reshape2::dcast(Mumtrans, SNP.name ~ maternal)

  #  # Mumtrans_sub<-Mumtrans_sub %>%
  #  #   group_by(SNP.name, maternal)%>%
  #  #   summarise(freq=n())%>%
  #  #   na.omit
  #  # Mumtrans_sub <- dcast(Mumtrans_sub, SNP.name ~ maternal)


  Dadtrans<-subset(phase_chr_random_clean, informative=="yes_pat"|informative=="yes_both")
  # #Dadtrans_sub<-subset(Dadtrans, SNP.name %in% MAF_compare_sub$SNP.name)

  Dadtrans<-Dadtrans %>%
    group_by(SNP.name, paternal)%>%
    summarise(freq=n())%>%
    na.omit

  Dadtrans <-reshape2::dcast(Dadtrans, SNP.name ~ paternal)

  #  # Dadtrans_sub<-Dadtrans_sub %>%
  #  #   group_by(SNP.name, paternal)%>%
  #  #   summarise(freq=n())%>%
  #  #   na.omit
  #  # 
  #  # Dadtrans_sub <- dcast(Dadtrans_sub, SNP.name ~ paternal)                         


  #summarise transmission frequencies by transmitting parent and offspring sex
  Mum_off_sextrans<-subset(phase_chr_random_clean, (informative=="yes_mat"|informative=="yes_both") & !sex=="Z")
  Mum_off_sextrans<-subset(Mum_off_sextrans, !sex=="unknown")
  # #Mum_off_sextrans_sub<-subset(Mum_off_sextrans, SNP.name %in% MAF_compare_sub$SNP.name)

  Mum_off_sextrans<-Mum_off_sextrans %>%
    group_by(SNP.name, maternal, sex) %>%
    summarise(freq=n()) %>%
    na.omit

  Mum_off_sextrans<-reshape2::dcast(Mum_off_sextrans, SNP.name+sex ~ maternal)

  #  # Mum_off_sextrans_sub<-Mum_off_sextrans_sub %>%
  #  #   group_by(SNP.name, maternal, sex) %>%
  #  #   summarise(freq=n()) %>%
  #  #   na.omit
  #  # 
  #  # Mum_off_sextrans_sub<-dcast(Mum_off_sextrans_sub, SNP.name+sex ~ maternal)


  Dad_off_sextrans<-subset(phase_chr_random_clean, (informative=="yes_pat"|informative=="yes_both") & !sex=="Z")
  Dad_off_sextrans<-subset(Dad_off_sextrans, !sex=="unknown")
  # #Dad_off_sextrans_sub<-subset(Dad_off_sextrans, SNP.name %in% MAF_compare_sub$SNP.name)

  Dad_off_sextrans<-Dad_off_sextrans %>%
    group_by(SNP.name, paternal, sex) %>%
    summarise(freq=n()) %>%
    na.omit

  Dad_off_sextrans<- reshape2::dcast(Dad_off_sextrans, SNP.name+sex ~ paternal)

  #  # Dad_off_sextrans_sub<-Dad_off_sextrans_sub %>%
  #  #   group_by(SNP.name, paternal, sex) %>%
  #  #   summarise(freq=n()) %>%
  #  #   na.omit
  #  # 
  #  # Dad_off_sextrans_sub<-dcast(Dad_off_sextrans_sub, SNP.name+sex ~ paternal)


  #offspring TD - TD by offspring sex (Mum and Dad allele togehter)
  otab<-phase_chr_random_clean[c("SNP.name","ID", "sex", "maternal", "paternal", "informative")]

  otab_pat<-subset(otab, informative=="yes_pat")
  otab_pat<-otab_pat[c("SNP.name", "ID", "sex", "paternal")]
  names(otab_pat)[names(otab_pat) == "paternal"] <- "value"

  otab_mat<-subset(otab, informative=="yes_mat")
  otab_mat<-otab_mat[c("SNP.name", "ID", "sex", "maternal")]
  names(otab_mat)[names(otab_mat) == "maternal"] <- "value"

  otab_both<-subset(otab, informative=="yes_both")
  otab_both<-otab_both[c("SNP.name", "ID", "sex", "maternal", "paternal")]
  otab_both <- reshape2::melt(otab_both, id.vars = c("SNP.name","ID", "sex"))
  otab_both<-otab_both[c("SNP.name", "ID", "sex", "value")]

  offspringtrans<-rbind(otab_both, otab_mat)
  offspringtrans<-rbind(offspringtrans, otab_pat)
  offspringtrans<-subset(offspringtrans, !sex=="Z")

  # #offspringtrans_sub<-subset(offspringtrans, SNP.name %in% MAF_compare_sub$SNP.name)

  #offspring transmission
  offspringtrans <- offspringtrans %>%
    group_by(SNP.name, sex, value) %>%
    summarize(freq=n()) %>%
    na.omit


  names(offspringtrans)[names(offspringtrans) == "value"] <- "allele"
  offspringtrans <- reshape2::dcast(offspringtrans, SNP.name + sex ~ allele) #dcast = oppposite of melt

  #  # offspringtrans_sub <- offspringtrans_sub %>%
  #  #   group_by(SNP.name, sex, value) %>%
  #  #   summarize(freq=n()) %>%
  #  #   na.omit
  #  # 
  #  # 
  #  # names(offspringtrans_sub)[names(offspringtrans_sub) == "value"] <- "allele"
  #  # offspringtrans_sub <- dcast(offspringtrans_sub, SNP.name + sex ~ allele) #dcast = oppposite of melt 


  #add chromosome minor allele and mendelian error information
  phase_chr_random_clean_single<-subset(phase_chr_random_clean, !duplicated(SNP.name))
  Mumtrans<-join(Mumtrans, phase_chr_random_clean_single[c("SNP.name", "Chromosome", "minor.allele")], by="SNP.name", type="inner")
  Mumtrans<-join(Mumtrans, SNP_error, by="SNP.name")
  Dadtrans<-join(Dadtrans, phase_chr_random_clean_single[c("SNP.name", "Chromosome", "minor.allele")], by="SNP.name")
  Dadtrans<-join(Dadtrans, SNP_error, by="SNP.name")
  offspringtrans<-join(offspringtrans, phase_chr_random_clean_single[c("SNP.name", "Chromosome", "minor.allele")], by="SNP.name")
  offspringtrans<-join(offspringtrans, SNP_error, by="SNP.name")
  Mum_off_sextrans<-join(Mum_off_sextrans, phase_chr_random_clean_single[c("SNP.name", "Chromosome", "minor.allele")], by="SNP.name")
  Mum_off_sextrans<-join(Mum_off_sextrans, SNP_error, by="SNP.name")
  Dad_off_sextrans<-join(Dad_off_sextrans, phase_chr_random_clean_single[c("SNP.name", "Chromosome", "minor.allele")], by="SNP.name")
  Dad_off_sextrans<-join(Dad_off_sextrans, SNP_error, by="SNP.name")

  # #subsets (filtered for phasing rate)
  # # Mumtrans_sub<-join(Mumtrans_sub, phase_chr_random_clean_single[c("SNP.name", "Chromosome", "minor.allele")], by="SNP.name", type="inner")
  # # Mumtrans_sub<-join(Mumtrans_sub, SNP_error, by="SNP.name")
  # # Dadtrans_sub<-join(Dadtrans_sub, phase_chr_random_clean_single[c("SNP.name", "Chromosome", "minor.allele")], by="SNP.name")
  # # Dadtrans_sub<-join(Dadtrans_sub, SNP_error, by="SNP.name")
  # # offspringtrans_sub<-join(offspringtrans_sub, phase_chr_random_clean_single[c("SNP.name", "Chromosome", "minor.allele")], by="SNP.name")
  # # offspringtrans_sub<-join(offspringtrans_sub, SNP_error, by="SNP.name")
  # # Mum_off_sextrans_sub<-join(Mum_off_sextrans_sub, phase_chr_random_clean_single[c("SNP.name", "Chromosome", "minor.allele")], by="SNP.name")
  # # Mum_off_sextrans_sub<-join(Mum_off_sextrans_sub, SNP_error, by="SNP.name")
  # # Dad_off_sextrans_sub<-join(Dad_off_sextrans_sub, phase_chr_random_clean_single[c("SNP.name", "Chromosome", "minor.allele")], by="SNP.name")
  # # Dad_off_sextrans_sub<-join(Dad_off_sextrans_sub, SNP_error, by="SNP.name")


  #append to data lists
  master_tab_list[[chr]]<-phase_chr_random_clean
  Mumtrans_list[[chr]]<-Mumtrans
  Dadtrans_list[[chr]]<-Dadtrans
  offspringtrans_list[[chr]]<-offspringtrans
  Mum_off_sextrans_list[[chr]]<-Mum_off_sextrans
  Dad_off_sextrans_list[[chr]]<-Dad_off_sextrans

  #sublists
  # Mumtrans_sub_list[[chr]]<-Mumtrans_sub
  # Dadtrans_sub_list[[chr]]<-Dadtrans_sub
  # offspringtrans_sub_list[[chr]]<-offspringtrans_sub
  # Mum_off_sextrans_sub_list[[chr]]<-Mum_off_sextrans_sub
  # Dad_off_sextrans_sub_list[[chr]]<-Dad_off_sextrans_sub
  
  
}

#save TD table lists
save(master_tab_list, file="TD_master_table_list_phased_AlphaPeel.RData")
save(Mumtrans_list, file="TD_mother_allele_list_phased_AlphaPeel.RData")
save(Dadtrans_list, file="TD_father_allele_list_phased_AlphaPeel.RData")
save(offspringtrans_list, file="TD_offspring_sex_list_phased_AlphaPeel.RData")
save(Mum_off_sextrans_list, file="TD_mother_allele_offspring_sex_list_phased_AlphaPeel.RData")
save(Dad_off_sextrans_list, file="TD_father_allele_offspring_sex_list_phased_AlphaPeel.RData")

#save sublists
save(Mumtrans_sub_list, file="TD_mother_allele_sub_list_phased_reclustered.RData")
save(Dadtrans_sub_list, file="TD_father_allele_sub_list_phased_reclustered.RData")
save(offspringtrans_sub_list, file="TD_offspring_sex_sub_list_phased_reclustered.RData")
save(Mum_off_sextrans_sub_list, file="TD_mother_allele_offspring_sex_sub_list_phased_reclustered.RData")
save(Dad_off_sextrans_sub_list, file="TD_father_allele_offspring_sex_sub_list_phased_reclustered.RData")

#save Mendelian error count table
#AlphaPeel run
write.table(SNP_error_all, file="Mendelian_error_count_AlphaPeel_run.txt", sep="\t", row.names = F, quote=F)
hist(SNP_error_all$error.count)


load("TD_master_table_list_phased_reclustered.RData")
df_master_table<-master_tab_list[[1]]

#load transmitted allele data frame list files and apply function to calculate Fisher's exact test per SNP for alleles transmitted
#by the mother, the father and for each the mother or father transmitted allele also split by offspring sex
source("scripts/TD_Fishers_test.R")

#alleles transmitted by mother
#load("TD_mother_allele_list_phased.RData")
#load("TD_mother_allele_list_phased_reclustered.RData")
#load("TD_mother_allele_sub_list_phased_reclustered.RData")
load("TD_mother_allele_list_phased_AlphaPeel.RData")
Mumtrans_result<-lapply(Mumtrans_list, TD_Fishers_test)
#Mumtrans_result<-lapply(Mumtrans_sub_list, TD_Fishers_test)
Mumtrans_result_df<-Mumtrans_result[[1]]
#save(Mumtrans_result, file="TD_mother_allele_results_list_phased.RData")
#save(Mumtrans_result, file="TD_mother_allele_results_list_phased_reclustered.RData")
save(Mumtrans_result, file="results/TD_mother_allele_results_list_phased_AlphaPeel.RData")
#save(Mumtrans_result, file="TD_mother_allele_results_sub_list_phased_reclustered.RData")


#alleles transmitted by mother split by offspring sex
#load("TD_mother_allele_offspring_sex_list_phased.RData")
#load("TD_mother_allele_offspring_sex_list_phased_reclustered.RData")
#load("TD_mother_allele_offspring_sex_sub_list_phased_reclustered.RData")
load("TD_mother_allele_offspring_sex_list_phased_AlphaPeel.RData")
Mum_off_sextrans_result<-lapply(Mum_off_sextrans_list, TD_Fishers_test)
#Mum_off_sextrans_result<-lapply(Mum_off_sextrans_sub_list, TD_Fishers_test)
Mum_off_sextrans_result_df<-Mum_off_sextrans_result[[1]]
#save(Mum_off_sextrans_result, file="TD_mother_allele_offspring_sex_results_list_phased.RData")
#save(Mum_off_sextrans_result, file="TD_mother_allele_offspring_sex_results_list_phased_reclustered.RData")
#save(Mum_off_sextrans_result, file="TD_mother_allele_offspring_sex_results_sub_list_phased_reclustered.RData")
save(Mum_off_sextrans_result, file="results/TD_mother_allele_offspring_sex_results_list_phased_AlphaPeel.RData")


#alleles transmitted by father
#load("TD_father_allele_list_phased.RData")
#load("TD_father_allele_list_phased_reclustered.RData")
#load("TD_father_allele_sub_list_phased_reclustered.RData")
load("TD_father_allele_list_phased_AlphaPeel.RData")
Dadtrans_result<-lapply(Dadtrans_list, TD_Fishers_test)
#Dadtrans_result<-lapply(Dadtrans_sub_list, TD_Fishers_test)
Dadtrans_result_df<-Dadtrans_result[[1]]
#save(Dadtrans_result, file="TD_father_allele_results_list_phased.RData")
#save(Dadtrans_result, file="TD_father_allele_results_list_phased_reclustered.RData")
#save(Dadtrans_result, file="TD_father_allele_results_sub_list_phased_reclustered.RData")
save(Dadtrans_result, file="results/TD_father_allele_results_list_phased_AlphaPeel.RData")


#alleles transmitted by father split by offspring sex
#load("TD_father_allele_offspring_sex_list_phased.RData")
#load("TD_father_allele_offspring_sex_list_phased_reclustered.RData")
#load("TD_father_allele_offspring_sex_sub_list_phased_reclustered.RData")
load("TD_father_allele_offspring_sex_list_phased_AlphaPeel.RData")
Dad_off_sextrans_result<-lapply(Dad_off_sextrans_list, TD_Fishers_test)
#Dad_off_sextrans_result<-lapply(Dad_off_sextrans_sub_list, TD_Fishers_test)
Dad_off_sextrans_result_df<-Dad_off_sextrans_result[[1]]
#save(Dad_off_sextrans_result, file="TD_father_allele_offspring_sex_results_list_phased.RData")
#save(Dad_off_sextrans_result, file="TD_father_allele_offspring_sex_results_list_phased_reclustered.RData")
#save(Dad_off_sextrans_result, file="TD_father_allele_offspring_sex_results_sub_list_phased_reclustered.RData")
save(Dad_off_sextrans_result, file="results/TD_father_allele_offspring_sex_results_list_phased_AlphaPeel.RData")


#alleles transmitted by mother or father split by offspring sex
#load("TD_offspring_sex_list_phased.RData")
#load("TD_offspring_sex_list_phased_reclustered.RData")
#load("TD_offspring_sex_sub_list_phased_reclustered.RData")
load("TD_offspring_sex_list_phased_AlphaPeel.RData")
offspringtrans_result<-lapply(offspringtrans_list, TD_Fishers_test)
#offspringtrans_result<-lapply(offspringtrans_sub_list, TD_Fishers_test)
offspringtrans_result_df<-offspringtrans_result[[1]]
#save(offspringtrans_result, file="TD_offspring_sex_results_list_phased.RData")
#save(offspringtrans_result, file="TD_offspring_sex_results_list_phased_reclustered.RData")
#save(offspringtrans_result, file="TD_offspring_sex_results_sub_list_phased_reclustered.RData")
save(offspringtrans_result, file="results/TD_offspring_sex_results_list_phased_AlphaPeel.RData")


  # load("TD_mother_allele_offspring_sex_results_list_phased.RData")
  # load("TD_mother_allele_results_list_phased.RData")
  # load("TD_father_allele_offspring_sex_results_list_phased.RData")
  # load("TD_father_allele_results_list_phased.RData")
  # load("TD_offspring_sex_results_list_phased.RData")

# load("TD_mother_allele_offspring_sex_results_list_phased_reclustered.RData")
# load("TD_mother_allele_results_list_phased_reclustered.RData")
# load("TD_father_allele_offspring_sex_results_list_phased_reclustered.RData")
# load("TD_father_allele_results_list_phased_reclustered.RData")
# load("TD_offspring_sex_results_list_phased_reclustered.RData")


# load("TD_mother_allele_offspring_sex_results_sub_list_phased_reclustered.RData")
# load("TD_mother_allele_results_sub_list_phased_reclustered.RData")
# load("TD_father_allele_offspring_sex_results_sub_list_phased_reclustered.RData")
# load("TD_father_allele_results_sub_list_phased_reclustered.RData")
# load("TD_offspring_sex_results_sub_list_phased_reclustered.RData")


load("results/TD_mother_allele_offspring_sex_results_list_phased_AlphaPeel.RData")
load("results/TD_mother_allele_results_list_phased_AlphaPeel.RData")
load("results/TD_father_allele_offspring_sex_results_list_phased_AlphaPeel.RData")
load("results/TD_father_allele_results_list_phased_AlphaPeel.RData")
load("results/TD_offspring_sex_results_list_phased_AlphaPeel.RData")

results_superlist<-list("Mum_off_sextrans" = Mum_off_sextrans_result, "Mumtrans" = Mumtrans_result ,"Dad_off_sextrans"=Dad_off_sextrans_result, 
                        "Dadtrans" = Dadtrans_result,"offspringtrans"= offspringtrans_result)

list_names<-c("Mum_off_sextrans", "Mumtrans", "Dad_off_sextrans", "Dadtrans", "offspringtrans")


for (l.name in list_names){
  res_l<-results_superlist[paste0(l.name)]
  df_all<-data.frame()
  for (chr in c(1:33)) {
    df<-res_l[[1]]
    df<-df[[chr]]
    df_all<-rbind(df_all, df)
  }
  assign(paste0(l.name, "_results_all_chr"), df_all)
}

#plot results as manhattan plots
link_map<-read.table("data/Cervus_elaphus_linkage_map_data.ordered.txt", header=T)
link_map_part<-link_map[, c("SNP.Name", "CEL.order")]
source("scripts/ggplot_manhattan.R")

#plot results for TD by offspring sex (disregarding parent sex)
offspringtrans_results_all_chr<-subset(offspringtrans_results_all_chr, !sex=="unknown")


p_offspringtrans<-manhattan_plot(df_plot = offspringtrans_results_all_chr, link_order = link_map_part, chr_no = 33, 
               p_threshold = 0.02, ylim_max=3,xlab_text = "CEL linkage group")

# ggsave("offspringtrans_sex_split_manhattan_plot_phased_reclustered.png", p_offspringtrans, width = 25,
#        height = 25, units = "cm")
 
ggsave("plots/offspringtrans_sex_split_manhattan_plot_AlphaPeel.png", p_offspringtrans, width = 25,
        height = 25, units = "cm")
 
# ggsave("offspringtrans_sex_split_manhattan_plot_phased_reclustered_sub.png", p_offspringtrans, width = 25,
#        height = 25, units = "cm")

# ggsave("offspringtrans_sex_split_manhattan_plot_phased.png", p_offspringtrans, width = 25,
#        height = 25, units = "cm")



# write.table(offspringtrans_results_all_chr, file="offspringtrans_TD_results_phased_reclustered.txt", sep="\t", 
#             col.names = T, row.names = F, quote=F)
 
 write.table(offspringtrans_results_all_chr, file="results/offspringtrans_TD_results_AlphaPeel.txt", sep="\t", 
             col.names = T, row.names = F, quote=F)
 

# write.table(offspringtrans_df.plot, file="offspringtrans_TD_results_phased_reclustered_sub.txt", sep="\t", 
#             col.names = T, row.names = F, quote=F)

 # write.table(offspringtrans_results_all_chr, file="offspringtrans_TD_results_phased.txt", sep="\t", 
 #                          col.names = T, row.names = F, quote=F)

#filter for top SNPs (at p_value , 0.02 based on genedrop simulations) -reclustered phased data (AlphaPeel)

offspringtrans_df_top<-subset(offspringtrans_results_all_chr, p_value<0.02)
write.table(offspringtrans_df_top, file="results/offspringtrans_TD_top_SNPs_AlphaPeel.txt.", sep="\t",
            row.names = F, col.names = T, quote=F)


#compare phasing rate to p_values
compare_MAF_tbl<-read.table("results/Compare_MAF_pre_post_phasing_AlphaPeel.txt", header = T)
offspringtrans_df.plot<-join(offspringtrans_results_all_chr, compare_MAF_tbl)
offspringtrans_df.plot[is.na(offspringtrans_df.plot)]<-0

p_phasing_p_value<-ggplot(offspringtrans_df.plot, aes((1-prop.unphased_peel), -log10(p_value)))+geom_point()+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))+
  xlab("phasing rate")+ylab(expression(-log[10](p)))

ggsave("plots/Phasing_vs_p_value_offspringtrans_phased.AlphaPeel.png", p_phasing_p_value, width = 25, height = 25, units = "cm")


#deviation plot

#offspringtrans_results_all_chr<-read.table("offspringtrans_TD_results_phased_reclustered.txt", header=T)
offspringtrans_results_all_chr<-read.table("results/offspringtrans_TD_results_AlphaPeel.txt", header=T)

offspringtrans_results_all_chr<-subset(offspringtrans_results_all_chr, !sex=="unknown")

#exclude SNPs with extremly low allele counts
freq.cutoff<-data.frame(freq.quantile=as.numeric(quantile(offspringtrans_results_all_chr$expected_freq, probs = 0.05)))
offspringtrans_results_all_chr<-subset(offspringtrans_results_all_chr, expected_freq > freq.cutoff$freq.quantile)

#exlude SNPs where allele counts are only available for one sex 
#(want to plot male vs female so need information for both sex per SNP)

offspringtrans_results_all_chr_sub<-subset(offspringtrans_results_all_chr, sex=="female")
offspringtrans_results_all_chr_sub2<-subset(offspringtrans_results_all_chr, sex=="male")

SNPs_uneven_sex<-subset(offspringtrans_results_all_chr_sub2, !SNP.name %in% offspringtrans_results_all_chr_sub$SNP.name)
offspringtrans_results_all_chr<-subset(offspringtrans_results_all_chr, !SNP.name %in% SNPs_uneven_sex$SNP.name)
offspringtrans_results_all_chr<-droplevels(offspringtrans_results_all_chr)

SNPs_uneven_sex2<-subset(offspringtrans_results_all_chr_sub, !SNP.name %in% offspringtrans_results_all_chr_sub2$SNP.name)
offspringtrans_results_all_chr<-subset(offspringtrans_results_all_chr, !SNP.name %in% SNPs_uneven_sex2$SNP.name)
offspringtrans_results_all_chr<-droplevels(offspringtrans_results_all_chr)

#AlphaPeel data after null distribution simulation
offspringtrans_results_all_chr[["significance"]]<-ifelse(offspringtrans_results_all_chr[["p_value"]]<0.02 , "sig", "non.sig")

table(offspringtrans_results_all_chr$significance)

offspringtrans_results_all_chr[["Deviation"]] <- (offspringtrans_results_all_chr[["A"]]- offspringtrans_results_all_chr$expected_freq)/offspringtrans_results_all_chr$expected_freq

offspringtrans_results_deviations <- reshape2::dcast(offspringtrans_results_all_chr[,c("SNP.name", "sex", "Deviation")], formula = SNP.name ~ sex)
offspringtrans_results_status <- reshape2::dcast(offspringtrans_results_all_chr[,c("SNP.name", "sex", "significance")], formula = SNP.name ~ sex)

names(offspringtrans_results_deviations) <- c("SNP", "Female.Deviation", "Male.Deviation")
names(offspringtrans_results_status) <- c("SNP", "Female.Status", "Male.Status")


offspringtrans_results_dev_plot<-join(offspringtrans_results_deviations, offspringtrans_results_status)

offspringtrans_results_dev_plot[["Label"]]<-ifelse(offspringtrans_results_dev_plot[["Female.Status"]]=="sig" & offspringtrans_results_dev_plot[["Male.Status"]]=="sig", "sig",
                                         ifelse(offspringtrans_results_dev_plot[["Female.Status"]]=="sig" &
                                                  offspringtrans_results_dev_plot[["Male.Status"]]=="non.sig", "sig.female",
                                                ifelse(offspringtrans_results_dev_plot[["Male.Status"]]=="sig" &
                                                         offspringtrans_results_dev_plot[["Female.Status"]]=="non.sig", "sig.male", "non.sig")))


offspringtrans_results_dev_plot[["Alpha"]] <- ifelse(offspringtrans_results_dev_plot$Label == "non.Sig.", 0.1, 0.7)

#what are the minimum and maximum deviations of significant SNPs for males and females?
offspringtrans_deviations_sig_fem<-subset(offspringtrans_results_dev_plot, Label == "sig" | Label == "sig.female")
min.female<-min(abs(offspringtrans_deviations_sig_fem$Female.Deviation)) #-0.32773 or absolute 0.10956
max.female<-max(abs(offspringtrans_deviations_sig_fem$Female.Deviation)) #0.33333333

offspringtrans_deviations_sig_male<-subset(offspringtrans_results_dev_plot, Label == "sig" | Label == "sig.male")  

min.male<-min(abs(offspringtrans_deviations_sig_male$Male.Deviation)) #-0.32727272 or absolute 0.107754
max.male<-max(abs(offspringtrans_deviations_sig_male$Male.Deviation)) #0.18734793 or absolute 0.327272

#save offspringtrans TD table with deviation information
write.table(offspringtrans_results_all_chr, file = "results/offspringtrans_TD_results_AlphaPeel.txt", sep="\t", 
            col.names = T, row.names = F, quote=F)

  
p_var_dev_off<-ggplot(offspringtrans_results_dev_plot, aes(Female.Deviation, Male.Deviation, colour = Label, alpha = Alpha)) +
  geom_point(size=2) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Female Deviation", y = "Male Deviation") +
  scale_colour_manual(values = c("grey70", "red", "blue"), labels=c("Not significant", "Females only", "Males only"))+
  guides(alpha=FALSE)+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=16, color = "black"), legend.background = element_rect(size=0.5, 
                                                                                               linetype = "solid", color="black"),
        axis.text = element_text(size=22, color="black"), axis.title = element_text(size=26),
        axis.ticks=element_line(size=1))+
  scale_y_continuous(labels=function(x){sprintf("%.2f", x)}, limits = c(-0.50, 0.50))+
  scale_x_continuous(labels=function(x){sprintf("%.2f", x)}, limits = c(-0.50, 0.50))

ggsave("plots/offspringtrans_deviation_plot_AlphaPeel.png", p_var_dev_off, width=25,
       height=25, units="cm")


#plot results for TD from mother by offspring sex 
Mum_off_sextrans_results_all_chr<-subset(Mum_off_sextrans_results_all_chr, !sex=="unknown")
p_Mum_off_sextrans<- manhattan_plot(df_plot = Mum_off_sextrans_results_all_chr, link_order = link_map_part, chr_no = 33, 
                                    p_threshold = 0.02, ylim_max=3, xlab_text = "CEL linkage group")

 # ggsave("Mumtrans_sex_split_manhattan_plot_phased_reclustered.png", p_Mum_off_sextrans, width = 25,
 #        height = 25, units = "cm")

# ggsave("Mumtrans_sex_split_manhattan_plot_phased_reclustered_sub.png", p_Mum_off_sextrans, width = 25,
#        height = 25, units = "cm")

ggsave("plots/Mumtrans_sex_split_manhattan_plot_AlphaPeel.png", p_Mum_off_sextrans, width = 25,
       height = 25, units = "cm")

#ggsave("Mumtrans_sex_split_manhattan_plot_phased.png", p_Mum_off_sextrans, width = 25,
       #height = 25, units = "cm")


# write.table(Mum_off_sextrans_results_all_chr, file="Mumtrans_offsex_TD_results_phased_reclustered.txt", sep="\t",
#             col.names = T, row.names = F, quote=F)

# write.table(Mum_off_sextrans_df.plot, file="Mumtrans_offsex_TD_results_phased_reclustered_sub.txt", sep="\t", 
#             col.names = T, row.names = F, quote=F)

 write.table(Mum_off_sextrans_results_all_chr, file="results/Mumtrans_offsex_TD_results_AlphaPeel.txt", sep="\t",
             col.names = T, row.names = F, quote=F) 
 
# write.table(Mum_off_sextrans_results_all_chr, file="Mumtrans_offsex_TD_results_phased.txt", sep="\t", 
#             col.names = T, row.names = F, quote=F)

Mum_off_sextrans_df_top<-subset(Mum_off_sextrans_results_all_chr, p_value<0.02)

write.table(Mum_off_sextrans_df_top, file="results/Mumtrans_offsex_TD_top_SNPs_AlphaPeel.txt",
            sep="\t", col.names = T, row.names = F, quote=F)

#plot results for TD from mother 
Mumtrans_results_all_chr$trans_parent<-NULL
p_Mumtrans<- manhattan_plot(df_plot = Mumtrans_results_all_chr, link_order = link_map_part, chr_no = 33, 
                                      p_threshold = 0.02, ylim_max = 3,xlab_text = "CEL linkage group")

# ggsave("Mumtrans_manhattan_plot_phased_reclustered.png", p_Mumtrans, width = 25,
#        height = 25, units = "cm")

# ggsave("Mumtrans_manhattan_plot_phased_reclustered_sub.png", p_Mumtrans, width = 25,
#        height = 25, units = "cm")

ggsave("plots/Mumtrans_manhattan_plot_AlphaPeel.png", p_Mumtrans, width = 25,
       height = 25, units = "cm")

 # ggsave("Mumtrans_manhattan_plot_phased.png", p_Mumtrans, width = 25,
 #        height = 25, units = "cm")

# write.table(Mumtrans_results_all_chr, file="Mumtrans_TD_results_phased_reclustered.txt", sep="\t",
#             col.names = T, row.names = F, quote=F)

# write.table(Mumtrans_df.plot, file="Mumtrans_TD_results_phased_reclustered_sub.txt", sep="\t", 
#             col.names = T, row.names = F, quote=F)

write.table(Mumtrans_results_all_chr, file="results/Mumtrans_TD_results_AlphaPeel.txt", sep="\t",
            col.names = T, row.names = F, quote=F)

 # write.table(Mumtrans_df.plot, file="Mumtrans_TD_results_phased.txt", sep="\t", 
 #            col.names = T, row.names = F, quote=F)



#filter for top SNPs (arbitrary at -log(p_value)f>3 or now) -reclustered phased data
Mumtrans_df_top<-subset(Mumtrans_results_all_chr, p_value<0.02)

write.table(Mumtrans_df_top, file="results/Mumtrans_TD_top_SNPs_AlphaPeel.txt", sep="\t",
            col.names = T, row.names = F, quote=F)


#plot results for TD from father by offspring sex 
Dad_off_sextrans_results_all_chr<-subset(Dad_off_sextrans_results_all_chr, !sex=="unknown")

p_Dad_off_sextrans<-manhattan_plot(df_plot = Dad_off_sextrans_results_all_chr, link_order = link_map_part, chr_no = 33, 
                                   p_threshold = 0.02, ylim_max = 3, xlab_text = "CEL linkage group")
  
  
# ggsave("Dadtrans_sex_split_manhattan_plot_phased_reclustered.png", p_Dad_off_sextrans, width = 25,
#       height = 25, units = "cm")

# ggsave("Dadtrans_sex_split_manhattan_plot_phased_reclustered_sub.png", p_Dad_off_sextrans, width = 25,
#        height = 25, units = "cm")

ggsave("plots/Dadtrans_sex_split_manhattan_plot_AlphaPeel.png", p_Dad_off_sextrans, width = 25,
              height = 25, units = "cm")

 
 # ggsave("Dadtrans_sex_split_manhattan_plot_phased.png", p_Dad_off_sextrans, width = 25,
 #        height = 25, units = "cm")


# write.table(Dad_off_sextrans_results_all_chr, file="Dadtrans_offsex_TD_results_phased_reclustered.txt", sep="\t",
#             col.names = T, row.names = F, quote=F)

# write.table(Dad_off_sextrans_df.plot, file="Dadtrans_offsex_TD_results_phased_reclustered_sub.txt", sep="\t", 
#             col.names = T, row.names = F, quote=F)

write.table(Dad_off_sextrans_results_all_chr, file="results/Dadtrans_offsex_TD_results_AlphaPeel.txt", sep="\t",
            col.names = T, row.names = F, quote=F)

# write.table(Dad_off_sextrans_df.plot, file="Dadtrans_offsex_TD_results_phased.txt", sep="\t", 
#             col.names = T, row.names = F, quote=F)

Dad_off_sextrans_df_top<-subset(Dad_off_sextrans_results_all_chr, p_value<0.02)

write.table(Dad_off_sextrans_df_top, file="results/Dadtrans_offsex_TD_top_SNPs_AlphaPeel.txt", sep="\t",
            col.names = T, row.names = F, quote=F)

#plot results for TD from father 
Dadtrans_results_all_chr$trans_parent<-NULL
p_Dadtrans<- manhattan_plot(df_plot = Dadtrans_results_all_chr, link_order = link_map_part, chr_no = 33, 
                            p_threshold = 0.02, ylim_max = 3, xlab_text = "CEL linkage group")

# ggsave("Dadtrans_manhattan_plot_phased_reclustered.png", p_Dadtrans, width = 25,
#       height = 25, units = "cm")

# ggsave("Dadtrans_manhattan_plot_phased_reclustered_sub.png", p_Dadtrans, width = 25,
#        height = 25, units = "cm")

ggsave("plots/Dadtrans_manhattan_plot_AlphaPhase.png", p_Dadtrans, width = 25,
       height = 25, units = "cm")
 
 # ggsave("Dadtrans_manhattan_plot_phased.png", p_Dadtrans, width = 25,
 #        height = 25, units = "cm")
 # 
# write.table(Dadtrans_results_all_chr, file="Dadtrans_TD_results_phased_reclustered.txt", sep="\t",
#             col.names = T, row.names = F, quote=F)

# write.table(Dadtrans_df.plot, file="Dadtrans_TD_results_phased_reclustered_sub.txt", sep="\t", 
#             col.names = T, row.names = F, quote=F)

write.table(Dadtrans_results_all_chr, file="results/Dadtrans_TD_results_AlphaPeel.txt", sep="\t",
            col.names = T, row.names = F, quote=F)

 # write.table(Dadtrans_df.plot, file="Dadtrans_TD_results_phased.txt", sep="\t", 
 #             col.names = T, row.names = F, quote=F)
 # 

Dadtrans_df_top<-subset(Dadtrans_results_all_chr, p_value<0.02)

write.table(Dadtrans_df_top, file ="results/Dadtrans_TD_top_SNPs_AlphaPeel.txt", sep="\t",
            col.names = T, row.names = F, quote = F)


##Mum and Dad trans (not by offspring sex) together#
# Mumtrans_results_all_chr<-read.table("Mumtrans_TD_results_phased_reclustered.txt", header=T)
# Dadtrans_results_all_chr<-read.table("Dadtrans_TD_results_phased_reclustered.txt", header=T)

Mumtrans_results_all_chr<-read.table("results/Mumtrans_TD_results_AlphaPeel.txt", header=T)
Dadtrans_results_all_chr<-read.table("results/Dadtrans_TD_results_AlphaPeel.txt", header=T)

#join Mum and Dad trans tables to make wrapped plot
Mumtrans_results_all_chr<-mutate(Mumtrans_results_all_chr, trans_parent="mother")
Dadtrans_results_all_chr<-mutate(Dadtrans_results_all_chr, trans_parent="father")

Mum_Dad_trans_results<-rbind(Mumtrans_results_all_chr, Dadtrans_results_all_chr)

p_Mum_Dadtrans<-manhattan_plot(df_plot = Mum_Dad_trans_results, link_order = link_map_part, chr_no = 33, 
               p_threshold = 0.02, ylim_max = 3, xlab_text = "CEL linkage group")



#ggsave("Mum_Dadtrans_manhattan_plot_phased.png", p_Mum_Dadtrans, width=25, height=25, units = "cm")

#ggsave("Mum_Dadtrans_manhattan_plot_phased_reclustered.png", p_Mum_Dadtrans, width=25, height=25, units = "cm")

ggsave("plots/Mum_Dadtrans_manhattan_plot_AlphaPeel.png", p_Mum_Dadtrans, width=25, height=25, units = "cm")


#deviation plot

#exclude SNPs with extremly low allele counts
freq.cutoff<-data.frame(freq.quantile=as.numeric(quantile(Mum_Dad_trans_results$expected_freq, probs = 0.05)))
Mum_Dad_trans_results<-subset(Mum_Dad_trans_results, expected_freq > freq.cutoff$freq.quantile)

#exlude SNPs where allele counts are only available for one parent sex 
#(want to plot male vs female so need information for both parents per SNP)
Mum_Dad_trans_results_sub<-subset(Mum_Dad_trans_results, trans_parent=="mother")
Mum_Dad_trans_results_sub2<-subset(Mum_Dad_trans_results, trans_parent=="father")

SNPs_uneven_sex<-subset(Mum_Dad_trans_results_sub2, !SNP.name %in% Mum_Dad_trans_results_sub$SNP.name)
Mum_Dad_trans_results<-subset(Mum_Dad_trans_results, !SNP.name %in% SNPs_uneven_sex$SNP.name)
Mum_Dad_trans_results<-droplevels(Mum_Dad_trans_results)

SNPs_uneven_sex2<-subset(Mum_Dad_trans_results_sub, !SNP.name %in% Mum_Dad_trans_results_sub2$SNP.name)
Mum_Dad_trans_results<-subset(Mum_Dad_trans_results, !SNP.name %in% SNPs_uneven_sex2$SNP.name)
Mum_Dad_trans_results<-droplevels(Mum_Dad_trans_results)

#add signifcance status to top SNPs
Mum_Dad_trans_results[["significance"]]<-ifelse(Mum_Dad_trans_results[["p_value"]]<0.02 , "sig", "non.sig")
table(Mum_Dad_trans_results$significance)

Mum_Dad_trans_results[["Deviation"]] <- (Mum_Dad_trans_results[["A"]]- Mum_Dad_trans_results$expected_freq)/Mum_Dad_trans_results$expected_freq

Mum_Dad_trans_results_deviations <- reshape2::dcast(Mum_Dad_trans_results[,c("SNP.name", "trans_parent", "Deviation")], formula = SNP.name ~ trans_parent)
Mum_Dad_trans_results_status <- reshape2::dcast(Mum_Dad_trans_results[,c("SNP.name", "trans_parent", "significance")], formula = SNP.name ~ trans_parent)

names(Mum_Dad_trans_results_deviations) <- c("SNP", "Male.Deviation", "Female.Deviation")
names(Mum_Dad_trans_results_status) <- c("SNP",  "Male.Status", "Female.Status")


Mum_Dad_trans_results_dev_plot<-join(Mum_Dad_trans_results_deviations, Mum_Dad_trans_results_status)

Mum_Dad_trans_results_dev_plot[["Label"]]<-ifelse(Mum_Dad_trans_results_dev_plot[["Female.Status"]]=="sig" & Mum_Dad_trans_results_dev_plot[["Male.Status"]]=="sig", "sig",
                                                   ifelse(Mum_Dad_trans_results_dev_plot[["Female.Status"]]=="sig" &
                                                            Mum_Dad_trans_results_dev_plot[["Male.Status"]]=="non.sig", "sig.female",
                                                          ifelse(Mum_Dad_trans_results_dev_plot[["Male.Status"]]=="sig" &
                                                                   Mum_Dad_trans_results_dev_plot[["Female.Status"]]=="non.sig", "sig.male", "non.sig")))
table(Mum_Dad_trans_results_dev_plot$Label)

Mum_Dad_trans_results_dev_plot[["Alpha"]] <- ifelse(Mum_Dad_trans_results_dev_plot$Label == "non.Sig.", 0.1, 0.7)


#what are the minimum and maximum deviations of significant SNPs for males and females?
Mum_Dad_trans_deviations_sig_fem<-subset(Mum_Dad_trans_results_dev_plot, Label == "sig" | Label == "sig.female")
min.mothers<-min(abs(Mum_Dad_trans_deviations_sig_fem$Female.Deviation)) #-0.25984 or absolute 0.10853
max.mothers<-max(abs(Mum_Dad_trans_deviations_sig_fem$Female.Deviation)) #0.1700 or absolute 0.2598

Mum_Dad_trans_deviations_sig_male<-subset(Mum_Dad_trans_results_dev_plot, Label == "sig" | Label == "sig.male")  
min.fathers<-min(abs(Mum_Dad_trans_deviations_sig_male$Male.Deviation)) #-0.25833 or absolute 0.11416
max.fathers<-max(abs(Mum_Dad_trans_deviations_sig_male$Male.Deviation)) #0.33846

#save Mum_Dad_trans TD table with deviation information
write.table(Mum_Dad_trans_results, file = "results/Mum_Dad_trans_TD_results_AlphaPeel.txt", sep="\t", 
            col.names = T, row.names = F, quote=F)


p_var_dev_mumdad<-ggplot(Mum_Dad_trans_results_dev_plot, aes(Female.Deviation, Male.Deviation, colour = Label, alpha = Alpha)) +
  geom_point(size=2) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Female Deviation", y = "Male Deviation") +
  scale_colour_manual(values = c("grey70", "red", "blue"), labels=c("Not significant", "Females only", "Males only"))+
  guides(alpha=FALSE)+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=16, color = "black"), legend.background = element_rect(size=0.5, 
                                                                                               linetype = "solid", color="black"),
        axis.text = element_text(size=22, color="black"), axis.title = element_text(size=26),
        axis.ticks=element_line(size=1))+
  scale_y_continuous(labels=function(x){sprintf("%.2f", x)}, limits = c(-0.50, 0.50))+
  scale_x_continuous(labels=function(x){sprintf("%.2f", x)}, limits = c(-0.50, 0.50))

ggsave("plots/MumDadtrans_deviation_plot_AlphaPeel.png", p_var_dev_mumdad, width=25,
       height=25, units="cm")


#make summary table of significant SNP TD deviation values
#all values based on absolute deviations, since sign is not important (randomly allocated allele label)
TD_deviations_sig<- data.frame(Status=rep(c("parent", "offspring"), 2),
                               Sex=rep(c("female", "male"), each=2), min=c(min.mothers, min.female, min.fathers, min.male),
                               max=c(max.mothers, max.female, max.fathers, max.male),
                               mean=c(mean(abs(Mum_Dad_trans_deviations_sig_fem$Female.Deviation)),
                                      mean(abs(offspringtrans_deviations_sig_fem$Female.Deviation)),
                                      mean(abs(Mum_Dad_trans_deviations_sig_male$Male.Deviation)),
                                      mean(abs(offspringtrans_deviations_sig_male$Male.Deviation))),
                               std=c(sd(abs(Mum_Dad_trans_deviations_sig_fem$Female.Deviation)),
                                      sd(abs(offspringtrans_deviations_sig_fem$Female.Deviation)),
                                      sd(abs(Mum_Dad_trans_deviations_sig_male$Male.Deviation)),
                                      sd(abs(offspringtrans_deviations_sig_male$Male.Deviation))))



print(xtable(TD_deviations_sig, type = "latex", 
             digits= c(1,1,1,4,4,4,4) ),
             file = "results/TD_significant_SNPs_effect_summary_tbl.tex")


#combine all TD data, including deviations
Mum_Dad_trans_results_all_chr<-read.table("results/Mum_Dad_trans_TD_results_AlphaPeel.txt", header=T)
offspringtrans_results_all_chr<-read.table("results/offspringtrans_TD_results_AlphaPeel.txt", header=T)

head(Mum_Dad_trans_results_all_chr)
head(offspringtrans_results_all_chr)

names(Mum_Dad_trans_results_all_chr)[names(Mum_Dad_trans_results_all_chr) == "trans_parent"]<-"sex"
Mum_Dad_trans_results_all_chr$sex<-ifelse(Mum_Dad_trans_results_all_chr$sex == "mother", "female", "male")
Mum_Dad_trans_results_all_chr$status<-"parent"
Mum_Dad_trans_results_all_chr<-Mum_Dad_trans_results_all_chr[, c("SNP.name", "sex", "status", "A", "B","Chromosome", "minor.allele",
                                                                 "p_value", "expected_freq", "Deviation")]

offspringtrans_results_all_chr$status<-"offspring"
offspringtrans_results_all_chr<-offspringtrans_results_all_chr[, c("SNP.name", "sex", "status", "A", "B","Chromosome", "minor.allele",
                                                                   "p_value", "expected_freq", "Deviation")]

TD_results_all<-rbind(Mum_Dad_trans_results_all_chr, offspringtrans_results_all_chr)

#add information about cow position and order of SNPs, plus info on MAF

#read in linkage map with new cattle genome positions 
linkmap_cow_pos<-read.table("data/Cervus_elaphus_linkage_map_additional_BT_pos.txt", header=T)
head(linkmap_cow_pos)
names(linkmap_cow_pos)[names(linkmap_cow_pos)=="SNP.Name"]<-"SNP.name"

TD_results_all<-join(TD_results_all, linkmap_cow_pos[, c("SNP.name", "CEL.order", "cow_ARS_UCD1.2_Chr",
                                                                   "cow_ARS_UCD1.2_pos")])
#read in MAF info table
MAF_info<-read.table("results/Compare_MAF_pre_post_phasing_AlphaPeel.txt", header=T)
head(MAF_info)

TD_results_all<-join(TD_results_all, MAF_info[, c("SNP.name", "MAF_unphased")])
names(TD_results_all)[names(TD_results_all) == "Chromosome"]<-"CEL.LG"
TD_results_all<-TD_results_all[, c("SNP.name", "sex", "status", "A", "B","expected_freq", "minor.allele", "Deviation",
                                            "p_value", "CEL.LG", "CEL.order", "cow_ARS_UCD1.2_Chr",
                                             "cow_ARS_UCD1.2_pos", "MAF_unphased")]

TD_results_all<-arrange(TD_results_all, CEL.LG, CEL.order)

write.table(TD_results_all, file="results/TD-all_results_info_complete.txt", sep="\t", 
            row.names = F, col.names = T, quote = F)


#make manhattan plot split by sex and status (offspring/parent)
TD_results_all<-read.table("results/TD-all_results_info_complete.txt", header = T)

TD_results_all[ ,"CEL.order"]<-as.numeric(TD_results_all[ ,"CEL.order"])
TD_results_all[ ,"CEL.LG"]<-as.numeric(TD_results_all[ ,"CEL.LG"])
TD_results_all<-arrange(TD_results_all, CEL.LG, CEL.order)
cum<-nrow(TD_results_all)
TD_results_all[ ,"pos_cuml"]<-c(1:cum)

#make separate x axis for chromosomes and add title variable to label plot with a facet title including an expression
axisdf <- plyr::ddply(TD_results_all,.(CEL.LG), summarize, center=(max(pos_cuml) + min(pos_cuml))/2)

label_seq<-axisdf$CEL.LG[seq(1,length(axisdf$CEL.LG),by=2)]
breaks_seq<-axisdf$center[seq_along(axisdf$center)%% 2 > 0]

TD_results_all<-subset(TD_results_all, !sex=="unknown")

#make custom facet labels
sex.labs<-c("Females", "Males")
names(sex.labs)<-c("female", "male")
status.labs<-c("Parent", "Offspring")
names(status.labs)<-c("parent", "offspring")

chr_no = 33
p_threshold = 0.02
xlab_text = "CEL Linkage Group"
ylim_max= 3.5
x_text_size = 16 
y_text_size = 16
point_size = 3
alpha_value = 0.5
chr_colours = c("steelblue3", "red3")

  p_manhattan<-ggplot(TD_results_all, aes(pos_cuml, -log10(p_value))) +
    geom_point(aes(color=as.factor(CEL.LG)), size=point_size, alpha=alpha_value)+
    scale_colour_manual(values = rep(chr_colours, chr_no))+
    scale_x_continuous(label= label_seq, breaks = breaks_seq)+ 
    scale_y_continuous(expand = c(0,0))+
    ylim(0,ylim_max)+ 
    geom_hline(yintercept =  -log10(p_threshold), linetype="dashed", color="black", size=1)+
    xlab(xlab_text)+ylab(expression(-log[10](p)))+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(size=x_text_size, colour = "black"),
          axis.text.y = element_text(size=y_text_size, colour = "black"), 
          axis.title.x = element_text(size=30, margin = margin(t=20)), 
          axis.title.y = element_text(size=30, margin = margin(r=20)),
          axis.line = element_line(colour = "black"),
          axis.ticks=element_line(size=1), strip.text = element_text(size=22), #, hjust=0.01
          strip.background = element_rect(colour="black", fill="lightgray"))+
    facet_grid(rows = vars(status), cols=vars(sex), labeller = labeller(sex=sex.labs, status=status.labs))

  
  
  ggsave("plots/TD_manhattan_plot_all_AlphaPeel.png", p_manhattan, width=35,
         height=25, units="cm")  
  
  
#make wrapped deviation plots
#join data frames for parent and offspring deviatio plots

Mum_Dad_trans_results_dev_plot$status<-"parent"
offspringtrans_results_dev_plot$status<-"offspring"  

TD_results_all_dev_plot<-rbind(Mum_Dad_trans_results_dev_plot, offspringtrans_results_dev_plot)
  
  p_dev<-ggplot(TD_results_all_dev_plot, aes(Female.Deviation, Male.Deviation, colour = Label, alpha = Alpha)) +
    geom_point(size=2) +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "Female Deviation", y = "Male Deviation") +
    scale_colour_manual(values = c("grey70", "red", "blue"), labels=c("Not significant", "Females only", "Males only"))+
    guides(alpha=FALSE)+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=16, color = "black"), legend.background = element_rect(size=0.5, 
                                                                 linetype = "solid", color="black"),
          axis.text = element_text(size=22, color="black"), axis.title = element_text(size=30),
          axis.ticks=element_line(size=1), strip.text = element_text(size=22), #, hjust=0.01
          strip.background = element_rect(colour="black", fill="lightgray"))+
    scale_y_continuous(labels=function(x){sprintf("%.2f", x)}, limits = c(-0.50, 0.50))+
    scale_x_continuous(labels=function(x){sprintf("%.2f", x)}, limits = c(-0.50, 0.50))+
    facet_wrap(vars(status),  nrow = 2, ncol = 1, labeller = labeller(status=status.labs))
  
ggsave("plots/TD_deviation_plot_all_AlphaPeel.png", p_dev, width=25,
       height=25, units="cm")   
  
#combine significant TD results in one table
Mumtrans_SNPs_sig<-read.table("results/Mumtrans_TD_top_SNPs_AlphaPeel.txt", header=T, stringsAsFactors = F)
Dadtrans_SNPs_sig<-read.table("results/Dadtrans_TD_top_SNPs_AlphaPeel.txt", header=T, stringsAsFactors = F)
offspringtrans_SNPs_sig<-read.table("results/offspringtrans_TD_top_SNPs_AlphaPeel.txt", header=T, stringsAsFactors = F)

#add sex and status (parent or offspring) information where needed
Mumtrans_SNPs_sig<-Mumtrans_SNPs_sig %>%
  mutate(sex="female")

Mumtrans_SNPs_sig<-Mumtrans_SNPs_sig %>%
  mutate(status="parent")

Mumtrans_SNPs_sig<-Mumtrans_SNPs_sig[, c("SNP.name", "sex", "status", "A", "B", "Chromosome", 
                                         "minor.allele", "error.count", "p_value", "expected_freq")]

Dadtrans_SNPs_sig<-Dadtrans_SNPs_sig %>%
  mutate(sex="male")

Dadtrans_SNPs_sig<-Dadtrans_SNPs_sig %>%
  mutate(status="parent")

Dadtrans_SNPs_sig<-Dadtrans_SNPs_sig[, c("SNP.name", "sex", "status", "A", "B", "Chromosome", 
                                         "minor.allele", "error.count", "p_value", "expected_freq")]

offspringtrans_SNPs_sig<-offspringtrans_SNPs_sig %>%
  mutate(status="offspring")

offspringtrans_SNPs_sig<-offspringtrans_SNPs_sig[, c("SNP.name", "sex", "status", "A", "B", "Chromosome", 
                                                     "minor.allele", "error.count", "p_value", "expected_freq")]

TD_significant_SNPs<-rbind(offspringtrans_SNPs_sig, Mumtrans_SNPs_sig)
TD_significant_SNPs<-rbind(TD_significant_SNPs, Dadtrans_SNPs_sig)
names(TD_significant_SNPs)[names(TD_significant_SNPs) == "Chromosome"]<- "CEL.LG"

#add information about cow position and order of SNPs, plus info on MAF

#read in linkage map with new cattle genome positions 
linkmap_cow_pos<-read.table("data/Cervus_elaphus_linkage_map_additional_BT_pos.txt", header=T)
head(linkmap_cow_pos)
names(linkmap_cow_pos)[names(linkmap_cow_pos)=="SNP.Name"]<-"SNP.name"

TD_significant_SNPs<-join(TD_significant_SNPs, linkmap_cow_pos[, c("SNP.name", "CEL.order", "cow_ARS_UCD1.2_Chr",
                                                                   "cow_ARS_UCD1.2_pos")])
#read in MAF info table
MAF_info<-read.table("results/Compare_MAF_pre_post_phasing_AlphaPeel.txt", header=T)
head(MAF_info)

TD_significant_SNPs<-join(TD_significant_SNPs, MAF_info[, c("SNP.name", "MAF_unphased")])

TD_significant_SNPs<-TD_significant_SNPs[, c("SNP.name", "sex", "status", "A", "B", "minor.allele", "expected_freq",
                                             "error.count", "p_value", "CEL.LG", "CEL.order", "cow_ARS_UCD1.2_Chr",
                                             "cow_ARS_UCD1.2_pos", "MAF_unphased")]

TD_significant_SNPs<-arrange(TD_significant_SNPs, CEL.LG, CEL.order)
rownames(TD_significant_SNPs)<-NULL

write.table(TD_significant_SNPs, file = "results/TD_significant_SNPs_parent.off_info_tbl.txt", sep = "\t",
            col.names = T, row.names = F, quote = F)

print(xtable(TD_significant_SNPs[ , -8], type = "latex", 
             display = c("d", "s", "s", "s", "fg", "fg", "s", "f", "fg", "d", "d", "d", "d", "fg")), 
      file = "results/TD_significant_SNPs_parent.off_info_tbl.tex")












#compare LINKPHASE and AlphaPeel results

offspringtrans_results_AlphaPeel<-read.table("results/offspringtrans_TD_results_AlphaPeel.txt", header=T)
offspringtrans_results_LINK<-read.table("results/offspringtrans_TD_results_phased_reclustered.txt", header = T)

#exclude SNPs missing in one of two tables
offspringtrans_results_AlphaPeel<-subset(offspringtrans_results_AlphaPeel, SNP.name %in% offspringtrans_results_LINK$SNP.name)
offspringtrans_results_LINK<-subset(offspringtrans_results_LINK, SNP.name %in% offspringtrans_results_AlphaPeel$SNP.name)

names(offspringtrans_results_AlphaPeel)[names(offspringtrans_results_AlphaPeel) == "expected_freq"] <- "AlphaPeel_efreq"
names(offspringtrans_results_AlphaPeel)[names(offspringtrans_results_AlphaPeel) == "p_value"] <- "AlphaPeel_p_value"

names(offspringtrans_results_LINK)[names(offspringtrans_results_LINK) == "expected_freq"] <- "LINKPHASE_efreq"
names(offspringtrans_results_LINK)[names(offspringtrans_results_LINK) == "p_value"] <- "LINKPHASE_p_value"

offspringtrans_results_AlphaPeel<-subset(offspringtrans_results_AlphaPeel, sex=="female")
offspringtrans_results_LINK<-subset(offspringtrans_results_LINK, sex=="female")

offspringtrans_results_phase_compare<-join(offspringtrans_results_AlphaPeel[, c("SNP.name", "AlphaPeel_efreq", "AlphaPeel_p_value")], 
                                           offspringtrans_results_LINK[, c("SNP.name", "LINKPHASE_efreq", "LINKPHASE_p_value" )])

p_phase_freq_compare<-ggplot(offspringtrans_results_phase_compare, aes(AlphaPeel_efreq, LINKPHASE_efreq))+geom_point()+stat_smooth(method = "lm")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size=25),
        axis.line = element_line(colour="black"))+
  xlab("expected count AlphaPeel")+ylab("expected count LINKPHASE")

ggsave("Compare_phasing_sample_size.png", p_phase_freq_compare, width=25,
       height=25, units="cm")

cor.test(offspringtrans_results_phase_compare$AlphaPeel_efreq, offspringtrans_results_phase_compare$LINKPHASE_efreq)

ggplot(offspringtrans_results_phase_compare, aes(AlphaPeel_p_value, LINKPHASE_p_value))+geom_point()+stat_smooth(method = "lm")+
  theme(axis.text = element_text(size=20, colour="black"), axis.title = element_text(size=25),
        axis.line = element_line(colour="black"))+
  xlab("p value AlphaPeel")+ylab("p value LINKPHASE")

#subset according to top AlphaPeel SNPs
#AlphaPeel p-value
#p_threshold<-as.numeric(quantile(offspringtrans_results_phase_compare$AlphaPeel_p_value, probs = 0.001))

offspringtrans_results_phase_compare_sub<-subset(offspringtrans_results_phase_compare, 
                                                 AlphaPeel_p_value<0.01)

p_phase_freq_compare_sub<-ggplot(offspringtrans_results_phase_compare_sub, aes(AlphaPeel_efreq, LINKPHASE_efreq))+geom_point()+stat_smooth(method = "lm")+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))+
  xlab("expected count AlphaPeel")+ylab("expected count LINKPHASE")

ggsave("Compare_phasing_sample_size_top_AlphaPeel.png", p_phase_freq_compare_sub, width=25,
       height=25, units="cm")

p_phase_pv_compare_sub<-ggplot(offspringtrans_results_phase_compare_sub, aes(AlphaPeel_p_value, LINKPHASE_p_value))+geom_point()+stat_smooth(method = "lm")+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))+
  xlab("p value AlphaPeel")+ylab("p value LINKPHASE")

ggsave("Compare_phasing_p_value_top_AlphaPeel.png", p_phase_pv_compare_sub, width=25,
       height=25, units="cm")

#subset according to top LINKPHASE SNPs
offspringtrans_results_phase_compare_sub2<-subset(offspringtrans_results_phase_compare, 
                                                 -log10(LINKPHASE_p_value)>3)

p_phase_freq_compare_sub2<-ggplot(offspringtrans_results_phase_compare_sub2, aes(AlphaPeel_efreq, LINKPHASE_efreq))+geom_point()+stat_smooth(method = "lm")+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))+
  xlab("expected count AlphaPeel")+ylab("expected count LINKPHASE")

ggsave("Compare_phasing_sample_size_top_LINKPHASE.png", p_phase_freq_compare_sub2, width=25,
       height=25, units="cm")

p_phase_pv_compare_sub2<-ggplot(offspringtrans_results_phase_compare_sub2, aes(AlphaPeel_p_value, LINKPHASE_p_value))+geom_point()+stat_smooth(method = "lm")+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))+
  xlab("p value AlphaPeel")+ylab("p value LINKPHASE")

ggsave("Compare_phasing_p_value_top_LINKPHASE.png", p_phase_pv_compare_sub2, width=25,
       height=25, units="cm")





