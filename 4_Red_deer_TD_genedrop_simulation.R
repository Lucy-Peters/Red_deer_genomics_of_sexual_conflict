#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Red deer: Genedrop simulation of Mendelian inheritance patterns using founder haplotypes #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(GeneDrop)
library(parallel)
library(data.table)

##Below this is all run on a linux server (ashworth)##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define output directory (scratch)
args <- commandArgs(trailingOnly = TRUE)
resultpath <- as.character(args[1])
print(resultpath)

#read in phased genotype matrix
phase_all_chr<-read.table("data/Deer_AlphaPeel_out_all_chr.txt", header=F)
colnames(phase_all_chr)[1]<-"ID"

#read in SNP marker info (SNP name, order, linkage group)
marker_map_all_chr<-read.table("data/Deer_AlphaPeel_marker_file_all_chr.txt", header=T, stringsAsFactors=F)

#list of genotyped IDs
IDs_genotyped<-read.table("data/IDs_genotyped_reclustered_SNPs.txt", header=T, stringsAsFactors=F)

#pedigree
ped<-read.table("data/Pedigree_deer_2017_19_merged_genedrop.txt", header=T, stringsAsFactors=F)

#for GeneDrop, make males 1 and females 2, unknown NA
ped$Sex<-ifelse(ped$Sex==1, 2, 
                ifelse(ped$Sex==2, 1, NA))

#exclude individuals with unknown sex
ped<-subset(ped, !is.na(Sex))

colnames(ped)<-c("ID", "Dam", "Sire", "Sex", "Cohort")

#remove everything from before 1980 just now
ped <- subset(ped, Cohort > 1980)
ped$Sire[which(!ped$Sire %in% ped$ID)] <- 0
ped$Dam [which(!ped$Dam  %in% ped$ID)] <- 0

#remove singletons
#ped<- ped[-which(ped$Sire == 0 & ped$Dam == 0 & !ped$ID %in% ped$Sire & !ped$ID %in% ped$Dam),]

#identify founders
ped <- arrange(ped, Cohort)
ped$Founder <- ifelse(ped$Sire == 0 & ped$Dam == 0, "1", "2")
ped <- arrange(ped, Founder)
founders <- ped$ID[ped$Founder == 1] %>% data.frame
colnames(founders)<-"ID"
ped<-arrange(ped, Founder)
ped$Founder<-NULL

founders_unknown<-subset(founders, !ID %in% IDs_genotyped$ID)
founders<-subset(founders, !ID %in% founders_unknown$ID)

founders_vec<-founders$ID
founders_unk_vec<-founders_unknown$ID

ped_completed <-complete_ped_links(ped, founders_vec,
                   founders_unk = founders_unk_vec, rep_years_sire = c(5,15),
                   rep_years_dam = c(3, 12),
                   add_dummy_parents = TRUE,
                   limit_cohort_low = TRUE) %>% data.frame()

 
dummy_founders<-data.frame(ID=grep("^Dummy.*", ped_completed$ID, value=T))

#fix pedigree
ped_fix<-fix_pedigree(ped_completed)%>%data.frame()

#re-define known and unknown founders in same order as completed pedigree
founders_all_ordered<-subset(ped_fix, Sire==0 & Dam==0)
founders_all_ordered$order<-c(1:nrow(founders_all_ordered))

# #if dummy founders will have simulated haplotypes
# founders<-join(founders, founders_all_ordered, by="ID")
# founders<-arrange(founders, order)
# founders<-founders[, "ID"]%>% data.frame()
# colnames(founders)<-"ID"
# 
# founders_unknown<-subset(founders_all_ordered, !ID %in% IDs_genotyped$ID)
# founders_unknown<-join(founders_unknown, founders_all_ordered)
# founders_unknown<-arrange(founders_unknown, order)
# founders_unknown<-founders_unknown[, "ID"]%>% data.frame()
# colnames(founders_unknown)<-"ID"
# 

#if dummy parents will have empty haplotypes
founders<-rbind(founders, dummy_founders)
founders<-join(founders, founders_all_ordered, by="ID")
founders<-arrange(founders, order)
founders<-founders[, "ID"]%>% data.frame()
colnames(founders)<-"ID"


founders_unknown<-subset(founders_all_ordered, !ID %in% founders$ID)
founders_unknown<-join(founders_unknown, founders_all_ordered)
founders_unknown<-arrange(founders_unknown, order)
founders_unknown<-founders_unknown[, "ID"]%>% data.frame()
colnames(founders_unknown)<-"ID"


#make pedigree matrix
ped_fix<-as.matrix(ped_fix)


#linkage map for recombination frequency information
linkmap<-read.table("data/TableS1_CervusElaphus_Final_Linkage_Map.txt", header=T, stringsAsFactors=F)
linkmap<-linkmap[,c("SNP.Name", "CEL.LG", "CEL.order", "cMPosition.SexAveraged")]
linkmap$order<-c(1:nrow(linkmap))

founder_haplos<-subset(phase_all_chr, ID %in% founders$ID)

#add empty haplotypes for dummy founders

#repeat IDs
dummy_haplos_IDs<-data.frame(ID=rep(dummy_founders$ID, 2))
dummy_haplos_IDs<-arrange(dummy_haplos_IDs, ID)

#create empty (all 9) haplotypes for dummy IDs
dummy_haplos<-phase_all_chr[-1][c(1:nrow(dummy_haplos_IDs)), ]
dummy_haplos[, ]<-sapply(dummy_haplos[, ], as.character)
dummy_haplos [dummy_haplos == "0"]<-"9"
dummy_haplos [dummy_haplos == "1"]<-"9"

dummy_haplos<-cbind(dummy_haplos_IDs, dummy_haplos)

#add dummy haplotypes to founder haplotypes
founder_haplos<-rbind(founder_haplos, dummy_haplos)

founders<-founders%>%
  mutate(order=c(1:nrow(founders)))

founder_haplos<-join(founder_haplos, founders, by="ID")
founder_haplos<-arrange(founder_haplos, order)
founder_haplos$order<-NULL
founder_haplos$ID<-NULL

founder_haplos[, ]<-sapply(founder_haplos[, ], as.numeric)

founder_haplos_ma<-as.matrix(founder_haplos)
dim(founder_haplos_ma)

#get number of SNPs and map distances
linkmap_sub<-subset(linkmap, SNP.Name %in% marker_map_all_chr$SNP.Name)
linkmap_sub<-arrange(linkmap_sub, order)

#map distances
cMdiff<-diff(linkmap_sub$cMPosition.SexAveraged)
cMdiff<-append(cMdiff, 0)
cMdiff[cMdiff < 0]<- 0

#no of SNPs
loci_no<-linkmap_sub %>%
  group_by(CEL.LG)%>%
  summarise(no.loci=n())
loci_no<-loci_no$no.loci

#founders as vectors
founders_unknown_vec<-founders_unknown$ID

#Genedrop run

#make genedrop function to run it i times
genedrop_sim<- function(i){
  
  genedrop(pedigree = ped_fix, map_dist = cMdiff, chr_loci_num = loci_no, found_hap = founder_haplos_ma,
           founders_unk = founders_unknown_vec,  to_raw = TRUE, sample_hap = FALSE)
  
}

genedrop_list<-mclapply(1:25, genedrop_sim, mc.cores = 10)

save(genedrop_list, file=paste0(resultpath,"Deer_genedrop_objects_sim_", run, ".RData"))

#make MAF per cohort tables
source("scripts/MAF_per_cohort_function.R")

MAF_per_cohort_all_sim<-mclapply(genedrop_list, MAF_per_cohort, marker_map = marker_map_all_chr, 
                                 IDs_keep = IDs_genotyped , mc.cores = 10)

save(MAF_per_cohort_all_sim, file=paste0(resultpath, "Genedrop_sim_MAF_per_cohort_", run, ".RData"))
   
#make allele transmission tables 
source("scripts/TD_tables_phased_by_chr_function.R")

print("Start making TD tables.")
    
TD_table_objects_list<-mclapply(genedrop_list, TD_tables_phased_by_chr, marker_map = marker_map_all_chr, 
                                IDs_keep = IDs_genotyped, chr_range=c(1:33), mc.cores=10)   

print("Done making TD tables.")

save(TD_table_objects_list, file=paste0(resultpath, "Genedrop_sim_TD_tables", run, ".RData"))

#combine all TD tables from simulations by type (e.g Mumtrans, Dadtrans, offspringtrans...)

Mumtrans_all_sim<-data.frame()
Mum_off_sextrans_all_sim<-data.frame()
Dadtrans_all_sim<-data.frame()
Dad_off_sextrans_all_sim<-data.frame()
offspringtrans_all_sim<-data.frame()   
SNP_error_all_sim<-data.frame()

for (TD_tbl_obj in TD_table_objects_list) {
  
  Mumtrans_part<-TD_tbl_obj@Mumtrans_tbl
  Mumtrans_all_sim<-rbind(Mumtrans_all_sim, Mumtrans_part)
  
  Mum_off_sextrans_part<-TD_tbl_obj@Mum_off_sextrans_tbl
  Mum_off_sextrans_all_sim<-rbind(Mum_off_sextrans_all_sim, Mum_off_sextrans_part)
  
  Dadtrans_part<-TD_tbl_obj@Dadtrans_tbl
  Dadtrans_all_sim<-rbind(Dadtrans_all_sim, Dadtrans_part)
  
  Dad_off_sextrans_part<-TD_tbl_obj@Dad_off_sextrans_tbl
  Dad_off_sextrans_all_sim<-rbind(Dad_off_sextrans_all_sim, Dad_off_sextrans_part)
  
  offspringtrans_part<-TD_tbl_obj@offspringtrans_tbl
  offspringtrans_all_sim<-rbind(offspringtrans_all_sim, offspringtrans_part)
  
  SNP_error_part<-TD_tbl_obj@SNP_error_tbl
  SNP_error_all_sim<-rbind(SNP_error_all_sim, SNP_error_part)
  
  
}

  
#Fisher's exact test
source("scripts/TD_Fishers_test.R")

print("Starting Fisher's exact tests on all simulations ...")

#TD from mother - disregarding offspring sex
Mumtrans_results_all<-TD_Fishers_test(Mumtrans_all_sim)

#TD from mother - split by offspring sex
Mum_off_sextrans_results_all<-TD_Fishers_test(Mum_off_sextrans_all_sim)

#TD from father - disregarding offspring sex
Dadtrans_results_all<-TD_Fishers_test(Dadtrans_all_sim)

#TD from father - split by offspring sex
Dad_off_sextrans_results_all<-TD_Fishers_test(Dad_off_sextrans_all_sim)

#alleles transmitted by mother or father split by offspring sex
offspringtrans_results_all<-TD_Fishers_test(offspringtrans_all_sim)


write.table(Mumtrans_results_all, file=paste0(resultpath, "Mumtrans_genedrop_sim_TD_tables", run, ".txt"), sep = "\t",
            col.names = T, row.names = F, quote = F)

write.table(Mum_off_sextrans_results_all, file=paste0(resultpath, "Mum_off_sextrans_genedrop_sim_TD_tables", run, ".txt"), sep = "\t",
            col.names = T, row.names = F, quote = F)

write.table(Dadtrans_results_all, file=paste0(resultpath, "Dadtrans_genedrop_sim_TD_tables", run, ".txt"), sep = "\t",
            col.names = T, row.names = F, quote = F)

write.table(Dad_off_sextrans_results_all, file=paste0(resultpath, "Dad_off_sextrans_genedrop_sim_TD_tables", run, ".txt"), sep = "\t",
            col.names = T, row.names = F, quote = F)

write.table(offspringtrans_results_all, file=paste0(resultpath, "offspringtranstrans_genedrop_sim_TD_tables", run, ".txt"), sep = "\t",
            col.names = T, row.names = F, quote = F)

