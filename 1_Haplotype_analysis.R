##File prep for phasing haplotypes##

library(plyr)
library(dplyr)
library(tidyr)
library(kinship2)
library(reshape2)
library(RODBC)
library(ggplot2)

#pedigree file

ped<-read.csv("Pedigree_2019-07-02_DB.csv", header=T)
head(ped)
ped<-ped[, c(1:3)]
colnames(ped)[1]<-"ID"
ped_dup<-subset(ped, duplicated(ID))
ped<-subset(ped, !duplicated(ID))


#check new pedigree (2019) for inconsistencies with old pedigree
ped_old<-read.table("Pedigree_Deer_2017-12.txt", header=T)
head(ped_old)
colnames(ped_old)<-c("ID", "Mother_old", "Father_old")

ped_sub_comp<-subset(ped, ID %in% ped_old$ID)
ped_old_sub<-subset(ped_old, ID %in% ped_sub_comp$ID)

ped_old_new_comp<-join(ped_sub_comp, ped_old_sub, by="ID")
ped_old_new_comp$Mother<-as.character(ped_old_new_comp$Mother)
ped_old_new_comp$Mother_old<-as.character(ped_old_new_comp$Mother_old)
ped_old_new_comp$Father<-as.character(ped_old_new_comp$Father)
ped_old_new_comp$Father_old<-as.character(ped_old_new_comp$Father_old)

ped_old_new_comp[is.na(ped_old_new_comp)]<-"unknown"


#categorise pedigree inconsistencies - do mother and father assignments agree? Is the mother/father
#present in one version of the pedigree but not the other?

ped_old_new_comp$ped_error<-ifelse(ped_old_new_comp$Mother!=ped_old_new_comp$Mother_old, "mat_error",
                                   ifelse(ped_old_new_comp$Father!=ped_old_new_comp$Father_old, "pat_error", 
                                          ifelse(ped_old_new_comp$Mother!=ped_old_new_comp$Mother_old & 
                                                   ped_old_new_comp$Father!=ped_old_new_comp$Father_old, "mat_pat_error", "no_error")))


ped_old_new_comp$ped_error[which(ped_old_new_comp$Mother=="unknown" & ped_old_new_comp$Mother_old!="unknown")]<-"Mum_new_missing"                                                                   
ped_old_new_comp$ped_error[which(ped_old_new_comp$Father=="unknown" & ped_old_new_comp$Father_old!="unknown")]<-"Dad_new_missing"
ped_old_new_comp$ped_error[which(ped_old_new_comp$Father=="unknown" & ped_old_new_comp$Father_old!="unknown" & ped_old_new_comp$Mother=="unknown" & ped_old_new_comp$Mother_old!="unknown")]<-"Both_new_missing"

ped_old_new_comp$ped_error[which(ped_old_new_comp$Mother!="unknown" & ped_old_new_comp$Mother_old=="unknown")]<-"Mum_old_missing"
ped_old_new_comp$ped_error[which(ped_old_new_comp$Father!="unknown" & ped_old_new_comp$Father_old=="unknown")]<-"Dad_old_missing"
ped_old_new_comp$ped_error[which(ped_old_new_comp$Father!="unknown" & ped_old_new_comp$Father_old=="unknown" & ped_old_new_comp$Mother!="unknown" & ped_old_new_comp$Mother_old=="unknown")]<-"Both_old_missing"

ped_old_new_comp$ped_error<-as.factor(ped_old_new_comp$ped_error)              
ped_old_new_comp$ped_error<-factor(ped_old_new_comp$ped_error, 
                                   levels = c("no_error","mat_error", "Mum_new_missing", "Mum_old_missing", "pat_error",
                                              "Dad_new_missing","Dad_old_missing", "Both_new_missing", "Both_old_missing"))
ped_old_new_comp_error<-subset(ped_old_new_comp, !ped_error=="no_error")
library(ggplot2)
 
ggplot(ped_old_new_comp_error, aes(ped_error))+geom_bar(position = "dodge")+
  scale_x_discrete(labels=c("maternity inconsistency", "missing dam 2019", "missing dam 2017", "paternity inconsistency",
                            "sire missing 2019", "sire missing 2017", "2019 both missing", "2017 both missing"))+
  theme(axis.text.x = element_text(size=16, colour="black", angle=90, vjust = 0.5), 
        axis.text.y = element_text(size=16, colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

#paternity and maternity inconsistencies are due to 'dummy' IDs assigned by Sequoia, i.e we don't know mother or father 
#but dummy IDs have multiple links in pedigree so need to act as placeholders to figure out relatedness among other individuals
#dummy IDs are not consistent between runs, i.e between pedigrees
#some dummy paternities have changed to proper IDs in new pedigree - so it's fine to use new pedigree paternities
#there is 1 true maternity error SCI87: old Mum BLO80 is true, new Mum BLOSS is false (grandmother); fixed below

#fix issue of missing parent in new pedigree that was present in old pedigree - use old parent ID  

#separte out errors of missing parents that should be known
#exclude maternal/paternal inconsistency entries
ped_old_new_comp_missing_error<-subset(ped_old_new_comp_error, !ped_error=="mat_error"&!ped_error=="pat_error")
#exclude errors where old parent is missing - is updated/fine in new pedigree
ped_old_new_comp_missing_error<-subset(ped_old_new_comp_missing_error, !ped_error=="Mum_old_missing" &!ped_error=="Dad_old_missing")

#exlude error entries based on dummy variable inconsitiencies (M0 or F0 IDs)
#exclude dummy varibales in ID column, old father and old mother ID column 
#(not real error as old fummy IDs are labelled randomly & not consistent with new pedigree dummy IDs)
exclude_ID<-grep("^M0|^F0", ped_old_new_comp_missing_error$ID)
ped_old_new_comp_missing_error<-ped_old_new_comp_missing_error[-c(exclude_ID), ]
exclude_Father_old<-grep("^M0", ped_old_new_comp_missing_error$Father_old)
ped_old_new_comp_missing_error<-ped_old_new_comp_missing_error[-c(exclude_Father_old), ]
exclude_Mother_old<-grep("^F0", ped_old_new_comp_missing_error$Mother_old)
ped_old_new_comp_missing_error<-ped_old_new_comp_missing_error[-c(exclude_Mother_old), ]

mother_dummy_IDs<-data.frame(ID=grep("^F0", ped_old_new_comp_missing_error$Mother, value=T))
father_dummy_IDs<-data.frame(ID=grep("^M0", ped_old_new_comp_missing_error$Father, value=T))
ped_old_new_comp_missing_error$ped_error<-as.character(ped_old_new_comp_missing_error$ped_error)
ped_old_new_comp_missing_error_fixed<-ped_old_new_comp_missing_error

ped_old_new_comp_missing_error_fixed$ped_error<-ifelse(ped_old_new_comp_missing_error_fixed$ped_error=="Both_old_missing" &
                                                   (ped_old_new_comp_missing_error_fixed$Mother %in% mother_dummy_IDs$ID)==T & 
                                                   (ped_old_new_comp_missing_error_fixed$Father %in% father_dummy_IDs$ID)==F, "Dad_new_missing",
                                                 ifelse(ped_old_new_comp_missing_error_fixed$ped_error=="Both_old_missing" &
                                                          (ped_old_new_comp_missing_error_fixed$Mother %in% mother_dummy_IDs$ID)==F & 
                                                           (ped_old_new_comp_missing_error_fixed$Father %in% father_dummy_IDs$ID)==T, "Mum_new_missing",
                                                        ifelse((ped_old_new_comp_missing_error_fixed$Father %in% father_dummy_IDs$ID)==T &
                                                                 (ped_old_new_comp_missing_error_fixed$Mother %in% mother_dummy_IDs$ID)==T, "no_error",
                                                               ped_old_new_comp_missing_error$ped_error)))

ped_old_new_comp_missing_error_fixed<-subset(ped_old_new_comp_missing_error_fixed, !ped_error=="no_error")

write.table(ped_old_new_comp_missing_error_fixed, "Pedigree_deer_missing_parent_2017_19.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

ped_old_new_cleaned<-ped_old_new_comp_error

ped_old_new_cleaned$Mother<-ifelse(ped_old_new_cleaned$ped_error=="Mum_new_missing"|ped_old_new_cleaned$ped_error=="mat_error"|
                                     ped_old_new_cleaned$ped_error=="Both_new_missing", 
                                   ped_old_new_cleaned$Mother_old, ped_old_new_cleaned$Mother)

ped_old_new_cleaned$Father<-ifelse(ped_old_new_cleaned$ped_error=="Dad_new_missing"| ped_old_new_cleaned$ped_error=="Both_new_missing", ped_old_new_cleaned$Father_old,
                                      ped_old_new_cleaned$Father)

ped_old_new_cleaned<-ped_old_new_cleaned[, c("ID", "Mother", "Father")]

ped_ok<-subset(ped, !ID %in% ped_old_new_cleaned$ID)
ped_updated<-rbind(ped_ok, ped_old_new_cleaned)


#check if missing Mums are due to missing assignment of 2018 (new pedigree) calves
ped_updated_mum_na<-subset(ped_updated, is.na(ped_updated$Mother))
calves_2018_no_mum<-grep(".*18$", ped_updated_mum_na$ID, value = T) #CML18 and PKN18 are born to study area Mums
                                                                    #need to update their maternity assignment

ped_updated$Mother<-as.character(ped_updated$Mother)
ped_updated$ID<-as.character(ped_updated$ID)
ped_updated$Father<-as.character(ped_updated$Father)

ped_updated$Mother[which(ped_updated$ID=="CML18")]<-"SLA06"
ped_updated$Mother[which(ped_updated$ID=="PKN18")]<-"CNK04"

#try to integrate as many IDs that only appear in old pedigree (might be in genotype data - excluded quite a few from genotype data
#because they weren't in the new pedigree)
# ped_old$links<-ifelse((ped_old$ID %in% ped_old$Mother_old==F) & 
#                              (ped_old$ID %in% ped_old$Father_old)==F, "no_link", "has_links")
ped_old_only<-subset(ped_old, !ID %in% ped_updated$ID)

#join IDs from old pedigree to updated pedigree if they have links (i.e are offsprin and parent)
#ped_old_only<-subset(ped_old_only, links=="has_links")
colnames(ped_old_only)[c(1:3)]<-c("ID", "Mother", "Father")
ped_updated<-rbind(ped_updated, ped_old_only[, c(1:3)])

#sort pedigree (ancestors before offspring)
library(pedigreeR)
ped_updated<- editPed(sire=ped_updated$Father, dam= ped_updated$Mother, label=ped_updated$ID) 
ped_updated<-ped_updated[, c("label", "dam", "sire")]
colnames(ped_updated)<-c("ID", "Mother", "Father")

#exclude IDs that are not in databas - not study area animals or ID is questioned
exclude_IDs<-data.frame(ID=grep(".*_X|.*_Q", ped_updated$ID, value = T))
ped_updated<-subset(ped_updated, !ID %in% exclude_IDs$ID)

#make sure all Mother and Father IDs appear in ID column

ped_updated$ped_status<-ifelse(!ped_updated$Mother%in%ped_updated$ID & is.na(ped_updated$Mother==F), "missing_ID_Mum",
                                ifelse(!ped_updated$Father%in%ped_updated$ID & is.na(ped_updated$Father==F), "missing_ID_Dad",
                                       "complete"))
                                 
Dads_missing<-subset(ped_updated, ped_status=="missing_ID_Dad")
Dads_missing<-as.data.frame(Dads_missing$Father)
Dads_missing<-na.omit(Dads_missing) #all missing Dad IDs NA -don't need to add

Mums_missing<-subset(ped_updated, ped_status=="missing_ID_Mum")
Mums_missing<-as.data.frame(Mums_missing$Mother)
Mums_missing<-na.omit(Mums_missing) #all missing Mum IDs NA -don't need to add

ped_updated<-ped_updated[, c("ID", "Mother", "Father")]
ped_updated$Mother[which(ped_updated$Mother=="unknown")]<-NA
ped_updated$Father[which(ped_updated$Father=="unknown")]<-NA

write.table(ped_updated, file = "Pedigree_deer_2017_19_merged.txt", sep = "\t",
            col.names = T, row.names = F, quote = F)

#LINKPHASE analysis: get haplotypes from SNP genotypes using the pedigree

# #read in new merged pedigree
# ped<-read.table("Pedigree_deer_2017_19_merged.txt", header=T)
# 
# #sort pedigree (ancestors before offspring) and make sure all sire/dam IDs appear in id column
# library(pedigreeR)
# ped<- editPed(sire=ped$Father, dam= ped$Mother, label=ped$ID) 
# ped<-ped[, c("label", "sire", "dam")]
# colnames(ped)<-c("ID", "Father", "Mother")
# 
# #make new ID column with IDs made up of numbers (integers) only
# ped<-ped %>%
#   mutate(ID_new=c(1:nrow(ped)))
# head(ped)
# 
# ped_temp<-ped[, c("ID", "ID_new")]
# head(ped_temp)
# 
# #make new Mum ID and join to ped data frame
# colnames(ped_temp)<- c("Mother", "Mother_ID")
# ped<-join(ped, ped_temp, by="Mother")
# 
# #make new Dad ID and join to ped data frame
# colnames(ped_temp)<- c("Father", "Father_ID")
# ped<-join(ped, ped_temp, by="Father")
# 
# #set NA/missing mother and father ID to 0
# ped$Mother<-as.character(ped$Mother)
# ped$Father<-as.character(ped$Father)
# ped$ID<-as.character(ped$ID)
# ped[is.na(ped)]<-0
# 
# length(unique(ped$ID))
# ID_index_linkphase<-ped[, c("ID", "ID_new")]
# 
# ped_sub<-ped[c("ID_new",  "Father_ID", "Mother_ID")]
# 
# 
# write.table(ped_sub, file="Deer_LINKPHASE_pedigree_file_new.txt", sep = "\t", col.names = F, row.names = F,
#             quote=F)
# write.table(ID_index_linkphase, file="ID_conversion_linkphase_new_ped.txt", sep="\t", col.names = T, row.names = F, quote=F)
# 
# 
# #genotype file
# 
# #geno_tab<-read.delim("Plate1-34_v2_sub2_recoded.ped", header=F)
# geno_tab<-read.delim("20200326_Deer_1to35_SNPQC_CEL_recoded.ped", header=F)
# 
# geno_tab_part<-geno_tab[, c(1:13)]
# head(geno_tab_part)
# # geno_tab_part<-geno_tab_part[, c(2, 7)]
# # colnames(geno_tab_part)[1]<-"ID"
# # table(geno_tab_part$ID)
# # geno_tab_part<-subset(geno_tab_part, ID %in% ped$ID)
# 
# geno_tab<-geno_tab[, c(2, 7:ncol(geno_tab))]
# colnames(geno_tab)[1]<-"ID"
# #table(geno_tab$ID)
# geno_tab<-subset(geno_tab, ID %in% ped$ID)
# genotyped_IDs<-data.frame(ID=geno_tab$ID)
# 
# write.table(genotyped_IDs, file = "IDs_genotyped_reclustered_SNPs.txt", sep="\t", col.names = T, 
#             row.names = F, quote = F)
# 
# geno_tab<-join(geno_tab, ped)
# colnames(geno_tab)[c(39588:39593)]
# geno_tab<-geno_tab[, c(39591, 2:39588)]
# colnames(geno_tab)[c(1:3)]
# 
# #read in map file data with SNP names and add them as column names to geno_tab data frame
# #geno_map<-read.table("Plate1-34_v2_sub2_recoded.map", header=F)
# geno_map<-read.delim("20200326_Deer_1to35_SNPQC_CEL_recoded.map", header=F)
# colnames(geno_tab)[-1]<-as.character(geno_map$V2)
# 
# #marker file
# #read in deer linkage map for SNP positions (and to split the genotype file per chromosome)
# link_map<-read.table("Cervus_elaphus_linkage_map_data.ordered.txt", header=T)
# head(link_map)
# 
# #only use SNPs from linkage map that are in the reclustered marker/map file
# link_map<-subset(link_map, SNP.Name %in% geno_map$V2)
# 
# #subset marker file and then genotype file per chromosome (using SNP names per chromosome from linkage map)
# library(data.table)
# chr=1
# for (chr in c(1:33)) {
#   link_map_sub<-subset(link_map, CEL.LG==chr)
#   #head(link_map_sub)
#   link_map_sub<-link_map_sub[c("CEL.order", "SNP.Name", "Estimated.Mb.Position")]
#   link_map_sub_dup<-subset(link_map_sub, duplicated(Estimated.Mb.Position))
#   link_map_sub<-subset(link_map_sub, !SNP.Name %in% link_map_sub_dup$SNP.Name)
#   
#   write.table(link_map_sub, file=paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/Linkphase_run/Input_files/Deer_LINKPHASE_marker_file_new_chr", chr,".txt"), 
#               sep="\t", col.names = F, row.names = F, quote = F)
#   
#   SNP_names<-as.character(droplevels(link_map_sub$SNP.Name))
#   geno_tab_sub<-geno_tab[c("ID_new",SNP_names)]
#   #order genotype columns according to positions in linkage map
#   geno_tab_dt<-setDT(geno_tab_sub[-1])
#   geno_tab_dt<-setcolorder(geno_tab_dt, as.character(link_map_sub$SNP.Name))
#   geno_tab_dt<-as.data.frame(geno_tab_dt)
#   geno_tab_dt<-cbind(geno_tab_sub[, 1], geno_tab_dt)
#   colnames(geno_tab_dt)[1]<-"ID"
#   
#   write.table(geno_tab_dt, file=paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/Linkphase_run/Input_files/Deer_LINKPHASE_genotype_file_new_chr", chr,".txt"), 
#               sep="\t", col.names = F, row.names = F, quote = F)
#   
# }
# 
# ID<-geno_tab_dt[, 1]
# 
# #make linkphase specification file for each chromosome
# link_spec<-read.delim("linkin.txt", header=F)
# link_spec$V1<-as.character(link_spec$V1)
# 
# for (chr in c(1:33)) {
#   link_spec_temp<-link_spec
#   link_spec_temp[4,1]<-paste0("Deer_LINKPHASE_genotype_file_new_chr", chr, ".txt")
#   link_spec_temp[6,1]<-paste0("Deer_LINKPHASE_marker_file_new_chr", chr, ".txt")
#   
#   write.table(link_spec_temp, file=paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/Linkphase_run/Input_files/linkin_chr", chr,".txt"), 
#               col.names = F, row.names = F, quote = F)
# }
# 
# 




##### AlphaPeel input files preparation (phasing & haplotype identification -instead of linkphase) - input files same as for AlphaImpute/Phase ####
library(plyr)
library(dplyr)
#genotype file
geno_tab<-read.delim("data/20200326_Deer_1to35_SNPQC_CEL_recoded.ped", header=F)
geno_tab_part<-geno_tab[, c(2, 1:2006)]
colnames(geno_tab_part)[1]<-"ID"
head(geno_tab_part)[c(1:10)]

#recode genotypes from 1 1 / 1 2 / 2 2 - where 1 represents the minor allele and 2 the major allele -
#to 0 / 1 / 2, missing genotypes coded as 0 will be set to 9 
#(key: 0 = homozygous minor allele, 1 = heterozygous, 2 = homozygous major allele and 9 = missing genotype)

###data frame with part of the SNP data for testing#####
# geno_tab_part[, ]<-sapply(geno_tab_part[, ], as.character)
# geno_tab_part[ geno_tab_part=="0 0" ]<-9
# geno_tab_part[ geno_tab_part=="1 1" ]<-0
# geno_tab_part[ geno_tab_part=="1 2" ]<-1
# geno_tab_part[ geno_tab_part=="2 2" ]<-2

###full data set###
#select columns with ID and genotypes (two columns per locus one for each allele) only
geno_tab<-geno_tab[, c(2, 7:ncol(geno_tab))]
colnames(geno_tab)[1]<-"ID"
head(geno_tab)[c(1:10)]
#recode
geno_tab[, ]<-sapply(geno_tab[, ], as.character)
geno_tab[ geno_tab=="0 0" ]<-9
geno_tab[ geno_tab=="1 1" ]<-0
geno_tab[ geno_tab=="1 2" ]<-1
geno_tab[ geno_tab=="2 2" ]<-2



#pedigree file
ped<-read.table("Pedigree_deer_2017_19_merged.txt", header=T, stringsAsFactors = F)
head(ped)
ped<-ped[, c(1:3)]
colnames(ped)[1]<-"ID"
ped[is.na(ped)]<-0
#ped<-subset(ped, !Mother==0)
ped<-ped[, c("ID", "Father", "Mother")]
ped<-subset(ped, !duplicated(ID))
ped<-subset(ped, ID %in% geno_tab$ID)

  geno_tab_part<-subset(geno_tab_part, ID %in% ped$ID)

write.table(geno_tab_part, file="Deer_AlphaPhase_genotype_file_part.txt", sep = " ",
            col.names = F, row.names = F, quote = F)

write.table(ped, file="Deer_AlphaPhase_pedigree_file.txt", sep = " ",
            col.names = F, row.names = F, quote = F)

ped<-read.table("results/Deer_AlphaPhase_pedigree_file.txt", header = F)

geno_tab<-subset(geno_tab, ID %in% ped$V1)

#read in map file data with SNP names and add them as column names to geno_tab data frame
geno_map<-read.table("data/20200326_Deer_1to35_SNPQC_CEL_recoded.map", header=F)
colnames(geno_tab)[-1]<-as.character(geno_map$V2)
#subset genotypes to only contain SNPs in linkage map
link_map<-read.table("data/Cervus_elaphus_linkage_map_data.ordered.txt", header = T)
head(link_map)
snps.use<-names(geno_tab)[names(geno_tab) %in% link_map$SNP.Name]

geno_tab_link<-geno_tab[, snps.use] 

#order genotype columns according to positions in linkage map
library(data.table)
geno_tab_dt<-setDT(geno_tab_link)
#subset link map to only contain SNPs in geno_tab (snps.use)
snps.use.df<-data.frame(SNP.Name=snps.use)
link_map_sub<-subset(link_map, SNP.Name %in% snps.use.df$SNP.Name)



#make loop to create genotype files for each chromosome separately
chr=1

for (chr in c(1:33) ) {
  geno_sub<-subset(link_map_sub, CEL.LG==chr)
  #geno_sub_dup<-subset(geno_sub, duplicated(Estimated.Mb.Position))
  #geno_sub<-subset(geno_sub, !SNP.Name %in% geno_sub_dup$SNP.Name)
  marker_map<-geno_sub[, c("CEL.order", "SNP.Name", "Estimated.Mb.Position")]
  
  write.table(marker_map, 
              file=paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaPeel_run/Input_files/Deer_AlphaPeel_marker_file_chr",chr,".txt"), sep = " ",
              col.names = T, row.names = F, quote = F)
  
  # geno_tab_temp<-setcolorder(geno_tab_dt, as.character(geno_sub$SNP.Name))
  # geno_tab_temp<-as.data.frame(geno_tab_temp)
  # geno_tab_temp<-cbind(geno_tab[, 1], geno_tab_temp)
  # colnames(geno_tab_dt)[1]<-"ID"
  # 
  # write.table(geno_tab_temp, 
  #             file=paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaPeel_run/Input_files/Deer_AlphaPeel_genotype_file_chr",chr,".txt"), sep = " ",
  #             col.names = F, row.names = F, quote = F)
}


#input file containing all chromosomes
#exclude sex chromosome snps
link_map_sub<-subset(link_map_sub, !CEL.LG==34)
snps.use_sub<-link_map_sub$SNP.Name
geno_tab_link<-geno_tab[, snps.use_sub] 
geno_tab_dt<-setDT(geno_tab_link)

geno_tab_temp<-setcolorder(geno_tab_dt, as.character(link_map_sub$SNP.Name))
geno_tab_temp<-as.data.frame(geno_tab_temp)
geno_tab_temp<-cbind(geno_tab[, 1], geno_tab_temp)
head(geno_tab_temp)[c(1:10)]

write.table(geno_tab_temp, 
            file=paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaImpute_run/Input_files/Deer_AlphaImpute_genotype_file_all.txt"),
            sep = " ",col.names = F, row.names = F, quote = F)
write.table(link_map_sub, 
            file=paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaImpute_run/Input_files/Deer_AlphaImpute_marker_file_all.txt"), 
            sep = " ", col.names = T, row.names = F, quote = F)

#get list of SNPs per chromosome (for AlphaPhase/AlphaImpute spec file)
SNP_info<-data.frame()

for (chr in c(1:33)){
  geno_files<-list.files(path="C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaImpute_run/Input_files",
                         pattern=paste0("^Deer_AlphaPhase_genotype_file_chr",chr, ".txt"))
  geno_files_list<-as.list(geno_files)
  for (gf in geno_files_list){
    assign("geno_df", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaImpute_run/Input_files/", 
                                        gf), header = F))
    SNP_info_temp<-data.frame(Chromosome=chr, N_SNPs=(ncol(geno_df)-1))
    SNP_info<-rbind(SNP_info, SNP_info_temp)
  }
} 

write.table(SNP_info, file= "results/SNPs_per_chromosome_info.txt", sep="\t", col.names = T, row.names = F)


#make spec files for AlphaImpute/AlphaPhase - specify chromosome specific genotype file and number of SNPs
#AlphaImpute
spec_file<-read.delim("data/AlphaImputeSpec.txt", header = F)
spec_file_sep<-data.frame(do.call('rbind', strsplit(as.character(spec_file$V1), ',', fixed=T)))

#add more tail and core length options
spec_file_sep[7,2]<-paste0(100,",",200,",",300,",",400,",",500)
spec_file_sep[8,2]<-paste0(50,",",150,",",250,",",350,",",450)

SNP_info<-read.table("results/SNPs_per_chromosome_info.txt", header=T)

for (chr in c(1:33)) {
  SNP_info_sub<-subset(SNP_info, Chromosome==chr)
  spec_file_temp<-data.frame(X1=c("PedigreeFile", "GenotypeFile", "SexChrom", "NumberSnp"), 
                             X2=c("DeerAlphaPhasePedigreeFile.txt", paste0("Deer_AlphaPhase_genotype_file_chr", chr, ".txt"),
                                  "No", SNP_info_sub$N_SNPs))
  spec_file_temp<-rbind(spec_file_temp, spec_file_sep[c(5:20), ])
  
  
  write.table(spec_file_temp, file = paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaImpute_run/Input_files/AlphaImputeSpec_chr", chr, ".txt"), 
              sep=",", col.names = F, row.names = F, quote = F)
  
}





#get information about kindepth, phased founders and phasing rate to define phased founder haplotypes


#read in merged pedigree file
ped<-read.table("Pedigree_deer_2017_19_merged.txt", header=T)
head(ped)

kindepth_ped<-kindepth(ped$ID, ped$Father, ped$Mother)
ped$kindepth<-kindepth_ped
max(ped$kindepth)
min(ped$kindepth)
hist(ped$kindepth)

#subset pedigree to kindepth of at least 2
ped_sub<-subset(ped, kindepth>=2)


#get sex and birthyear data
db<-"C:\\Users\\s1767711\\Documents\\Red_deer\\database\\current\\RedDeer1.87.accdb"
con<-odbcConnectAccess2007(db)
queries<-sqlTables(con, tableType = "VIEW")
by_sex_data<-sqlFetch(con, "Birthyear_and_sex_data")
odbcClose(con)

head(by_sex_data)
colnames(by_sex_data)<-c("ID", "sex", "birthyear") #sex in database file (2=male, 1=female)

ped<-join(ped, by_sex_data, by="ID")

ped_meankin<-ped %>%
  group_by(birthyear)%>%
  summarise(mean_kindepth=mean(kindepth))

ped_meankin<-as.data.frame(ped_meankin)
max(ped_meankin$mean_kindepth)
ped_meankin<-na.omit(ped_meankin)

p_kindepth<-ggplot(ped_meankin, aes(as.numeric(birthyear), mean_kindepth))+geom_point(stat="identity")+
  stat_smooth()+
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y=element_text(size = 16, colour = "black"),
        axis.title = element_text(size=20))+xlab("cohort")+ylab("mean kindepth")

ggsave("Kindepth_per_cohort.png", p_kindepth, width = 25, height = 25, units = "cm" )


#filter for true founders - kindepth 0 but have offspring in pedigree
ped$true_founder<-ifelse(ped$kindepth==0 & ped$ID %in% ped$Mother, "yes",
                         ifelse(ped$kindepth==0 & ped$ID %in% ped$Father, "yes", "no"))
table(ped$true_founder)

ped_founders<-subset(ped, true_founder=="yes")

ped_founders<-ped_founders %>%
  group_by(birthyear)%>%
  summarise(N_founders=n())

ped_founders<-na.omit(ped_founders)



p_founder<-ggplot(ped_founders, aes(as.numeric(birthyear), N_founders))+geom_point(stat="identity")+
  stat_smooth(method="loess")+
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y=element_text(size = 16, colour = "black"),
        axis.title = element_text(size=20))+xlab("cohort")+ylab("N founders")

ggsave("Founders_per_cohort.png", p_founder, width = 25, height = 25, units = "cm")

#how many individuals are genotyped in each cohort?
genotyped_IDs<-read.table("IDs_genotyped_reclustered_SNPs.txt", header=T)

geno_cohort<-join(genotyped_IDs, ped, by="ID")
cohorts<-unique(geno_cohort$birthyear)
cohorts<-sort(cohorts)

co<-cohorts[1]
geno_cohort_start<-subset(geno_cohort, birthyear==co)
geno_cohort_cumsum_all<-geno_cohort_start%>%
  summarise(ID_cumsum=n())
geno_cohort_cumsum_all$cohort<-cohorts[1]

for (c in c(2:length(cohorts))) {
  co<-cohorts[c]
  geno_cohort_sub<-subset(geno_cohort, birthyear==co)
  geno_cohort_sum<-geno_cohort_sub%>%
    summarise(N_geno_ID=n())
  co_pre<-data.frame(previous=cohorts[c(1:(c-1))])
  geno_cohort_pre<-subset(geno_cohort, birthyear %in% co_pre$previous)
  geno_pre_sum<-geno_cohort_pre%>%
    summarise(N_geno_ID=n())
  geno_cohort_cumsum<-data.frame(ID_cumsum=(geno_cohort_sum$N_geno_ID+geno_pre_sum$N_geno_ID),
                                 cohort=geno_cohort_sub$birthyear[1])
  geno_cohort_cumsum_all<-rbind(geno_cohort_cumsum_all, geno_cohort_cumsum)
}


p_geno_cumsum<-ggplot(geno_cohort_cumsum_all, aes(as.numeric(cohort), ID_cumsum))+geom_point(stat="identity")+
  stat_smooth()+
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y=element_text(size = 16, colour = "black"),
        axis.title = element_text(size=20))+xlab("cohort")+ylab("cumsum genotyped IDs")

ggsave("Cumsum_genotyped_IDs_per_cohort.png", p_geno_cumsum, width = 25, height=25, units = "cm")


#compare phased genotype concordance of LINKPHASE and AlphaPeel

#ID_index_link<-read.table("ID_conversion_linkphase.txt", header=T)
ID_index_link<-read.table("results/ID_conversion_linkphase_new_ped.txt", header=T)

phase_compare_ID_all<-data.frame()
phase_compare_SNP_all<-data.frame()

for (chr in c(1:33)) {
  
   #LINKPHASE
   assign("phase_chr_link", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/Linkphase_run/Results/phasesLINK_chr",
                                         chr), header=F))
   assign("marker_map_chr_link", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/Linkphase_run/Input_files/Deer_LINKPHASE_marker_file_new_chr",
                                              chr,".txt"), header=F))

  #AlphaPeel phase files
  assign("phase_chr_peel", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaPeel_run/Results/AlphaPeel_out_chr", 
                                        chr, ".called_phase.0.95"), header=F))
  assign("marker_map_chr_peel", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaPeel_run/Input_files/Deer_AlphaPeel_marker_file_chr", 
                                             chr,".txt"), header=T))
  
  
  #LINKPHASE
  colnames(marker_map_chr_link)<-c("marker.order", "SNP.name", "Est.Mb.pos")
  colnames(phase_chr_link)<-c("ID_new", "hap.origin", as.character(marker_map_chr_link$SNP.name))
  SNP.names_link<-c(as.character(marker_map_chr_link$SNP.name))
  phase_chr_link<-join(phase_chr_link, ID_index_link, by="ID_new")
  phase_chr_link<-phase_chr_link[, c("ID", "hap.origin", SNP.names_link)]
  phase_chr_link$hap.origin<-ifelse(phase_chr_link$hap.origin == 1, "paternal", "maternal")
  phase_chr_melt_link<-melt(phase_chr_link, id.vars = c("ID", "hap.origin"))
  names(phase_chr_melt_link)[names(phase_chr_melt_link)=="variable"]<-"SNP.name"
  names(phase_chr_melt_link)[names(phase_chr_melt_link)=="value"]<-"allele_link"
  
  #AlphaPeel
  colnames(marker_map_chr_peel)<-c("marker.order", "SNP.name", "Est.Mb.pos")
  colnames(phase_chr_peel)<-c("ID", as.character(marker_map_chr_peel$SNP.name))
  SNP.names_peel<-c(as.character(marker_map_chr_peel$SNP.name))
  IDs_discard<-data.frame(ID=grep("^FatherOf|^MotherOf", phase_chr_peel$ID, value = T)) #get rid of dummy IDs
  IDs_discard2<-subset(phase_chr_peel, !ID %in% phase_chr_melt_link$ID)
  IDs_discard2<-data.frame(ID=IDs_discard2[, "ID"])
  IDs_discard<-rbind(IDs_discard, IDs_discard2)
  phase_chr_peel<-subset(phase_chr_peel, !ID %in% IDs_discard$ID)
  hap.origin.names<-as.vector(c("paternal", "maternal"))
  phase_chr_peel$hap.origin<-rep(hap.origin.names, (nrow(phase_chr_peel)/2))
  phase_chr_melt_peel<- reshape2::melt(phase_chr_peel, id.vars = c("ID", "hap.origin"))
  names(phase_chr_melt_peel)[names(phase_chr_melt_peel)=="variable"]<-"SNP.name"
  names(phase_chr_melt_peel)[names(phase_chr_melt_peel)=="value"]<-"allele_peel"
  
  #code AlphaPeel alleles the same as in LINKPHASE
  phase_chr_melt_peel[phase_chr_melt_peel=="1"]<-3 #dummy name for major allele before replacing it with 2
  phase_chr_melt_peel[phase_chr_melt_peel=="0"]<-1
  phase_chr_melt_peel[phase_chr_melt_peel=="3"]<-2
  phase_chr_melt_peel[phase_chr_melt_peel=="9"]<-0
  
  phase_chr_melt_peel<-subset(phase_chr_melt_peel, SNP.name %in% phase_chr_melt_link$SNP.name)
  
  phase_chr_comp<-join(phase_chr_melt_link, phase_chr_melt_peel)
  
  phase_chr_comp$phase_error<-ifelse(phase_chr_comp$allele_link==0 & phase_chr_comp$allele_peel!=0, 
                                     "missing_link_only",
                                     ifelse(phase_chr_comp$allele_link!=0 & phase_chr_comp$allele_peel==0,
                                            "missing_peel_only",
                                            ifelse(phase_chr_comp$allele_link!=phase_chr_comp$allele_peel, "phase_dis",
                                                   ifelse(phase_chr_comp$allele_link==0 & phase_chr_comp$allele_peel==0, "missing_both","phase_agree"))))
 
  phase_chr_comp_ID<-phase_chr_comp %>%
    group_by(ID, phase_error) %>%
    summarize(count=n())
  
  phase_chr_comp_ID<-reshape2::dcast(phase_chr_comp_ID, ID ~ phase_error, value.var = "count")
  
  phase_chr_comp_ID<-phase_chr_comp_ID %>%
    mutate(Chromosome=chr)
  
  phase_chr_comp_SNP<-phase_chr_comp %>%
    group_by(SNP.name, phase_error) %>%
    summarize(count=n())
  
  phase_chr_comp_SNP<-reshape2::dcast(phase_chr_comp_SNP, SNP.name ~ phase_error, value.var = "count")
  
  phase_chr_comp_SNP<-phase_chr_comp_SNP %>%
    mutate(Chromosome=chr)
  
  phase_compare_SNP_all<-rbind(phase_compare_SNP_all, phase_chr_comp_SNP)
  
  phase_compare_ID_all<-rbind(phase_compare_ID_all, phase_chr_comp_ID)
  
}

head(phase_compare_ID_all)

#sum phasing information per individual over all chromosomes
phase_compare_ID_sum<-phase_compare_ID_all %>%
  group_by(ID) %>%
  summarise(sum(missing_both), sum(missing_link_only),
            sum(missing_peel_only), sum(phase_agree), sum(phase_dis))

#get total count of unphased SNPs/ID and proportion (unphased)
phase_compare_ID_sum[is.na(phase_compare_ID_sum)]<-0
phase_compare_ID_sum$phase_na_peel<-phase_compare_ID_sum$`sum(missing_both)`+phase_compare_ID_sum$`sum(missing_peel_only)`
phase_compare_ID_sum$phase_na_link<-phase_compare_ID_sum$`sum(missing_both)`+phase_compare_ID_sum$`sum(missing_link_only)`

phase_compare_ID_prop<-phase_compare_ID_sum[, -1]/(nrow(phase_compare_SNP_all)*2)
phase_compare_ID_prop<-cbind(phase_compare_ID_sum[, 1], phase_compare_ID_prop)
colnames(phase_compare_ID_prop)<-c("ID", "prop.missing_both", "prop.missing_link_only", "prop.missing_peel_only",
                                   "prop.phase_agree", "prop.phase_dis", "prop.unphased_peel", "prop.unphased_link")

#add birthyear, sex, kindepth and founder information to SNP/ID phasing rate data

#pedigree for kindepth information
ped<-read.table("results/Pedigree_deer_2017_19_merged.txt", header=T)
kindepth_ped<-kindepth(ped$ID, ped$Father, ped$Mother)
ped$kindepth<-kindepth_ped  


#get sex and birthyear data
db<-"C:\\Users\\s1767711\\Documents\\Red_deer\\database\\current\\RedDeer1.93.accdb"
con<-odbcConnectAccess2007(db)
queries<-sqlTables(con, tableType = "VIEW")
by_sex_data<-sqlFetch(con, "Birthyear_and_sex_data")
odbcClose(con)

head(by_sex_data)
colnames(by_sex_data)<-c("ID", "sex", "birthyear") #sex in database file (2=male, 1=female)
by_sex_data$sex<-ifelse(by_sex_data$sex==1, "female", ifelse(by_sex_data$sex==2, "male", NA))
ped<-join(ped, by_sex_data, by="ID")

#filter for true founders - kindepth 0 but have offspring in pedigree
ped$true_founder<-ifelse(ped$kindepth==0 & ped$ID %in% ped$Mother, "yes",
                         ifelse(ped$kindepth==0 & ped$ID %in% ped$Father, "yes", "no"))

#join phasing and ped tables
phase_compare_ID_prop<-join(phase_compare_ID_prop, ped[, c("ID", "sex", "birthyear", "kindepth", "true_founder")])
write.table(phase_compare_ID_prop, file="Phasing_rate_link.peel_kindepth_info.txt", sep="\t", col.names = T,
            row.names = F, quote = F)

##compare phasing rate against kindepth##

#LINKPHASE
phase_na_mean_link<-phase_compare_ID_prop %>%
  group_by(kindepth) %>%
  summarise(mean_phasing_rate=mean((1-prop.unphased_link)))

p_phasing<-ggplot(phase_na_mean_link, aes(kindepth, mean_phasing_rate))+geom_point(stat="identity")+
  stat_smooth()+
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y=element_text(size = 16, colour = "black"),
        axis.title = element_text(size=20))+xlab("kindepth")+ylab("mean phasing rate")
  
ggsave("Phasing_rate_vs_kindepth.png", p_phasing, width=25, height = 25, units = "cm")

#AlphaPeel
phase_na_mean_peel<-phase_compare_ID_prop %>%
  group_by(kindepth) %>%
  summarise(mean_phasing_rate=mean((1-prop.unphased_peel)))

phase_na_mean_peel<-phase_compare_ID_prop %>%
  group_by(birthyear) %>%
  summarise(mean_phasing_rate=mean((1-prop.unphased_peel)))%>%
  na.omit()


p_phasing_peel<-ggplot(phase_na_mean_peel, aes(birthyear, mean_phasing_rate))+ #kindepth
  geom_point(stat="identity")+stat_smooth()+
  theme_bw()+
  theme(axis.text.x = element_text(size=20, colour = "black"), axis.ticks = element_line(size=1),
        axis.text.y=element_text(size = 20, colour = "black"),
        axis.title = element_text(size=26))+
  xlab("Cohort")+ylab("Mean phasing rate")+
  scale_x_continuous(breaks = c(1960, 1965, 1970, 1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015))


ggsave("Phasing_rate_vs_kindepth_AlphaPeel.png", p_phasing_peel, width=25, height = 25, units = "cm")
ggsave("plots/Phasing_rate_vs_cohort_AlphaPeel.png", p_phasing_peel, width=25, height = 25, units = "cm")

#compare kindepth vs phasing rate for LINKPHASE vs AlphaPeel
phase_na_mean_link<-phase_na_mean_link %>%
  mutate(program="LINKPHASE")

phase_na_mean_peel<-phase_na_mean_peel %>%
  mutate(program="AlphaPeel")

phase_na_mean_link_peel<-rbind(phase_na_mean_link, phase_na_mean_peel)
library(RColorBrewer)
my_colors <- RColorBrewer::brewer.pal(9, "Blues")[1:9]
my_colors_sub<-my_colors[c(5,8)]

p_phasing_linkpeel<-ggplot(phase_na_mean_link_peel, aes(kindepth, mean_phasing_rate, color=program))+
  geom_point(stat="identity")+stat_smooth()+ 
  scale_colour_manual(values = my_colors_sub )+
  theme(axis.text.x = element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size = 16, colour = "black"),
        axis.line = element_line(colour="black"),
        legend.text = element_text(size=16), legend.title = element_text(size = 18),
        axis.title.y= element_text(size=20, margin=margin(r=20)),
        axis.title.x= element_text(size=20, margin=margin(t=20)))+xlab("kindepth")+ylab("mean phasing rate")

ggsave("Phasing_rate_vs_kindepth_link_peel.png", p_phasing_linkpeel, width=25, height = 20, units = "cm")


##compare mean phasing rate of LINKPHASE and AlphaPeel and mean kindepth between the two sexes##
phase_compare_ID_prop<-subset(phase_compare_ID_prop, !is.na(sex) )

phase_na_sex_mean<-phase_compare_ID_prop %>%
  group_by(sex) %>%
  summarise(mean_phasing_rate=mean(1-prop.unphased_link))

phase_kindepth_sex<-phase_compare_ID_prop %>%
  group_by(sex) %>%
  summarise(mean_kindepth=mean(kindepth))

phase_na_sex_mean$variable<-"phasing_rate_link"
names(phase_na_sex_mean)[names(phase_na_sex_mean)=="mean_phasing_rate"]<-"value"
phase_na_sex_mean$value<-phase_na_sex_mean$value*7

phase_kindepth_sex$variable<-"kindepth"
names(phase_kindepth_sex)[names(phase_kindepth_sex)=="mean_kindepth"]<-"value"

phase_kin_rate_mean<-rbind(phase_kindepth_sex, phase_na_sex_mean)

phase_na_sex_mean_peel<-phase_compare_ID_prop %>%
  group_by(sex) %>%
  summarise(mean_phasing_rate=mean(1-prop.unphased_peel))

phase_na_sex_mean_peel$variable<-"phasing_rate_peel"
names(phase_na_sex_mean_peel)[names(phase_na_sex_mean_peel)=="mean_phasing_rate"]<-"value"
phase_na_sex_mean_peel$value<-phase_na_sex_mean_peel$value*7

phase_kin_rate_mean<-rbind(phase_kin_rate_mean, phase_na_sex_mean_peel)

p<-ggplot(phase_kin_rate_mean, aes(x=sex, y=value))+
  geom_bar(aes(fill=variable), stat="identity", position = "dodge")+xlab("sex")+ylab("mean kindepth")

p2<-p+ scale_y_continuous(sec.axis = sec_axis(~./7, name = "mean phasing rate"))+
    scale_fill_brewer(type="qual",palette = "Accent", labels=c("kindepth", "LINKPHASE", "AlphaPeel"))+
    theme(axis.text.x = element_text(size=16, colour = "black"), 
          axis.text.y=element_text(size = 16, colour = "black"),
          axis.title.y = element_text(size=20, margin = margin(t = 20, r = 20, b = 0, l = 0)),
          axis.title.y.right  = element_text(size=20, margin = margin(t = 0, r = 0, b = 0, l = 20)),
          axis.title.x = element_text(size=20, margin = margin(t = 20, r = 0, b = 0, l = 0)),
          legend.title = element_blank(), legend.text = element_text(size=12, colour="black"))
        
ggsave("Phasing_rate_linkpeel_kindepth_by_sex.png", p2, width = 25, height = 25, units = "cm")


#compare ID phasing rate between LINKPHASE and AlphaPeel

p_ID_unphased_linkpeel<-ggplot(phase_compare_ID_prop, aes(prop.unphased_peel, prop.unphased_link))+geom_point()+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))+
  xlab("Prop. unphased SNPs/ID AlphaPeel")+ylab("Prop. unphased SNPs/ID LINKPHASE")

ggsave("ID_phasing_rate_link_vs_peel.png", p_ID_unphased_linkpeel, width = 25, height = 25, units = "cm")

head(phase_compare_SNP_all)
#get total count of unphased IDs/SNP and proportion (unphased)
phase_compare_SNP_all$phase_na_peel<-phase_compare_SNP_all$missing_both+phase_compare_SNP_all$missing_peel_only
phase_compare_SNP_all$phase_na_link<-phase_compare_SNP_all$missing_both+phase_compare_SNP_all$missing_link_only
phase_compare_SNP_all<-phase_compare_SNP_all[, c(1,7, 2:6, 8:9)]
phase_compare_SNP_prop<-phase_compare_SNP_all[, -c(1,2)]/(nrow(phase_compare_ID_sum)*2)
phase_compare_SNP_prop<-cbind(phase_compare_SNP_all[, c(1,2)], phase_compare_SNP_prop)
colnames(phase_compare_SNP_prop)<-c("SNP.name", "Chromosome", "prop.missing_both", "prop.missing_link_only", "prop.missing_peel_only",
                                    "prop.phase_agree", "prop.phase_dis", "prop.unphased_peel", "prop.unphased_link")

#check MAF before and after phasing
IDs_genotyped<-read.table("results/IDs_genotyped_reclustered_SNPs.txt", header=T)

phase_chr_MAF_all<-data.frame()

for (chr in c(1:33)) {
  
  # #LINKPHASE
  # assign("phase_chr", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/Linkphase_run/Results/phasesLINK_chr", 
  #                                       chr), header=F))
  # assign("marker_map_chr", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/Linkphase_run/Input_files/Deer_LINKPHASE_marker_file_new_chr", 
  #                                            chr,".txt"), header=F))
  
  #AlphaPeel phase files
  assign("phase_chr", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaPeel_run/Results/AlphaPeel_out_chr", 
                                        chr, ".called_phase.0.95"), header=F))
  assign("marker_map_chr", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaPeel_run/Input_files/Deer_AlphaPeel_marker_file_chr", 
                                             chr,".txt"), header=T))
  
  colnames(marker_map_chr)<-c("marker.order", "SNP.name", "Est.Mb.pos")
  
  #LINKPHASE
  # colnames(phase_chr)<-c("ID_new", "hap.origin", as.character(marker_map_chr$SNP.name))
  # SNP.names<-c(as.character(marker_map_chr$SNP.name))
  # phase_chr$hap.origin<-ifelse(phase_chr$hap.origin == 1, "paternal", "maternal")
  # phase_chr_melt<-melt(phase_chr, id.vars = c("ID_new", "hap.origin"))
  # phase_chr_melt_na<-subset(phase_chr_melt, value==0)
  # phase_chr_melt<-subset(phase_chr_melt, !value==0)
  
  #AlphaPeel
  colnames(phase_chr)<-c("ID", as.character(marker_map_chr$SNP.name))
  SNP.names<-c(as.character(marker_map_chr$SNP.name))
  phase_chr<-subset(phase_chr, ID %in% IDs_genotyped$ID)
  hap.origin.names<-as.vector(c("paternal", "maternal"))
  phase_chr$hap.origin<-rep(hap.origin.names, (nrow(phase_chr)/2))
  phase_chr_melt<- reshape2::melt(phase_chr, id.vars = c("ID", "hap.origin"))
  names(phase_chr_melt)[names(phase_chr_melt)=="variable"]<-"SNP.name"
  
  #exclude entries with missing allele
  phase_chr_melt<-subset(phase_chr_melt, !value==9)
  
  phase_chr_freq<-phase_chr_melt%>%
    group_by(SNP.name, value)%>%
    summarise(N=n())
  
  phase_chr_freq<-phase_chr_freq%>%
    group_by(SNP.name)%>%
    summarise(total=sum(N))%>%
    left_join(phase_chr_freq, ., by="SNP.name")
  
  #LINKPHASE
  # phase_chr_MAF<-phase_chr_freq%>%
  #   mutate(fraction=N/total)%>%
  #   filter(value==1)
  
  #AlphaPeel
  phase_chr_MAF<-phase_chr_freq%>%
    mutate(fraction=N/total)%>%
    filter(value==0)
  
  names(phase_chr_MAF)[names(phase_chr_MAF)=="fraction"]<-"MAF"
  names(phase_chr_MAF)[names(phase_chr_MAF)=="N"]<-"N_minor"
  phase_chr_MAF<-as.data.frame(phase_chr_MAF)
  
  phase_chr_MAF_all<-rbind(phase_chr_MAF_all, phase_chr_MAF)
  
}



#get MAF information of new reclustered (unphased) data (MAF calculated in plink)
names(phase_chr_MAF_all)[names(phase_chr_MAF_all)=="MAF"]<-"MAF_phased"
geno_MAF_reclustered<-read.table("data/20200326_Deer_1to35_SNPQC_CEL_recoded.frq", header=T)
names(geno_MAF_reclustered)[names(geno_MAF_reclustered)=="SNP"]<-"SNP.name"
names(geno_MAF_reclustered)[names(geno_MAF_reclustered)=="MAF"]<-"MAF_unphased"

MAF_compare<-join(phase_chr_MAF_all, geno_MAF_reclustered[, c("SNP.name", "MAF_unphased")], by="SNP.name")
MAF_compare<-join(MAF_compare, phase_compare_SNP_prop[, c("SNP.name", "Chromosome", "prop.unphased_peel")], by="SNP.name")
MAF_compare_sub<-subset(MAF_compare, prop.unphased_peel<0.10)

p_maf_compare<-ggplot(MAF_compare, aes(MAF_phased, MAF_unphased))+geom_point()+stat_smooth(method = "lm")+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))+
  xlab("MAF after phasing")+ylab("MAF before phasing")

#ggsave("Unphased_reclustered_vs_phased_MAF.png", p_maf_compare, width = 25, height = 25, units = "cm")
ggsave("plots/Unphased_reclustered_vs_phased_MAF_AlphaPeel.png", p_maf_compare, width = 25, height = 25, units = "cm")


p_maf_vs_unphased<-ggplot(MAF_compare, aes(prop.unphased_peel, MAF_phased))+geom_point()+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))+
  xlab("Proportion unphased IDs/SNP")+ylab("MAF after phasing")

#ggsave("UnphasedIDs_prop_vs_unphasedMAF_reclustered.png", p_maf_vs_unphased, width = 25, height = 25, units = "cm")
ggsave("UnphasedIDs_prop_vs_unphasedMAF_AlphaPeel.png", p_maf_vs_unphased, width = 25, height = 25, units = "cm")

+stat_smooth(method = "lm")

# write.table(MAF_compare, file="Compare_MAF_pre_post_phasing_reclustered.txt", sep="\t", col.names = T,
#             row.names = F, quote=F)
write.table(MAF_compare, file="results/Compare_MAF_pre_post_phasing_AlphaPeel.txt", sep="\t", col.names = T,
            row.names = F, quote=F)

#compare SNP phasing rate between LINKPHASE and AlphaPeel

p_SNP_unphased_linkpeel<-ggplot(phase_compare_SNP_prop, aes(prop.unphased_peel, prop.unphased_link))+geom_point()+stat_smooth(method = "lm")+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20))+
  xlab("Prop. unphased IDs/SNP AlphaPeel")+ylab("Prop. unphased IDs/SNP LINKPHASE")


ggsave("SNP_phasing_rate_link_vs_peel.png", p_SNP_unphased_linkpeel, width = 25, height = 25, units = "cm")

write.table(phase_compare_SNP_prop, "Phasing_rate_SNP_linkpeel_compare.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)



#compare MAF (unphased) between sexes

#read in sex data
by_sex_data<-read.table("data/Birthyear_sex_data_2019.txt", header=TRUE)

#list of genotyped IDs
IDs_genotyped<-read.table("results/IDs_genotyped_reclustered_SNPs.txt", header=T)

IDs_sexed<-join(IDs_genotyped, by_sex_data, by="ID")
IDs_sexed$fam_id<-1

#split IDs by sex
IDs_female<-subset(IDs_sexed, sex==1)
IDs_female<-IDs_female[, c("fam_id", "ID")]
IDs_male<-subset(IDs_sexed, sex==2)
IDs_male<-IDs_male[, c("fam_id", "ID")]

#write out
write.table(IDs_female, file = "results/IDs_genotyped_female.txt", sep="\t", col.names = F, row.names = F, quote = F)
write.table(IDs_male, file = "results/IDs_genotyped_male.txt", sep="\t", col.names = F, row.names = F, quote = F)


#compare sex specific MAF (unphased) - MAF calculation done in plink
geno_MAF_fem<-read.table("data/20200326_Deer_1to35_SNPQC_CEL_recoded_fem.frq", header=T)
names(geno_MAF_fem)[names(geno_MAF_fem)=="MAF"]<-"MAF_fem"

geno_MAF_male<-read.table("data/20200326_Deer_1to35_SNPQC_CEL_recoded_male.frq", header=T)
names(geno_MAF_male)[names(geno_MAF_male)=="MAF"]<-"MAF_male"

geno_MAF_sex<-join(geno_MAF_fem, geno_MAF_male, by=c("SNP", "CHR"))
geno_MAF_sex$NCHROBS<-NULL
geno_MAF_sex$A1<-NULL
geno_MAF_sex$A2<-NULL

#only phased SNPs
marker_map_all_chr<-read.table("results/Deer_AlphaPeel_marker_file_all_chr.txt", header=T)
geno_MAF_sex_sub<-subset(geno_MAF_sex, SNP %in% marker_map_all_chr$SNP.Name)

p_MAF_sex<-ggplot(geno_MAF_sex_sub, aes(MAF_fem, MAF_male))+geom_point()+stat_smooth(method = "lm")+
  theme(axis.text = element_text(size=16, colour = "black"), axis.title = element_text(size=20))+
  xlab("female specific MAF")+ylab("male specific MAF")

cor.test(geno_MAF_sex$MAF_fem, geno_MAF_sex$MAF_male)

ggsave("plots/Sex_specific_MAF_phased_SNPs.png", p_MAF_sex, width = 25, height = 25, units = "cm" )

