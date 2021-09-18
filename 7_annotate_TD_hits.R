#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Red deer: annotate transmission distortion hits - using phased (AlphaPeel) and data based on genedrop simulations #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(tidyr)
library(plyr)
library(dplyr)
library(reshape2)
library(biomaRt)
library(xtable)

#read in linkage map with new cattle genome positions and MAF info table
linkmap_cow_pos<-read.table("data/Cervus_elaphus_linkage_map_additional_BT_pos.txt", header=T)
head(linkmap_cow_pos)
names(linkmap_cow_pos)[names(linkmap_cow_pos)=="SNP.Name"]<-"SNP.name"

#what is the distance between SNPs (based on cow positions)?

linkmap_cow_pos_dist<-data.frame()

for (chr in c(1:29)) {
  
  linkmap_cow_pos_chr<-subset(linkmap_cow_pos, cow_ARS_UCD1.2_Chr==chr)
  linkmap_cow_pos_chr<-subset(linkmap_cow_pos_chr, !is.na(cow_ARS_UCD1.2_pos))
  linkmap_cow_pos_chr<-arrange(linkmap_cow_pos_chr, cow_ARS_UCD1.2_pos)
  
  marker.dist<-diff(linkmap_cow_pos_chr$cow_ARS_UCD1.2_pos)
  marker.dist<-append(marker.dist, NA)
  linkmap_cow_pos_chr$marker.dist_cow1.2<-marker.dist
  
  linkmap_cow_pos_dist<-rbind(linkmap_cow_pos_dist, linkmap_cow_pos_chr)
  
}

#take out 0 distances (sika and red deer markers with same positions)
linkmap_cow_pos_dist<-subset(linkmap_cow_pos_dist, !marker.dist_cow1.2==0)

#plot marker distances
hist(linkmap_cow_pos_dist$marker.dist_cow1.2, breaks=100)
mean.dist<-mean(linkmap_cow_pos_dist$marker.dist_cow1.2, na.rm = T) #70kbp
min(linkmap_cow_pos_dist$marker.dist_cow1.2, na.rm = T)

MAF_info<-read.table("results/Compare_MAF_pre_post_phasing_AlphaPeel.txt", header=T)
head(MAF_info)

#SNPs showing TD by offspring sex
offspringtrans_SNPs_sig<-read.table("results/offspringtrans_TD_top_SNPs_AlphaPeel.txt", header=T)
head(offspringtrans_SNPs_sig)

#join MAF and cow position info to sig SNPs
offspringtrans_SNPs_sig<-join(offspringtrans_SNPs_sig, linkmap_cow_pos[, c("SNP.name", "cow_ARS_UCD1.2_pos", "cow_ARS_UCD1.2_Chr")])
offspringtrans_SNPs_sig<-join(offspringtrans_SNPs_sig, MAF_info[, c("SNP.name", "MAF_unphased")])

#filter out SNP with low MAF and no matching position in new cattle genome assembly
offspringtrans_SNPs_sig<-subset(offspringtrans_SNPs_sig, MAF_unphased>0.05)
offspringtrans_SNPs_sig<-subset(offspringtrans_SNPs_sig, !is.na(cow_ARS_UCD1.2_pos))

ensembl_btaurus    = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", 
                                dataset=c("btaurus_gene_ensembl"))

listMarts()
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")
datasets<-listDatasets(ensembl) #btaurus is version ARS-UCD1.2


#define region to be searched for genes covering the area of all significant windows
offspringtrans_SNPs_sig<-arrange(offspringtrans_SNPs_sig, cow_ARS_UCD1.2_Chr)
offspringtrans_SNPs_sig<- offspringtrans_SNPs_sig %>%
  mutate(region_full=paste0(cow_ARS_UCD1.2_Chr, ":",(cow_ARS_UCD1.2_pos-70000), ":", (cow_ARS_UCD1.2_pos+70000) ))


#extract gene names and description, transcript start and end points etc from ensembl database (per SNP region)
#using the cow genome


offspringtrans_SNPs_info_full<-data.frame()
#n<-1

for (n in c(1:nrow(offspringtrans_SNPs_sig))) {
  
  offspringtrans_SNPs_sig_sub<-offspringtrans_SNPs_sig[n, ]
  
  offspringtrans_SNps_info<-getBM(attributes=c('external_gene_name','description', 
                                             'phenotype_description', 'go_id', 'ensembl_gene_id' ,'name_1006', 'definition_1006', 
                                             'transcript_start', 'transcript_end', 'start_position', 'end_position'),
                                filters = "chromosomal_region",
                                values = offspringtrans_SNPs_sig_sub[["region_full"]], mart = 
                                  ensembl_btaurus)
  offspringtrans_SNps_info<-offspringtrans_SNps_info %>%
    mutate(SNP=offspringtrans_SNPs_sig_sub[["SNP.name"]],sex=offspringtrans_SNPs_sig_sub[["sex"]])
  
  
  offspringtrans_SNPs_info_full<-rbind(offspringtrans_SNPs_info_full, offspringtrans_SNps_info)
  
}

offspringtrans_SNPs_info_full<-subset(offspringtrans_SNPs_info_full, !description=="")

offspringtrans_SNPs_info_full<-offspringtrans_SNPs_info_full %>%
  mutate(species="btaurus")

##put in gene name in description column as external gene name if latter is missing
offspringtrans_SNPs_info_full$external_gene_name<-ifelse(offspringtrans_SNPs_info_full$external_gene_name=="", 
                                                         offspringtrans_SNPs_info_full$description,
                                                         offspringtrans_SNPs_info_full$external_gene_name)


#find homologs from mice and human genomes

offspringtrans_SNPs_homologs_full<-data.frame()

for (n in c(1:nrow(offspringtrans_SNPs_sig))) {
  
  offspringtrans_SNPs_sig_sub<-offspringtrans_SNPs_sig[n, ]
  
  offspringtrans_SNPs_homologs<-getBM(attributes=c('mmusculus_homolog_ensembl_gene', 'mmusculus_homolog_associated_gene_name', 'mmusculus_homolog_orthology_type',
                                                   'mmusculus_homolog_orthology_confidence','hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name', 'hsapiens_homolog_orthology_type',
                                                   'hsapiens_homolog_orthology_confidence'),
                                  filters = "chromosomal_region",
                                  values = offspringtrans_SNPs_sig_sub[["region_full"]], mart = ensembl_btaurus)
  
  offspringtrans_SNPs_homologs<-offspringtrans_SNPs_homologs %>%
    mutate(SNP=offspringtrans_SNPs_sig_sub[["SNP.name"]],sex=offspringtrans_SNPs_sig_sub[["sex"]])
  
  
  offspringtrans_SNPs_homologs_full<-rbind(offspringtrans_SNPs_homologs_full, offspringtrans_SNPs_homologs)
  
}


offspringtrans_SNPs_homologs_full<-subset(offspringtrans_SNPs_homologs_full, !mmusculus_homolog_ensembl_gene =="" & !hsapiens_homolog_ensembl_gene=="")
offspringtrans_SNPs_homologs_full<-arrange(offspringtrans_SNPs_homologs_full, SNP, sex)

#add homolog info to gene table
#make sure gene table only contains distinct SNP, sex gene pairings - get rid of all GO term information
offspringtrans_SNPs_info_genes_distinct<-offspringtrans_SNPs_info_full %>%
  distinct(SNP, sex, ensembl_gene_id, external_gene_name, description, .keep_all=T)

offspringtrans_SNPs_info_genes_distinct<-offspringtrans_SNPs_info_genes_distinct[, c("SNP", "sex", "ensembl_gene_id",
                                                                                     "external_gene_name", "description")]

offspringtrans_SNPs_info_genes_distinct<-arrange(offspringtrans_SNPs_info_genes_distinct, SNP, sex)

offspringtrans_SNPs_homologs_full_sub<-subset(offspringtrans_SNPs_homologs_full, SNP %in% offspringtrans_SNPs_info_genes_distinct$SNP)

offspringtrans_SNPs_info_genes_homologs<-join(offspringtrans_SNPs_info_genes_distinct, offspringtrans_SNPs_homologs_full_sub)


#get rid of duplicated rows by matching mouse and human orthologs to cow ensembl genes

for (nr in c(1:nrow(offspringtrans_SNPs_info_genes_homologs))) {
  
  offspringtrans_SNPs_info_genes_homologs$mouse_match[nr]<-grepl(offspringtrans_SNPs_info_genes_homologs$external_gene_name[nr], 
  x = offspringtrans_SNPs_info_genes_homologs$mmusculus_homolog_associated_gene_name[nr], ignore.case=TRUE)
  
}


for (nr in c(1:nrow(offspringtrans_SNPs_info_genes_homologs))) {
  
  offspringtrans_SNPs_info_genes_homologs$human_match[nr]<-grepl(offspringtrans_SNPs_info_genes_homologs$external_gene_name[nr], 
                                                                 x = offspringtrans_SNPs_info_genes_homologs$hsapiens_homolog_associated_gene_name[nr], ignore.case=TRUE)
  
}

#genes with clear orthologs in humans or mice
offspringtrans_SNPs_homologs<-subset(offspringtrans_SNPs_info_genes_homologs, !(mouse_match==FALSE & human_match==FALSE))

#make sure there is no more duplicated rows 
offspringtrans_SNPs_homologs<- offspringtrans_SNPs_homologs %>%
  distinct(SNP, sex, ensembl_gene_id, mmusculus_homolog_associated_gene_name,
           hsapiens_homolog_associated_gene_name, .keep_all=T)

#genes with no clear orthologs in humans or mice 
offspringtrans_SNPs_no_homologs<-subset(offspringtrans_SNPs_info_genes_homologs, 
                                        !external_gene_name %in% offspringtrans_SNPs_homologs$external_gene_name)

#make sure there is no more duplicated rows
offspringtrans_SNPs_no_homologs<- offspringtrans_SNPs_no_homologs %>%
  distinct(SNP, sex, ensembl_gene_id, .keep_all=T)

offspringtrans_SNPs_no_homologs[, c("mmusculus_homolog_ensembl_gene","mmusculus_homolog_associated_gene_name",
                                    "mmusculus_homolog_orthology_type", "mmusculus_homolog_orthology_confidence",
                                    "hsapiens_homolog_ensembl_gene","hsapiens_homolog_associated_gene_name","hsapiens_homolog_orthology_type",       
                                    "hsapiens_homolog_orthology_confidence")]<-NA

offspringtrans_SNPs_ortholog_info<-rbind(offspringtrans_SNPs_homologs, offspringtrans_SNPs_no_homologs)
offspringtrans_SNPs_ortholog_info<-arrange(offspringtrans_SNPs_ortholog_info, sex, SNP)

#add transcript positions of genes and SNP position information
offspringtrans_SNPs_ortholog_info<-join(offspringtrans_SNPs_ortholog_info, 
                                        offspringtrans_SNPs_info_full[, c("ensembl_gene_id", "start_position", "end_position")])

names(offspringtrans_SNPs_ortholog_info)[names(offspringtrans_SNPs_ortholog_info) =="SNP"]<-"SNP.name"

offspringtrans_SNPs_ortholog_info<-join(offspringtrans_SNPs_ortholog_info, offspringtrans_SNPs_sig[, c("SNP.name", "cow_ARS_UCD1.2_pos", "cow_ARS_UCD1.2_Chr")])

#make sure there is no more duplicated rows
offspringtrans_SNPs_ortholog_info<- offspringtrans_SNPs_ortholog_info %>%
  distinct(SNP.name, sex, ensembl_gene_id, mmusculus_homolog_associated_gene_name,
           hsapiens_homolog_associated_gene_name, .keep_all=T)


offspringtrans_SNPs_ortholog_info$gene_hit<-ifelse(offspringtrans_SNPs_ortholog_info$start_position <= offspringtrans_SNPs_ortholog_info$cow_ARS_UCD1.2_pos &
                                                     offspringtrans_SNPs_ortholog_info$end_position >= offspringtrans_SNPs_ortholog_info$cow_ARS_UCD1.2_pos, "yes", "no" )


offspringtrans_SNPs_info_full<-offspringtrans_SNPs_info_full[, c("SNP", "sex", "external_gene_name", "description", 
                                                                 "ensembl_gene_id", "go_id", "name_1006", "species" )]


write.table(offspringtrans_SNPs_info_full, file="results/offspringtrans_sig_SNPs_geneTbl_AlphaPeel.txt", sep="\t", row.names = F,
            col.names = T, quote = F)             

write.table(offspringtrans_SNPs_ortholog_info, file="results/offspringtrans_sig_SNPs_orthologs_transcriptInfo_AlphaPeel.txt", sep="\t", 
            row.names = F, col.names = T, quote = F)

#how many unique genes are present in the region?
offspringtrans_SNPs_info_full_sub<-offspringtrans_SNPs_info_full %>%
  distinct(ensembl_gene_id, external_gene_name, description, sex,  .keep_all = T) #76 unique genes (males and females together; female:25, males:51) )

#shared genes between males and females
offspringtrans_SNPs_info_full_sub2<-subset(offspringtrans_SNPs_info_full_sub, duplicated(ensembl_gene_id))
offspringtrans_SNPs_info_shared<-subset(offspringtrans_SNPs_info_full, 
                                        external_gene_name %in% offspringtrans_SNPs_info_full_sub2$external_gene_name)


#parent TD


#SNPs showing TD by parent sex
Mumtrans_SNPs_sig<-read.table("results/Mumtrans_TD_top_SNPs_AlphaPeel.txt", header=T)
Mumtrans_SNPs_sig$sex<-"female"
Dadtrans_SNPs_sig<-read.table("results/Dadtrans_TD_top_SNPs_AlphaPeel.txt", header=T)
Dadtrans_SNPs_sig$sex<-"male"

MumDadtrans_SNPs_sig<-rbind(Mumtrans_SNPs_sig, Dadtrans_SNPs_sig )

#join MAF and cow position info to sig SNPs
MumDadtrans_SNPs_sig<-join(MumDadtrans_SNPs_sig, linkmap_cow_pos[, c("SNP.name", "cow_ARS_UCD1.2_pos", "cow_ARS_UCD1.2_Chr")])
MumDadtrans_SNPs_sig<-join(MumDadtrans_SNPs_sig, MAF_info[, c("SNP.name", "MAF_unphased")])

#filter out SNP with low MAF and no matching position in new cattle genome assembly
MumDadtrans_SNPs_sig<-subset(MumDadtrans_SNPs_sig, MAF_unphased>0.05)
MumDadtrans_SNPs_sig<-subset(MumDadtrans_SNPs_sig, !is.na(cow_ARS_UCD1.2_pos))


ensembl_btaurus    = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", 
                                dataset=c("btaurus_gene_ensembl"))

listMarts()
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")
datasets<-listDatasets(ensembl) #btaurus is version ARS-UCD1.2


#define region to be searched for genes covering the area of all significant windows
MumDadtrans_SNPs_sig<-arrange(MumDadtrans_SNPs_sig, cow_ARS_UCD1.2_Chr)
MumDadtrans_SNPs_sig<- MumDadtrans_SNPs_sig %>%
  mutate(region_full=paste0(cow_ARS_UCD1.2_Chr, ":",(cow_ARS_UCD1.2_pos-70000), ":", (cow_ARS_UCD1.2_pos+70000) ))

#SNPs that are significant in parents and offspring
parent_offspring_SNPs<-subset(MumDadtrans_SNPs_sig, SNP.name %in% offspringtrans_SNPs_sig$SNP.name)
parent_offspring_SNPs$status<-"parent"
offspringtrans_SNPs_sig$status<-"offspring"
offspring_parent_SNPs<-subset(offspringtrans_SNPs_sig, SNP.name %in% MumDadtrans_SNPs_sig$SNP.name)

parent_offspring_SNPs_sig<-rbind(parent_offspring_SNPs, offspring_parent_SNPs)
parent_offspring_SNPs_sig<-parent_offspring_SNPs_sig[, c(1:9, 14, 10:13)]
parent_offspring_SNPs_sig<-arrange(parent_offspring_SNPs_sig, SNP.name)


write.table(parent_offspring_SNPs_sig, file="results/parent_offspring_sig_SNPs_shared_AlphaPeel.txt", sep="\t", row.names = F,
            col.names = T, quote = F)     

parent_offspring_SNPs_sig<-parent_offspring_SNPs_sig[, c(1:5, 7, 9:10, 13)]
parent_offspring_SNPs_sig<-parent_offspring_SNPs_sig[, c(1,7,8,2,3,5,6,4,9)]

#latex table
print(xtable(parent_offspring_SNPs_sig, type = "latex", digits = c(1,1,1,1,1,1,1,4,1,2),
             display = c("d", "s", "s", "s", "d", "d", "s", "f", "d", "f") ),
      file="results/parent_offspring_sig_SNPs_shared_AlphaPeel_tbl.tex")


#extract gene names and description, transcript start and end points etc from ensembl database (per SNP region)
#using the cow genome


MumDadtrans_SNPs_info_full<-data.frame()
#n<-1

for (n in c(1:nrow(MumDadtrans_SNPs_sig))) {
  
  MumDadtrans_SNPs_sig_sub<-MumDadtrans_SNPs_sig[n, ]
  
  MumDadtrans_SNPs_info<-getBM(attributes=c('external_gene_name','description', 
                                            'phenotype_description', 'go_id', 'ensembl_gene_id' ,'name_1006', 'definition_1006', 
                                            'transcript_start', 'transcript_end', 'start_position', 'end_position'),
                                  filters = "chromosomal_region",
                                  values = MumDadtrans_SNPs_sig_sub[["region_full"]], mart = 
                                    ensembl_btaurus)
  MumDadtrans_SNPs_info<-MumDadtrans_SNPs_info %>%
    mutate(SNP=MumDadtrans_SNPs_sig_sub[["SNP.name"]],sex=MumDadtrans_SNPs_sig_sub[["sex"]])
  
  
  MumDadtrans_SNPs_info_full<-rbind(MumDadtrans_SNPs_info_full, MumDadtrans_SNPs_info)
  
}


MumDadtrans_SNPs_info_full<-subset(MumDadtrans_SNPs_info_full, !description=="")

MumDadtrans_SNPs_info_full<-MumDadtrans_SNPs_info_full %>%
  mutate(species="btaurus")

##put in gene name in description column as external gene name if latter is missing
MumDadtrans_SNPs_info_full$external_gene_name<-ifelse(MumDadtrans_SNPs_info_full$external_gene_name=="", 
                                                      MumDadtrans_SNPs_info_full$description,
                                                      MumDadtrans_SNPs_info_full$external_gene_name)


#find homologs from mice and human genomes

MumDadtrans_SNPs_homologs_full<-data.frame()

for (n in c(1:nrow(MumDadtrans_SNPs_sig))) {
  
  MumDadtrans_SNPs_sig_sub<-MumDadtrans_SNPs_sig[n, ]
  
  MumDadtrans_SNPs_homologs<-getBM(attributes=c('mmusculus_homolog_ensembl_gene', 'mmusculus_homolog_associated_gene_name', 'mmusculus_homolog_orthology_type',
                                                'mmusculus_homolog_orthology_confidence','hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name', 'hsapiens_homolog_orthology_type',
                                                'hsapiens_homolog_orthology_confidence'),
                                   filters = "chromosomal_region",
                                   values = MumDadtrans_SNPs_sig_sub[["region_full"]], mart = ensembl_btaurus)
  
  MumDadtrans_SNPs_homologs<-MumDadtrans_SNPs_homologs %>%
    mutate(SNP=MumDadtrans_SNPs_sig_sub[["SNP.name"]],sex=MumDadtrans_SNPs_sig_sub[["sex"]])
  
  
  MumDadtrans_SNPs_homologs_full<-rbind(MumDadtrans_SNPs_homologs_full, MumDadtrans_SNPs_homologs)
  
}


MumDadtrans_SNPs_homologs_full<-subset(MumDadtrans_SNPs_homologs_full, !mmusculus_homolog_ensembl_gene =="" & !hsapiens_homolog_ensembl_gene=="")
MumDadtrans_SNPs_homologs_full<-arrange(MumDadtrans_SNPs_homologs_full, SNP, sex)

#add homolog info to gene table
#make sure gene table only contains distinct SNP, sex gene pairings - get rid of all GO term information
MumDadtrans_SNPs_info_genes_distinct<-MumDadtrans_SNPs_info_full %>%
  distinct(SNP, sex, ensembl_gene_id, external_gene_name, description, .keep_all=T)

MumDadtrans_SNPs_info_genes_distinct<-MumDadtrans_SNPs_info_genes_distinct[, c("SNP", "sex", "ensembl_gene_id",
                                                                               "external_gene_name", "description")]

MumDadtrans_SNPs_info_genes_distinct<-arrange(MumDadtrans_SNPs_info_genes_distinct, SNP, sex)

MumDadtrans_SNPs_homologs_full_sub<-subset(MumDadtrans_SNPs_homologs_full, SNP %in% MumDadtrans_SNPs_info_genes_distinct$SNP)

MumDadtrans_SNPs_info_genes_homologs<-join(MumDadtrans_SNPs_info_genes_distinct, MumDadtrans_SNPs_homologs_full_sub)


#get rid of duplicated rows by matching mouse and human orthologs to cow ensembl genes

for (nr in c(1:nrow(MumDadtrans_SNPs_info_genes_homologs))) {
  
  MumDadtrans_SNPs_info_genes_homologs$mouse_match[nr]<-grepl(MumDadtrans_SNPs_info_genes_homologs$external_gene_name[nr], 
                                                              x = MumDadtrans_SNPs_info_genes_homologs$mmusculus_homolog_associated_gene_name[nr], ignore.case=TRUE)
  
}


for (nr in c(1:nrow(MumDadtrans_SNPs_info_genes_homologs))) {
  
  MumDadtrans_SNPs_info_genes_homologs$human_match[nr]<-grepl(MumDadtrans_SNPs_info_genes_homologs$external_gene_name[nr], 
                                                              x = MumDadtrans_SNPs_info_genes_homologs$hsapiens_homolog_associated_gene_name[nr], ignore.case=TRUE)
  
}

#genes with clear orthologs in humans or mice
MumDadtrans_SNPs_homologs<-subset(MumDadtrans_SNPs_info_genes_homologs, !(mouse_match==FALSE & human_match==FALSE))

#make sure there is no more duplicated rows 
MumDadtrans_SNPs_homologs<- MumDadtrans_SNPs_homologs %>%
  distinct(SNP, sex, ensembl_gene_id, mmusculus_homolog_associated_gene_name,
           hsapiens_homolog_associated_gene_name, .keep_all=T)

#genes with no clear orthologs in humans or mice 
MumDadtrans_SNPs_no_homologs<-subset(MumDadtrans_SNPs_info_genes_homologs, 
                                     !external_gene_name %in% MumDadtrans_SNPs_homologs$external_gene_name)

#make sure there is no more duplicated rows
MumDadtrans_SNPs_no_homologs<- MumDadtrans_SNPs_no_homologs %>%
  distinct(SNP, sex, ensembl_gene_id, .keep_all=T)

MumDadtrans_SNPs_no_homologs[, c("mmusculus_homolog_ensembl_gene","mmusculus_homolog_associated_gene_name",
                                 "mmusculus_homolog_orthology_type", "mmusculus_homolog_orthology_confidence",
                                 "hsapiens_homolog_ensembl_gene","hsapiens_homolog_associated_gene_name","hsapiens_homolog_orthology_type",       
                                 "hsapiens_homolog_orthology_confidence")]<-NA

MumDadtrans_SNPs_ortholog_info<-rbind(MumDadtrans_SNPs_homologs, MumDadtrans_SNPs_no_homologs)
MumDadtrans_SNPs_ortholog_info<-arrange(MumDadtrans_SNPs_ortholog_info, sex, SNP)

#add transcript positions of genes and SNP position information
MumDadtrans_SNPs_ortholog_info<-join(MumDadtrans_SNPs_ortholog_info, 
                                     MumDadtrans_SNPs_info_full[, c("ensembl_gene_id", "start_position", "end_position")])

names(MumDadtrans_SNPs_ortholog_info)[names(MumDadtrans_SNPs_ortholog_info) =="SNP"]<-"SNP.name"

MumDadtrans_SNPs_ortholog_info<-join(MumDadtrans_SNPs_ortholog_info, MumDadtrans_SNPs_sig[, c("SNP.name", "cow_ARS_UCD1.2_pos", "cow_ARS_UCD1.2_Chr")])

#make sure there is no more duplicated rows
MumDadtrans_SNPs_ortholog_info<- MumDadtrans_SNPs_ortholog_info %>%
  distinct(SNP.name, sex, ensembl_gene_id, mmusculus_homolog_associated_gene_name,
           hsapiens_homolog_associated_gene_name, .keep_all=T)


MumDadtrans_SNPs_ortholog_info$gene_hit<-ifelse(MumDadtrans_SNPs_ortholog_info$start_position <= MumDadtrans_SNPs_ortholog_info$cow_ARS_UCD1.2_pos &
                                                  MumDadtrans_SNPs_ortholog_info$end_position >= MumDadtrans_SNPs_ortholog_info$cow_ARS_UCD1.2_pos, "yes", "no" )


MumDadtrans_SNPs_info_full<-MumDadtrans_SNPs_info_full[, c("SNP", "sex", "external_gene_name", "description", 
                                                           "ensembl_gene_id", "go_id", "name_1006", "species" )]


write.table(MumDadtrans_SNPs_info_full, file="results/MumDadtrans_sig_SNPs_geneTbl_AlphaPeel.txt", sep="\t", row.names = F,
            col.names = T, quote = F)             

write.table(MumDadtrans_SNPs_ortholog_info, file="results/MumDadtrans_sig_SNPs_orthologs_transcriptInfo_AlphaPeel.txt", sep="\t", 
            row.names = F, col.names = T, quote = F)

#how many unique genes are present in the region?
MumDadtrans_SNPs_info_full_sub<-MumDadtrans_SNPs_info_full %>%
  distinct(ensembl_gene_id, external_gene_name, description, sex,  .keep_all = T) #66; male: 42, female: 24

#shared genes between males and females
MumDadtrans_SNPs_info_full_sub2<-subset(MumDadtrans_SNPs_info_full_sub, duplicated(ensembl_gene_id))
MumDadtrans_SNPs_info_shared<-subset(MumDadtrans_SNPs_info_full, 
                                     ensembl_gene_id %in% MumDadtrans_SNPs_info_full_sub2$ensembl_gene_id)

#combine gene annotation information for offspring and parents
MumDadtrans_SNPs_ortholog_info$status<-"parent"
offspringtrans_SNPs_ortholog_info$status<-"offspring"
gene_info_orthologs_all<-rbind(MumDadtrans_SNPs_ortholog_info, offspringtrans_SNPs_ortholog_info)

gene.hit_info_orthologs<-subset(gene_info_orthologs_all, gene_hit=="yes")

#how many genes have SNp its within coding region
gene.hit_info_orthologs_sub<-gene.hit_info_orthologs %>%
  distinct(ensembl_gene_id) #42


#gene hits by sex and status
gene.hit_info_orthologssub2<-gene.hit_info_orthologs %>%
  distinct(ensembl_gene_id, sex, status) %>%
  group_by(sex, status)%>%
  summarise(count=n())  #offspring: female 11, male 13; parents: mothers 7, fathers 14
  

#make table for chapter
colnames(gene_info_orthologs_all)
gene_info_orthologs_all<-gene_info_orthologs_all[, c("SNP.name" , "sex", "status", "ensembl_gene_id", "external_gene_name",
                                                     "description", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name",
                                                     "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name", 
                                                     "start_position", "end_position", "cow_ARS_UCD1.2_pos", "cow_ARS_UCD1.2_Chr")]



write.table(gene_info_orthologs_all, file="results/All_sig_SNPs_orthologs_transcriptInfo_AlphaPeel.txt", sep="\t",
            row.names = F, col.names = T, quote=F)

#add GO terms to table with all transcript info
gene_info_orthologs_all<-read.table("results/All_sig_SNPs_orthologs_transcriptInfo_AlphaPeel.txt", header=T, sep="\t")

MumDadtrans_SNPs_info_full<-read.table("results/MumDadtrans_sig_SNPs_geneTbl_AlphaPeel.txt", sep="\t", header=T)
colnames(MumDadtrans_SNPs_info_full)[1]<-"SNP.name"

offspringtrans_SNPs_info_full<-read.table("results/offspringtrans_sig_SNPs_geneTbl_AlphaPeel.txt", sep="\t", header=T)
colnames(offspringtrans_SNPs_info_full)[1]<-"SNP.name"

gene_info_orthologs_all.1<-subset(gene_info_orthologs_all, SNP.name %in% MumDadtrans_SNPs_info_full$SNP.name)
gene_info_orthologs_all.1<-join(MumDadtrans_SNPs_info_full[, c("SNP.name", "go_id", "ensembl_gene_id", "name_1006")], gene_info_orthologs_all.1, 
                                by=c("SNP.name", "ensembl_gene_id"))

gene_info_orthologs_all.2<-subset(gene_info_orthologs_all, SNP.name %in% offspringtrans_SNPs_info_full$SNP.name)
gene_info_orthologs_all.2<-join(offspringtrans_SNPs_info_full[, c("SNP.name", "go_id", "ensembl_gene_id", "name_1006")], gene_info_orthologs_all.2, 
                                by=c("SNP.name", "ensembl_gene_id"))

gene_info_orthologs_all<-rbind(gene_info_orthologs_all.1, gene_info_orthologs_all.2)
gene_info_orthologs_all<-gene_info_orthologs_all %>%
                         distinct(SNP.name, sex, status, ensembl_gene_id, go_id, .keep_all =T)


#make subset of gene info table SNPs significant in both a parent and offspring and SNPs under balancing selection
parent_offspring_SNPs_sig<-read.table("results/parent_offspring_sig_SNPs_shared_AlphaPeel.txt", header=T)
all_SNPs_empirical_p_sub_sig<-read.table("results/TD_SNPs_sig_selection.txt", header=T)
colnames(all_SNPs_empirical_p_sub_sig)[1]<-"SNP.name"

gene_info_orthologs_sub<-subset(gene_info_orthologs_all, SNP.name %in% parent_offspring_SNPs_sig$SNP.name)
gene_info_orthologs_sub<-gene_info_orthologs_sub %>%
  distinct(SNP.name, ensembl_gene_id, go_id, .keep_all =T)

gene_info_orthologs_sub2<-subset(gene_info_orthologs_all, SNP.name %in% all_SNPs_empirical_p_sub_sig$SNP.name)

gene_info_orthologs_sub<-rbind(gene_info_orthologs_sub, gene_info_orthologs_sub2)

gene_info_orthologs_sub_merged<-gene_info_orthologs_sub %>%
  group_by(SNP.name, external_gene_name) %>%
  summarise(go_ids = paste0(go_id, collapse = ","), 
            go_terms = paste0(name_1006, collapse = ",")) 

#latex table
print(xtable(gene_info_orthologs_sub_merged, type = "latex",
             display = c("d", "s", "s", "s", "s") ),
      file="results/transgenerational_plus_selected_TD_SNPs_genetbl.tex")


#Gowinda gene enrichment analysis

#make SNP input files
head(linkmap_cow_pos)
marker_map_all_chr<-read.table("results/Deer_AlphaPeel_marker_file_all_chr.txt", header=T)
linkmap_cow_pos_sub<-subset(linkmap_cow_pos, SNP.name %in% marker_map_all_chr$SNP.Name)

#all SNPs
gowinda_all_SNPs<-linkmap_cow_pos_sub[, c("cow_ARS_UCD1.2_Chr", "cow_ARS_UCD1.2_pos")]
gowinda_all_SNPs<-subset(gowinda_all_SNPs, !is.na(cow_ARS_UCD1.2_pos))

#SNPs sig in female offspring
gowinda_SNPs_female<-subset(offspringtrans_SNPs_sig, sex=="female")
gowinda_SNPs_female<-gowinda_SNPs_female[, c("cow_ARS_UCD1.2_Chr", "cow_ARS_UCD1.2_pos")]

#SNPs sig in male offspring
gowinda_SNPs_male<-subset(offspringtrans_SNPs_sig, sex=="male")
gowinda_SNPs_male<-gowinda_SNPs_male[, c("cow_ARS_UCD1.2_Chr", "cow_ARS_UCD1.2_pos")]

#SNPs sig in mothers
gowinda_SNPs_mother<-subset(MumDadtrans_SNPs_sig, sex=="female")
gowinda_SNPs_mother<-gowinda_SNPs_mother[, c("cow_ARS_UCD1.2_Chr", "cow_ARS_UCD1.2_pos")]

#SNPs sig in fathers
gowinda_SNPs_father<-subset(MumDadtrans_SNPs_sig, sex=="male")
gowinda_SNPs_father<-gowinda_SNPs_father[, c("cow_ARS_UCD1.2_Chr", "cow_ARS_UCD1.2_pos")]


write.table(gowinda_all_SNPs, "results/Gowinda_total_SNP_file.txt", sep="\t", col.names = F, row.names = F, quote = F)
write.table(gowinda_SNPs_male, "results/Gowinda_candidate_SNP_file_males.txt", sep="\t", col.names = F, row.names = F, quote = F)
write.table(gowinda_SNPs_female, "results/Gowinda_candidate_SNP_file_females.txt", sep="\t", col.names = F, row.names = F, quote = F)
write.table(gowinda_SNPs_mother, "results/Gowinda_candidate_SNP_file_mothers.txt", sep="\t", col.names = F, row.names = F, quote = F)
write.table(gowinda_SNPs_father, "results/Gowinda_candidate_SNP_file_fathers.txt", sep="\t", col.names = F, row.names = F, quote = F)

#get functional annotation file for enrichment analysis

functional_annotation<-getBM(attributes=c('go_id', 'name_1006','ensembl_gene_id'), mart = ensembl_btaurus)
functional_annotation<-subset(functional_annotation, !go_id=="")

functional_annotation_all<-data.frame()

go_ids_all<-unique(functional_annotation$go_id)

for (go in go_ids_all) {
  
  functional_annotation_sub<-subset(functional_annotation, go_id==go)
  
  functional_genes<-data.frame(gene_id=functional_annotation_sub[, c("ensembl_gene_id")])  
  functional_genes<-subset(functional_genes, !duplicated(gene_id))
  functional_genes<-data.frame(gene_id=t(functional_genes))
  colnames(functional_genes)[c(1, ncol(functional_genes))]<-c("first", "last")
  
  if (length(functional_genes)>1){
    functional_genes_merge<- functional_genes %>%
      unite(gene_id_merged, first:last, sep=" ")%>%
      select(gene_id_merged)
    
  } else {
    functional_genes_merge<-functional_genes
    colnames(functional_genes_merge)<-"gene_id_merged"
  }
  
  
  functionaL_annotation_temp<-data.frame(go_id=functional_annotation_sub$go_id[1], description=functional_annotation_sub$name_1006[1],
                                         gene_id=functional_genes_merge$gene_id_merged)
  
  functional_annotation_all<-rbind(functional_annotation_all, functionaL_annotation_temp)
  
  
}


write.table(functional_annotation_all, file="results/Gowinda_go_annotation_file_ensemblIDs_b.taurus.txt", sep="\t", col.names = F, 
            row.names = F, quote = F)




