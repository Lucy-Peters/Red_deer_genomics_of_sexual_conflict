#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Red deer: allele frequency change over time after genedrop simulation  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(plyr)
library(dplyr)
library(ggplot2)
library(GeneDrop)
library(purrr)
library(furrr)

# #empirical MAF per cohort
# IDs_genotyped<-read.table("results/IDs_genotyped_reclustered_SNPs.txt", header=T)
# ped<-read.table("results/Pedigree_deer_2017_19_merged_genedrop.txt", header=T)
# 
# #make sure allele change is comparable to simulation by excluding alleles inherited by dummy parents in simulation from 
# #empirical data set
# load("results/Deer_genedrop_objects_sim_part1.RData")
# genedrop_obj<-genedrop_list[[1]]
# 
# ID_cohorts<-ped[, c("ID", "Cohort")] %>%data.frame()
# 
# phase_chr_MAF_all<-data.frame()
# 
# for (chr in c(1:33)) {
#   
#   #AlphaPeel phase files
#   assign("phase_chr", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaPeel_run/Results/AlphaPeel_out_chr", 
#                                         chr, ".called_phase.0.95"), header=F))
#   assign("marker_map_chr", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaPeel_run/Input_files/Deer_AlphaPeel_marker_file_chr", 
#                                              chr,".txt"), header=T))
#   
#   #calculate allele frequencies per cohort
#   genedrop_ped<-genedrop_obj@pedigree %>%data.frame()
#   colnames(genedrop_ped)<-c("ID", "Father", "Mother", "sex", "Cohort")
#   
#   dummy_founders<-data.frame(ID=grep("^Dummy.*", genedrop_ped$ID, value=T))
#   
#   genedrop_ped$dummy_off<-ifelse(genedrop_ped$Father %in% dummy_founders$ID, "yes_sire", 
#                                  ifelse(genedrop_ped$Mother %in% dummy_founders$ID, "yes_dam", "no"))
#   
#   dummy_offspring<-subset(genedrop_ped, !dummy_off=="no")
#   
#   colnames(marker_map_chr)<-c("marker.order", "SNP.name", "Est.Mb.pos")
#   
#   marker_map_chr<-marker_map_chr %>%
#     mutate(CEL.LG=chr)
#   
#   colnames(phase_chr)<-c("ID", as.character(marker_map_chr$SNP.name))
#   SNP.names<-c(as.character(marker_map_chr$SNP.name))
#   phase_chr<-join(phase_chr, ID_cohorts, by="ID")
#   
#   #make sure IDs are the same as in simulation
#   phase_chr<-subset(phase_chr, ID %in% genedrop_ped$ID)
#   phase_chr<-subset(phase_chr, ID %in% IDs_genotyped$ID)
#   
#   cohorts<-c(1981:2018)
#   
#   phase_chr_MAF_all_co<-data.frame()
#   
#   for (co in cohorts) {
#     
#     phase_chr_sub<-subset(phase_chr, Cohort==co)
#     
#     hap.origin.names<-as.vector(c("paternal", "maternal"))
#     phase_chr_sub$hap.origin<-rep(hap.origin.names, (nrow(phase_chr_sub)/2))
#     phase_chr_melt<- reshape2::melt(phase_chr_sub, id.vars = c("ID", "hap.origin"))
#     names(phase_chr_melt)[names(phase_chr_melt)=="variable"]<-"SNP.name"
#     
#     phase_chr_melt<-join(phase_chr_melt, genedrop_ped, by="ID") 
#     
#     phase_chr_melt$value<-ifelse(phase_chr_melt$dummy_off =="yes_sire" & phase_chr_melt$hap.origin=="paternal", 9,
#                                  ifelse(phase_chr_melt$dummy_off =="yes_dam" & phase_chr_melt$hap.origin=="maternal", 9,
#                                         phase_chr_melt$value))                          
#     
#     
#     
#     #exclude entries with missing allele
#     phase_chr_melt<-subset(phase_chr_melt, !value==9)
#     
#     phase_chr_freq<-phase_chr_melt%>%
#       group_by(SNP.name, value)%>%
#       summarise(N=n())
#     
#     phase_chr_freq<-phase_chr_freq%>%
#       group_by(SNP.name)%>%
#       summarise(total=sum(N))%>%
#       left_join(phase_chr_freq, ., by="SNP.name")
#     
#     #AlphaPeel
#     phase_chr_MAF<-phase_chr_freq%>%
#       mutate(fraction=N/total)%>%
#       filter(value==0)
#     
#     names(phase_chr_MAF)[names(phase_chr_MAF)=="fraction"]<-"MAF"
#     names(phase_chr_MAF)[names(phase_chr_MAF)=="N"]<-"N_minor"
#     phase_chr_MAF<-as.data.frame(phase_chr_MAF)
#     
#     SNPs_MAF_miss<-subset(marker_map_chr, !SNP.name %in% phase_chr_MAF$SNP.name)
#     
#     if (nrow(SNPs_MAF_miss)!=0){
#       phase_chr_frq_sub<-subset(phase_chr_freq, SNP.name %in% SNPs_MAF_miss$SNP.name)
#       phase_chr_MAF_fixed<-data.frame(SNP.name=unique(phase_chr_frq_sub$SNP.name),
#                                       value=0, N_minor=0, total=phase_chr_freq$total[1],
#                                       MAF=1)
#       phase_chr_MAF<-rbind(phase_chr_MAF, phase_chr_MAF_fixed)
#     }
#     
#     #keep SNP name, cohort, no of individuals and MAF
#     phase_chr_MAF$N_IDs<-phase_chr_MAF$total/2
#     phase_chr_MAF$Cohort<-co
#     phase_chr_MAF<-join(phase_chr_MAF, marker_map_chr[c("SNP.name", "CEL.LG")])
#     phase_chr_MAF<-phase_chr_MAF[c("SNP.name", "CEL.LG","MAF", "N_IDs", "Cohort")]
#     
#     phase_chr_MAF_all_co<-rbind(phase_chr_MAF_all_co, phase_chr_MAF)
#     
#   }
#   
#   phase_chr_MAF_all<-rbind(phase_chr_MAF_all, phase_chr_MAF_all_co)
# }
# 
# 
# 
# write.table(phase_chr_MAF_all, file="results/MAF_per_cohort_AlphaPeel.txt", sep="\t", row.names = F, 
#             col.names = T, quote = F)
# 


#empirical p_value calculation run on ashworth server

#define output directory (scratch)
args <- commandArgs(trailingOnly = TRUE)
resultpath <- as.character(args[1])


#empirical MAF per cohort
phase_chr_MAF_all<-read.table("results/MAF_per_cohort_AlphaPeel.txt", header=T, stringsAsFactors = F)


#determine if simulated MAF data has equal or more extreme values for indicators of directional or balancing selection and determine 
#empirical p_value

#list of SNPs significant for TD
Mumtrans_SNPs_sig<-read.table("results/Mumtrans_TD_top_SNPs_AlphaPeel.txt", header=T, stringsAsFactors = F)
Dadtrans_SNPs_sig<-read.table("results/Dadtrans_TD_top_SNPs_AlphaPeel.txt", header=T, stringsAsFactors = F)
offspringtrans_SNPs_sig<-read.table("results/offspringtrans_TD_top_SNPs_AlphaPeel.txt", header=T, stringsAsFactors = F)

TD_snps<-unique(Mumtrans_SNPs_sig$SNP.name)
TD_snps<-append(TD_snps, unique(Dadtrans_SNPs_sig$SNP.name))
TD_snps<-append(TD_snps, unique(offspringtrans_SNPs_sig$SNP.name))


source("scripts/empirical_p_from_simulation_function.R")

#this is how you would get empirical_p_all_SNPs data frame with function (but I did it without function first)

#run function in safe mode
allele_frq_change_test_safely<-purrr::safely(allele_frq_change_test)

# #run function in parallel
# availableCores()
# #specify number of cores to parallelise over
# plan(multiprocess, workers = 2)
# 
# options(future.globals.maxSize = +Inf)
# 
# empirical_p_all_SNPs<-future_map(TD_snps, allele_frq_change_test_safely, output_option = "empirical_p", 
#                       path_out = resultpath, empirical_MAF_df = phase_chr_MAF_all)
                     
                      
extract_results<-function(l, object, type = "tbl"){
  
  if(type == "tbl" & class(object[[l]]$result)!="NULL"){
    
    result_extracted<-object[[l]]$result@empirical_p_all_SNPs_tbl
    
  } else if ( type == "plot" & class(object[[l]]$result)!="NULL"){
    
    result_extracted<-object[[l]]$result@p_MAF_SNP_plot
    
  } else if (type == "plot_only" & class(object[[l]]$result)!="NULL"){
    
    result_extracted<-object[["result"]]@p_MAF_SNP_plot
    
  } else {
    
    result_extracted<-data.frame()
    
  }
  
  result_extracted
}


# results_length<-c(1:length(empirical_p_all_SNPs))
# 
# empirical_p_all_SNPs_tbl<-map(results_length, extract_results, object = empirical_p_all_SNPs, type = "tbl") %>%
#                          bind_rows()

# write.table(empirical_p_all_SNPs_tbl, file=paste0(resultpath, "TD_sig_SNPs_allele_frq_change_p_values_noMAF1.txt"), sep="\t",
#             col.names = T, row.names = F, quote=F)

#make plots of significant TD SNPs showing significant allele frequency change signal
#read in empirical_p_value table
empirical_p_all_SNPs<-read.table("results/TD_sig_SNPs_allele_frq_change_p_values_noMAF1.txt", header=T, stringsAsFactors = F)
empirical_p_all_SNPs_sub<-empirical_p_all_SNPs[c(1:2), ]

TD_selected_snp_plots<- allele_frq_change_test_safely(output_option = "plotting", input_object_range = c(1:40), path_out = resultpath, 
                                                      plot_name = "TD", empirical_p_in = empirical_p_all_SNPs, empirical_MAF_df = phase_chr_MAF_all)

save(TD_selected_snp_plots, file = paste0(resultpath, "TD_sig_SNPs_allele_frq_change_plots.RData"))


#apply FDR correction to empirical p-values for significant TD SNPs
empirical_p_all_SNPs_tbl<-read.table("results/TD_sig_SNPs_allele_frq_change_p_values_noMAF1.txt", header=T)
empirical_p_all_SNPs_tbl<-subset(empirical_p_all_SNPs_tbl, !duplicated(SNP.Name))
empirical_p_all_SNPs_fdr<-apply(empirical_p_all_SNPs_tbl[, c(8:10)], 2, p.adjust, method="fdr")
colnames(empirical_p_all_SNPs_fdr)<-c("p_direct.fdr", "p_balance.fdr", "p_slope.fdr")
empirical_p_all_SNPs_tbl<-cbind(empirical_p_all_SNPs_tbl, empirical_p_all_SNPs_fdr)

#add sex and status to table
empirical_p_all_SNPs_tbl$status<-ifelse(empirical_p_all_SNPs_tbl$SNP.Name %in% Mumtrans_SNPs_sig$SNP.name | 
                                          empirical_p_all_SNPs_tbl$SNP.Name %in% Dadtrans_SNPs_sig$SNP.name, "parent", "offspring")

empirical_p_all_SNPs_tbl$sex<-ifelse(empirical_p_all_SNPs_tbl$SNP.Name %in% Mumtrans_SNPs_sig$SNP.name, "female",
                                     ifelse(empirical_p_all_SNPs_tbl$SNP.Name %in% Dadtrans_SNPs_sig$SNP.name, "male",
                                            ifelse(empirical_p_all_SNPs_tbl$SNP.Name %in% offspringtrans_SNPs_sig$SNP.name,
                                                   offspringtrans_SNPs_sig$sex, NA)))


empirical_p_all_SNPs_tbl<-empirical_p_all_SNPs_tbl[, c("SNP.Name", "status", "sex", "total_MAF_change", "count_direct", "p_direct", "p_direct.fdr",
                                                       "total_MAF_dist", "count_balance", "p_balance", "p_balance.fdr",
                                                       "MAF_slope",  "count_slope", "p_slope", "p_slope.fdr")]

write.table(empirical_p_all_SNPs_tbl, file="results/TD_sig_SNPs_allele_frq_change_p_values_noMAF1.txt", sep="\t",
            col.names = T, row.names = F, quote=F)

#make table with SNPs under significant selection
empirical_p_all_SNPs_sig<-subset(empirical_p_all_SNPs_tbl, p_direct.fdr<0.05 | p_balance.fdr<0.05 | p_slope.fdr<0.05)

#add sex and status to table
empirical_p_all_SNPs_sig$status<-ifelse(empirical_p_all_SNPs_sig$SNP.Name %in% Mumtrans_SNPs_sig$SNP.name | 
                                          empirical_p_all_SNPs_sig$SNP.Name %in% Dadtrans_SNPs_sig$SNP.name, "parent", "offspring")

empirical_p_all_SNPs_sig$sex<-ifelse(empirical_p_all_SNPs_sig$SNP.Name %in% Mumtrans_SNPs_sig$SNP.name, "female",
                                     ifelse(empirical_p_all_SNPs_sig$SNP.Name %in% Dadtrans_SNPs_sig$SNP.name, "male",
                                            ifelse(empirical_p_all_SNPs_sig$SNP.Name %in% offspringtrans_SNPs_sig$SNP.name,
                                                   offspringtrans_SNPs_sig$sex, NA)))

empirical_p_all_SNPs_sig<-empirical_p_all_SNPs_sig[, c("SNP.Name", "status", "sex", "total_MAF_dist", "p_balance.fdr", 
                                                       "total_MAF_change", "p_direct.fdr", "MAF_slope", "p_slope.fdr")]

library(xtable)
print(xtable(empirical_p_all_SNPs_sig, type = "latex", 
             digits= c(1,1,1,1,4,4,4,4,4,4) ),
      file = "results/TD_significant_selection_SNPs_tbl.tex")

#plot signifcant SNPs in one plot
snps.sig<-empirical_p_all_SNPs_sig$SNP.Name
snp<-"cela1_red_15_32873496"

TD_SNPs_sigMAF_df<-data.frame()

for (snp in snps.sig) {
  command<-paste0("TD_selected_snp_plots[[", "'", "result", "'", "]]@p_MAF_SNP_plot$", snp)
  plot<-eval(parse(text = command))
  df<-plot$data
  
  TD_SNPs_sigMAF_df<-rbind(TD_SNPs_sigMAF_df, df)
  
}

p_MAF_SNP<-ggplot(TD_SNPs_sigMAF_df, aes(Cohort, MAF, group=iteration))+
  geom_line(aes(color=source, size=source))+
  scale_color_manual(values = c("red", "black"))+
  scale_size_manual(values= c(rep(0.8, (max(TD_SNPs_sigMAF_df$iteration)-1)), 2)) +
  scale_x_discrete(breaks = TD_SNPs_sigMAF_df$Cohort[seq(1,length(TD_SNPs_sigMAF_df$Cohort), by=2)]) +
  ylim(0, (max(TD_SNPs_sigMAF_df$MAF)+0.05))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=16, colour = "black", angle = -45, vjust=-0.4, margin = margin(t=15)),
        axis.text.y = element_text(size=16, colour = "black"),
        axis.title.x = element_text(size=35, margin = margin(t=20)),
        axis.title.y = element_text(size=35, margin = margin(r=20)),
        axis.line = element_line(colour = "black"), axis.ticks=element_line(size=1), strip.text = element_text(size=20) )+
  facet_wrap(vars(SNP.Name), nrow = 7, ncol = 2)


ggsave("plots/TD_sigSNPs_MAF_change_sim.png", p_MAF_SNP, width = 30, height = 40, units = "cm")

#run empirical p-value test for allele frequency change for all SNPs
#get list of SNPs
linkmap<-read.table("data/TableS1_CervusElaphus_Final_Linkage_Map.txt", header=T, stringsAsFactors=F)
SNP_list<-c(linkmap$SNP.Name)

#run function in parallel
plan(multiprocess, workers = 10)

options(future.globals.maxSize = +Inf)

all_SNPs_empirical_p<-future_map(SNP_list, allele_frq_change_test_safely, output_option = "both", input_object_range = c(1:40),
                                 path_out = resultpath, plot_name = "All_SNPs", empirical_MAF_df = phase_chr_MAF_all)

load("results/All_SNPs_allele_frq_change_empirical_p.RData")
all_SNPs_empirical_p_sub<-all_SNPs_empirical_p[c(1:100)]

results_length<-as.list(c(1:length(all_SNPs_empirical_p)))

all_SNPs_empirical_p_tbl<-map(results_length, extract_results, object = all_SNPs_empirical_p, type = "tbl")%>%
  bind_rows()

write.table(all_SNPs_empirical_p_tbl, file=paste0(resultpath, "All_SNPs_allele_frq_change_p_values_noMAF1.txt"), sep="\t",
            col.names = T, row.names = F, quote=F)

#all_SNPs_empirical_p_plots<-map(results_length, extract_results, object = all_SNPs_empirical_p, type = "plot") 
# 

#             
# save(all_SNPs_empirical_p_plots, file = paste0(resultpath, "All_SNPs_allele_frq_change_plots.RData"))

all_SNPs_empirical_p_tbl<-read.table("results/All_SNPs_allele_frq_change_p_values_noMAF1.txt", header = T)

#apply FDR correction to all SNPs table
all_SNPs_empirical_p_tbl<-subset(all_SNPs_empirical_p_tbl, !duplicated(SNP.Name))
all_SNPs_empirical_p_tbl_fdr<-apply(all_SNPs_empirical_p_tbl[, c(8:10)], 2, p.adjust, method="fdr")
colnames(all_SNPs_empirical_p_tbl_fdr)<-c("p_direct.fdr", "p_balance.fdr", "p_slope.fdr")
all_SNPs_empirical_p_tbl<-cbind(all_SNPs_empirical_p_tbl, all_SNPs_empirical_p_tbl_fdr)

#get total number of iteration used for empiricl p-value calculation (some filtered out bc MAF went to 1)
all_SNPs_empirical_p_tbl$iterations<-(all_SNPs_empirical_p_tbl$count_direct+1)/all_SNPs_empirical_p_tbl$p_direct
all_SNPs_empirical_p_tbl$iterations<-all_SNPs_empirical_p_tbl$iterations-1

#subset the all SNP empirical p value table to significant TD SNPs 

#list of SNPs significant for TD
Mumtrans_SNPs_sig<-read.table("results/Mumtrans_TD_top_SNPs_AlphaPeel.txt", header=T, stringsAsFactors = F)
Dadtrans_SNPs_sig<-read.table("results/Dadtrans_TD_top_SNPs_AlphaPeel.txt", header=T, stringsAsFactors = F)
offspringtrans_SNPs_sig<-read.table("results/offspringtrans_TD_top_SNPs_AlphaPeel.txt", header=T, stringsAsFactors = F)

TD_snps<-unique(Mumtrans_SNPs_sig$SNP.name)
TD_snps<-append(TD_snps, unique(Dadtrans_SNPs_sig$SNP.name))
TD_snps<-append(TD_snps, unique(offspringtrans_SNPs_sig$SNP.name))
TD_snps<-data.frame(SNP.Name=TD_snps)

all_SNPs_empirical_p_tbl_sub<-subset(all_SNPs_empirical_p_tbl, SNP.Name %in% TD_snps$SNP.Name)


#add sex and status to table
all_SNPs_empirical_p_tbl_sub$status<-ifelse(all_SNPs_empirical_p_tbl_sub$SNP.Name %in% Mumtrans_SNPs_sig$SNP.name | 
                                              all_SNPs_empirical_p_tbl_sub$SNP.Name %in% Dadtrans_SNPs_sig$SNP.name, "parent", "offspring")

all_SNPs_empirical_p_tbl_sub$sex<-ifelse(all_SNPs_empirical_p_tbl_sub$SNP.Name %in% Mumtrans_SNPs_sig$SNP.name, "female",
                                     ifelse(all_SNPs_empirical_p_tbl_sub$SNP.Name %in% Dadtrans_SNPs_sig$SNP.name, "male",
                                            ifelse(all_SNPs_empirical_p_tbl_sub$SNP.Name %in% offspringtrans_SNPs_sig$SNP.name,
                                                   offspringtrans_SNPs_sig$sex, NA)))


all_SNPs_empirical_p_tbl_sub<-all_SNPs_empirical_p_tbl_sub[, c("SNP.Name", "status", "sex", "total_MAF_change", "count_direct", "p_direct", "p_direct.fdr",
                                                       "total_MAF_dist", "count_balance", "p_balance", "p_balance.fdr",
                                                       "MAF_slope",  "count_slope", "p_slope", "p_slope.fdr")]

write.table(all_SNPs_empirical_p_tbl_sub, file="results/TD_sig_SNPs_allele_frq_change_p_values_noMAF1.txt", sep="\t",
            col.names = T, row.names = F, quote=F)

#make table with SNPs under significant selection
all_SNPs_empirical_p_sub_sig<-subset(all_SNPs_empirical_p_tbl_sub, p_direct.fdr<0.05 | p_balance.fdr<0.05 | p_slope.fdr<0.05)

all_SNPs_empirical_p_sub_sig<-all_SNPs_empirical_p_sub_sig[, c("SNP.Name", "status", "sex", "total_MAF_dist", "p_balance.fdr", 
                                                               "total_MAF_change", "p_direct.fdr", "MAF_slope", "p_slope.fdr")]

write.table(all_SNPs_empirical_p_sub_sig, file="results/TD_SNPs_sig_selection.txt", sep="\t",
            col.names = T, row.names = F, quote=F)

library(xtable)
print(xtable(all_SNPs_empirical_p_sub_sig, type = "latex", 
             digits= c(1,1,1,1,4,4,4,4,4,4) ),
      file = "results/TD_significant_selection_SNPs_tbl.tex")

#plot signifcant SNPs in one plot
load("results/TD_sig_SNPs_allele_frq_change_plots.RData")

snps.sig<-all_SNPs_empirical_p_sub_sig$SNP.Name
snp<-"cela1_red_15_32873496"

TD_SNPs_sigMAF_df<-data.frame()

for (snp in snps.sig) {
  command<-paste0("TD_selected_snp_plots[[", "'", "result", "'", "]]@p_MAF_SNP_plot$", snp)
  plot<-eval(parse(text = command))
  df<-plot$data
  
  TD_SNPs_sigMAF_df<-rbind(TD_SNPs_sigMAF_df, df)
  
}

alpha<-ifelse(TD_SNPs_sigMAF_df$source == "simulation", 0.001, 0.05)

p_MAF_SNP<-ggplot(TD_SNPs_sigMAF_df, aes(Cohort, MAF, group=iteration))+
  geom_line(aes(color=source, size=source, alpha=alpha))+
  scale_color_manual(values = c("red", "black"))+
  scale_size_manual(values= c(rep(0.8, (max(TD_SNPs_sigMAF_df$iteration)-1)), 2)) +
  scale_x_discrete(breaks = TD_SNPs_sigMAF_df$Cohort[seq(1,length(TD_SNPs_sigMAF_df$Cohort), by=2)]) +
  ylim(0, (max(TD_SNPs_sigMAF_df$MAF)+0.05))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=16, colour = "black", angle = -45, vjust=-0.4, margin = margin(t=15)),
        axis.text.y = element_text(size=16, colour = "black"),
        axis.title.x = element_text(size=35, margin = margin(t=20)),
        axis.title.y = element_text(size=35, margin = margin(r=20)),
        axis.line = element_line(colour = "black"), axis.ticks=element_line(size=1), strip.text = element_text(size=20) )+
  facet_wrap(vars(SNP.Name), nrow = 7, ncol = 2)


ggsave("plots/TD_sigSNPs_MAF_change_sim.png", p_MAF_SNP, width = 30, height = 40, units = "cm")



