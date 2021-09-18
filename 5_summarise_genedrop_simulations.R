#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Red deer: summarising all genedrop simulations #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)

#define output directory (scratch)
args <- commandArgs(trailingOnly = TRUE)
resultpath <- as.character(args[1])
print(resultpath)

#combine transmission tables from all simulation runs by type
table_types<-c("Mumtrans","Mum_off_sextrans","Dadtrans", "Dad_off_sextrans","offspringtranstrans")

for (type in table_types) {
  df_all<-data.frame()
  for (run in c(1:19, 21:22, 24:30, 32:33, 35:37, 39:40)) {
    assign("df_part", read.table(paste0("results/", type,"_genedrop_sim_TD_tablespart", run, ".txt"), header = TRUE)) 
    df_all<-rbind(df_all, df_part)
  }
  assign(paste0(type, "_results_all_sim"), df_all)
}


# #p-value thresholds - before and after cut off of lowest 5% sample size SNPs
# #mother TD
# Mumtrans_p_value<-as.numeric(quantile(Mumtrans_results_all_sim$p_value, probs = 0.001))
# 
# Mumtrans_p_value_frq.cut<- Mumtrans_results_all_sim%>%
#   filter(expected_freq>quantile(Mumtrans_results_all_sim$expected_freq, probs = 0.05))
# Mumtrans_p_value_frq.cut<-as.numeric(quantile(Mumtrans_p_value_frq.cut$p_value, probs = 0.001))
# 
# #mother by offspring sex TD  
# Mum_off_sextrans_p_value<-as.numeric(quantile(Mum_off_sextrans_results_all_sim$p_value, probs = 0.001))
# 
# Mumtrans_off_sex_p_value_frq.cut<- Mum_off_sextrans_results_all_sim%>%
#   filter(expected_freq>quantile(Mum_off_sextrans_results_all_sim$expected_freq, probs = 0.05))
# Mumtrans_off_sex_p_value_frq.cut<-as.numeric(quantile(Mumtrans_off_sex_p_value_frq.cut$p_value, probs = 0.001))
# 
# #father TD
# Dadtrans_p_value<-as.numeric(quantile(Dadtrans_results_all_sim$p_value, probs = 0.001))
# 
# Dadtrans_p_value_frq.cut<- Dadtrans_results_all_sim%>%
#   filter(expected_freq>quantile(Dadtrans_results_all_sim$expected_freq, probs = 0.05))
# Dadtrans_p_value_frq.cut<-as.numeric(quantile(Dadtrans_p_value_frq.cut$p_value, probs = 0.001))
# 
# #father by offspring sex TD
# Dad_off_sextrans_p_value<-as.numeric(quantile(Dad_off_sextrans_results_all_sim$p_value, probs = 0.001))
# 
# Dadtrans_off_sex_p_value_frq.cut<- Dad_off_sextrans_results_all_sim%>%
#   filter(expected_freq>quantile(Dad_off_sextrans_results_all_sim$expected_freq, probs = 0.05))
# Dadtrans_off_sex_p_value_frq.cut<-as.numeric(quantile(Dadtrans_off_sex_p_value_frq.cut$p_value, probs = 0.001))
# 
# #TD by offspring sex
# offspringtrans_p_value<-as.numeric(quantile(offspringtranstrans_results_all_sim$p_value, probs = 0.001))
# 
# offspringtrans_p_value_frq.cut<- offspringtranstrans_results_all_sim%>%
#   filter(expected_freq>quantile(offspringtranstrans_results_all_sim$expected_freq, probs = 0.05))
# offspringtrans_p_value_frq.cut<-as.numeric(quantile(offspringtrans_p_value_frq.cut$p_value, probs = 0.001))
# 
# #make summary table with p_values
# print("Making summary table with p_values ...")
# 
# genedrop_sim_p_value_thresholds<-data.frame(TD_category=c( "Mumtrans", "Mum_off_sextrans", "Dadtrans",
#                                                            "Dad_off_sextrans", "offspringtrans"),
#                                             p_value=c(Mumtrans_p_value, Mum_off_sextrans_p_value, Dadtrans_p_value,
#                                                       Dad_off_sextrans_p_value, offspringtrans_p_value),
#                                             p_value_frq.cut=c(Mumtrans_p_value_frq.cut, Mumtrans_off_sex_p_value_frq.cut,
#                                                               Dadtrans_p_value_frq.cut, Dadtrans_off_sex_p_value_frq.cut,
#                                                               offspringtrans_p_value_frq.cut))
# 
# 
# 
# write.table(genedrop_sim_p_value_thresholds, file=paste0(resultpath, "Genedrop_sim_TD_pvalue_thr.txt"), sep = "\t",
#             col.names = T, row.names = F, quote = F)


#minimum p-values
Mumtrans_min<-as.numeric(min(Mumtrans_results_all_sim$p_value))
Mum_off_sextran_min<-as.numeric(min(Mum_off_sextrans_results_all_sim$p_value))

Dadtrans_min<-as.numeric(min(Dadtrans_results_all_sim$p_value))
Dad_off_sextran_min<-as.numeric(min(Dad_off_sextrans_results_all_sim$p_value))

offspringtrans_min<-as.numeric(min(offspringtranstrans_results_all_sim$p_value))

genedrop_sim_p_value_min<-data.frame(TD_category=c( "Mumtrans", "Mum_off_sextrans", "Dadtrans",
                                                    "Dad_off_sextrans", "offspringtrans"),
                                     p_value=c(Mumtrans_min, Mum_off_sextrans_min, Dadtrans_min,
                                               Dad_off_sextrans_min, offspringtrans_min))

write.table(genedrop_sim_p_value_min, file=paste0(resultpath, "Genedrop_sim_TD_pvalue_min.txt"), sep = "\t",
            col.names = T, row.names = F, quote = F)
