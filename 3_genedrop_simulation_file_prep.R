#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Red deer: file preparation to run genedrop simulation  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(plyr)
library(dplyr)
library(ggplot2)

#make large tables with genotype and marker information of all chromosomes

phase_all_chr<-data.frame()
marker_map_all_chr<-data.frame()

for (chr in c(1:33)) {

  #AlphaPeel phase files
  assign("phase_chr_in", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaPeel_run/Results/AlphaPeel_out_chr",
                                           chr, ".called_phase.0.95"), header=F))

  #AlphaPeel marker files not required for AlphaPeel but needed to assign SNP names
  assign("marker_map_chr", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/sexual_conflict/Phasing_haplotypes/AlphaPeel_run/Input_files/Deer_AlphaPeel_marker_file_chr",
                                             chr,".txt"), header=T))

  if(chr==1){

    phase_all_chr<-phase_chr_in

  } else {

   phase_all_chr<-cbind(phase_all_chr, phase_chr_in[, -1])

  }

  marker_map_all_chr<-rbind(marker_map_all_chr, marker_map_chr)

}



#add chromosome information to marker map
linkmap<-read.table("data/TableS1_CervusElaphus_Final_Linkage_Map.txt", header=T)
linkmap<-linkmap[,c("SNP.Name", "CEL.LG", "CEL.order", "cMPosition.SexAveraged")]

names(marker_map_all_chr)[which(names(marker_map_all_chr)=="SNP.name")]<-"SNP.Name"
marker_map_all_chr<-join(marker_map_all_chr, linkmap[c("SNP.Name", "CEL.LG")])

write.table(phase_all_chr, file = "results/AlphaPeel_out_all_chr.txt", sep="\t", col.names = F, row.names = F,
            quote=F)
write.table(marker_map_all_chr, file="results/Deer_AlphaPeel_marker_file_all_chr.txt", sep="\t", col.names = T,
            row.names = F, quote=F)


#need to split iterations for simulation because 1000 iterations creates vectors that are too big
#parallelise simulations additionally to mclapply by splitting iterations into different scripts

#write separate scripts for 40 runs of 25 simulation iterations

x<-readLines("scripts/Red_deer_TD_genedrop_simulation.R")

for (i in c(1:40)){
  write(paste0("run =" , "'", "part", i, "'"), paste0("scripts/Red_deer_TD_genedrop_simulation_run", i,".R"))
  write(x, paste0("scripts/Red_deer_TD_genedrop_simulation_run", i,".R"), append = T)
}



