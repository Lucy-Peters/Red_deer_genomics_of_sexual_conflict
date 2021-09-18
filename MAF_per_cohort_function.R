#function to calculate MAF per cohort for a genedrop object

MAF_per_cohort<-function(genedrop_obj, marker_map, IDs_keep){
  
  haplos_sim<-extract_genotype_mat(genedrop_obj) %>%data.frame()
  
  #calculate allele frequencies per cohort
  genedrop_ped<-genedrop_obj@pedigree %>%data.frame()
  colnames(genedrop_ped)<-c("ID", "Father", "Mother", "sex", "Cohort")
  ID_cohorts<-genedrop_ped[, c("ID", "Cohort")] %>%data.frame()
  ID_cohorts<-data.frame(ID=rep(ID_cohorts$ID, each=2), Cohorts=rep(ID_cohorts$Cohort, each=2))
  
  haplos_sim<-cbind(ID_cohorts, haplos_sim)
  
  
  colnames(haplos_sim)[c(-1, -2)]<-c(as.character(marker_map$SNP.Name))
  haplos_sim<-subset(haplos_sim, ID %in% IDs_keep$ID)
  cohorts<-unique(haplos_sim$Cohorts)
  #dummy_founders<-data.frame(ID=grep("^Dummy.*", haplos_sim$ID, value=T))
  

  haplos_sim_MAF_all<-data.frame()
  
  for (co in cohorts) {
    haplos_sim_sub<-subset(haplos_sim, Cohorts==co)
    #haplos_sim_sub<-subset(haplos_sim_sub, !ID %in% dummy_founders$ID)
    hap.origin.names<-as.vector(c("paternal", "maternal"))
    haplos_sim_sub$hap.origin<-rep(hap.origin.names, (nrow(haplos_sim_sub)/2))
    haplos_sim_sub$Cohorts<-NULL
    haplos_sim_melt<- reshape2::melt(haplos_sim_sub, id.vars = c("ID", "hap.origin"), variable.name="SNP.Name")
    
    print(paste0(co, " ", length(unique(haplos_sim_melt$ID))))
    
    #exclude missing allele entries
    haplos_sim_melt<-subset(haplos_sim_melt, !value==9)
    
    print(paste0("No IDs after excluding missing alleles: ", length(unique(haplos_sim_melt$ID))))
    
    haplos_sim_freq<-haplos_sim_melt%>%
      group_by(SNP.Name, value)%>%
      summarise(N=n())
    
    haplos_sim_freq<-haplos_sim_freq%>%
      group_by(SNP.Name)%>%
      summarise(total=sum(N))%>%
      left_join(haplos_sim_freq, ., by="SNP.Name")
    
    haplos_sim_MAF<-haplos_sim_freq%>%
      mutate(fraction=N/total)%>%
      filter(value==0)
    
    names(haplos_sim_MAF)[names(haplos_sim_MAF)=="fraction"]<-"MAF"
    names(haplos_sim_MAF)[names(haplos_sim_MAF)=="N"]<-"N_minor"
    haplos_sim_MAF<-as.data.frame(haplos_sim_MAF)
    
    SNPs_MAF_miss<-subset(marker_map, !SNP.Name %in% haplos_sim_MAF$SNP.Name)
    
    if (nrow(SNPs_MAF_miss)!=0){
      haplos_sim_freq_sub<-subset(haplos_sim_freq, SNP.Name %in% SNPs_MAF_miss$SNP.Name)
      haplos_sim_MAF_fixed<-data.frame(SNP.Name=unique(haplos_sim_freq_sub$SNP.Name),
                                       value=0, N_minor=0, total=haplos_sim_freq$total[1],
                                       MAF=1)
      haplos_sim_MAF<-rbind(haplos_sim_MAF, haplos_sim_MAF_fixed)
    }
    
    #keep SNP name, cohort, no of individuals and MAF
    haplos_sim_MAF$N_IDs<-haplos_sim_MAF$total/2
    haplos_sim_MAF$Cohort<-co
    haplos_sim_MAF.dt<-data.table(haplos_sim_MAF)
    marker_map.dt<-data.table(marker_map[c("SNP.Name", "CEL.LG")])
    haplos_sim_MAF<-as.data.frame(marker_map.dt[haplos_sim_MAF.dt, on = c(SNP.Name = "SNP.Name")])
    haplos_sim_MAF<-haplos_sim_MAF[c("SNP.Name", "CEL.LG","MAF", "N_IDs", "Cohort")]
    
    haplos_sim_MAF_all<-rbind(haplos_sim_MAF_all, haplos_sim_MAF)
    
  }
  
  haplos_sim_MAF_all
  
  
  
}


