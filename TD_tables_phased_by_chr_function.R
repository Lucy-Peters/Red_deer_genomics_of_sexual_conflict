#function to make allele transmission tables
TD_tables_phased_by_chr<-function(genedrop_obj, marker_map, IDs_keep, chr_range){
  
  #define output object
  TD_table_object <- setClass("TD_table_object",
                              slots = c(Mumtrans_tbl = "data.frame",
                                        Dadtrans_tbl = "data.frame",
                                        offspringtrans_tbl = "data.frame",
                                        Mum_off_sextrans_tbl = "data.frame",
                                        Dad_off_sextrans_tbl = "data.frame",
                                        SNP_error_tbl = "data.frame"))
  
  #get genotype matrix
  haplos_sim<-extract_genotype_mat(genedrop_obj) %>%data.frame()
  
  #extract pedigree, cohort and ID information
  genedrop_ped<-genedrop_obj@pedigree %>%data.frame()
  colnames(genedrop_ped)<-c("ID", "Father", "Mother", "sex", "Cohort")
  ID_cohorts<-genedrop_ped[, c("ID", "Cohort")] %>%data.frame()
  ID_cohorts<-data.frame(ID=rep(ID_cohorts$ID, each=2), Cohorts=rep(ID_cohorts$Cohort, each=2))
  
  haplos_sim<-cbind(ID_cohorts["ID"], haplos_sim)
  
  colnames(haplos_sim)[-1]<-c(as.character(marker_map$SNP.Name))
  haplos_sim<-subset(haplos_sim, ID %in% IDs_keep$ID)
  
  Mumtrans_all<-data.frame()
  Mum_off_sextrans_all<-data.frame()
  Dadtrans_all<-data.frame()
  Dad_off_sextrans_all<-data.frame()
  offspringtrans_all<-data.frame()   
  SNP_error_all<-data.frame()
  
  for (chr in chr_range) {
    
    #split map by chromosome
    marker_map_sub<-subset(marker_map, CEL.LG==chr)
    #snps on chromosome
    SNP.names<-c(as.character(marker_map_sub$SNP.Name))
    SNP.names.df<-data.frame(SNP.Name=SNP.names)
    #split phased data by chromosome
    haplos_sim_sub<-haplos_sim[, SNP.names] %>% data.frame()

    #correct SNP names (remove X at start of name in the sheep)
    snp.names_haplos_sim_sub<-colnames(haplos_sim_sub) %>% data.frame
    colnames(snp.names_haplos_sim_sub)<-"SNP.Name"
    snp_matches_sub<-subset(snp.names_haplos_sim_sub, !SNP.Name %in% SNP.names.df$SNP.Name)
    names(haplos_sim_sub)[!names(haplos_sim_sub) %in% SNP.names.df$SNP.Name]<- substring(snp_matches_sub$SNP.Name, 2)
    
    haplos_sim_sub<-cbind(haplos_sim["ID"], haplos_sim_sub)
    
    #treat simulated haplotypes as empirical data set and create TD tables
    phase_chr<-haplos_sim_sub
    
    hap.origin.names<-as.vector(c("paternal", "maternal"))
    phase_chr$hap.origin<-rep(hap.origin.names, (nrow(phase_chr)/2))
    
    phase_chr_melt<- reshape2::melt(phase_chr, id.vars = c("ID", "hap.origin"), variable.name="SNP.Name")
    phase_chr_cast<- reshape2::dcast(phase_chr_melt, SNP.Name+ID~hap.origin, value.var = "value" )
    
    phase_chr_cast.dt<-data.table(phase_chr_cast)
    genedrop_ped.dt<-data.table(genedrop_ped)
    marker_map_sub.dt<-data.table(marker_map_sub[c("SNP.Name", "CEL.LG")])
    phase_chr_cast.dt<-genedrop_ped.dt[phase_chr_cast.dt, on = c(ID="ID")]
    phase_chr_cast<-as.data.frame(marker_map_sub.dt[phase_chr_cast.dt, on = c(SNP.Name="SNP.Name")])
    
    phase_chr_cast$sex<-ifelse(phase_chr_cast$sex==2, "female", 
                               ifelse(phase_chr_cast$sex==1,"male", "unknown"))
    
    #allocate random allele labels 
    
    #phase_chr_random<-data.frame()
    
    #for (snp in SNP.names) {
      #phase_chr_cast_sub<-subset(phase_chr_cast, SNP.Name==snp)
      
      #draw<-rbinom(1, 1, 0.5)
      #if(draw==1){
        
        #AlphaPeel
        #allele 0 = A, allele 1 = B
        #phase_chr_cast_sub$maternal<-ifelse(phase_chr_cast_sub$maternal==0, "A",
                                            #ifelse(phase_chr_cast_sub$maternal==1, "B", NA))
        
        #phase_chr_cast_sub$paternal<-ifelse(phase_chr_cast_sub$paternal==0, "A",
                                            #ifelse(phase_chr_cast_sub$paternal==1, "B", NA))
        
        #phase_chr_cast_sub$minor.allele<-"A"
        
        
        
      #} else {
        
        #AlphaPeel
        #allele 0 = B, allele 1 = A
        #phase_chr_cast_sub$maternal<-ifelse(phase_chr_cast_sub$maternal==0, "B",
                                            #ifelse(phase_chr_cast_sub$maternal==1, "A", NA))
        
        #phase_chr_cast_sub$paternal<-ifelse(phase_chr_cast_sub$paternal==0, "B",
                                            #ifelse(phase_chr_cast_sub$paternal==1, "A", NA))
        
        #phase_chr_cast_sub$minor.allele<-"B"
        
      #}
      
      #phase_chr_random<-rbind(phase_chr_random, phase_chr_cast_sub)
    #}
    
    #assign allele label non-randomly
    #allele 0 = B, allele 1 = A
    phase_chr_cast$maternal<-ifelse(phase_chr_cast$maternal==0, "B",
                                            ifelse(phase_chr_cast$maternal==1, "A", NA))
    phase_chr_cast$paternal<-ifelse(phase_chr_cast$paternal==0, "B",
                                            ifelse(phase_chr_cast$paternal==1, "A", NA))
                
    phase_chr_cast$minor.allele<-"B"    
    
    phase_chr_random<-phase_chr_cast


    phase_chr_random_mum<-phase_chr_random[, c("SNP.Name", "ID", "maternal", "paternal")]
    colnames(phase_chr_random_mum)<-c("SNP.Name", "Mother", "maternal_mother", "paternal_mother")
    
    phase_chr_random.dt<-data.table(phase_chr_random, key=c("SNP.Name", "Mother"))
    phase_chr_random_mum.dt<-data.table(phase_chr_random_mum, key=c("SNP.Name", "Mother"))
    phase_chr_random<-as.data.frame(phase_chr_random_mum.dt[phase_chr_random.dt])
    
    phase_chr_random_dad<-phase_chr_random[, c("SNP.Name", "ID", "maternal", "paternal")]
    colnames(phase_chr_random_dad)<-c("SNP.Name", "Father", "maternal_father", "paternal_father")
    
    phase_chr_random.dt<-data.table(phase_chr_random, key=c("SNP.Name", "Father"))
    phase_chr_random_dad.dt<-data.table(phase_chr_random_dad, key=c("SNP.Name", "Father"))
    phase_chr_random<-as.data.frame(phase_chr_random_dad.dt[phase_chr_random.dt])
    
    
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
      group_by(SNP.Name)%>%
      summarize(error.count=n()) %>%
      na.omit
    
    SNP_error<-as.data.frame(SNP_error)
    
    #exclude any transmissions with Mendelian error
    phase_chr_random_clean<-subset(phase_chr_random_clean, !mend_error=="mat_error" & !mend_error=="pat_error")
    
    #summarise transmission frequencies by transmitting parent and convert into format that can be used for
    #Fisher's exact test (dcast - 1 row per SNP, columns are alleles A and B, values are frequencies )
    
    Mumtrans<-subset(phase_chr_random_clean, informative=="yes_mat"|informative=="yes_both")
    
    Mumtrans<-Mumtrans %>%
      group_by(SNP.Name, maternal)%>%
      summarise(freq=n())%>%
      na.omit
    
    
    Mumtrans <- reshape2::dcast(Mumtrans, SNP.Name ~ maternal)
    
    
    Dadtrans<-subset(phase_chr_random_clean, informative=="yes_pat"|informative=="yes_both")
    
    Dadtrans<-Dadtrans %>%
      group_by(SNP.Name, paternal)%>%
      summarise(freq=n())%>%
      na.omit
    
    Dadtrans <-reshape2::dcast(Dadtrans, SNP.Name ~ paternal)
    
    
    #summarise transmission frequencies by transmitting parent and offspring sex
    Mum_off_sextrans<-subset(phase_chr_random_clean, (informative=="yes_mat"|informative=="yes_both") & !sex=="Z")
    Mum_off_sextrans<-subset(Mum_off_sextrans, !sex=="unknown")
    
    Mum_off_sextrans<-Mum_off_sextrans %>%
      group_by(SNP.Name, maternal, sex) %>%
      summarise(freq=n()) %>%
      na.omit
    
    Mum_off_sextrans<-reshape2::dcast(Mum_off_sextrans, SNP.Name+sex ~ maternal)
    
    Dad_off_sextrans<-subset(phase_chr_random_clean, (informative=="yes_pat"|informative=="yes_both") & !sex=="Z")
    Dad_off_sextrans<-subset(Dad_off_sextrans, !sex=="unknown")
    
    Dad_off_sextrans<-Dad_off_sextrans %>%
      group_by(SNP.Name, paternal, sex) %>%
      summarise(freq=n()) %>%
      na.omit
    
    Dad_off_sextrans<- reshape2::dcast(Dad_off_sextrans, SNP.Name+sex ~ paternal)
    
    
    #offspring TD - TD by offspring sex (Mum and Dad allele togehter)
    otab<-phase_chr_random_clean[c("SNP.Name","ID", "sex", "maternal", "paternal", "informative")]
    
    otab_pat<-subset(otab, informative=="yes_pat")
    otab_pat<-otab_pat[c("SNP.Name", "ID", "sex", "paternal")]
    names(otab_pat)[names(otab_pat) == "paternal"] <- "value"
    
    otab_mat<-subset(otab, informative=="yes_mat")
    otab_mat<-otab_mat[c("SNP.Name", "ID", "sex", "maternal")]
    names(otab_mat)[names(otab_mat) == "maternal"] <- "value"
    
    otab_both<-subset(otab, informative=="yes_both")
    otab_both<-otab_both[c("SNP.Name", "ID", "sex", "maternal", "paternal")]
    otab_both <- reshape2::melt(otab_both, id.vars = c("SNP.Name","ID", "sex"))
    otab_both<-otab_both[c("SNP.Name", "ID", "sex", "value")]
    
    offspringtrans<-rbind(otab_both, otab_mat)
    offspringtrans<-rbind(offspringtrans, otab_pat)
    offspringtrans<-subset(offspringtrans, !sex=="Z")
    
    #offspring transmission
    offspringtrans <- offspringtrans %>%
      group_by(SNP.Name, sex, value) %>%
      summarize(freq=n()) %>%
      na.omit
    
    
    names(offspringtrans)[names(offspringtrans) == "value"] <- "allele"
    offspringtrans <- reshape2::dcast(offspringtrans, SNP.Name + sex ~ allele) #dcast = opposite of melt
    
    
    #add chromosome minor allele and mendelian error information
    phase_chr_random_clean_single<-subset(phase_chr_random_clean, !duplicated(SNP.Name))
    Mumtrans<-join(Mumtrans, phase_chr_random_clean_single[c("SNP.Name", "CEL.LG", "minor.allele")], by="SNP.Name", type="inner")
    Mumtrans<-join(Mumtrans, SNP_error, by="SNP.Name")
    Dadtrans<-join(Dadtrans, phase_chr_random_clean_single[c("SNP.Name", "CEL.LG", "minor.allele")], by="SNP.Name")
    Dadtrans<-join(Dadtrans, SNP_error, by="SNP.Name")
    offspringtrans<-join(offspringtrans, phase_chr_random_clean_single[c("SNP.Name", "CEL.LG", "minor.allele")], by="SNP.Name")
    offspringtrans<-join(offspringtrans, SNP_error, by="SNP.Name")
    Mum_off_sextrans<-join(Mum_off_sextrans, phase_chr_random_clean_single[c("SNP.Name", "CEL.LG", "minor.allele")], by="SNP.Name")
    Mum_off_sextrans<-join(Mum_off_sextrans, SNP_error, by="SNP.Name")
    Dad_off_sextrans<-join(Dad_off_sextrans, phase_chr_random_clean_single[c("SNP.Name", "CEL.LG", "minor.allele")], by="SNP.Name")
    Dad_off_sextrans<-join(Dad_off_sextrans, SNP_error, by="SNP.Name")
    
    
    Mumtrans_all<-rbind(Mumtrans_all, Mumtrans)
    Mum_off_sextrans_all<-rbind(Mum_off_sextrans_all, Mum_off_sextrans)
    Dadtrans_all<-rbind(Dadtrans_all, Dadtrans)
    Dad_off_sextrans_all<-rbind(Dad_off_sextrans_all, Dad_off_sextrans)
    offspringtrans_all<-rbind(offspringtrans_all, offspringtrans)   
    SNP_error_all<-rbind(SNP_error_all, SNP_error)
    
    print(paste0("Done tables for chromosome ", chr))
    
  }
  
  
  
  TD_tables_out<-new("TD_table_object",
                     Mumtrans_tbl = Mumtrans_all,
                     Dadtrans_tbl = Dadtrans_all,
                     offspringtrans_tbl = offspringtrans_all,
                     Mum_off_sextrans_tbl = Mum_off_sextrans_all,
                     Dad_off_sextrans_tbl = Dad_off_sextrans_all,
                     SNP_error_tbl = SNP_error_all)
  
}

