#function to calculate empirical p-values based on genedrop simulation testing for signals
#of directional and balancing selection and plot trajectories of significant SNPs

#snp = character, name of SNP to be tested for selection; path = path to genedrop output containing simulated MAF per cohort; 
#input_object_range = number of genedrop objects each containing a subset of simulation outputs; 
#output_option = choose between both = generating both table with empirical p_values and plot SNPs under selection;
#path_out = output path for plots;
#empirical_p = only calculate table with empirical p_values or plotting = only generate plots of SNPs under selection;
#empirical_MAF_df = data frame containing MAF per cohort for the empirical data set

allele_frq_change_test<-function(snp, path = paste0("results/Genedrop_sim_MAF_per_cohort_part"), input_object_range = c(1:40),
         output_option = "both", path_out = resultpath, plot_name, empirical_p_in = NULL, empirical_MAF_df = phase_chr_MAF_all){
  
  #define output object
  TD_allele_frq_object <- setClass("TD_allele_frq_object",
                              slots = c(empirical_p_all_SNPs_tbl = "data.frame",
                                        p_MAF_SNP_plot = "list"))
  
  if(output_option == "both" || output_option == "empirical_p" ){
    
    #results df
    empirical_p_all_SNPs<-data.frame()
    
    print(paste0("name of SNP ", snp))
    
    MAF_per_cohort_SNP_all_sim<-data.frame()
    i=1
    
    for (p in input_object_range) {
      
      print(paste0("number of Robject containing multiple iterations each: ", p))
      
      load(paste0(path, p, ".RData"))
      
      MAF_per_cohort_SNP_all<-data.frame()
      
      for (df in MAF_per_cohort_all_sim) {
        
        print(paste0("head of df in MAF_per_cohort_all_sim")) 
        print(head(df))
        
        if (length(df)!=0 & class(df)=="data.frame"){
          
          df<-df %>%
            mutate(iteration=i)
          
          print(paste0("iteration ", i))
          
          MAF_per_cohort_SNP<- df %>%
            filter(SNP.Name == snp)
          
          #note whether there are MAFs of 1 in simulated data
          MAF1_match<-ifelse(MAF_per_cohort_SNP$MAF==1, 1, 0)
          MAF1_match<-sum(MAF1_match)
          
          print(paste0("sum of MAF1 = ", MAF1_match))
          
          MAF_per_cohort_SNP$MAF1<-ifelse(MAF1_match!=0, "yes", "no")
          
          MAF_per_cohort_SNP_all<-rbind(MAF_per_cohort_SNP_all, MAF_per_cohort_SNP)
          
        }
        
        i=i+1  
        
      }  
      
      
      #remove iterations with MAF of 1 data points
      MAF_per_cohort_SNP_all<-subset(MAF_per_cohort_SNP_all, !MAF1 == "yes")
      
      MAF_per_cohort_SNP_all_sim<-rbind(MAF_per_cohort_SNP_all_sim, MAF_per_cohort_SNP_all)
      
    }
    
    
    #subset empirical MAF per cohort df by SNP of interest
    phase_chr_MAF_SNP<- empirical_MAF_df %>%
      filter(SNP.name==snp)
    
    phase_chr_MAF_SNP<-phase_chr_MAF_SNP %>%
      mutate(iteration=(max(MAF_per_cohort_SNP_all_sim$iteration)+1), source="empirical")
    
    names(phase_chr_MAF_SNP)[which(names(phase_chr_MAF_SNP)=="SNP.name")]<-"SNP.Name"
    
    MAF_per_cohort_SNP_all_sim<- MAF_per_cohort_SNP_all_sim %>%
      mutate(source="simulation")
    
    n_direct<-0
    n_slope<-0
    n_balance<-0
    
    #make sure data frames are ordered by time point
    phase_chr_MAF_SNP<-arrange(phase_chr_MAF_SNP, Cohort)
    
    #calculate total change over whole time period for empirical data
    phase_change<-abs(phase_chr_MAF_SNP$MAF[nrow(phase_chr_MAF_SNP)]-phase_chr_MAF_SNP$MAF[1])
    
    #calculate magnitude of slope for allele frequency change over whole time period for empirical data
    reg_model<-lm(MAF ~ Cohort, data = phase_chr_MAF_SNP)
    phase_slope<-abs(summary(reg_model)$coefficients[2])
    
    #calculate total covered distances for empirical time series
    phase_dist<-sum(abs(diff(phase_chr_MAF_SNP$MAF)))
    
    
    for(it in unique(MAF_per_cohort_SNP_all_sim$iteration)){
      
      #subset simulation MAF estimates by iteration
      MAF_per_cohort_SNP_sim_sub<-subset(MAF_per_cohort_SNP_all_sim, iteration==it)
      
      #make sure data frames are ordered by time point
      MAF_per_cohort_SNP_sim_sub<-arrange(MAF_per_cohort_SNP_sim_sub, Cohort)
      
      #calculate total change over whole time period for simulated data
      simulated_change<-abs(MAF_per_cohort_SNP_sim_sub$MAF[nrow(MAF_per_cohort_SNP_sim_sub)]-MAF_per_cohort_SNP_sim_sub$MAF[1])
      
      #calculate magnitude of slope for allele frequency change over whole time period for simulated data
      reg_model_sim<-lm(MAF ~ Cohort, data = MAF_per_cohort_SNP_sim_sub)
      simulated_slope<-abs(summary(reg_model_sim)$coefficients[2])
      
      #calculate total covered distances for simulated time series
      simulated_dist<-sum(abs(diff(MAF_per_cohort_SNP_sim_sub$MAF)))
      
      if(simulated_change >= phase_change){
        n_direct<-n_direct +1
      } else if (simulated_dist <= phase_dist) {
        n_balance<-n_balance+1
      } else if (simulated_slope >= phase_slope) {
        n_slope<-n_slope+1
      }
      
    }
    
    
    empirical_p_df<-data.frame( SNP.Name = snp, count_direct = n_direct, count_balance=n_balance, count_slope=n_slope,
                                total_MAF_change=phase_change, total_MAF_dist=phase_dist, MAF_slope=phase_slope,
                                p_direct = (n_direct+1)/(length(unique(MAF_per_cohort_SNP_all_sim$iteration))+1), 
                                p_balance = (n_balance+1)/(length(unique(MAF_per_cohort_SNP_all_sim$iteration))+1),
                                p_slope = (n_slope+1)/(length(unique(MAF_per_cohort_SNP_all_sim$iteration))+1) )
    
    
    empirical_p_all_SNPs<-rbind(empirical_p_all_SNPs, empirical_p_df)
    
   
  }
  
  
  if(is.null(empirical_p_in) == TRUE){
    
    SNPs_empirical_p_sig<- as.data.frame(empirical_p_all_SNPs %>%
                                           filter(p_direct<0.05 | p_slope<0.05 | p_balance<0.05))
                                           
    SNPs_empirical_p_sig<-SNPs_empirical_p_sig$SNP.Name
    
  } else {
    
    SNPs_empirical_p_sig<- as.data.frame(empirical_p_in %>%
                                           filter(p_direct<0.05 | p_slope<0.05 | p_balance<0.05))
    
    SNPs_empirical_p_sig<-SNPs_empirical_p_sig$SNP.Name  
   
    
  }
  
  plot_list<-list()
  
  for (snp_sig in SNPs_empirical_p_sig) { 
    
    if (output_option == "both" || output_option == "plotting"){
      
      #subset empirical MAF per cohort df by SNP of interest
      phase_chr_MAF_SNP<- empirical_MAF_df %>%
        filter(SNP.name==snp_sig)
      
      MAF_per_cohort_SNP_all_sim2<-data.frame()
      j=1
      
      for (pa in input_object_range) {
        
        print(paste0("for plots: number of Robjects containing multiple iterations each: ", pa))
        
        load(paste0(path, pa, ".RData"))
        
        MAF_per_cohort_SNP_all2<-data.frame()
        
        for (df2 in MAF_per_cohort_all_sim) {
          
          if (length(df2)!=0 & class(df2)=="data.frame"){
            
            df2<-df2 %>%
              mutate(iteration=j)
            
            print(paste0("for plots:iteration ", j))
            
            MAF_per_cohort_SNP2<- df2 %>%
              filter(SNP.Name == snp_sig)
            
            #note whether there are MAFs of 1 in simulated data
            MAF1_match2<-ifelse(MAF_per_cohort_SNP2$MAF==1, 1, 0)
            MAF1_match2<-sum(MAF1_match2)
            
            print(paste0("sum of MAF1 = ", MAF1_match2))
            
            MAF_per_cohort_SNP2$MAF1<-ifelse(MAF1_match2!=0, "yes", "no")
            
            MAF_per_cohort_SNP_all2<-rbind(MAF_per_cohort_SNP_all2, MAF_per_cohort_SNP2)
            
          }
          
          j=j+1  
          
        }  
        
        
        #remove iterations with MAF of 1 data points
        MAF_per_cohort_SNP_all2<-subset(MAF_per_cohort_SNP_all2, !MAF1 == "yes")
        
        MAF_per_cohort_SNP_all_sim2<-rbind(MAF_per_cohort_SNP_all_sim2, MAF_per_cohort_SNP_all2)
        
      }
      
       phase_chr_MAF_SNP<-phase_chr_MAF_SNP %>%
         mutate(iteration=(max(MAF_per_cohort_SNP_all_sim2$iteration, na.rm = T)+1), source="empirical")
      
      names(phase_chr_MAF_SNP)[which(names(phase_chr_MAF_SNP)=="SNP.name")]<-"SNP.Name" 
       
      MAF_per_cohort_SNP_all_sim2<- MAF_per_cohort_SNP_all_sim2 %>%
        mutate(source="simulation")
      
      MAF_per_cohort_SNP_all_sim2$MAF1<-NULL
      
      MAF_per_cohort_SNP_all_sim2<-rbind(MAF_per_cohort_SNP_all_sim2, phase_chr_MAF_SNP)
      
      p_MAF_SNP<-ggplot(MAF_per_cohort_SNP_all_sim2, aes(Cohort, MAF, group=iteration))+
        geom_line(aes(color=source, size=source))+
        scale_color_manual(values = c("red", "black"))+
        scale_size_manual(values= c(rep(0.8, (max(MAF_per_cohort_SNP_all_sim2$iteration)-1)), 2)) +
        scale_x_discrete(breaks = MAF_per_cohort_SNP_all_sim2$Cohort[seq(1,length(MAF_per_cohort_SNP_all_sim2$Cohort), by=2)]) +
        ylim(0, (max(MAF_per_cohort_SNP_all_sim2$MAF)+0.05))+
        theme_bw()+
        theme(legend.position = "none",
              axis.text.x = element_text(size=16, colour = "black", angle = -45, vjust=-0.4, margin = margin(t=15)),
              axis.text.y = element_text(size=20, colour = "black"),
              axis.title.x = element_text(size=35, margin = margin(t=20)),
              axis.title.y = element_text(size=35, margin = margin(r=20)),
              axis.line = element_line(colour = "black"),
              axis.ticks=element_line(size=1))
      
      ggsave(paste0(path_out, plot_name, "_allele_frq_change_sig_", snp_sig, ".png"), p_MAF_SNP, 
             width = 25, height = 25, units = "cm")
      
      plot_list[[snp_sig]]<-p_MAF_SNP 
      
    }
    
    
  }

  
  if(output_option == "both"){
    
    allele_frq_change_out<-new("TD_allele_frq_object",
                               empirical_p_all_SNPs_tbl = empirical_p_all_SNPs,
                               p_MAF_SNP_plot = plot_list)
    
  } else if(output_option == "empirical_p"){
    
    plot_list<-list()
    
    allele_frq_change_out<-new("TD_allele_frq_object",
                               empirical_p_all_SNPs_tbl = empirical_p_all_SNPs,
                               p_MAF_SNP_plot = plot_list) 
    
  } else {
    
    empirical_p_all_SNPs<-data.frame()
    
    allele_frq_change_out<-new("TD_allele_frq_object",
                               empirical_p_all_SNPs_tbl = empirical_p_all_SNPs,
                               p_MAF_SNP_plot = plot_list) 
    
  }
  
  
  
}
