#function for manhattan plots using ggplot; option to split by sex
#need data frame with columns SNP.name, Chromosome and p_value; variable link_order needs to contain SNP names and SNP order only
#if plot is split by parent (for transmission distortion) column name must be trans_parent

manhattan_plot<-function(df_plot, link_order, y.var = "p_value", chr_no, p_threshold, y_scales = NULL, xlab_text = "Chromosome", ylim_max=6, x_text_size=18, 
                         y_text_size=22, point_size=3, alpha_value=0.5, chr_colours=c("steelblue3", "red3")){
  
  colnames(link_order)<-c("SNP.name", "order")
  #freq.cutoff<-data.frame(freq.quantile=as.numeric(quantile(df_plot$expected_freq, probs = 0.05)))
  #df_plot<-subset(df_plot, expected_freq > freq.cutoff$freq.quantile)
  df_plot<-join(df_plot, link_order)
  df_plot[ ,"order"]<-as.numeric(df_plot[ ,"order"])
  df_plot[ ,"Chromosome"]<-as.numeric(df_plot[ ,"Chromosome"])
  df_plot<-arrange(df_plot, Chromosome, order)
  cum<-nrow(df_plot)
  df_plot[ ,"pos_cuml"]<-c(1:cum)
  
  #make separate x axis for chromosomes and add title variable to label plot with a facet title including an expression
  axisdf <- plyr::ddply(df_plot,.(Chromosome), summarize, center=(max(pos_cuml) + min(pos_cuml))/2)
  
  label_seq<-axisdf$Chromosome[seq(1,length(axisdf$Chromosome),by=2)]
  breaks_seq<-axisdf$center[seq_along(axisdf$center)%% 2 > 0]
  
  if(length(grep("Sex|sex", colnames(df_plot)))!=0){
    
    df_plot<-subset(df_plot, !sex=="unknown")
    
    #make custom facet labels
    sex.labs<-c("A. Females", "B. Males")
    names(sex.labs)<-c("female", "male")
    
    if (y.var == "p_value"){
    
    p_manhattan<-ggplot2::ggplot(df_plot, aes(pos_cuml, -log10(p_value))) +
      geom_point(aes(color=as.factor(Chromosome)), size=point_size, alpha=alpha_value)+
      scale_colour_manual(values = rep(chr_colours, chr_no))+
      scale_x_continuous(label= label_seq, breaks = breaks_seq)+ 
      scale_y_continuous(expand = c(0,0))+
      ylim(0,ylim_max)+ 
      geom_hline(yintercept =  -log10(p_threshold), linetype="dashed", color="black", size=1)+
      xlab(xlab_text)+ylab(expression(-log[10](p)))+
      theme_bw()+
      theme(legend.position = "none", axis.text.x = element_text(size=x_text_size, colour = "black"),
            axis.text.y = element_text(size=y_text_size, colour = "black"), 
            axis.title.x = element_text(size=35, margin = margin(t=20)), 
            axis.title.y = element_text(size=35, margin = margin(r=20)),
            axis.line = element_line(colour = "black"),
            axis.ticks=element_line(size=1), strip.text = element_text(size=22, hjust=0.01),
            strip.background = element_rect(colour="black", fill="lightgray"))+
      facet_wrap(vars(sex), nrow = 2, ncol = 1,labeller = labeller(sex=sex.labs))
    
    } else if (y.var == "effect_size"){
      
      p_manhattan<-ggplot2::ggplot(df_plot, aes(pos_cuml, effect_size)) +
        geom_point(aes(color=as.factor(Chromosome)), size=point_size, alpha=alpha_value)+
        scale_colour_manual(values = rep(chr_colours, chr_no))+
        scale_x_continuous(label= label_seq, breaks = breaks_seq)+ 
        scale_y_continuous(limits = c(0, NA))+
        xlab(xlab_text)+ylab("Effect size (probability)")+
        theme_bw()+
        theme(legend.position = "none", axis.text.x = element_text(size=x_text_size, colour = "black"),
              axis.text.y = element_text(size=y_text_size, colour = "black"), 
              axis.title.x = element_text(size=30, margin = margin(t=20)), 
              axis.title.y = element_text(size=30, margin = margin(r=20)),
              axis.line = element_line(colour = "black"),
              axis.ticks=element_line(size=1), strip.text = element_text(size=22, hjust=0.01),
              strip.background = element_rect(colour="black", fill="lightgray"))+
        facet_wrap(vars(sex), nrow = 2, ncol = 1,labeller = labeller(sex=sex.labs), scales = "free_y")
        #facet_grid_sc(rows = vars(sex), scales = list(y = scales_y), labeller = labeller(sex=sex.labs))
      
    }
    
    
  } else if (length(grep("trans_parent", colnames(df_plot)))!=0){
    
    #make custom facet labels
    parent.labs<-c("A. Females", "B. Males")
    names(parent.labs)<-c("mother", "father")
    
    parent_levels<-c("mother", "father")
    df_plot$parent.ordered<-factor(df_plot$trans_parent, levels = parent_levels)
    
    p_manhattan<-ggplot2::ggplot(df_plot, aes(pos_cuml, -log10(p_value))) +
      geom_point(aes(color=as.factor(Chromosome)), size=point_size, alpha=alpha_value)+
      scale_colour_manual(values = rep(chr_colours, chr_no))+
      scale_x_continuous(label= label_seq, breaks = breaks_seq)+ 
      scale_y_continuous(expand = c(0,0))+
      ylim(0,ylim_max)+ 
      geom_hline(yintercept =  -log10(p_threshold), linetype="dashed", color="black", size=1)+
      xlab(xlab_text)+ylab(expression(-log[10](p)))+
      theme_bw()+
      theme(legend.position = "none", axis.text.x = element_text(size=x_text_size, colour = "black"),
            axis.text.y = element_text(size=y_text_size, colour = "black"), 
            axis.title.x = element_text(size=35, margin = margin(t=20)), 
            axis.title.y = element_text(size=35, margin = margin(r=20)),
            axis.line = element_line(colour = "black"),
            axis.ticks=element_line(size=1), strip.text = element_text(size=22, hjust=0.01),
            strip.background = element_rect(colour="black", fill="lightgray"))+
      facet_wrap(vars(parent.ordered), nrow = 2, ncol = 1, labeller=labeller(parent.ordered=parent.labs))
    
    
  } else {
    
    if (y.var == "p_value"){
    
    p_manhattan<-ggplot2::ggplot(df_plot, aes(pos_cuml, -log10(p_value))) +
      geom_point(aes(color=as.factor(Chromosome)), size=point_size, alpha=alpha_value)+
      scale_colour_manual(values = rep(chr_colours, chr_no))+
      scale_x_continuous(label= label_seq, breaks = breaks_seq)+ 
      scale_y_continuous(expand = c(0,0))+
      ylim(0,ylim_max)+ 
      geom_hline(yintercept =  -log10(p_threshold), linetype="dashed", color="black", size=1)+
      xlab(xlab_text)+ylab(expression(-log[10](p)))+
      theme_bw()+
      theme(legend.position = "none", axis.text.x = element_text(size=x_text_size, colour = "black"),
            axis.text.y = element_text(size=y_text_size, colour = "black"), 
            axis.title.x = element_text(size=35, margin = margin(t=20)), 
            axis.title.y = element_text(size=35, margin = margin(r=20)),
            axis.line = element_line(colour = "black"),
            axis.ticks=element_line(size=1))
    
    } else if (y.var == "effect_size"){
      
      p_manhattan<-ggplot2::ggplot(df_plot, aes(pos_cuml, effect_size)) +
        geom_point(aes(color=as.factor(Chromosome)), size=point_size, alpha=alpha_value)+
        scale_colour_manual(values = rep(chr_colours, chr_no))+
        scale_x_continuous(label= label_seq, breaks = breaks_seq)+ 
        scale_y_continuous(expand = c(0,0))+
        ylim(0,ylim_max)+ 
        xlab(xlab_text)+ylab(expression(-log[10](p)))+
        theme_bw()+
        theme(legend.position = "none", axis.text.x = element_text(size=x_text_size, colour = "black"),
              axis.text.y = element_text(size=y_text_size, colour = "black"), 
              axis.title.x = element_text(size=35, margin = margin(t=20)), 
              axis.title.y = element_text(size=35, margin = margin(r=20)),
              axis.line = element_line(colour = "black"),
              axis.ticks=element_line(size=1))
      
    }
    
    
  }
  
 
  p_manhattan
  
}