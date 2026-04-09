source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)),'/_globalStuff.R'))
###############################################################################
# Description
################################################################################
script_title='Figure_2_5_alpha_diversity_stats'
set_output(script_title)
options(scipen = 999) # don;t use scientific notation
################################################################################
# Running info and switches
################################################################################

make_new_data_files=testing='no'

make_new_data_files='yes'

# testing='yes';testing_r=7
# testing='yes';testing_r=8
starting_r=7
ending_r=8
################################################################################
# Extra Script Variables and functions
################################################################################

metric_order=c('percent_Muri','percent_Lachno','percent_Propi','percent_Archaea','richness', 'evenness', 'shannon','n')

################################################################################
# Extra Script Variables and functions
################################################################################
################################################################################
# plotting settings
################################################################################
strip_text_size=12
x_text_size=10
margin_size=10

pool_height=10
pool_width=7
ind_height=10
ind_width=8

p_value_size=3
geom_text_size=3
ind_text_size=3
common_theme <- theme(
 # legend.position = "bottom",
  plot.margin = margin(margin_size, margin_size, margin_size, margin_size),
  plot.title = element_blank(),
  axis.title = element_blank(),
  plot.background = element_rect(fill = color_palette['general_background']),
  axis.text.y=element_text(face='bold'),
  strip.text = element_text(size = 11)
  
)

gPlot <- function(p) {
  p=p+
    scale_color_manual(values = color_palette,labels=palette_label)+
    common_theme+
    labs(
      color='',
    #  title = taxa_levs,
      title = '',
      x = '',
      y = '',
      #   caption = ''
    )
  #print(p)
}

################################################################################
# loop through all matrix files 
################################################################################

if(testing=='yes'){ starting_r=ending_r=testing_r}
if(!exists('ending_r')){ending_r=n_rows=nrow(matrix_names)}
if(!exists('starting_r')){starting_r=1}
# load matrix from matrix_names data
################################################################################
for (r in starting_r:ending_r){ # main loop for different data sets
  ################################################################################

  
  
  taxa_levs=matrix_names[r,'taxa_levs']
  primer_set=matrix_names[r,'primer_sets']
  p_title=paste0(taxa_levs,'_',script_title,'_',primer_set)
  
  
  if(primer_set=='Ind'){  geom_text_size=ind_text_size; sample_types=sample_type_order}
  if(primer_set=='Pool'){  geom_text_size=ind_text_size; sample_types=c('Skin','Soil')}
  ################################################################################
  if (make_new_data_files=='yes'){ # create data files for plotting
    ################################################################################
    matrix_df=read.csv(matrix_names[r,'file_path'],check.names=FALSE,row.names=1)%>%
    #  matrix_df=read.csv(matrix_names[r,'rarefied_path'],check.names=FALSE,row.names=1)%>%
      t()%>%as.data.frame()
    
    # define variables for loaded matrix
    feature_names=names(matrix_df)

    ################################################################################
    # get percent Archaea and Propi/cutibacterium
    ################################################################################
    
    df2=matrix_df%>%
      xPlode_sample_name()%>%
      pivot_longer(col=all_of(feature_names),names_to='feature',values_to='counts')%>%
      group_by(feature,Type)%>%
      mutate(feature_Type_counts=sum(counts))%>%
      filter(feature_Type_counts>0)%>%
      mutate(taxa=feature)%>%
      group_by(sample_name) %>%
      mutate(sample_counts=sum(counts))%>%
     
     #filter(Type=='Fecal')%>%
      ungroup()
    
    if(matrix_names[r,"taxa_levs"]=='ASV'){
      
      df2=df2%>%
        rename(ASV=taxa)%>%
        left_join(ASV_taxa%>%select(ASV,taxa,Family))%>%
      ungroup()
    }
    names(df2)
    df_frequency=df2 %>% 
      select(sample_name,primer_label,Type,Exp,feature,counts) %>% 
      filter(counts >= 1 & counts <= 20)%>%
      mutate(values=counts) %>% 
      group_by(primer_label,Type,Exp) %>%
      dplyr::count(values)%>%
      rename(frequency = n) %>%
      ungroup()
    names(df_frequency)
    
    ggplot(df_frequency,aes(x = values,y=frequency))+
      facet_grid2(Exp+Type~ primer_label,scales='free')+
      geom_col()
    
    ggsave(paste0(output_plot,p_title,'_frequency.png'))
    
    Domain_df <- df2 %>%
      mutate(Domain = sapply(taxa, extract_domain)) %>%
      group_by(sample_name, Domain) %>%
      summarise(Domain_counts = sum(counts),
                sample_counts = mean(sample_counts),
                .groups = 'drop')%>%
      mutate(Domain_percentage=100*Domain_counts/sample_counts)%>%
      ungroup()%>%
      filter(Domain=='Archaea')%>%
      rename(smp_percent_Archaea=Domain_percentage)%>%
      select(sample_name,smp_percent_Archaea)
    Propi_df <- df2 %>%
      mutate(Propi = sapply(taxa, extract_propi)) %>%
      group_by(sample_name, Propi) %>%
      summarise(Propi_counts = sum(counts),
                sample_counts = mean(sample_counts),
                .groups = 'drop')%>%
      mutate(Propi_percentage=100*Propi_counts/sample_counts)%>%
      ungroup()%>%
      filter(Propi=='Propi')%>%
      rename(smp_percent_Propi=Propi_percentage)%>%
      select(sample_name,smp_percent_Propi)
    Lachno_df <- df2 %>%
      mutate(Lachno = sapply(taxa, extract_lachno)) %>%
      group_by(sample_name, Lachno) %>%
      summarise(Lachno_counts = sum(counts),
                sample_counts = mean(sample_counts),
                .groups = 'drop')%>%
      mutate(Lachno_percentage=100*Lachno_counts/sample_counts)%>%
      ungroup()%>%
      filter(Lachno=='Lachno')%>%
      rename(smp_percent_Lachno=Lachno_percentage)%>%
      select(sample_name,smp_percent_Lachno)
    
    Muri_df <- df2 %>%
      mutate(Muri = sapply(taxa, extract_muri)) %>%
      group_by(sample_name, Muri) %>%
      summarise(Muri_counts = sum(counts),
                sample_counts = mean(sample_counts),
                .groups = 'drop')%>%
      mutate(Muri_percentage=100*Muri_counts/sample_counts)%>%
      ungroup()%>%
      filter(Muri=='Muri')%>%
      rename(smp_percent_Muri=Muri_percentage)%>%
      select(sample_name,smp_percent_Muri)
    
    ################################################################################
    # Generate alpha diversity metrics
    ################################################################################
    
    end_count=1
    
    metrics_df <- matrix_df %>%
      as.data.frame() %>%
      xPlode_sample_name()%>%
      rowwise()%>%
      mutate(
        mn = mean(c_across(any_of(feature_names))),
        sd=sd(c_across(any_of(feature_names))),
        smp_n=sum(c_across(any_of(feature_names))),
        smp_shannon=qShannon(c_across(any_of(feature_names))),
        smp_evenness=qEvenness(c_across(any_of(feature_names))),
        smp_richness = qRichness(c_across(any_of(feature_names))),
      )%>%
      select(-any_of(feature_names))%>%
      full_join(Domain_df)%>%
      full_join(Propi_df)%>%
      full_join(Lachno_df)%>%
      full_join(Muri_df)%>%
      group_by(Type,Primer)%>%
      mutate(across(starts_with("smp_"), ~mean(., na.rm = TRUE), .names = "{sub('smp_', 'mean_', .col)}"))%>%
      rename_with(~gsub("smp_", "", .), starts_with("ratio_smp_"))%>%
      ungroup()
    
    write.csv(metrics_df,paste0(output_data,p_title,'.csv'),row.names = FALSE)
    
    ################################################################################
  } # skipping above if (make_new_data_files!='yes')
  ################################################################################
  
  ################################################################################
  print('prepare for plotting') # beginning plotting section
  ################################################################################
  
  ################################################################################ 
  # read data files from save
  if (make_new_data_files!='yes'){metrics_df=read.csv(paste0(output_data,p_title,'.csv'),check.names = FALSE)}
  ################################################################################

  #################################################################################
  print('# changing to long form')
  #################################################################################
  #stop()
  result_df=data_df=metrics_df%>%
    pivot_longer(cols = starts_with("smp_"), names_to = "metric")%>%
    filter(
      (Type %in% c('Skin') & metric != 'smp_percent_Muri'& metric != 'smp_percent_Archaea') |
        (Type %in% c('Soil', 'Wastewater') & metric != 'smp_percent_Propi'& metric != 'smp_percent_Muri') |
        (Type %in% c('Feces') & metric != 'smp_percent_Propi'& metric != 'smp_percent_Archaea') 
    )%>%
    filter(metric!='smp_percent_Lachno')%>%
    group_by(Type,Primer,metric)%>%
    mutate(mean_value=mean(value))%>%
    ungroup()
  
  #################################################################################
  print('# T-test if comparing only 2 primers across combinations of Type and metric')
  #################################################################################
  unique_primers=as.character(unique(data_df$Primer))
  qPrint(unique_primers)
  qPrint(length(unique_primers))
  
  if(length(unique_primers)==2){
    t_test_results =data_df%>%
      #filter((Type=='Skin' & metric!='smp_percent_Archaea')|(Type=='Soil' & metric!='smp_percent_Propi')|Type %in% c('Fecal','Water'))%>%
      group_by(Type, metric) %>%
      # Perform a t-test within each group
      summarise(t_test_summary = {
        if (n_distinct(value) == 1) {
          list(data.frame(p.value = NA, statistic = NA))
        } else {
          list(broom::tidy(wilcox.test(value ~ Primer)))
          #list(broom::tidy(wilcox.test(value ~ Primer, exact = FALSE)))
        }
      },.groups = 'drop') %>%
      # Unnest the t-test summary for clearer results
      unnest(t_test_summary)
    t_test_df=t_test_results%>%
      select(Type,metric,p.value)
    
    # Join t-test results back to the original data
    result_df <- data_df %>%
      left_join(t_test_df, by = c("Type", "metric"))%>%
      #result_df <-t_test_results%>%
      rename(p_value='p.value')%>%
      mutate(metric = factor(metric, levels = rev(paste0("smp_", metric_order))))
  }
  
  ################################################################################
  plot_list <- list()
  
  for (t in sample_types){ # setting up separate plots by Type 
    ################################################################################ 
    
    qPrint(sample_types)
    qPrint(t)
  
    df=result_df%>%
      filter(Type==t)%>%
      mutate(metric = factor(metric, levels = rev(paste0("smp_", metric_order))))
    
    ###############################################################################
    print('# strip backgrounds with outlines')
    ################################################################################
    unique_x_labels=unique(df$Type)
    
    outline_color_x='black'
    
    s=paste0('backgrounds_x <- list(element_rect(color=outline_color_x,fill = color_palette["',paste0(unique_x_labels,collapse='"]),element_rect(color=outline_color_x,fill = color_palette["'),'"]))')
    print(s)
    eval(parse(t=s))
    backgrounds_x
    
    unique_y_labels =c(rep(unique_x_labels,(length(unique(df$metric))/length(unique(df$Type)))))
    
    outline_color_y='black'
    
    s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = color_palette["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = color_palette["'),'"]))')
    eval(parse(t=s))
    
    face_y=face=c(rep('bold',5))
    if(t=='Skin'){face_y=face=c(rep('bold',4),'bold.italic')}
    ################################################################################
    print('# Changing strip labels')
    ################################################################################
    print(metric_order)
    new_metric=c('% Muribaculaceae','% Lachnospiraceae','% Cutibacterium','% Archaea','Richness', 'Evenness', 'Shannon','Reads')
    names(new_metric)=metric_order
    # new_Type=c('Feces','Skin','Soil', 'Wastewater')
    # names(new_Type)=Type_order
   # mutate(Type = case_when(Type %in% names(new_Type) ~ new_Type[Type],TRUE ~ Type))%>%
    
    ################################################################################
    print('# Plotting')
    ################################################################################ 
    
    if(primer_set=='Ind'){
      
      p=df%>%
        mutate(metric = gsub("smp_", "", metric))%>%
        #mutate(Type = case_when(Type %in% names(new_Type) ~ new_Type[Type],TRUE ~ Type))%>%
        mutate(metric = case_when(metric %in% names(new_metric) ~ new_metric[metric],TRUE ~ metric))%>%
        mutate(metric = factor(metric, levels = rev(new_metric)))%>%
        mutate(primer_label = factor(primer_label, levels = primer_label_order))%>%
        mutate(mean_value_label=signif(mean_value,3))%>%
        mutate(mean_value_label=if_else(mean_value_label>1000,round(mean_value,0),mean_value_label))%>%
        
        ggplot(aes(x = primer_label, y = value, color = primer_label)) +
        geom_boxplot(width=.6,outlier.size = .5)+
        geom_text(aes(label = paste0('  ',mean_value_label)),y=0,vjust=.5,hjust=0,color='black', show.legend=FALSE,size=geom_text_size,angle=90)+
        facet_grid2(metric~ Type,
                    strip = strip_themed(background_x = backgrounds_x,background_y = backgrounds_y,
                                         text_x=elem_list_text(
                                           face=c(rep('bold',length(unique_x_labels))),
                                           size=c(rep(strip_text_size,length(unique_x_labels)))),
                                         text_y=elem_list_text(
                                           face=face_y,
                                           size=c(rep(strip_text_size,length(unique_y_labels))))),
                    scale='free')+
        scale_y_continuous(limits=c(0,NA))+
        scale_x_discrete(labels = function(x) {
          x_label <- df$primer_label[match(x, df$primer_label)]
          color <- color_palette[x_label]
          paste0("<span style='color:", color, "'>", x_label, "</span>")
        }) +
        theme(legend.position='none',
          axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5,face='bold',size=x_text_size)  # Rotate and enable markdown for x-axis labels
        )+
        ggtitle(taxa_levs)+
        
        labs(
          title = taxa_levs,
          x = '',
          y = '',
        )
      
      plot_list[[t]] <- gPlot(p)+ theme(legend.position='none')
      
    }
    ############################################################################
    
    if(length(unique_primers)==2){
      
      ############################################################################
      x_labels <- setNames(df$primer_label, df$primer_label)
      breaks=as.character(unique(df$primer_label))
      breaks
      break_labels=palette_label[breaks]
      
      p=df%>%
        mutate(metric = gsub("smp_", "", metric))%>%
        mutate(metric = case_when(metric %in% names(new_metric) ~ new_metric[metric],TRUE ~ metric))%>%
        mutate(metric = factor(metric, levels = rev(new_metric)))%>%
        mutate(p_label=paste0('p-value = ',signif(p_value, 3)))%>%
        mutate(p_label = paste0('p-value = ', format(p_value, scientific = TRUE, digits = 3)))%>%
        mutate(significance=ifelse(p_value<.05,'sig','ns'))%>%
        mutate(p_label = ifelse(metric=='Reads', "", p_label))%>%
        mutate(mean_value_label=signif(mean_value,3))%>%
        mutate(mean_value_label=if_else(mean_value_label>1000,round(mean_value,0),mean_value_label))%>%
        ggplot(aes(x = primer_label, y = value, color = primer_label)) +
        geom_boxplot(width=.6,outlier.size = .5)+
        facet_grid2(metric~ Type,
                    strip = strip_themed(background_x = backgrounds_x,background_y = backgrounds_y,
                                         text_x=elem_list_text(
                                           face=c(rep('bold',length(unique_x_labels))),
                                           size=c(rep(strip_text_size,length(unique_x_labels)))),
                                         text_y=elem_list_text(
                                           face=face_y,
                                           size=c(rep(strip_text_size,length(unique_y_labels))))),
                    scale='free')+
        geom_text(aes(label = p_label, x = 1.5, y = 0 ,color=significance), hjust = .5,vjust = -10.5,size=p_value_size, fontface='bold',show.legend = FALSE)+
        geom_text(aes(y = 0, label = mean_value_label),color='black', vjust = -.5,show.legend=FALSE,size=geom_text_size)+
        scale_y_continuous(limits=c(0,NA))+
        scale_x_discrete(labels = function(x) {
          x_label <- df$primer_label[match(x, df$primer_label)]
          color <- color_palette[x_label]
          paste0("<span style='color:", color, "'>", x_label, "</span>")
        }) +
        common_theme+
        theme(
          axis.text.x=element_blank(),
          axis.ticks.x =element_blank() 
          #axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5,face='bold',size=x_text_size)  # Rotate and enable markdown for x-axis labels
          
        )+
        scale_color_manual(
          values = color_palette,
          labels = palette_label,
          breaks = as.character(unique(df$primer_label))
        )+
        
        labs(
          color='',
          title = taxa_levs,
          x = '',
          y = '',
          #   caption = ''
        )
      p
      plot_list[[t]] <- p
      
    }
  }
  

  ############################################################################
  print('# save plots')
  ############################################################################
  type_plots=grid.arrange(plot_list[['Skin']],plot_list[['Soil']],nrow=1,widths = c(1, 1)) 

  
  if(primer_set=='Pool'){
    type_plots <- wrap_plots(
      plot_list[['Skin']] + theme(legend.position = "bottom"),
      plot_list[['Soil']] + theme(legend.position = "bottom"),
      nrow = 1,
      widths = c(1, 1)
    ) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    ggsave(paste0(output_plot,p_title,'.png'),type_plots,width=pool_width,height=pool_height)
  }else{
    ggsave(paste0(output_plot,p_title,'.png'),type_plots,width=ind_width,height=ind_height)
    type_plots=grid.arrange(plot_list[['Feces']],plot_list[['Wastewater']],nrow=1,widths = c(1, 1)) 
    ggsave(paste0(output_plot,p_title,'supplemental.png'),type_plots,width=ind_width,height=ind_height)
  }
}
############################################################################
print('# End')
############################################################################
