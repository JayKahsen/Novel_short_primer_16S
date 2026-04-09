get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)))
  }
  getwd()
}

source(file.path(get_script_dir(), "_globalStuff.R"))
################################################################################
# Description
################################################################################
script_title='Figure_4A_4B_6A_6B_abundance'
set_output(script_title)
################################################################################
# Running info and switches
################################################################################

make_new_data_files=testing='no'

make_new_data_files='yes'
#use_full_data='yes'

number_size=1

starting_r=1
ending_r=2

#testing='yes';testing_r=2;number_size=3
#testing='yes';testing_r=1;number_size=1

################################################################################
# function variables
################################################################################
big_constant_number=1000
sc=1000
x_min=2/sc;x_alt=10

comparison_set<-c('primer_label') 
################################################################################
# common Theme and gPlot
################################################################################
sc_log_labels <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    # } else if (val == 1*sc) {
    #   return(expression(10^0))
    } else {val=val/sc
   # return(as.expression(bquote(10^.(round(log10(val))))))
    return(paste0("<b>10<sup>", (round(log10(val))), "</sup></b>"))
    }
  })
}

log2_labels_bold <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    } else if (val == 1) {
      return(expression(2^0))
    } else {val=val
    return(as.expression(bquote(2^.(round(log2(val))))))
    }
  })
}

log2_labels_bold <- function(x) {
  # Requires theme(axis.text.y = element_markdown())
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    } else {
      exponent <- round(log2(val))
      return(paste0("<b>2<sup>", exponent, "</sup></b>"))  # Apply bold to 10 and superscript to exponent
    }
  })
}

# annotate_text_size=6
# text_size=16

magnify=1

margin_size=10
annotate_text_size=6*magnify
axis_text_size=16*magnify
axis_title_text_size=18*magnify
strip_text_size=20*magnify
legend_text_size=20*magnify
rel_x_text_size=1.7


common_theme=  theme(
  strip.text = element_text(face='bold',size = strip_text_size),
  strip.text.y = element_text(face='bold',size = strip_text_size),
  plot.margin = margin(margin_size, margin_size, margin_size, margin_size),
  plot.background = element_rect(fill = color_palette['general_background']),
  legend.background = element_rect(fill = color_palette['general_background']),
  legend.title = element_blank(),
  legend.text = element_text(size = legend_text_size),
  axis.title.y = element_text(face='bold',size = axis_title_text_size,margin = margin(r = 10)),
  axis.title.x = element_text(face='bold',size = axis_title_text_size,margin = margin(t = 10)),
  plot.title = element_blank(),
  axis.text.x=element_markdown(face='bold',size=rel(rel_x_text_size)),
  #axis.text.y.right=element_markdown(face='bold',size=rel(rel_x_text_size)),
  plot.caption = element_text(hjust = 0.5))


gPlot <- function(p) {
  p=p+
    common_theme +
    scale_fill_manual(values = color_palette,guide = guide_legend(reverse = TRUE)) +
    scale_color_manual(values = color_palette,guide = guide_legend(reverse = TRUE))+
    labs(
      title = '',
      x = "",
      y = '',
      caption = ''  
    )
  print(p)
  p  
}



custom_strip_text <- function(strip_colors) {
  theme(
    strip.text = element_text(color = function(label) strip_colors[label])
  )
}

adjust_labels <- function(x) {
  scales::trans_format("log", math_format(exp(.x)))(x)
}
custom_breaks <- function(x) {
  max_number=max(as.numeric(x))
  breaks=unique(x)
  return(breaks)
}
custom_labels <- function(x) {
  labels=c((rep('',(length(x)-1))),paste(length(x),taxa_plural))
  
  return(labels)
}

custom_breaks_labels <- function(x) {
  unique_x <- unique(x)  # Ensure breaks are unique values
  breaks <- tail(unique_x, 1)  # Get the last unique value if needed
  
  qPrint(breaks)
  qPrint(labels)
  
  labels <- length(unique_x)
  assign("labels_from_breaks", labels, envir = .GlobalEnv)
  qPrint(breaks)
  qPrint(labels)
  qPrint(labels_from_breaks)
  return(breaks)
}
################################################################################
# prepping matrix
################################################################################

prepare_df=function(df=matrix_df){
  df1=df%>%
    as.data.frame()%>%
    xPlode_sample_name()%>%
    filter(Type %in% sample_types)%>%
    pivot_longer(col=any_of(feature_names),names_to='feature',values_to='counts')%>%
    left_join(feature_labels_df %>% select(feature, taxa))%>%
    ungroup()%>%
    mutate(Domain = sapply(taxa, extract_domain))%>%
    # calculate differential / relative abundance
    group_by(sample_name)%>%
    mutate(sample_counts = sum(counts)) %>%
    group_by(sample_name,feature)%>%
    mutate(feature_counts = sum(counts)) %>%
    group_by(sample_name,Domain)%>%
    mutate(Domain_counts = sum(counts)) %>%
    ungroup()%>%   
    mutate(rel_abun=counts/sample_counts)%>%
    mutate(domain_diff_abun=Domain_counts/sample_counts)%>%
    mutate(feature_diff=feature_counts/sample_counts)%>%
    group_by(Type,feature)%>%
    mutate(Type_feature_counts=sum(counts))%>%
    filter(Type_feature_counts>0)%>%
    ungroup()
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
  
  if (r==1){number_size=1}
  if (r==2){number_size=3}

  
  taxa_levs=matrix_names[r,'taxa_levs']
  primer_set=matrix_names[r,'primer_sets']
  plot_set=paste0(taxa_levs,'_',primer_set)
  taxa_plural=matrix_names[r,'taxa_plural']
  
  p_title=paste0(script_title,'_',taxa_levs,'_',primer_set)
  
  qPrint(p_title)
  qPrint(primer_set)
  
  
  matrix_df=read.csv(matrix_names[r,'file_path'],check.names=FALSE,row.names=1)%>%
    t()%>%as.data.frame()%>%
    xPlode_sample_name()%>%
    filter(Type %in% sample_types)%>%
    imPlode_sample_name()%>%
    mutate_all(as.numeric)%>%
    select(which(colSums(.) > 0))
  
  ##############################################################################
  # Make data files
  
  ##############################################################################
  # get names ofall features in data set
  feature_names=names(matrix_df)
  
  # time check
  start_time <- Sys.time()
  print('start_time')
  print(start_time)
  
  # compare features
  kw_group_results=Kruskal_Mann_Whitney_Test(matrix_df,comparison_set,'Type')
  
  
  
  ################################################################################
  # create dataframes
  ################################################################################
  
  feature_labels_df=feature_labels(df = matrix_df)
  
  #normed_df=matrix_df/rowSums(matrix_df)
  
  standard_df=prepare_df(matrix_df)
  
  ################################################################################
  print('plotting preperation')
  ################################################################################
  
  feature_labels_df=feature_labels('Type', df = matrix_df)
  
  st_df1=standard_df%>%
    select(-setdiff(meta_data, "sample_name"))%>%
    left_join(meta)%>%
    group_by(feature,Type)%>%
    mutate(Type_count=sum(counts))%>%
    ungroup()%>%
    mutate(Type=factor(Type, levels = rev(sample_type_order)))%>%
    mutate(primer_label=factor(primer_label, levels = primer_label_order))%>%
    mutate(Type=factor(Type, levels = (sample_type_order)))%>%
    left_join(feature_labels_df%>%select(-c(taxa,Domain)), by = c("feature", "Type"))%>%
    group_by(feature, Domain,Type,primer_label) %>%
    summarize(group_mean_counts = mean(counts),
              mean_rel_abun = mean(rel_abun), .groups = 'drop')%>%
    ungroup%>%
    left_join(feature_labels_df%>%select(-c(taxa,Domain)), by = c("feature", "Type"))%>%
    mutate(Type=factor(Type, levels = (sample_type_order)))%>%
    mutate(primer_label=factor(primer_label, levels = primer_label_order))%>%
    ungroup()
  
  mN_df1=st_df1%>%
    group_by(Type,feature)%>%
    mutate(mean_counts=mean(group_mean_counts))%>%
    ungroup()%>%
    mutate(mean_normalized_counts=group_mean_counts/mean_counts)%>%
    #select(ASV,Domain,Type,mean_normalized,plot_order,ASV_label2)%>%
    left_join(kw_group_results, by = c("feature", "Type")) %>%
    ungroup()
  
  set_name=''
  for(i in 1:2){
    if(i==1){
      plotting_set=c('Skin','Soil')
      set_name='_SS'
    }else{
      if(primer_set=='Pool'){next}
      plotting_set=c('Feces','Wastewater')
      set_name='_FW'
    }
    st_df=st_df1%>%
      mutate(Type=factor(Type, levels = (sample_type_order)))%>%
      mutate(primer_label=factor(primer_label, levels = rev(primer_label_order)))%>%
      filter(Type%in% plotting_set)%>%
      arrange(Type)
    mN_df=mN_df1%>%
      mutate(Type=factor(Type, levels = (sample_type_order)))%>%
      mutate(primer_label=factor(primer_label, levels = rev(primer_label_order)))%>%
      filter(Type%in% plotting_set)%>%
      arrange(Type)
    
    ################################################################################
    print('begin plotting')
    #########################################################################################
    ###############################################################################
    # strip backgrounds with outlines
    ################################################################################
    
    unique_y_labels=unique(st_df$Type)
    
    outline_color_y='black'
    
    s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = color_palette["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = color_palette["'),'"]))')
    eval(parse(t=s))
    
    ################################################################################
    print('bar plots')
    ################################################################################
    
    st_plot=ggplot(st_df, aes(x = group_mean_counts*sc, y = plot_order)) +
      geom_col(aes(fill = primer_label),width=.8,position = position_dodge(width = 0.8))+
      #geom_boxplot(aes(color=Type))+
      scale_y_discrete(breaks = custom_breaks, labels = custom_labels)+
      facet_grid2(Type~.,
                  strip = strip_themed(
                    background_y = backgrounds_y,
                    text_y=elem_list_text(
                      face=c(rep('bold',length(unique_y_labels))),
                      size=c(rep(strip_text_size,length(unique_y_labels))))),
                  scale='free', space = "free")+
      scale_x_log10(labels = sc_log_labels,breaks = log_breaks)+
      labs(
        title = '',
        x = "",
        y = paste(taxa_plural,'ordered by average relative abundance'),
        caption = ''  )+
      ggtitle(taxa_levs)
    st_plot
    
    mN_plot=ggplot(mN_df, aes(x = mean_normalized_counts, y = plot_order, fill = Domain)) +
      geom_col(aes(fill = primer_label),width=.8,position = position_dodge(width = 0.8))+
      geom_vline(aes(xintercept=1),color='magenta',linetype='dashed')+
      facet_grid2(Type~.,
                  strip = strip_themed(
                    background_y = backgrounds_y,
                    text_y=elem_list_text(
                      face=c(rep('bold',length(unique_y_labels))),
                      size=c(rep(strip_text_size,length(unique_y_labels))))),
                  scale='free', space = "free")+
      scale_x_continuous(trans = "log2",labels = log2_labels_bold,expand = expansion(mult = c(0, 0.2))) +
      scale_y_discrete(breaks = custom_breaks, labels = custom_labels,position='right')+
      common_theme +
      scale_fill_manual(values = color_palette) +
      scale_color_manual(values = color_palette)+
      labs(
        title = '',
        x = "",
        y = paste(taxa_plural,'ordered by average relative abundance'),
        caption = ''  )
    mN_plot
    
    if(taxa_levs=='Phylum'){
      
      y_labels <- setNames(st_df$feature_label2, st_df$plot_order)
      
      st_plot=st_plot+
        geom_text(
          aes(
            label = signif(group_mean_counts, 3),
            x = ifelse(group_mean_counts < x_min, x_alt, group_mean_counts * sc),
            color = ifelse(group_mean_counts < x_min, 'black', 'white'),
            group = primer_label
          ),
          size = number_size,
          hjust = 1.1,
          fontface = 'bold',
          show.legend = FALSE,
          position = position_dodge(width = 0.8))+
        scale_y_discrete(labels = y_labels)+ 
        theme(
          axis.text.y = element_markdown(face='bold',size=axis_text_size))+
        labs(
          title = '',
          x = "mean counts",
          #y = paste(taxa_plural,'ordered by average relative abundance'),
          caption = ''  )
      
      y_labels <- setNames(mN_df$feature_label2, mN_df$plot_order)
      
      mN_plot=mN_plot+
        geom_text(aes(x=Inf,label=adjusted_significance,color=significance),color='gray40',hjust=1.1,size=annotate_text_size)+
        scale_y_discrete(labels = y_labels,position='right')+ 
        theme(
          axis.text.y.right = element_markdown(face='bold',size=axis_text_size),
          legend.position = "none")+  
        labs(
          title = '',
          x = "RATIO: normalized by the /mean",
         # y = paste(taxa_plural,'ordered by average relative abundance'),
          caption = ''  )
      
    }else{
      mN_plot=mN_plot+
        scale_y_discrete(breaks = custom_breaks, labels = custom_labels,position='right')+ 
        theme(
           axis.text.y.right = element_markdown('bold',size=axis_text_size),
          legend.position = "none")  
    }
    
    st_plot=gPlot(st_plot)
    mN_plot=gPlot(mN_plot)
    combined_plot=grid.arrange(st_plot,mN_plot,nrow=1,widths=c(1.2,1))
    
    plot_height=48
    if(taxa_levs=='Phylum'){
      ggsave(paste0(output_plot,taxa_levs,'_',script_title,'_',primer_set,set_name,'.png'),combined_plot,width=16,height=18)
    }else{
      ggsave(paste0(output_plot,taxa_levs,'_',script_title,'_',primer_set,set_name,'.png'),combined_plot,width=28,height=48)
    }
    
    #if(primer_set=='Pool'){next}
  }
}
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################
# Citation for R
citation()

# Citation for compositions
citation("compositions")

# Citation for vegan
citation("vegan")

# Citation for tidyverse
citation("tidyverse")

# Citation for ggplot2
citation("ggplot2")

# Citation for ggtext
citation("ggtext")





