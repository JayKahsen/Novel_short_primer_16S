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
script_title='Figure_4C_6C_ordination'
script_description <- "
Creates the ordination figure panels and their PERMANOVA and PERMDISP annotations.
This repository preserves the published statistical workflow by calculating those
tests on distances from the plotted PCA score space rather than the full CLR matrix.
"

set_output(script_title)
################################################################################
# Running info and switches
################################################################################

testing=make_new_data_files='no'
make_new_data_files='yes'

# testing='yes';testing_r=7
# testing='yes';testing_r=8
starting_r=7
ending_r=8

#testing='yes';testing_r=5;starting_r=1
included_Primer=included_Primers=primer_order

number_of_permutations=1000
#number_of_permutations=2

################################################################################
################################################################################
# common Theme and gPlot
################################################################################
annotate_text_size=6
vjust_margin=1.5
text_size=18
title_size=text_size+2
strip_text_size=title_size


margin_size=20
loading_text_size=2.1

common_theme=  theme(
  plot.title = element_blank(),
  plot.background = element_rect(fill = color_palette['general_background']),
  legend.background = element_rect(fill = color_palette['general_background']),
  legend.title = element_blank(),
  legend.text = element_text(size = title_size),
  axis.title.x = element_markdown(size = title_size,margin = margin(t = 10)),
  axis.title.y = element_markdown(size = title_size),
  axis.text = element_text(size = text_size),
  plot.caption = element_text(hjust = 0.5))



gPlot <- function(p) {
  p=p+
    scale_color_manual(values = color_palette)+
    common_theme 
  p  
}
################################################################################
# generic information and functions
################################################################################
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
  matrix_df=read.csv(matrix_names[r,'file_path'],check.names=FALSE,row.names=1)%>%
    t()%>%as.data.frame()
  
  # define variables for loaded matrix
  feature_names=names(matrix_df)
  taxa_levs=matrix_names[r,'taxa_levs']
  primer_set=matrix_names[r,'primer_sets']
  p_title=paste0(script_title,'_',taxa_levs,'_',primer_set)
  
  # time check
  start_time <- Sys.time()
  print('start_time')
  print(p_title)
  print(start_time)
  
  ################################################################################
  # prepping matrix
  ################################################################################
  
  prep_df=matrix_df%>%
    as.data.frame()%>%
    xPlode_sample_name()%>%
    filter(Type %in% sample_types)%>%
    ungroup()
  
  #initialize
  pcs_df=loadings_df=percent_explained_df=data.frame()
  
  ################################################################################
  # looping by type
  ################################################################################
  permANOVA=permDisp=plot_list <- list()
  Type=unique(prep_df$Type)
  p_vector=1:length(Type)
  names(p_vector)=Type
  f_vector=disp_p_vector=disp_f_vector=p_vector
  
  for (t in Type){
    qPrint(Type)
    qPrint(t)
    
    Type_df=prep_df%>%
      filter(Type==t)%>%
      imPlode_sample_name()%>%
      select(where(~ sum(.) > 0)) 
    
    meta_df=Type_df%>%
      xPlode_sample_name()
    ################################################################################
    print('# distance matrix,mds, points, percentage explained, assemble types')
    ################################################################################
    
    pseudo_count= 0.5*(min(Type_df[Type_df>0]))
    pseudo_count_df=Type_df+pseudo_count
    normalized_df=pseudo_count_df/rowSums(pseudo_count_df)
    data=as.data.frame(clr(normalized_df))
    
    # Perform PCA
    results=pca_result <- prcomp(data, scale. = FALSE)
    summary(results)
    
    pca_result_df=as.data.frame(pca_result$x)%>%
      ungroup
    pca_result_df=as.data.frame(pca_result$rotation)%>%
      ungroup
    
    # Extract the principal components
    pcs <- as.data.frame(pca_result$x)%>%
      xPlode_sample_name()%>%
      rename(x='PC1',y='PC2')%>%
      select(any_of(meta_data),x,y)
    
    # percent explained
    percent_explained=format(round(results$sdev^2 / sum(results$sdev^2)*100,1),nsmall=1,trim=TRUE)%>%
      as.data.frame()%>%
      slice(1:2) %>%
      t()%>%
      as.data.frame() %>%
      `colnames<-`(c('x_lab', 'y_lab'))%>%
      mutate(Type=t)
    
    percent_explained_df=rbind(percent_explained_df,percent_explained)
    # loadings_df=rbind(loadings_df,loadings)
    pcs_df=rbind(pcs_df,pcs)
    
    
    #############################################################################
    print('permANOVA')
    #############################################################################
    
    meta_df=Type_df%>%
      xPlode_sample_name()
    
    qPrint(Type)
    qPrint(t)
    
    # The published version tested dispersion and PERMANOVA on the plotted PCA
    # score space rather than the full CLR matrix.
    pca_scores <- pca_result$x[, 1:2]
    
    dist_matrix <- dist(pca_scores)
    
    # Perform PERMANOVA
    result <- adonis2(dist_matrix ~ Primer,data=meta_df,permutations = number_of_permutations)
    print(result)
    p_value=result$'Pr(>F)'[1]
    f_statistic <- result$F[1]
    qPrint(p_value)
    qPrint(f_statistic)
    
    p_vector[t]=signif(p_value,2)
    f_vector[t]=signif(f_statistic,4)
    qPrint(t)
    qPrint(p_vector)
    
    # View results
    print(result)
    permANOVA[[t]]=result
    #############################################################################
    print('permDISP')
    #############################################################################
    
    permdisp_res <- betadisper(dist_matrix, group = meta_df$Primer)
    
    # View PERMDISP results
    print(permdisp_res)
    
    # PERMDISP annotation is based on the permutation-test result for the same
    # distance matrix used in the ordination and PERMANOVA sections.
    perm_test <- permutest(permdisp_res, permutations = 999)
    print(perm_test)
    
    perm_p<- perm_test$tab[1, "Pr(>F)"]
    print(perm_p)
    perm_f <- perm_test$tab[1, "F"]
    print(perm_f)
    
    # Plot the results
    plot(permdisp_res)
    
    disp_p_vector[t]=signif(perm_p,2)
    disp_f_vector[t]=signif(perm_f,4)
    permDisp[[t]]=perm_test
    
  }
  #############################################################################
  print('plotting begins')
  ##############################################################################
  disp_p_label=disp_f_lapel=p_label=f_label=''
  
  for (t in Type){print(t)
    
    p_label=paste0('PERMANOVA  p = ' ,p_vector[t])
    disp_p_label=paste0('PERMDISP  p = ' ,disp_p_vector[t])
    
    # f_label=paste0('PERMANOVA F = ' ,f_vector[t])
    # disp_f_label=paste0('PERMDISP  F = ' ,disp_f_vector[t])
    
    x_lab=percent_explained_df%>%
      filter(Type==t)%>%
      pull(x_lab)
    
    y_lab=percent_explained_df%>%
      filter(Type==t)%>%
      pull(y_lab)
    
    df=pcs_df%>%
      filter(Type==t)
    
    df$Primer=factor(df$Primer,levels=primer_order)
    
    df_matrix=matrix_df
    Type_feature_labels=feature_labels('Type')%>%
      filter(Type==t)%>%
      select(-Type)
    
    lim_adj=2
    ############################################################################
    
    Group=unique(df$Group)
    
    breaks_asv <- as.character(1:10)
    breaks_Primer = labels_Primer=c(primer_order,Group) 
    breaks_Primer = labels_Primer=c(primer_order,Group) 
    
    shape_breaks=shape_labels=Group
    
    breaks <- c(breaks_asv,primer_order,Group)
    #labels <- c(labels_asv,primer_order,Group)
    ############################################################################
    
    ###############################################################################
    # strip backgrounds with outlines
    ################################################################################
    unique_x_labels=  unique(df$Type)
    
    outline_color_x='black'
    
    s=paste0('backgrounds_x <- list(element_rect(color=outline_color_x,fill = color_palette["',paste0(unique_x_labels,collapse='"]),element_rect(color=outline_color_x,fill = color_palette["'),'"]))')
    print(s)
    eval(parse(t=s))
    backgrounds_x
    
    unique_y_labels =rev(unique(df$Type))
    unique_y_labels
    
    outline_color_y='black'
    
    s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = color_palette["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = color_palette["'),'"]))')
    eval(parse(t=s))
    
    #colors <- color_palette[match(intersect(primer_order,unique(df$Primer)), names(color_palette))]
    
    # unique(df$primer_label)
    
    #############################################################################
    #if(t=='Skin'){stop()}
    
    # df=df%>%
    #   mutate(Type = case_when(Type %in% names(new_Type) ~ new_Type[Type],TRUE ~ Type))%>%
    #   ungroup()
    
    spacer=' '
    if(t %in% c('Feces','Wastewater')){
      spacer='                                                                  '
    }
    
    p=ggplot(df, aes(x = x, y = y)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed")+
      geom_point(aes(color=primer_label))+
      xlab(paste('PCA 1',x_lab,'% Explained'))+
      ylab(paste('PCA 2',y_lab,'% Explained'))+
      stat_ellipse(aes(group = primer_label, color = primer_label),linewidth=.7, level = 0.95, type = "norm") +
      stat_ellipse(data=df%>%filter(primer_label %in% c('StandardO','TruncatedFO','Standard','Truncated','St','Tr')),aes(group = primer_label,  color = primer_label),linewidth=1.4, level = 0.95, type = "norm") +
      
      
      facet_wrap2(~Type, strip.position = "right",
                  strip = strip_themed(
                    background_x = backgrounds_x, background_y =backgrounds_y,
                    text_x=elem_list_text(face=c('bold'),
                                          size=c(rep(strip_text_size,length( unique_x_labels)))),
                    text_y=elem_list_text(face=c('bold'),
                                          size=c(rep(strip_text_size,length( unique_y_labels))))),
                  scale = 'free') +
      ggtitle(paste(t,p_title,'biplot'))+
      theme(legend.position = "bottom",
            plot.margin = margin(margin_size, margin_size, margin_size/2, margin_size))+
      scale_color_manual(values = color_palette)+
      guides(color = guide_legend(override.aes = list(shape = 22,linewidth=8)))+
      
      annotate("text", x = -Inf, y = Inf, label = paste0(spacer,p_label), hjust = 0,vjust=vjust_margin*1,color='orangered',size=annotate_text_size)+
      #annotate("text", x = Inf, y = Inf, label = disp_anova, hjust = 1.1,vjust=vjust_margin*2,color='orangered',size=annotate_text_size)+
      annotate("text", x = -Inf, y = Inf, label = paste0(spacer,disp_p_label), hjust = 0,vjust=vjust_margin*2,color='orangered',size=annotate_text_size)+
      
      # annotate("text", x = -Inf, y = Inf, label = paste0(spacer,f_label), hjust = 0,vjust=vjust_margin*3,color='blue',size=annotate_text_size)+
      # annotate("text", x = -Inf, y = Inf, label = paste0(spacer,disp_f_label), hjust = 0,vjust=vjust_margin*4,color='blue',size=annotate_text_size)+
      
      labs(color =NULL,fill=NULL,shape =NULL,size=NULL)+
      ggtitle(taxa_levs)
    p
    
    gPlot(p)
    # df %>% filter(primer_label == 'Tr') %>% select(x, y) %>% summary()
    
    #############################################################################
    # df_loadings$ASV_label
    ################################################################################
    #############################################################################
    #library(grid)
    
    
    if( t%in% c('Soil','Wastewater')){p=p+theme(legend.position = "none",
                                                plot.margin = margin(margin_size/2, margin_size, margin_size, margin_size))
    }
    plot_list[[t]] <- gPlot(p)
    
  }
  
  
  # ############################################################################
  # print('# save plots')
  # ############################################################################
  plot_width=9
  plot_height=18
  type_plots=grid.arrange(plot_list[['Skin']],plot_list[['Soil']],ncol=1,heights = c(1.05, 1))
  
  if(primer_set=='Pool'){
    ggsave(paste0(output_plot,p_title,'.png'),type_plots,width=plot_width,height=plot_height)
  }else{
    ggsave(paste0(output_plot,p_title,'.png'),type_plots,width=plot_width,height=plot_height)
    type_plots=grid.arrange(plot_list[['Feces']],plot_list[['Wastewater']],ncol=1,heights = c(1.05, 1)) 
    ggsave(paste0(output_plot,p_title,'supplemental.png'),type_plots,width=plot_width,height=plot_height)
  }
  
  
  
}
################################################################################


################################################################################
permANOVA
permDisp
