
#' @param graph_mtx The graph matrix --> Inv Cov
#' @param thres The thereshold for showing edges according to partial correlation matrix
#' @param regionname This is the Patient ID for plot
#' @return A lot of things
#' 
#' 
#' 

pacman::p_load(spdep,fpc,GWmodel,sp,ppcor,spatstat,reshape2, dplyr,
               tmap,ggplot2,tidyverse, Matrix, igraph, ggraph, matrixcalc)

Network_VisualizationSB <- function( graph_mtx, dinverse, thres, regionname){
  # thres<- 0.1
  # regionname <- 'PATID_201_SLIDE_1'
  # graph_mtx<-PATID_201_SLIDE_1_ICOV[['inv']]
  # rmyname<-colnames(analyzeddata[[1]]$Shared_Net)
  # graph_mtx = my_icov;
  # thres = 0.1;
  # regionname = mypatID ;
  
  
  set.seed(1234)
  ggplotnames <- colnames(graph_mtx);
  
  # change to parcor
  graph_mtx1 <- abs(-diag(1 / sqrt(diag(graph_mtx)))  %*% graph_mtx %*% diag(1 / sqrt(diag(graph_mtx))) )
  diag(graph_mtx)<-0

  # processing the inverse covariance matrix to plot Graph
  # upper triangular matrix
  cor_list <- reshape2:::melt.matrix(graph_mtx) %>%  filter(Var1 != Var2)
  feature_f <- unique(cor_list$Var1)
  feature_f <- feature_f[order(feature_f)]
  cor_list$Var1 <- factor(cor_list$Var1, levels = feature_f)
  cor_list$Var2 <- factor(cor_list$Var2, levels = feature_f)
  
  # [1] added 5 Nov 2026: current dplyr requires one value per group from summarise(),
  # and in a symmetric matrix both entries of a pair tie
  abs_max_fcn <- function(x) { return(x[which(abs(x) == max(abs(x)))[1]]) }
  abs_min_fcn <- function(x) { return(x[which(abs(x) == min(abs(x)))[1]]) }
  
  cor_list <- cor_list %>% mutate(ind_max = pmax(as.numeric(Var1),as.numeric(Var2)),
                                  ind_min = pmin(as.numeric(Var1),as.numeric(Var2))) %>%
              group_by(ind_max, ind_min) %>%
              summarise(cor_max = abs_max_fcn(value),
              cor_min = abs_min_fcn(value),
              feature1 = Var1[which.min(as.numeric(Var1))],
              feature2 = Var2[which.max(as.numeric(Var2))],
              .groups = "drop") %>% distinct(feature1, feature2, .keep_all = TRUE) %>% suppressWarnings()
  
  cor_phi <- select(cor_list,-ind_max,-ind_min,-cor_min)
  
  cor_phi$cor_max <- ifelse(cor_phi$cor_max > 1,1,cor_phi$cor_max)
  cor_phi$cor_max <- ifelse(cor_phi$cor_max < -1,-1,cor_phi$cor_max)
  colnames(cor_phi) <- c("Correlation", "feature1", "feature2" )
  
  cor_phi_list <- cor_phi %>% filter(abs(Correlation) >= thres)
  
  sum_f1 <- cor_phi_list %>% group_by(feature1) %>% summarise(cor_sum1 = sum(abs(Correlation)))
  sum_f2 <- cor_phi_list %>% group_by(feature2) %>% summarise(cor_sum2 = sum(abs(Correlation)))
  colnames(sum_f1)[1] = colnames(sum_f2)[1] = "feature"
  sum_all <- full_join(sum_f1,sum_f2)
  sum_all[is.na(sum_all)] <- 0
  sum_all <- mutate(sum_all, node_size = cor_sum1 + cor_sum2) %>% select(feature,node_size) %>%
    mutate(id = 1:nrow(sum_all))
  
  cor_phi_list <- select(cor_phi_list, "feature1","feature2","Correlation")
  colnames(cor_phi_list) <- c("from","to","Correlation")
  
  g <- igraph::graph_from_data_frame(cor_phi_list, directed=FALSE, vertices=sum_all)
  node_range = c(min(sum_all$node_size), max(sum_all$node_size))
  edge_abs_range = c(min(abs(cor_phi_list$Correlation)), max(abs(cor_phi_list$Correlation)))
  node_size <- sum_all$node_size
  
 
  p3 <- ggraph::ggraph(g, layout = "circle") +
    ggraph::geom_edge_arc(aes(color = factor(sign(Correlation), levels = c(1,-1)),
                              width = abs(Correlation)),strength = 0.2,alpha = 0.75)+
    scale_edge_width(name = "Magnitude", breaks = c(0.2,0.4,0.6,0.8),
                     range = 1*edge_abs_range)+
    scale_edge_color_discrete(name = "Sign", labels = c("Positive", "Negative")) +
    #scale_edge_color_gradient2(   low = "grey", mid = "goldenrod", high = "red", midpoint = 1.5) +
    geom_node_point(aes(size=node_size), color= 'blue') +
    scale_size(name = "Connectivity Degree", range = 2*node_range) +
    theme_void() +  coord_fixed() +
    geom_node_text(aes(label=name,
                       x = x * 1.15, y = y* 1.15,angle = ifelse(atan(-(x/y))*(180/pi) < 0,
                                                                90 + atan(-(x/y))*(180/pi),270 + atan(-x/y)*(180/pi)),
                       hjust = ifelse(x > 0, 0 ,1)), size=2.6, alpha=16, fontface = "bold")+
    expand_limits(x = c(-2, 2), y = c(-2, 2))+
    guides(size = guide_legend(order = 1))+
    theme(text = element_text(family = "Calibri", #size = 14,
                              face = "bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0.5),
          legend.title=element_text(size = 6,face = "bold"),
          legend.text = element_text(size = 5),legend.position = "bottom",
          legend.box="vertical",legend.margin = margin(0, 0, -0.5, -0.5),
          legend.spacing.x = unit(0, "mm"), legend.spacing.y = unit(0, "mm")) +
    labs(title = " Slide Specific Networks",
         subtitle = paste0(' ', regionname, ' ')) 
  # labs(caption = paste0("Gene connectivity with thres", thres)) 
  p3
  
  return(list('graph'=g, 'node_range' = node_range, 'edge_abs_range' = edge_abs_range, 
              'node_size' = node_size, 'cor_phi' = cor_phi, 'cor_list' = cor_list, 'myplot' = p3))
  
}
