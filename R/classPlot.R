classPlot <- function(dfc,
                      coordsgrps,
                      metric,
                      metric2,
                      samp){

  metric <- sym(metric)
  metric2 <- sym(metric2)
  
  #--- sample cells based on 'samp' parameter ---#
  
  dfc <- dfc %>%
    na.omit %>%
    group_by(class) %>%
    slice_sample(prop = samp)
  
  #--- plot individual cells coloured by associated class ---#

  p <- ggplot(mapping = aes(x = !!metric, y = !!metric2))+
    geom_point(data = dfc, alpha = 0.3, aes(color = as.factor(class)))+
    theme_bw() +
    theme(legend.position = "none")
  
  #--- add class boundary boxes to delineate class extents ---#

  for(i in 1:nrow(coordsgrps)){
    
    data <- coordsgrps$data[[i]]
    
    p <- p + geom_rect(data = data,
                       aes(xmin = min(!!metric),
                           xmax = max(!!metric),
                           ymin = min(!!metric2),
                           ymax = max(!!metric2)),
                       colour = "black",
                       fill = NA)
  }

  return(p)
}
