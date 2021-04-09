classPlot <- function(dfc,
                      coordsgrps,
                      var1,
                      var2,
                      samp){

  var1 <- sym(var1)
  var2 <- sym(var2)
  
  #--- sample cells based on 'samp' parameter ---#
  
  dfc <- dfc %>%
    na.omit %>%
    group_by(class) %>%
    slice_sample(prop = samp)
  
  #--- plot individual cells coloured by associated class ---#

  p <- ggplot(mapping = aes(x = !!var1, y = !!var2))+
    geom_point(data = dfc, alpha = 0.3, aes(color = as.factor(class)))+
    theme_bw() +
    theme(legend.position = "none")
  
  #--- add class boundary boxes to delineate class extents ---#

  for(i in 1:nrow(coordsgrps)){
    
    data <- coordsgrps$data[[i]]
    
    p <- p + geom_rect(data = data,
                       aes(xmin = min(!!var1),
                           xmax = max(!!var1),
                           ymin = min(!!var2),
                           ymax = max(!!var2)),
                       colour = "black",
                       fill = NA)
  }

  return(p)
}
