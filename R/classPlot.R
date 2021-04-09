classPlot <- function(dfc,
                      coordsgrps,
                      var1,
                      var2,
                      samp){

  var1 <- sym(var1)
  var2 <- sym(var2)
  
  dfc <- dfc %>%
    na.omit %>%
    group_by(class) %>%
    slice_sample(prop = samp)

  p <- ggplot(mapping = aes(x=!!var1,y=!!var2))+
    geom_point(data = dfc, aes(color = as.factor(class)))+
    theme_bw() +
    theme(legend.position = "none")

  for(i in 1:nrow(coordsgrps)){
    p <- p + geom_rect(data=coordsgrps$data[[i]],
                       aes(xmin=min(!!var1),
                           xmax=max(!!var1),
                           ymin=min(!!var2),
                           ymax=max(!!var2)),
                       colour="black",
                       fill=NA)
  }

  return(p)
}
