# This function needs to be rewritten. A lot of redundand code!


#' Calculate bias, RMSE, corellation etc of two vectors.
#'
#' @description Calculates: corelation coeficient,
#'  bias (absolute and relative), RMSE (absolute and relative), and performes statistical test to check if the differences between \code{reference} and \code{estimate} are significant.
#'  Correlation coefficient as well as statistical test depend on the \code{dist.normal} value.
#'  Optionally statistics can be calculated for groups of data, specified with the \code{by} parameter.
#' @param reference a vector of reference values
#' @param estimate a vector of estimated values
#' @param by Optional grouping variable
#' @param noinfo Logical. Should the additional information on the calculations be displayed?
#' @param dist.normal Logical. Is the distribution normal? If TRUE t-test and Pearson correlation is calculated. If FALSE, Wilcoxon paired test and Spearman corelation is calculated.
#' @param short Logical. Calculate only a subset of summary statistics.
#' @return summary statistics of the differences between \code{reference} and \code{estimate}.
#' @export

calc.error<-function(reference,estimate,by=NULL,noinfo=TRUE,dist.normal=TRUE,short=FALSE) {

  #check input variables
  if (!is.numeric(reference) | !is.numeric(estimate)) {
    stop("Input values not numeric!")
  }

  #check if input values are equal length
  if(length(reference) != length(estimate)) {
    stop("Input variable not of equal length")
  }

  if(is.null(by) == FALSE) {
    if(length(reference) != length(by) | length(reference) != length(by)) {
      stop("Input variable not of equal length")
    }
  }


  if(noinfo==FALSE) { #display information about how the differnce is calculated
    message("The difference is calculated with following equations:") ;
    message(paste("bias =",deparse(substitute(estimate)),"-",deparse(substitute(reference))));
    message(paste("bias% = (",deparse(substitute(estimate)),"-",deparse(substitute(reference)),") / ",deparse(substitute(reference)),"* 100"))
    if (dist.normal == T) {
      message("Correlation method: Pearson. Stat test: t-test (paired)")
    } else {
      message("Correlation method: Spearman Stat test: Wilcoxon (paired)")
    }
  }
  d<-estimate-reference

  #stats in groups
  if(is.null(by) == FALSE) {
    df <- data.frame(reference,estimate,by)  #create temporary data frame
    f <- levels(as.factor(df$by))

    #create placeholders for results
    count<-cor_coeff<-bias<-RMSE<-bias_perc<-RMSE_perc<-p_val<-stat<-c()

    for (i in 1:length(f)) {
      subs <- subset(df,df$by==f[i])
      a <- calc_error_internal(subs,dist.normal = dist.normal, short = short)

      count[i] <- a[2]
      cor_coeff[i] <- a[1]
      bias[i] <- a[3]
      bias_perc[i] <- a[4]
      RMSE[i] <- a[5]
      RMSE_perc[i] <- a[6]
      stat[i] <- a[7]
      p_val[i] <- a[8]
    }

    stat_for_all<-calc.error(reference,estimate,noinfo=TRUE,dist.normal=dist.normal) #recursive!:) #stats for all data to add as the last row
    out<-data.frame(factor=f,count,cor_coeff,bias,bias_perc,RMSE,RMSE_perc,stat,p_val)
    out<-rbind(out,stat_for_all)
    return(out)

  } else {  #calculations for ungrouped data
    f<-"all"
    df <- data.frame(reference,estimate)
    a <- calc_error_internal(df,dist.normal = dist.normal, short = short)

    count <- a[2]
    cor_coeff <- a[1]
    bias <- a[3]
    bias_perc <- a[4]
    RMSE <- a[5]
    RMSE_perc <- a[6]
    stat <- a[7]
    p_val <- a[8]

    return(data.frame(factor=f,count,cor_coeff,bias,bias_perc,RMSE,RMSE_perc,stat,p_val))
  }
}


calc_error_internal <- function(dfs, dist.normal = dist.normal, short = short) {
  count_int     <- length(dfs$reference)
  bias_int      <- mean(dfs$estimate - dfs$reference,na.rm=T)
  RMSE_int      <- rmse(dfs$estimate,dfs$reference)
  bias_perc_int <- bias_int / mean(dfs$reference, na.rm=T)*100
  RMSE_perc_int <- RMSE_int / abs(mean(dfs$reference,na.rm=T))*100

  if (is.infinite(bias_perc_int)) bias_perc <- mean((dfs$estimate - dfs$reference) / ifelse(dfs$reference==0,NA,dfs$reference*100),na.rm=T)


  if (dist.normal == T) {
    if (short==FALSE) {test <- t.test(x = dfs$reference,y = dfs$estimate,paired=T)}
    if (short==FALSE) {p_val <- test$p.value; stat <- as.numeric(test$statistic)} else {p_val<-NA; stat<-NA}
    if (short==FALSE) {cor_coeff<-cor(dfs$reference,dfs$estimate,method="pearson",use="pairwise.complete.obs")} else {cor_coeff <- NA}
  } else {
    if (short==FALSE) {test <-wilcox.test(dfs$reference,dfs$estimate,paired=T)}
    if (short==FALSE) {p_val<-test$p.value; stat<-as.numeric(test$statistic)} else {p_val<-NA; stat<-NA}
    if (short==FALSE) {cor_coeff<-cor(dfs$reference,dfs$estimate,method="spearman",use="pairwise.complete.obs")} else {cor_coeff <- NA}
  }

  return(c(cor_coeff,count_int, bias_int, bias_perc_int, RMSE_int,RMSE_perc_int, stat,p_val))
}


#' round a number to nearest given integer
#'
#' @param x numeric value, can be a vector.
#' @param base integer to round to. Default = 10
#' @return \code{x} rounded to nearest \code{base}.
#' @examples
#' mround(11.32)
#' mround(seq(1,10,by=0.5),2)
#' @export

mround <- function(x,base=10){
  base*round(x/base)
}


#' Calculate root mean square error
#'
#' @param x a single value or vector of observed (reference) values
#' @param y a single value of vector of the same length as \code{x} of predicted values
#' @return root meas square error
#' @examples
#' rmse(x,y)
#' @export

rmse <- function(x,y) {
  sqrt(mean((x - y)^2, na.rm = TRUE))
}


#' Calculate most frequent value of a vector
#'
#' @param x a vector
#' @return mode value of \code{x}
#' @examples
#' mode(rnorm(20))
#' @export


mode <- function(x) {
  if(is.numeric(x)) {
    x <- round(x,2)
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#' Calculate R2
#'
#' @param obs a vector of reference value
#' @param pred a vector of estimated (predicted) values
#' @return R2 value
#' @export

r_square <- function(obs,pred) {
  resid <- obs - pred
  1 - (var(resid, na.rm = TRUE) / var(obs, na.rm = TRUE))
}

#' Custom ggplot template
#'
#' @param base_size Optional font size.
#' @param base_family Optional font name.
#' @return A ggplot template.
#' @aliases bw
#' @export


theme_pt <- function(base_size = 12, base_family = "") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) + #%+replace%
    ggplot2::theme(
      #panel.border = ggplot2::element_rect(colour = "black", fill = F, size = 1),
      axis.text = ggplot2::element_text(margin = ggplot2::margin(10,10,10,10)),
      plot.margin = grid::unit(c(1.2, 1.2, 1.2, 1.2), "lines"),
      axis.ticks.length= grid::unit(0.15,"cm"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = F, size = 1),
      legend.key = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(10,0,0,0)),#element_text(hjust=0.5,vjust=0.5),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(0,10,0,0),angle=90),#(hjust=0.5,vjust=1.5,angle=90),
      #strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(lineheight=1.5),
      strip.background = ggplot2::element_blank(),#no border for facet titles
      legend.position="bottom", # legend on bottom
      legend.title = ggplot2::element_blank () #no title for legend
    )
}

bw<-scatter:::theme_pt()
