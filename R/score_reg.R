#' Regression of continuous response on the score from multivariate analysis
#'
#' This function performs regression of continuous response on the score from multivariate analysis like (s)3CA and (s)cPCA.
#' @import ggplot2
#' @param x score from multivariate analysis.
#' @param y continous external response.
#' @param plot whether or not plot the regression fit, default is FALSE.
#' @param title title for regression plot.
#' @export
#' @return regression fit, i.e., lm object.




score_reg <- function(x, y, plot = FALSE, title = ""){
  y <- y - mean(y)
  fit <- lm(y~., data = x)
  dat <- as.data.frame(cbind(x, y))
  if(plot){
    p <- ggplot(data = dat)+
      geom_point(aes(x = dat[,1], y = dat[,2], color = dat[,3]), show.legend = F)+
      labs(x = colnames(dat)[1], y = colnames(dat)[2], title = title)
    theme_bw()
    print(p)
  }
  return(fit)
}
