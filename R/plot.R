#' plot.cs_sampling
#'
#' Produce pairwise plots comparing each parameter distribution before and after adjustment.
#' Takes the output list from cs_sampling as input.
#'
#' @param x - object of type cs_sampling
#' @param varnames - optional vector of names of subset of variable for pairs plotting
#' @param ... - generic parameter
#' @return The output of \code{\link[GGally]{ggpairs}}.
#' @import GGally
#' @import ggplot2
#' @method plot cs_sampling
#' @rdname plot.cs_sampling
#' @export
plot.cs_sampling <- function(x, varnames = NULL, ...) {

  datpl <- data.frame(rbind(as.matrix(x$sampled_parms), as.matrix(x$adjusted_parms))
                      , as.factor(c(rep("NO", dim(x$sampled_parms)[1]), rep("YES", dim(x$adjusted_parms)[1]))))
  names(datpl)[dim(x$sampled_parms)[2]+1] <- c("Adjust")
  rownames(datpl) <- NULL

  #subset to varnames
  if(!is.null(varnames)){datpl <- datpl[, c(varnames, "Adjust")]}

  my_ellipse <- function(data, mapping){
    ggplot2::ggplot(data = data, mapping = mapping) +
      ggplot2::geom_point()+
      ggplot2::stat_ellipse(level = 0.90, type = "norm", size = 2)
  }

  my_violin <- function(data, mapping){
    ggplot2::ggplot(data = data, mapping = mapping) +
      ggplot2::geom_violin(trim=TRUE,draw_quantiles = c(0.05, 0.5, 0.95),alpha=0.5, size = 1.5)
  }

  p1 <- GGally::ggpairs(datpl, mapping = aes(color = Adjust, alpha = 0.5), columns = c(1:(dim(datpl)[2]-1)),
                lower = list(continuous = my_ellipse))
  return(p1)
}
