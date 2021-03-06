#' Set parameters for LOESS and K-fold cross-validation.
#' 
#' Set parameters for LOESS and the K-fold cross-validation procedure in tafa.
#' @param degree The degree of the locally-fitted polynomial to use. Default is 1 (i.e., linear)
#' @param folds The number of random partitions (folds) to use in K-fold CV. Default it 10.
#' @param iteration The number of repetitions of K-fold CV to perform. Default is 10.
#' @param loss The loss function (prediction error metric from the fit) to evaluate. Default is mean absolute deviation (MAD); alternative is mean squared error (MSE).
#' @param setspan Optional; can define the smoothing parameter for LOESS to use and bypass the K-fold CV. Values can be greater than 0 and <= 1.
#' @export
#' @seealso \code{\link[stats]{loess}}, \code{\link{tafa}}
#' 
#' @details 
#' This control function is akin to how \code{\link[stats]{loess.wrapper}} is
#' used with \code{\link[stats]{loess}}. Here, certain parameters of \code{\link[stats]{loess}}
#' are set as well as specifics of the K-fold cross-validation (CV) procedure.
#' Repeated K-fold CV is used by tafa to estimate f_opt, the optimized
#' smoothing parameter for LOESS. For details on K-fold CV, I recommend
#' Kohavi (1995). 
#' 
#' 10 iterations of 10-fold CV is the default setting for tafa, which performs
#' well with water quality data (Simpson and Haggard, 2016). The user has the 
#' option to specify mean squared error (MSE) as the loss function of prediction
#'  though it may too harshly penalize certain fits due to the squared error
#'  term. Mean absolute deviation (MAD) seems to perform better (and so it's 
#' the default), though it has not been fully studied (there is extensive 
#' discussion in the literature on when various loss functions are more 
#' appropriate to use). 
#' 
#' The user may also specify the smoothing parameter for LOESS to use via 
#' the \code{setspan} argument. Past water quality studies have defaulted
#'  to using a value of 0.5, which generally performs well (simpson and
#'  Haggard, 2016).
#' 
#' @references 
#' Kohavi, R. 1995. A study of cross-validation and bootstrap for accuracy
#' estimation and model selection. International Joint Conference on 
#' Artificial Intelligence 14(2): 1137-1145.
#' 
#' Simpson, Z.P. and B.E. Haggard. 2016. An optimized procedure for 
#' flow-adjustment of constituent concentrations for trend analysis.
#' In preparation.



tafa.control <- function(degree = 1, folds = 10, iteration = 10, loss = c("MAD","MSE"), setspan = NULL){
  list(degree = degree, folds = folds, iteration = iteration, loss=match.arg(loss), setspan = setspan)
}