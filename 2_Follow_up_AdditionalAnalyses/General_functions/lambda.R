
#' Calculate genomic inflation factor lambda
#'
#' @param pvalues vector of p values
#'
#' @return lambda value
#' @export
lambda <- function(pvalues){
  chisq <- qchisq(1-pvalues, 1)
  lambda <- median(chisq)/qchisq(0.5, 1)
  return(lambda)
}