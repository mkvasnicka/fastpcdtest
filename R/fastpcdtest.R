#' fastpcdtest: A package for fast and robust calculation of Pesaran's CD test.
#'
#' The fastpcdtest package calculates Pesaran's CD test. It substitutes
#' for pcdtest in plm package in two situations:
#' 1) when the panel is huges (since it is faster than pcdtest in plm package) and
#' 2) when the panel is highly unbalanced (since is calculates CD test even when
#' some observations do not have overlap in their time dimension).
#'
#' @section fastpcdtest functions:
#' The package provides only one function, \code{pcdtest}, which overshadows
#' the test of the same name in \bold{plm} package. However, internally, there are
#' two Rcpp functions, one for balanced panels and other for unbalanced panels.
#' R wrapper function \code{pcdtest} automatically selects the proper function.
#'
#' @docType package
#' @name fastpcdtest
NULL




# fastpcdtest(mod) ... performs Pesaran CD test for cross-sectional dependence in panels;
# should work for both balanced and unbalanced panels
# (for balanced panels, switching to balanced cpp function could be faster)
#    mod ... plm model or data.frame with specified index
#    index ... character vector of length three containing names of columns
#              with cross-section unit (1st position), time unit (2nd position) and
#              residuals (3rd position).
# number of skipped correlations is in htest slot bad.obs

#' Calculate Pesaran's CD test.
#'
#' description
#'
#' details
#'
#' @param mod Either a plm model or data.frame with specified index.
#' @param index Character vector of length three containing names of columns with
#'              cross-section unit (1st position), time unit (2nd position) and
#'              residuals (3rd position).
#' @param min_common
#' @return Object of htest class with the results of the test.
#' @examples
#' #to be added
#' @export
pcdtest <- function(mod, index = NULL, min_common = 2L){
  # check the input
  if(!any((c("plm","data.frame") %in% class(mod)))){
    stop("The argument must be estimated plm model or data.frame with index specified!")
  }
  if("data.frame" %in% class(mod)){
    if(is.null(index)){stop("The argument index must be specified for a data.frame!")}
    if(length(index)!=3 | !is.character(index)){
      stop("The argument index must be character vector of length three!")
    }
    if(!all(index %in% names(mod))){
      stop("Index does not contain valid column names!")
    }
  }
  if(!(is.integer(min_common) & (min_common >= 1)))
    stop("Minimal number of common observations for any pair of individuals must be integer higher than one!")
  # require what is needed
  require(plm)
  require(tidyr)
  require(dplyr)
  require(magrittr)
  # create the matrix of residuals
  if("plm" %in% class(mod)){
    err  <- resid(mod)
    nerr <- length(err)
    merr <- cbind(resid(mod) %>% attr(.,"index") %>% as.data.frame(),
                  resid(mod)  %>% as.vector()) %>% magrittr::set_names(c("N","Tu","e"))
  }else{
    merr <- mod[,index] %>% magrittr::set_names(c("N","Tu","e"))
    nerr <- nrow(merr)
  }
  merr <- merr %>% tidyr::spread(Tu,e,fill=NA) %>% dplyr::select(-N) %>%
    as.matrix
  # calculate CD and z stratistics
  if(prod(dim(merr))>nerr){
    # for unbalanced panel
    result <- fastpcdtest_unbalanced_(merr, min_common)
    z <- result[1]
    invalid_cases <- result[2]
    invalid_cases_percent <- round(100 * invalid_cases / (nrow(merr) * (nrow(merr) - 1)), 2)
    if(invalid_cases > 0)
      warning(paste("Highly unbalanced panel:",
                    invalid_cases, "(i.e.", invalid_cases_percent,
                    "%) of correlation pairs could not been calculated."))
    if(invalid_cases_percent > 0.25)
      warning("Results of the test may be unreliable since too many correlation pair skipped.")
  }else{
    # for balanced panel
    z <- fastpcdtest_balanced_(merr)
    invalid_cases <- 0
    invalid_cases_percent <- 0
  }

  # calulate p-value
  pval <- 2 * pnorm(abs(z), lower.tail=FALSE)

  # set names
  names(z) <- "z"
  names(pval) <- "p-value"

  # output htest object
  structure(list(statistic = z,
                 # null.value = 0,
                 bad.obs = invalid_cases,
                 parameter = NULL,
                 method = "Pesaran CD test for cross-sectional dependence in panels",
                 data.name=deparse(substitute(mod)),
                 alternative = "cross-sectional dependence",
                 p.value = pval),
            class = "htest")
}



.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to fastpcdtest package.",
                        "Jak citovat.")
}
