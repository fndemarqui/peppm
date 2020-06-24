#' GTE Corporation telecommunication systems
#'
#' @name telecom
#' @docType data
#' @author Fabio N. Demarqui \email{fndemarqui@est.ufmg.br}
#' @keywords datasets
#' @description  Failure times (in days) of 125 installed telecommunication systems installed by GTE Corporation.
#' @format A data frame with 125 rows and 2 variables:
#' \itemize{
#'   \item time: vector of failure times (in days)
#'   \item status: vector of failure indicator
#' }
#' @examples 
#' library(peppm)
#' data(telecom)
#' fit1 <- with(telecom, peppm(time, status, cohesion=1, nburnin=0, nlag=1, npost=100))
#' fit2 <- with(telecom, peppm(time, status, cohesion=2, nburnin=0, nlag=1, npost=100))
#' fit3 <- with(telecom, peppm(time, status, cohesion=3, nburnin=0, nlag=1, npost=100))
#' fit4 <- with(telecom, peppm(time, status, cohesion=4, nburnin=0, nlag=1, npost=100))
#' # time grid associated with the first line of the matrix U:
#' 
#' @references  Piecewise Exponential Estimator for the Survival Function. J. S. Kim and F. Proschan.
#' IEEE TRANSACTIONS ON RELIABILITY, VOL. 40, NO. 2, 1991.
#' 
NULL