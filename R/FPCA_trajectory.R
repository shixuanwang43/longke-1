#' @title FPCA_trajectory
#'
#' @description Function used to perform functional principal component analysis (FPCA)
#'                for a single functional variable
#'
#' @usage FPCA_trajectory(data,...)
#'
#' @param data A long format data matrix containing 3 columns ordered by
#'               time, subject ID, variable where the measurement time of
#'               the longitudinal data should be discretized
#' @param ... Arguments to be passed to fdapace::FPCA
#'
#' @return A list containing two elements
#' \item{fpca_target}{A FPCA object}
#' \item{target_fit}{A num.t x num.sub matrix containing the imputated longitudinal trajectories
#'                    where num.t is the total number of the discrete measurement time and
#'                    num.sub is the total number of subjects}
#' @examples
#' t_all = 1:50
#' data = datagen(ntotal=350,ntest=50,t_all=t_all,t_split=25,seed=1)
#' data.sample = data$test[,c(1,2,3)]
#' # In this case, num.t=50 and num.sub=50 since we only used 50 subjects in the testing data
#' data.FPCA = FPCA_trajectory(data.sample,list(dataType='Sparse',
#'                 error=FALSE, kernel='gauss', verbose=FALSE, nRegGrid=length(t_all)))
#' data.FPCA$target_fit
#'
#' @references
#' \cite{Carroll, C., Gajardo, A., Chen, Y., Dai, X., Fan, J., Hadjipantelis, P. Z., ... & Wang, J. L. (2020).
#' fdapace: Functional data analysis and empirical dynamics. R package version 0.5, 4.}
#'
#' \cite{Yao, F., Müller, H. G., & Wang, J. L. (2005).
#' Functional data analysis for sparse longitudinal data. Journal of the American statistical association, 100(470), 577-590.}
#'
#' @export
FPCA_trajectory = function(data,...) {

  ## data: a long format Nx3 data matrix containing 3 columns:
  ##       time, subject ID, variable. N is the total sample size.
  ##       N is the number of observations instead of number of subjects.
  ##       time should be discretized to integer.

  ## data transformation
  uid = unique(data[,2]) # unique id
  n_subject = length(uid) # number of subjects
  target = time = list() # target and time are list that stores
                         # observations and time point of each subject
  for (i in 1:n_subject){
    target[i] =  list((data[data[,2] == uid[i],])[,3])
    time[i] = list((data[data[,2] == uid[i],])[,1])
  }

  #library(fdapace)
  #t_all = seq(min(data[,2]),max(data[,2]),1) # discrete time point i.e. 1,2,3,4 !!
  df = tibble(uid,target,time) # tibble format to fit FPCA from fdapace
  fpca_target = fdapace::FPCA(df$target, df$time,...)
  target_fit = fitted(fpca_target)
  return(list(target_fit=target_fit,
              fpca_target=fpca_target))

}
