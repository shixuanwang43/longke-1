#' @title KE_trajSCB
#'
#' @description Function used to derive simultaneous confidence band (SCB) for the predicted response trajectory
#'
#' @usage KE_trajSCB(KE.fit.object,nboot=500,alpha=0.05)
#'
#' @param KE.fit.object An object whose class is KE (you can get it by letting ke = KE.fit())
#' @param nboot Number of bootstrap sample size to construct SCB
#' @param alpha Confidence level for bootstrap SCB of predicted response trajectory
#'
#' @return A list containing num.test elements where the num.test represents the number of testing subjects.
#'         Within each element, there is a list containing 3 elements:
#' \item{se}{A vector containing standard errors at each discrete measurement time}
#' \item{traj.upper}{A vector containing upper bound of the testing subject at each measurement time}
#' \item{traj.lower}{A vector containing lower bound of the testing subject at each measurement time}
#'
#' @seealso \code{\link{KE_fit}}
#' @examples
#' \donttest{t_all = 1:50}
#' \donttest{T1=c(1,25);T2=c(26,50)}
#' \donttest{data = datagen(ntotal=350,ntest=50,t_all=t_all,t_split=25,seed=1)}
#' \donttest{train = data$train}
#' \donttest{test = data$test}
#' \donttest{ke.fit = KE_fit(train=train,test=test,T1=T1,T2=T2,bw_time=1,bw_subj = 0.2)}
#' \donttest{ketraj.toy = KE_trajSCB(KE.fit.object = ke.fit,
#'             nboot=10,alpha=0.05)}
#' @references
#' \cite{Wang S, Kim S, Cho H, Chang W.
#'       Nonparametric predictive model for sparse and irregular longitudinal data. (2023+)}
#'
#' \cite{Kim, S., Ryan Cho, H., & Kim, M. O. (2021).
#' Predictive generalized varying‐coefficient longitudinal model. Statistics in Medicine, 40(28), 6243-6259.}
#' @export



KE_trajSCB = function(KE.fit.object,nboot=500,alpha=0.05) {



  stopifnot(class(KE.fit.object)=="KE")
  ## KE.fit.object: This should be a KE object from function KE.fit()
  ##     nboot    : number of bootstrap sample replicates to construct
  ##                simultaneous confidence band (SCB) for the
  ##                predicted trajectory over T2
  ##     alpha    : confidence level for bootstrap SCB of estimated response trajectory

  ## Read input of KE.fit()
  train = KE.fit.object$input.list$train
  test = KE.fit.object$input.list$test
  T1 = KE.fit.object$input.list$T1
  T2 = KE.fit.object$input.list$T2
  w.hat = KE.fit.object$w.hat
  fpca.fit = KE.fit.object$fpca.fit
  proxycoeff = KE.fit.object$proxycoeff
  bw_time = KE.fit.object$input.list$bw_time
  bw_subj = KE.fit.object$input.list$bw_subj
  testTraj = KE.fit.object$testTraj


  n.train = length(unique(train[,2])) # number of training subject
  n.test = length(unique(test[,2])) # number of testing subject


  ## Construct SCB for ith subject
  T2.all = seq(T2[1],T2[2]) # discrete time point over T2 i.e. 1,2,3,4 !!
  out.CI = list() # output SCB for all testing subjects
  for (i in 1:n.test) { # for each testing subject


    boot.new.pred.tra = matrix(NA,nrow=nboot,ncol=T2[2]-T2[1]+1) # predicted response trajectory
                                                                 # for each bootstrap
    for (j in 1:nboot) { # for each time of bootstrap

      set.seed(i*1000+j)
      boot.id = sample(unique(train[,2]),replace = TRUE) # get bootstrap sample ID

      new.pred.tra = c() # predicted subject over T2 for a testing subject
      for(tt in T2.all) { # for every measurement time of interest

        train.T2 = train[train[,1]>=T2[1],c(1,2,3)]# training response observations over T2

        # Compute weighted obs
        real_i = rep(NA,n.train)# sum(Y_ik*K_ik(s))
        real_i_ker = rep(NA,n.train)# K_w(w_hat)
        for (k in 1:n.train) { # for kth training subject

          train.single.subj = train.T2[train.T2[,2]==boot.id[k],]# observations for kth training subject
          mi = nrow(train.single.subj)# number of observations for kth training subject

          distance_r = train.single.subj[,1]-tt# time difference
          kernel_r = dnorm(distance_r/bw_time)# K_ik(s)
          real_i[k] = sum(kernel_r*train.single.subj[,3])/mi # sum(Y_ik*K_ik(s))/n_i
          real_i_ker[k] = sum(kernel_r)/mi # sum(K_ik(s))/n_i
        }
        kernel = dnorm(sqrt(w.hat[[i]][boot.id]-proxycoeff[1])/bw_subj)   # K_w(w_hat)
        new.pred.tra[tt-T2[1]+1] = sum(real_i * kernel)/sum(real_i_ker * kernel)# mu_hat(s,0)

      }
      boot.new.pred.tra[j,] = new.pred.tra
    }

    Y_pred = testTraj[i,] # ith predicted response trajectory over T2
    boot_data = boot.new.pred.tra
    s_boot = apply(boot_data,2,sd) # standard deviation at each measurement time
    L = c() #
    for (ii in 1:nrow(boot_data)) {
      L[ii] = max(abs(boot_data[ii,]-Y_pred)/s_boot) # supremum of standardized deviation in absolute value
    }

    out.CI[[i]] = list(se = s_boot*quantile(L,1-alpha),
                       traj.upper = Y_pred+s_boot*quantile(L,1-alpha),
                       traj.lower = Y_pred-s_boot*quantile(L,1-alpha))

    cat(paste0(i,"th subject completed \n"))

  }

  return(out.CI)

}
