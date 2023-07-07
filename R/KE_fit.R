#' @title KE_fit
#'
#' @description Function used to predict response trajectory by nonparametric kernel estimator
#'
#' @usage KE_fit(train,test,T1,T2,bw_time,bw_subj,alpha=0.05,seed=1,coefCI=FALSE)
#'
#' @param train A long format data matrix containing columns ordered by
#'              time, subject ID, response, predictor1, predictor2, ...
#'              where the measurement time of the longitudinal data should be discretized within T1.
#' @param test A long format data matrix containing columns ordered by
#'             time, subject ID, response, predictor1, predictor2, ...
#'             where the measurement time of the longitudinal data should be discretized within T2.
#' @param T1 A measurement time domain where the functional predictors are measured within
#' @param T2 A measurement time domain where the functional response is of interest to predict
#'
#'
#' @param bw_time (optimal) time bandwidth
#' @param bw_subj (optimal) trajectory/subject bandwidth
#' @param alpha confidence level for bootstrap CI of alpha_0, alpha_1, ...
#' @param seed A random seed fo producing replicable bootstrap CI of alpha_0, alpha_1, ...
#' @param coefCI Logical statement: TRUE to derive bootstrap CI of alpha0, alpha1, ... default is FALSE
#'
#' @return A list containing 6 elements
#' \item{testTraj}{A num.test x num.T2 matrix containing num.test subjects' trajectories where num.T2
#'                 is the total number of the discrete measurement time over T2}
#' \item{proxycoeff}{Coefficient estimation for the non-negative least square regression. From left to
#'                   right they are alpha_0, alpha_1, ...}
#' \item{fpca.fit}{A list containing FPCA fit for the functional predictors and the functional response}
#' \item{w.hat}{A list containing num.test elements where ith element contains the proxy distance/similarity
#'                between ith testing subject and other training subjects}
#' \item{bootCI.mean}{Bootstrap confidence interval of alpha_0, alpha_1, ...}
#' \item{input.list}{A list containing the input arguments}
#'
#'
#' @examples
#' \donttest{t_all = 1:50}
#' \donttest{T1=c(1,25);T2=c(26,50)}
#' \donttest{data = datagen(ntotal=350,ntest=50,t_all=t_all,t_split=25,seed=1)}
#' \donttest{train = data$train}
#' \donttest{test = data$test}
#' \donttest{ke.fit = KE_fit(train=train,test=test,T1=T1,T2=T2,bw_time=1,bw_subj = 0.2)}
#'
#' @references
#' \cite{Wang S, Kim S, Cho H, Chang W.
#'       Nonparametric predictive model for sparse and irregular longitudinal data. (2023+)}
#' @export

KE_fit = function(train,
                  test,
                  T1,T2,
                  bw_time,bw_subj,
                  alpha=0.05,
                  seed=1,
                  coefCI = FALSE) {


  stopifnot( "matrix" %in% class(test));stopifnot( "matrix" %in% class(train))
  stopifnot(length(T1)==2);stopifnot(length(T2)<=2)
  #source("FPCA_trajectory.R")

  #library(bvls)
  ##   train: a long format data matrix containing columns:
  ##          time, subject ID, response, predictor1, predictor2, ... over T1+T2
  ##          the order has to be like this!!!!
  ##    test: a long format data matrix containing columns:
  ##          time, subject ID, response predictor1, predictor2, ... over T1
  ##          Note the response over T1 is also considered as a predictor but here
  ##          I do not list it as predictor1 in case of confusion
  ##      T1: discrete time interval where the functional predictor and response are both known
  ##      T2: discrete time interval where the functional response is of interest to predict
  ## bw_time: (optimal) time bandwidth
  ## bw_traj: (optimal) trajectory/subject bandwidth
  ##          The above tuning parameters could be derived by KE_selection()
  ##   alpha: confidence level for bootstrap CI of alpha_0, alpha_1, ...
  ##    seed: produce replicable bootstrap CI of alpha_0, alpha_1, ...
  ##  coefCI: a logical statement: TRUE to derive bootstrap CI of alpha0, alpha1, ...
  ##                               default is FALSE

  t_all = seq(min(train[,1]),max(train[,1]),1) # discrete time point i.e. 1,2,3,4 !!
  nfeature = ncol(train)-2 # number of features (response + predictor)



  ## Fitting FPCA to impute the sparse and irregular longitudinal data
  fpca.tra = list() # list of multiple matrix nxm:
  #                      n: number of subject; m: span of measurement time
  fpca.fit = list() # fpca object used to predict testing data
  for (i in 1:nfeature) {

    dat.temp = train[,c(1,2,(i+2))]
    fit.temp = FPCA_trajectory(dat.temp,list(dataType='Sparse',
                                             error=FALSE, kernel='gauss', verbose=FALSE,
                                             nRegGrid=length(t_all)))
    fpca.tra[[i]] = fit.temp$target_fit # imputed trajectory in a matrix format
    fpca.fit[[i]] = fit.temp$fpca_target # fpca model object
  }

  tmin=min(train[,1]);tmax=max(train[,1]); # smallest and largest measurement time
  aa=T1[1]-tmin+1;bb=T1[2]-tmin+1;cc=T2[2]-tmin+1 # measurement time index of T1 and T2
  id_num=length(unique(train[,2])) # total number of training subjects


  ## Compute training subject similarity over T1 and response similarity over T2
  distance.t1.joint = distance.t2.response = c()
  for (nf in 1:nfeature) {
    k=1
    distance.t1.temp=c()
    for(i in 1:(id_num-1)) {
      for(j in (i+1):id_num) {
        if (nf==1) { # get response similarity over T2
          distance.t2.response[k] = mean((fpca.tra[[1]][i,(bb+1):cc]-fpca.tra[[1]][j,(bb+1):cc])^2)
        }
        distance.t1.temp[k] = mean((fpca.tra[[nf]][i,(aa:bb)]-fpca.tra[[nf]][j,(aa:bb)])^2)
        k=k+1
      }
    }
    distance.t1.joint=cbind(distance.t1.joint,distance.t1.temp)
  }
  distance.comb = cbind(distance.t2.response,distance.t1.joint) # id_num*(id_num-1)/2 by nfeature matrix
                                                                # that stores the similarity between
                                                                # ith subject and jth subject (i neq j)


  ## Fitting constraint linear regression to predict the response similarity over T2 by predictor similarity over T1
  nnls_distance = bvls::bvls(A=cbind(rep(1,nrow(distance.t1.joint)),distance.t1.joint), b=distance.t2.response,
                       bl = rep(0,nfeature+1), bu = rep(Inf, nfeature+1))



  ## Bootstrap confidence interval for coefficient in the constraint linear regression
  ## they are alpha_0, alpha_1, ...
  if (coefCI) {

    set.seed(seed)
    boot_sample_size = 1000 # bootstrap sample size
    warehouse = matrix(NA,boot_sample_size,nfeature+1) # storing each bootstrapped alpha_0, alpha_1, ...
    for (i in 1:boot_sample_size) {

      boot.index = sample(1:nrow(distance.t1.joint), size=nrow(distance.t1.joint), replace = TRUE, prob = NULL)
      boot.distance.t2.response = distance.t2.response[boot.index]
      boot.distance.t1.joint = distance.t1.joint[boot.index,]
      boot_nnls = bvls::bvls(A=cbind(rep(1,nrow(boot.distance.t1.joint)),boot.distance.t1.joint), b=boot.distance.t2.response,
                       bl = rep(0,nfeature+1), bu = rep(Inf, nfeature+1))
      warehouse[i,] = boot_nnls$x
    }
    bootCI = t(apply(warehouse,2,quantile,c(alpha/2,1-alpha/2)))
    rownames(bootCI)=paste0("alpha_",0:nfeature)
    bootCI.mean = cbind(nnls_distance$x,bootCI) # bootstrap CI and coefficient estimation for alpha_0, alpha_1, ...
    colnames(bootCI.mean)[1]="Est. Coefficient"
  } else {bootCI.mean="NULL"}





  ## here we use proposed kernel estimator to predict the response trajectory over T2
  sub.test = unique(test[,2]) # testing subject ID
  n.test = length(sub.test) # number of testing subject
  w.hat = list() # proxy
  new.pred.tra.joint = matrix(NA,nrow=n.test,ncol=T2[2]-T2[1]+1) # predicted response trajectory over T2 for all testing subjects
  for (i in 1:n.test) {


    test.temp = test[test[,2]==sub.test[i],]

    distance.t1.test = matrix(1,nrow=id_num,ncol=1+nfeature)
    ## Predict the testing trajectory over T1 using FPCA
    for (j in 1:nfeature) {
      # fpca.test.t1 is imputed predictors trajectory over T1 by FPCA for testing subjects
      fpca.test.t1 = predict(fpca.fit[[j]], newLy=list(test.temp[,j+2]), newLt=list(test.temp[,1]))$predCurves[,aa:bb]
      # distance.t1.test computes similarity between testing subjects over T1 by L2 norm
      distance.t1.test[,j+1] = apply(fpca.tra[[j]][,aa:bb], 1, function(x) mean((x-fpca.test.t1)^2))
    }
    w.hat[[i]] = distance.t1.test%*%as.matrix(nnls_distance$x) # proxy


    t_T2 = seq(T2[1],T2[2],1) # discrete time point over T2 i.e. 1,2,3,4 !!
    new.pred.tra=c() # predicted subject over T2 for a testing subject
    for (j in t_T2) {

      train.T2 = train[train[,1]>=T2[1],c(1,2,3)] # training response observations over T2
      train.id = unique(train[,2]) # training subject ID

      # Compute weighted obs
      real_i = rep(NA,id_num) # sum(Y_ik*K_ik(s))
      real_i_ker = rep(NA,id_num) # K_w(w_hat)
      for (k in 1:id_num) { # for kth training subject

        train.single.subj = train.T2[train.T2[,2]==train.id[k],] # observations for kth training subject
        mi = nrow(train.single.subj) # number of observations for kth training subject

        distance_r = train.single.subj[,1]-j # time difference
        kernel_r = dnorm(distance_r/bw_time) # K_ik(s)
        real_i[k] = sum(kernel_r*train.single.subj[,3])/mi # sum(Y_ik*K_ik(s))/n_i
        real_i_ker[k] = sum(kernel_r)/mi # sum(K_ik(s))/n_i

      }
      kernel = dnorm(sqrt(w.hat[[i]]-nnls_distance$x[1])/bw_subj)  # K_w(w_hat)
      new.pred.tra[j-T2[1]+1] = sum(real_i * kernel)/sum(real_i_ker * kernel) # mu_hat(s,0)
    }

    new.pred.tra.joint[i,] = new.pred.tra

  }

  # output the input element
  input.list = list(train=train,
                    test=test,
                    T1=T1,T2=T2,
                    bw_time=bw_time,bw_subj=bw_subj)

  KE = list(testTraj = new.pred.tra.joint, # testing trajectory by KE
            proxycoeff = nnls_distance$x, # coefficient estimation of constraint linear regression
            fpca.fit = fpca.fit, # fitted FPCA
            w.hat = w.hat, # proxy
            bootCI.mean = bootCI.mean, # bootstrap CI of alpha0, alpha1, ...
            input.list = input.list)
  class(KE) = "KE"
  KE

}
