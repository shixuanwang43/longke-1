#' @title KE_bwselection
#'
#' @description Function used to perform leave-one-subject-out cross validation to select
#'              optimal time bandwidth (b_s) and trajectory bandwidth (b_w)
#'
#' @usage KE_bwselection(data,bw_time,bw_subj,T1,T2)
#'
#' @param data A long format data matrix containing columns ordered by
#'                 time, subject ID, response, predictor1, predictor2, ...
#'                 where the measurement time of the longitudinal data should be discretized
#' @param bw_time A numeric vector that contains the candidate time bandwidths
#' @param bw_subj A numeric vector that contains the candidate trajectory bandwidths
#' @param T1 A measurement time domain where the functional predictors are measured within
#' @param T2 A measurement time domain where the functional response is of interest to predict
#'
#' @return A list containing 3 elements
#' \item{BWSelecStep}{Total SSE for each bandwidth combination}
#' \item{optimalBW}{A vector containing the optimal time/trajectory bandwidth}
#' \item{RunningTime}{Running time of the bandwidth selection}
#'
#' @examples
#' \donttest{t_all = 1:50}
#' \donttest{T1=c(1,25);T2=c(26,50)}
#' \donttest{data = datagen(ntotal=350,ntest=50,t_all=t_all,t_split=25,seed=1)}
#' \donttest{data.sample = data$train}
#' \donttest{bwsele.toy = KE_bwselection(data=data.sample,}
#' \donttest{bw_time=c(1,2),bw_subj=c(0.1,0.5),T1=T1,T2=T2)}
#' \donttest{bwsele.toy$optimalBW}
#'
#' @references
#' \cite{Wang S, Kim S, Cho H, Chang W.
#'       Nonparametric predictive model for sparse and irregular longitudinal data. (2023+)}
#' @export
KE_bwselection = function(data,
                          bw_time,
                          bw_subj,
                          T1,T2) {

  bw_time=as.vector(bw_time);bw_subj=as.vector(bw_subj)
  stopifnot( "matrix" %in% class(data))
  stopifnot(class(bw_subj)=="numeric")
  stopifnot(class(bw_time)=="numeric")
  stopifnot(length(T1)==2)
  stopifnot(length(T2)<=2)
  #source("FPCA_trajectory.R")
  #library(bvls)
  ##    data : a long format data matrix containing columns:
  ##           time, subject ID, response, predictor1, predictor2, ... over T1+T2
  ##           the order has to be like this!!!!
  ## bw_time : a vector that contains the candidate time bandwidths
  ## bw_subj : a vector that contains the candidate trajectory bandwidths
  ##      T1 : discrete time interval where the functional predictor and response are both known
  ##      T2 : discrete time interval where the functional response is of interest to predict

  t_all = seq(min(data[,1]),max(data[,1]),1) # discrete time points i.e. 1,2,3,4 !!
  nfeature = ncol(data)-2 # number of features (response + predictor)


  fpca.tra = list() # list of multiple matrix nxm:
                    # n: number of subject; m: range of all measurement time
  fpca.fit = list() # fpca object used to predict testing data
  for (i in 1:nfeature) {

    dat.temp = data[,c(1,2,(i+2))]
    fit.temp = FPCA_trajectory(dat.temp,list(dataType='Sparse',
                                             error=FALSE, kernel='gauss', verbose=FALSE,
                                             nRegGrid=length(t_all)))
    fpca.tra[[i]] = fit.temp$target_fit # imputed trajectory in a matrix format
    fpca.fit[[i]] = fit.temp$fpca_target # fpca model object
  }

  tmin=min(data[,1]);tmax=max(data[,1]); # smallest and largest measurement time
  aa=T1[1]-tmin+1;bb=T1[2]-tmin+1;cc=T2[2]-tmin+1 # measurement time index of T1 and T2
  id_num=length(unique(data[,2])) # total number of subjects

  ## Compute training subject similarity over T1 and response similarity over T2
  distance.t1.joint = distance.t2.response = c()
  for (nf in 1:nfeature) {
    k=1
    distance.t1.temp=c()
    for(i in 1:id_num) {
      for(j in 1:id_num) {
        if (nf==1) { # get response similarity over T2
          distance.t2.response[k] = mean((fpca.tra[[1]][i,(bb+1):cc]-fpca.tra[[1]][j,(bb+1):cc])^2)
        }
        distance.t1.temp[k] = mean((fpca.tra[[nf]][i,(aa:bb)]-fpca.tra[[nf]][j,(aa:bb)])^2)
        k=k+1
      }
    }
    distance.t1.joint=cbind(distance.t1.joint,distance.t1.temp)
  }
  distance.comb = cbind(distance.t2.response,distance.t1.joint)# id_num*id_num by nfeature matrix
                                                               # that stores the similarity between
                                                               # ith subject and jth subject (i=j is possible)
  distance.comb = distance.comb[distance.comb[,1]!=0,] # remove similarity between ith and ith subject


  ## Fitting constraint linear regression to predict the response similarity over T2 by predictor similarity over T1
  nnls_distance = bvls::bvls(A=cbind(rep(1,nrow(distance.comb)),distance.comb[,2:(nfeature+1)]), b=distance.comb[,1],
                       bl = rep(0,nfeature+1), bu = rep(Inf, nfeature+1))



  ## Compute distance/similarity between ith subject and jth subject (i neq j) over T1
  distance.t1.joint.input = cbind(rep(1,nrow(distance.comb)),distance.comb[,2:(nfeature+1)])
  sub.d.hat = matrix(rep(NA,id_num*(id_num-1)),nrow = (id_num-1),ncol = id_num)
              # suppose we have 300 training subject, sub.d.hat is a 299x300 matrix
              # the ith column vector measures the similarity between ith subject and the rest 299 subjects
  for (i in 1:id_num) {
    sub.temp = distance.t1.joint.input[(1+(id_num-1)*(i-1)):((id_num-1)*i),]
    sub.d.hat[,i] = as.matrix(sub.temp)%*% as.matrix(nnls_distance$x)
  }





  ##  This function was designed to estimate the pointwise response trajectory for training data
  ##  Note that the variable in the global environment is also used in this function
  ##  ii: ii_th subject in the training data
  ##  hh: trajectory bandwidth
  ##  hhr: time bandwidth
  ##  tt: time point
  Single_Subject_Tra_Diff = function(ii,hh,hhr,tt) {


    train_data = data[data[,1]>T1[2],1:3] # training response observations over T2
    train_id = unique(train_data[,2])# training subject ID
    real_i = c()# sum(Y_ik*K_ik(s))
    real_i_ker = c() # K_w(w_hat)
    for (i in 1:(id_num-1)) {# for ith training subject

      valid.subj = train_data[train_data[,2]==train_id[i],]# observations for kth training subject
      mi = nrow(valid.subj) # number of observations for ith subject

      # Weighted avg with time dependency/kernel
      distance_r = valid.subj[,1]-tt # time difference
      kernel_r = dnorm(distance_r/hhr) # K_ik(s)
      real_i[i] = sum(kernel_r*valid.subj[,3])/mi# sum(Y_ik*K_ik(s))/n_i
      real_i_ker[i] = sum(kernel_r)/mi # sum(K_ik(s))/n_i

    }

    # Weighted avg with subject similarity/kernel
    distance = sub.d.hat[,ii]
    kernel = dnorm(sqrt(distance-nnls_distance$x[1])/hh)  # K_w(w_hat)

    # Point estimation at tt for subject ii
    ab = sum(real_i * kernel)/sum(real_i_ker * kernel) # mu_hat(s,0)
    return(ab)
  }


  ## Cross-validation to select optimal bandwidth
  ## The optimal combination of bandwidths is selected
  ##     by locating the smallest total_sse_new
  id.subj = unique(data[,2])
  total_sse_new = c() # Total sum of squared error for all training subjects
  sse = c() # sum of sqaured error for a single subject

  cnt = 1
  start.time = proc.time()
  for (h in bw_time) { # for each time bandwidth

    for (v in bw_subj) { # for each trajectory bandwidth


      for (j in id.subj) { # Error computation

        test_data = data[data[,2]==j,]

        # Point estimation for the trajectory
        pe.tra = c()
        for (i in test_data[,1]) {
          pe.tra[i-T1[2]] = Single_Subject_Tra_Diff(j,v,h,i)
        }
        # Total squared error for j^th subject between true measurement and trajectory point estimation
        sse[j] = sum((na.omit(pe.tra)-test_data[,3])^2)
      }
      # Total squared error for all subjects
      total_sse_new[cnt] = sum(sse)
      cat("\n","--------------------------------", "\n",
          paste0(cnt,"th bandwidth combination,"), "\n",
          "time bandwidth: ", h, "\n",
          "trajectory bandwidth: ", v, "\n",
          "Total SSE: ", sum(sse),"\n",
          "--------------------------------", "\n")
      cnt = cnt + 1
    }
  }
  end.time = proc.time() - start.time
  bw.comb = crossing(bw_time,bw_subj)
  bw.output = cbind(bw.comb,total_SSE=total_sse_new)
  bw.optim = bw.output[which.min(bw.output[,3]),]

  return(list(BWSelecStep = bw.output, # output total sse for every bandwidth combination
              optimalBW = bw.optim, # optimal bandwidth
              RunningTime = end.time)) # running time
}
