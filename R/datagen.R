#' @title Simulate longitudinal data
#'
#' @description Function used to simulate sample sparse and irregular longitudinal data
#'
#' @usage datagen(ntotal,ntest,t_all,t_split,seed)
#'
#' @param ntotal Number of total longitudinal subjects
#' @param ntest Number of total longitudinal subjects in the testing set
#' @param t_all Vector of discrete measurement time (i.e 1,2,3,4,...)
#' @param t_split A measurement time where the longitudinal response is of interest to predict
#'                  after this t_split
#' @param seed Seed to derive replicable data
#'
#' @return A list containing two elements
#' \describe{
#' \item{train}{A long format data matrix containing one functional response (yy) and
#'              two functional predictors (xx,zz) with (ntotal-ntest) subjects}
#' \item{test}{A long format data matrix containing one functional response (yy) and
#'              two functional predictors (xx,zz) with (ntest) subjects}
#' }
#'
#' @examples
#' data = datagen(ntotal=350,ntest=50,t_all=1:50,t_split=25,seed=1)
#' data$test
#' data$train
#' @export

datagen = function(ntotal,
                   ntest,
                   t_all,
                   t_split,
                   seed) {


  stopifnot(ntotal>ntest)
  stopifnot((t_split>min(t_all))&(t_split<max(t_all)))
  #library(tidyverse)
  #library(mvtnorm)

  ss=seed # specification of random seed to derive replicable simulation result
  set.seed(ss)
  ### The number of subjects in the training and testing set
  n=ntotal
  ### The number of subjects in the testing set
  n_test=ntest
  ### Total time period of interest, that is the union of T1 and T2
  t_all = t_all
  t_len = length(t_all)

  ### Generate the training and testing sets
  ######## X,Z ########
  ni = length(t_all)
  X = matrix(rep(NA,length(t_all)*n),nrow = length(t_all))
  Z = matrix(rep(NA,length(t_all)*n),nrow = length(t_all))

  ### Hyperparameter of Covariance matrix in Gaussian process (length scale, variance)
  ### Squared exponential covariance function is used.
  hyperp = c(50,1) # lengthscale = 50 and variance = 1
  for (i in 1:n) {

    l = abs(outer(t_all,t_all,"-"))
    Sigma_SE = exp(-l^2/(2*hyperp[1]^2))
    sss = Sigma_SE*hyperp[2] # GP covariance matrix
    a0=runif(1,-3.5,3.5) # random subject effect

    mu.x = t_all/5 # marginal mean function of x
    X[,i] = mu.x+a0+mvtnorm::rmvnorm(1,rep(0,t_len),sss)

    mu.z =  5*cos(2*pi*t_all/50) # marginal mean function of z
    Z[,i] = mu.z+a0+mvtnorm::rmvnorm(1,rep(0,t_len),sss)

  }
  ######## Y ########
  Y = matrix(rep(NA,length(t_all)*n),nrow = length(t_all))
  hyperp = c(50,1)
  ni = length(t_all)
  for (i in 1:n) {

    l = abs(outer(t_all,t_all,"-"))
    #lx = abs(outer(X[,i],X[,i],"-"))
    #lz = abs(outer(Z[,i],Z[,i],"-"))
    Sigma_SE = exp(-l^2/(2*hyperp[1]^2))
    sss = Sigma_SE*hyperp[2]

    mu.v = rep(0,length(t_all))# marginal mean function of y

    # Linear Model
    Y[,i] = X[,i]/5+Z[,i]+mvtnorm::rmvnorm(1,rep(0,t_len),sss)

    # Quadratic Model
    #Y[,i] = X[,i]/5+Z[,i]+X[,i]^2/50 +mvtnorm::rmvnorm(1,rep(0,t_len),sss)

    # Interaction Model
    #Y[,i] = X[,i]/5+Z[,i]+X[,i]^2/50+X[,i]*Z[,i]/50 +mvtnorm::rmvnorm(1,rep(0,t_len),sss)

  }

  ### Discrete measurement times
  index_Y = seq(1,length(t_all),1)
  index_X = seq(1,length(t_all),1)
  index_Z = seq(1,length(t_all),1)
  ## Derive Sparse and irregular longitudinal index
  set.seed(ss)
  mat_Z=mat_Y=mat_X=ni_index_Z=ni_index_Y = ni_index_X = ni_index_complement_Z=ni_index_complement_Y = ni_index_complement_X = list(c())
  for (i in 1:n) {

    tenth = trunc(t_len/10)*10
    kk.t = tenth/5
    kk=(5*(0:(kk.t-1))+purrr::rdunif(kk.t,1,5)) # target index of measurement time (that would be picked)
    ni_index_Y[[i]] = kk
    ni_index_X[[i]] = kk
    ni_index_Z[[i]] = kk

    mat_Y[[i]] = match(ni_index_Y[[i]],t_all)
    mat_X[[i]] = match(ni_index_X[[i]],t_all)
    mat_Z[[i]] = match(ni_index_Z[[i]],t_all)

    ni_index_complement_Y[[i]] = index_Y[-mat_Y[[i]]] # complement target index of measurement time
    ni_index_complement_X[[i]] = index_X[-mat_X[[i]]] # (that would be removed)
    ni_index_complement_Z[[i]] = index_Z[-mat_Z[[i]]]
  }
  ### Derive sparse and irregular longitudinal features (Y_i(t),X_i(t),Z_i(t)) for ith subject
  Y_sparse = Y
  X_sparse = X
  Z_sparse = Z
  for (i in 1:n) {
    Y_sparse[,i][ni_index_complement_Y[[i]]] = NA
    X_sparse[,i][ni_index_complement_X[[i]]] = NA
    Z_sparse[,i][ni_index_complement_Z[[i]]] = NA
  }
  set.seed(ss)
  test_index = sample(x=1:ncol(Y),size=n_test) # randomly select t_len testing subjects
  test_Y = Y[,test_index];test_X = X[,test_index];test_Z = Z[,test_index]
  test_Y_sparse =  Y_sparse[,test_index]; test_X_sparse = X_sparse[,test_index];test_Z_sparse = Z_sparse[,test_index]
  Y_sparse = Y_sparse[,-test_index]; X_sparse = X_sparse[,-test_index];Z_sparse = Z_sparse[,-test_index]
  Y = Y[,-test_index]; X = X[,-test_index];Z = Z[,-test_index]
  ###  Transform to long format data matrix for training and testing data set
  n = ncol(Y) # number of training subjects
  y_flm = Y_sparse %>%
    t() %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "fsd", values_to = "vvv") %>%
    mutate("age" = rep(seq(1:length(t_all)),n), "ID" = sort(rep(seq(1,n),length(t_all)))) %>%
    drop_na() %>%
    select(-"fsd") %>% rename("yy" = "vvv")
  x_flm = X_sparse %>%
    t() %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "fsd", values_to = "vvv") %>%
    mutate("age" = rep(seq(1:length(t_all)),n), "ID" = sort(rep(seq(1,n),length(t_all)))) %>%
    drop_na() %>%
    select(-"fsd") %>% rename("xx" = "vvv")
  z_flm = Z_sparse %>%
    t() %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "fsd", values_to = "vvv") %>%
    mutate("age" = rep(seq(1:length(t_all)),n), "ID" = sort(rep(seq(1,n),length(t_all)))) %>%
    drop_na() %>%
    select(-"fsd") %>% rename("zz" = "vvv")

  y_flm_test = test_Y_sparse %>%
    t() %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "fsd", values_to = "vvv") %>%
    mutate("age" = rep(seq(1:length(t_all)),n_test), "ID" = sort(rep(seq(1,n_test),length(t_all)))) %>%
    drop_na() %>%
    select(-"fsd") %>% rename("yy" = "vvv")
  x_flm_test = test_X_sparse %>%
    t() %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "fsd", values_to = "vvv") %>%
    mutate("age" = rep(seq(1:length(t_all)),n_test), "ID" = sort(rep(seq(1,n_test),length(t_all)))) %>%
    drop_na() %>%
    select(-"fsd") %>% rename("xx" = "vvv")
  z_flm_test = test_Z_sparse %>%
    t() %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "fsd", values_to = "vvv") %>%
    mutate("age" = rep(seq(1:length(t_all)),n_test), "ID" = sort(rep(seq(1,n_test),length(t_all)))) %>%
    drop_na() %>%
    select(-"fsd") %>% rename("zz" = "vvv")

  ## Split data into training and testing
  joint.temp1 = full_join(z_flm,x_flm,by = c("age", "ID"))
  joint.temp2 = full_join(z_flm_test,x_flm_test,by = c("age", "ID"))
  train = full_join(y_flm,joint.temp1,by = c("age", "ID")) %>% select("age","ID","yy","xx","zz") %>%as.matrix()
  test.temp =  full_join(y_flm_test,joint.temp2,by = c("age", "ID")) %>% select("age","ID","yy","xx","zz")
  test = as.matrix(test.temp[test.temp$age<=t_split,])
  return(list(train=train,
              test=test))
}

