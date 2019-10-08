#########################################################################################
#########################################################################################
##
## Simulation study for evaluating the bivariate probit model with sample selection
## and random intercept
##
## Consideration of different missing-data mechanisms: 
## -> MAR, MNAR selection model, MNAR non-selection model
##
## @author Angelina Hammon
##
## 05.09.2019
##
#########################################################################################
#########################################################################################


# install package for model calculation: 
install.packages("bivprob.quad_1.0.tar.gz", repos = NULL, type="source")

library(mvtnorm)
library(mice)
library(lme4)
library(PanelCount)
library(doSNOW)
library(Rcpp)
library(bivprob.quad)



### true values ### 

model_bd <- vector("numeric",length=4)

rhos <- c(0,0.3,0.6,0.9)
taus <- c(0,0.1,0.3,0.5)

# Non-Heckman MNAR & MAR
rho <- rhos[1]
tau <- taus[1]

# Heckman MNAR:
rho <- rhos[4]  # change depending on scenario
tau <- taus[4] 

for(i in 1:100000) {
  n <- 2500
  m <- 50
  nj <- n/m
  
  x1 <- rnorm(n,0,0.3)
  x2 <- rnorm(n,0,0.8)
  x3 <- rnorm(n,0,4)
  
  vc <- diag(2)
  vc[2,1] <- vc[1,2] <- rho
  eps <- rmvnorm(n,rep(0,2),vc)
  eps_o <- eps[,1]
  eps_s <- eps[,2]
  
  vc_re <- matrix(1,2,2)
  vc_re[1,1] <- 0.5
  vc_re[2,2] <- 0.9
  vc_re[2,1] <- vc_re[1,2] <- sqrt(vc_re[1,1])*sqrt(vc_re[2,2])*tau
  
  alpha <- rmvnorm(m,rep(0,2),vc_re) 
  alpha_s <- rep(alpha[,1], each=nj)  
  alpha_o <- rep(alpha[,2], each=nj)  
  
  y <- 0.25+1*x1+0.5*x2+eps_o+alpha_o
  y <- ifelse(y<0,0,1)
  
  r <- 0.5+1.5*x1-0.25*x2+0.1*x3+eps_s+alpha_s
  r <- ifelse(r<0,0,1)
  
  group <- rep(1:m, each=nj) 
  
  data <- data.frame(y,r,x1,x2,x3,group)
  
  model <- glmer(y~x1+x2+(1|group),data=data,family=binomial("probit"),nAGQ=5)
  
  model_bd <- model_bd + c(model@beta,model@theta^2)
  
}

model_bd <- model_bd/100000



### Functions ###

coverage <- function(value, CI.low, CI.upper){
  ifelse(CI.low <= value && CI.upper >= value,1,0)
}

rel.bias <- function(value,est) {
  rel_bias <- (est-value)/value
  return(rel_bias)
}

RMSE <- function(value,est){
  mse <- (est-value)^2
  rmse <- sqrt(mse)
  return(mse)
}





### Data Generation ###

cl <- makeSOCKcluster(10)
registerDoSNOW(cl)
clusterSetupRNG(cl)
clusterExport(cl,"rel.bias")
clusterExport(cl,"coverage")
clusterExport(cl,"RMSE")
clusterExport(cl,"model_bd")
clusterExport(cl,"rhos")
clusterExport(cl,"taus")
clusterEvalQ(cl, {
  library(mice)
  library(mvtnorm)
  library(doParallel)
  library(lme4)
  library(Hmisc)
  library(Rcpp)
  library(RcppArmadillo)
  library(corpcor)
  library(compiler)
  library(nloptr)
  library(micemd)
  library(miceMNAR)
  library(bivprob.quad)
  source("function_MI_analysis.R")
  source("binom.2l.r")
  source("mice_heckman_binary_re1.R")
  source("mice_heckman_binary_re2.R")
})


start_time <- Sys.time()
sim <- 
  parLapply(cl = cl, 1:1000, fun = function(no, mech="non-heckman",rho=rhos[1],tau=taus[1],M=10){
    
    n <- 2500
    
    m <- 50
    nj <- n/m
    
    x1 <- rnorm(n,0,0.3)
    x2 <- rnorm(n,0,0.8)
    x3 <- rnorm(n,0,4)
    
    
    ## Missing-data mechanisms ##
    
    if(mech=="heckman") {
      
      if(rho == 0) {stop("rho is equal to 0")}
      
      vc <- diag(2)
      vc[2,1] <- vc[1,2] <- rho
      eps <- rmvnorm(n,rep(0,2),vc)
      eps_o <- eps[,1]
      eps_s <- eps[,2]
      
      vc_re <- matrix(1,2,2)
      vc_re[1,1] <- 0.5
      vc_re[2,2] <- 0.9
      vc_re[2,1] <- vc_re[1,2] <- sqrt(vc_re[1,1])*sqrt(vc_re[2,2])*tau
      
      alpha <- rmvnorm(m,rep(0,2),vc_re) 
      alpha_s <- rep(alpha[,1], each=nj)  
      alpha_o <- rep(alpha[,2], each=nj)  
      
      y <- 0.25+1*x1+0.5*x2+eps_o+alpha_o
      y <- ifelse(y<0,0,1)
      
      r <- 0.5+1.5*x1-0.25*x2+0.1*x3+eps_s+alpha_s
      r <- ifelse(r<0,0,1)
      
      
    } else if(mech=="non-heckman"){
      
      if(rho != 0) {stop("rho is unequal to 0")}
      
      vc <- diag(2)
      vc[2,1] <- vc[1,2] <- rho
      eps <- rmvnorm(n,rep(0,2),vc)
      eps_o <- eps[,1]
      eps_s <- eps[,2]
      
      vc_re <- matrix(1,2,2)
      vc_re[1,1] <- 0.5
      vc_re[2,2] <- 0.9
      vc_re[2,1] <- vc_re[1,2] <- sqrt(vc_re[1,1])*sqrt(vc_re[2,2])*tau 
      
      alpha <- rmvnorm(m,rep(0,2),vc_re) 
      alpha_s <- rep(alpha[,1], each=nj)  
      alpha_o <- rep(alpha[,2], each=nj)  
      
      y <- 0.25+1*x1+0.5*x2+eps_o+alpha_o
      
      prob <- pnorm(0.8+1.75*y+1.5*x1-2.5*x2+alpha_s)
      r <- rbinom(n,1,p=prob)
      
      y <- ifelse(y<0,0,1)
      
      
    } else{ # MAR
      if(rho != 0) {stop("rho is unequal to 0")}
      
      vc <- diag(2)
      vc[2,1] <- vc[1,2] <- rho
      eps <- rmvnorm(n,rep(0,2),vc)
      eps_o <- eps[,1]
      eps_s <- eps[,2]
      
      vc_re <- matrix(1,2,2)
      vc_re[1,1] <- 0.5
      vc_re[2,2] <- 0.9
      vc_re[2,1] <- vc_re[1,2] <- sqrt(vc_re[1,1])*sqrt(vc_re[2,2])*tau
      
      alpha <- rmvnorm(m,rep(0,2),vc_re) 
      alpha_s <- rep(alpha[,1], each=nj)  
      alpha_o <- rep(alpha[,2], each=nj)  
      
      y <- 0.25+1*x1+0.5*x2+eps_o+alpha_o
      y <- ifelse(y<0,0,1)
      
      r <- 0.5+1.5*x1-0.25*x2+0.1*x3+eps_s+alpha_s
      r <- ifelse(r<0,0,1)
      
      
    }
    
    
    group <- rep(1:m, each=nj) 
    
    
    ## BD ##
    data_comp <- data.frame(y,r,x1,x2,x3,group)
    
    bd <- glmer(y~x1+x2+(1|group),data=data_comp,family=binomial("probit"),nAGQ=5)
    
    CI_bd <- confint(bd,method="Wald")[2:(length(bd@beta)+1),]

    cov_bd <- NULL
    for(j in 1:(length(model_bd)-1)){
      cov_bd <- c(cov_bd,coverage(model_bd[j],CI_bd[j,1],CI_bd[j,2]))
    }
    
    rmse.est_bd <- RMSE(model_bd,c(bd@beta,bd@theta^2))
    
    var_bd <- c((coef(summary(bd))[,2])^2,0)
    
    
    ## Generation of missing values ##
    
    # based on upper specified mechanism in y:
    y[r==0] <- NA
    
    data <- data.frame(y,r,x1,x2,x3,group)
    
    
    ### Imputation ###
    
    ini <- mice(data,m=1,maxit=0) #,visitSequence="arabic")
    
    ## MNAR Heckman AGHQ ##
    pred_MNAR <- ini$pred
    pred_MNAR["y","r"] <- 0
    pred_MNAR["y","group"] <- -2
    
    imp_MNAR_aghq <- mice(data,m=M,maxit=1,method=c("2l.heckman1step_re_aghq","","","","",""),pred=pred_MNAR,print=F,QP=rep(10,2),draw=T,excl="x3")
    
    glm_MNAR_aghq <- with(imp_MNAR_aghq, glmer(y~x1+x2 + (1|group),family=binomial("probit"),nAGQ=5)) 
    est <- t(sapply(1:M,function(x) glm_MNAR_aghq$analyses[[x]]@beta))
    var <- t(sapply(1:M,function(x) diag(vcov(glm_MNAR_aghq$analyses[[x]]))))
    re <- t(sapply(1:M,function(x)  glm_MNAR_aghq$analyses[[x]]@theta^2))
    pool_glm_MNAR_aghq <- MI.analysis(est,var,m=M)
    pool_glm_MNAR_aghq <- rbind(pool_glm_MNAR_aghq,c(mean(re),rep(0,3)))
    
    cov_MNAR_aghq <- NULL
    for(j in 1:ncol(est)){
      cov_MNAR_aghq <- c(cov_MNAR_aghq,coverage(model_bd[j],pool_glm_MNAR_aghq[j,2],pool_glm_MNAR_aghq[j,3]))
    }
    
    rmse.est_MNAR_aghq <- RMSE(model_bd,pool_glm_MNAR_aghq[,1])
    
    var_MNAR_aghq <- pool_glm_MNAR_aghq[,4]
    
    
    ## MNAR Heckman GHQ ##
    imp_MNAR_ghq <- mice(data,m=M,maxit=1,method=c("2l.heckman1step_re_ghq","","","","",""),pred=pred_MNAR,print=F,QP=rep(10,2),draw=T,excl="x3")
    
    glm_MNAR_ghq <- with(imp_MNAR_ghq, glmer(y~x1+x2 + (1|group),family=binomial("probit"),nAGQ=5)) 
    est <- t(sapply(1:M,function(x) glm_MNAR_ghq$analyses[[x]]@beta))
    var <- t(sapply(1:M,function(x) diag(vcov(glm_MNAR_ghq$analyses[[x]]))))
    re <- t(sapply(1:M,function(x)  glm_MNAR_ghq$analyses[[x]]@theta^2))
    pool_glm_MNAR_ghq <- MI.analysis(est,var,m=M)
    pool_glm_MNAR_ghq <- rbind(pool_glm_MNAR_ghq,c(mean(re),rep(0,3)))
    
    cov_MNAR_ghq <- NULL
    for(j in 1:ncol(est)){
      cov_MNAR_ghq <- c(cov_MNAR_ghq,coverage(model_bd[j],pool_glm_MNAR_ghq[j,2],pool_glm_MNAR_ghq[j,3]))
    }
    
    rmse.est_MNAR_ghq <- RMSE(model_bd,pool_glm_MNAR_ghq[,1])
    
    var_MNAR_ghq <- pool_glm_MNAR_ghq[,4]
    
    
    ## MAR ##     pred_MAR <- ini$pred
    pred_MAR <- ini$pred
    pred_MAR[,"r"] <- 0
    pred_MAR["y","group"] <- -2
    
    imp_MAR <- mice(data,m=M,maxit=1,method=c("2l.binary","","","","",""),pred=pred_MAR,print=F)
    
    glm_MAR <- with(imp_MAR, glmer(y~x1+x2 + (1|group),family=binomial("probit"),nAGQ=5)) 
    est <- t(sapply(1:M,function(x) glm_MAR$analyses[[x]]@beta))
    var <- t(sapply(1:M,function(x) diag(vcov(glm_MAR$analyses[[x]]))))
    re <- t(sapply(1:M,function(x)  glm_MAR$analyses[[x]]@theta^2))
    pool_glm_MAR <- MI.analysis(est,var,m=M)
    pool_glm_MAR <- rbind(pool_glm_MAR,c(mean(re),rep(0,3)))
    
    cov_MAR <- NULL
    for(j in 1:ncol(est)){
      cov_MAR <- c(cov_MAR,coverage(model_bd[j],pool_glm_MAR[j,2],pool_glm_MAR[j,3]))
    }
    
    rmse.est_MAR <- RMSE(model_bd,pool_glm_MAR[,1])
    
    var_MAR <- pool_glm_MAR[,4]
    
    
    
    ## two stage method of package micemd ##
    
    data2 <- data
    data2$y <- as.factor(data2$y)
    
    imp_MAR_twostage <- mice(data2,m=M,maxit=1,method=c("2l.2stage.bin","","","","",""),pred=pred_MAR,print=F)
    
    glm_MAR_twostage <- with(imp_MAR_twostage, glmer(y~x1+x2 + (1|group),family=binomial("probit"),nAGQ=5)) 
    est <- t(sapply(1:M,function(x) glm_MAR_twostage$analyses[[x]]@beta))
    var <- t(sapply(1:M,function(x) diag(vcov(glm_MAR_twostage$analyses[[x]]))))
    re <- t(sapply(1:M,function(x)  glm_MAR_twostage$analyses[[x]]@theta^2))
    pool_glm_MAR_twostage <- MI.analysis(est,var,m=M)
    pool_glm_MAR_twostage  <- rbind(pool_glm_MAR_twostage,c(mean(re),rep(0,3)))
    
    cov_MAR_twostage <- NULL
    for(j in 1:ncol(est)){
      cov_MAR_twostage <- c(cov_MAR_twostage,coverage(model_bd[j],pool_glm_MAR_twostage[j,2],pool_glm_MAR_twostage[j,3]))
    }
    
    rmse.est_MAR_twostage <- RMSE(model_bd,pool_glm_MAR_twostage[,1])
    
    var_MAR_twostage <- pool_glm_MAR_twostage[,4]
    
    
    ## MNAR model with ML structure (Galimard) ##
    
    JointModelEq <- generate_JointModelEq(data=data2,varMNAR = "y")
    
    JointModelEq[,"y_var_sel"] <- c(0,0,1,1,1,0)
    JointModelEq[,"y_var_out"] <- c(0,0,1,1,0,0)
    
    arg <- MNARargument(data=data2,varMNAR="y",JointModelEq=JointModelEq)
    
    arg$predictorMatrix[,c("group","r")] <- 0 
    arg$predictorMatrix[c("r","x1","x2","x3","group"),] <- 0
    
    
    imp_MNAR_gal <- mice(data = arg$data_mod, method = arg$method,predictorMatrix = arg$predictorMatrix,
                         JointModelEq=arg$JointModelEq,control=arg$control,m=M)
    
    glm_MNAR_gal <- with(imp_MNAR_gal, glmer(y~x1+x2 + (1|group),family=binomial("probit"),nAGQ=5)) 
    est <- t(sapply(1:M,function(x) glm_MNAR_gal$analyses[[x]]@beta))
    var <- t(sapply(1:M,function(x) diag(vcov(glm_MNAR_gal$analyses[[x]]))))
    re <- t(sapply(1:M,function(x)  glm_MNAR_gal$analyses[[x]]@theta^2))
    pool_glm_MNAR_gal <- MI.analysis(est,var,m=M)
    pool_glm_MNAR_gal  <- rbind(pool_glm_MNAR_gal,c(mean(re),rep(0,3)))
    
    cov_MNAR_gal <- NULL
    for(j in 1:ncol(est)){
      cov_MNAR_gal <- c(cov_MNAR_gal,coverage(model_bd[j],pool_glm_MNAR_gal[j,2],pool_glm_MNAR_gal[j,3]))
    }
    
    rmse.est_MNAR_gal <- RMSE(model_bd,pool_glm_MNAR_gal[,1])
    
    var_MNAR_gal <- pool_glm_MNAR_gal[,4]
    
    
    ## CC ##
    cc <- glmer(y~x1+x2+(1|group),data=data,family=binomial("probit"),nAGQ=5)
    
    CI_cc <- confint(cc,method="Wald")[2:(length(cc@beta)+1),]

    cov_CC <- NULL
    for(j in 1:(length(model_bd)-1)){
      cov_CC <- c(cov_CC,coverage(model_bd[j],CI_cc[j,1],CI_cc[j,2]))
    }
    
    rmse.est_cc <- RMSE(model_bd,c(cc@beta,cc@theta^2))
    
    var_cc <- c((coef(summary(cc))[,2])^2,0)
    
    
    results <- cbind(pool_glm_MAR[,1], pool_glm_MAR_twostage[,1],
                     pool_glm_MNAR_aghq[,1],
                     pool_glm_MNAR_ghq[,1],
                     pool_glm_MNAR_gal[,1],
                     c(cc@beta,cc@theta^2), c(bd@beta,bd@theta^2),
                     rel.bias(model_bd,pool_glm_MAR[,1]), rel.bias(model_bd,pool_glm_MAR_twostage[,1]),
                     rel.bias(model_bd,pool_glm_MNAR_aghq[,1]), 
                     rel.bias(model_bd,pool_glm_MNAR_ghq[,1]),  
                     rel.bias(model_bd,pool_glm_MNAR_gal[,1]),  
                     rel.bias(model_bd, c(cc@beta,cc@theta^2)), rel.bias(model_bd,c(bd@beta,bd@theta^2)),
                     cov_MAR,cov_MAR_twostage,  
                     cov_MNAR_aghq,cov_MNAR_ghq,cov_MNAR_gal,
                     cov_CC,cov_bd,
                     rmse.est_MAR,rmse.est_MAR_twostage,rmse.est_MNAR_aghq,rmse.est_MNAR_ghq,rmse.est_MNAR_gal,
                     rmse.est_cc,rmse.est_bd,
                     var_MAR,var_MAR_twostage,var_MNAR_aghq,var_MNAR_ghq, var_MNAR_gal,
                     var_cc,var_bd)
    
    results[4,c(15:21)] <- 0
    
    rm(data,data_comp,imp_MAR,imp_MNAR_aghq,imp_MNAR_ghq,data2,imp_MAR_twostage,imp_MNAR_gal) 
    
    return(results)
    
  })
end_time <- Sys.time() 
end_time - start_time 

stopCluster(cl)


## Diagnostics ##

est_MAR <- rowMeans(sapply(sim,function(x) x[,1]))
est_MAR_twostage <- rowMeans(sapply(sim,function(x) x[,2]))
est_MNAR_aghq <- rowMeans(sapply(sim,function(x) x[,3]))
est_MNAR_ghq <- rowMeans(sapply(sim,function(x) x[,4]))
est_MNAR_gal <- rowMeans(sapply(sim,function(x) x[,5]))
est_CC <- rowMeans(sapply(sim,function(x) x[,6]))   
est_bd <- rowMeans(sapply(sim,function(x) x[,7]))   

bias_MAR <- rowMeans(sapply(sim,function(x) x[,8]))
bias_MAR_twostage <- rowMeans(sapply(sim,function(x) x[,9]))
bias_MNAR_aghq <- rowMeans(sapply(sim,function(x) x[,10]))
bias_MNAR_ghq <- rowMeans(sapply(sim,function(x) x[,11]))
bias_MNAR_gal <- rowMeans(sapply(sim,function(x) x[,12]))
bias_CC <- rowMeans(sapply(sim,function(x) x[,13]))
bias_bd <- rowMeans(sapply(sim,function(x) x[,14]))

cov_MAR <- rowMeans(sapply(sim,function(x) x[,15]))
cov_MAR_twostage <- rowMeans(sapply(sim,function(x) x[,16]))
cov_MNAR_aghq <- rowMeans(sapply(sim,function(x) x[,17]))
cov_MNAR_ghq <- rowMeans(sapply(sim,function(x) x[,18]))
cov_MNAR_gal <- rowMeans(sapply(sim,function(x) x[,19]))
cov_CC <- rowMeans(sapply(sim,function(x) x[,20]))
cov_bd <- rowMeans(sapply(sim,function(x) x[,21]))

rmse.est_MAR <- rowMeans(sapply(sim,function(x) x[,22]))
rmse.est_MAR_twostage <- rowMeans(sapply(sim,function(x) x[,23]))
rmse.est_MNAR_aghq <- rowMeans(sapply(sim,function(x) x[,24]))
rmse.est_MNAR_ghq <- rowMeans(sapply(sim,function(x) x[,25]))
rmse.est_MNAR_gal <- rowMeans(sapply(sim,function(x) x[,26]))
rmse.est_CC <- rowMeans(sapply(sim,function(x) x[,27]))
rmse.est_bd <- rowMeans(sapply(sim,function(x) x[,28]))

it <- length(sim)
emp.se_MAR <- sqrt(rowSums(sapply(sim,function(x) (x[,1]- est_MAR)^2))*(1/(it-1)))
emp.se_MAR_twostage <- sqrt(rowSums(sapply(sim,function(x) (x[,2]- est_MAR_twostage)^2))*(1/(it-1)))
emp.se_MNAR_aghq <- sqrt(rowSums(sapply(sim,function(x) (x[,3]- est_MNAR_aghq)^2))*(1/(it-1)))
emp.se_MNAR_ghq <- sqrt(rowSums(sapply(sim,function(x) (x[,4]- est_MNAR_ghq)^2))*(1/(it-1)))
emp.se_MNAR_gal <- sqrt(rowSums(sapply(sim,function(x) (x[,5]- est_MNAR_gal)^2))*(1/(it-1)))
emp.se_CC <- sqrt(rowSums(sapply(sim,function(x) (x[,6]- est_CC)^2))*(1/(it-1)))
emp.se_bd <- sqrt(rowSums(sapply(sim,function(x) (x[,7]- est_bd)^2))*(1/(it-1)))

av.se_MAR <- sqrt(rowMeans(sapply(sim,function(x) x[,29])))
av.se_MAR_twostage <- sqrt(rowMeans(sapply(sim,function(x) x[,30])))
av.se_MNAR_aghq <- sqrt(rowMeans(sapply(sim,function(x) x[,31])))
av.se_MNAR_ghq <- sqrt(rowMeans(sapply(sim,function(x) x[,32])))
av.se_MNAR_gal <- sqrt(rowMeans(sapply(sim,function(x) x[,33])))
av.se_CC <- sqrt(rowMeans(sapply(sim,function(x) x[,34])))
av.se_bd <- sqrt(rowMeans(sapply(sim,function(x) x[,35])))



## final results for chosen scenario ## 
res <- cbind(est_MAR, est_MAR_twostage, est_MNAR_aghq, est_MNAR_ghq,est_MNAR_gal, est_CC, est_bd, 
             bias_MAR, bias_MAR_twostage,bias_MNAR_aghq, bias_MNAR_ghq ,bias_MNAR_gal,  bias_CC, bias_bd, 
             cov_MAR, cov_MAR_twostage, cov_MNAR_aghq, cov_MNAR_ghq,cov_MNAR_gal, cov_CC, cov_bd,
             rmse.est_MAR, rmse.est_MAR_twostage, rmse.est_MNAR_aghq, rmse.est_MNAR_ghq, rmse.est_MNAR_gal,rmse.est_CC, rmse.est_bd,
             emp.se_MAR, emp.se_MAR_twostage, emp.se_MNAR_aghq, emp.se_MNAR_ghq,emp.se_MNAR_gal, emp.se_CC, emp.se_bd,
             av.se_MAR, av.se_MAR_twostage,av.se_MNAR_aghq,av.se_MNAR_ghq,av.se_MNAR_gal,av.se_CC,av.se_bd) 

