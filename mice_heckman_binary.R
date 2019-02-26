########################################################################
########################################################################
##
## MI function for imputing MNAR binary data 
## (usable with R package mice)
##
## @author Angelina Hammon
##
## 26.02.2019
##
########################################################################
########################################################################


# loading necessary packages:
require(MASS)
require(mvtnorm)
require(corpcor)
require(maxLik)
require(compiler)
require(Matrix)
require(VGAM)

# censored bivariate probit model:
mice.impute.heckman1step <- function(y, ry, x, type, draw=T, excl="",...) { 

  nam <- colnames(x) 
  
  pos.exc <- grep(excl,nam)  
  
  ind_s <- length(nam)+1   
  ind_y <- length(nam)-length(pos.exc)+1

  ry <- as.numeric(ry)
  data <- data.frame(x,y,ry)
  data$y <- as.numeric(as.character(data$y))
  colnames(data) <- c(colnames(x),"y","ry")
  
  
  ## Calculate one-step heckman model ##
  exc <- pos.exc 
  
  sel <- paste(nam, collapse="+")
  sel <- as.formula(paste("ry","~",sel,sep=""))
  out <- paste(nam[-exc], collapse="+")
  out <- as.formula(paste("y","~",out,sep=""))
  
  XS <- data.matrix(cbind(rep(1,length(y)),x))
  XO <- data.matrix(cbind(rep(1,length(y)),x[,-exc]))
  
  # estimating the model:
  heckit <- loglik_bivprob_est(data=data,yS=ry, yO=y, XO=XO, XS=XS, kS=ind_s, kO=ind_y,sel = sel, out=out) 
  
  
  ## Get parameter estimates of the heckman model ## 
  phi <- as.matrix(nearPD(ginv(heckit$hessian))$mat)
  q_hat <- heckit$estimates
  

  ## Draw ALL parameter values ##
  
  if(draw==T){
    
    q_star <- rmvnorm(1,q_hat,phi,method = "chol")
    beta_s_star <- q_star[1:ind_s]          
    beta_star <- q_star[(ind_s+1):(ind_s+ind_y)] 
    rho_star <- tanh(q_star[length(q_star)])
    
    cat("rho gezogen: ", rho_star, "  ","rho: ",tanh(q_hat[ind_s+ind_y+1]),"\n")
    
  }else{ ## just draw betas ##
    
    q_star <- rmvnorm(1,q_hat,phi,method = "chol")
    beta_s_star <- q_star[1:ind_s]          
    beta_star <- q_star[(ind_s+1):(ind_s+ind_y)] 
    rho_star <- tanh(q_hat[length(q_hat)])
    
  }
  

  ## Calculate prob based on the new parameter values ## 
  lin.predS <- XS %*% beta_s_star
  lin.predO <- XO %*% beta_star 
  
  lin.predS <- as.numeric(lin.predS)[!ry]
  lin.predO <- as.numeric(lin.predO)[!ry]
  
  
  Sigma <- matrix(-rho_star,ncol=2,nrow=2) 
  diag(Sigma) <- 1
  
  p <- sapply(1:sum(!ry),function(i) pmvnorm(c(-Inf,-Inf),c(lin.predO[i],-lin.predS[i]),mean=c(0,0),corr=Sigma)/pnorm(-lin.predS[i]))
  
  
  ## Draw Ymis based on new probability ##
  vec <- (runif(length(p)) <= p)
  vec[vec] <- 1
  if (is.factor(y)) {
    vec <- factor(vec, c(0, 1), levels(y))
  }
  
  print(table(vec))
  return(vec)
}

mice.impute.heckman1step <- compiler::cmpfun(mice.impute.heckman1step)



# estimation function: 
loglik_bivprob_est <- function(data, yS, yO, XO, XS, kS, kO, sel, out, ...) {
  
  r_start <- glm(sel,data=data,family=binomial(link=probit))
  
  y_start <- glm(out,data=data,family=binomial(link=probit))
  
  start <- c(coef(r_start),coef(y_start), atanh(0.5))
  
  opt <- optim(par=start, fn=loglik_bivprob, gr=loglik_bivprob_grad, method=c("BFGS"), yS, yO, XO, XS, kS, kO)
  
  hess <- numericHessian(loglik_bivprob,loglik_bivprob_grad,t0=opt$par, yS=yS, yO=yO, XO=XO, XS=XS, kS=kS, kO=kO)
  
  res <- list(estimates=opt$par,hessian=hess)
  
  return(res)
}


# likelihood: 
loglik_bivprob <- function(param, yS, yO, XO, XS, kS, kO) {
  
  betaS    <- param[1:kS]
  betaO    <- param[(kS+1):(kS+kO)]
  rho      <- tanh(param[kS+kO+1])
  
  
  # linear predictors:
  lin.predS <- XS %*% betaS
  lin.predO <- XO %*% betaO
  
  
  # Probabilities entering the likelihood function:
  p0 <- ifelse(yS==0,pnorm(-lin.predS),0)                                          ## P(yS=0)
  p11 <- ifelse(yS==1 & yO==1,pbinorm(lin.predS,lin.predO, cov12 = rho),0)         ## P(yS=1 & yO=1)
  p10 <- ifelse(yS==1 & yO==0,pbinorm(lin.predS,-lin.predO, cov12 = -rho),0)       ## P(yS=1 & yO=0)
  
  p10 <- ifelse(p10==0,1.1111e-15,p10)
  p11 <- ifelse(p11==0,1.1111e-15,p11)
  p0 <- ifelse(p0==0,1.1111e-15,p0)
  
  # Log-Likelihood function for every unit i:
  loglik <- (yS==0)*log(p0) + (yS==1 & yO==1)*log(p11) + (yS==1 & yO==0)*log(p10)
  
  
  # Total Log-Likelihood:
  return(-sum(loglik))
  
}


# analytic gradients: 
loglik_bivprob_grad <- function(param, yS, yO, XO, XS, kS, kO) {
  
  betaS    <- param[1:kS]
  betaO    <- param[(kS+1):(kS+kO)]
  rho      <- tanh(param[kS+kO+1])
  
  
  # linear predictors:
  lin.predS <- XS %*% betaS
  lin.predO <- XO %*% betaO
  
  p0 <- ifelse(yS==0,pnorm(-lin.predS),0)                                  ## P(yS=0)
  p11 <- ifelse(yS==1 & yO==1,pbinorm(lin.predS,lin.predO, cov12 = rho),0)         ## P(yS=1 & yO=1)
  p10 <- ifelse(yS==1 & yO==0,pbinorm(lin.predS,-lin.predO, cov12 = -rho),0)       ## P(yS=1 & yO=0)
  
  p10 <- ifelse(p10==0,1.1111e-15,p10)
  p11 <- ifelse(p11==0,1.1111e-15,p11)
  p0 <- ifelse(p0==0,1.1111e-15,p0)
  
  ### Gradients ###
  
  grad <- vector("numeric",length=8)
  
  r <- sqrt(1-rho^2)
  
  ## betaS:
  dp0 <-  apply(XS,2, function(x) x* (-dnorm(-(lin.predS))))
  dp10 <- apply(XS,2, function(x) x* dnorm(lin.predS) * pnorm((-(lin.predO)+rho*lin.predS)/r))
  dp11 <- apply(XS,2, function(x) x* dnorm(lin.predS) * pnorm((lin.predO-rho*lin.predS)/r))
  
  dprod <- (yS==0)*apply(dp0,2,function(x) x/p0) + (yS==1 & yO==1)*apply(dp11,2, function(x) x/p11) + (yS==1 & yO==0)*apply(dp10,2, function(x) x/p10) 
  
  dbetaS <- colSums(dprod)      
  grad[1:kS] <- dbetaS
  
  
  ## betaO:
  dp10 <- apply(XO,2, function(x) x* (-dnorm(-(lin.predO))) * pnorm((lin.predS - rho*lin.predO)/r))
  dp11 <- apply(XO,2, function(x) x* dnorm(lin.predO) * pnorm((lin.predS - rho*lin.predO)/r))
  
  dprod <- (yS==1 & yO==1)*apply(dp11,2,function(x) x/p11) + (yS==1 & yO==0)*apply(dp10,2, function(x) x/p10) 
  
  dbetaO <- colSums(dprod)    
  
  grad[(kS+1):(kS+kO)] <- dbetaO
  
  
  ## rho: 
  dr10 <- -(1-rho^2)*dbinorm(lin.predS,-(lin.predO), mean1 = 0, mean2 = 0, cov12 = -rho)
  dr11 <- (1-rho^2)*dbinorm(lin.predS,lin.predO,  mean1 = 0, mean2 = 0,cov12 = rho)
  
  dprod <- (yS==1 & yO==1)*(dr11/p11) + (yS==1 & yO==0)*(dr10/p10) 
  
  drho <- sum(dprod)
  
  grad[(kS+kO+1)] <- drho
  
  # return gradient vector:
  return(-grad)
  
}
