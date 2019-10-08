# main function for multilevel imputation with lme4 which
# is just wrapper by methods "2l.continuous", "2l.binary" and "2l.pmm"
mice.impute.2l.lmer <- function(y, ry, x, type, intercept=TRUE,
                                groupcenter.slope=FALSE, draw.fixed=TRUE,
                                random.effects.shrinkage=1E-6,
                                glmer.warnings=TRUE, 
                                model = "continuous" , donors = 3 ,
                                match_sampled_pars = FALSE ,
                                blme_use = FALSE ,
                                blme_args = NULL , 
                                ...){
  
  # *** ...............................
  # preliminary calculations
  clus <- x[,type==-2]
  clus_unique <- unique(clus)
  ngr <- length(clus_unique)
  clus_name <- colnames(x)[type==-2]  # name of cluster identifier
  
  #zz0 <- Sys.time()	
  
  # arguments for lmer model	
  if ( model == "binary"){	
    lmer_family <- stats::binomial(link="logit")
    if ( blme_use){
      lmer_function <- blme::bglmer						
    } else {
      lmer_function <- lme4::glmer
    }
  }
  if ( model %in% c("continuous","pmm") ){	
    if ( blme_use){
      lmer_function <- blme::blmer						
    } else {
      lmer_function <- lme4::lmer
    }		
  }	
  
  #--- add group means (if needed)
  res <- mice_multilevel_add_groupmeans( y=y , ry=ry , x=x , type = type ,
                                         groupcenter.slope = groupcenter.slope)
  x <- res$x
  type <- res$type
  #--- create formulas for lme4
  rhs.f <- mice_multilevel_create_formula( 
    variables = colnames(x)[type %in% c(1,2)] , 
    include_intercept = intercept )
  rhs.r <- mice_multilevel_create_formula( 
    variables = colnames(x)[type == 2] , 
    include_intercept = TRUE )	
  # combine formulas
  fml <- paste0( "dv._lmer~", rhs.f, "+(", rhs.r,"|", clus_name ,")" )
  
  #*** prepare arguments for lmer estimation
  y1 <- y
  y1[!ry] <- NA
  dat_lme4 <- data.frame(dv._lmer=y1, x)
  lmer_args <- list( formula = fml , data = dat_lme4 ,
                     na.action = "na.omit"  )
  if ( model == "binary"){
    lmer_args$family <- lmer_family
  }
  
  # apply blme arguments if provided
  lmer_args <- mice_multilevel_imputation_blme_args(
    lmer_args = lmer_args , blme_args = blme_args )
  
  # fit based on observed y
  fit <- mice_multilevel_doCall_suppressWarnings( 
    what = lmer_function , args = lmer_args , 
    warnings = glmer.warnings )
  
  # clusters without missing values
  clus0 <- clus[!ry]	
  
  #--- draw fixed effects
  b.est <- b.star <- lme4::fixef(fit)
  if( draw.fixed ){     # posterior draw for fixed effects		
    b.star <- mice_multilevel_draw_rnorm1( mu = b.star , Sigma = vcov(fit) )
  } 
  
  #--- extract posterior distribution of random effects
  fl <- lme4::getME(fit, "flist")[[1]]	
  # ind <- match( clus0, clus_unique)
  index_clus <- match( clus, clus_unique)
  # clusters with at least one observation
  clus_obs <- match( unique(fl), clus_unique)
  fit_VarCorr <- lme4::VarCorr(fit)	
  vu <- fit_VarCorr[[1]][,,drop=FALSE ]     # random effects (co-)variance
  # extract random effects
  re0 <- lme4::ranef(fit, condVar=TRUE)[[1]]
  NR <- ncol(re0)
  re <- matrix(0, nrow=ngr, ncol=NR)     # re: 0 if fully unobserved
  re[clus_obs,] <- as.matrix(re0)        # re: EAP if partially observed
  
  pv0 <- attr(re0, "postVar")
  pv <- array(0, dim=c(NR,NR,ngr))
  pv[,,clus_obs] <- pv0                # pv: post. variance if partially observed
  pv[,,-clus_obs] <- vu                # pv: random effects cov. if fully unobserved
  
  #--- draw random effects
  u <- mice_multilevel_imputation_draw_random_effects( mu = re , Sigma = pv  ,
                                                       ridge = random.effects.shrinkage )
  
  #--- x and z for prediction
  x0 <- as.matrix( x[ ,type>=1,drop=FALSE ] )
  z0 <- as.matrix( x[ ,type==2,drop=FALSE ] )
  if(intercept){
    x0 <- cbind(1,x0)
    z0 <- cbind(1,z0)
  }	
  
  #--- compute predicted values including fixed and random part
  predicted <- x0 %*% b.star + rowSums( z0 * u[index_clus ,1:NR,drop=FALSE]) 
  # predicted values for cases with missing data
  predicted0 <- predicted[ !ry ]
  # predicted values for cases with observed data
  if ( model == "pmm"){
    # use non-sampled values here, see corresponding mice approach
    # after this function
    if (match_sampled_pars){
      # non-mice approach		
      pred <- x0 %*% b.star + rowSums( z0 * u[index_clus ,1:NR,drop=FALSE] )
    } else {
      # "the mice approach"
      pred <- x0 %*% b.est + rowSums( z0 * re[index_clus ,1:NR,drop=FALSE] )		
    }				
    predicted1 <- pred[ ry ]	
  }
  
  #---- draw imputations
  if ( model == "binary"){	
    imp <- mice_multilevel_draw_binomial( probs = antilogit(predicted0) )
  }
  if ( model == "continuous"){	
    sigma <- attr( fit_VarCorr ,"sc")
    imp <- mice_multilevel_imputation_draw_residuals(
      predicted = predicted0 , sigma = sigma  )		
  }
  if ( model == "pmm"){
    imp <- mice_multilevel_imputation_pmm5(y=y, ry=ry, x, 
                                           yhatobs= predicted1 , yhatmis = predicted0, 
                                           donors=donors , noise = 1E5 , ...)		
  }
  #--- output imputed values
  return(imp)
}

mice_multilevel_add_groupmeans <- function( y , ry , x , type ,
                                            groupcenter.slope ){
  # add groupmeans in the regression model
  if ( any( type %in% c(3,4) ) ){ 
    x0 <- as.matrix(cbind( x[,type==-2], x[,type %in% c(3,4)] ))
    colnames(x0) <- c( colnames(x)[type==-2], colnames(x)[type %in% c(3,4)] )
    type0 <- c( -2, rep(1,ncol(x0)-1) )
    x0.aggr <- as.matrix( mice_multilevel_impute_groupmean(y=y, ry=ry, x=x0, 
                                                           type=type0, grmeanwarning=FALSE ))
    colnames(x0.aggr) <- paste0("M._", colnames(x0)[-1])
    # group mean centering
    if ( groupcenter.slope ){ 
      x0.aggr1 <- as.matrix(x0.aggr)
      colnames(x0.aggr1) <- colnames(x0)[-1]
      x0cent <- x0[,-1] - x0.aggr1
      x[ , colnames(x0cent) ] <- x0cent
    }
    # combine covariate matrix
    x <- cbind( x , x0.aggr )
    # add type
    type1 <- c( type , rep(1 , ncol(x0.aggr) ) )
    names(type1) <- c( names(type) , colnames(x0.aggr) )   
    type1[ type1 == 3 ] <- 1
    type1[ type1 == 4 ] <- 2
    type <- type1
  }		
  res <- list( "x" = x , "type" = type)	
  return(res)	
}	


mice_multilevel_create_formula <- function( 
  variables , include_intercept ){
  #---
  intercept_code <- if ( include_intercept ){ 1 } else { 0 }
  fml <- paste0( c( intercept_code , variables ), 
                 collapse="+" )
  return(fml)
}

mice_multilevel_imputation_blme_args <- function(lmer_args , blme_args){
  if ( ! is.null( blme_args) ){
    NL <- length(blme_args)
    for (nn in 1:NL){
      name_nn <- names(blme_args)[nn]
      lmer_args[[ name_nn ]] <- blme_args[[ nn ]]	
    }
  }
  return(lmer_args)
}

mice_multilevel_doCall_suppressWarnings <- function( what , args , warnings = TRUE){
  if (warnings){
    res <- do.call( what = what , args = args)
  } else {
    suppressWarnings(
      res <- do.call( what = what , args = args)
    )
  }
  return(res)
}	

mice_multilevel_draw_rnorm1 <- function( mu  , Sigma){	
  #----	
  #b.star <- b.star + as.vector( t(chol(vcov(fit))) %*% rnorm(length(b.star)) )		
  NP <- length(mu)	
  res <- mu + as.vector( t( chol(Sigma) %*% stats::rnorm(NP) ) )	
  return(res)
}

mice_multilevel_imputation_draw_random_effects <- function( mu , Sigma ,
                                                            ridge = 1E-50 ){
  
  dim_Sigma <- dim(Sigma)
  ngr <- dim_Sigma[3]
  NR <- dim_Sigma[1]
  # draw random effects
  u <- matrix(0, nrow=ngr, ncol=NR)
  if (NR==1){
    u[,1] <- stats::rnorm(ngr, mean=mu[,1], sd= sqrt(Sigma[1,1,]) )
  } else {
    for(i in 1:ngr){
      #-- compute covariance matrix with ridge		
      Sigma1 <- Sigma[,,i] + diag(ridge,NR)
      # Cholesky matrix of Sigma1
      Sigma_chol <- chol(Sigma1)
      # simulation
      rn <- stats::rnorm(NR, mean=0 , sd=1)
      u[i,] <- mu[i,] + as.vector( t( Sigma_chol ) %*% rn )
    }
  }
  return(u)
}	

mice_multilevel_draw_binomial <- function( probs ){	
  N <- length(probs)
  rn <- stats::runif(N, 0, 1)
  res <- 1*(rn < probs)
  return(res)
}

antilogit <- function(p){
  return( 1 / ( 1 + exp(-p) ) )
}


mice.impute.2l.binary <- function (y, ry, x, type, intercept = TRUE, groupcenter.slope = FALSE, 
          draw.fixed = TRUE, random.effects.shrinkage = 1e-06, glmer.warnings = TRUE, 
          blme_use = FALSE, blme_args = NULL, ...) {
  imp <- mice.impute.2l.lmer(y = y, ry = ry, x = x, type = type, 
                             intercept = intercept, groupcenter.slope = groupcenter.slope, 
                             draw.fixed = draw.fixed, random.effects.shrinkage = random.effects.shrinkage, 
                             glmer.warnings = glmer.warnings, blme_use = blme_use, 
                             blme_args = blme_args, model = "binary", ...)
  return(imp)
}