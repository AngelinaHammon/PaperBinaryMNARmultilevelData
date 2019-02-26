
#########################################################################################
#########################################################################################
##
## Application Bivariate Probit Imputation
## Correlation grade maths (binary: 1,2,3 vs. 4,5,6 ) and compet/reason/sex/mig/asp 
##
## @author Sabine Zinn & Angelina Hammon
##
## 26.02.2019
##
#########################################################################################
#########################################################################################



 
rm(list=ls())

setwd("Z:\\Abteilung_3_(03)\\Statistik_(04)\\MA_(02)\\_StGeIm\\AH\\MNAR_application")


library(readstata13)
library(ICCbin)
library(oglmx)
library(bivprob.quad)
library(mice)
library(mitml)
library(miceadds)
library(miceMNAR)
library(lme4)
library(mitools)
library(mitml)
library(BaylorEdPsych)
library(mvnmle)

source("mice_heckman_binary_re2.R")
source("mice_heckman_binary.R")


# ----------
# Read Data
# ----------
dat <- read.dta13("_datAnaExp_BivImp.dta")

# -----------------------------
# Complete Case Analysis (MCAR)
# -----------------------------
modICC <- glmer(aspH ~ (1|ID_i), family= binomial(link = "logit"), data=dat)
summary(modICC)

modCC_2L <- glmer(aspH ~ MN + DN + sex + 
                    + as.factor(lesC) + as.factor(mathC) + 
                    + UniAdm + mig + + UniAdmAv + 
                    + (1|ID_i), 
              family= binomial(link = "logit"), data=dat,
              glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
summary(modCC_2L) 

# --------------------------------------------------------------------
# Little's MCAR test
# --------------------------------------------------------------------
datTest <- dat[,c("MN", "DN", "sex", "aspH", "UniAdm", "mig", "lesC", "mathC")]    
littleTest <- LittleMCAR(datTest)
littleTest$p.value
# small p-value: much evidence against MCAR (H_0)

# --------------------------------------------------------------------
# Assess Exclusion criterion
# --------------------------------------------------------------------
# for aspiration
selYN <- ifelse(is.na(dat$aspH),1,0) 
tbl <- table(selYN, dat$indFeld)
chisq.test(tbl, correct=F) # very high correlation (measured by p-value), around 96%
tbl <- table(dat$aspH, dat$indFeld)
chisq.test(tbl, correct=F) # very low correlation (measured by p-value), nearly zero
# for maternal education
selYN <- ifelse(is.na(dat$UniAdm),1,0) 
tbl <- table(selYN, dat$indFeld)
chisq.test(tbl, correct=F) # moderate correlation (measured by p-value), around 46%
tbl <- table(dat$UniAdm, dat$indFeld)
chisq.test(tbl, correct=F) # low correlation (measured by p-value), around 17%

# ---------------
# Missing pattern
# ---------------
md.pattern(datTest)
tab <- md.pattern(datTest)[nrow(md.pattern(datTest)),]
round(tab / nrow(dat),2) # 54% in uni admission (uniAdm: Yes=1/No=0)

propNA <- function(vec){
  tab <- table(is.na(vec))
  ind <- which(names(tab)==TRUE)
  le <- length(vec)
  res <- ifelse(length(ind)==0,0,tab[ind]/le)
  return(res)
}
A <- aggregate(dat$UniAdm, list(ID_i=dat$ID_i), propNA)

namImp <- c("ID_t", "ID_i", "sex", "lesC", "mathC", "MN", "DN", "UniAdm", "aspH", "mig", "indFeld")
datI <- dat[,colnames(dat) %in% namImp]
md.pattern(datI) # missing values in aspH, mig, DN, MN, mathC, lesC, UniAdm

colnames(datI)[colnames(datI) %in% "ID_i"] <- "group"
datI$sex <- as.factor(datI$sex)
datI$UniAdm <- as.factor(datI$UniAdm)
datI$mathC <- as.factor(datI$mathC)
datI$lesC <- as.factor(datI$lesC)

# --------
# with MAR
# --------
ini <- mice(datI,m=1,maxit=0)
ini 
pred_MAR <- ini$pred
pred_MAR[,"ID_t"] <- 0
pred_MAR[,"group"] <- 0 
pred_MAR["aspH","group"] <- -2 # cluster id in the data set has to be assigned in the predictor matrix as -2
pred_MAR["mig","group"] <- -2 
meth <- ini$method
meth["aspH"] <- "2l.binary"
meth["mig"] <- "2l.binary"
meth["sex"] <- "logreg"
meth["UniAdm"] <- "logreg"
meth["MN"] <- "pmm"
meth["DN"] <- "pmm"
meth["lesC"] <- "polyreg"
meth["mathC"] <- "polyreg"

imp <- mice(datI,m=20,maxit=50,method=meth,pred=pred_MAR,seed=1234)

cc <- 0
getGlmer <- function(){
  cc <<- cc+1
  datIm <- complete(imp, action = cc)
  Z <- aggregate(as.numeric(as.character(datIm$UniAdm)), list(group=datIm$group), mean)
  colnames(Z)[2] <- "UniAdmAv"
  datIm <- merge(datIm, Z, by="group")
  model <- "aspH ~ MN + DN + sex + UniAdm + mig + as.factor(lesC) + as.factor(mathC) + UniAdmAv + (1|group)"
  return(lme4::glmer(formula=model, data=datIm, family = binomial(link="logit"), 
                       glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))))
}
fit <- with(imp, getGlmer())
summary(pool(fit))
testEstimates(as.mitml.result(fit), var.comp = TRUE)$var.comp

# -----
# with MNAR
# ---------

datI$r1 <- ifelse(is.na(datI$aspH),0,1) # indicator showing whether a values in MNAR variable has been observed (=1), or not (=0)
datI$r2 <- ifelse(is.na(datI$UniAdm),0,1) # indicator showing whether a values in MNAR variable has been observed (=1), or not (=0)
ini <- mice(datI,m=1,maxit=0) 
ini 
pred_MNAR <- ini$pred
pred_MNAR[,"ID_t"] <- 0
pred_MNAR[,"group"] <- 0
pred_MNAR["aspH","group"] <- -2  # cluster id in the data set has to be assigned in the predictor matrix as -2
pred_MNAR["mig","group"] <- -2 
# exclusion criterion has to be assigned in the predictor matrix as 1 and 
# its variable name has to be specified with the excl="..." option in the mice function
pred_MNAR["aspH","indFeld"] <- 1
pred_MNAR["UniAdm","indFeld"] <- 1
# for multivariate missing data use missing indicator of variable 
# that is supposed to be MNAR as predictor for other incomplete variables
pred_MNAR["aspH","r1"] <- 0
pred_MNAR["UniAdm","r1"] <- 1 
pred_MNAR["mig","r1"] <- 1
pred_MNAR["MN","r1"] <- 1
pred_MNAR["DN","r1"] <- 1
pred_MNAR["sex","r1"] <- 1
pred_MNAR["lesC","r1"] <- 1
pred_MNAR["mathC","r1"] <- 1
pred_MNAR["UniAdm","r2"] <- 0
pred_MNAR["aspH","r2"] <- 1 
pred_MNAR["mig","r2"] <- 1
pred_MNAR["MN","r2"] <- 1
pred_MNAR["DN","r2"] <- 1
pred_MNAR["sex","r2"] <- 1
pred_MNAR["lesC","r2"] <- 1
pred_MNAR["mathC","r2"] <- 1

meth <- ini$method
meth["UniAdm"] <- "heckman1step"
meth["aspH"] <- "2l.heckman1step_re_aghq"
meth["mig"] <- "2l.binary"
meth["sex"] <- "logreg"
meth["MN"] <- "pmm"
meth["DN"] <- "pmm"
meth["lesC"] <- "polyreg"
meth["mathC"] <- "polyreg"

# Input parameters of biv-probit imputation model:
# - draw: set TRUE if parameters should be drawn from approximate normal posterior
# - H: number of quadrature points for the two random intercepts (needs to be a vector of length 2)
# - excl: give variable name of exclusion criterion as a string 

imp <- mice(datI,m=20,maxit=50,method=meth,pred=pred_MNAR,draw=T, excl="indFeld", seed=1234)

cc <- 0
getGlmer <- function(){
  cc <<- cc+1
  datIm <- complete(imp, action = cc)
  Z <- aggregate(as.numeric(as.character(datIm$UniAdm)), list(group=datIm[,"group"]), mean)
  colnames(Z)[2] <- "UniAdmAv"
  datIm <- merge(datIm, Z, by="group")
  model <- "aspH ~ MN + DN + as.factor(mathC) + as.factor(lesC) + sex + UniAdm + mig + UniAdmAv + (1|group)"
  return(lme4::glmer(formula=model, data=datIm, family = binomial(link="logit"), 
                     glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))))
}
fit <- with(imp, getGlmer())
summary(pool(fit))
testEstimates(as.mitml.result(fit), var.comp = TRUE)$var.comp
