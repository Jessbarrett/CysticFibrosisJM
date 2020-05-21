




rm(list=ls())

# Load required packages
library(nlme)
library(mnormt)
library(plyr)
library(tensor)
library(numDeriv)



################## FIT THE JOINT MODEL TO A RANDOM SAMPLE #######################
## Get initial values for full sample by fitting the model to a random sample ###

### Get the random sample
# The file CF_JM_functions.R contains all the functions needed for setting up the data 
# and fitting the model 
source("CF_JM_functions.R")
# getdat is a function which creates the R data frames
temp <- getdat(whiteonly=T, rs=T, fortyplus=F, datafile="JM_dataset3.RData")
cfdat.rs <- temp$data
cfsurv.rs <- temp$survdata
rm(temp)

# Set up the data and formulae for model fit
dat.y <- with(cfdat.rs,data.frame(id,fev,ageatvisit,birthyear_c,sex,F508.0.1.NT))
dat.y$ageatvisitsq <- dat.y$ageatvisit^2
# Formula for longitudinal model
fmla.y <- (fev ~ ageatvisit + ageatvisitsq + birthyear_c + sex + F508.0.1.NT + sex:ageatvisit + F508.0.1.NT:ageatvisit
           + sex:birthyear_c + F508.0.1.NT:birthyear_c )
dat.t <- with(cfsurv.rs,data.frame(id,status,time,startage,birthyear_c,sex,F508.0.1.NT))

# Remove 3 individuals with time=0 (these are individuals whose death or censoring date occurs before their first visit date)
ids.to.remove <- dat.t$id[dat.t$time==0]
dat.t <- dat.t[!(dat.t$id %in% ids.to.remove),]
dat.y <- dat.y[!(dat.y$id %in% ids.to.remove),]
dat.t$id <- factor(as.character(dat.t$id))
dat.y$id <- factor(as.character(dat.y$id))

# fitmodel fits the joint model
cffit.rs <- fitmodel(Dat=dat.y, longformula=fmla.y, Survdat=dat.t, intervaltimes=((0:19)+0.5), printpars=T, method.lme="optim", maxiter=10000)
save(list=ls(),file="cffit_RS.RData")
summary.modelfit(cffit.rs.epinet, digits=3)



############## FIT THE JOINT MODEL TO THE FULL DATASET #############################
############## USING RANDOM SAMPLE ESTIMATES AS INITIAL VALUES #####################

rm(list=ls())

# Use fit from old data for initial standard errors
# !!!!! Put the values in here
load(file=paste(path_output_olddata,"cffit_9_new.Rdata",sep=""))
initses <- cffit.9$se
rm(list=setdiff(ls(),c("initses")))
# Use fit to random sample for intial parameter values
load(file="cffit_RS.RData")
inits <- cffit.rs$estimate
rm(list=setdiff(ls(),c("inits","initses")))

# Now load the data
source("CF_JM_functions.R")
temp <- getdat(whiteonly=T, rs=F, fortyplus=F, datafile="JM_dataset3.RData")
cfdat <- temp$data
cfsurv <- temp$survdata
rm(temp)

# Set up the data and formulae for model fit
dat.y <- with(cfdat,data.frame(id,fev,ageatvisit,birthyear_c,sex,F508.0.1.NT))
dat.y$ageatvisitsq <- dat.y$ageatvisit^2
# Formula for longitudinal model
fmla.y <- (fev ~ ageatvisit + ageatvisitsq + birthyear_c + sex + F508.0.1.NT + sex:ageatvisit + F508.0.1.NT:ageatvisit
                       + sex:birthyear_c + F508.0.1.NT:birthyear_c)
dat.t <- with(cfsurv,data.frame(id,status,time,startage,birthyear_c,sex,F508.0.1.NT))

# Remove 13 individuals with time=0 (these are individuals whose death or censoring date occurs before their first visit date)
ids.to.remove <- dat.t$id[dat.t$time==0]
dat.t <- dat.t[!(dat.t$id %in% ids.to.remove),]
dat.y <- dat.y[!(dat.y$id %in% ids.to.remove),]
dat.t$id <- factor(as.character(dat.t$id))
dat.y$id <- factor(as.character(dat.y$id))

# fitmodel fits the joint model
system.time(cffit <- fitmodel(Dat=dat.y, longformula=fmla.y, Survdat=dat.t, intervaltimes=((0:19)+0.5), printpars=T, pmnormtol=1e-4, inits=inits, initses=initses))
save(list=ls(),file="cffit.Rdata")

# Display results
summary.modelfit(cffit, digits=3)



############### USE SMALLER PMNORM TOLERANCE TO GET STANDARD ERRORS ###################

rm(list=ls())
# Load data
load(file="cffit.Rdata")

# Scaling of standard errors used for likelihood maximisation
scaleses <- cffit.epinet$estimate/cffit.epinet$fit$par
# pmnorm tolerance to use for standard error estimation
pmnormtolforse <- 1e-8
# Get the Hessian
hess <- hessian(func=loglik, x=cffit.epinet$fit$par, dat=dat.y, fmla=fmla.y, surv=dat.t,  ses=scaleses,intervaltimes=((0:19)+0.5),
                      printpars=F, pmnormtol=pmnormtolforse)
# Use the nearPD function to get a positive defininte estimate for the hessian 
library(Matrix)
hessinv <- solve(hess)
nearpdfit <- nearPD(hessinv)
hessinvnew <- nearpdfit$mat
se <- sqrt(diag(hessinvnew))*scaleses
# Save standard error estimates in cffit
cffit.epinet$se.new <- se.test
# Save results
save(list=ls(),file="cffit.Rdata")

## Display results
summary.modelfit(cffit.epinet, digits=3) 



########### CHECK RESULT BY RE-RUNNING USING DIFFERENT INITIAL VALUES ########################
###### Starting at different initial values gives similar results                      #######

# Load the data
source("CF_JM_functions.R")
temp <- getdat(whiteonly=T, rs=F, fortyplus=F, datafile="JM_dataset3.RData")
cfdat <- temp$data
cfsurv <- temp$survdata
rm(temp)

# Set up the data and formulae for model fit
dat.y <- with(cfdat,data.frame(id,fev,ageatvisit,birthyear_c,sex,F508.0.1.NT))
dat.y$ageatvisitsq <- dat.y$ageatvisit^2
# Formula for longitudinal model
fmla.y <- (fev ~ ageatvisit + ageatvisitsq + birthyear_c + sex + F508.0.1.NT + sex:ageatvisit + F508.0.1.NT:ageatvisit
           + sex:birthyear_c + F508.0.1.NT:birthyear_c)
dat.t <- with(cfsurv,data.frame(id,status,time,startage,birthyear_c,sex,F508.0.1.NT))

# Remove 13 individuals with time=0 (these are individuals whose death or censoring date occurs before their first visit date)
ids.to.remove <- dat.t$id[dat.t$time==0]
dat.t <- dat.t[!(dat.t$id %in% ids.to.remove),]
dat.y <- dat.y[!(dat.y$id %in% ids.to.remove),]
dat.t$id <- factor(as.character(dat.t$id))
dat.y$id <- factor(as.character(dat.y$id))

# Joint model using initial values estimated by fitmodel
system.time(cffitnew <- fitmodel(Dat=dat.y, longformula=fmla.y, Survdat=dat.t, intervaltimes=((0:19)+0.5), printpars=T, pmnormtol=1e-4))
save(list=ls(),file="cffitnew.Rdata")


################ FITTING INDEPENDENT SURVIVAL MODELS ###########################


rm(list=ls())

# Load the data
source("CF_JM_functions.R")
temp <- getdat(whiteonly=T, rs=F, fortyplus=F, datafile="JM_dataset3.RData")
cfdat <- temp$data
cfsurv <- temp$survdata
rm(temp)

# Set up the data and formulae for model fit
dat.y <- with(cfdat,data.frame(id,fev,ageatvisit,birthyear_c,sex,F508.0.1.NT))
dat.y$ageatvisitsq <- dat.y$ageatvisit^2
dat.t <- with(cfsurv,data.frame(id,status,time,startage,birthyear_c,sex,F508.0.1.NT))

# Remove 13 individuals with time=0 (these are individuals whose death or censoring date occurs before their first visit date)
ids.to.remove <- dat.t$id[dat.t$time==0]
dat.t <- dat.t[!(dat.t$id %in% ids.to.remove),]
dat.y <- dat.y[!(dat.y$id %in% ids.to.remove),]
dat.t$id <- factor(as.character(dat.t$id))
dat.y$id <- factor(as.character(dat.y$id))

# Mid-points of time intervals 
intervaltimes=((0:19)+0.5)
# Add baseline FEV to survival data frame
indfev <- match(dat.t$id,dat.y$id)
basefev <- dat.y$fev[indfev]
dat.t <- cbind(dat.t, basefev)
# survival dataset with one row per year interval
ind  <- rep(1:nrow(dat.t),dat.t$time)
dat.tnew <- dat.t[ind,]
# generate a time interval variable
temp <- rep(1,nrow(dat.tnew))
dat.tnew$interval <- unlist(tapply(temp,dat.tnew$id,cumsum))
dat.tnew$intervalage <- intervaltimes[dat.tnew$interval]+dat.tnew$startage
# generate the outcome for probit analysis of survival data
dat.tnew$Yes <- 1
dat.tnew$Yes[dat.tnew$status==0 & dat.tnew$interval==dat.tnew$time] <- 0

# Survival model without FEV as a predictor
survfit1 <- glm(Yes ~ intervalage + birthyear_c + sex + F508.0.1.NT, family=binomial(link="probit"), data=dat.tnew)
# Survival model with baseline FEV as a predictor
survfit2 <- glm(Yes ~ basefev + intervalage + birthyear_c + sex + F508.0.1.NT, family=binomial(link="probit"), data=dat.tnew)

## adjust also for FEV as time-dependent covariate
# Get the time interval of each visit in the FEV data
reps <- table(dat.y$id)
startage <- tapply(dat.y$ageatvisit,dat.y$id,min)
visittime <- dat.y$ageatvisit - rep(startage,reps)
timeinterval.y <- as.integer(floor(visittime)+1)
# Gives 9741x19 matrix of mean FEV in each time interval for each id, NA if no FEV measurement taken
meanfevbyinterval <- tapply(dat.y$fev,data.frame(dat.y$id,timeinterval.y),mean)
# Add an extra column of NA's for time interval 20
meanfevbyinterval <- cbind(meanfevbyinterval,NA)
# create vector of FEV means from meanfevbyinterval, truncated by survival time for each id
meanfev <- NULL
for(i in 1:nrow(meanfevbyinterval)){
  meanfev <- c(meanfev,meanfevbyinterval[i,1:dat.t$time[i]])
}
# replace NA's with Last Observation Carried Forward
fevlocf <- meanfev
while(any(is.na(fevlocf))==T){
  indna <- which(is.na(fevlocf))
  fevlocf[indna] <- fevlocf[indna-1]
}
dat.tnew <- cbind(dat.tnew,fevlocf)
survfit3 <- glm(Yes ~ fevlocf + intervalage + birthyear_c + sex + F508.0.1.NT, family=binomial(link="probit"), data=dat.tnew)

summary(survfit1)
summary(survfit2)
summary(survfit3)
confint(survfit1)
confint(survfit2)
confint(survfit3)


################### FITTING A TWO-STAGE MODEL ###########################


rm(list=ls())

# Load the data
source("CF_JM_functions.R")
temp <- getdat(whiteonly=T, rs=F, fortyplus=F, datafile="JM_dataset3.RData")
cfdat <- temp$data
cfsurv <- temp$survdata
rm(temp)

# Set up the data and formulae for model fit
dat.y <- with(cfdat,data.frame(id,fev,ageatvisit,birthyear_c,sex,F508.0.1.NT))
dat.y$ageatvisitsq <- dat.y$ageatvisit^2
dat.t <- with(cfsurv,data.frame(id,status,time,startage,birthyear_c,sex,F508.0.1.NT))

# Remove 13 individuals with time=0 (these are individuals whose death or censoring date occures before their first visit date)
ids.to.remove <- dat.t$id[dat.t$time==0]
dat.t <- dat.t[!(dat.t$id %in% ids.to.remove),]
dat.y <- dat.y[!(dat.y$id %in% ids.to.remove),]
dat.t$id <- factor(as.character(dat.t$id))
dat.y$id <- factor(as.character(dat.y$id))

intervaltimes=((0:19)+0.5)
# survival dataset with one row per year interval
ind  <- rep(1:nrow(dat.t),dat.t$time)
dat.t.2s <- dat.t[ind,]
# generate a time interval variable
temp <- rep(1,nrow(dat.t.2s))
dat.t.2s$interval <- unlist(tapply(temp,dat.t.2s$id,cumsum))
dat.t.2s$intervalage <- intervaltimes[dat.t.2s$interval]+dat.t.2s$startage
# generate a ageatvisit variable equal to intervalage
dat.t.2s$ageatvisit <- dat.t.2s$intervalage
dat.t.2s$ageatvisitsq <- dat.t.2s$ageatvisit^2

# generate the outcome for probit analysis of survival data
dat.t.2s$Yes <- 1
dat.t.2s$Yes[dat.t.2s$status==0 & dat.t.2s$interval==dat.t.2s$time] <- 0

fmla.y <- (fev ~ ageatvisit + ageatvisitsq + birthyear_c + sex + F508.0.1.NT + sex:ageatvisit + F508.0.1.NT:ageatvisit
                                   + sex:birthyear_c + F508.0.1.NT:birthyear_c)
longfit <- lme(fmla.y, random=~ageatvisit|id, data=dat.y, control=list(opt="optim"))
dat.t.2s$fevpred <- predict(longfit, newdata=dat.t.2s)
survfit4 <- glm(Yes ~ fevpred + intervalage + birthyear_c + sex + F508.0.1.NT, family=binomial(link="probit"), data=dat.t.2s)
summary(survfit4)
summary(longfit)


