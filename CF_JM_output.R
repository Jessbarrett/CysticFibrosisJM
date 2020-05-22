#######################################################################################
##################    POPULATION CHARACTERISTICS TABLE   ##############################

rm(list=ls())

# Load the data
source("CF_JM_functions.R")
temp <- getdat(whiteonly=T, rs=F, fortyplus=F, datafile="JM_dataset3.RData")
cfdat <- temp$data
cfsurv <- temp$survdata
rm(temp)

# Dichotomise IMD
cfsurv$imd_bin <- as.numeric(cfsurv$imd<0)
cfdat$imd_bin <- as.numeric(cfdat$imd<0)

# Get age at diagnosis data from separate file 
load(paste(path_epidata,"IdsDiagnosisAge.RData",sep=""))
# Select the correct dataframe from the list dataList 
ageatdiagdata <- dataList[[1]]
hist(ageatdiagdata$dmg_ageatdiagnosis)
#### Note 87 have missing age at diagnosis
sum(is.na(ageatdiagdata$dmg_ageatdiagnosis))
# Add age at diagnosis to longitudinal dataset
cfdat$ageatdiag <- NA
cfdat$ageatdiag <- ageatdiagdata$dmg_ageatdiagnosis[match(cfdat$id,ageatdiagdata$S01CaseId_Original)]
# Add age at diagnosis to survival dataset
cfsurv$ageatdiag <- NA
cfsurv$ageatdiag <- ageatdiagdata$dmg_ageatdiagnosis[match(cfsurv$id,ageatdiagdata$S01CaseId_Original)]
rm(dataList, ageatdiagdata)

# Remove 13 individuals with time=0 (these are individuals whose death or censoring date occurs before their first visit date)
ids.to.remove <- cfsurv$id[cfsurv$time==0]
cfsurv <- cfsurv[!(cfsurv$id %in% ids.to.remove),]
cfdat <- cfdat[!(cfdat$id %in% ids.to.remove),]
cfsurv$id <- factor(as.character(cfsurv$id))
cfdat$id <- factor(as.character(cfdat$id))

dat.t <- with(cfsurv,data.frame(id,status,time,startage,birthyear_c,sex,F508.0.1.NT,imd,imd_bin,ageatdiag))
dat.y <- with(cfdat,data.frame(id,fev,ageatvisit,birthyear_c,sex,F508.0.1.NT))

# No. patients
dim(dat.t)
length(unique(dat.y$id))
# No. fev measurements
dim(dat.y)
any(is.na(dat.y$fev))
# No. fev measurements per person
nfev <- tapply(dat.y$fev,dat.y$id,length)
median(nfev)
min(nfev)
max(nfev)
mean(nfev)
sd(nfev)
# Percent. with X or more measurements
X <- 3
sum(nfev >= X)*100/length(nfev)
# person-years of follow-up     
sum(cfsurv$exacttime)
# Baseline characteristics by sex
table(dat.t$sex)    # 0 = Female, 1 = Male
prop.table(table(dat.t$sex))*100
table(dat.y$sex)
prop.table(table(dat.y$sex))*100
# No FEV measurements per person by sex
tapply(nfev,dat.t$sex,mean)
tapply(nfev,dat.t$sex,sd)
# Deaths
table(dat.t$status)    # 0 = Dead, 1 = Censored
table(dat.t$status,dat.t$sex,dnn=c("Status","Sex"))
prop.table(table(dat.t$status,dat.t$sex,dnn=c("Status","Sex")),margin=1)*100
# F508 alleles              
table(cfsurv$F508_class)
prop.table(table(cfsurv$F508_class))*100
table(cfsurv$F508_class,cfsurv$sex,dnn=c("F508","Sex"))
prop.table(table(cfsurv$F508_class,cfsurv$sex,dnn=c("F508","Sex")),margin=2)*100
# Birth cohort
table(cfsurv$birthcohort)               
prop.table(table(cfsurv$birthcohort))*100
table(cfsurv$birthcohort,cfsurv$sex)
prop.table(table(cfsurv$birthcohort,cfsurv$sex),margin=2)*100
# Age at entry into study
mean((cfsurv$startage+5))
sd(cfsurv$startage)
median(cfsurv$startage+5)
tapply(cfsurv$startage+5,cfsurv$sex,mean)
tapply(cfsurv$startage+5,cfsurv$sex,sd)
# IMD Z score
sum(is.na(dat.t$imd))
mean(dat.t$imd,na.rm=T)
sd(dat.t$imd,na.rm=T)
tapply(cfsurv$imd,cfsurv$sex,mean,na.rm=T)
tapply(cfsurv$imd,cfsurv$sex,sd,na.rm=T)
# Low IMD Z score
table(dat.t$imd_bin)
prop.table(table(dat.t$imd_bin))*100
table(cfsurv$imd_bin,cfsurv$sex,dnn=c("IMD","sex"))
prop.table(table(cfsurv$imd_bin,cfsurv$sex,dnn=c("IMD","sex")),margin=2)*100
# Age at diagnosis
sum(is.na(dat.t$ageatdiag))
mean(dat.t$ageatdiag,na.rm=T)
sd(dat.t$ageatdiag,na.rm=T)
median(dat.t$ageatdiag,na.rm=T)
tapply(cfsurv$ageatdiag,cfsurv$sex,mean,na.rm=T)
tapply(cfsurv$ageatdiag,cfsurv$sex,sd,na.rm=T)

# Histogram of age at entry into study
hist(dat.t$startage, main="", xlab="Age at study entry")


##################################################################################################
################			PROBABILITY RATIOS AND STANDARD ERRORS		##########################

rm(list=ls())
load(file="cffit.Rdata")

### Probability ratio for difference in current level of FEV

# Parameter estimates and covariance matrix for probit linear predictor
alpha <- cffit.epinet$estimate[c(11:15,20:21)]
Vfull <- scaleses%*%t(scaleses) * hess.inv.test
Valpha <- Vfull[c(11:15,20:21),c(11:15,20:21)]

# Covariate values for 20-year old male born in 1980 with 2 F508 alleles 
age <- 15.5
byr <- -6
sex <- 1
f508 <- 1
X = c(1,age,byr,sex,f508)

# Current FEV, given FEV intercept = pop average, FEV slope = pop average
Xlong = c(1,age,age^2,byr,sex,f508,age*sex,age*f508,byr*sex,byr*f508)
beta <- cffit.epinet$estimate[1:10]
currentfev1 <- sum(Xlong*beta)
# Current FEV, given FEV intercept = pop average-10, FEV slope = pop average
currentfev2 <- currentfev1-10
# FEV slope
dXlong = c(0,1,2*age,0,0,0,sex,f508,0,0)
fevslope <- sum(dXlong*beta)

# Probabilities of one-year survival
Xsurv1 <- c(X,currentfev1,fevslope)
Xsurv2 <- c(X,currentfev2,fevslope)
prob1 <- pnorm(sum(Xsurv1*alpha))
prob2 <- pnorm(sum(Xsurv2*alpha))

# Relative risk
rr <- (1-prob2)/(1-prob1)
# Standard errors for probit linear predictors
xb1var <- as.numeric(Xsurv1%*%Valpha%*%Xsurv1)
xb2var <- as.numeric(Xsurv2%*%Valpha%*%Xsurv2)
# Covariance between xb1 and xb2
xb1xb2cov <- as.numeric(Xsurv1%*%Valpha%*%Xsurv2)
# Use delta method to get standard error for the log relative risk (note derivative of pnorm is dnorm)
xbV <- matrix(c(xb1var,xb1xb2cov,xb1xb2cov,xb2var),nrow=2)
gradrr1 <- dnorm(sum(Xsurv1*alpha))/(1-prob1)
gradrr2 <- -dnorm(sum(Xsurv2*alpha))/(1-prob2)
gradrr <- c(gradrr1,gradrr2)
rrvar <- gradrr%*%xbV%*%gradrr
# 95% Confidence interval
logrrlower <- log(rr)-qnorm(0.975)*sqrt(rrvar)
rrlower <- exp(logrrlower)
logrrupper <- log(rr)+qnorm(0.975)*sqrt(rrvar)
rrupper <- exp(logrrupper)
rr
rrlower
rrupper

# To get ratio of survival probabilities
# Relative risk
rr <- prob2/prob1
# Standard errors for probit linear predictors
xb1var <- as.numeric(Xsurv1%*%Valpha%*%Xsurv1)
xb2var <- as.numeric(Xsurv2%*%Valpha%*%Xsurv2)
# Covariance between xb1 and xb2
xb1xb2cov <- as.numeric(Xsurv1%*%Valpha%*%Xsurv2)
# Use delta method to get standard error for the log relative risk (note derivative of pnorm is dnorm)
xbV <- matrix(c(xb1var,xb1xb2cov,xb1xb2cov,xb2var),nrow=2)
gradrr1 <- -dnorm(sum(Xsurv1*alpha))/prob1
gradrr2 <- dnorm(sum(Xsurv2*alpha))/(prob2)
gradrr <- c(gradrr1,gradrr2)
rrvar <- gradrr%*%xbV%*%gradrr
# 95% Confidence interval
logrrlower <- log(rr)-qnorm(0.975)*sqrt(rrvar)
rrlower <- exp(logrrlower)
logrrupper <- log(rr)+qnorm(0.975)*sqrt(rrvar)
rrupper <- exp(logrrupper)

### Probability ratio for difference of -1 in FEV slope
fevslope2 <- fevslope - 1
Xsurv4 <- c(X,currentfev1,fevslope2)
prob4 <- pnorm(sum(Xsurv4*alpha))
rr <- (1-prob4)/(1-prob1)

# Standard errors for probit linear predictors
xb1var <- as.numeric(Xsurv1%*%Valpha%*%Xsurv1)
xb4var <- as.numeric(Xsurv4%*%Valpha%*%Xsurv4)
# Covariance between xb1 and xb4
xb1xb4cov <- as.numeric(Xsurv1%*%Valpha%*%Xsurv4)
# Use delta method to get standard error for the log relative risk (note derivative of pnorm is dnorm)
xbV <- matrix(c(xb1var,xb1xb4cov,xb1xb4cov,xb4var),nrow=2)
gradrr1 <- dnorm(sum(Xsurv1*alpha))/(1-prob1)
gradrr2 <- -dnorm(sum(Xsurv4*alpha))/(1-prob4)
gradrr <- c(gradrr1,gradrr2)
rrvar <- gradrr%*%xbV%*%gradrr
# 95% Confidence interval
logrrlower <- log(rr)-qnorm(0.975)*sqrt(rrvar)
rrlower <- exp(logrrlower)
logrrupper <- log(rr)+qnorm(0.975)*sqrt(rrvar)
rrupper <- exp(logrrupper)
rr
rrlower
rrupper

###To get ratio of survival probabilities
### Probability ratio for difference of -1 in FEV slope
fevslope2 <- fevslope - 1
Xsurv4 <- c(X,currentfev1,fevslope2)
prob4 <- pnorm(sum(Xsurv4*alpha))
rr <- (prob4)/(prob1)

# Standard errors for probit linear predictors
xb1var <- as.numeric(Xsurv1%*%Valpha%*%Xsurv1)
xb4var <- as.numeric(Xsurv4%*%Valpha%*%Xsurv4)
# Covariance between xb1 and xb4
xb1xb4cov <- as.numeric(Xsurv1%*%Valpha%*%Xsurv4)
# Use delta method to get standard error for the log relative risk (note derivative of pnorm is dnorm)
xbV <- matrix(c(xb1var,xb1xb4cov,xb1xb4cov,xb4var),nrow=2)
gradrr1 <- -dnorm(sum(Xsurv1*alpha))/(prob1)
gradrr2 <- dnorm(sum(Xsurv4*alpha))/(prob4)
gradrr <- c(gradrr1,gradrr2)
rrvar <- gradrr%*%xbV%*%gradrr
# 95% Confidence interval
logrrlower <- log(rr)-qnorm(0.975)*sqrt(rrvar)
rrlower <- exp(logrrlower)
logrrupper <- log(rr)+qnorm(0.975)*sqrt(rrvar)
rrupper <- exp(logrrupper)
rr
rrlower
rrupper


############################################################################################################
######          FIGURE 2: INTERPRETING JOINT MODEL RESULTS: DIFFERENT VALUES OF RANDOM EFFECTS         #####

rm(list=ls())
par(mfrow=c(1,2))
load(file="cffit.Rdata",sep=""))
est <- cffit.epinet$estimate

RE1 <- c(20,1,0)    # No quadratic age RE in the model, so set to 0 for everyone
RE2 <- c(0,0,0)
RE3 <- c(-20,-1,0)

# Fixed effects (intercept, slope and quadrtic term) for a 20-year old man (age in data is age - 5) born in 1980 (mean birth year=1986, birthyear_c=-6), two F508 alleles
byr <- -6
f508 <- 0     
sex <- 1
X1.int <- c(1,byr,sex,f508,byr*sex,byr*f508)
X1.slope <- c(1,sex,f508)
X1.quad <- 1
beta.int <- est[c(1,4:6,9:10)]
beta.slope <- est[c(2,7,8)]
beta.quad <- est[3]
intslope.fixed <- c(sum(X1.int*beta.int),sum(X1.slope*beta.slope),X1.quad*beta.quad)
# Add in the random effects
intslope1 <- intslope.fixed+RE1
intslope2 <- intslope.fixed+RE2
intslope3 <- intslope.fixed+RE3

f1 <- function(age){ intslope1[1] + intslope1[2]*(age-5) + intslope1[3]*(age-5)^2}
f2 <- function(age){ intslope2[1] + intslope2[2]*(age-5) + intslope2[3]*(age-5)^2}
f3 <- function(age){ intslope3[1] + intslope3[2]*(age-5) + intslope3[3]*(age-5)^2}

curve(f1, xlim=c(20,30), ylim=c(0,120), xlab="Age (yrs)", ylab="%FEV1", bty="l", lwd=2, lty=1, main="Longitudinal trajectories")
curve(f2, xlim=c(20,30), ylim=c(0,120), xlab="Age (yrs)", ylab="%FEV1", bty="l", lwd=2, add=T,lty=2)
curve(f3, xlim=c(20,30), ylim=c(0,120), xlab="Age (yrs)", ylab="%FEV1", bty="l", lwd=2, add=T, lty=3)
legend("topright", lwd=2, lty=1:3, legend=c("High intercept, shallow slope", "Pop. intercept, pop. slope", 
                                            "Low intercept, steep slope"), cex=0.7, bty="n")

prob1 <- 1
prob2 <- 1
prob3 <- 1
beta.surv <- est[c(11:15,20:21)]
beta.long <- est[1:10]
beta.slope <- est[c(2,7,8)]
X1.slope <- c(1,sex,f508)
for(i in 1:10){
  # Age at midpoint of current interval
  currentage <- (15+i-0.5)
  X1.long <- c(1,currentage,currentage^2,byr,sex,f508,currentage*sex,currentage*f508,byr*sex,byr*f508)
  # Current FEV values
  currentfev1 <- sum(X1.long*beta.long)+RE1[1]+RE1[2]*currentage
  currentfev2 <- sum(X1.long*beta.long)+RE2[1]+RE2[2]*currentage
  currentfev3 <- sum(X1.long*beta.long)+RE3[1]+RE3[2]*currentage
  currentslope1 <- sum(X1.slope*beta.slope) + 2*beta.quad*currentage + RE1[2]
  currentslope2 <- sum(X1.slope*beta.slope) + 2*beta.quad*currentage + RE2[2]
  currentslope3 <- sum(X1.slope*beta.slope) + 2*beta.quad*currentage + RE3[2]
  # probit linear predictors and survival probabilities
  Xsurv1 <- c(1,currentage,byr,sex,f508,currentfev1,currentslope1)
  Xsurv2 <- c(1,currentage,byr,sex,f508,currentfev2,currentslope2)
  Xsurv3 <- c(1,currentage,byr,sex,f508,currentfev3,currentslope3)
  mean1.i <- sum(Xsurv1*beta.surv)
  prob1.i <- pnorm(mean1.i)
  prob1 <- c(prob1,prob1[length(prob1)]*prob1.i)
  mean2.i <- sum(Xsurv2*beta.surv)
  prob2.i <- pnorm(mean2.i)
  prob2 <- c(prob2,prob2[length(prob2)]*prob2.i)
  mean3.i <- sum(Xsurv3*beta.surv)
  prob3.i <- pnorm(mean3.i)
  prob3 <- c(prob3,prob3[length(prob3)]*prob3.i)
}

x <- rep(20:30,c(1,rep(2,length(prob1)-1)))
y1 <- rep(prob1,c(rep(2,length(prob1)-1),1))
y2 <- rep(prob2,c(rep(2,length(prob2)-1),1))
y3 <- rep(prob3,c(rep(2,length(prob3)-1),1))

plot(x,y1, type="l", ylim=c(0,1), lty=1, main="Survival probabilities", xlab="Age (yrs)", ylab="Proportion alive", bty="l")
par(new=T)
plot(x,y2, type="l", ylim=c(0,1), lty=2, xlab="", ylab="", bty="l")
par(new=T)
plot(x,y3, type="l", ylim=c(0,1), lty=3, xlab="", ylab="", bty="l")
par(mfrow=c(1,1))


##################################################################################################
########      FIG 3: INTERPRETING JOINT MODEL RESULTS: COMPARING MALES AND FEMALES       ############

rm(list=ls())
load(file="cffit.Rdata")
est <- cffit.epinet$estimate

# Fixed effects (intercept and slope) for a 20-year old man (age in data is age - 5) born in 1980 (mean birth year=1986, birthyear_c=-6), 
# two F508 alleles
# Fixed effects (intercept, slope and quadrtic term) for a 20-year old man (age in data is age - 5) born in 1980 (mean birth year=1986, birthyear_c=-6), two F508 alleles
byr <- -6
f508 <- 0
sex1 <- 1
X1.int <- c(1,byr,sex1,f508,byr*sex1,byr*f508)
X1.slope <- c(1,sex1,f508)
X1.quad <- 1

# Fixed effects for a similar 20-year old woman
sex2 <- 0
X2.int <- c(1,byr,sex2,f508,byr*sex2,byr*f508)
X2.slope <- c(1,sex2,f508)
X2.quad <- 1

# Longitudinal linear predictors
beta.int <- est[c(1,4:6,9:10)]
beta.slope <- est[c(2,7,8)]
beta.quad <- est[3]
intslope1 <- c(sum(X1.int*beta.int),sum(X1.slope*beta.slope),X1.quad*beta.quad)
intslope2 <- c(sum(X2.int*beta.int),sum(X2.slope*beta.slope),X2.quad*beta.quad)

f1 <- function(age){ intslope1[1] + intslope1[2]*(age-5) + intslope1[3]*(age-5)^2}
f2 <- function(age){ intslope2[1] + intslope2[2]*(age-5) + intslope2[3]*(age-5)^2}

curve(f1, xlim=c(20,30), ylim=c(0,100), xlab="Age (yrs)", ylab="%FEV1", bty="l", lwd=2, lty=1, main="Longitudinal trajectories")
curve(f2, xlim=c(20,30), ylim=c(0,100), xlab="Age (yrs)", ylab="%FEV1", bty="l", lwd=2, lty=2, add=T)
legend("topright", lty=1:2, lwd=2, legend=c("Males", "Females"), cex=0.7, bty="n")

prob1 <- 1
prob2 <- 1
beta.surv <- est[c(11:15,20:21)]
beta.long <- est[1:10]
beta.slope <- est[c(2,7,8)]
X1.slope <- c(1,sex1,f508)
X2.slope <- c(1,sex2,f508)
for(i in 1:10){
# Age at midpoint of current interval
	currentage <- (15+i-0.5)
	X1.long <- c(1,currentage,currentage^2,byr,sex1,f508,currentage*sex1,currentage*f508,byr*sex1,byr*f508)
	X2.long <- c(1,currentage,currentage^2,byr,sex2,f508,currentage*sex2,currentage*f508,byr*sex2,byr*f508)
	
# Current FEV and slope values
	currentfev1 <- sum(X1.long*beta.long)
	currentfev2 <- sum(X2.long*beta.long)
	currentslope1 <- sum(X1.slope*beta.slope) + 2*beta.quad*currentage 
	currentslope2 <- sum(X2.slope*beta.slope) + 2*beta.quad*currentage
	
# probit linear predictors and survival probabilities
	Xsurv1 <- c(1,currentage,byr,sex1,f508,currentfev1,currentslope1)
	Xsurv2 <- c(1,currentage,byr,sex2,f508,currentfev2,currentslope2)
	mean1.i <- sum(Xsurv1*beta.surv) 
	prob1.i <- pnorm(mean1.i)
 	prob1 <- c(prob1,prob1[length(prob1)]*prob1.i)
	mean2.i <- sum(Xsurv2*beta.surv)
	prob2.i <- pnorm(mean2.i)
 	prob2 <- c(prob2,prob2[length(prob2)]*prob2.i)
}

x <- rep(20:30,c(1,rep(2,length(prob1)-1)))
y1 <- rep(prob1,c(rep(2,length(prob1)-1),1))
y2 <- rep(prob2,c(rep(2,length(prob2)-1),1))

plot(x,y1, type="l", ylim=c(0,1), lty=1, main="Survival probabilities", xlab="Age (yrs)", ylab="Proportion alive", bty="l")
par(new=T)
plot(x,y2, type="l", ylim=c(0,1), lty=2, xlab="", ylab="", bty="l")


#################################################################################################################
######          INTERPRETING JOINT MODEL RESULTS: DEMONSTRATE EFFECT OF VARIATION IN RANDOM EFFECTS         #####

library(MASS)

# Load results of model fit
rm(list=ls())
load(file="cffit.Rdata")
est <- cffit.epinet$estimate

frhoinv <- function(x){
  x <- exp(x)/(1+exp(x))
  2*x-1
}

# Function to plot longitudinal trajectories and survival curves
heterogen_plot <- function(){
  par(mfrow=c(1,2))
  
  # Draw N samples from random effects distribution
  N <- 50
  sigma0 <- exp(est[names(est)=="sigma0"])
  sigma1 <- exp(est[names(est)=="sigma1"])
  rho <- frhoinv(est[names(est)=="rho"])
  REcov <- sigma0*sigma1*rho
  Sigma <- matrix(c(sigma0^2,REcov,REcov,sigma1^2),2,2)
  set.seed(10896)
  RE <- mvrnorm(n=N, rep(0, 2), Sigma)
  RE <- cbind(RE,rep(0,nrow(RE)))               # Add a 0 for quadratic term which is treated as fixed effect only
  
  # Fixed effects (intercept and slope) for a 20-year old man (age in data is age - 5) born in 1980 (mean birth year=1986, birthyear_c=-6), 
  # two F508 alleles
  # Fixed effects (intercept, slope and quadrtic term) for a 20-year old man (age in data is age - 5) born in 1980 (mean birth year=1986, birthyear_c=-6), two F508 alleles
  byr <- -6
  f508 <- 0     
  sex <- 1
  X1.int <- c(1,byr,sex,f508,byr*sex,byr*f508)
  X1.slope <- c(1,sex,f508)
  X1.quad <- 1
  beta.int <- est[c(1,4:6,9:10)]
  beta.slope <- est[c(2,7,8)]
  beta.quad <- est[3]
  intslope.fixed <- c(sum(X1.int*beta.int),sum(X1.slope*beta.slope),X1.quad*beta.quad)
  
  # Plot individual trajectories
  for(j in 1:N){
    intslope <- intslope.fixed+RE[j,]
    fRE <- function(age){ intslope[1] + intslope[2]*(age-5) + intslope[3]*(age-5)^2}
    curve(fRE, xlim=c(20,30), ylim=c(0,120), bty="l", lwd=1, col="grey", xlab="", ylab="")
    par(new=T)
  }
  
  # Plot population trajectory
  f1 <- function(age){ intslope.fixed[1] + intslope.fixed[2]*(age-5) + intslope.fixed[3]*(age-5)^2}
  curve(f1, xlim=c(20,30), ylim=c(0,120), xlab="Age (yrs)", ylab="%FEV1", bty="l", lwd=2, main="Longitudinal trajectories")
  
  # Plot individual survival curves
  beta.surv <- est[c(11:15,20:21)]
  beta.long <- est[1:10]
  X.slope <- c(1,sex,f508)
  for(j in 1:N){
    prob <- 1
    for(i in 1:10){
      # Age at midpoint of current interval
      currentage <- (15+i-0.5)
      X.long <- c(1,currentage,currentage^2,byr,sex,f508,currentage*sex,currentage*f508,byr*sex,byr*f508)
      # Current FEV values
      currentfev <- sum(X.long*beta.long)+RE[j,1]+RE[j,2]*currentage
      currentslope <- sum(X.slope*beta.slope) + 2*beta.quad*currentage + RE[j,2]
      # probit linear predictors and survival probabilities
      Xsurv <- c(1,currentage,byr,sex,f508,currentfev,currentslope)
      
      mean.i <- sum(Xsurv*beta.surv) 
      prob.i <- pnorm(mean.i)
      prob <- c(prob,prob[length(prob)]*prob.i)
    }
    
    x <- rep(20:30,c(1,rep(2,length(prob)-1)))
    y <- rep(prob,c(rep(2,length(prob)-1),1))
    
    plot(x,y, type="l", ylim=c(0,1), main="Survival probabilities", xlab="Age (yrs)", ylab="Proportion alive", bty="l", col="grey")
    par(new=T)
  }
  
  # Plot population survival
  prob <- 1
  for(i in 1:10){
    # Age at midpoint of current interval
    currentage <- (15+i-0.5)
    X.long <- c(1,currentage,currentage^2,byr,sex,f508,currentage*sex,currentage*f508,byr*sex,byr*f508)
    # Current FEV values
    currentfev <- sum(X.long*beta.long)
    currentslope <- sum(X.slope*beta.slope) + 2*beta.quad*currentage
    # probit linear predictors and survival probabilities
    Xsurv <- c(1,currentage,byr,sex,f508,currentfev,currentslope)
    mean.i <- sum(Xsurv*beta.surv) 
    prob.i <- pnorm(mean.i)
    prob <- c(prob,prob[length(prob)]*prob.i)
  }
  
  x <- rep(20:30,c(1,rep(2,length(prob)-1)))
  y <- rep(prob,c(rep(2,length(prob)-1),1))
  
  plot(x,y, type="l", ylim=c(0,1), main="Survival probabilities", xlab="Age (yrs)", ylab="Proportion alive", bty="l", lwd=2)
  
  par(mfrow=c(1,1))
}


######################################################################################
######################## LONGITUDINAL MODEL DIAGNOSTIC PLOTS #########################

rm(list=ls())

# Load the data
source("CF_JM_functions.R")
temp <- getdat(whiteonly=T, rs=F, fortyplus=F, datafile="JM_dataset3.RData")
cfdat <- temp$data
cfsurv <- temp$survdata
rm(temp)

# Remove 13 individuals with time=0 (these are individuals whose death or censoring date occures before their first visit date)
ids.to.remove <- cfsurv$id[cfsurv$time==0]
cfsurv <- cfsurv[!(cfsurv$id %in% ids.to.remove),]
cfdat <- cfdat[!(cfdat$id %in% ids.to.remove),]
cfsurv$id <- factor(as.character(cfsurv$id))
cfdat$id <- factor(as.character(cfdat$id))

# Set up the data and formulae for model fit
dat.y <- with(cfdat,data.frame(id,fev,ageatvisit,birthyear_c,sex,F508.0.1.NT))
dat.y$ageatvisitsq <- dat.y$ageatvisit^2
fmla.y <- (fev ~ ageatvisit + ageatvisitsq + birthyear_c + sex + F508.0.1.NT + sex:ageatvisit + sex:birthyear_c + F508.0.1.NT:ageatvisit + F508.0.1.NT:birthyear_c)
longfit1 <- lme(fmla.y, random=~ageatvisit|id, dat=dat.y, control=list(opt="optim"))
plot(longfit1)

## Try it with log-transformed FEV
dat.y$logfev <- log(dat.y$fev)
fmla.y2 <- (logfev ~ ageatvisit + ageatvisitsq + birthyear_c + sex + F508.0.1.NT + sex:ageatvisit + sex:birthyear_c + F508.0.1.NT:ageatvisit + F508.0.1.NT:birthyear_c)
longfit2 <- lme(fmla.y2, random=~ageatvisit|id, dat=dat.y)
plot(longfit2)


################################################################################################
##########################    SURVIVAL MODEL DIAGNOSTIC PLOTS    ###############################

rm(list=ls())
load("CF_survresults.Rdata")

library(DHARMa)

set.seed(1411)
res = simulateResiduals(survfit3)
plot(res)


######################################################################################################
##############################     SPAGHETTI PLOT WITH MEAN SMOOTHER     #############################

rm(list=ls())

# Load the data
source("CF_JM_functions.R")
temp <- getdat(whiteonly=T, rs=F, fortyplus=F, datafile="JM_dataset3.RData")
cfdat <- temp$data
cfsurv <- temp$survdata
rm(temp)

# Remove 13 individuals with time=0 (these are individuals whose death or censoring date occures before their first visit date)
ids.to.remove <- cfsurv$id[cfsurv$time==0]
cfsurv <- cfsurv[!(cfsurv$id %in% ids.to.remove),]
cfdat <- cfdat[!(cfdat$id %in% ids.to.remove),]
cfsurv$id <- factor(as.character(cfsurv$id))
cfdat$id <- factor(as.character(cfdat$id))

cfdat$trueage <- cfdat$ageatvisit+5
cfdat <- as.data.frame(cfdat)

library(ggplot2)
p <- ggplot(data = as.data.frame(cfdat), aes(x = trueage, y = fev, group = id))
p + geom_point() + stat_smooth(aes(group = 1))
p + geom_line()
p + geom_line() + stat_smooth(aes(group = 1)) + labs(x="Age (yrs)", y = "%FEV1")

set.seed(154873)
sam <- sample(cfsurv$id, size=200)
cfdat.sam <- cfdat[cfdat$id %in% sam,]
cfdat.sam$id <- factor(as.character(cfdat.sam$id))

psam <- ggplot(data = as.data.frame(cfdat.sam), aes(x = trueage, y = fev, group = id))
psam + geom_line() + stat_smooth(aes(group = 1)) + labs(x="Age (yrs)", y = "%FEV1")


