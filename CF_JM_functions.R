#########################################################################################################
#################################### GENERATES THE CF DATAFRAMES ########################################


getdat <- function(whiteonly=T, rs=F, fortyplus=F, datafile){
# Reads in the CF data.
# whiteonly=T to include only white ethnicities
# fortyplus=F to censor at age 40
# rs=T to select a randomsample, 20% of patients  
# datafile is the path to the .Rdata file containing the dataset 
	temp <- setdat(datafile=datafile,whiteonly=whiteonly,fortyplus=fortyplus)
	cfdat <- temp$data
	cfsurv <- temp$surv
	rm(temp)

# Take a random sample of the data
	if(rs){
		prop <- 0.2
		set.seed(154873)
		sam <- sample(cfsurv$id, size=round(0.2*nrow(cfsurv)))

		cfdat <- cfdat[cfdat$id %in% sam,]
		cfsurv <- cfsurv[cfsurv$id %in% sam,]
		cfdat$id <- factor(as.character(cfdat$id))
		cfsurv$id <- factor(as.character(cfsurv$id))
	}

	out <- list(data=cfdat, survdata=cfsurv)
}


setdat <- function(datafile="",whiteonly,fortyplus){
# Creates longitudinal and survival dataframes.

# Load the data with the following variables (in long format)
# id = patient id
# fev = FEV_1
# sex
# nonwhite = 0 if white, 1 otherwise
# imd = index of multiple deprivation score
# birthyear
# ageatvisit = age at this row's visit
# F508_class = class of F508 alleles (!!! check classes)
# time_sincebaseline_next is time of next visit, censored by death or 31st December 2015
# d=1 if death occurred in this row's time interval
	load(datafile)
  
# Check number of observations over age 40
  percentover40 <- sum(dataset5$ageatvisit>40)*100/nrow(dataset5)
    
# Variable for homozygosity in F508 alleles
  dataset5$F508.0.1.NT <- as.numeric(dataset5$F508_class!="Homoz")
  
# order data by id and visitage
  dataset5 <- dataset5[order(dataset5$id,dataset5$ageatvisit),] 

# Create a variable for age at entry into dataset    
  startagebyid <- tapply(dataset5$ageatvisit,dataset5$id,min)
  dataset5$startage <- rep(startagebyid,table(dataset5$id))
  
# Birth cohort
  dataset5$birthcohort <- cut(dataset5$birthyear,breaks=c(1920,1959.5,1964.5,1969.5,1974.5,1979.4,1984.5,1989.5,1994.5,2020),
                     labels=c("<1960","1960-1964","1965-1969","1970-1974","1975-1979","1980-1984","1985-1989","1990-1994",">1995"))
  
# Create a variable for death/censoring time
  timebyid <- tapply(dataset5$time_sincebaseline_next,dataset5$id,max)
  dataset5$time <- rep(timebyid,table(dataset5$id))
  
# Create an indicator variable for death/censoring
  diedbyid <- tapply(dataset5$d,dataset5$id,max)
  dataset5$died <- rep(diedbyid,table(dataset5$id))
  
	FEV <- dataset5[, c("id", "fev", "ageatvisit", "sex", "F508.0.1.NT", "imd", "birthyear", "nonwhite")]
	SUR <- dataset5[, c("id", "time", "died", "startage", "sex", "F508.0.1.NT", "F508_class", "imd", "birthyear", "birthcohort","nonwhite")]
	SUR <- SUR[!duplicated(SUR),]
	SUR <- SUR[match(unique(SUR$id),SUR$id),]
	
# if whiteonly=T restrict to white population
	if(whiteonly){
		FEV <- FEV[FEV$nonwhite==0,]
		SUR <- SUR[SUR$nonwhite==0,]
	}

# if fortyplus=F then truncate longitudinal data at age 40 and remove ids who have inital age > 40
	if(!fortyplus){
		FEV <- FEV[FEV$ageatvisit<40,]
		SUR <- SUR[SUR$startage<40,]		
	}

# factor IDs
	SUR$id <- as.character(SUR$id)
	FEV$id <- as.character(FEV$id)
	idall <- unique(c(SUR$id,FEV$id))
	SUR$id <- factor(SUR$id,levels=idall)
	FEV$id <- factor(FEV$id,levels=idall)

# Set up variables for longitudinal dataset
# "ageatvisit" is the age at visit minus 5 years, because lung function can't be measured before 5 years
	FEV$ageatvisit <- FEV$ageatvisit-5
# birthyear will be centred below using the mean from the survival dataset
	FEV$birthyear_c <- FEV$birthyear	
	fevdata <- FEV[,c("id", "ageatvisit", "birthyear_c", "fev", "sex", "F508.0.1.NT", "imd")]

# Set up variables for survival dataset
	SUR$status <- 1-SUR$died     # Status is 0=died, 1=censored
# If fortyplus=F then censor survival data at age 40
	if(!fortyplus){
	  deathageyrs <- ifelse(SUR$died==1, SUR$time+SUR$startage, 0)
		SUR$status[deathageyrs > 40] <- 1
		SUR$time[deathageyrs > 40] <- (40-SUR$startage[deathageyrs>40])
	}
	SUR$exacttime <- SUR$time
	SUR$time <- as.integer(ceiling(SUR$time))
	mean.byr <- mean(SUR$birthyear)
	SUR$birthyear_c <- SUR$birthyear - mean.byr
	surdata <- SUR[,c("id", "status", "time", "exacttime","startage", "birthyear_c", "birthcohort", "sex", "F508.0.1.NT", "F508_class", "imd")]
	
# centre longitudinal birth year
	fevdata$birthyear_c <- fevdata$birthyear_c - mean.byr
	cat("Mean birth year is ",mean.byr,"\n")
	
	out <- list(data=fevdata, survdata=surdata)
}



##################################################################################################
#################################### FITS THE MODEL ##############################################

fitmodel <- function(Dat, longformula, Survdat, inits=NA, intervaltimes=((0:15)+0.5), initses=rep(1,length(inits)), printpars=F, maxiter=10000, pmnormtol=1e-3,
		method.lme=c("optim","nlminb")){
# fits the random intercept/slope model
# Dat is longitudinal dataframe
# longformula is formula to specify longitudinal model fixed effects
# Survdat is survival dataframe
# use inits to specify initial parameter values for optimisation routine
# intervaltimes specifies the mid-points of the time intervals
# initses is initial estimate of standard errors, used to scale parameters for improved performance of opitmisatin routine
# printpars and maxiter are arguments passed to optimisation routine 
# pmnormtol is tolerance in calculation of normal CDFs
# method.lme specifies optimisation for method in lme when generating initial parameters, either "optim" or "nlminb"

# Check all survival times are > 0
	if(any(Survdat$time <= 0)){
		print("Error: survival data contains times <= 0")
		return(NULL)
	}

# if inits not specified fit separate longitudinal and survival models to get initial parameter values
	if(any(is.na(inits))){
# initial parameter values for longitudinal data
		if(method.lme == "optim"){								
# use do.call, otherwise predict below can't find parsed formula
			longfit <- do.call(lme,list(fixed=longformula,random=~ageatvisit|id,data=Dat,control=list(opt="optim")))
		}
		if(method.lme == "nlminb"){
			longfit <- do.call(lme,list(fixed=longformula,random=~ageatvisit|id,data=Dat))
		}
		meanpars.y <- fixef(longfit)
		varpars <- c(longfit$sigma,intervals(longfit)$reStruct$id[,2])
		se.y <- sqrt(diag(longfit$varFix))
		se.var <- sqrt(diag(longfit$apVar))[c(4,1:3)] # N.B. apVar already gives transformed standard errors
		varpars <- c(log(varpars[1:3]),frho(varpars[4]))
		
# initial parameter values for survival data
# survival data in long format
		ind  <- rep(1:nrow(Survdat),Survdat$time)
		Survdatnew <- Survdat[ind,]
# generate a time interval variable
		temp <- rep(1,nrow(Survdatnew))
		Survdatnew$interval <- unlist(tapply(temp,Survdatnew$id,cumsum))
		Survdatnew$intervalage <- intervaltimes[Survdatnew$interval]+Survdatnew$startage
# To include current fev value and slope in survival model, get y covariate matrix in long survival format (i.e. one row per time interval)
		ind.y <- match(Survdatnew$id,Dat$id)
		Datnew <- Dat[ind.y,-which(names(Dat) %in% c("ageatvisit","ageatvisitsq"))]
		Datnew$ageatvisit <- Survdatnew$intervalage
		Datnew$ageatvisitsq <- Survdatnew$intervalage^2
		Xyt <- model.matrix(longformula, data=Datnew) 	
		Survdatnew$currentfev <- as.vector(Xyt %*% meanpars.y) + ranef(longfit)[ind,1] + ranef(longfit)[ind,2]*Datnew$ageatvisit
		slopeind <- c(match("ageatvisit",names(meanpars.y)),match("ageatvisitsq",names(meanpars.y)),grep("ageatvisit:",names(meanpars.y)),grep(":ageatvisit",names(meanpars.y)))
		Xslope <- Xyt[,slopeind]
		Xslope <- Xslope/Xslope[,1]     # Divide by age at visit
		Xslope[,which(colnames(Xslope)=="ageatvisitsq")] <- 2*Xslope[,which(colnames(Xslope)=="ageatvisitsq")]
		Survdatnew$fevslope <- Xslope %*% meanpars.y[slopeind] + ranef(longfit)[ind,2]  
# generate the outcome for probit analysis of survival data
		Survdatnew$Yes <- 1
		Survdatnew$Yes[Survdatnew$status==0 & Survdatnew$interval==Survdatnew$time] <- 0
# formula for survival analysis
		xnamsurv <- names(Survdat)[-(1:4)]
		xnamsurv <- c(xnamsurv,"currentfev","fevslope")
		fmla.surv <- as.formula(paste("Yes ~ intervalage + ", paste(xnamsurv, collapse= "+")))
		survfit <- glm(fmla.surv, family=binomial(link="probit"), data=Survdatnew)
		pt <- length(coef(survfit))-2
		meanpars.t <- coef(survfit)[1:pt]
		gampars <- coef(survfit)[(pt+1):(pt+2)]
		se.t <- sqrt(diag(vcov(survfit)))[1:pt]
		se.gam <- sqrt(diag(vcov(survfit)))[(pt+1):(pt+2)]

		initsnew <- c(meanpars.y,meanpars.t,varpars,gampars)
		initsesnew <- c(se.y,se.t,se.var,se.gam)
# To allow for some initial values to be specified via inits
		if(any(!is.na(inits))){
			ind.na <- which(is.na(inits))
			inits[ind.na] <- initsnew[ind.na]
 			initses[ind.na] <- initsesnew[ind.na]
		} else {
			inits <- initsnew
			initses <- initsesnew
		}
	}

# Scale initial parameter values by standard errors
	inits <- inits/initses
		
# fit the model with Hessian
	fit <- optim(par=inits, fn=loglik, dat=Dat, fmla=as.formula(longformula), surv=Survdat, intervaltimes=intervaltimes, printpars=printpars, ses=initses, 
				pmnormtol=pmnormtol, hessian=T, method="BFGS", control=list(maxit=maxiter, reltol=1e-9))	
	fithess <- fit$hessian

# Use the Hessian to estimate standard errors
	estimate <- fit$par*initses
	hessinv <- try(solve(fithess))
	if(class(hessinv)=="try-error"){
		fitses <- NULL
	} else {
		fitses <- sqrt(diag(hessinv))*initses
	}

# Name estimates and standard errors from optimisation output
	longnames <- c("(Intercept)",attr(terms(longformula),"term.labels"))
	survnames <- c("(Intercept)","Age",names(Survdat)[-(1:4)])
	names(estimate) <- c(longnames,survnames,"sigmay","sigma0","sigma1","rho","gamma0","gamma1")
	if(class(hessinv)!="try-error"){
		names(fitses) <- names(estimate)
	}	

	out <- list(estimate=estimate, se=fitses, fit=fit, fithess=fithess, longformula=longformula)
	out
}







##############################################################################################################
############################### FUNCTIONS TO USE IN LOGLIK AND GETDAT ########################################


 
Sigma.is <- function(upars){
# generates random effects covariance matrix for intercept/slope model
	sig0 <- upars[1]
	sig1 <- upars[2]
	rho <- upars[3]
    	mat <- matrix(c(sig0^2,rho*sig0*sig1,rho*sig0*sig1,sig1^2),nrow=2,ncol=2)
    	mat
}


L.cv <- function(p=nt,gampars){
# specifies how intercept/slope random effects enter into survival model 
# current value model
  gam0 <- gampars[1]
  gam1 <- gampars[2]
  mat <- matrix(0,nrow=p,ncol=2)
  mat[,1] <- gam0
  mat[,2] <- gam1
  mat[,2] <- mat[,2]+(0:(p-1)+0.5)*gam0
  mat
}

frho <- function(p){
# transform correlation range (-1,1) to entire real line  
	p <- (p+1)/2 
	log(p/(1-p))
}

frhodiff <- function(p){
# differential of frho
	1/(1-p^2)
}


################################################################################################################
#################################### CALCULATES THE LOGLIKELIHOOD ##############################################


loglik <- function(pars,  ses=NULL, dat, fmla, surv, intervaltimes=((0:15)+0.5), pmnormtol=1e-3, printpars=F, timeinteraction=F){
# likelihood for random intercept/slope model
# parameters
	if(!is.null(ses))	pars <- pars*ses
  	if(printpars) print(pars)
  	n <- dim(surv)[1]
  	ny <- dim(dat)[1]
  	nt <- max(surv$time)
  	
	tf <- terms(fmla)
	py <- length(attr(tf,"term.labels"))+attr(tf,"intercept")
	pt <- ncol(surv)-2
	pyt <- py+pt
	by <- pars[1:py]
	bt <- pars[(py+1):pyt]         
	sig.y <- exp(pars[(pyt+1)])

# calculate random effects covariance Sig.u and L, which gives random effects terms in survival model   	
	rho <- 2*exp(pars[pyt+4])/(1+exp(pars[pyt+4]))-1
  	upars <- c(exp(pars[(pyt+2):(pyt+3)]),rho)
  	gampars <- pars[(pyt+5):(pyt+6)]
  	Sig.u <- Sigma.is(upars)
  	L <- L.cv(p=nt, gampars)
	nu <- 2
	  
# longitudinal data
  	id <- dat$id
  	y <- dat$fev
	time.y <- dat$ageatvisit
  	X.y <- model.matrix(fmla, data=dat)

# survival data
	status <- surv$status
	time <- surv$time
# transform survival data to long format (one row per time interval)
	intervaltime <- rep(intervaltimes, n)
	ind.t <- rep(1:n,rep(nt,n))
	X.t <- surv[,-(1:3)]
	X.t <- cbind(1,X.t[ind.t,])
	X.t$startage <- X.t$startage+intervaltime
	names(X.t)[2] <- "intervalage"
# To include current fev value and slope in survival model, get y covariate matrix in long survival format (i.e. one row per time interval)
	ind.y <- match(surv$id[ind.t],dat$id)
	datnew <- dat[ind.y,-which(names(dat) %in% c("ageatvisit","ageatvisitsq"))]
	datnew$ageatvisit <- X.t$intervalage
	datnew$ageatvisitsq <- X.t$intervalage^2
	X.y.t <- model.matrix(fmla, data=datnew)
	names(by) <- colnames(X.y.t) 	
	X.t$currentfev <- X.y.t %*% by
	slopeind <- c(match("ageatvisit",names(by)),match("ageatvisitsq",names(by)),grep("ageatvisit:",names(by)),grep(":ageatvisit",names(by)))
	X.slope <- X.y.t[,slopeind]
# remove time from slope covariate
	X.slope <- X.slope/X.slope[,1]
	X.slope[,which(colnames(X.slope)=="ageatvisitsq")] <- 2*X.slope[,which(colnames(X.slope)=="ageatvisitsq")]
	X.t$fevslope <- X.slope %*% by[slopeind]  
	

# calculate H - my H is different from the one in the paper, it's A^T A / \nu^2  
# H is a square matrix which is different for each person, so put in a 3-dimensional array
  	H <- array(0,dim=c(2,2,n))
  	H[1,1,] <- table(factor(id,levels=surv$id)) 
  	H[1,2,] <- tapply(time.y,factor(id,levels=surv$id),sum)
  	H[2,1,] <- H[1,2,]
  	H[2,2,] <- tapply(time.y^2,factor(id,levels=surv$id),sum)
  	H <- ifelse(is.na(H),0,H)
  	H <- H / sig.y^2
	
# Calculate V - in the paper notation it's H^{-1} = (A^T A / \nu^2 + \Sigma^{-1})^{-1} (again different for each person)  
  	Sig.u.inv <- solve(Sig.u)
  	getV <- function(M){solve(M+Sig.u.inv)}
  	V <- aaply(H,3,getV)
  	V <- aperm(V,perm=c(2,3,1))

# remove fixed effects from y
	yran <- y - as.matrix(X.y) %*% as.matrix(by)

# calculate h - in the paper notation it's A^T (y - X \beta) / \nu^2
   	h <- matrix(nrow=n,ncol=2)
 	h[,1] <- tapply(yran,factor(id,levels=surv$id),sum) /sig.y^2
  	h[,2] <- tapply(yran*time.y,factor(id,levels=surv$id),sum) /sig.y^2 
  	h <- ifelse(is.na(h),0,h)
	
# calculate q.Phi, arguments of multivariate normal CDFs (each person has one for each time-point) 
  	Vh <- matrix(nrow=nu,ncol=n)
  	for(j in 1:nu){
    		tmp <- V[j,,] * t(h)
    		Vh[j,] <- apply(tmp,2,sum)      
  	}

	betaX <- as.matrix(X.t) %*% as.matrix(c(bt,gampars))
	betaX <- matrix(betaX, nrow=nt, ncol=n)
	q.Phi <- betaX + L%*%Vh

# calculate the multivariate normal CDFs and put results in Phi, Phi[,i] is the probability of surviving to the start of the ith interval
	Phi <- matrix(0,nrow=n,ncol=(nt+1))
# First column of Phi is the probability of surviving to time 0
  	Phi[,1] <- 1
	I <- array(diag(nt),dim=c(nt,nt,n))
  	for(i in 1:nt){
  		Inew <- I[1:i,1:i,]
  	
   		Lnew <- L[1:i,]
  		if(i==1){
  			tmp <- tensor(as.matrix(Lnew),V,1,1)
  			LVL <- tensor(tmp,as.matrix(Lnew),2,1)
  		} else {
  			tmp <- tensor(Lnew,V,2,1)
  			LVL <- tensor(Lnew,tmp,2,2)
  		}
     		Sig.Phi <- Inew + LVL
    		mu.Phi <- rep(0,i)
 
# only calculate the components of Phi that are actually going to be used in the likelihood 
# (i.e. time for censored individuals and time and (time - 1) for deaths)   
    		ids1 <- c(which((time==i)&status==0),which((time==(i+1))&status==0))
    		ids2 <- which((time==i)&status==1)
    		ids <- c(ids1,ids2)
   
    		p.Phi <- rep(0,n)
    		if(i==1){
      		p.Phi[ids] <- pnorm(q.Phi[1,ids],mean=0,sd=sqrt(Sig.Phi[ids]))
    		} else {
      		for(k in ids){
        			p.Phi[k] <- pmnorm(q.Phi[1:i,k],mu.Phi,Sig.Phi[,,k],abseps=pmnormtol)
      		}
    		}
    		Phi[,(i+1)] <- p.Phi
  	}

# Psi[,i] is the probability of dying in time interval i 
	Psi <- matrix(nrow=n,ncol=nt)
  	for(i in 1:nt){
    		Psi[,i] <- Phi[,i] - Phi[,(i+1)]
  	}

# calculate all the terms in the likelihood
  	term1 <- 0
  	for(i in 1:nt){
    		ids <- which((time==i)&status==0)
    		term1 <- term1 + sum(log(Psi[ids,i]))
  	}
  	term2 <- 0
  	for(i in 1:nt){
    		ids <- which((time==i)&status==1)
    		term2 <- term2 + sum(log(Phi[ids,(i+1)]))
  	}
	term3 <- - sum(yran^2)/(2*sig.y^2)  
    	tmp <- matrix(nrow=nu,ncol=nu)
  	for(i in 1:nu){
    		for(j in 1:nu){
      		tmp[i,j] <- sum(h[,i] * V[i,j,] * h[,j])
    		}	
  	}
  	term4 <- 0.5 * sum(tmp)
  	term5 <- -ny*log(sig.y)
 	term6 <- - 0.5*n*log(det(Sig.u))
   	term7 <- 0.5*sum(log(apply(V,3,det)))

  	loglik <- term1 + term2 + term3 + term4 + term5 + term6 + term7 
  	if (is.finite(loglik)) {
      	return(-loglik)
    	} else {
        	return(1e+09)
    	}

}



#################################################################################################################
####################################        PRINT RESULTS     ###################################################

frhoinv <- function(x){
	x <- exp(x)/(1+exp(x))
	2*x-1
}

frhoinvdiff <- function(x){
	2*exp(x)/(1+exp(x))^2
}

summary.modelfit <- function(modelfit, digits=3){

	if(modelfit$fit$convergence == 0){ cat("\n Optimisation routine converged.\n\n")
	} else { cat("\n Warning: optimisation routine didn't converge (", modelfit$fit$message,")\n\n") }

	est <- modelfit$estimate
	se <- modelfit$se
	if(any(is.na(se)))  se <- modelfit$se.new 
	if(all(is.na(se)))  se <- modelfit$se

	ind.sigy <- which(names(est)=="sigmay")
	ind.sig1 <- which(names(est)=="sigma1")
	ind.rho <- which(names(est)=="rho")
	ind.gam0 <- which(names(est)=="gamma0")
	ind.gam1 <- which(names(est)=="gamma1")
	
# 95% Confidence intervals
	lci <- est-qnorm(0.975)*se
	uci <- est+qnorm(0.975)*se
# P-values
	pval <- 2*(1-pnorm(abs(est)/se))

	estnew <- c(est[1:(ind.sigy-1)],exp(est[ind.sigy:ind.sig1]),frhoinv(est[ind.rho]),est[ind.gam0:ind.gam1])
	senew <- c(se[1:(ind.sigy-1)],exp(est[ind.sigy:ind.sig1])*se[ind.sigy:ind.sig1],frhoinvdiff(est[ind.rho])*se[ind.rho],
			se[ind.gam0:ind.gam1])
	lcinew <- c(lci[1:(ind.sigy-1)],exp(lci[ind.sigy:ind.sig1]),frhoinv(lci[ind.rho]),lci[ind.gam0:ind.gam1])
	ucinew <- c(uci[1:(ind.sigy-1)],exp(uci[ind.sigy:ind.sig1]),frhoinv(uci[ind.rho]),uci[ind.gam0:ind.gam1])
	loglikelihood <- -modelfit$fit$value
	AIC <- 2*length(modelfit$est) + 2*modelfit$fit$value

	longformula = modelfit$longformula
	tf <- terms(longformula)
	py <- length(attr(tf,"term.labels"))+attr(tf,"intercept")
	longtable <- matrix(ncol=5, nrow=py) 
	longtable[,1] <- estnew[1:py]
	longtable[,2] <- senew[1:py]
	longtable[,3] <- lcinew[1:py]
	longtable[,4] <- ucinew[1:py]
	longtable[,5] <- pval[1:py]
	colnames(longtable) <- c("Value","Std.Error","Lower 95% CI", "Upper 95% CI", "P-value")
	rownames(longtable) <- names(estnew)[1:py]

	pt <- ind.sigy - py - 1
	survtable <- matrix(ncol=5, nrow=pt)
	survtable[,1] <- estnew[(py+1):(ind.sigy-1)]
	survtable[,2] <- senew[(py+1):(ind.sigy-1)]
	survtable[,3] <- lcinew[(py+1):(ind.sigy-1)]
	survtable[,4] <- ucinew[(py+1):(ind.sigy-1)]
	survtable[,5] <- pval[(py+1):(ind.sigy-1)]
	colnames(survtable) <- c("Value","Std.Error","Lower 95% CI", "Upper 95% CI", "P-value")
	rownames(survtable) <- names(estnew)[(py+1):(ind.sigy-1)]

	gamtable <- matrix(ncol=5, nrow=2)
	gamtable[,1] <- estnew[ind.gam0:ind.gam1]
	gamtable[,2] <- senew[ind.gam0:ind.gam1]
	gamtable[,3] <- lcinew[ind.gam0:ind.gam1]
	gamtable[,4] <- ucinew[ind.gam0:ind.gam1]
	gamtable[,5] <- pval[ind.gam0:ind.gam1]
	colnames(gamtable) <- c("Value","Std.Error","Lower 95% CI", "Upper 95% CI","P-value")
	rownames(gamtable) <- c("Current value","Slope")

	retable <- matrix(ncol=5, nrow=4)
	retable[,1] <- estnew[c((ind.sigy+1):ind.rho,ind.sigy)]
	retable[,2] <- senew[c((ind.sigy+1):ind.rho,ind.sigy)]
	retable[,3] <- lcinew[c((ind.sigy+1):ind.rho,ind.sigy)]
	retable[,4] <- ucinew[c((ind.sigy+1):ind.rho,ind.sigy)]
	retable[,5] <- pval[c((ind.sigy+1):ind.rho,ind.sigy)]
	colnames(retable) <- c("Value","Std.Error","Lower 95% CI", "Upper 95% CI", "P-value")
	rownames(retable) <- c("Intercept SD", "Slope SD", "Intercept/slope corr.", "Residual SD")
	
	cat("Longitudinal model estimates \n")
	longtable <- format(round(longtable,digits), nsmall=digits)
	print(longtable, quote=F)
	cat("\nSurvival model estimates \n")
	survtable <- format(round(survtable,digits), nsmall=digits)
	print(survtable, quote=F)
	cat("Association parameters \n")
	gamtable <- format(round(gamtable,digits), nsmall=digits)
	print(gamtable, quote=F)
	
	cat("\nRandom effects parameters\n")
	retable <- format(round(retable,digits), nsmall=digits)
	print(retable, quote=F)

	cat("\nloglikelihood = ",loglikelihood)	
	cat("\nAIC = ",AIC,"\n")
	
	out=list(longtable=as.table(longtable), survtable=as.table(survtable), gamtable=as.table(gamtable), retable=as.table(retable))
}


