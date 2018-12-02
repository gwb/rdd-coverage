rm(list=ls())

# function gen.2RRDD
# last update: 4/2/2013 Kristin 
# summary: generates multiple samples from a 2 rating RDD design
#
# inputs:
# 	n = sample size
#     s = number of samples
# 	mu1 = mean of rating 1
#	mu2 = mean of rating 2
#	sigma.r1 = sd of true rating 1
#	sigma.r2 = sd of true rating 2
#	rho.r1r2 = correlation between true ratings
#	cut1.quantile = quantile for cutting rating 1
#	cut2.quantile = quantile for cutting rating 2
# 	cut1.value = value for cutting rating 1
#	cut2.value = value for cutting rating 2
#	param.Y0 = parameters for fn relating ratings and Y0 
#	param.eff.T2 = parameters for fn relating ratings and effect of failing both
#	param.eff.T1 = parameters for fn relating ratings and effect of failing r1
#	param.eff.T3 = parameters for fn relating ratings and effect of failing r2
#	param.eff.T = parameters for fn relating ratings and effect for any way of rec treatment
#		# for full model with terms r1,r2,r1*r2,r1^2,r2^2
#	sigma.noise<-sd of error term in surface response model
#	sigma.E<-sd of random error in effects (assumed same for all params)
#	data.type = "full","observed",or "none" (no data returned - just parameter values)
#
# outputs:
#	simdat = either full or observed dataset
#	add some descriptive info - parameter values for observed data?


gen.2RRDD<-function(n,s,mu1,mu2,sigma.r1,sigma.r2=sigma.r1,rho.r1r2,rho.r1Y0,rho.r2Y0=rho.r1Y0,
lambda.r1=1,lambda.r2=1,lambda.Y0=1,cut1.quantile,cut2.quantile,cut1.value,cut2.value,
param.Y0,param.eff.T2,param.eff.T1,param.eff.T3,param.eff.T,
sigma.E,sigma.noise,data.type="observed") {

# load needed libraries
require(MASS) # need for mvrnorm()

# sample id
sid<-rep(1:s,each=n)

# generate 2 ratings and counterfactual outcomes with mv normal distribution
    # true ratings 
	cov.r1r2<-rho.r1r2*sigma.r1*sigma.r2
	Sigma=matrix(c(sigma.r1^2,cov.r1r2,cov.r1r2,sigma.r2^2),c(2,2))
	r.true<-mvrnorm(n*s,mu=c(mu1,mu2),Sigma=Sigma)
	r1.true<-r.true[,1]
	r2.true<-r.true[,2]

    # observed ratings 
	var.r1<-var(r1.true)/lambda.r1-var(r1.true)
	var.r2<-var(r2.true)/lambda.r2-var(r2.true)
	r1.obs<-r1.true+rnorm(n*s,0,sqrt(var.r1))
	r2.obs<-r2.true+rnorm(n*s,0,sqrt(var.r2))

# assign treatment based on specified cut-point in the observed ratings
	if (is.null(cut1.quantile)==FALSE) cut1<-quantile(r1.obs,probs=cut1.quantile)
	if (is.null(cut2.quantile)==FALSE) cut2<-quantile(r2.obs,probs=cut2.quantile)
	if (is.null(cut1.value)==FALSE) cut1<-cut1.value
	if (is.null(cut2.value)==FALSE) cut2<-cut2.value
	r1fail<-1*(r1.obs<cut1)	
	r2fail<-1*(r2.obs<cut2)	           
	T1<-1*(r1fail==1 & r2fail==0)	# upper left quad
	T3<-1*(r2fail==1 & r1fail==0)	# lower rt quad
	T2<-1*(r1fail==1 & r2fail==1) # lower left quad
	T<-1*(T1==1 | T2==1 | T3==1)

# center ratings around cut
	r1<-r1.obs-cut1
	r2<-r2.obs-cut2
#	r1<-(r1.obs-mean(r1.obs))
#	r2<-(r2.obs-mean(r2.obs))
	cut1.old<-cut1
	cut2.old<-cut2
	cut1<-0
	cut2<-0

# generate counterfactuals
    # true counterfactual
	dm<-data.frame(intercept=1,r1=r1,r2=r2,r1r2=r1*r2,r1.2=r1^2,r2.2=r2^2)
	Y0<-apply(sweep(dm,MARGIN=2,param.Y0,'*'),1,sum) + rnorm(n,0,sigma.noise)    # observed ratings and counterfactual
	var.Y0<-var(Y0)/lambda.Y0-var(Y0)
	Y0.obs<-Y0+rnorm(n*s,0,sqrt(var.Y0))

# generate impacts
	eff.err<-rnorm(n*s,0,sigma.E)  
         # effect  
	E.true.T<- apply(sweep(dm,MARGIN=2,param.eff.T,'*'),1,sum) + eff.err
	E.true.T2<-apply(sweep(dm,MARGIN=2,param.eff.T2,'*'),1,sum) + eff.err
	E.true.T1<-apply(sweep(dm,MARGIN=2,param.eff.T1,'*'),1,sum) + eff.err 
       E.true.T3<-apply(sweep(dm,MARGIN=2,param.eff.T3,'*'),1,sum) + eff.err       
                                                               
# generate potential outcomes under treatment
	Y1.T2<-Y0+E.true.T2
	Y1.T1<-Y0+E.true.T1  
	Y1.T3<-Y0+E.true.T3
	
 	var.Y1.T2<-var(Y1.T2)/lambda.Y0-var(Y1.T2)
 	var.Y1.T1<-var(Y1.T1)/lambda.Y0-var(Y1.T1)	                                 
 	var.Y1.T3<-var(Y1.T3)/lambda.Y0-var(Y1.T3)
  
	Y1.obs.T2 <- Y1.T2+rnorm(n*s,0,sqrt(var.Y1.T2))
	Y1.obs.T1 <- Y1.T1+rnorm(n*s,0,sqrt(var.Y1.T1))
	Y1.obs.T3 <- Y1.T3+rnorm(n*s,0,sqrt(var.Y1.T3))                            

# assign potential outcome based on treatment status
  Y.obs<-Y0.obs
  Y.obs[T1==1]<-Y1.obs.T1[T1==1] 
  Y.obs[T3==1]<-Y1.obs.T3[T3==1]
  Y.obs[T2==1]<-Y1.obs.T2[T2==1]

# create dataframes (full and observed)
fulldat <-data.frame(
sampleid=sid,rating1=r1,rating1.obs=r1.obs,rating2=r2,rating2.obs=r2.obs,
r1fail,r2fail,cut1,cut2,T1,T2,T3,T,
Y0,Y0.obs,Y1.T2,Y1.T1,Y1.T3,Y1.obs.T2,
Y1.obs.T2,Y1.obs.T3,Y.obs)

obsdat<-data.frame(sid,r1,r2,T,T1,T2,T3,Y.obs)

colnames(obsdat)<-c("sampleID","rating1","rating2","T","T1","T2","T3","Y")


if (data.type=="full") simdat<-fulldat
if (data.type=="observed") simdat<-obsdat

# parameter values
	# sigma
obs.sigma.r1<-sd(r1.obs)
obs.sigma.r2<-sd(r2.obs)
obs.sigma.Y<-sd(Y.obs)
	# rho
obs.rho.r1r2<-cor(r1.obs,r2.obs)
obs.rho.r1Y0<-cor(r1.obs,Y0.obs)
obs.rho.r2Y0<-cor(r2.obs,Y0.obs)
	# lambda
obs.lambda.r1<-summary(lm(r1.obs~r1))$r.squared
obs.lambda.r2<-summary(lm(r2.obs~r2))$r.squared
obs.lambda.Y0<-summary(lm(Y0.obs~Y0))$r.squared
	# effects
	# design matrix for full model 

#data frames with ratings replaced by cutpoint value:
dm.T.r1cut<-data.frame(intercept=dm$intercept,r1=cut1,r2=r2,r1r2=cut1*r2,r1.2=cut1^2,r2.2=r2^2)
dm.T.r2cut<-data.frame(intercept=dm$intercept,r1=r1,r2=cut2,r1r2=r1*cut2,r1.2=r1^2,r2.2=cut2^2)
dm.T.bcut<-data.frame(intercept=dm$intercept,r1=cut1,r2=cut2,r1r2=cut1*cut2,r1.2=cut1^2,r2.2=cut2^2)

# now get individuals' effects IF at frontiers
# note notation variablename.X.Y where X denotes treatment and Y denotes frontier

  # impacts at r1 cut
E.T1.r1<-apply(sweep(dm.T.r1cut,MARGIN=2,param.eff.T1,'*'),1,sum)
  # impacts at r2 cut
E.T3.r2<-apply(sweep(dm.T.r2cut,MARGIN=2,param.eff.T3,'*'),1,sum)

# parameter values
  tau.T <- mean(E.true.T) # marginal for whole pop
  tau.T1<- mean(E.true.T1[T1==1]) # marginal for those in Q1
  tau.T2<- mean(E.true.T2[T2==1]) # marginal for those in Q2
  tau.T3<- mean(E.true.T3[T3==1]) # marginal for those in Q3
  tau.T.r2fail<-mean(E.true.T[r2fail==0]) # marginal for those passing r2
  tau.T.r1fail<-mean(E.true.T[r1fail==0]) # marginal for those passing r1
  # impacts at frontiers...
   	    # means and sds of conditional density
       mu.r1gr2e0<- mean(r1)+rho.r1r2*sd(r1)*(0-mean(r2))/sd(r2)
       mu.r2gr1e0<- mean(r2)+rho.r1r2*sd(r2)*(0-mean(r1))/sd(r1)
	 sd.r1gr2e0<- sd(r1)*sqrt(1-rho.r1r2^2)
	 sd.r2gr1e0<- sd(r2)*sqrt(1-rho.r1r2^2)
	    # normal density function
	 f.gaus<-function(...,mu,sigma) {(1/sqrt(2*pi*sigma^2))*exp(-(...-mu)^2 / (2*sigma^2))} 
	    # conditional effect at r1=0 given r2, and at r2=0 given r1
	 Eatr1 <- function(...) {param.eff.T1[1]+param.eff.T1[3]*...+param.eff.T1[6]*...*...} # fn of r2
	 Eatr2 <- function(...) {param.eff.T3[1]+param.eff.T3[2]*...+param.eff.T3[5]*...*...} # fn of r1
	    # integrate conditional effect * conditional density
	 num.integrand.r1<-function(...) {Eatr1(...)*f.gaus(...,mu=mu.r2gr1e0,sigma=sd.r2gr1e0)}
	 num.r1<-integrate(num.integrand.r1,lower=0,upper=Inf)
	 num.integrand.r2<-function(...) {Eatr2(...)*f.gaus(...,mu=mu.r1gr2e0,sigma=sd.r1gr2e0)}
	 num.r2<-integrate(num.integrand.r2,lower=0,upper=Inf)
	    # integrate conditional density
	 den.integrand.r1<-function(...) {f.gaus(...,mu=mu.r2gr1e0,sigma=sd.r2gr1e0)}
	 den.r1<-integrate(den.integrand.r1,lower=0,upper=Inf)
	 den.integrand.r2<-function(...) {f.gaus(...,mu=mu.r1gr2e0,sigma=sd.r1gr2e0)}
	 den.r2<-integrate(den.integrand.r2,lower=0,upper=Inf)
	    # mean effect at r1
	 tau.T1.r1<-num.r1$value/den.r1$value
	    # mean efect at r2
	 tau.T3.r2<-num.r2$value/den.r2$value
	    # weights for effect among those at either r1 or r2
	 wt.denom<-den.r1$value+den.r2$value
	 wt.r1<-den.r1$value/wt.denom
	 wt.r2<-den.r2$value/wt.denom
	 tau.T.b<-wt.r1*tau.T1.r1+wt.r2*tau.T3.r2

 effects<-c(tau.T,tau.T1,tau.T2,tau.T3,tau.T.r2fail,tau.T.r1fail,tau.T1.r1,tau.T3.r2,tau.T.b)
  names(effects)<-c("ATE","ATE in T1","ATE in T2","ATE in T3","ATE r2pass","ATE r2pass",
"AFE T1","AFE T3","AFE T1 or T3")

# collect all parameter values for outputting
observed<-data.frame(cut1,cut2,obs.sigma.r1,obs.sigma.r2,obs.sigma.Y,
obs.rho.r1r2,obs.rho.r1Y0,obs.rho.r2Y0,
obs.lambda.r1,obs.lambda.r2,obs.lambda.Y0) 

parameters<-list(effects,observed)
names(parameters)<-c("effects","observed")

# stuff to be returned
out<-list(simdat,parameters)
names(out)<-c("data","parameters")
return(out)

} #end function gen.2RRDD


#######################################
#for(sim in c(2,4)) {	
sim<-4									# *** changes with sim
 	n <- 1000 # = sample size
      s <- 500 # = number of samples
 	mu1 <-0 # = mean of rating 1
	mu2 <-0 #= mean of rating 2
	sigma.r1 <- 1 #= sd of true rating 1				# * may change for some runs
	sigma.r2 <- 1 #= sd of true rating 2 (always = ratings)	# * may change for some runs
	rho.r1r2 <- .2 # = correlation between true ratings            # *** changes with sim
#	for (rho.r1r2 in c(.2,.5,.9)){
	cut1.quantile <- 0.30 # = quantile for cutting rating 1       # *** changes with sim
	cut2.quantile <- 0.70 # = quantile for cutting rating 2	# *** changes with sim
	cut1.value <- NULL # = value for cutting rating 1
	cut2.value <- NULL # = value for cutting rating 2
	sigma.noise <- 1  							# * may change for some runs
	sigma.E<- 0 # =sd of random error in effects 
       ddir<-     "C:/Users/shane/MAV/Data/"                           # *** change with user
			#"C:/Users/porterk/Documents/"	
			#"C:/Users/shane/MAV/Data/"
     

if (sim==1) {
  sim.type<-"hom"
  sim.folder<- "Yeq0p4T_0p5R1_R2/" 
  param.Y0<-c(0,0.5,1,0,0,0)        # FIRST VALUE ALWAYS 0, COEFF ON R1, COEFF ON R2, ALWAYS 0, COEFF ON R1sq, COEFF ON R2sq    
  param.eff.T2<-c(0.4,0,0,0,0,0)
  param.eff.T1<-c(0.4,0,0,0,0,0)
  param.eff.T3<-c(0.4,0,0,0,0,0)
  param.eff.T<-c(0.4,0,0,0,0,0)
}

if (sim==2) {
  sim.type<-"hom"
  sim.folder<- "Yeq0p4T_0p5R1_R2_2R1sq_1R2sq/" 
  param.Y0<-c(0,0.5,1,0,2,1) 
  param.eff.T2<-c(0.4,0,0,0,0,0)
  param.eff.T1<-c(0.4,0,0,0,0,0)
  param.eff.T3<-c(0.4,0,0,0,0,0)
  param.eff.T<-c(0.4,0,0,0,0,0)
}

if (sim==3) {
  sim.type<-"htA"
  sim.folder<- "Yeq0p4T_0p5R1_R2_neg0p1TR1_neg0p2TR2/" 
  param.Y0<-c(0,0.5,1,0,0,0)
  param.eff.T2<-c(0.4,-0.1,-0.2,0,0,0)            # COEFF ON T, COEFF OF INTERACTION ON T AND R1, COEFF OF INTERACTION ON T AND R2, ALWAYS 0,
  param.eff.T1<-c(0.4,-0.1,-0.2,0,0,0)	          # COEFF OF INTERACTION ON T AND R1sq, COEFF OF INTERACTION ON T AND R2sq, 
  param.eff.T3<-c(0.4,-0.1,-0.2,0,0,0)            # ALL T VALUES ARE THE SAME – COPY AND PASTE
  param.eff.T<- c(0.4,-0.1,-0.2,0,0,0)
}
if (sim==3.5) {
  sim.type<-"htA"
  sim.folder<- "Yeq0p4T_0p5R1_R2_neg0p4TR1_neg0p4TR2/" 
  param.Y0<-c(0,0.5,1,0,0,0)
  param.eff.T2<-c(0.4,-0.4,-0.4,0,0,0)           # COEFF ON T, COEFF OF INTERACTION ON T AND R1, COEFF OF INTERACTION ON T AND R2, ALWAYS 0,
  param.eff.T1<-c(0.4,-0.4,-0.4,0,0,0)	      # COEFF OF INTERACTION ON T AND R1sq, COEFF OF INTERACTION ON T AND R2sq, 
  param.eff.T3<-c(0.4,-0.4,-0.4,0,0,0)           # ALL T VALUES ARE THE SAME – COPY AND PASTE
  param.eff.T<- c(0.4,-0.4,-0.4,0,0,0)
}
if (sim==4) {
  sim.type<-"htA"
  sim.folder<- "Yeq0p4T_0p5R1_R2_neg0p1TR1_neg0p2TR2_2R1sq_1R2sq_neg0p8TR1sq_neg0p8TR2sq/" 
  param.Y0<-c(0,0.5,1,0,2,1)
  param.eff.T2<-c(0.4,-0.1,-0.2,0,-0.8,-0.8)           # COEFF ON T, COEFF OF INTERACTION ON T AND R1, COEFF OF INTERACTION ON T AND R2, ALWAYS 0,
  param.eff.T1<-c(0.4,-0.1,-0.2,0,-0.8,-0.8)	      # COEFF OF INTERACTION ON T AND R1sq, COEFF OF INTERACTION ON T AND R2sq, 
  param.eff.T3<-c(0.4,-0.1,-0.2,0,-0.8,-0.8)           # ALL T VALUES ARE THE SAME – COPY AND PASTE
  param.eff.T<- c(0.4,-0.1,-0.2,0,-0.8,-0.8)
}

gendat <- gen.2RRDD(n=n,s=s,mu1=mu1,mu2=mu2,sigma.r1=sigma.r1,sigma.r2=sigma.r2,sigma.E=sigma.E,sigma.noise=sigma.noise,rho.r1r2=rho.r1r2,lambda.r1=1,lambda.r2=1,lambda.Y0=1,
cut1.quantile=cut1.quantile,cut2.quantile=cut2.quantile,cut1.value=cut1.value,cut2.value=cut2.value,

param.Y0=param.Y0,param.eff.T2=param.eff.T2,param.eff.T1=param.eff.T1,param.eff.T3=param.eff.T3,param.eff.T=param.eff.T,
data.type="observed") 
dat<-gendat$data
effect_nm <- c("ATE","ATE in T1","ATE in T2","ATE in T3","ATE r1pass", "ATE r2pass","AFE T1","AFE T3","AFE T1 or T2")
param<-data.frame(gendat$parameters,effect_nm)

fname1 <- paste(sim.type,mu1,mu2,sigma.r1,sigma.r2,"_",sigma.E,sigma.noise,cut1.quantile*100,cut2.quantile*100,rho.r1r2*100,round(n/1000),"k",round(s/100),"h",sep="")
fname2 <- paste(sim.type,mu1,mu2,sigma.r1,sigma.r2,"_",sigma.E,sigma.noise,cut1.quantile*100,cut2.quantile*100,rho.r1r2*100,round(n/1000),"k",round(s/100),"h","_p",sep="")

write.csv(dat,file=paste(ddir,sim.folder,fname1,".csv",sep=""))
write.csv(param,file=paste(ddir,sim.folder,fname2,".csv",sep=""))
#}#}


