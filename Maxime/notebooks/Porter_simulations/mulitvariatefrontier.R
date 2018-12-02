#### FUNCTION mf: does estimation for both 2-RDD parameters using either the multivariate or frontier approach
# -- currently only coded for global analysis

# PARAMETERS:
# analysis.type = "global" or "local"
# method = "surface" (i.e. multivariate) or "frontier" (i.e. subset)
# het = TRUE if effects vary with the ratings
# dat = dataset with all samples
# s = number of samples
# model.r1 = surface response model functional form for estimating param at rating 1 frontier
# model.r2 = "										  " rating 2 froniter
# num.range = range for banwidth selection (NULL if analysis.type=="global")
# bandwidths = vector of percentiles for candidate bandwidths (NULL if analysis.type=="global")
# models = vector of candidate models (NULL if analysis.type="global")
# fn.out = filename for writing out results
# odir = output directory
# sbonly = TRUE if only running bw selection without estimation
# estonly = TRUE if reading in previous bw selection results for estimation
# incl.var = TRUE if including covariates IN ESTIMATION MODEL - NOT BW SELECT MODELS

mf<-function(analysis.type,method,het,dat,s,model.r1,model.r2,num.range,bandwidths,models,fn.out,fn.bw,odir,sbonly,estonly,incl.covar) {

	if (analysis.type=="local" & estonly==FALSE) {
		source(paste(pdir,"select.bandwidth_0507.R",sep="")) 
		sb.r1<-select.bandwidth(method=method,rating=1,num.range=num.range,sampleID=dat$sampleID,Y=dat$Y,T=dat$T,
				  threshold.rating=dat$rating1,other.rating=dat$rating2,
				  bandwidths=bandwidths,models=models,odir=odir,fn=fn.bw)
		sb.r2<-select.bandwidth(method=method,rating=2,num.range=num.range,sampleID=dat$sampleID,Y=dat$Y,T=dat$T,
				  threshold.rating=dat$rating2,other.rating=dat$rating1,
				  bandwidths=bandwidths,models=models,odir=odir,fn=fn.bw)
	}
	if (analysis.type=="local" & estonly==TRUE) {
		sb.r1<-read.csv(paste(odir,fn.bw,"_bw_r",1,".csv",sep=""))
		sb.r2<-read.csv(paste(odir,fn.bw,"_bw_r",2,".csv",sep=""))
	}
     
	if (sbonly==FALSE) {
	out1<-matrix(NA,s,5)
      colnames(out1)<-c("sampleID","estimate","se from glm","rho.hat","sampsize")
	out2<-out1

      for (i in 1:s) {
		print(i)
         sampdat<-dat[dat$sampleID==i,]
	   if (analysis.type=="global") ww.r1<-ww.r2<-999
         if (analysis.type=="local") {
		ww.r1<-as.numeric(sb.r1[i,"Opt2.Stage2.h"])
		ww.r2<-as.numeric(sb.r2[i,"Opt2.Stage2.h"])
	   }
	 if (method=="frontier") {
	   temp.r1<-sampdat[sampdat$rating2>0,]
	   dat.r1<-temp.r1[temp.r1$rating1>-ww.r1 & temp.r1$rating1<ww.r1,]
	   temp.r2<-sampdat[sampdat$rating1>0,]
	   dat.r2<-temp.r2[temp.r2$rating2>-ww.r2 & temp.r2$rating2<ww.r2,]
	   }

			### plot to check
#			plot(sampdat$rating1,sampdat$rating2);abline(v=0);abline(h=0)
#			points(temp.r1$rating1,temp.r1$rating2,col="blue")
#			points(dat.r1$rating1,dat.r1$rating2,col="red")
#			points(temp.r2$rating1,temp.r2$rating2,col="blue")
#			points(dat.r2$rating1,dat.r2$rating2,col="red")

	 if (method=="surface") dat.r1<-dat.r2<-sampdat

      # fit model of response surface
		fitm1<-glm(model.r1,data=dat.r1)
		fitm2<-glm(model.r2,data=dat.r2)

         if (het==FALSE | analysis.type=="local") {
		tau.r1<-coef(fitm1)["T"]
		tau.r2<-coef(fitm2)["T"]
		se.r1<-summary(fitm1)$coef["T","Std. Error"]
		se.r2<-summary(fitm2)$coef["T","Std. Error"]
		}
         if (het==TRUE & analysis.type=="global") {
		eff.cond.r1<-function(...) {coef(fitm1)["T"]+coef(fitm1)["T:rating2"]*...}
		eff.cond.r2<-function(...) {coef(fitm2)["T"]+coef(fitm2)["T:rating1"]*...}
	     if (length(coef(fitm1))>8) {
		eff.cond.r1<-function(...) {coef(fitm1)["T"]+coef(fitm1)["T:rating2"]*...+coef(fitm1)["T:I(rating2^2)"]*...*...}
		eff.cond.r2<-function(...) {coef(fitm2)["T"]+coef(fitm2)["T:rating1"]*...+coef(fitm2)["T:I(rating1^2)"]*...*...}
	     }
		se.r1<-se.r2<-NA
		}
      #  get estimates for conditional densities
		rhohat<-cor(sampdat$rating1,sampdat$rating2)
            mu.r1gr2e0<- mean(sampdat$rating1)+rhohat*sd(sampdat$rating1)*(0-mean(sampdat$rating2))/sd(sampdat$rating2)
            mu.r2gr1e0<- mean(sampdat$rating2)+rhohat*sd(sampdat$rating2)*(0-mean(sampdat$rating1))/sd(sampdat$rating1)
            sd.r1gr2e0<- sd(sampdat$rating1)*sqrt(1-rhohat^2)
            sd.r2gr1e0<- sd(sampdat$rating2)*sqrt(1-rhohat^2)
                
	 # normal density function
            f.gaus<-function(...,mu,sigma) {(1/sqrt(2*pi*sigma^2))*exp(-(...-mu)^2 / (2*sigma^2))} 
               
	if (het==TRUE & analysis.type=="global") {
       # integrate conditional effect * conditional density
             num.integrand.r1<-function(...) {eff.cond.r1(...)*f.gaus(...,mu=mu.r2gr1e0,sigma=sd.r2gr1e0)}
             num.r1<-integrate(num.integrand.r1,lower=0,upper=Inf)
             num.integrand.r2<-function(...) {eff.cond.r2(...)*f.gaus(...,mu=mu.r1gr2e0,sigma=sd.r1gr2e0)}
             num.r2<-integrate(num.integrand.r2,lower=0,upper=Inf)
                    
	# integrate conditional density
		den.integrand.r1<-function(...) {f.gaus(...,mu=mu.r2gr1e0,sigma=sd.r2gr1e0)}
      	den.r1<-integrate(den.integrand.r1,lower=0,upper=Inf)
      	den.integrand.r2<-function(...) {f.gaus(...,mu=mu.r1gr2e0,sigma=sd.r1gr2e0)}
      	den.r2<-integrate(den.integrand.r2,lower=0,upper=Inf)

	# mean effect at r1
      	tau.r1<-num.r1$value/den.r1$value
      # mean efect at r2
		tau.r2<-num.r2$value/den.r2$value
 	} # end het loop
		out1[i,]<-c(i,tau.r1,se.r1,rhohat,dim(dat.r1)[1])
		out2[i,]<-c(i,tau.r2,se.r2,rhohat,dim(dat.r2)[1])
	
	} # end loop thru samples

	both<-rbind(out1,out2)
	parameter<-c(rep("psi.r1",s),rep("psi.r2",s))
	out<-data.frame(parameter,both)

      save(out,file=paste(odir,fn.out,".Rdata",sep=""),compress=FALSE)
  	write.csv(out,file=paste(odir,fn.out,".csv",sep=""))
	return(out)
	} # end sbonly==FALSE

} # end function mf


#########################################################

analysis.type<-"local"						# **
method<-"frontier"  # ** options: "frontier" or "surface"
#for (rho in c("20","50", "90")) {
rho<- 50  								# ** options: "20", "50", "90"
cuts <- "5050" 							# ** options: "5050" "3030" "3070"
#for (cuts in c("3030", "3070")) {
s<-500 # change only if we change num of samples7
#for (sim in c(2,4)) {
sim<-4 								# **
sbonly<-FALSE							# **
estonly<-FALSE
#for (incl.covar in c("TRUE","FALSE")){						# **
incl.covar<-FALSE							# **
if (analysis.type=="global") incl.covar<-TRUE         

### specifcations for bw selection
if (analysis.type=="local") {
	pre.type<-4 # prefix on data folder for local analysis (4 datafiles per sim)
	num.range<-100
	bandwidths<-c(0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.75,0.9)
	models<-c("Y~threshold.rating")
#	if (incl.covar==TRUE) models<-c("Y~threshold.rating+other.rating")
#        IF incl.covar==TRUE INCLUDE COVAR IN FINAL ESTIMATION MODEL ONLY 
}
if (analysis.type=="global") {
	pre.type<-""  
	num.range<-bandwidths<-models<-NULL
}

if (sim==1) {
pre<-"hom"
het<-FALSE
sim.folder<- paste(pre.type,"Yeq0p4T_0p5R1_R2/",sep="")
model.r1<-as.formula("Y~T + rating1 + rating2")
model.r2<-model.r1 }

if (sim==2) {
pre<-"hom"
het<-FALSE
sim.folder<- paste(pre.type,"Yeq0p4T_0p5R1_R2_2R1sq_1R2sq/",sep="")
model.r1<-as.formula("Y~T + rating1 + rating2 + I(rating1^2) + I(rating2^2)") #I(rating1^2) + I(rating2^2) rating1:rating1 + rating2:rating2
model.r2<-model.r1 }

if (sim==3) {
pre<-"htA"
het<-TRUE
sim.folder<- paste(pre.type,"Yeq0p4T_0p5R1_R2_neg0p1TR1_neg0p2TR2/",sep="")
model.r1<-as.formula("Y~T + rating1 + rating2 + T:rating1 + T:rating2")
model.r2<-model.r1 }

if (sim==3.5) {
pre<-"htA"
het<-TRUE
sim.folder<- paste(pre.type,"Yeq0p4T_0p5R1_R2_neg0p4TR1_neg0p4TR2/",sep="")
model.r1<-as.formula("Y~T + rating1 + rating2 + T:rating1 + T:rating2")
model.r2<-model.r1 }

if (sim==4) {
pre<-"htA"
het<-TRUE
#sim.folder<- paste(pre.type,"Yeq0p4T_0p5R1_R2_neg0p1TR1_neg0p2TR2_2R1sq_1R2sq_neg0p08TR1sq_neg0p08TR2sq/",sep="")
#sim.folder<- paste(pre.type,"Yeq0p4T_0p5R1_R2_neg0p1TR1_neg0p2TR2_0p02R1sq_0p01R2sq_neg0p008TR1sq_neg0p008TR2sq/",sep="")
sim.folder<- paste(pre.type,"Yeq0p4T_0p5R1_R2_neg0p1TR1_neg0p2TR2_0p2R1sq_0p1R2sq_neg0p008TR1sq_neg0p008TR2sq/",sep="")
model.r1<-as.formula("Y~T + rating1 + rating2 + I(rating1^2) + I(rating2^2) + rating1:T + rating2:T + I(rating1^2):T + I(rating2^2):T")
model.r2<-model.r1 }

  ## overwrite models if local analysis - same for all sims
if (analysis.type=="local") {
  if (incl.covar==TRUE) {
	model.r1<-as.formula("Y~T+rating1+T:rating1+rating2")
	model.r2<-as.formula("Y~T+rating2+T:rating2+rating1")
  }	
  if (incl.covar==FALSE) {
	model.r1<-as.formula("Y~T+rating1+T:rating1")
	model.r2<-as.formula("Y~T+rating2+T:rating2")
  }
}

pdir<- "H:/K12/MAV/Programs/Current_Programs_for_Analyses/"
ddir<- #"C:/Users/shane/MAV/Data/"   
	 #"H:/K12/MAV/Data/" 
	 "C:/Users/porterk/Documents/"
odir<- paste("H:/K12/MAV/Output/",method,"/",sim.folder,sep="")
fn.out<-paste(pre,"cov",incl.covar,"0011_01",cuts,rho,"1k5h_est",sep="") 
fn.bw<-paste(pre,"0011_01",cuts,rho,"1k5h_est",sep="") 

if (analysis.type=="global") {
	dat<-read.csv(paste(ddir,sim.folder,pre,"0011_01",cuts,rho,"1k5h.csv",sep=""))
}					
if (analysis.type=="local") {
	dat1<-read.csv(paste(ddir,sim.folder,pre,"0011_01",cuts,rho,"5k1h1.csv",sep=""))
	dat2<-read.csv(paste(ddir,sim.folder,pre,"0011_01",cuts,rho,"5k1h2.csv",sep=""))
	dat3<-read.csv(paste(ddir,sim.folder,pre,"0011_01",cuts,rho,"5k1h3.csv",sep=""))
	dat4<-read.csv(paste(ddir,sim.folder,pre,"0011_01",cuts,rho,"5k1h4.csv",sep=""))
	dat<-rbind(dat1,dat2,dat3,dat4)
	dat$sampleID<-rep(1:s,each=5000)
}

runfun<-mf(analysis.type,method,het,dat,s,model.r1,model.r2,num.range,bandwidths,models,fn.out,fn.bw,odir,sbonly,estonly,incl.covar)

summary(runfun)

if (analysis.type=="global") params<-read.csv((paste(ddir,sim.folder,pre,"0011_01",cuts,rho,"1k5h_p.csv",sep=""))) # ****change
if (analysis.type=="local") params<-read.csv((paste(ddir,sim.folder,pre,"0011_01",cuts,rho,"5k1h1_p.csv",sep=""))) # ****change
w.p<-params[,"effect_nm"]
truth.r1<-params[w.p=="AFE T1","effects"]
truth.r2<-params[w.p=="AFE T3","effects"]
tt1<-t.test(runfun[runfun[,1]=="psi.r1","estimate"],mu=truth.r1)
tt2<-t.test(runfun[runfun[,1]=="psi.r2","estimate"],mu=truth.r2)
tt1
tt2

E.r1<-mean(runfun[runfun[,1]=="psi.r1","estimate"])
E.r2<-mean(runfun[runfun[,1]=="psi.r2","estimate"])
bias.r1<-E.r1-truth.r1
bias.r2<-E.r2-truth.r2
var.r1<-var(runfun[runfun[,1]=="psi.r1","estimate"])
var.r2<-var(runfun[runfun[,1]=="psi.r2","estimate"])
MSE.r1<-bias.r1^2+var.r1
MSE.r2<-bias.r2^2+var.r2
var.glm.r1<-mean(runfun[runfun[,1]=="psi.r1","se.from.glm"]^2)
var.glm.r2<-mean(runfun[runfun[,1]=="psi.r2","se.from.glm"]^2)
  # coverage using estimated se
ci.l.r1<-runfun[runfun[,1]=="psi.r1","estimate"]-1.96*runfun[runfun[,1]=="psi.r1","se.from.glm"]
ci.u.r1<-runfun[runfun[,1]=="psi.r1","estimate"]+1.96*runfun[runfun[,1]=="psi.r1","se.from.glm"]
cp1<-1*(truth.r1>ci.l.r1 & truth.r1<ci.u.r1)
ci.l.r2<-runfun[runfun[,1]=="psi.r2","estimate"]-1.96*runfun[runfun[,1]=="psi.r2","se.from.glm"]
ci.u.r2<-runfun[runfun[,1]=="psi.r2","estimate"]+1.96*runfun[runfun[,1]=="psi.r2","se.from.glm"]
cp2<-1*(truth.r2>ci.l.r2 & truth.r2<ci.u.r2)
  # coverage if used true var
cci.l.r1<-runfun[runfun[,1]=="psi.r1","estimate"]-1.96*sqrt(var.r1)
cci.u.r1<-runfun[runfun[,1]=="psi.r1","estimate"]+1.96*sqrt(var.r1)
ccp1<-1*(truth.r1>cci.l.r1 & truth.r1<cci.u.r1)
cci.l.r2<-runfun[runfun[,1]=="psi.r2","estimate"]-1.96*sqrt(var.r2)
cci.u.r2<-runfun[runfun[,1]=="psi.r2","estimate"]+1.96*sqrt(var.r2)
ccp2<-1*(truth.r2>cci.l.r2 & truth.r2<cci.u.r2)
  # quantiles of estimates
qq1<-quantile(runfun[runfun[,1]=="psi.r1","estimate"],probs=c(0.025,0.975))
qq2<-quantile(runfun[runfun[,1]=="psi.r2","estimate"],probs=c(0.025,0.975))
res.r1<-c(E.r1,bias.r1,tt1$p.value,var.r1,var.glm.r1,MSE.r1,sqrt(MSE.r1),mean(cp1),mean(ccp1),qq1)
res.r2<-c(E.r2,bias.r2,tt2$p.value,var.r2,var.glm.r2,MSE.r2,sqrt(MSE.r2),mean(cp2),mean(ccp2),qq2)
res<-rbind(res.r1,res.r2)
colnames(res)<-c("Mean","Bias","p.value","Var","Mn Samp Var","MSE","RMSE","cp w est se","cp with true var","2.5%","97.5%")
res
save(res,file=paste(odir,fn.out,"_res.Rdata",sep=""))
write.csv(res,file=paste(odir,fn.out,"_res.csv",sep=""))
#}#}#}

