#method = estimation method: "frontier", "fuzzyIV", or "binding"
#rating = rating of param of interest (1 or 2)
#num.range = range for bandwidth selection (number of points on each side of threshold)
#sampleID = vector of sample id's
#threshold.rating = vector of rating 1 values, for all samples
#other.rating = vector of rating 2 values, for all samples
#T = vector of treatment indicators, for all samples
#Y = vector of outcomes, for all samples
#threshold.rating = threshold for parameter of interest (=threshold.rating or other.rating)
#bandwidths = percentiles that define candidate bandwidths
#models = vector of reg models (e.g. c("Y~threshold.rating","Y~threshold.rating+other.rating")
#odir = path for writing out saved results
#fn = filename for writing out bandwidth/model selection results

mn.fn<-function(...) {mean(...,na.rm=TRUE)}

select.bandwidth<-function(method,rating,num.range,sampleID,T,Y,threshold.rating,other.rating,bandwidths,models,odir,fn,plot=FALSE) {

  # load Hmisc packagage for tabling
require(Hmisc)

  # to hold results
ns<-max(sampleID)
MSE.by.sample<-MSEr.by.sample<-MSE.T.by.sample<-vector("list",ns)
bw.select.option1<-bw.select.option2<-matrix(NA,ns,4)
colnames(bw.select.option1)<-c("Opt1.Stage1.p","Opt1.Stage1.h","Opt1.Stage2.p","Opt1.Stage2.h")
colnames(bw.select.option2)<-c("Opt2.Stage1.p","Opt2.Stage1.h","Opt2.Stage2.p","Opt2.Stage2.h")
model.select<-numeric(ns)

for (s in 1:max(sampleID)) {
  print(paste("bw selection for sample=",s))
  dat.samp<-data.frame(Y,T,threshold.rating,other.rating)[sampleID==s,]

  # sample data, which passed other rating
 
  if (method=="frontier") pass.samp<-dat.samp[dat.samp$other.rating>0,]
  if (method=="binding") {
    ratings<-dat.samp[,c("threshold.rating","other.rating")]
    r.c<-apply(ratings,1,min)
    T.c<-1*(r.c<0)
    pass.samp<-data.frame(Y=dat.samp$Y,T.c,threshold.rating=r.c,other.rating1=ratings[,1],other.rating2=ratings[,2])
  }
  if (method=="fuzzyIV") pass.samp<-dat.samp

		### PLOT TO CHECK
	if (plot==TRUE) {
		pass.samp[1:5,]
		plot(pass.samp$threshold.rating,pass.samp$other.rating)
		abline(v=0)
	}

  # identify which observations fall within range
w.r<-which(pass.samp$threshold.rating>0)
w.l<-which(pass.samp$threshold.rating<0)
s.r<-sort.int(pass.samp$threshold.rating[w.r],index.return=TRUE)
s.l<-sort.int(pass.samp$threshold.rating[w.l],decreasing=TRUE,index.return=TRUE)
range.l<-s.l[[2]][1:num.range] # index of first num.range points to left or below cut
range.r<-s.r[[2]][1:num.range] # index of first num.range points to right or above cut

		### PLOT TO CHECK
	if (plot==TRUE) {
		pass.samp[1:5,]
		  # subsetted data
		plot(pass.samp$threshold.rating,pass.samp$other.rating)
		abline(v=0)
		  # points to right
		points(pass.samp$threshold.rating[w.r],pass.samp$other.rating[w.r],col="blue")
		  # points in range
		points(pass.samp$threshold.rating[w.r][range.r],pass.samp$other.rating[w.r][range.r],col="red")

		### for binding score
		plot(pass.samp$threshold.rating,pass.samp$Y)
		abline(v=0)
		  # points to right
		points(pass.samp$threshold.rating[w.r],pass.samp$Y[w.r],col="blue")
		  # points in range
		points(pass.samp$threshold.rating[w.r][range.r],pass.samp$Y[w.r][range.r],col="red")
	}

  # identify values associated with various percentiles to define bounds of bandwidths  
	# using the side of the cut with the fewest points then using same values on other side
left.min<-min(length(w.l),length(w.r))==length(w.l)
ifelse (left.min,
	  bw<- -quantile(pass.samp$threshold.rating[w.l],probs=1-bandwidths),
	  bw<- quantile(pass.samp$threshold.rating[w.r],probs=bandwidths)
	  )
bw.r<- bw
bw.l<- -bw

  ss<-matrix(NA,(num.range),length(bw.l)) # to hold sample sizes for one sample
  MSE<-MSEr<-MSE.T<-matrix(NA,length(bw.l),length(models)) 

  # loop through bandwidths
for (b in 1:length(bw)) {
 # print(paste("bandwidth",b))
  num.mod<-length(models)
  predYr<-predYl<-predTr<-predTl<-matrix(NA,nrow=num.range,ncol=num.mod)

  # within bandwidth loop thru points in range
for (i in 1:num.range) {
 #   print(paste("range pt #",i))
   # rating of point
	rr.i<-pass.samp$threshold.rating[w.r][range.r[i]]
	ll.i<-pass.samp$threshold.rating[w.l][range.l[i]]
   # get data on right or left side of cut
	dat.r<-pass.samp[w.r,]
	dat.l<-pass.samp[w.l,]
 
		### PLOT TO CHECK
	if (plot==TRUE) {
		  # dat.r is data to the right of cut
		plot(dat.r$threshold.rating,dat.r$other.rating,xlim=c(0,10))
		abline(v=0)
		  # range
		points(dat.r$threshold[range.r],dat.r$other.rating[range.r],col="red")
		  # rating of point of focus 
		abline(v=rr.i,col="green")
	}

   # create squared terms
#	threshold.rating.2.r<-dat.r$threshold.rating^2
#	other.rating.2.r<-dat.r$other.rating^2
#	threshold.rating.2.l<-dat.l$threshold.rating^2
#	other.rating.2.l<-dat.l$other.rating^2
   # create indicator that threshold.rating is greater than cut
#      I1.r.threshold.rating<-1*(dat.r$threshold.rating>0)
#      I1.l.threshold.rating<-1*(dat.l$threshold.rating>0)
   # add squared terms and indicator to datasets
	 dat.r2<-dat.r
	 dat.l2<-dat.l
#	dat.r2<-data.frame(Y=dat.r$Y,T=dat.r$T,threshold.rating=dat.r$threshold.rating,other.rating=dat.r$other.rating,threshold.rating.2=threshold.rating.2.r,other.rating.2=other.rating.2.r,I1.threshold.rating=I1.r.threshold.rating)
#	dat.l2<-data.frame(Y=dat.l$Y,T=dat.l$T,threshold.rating=dat.l$threshold.rating,other.rating=dat.l$other.rating,threshold.rating.2=threshold.rating.2.l,other.rating.2=other.rating.2.l,I1.threshold.rating=I1.l.threshold.rating)
   # get data to right or left of point
	dat.r.i<-dat.r2[dat.r2$threshold.rating>rr.i,]
	dat.l.i<-dat.l2[dat.l2$threshold.rating<ll.i,]
   # limit data by bandwidth
	dat.r.bw<-dat.r.i[dat.r.i$threshold.rating<(rr.i+bw.r[b]),]
	dat.l.bw<-dat.l.i[dat.l.i$threshold.rating>(ll.i+bw.l[b]),]
      #ss.l[i,b]<-dim(dat.l.bw)[1]

 		### PLOT TO CHECK
	if (plot==TRUE) {
		  # dat.r2 (same as dat.r) is data to the right of cut
		plot(dat.r2$threshold.rating,dat.r2$other.rating,xlim=c(0,10))
		abline(v=0)
		  # range
		points(dat.r$threshold[range.r],dat.r$other.rating[range.r],col="red")
		  # rating of point of focus 
		abline(v=rr.i,col="green")
		  # data to right of point of focus
		points(dat.r.i$threshold.rating,dat.r.i$other.rating,col="green")
		  # limiting data to bw
		points(dat.r.bw$threshold.rating,dat.r.bw$other.rating,col="green")
		  # change b to 5 and rerun dat.r.bw
		points(dat.r.bw$threshold.rating,dat.r.bw$other.rating,col="green")
	}

   # fit regs
     for (m in 1:num.mod) {

	   form<-as.formula(models[m])

       if (method%in%c("frontier","binding")) {
	   fit.r<-try(fastLm(form,data=dat.r.bw),silent=TRUE)
	   fit.l<-try(fastLm(form,data=dat.l.bw),silent=TRUE)
	   newdat.r<-dat.r2[dat.r2$threshold.rating==rr.i,]
	   newdat.l<-dat.l2[dat.l2$threshold.rating==ll.i,]

 			### PLOT TO CHECK
		if (plot==TRUE) {	
	  	 	plot(dat.r$threshold.rating,dat.r$Y)
	    	 	 # observed data at rr.i
	  	 	points(newdat.r$threshold.rating,newdat.r$Y,col="green",pch=19)
		}

  	   ifelse (class(fit.r)!="try-error", predYr[i,m]<-predict(fit.r,newdata=newdat.r),predYr[i,m]<-NA)
	   ifelse (class(fit.l)!="try-error", predYl[i,m]<-predict(fit.l,newdata=newdat.l),predYl[i,m]<-NA)

 			### PLOT TO CHECK
		if (plot==TRUE) {
			 # predicted value at rr.i
	   		points(rr.i,predYr[i,m],col="red",pch=19)
		}

        } # end if method frontier or binding loop

      if (method=="fuzzyIV") { # only works for one (linear) model
	   # stage 1
         fit.s1.r<-fastLm(T~threshold.rating+other.rating,data=dat.r.bw) # note other.rating is NULL if incl.covar=FALSE
		# don't fit on left side (probT=1)
	   # stage 2
         fit.s2.l<-fastLm(form,data=dat.l.bw)
         fit.s2.r<-fastLm(form,data=dat.r.bw)

         ww.r<-which(dat.r2$threshold.rating==rr.i)
	   newdat.r<-dat.r2[ww.r,]
  	   predTr[i,m]<-predict(fit.s1.r,newdata=newdat.r)
         predYr[i,m]<-predict(fit.s2.r,newdata=newdat.r) 

			### PLOT TO CHECK - STAGE 1
		if (plot==TRUE) {
	  	 	plot(dat.r$threshold.rating,dat.r$T)
	    	 	 # observed data at rr.i
	  	 	points(newdat.r$threshold.rating,newdat.r$T,col="green",pch=19)
			 # predicted value at rr.i
			points(rr.i,predTr[i,m],col="red",pch=19)
			
			### PLOT TO CHECK - STAGE 2
	  	 	plot(dat.r$threshold.rating,dat.r$Y)
	    	 	 # observed data at rr.i
	  	 	points(newdat.r$threshold.rating,newdat.r$Y,col="green",pch=19)
			 # predicted value at rr.i
			points(rr.i,predYr[i,m],col="red",pch=19)
       	}

        } # end if method fuzzyIV loop
	} # end m loop (loop thru models)
} # end loop for points

  # MSE for outcome
if (method%in%c("frontier","binding")) {

		### SHOW THESE ARE SAME NUMBERS...
 	 	#print(cbind(predYr,pass.samp$Y[w.r][range.r])[i,])

  predYrmTruth<-predYr-pass.samp$Y[w.r][range.r]
  predYlmTruth<-predYl-pass.samp$Y[w.l][range.l]
  predYmTruth<-rbind(predYrmTruth,predYlmTruth)
  predYmTruth.Sqr<-predYmTruth^2
  MSE[b,]<-apply(predYmTruth.Sqr,2,mn.fn)
}

  # if fuzzyIV method, MSE for treatment fit too (and only one one side)
if (method=="fuzzyIV") {
    # Stage 1
  predYrmTruth<-predYr-pass.samp$Y[w.r][range.r]
  predYlmTruth<-predYr-pass.samp$Y[w.l][range.l]
  predYrmTruth.Sqr<-predYrmTruth^2 # for when num.bw=1
  predYmTruth<-rbind(predYrmTruth,predYlmTruth)
  predYmTruth.Sqr<-predYmTruth^2
  MSEr[b,]<-apply(predYrmTruth.Sqr,2,mean) # right side only
  MSE[b,]<-apply(predYmTruth.Sqr,2,mean) # using both sides

    # Stage 2
  predTrmTruth<-predTr-pass.samp$T[w.r][range.r]
  predTmTruth<-predTrmTruth
  predTmTruth.Sqr<-predTmTruth^2
  MSE.T[b,]<-apply(predTmTruth.Sqr,2,mean)
}

} # end loop through bandwidths

colnames(MSE)<-colnames(MSEr)<-models
rownames(MSE)<-rownames(MSEr)<-rownames(MSE.T)<-bw.r

 # using both sides of frontier for Stage 2
MSE.by.sample[[s]]<-MSE
minMSE.ind<-1*(MSE==min(MSE))
bw.which<-which(apply(minMSE.ind,1,sum)>=1)
 # using just right side for Stage 2 (same code as above using both sides)
MSEr.by.sample[[s]]<-MSEr
minMSEr.ind<-1*(MSEr==min(MSEr))
bw.which.r<-which(apply(minMSEr.ind,1,sum)>=1)
 # Stage 1 is NA except for fuzzyIV
bw.T.which<-999

  if (method=="fuzzyIV") {
MSE.T.by.sample[[s]]<-MSE.T
minMSE.T.ind<-1*(MSE.T==min(MSE.T))
bw.T.which<-which(apply(minMSE.T.ind,1,sum)>=1)
  }

  # Option 1: 
  # FOR FRONTIER AND BINDING SCORE: min MSE for Y model
  # FUZZY IV: min(MSE for Y model,MSE for T model)- both using right side only
  # same bw selected for both stages

bw.o1<-min(bw.which.r,bw.T.which)

  # Option 2: 
  # FOR FRONTIER AND BINDING SCORE: same as Option 1
  # FOR FUZZY IV: select h for each stage

bw.o2.S1<-bw.T.which
bw.o2.S2<-bw.which

bw.select.option1[s,]<-c(bandwidths[bw.o1],bw.r[bw.o1],bandwidths[bw.o1],bw.r[bw.o1])
bw.select.option2[s,]<-c(bandwidths[bw.o2.S1],bw.r[bw.o2.S1],bandwidths[bw.o2.S2],bw.r[bw.o2.S2])
model.select[s]<- models[apply(minMSE.ind,2,sum)>=1]
} #end loop through samples	

# write out csv file of selection results
bw.mod<-data.frame(bw.select.option1,bw.select.option2,model.select)
#write.csv(bw.mod,file=paste(odir,fn,"_bw_r",rating,".csv",sep=""))
#save(bw.mod,file=paste(odir,fn,"_bw_r",rating,".Rdata",sep=""))
return(bw.mod)

#out<-list(MSE.T.by.sample,MSE.by.sample,bw.select.option1,bw.select.option2,model.select)
#names(out)<-c("T MSE by sample","Y MSE by sample","bw selection option 1","bw selection option 2","model")
#return(out)
} # end function		

