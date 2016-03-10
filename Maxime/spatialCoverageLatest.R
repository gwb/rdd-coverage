##Zach Branson
##Extrapolating with Directional Derivatives

require("MASS")
require("rstan")
#We want to simulate a GP; we'll simulate from the model y = shift + X*B + e
#x will be the locations, and y will be the response surface
#The shift will be zero to the left of the boundary, and
#a constant to the right of the boundary.
#e is from a GP with a squared exponential covariance and nuggent tausq

#First we'll generate data uniformly within the region (0, 10)
#These will be our spatial observations
set.seed(123)
#number of datasets
numberDatasets = 1000
#The rows of x1.star correspond to each of the 1,000 datasets
#the columns correspond to the 20 datapoints in the region (0,5)
x1.star = matrix(runif(20*numberDatasets, 0, 5), nrow = numberDatasets, 20)
#x2.star are the 20 datapoints in the region (5, 10)
x2.star = matrix(runif(20*numberDatasets, 5, 10), nrow = numberDatasets, 20)
#Order each row of x.star,
#so the left-most values are first, the right-most values last:
x1.star = t(apply(x1.star, MARGIN = 1, FUN = sort))
x2.star = t(apply(x2.star, MARGIN = 1, FUN = sort))

#Now we'll combine x1.star and x2.star to just make one matrix of the
#locations, where each row corresponds to a particular data set.
x.star = cbind(x1.star, x2.star)


#Following the Banerjee (2003) Section 5 example, we'll set
#sigmasq = 1, phi = 1.05
#We'll also set beta = 1, tausq = 0.1
sigmasq.star = 1; phi.star = 1.05; beta.star = 1; tausq.star = 0.1
#Now we'll determine the distances between Xs and the corresponding covariance matrix
distances.x.star = list(); covMat.star = list()
#and we'll also record the responses
y.star = matrix(nrow = nrow(x.star), ncol = ncol(x.star))
#Our model will be
# y = mu + X*beta + epsilon
#where mu = 3 for X > 5 and mu = 0 for X < 5
#where epsilon ~ GP(0, covMat.star)
#So, let's record the mu and epsilon
mu.star = ifelse(x.star < 5, 0, 3)
epsilon.star = matrix(nrow = nrow(x.star), ncol = ncol(x.star))

#Creating the covariance matrix for each dataset,
#as well as the error terms (epsilon.star) and response surface (y.star)
for(i in 1:nrow(x.star)){
  distances.x.star[[i]] = as.matrix(dist(x.star[i,]))
  covMat.star[[i]] = sigmasq.star*exp(-phi.star*distances.x.star[[i]]^2) +
    diag(tausq.star, ncol(x.star))
  #Again, epsilon ~ GP(0, covMat.star), so
  epsilon.star[i,] = mvrnorm(n = 1, mu = rep(0, ncol(x.star)),
                         Sigma = covMat.star[[i]])
  #And thus our responses are
  y.star[i,] = mu.star[i,] + x.star[i,]*beta.star + epsilon.star[i,]
}

#Now we'll divide up the response by region.
y1.star = y.star[, 1:20]; y2.star = y.star[, 21:40]

#We have several options for estimating the treatment effect:
# 1) We can fit a shifted GP model
# 2) We can fit two GPs (one in each region) and extrapolate
# 3) We can fit a derivative process in each region and extrapolate

##########################
#####SHIFTED GP MODEL#####
##########################

#In the fitGPShift() function, we use a "Stan fit object" called gpShiftFitTest
#so let's create gpShiftFitTest.
#gpShiftFitTest runs the C++ code to fit the model, and then
#the fitGPShift() function points to gpShiftFitTest so it doesn't have to
#run the C++ code again (which takes some time).
gpShiftFitTest = stan(file = "shiftedGPFit.stan",
                      data =list(x=x.star[1,],y=y.star[1,],N=ncol(x.star)),
                      iter=1000,chains=1)
#Note: Currently only runs one chain with 1,000 iterations;
#this has worked fine so far, but for more complicated models might want to do more
#chains/iterations.

#Function to fit the shifted GP model
#x1 are the values in the (0,5) region; x2 are values in the (5, 10) region
#same goes for y1, y2.
fitGPShift = function(x1, x2, y1, y2){
  data.listShift = list(x=c(x1,x2),
                        y=c(y1,y2),
                        N=length(x1)+length(x2))
  fit = stan(fit = gpShiftFitTest,
             data = data.listShift,
             iter = 1000, chains = 1)
  print(fit)
  fit.extract = extract(fit, permuted=TRUE)
  return(fit.extract)
}

#Create a list which contains the 1,000 fits for the shifted GP model
#(i.e., the fit of the model on each of the 1,000 datasets)
#This takes some time; created a list and wrote a for-loop to add to the list
#so you can stop the code and add to the list later, if necessary.
gpShiftFits = list()
for(i in 1:1000){
  gpShiftFits[[i]] = fitGPShift(x1=x1.star[i,], x2 = x2.star[i,],
                                y1 = y1.star[i,], y2 = y2.star[i,])
}

#A 1000 x 3 matrix which contains the different estimates for the shift parameter
#and the corresponding credible intervals.
#The first column are the estimates (posterior median);
#the second column is the lower bound; the third column is the upper bound
gpShiftFitsMatrix = matrix(nrow=1000,ncol=3)
for(i in 1:1000){
  gpShiftFitsMatrix[i,1] = as.numeric(quantile(gpShiftFits[[i]]$shift,
                                               probs = 0.5))
  gpShiftFitsMatrix[i,2] = as.numeric(quantile(gpShiftFits[[i]]$shift,
                                               probs = c(0.025,0.975))[1])
  gpShiftFitsMatrix[i,3] = as.numeric(quantile(gpShiftFits[[i]]$shift,
                                               probs = c(0.025,0.975))[2])
}

#The true treatment effect is 3; which CIs contain 3?
length(which(gpShiftFitsMatrix[,2] <= 3 & gpShiftFitsMatrix[,3] >= 3))
#mean CI length
mean(gpShiftFitsMatrix[,3]-gpShiftFitsMatrix[,2])

#####################
#####FIT TWO GPS#####
#####################

##Now we'll fit two GPs (one in each region) and extrapolate
#Similar to fitting the shifted GP model,
#we'll run the C++ code once,
#and then create functions that point to these "Stan fit objects".
gpFit1.test = stan(file = "gp2fit.stan",
                   data=list(x=x1.star[1,],
                             y=y1.star[1,],
                             N = length(x1.star[1,])),
                   iter=1000,chains=1)
gpFit2.test = stan(file = "gp2fit.stan",
                   data=list(x=x2.star[1,],
                             y=y2.star[1,],
                             N = length(x2.star[1,])),
                   iter=1000,chains=1)

#Function for fitting a GP in the left-hand (0, 5) region
fitGP1 = function(x1, y1){
  data.list1 = list(x=x1,
                    y=y1,
                    N=length(x1))
  
  gpFit1 <- stan(fit=gpFit1.test,
                 data=data.list1,
                 iter=1000, chains = 1)
  #Note: permuted = TRUE means the draws are taken AFTER burn-in
  #could set permuted = FALSE to keep ALL draws.
  gpFit1.extract = extract(gpFit1, permuted=TRUE)
  return(gpFit1.extract)
}
#Function for fitting a GP in the right-hand (5, 10) region
fitGP2 = function(x2, y2){
  data.list2 = list(x=x2,
                    y=y2,
                    N=length(x2))
  gpFit2 <- stan(fit=gpFit2.test,
                 data=data.list2,
                 iter=1000, chains = 1)
  gpFit2.extract = extract(gpFit2, permuted=TRUE)
  return(gpFit2.extract)
}

#Create a list, where each item corresponds to fitting a GP to a different dataset
#We create two different lists, one for each region, so we can analyze any noticeable
#differences between the left-hand and right-hand GPs.
gpFit1s = list(); gpFit2s = list()
for(i in 1:nrow(x1.star)){
  gpFit1s[[i]] = fitGP1(x1.star[i,], y1.star[i,])
  gpFit2s[[i]] = fitGP2(x2.star[i,], y2.star[i,])
}

#Draw from the posterior of each GP
#postPredY1Draws are posterior predictive draws for the left-hand GP;
#postPredY2Draws are posterior predictive draws for the right-hand GP.
postPredY1Draws = postPredY2Draws = replicate(length(gpFit1s),
                                      matrix(nrow = length(gpFit1s[[1]]$mu),
                                             ncol = ncol(x1.star)),
                                      simplify = FALSE)
#The index m corresponds to each of the 500 posterior draws for each parameter
#(i.e., mean of the GP, sigmasq, and tausq).
#Thus, for each draw from a parameter's posterior, we are drawing from
#the posterior predictive distribution for the response surface.
for(i in 1:length(postPredY1Draws)){
  for(m in 1:nrow(postPredY1Draws[[1]])){
    postPredY1Draws[[i]][m,] = mvrnorm(1, mu = gpFit1s[[i]]$mu_vec[m,],
                               Sigma = diag(gpFit1s[[i]]$sigmasq[m] + gpFit1s[[i]]$tausq[m],
                                            ncol(x1.star)))
    postPredY2Draws[[i]][m,] = mvrnorm(1, mu = gpFit2s[[i]]$mu_vec[m,],
                               Sigma = diag(gpFit2s[[i]]$sigmasq[m] + gpFit2s[[i]]$tausq[m],
                                            ncol(x2.star)))
  }
}




#joint kriging

#  mu, beta, sigmasq, tausq, and phi are all parameters, and 500 draws
#of each are already obtained from the Stan model.
#  X is spatial data
#  Y is response surface data
#  predX are new spatial data (that we have to do kriging on, because we
#don't have the corresponding Y values for each predX value).
#For example, if our right-most observation in the LEFT-HAND region is 4.5,
#predX will be values between 4.5 and 5.
drawKrigingY = function(n, mu, beta, sigmasq, tausq, phi, X, Y, predX){
  
  #First we need distances among the X values
  distances = matrix(nrow=length(X), ncol=length(X))
  for(i in
      1:length(X)){
    for(j in 1:length(X)){
      distances[i,j] = X[i] - X[j]
    }
  }
  #and now we need the distances among the *new* X values
  predDistances = matrix(nrow=length(predX), ncol=length(predX))
  for(i in 1:length(predX)){
    for(j in 1:length(predX)){
      predDistances[i,j] = predX[i] - predX[j]
    }
  }
  #and also we need the distances between each of the X values and the *new* X values
  predXDistances = t(sapply(X, function(x) x - predX))
  
  #Recall that mu is one of our paramters whose posterior we already drew from
  #using the Stan model.
  #Thus, for each posterior draw of the parameters, we'll draw from the posterior
  #predictive distribution for each new X (i.e., each predX)
  #Thus, we'll have length(mu)-many draws from the posterior predictive distribution.
  posteriorPredictiveDraws = matrix(nrow = length(mu), ncol = length(predX))
  for(m in 1:length(mu)){
    
    #the covariance matrix for the Xs is
    estCovMat = sigmasq[m]*exp(-phi[m]*distances^2) + diag(tausq[m], length(X))
    #and the covariance matrix for the new Xs is
    predCovMat = sigmasq[m]*exp(-phi[m]*predDistances^2)
    #and the covariance matrix between the Xs and new Xs is
    predCrossCovMat = sigmasq[m]*exp(-phi[m]*predXDistances^2) 
    
    #using conditional MVN theory, we can find the distribution of
    #p(predX | X, Y, theta)
    #the mean of this distribution is
    predYMean = mu[m] + predX*beta[m] + t(predCrossCovMat)%*%solve(estCovMat)%*%(Y - (mu[m] + X*beta[m]))
    #and the covariance matrix of this distribution is
    predYCovMat = predCovMat - t(predCrossCovMat)%*%solve(estCovMat)%*%predCrossCovMat
    #Therefore, using the above mean and covariance, we can draw from the posterior
    #predictive distribution for each predX.
    posteriorPredictiveDraws[m,] = mvrnorm(n = n, mu = predYMean, Sigma = predYCovMat)
  }
  return(posteriorPredictiveDraws)
} 


#draws from the posterior predictive distribution of both GPs (one in each region)
krigingDrawsY1=krigingDrawsY2=list() 
predX1 = predX2 = matrix(nrow = nrow(x1.star), ncol = 50)

#predX1 is the right-most x value in the LEFT-HAND region to 5
#predX2 is 5 up to the left-most x value in the  RIGHT-HAND region.
#Thus, krigingDrawsY1 are the y values that we would predict if
#we were given, say, x values between 4.5 and 5
#And krigingDraws Y2 are the y values that we would predict if
#we were given, say, x values between 5 and 5.5
for(i in 1:length(gpFit1s)){
  predX1[i,] = seq(x1.star[i,ncol(x1.star)], 5, length = 50)
  predX2[i,] = seq(5, x2.star[i, 1], length = 50)
  krigingDrawsY1[[i]] = drawKrigingY(n=1,gpFit1s[[i]]$mu,gpFit1s[[i]]$beta,
                                             gpFit1s[[i]]$sigmasq,gpFit1s[[i]]$tausq,
                                             gpFit1s[[i]]$phi,
                                             X=x1.star[i,],Y=y1.star[i,],predX=predX1[i,])
  krigingDrawsY2[[i]] = drawKrigingY(n=1,gpFit2s[[i]]$mu,gpFit2s[[i]]$beta,
                                             gpFit2s[[i]]$sigmasq,gpFit2s[[i]]$tausq,
                                             gpFit2s[[i]]$phi,
                                             X=x2.star[i,],Y=y2.star[i,],predX=predX2[i,])
}

#Averaging over the kriging draws so we can have an estimated y for each predX
krigingDrawsY1Mean = krigingDrawsY2Mean = 
  matrix(nrow=length(krigingDrawsY1), ncol = length(predX1))
for(i in 1:length(krigingDrawsY1)){
  krigingDrawsY1Mean[i,] = colMeans(krigingDrawsY1[[i]])
  krigingDrawsY2Mean[i,] = colMeans(krigingDrawsY2[[i]])
}

##PLOT SHOWING EXTRAPOLATION VIA TWO GPs
plot(c(x1,x2), c(y1,y2))
lines(x1, colMeans(postPredY1Draws))
lines(x2, colMeans(postPredY2Draws))
abline(v = 5)
#interal lines
for(i in 1:nrow(krigingDrawsY1)){
  lines(predX1, krigingDrawsY1[i,], col = "gray")
  lines(predX2, krigingDrawsY2[i,], col = "gray")
}
#mean lines
lines(predX1, krigingDrawsY1Mean, col = "purple")
lines(predX2, krigingDrawsY2Mean, col = "purple")

#Plot showing a few extrapolations via 2 GPs
par(mfrow=c(3,3))
for(i in 1:9){
  plot(x.star[i,],y.star[i,])
  lines(predX1, colMeans(krigingDrawsY1[[i]]), col = "purple")
  lines(predX2, colMeans(krigingDrawsY2[[i]]), col = "purple")
  abline(v=5)
}

#Our estimated treatment effect will be the difference between two values:
#   1) The value we would expect at x = 5 for the left-hand GP
#   2) The value we would expect at x = 5 for the right-hand GP

#estEffect2GPs contains the estimated treatment effect via the 2 GPs method.
estEffect2GPs = matrix(nrow = length(krigingDrawsY1),
                             ncol = nrow(krigingDrawsY1[[1]]))
for(i in 1:length(krigingDrawsY1)){
  for(j in 1:nrow(krigingDrawsY1[[1]])){
    estEffect2GPs[i,j] = krigingDrawsY2[[i]][j,1] - krigingDrawsY1[[i]][j,length(predX1[i,])]
  }
}

#Now we'll create a matrix that contains the estimated effect and 95% CI for
#each simulated data set.
#First column corresponds to estimated treatment effects;
#second and third columns are lower and upper bounds of CIs
estEffect2GPsMatrix = matrix(nrow=length(krigingDrawsY1), ncol = 3)
for(i in 1:length(krigingDrawsY1)){
  estEffect2GPsMatrix[i,1] = as.numeric(quantile(estEffect2GPs[i,],
                                                 probs = 0.5))
  estEffect2GPsMatrix[i,2] = as.numeric(quantile(estEffect2GPs[i,],
                                                       probs = c(0.025,0.975)))[1]
  estEffect2GPsMatrix[i,3] = as.numeric(quantile(estEffect2GPs[i,],
                                                       probs = c(0.025,0.975)))[2]
}
#How many CIs contain the true treatment effect of 3?
length(which(estEffect2GPsMatrix[,2] <= 3 & estEffect2GPsMatrix[,3] >= 3))
#Mean CI length
mean(estEffect2GPsMatrix[,3] - estEffect2GPsMatrix[,2])

#Can think about which datasets contain observations closest to the boundary:
#Which datasets have the largest gaps between regions?
xGapDistances = x2.star[,1]-x1.star[,ncol(x1.star)]
xGapCutOff = as.numeric(quantile(xGapDistances, probs = 0.5))
smallGapIndex = which(xGapDistances < xGapCutOff)
largeGapIndex = which(xGapDistances >= xGapCutOff)
#Is the coverage for datasets with a small gap better than
#that of datasets with a large gap?
length(which(estEffect2GPsMatrix[smallGapIndex,2] <= 3 & estEffect2GPsMatrix[smallGapIndex,3] >= 3))
length(which(estEffect2GPsMatrix[largeGapIndex,2] <= 3 & estEffect2GPsMatrix[largeGapIndex,3] >= 3))

##################################
##DIRECTIONAL DERIVATIVE PROCESS##
##################################

#The easiest way to estimate a directional derivative is to estimate the directive.
#We do this by performing kriging on our observations closest to the boundary,
#and then perform kriging on an x value very close to this observation,
#and then calculate the slope of the estimated y for each of these points.
#Banerjee (2003) shows that this is a close approximation
#to the "true directional derivative"

#h is how far away we perform kriging from our observation closest to the boundary
#We want h to be small (Banerjee 2003 suggests h = 0.01)
#u is the "direction"; because we are in one dimension, u = 1, but u can be
#more complicated for high dimensions.
#For example, if we wanted to estimate the derivative to the northeast in two-dimensions,
#u = (0.5, 0.5).
h = 0.01; u = 1
#Observation in (0,5) closest to the boundary
predX1Closest = x1.star[1,ncol(x1.star)]
#Observation in (5, 10) closest to the boundary
predX2Closest = x2.star[1,1]

#The x values that we will perform kriging on
#x1H and x2H are simply perturbations around the observations closest to the boundary
x1H = predX1Closest - h*u
x2H = predX2Closest + h*u

#Now we'll obtain kriging draws for these perturbed values
predY1andYH=predY2andYH=replicate(1000, matrix(nrow=500,ncol=2),simplify=FALSE)
#and also estimate the derivative in the left-hand and right-hand regions
estDeriv1=estDeriv2=
  y1EstDeriv=y2EstDeriv=
  estDerivEffect=matrix(nrow=1000,ncol=500)
#We'll do this for each of the 1,000 datasets
for(i in 1:1000){
  #These are our observed values closest to the boundary
  predX1Closest = x1.star[i,ncol(x1.star)]
  predX2Closest = x2.star[i,1]
  #These are our perturbations around these closest observed values
  x1H = predX1Closest - h*u
  x2H = predX2Closest + h*u
  #Perform kriging on the observations closest to the boundary and
  #the perturbed values

  #left-hand region
  predY1andYH[[i]] = drawKrigingY(n=1, gpFit1s[[i]]$mu,
                                       gpFit1s[[i]]$beta,
                                       gpFit1s[[i]]$sigmasq,
                                       gpFit1s[[i]]$tausq,
                                       gpFit1s[[i]]$phi,
                                       X = x1.star[i,], Y = y1.star[i,],
                                       predX = c(x1H, predX1Closest))
  #right-hand region
  predY2andYH[[i]] = drawKrigingY(n=1, gpFit2s[[i]]$mu,
                                       gpFit2s[[i]]$beta,
                                       gpFit2s[[i]]$sigmasq,
                                       gpFit2s[[i]]$tausq,
                                       gpFit2s[[i]]$phi,
                                       X = x2.star[i,], Y = y2.star[i,],
                                       predX = c(predX2Closest, x2H))
  
  #The estimated derivative on both sides
  estDeriv1[i,] = (predY1andYH[[i]][,2] - predY1andYH[[i]][,1])/h
  estDeriv2[i,] = (predY2andYH[[i]][,1] - predY2andYH[[i]][,2])/h

  #Now that we've estimated the derivative in the left-hand and right-hand regions,
  #we extrapolate to the boundary via the estimated derivative 
  #In other words, we are again estimating what the y value should be at x = 5
  #given the estimated derivative for the left-hand side, and doing the
  #same thing for the right-hand side.
  y1EstDeriv[i,] = predY1andYH[[i]][,2] + (5-predX1Closest)*estDeriv1[i,]
  y2EstDeriv[i,] = predY2andYH[[i]][,1] + (predX2Closest-5)*estDeriv2[i,]
  
  #The estimated treatment effect is the difference between the value we estimate
  #at x = 5 for the left-hand side and that of the right-hand side.
  estDerivEffect[i,] = y2EstDeriv[i,]-y1EstDeriv[i,]
}

#1000 x 3 matrix, where the first column contains the estimated treatment effect
#via the directional derivative process, and the second and third columns
#contain the lower and upper bounds of the corresponding credible intervals.
estDerivEffectMatrix = matrix(nrow=1000,ncol=3)
for(i in 1:1000){
  estDerivEffectMatrix[i,1] = mean(estDerivEffect[i,])
  estDerivEffectMatrix[i,2] = as.numeric(quantile(estDerivEffect[i,],
                                                      prob = c(0.025,0.975))[1])
  estDerivEffectMatrix[i,3] = as.numeric(quantile(estDerivEffect[i,],
                                                      prob = c(0.025,0.975))[2])
}

#Which CIs contain the true treatment effect of 3?
length(which(estDerivEffectMatrix[,2] <= 3 & estDerivEffectMatrix[,3] >= 3))
#average length of CIs
mean(estDerivEffectMatrix[,3] - estDerivEffectMatrix[,2])
#average length of CIs for datasets with observations far from the boundary
#as well as for datasets with observations close to the boundary.
mean(estDerivEffectMatrix[largeGapIndex,3]-estDerivEffectMatrix[largeGapIndex,2])
mean(estDerivEffectMatrix[smallGapIndex,3]-estDerivEffectMatrix[smallGapIndex,2])

#Plots showing the directional derivative process for some datasets
par(mfrow=c(3,3))
for(i in 1:9){
  plot(c(x1.star[i,],x2.star[i,]),c(y1.star[i,],y2.star[i,]),
       xlab = "locations", ylab = "response")
  lines(x1.star[i,], colMeans(postPredY1Draws[[i]]), col = "purple")
  lines(x2.star[i,], colMeans(postPredY2Draws[[i]]), col = "purple")
  lines(predX1[i,], colMeans(krigingDrawsY1[[i]]), col = "purple")
  lines(predX2[i,], colMeans(krigingDrawsY2[[i]]), col = "purple")
  lines(c(x1.star[i,20], 5), c(mean(predY1andYH[[i]][,2]), mean(y1EstDeriv[i,])),
        col="green")
  lines(c(5, x2.star[i,1]), c(mean(y2EstDeriv[i,]), mean(predY2andYH[[i]][,2])),
        col="green")
  abline(v=5)
}
