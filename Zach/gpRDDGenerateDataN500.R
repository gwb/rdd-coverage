#FUNCTIONS PROVIDED BY GUILLAUME BASSE

#generates x data
get.x <- function(n=500){
    return(2 * rbeta(n, 2, 4) - 1)
}

#functions that generate y data from various papers
m.lee <- function(x, sig=0.1295){
    m <- ifelse(x < 0,
                0.48 + 1.27 * x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5,
                0.52 + 0.84 * x - 3.00*x^2 + 7.99*x^3  - 9.01*x^4 + 3.56*x^5)
    return(rnorm(length(x), m, sd=sig))
}
m.quad <- function(x, sig=0.1295){
    m <- ifelse(x < 0, 3*x^2, 4*x^2)
    return(rnorm(length(x), m, sd=sig))
}
m.cubic <- function(x, sig=0.1295){
    m <- ifelse(x < 0, 3*x^3, 4*x^3)
    return(rnorm(length(x), m, sd=sig))
}
m.cate.1 <- function(x, sig=0.1295){
    m <- 0.42 + 0.84 * x - 3.00 * x^2 + 7.99 * x^3 - 9.01 * x^4 + 3.56 * x^5 +
        ifelse(x < 0, 0, 0.1)
    return(rnorm(length(x), m, sd=sig))
}
m.cate.2 <- function(x, sig=0.1295){
    m <- 0.42 + 0.84 * x + 7.99 * x^3 - 9.01 * x^4 + 3.56 * x^5 +
        ifelse(x < 0, 0, 0.1)
    return(rnorm(length(x), m, sd=sig))
}
m.ludwig <- function(x, sig=0.1295){
    m <- ifelse(x < 0,
                3.71 + 2.30*x + 3.28*x^2 + 1.45*x^3 + 0.23*x^4 + 0.03*x^5,
                0.36 + 18.49*x - 54.81*x^2 + 74.30*x^3 - 45.02*x^4 + 9.83*x^5)
    return(rnorm(length(x), m, sd=sig))
}
m.curvature <- function(x, sig=0.1295){
    m <- ifelse(x < 0,
                0.48 + 1.27 * x - 0.5*7.18*x^2 + 0.7*20.21*x^3 + 1.1*21.54*x^4 + 1.5*7.33*x^5,
                0.52 + 0.84 * x - 0.1*3.00*x^2 - 0.3*7.99*x^3  - 0.1*9.01*x^4 + 3.56*x^5)
    return(rnorm(length(x), m, sd=sig))
}

#simulate quadratic data
simulateQuadData = function(n, gapSize = 0, sig = 0.1295){
	xs = get.x(n=n)
	ys = m.quad(x = xs, sig = sig)
	sortIndex = sort.int(xs,index.return=TRUE)$ix
	xs = xs[sortIndex]; ys = ys[sortIndex]
	data = data.frame(x = xs, y = ys)
	return(data)
}
simulateCubicData = function(n, gapSize = 0, sig = 0.1295){
  xs = get.x(n=n)
  ys = m.cubic(x = xs, sig = sig)
  sortIndex = sort.int(xs,index.return=TRUE)$ix
  xs = xs[sortIndex]; ys = ys[sortIndex]
  data = data.frame(x = xs, y = ys)
  return(data)
}
#simulates data using the function m.lee
simulateLeeData = function(n, gapSize = 0, sig = 0.1295){
	xs = get.x(n=n)
	ys = m.lee(x = xs, sig = sig)
	sortIndex = sort.int(xs,index.return=TRUE)$ix
	xs = xs[sortIndex]; ys = ys[sortIndex]
	data = data.frame(x = xs, y = ys)
	return(data)
}
#If c = 1, uses m.cate.1; if c = 2, uses m.cate.2
simulateCateData = function(n, c = 1, gapSize = 0, sig = 0.1295){
	xs = get.x(n=n)
	if(c == 1){
		ys = m.cate.1(x = xs, sig = sig)
	}
	if(c == 2){
		ys = m.cate.2(x = xs, sig = sig)
	}
	sortIndex = sort.int(xs,index.return=TRUE)$ix
	xs = xs[sortIndex]; ys = ys[sortIndex]
	data = data.frame(x = xs, y = ys)
	return(data)
}
#simulate ludwig data
simulateLudwigData = function(n, gapSize = 0){
  xs = get.x(n=n)
  ys = m.ludwig(x = xs)
  sortIndex = sort.int(xs,index.return=TRUE)$ix
  xs = xs[sortIndex]; ys = ys[sortIndex]
  data = data.frame(x = xs, y = ys)
  return(data)
}
#simulate curvature data
simulateCurvatureData = function(n, gapSize = 0){
  xs = get.x(n=n)
  ys = m.curvature(x = xs)
  sortIndex = sort.int(xs,index.return=TRUE)$ix
  xs = xs[sortIndex]; ys = ys[sortIndex]
  data = data.frame(x = xs, y = ys)
  return(data)
}

##GENERATE DATASETS##
nsims=1000
set.seed(123)
leeDatasets = replicate(nsims, simulateLeeData(n=500), simplify = FALSE)
for (i in 1:nsims){
    write.csv(leeDatasets[i], file=sprintf("saved_simData/lee_%d.csv", i))
}
set.seed(123)
quadDatasets = replicate(nsims, simulateQuadData(n=500), simplify = FALSE)
for (i in 1:nsims){
    write.csv(quadDatasets[i], file=sprintf("saved_simData/quad_%d.csv", i))
}
set.seed(123)
cubicDatasets = replicate(nsims, simulateCubicData(n=500), simplify = FALSE)
for (i in 1:nsims){
    write.csv(cubicDatasets[i], file=sprintf("saved_simData/cubic_%d.csv", i))
}
set.seed(123)
cate1Datasets = replicate(nsims, simulateCateData(n=500, c = 1), simplify = FALSE)
for (i in 1:nsims){
    write.csv(cate1Datasets[i], file=sprintf("saved_simData/cate1_%d.csv", i))
}
set.seed(123)
cate2Datasets = replicate(nsims, simulateCateData(n=500, c = 2), simplify = FALSE)
for (i in 1:nsims){
    write.csv(cate2Datasets[i], file=sprintf("saved_simData/cate2_%d.csv", i))
}
set.seed(123)
ludwigDatasets = replicate(nsims, simulateLudwigData(n=500), simplify = FALSE)
for (i in 1:nsims){
    write.csv(ludwigDatasets[i], file=sprintf("saved_simData/ludwig_%d.csv", i))
}
set.seed(123)
curvatureDatasets = replicate(nsims, simulateCurvatureData(n=500), simplify = FALSE)
for (i in 1:nsims){
    write.csv(curvatureDatasets[i], file=sprintf("saved_simData/curvature_%d.csv", i))
}

##PERFORM KRIGING##
performKrigingSameParams = function(n = 1, data, stanFit, boundary = 0, length = 50){
  x1 = data$x[which(data$x < boundary)]
  y1 = data$y[which(data$x < boundary)]
  x2 = data$x[which(data$x >= boundary)]
  y2 = data$y[which(data$x >= boundary)]
  
  predX1 = getPredX1(data = data, boundary = boundary, length = length)
  predX2 = getPredX2(data = data, boundary = boundary, length = length)
  
  #first, perform kriging on the left side
  krigingLeft = drawPostPredY(n = n,
                mu = stanFit$mu1,
                beta = stanFit$beta1,
                  sigmasq = stanFit$sigmasq,
                              phi = stanFit$phi,
                              etasq = stanFit$etasq,
                              X = x1,
                              Y = y1,
                              predX = predX1)
  #now perform kriging on the right side
  krigingRight = drawPostPredY(n = n,
                  mu = stanFit$mu2,
                  beta = stanFit$beta2,
                  sigmasq = stanFit$sigmasq,
                              phi = stanFit$phi,
                              etasq = stanFit$etasq,
                              X = x2,
                              Y = y2,
                              predX = predX2)
  
  #krigingLeft and krigingRight are matrices
  #We'll cbind these matrices, so that
  #The first 100 columns are for the left,
  #the next 100 columns for the right.
  krigingMatrix = cbind(krigingLeft, krigingRight)
  return(krigingMatrix)
}
drawPostPredY = function(n = 1, mu, sigmasq, phi, etasq, X, Y, predX,beta){
  
  #First we need distances among the X values
  distances = matrix(nrow=length(X), ncol=length(X))
  for(i in 1:length(X)){
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
  
  #sigmasq and etasq are paramters whose posterior we already drew from
  #using the Stan model.
  #Thus, for each posterior draw of the parameters, we'll draw from the posterior
  #predictive distribution for each new X (i.e., each predX)
  #Thus, we'll have length(sigmasq)-many draws from the posterior predictive distribution.
  posteriorPredictiveDraws = list(length = length(sigmasq))
  for(m in 1:length(sigmasq)){
    
    #the covariance matrix for the Xs is
    estCovMat = sigmasq[m]*exp(-phi[m]*distances^2) + diag(etasq[m], length(X))
    #and the covariance matrix for the new Xs is
    predCovMat = sigmasq[m]*exp(-phi[m]*predDistances^2)
    #and the covariance matrix between the Xs and new Xs is
    predCrossCovMat = sigmasq[m]*exp(-phi[m]*predXDistances^2)
    
    #using conditional MVN theory, we can find the distribution of
    #p(predX | X, Y, theta)
    #the mean of this distribution is
    predYMean = (mu[m] + predX*beta[m]) + t(predCrossCovMat)%*%solve(estCovMat)%*%(Y - mu[m] - X*beta[m] )
    #and the covariance matrix of this distribution is
    predYCovMat = predCovMat - t(predCrossCovMat)%*%solve(estCovMat)%*%predCrossCovMat
    #Therefore, using the above mean and covariance, we can draw from the posterior
    #predictive distribution for each predX.
    posteriorPredictiveDraws[[m]] = mvrnorm(n = n, mu = predYMean, Sigma = predYCovMat)
  }
  #right now, posteriorPredictiveDraws is a list of matrices,
  #and we'd like to collapse this into one big matrix
  posteriorPredictiveDrawsMatrix = do.call(rbind, posteriorPredictiveDraws)
  return(posteriorPredictiveDrawsMatrix)
} 
getPredX1 = function(data, boundary, length = 50){
  x1 = data$x[which(data$x < boundary)]
  predX1 = c(x1, seq(max(x1), boundary, length = length))
  return(predX1)
}
getPredX2 = function(data, boundary, length = 50){
  x2 = data$x[which(data$x >= boundary)]
  predX2 = c(seq(boundary, min(x2), length = length), x2)
  return(predX2)
}
