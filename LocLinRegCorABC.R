################
## LocLinRegCorABC.R
##
## @authors Martyna Lukaszewicz, Erkan Buzbas, Paul Hohenlohe
## @contact martyna@uidaho.edu
##
## @license This is free software distributed without any warranty, without
##          even the implied warranty of merchantability or fitness for a particular
##          purpose. You are welcome to use, modify, or distribute it 
##          under the terms of the GNU General Public License,
##          version 3, as published by the Free Software Foundation.
##          See http://www.r-project.org/Licenses/GPL-3
##
################
## @description Do local linear regression correction
##              (implement Beaumont et al. 2002) for 
##              post-sampling adjustment of 
##              posterior samples obtained by ABC rejection
##              
##
## @param posteriorsamples : List of 3, vector of (n x 1) dist, matrix of (n x nstats)
##                           of summary statistics corresponding to accepted values 
##                           of parameters, and a matrix of (n x npar) accepted values of
##                           parameters as from the posterior in the rejection
##
## @param posteriorsamples[[1]] : Euclidean distances of accepted datasets, 1 x nsample vector
## @param posteriorsamples[[2]] : Proposed scaled accepted summary statistics, nsample x nstats matrix
## @param posteriorsamples[[3]] : Proposed scaled accepted parameters, nsample x nstats matrix
## @param posteriorsamples[[1]] : Summary statistics from observed data,  1 x s vector 
##
## @return posteriorsamples.adj : Posterior samples after local linear regression 
##                                adjustment is performed. Ordered values of model 
##                                parameters sampled in posterior, n x npar matrix
##
##
## @lastChange 2021-04-14
##
## @changes
##  
##  Updated posteriorsamples.adj with kernel=epanechnokov method linear regression correction 
##
## Example Usage:
## posteriorsamples <- vector ("list",4)
## npar <- 2
## nstat <- 4
## nsample <- 10
## posteriorsamples[[1]] <- runif(nsample)
## posteriorsamples[[2]] <- matrix(runif(nsample*nstat), ncol=nstats, nrow=nsample)
## posteriorsamples[[3]] <- matrix(runif(nsample*npar), ncol=npar, nrow=nsample)
## posteriorsamples[[4]] <- runif(nstat)
## posterior.adj <- LocLinregCorABC(posteriorsamples)

################
LocLinRegCorABC <- function(posteriorsamples,kernel="epanechnikov"){
  sstat<- posteriorsamples[[2]]   
  par <- posteriorsamples[[3]]
  nsample <-  dim(par)[1]
  npar <- dim(par)[2]
  ##
  if (kernel == "epanechnikov"){
    ## Epanechnikov kernel weights for local linear regression
    regweights <- 1-(posteriorsamples[[1]]/posteriorsamples[[1]][nsample])^2 
  }
  if (kernel == "gaussian"){
    ## Gaussian kernel weights for local linear regression
    regweights <- 1/sqrt(2 * pi) * exp(-0.5 * (posteriorsamples[[1]]/(0.5*posteriorsamples[[1]][nsample]) )^2)
  }
  
  locreg <- lsfit(sstat,par,wt=regweights)
  pred <- matrix(t(locreg$coeff) %*% c(1,posteriorsamples[[4]]), ncol=npar, nrow=nsample, byrow=TRUE)  ## pred = betas x X 
  res <- locreg$residuals  ## epsilon
  resmean <- apply(res,FUN=mean,2)
  res <- sapply(1:npar,FUN=function(x){res[,x]-resmean[x]})
  pred <- sapply(1:npar,FUN=function(x){pred[,x]+resmean[x]})
  posteriorsamples.adj <- pred + res
  return(posteriorsamples.adj)
}