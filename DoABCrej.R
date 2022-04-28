################
## DoABCrej.R
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
## @description Do ABC with rejection
##
## @param ssim : Summary statistics from simulated samples, m x s matrix, m>1
## @param sobs : Summary statistics from observed data,  1 x s vector 
## @param parprop  : Proposed values of model parameters, m x p matrix
## @ param tolerance : Proportion of ssim values to be considered for 
##                local linear regression 
##                and ABC (posterior sample size to be returned
##                is tolerance x length(ssim))
##
## rows of ssim and parprop match
##
## @return posteriorsamples : Ordered values of model parameters sampled in posterior 
##                            and their corresponding errrors n(epsilon) x (p+1) vector. 
##                            p dimensions for parameters and 1 for error
##
##
## @lastChange 2021-04-12
##
## @changes
##   added scaling() (2021-04-12)
##   updated 'dist' and 'eucdist' (2021-04-07)
##   Included an option for returning constant number of samples (n) with varying maximum error (epsilon) (2018-06-25)
##
##
## Example Usage:
## nsim = 100
## parprop = runif(nsim,1,10)
## ssim = matrix(NA,nsim,1) 
## for ( i in 1: length(parprop)){
## ssim[i,] = mean(rnorm(nsim,parprop[i],1))
## }
## sobs = mean(rnorm(nsim,parprop[1],1))
## posterior <- DoABCrej(ssim, sobs, parprop, tolerance)
################
DoABCrej <- function(ssim, sobs, parprop, tolerance){
  if (!is.vector(sobs)) 
    stop("'sobs' has to be a (nsim*tolerance) vector.")
#
  nsim = dim(ssim)[1]
  nstats = dim(ssim)[2]
  n = round(nsim*tolerance)

# Scaling of summary statistics
  ssim.scaled <- ssim
  sobs.scaled <- sobs
  for(j in 1:nstats){
    scaled <- scaling(ssim[,j])  ## calculate MAD for each of 'nstats' summary statistics
                           ## default constant for MAD is 1.4826 for uniform distribution
    ssim.scaled[,j] <- ssim[,j]/scaled
    sobs.scaled[j] <- sobs[j]/scaled
  }
  
  ## calculate euclidean distance of sobs.scaled from each instance of ssim
  
  dist = (ssim.scaled - matrix(sobs.scaled, nrow = nsim, ncol = nstats, byrow = TRUE))

  ## calculate euclidean distance
  eucdist = rowSums(dist^2)
  
  ord = order(eucdist)[1:n] # find smallest n distances 
  posteriorsamples <- vector("list", 3)
  posteriorsamples[[1]] <- eucdist[ord]
  posteriorsamples[[2]] <- ssim.scaled[ord,]
  posteriorsamples[[3]] <- parprop[ord,]
  return(posteriorsamples)
}  



