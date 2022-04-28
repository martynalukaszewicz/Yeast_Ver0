################
## calc_signFst.R
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
## @description calculate signFst. Fixation index (Fst) calculation based on Hartl and Clark 2007.
##              If (pX-pY) < 0, signFst=-Fst, else signFst=Fst
##
## @param pX : Proportion of major allele (allele 1) in population X, nsim x nstat matrix
## @param pY : Proportion of major allele (allele 1) in population Y, nsim x nstat matrix
## @param nsim: Number of simulations of Simulator
## @param nstat: Number of summary statistics 
##
##
## @return signFst : Updated fixation index, with negative sign when (pX-pY) < 0.
##                             
##
## @lastChange
##
## @changes
##
##
## Example Usage:
## pX = matrix(runif(nsim*nstat,0,1),nrow=nsim,ncol=nstat)
## pY = matrix(runif(nsim*nstat,0,1),nrow=nsim,ncol=nstat)
## signFst<- calc_signFst(pX,pY)
################
## dependencies:
require(dplyr) ## %>%
#####################


calc_signFst <- function(pX,pY){
  ## dependencies:
  require(dplyr) ## %>%
  #####################
  ## @param pX : Proportion of major allele (allele 1) in population X, nsim x nstat matrix
  ## @param pY : Proportion of major allele (allele 1) in population Y, nsim x nstat matrix
  ## @param nsim: Number of simulations of Simulator
  ## @param nstat: Number of summary statistics 
  
  pbar <- 0.5*(pX+pY)  ## mean proportion of major allele
  sigmasq <- 0.5*(pX-pbar)^2+0.5*(pY-pbar)^2  ## variance
  Fst <- sigmasq/(pbar*(1-pbar))  ## Fst defined in Hartl and Clark 2007 for bi-alleleic case
  Fst[is.na(Fst)] <- 0  ## if pbar (or 1-pbar) is 0, allele is fixed in both populations, Fst=0
  nstat <- Fst %>% ncol()
  
  for (i in 1:nstat){
    check <- (pX[,i]-pY[,i])<0
    Fst[check,i] <- (-1)*Fst[check,i]
  }
  
  return(Fst)
}
