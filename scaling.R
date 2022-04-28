################
## scaling.R
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
## @description scaling of parameters 
##
## @param x: perform mad() on x
##
##
## @lastChange 2021-04-12
##
## @changes
##   
################

scaling <- function(x){
  if(mad(x) != 0)
    return (mad(x))
  else
    return (1)  ## no scaling
}
