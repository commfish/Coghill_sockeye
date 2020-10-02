#State Space Model Function with autocorrelation
# code after Fair et al. 2012 pg. 90
mod=function(){
  
  #PRIORS
  #lnalpha ~ dnorm(0,1.0E-6)%_%T(0,10) #uninformative
  #beta ~ dnorm(0,1.0E-6)%_%T(0,10)  #uninformative, normal distribution, constrained to be >0
  #phi ~ dnorm(0,1.0E-6)%_%T(-0.98,0.98) #AR(1) model so phi IS included and does not = zero. uninformative btwn -1 & 1
  lnalpha ~ dnorm(0,1.0E-6)%_%T(1.0E-6, )#uninformative;*#normal distribution with mean 0 and large variance; constrained to be >0 since more biologically conservative (pg. 406)
  beta ~ dnorm(0,1.0E-6)%_%T(1.0E-6,)   #uninformative; normal distrib; constrained to be >0 *          
  phi ~ dnorm(0,1.0E-6)%_%T(-1, 1)#uninformative; btw -1 and 1
  
  resid.red.0 ~ dnorm(0,tau.red)
  sigma.white ~ dunif(0,10)
  
  
  for(y in 1:n) {lnRS[y] ~ dnorm(mean2.lnRS[y],tau.white) }  #Is this a prior or not????
  
  mean2.lnRS[1] <- mean1.lnRS[1] + phi * resid.red.0  
  for (y in 2:n) { mean2.lnRS[y] <- mean1.lnRS[y] + phi * resid.red[y-1] }   #AR1
  
  for(y in 1:n) {  mean1.lnRS[y] <- lnalpha - beta * S[y]  } #This is the Ricker model
  for(y in 1:n) {  resid.red[y]     <- lnRS[y] - mean1.lnRS[y]  }
  for(y in 1:n) {  resid.white[y] <- lnRS[y] - mean2.lnRS[y]  }
  
  alpha <- exp(lnalpha.c) #exponentiate to solve for alpha
  sigma.red <- 1 / sqrt(tau.red)
  tau.white <- 1 / sigma.white / sigma.white
  tau.red <- tau.white * (1-phi*phi)

  lnalpha.c <- lnalpha + (sigma.red * sigma.red / 2)  #adjust for calculating means of R.msy, S.msy etc.for the AR model
  
  S.max <- 1 / beta
  S.eq <- S.max * lnalpha.c 
  
  S.msy <- S.eq * (0.5 - 0.07*lnalpha.c)  #Hilborn approximation to calculate Smsy
  U.msy <- lnalpha.c * (0.5 - 0.07*lnalpha.c)  #Hilborn approximation of U.msy
  R.msy <- S.msy * exp(lnalpha.c - beta * S.msy) #Xinxian's calculation of R.msy
  MSY <- step(R.msy-S.msy)*(R.msy-S.msy) #if R.msy < S.msy then MSY=0.
}


