\name{rcw1}
\alias{rcw1}
\title{
Parameter Estimation in LTM using Regression Calibration with Interval-censored Data and Measurement Error
}
\description{
This function should be used when only a single replication of the error-prone covariate is available. If multiple replications are available please see rc() function. This function estimates the parameters of the semiparametric linear transformation model using regression calibration method when time-to-event is interval-censored and a covariate is measured with error. Estimated standard errors of the model parameters are also provided. A right censoring indicator must be available.
}
\usage{
rcw1(datamat, wmat, rfix, gridlen, ntimp, sigma2u)
}
\arguments{
  \item{datamat}{
A data matrix with as many rows as there are subjects (n) in the study. The first column contains left points of the intervals, second column contains right points of the intervals, third column is the right censoring indicator and the fourth (and final) column is the error-free covariate.
}
  \item{wmat}{
A matrix of surrogate measurements for the covariate with measurement error.
}
  \item{rfix}{
A value that characterizes the error density in the linear transformation model. The value 0 corresponds to the Cox PH model and 1 corresponds to Proportional Odds model.
}
  \item{gridlen}{
A non-negative value representing the desired grid length used to divide the observed failure-time intervals.
}
  \item{ntimp}{
Number of failure time imputations required.
}
  \item{sigma2u}{
A known value for the measurement error variance.
}
}
\value{
Parameter estimates of the semiparamteric linear transformation model and their corresponding standard errors.
}
\references{
Mandal, S., Wang, S. and Sinha, S. (2019+). Analysis of Linear Transformation Models with Covariate Measurement Error and Interval Censoring. (accepted, Statistics In Medicine)
}
\author{
Soutrik Mandal, Suojin Wang and Samiran Sinha
}
\note{
Larger sample sizes or imputation numbers will result in longer run times.
}


\examples{
## this function is used in generating epsilon from its CDF
  rsep= function(u,r)
        {
          if(r==0)
            return(  log(-log(1-u))  )
          else
            return(  log((exp(-r*log(1-u))-1)/r)  )
        }


  n= 30 #200     # sample size; small number used for quick demonstration only
  rfix= 1.0
  sigma2u= 0.5
  
  m= 3 #10       # imputed datasets for failure time; small number used for quick demonstration only
  
  nrep= 1        # number of repeated measurement of error prone covariate

  gridlen= 0.1

  result= NULL
  ah=1

  set.seed(ah)

#  z1= rnorm(n, mean= 0, sd=1)
  z1= (rgamma(n,shape=2,scale=2)-4)/sqrt(8)
  z2= rbinom(n,1,0.5)

  ugen= runif(n)
  ep= rsep(ugen,rfix)

  truebeta= c(-1,1)

  logt= -truebeta[1]*z1 -truebeta[2]*z2 + ep + 3
  ttime= exp(logt)

  cen= runif(n,0,0.0001)

## creating tau matrix to locate the actual times in each rows and form the corresponding intervals
  len= 0.15
  taumat= matrix(0,n,10)
  taumat[,10]= 9000000000             # if you're changing this, change rcpos below
  taumat[,2]= cen
  for(i1 in 3:9)
    taumat[,i1]= taumat[,2]+(i1-1)*len

## now forming the intervals
  right1= rep(0,n)
  left= rep(0,n)
  for(i2 in 1:n)
  {
    lenleft= length(which(taumat[i2,2:9]<ttime[i2]))
    leftpot= rep(0,8)
    leftpot[1:lenleft]= 1
    missvec1= c(rbinom(4,1,0.7),rbinom(4,1,0.5))
    left[i2]= max(leftpot*missvec1*taumat[i2,2:9])
    
    if(left[i2]==0)
      left[i2]= taumat[i2,2]
    
    lenright= length(which(taumat[i2,2:10]>ttime[i2]))
    rightpot= rep(0,9)
    rightpot[(9-lenright+1):9]= 1
    missvec2= c(rbinom(4,1,0.7),rbinom(4,1,0.5),1)    
    right1temp= rightpot*missvec2*taumat[i2,2:10]
    right1[i2]= min(right1temp[right1temp!=0])    
  }

  rcpos= which(right1==9000000000)
  lrcpos= length(rcpos)
  notrcpos= (1:n)[-rcpos]
  
  delta_temp= rep(1,n) # del=1 are uncensored observations
  delta_temp[rcpos]= 0
  k= sum(delta_temp)   # this is just the number of data points that are not right censored
  
## measurement error generation
  umat= matrix((rgamma(n*nrep,shape=2,scale=2)-4)*0.5/sqrt(8),n,nrep)

  wmat= z1+umat


  datamat= cbind(left,right1,delta_temp,z2)
  ntimp= m

  library(icemelt)
  out_rcw1= rcw1(datamat, wmat, rfix, gridlen, ntimp, sigma2u)
  out_rcw1

}

