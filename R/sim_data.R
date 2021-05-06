#' Simulated data
#' @description  Function that generates simulated data with five covariates: three confounders,one instrumental variable and one predictor of the outcome (scenario from the paper of Wang) 
#' @param N sample size 
#' @param beta exposure marginal effect on outcome 
#' @return a list with outcome y,exposure x and covariates U
#' @rdname sim_data
#' @export

sim_data <- function(N,beta){

  Coef.X.true=c(1,1,.1,0,.7)
  #beta=.1
  Coef.Y.true= c(1,.1,1,1,0,beta)
  p = 5
  U.mean=rep(0,p)
  U.sigma=diag(1,p)
  
  expit=function(a){exp(a)/(1+exp(a))}
  
  Utemp= MASS::mvrnorm(N,U.mean,U.sigma)

  U=Utemp
  logitpx=U%*%Coef.X.true
  px=expit(logitpx)
  x=rbinom(N,1,px)
  
  logitpy0=U%*%Coef.Y.true[1:p]
  logitpy1=logitpy0+Coef.Y.true[p+1]
  y0=stats::rbinom(N,1,expit(logitpy0))
  y1=stats::rbinom(N,1,expit(logitpy1))
  y=rep(NA,N)
  y[x==0]=y0[x==0]
  y[x==1]=y1[x==1]
  y=as.matrix(y)
  x=as.matrix(x)
  U = data.frame(U)
  return(list(y,x,U))
  
}

