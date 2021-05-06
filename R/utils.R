#' Internal BACps functions
#' @description Internal BACps helper functions for the 1st Stage
#' @details These functions are not intended for use by users.
#' @name BACps-internal
NULL

##### STAGE 1 ############

########################################
###### Log Likelihood functions  #######
########################################
#' @rdname BACps-internal

logL=function(N,y,U,alpha,delta){
  ones=rbind(rep(1,N))
  logL=t(y)%*%t(delta[1]+(delta[-1])%*%t(U[,alpha!=0]))- ones%*%t(log(ones+exp(delta[1]+(delta[-1])%*%t(U[,alpha!=0])))) # with intercept given by 1
  logL
}

#' @rdname BACps-internal
##############################
###### MH functions  #######
###########################

MH.delta=function(N,y,U,delta0,alpha0,delta.mle.model.alpha0,delta.cov.model.alpha0){
  logL.delta0=logL(N,y,U,alpha0,delta0)
  accept.delta=0
  #a Proposal is a combination of an independent chain (not necessary symmentric) and dependent chain (symmetric)
  thres.prop=runif(1)   
  if(thres.prop>=0.5){
    delta.temp=MASS::mvrnorm(n=1,delta.mle.model.alpha0,delta.cov.model.alpha0) # g(delta)=N(delta.mle, var-cov)
  } else { 
    delta.temp=delta0+MASS::mvrnorm(n=1,rep(0,length(delta0)),diag(10,nrow=length(delta0),ncol=length(delta0))) # delta1=delta0+error so g(delta'|delta*)=N(delta*,var-cov)
  } 
  delta1=delta.temp
  logL.delta1=logL(N,y,U,alpha0,delta1)
  extra.term.delta0=c(delta0)%*%solve(diag(10,nrow=length(delta0),ncol=length(delta0)))%*%cbind(c(delta0))
  extra.term.delta1=c(delta1)%*%solve(diag(10,nrow=length(delta0),ncol=length(delta0)))%*%cbind(c(delta1)) 
  
  if(thres.prop>=0.5){
    ratio.prop1=mvtnorm::dmvnorm(delta0,delta.mle.model.alpha0,delta.cov.model.alpha0)
    ratio.prop2=mvtnorm::dmvnorm(delta1,delta.mle.model.alpha0,delta.cov.model.alpha0)
    BF.delta=exp(logL.delta1-logL.delta0-1/2*extra.term.delta1+1/2*extra.term.delta0)*(ratio.prop1/ratio.prop2)
  } else{
    BF.delta=exp(logL.delta1-logL.delta0-1/2*extra.term.delta1+1/2*extra.term.delta0)
  }  # if the proposal density is not symmetric i should have another term there
  if(BF.delta>=1 || (BF.delta<1 & runif(1)<BF.delta)){
    delta0=delta1  
    accept.delta=1
  }
  list(delta0=delta0,accept.delta=accept.delta)
}

#' @rdname BACps-internal

MH.alpha.s1_a=function(N,y,U,alpha0,alpha1,delta0,delta_t_1,q_delta_j){
  accept.alpha=0
  logL.alpha0=logL(N,y,U,alpha0,delta0)
  logL.alpha1=logL(N,y,U,alpha1,delta_t_1)
  ratio.prior.delta0=exp((-1/2)*c(delta0)%*%diag(1/10^10,nrow=length(delta0),ncol=length(delta0))%*%c(delta0))
  ratio.prior.delta1=exp(-1/2*c(delta_t_1)%*%diag(1/10^10,nrow=length(delta_t_1),ncol=length(delta_t_1))%*%c(delta_t_1))
  BF.alpha=exp(logL.alpha1-logL.alpha0)*(ratio.prior.delta1/ratio.prior.delta0)*(1/q_delta_j) # the ratio of the prior is very close to one!!!!
  if(BF.alpha>=1 || (BF.alpha<1 & runif(1)<BF.alpha)){
    alpha0=alpha1
    delta0=delta_t_1
    logL.alpha0=logL.alpha1
    accept.alpha=1
  }
  list(alpha0=alpha0,delta0=delta0,logL.alpha0=logL.alpha0,accept.alpha=accept.alpha)
}

#' @rdname BACps-internal

MH.alpha.s1_b=function(N,y,U,alpha0,alpha1,delta0,delta_t_1,q_delta_j){
  accept.alpha=0
  logL.alpha0=logL(N,y,U,alpha0,delta0)
  logL.alpha1=logL(N,y,U,alpha1,delta_t_1)
  
  #ratio.prior.delta0=dmvnorm(delta0,rep(0,length(delta0)),unit.info.prior.var_0)
  #ratio.prior.delta1=dmvnorm(delta_t_1,rep(0,length(delta_t_1)),unit.info.prior.var_1)
  ratio.prior.delta0=exp((-1/2)*c(delta0)%*%diag(1/10^10,nrow=length(delta0),ncol=length(delta0))%*%c(delta0))
  ratio.prior.delta1=exp(-1/2*c(delta_t_1)%*%diag(1/10^10,nrow=length(delta_t_1),ncol=length(delta_t_1))%*%c(delta_t_1))
  
  BF.alpha=exp(logL.alpha1-logL.alpha0)*(ratio.prior.delta1/ratio.prior.delta0)*(q_delta_j) # the ratio of the prior is very close to one!!!!
  if(BF.alpha>=1 || (BF.alpha<1 & runif(1)<BF.alpha)){
    alpha0=alpha1
    delta0=delta_t_1
    logL.alpha0=logL.alpha1
    accept.alpha=1
  }
  list(alpha0=alpha0,delta0=delta0,logL.alpha0=logL.alpha0,accept.alpha=accept.alpha)
}

####### MH intercept ##################
#' @rdname BACps-internal

MH.delta.intercept=function(N,y,U,alpha0,delta0){
  U.temp=matrix(rep(alpha0,N),nrow=N,byrow=TRUE)*as.matrix(U)
  fit.model.alpha0=summary(stats::glm(y~U.temp,family="binomial"))
  delta.mle.model.alpha0=fit.model.alpha0$coef[1,1]
  delta.cov.model.alpha0=fit.model.alpha0$cov.unscaled*10
  
  thres.prop=runif(1)   
  if(thres.prop>=0.5){
    delta.interc=MASS::mvrnorm(n=1,delta.mle.model.alpha0,delta.cov.model.alpha0[1,1]) # g(delta)=N(delta.mle, var-cov)
  } else { 
    delta.interc=delta0[1]+MASS::mvrnorm(n=1,0,10) #putting delta0 in the same length of alpha
    # delta1=delta0+error so g(delta'|delta*)=N(delta*,var-cov)
  } 
  logL.interc0=logL(N,y,U,alpha0,delta0)
  logL.interc1=logL(N,y,U,alpha0,c(delta.interc,delta0[-1]))
  ratio.prior.delta0=exp((-1/2)*c(delta0[1])%*%solve(unit.info.prior.var_0[1,1])%*%c(delta0[1]))
  ratio.prior.delta1=exp(-1/2*c(delta.interc)%*%solve(unit.info.prior.var_0[1,1])%*%c(delta.interc))
  if(thres.prop>=0.5){
    ratio.prop1=stats::dnorm(delta0[1],delta.mle.model.alpha0,delta.cov.model.alpha0[1,1])
    ratio.prop2=stats::dnorm(delta.interc,delta.mle.model.alpha0,delta.cov.model.alpha0[1,1])
    BF.interc=exp(logL.interc1-logL.interc0)*(ratio.prior.delta1/ratio.prior.delta0)*(ratio.prop1/ratio.prop2)
  } else{
    BF.interc=exp(logL.interc1-logL.interc0)*(ratio.prior.delta1/ratio.prior.delta0) # this is ratio logLik times ratio of prior
  } 
  if(BF.interc>=1 || (BF.interc<1 & runif(1)<BF.interc)){
    
    delta0[1]=delta.interc
  }
  list(delta0=delta0)
}

##############################
#' @rdname BACps-internal

countpattern=function (x, matching = FALSE) 
{
  nvar <- dim(x)[2]
  n <- dim(x)[1]
  b <- matrix(0, 2^nvar, nvar)
  for (i in 1:nvar) b[, nvar + 1 - i] <- rep(rep(c(0, 1), c(2^(i - 
                                                                 1), 2^(i - 1))), 2^(nvar - i))
  namespat <- b[, 1]
  for (i in 2:nvar) namespat <- paste(namespat, b[, i], sep = "")
  xpat <- x[, 1]
  for (i in 2:nvar) xpat <- 2 * xpat + x[, i]
  xpat <- xpat + 1
  pat <- tabulate(xpat, nbins = 2^nvar)
  names(pat) <- namespat
  if (matching) 
    return(list(pat = pat, matching = xpat))
  else return(pat)
}

###################2nd Stage ###########
###### Log Likelihood functions  #######


#' @rdname BACps-internal

logL_PS=function(x, U, alphaX,gamma){
  l_PS=t(x)%*%t(gamma[1]+(gamma[-1])%*%t(U[,alphaX!=0]))- sum(t(log(1+exp(gamma[1]+gamma[-1]%*%t(U[alphaX!=0]))))) # with intercept given by 1
  l_PS
}


#' @rdname BACps-internal

logL_Outcome=function(y,x,U, alphaX,gamma,Yparams){ 
  expit_temp = exp(gamma[1]+(gamma[-1]%*%t(U[,alphaX!=0])))
  ps= expit_temp/(1+expit_temp)
  qs_ps_hat_t=quantile(ps,c(.2,.4,.6,.8))
  strat_ps_hat_t=ifelse(ps<=qs_ps_hat_t[1],1,ifelse(ps<=qs_ps_hat_t[2],2,ifelse(ps<=qs_ps_hat_t[3],3,ifelse(ps<=qs_ps_hat_t[4],4,5))))
  h_ps=t(as.matrix(rbind(strat_ps_hat_t==2,strat_ps_hat_t==3,strat_ps_hat_t==4,strat_ps_hat_t==5)))
  Z=cbind(1,x,h_ps,U[,alphaX!=0])
  l_Y=sum(y*(as.matrix(Z)%*%Yparams))-sum(log(1+exp(as.matrix(Z)%*%Yparams))) 
  #l_Y=t(y)%*%(Yparams[1]+Yparams[2]*x+h_ps%*%Yparams[3:6]%*%+U[,alphaX!=0]%*%Yparams[-(1:6)])- sumt(log(1+exp(Yparams[1]+Yparams[2]*x+h_ps%*%Yparams[3:6]%*%+U[,alphaX!=0]%*%Yparams[-(1:6)]))) # with intercept given by 1
  l_Y
}

##################### Updating functions #####################
#' @rdname BACps-internal

update.params=function(y,x,U,gamma_mle,gamma_cov,outcome_mle,outcome_cov,nparams,ngamma_t,params_t,alphaX,likhood_t){
  #actualizo primero los gamma
  for(i in 1:length(params_t[1:ngamma_t])){
    oldparam <- params_t[1:ngamma_t][i]
    thres_1=runif(1)
    if(thres_1>=0.5){
      gamma_prop=MASS::mvrnorm(n=1,gamma_mle[i],gamma_cov[i,i]) # g(delta)=N(delta.mle, var-cov) , este vector como viene del fitted model solamente incluye las cov en el residual adj indicadas por alphaX
    } else { 
      gamma_prop=MASS::mvrnorm(n=1,oldparam,gamma_cov[i,i]) # delta1=delta0+error so g(delta'|delta*)=N(delta*,var-cov) esta es simetrica
    } 
    params_t[1:ngamma_t][i]=gamma_prop
    # To calculate acceptance probability
    newlikhood <- logL_PS(x,U,alphaX,params_t[1:ngamma_t])+logL_Outcome(y,x,U,alphaX,params_t[1:ngamma_t],params_t[-(1:ngamma_t)])
    if(thres_1>=0.5){
      num <- newlikhood + log(dnorm(params_t[1:ngamma_t][i],0,5))-log(dnorm(gamma_prop,gamma_mle[i],gamma_cov[i,i]))
      den <- likhood_t + log(dnorm(oldparam,0,5))-log(dnorm(oldparam,gamma_mle[i],gamma_cov[i,i]))
    }else{
      num <- newlikhood + log(dnorm(params_t[1:ngamma_t][i],0,5))
      den <- likhood_t + log(dnorm(oldparam,0,5))
    }
    A <- min(1,exp(num-den))
    # Accept/reject step:
    u <- runif(1)
    if (u <= A) { 
      likhood_t <- newlikhood
    }else{params_t[1:ngamma_t][i] <- oldparam }
  }
  #actualizo outcome
  for(i in 1:length(params_t[-(1:ngamma_t)])){
    oldparam <- params_t[-(1:ngamma_t)][i]
    #actualizo primero los gamma
    thres_2=runif(1) 
    if(thres_2>=0.5){
      outcome_prop=MASS::mvrnorm(n=1,outcome_mle[i],outcome_cov[i,i]) # g(delta)=N(delta.mle, var-cov) , este vector como viene del fitted model solamente incluye las cov en el residual adj indicadas por alphaX
    } else { 
      outcome_prop=MASS::mvrnorm(n=1,oldparam,outcome_cov[i,i]) # delta1=delta0+error so g(delta'|delta*)=N(delta*,var-cov)
    } 
    params_t[-(1:ngamma_t)][i]=outcome_prop
    # To calculate acceptance probability
    newlikhood <- logL_PS(x,U,alphaX,params_t[1:ngamma_t])+logL_Outcome(y,x,U,alphaX,params_t[1:ngamma_t],params_t[-(1:ngamma_t)])
    if(thres_2>=0.5){
      num <- newlikhood + log(dnorm(params_t[-(1:ngamma_t)][i],0,5))-log(dnorm(outcome_prop,outcome_mle[i],outcome_cov[i,i]))
      den <- likhood_t + log(dnorm(oldparam,0,5))-log(dnorm(oldparam,outcome_mle[i],outcome_cov[i,i])) 
    }else{
      num <- newlikhood + log(dnorm(params_t[-(1:ngamma_t)][i],0,5))
      den <- likhood_t + log(dnorm(oldparam,0,5))
    }
    A <- min(1,exp(num-den))
    # Accept/reject step:
    u <- runif(1)
    if (u <= A) { 
      likhood_t <- newlikhood
    }else {params_t[-(1:ngamma_t)][i] <- oldparam }
  }
  
  output <- c(params_t, likhood_t)
  output
}
######################
#' @rdname BACps-internal

update.params_when_alpha_0=function(y,x,U,nparams,ngamma_t,params_t,alphaX,likhood_t){
  for (i in 1:nparams) {
    oldparam <- params_t[i]
    params_t[i] <- runif(1, params_t[i]-.1, params_t[i]+.1) #uniform random walk MH
    # To calculate acceptance probability
    newlikhood <- logL_PS(x,U,alphaX,params_t[1:ngamma_t])+logL_Outcome(y,x,U,alphaX,params_t[1:ngamma_t],params_t[-(1:ngamma_t)])
    num <- newlikhood + log(dnorm(params_t[i],0,5))
    den <- likhood_t + log(dnorm(oldparam,0,5))
    A <- min(1,exp(num-den))
    # Accept/reject step:
    u <- runif(1)
    if (u <= A) { 
      likhood_t <- newlikhood
    }
    else { params_t[i] <- oldparam }
  }
  output <- c(params_t, likhood_t)
  output
}
######################

### Function for updating the model
#' @rdname BACps-internal

updatemodel=function(N,p,y,x,U,ngamma_t,params_t,alphaX_t,likhood_t,MX.nsims){
  
  repeat{
    alphaX_new=alphaX_t
    j=sample(p,1)
    alphaX_new[j]=1-alphaX_t[j]
    if(sum(alphaX_new)>0) break
  }
  
  index_alphaX_t=which(sapply(1:dim(MX.nsims)[1],function(i) identical(as.integer(alphaX_t),as.integer(MX.nsims[i,1:p])))==TRUE)
  index_alphaX_new=which(sapply(1:dim(MX.nsims)[1],function(i) identical(as.integer(alphaX_new),as.integer(MX.nsims[i,1:p])))==TRUE)  
  
  accept.alphaX=0
  # If covariate is not in the current model propose a candidate value
  if(alphaX_t[j]==0) {
    U.temp=matrix(rep(alphaX_new,N),nrow=N,byrow=TRUE)*as.matrix(U)
    gamma_mle_new=glm(x~U.temp,family="binomial")$coef
    phi_mle_new=glm(y~x+U.temp,family="binomial")$coef
    if( sum(alphaX_new)==0){
      gamma_j=rnorm(n=1,mean=gamma_mle_new[1],sd=.5)
      phi_j=0   #rnorm(n=1,mean=phi_mle_new[2],sd=.5)
    }else{
      gamma_j=rnorm(n=1,mean=gamma_mle_new[j+1],sd=.5)
      phi_j=rnorm(n=1,mean=phi_mle_new[j+2],sd=.5) 
    }
    gamma_temp=rep(NA,p+1)
    gamma_temp[1]=params_t[1:ngamma_t][1]
    gamma_temp[-1][alphaX_t!=0]=params_t[1:ngamma_t][-1]
    gamma_temp[j+1]=gamma_j
    gamma_new=gamma_temp[is.na(gamma_temp)==FALSE]
    ngamma_new=length(gamma_new)
    phi_temp=rep(NA,p)
    phi_temp[alphaX_t!=0]=params_t[-(1:ngamma_t)][-(1:6)]
    phi_temp[j]=phi_j
    phi_new=phi_temp[is.na(phi_temp)==FALSE]
    nphi_new=length(phi_new)
    params_new=c(gamma_new,params_t[-(1:ngamma_t)][(1:6)],phi_new)
  } else {
    gamma_temp=rep(NA,p+1)
    gamma_temp[1]=params_t[1:ngamma_t][1]
    gamma_temp[-1][alphaX_t!=0]=params_t[1:ngamma_t][-1]
    gamma_j=gamma_temp[-1][(alphaX_t-alphaX_new)!=0]
    gamma_temp[j+1]=NA
    gamma_new=gamma_temp[is.na(gamma_temp)==FALSE]
    ngamma_new=length(gamma_new)
    phi_temp=rep(NA,p)
    phi_temp[alphaX_t!=0]=params_t[-(1:ngamma_t)][-(1:6)]
    phi_j=phi_temp[(alphaX_t-alphaX_new)!=0]
    phi_temp[j]=NA
    phi_new=phi_temp[is.na(phi_temp)==FALSE]
    nphi_new=length(phi_new)
    params_new=c(gamma_new,params_t[-(1:ngamma_t)][(1:6)],phi_new)
    
  }
  # calculate the acceptance probability
  e_new=MX.nsims[index_alphaX_new,dim(MX.nsims)[2]]
  e_t=MX.nsims[index_alphaX_t,dim(MX.nsims)[2]]
  if(alphaX_t[j]==0) { 
    #newlikhood=logL_PS(alphaX_new,params_new[1:ngamma_new])-2*log(sum(alphaX_t)*N)*(e_t/e_new)+logL_Outcome(alphaX_new,params_new[1:ngamma_new],params_new[-(1:ngamma_new)])
    newlikhood=logL_PS(x,U,alphaX_new,params_new[1:ngamma_new])-2*log(N)*(e_t/e_new)+logL_Outcome(y,x,U,alphaX_new,params_new[1:ngamma_new],params_new[-(1:ngamma_new)])
    likhood_t=logL_PS(x, U,alphaX_t,params_t[1:ngamma_t])+logL_Outcome(y,x,U,alphaX_t,params_t[1:ngamma_t],params_t[-(1:ngamma_t)])
    
  }else{
    newlikhood=logL_PS(x,U,alphaX_new,params_new[1:ngamma_new])+logL_Outcome(y,x,U,alphaX_new,params_new[1:ngamma_new],params_new[-(1:ngamma_new)])
    likhood_t=logL_PS(x, U,alphaX_t,params_t[1:ngamma_t])-2*log(N)*(e_new/e_t)+logL_Outcome(y,x,U,alphaX_t,params_t[1:ngamma_t],params_t[-(1:ngamma_t)])
    
  }
  if(sum(alphaX_t)==0){
    num=newlikhood+log(mvtnorm::dmvnorm(gamma_new,rep(0,ngamma_new),diag(5,ngamma_new,ngamma_new)))+log(mvtnorm::dmvnorm(phi_new,rep(0,nphi_new),diag(5,nphi_new,nphi_new)))
    den=likhood_t+log(mvtnorm::dmvnorm(params_t[1:ngamma_t],rep(0,ngamma_t),diag(5,ngamma_t,ngamma_t)))
    
    num=num+log(MX.nsims[index_alphaX_new,dim(MX.nsims)[2]])
    den=den+log(MX.nsims[index_alphaX_t,dim(MX.nsims)[2]])+log(dnorm(gamma_j,0,.5))+log(dnorm(phi_j,0,.5))
    
  }else{
    if(alphaX_t[j]==0){
      num=newlikhood+log(mvtnorm::dmvnorm(gamma_new,rep(0,ngamma_new),diag(5,ngamma_new,ngamma_new)))+log(mvtnorm::dmvnorm(phi_new,rep(0,nphi_new),diag(5,nphi_new,nphi_new)))
      den=likhood_t+log(mvtnorm::dmvnorm(params_t[1:ngamma_t],rep(0,ngamma_t),diag(5,ngamma_t,ngamma_t)))+log(mvtnorm::dmvnorm(params_t[-(1:ngamma_t)][-(1:6)],rep(0,(ngamma_t-1)),diag(5,(ngamma_t-1),(ngamma_t-1))))
      
      num=num+log(MX.nsims[index_alphaX_new,dim(MX.nsims)[2]])
      den=den+log(MX.nsims[index_alphaX_t,dim(MX.nsims)[2]])+log(dnorm(gamma_j,mean=gamma_mle_new[j+1],sd=.5))+log(dnorm(phi_j,mean=phi_mle_new[j+2],sd=.5))
    } else {
      if(sum(alphaX_new)==0){
        num=newlikhood+log(mvtnorm::dmvnorm(gamma_new,rep(0,ngamma_new),diag(5,ngamma_new,ngamma_new)))
        den=likhood_t+log(mvtnorm::dmvnorm(params_t[1:ngamma_t],rep(0,ngamma_t),diag(5,ngamma_t,ngamma_t)))+log(mvtnorm::dmvnorm(params_t[-(1:ngamma_t)][-(1:6)],rep(0,(ngamma_t-1)),diag(5,(ngamma_t-1),(ngamma_t-1))))
        
        num=num+log(MX.nsims[index_alphaX_new,dim(MX.nsims)[2]])+log(dnorm(gamma_j,0,.5))
        den=den+log(MX.nsims[index_alphaX_t,dim(MX.nsims)[2]])
        
      }else{
        num=newlikhood+log(mvtnorm::dmvnorm(gamma_new,rep(0,ngamma_new),diag(5,ngamma_new,ngamma_new)))+log(mvtnorm::dmvnorm(phi_new,rep(0,nphi_new),diag(5,nphi_new,nphi_new)))
        den=likhood_t+log(mvtnorm::dmvnorm(params_t[1:ngamma_t],rep(0,ngamma_t),diag(5,ngamma_t,ngamma_t)))+log(mvtnorm::dmvnorm(params_t[-(1:ngamma_t)][-(1:6)],rep(0,(ngamma_t-1)),diag(5,(ngamma_t-1),(ngamma_t-1))))
        
        num=num+log(MX.nsims[index_alphaX_new,dim(MX.nsims)[2]])+log(dnorm(gamma_j,0,.5))+log(dnorm(phi_j,0,.5))
        den=den+log(MX.nsims[index_alphaX_t,dim(MX.nsims)[2]])
      }
    }
  }
  
  #Acceptance probability of RJ step:
  A <- min(1,exp(num-den))
  u <- runif(1)
  if (u <= A) { 
    likhood_t <- newlikhood
    alphaX_t=alphaX_new
    params_t=params_new
    accept.alphaX=1
  }
  
  output <- c(alphaX_t,params_t,likhood_t,accept.alphaX)
}
