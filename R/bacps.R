#' Bayesian Adjustment for Confounding
#' @description  Main Algorithm that performs Bayesian Adjusting for confounding carrying out 1st and 2nd stage.
#' @param w a data.frame with at least the following variables: 
#' @param x a matrix with the exposure (binary) variable  
#' @param y a matrix with the response (binary) variable
#' @param U a data frame with the covariates
#' @param niter_1st iteration MCMC for the 1st stage
#' @param niter_2nd iteration MCMC for the 2nd stage
#' @param burnin_1st burn-in MCMC for the 2nd stage
#' @param burnin_2nd burn-in MCMC for the 2nd stage
#' @return a list  
#' @rdname bacps
#' @export 


bacps <- function(w,x,y,U, niter_1st, niter_2nd, burnin_1st, burnin_2nd){
  
  #################################
  ######## C++ Functions ##############
  ################################
  
  src_bin='
  double iter=Rcpp::as<double>(n);
  Rcpp::NumericVector xx=Rcpp::NumericVector (Rcpp::Dimension(iter,1,1));
  Rcpp::NumericMatrix mat(mx);
  double bin=1;
  int p = mat.ncol();
  int k = mat.nrow();
  for (int i=0; i<k;i++){
  double bin=1;
  for(int j=0; j<p;j++){
  if (mat(i,j)==1){
  bin += pow(2,j);
  }
  }
  xx[i]=bin;
  }
  return Rcpp::wrap(xx);
  '
  fool_bin=inline::cxxfunction(signature(n="integer",mx="numeric"),src_bin,plugin="Rcpp")
  ###########################
  
  src_postY='
  Rcpp::NumericMatrix index(mx);
  Rcpp::NumericMatrix binalpha(my);
  Rcpp::NumericMatrix final(mf);
  int k = index.nrow();
  int p = binalpha.ncol();
  int s = binalpha.nrow();
  for (int i=0; i<k;i++){
  for (int j=0; j<s;j++){
  if (index[i]==binalpha(j,0)){
  final(i,_)=binalpha.row(j);
  }
  }
  }
  return Rcpp::wrap(final);
  '
  fool_postY=inline::cxxfunction(signature(mx="numeric",my="numeric",mf="numeric"),src_postY,plugin="Rcpp")
  ##################
  
  src_Xnot_inY='
  Rcpp::NumericMatrix X(mx);
  Rcpp::NumericMatrix Y(my);
  Rcpp::NumericMatrix XnotY(mf);
  int m = X.nrow();
  int p = X.ncol();
  int k = Y.nrow();
  
  for (int i=0; i<m;i++){
  for (int j=0; j<k;j++){
  for (int s=0; s<p;s++){
  if (X(i,s)!=Y(j,s)){
  XnotY(i,j) +=1;
  }
  }
  
  }
  }
  return Rcpp::wrap(XnotY);
  '
  fool_xnoty=inline::cxxfunction(signature(mx="numeric",my="numeric",mf="numeric"),src_Xnot_inY,plugin="Rcpp")
  ################# end functions #####################################################################
  
  
  p = dim(U)[2]
  # constraint
  x0_y0=w/(1+w)
  x1_y0=1/(1+w)
  x0_y1=0.5
  
  AB <- matrix(0,2^p, p)
  for (i in 1:p) AB[, p + 1 - i] <- rep(rep(c(0, 1), c(2^(i - 1), 2^(i - 1))), 2^(p - i))
  AB=AB[-1,]
  
  A=matrix(NA,2^p-1,2^p-1)
  for(i in 1:((2^p)-1)){
    q=AB[i,]
    for(j in 1:(2^p-1)){
      notYnotX=0
      for(m in 1:p){
        if(AB[j,m]==0 & q[m]==0) notYnotX=notYnotX+1
      }
      A[i,j]=notYnotX
    }
  }
  
  B=matrix(NA,2^p-1,2^p-1)
  
  for(i in 1:(2^p-1)){
    q=AB[i,]
    for(j in 1:(2^p-1)){
      notYinX=0
      for(m in 1:p){
        if(AB[j,m]==0 & q[m]==1) notYinX=notYinX+1
      }
      B[i,j]=notYinX
    }
  }
  
  inY=apply(AB,1,sum)
  C=matrix(x0_y1,2^p-1,2^p-1)
  for(i in 1:(2^p-1)){
    C[i,]=C[i,]^inY
  }
  
  p <- dim(U)[2]
  N <- length(y)
  y0=y[x==0]
  U0=U[x==0,]
  N0=length(y0)
  
  expit=function(a){exp(a)/(1+exp(a))}
  
  full.y.model=stats::glm(y0~as.matrix(U0),family="binomial")
  delta.mle.full.model=full.y.model$coef
  cov.model.full=summary(full.y.model)$cov.unscaled*10
  ############################
  
  #####  initial values  #####
  alpha0= stats::rbinom(p,1,.5) #initial values alpha
  U.temp=matrix(rep(alpha0,N0),nrow=N0,byrow=TRUE)*as.matrix(U0)
  fit.model.alpha0=summary(glm(y0~U.temp,family="binomial"))
  delta.mle.model.alpha0=fit.model.alpha0$coef[,1]
  delta.cov.model.alpha0=fit.model.alpha0$cov.unscaled*10
  delta0=MASS::mvrnorm(n=1,delta.mle.model.alpha0,delta.cov.model.alpha0) #initial values delta
  delta_all=rep(NA,p+1)
  delta_all[1]=delta0[1]
  delta_all[-1][alpha0!=0]=delta0[-1]
  
  niter= niter_1st  #10000
  burnin= burnin_1st   #2000
  chain.alpha=matrix(NA,niter,p)
  chain.delta=matrix(NA,niter,p+1)
  chain.alpha[1,]=alpha0
  chain.delta[1,]=delta_all
  accept.alpha=0
  accept.delta=0
  k=1
  # Start the clock!
  ptm_rjmcmc <- proc.time()
  
  repeat {
    #step 1) update the parameters in the current model
    U.temp=matrix(rep(alpha0,N0),nrow=N0,byrow=TRUE)*as.matrix(U0)
    fit.model.alpha0=summary(stats::glm(y0~U.temp,family="binomial"))
    delta.mle.model.alpha0=fit.model.alpha0$coef[,1]
    delta.cov.model.alpha0=fit.model.alpha0$cov.unscaled*10
    delta=MH.delta(N,y,U,delta0,alpha0,delta.mle.model.alpha0,delta.cov.model.alpha0)
    delta0=delta$delta0
    #step 2.a) generate a proposed variable j  
    alpha1=alpha0
    j=sample(p,1)
    alpha1[j]=1-alpha0[j]
    delta_j=rnorm(n=1,mean=delta.mle.full.model[1+j],sd=cov.model.full[1+j,1+j])
    q_delta_j=dnorm(delta_j,mean=delta.mle.full.model[1+j],sd=cov.model.full[[1+j,1+j]])
    #unit information prior variance, equation 9 paper Dellaportas et al. 2003 
    #unit.info.prior.var_0=4*solve(t(cbind(rep(1,N[1]),as.matrix(U[,alpha0!=0])))%*%cbind(rep(1,N[1]),as.matrix(U[,alpha0!=0])))
    #unit.info.prior.var_1=4*solve(t(cbind(rep(1,N[1]),as.matrix(U[,alpha1!=0])))%*%cbind(rep(1,N[1]),as.matrix(U[,alpha1!=0])))
    
    #step 2.b)
    if(alpha0[j]==0){
      delta_t_1.temp=rep(NA,p+1)
      delta_t_1.temp[1]=delta0[1]
      delta_t_1.temp[-1][alpha0!=0]=delta0[-1]
      delta_t_1.temp[j+1]=delta_j
      delta_t_1=delta_t_1.temp[is.na(delta_t_1.temp)==FALSE]
      alpha=MH.alpha.s1_a(N,y,U,alpha0,alpha1,delta0,delta_t_1,q_delta_j) 
    }
    #step 2.c)
    else{
      delta_t_1.temp=rep(NA,p+1)
      delta_t_1.temp[1]=delta0[1]
      delta_t_1.temp[-1][alpha0!=0]=delta0[-1]
      delta_t_1.temp=delta_t_1.temp[-(1+j)]
      delta_t_1=delta_t_1.temp[is.na(delta_t_1.temp)==FALSE]
      alpha=MH.alpha.s1_b(N,y,U,alpha0,alpha1,delta0,delta_t_1,q_delta_j)
    }
    accept.alpha=accept.alpha+alpha$accept.alpha
    alpha0=alpha$alpha0
    delta0=alpha$delta0
    
    k=k+1
    chain.alpha[k,]=alpha0
    delta_all_0=rep(NA,p+1)
    delta_all_0[1]=delta0[1]
    delta_all_0[-1][alpha0!=0]=delta0[-1]
    chain.delta[k,]=delta_all_0
    if(k==niter) break
  }
  
  
  #acceptation ratio of parameter alpha
  print('acceptation ratio of parameter alpha:')
  print(accept.alpha/niter)
  
  #posterior distribution of alphaY
  post.alphaY=countpattern(chain.alpha[burnin:niter,])/(niter-burnin+1)
  
  post.alphaX_Y_w=(x0_y0^A*x1_y0^B*C)%*%as.matrix(cbind(post.alphaY))[-1]
  post.alphaX_Y_w=post.alphaX_Y_w/sum(post.alphaX_Y_w)
  
  MX_w=cbind(AB,post.alphaX_Y_w)
  colnames(MX_w)[p+1]="post.alphaX|Y"
  
  print('First Stage Done!')
  print(proc.time() - ptm_rjmcmc)
  
  ######     2nd Stage    ##########
  
  src_indexMX='
  Rcpp::NumericVector alpha(mx);
  Rcpp::NumericMatrix MX(my);
  int m=MX.nrow();
  int s=0;
  Rcpp::NumericVector alpha1=Rcpp::NumericVector (Rcpp::Dimension(1,25,1));
  for (int i=0; i<m;i++){
  alpha1=MX.row(i);
  if( std::equal(alpha.begin(), alpha.end(), alpha1.begin()) ){
  s=1;
  }
  }
  return Rcpp::wrap(s);
  '
  fool_indexMX=inline::cxxfunction(signature(mx="numeric",my="numeric"),src_indexMX,plugin="Rcpp")
  
  src_indexmx='
  Rcpp::NumericVector alpha(mx);
  Rcpp::NumericMatrix MX(my);
  int m=MX.nrow();
  Rcpp::NumericVector alpha1=Rcpp::NumericVector (Rcpp::Dimension(1,25,1));
  for (int i=0; i<m;i++){
  alpha1=MX.row(i);
  if( std::equal(alpha.begin(), alpha.end(), alpha1.begin()) ){
  return Rcpp::wrap(i);
  }
  }
  '
  fool_indexmx=inline::cxxfunction(signature(mx="numeric",my="numeric"),src_indexmx,plugin="Rcpp")
  
  src_bin='
  double iter=Rcpp::as<double>(n);
  Rcpp::NumericVector xx=Rcpp::NumericVector (Rcpp::Dimension(iter,1,1));
  Rcpp::NumericMatrix mat(mx);
  int p = mat.ncol();
  int k = mat.nrow();
  for (int i=0; i<k;i++){
  double bin=1;
  for(int j=0; j<p;j++){
  if (mat(i,j)==1){
  bin += pow(2,j);
  }
  }
  xx[i]=bin;
  }
  return Rcpp::wrap(xx);
  '
  fool_bin=inline::cxxfunction(signature(n="integer",mx="numeric"),src_bin,plugin="Rcpp")
  
  src_postY='
  Rcpp::NumericMatrix index(mx);
  Rcpp::NumericMatrix binalpha(my);
  Rcpp::NumericMatrix final(mf);
  int k = index.nrow();
  int s = binalpha.nrow();
  for (int i=0; i<k;i++){
  for (int j=0; j<s;j++){
  if (index[i]==binalpha(j,0)){
  final(i,_)=binalpha.row(j);
  }
  }
  }
  return Rcpp::wrap(final);
  '
  fool_postY=inline::cxxfunction(signature(mx="numeric",my="numeric",mf="numeric"),src_postY,plugin="Rcpp")
  ######  end C++ functions ####################
  ncov=p
  probmodel=array(0,2^ncov)         # to be used to compute the posterior probability for each model
  model=array(0,dim=c(2^ncov,ncov)) # model identifier
  nparams=p+1+2+4+p# total parameters
  niter= niter_2nd #150000
  burnin= burnin_2nd #10000
  #########################################################
  # Read priors and proposals parameter values for MH and RJ steps
  
  # Normal priors on regression coefficients
  
  mu_logitps=rep(0,p+1)
  mu_logity=rep(0,2+4+p)
  sig2_logitps=rep(10,p+1)
  sig2_logity=rep(0,2+4+p)
  
  # Parameters for M-H updates 
  
  # Parameters for RJ updates for gamma and phi
  rjmu=c(rep(0,p),rep(0,p))
  rjsig2=c(rep(.5,p),rep(.5,p))
  rjsig=sqrt(rjsig2)
  
  # Prior Model (it would be given by MX)
  #prob_alphaX
  
  #############################################
  
  # Set initial parameter values and calculate log-likelihood
  probparam=array(0,nparams) # it will contain posterior marginal prob that each parameter is in the model
  sum_gamma=array(0,p+1)
  sum_gamma2=array(0,p+1)# model averaged posterior mean
  sum_y=array(0,2+4+p) # standard deviation
  sum_y2=array(0,2+4+p)
  
  ptm_rjmcmc <- proc.time()
  
  MX=MX_w
  #### alphaX model  
  alphaX_t=as.numeric(MX[sample(dim(MX)[1],1),-(p+1)])
  U_t=matrix(rep(alphaX_t,N),nrow=N,byrow=TRUE)*as.matrix(U)
  fit.PS.init=stats::glm(x~U_t,family="binomial")
  ps_hat_t=fit.PS.init$fitted.values
  qs_ps_hat_t=quantile(ps_hat_t,c(.2,.4,.6,.8))
  strat_ps_hat_t=ifelse(ps_hat_t<=qs_ps_hat_t[1],1,ifelse(ps_hat_t<=qs_ps_hat_t[2],2,ifelse(ps_hat_t<=qs_ps_hat_t[3],3,ifelse(ps_hat_t<=qs_ps_hat_t[4],4,5))))
  h_ps_hat_t=t(as.matrix(rbind(strat_ps_hat_t==2,strat_ps_hat_t==3,strat_ps_hat_t==4,strat_ps_hat_t==5)))
  Z_t=cbind(x,h_ps_hat_t,U_t) #Z=[1,X,quantiles PS, U]
  fit.Y.init=stats::glm(y~Z_t,family="binomial")
  
  gamma_t=summary(fit.PS.init)$coef[,1] # initial values of gamma
  gamma.cov.model_t=summary(fit.PS.init)$cov.unscaled*10
  outcome.params_t=summary(fit.Y.init)$coef[,1]  # initial values of outcome params
  outcome.cov_t=summary(fit.Y.init)$cov.unscaled*10
  
  gamma_all.init=rep(NA,p+1)
  gamma_all.init[1]=gamma_t[1]
  gamma_all.init[-1][alphaX_t!=0]=gamma_t[-1]
  outcome.params_all.init=rep(NA,2+4+p)
  outcome.params_all.init[1:6]=outcome.params_t[1:6]
  outcome.params_all.init[-(1:6)][alphaX_t!=0]=outcome.params_t[-(1:6)]
  
  lik.ps.init=logL_PS(x,U, alphaX_t,gamma_t)
  lik.y.init=logL_Outcome(y,x,U,alphaX_t,gamma_t,outcome.params_t)
  likhood_t=lik.ps.init+lik.y.init
  
  params_t=c(gamma_t,outcome.params_t)
  ngamma_t=length(gamma_t)
  nparams_t=length(params_t)
  
  accept.alphaX=0
  k=1
  chain.alphaX=matrix(NA,niter,p)
  chain.gamma=matrix(NA,niter,p+1) #plus intercept
  chain.outcome=matrix(NA,niter,(2+4+p)) # (intercept and betaX)
  chain.alphaX[1,]=alphaX_t
  chain.gamma[1,]=gamma_all.init
  chain.outcome[1,]=outcome.params_all.init
  ######## RJ MCMC 
  repeat{
    if(sum(alphaX_t)>0){
      U_t=matrix(rep(alphaX_t,N),nrow=N,byrow=TRUE)*as.matrix(U)
      fit.PS.init=stats::glm(x~U_t,family="binomial")
      ps_hat_t=fit.PS.init$fitted.values
      qs_ps_hat_t=quantile(ps_hat_t,c(.2,.4,.6,.8))
      strat_ps_hat_t=ifelse(ps_hat_t<=qs_ps_hat_t[1],1,ifelse(ps_hat_t<=qs_ps_hat_t[2],2,ifelse(ps_hat_t<=qs_ps_hat_t[3],3,ifelse(ps_hat_t<=qs_ps_hat_t[4],4,5))))
      h_ps_hat_t=t(as.matrix(rbind(strat_ps_hat_t==2,strat_ps_hat_t==3,strat_ps_hat_t==4,strat_ps_hat_t==5)))
      Z_t=cbind(x,h_ps_hat_t,U_t) #Z=[1,X,quantiles PS, U]
      fit.Y.init=stats::glm(y~Z_t,family="binomial")
      gamma_mle=summary(fit.PS.init)$coef[,1] 
      gamma_cov=summary(fit.PS.init)$cov.unscaled*10
      outcome_mle=summary(fit.Y.init)$coef[,1]  
      outcome_cov=summary(fit.Y.init)$cov.unscaled*10
      
      # Update the parameters presents in the model
      output_1=update.params(y,x,U, gamma_mle,gamma_cov,outcome_mle,outcome_cov,nparams_t,ngamma_t,params_t,alphaX_t,likhood_t)
    }else{
      output_1=update.params_when_alpha_0(y,x,U,nparams_t,ngamma_t,params_t,alphaX_t,likhood_t)
    }
    params_t=output_1[1:nparams_t]
    likhood_t=output_1[nparams_t+1]
    #aca itero y me quedo con la ultima iteracion
    
    # Update the covariates in the model
    output_2=updatemodel(N,p,y,x,U,ngamma_t,params_t,alphaX_t,likhood_t,MX)
    noutput=length(output_2)
    alphaX_t=output_2[1:p]
    params_t=output_2[(p+1):(noutput-2)]
    likhood_t=output_2[(noutput-1)]
    accept.alphaX=accept.alphaX+output_2[noutput]
    ngamma_t=sum(alphaX_t)+1
    nparams_t=length(params_t)
    
    #Record the parameters
    k=k+1
    chain.alphaX[k,]=alphaX_t
    gamma_all_t=rep(NA,p+1)
    gamma_all_t[1]=params_t[1]
    gamma_all_t[-1][alphaX_t!=0]=params_t[1:ngamma_t][-1]
    outcome.params_all_t=rep(NA,2+4+p)
    outcome.params_all_t[1:6]=params_t[-(1:ngamma_t)][1:6]
    outcome.params_all_t[-(1:6)][alphaX_t!=0]=params_t[-(1:ngamma_t)][-(1:6)]
    chain.gamma[k,]=gamma_all_t
    chain.outcome[k,]=outcome.params_all_t
    
    if(k==niter) break
  }
  
  #acceptation ratio of parameters 
  print(accept.alphaX/niter)
  
  #Posterior Model probability distribution of alphaX
  
  BIN=fool_bin(niter-burnin+1,chain.alphaX[burnin:niter,])
  models=as.matrix(as.numeric(names(table(BIN))))
  temp=matrix(0,length(models),p+1)
  ident_each_selectedmodel=fool_postY(models,cbind(BIN,chain.alphaX[burnin:niter,]),temp)
  post.alphaX=table(BIN)/(niter-burnin+1)
  
  stage2_models=cbind(chain.outcome[burnin:niter,],chain.gamma[burnin:niter,],chain.alphaX[burnin:niter,])
  stage2_post_alphaX=cbind(ident_each_selectedmodel,post.alphaX) 
  
  print('2nd Stage Done!')
  print(proc.time() - ptm_rjmcmc)
  
  #####################################   
  niter_burnin=niter-burnin+1
  nY=p+2+4 # covariates + intercept + exposure + quintiles PS 
  nPS=1+p
  chains=stage2_models
  k=1
  triangle=matrix(NA,niter_burnin,1)
  for(k in 1:niter_burnin){
    chains[k,][is.na(chains[k,])]=0
    alphaX=chains[k,-(1:(nY+nPS))]
    gamma=chains[k,(nY+1):(nY+nPS)]
    ps_hat_t=expit(gamma[1]+(gamma[-1][alphaX!=0]%*%t(U[,alphaX!=0])))
    qs_ps_hat_t=quantile(ps_hat_t,c(.2,.4,.6,.8))
    strat_ps_hat_t=ifelse(ps_hat_t<=qs_ps_hat_t[1],1,ifelse(ps_hat_t<=qs_ps_hat_t[2],2,ifelse(ps_hat_t<=qs_ps_hat_t[3],3,ifelse(ps_hat_t<=qs_ps_hat_t[4],4,5))))
    h_ps=t(as.matrix(rbind(strat_ps_hat_t==2,strat_ps_hat_t==3,strat_ps_hat_t==4,strat_ps_hat_t==5)))
    Z=cbind(1,x,h_ps,U[,alphaX!=0])
    Yparams=c(chains[k,1:(2+4)],chains[k,7:(p+6)][alphaX!=0])
    Prob_y_x_c=expit(as.matrix(Z)%*%Yparams)
    triangle[k]=(sum(Prob_y_x_c[x==1])/sum(x==1))-(sum(Prob_y_x_c[x==0])/sum(x==0))
  }
  
  sum_triangle_alpha=matrix(NA,dim(stage2_post_alphaX)[1],1)
  sigma2_triangle=matrix(NA,dim(stage2_post_alphaX)[1],1)
  l=1
  for(l in 1:(dim(stage2_post_alphaX)[1])){
    pp1=sapply(1:niter_burnin,function(i) identical(as.integer(stage2_post_alphaX[l,2:(p+1)]),as.integer(stage2_models[i,(nY+nPS+1):ncol(stage2_models)])))==TRUE
    
    sum_triangle_alpha[l]=sum(triangle[pp1])/sum(pp1)
    sigma2_triangle[l]=sum((triangle[pp1]-sum_triangle_alpha[l])^2)/sum(pp1)
  }
  #mse=sum(stage2_post_alphaX[,p+2]*(sum_triangle_alpha-.1)^2)
  ace=sum(sum_triangle_alpha*stage2_post_alphaX[,p+2])
  
  return( list(MX_w, post.alphaY,stage2_models, stage2_post_alphaX, ace))
  
  
}






