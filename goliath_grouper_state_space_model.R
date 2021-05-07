rm(list=ls())
library(rstan);library(loo);library(bayesplot);library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
library(here)
here()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###Functions####
ts_reef = function(X){
    occ_by_year<- X %>% group_by(year) %>% summarize(n.occ=sum(occ),n.surv=n(),p.occ=n.occ/n.surv,sp=sp$commonname[match(unique(sp_x$speciesid),sp$speciesid)])
    abun_by_year<- X %>% group_by(year,geogr) %>% summarize(site_abun=mean(abund_trans),sd_abund=sd(abund_trans),n.surv=n()) %>% group_by(year) %>% summarize(mean_abund=mean(site_abun),n.survs=sum(n.surv),n.sites=n())
    total_sd<- X %>% group_by(year) %>% summarize(sd=sd(abund_trans))
    
    comb<- left_join(occ_by_year,abun_by_year)
    comb2<- left_join(comb,total_sd)
    ts_dat<- comb2
  return(ts_dat)
}

abund_tranfs<- function(probs){
  sum<- probs[1]*0+probs[2]*1+probs[3]*2+probs[4]*11+probs[5]*101
  return(sum)
}

ord_to_n<- function(x,c){
  abund_x<- numeric(length(x))
  p= matrix(ncol=5,nrow=length(x))
  if(ncol(c)==2){
      p[,1]=plogis(c[,1]-x)
      p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
      p[,3]=1-plogis(c[,2]-x)
      p[,4]=0
      p[,5]=0
      for(i in 1:length(x)){
        abund_x[i]=abund_tranfs(p[i,])  
      }
  }
  if(ncol(c)==3){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
    p[,4]=1-plogis(c[,3]-x)
    p[,5]=0
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  if(ncol(c)==4){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
    p[,4]=plogis(c[,4]-x[,i])-plogis(c[,3]-x[,i])
    p[,5]=1-plogis(c[,4]-x[,i])
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  return(abund_x)
}

TS_stan_state_only_plot_MARSS<- function(sp,GZ,params1,TT,ts){
  pdf(paste('./figures/',paste(sp,GZ,'state_space',sep='_'),'.pdf',sep=''),width=8,height=6)
  lambda_mat<- list()  
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,1200))
    
    if(ncol(params1$c)==2){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
        reef_coef[,3]<-1-plogis(params1$c[,2]-params1$a_yr[,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
        reef_coef[,8]<-1-plogis(params1$c[,2]-params1$x[,i])
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0

      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1$c)==3){

        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,i])-plogis(params1$c[,2]-params1$a_yr[,i])
        reef_coef[,4]<- 1-plogis(params1$c[,3]-params1$a_yr[,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i]-params1$a)-plogis(params1$c[,1]-params1$x[,i]-params1$a)
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i]-params1$a)-plogis(params1$c[,2]-params1$x[,i]-params1$a)
        reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,i]-params1$a)
        reef_coef[,10]<- 0
    
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1$c)==4){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr2[,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr2[,i])-plogis(params1$c[,1]-params1$a_yr2[,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr2[,i])-plogis(params1$c[,2]-params1$a_yr2[,i])
        reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr2[,i])-plogis(params1$c[,3]-params1$a_yr2[,i])
        reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr2[,i])
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i]-params1$a)
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i]-params1$a)-plogis(params1$c[,1]-params1$x[,i]-params1$a)
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i]-params1$a)-plogis(params1$c[,2]-params1$x[,i]-params1$a)
        reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,i]-params1$a)-plogis(params1$c[,3]-params1$x[,i]-params1$a)
        reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,i]-params1$a)
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    lambda_mat[[i]]=reef_coef
  }  
  
  
  x_mat<- data.frame(median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:TT){
    x_mat[i,1]=median(lambda_mat[[i]]$lambda.x)
    x_mat[i,2]=quantile(lambda_mat[[i]]$lambda.x,0.05)
    x_mat[i,3]=quantile(lambda_mat[[i]]$lambda.x,0.95)
  }
  
  y_mat<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  
  for(i in 1:TT){
    y_mat[i,2]=median(lambda_mat[[i]]$lambda.y)
    y_mat[i,3]=quantile(lambda_mat[[i]]$lambda.y,0.05)
    y_mat[i,4]=quantile(lambda_mat[[i]]$lambda.y,0.95)
  }
  
  par(xpd=T)
  plot(y_mat$median.reef~y_mat$year,type='n',ylim=c(min(x_mat),max(c(max(y_mat[,2]),max(x_mat)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=paste(sp,GZ,sep=', '))
  lines(x_mat[,1]~y_mat$year,lty=5,lwd=2,col='darkcyan')
  x<- c(y_mat$year, rev(y_mat$year))
  y1<- c(x_mat[,2], rev(x_mat[,3]))
  polygon(x, y1, col = adjustcolor('darkcyan', alpha = 0.1), border=NA) # Add uncertainty polygon

  lines(y_mat$median.reef~y_mat$year,col='dodgerblue4',lwd=2)
  points(y_mat$median.reef~y_mat$year,col='white',pch=21,bg='dodgerblue4',cex=1.5)
  text(y=rep(0,nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
  dev.off()
}


###REEF Data Filtering####
####Stan model - ordinal abundance####
SS_trend_ord<-"functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
  int<lower=1> N;//number of observations (REEF surveys)
  int y[N]; //abundance category for each survey
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix X
  int K; //ordinal levels
  int TT; // timespan
  int<lower=0> N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year
}
parameters {
  ordered[K-1] c; //cutpoints
  real x0; //initial popn size

  //deviations from intercept
  vector[Z] beta; //effort coefficients - RVC
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //variance on the deviance components
  real<lower = 0> sd_site;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  real<lower = 0> sd_r;
  real<lower = 0> sd_q;
  
  vector[TT] pro_dev; //process deviations
  vector[N_yr] obs_dev; //observation deviations 
}

transformed parameters{
  vector[TT] x;
  vector[N_yr] a_yr;

  x[1] = x0 + pro_dev[1];
   
  for(t in 2:TT){
    x[t] = x[t-1] + pro_dev[t];
  }
   
  for(i in 1:N_yr){
    a_yr[i] = x[yr_index[i]] + obs_dev[i]; 
  }
}  

model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta ~ normal(0,2); //covariates
  x0 ~ normal(0,5); //initial state

  //variance terms
  sd_q ~inv_gamma(2,0.25);
  sd_r ~ inv_gamma(2,0.25);
  sd_site ~ inv_gamma(2, 1);
  sd_dv ~ inv_gamma(3, 1);
  sd_dmy ~ inv_gamma(2, 1);
  
  //varying intercepts
  a_site ~ normal(0, sd_site);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);
  
  for(t in 1:TT){
    pro_dev[t] ~ normal(0, sd_q);
    obs_dev[t] ~ normal(0,sd_r);
  }

  y ~ ordered_logistic(a_yr[year_id]+a_site[site]+a_dv[diver]+a_dmy[dmy]+X*beta,c);
  
}
"
###Florida (All Regions) Population Trajectory####
gg_occs<- read.csv(here('data','Goliath_Grouper_surveys_fl.csv'))
gg_ts<- ts_reef(gg_occs)

X<- matrix(data=c(scale(as.numeric(gg_occs$btime)),scale(as.numeric(gg_occs$averagedepth)),scale(as.numeric(gg_occs$visibility)),scale(as.numeric(gg_occs$current)),gg_occs$exp_binary),ncol=5,nrow=nrow(gg_occs))

gg_SS<- rstan::stan(model_code = SS_trend_ord, data = list(y =gg_occs$abundance2,
                                                                N = nrow(gg_occs),
                                                                site=as.numeric(factor(gg_occs$geogr)),
                                                                N_site=length(unique(gg_occs$geogr)),
                                                                diver=as.numeric(factor(gg_occs$fish_memberid)),
                                                                N_dv=length(unique(gg_occs$fish_memberid)),
                                                                dmy=as.numeric(factor(gg_occs$site_dmy)),
                                                                N_dmy=length(unique(gg_occs$site_dmy)),
                                                                K=length(unique(gg_occs$abundance)),
                                                                X=X,
                                                                Z=ncol(X),
                                                                TT=27,
                                                                N_yr=length(unique(gg_occs$year)),
                                                                yr_index=sort(unique(as.numeric(factor(gg_occs$year)))),
                                                                year_id=as.numeric(factor(gg_occs$year))),
                   pars = c('c','sd_site','sd_dv','sd_dmy','sd_r','sd_q','x','a_yr','beta'),
                   control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 4, iter = 450, thin = 1)


gg_params<- rstan::extract(gg_SS)

diff=(ord_to_n(gg_params$x[,27],gg_params$c)/ord_to_n(gg_params$x[,17],gg_params$c))
median(diff)-1
quantile(diff,0.05)-1
quantile(diff,0.95)-1

TS_stan_state_only_plot_MARSS(sp='Goliath Grouper',GZ='Florida (All Regions)',params1=gg_params,TT=27,ts=gg_ts)
dev.off()

##Florida Keys to Dry Tortugas trajectory #####
gg_occs_fk<- gg_occs[which(substring(gg_occs$geogr4,1,2)==34),]
gg_ts_fk<- ts_reef(gg_occs_fk,sp=goliath)
length(unique(gg_occs_fk$site))

##FK model
X<- matrix(data=c(scale(as.numeric(gg_occs_fk$btime)),scale(as.numeric(gg_occs_fk$averagedepth)),scale(as.numeric(gg_occs_fk$visibility)),scale(as.numeric(gg_occs_fk$current)),gg_occs_fk$exp_binary),ncol=5,nrow=nrow(gg_occs_fk))

gg_SS_fk<- rstan::stan(model_code = SS_trend_ord, data = list(y =gg_occs_fk$abundance2,
                                                           N = nrow(gg_occs_fk),
                                                           site=as.numeric(factor(gg_occs_fk$geogr)),
                                                           N_site=length(unique(gg_occs_fk$geogr)),
                                                           diver=as.numeric(factor(gg_occs_fk$fish_memberid)),
                                                           N_dv=length(unique(gg_occs_fk$fish_memberid)),
                                                           dmy=as.numeric(factor(gg_occs_fk$site_dmy)),
                                                           N_dmy=length(unique(gg_occs_fk$site_dmy)),
                                                           K=length(unique(gg_occs_fk$abundance)),
                                                           X=X,
                                                           Z=ncol(X),
                                                           TT=27,
                                                           N_yr=length(unique(gg_occs_fk$year)),
                                                           yr_index=sort(unique(as.numeric(factor(gg_occs_fk$year)))),
                                                           year_id=as.numeric(factor(gg_occs_fk$year))),
                    pars = c('c','sd_site','sd_dv','sd_dmy','sd_r','sd_q','x','a_yr'),
                    control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 4, iter = 450, thin = 1)

gg_params_fk<- rstan::extract(gg_SS_fk)

TS_stan_state_only_plot_MARSS(sp='Goliath Grouper',GZ='Florida Keys',params1=gg_params_fk,TT=27,ts=gg_ts_fk)
dev.off()

diff=(ord_to_n(gg_params_fk$x[,27],gg_params_fk$c)/ord_to_n(gg_params_fk$x[,17],gg_params_fk$c))
median(diff)-1
quantile(diff,0.05)-1
quantile(diff,0.95)-1

##East Coast population trajectory ####
gg_occs_ec<- gg_occs[which(substring(gg_occs$geogr4,1,2)==33),]
gg_ts_ec<- ts_reef(gg_occs_ec,sp=goliath)
length(unique(gg_occs_ec$site))

X<- matrix(data=c(scale(as.numeric(gg_occs_ec$btime)),scale(as.numeric(gg_occs_ec$averagedepth)),scale(as.numeric(gg_occs_ec$visibility)),scale(as.numeric(gg_occs_ec$current)),gg_occs_ec$exp_binary),ncol=5,nrow=nrow(gg_occs_ec))

gg_SS_ec<- rstan::stan(model_code = SS_trend_ord, data = list(y =gg_occs_ec$abundance2,
                                                              N = nrow(gg_occs_ec),
                                                              site=as.numeric(factor(gg_occs_ec$geogr)),
                                                              N_site=length(unique(gg_occs_ec$geogr)),
                                                              diver=as.numeric(factor(gg_occs_ec$fish_memberid)),
                                                              N_dv=length(unique(gg_occs_ec$fish_memberid)),
                                                              dmy=as.numeric(factor(gg_occs_ec$site_dmy)),
                                                              N_dmy=length(unique(gg_occs_ec$site_dmy)),
                                                              K=length(unique(gg_occs_ec$abundance)),
                                                              X=X,
                                                              Z=ncol(X),
                                                              TT=27,
                                                              N_yr=length(unique(gg_occs_ec$year)),
                                                              yr_index=sort(unique(as.numeric(factor(gg_occs_ec$year)))),
                                                              year_id=as.numeric(factor(gg_occs_ec$year))),
                       pars = c('c','sd_site','sd_dv','sd_dmy','sd_r','sd_q','x','a_yr','beta'),
                       control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 4, iter = 450, thin = 1)

gg_params_ec<- rstan::extract(gg_SS_ec)

diff=(ord_to_n(gg_params_ec$x[,27],gg_params_ec$c)/ord_to_n(gg_params_ec$x[,17],gg_params_ec$c))
median(diff)-1
quantile(diff,0.05)-1
quantile(diff,0.95)-1

TS_stan_state_only_plot_MARSS(sp='Goliath Grouper',GZ='Florida (East Coast)',params1=gg_params_ec,TT=27,ts=gg_ts_ec)
dev.off()

##Gulf of Mexico - data summary###
##East Coast
gg_occs_mx<- gg_occs[which(substring(gg_occs$geogr4,1,1)==2),]
gg_occs_mx<- subset(gg_occs_mx,year>=1998)
gg_ts_mx<- ts_reef(gg_occs_ec,sp=goliath)
length(unique(gg_occs_mx$site))
