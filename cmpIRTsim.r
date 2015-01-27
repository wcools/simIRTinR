#
# study to compare various IRT estimation methods
#
# for simulated data different R packages are used
# note that this paradigm allows for comparing various designs
# for example, different types of imbalanced data and sample sizes

# of interest are
library(ltm)
library(mirt)
library(MCMCpackd)
library(rstan)
library(reshape)

# ------------- #
# simulate data	#
# ------------- #
# see simIRT.r

# PROFICIENCY (subject abilities)
# set number of subjects and from uniform equidistant intervals generate z-scores to represent standard normal distribution
nrSubjects <- 256
S <- round(qnorm(seq(.01,(1-.01),length.out=nrSubjects),0,1),3)
names(S) <- paste0('s',sprintf(eval(paste0("%0",nchar(nrSubjects),"d")),1:length(S)))

# DIFFICULTY & DISCRIMINATION (item properties)
# set number of items and from uniform equidistant intervals generate scores to represent uniform distribution
nrItems <- 64									# must be dividable by 4 and not equal to 10, 100, 1000, ... for discrimination
# difficulties
I <- round(qlogis(seq(.05,.95,length.out=nrItems)),3)
names(I) <- paste0('i',sprintf(eval(paste0("%0",nchar(nrItems),"d")),1:length(I)))
# discrimination
D <- c(seq(.2,2,length.out=nrItems/4),seq(2,.2,length.out=nrItems/4),seq(2,.2,length.out=nrItems/4),seq(.2,2,length.out=nrItems/4))
names(D) <- names(I)

# SUCCESS PROBABILITIES
# assume that each subject responds to each item
# assume a success probability that is inversely related to the logit of S-I
# reset the discrimination D as rep(1,length(D)) to go from 2PL to 1PL
P2 <- round(plogis(sweep(outer(S,I,'-'),MARGIN=2,D,`*`)),3)
P1 <- round(plogis(sweep(outer(S,I,'-'),MARGIN=2,rep(1,length(D)),`*`)),3)
set.seed(312)
O1 <- matrix(unlist(Map(rbinom,prob=P1,size=1,n=1)),ncol=length(I))
O2 <- matrix(unlist(Map(rbinom,prob=P2,size=1,n=1)),ncol=length(I))
# alternative way to generate the 1 and 0's from a binomial
# O2 <- matrix(mapply(function(r,c) rbinom(1,1,P2[r,c]), row(P2), col(P2) ),nrow=nrow(P2),ncol=ncol(P2))
# O1 <- matrix(mapply(function(r,c) rbinom(1,1,P1[r,c]), row(P1), col(P1) ),nrow=nrow(P1),ncol=ncol(P1))
dimnames(P1) <- dimnames(P2) <- dimnames(O1) <- dimnames(O2) <- list(names(S),names(I))

# ------------- #
# IRT analyses	#
# ------------- #

sDta <- O2
lDta <- melt(sDta)
names(lDta) <- c("user","item","score")
lDta <- lDta[!is.na(lDta$score),]

# LTM and mirt
rltm <- ltm(sDta~z1)
cltm <- coef(rltm)
rmrt <- mirt(sDta,1,itemtype="2PL")
cmrt <- do.call(rbind,coef(rmrt)[1:length(I)])[,2:1]
dimnames(cmrt)[[1]] <- dimnames(cltm)[[1]]
cmrt[,1] <- -1*cmrt[,1]/cmrt[,2]
rmci <- MCMCirtKd(sDta,1,burnin=1000,mcmc=5000,thin=5,store.item=T,store.ability=T,verbose=2500)

# MCMCpack
cmci <- cmrt
cmci[,1] <- apply(rmci[,substring(dimnames(rmci)[[2]],1,2)=="al"],2,mean)/apply(rmci[,substring(dimnames(rmci)[[2]],1,2)=="be"],2,mean)
cmci[,2] <- 1.7*apply(rmci[,substring(dimnames(rmci)[[2]],1,2)=="be"],2,mean)
# note that difficulty in mcmc is centered! (so diff mcmc = diff mirt - mean diff mirt)

# ELO
runElo <- function(response,ability,difficulty,step){
 change <- step*(response-plogis(ability-difficulty))
 return(list(ability=ability+change,difficulty=difficulty-change,change=change))
}
n <- length(unique(lDta$user)); k <- length(unique(lDta$item))
n2 <- lDta[sample(1:nrow(lDta)),]

n2s <- vector("numeric",length=n); names(n2s) <- unique(lDta$user)
tmp <- aggregate(list(score=lDta$score),list(user=lDta$user),mean)
n2s[] <- qlogis(tmp[order(match(tmp[,1],names(n2s))),"score"])

n2i <- vector("numeric",length=k); names(n2i) <- unique(lDta$item)
tmp <- aggregate(list(score=lDta$score),list(item=lDta$item),mean)
n2i[] <- qlogis(tmp[order(match(tmp[,1],names(n2i))),"score"])

for(it in 1:nrow(n2)){
	obs <- runElo(n2[it,'score'], n2s[as.character(n2[it,'user'])], n2i[as.character(n2[it,'item'])], .5*(nrow(n2)-it+1)/(it^1.2+nrow(n2))+.01)
	n2s[as.character(n2[it,'user'])] <- obs$ability; n2i[as.character(n2[it,'item'])] <- obs$difficulty
}
celor <- cbind(diff=n2i)
aelor <- cbind(able=n2s)

# rstan
my2pl <- '
data {
  int<lower=1> J;                // number of students
  int<lower=1> K;                // number of questions
  int<lower=1> N;                // number of observations
  int<lower=1,upper=J> jj[N];    // student for observation n
  int<lower=1,upper=K> kk[N];    // question for observation n
  int<lower=0,upper=1> y[N];     // correctness of observation n
}
parameters {    
  real delta;                    // mean student ability
  real alpha[J];                 // ability for student j - mean ability
  real beta[K];                  // difficulty for question k
  real log_gamma[K];             // discriminativeness for question k
  real<lower=0> sigma_alpha;     // sd of student abilities  
  real<lower=0> sigma_beta;      // sd of question difficulties 
  real<lower=0> sigma_gamma;     // sd of log question discriminativeness
}
model {
  alpha ~ normal(0,sigma_alpha); 
  beta ~ normal(0,sigma_beta);   
  log_gamma ~ normal(0,sigma_gamma);
  delta ~ cauchy(0,5);
  sigma_alpha ~ cauchy(0,5);
  sigma_beta ~ cauchy(0,5);
  sigma_gamma ~ cauchy(0,5);
  for (n in 1:N)
    y[n] ~ bernoulli_logit( exp(log_gamma[kk[n]])
                            * (alpha[jj[n]] - beta[kk[n]] + delta) );
}
'
J <- length(unique(lDta$user))
K <- length(unique(lDta$item))
N <- nrow(lDta)

y <- lDta$score
jj <- as.numeric(substring(lDta$user,2,3))
kk <- as.numeric(substring(lDta$item,2,3))

dat = list(K=K,J=J,N=N,y=y,jj=jj,kk=kk)
fit0 <- stan(model_name="s1024i32pl2",model_code = my2pl, data=dat, iter = 10, verbose = FALSE) 
# rstn <- stan(fit = fit0, data = dat, iter = 8000, warmup=2000, chains = 4)

cstn <- cmrt
rstnx <- extract(rstn,permuted=T)
cstn[,1] <- colMeans(rstnx$beta)
cstn[,2] <- exp(colMeans(rstnx$log_gamma))

# plot
tmpPar <- par()

par(bg="#202020",col="grey",col.axis="grey",col.lab="grey",col.main="grey",col.sub="grey")

plot(seq(.05,.95,length.out=length(I)),cltm[,1],xlim=c(0.05,.95),type='n',xaxt='n',yaxt='n',ylim=c(-5,7),xlab="Items ~ P(Success)",ylab="Estimated Difficulty & Discrimination",main="IRT calibration in R for 2PL items")
lines(plogis(I),I,lwd=1.5)
points(plogis(I),I,cex=D)

# ltm
lines(seq(.05,.95,length.out=length(I)),cltm[,1],lwd=2,col=7)
points(seq(.05,.95,length.out=length(I)),cltm[,1],cex=cltm[,2],col=7)
# mirt
lines(seq(.05,.95,length.out=length(I)),cmrt[,1],col=2,lwd=2)
points(seq(.05,.95,length.out=length(I)),cmrt[,1],cex=cmrt[,2],col=2)
# MCMCpack
lines(seq(.05,.95,length.out=length(I)),cmci[,1],col=3,lwd=2)
points(seq(.05,.95,length.out=length(I)),cmci[,1],cex=cmci[,2],col=3)
# Elo
lines(seq(.05,.95,length.out=length(I)),celor,col=4,lwd=2)
# rstan
# lines(seq(.05,.95,length.out=length(I)),1.4*(cstn[,1]+mean(cltm[,1])),col=4,lwd=2)
# points(seq(.05,.95,length.out=length(I)),1.4*(cstn[,1]+mean(cltm[,1])),cex=cstn[,2],col=4)
# legend("topleft",c("ltm","mirt (diff=-1*diff/disc)","MCMCpack (disc=1.7*disc)","rstan (diff=1.4*(diff+.4))","elo"),lty=1,lwd=2,col=1:5,bty='n')
legend("topleft",c("ltm","mirt (diff=-1*diff/disc)","MCMCpack (disc=1.7*disc)","Elo rating"),lty=1,lwd=2,col=c(7,2:4),bty='n')
points(seq(.05,.95,length.out=length(I)),rep(-5,length(I)),pch=19,col=8,cex=D)
axis(1,names(I),at=plogis(I),cex.axis=.5,col="grey")
axis(2,at=seq(-5,5,1),cex.axis=1,las=1,col="grey")

par(tmpPar)

setwd("C:/Users/u0032822/Documents/code/GitHubRepos/simIRTinR")
savePlot("cmpIRTsim.png",type="png")
save.image("cmpIRTsim.RData")
