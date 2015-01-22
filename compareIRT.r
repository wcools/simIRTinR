# study to compare various IRT estimation methods

# of interest are
	# ltm
	# mirt
	# MCMCpack
	# Elo
	# rstan

# 1. simulation of data
# using a grid of population values, data are generated with
# seeds 101:120
# number of subjects 50, 100, 200, 400, 800, 1600
# number of items 33, 65, 129	-> consider abilities
# patterns of missingness: adaptive testing, blocking, random (only for 33 items)
 
# the data are stord in an R workspace: simirt.RData in the following working directory
setwd("C:\\Users\\u0032822\\Documents\\code\\IRTestimation")

# use the population matrix for 50 subject proficiencies (normally distributed) that respond to 65 items with differing difficulty / discrimination (uniformly distributed)
load("sim_1pl_1024_0032.RData")
load("sim_2pl_1024_0032.RData")
ls()
popMx1 <- sim_1pl_1024_0032_01
popMx2 <- sim_2pl_1024_0032_01
dimnames(popMx1) <- list(paste0("s",1:nrow(popMx1)),paste0("i",1:ncol(popMx1)))
dimnames(popMx2) <- list(paste0("s",1:nrow(popMx2)),paste0("i",1:ncol(popMx2)))
popLf1 <- melt(popMx1)
popLf2 <- melt(popMx2)
names(popLf1) <- names(popLf2) <- c("user","item","score")

# LTM and mirt
rltm <- ltm(popMx2~z1)
cltm <- coef(rltm)
rmrt <- mirt(popMx2,1,itemtype="2PL")
cmrt <- do.call(rbind,coef(rmrt)[1:32])[,2:1]
dimnames(cmrt)[[1]] <- dimnames(cltm)[[1]]
cmrt[,1] <- -1*cmrt[,1]/cmrt[,2]
rmci <- MCMCirtKd(popMx2,1,burnin=1000,mcmc=5000,thin=5,store.item=T,store.ability=T,verbose=2500)

# MCMCpack
cmci <- cmrt
cmci[,1] <- -1*apply(rmci[,substring(dimnames(rmci)[[2]],1,2)=="al"],2,mean)/apply(rmci[,substring(dimnames(rmci)[[2]],1,2)=="be"],2,mean)
cmci[,2] <- -1.7*apply(rmci[,substring(dimnames(rmci)[[2]],1,2)=="be"],2,mean)
# note that difficulty in mcmc is centered! (so diff mcmc = diff mirt - mean diff mirt)

# ELO
runElo <- function(response,ability,difficulty,step){
 change <- step*(response-plogis(ability-difficulty))
 return(list(ability=ability+change,difficulty=difficulty-change,change=change))
}
n <- length(unique(popLf2$user)); k <- length(unique(popLf2$item))
n2 <- popLf2[sample(1:nrow(popLf2)),]

n2s <- vector("numeric",length=n); names(n2s) <- unique(popLf2$user)
tmp <- aggregate(list(score=popLf2$score),list(user=popLf2$user),mean)
n2s[] <- qlogis(tmp[order(match(tmp[,1],names(n2s))),"score"])

n2i <- vector("numeric",length=k); names(n2i) <- unique(popLf2$item)
tmp <- aggregate(list(score=popLf2$score),list(item=popLf2$item),mean)
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
J <- length(unique(popLf2$user))
K <- length(unique(popLf2$item))
N <- nrow(popLf2)

y <- popLf2$score
jj <- as.numeric(substring(popLf2$user,2,3))
kk <- as.numeric(substring(popLf2$item,2,3))

dat = list(K=K,J=J,N=N,y=y,jj=jj,kk=kk)
fit0 <- stan(model_name="s1024i32pl2",model_code = my2pl, data=dat, iter = 10, verbose = FALSE) 
# rstn <- stan(fit = fit0, data = dat, iter = 8000, warmup=2000, chains = 4)

cstn <- cmrt
rstnx <- extract(rstn,permuted=T)
cstn[,1] <- colMeans(rstnx$beta)
cstn[,2] <- exp(colMeans(rstnx$log_gamma))

# plot
png(file="cmpIRTinR.png")
plot(seq(.05,.95,length.out=32),cltm[,1],xlim=c(0,1),type='n',ylim=c(-3,3),xlab="P(Success)",ylab="Difficulty",main="IRT calibration in R for 2PL items")
lines(seq(.05,.95,length.out=32),cltm[,1],lwd=2,col=1)
points(seq(.05,.95,length.out=32),cltm[,1],cex=cltm[,2])
lines(seq(.05,.95,length.out=32),cmrt[,1],col=2,lwd=2)
points(seq(.05,.95,length.out=32),cmrt[,1],cex=cmrt[,2],col=2)
lines(seq(.05,.95,length.out=32),cmci[,1]+mean(cltm[,1]),col=3,lwd=2)
points(seq(.05,.95,length.out=32),cmci[,1]+mean(cltm[,1]),cex=cmci[,2],col=3)
lines(seq(.05,.95,length.out=32),1.4*(cstn[,1]+mean(cltm[,1])),col=4,lwd=2)
points(seq(.05,.95,length.out=32),1.4*(cstn[,1]+mean(cltm[,1])),cex=cstn[,2],col=4)
lines(seq(.05,.95,length.out=32),celor+mean(cltm[,1]),col=5,lwd=2)
legend("topleft",c("ltm","mirt (diff=-1*diff/disc)","MCMCpack (diff=-1*diff+.4; disc=-1.7*disc)","rstan (diff=1.4*(diff+.4))","elo"),lty=1,lwd=2,col=1:5)
points(seq(.05,.95,length.out=32),rep(-3,32),pch=19,col=8,cex=rep(c(c(1,1.5,.5,2.5,.05,1.5,.5,1),c(1,1.5,.5,2.5,.05,1.5,.5,1)[8:1]),32/16))
dev.off()
# savePlot("compIRTinR.png",type="png")
