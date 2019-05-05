library(mvtnorm);library(MASS);library(msm)

#############################################################################

mvnorm.dens = function(ind,z,rho){
	s = matrix(c(1,rho[ind],rho[ind],1),2,2)
	f = dmvnorm(z[ind,],rep(0,2),s,log=T)
return(f)}

rho.f = function(tau,x){
	rho = c(tanh(.5*(x%*%tau)))
	rho = ifelse(rho<1-10^(-3),rho,1-10^(-3))
	rho = ifelse(rho>-1+10^(-3),rho,-1+10^(-3))
return(rho)}

trunc.norm = function(y,mu,size,rho,z){
	lb = pnbinom(y-1,mu=mu,size=size); ub = pnbinom(y,mu=mu,size=size)
	lb = ifelse(lb<ub,lb,lb-10^(-6))
	z.draw = rtnorm(length(y),rho*z,sqrt(1-rho^2),qnorm(lb),qnorm(ub))
return(z.draw)}

Llike = function(y,mu,size,rho,z){
	ub = pnorm((qnorm(pnbinom(y,mu=mu,size=size))-rho*z)/sqrt(1-rho^2))
	ub = ifelse(ub<1,ub,ub-10^(-7))
	lb = pnorm((qnorm(pnbinom(y-1,mu=mu,size=size))-rho*z)/sqrt(1-rho^2))
	lb = ifelse(lb<1,lb,lb-10^(-6))
	lb = ifelse(lb<ub,lb,lb-10^(-6))
return(sum(log(ub-lb)))}

b.post = function(b,y,x,size,rho,z,prior.mean,prior.cov){
	Lprior = dmvnorm(b,prior.mean,prior.cov,log=T)
	mu = c(exp(x%*%b))
	Llike = Llike(y,mu,size,rho,z)
return(Llike+Lprior)}

g.post = function(g,y,x,mu,rho,z,prior.mean,prior.cov){
	Lprior = dmvnorm(g,prior.mean,prior.cov,log=T)
	size = c(exp(x%*%g))
	Llike = Llike(y,mu,size,rho,z)
return(Llike+Lprior)}

tau.post = function(tau,x,z,prior.mean,prior.cov){
	Lprior = dmvnorm(tau,prior.mean,prior.cov,log=T)
	ind = matrix(1:(dim(x)[1]),ncol=1)
	Llike = sum(apply(ind,1,mvnorm.dens,z,rho.f(tau,x)))
return(Lprior+Llike)}

b.update = function(b,y,x,size,rho,z,prior.mean,prior.cov,iter,V0,C,B=100){
	if(iter==2){b.old = b} else {b.old = b[iter-1,]}

	if(iter<=B){
		repeat{
		b.star = c(rmvnorm(1,mean=b.old,sigma=C*V0))
		r = b.post(b.star,y,x,size,rho,z,prior.mean,prior.cov) - 
			b.post(b.old,y,x,size,rho,z,prior.mean,prior.cov)
		if(is.na(r)==F) break}
	} else{
		repeat{
		b.star = c(rmvnorm(1,mean=b.old,sigma=C*(1-1/(iter-1))*cov(b)))
		r = b.post(b.star,y,x,size,rho,z,prior.mean,prior.cov) - 
			b.post(b.old,y,x,size,rho,z,prior.mean,prior.cov)
		if(is.na(r)==F) break}
	}
	if(r > log(runif(1))){
		return(list(b.draw=b.star,mu=c(exp(x%*%b.star)),accept=1))
	} else return(list(b.draw=b.old,mu=c(exp(x%*%b.old)),accept=0))
}

g.update = function(g,y,x,mu,rho,z,prior.mean,prior.cov,iter,V0,C,B=100){
	if(iter==2){g.old = g} else {g.old = g[iter-1,]}

	if(iter<=B){
		repeat{
		g.star = c(rmvnorm(1,mean=g.old,sigma=C*V0))
		r = g.post(g.star,y,x,mu,rho,z,prior.mean,prior.cov) - 
			g.post(g.old,y,x,mu,rho,z,prior.mean,prior.cov)
		if(is.na(r)==F) break}
	} else{
		repeat{
		g.star = c(rmvnorm(1,mean=g.old,sigma=C*(1-1/(iter-1))*cov(g)))
		r = g.post(g.star,y,x,mu,rho,z,prior.mean,prior.cov) - 
			g.post(g.old,y,x,mu,rho,z,prior.mean,prior.cov)
		if(is.na(r)==F) break}
	}
	if(r > log(runif(1))){
		return(list(g.draw=g.star,size=c(exp(x%*%g.star)),accept=1))
	} else return(list(g.draw=g.old,size=c(exp(x%*%g.old)),accept=0))
}

tau.update = function(tau,x,z,prior.mean,prior.cov,iter,V0,C,B=200){
	if(iter==2){tau.old = tau} else {tau.old = tau[iter-1,]}

	if(iter<=B){
		repeat{
		tau.star = c(rmvnorm(1,mean=tau.old,sigma=C*V0))
		r = tau.post(tau.star,x,z,prior.mean,prior.cov) - 
			tau.post(tau.old,x,z,prior.mean,prior.cov)
		if(is.na(r)==F) break}
	} else{
		repeat{
		tau.star = c(rmvnorm(1,mean=tau.old,
			sigma=C*(1-1/(iter-1))*cov(tau)))
		r = tau.post(tau.star,x,z,prior.mean,prior.cov) - 
			tau.post(tau.old,x,z,prior.mean,prior.cov)
		if(is.na(r)==F) break}
	}
	if(r > log(runif(1))){
		return(list(tau.draw=tau.star,rho=rho.f(tau.star,x)))
	} else return(list(tau.draw=tau.old,rho=rho.f(tau.old,x)))
}

#############################################################################
### g-prior type prior distribution for regression coefficients;
### V0 and C correspond to the covariance matrix and step size of the proposal
### distribution in the adaptive MCMC scheme (Haario et al. 2001)

n = length(x); x.mat = cbind(rep(1,n),x)
prior.mean = rep(0,2); prior.cov = n*solve(t(x.mat)%*%x.mat) 
V0 = solve(t(x.mat)%*%x.mat)

z = matrix(rnorm(2*n),ncol=2)

n.mcmc = 10000
burn = 2000; thin = 10
b1 = b2 = g1 = g2 = matrix(0,ncol=2,nrow=n.mcmc)
tau = matrix(0,ncol=2,nrow=n.mcmc)
g1.draw = g2.draw = list(g.draw=c(0,0),size=c(exp(x.mat%*%c(0,0))))
tau.draw = list(tau.draw=c(0,0),rho=rep(0,n))

for(i in 2:n.mcmc){
b1.draw = b.update(b1[1:(i-1),],y[,1],x.mat,g1.draw$size,tau.draw$rho,z[,2],
	prior.mean,prior.cov,iter=i,V0,C=.025)
b1[i,] = b1.draw$b.draw
g1.draw = g.update(g1[1:(i-1),],y[,1],x.mat,b1.draw$mu,tau.draw$rho,z[,2],
	prior.mean,prior.cov,iter=i,V0,C=.025)
g1[i,] = g1.draw$g.draw
z[,1] = trunc.norm(y[,1],b1.draw$mu,g1.draw$size,tau.draw$rho,z[,2])

b2.draw = b.update(b2[1:(i-1),],y[,2],x.mat,g2.draw$size,tau.draw$rho,z[,1],
	prior.mean,prior.cov,iter=i,V0,C=.025)
b2[i,] = b2.draw$b.draw
g2.draw = g.update(g1[1:(i-1),],y[,2],x.mat,b2.draw$mu,tau.draw$rho,z[,1],
	prior.mean,prior.cov,iter=i,V0,C=.025)
g2[i,] = g2.draw$g.draw
z[,2] = trunc.norm(y[,2],b2.draw$mu,g2.draw$size,tau.draw$rho,z[,1])

tau.draw = tau.update(tau[1:(i-1),],x.mat,z,prior.mean,prior.cov,iter=i,V0,C=.05)
tau[i,] = tau.draw$tau.draw
}