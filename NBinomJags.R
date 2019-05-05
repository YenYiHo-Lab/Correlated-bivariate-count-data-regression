library(rjags)

### Indirect Poisson-gamma mixture approach
# y is n-by-2 response matrix in integers; x is n-by-1 covariate

jags_data = list(N=n,Y=y,x=x)
jags_inits = list(beta0=c(0,0),beta1=c(0,0),gamma0=c(0,0),gamma1=c(0,0),
	tau=c(0,0))
jags_model = jags.model("indirect.txt",data=jags_data,n.adapt=0,
	inits=jags_inits)
update(jags_model,n.iter=2000)
coda_sample = coda.samples(jags_model,variable.names=c('tau'),
	n.iter=8000,thin=10)

### Direct approach

jags_data = list(N=n,Y=y,x=x)
jags_inits = list(beta=rep(0,4),gamma=rep(0,4),tau=rep(0,2))
jags_model = jags.model("direct.txt",data=jags_data,n.adapt=0,
	inits=jags_inits)
update(jags_model,n.iter=2000)
coda_sample = coda.samples(jags_model,variable.names=c('tau'),
	n.iter=8000,thin=10)