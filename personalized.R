library(pls); library(MASS); library(msm); library(MCMCpack); library(dlm);library(bindata)

log.tau <- function(tau,Z.tmp,Y_E_star,theta,beta_E,X_E)
	sum( -(Y_E_star+theta+X_E%*%beta_E)^2/(2*exp(2*tau*Z.tmp))-tau*Z.tmp )

mcmc.arms <- function(Y_T,Y_E,Z.tmp,tox.scores,eff.scores,n,n.tox,n.eff,N.post,N.burnin) {
	

	beta_T_t <- matrix(0,N.post,(n.tox+1))
	beta_E_t <- matrix(0,N.post,(n.eff+1))

	beta0_T <- rep(beta.prior.mean,(n.tox+1)) 
	beta0_E <- rep(beta.prior.mean,(n.eff+1)) 
	Sigma0_T <- diag(beta.prior.sd^2,(n.tox+1)) 
	Sigma0_E <- diag(beta.prior.sd^2,(n.eff+1)) 

	alpha_t <- rep(0,N.post)
	gamma_t <- alpha_t
	tau_t <- alpha_t
	sigma2_t <- alpha_t

	theta <- rnorm(n,0,1)
	gamma <- 1
	alpha <- 1
	beta_T <- rep(.1,(n.tox+1))
	beta_E <- rep(.1,(n.eff+1))
	tau <- .5
	sigma2 <- .2
	X_T <- cbind(1,tox.scores); X_E <- cbind(1,eff.scores)

	Y_T_star_1 <- ifelse(Y_T==1,-.5,Y_T)
	Y_T_star_2 <- ifelse(Y_T_star_1==2,1,Y_T_star_1)
	Y_T_star <- ifelse(Y_T_star_2==3,3,Y_T_star_2)

	Y_E_star_1 <- ifelse(Y_E==1,-.5,Y_E)
	Y_E_star_2 <- ifelse(Y_E_star_1==2,1,Y_E_star_1)
	Y_E_star <- ifelse(Y_E_star_2==3,3,Y_E_star_2)

	for (ite in 1:N.post) {
		min_gamma <- ifelse(sum(Y_T==2)==0,0,max(max(subset(Y_T_star,Y_T==2)),0))
		max_gamma <- ifelse(sum(Y_T==3)==0,U_gamma,min(min(subset(Y_T_star,Y_T==3)),U_gamma))
		gamma <- runif(1,min_gamma,max_gamma)
		gamma_t[ite] <- gamma
		gamma_cut <- c(-Inf,0,gamma,Inf)

		Y_T_star <- rtnorm(n,-theta-X_T%*%beta_T,rep(1,n),lower=gamma_cut[Y_T],upper=gamma_cut[Y_T+1])

		Sigma.betaT <- solve(solve(Sigma0_T)+t(X_T)%*%X_T)
		mean.betaT <- Sigma.betaT%*%( solve(Sigma0_T)%*%beta0_T-t(X_T)%*%(Y_T_star+theta) )
		beta_T <- mvrnorm(1,mean.betaT,Sigma.betaT)
		beta_T_t[ite,] <- t(beta_T)	

		min_alpha <- ifelse(sum(Y_E==2)==0,0,max(max(subset(Y_E_star,Y_E==2)),0))
		max_alpha <- ifelse(sum(Y_E==3)==0,U_alpha,min(min(subset(Y_E_star,Y_E==3)),U_alpha))
		alpha <- runif(1,min_alpha,max_alpha)
		alpha_t[ite] <- alpha
		alpha_cut <- c(-Inf,0,alpha,Inf)

		Y_E_star <- rtnorm(n,-theta-X_E%*%beta_E,exp(tau*Z.tmp),lower=alpha_cut[Y_E],upper=alpha_cut[Y_E+1])

		Sigma.Y.star <- diag(exp(2*tau*Z.tmp))
		Sigma.betaE <- solve(solve(Sigma0_E)+t(X_E)%*%solve(Sigma.Y.star)%*%X_E)
		mean.betaE <- Sigma.betaE%*%( solve(Sigma0_E)%*%beta0_E-t(X_E)%*%solve(Sigma.Y.star)%*%(Y_E_star+theta) )
		beta_E <- mvrnorm(1,mean.betaE,Sigma.betaE)	
		beta_E_t[ite,] <- t(beta_E)	
	
		var.theta <- 1/( 1+1/exp(2*tau*Z.tmp)+1/sigma2 )
		mean.theta <- -(Y_T_star+X_T%*%beta_T+(Y_E_star+X_E%*%beta_E)/exp(2*tau*Z.tmp))*var.theta 
		theta <- rnorm(n,mean.theta,sqrt(var.theta))
		
		sigma2 <- rinvgamma(1,n/2+eta,eta+sum(theta^2)/2)
		sigma2_t[ite] <- sigma2
		logden <- function(x) 
			log.tau(x,Z.tmp,Y_E_star,theta,beta_E,X_E)
		tau <- arms(tau,logden,function(x) ((x>-cc)*(x<cc)),1)
		tau_t[ite] <- tau	

	}
	list(beta_T=beta_T_t[(N.burnin+1):N.post,],beta_E=beta_E_t[(N.burnin+1):N.post,],alpha=alpha_t[(N.burnin+1):N.post],sigma2=sigma2_t[(N.burnin+1):N.post],gamma=gamma_t[(N.burnin+1):N.post],tau=tau_t[(N.burnin+1):N.post])	
}

marginal.prob.integrand <- function(theta,tox,cut,beta,dose,cov.vec,sigma2,tau,true.prob) {
	
	if (true.prob==0)
		Xmat <- c(1,cov.vec) else Xmat <- c(1,dose,cov.vec,dose*cov.vec)
	ifelse(tox==1,pnorm(theta+cut+t(beta)%*%Xmat),pnorm((theta+cut+t(beta)%*%Xmat)/exp(tau*dose)))*exp(-theta^2/(2*sigma2))/(sqrt(2*pi*sigma2))

}

marginal.prob <- function(input.vec,cov.vec,tox,p,true.prob) {
	
	beta <- input.vec[1:(p-3)]; cut <- input.vec[p-2]; sigma2 <- input.vec[p-1]; tau <- input.vec[p]
	
	cut.vec <- c(0,cut)
	prob <- matrix(0,n.dose,3)
	for (i in 1:n.dose) {
		if (true.prob==0) cov.vec0 <- cov.vec[i,] else cov.vec0 <- cov.vec
		prob[i,1] <- integrate(marginal.prob.integrand,-Inf,Inf,tox,cut=cut.vec[1],beta=beta,dose=doses[i],cov.vec=cov.vec0,sigma2=sigma2,tau,true.prob)$value
		cumu2 <- integrate(marginal.prob.integrand,-Inf,Inf,tox,cut=cut.vec[2],beta=beta,dose=doses[i],cov.vec=cov.vec0,sigma2=sigma2,tau,true.prob)$value
		prob[i,2] <- cumu2-prob[i,1]
		prob[i,3] <- 1-cumu2
	}
	return(as.vector(prob))
}

joint.prob.integrand <- function(theta,j,k,cutT,cutE,betaT,betaE,dose,pred.tox.i,pred.eff.i,sigma2,tau,true.prob) {
	if (true.prob==0) {
		XmatT <- c(1,pred.tox.i)
		XmatE <- c(1,pred.eff.i)
	} else {
		XmatT <- c(1,dose,pred.tox.i,dose*pred.tox.i)
		XmatE <- c(1,dose,pred.eff.i,dose*pred.eff.i)
	}
	cumuT.1 <- pnorm(theta+t(betaT)%*%XmatT)
	cumuT.2 <- pnorm(theta+cutT+t(betaT)%*%XmatT)
	cumuE.1 <- pnorm((theta+t(betaE)%*%XmatE)/exp(tau*dose))
	cumuE.2 <- pnorm((theta+cutE+t(betaE)%*%XmatE)/exp(tau*dose))

	if (j == 1) tox.part <- cumuT.1
	if (j == 2) tox.part <- cumuT.2-cumuT.1
	if (j == 3) tox.part <- 1-cumuT.2
	if (k == 1) eff.part <- cumuE.1
	if (k == 2) eff.part <- cumuE.2-cumuE.1
	if (k == 3) eff.part <- 1-cumuE.2


	return(tox.part*eff.part*exp(-theta^2/(2*sigma2))/(sqrt(2*pi*sigma2)))

}

pair.ind <- function(j,k) {
	if (j==1) pair <- k
	if (j==2) pair <- k+3
	if (j==3) pair <- k+6
	return(pair)
}

joint.prob <- function(input.vec,pred.tox,pred.eff,n.comp.tox,n.comp.eff) {
	betaT <- input.vec[1:(n.comp.tox+1)]; gamma <- input.vec[n.comp.tox+2]; 
	betaE <- input.vec[(n.comp.tox+3):(n.comp.tox+n.comp.eff+3)]; alpha <- input.vec[n.comp.tox+n.comp.eff+4]; 
	sigma2 <- input.vec[n.comp.tox+n.comp.eff+5]; tau <- input.vec[n.comp.tox+n.comp.eff+6]

	prob <- matrix(0,ncol=n.dose,nrow=9)
	for (i in 1:n.dose) {
		for (j in 1:3) {
			for (k in 1:3) {
				prob[pair.ind(j,k),i] <- integrate(joint.prob.integrand,-Inf,Inf,j,k,cutT=gamma,cutE=alpha,betaT=betaT,betaE=betaE,dose=doses[i],pred.tox[i,],pred.eff[i,],sigma2=sigma2,tau,0)$value
			}
		}
	}
	return(as.vector(prob))
}

true.joint.prob <- function(input.vec,pred.tox,pred.eff,p) {
	
	betaT <- input.vec[1:(p-16)]; gamma <- input.vec[p-15]; 
	betaE <- input.vec[(p-14):(p-3)]; alpha <- input.vec[p-2]; 
	sigma2 <- input.vec[p-1]; tau <- input.vec[p]

	
	prob <- matrix(0,ncol=n.dose,nrow=9)
	for (i in 1:n.dose) {
		for (j in 1:3) {
			for (k in 1:3) {
				prob[pair.ind(j,k),i] <- integrate(joint.prob.integrand,-Inf,Inf,j,k,cutT=gamma,cutE=alpha,betaT=betaT,betaE=betaE,dose=doses[i],pred.tox,pred.eff,sigma2=sigma2,tau,1)$value
			}
		}
	}
	return(prob)
}



p.E <- function(v,phi)
	sum(v>=phi)/(N.post-N.burnin)
p.T <- function(v)
	sum(v<=phi_T)/(N.post-N.burnin)

summ.mcmc <- function(mcmc.tmp,pred.tox,pred.eff,n.comp.tox,n.comp.eff) {
	tox.mcmc.matrix <- cbind(mcmc.tmp$beta_T,mcmc.tmp$gamma,mcmc.tmp$sigma2,mcmc.tmp$tau)
	eff.mcmc.matrix <- cbind(mcmc.tmp$beta_E,mcmc.tmp$alpha,mcmc.tmp$sigma2,mcmc.tmp$tau)
	summ.tox <- t(apply(tox.mcmc.matrix,1,marginal.prob,cov.vec=pred.tox,tox=1,p=4+n.comp.tox,0))
	summ.eff <- t(apply(eff.mcmc.matrix,1,marginal.prob,cov.vec=pred.eff,tox=0,p=4+n.comp.eff,0))
	
	summ.tox.3 <- apply(summ.tox[,(11:15)],2,p.T)
	summ.eff.3 <- apply(summ.eff[,(11:15)],2,p.E,phi_E3)
	summ.eff.2 <- apply(summ.eff[,(6:10)],2,p.E,phi_E2)
	accept.tox.ind <- which(summ.tox.3>delta_T)
	if (length(accept.tox.ind)==0) accept.ind <- NULL else {
		accept.ind <- NULL
		for (i in 1:length(accept.tox.ind)) {
			if (summ.eff.3[accept.tox.ind[i]]>=delta_E3) accept.ind <- c(accept.ind,accept.tox.ind[i])
			else if (summ.eff.2[accept.tox.ind[i]]>=delta_E2) accept.ind <- c(accept.ind,accept.tox.ind[i])
		}
	}
	
	input.matrix <- cbind(mcmc.tmp$beta_T,mcmc.tmp$gamma,mcmc.tmp$beta_E,mcmc.tmp$alpha,mcmc.tmp$sigma2,mcmc.tmp$tau)
	input.vec <- apply(input.matrix,2,mean)
	joint.probability <- matrix(joint.prob(input.vec,pred.tox,pred.eff,n.comp.tox,n.comp.eff),nrow=9)
	utility <- t(joint.probability)%*%as.vector(Utility)
	list(accept.ind=accept.ind,utility=utility)
}



which1 <- function(v)
	which(v==1)


	
n.marker <- 5
true.prob.5 <- c(.5,.55,.5,.7,.3) 
true.cor.12 <- .9; 
	
N.max <- 60; coh.size <- 3
doses <- c(.1,.3,.5,.7,.9); n.dose <- length(doses)

beta.prior.mean <- 0;beta.prior.sd <- 1.3

U_gamma <- 6; U_alpha <- 6; eta <- .1; cc <- 2
N.post <- 1300; N.burnin <- 300
phi_T <- .3; phi_E3 <- .2; phi_E2 <- .25
delta_T <- .05; delta_E3 <- .05; delta_E2 <- .05

Utility <- matrix(c(10,5,0,50,20,10,100,60,20),byrow=T,nrow=3)


patient.cov.simulation <- function(N.max,true.prob.5,true.cor.12) {
	x.12 <- rmvbin(N.max,margprob=true.prob.5[1:2],bincorr=matrix(c(1,true.cor.12,true.cor.12,1),ncol=2))
	x.3 <- rbinom(N.max,1,true.prob.5[3])
	x.4 <- rbinom(N.max,1,true.prob.5[4])
	x.5 <- rbinom(N.max,1,true.prob.5[5])
	
	x <- cbind(x.12,x.3,x.4,x.5)
	return(x)
}
patient.out <- function(Z.tmp,patient.index,patient.cov,true.betaT,true.betaE,true.sigma2,true.tau,true.gamma,true.alpha) {
	n.pat <- length(patient.index)
	cov <- patient.cov[patient.index,]
	theta <- rnorm(n.pat,0,sqrt(true.sigma2))
	if (n.pat==1) X <- c(1,Z.tmp,cov,Z.tmp*cov) else X <- cbind(1,rep(Z.tmp,n.pat),cov,Z.tmp*cov)
	(p.1.T <- pnorm(theta+X%*%true.betaT))
	cumu.2.T <- pnorm(theta+true.gamma+X%*%true.betaT)
	cumu.T <- cbind(p.1.T,cumu.2.T,1)

	(p.1.E <- pnorm((theta+X%*%true.betaE)/exp(true.tau*Z.tmp)))
	cumu.2.E <- pnorm((theta+true.alpha+X%*%true.betaE)/exp(true.tau*Z.tmp))
	cumu.E <- cbind(p.1.E,cumu.2.E,1)
	
	r.T <- runif(n.pat,0,1)
	r.E <- runif(n.pat,0,1)
	out.T <- rep(0,n.pat)
	out.E <- rep(0,n.pat)
	for (i in 1:n.pat) {
		out.T[i] <- min(which(r.T[i]<cumu.T[i,]))
		out.E[i] <- min(which(r.E[i]<cumu.E[i,]))
	}
	list(out.T=out.T,out.E=out.E)
}



random200 <- matrix(0,nrow=32,ncol=5)
random200[1,] <- c(1,1,0,0,0)
random200[2,] <- c(1,1,0,1,0)
random200[3,] <- c(1,1,0,0,1)
random200[4,] <- c(1,1,1,0,0)
random200[5,] <- c(1,1,1,1,0)
random200[6,] <- c(1,1,0,1,1)
random200[7,] <- c(1,1,1,0,1)
random200[8,] <- c(1,1,1,1,1)

random200[9:16,] <- random200[1:8,]
random200[9:16,1:2] <- 0

random200[17:24,] <- random200[1:8,]
random200[17:24,1] <- 0

random200[25:32,] <- random200[1:8,]
random200[25:32,2] <- 0

repre.pattern <- random200[1:16,]

patient.cov.simulation200 <- function(N.max,true.prob.7,true.cor.12,true.cor.23,true.cor.13,true.cor.45) {
	x.1to3 <- rmvbin(N.max,margprob=true.prob.7[1:3],bincorr=matrix(c(1,true.cor.12,true.cor.13,true.cor.12,1,true.cor.23,true.cor.13,true.cor.23,1),ncol=3))
	x.45 <- rmvbin(N.max,margprob=true.prob.7[4:5],bincorr=matrix(c(1,true.cor.45,true.cor.45,1),ncol=2))
	x.6 <- rbinom(N.max,1,true.prob.7[6])
	x.7 <- rbinom(N.max,1,true.prob.7[7])
	x.8 <- apply(t(rmultinom(N.max,1,c(1/3,1/3,1/3))),1,which1)-2
	x <- cbind(x.1to3,x.45,x.6,x.7,x.8)
	return(x)
}
count.rows <- 
   function(x) 
   { 
     order.x <- do.call(order,as.data.frame(x)) 
     equal.to.previous <- 
       rowSums(x[tail(order.x,-1),] != x[head(order.x,-1),])==0 
     tf.runs <- rle(equal.to.previous) 
     counts <- c(1, 
                 unlist(mapply( function(x,y) if (y) x+1 else (rep(1,x)), 
                               tf.runs$length, tf.runs$value ))) 
     counts <- counts[ c(diff(counts) <= 0, TRUE ) ] 
     unique.rows <- which( c(TRUE, !equal.to.previous ) ) 
     cbind( counts, x[order.x[ unique.rows ], ,drop=F] ) 
   } 

unique.random200 <- count.rows(random200)
n.unique.random200 <- nrow(unique.random200)
unique.rows.random200 <- unique.random200[,-1]
num.unique.random200 <- unique.random200[,1]



which.less <- function(v)
	which(v<=phi_T)
which.greater <- function(v,phi)
	which(v>=phi)
which.less200 <- function(v)
	which(v<=.4)

true.best.dose200 <- function(pattern,true.betaT,true.betaE,true.gamma,true.alpha,true.sigma2,true.tau) {
	
	true.marginal.tox3 <- matrix(0,nrow(pattern),5)
	true.marginal.eff3 <- matrix(0,nrow(pattern),5)
	true.marginal.eff2 <- matrix(0,nrow(pattern),5)
	for (i in 1:nrow(pattern)) {
		true.marginal.tox3[i,] <- matrix(marginal.prob(c(true.betaT,true.gamma,true.sigma2,true.tau),pattern[i,],1,2*n.marker+5,1),nrow=5)[,3]
		true.marginal.eff.i <- matrix(marginal.prob(c(true.betaE,true.alpha,true.sigma2,true.tau),pattern[i,],0,2*n.marker+5,1),nrow=5)
		true.marginal.eff3[i,] <- true.marginal.eff.i[,3]	
		true.marginal.eff2[i,] <- true.marginal.eff.i[,2]
	}


	accept.tox.ind <- apply(true.marginal.tox3,1,which.less200)
	accept.ind <- vector("list",nrow(pattern))
	uti <- matrix(0,nrow(pattern),5)
	best.dose.ind <- vector("list",nrow(pattern))
	for (i in 1:nrow(pattern)) {
		if (is.matrix(accept.tox.ind)==T) accept.tox <- accept.tox.ind[,i] else
		accept.tox <- accept.tox.ind[[i]]
		accept.eff3 <- true.marginal.eff3[i,accept.tox]
		accept.eff2 <- true.marginal.eff2[i,accept.tox]
		accept.eff3.ind <- accept.tox[which.greater(accept.eff3,.1)]
		accept.eff2.ind <- accept.tox[which.greater(accept.eff2,.15)]
		joint.proba <-  true.joint.prob(c(true.betaT,true.gamma,true.betaE,true.alpha,true.sigma2,true.tau),pattern[i,],pattern[i,],4*n.marker+8)
		uti[i,] <- t(joint.proba)%*%as.vector(Utility)
		if (length(unique(c(accept.eff3.ind,accept.eff2.ind)))!=0) {
			accept.ind[[i]] <- sort(unique(c(accept.eff3.ind,accept.eff2.ind)))
			uti.accept.i <- uti[i,accept.ind[[i]]]
			max.uti.accept.i <- max(uti.accept.i)
			best.dose.ind[[i]] <- accept.ind[[i]][which(uti.accept.i>=0.9*max.uti.accept.i)]
		}	 else {
				accept.ind[[i]] <- 0		
				best.dose.ind[[i]] <- 0
			} 
	}
	list(uti=uti,best.dose=best.dose.ind,accept.ind=accept.ind,true.tox3=true.marginal.tox3,true.eff2=true.marginal.eff2,true.eff3=true.marginal.eff3)
}
true.best.dose <- function(pattern,true.betaT,true.betaE,true.gamma,true.alpha,true.sigma2,true.tau) {
	
	true.marginal.tox3 <- matrix(0,nrow(pattern),5)
	true.marginal.eff3 <- matrix(0,nrow(pattern),5)
	true.marginal.eff2 <- matrix(0,nrow(pattern),5)
	for (i in 1:nrow(pattern)) {
		true.marginal.tox3[i,] <- matrix(marginal.prob(c(true.betaT,true.gamma,true.sigma2,true.tau),pattern[i,],1,2*n.marker+5,1),nrow=5)[,3]
		true.marginal.eff.i <- matrix(marginal.prob(c(true.betaE,true.alpha,true.sigma2,true.tau),pattern[i,],0,2*n.marker+5,1),nrow=5)
		true.marginal.eff3[i,] <- true.marginal.eff.i[,3]	
		true.marginal.eff2[i,] <- true.marginal.eff.i[,2]
	}


	accept.tox.ind <- apply(true.marginal.tox3,1,which.less)
	accept.ind <- vector("list",nrow(pattern))
	uti <- matrix(0,nrow(pattern),5)
	best.dose.ind <- vector("list",nrow(pattern))
	for (i in 1:nrow(pattern)) {
		if (is.matrix(accept.tox.ind)==T) accept.tox <- accept.tox.ind[,i] else
		accept.tox <- accept.tox.ind[[i]]
		accept.eff3 <- true.marginal.eff3[i,accept.tox]
		accept.eff2 <- true.marginal.eff2[i,accept.tox]
		accept.eff3.ind <- accept.tox[which.greater(accept.eff3,phi_E3)]
		accept.eff2.ind <- accept.tox[which.greater(accept.eff2,phi_E2)]
		joint.proba <-  true.joint.prob(c(true.betaT,true.gamma,true.betaE,true.alpha,true.sigma2,true.tau),pattern[i,],pattern[i,],4*n.marker+8)
		uti[i,] <- t(joint.proba)%*%as.vector(Utility)
		if (length(unique(c(accept.eff3.ind,accept.eff2.ind)))!=0) {
			accept.ind[[i]] <- sort(unique(c(accept.eff3.ind,accept.eff2.ind)))
			uti.accept.i <- uti[i,accept.ind[[i]]]
			max.uti.accept.i <- max(uti.accept.i)
			best.dose.ind[[i]] <- accept.ind[[i]][which(uti.accept.i>=0.9*max.uti.accept.i)]
		}	 else {
				accept.ind[[i]] <- 0		
				best.dose.ind[[i]] <- 0
			} 
	}
	list(uti=uti,best.dose=best.dose.ind,accept.ind=accept.ind,true.tox3=true.marginal.tox3,true.eff2=true.marginal.eff2,true.eff3=true.marginal.eff3)
}


p.deescalate <- 0.3585195
p.escalate <- 0.2364907

main <- function(scenario) {

	
	
	
 	true.betaT <- c(-.2,1,.5,-.45,.1,0,.1,-.7,-.6,0,-.1,-1)
	length(true.betaT)
	 true.betaE <- c(-.4,1.4,-.6,0,-1,0,0,1.1,-1,0,-2,-.1)
	length(true.betaE)
 	true.sigma2 <- 0.5
	 true.tau <- -1
 	true.gamma <- 2.5
 	true.alpha <- 1.3

	
	
	

	true.200 <- true.best.dose200(random200,true.betaT,true.betaE,true.gamma,true.alpha,true.sigma2,true.tau)
	
	

	patient.cov <- patient.cov.simulation(N.max,true.prob.5,true.cor.12)
	
	num.tox <- 0
	dose.ind <- 1
	n <- 0
	Z <- NULL
	Y_T <- NULL; Y_E <- NULL

	for (i in 1:7) {

		
		pat.out <- patient.out(doses[dose.ind],((n+1):(n+coh.size)),patient.cov,true.betaT,true.betaE,true.sigma2,true.tau,true.gamma,true.alpha)
		tox.out <- pat.out$out.T
		eff.out <- pat.out$out.E
		Y_T <- c(Y_T,tox.out)
		Y_E <- c(Y_E,eff.out)
		num.tox <- sum(tox.out==3)		
		Z <- c(Z,rep(doses[dose.ind],coh.size))	
		n <- n + coh.size

		p.hat <- sum((Y_T[Z==doses[dose.ind]])==3)/length(Y_T[Z==doses[dose.ind]])

		if ((dose.ind==1) & (p.hat>p.deescalate)) act <- 0 else {
			act <- 1
			if (p.hat>p.deescalate) dose.ind <- dose.ind-1
			if ((dose.ind<5) & (p.hat<p.escalate)) dose.ind <- dose.ind+1
			if ((dose.ind==5) & (p.hat<p.escalate)) dose.ind <- dose.ind
			if ((p.hat<p.deescalate) & (p.hat>p.escalate)) dose.ind <- dose.ind
		}
	}
	if (act==0) {
		dose.ind <- 0 
		act <- 0 } else {
	
		n <- n+1; n.0 <- 0
		while (n<=N.max) {
			pat.out <- patient.out(doses[dose.ind],n,patient.cov,true.betaT,true.betaE,true.sigma2,true.tau,true.gamma,true.alpha)
			Y_T <- c(Y_T,rep(0,n.0),pat.out$out.T)
			Y_E <- c(Y_E,rep(0,n.0),pat.out$out.E)
			Z <- c(Z,rep(0,n.0),doses[dose.ind])

			ind.0 <- which(Z!=0)
			Y_T.mcmc <- Y_T[ind.0]; Y_E.mcmc <- Y_E[ind.0]; Z.mcmc <- Z[ind.0]; n.mcmc <- length(ind.0)

			
			data.tox <- data.frame(Y_T=Y_T.mcmc,X=rep(0,length(ind.0)))

			tox.dummy <- matrix(0,length(Y_T.mcmc),3)
			for (i in 1:length(Y_T.mcmc)) {
				if (Y_T.mcmc[i]==1) tox.dummy[i,1] <- 1
				if (Y_T.mcmc[i]==2) tox.dummy[i,2] <- 1
				if (Y_T.mcmc[i]==3) tox.dummy[i,3] <- 1
			}
			data.tox$dummy <- tox.dummy

			data.tox$X <- cbind(Z.mcmc,patient.cov[ind.0,],Z.mcmc*patient.cov[ind.0,])
			tox.cppls <- cppls(dummy~X,2,data=data.tox)

			n.comp.tox <- 2
			tox.scores <- tox.cppls$scores
			

			data.eff <- data.frame(Y_E=Y_E.mcmc,X=rep(0,length(ind.0)))

			eff.dummy <- matrix(0,length(Y_E.mcmc),3)
			for (i in 1:length(Y_E.mcmc)) {
				if (Y_E.mcmc[i]==1) eff.dummy[i,1] <- 1
				if (Y_E.mcmc[i]==2) eff.dummy[i,2] <- 1
				if (Y_E.mcmc[i]==3) eff.dummy[i,3] <- 1
			}
			data.eff$dummy <- eff.dummy

			data.eff$X <- cbind(Z.mcmc,patient.cov[ind.0,],Z.mcmc*patient.cov[ind.0,])
			eff.cppls <- cppls(dummy~X,2,data=data.eff)
			n.comp.eff <- 2
			eff.scores <- eff.cppls$scores

			
			tox.scores0 <- .5*stdize(tox.scores) 
			sd.score.tox <- apply(tox.scores,2,sd)
			
			eff.scores0 <- .5*stdize(eff.scores) 
			sd.score.eff <- apply(eff.scores,2,sd)
			

			mcmc.tmp <- mcmc.arms(Y_T.mcmc,Y_E.mcmc,Z.mcmc,tox.scores0,eff.scores0,n.mcmc,n.comp.tox,n.comp.eff,N.post,N.burnin) 

			if (n<N.max) {
				dose.ind <- 0
				n.0 <- 0
				while (dose.ind==0)  {
					n <- n+1
					pat.cov <- t(as.matrix(patient.cov[n,]))
					pat.cov.dose <- matrix(0,nrow=n.dose,ncol=2*n.marker+1)
					

					for (i in 1:n.dose) 
						pat.cov.dose[i,] <- c(doses[i],pat.cov,doses[i]*pat.cov)
					pred.tox0 <- predict(tox.cppls,comps=1:n.comp.tox,newdata=pat.cov.dose,type="score")
					pred.eff0 <- predict(eff.cppls,comps=1:n.comp.eff,newdata=pat.cov.dose,type="score")

					pred.tox <- stdize(pred.tox0,center=F,scale=2*sd.score.tox)
					pred.eff <- stdize(pred.eff0,center=F,scale=2*sd.score.eff)

					mcmc.summary <- summ.mcmc(mcmc.tmp,pred.tox,pred.eff,n.comp.tox,n.comp.eff) 
					accept.ind <- mcmc.summary$accept.ind
					if (length(accept.ind)==0) dose.ind <- 0 else if (length(accept.ind)==1) dose.ind <- accept.ind else {
					
						uti <- mcmc.summary$utility[accept.ind]
						uti.normalize <- uti/sum(uti)
						cumu.uti <- uti
						for (i in 1:length(accept.ind)) 
							cumu.uti[i] <- sum(uti.normalize[1:i])
						r <- runif(1,0,1)
						dose.ind <- accept.ind[min(which(r<cumu.uti))]
						
						h.dose.ind <- (max(Z)*10+1)/2
						if (dose.ind>h.dose.ind) dose.ind <- h.dose.ind+1
					}
					if (dose.ind==0)	
						n.0 <- n.0+1
					if ((n.0>=10) | (n==N.max)) break
			
				}
			}
			if ((n.0>=10) | (length(Z)==N.max) | ((dose.ind==0) & (n==N.max))) break
		}	
		if ((dose.ind==0) & (length(Z)==(N.max-1))) {
			Z <- c(Z,0)
			Y_T <- c(Y_T,0)
			Y_E <- c(Y_E,0)
		}
		dose.recomm <- rep(0,nrow(repre.pattern))	
		
		for (ii in 1:nrow(repre.pattern)) {
			pat.cov.dose <- matrix(0,nrow=n.dose,ncol=2*n.marker+1)
			for (i in 1:n.dose) 
				pat.cov.dose[i,] <- c(doses[i],repre.pattern[ii,],doses[i]*repre.pattern[ii,])
			pred.tox0 <- predict(tox.cppls,comps=1:n.comp.tox,newdata=pat.cov.dose,type="score")
			pred.eff0 <- predict(eff.cppls,comps=1:n.comp.eff,newdata=pat.cov.dose,type="score")

			pred.tox <- stdize(pred.tox0,center=F,scale=2*sd.score.tox)
			pred.eff <- stdize(pred.eff0,center=F,scale=2*sd.score.eff)

			mcmc.summary <- summ.mcmc(mcmc.tmp,pred.tox,pred.eff,n.comp.tox,n.comp.eff) 
			accept.ind <- mcmc.summary$accept.ind
			if (length(accept.ind)==0) dose.recomm[ii] <- 0 else if (length(accept.ind)==1) dose.recomm[ii] <- accept.ind else {
					
				uti <- mcmc.summary$utility[accept.ind]
				dose.recomm[ii] <- accept.ind[which(uti==max(uti))]

			}
			
		}


		dose.recomm200 <- rep(0,nrow(random200))	
		true.uti.200 <- true.200$uti	
		true.accept.200 <- true.200$accept.ind
		for (ii in 1:nrow(random200)) {
			pat.cov.dose <- matrix(0,nrow=n.dose,ncol=2*n.marker+1)
			for (i in 1:n.dose) 
				pat.cov.dose[i,] <- c(doses[i],random200[ii,],doses[i]*random200[ii,])
			pred.tox0 <- predict(tox.cppls,comps=1:n.comp.tox,newdata=pat.cov.dose,type="score")
			pred.eff0 <- predict(eff.cppls,comps=1:n.comp.eff,newdata=pat.cov.dose,type="score")

			pred.tox <- stdize(pred.tox0,center=F,scale=2*sd.score.tox)
			pred.eff <- stdize(pred.eff0,center=F,scale=2*sd.score.eff)

			mcmc.summary200 <- summ.mcmc(mcmc.tmp,pred.tox,pred.eff,n.comp.tox,n.comp.eff) 
			
			accept <- mcmc.summary200$accept.ind
			uti <- mcmc.summary200$utility[accept]
			if (length(accept)!=0)
				dose.recomm200[ii] <- accept[which(uti==max(uti))]
			else dose.recomm200[ii] <- 0
			
		}
		
		
	}
		if (act==0)
			list(dose.recomm=rep(0,nrow(repre.pattern)),dose.recomm200=rep(0,200),true.200=true.200)
		else
			list(dose.recomm=dose.recomm,dose.recomm200=dose.recomm200,true.200=true.200)

}



n.sim <- 100
dose.recomm <- matrix(0,n.sim,nrow(repre.pattern))
dose.recomm200 <- matrix(0,n.sim,32)

for (i in 1:n.sim) {
	tmp <- main(scenario=1) 
	dose.recomm[i,] <- tmp$dose.recomm
	dose.recomm200[i,] <- tmp$dose.recomm200
		
}
true200 <- tmp$true.200

result <- rep(0,32)
for (i in 1:32) {
	best.dose.i <- true200$best.dose[[i]]
	if (length(best.dose.i)!=0) {
		n.best.i <- length(best.dose.i)
		select.i <- dose.recomm200[,i]
	
		for (j in 1:n.best.i)
			result[i] <- result[i]+sum(select.i==best.dose.i[j])
	} else result[i] <- -10
}
result <- result/n.sim








