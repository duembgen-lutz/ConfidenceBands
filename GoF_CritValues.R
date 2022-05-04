# Critical values for various goodness-of-fit tests
# based on the uniform empirical process

source("NoeRecursion.R")

# In what follows,
#    U[1] < U[2] < ... < U[n]
# are the order statistics of n independent random variables
# with uniform distribution on [0,1].

# Kolmogorov-Smirnov tests ----

KS_CDF <- function(n, kappav)
	# Computes the distribution function of the
	# kolmogorov-Smirnov test statistic at the points
	# specified by kappav.
{
	u0 <- (0:(n-1))/n
	u1 <- (1:n)/n
	nk <- length(kappav)
	CDF <- rep(1,nk)
	for (j in 1:nk){
		av <- pmax(u1 - kappav[j], 0)
		bv <- pmin(u0 + kappav[j], 1)
		CDF[j] <- NoeRecursion(av,bv)
	}
	return(CDF)
}

KS_Critval <- function(n, pv=c(0.5,0.9,0.95,0.99),digits=5,
					   kappa0=NULL)
	# Computes the pv-quantiles of the
	# Kolmogorov-Smirnov test statistic
	#    max( max(U[i] - (i-1)/n, i/n - U[i]) : i=1,...,n)
	# rounded up to multiples of 10^(-digits).
	# The input parameter pv is a number or a vector of
	# strictly increasing numbers in (0,1).
	# The output is a vector kappav containing the
	# pv-quantiles, rounded up to multiples of 10^(-digits).
	# The optional input parameter kappa0 is a first guess of
	# kappav[1].
{
	u0 <- (0:(n-1))/n
	u1 <- (1:n)/n
	np <- length(pv)
	MK <- matrix(0,np,2)
	MK[,2] <- Inf
	dimnames(MK)[[1]] <- pv
	MP <- matrix(0,np,2)
	MP[,2] <- 1
	# MK[j,1] and MK[j,2] are a lower and an upper bound,
	# respectively, for 10^digits*kappav[j], while
	# MP[j,1] and MP[j,2] are the corresponding values
	# of the distribution function.
	if (!is.null(kappa0)){
		k <- floor(10^digits*kappa0)
	}else{
		k <- ceiling(10^digits*sqrt(log(2/(1 - pv[1]))/2/n))
	}
	av <- pmax(u1 - 10^(-digits)*k, 0)
	bv <- pmin(u0 + 10^(-digits)*k, 1)
	q <- NoeRecursion(av,bv)
	tmp <- (pv > q)
	MK[tmp,1] <- pmax(MK[tmp,1],k)
	MK[!tmp,2] <- pmin(MK[!tmp,2],k)
	MP[tmp,1] <- pmax(MP[tmp,1],q)
	MP[!tmp,2] <- pmin(MP[!tmp,2],q)
	print(cbind(10^(-digits)*MK,MP))
	for (j in 1:np){
		while (MK[j,2] == Inf){
			k <- 2*k
			av <- pmax(u1 - 10^(-digits)*k, 0)
			bv <- pmin(u0 + 10^(-digits)*k, 1)
			q <- NoeRecursion(av,bv)
			tmp <- (pv > q)
			MK[tmp,1] <- pmax(MK[tmp,1],k)
			MK[!tmp,2] <- pmin(MK[!tmp,2],k)
			MP[tmp,1] <- pmax(MP[tmp,1],q)
			MP[!tmp,2] <- pmin(MP[!tmp,2],q)
			print(cbind(10^(-digits)*MK,MP))
		}
		while (MK[j,2] > MK[j,1] + 1){
			pdiff <- MP[j,2] - MP[j,1]
			if (pdiff > 0){
				lambda <- (pv[j] - MP[j,1])/pdiff
				lambda <- (1/2 + 7*lambda)/8
				k <- max(MK[j,1]+1,
						 min(MK[j,2]-1,
						 	floor((1-lambda)*MK[j,1] + lambda*MK[j,2])))
			}else{
				k <- floor((MK[j,1] + MK[j,2])/2)
			}
			av <- pmax(u1 - 10^(-digits)*k, 0)
			bv <- pmin(u0 + 10^(-digits)*k, 1)
			q <- NoeRecursion(av,bv)
			tmp <- (pv > q)
			MK[tmp,1] <- pmax(MK[tmp,1],k)
			MK[!tmp,2] <- pmin(MK[!tmp,2],k)
			MP[tmp,1] <- pmax(MP[tmp,1],q)
			MP[!tmp,2] <- pmin(MP[!tmp,2],q)
			print(cbind(10^(-digits)*MK,MP))
		}
	}
	return(10^(-digits)*MK[,2])
}

KS_TestCritval <- function(n,kappav,pv=c(0.5,0.9,0.95,0.99),
						   mcsim=10000)
	# Test empirically whether a given critical value
	# kappa could be the (1 - alpha)-quantile of the
	# Kolmogorov-Smirnov test statistic.
{
	u0 <- (0:(n-1))/n
	u1 <- (1:n)/n
	np <- length(pv)
	H <- rep(0,np)
	for (s in 1:mcsim){
		U <- sort(runif(n))
		H <- H + (max(max(U - u0),max(u1 - U)) <= kappav)
	}
	Table <- cbind('kappa'=kappav,
				   'p'=pv,
				   'p.emp'=H/mcsim,
				   'p-value'=2*pmin(pbinom(H,mcsim,pv),
				   				 pbinom(mcsim-H,mcsim,1-pv)))
	print(Table,quote=FALSE)
}

# Stepanova-Pavlenko tests (for confidence bands) ----

SP_CDF <- function(n, kappav,gamma=1)
	# Computes the distribution function of the
	# Stepanova-Pavlenko test statistic
	#    max( max{(U[i] - (i-1)/n)/q((i-1)/n,
	#             (i/n - U[i])/q(i/n)} : i=1,...,n)
	# at the points specified by kappav.
{
	u <- (1:(n-1))/n
	v <- u*(1-u)
	qv <- sqrt(v*log(log(gamma/v))/n)
	nk <- length(kappav)
	CDF <- rep(1,nk)
	for (j in 1:nk){
		av <- c(pmax(u - kappav[j]*qv, 0), 0)
		av[n] <- av[n-1]
		bv <- c(1,pmin(u + kappav[j]*qv, 1))
		bv[1] <- bv[2]
		CDF[j] <- NoeRecursion(av,bv)
	}
	return(CDF)
}

SP_Critval <- function(n, pv=c(0.5,0.9,0.95,0.99),gamma=1,
					   digits=5,
					   kappa0=NULL)
	# Computes the (1 - alpha)-quantile of the
	# Stepanova-Pavlenko test statistic
	#    max( max{(U[i] - (i-1)/n)/q((i-1)/n,
	#             (i/n - U[i])/q(i/n)} : i=1,...,n)
	# rounded up to multiples of 10^(-digits), where
	#    q(u) := sqrt(u*(1-u)*log(log(gamma/[u*(1-u)])))
	# for 0 < u < 1 and q(0) := q(1) := Inf.
{
	u <- (1:(n-1))/n
	v <- u*(1-u)
	qv <- sqrt(v*log(log(gamma/v))/n)
	np <- length(pv)
	MK <- matrix(0,np,2)
	MK[,2] <- Inf
	dimnames(MK)[[1]] <- pv
	MP <- matrix(0,np,2)
	MP[,2] <- 1
	# MK[j,1] and MK[j,2] are a lower and an upper bound,
	# respectively, for 10^digits*kappav[j], while
	# MP[j,1] and MP[j,2] are the corresponding values
	# of the distribution function.
	if (!is.null(kappa0)){
		k <- floor(10^digits*kappa0)
	}else{
		k <- ceiling(10^digits*qnorm((1 + pv[1])/2))
	}
	av <- c(pmax(u - 10^(-digits)*k*qv, 0), 0)
	av[n] <- av[n-1]
	bv <- c(1,pmin(u + 10^(-digits)*k*qv, 1))
	bv[1] <- bv[2]
	q <- NoeRecursion(av,bv)
	tmp <- (pv > q)
	MK[tmp,1] <- pmax(MK[tmp,1],k)
	MK[!tmp,2] <- pmin(MK[!tmp,2],k)
	MP[tmp,1] <- pmax(MP[tmp,1],q)
	MP[!tmp,2] <- pmin(MP[!tmp,2],q)
	print(cbind(10^(-digits)*MK,MP))
	for (j in 1:np){
		while (MK[j,2] == Inf){
			k <- 2*k
			av <- c(pmax(u - 10^(-digits)*k*qv, 0), 0)
			av[n] <- av[n-1]
			bv <- c(1,pmin(u + 10^(-digits)*k*qv, 1))
			bv[1] <- bv[2]
			q <- NoeRecursion(av,bv)
			tmp <- (pv > q)
			MK[tmp,1] <- pmax(MK[tmp,1],k)
			MK[!tmp,2] <- pmin(MK[!tmp,2],k)
			MP[tmp,1] <- pmax(MP[tmp,1],q)
			MP[!tmp,2] <- pmin(MP[!tmp,2],q)
			print(cbind(10^(-digits)*MK,MP))
		}
		while (MK[j,2] > MK[j,1] + 1){
			pdiff <- MP[j,2] - MP[j,1]
			if (pdiff > 0){
				lambda <- (pv[j] - MP[j,1])/pdiff
				lambda <- (1/2 + 7*lambda)/8
				k <- max(MK[j,1]+1,
						 min(MK[j,2]-1,
						 	floor((1-lambda)*MK[j,1] + lambda*MK[j,2])))
			}else{
				k <- floor((MK[j,1] + MK[j,2])/2)
			}
			av <- c(pmax(u - 10^(-digits)*k*qv, 0), 0)
			av[n] <- av[n-1]
			bv <- c(1,pmin(u + 10^(-digits)*k*qv, 1))
			bv[1] <- bv[2]
			q <- NoeRecursion(av,bv)
			tmp <- (pv > q)
			MK[tmp,1] <- pmax(MK[tmp,1],k)
			MK[!tmp,2] <- pmin(MK[!tmp,2],k)
			MP[tmp,1] <- pmax(MP[tmp,1],q)
			MP[!tmp,2] <- pmin(MP[!tmp,2],q)
			print(cbind(10^(-digits)*MK,MP))
		}
	}
	return(10^(-digits)*MK[,2])
}

SP_TestCritval <- function(n,kappav,pv=c(0.5,0.9,0.95,0.99),
						   gamma=1,
						   mcsim=10000)
	# Test empirically whether given critical values
	# kappav could be pv-quantiles of the
	# Stepanova-Pavlenko test statistic.
{
	u <- (1:(n-1))/n
	v <- u*(1-u)
	qv <- sqrt(v*log(log(gamma/v))/n)
	np <- length(pv)
	H <- rep(0,np)
	for (sim in 1:mcsim){
		U <- sort(runif(n))
		Ts <- max(max((U[2:n] - u)/qv),
				  max((u - U[1:(n-1)])/qv))
		H <- H + (Ts <= kappav)
	}
	Table <- cbind('kappa'=kappav,
				   'p'=pv,
				   'p.emp'=H/mcsim,
				   'p-value'=2*pmin(pbinom(H,mcsim,pv),
				   				 pbinom(mcsim-H,mcsim,1-pv)))
	print(Table,quote=FALSE)
}


# Berk-Jones tests ----

BJ_CDF <- function(n, s=1, kappav)
	# Computes the distribution function of the
	# Berk-Jones test statistic at the points specified
	# by kappav. The test statistic is given by
	#    sup_{t in (0,1)} n K_s(Ghat(t),t) ,
	# where Ghat is the empricial distribution function
	# of n independent random variables with uniform
	# distribution on [0,1].
{
	nk <- length(kappav)
	CDF <- rep(1,nk)
	for (j in 1:nk){
		AB <- BJ_band(n,s,kappav[j])
		CDF[j] <- NoeRecursion(AB[2:(n+1),1],AB[1:n,2])
	}
	return(CDF)
}

BJ_Critval <- function(n, s=1, pv=c(0.5,0.9,0.95,0.99),
					   digits=5,
					   kappa0=NULL)
	# Computes the (1 - alpha)-quantile of the
	# Berk-Jones test statistic
	#    sup_{t in (0,1)} n K_s(Ghat(t),t) ,
	# where Ghat is the empricial distribution function
	# of n independent random variables with uniform
	# distribution on [0,1].
{
	np <- length(pv)
	MK <- matrix(0,np,2)
	MK[,2] <- Inf
	dimnames(MK)[[1]] <- pv
	MP <- matrix(0,np,2)
	MP[,2] <- 1
	# MK[j,1] and MK[j,2] are a lower and an upper bound,
	# respectively, for 10^digits*kappav[j], while
	# MP[j,1] and MP[j,2] are the corresponding values
	# of the distribution function.
	if (!is.null(kappa0)){
		k <- floor(10^digits*kappa0)
	}else{
		k <- ceiling(10^digits*qchisq(pv[1],1)/2)
	}
	AB <- BJ_band(n,s,10^(-digits)*k)
	q <- NoeRecursion(AB[2:(n+1),1],AB[1:n,2])
	tmp <- (pv > q)
	MK[tmp,1] <- pmax(MK[tmp,1],k)
	MK[!tmp,2] <- pmin(MK[!tmp,2],k)
	MP[tmp,1] <- pmax(MP[tmp,1],q)
	MP[!tmp,2] <- pmin(MP[!tmp,2],q)
	print(cbind(10^(-digits)*MK,MP))
	for (j in 1:np){
		while (MK[j,2] == Inf){
			k <- 2*k
			AB <- BJ_band(n,s,10^(-digits)*k)
			q <- NoeRecursion(AB[2:(n+1),1],AB[1:n,2])
			tmp <- (pv > q)
			MK[tmp,1] <- pmax(MK[tmp,1],k)
			MK[!tmp,2] <- pmin(MK[!tmp,2],k)
			MP[tmp,1] <- pmax(MP[tmp,1],q)
			MP[!tmp,2] <- pmin(MP[!tmp,2],q)
			print(cbind(10^(-digits)*MK,MP))
		}
		while (MK[j,2] > MK[j,1] + 1){
			pdiff <- MP[j,2] - MP[j,1]
			if (pdiff > 0){
				lambda <- (pv[j] - MP[j,1])/pdiff
				lambda <- (1/2 + 7*lambda)/8
				k <- max(MK[j,1]+1,
						 min(MK[j,2]-1,
						 	floor((1-lambda)*MK[j,1] + lambda*MK[j,2])))
			}else{
				k <- floor((MK[j,1] + MK[j,2])/2)
			}
			AB <- BJ_band(n,s,10^(-digits)*k)
			q <- NoeRecursion(AB[2:(n+1),1],AB[1:n,2])
			tmp <- (pv > q)
			MK[tmp,1] <- pmax(MK[tmp,1],k)
			MK[!tmp,2] <- pmin(MK[!tmp,2],k)
			MP[tmp,1] <- pmax(MP[tmp,1],q)
			MP[!tmp,2] <- pmin(MP[!tmp,2],q)
			print(cbind(10^(-digits)*MK,MP))
		}
	}
	return(10^(-digits)*MK[,2])
}

BJ_TestCritval <- function(n,s=1,kappav,
						   pv=c(0.5,0.9,0.95,0.99),
						   mcsim=10000)
	# Test empirically whether given critical values
	# kappav could be pv-quantiles of the
	# Berk-Jones test statistic.
{
	u0 <- (0:(n-1))/n
	u1 <- (1:n)/n
	np <- length(pv)
	H <- rep(0,np)
	for (sim in 1:mcsim){
		U <- sort(runif(n))
		Ts <- n*max(max(KSUT(u0,U,s)),max(KSUT(u1,U,s)))
		H <- H + (Ts <= kappav)
	}
	Table <- cbind('kappa'=kappav,
				   'p'=pv,
				   'p.emp'=H/mcsim,
				   'p-value'=2*pmin(pbinom(H,mcsim,pv),
				   				 pbinom(mcsim-H,mcsim,1-pv)))
	print(Table,quote=FALSE)
}



# Duembgen-Wellner tests ----

DW_CDF <- function(n, s=1, nu=1, kappav)
	# Computes the distribution function of the
	# Duembgen-Wellner test statistic at the points specified
	# by kappav. The test statistic is given by
	#    sup_{t in (0,1)}
	#       (n K_s(Ghat(t), t) - C_nu(Ghat(t),t)),
	# where Ghat is the empricial distribution function
	# of n independent random variables with uniform
	# distribution on [0,1].
{
	nk <- length(kappav)
	CDF <- rep(1,nk)
	for (j in 1:nk){
		AB <- DW_band(n,s,nu,kappav[j])
		CDF[j] <- NoeRecursion(AB[2:(n+1),1],AB[1:n,2])
	}
	return(CDF)
}

DW_Critval <- function(n, s=1, nu=1,
					   pv=c(0.5,0.9,0.95,0.99),
					   digits=5,
					   kappa0=NULL)
	# Computes the (1 - alpha)-quantile of the
	# Berk-Jones test statistic
	#    sup_{t in (0,1)} n K_s(Ghat(t),t) ,
	# where Ghat is the empricial distribution function
	# of n independent random variables with uniform
	# distribution on [0,1].
{
	np <- length(pv)
	MK <- matrix(0,np,2)
	MK[,2] <- Inf
	dimnames(MK)[[1]] <- pv
	MP <- matrix(0,np,2)
	MP[,2] <- 1
	# MK[j,1] and MK[j,2] are a lower and an upper bound,
	# respectively, for 10^digits*kappav[j], while
	# MP[j,1] and MP[j,2] are the corresponding values
	# of the distribution function.
	if (!is.null(kappa0)){
		k <- floor(10^digits*kappa0)
	}else{
		k <- ceiling(10^digits*qchisq(pv[1],1)/2)
	}
	AB <- DW_band(n,s,nu,10^(-digits)*k)
	q <- NoeRecursion(AB[2:(n+1),1],AB[1:n,2])
	tmp <- (pv > q)
	MK[tmp,1] <- pmax(MK[tmp,1],k)
	MK[!tmp,2] <- pmin(MK[!tmp,2],k)
	MP[tmp,1] <- pmax(MP[tmp,1],q)
	MP[!tmp,2] <- pmin(MP[!tmp,2],q)
	print(cbind(10^(-digits)*MK,MP))
	for (j in 1:np){
		while (MK[j,2] == Inf){
			k <- 2*k
			AB <- DW_band(n,s,nu,10^(-digits)*k)
			q <- NoeRecursion(AB[2:(n+1),1],AB[1:n,2])
			tmp <- (pv > q)
			MK[tmp,1] <- pmax(MK[tmp,1],k)
			MK[!tmp,2] <- pmin(MK[!tmp,2],k)
			MP[tmp,1] <- pmax(MP[tmp,1],q)
			MP[!tmp,2] <- pmin(MP[!tmp,2],q)
			print(cbind(10^(-digits)*MK,MP))
		}
		while (MK[j,2] > MK[j,1] + 1){
			pdiff <- MP[j,2] - MP[j,1]
			if (pdiff > 0){
				lambda <- (pv[j] - MP[j,1])/pdiff
				lambda <- (1/2 + 7*lambda)/8
				k <- max(MK[j,1]+1,
						 min(MK[j,2]-1,
						 	floor((1-lambda)*MK[j,1] + lambda*MK[j,2])))
			}else{
				k <- floor((MK[j,1] + MK[j,2])/2)
			}
			AB <- DW_band(n,s,nu,10^(-digits)*k)
			q <- NoeRecursion(AB[2:(n+1),1],AB[1:n,2])
			tmp <- (pv > q)
			MK[tmp,1] <- pmax(MK[tmp,1],k)
			MK[!tmp,2] <- pmin(MK[!tmp,2],k)
			MP[tmp,1] <- pmax(MP[tmp,1],q)
			MP[!tmp,2] <- pmin(MP[!tmp,2],q)
			print(cbind(10^(-digits)*MK,MP))
		}
	}
	return(10^(-digits)*MK[,2])
}

DW_TestCritval <- function(n,s=1,nu=1,kappav,
						   pv=c(0.5,0.9,0.95,0.99),
						   mcsim=10000)
	# Test empirically whether given critical values
	# kappav could be the pv-quantiles of the
	# Duembgen-Wellner test statistic.
{
	u0 <- (0:(n-1))/n
	u1 <- (1:n)/n
	np <- length(pv)
	H <- rep(0,np)
	for (sim in 1:mcsim){
		U <- sort(runif(n))
		Ts <- max(max(n*KSUT(u0,U,s) - CNuUT(u0,U,nu)),
				  max(n*KSUT(u1,U,s) - CNuUT(u1,U,nu)))
		H <- H + (Ts <= kappav)
	}
	Table <- cbind('kappa'=kappav,
				   'p'=pv,
				   'p.emp'=H/mcsim,
				   'p-value'=2*pmin(pbinom(H,mcsim,pv),
				   				 pbinom(mcsim-H,mcsim,1-pv)))
	print(Table,quote=FALSE)
}


# Auxiliary functions ----

KSUT <- function(u,t,s)
	# u and t should be scalars or
	# vectors of the same length
{
	if (s == 0){
		kv <- rep(0,length(u))
		tmp <- (t==0)
		kv[tmp] <- -log(1 - u[tmp])
		tmp <- (t==1)
		kv[tmp] <- -log(u[tmp])
		tmp <- (t>0 & t<1)
		kv[tmp] <- t[tmp]*log(t[tmp]/u[tmp]) +
			(1 - t[tmp])*log((1 - t[tmp])/(1 - u[tmp]))
		return(kv)
	}
	if (s == 1){
		kv <- rep(0,length(t))
		tmp <- (u==0)
		kv[tmp] <- -log(1 - t[tmp])
		tmp <- (u==1)
		kv[tmp] <- -log(t[tmp])
		tmp <- (u>0 & u<1)
		kv[tmp] <- u[tmp]*log(u[tmp]/t[tmp]) +
			(1 - u[tmp])*log((1 - u[tmp])/(1 - t[tmp]))
		return(kv)
	}
	return((t^(1-s)*u^s +
				(1 - t)^(1 - s)*(1 - u)^s - 1)/s/(s-1))
}

CNuUT <- function(u,t,nu=1)
{
	C <- rep(0,max(length(u),length(t)))
	x <- pmax(u,t)
	tmp <- (x < 0.5)
	C[tmp] <- log(1 - log(1 - (2*x[tmp] - 1)^2))
	C[tmp] <- C[tmp] + nu*log(1 + C[tmp]^2)
	x <- pmin(u,t)
	tmp <- (x > 0.5)
	C[tmp] <- log(1 - log(1 - (2*x[tmp] - 1)^2))
	C[tmp] <- C[tmp] + nu*log(1 + C[tmp]^2)
	return(C)
}

BJ_band <- function(n,s=1,kappa,prec=10^(-7))
	# For a given parameter s in (0,2] and given
	# critical value kappa, this procedure computes
	# (approximately) the coefficients of the BJ-confidence
	# band resulting from the statistic
	#    T_s(F) := sup_{x in R} n K_s(Fhat(x), F(x)).
	# The inequality T_s(F) <= kappa is equivalent to
	#    AB[i+1,1] <= F(x) <= AB[i+1,2]
	# for 0 <= i <= n and X_{n;i} <= x < X_{n:i+1}.
	# In other words, the random vector U with components
	#    U[i] = F(X_{n:i})
	# satisfies
#    a[i] <= U[i] <= b[i]  for 1 <= i <= n,
# where a[i] := AB[i+1,1] and b[i] := AB[i,2].
{
	if (s <= 0){
		return("Sorry, the bands haven't been implemented for s <= 0.")
	}
	b <- rep(1,n+1)
	for (j in (n-1):0){
		b2 <- b[j+2]
		K2 <- n*KSUT(j/n,b2,s)
		b1 <- (j/n + b2)/2
		K1 <- n*KSUT(j/n,b1,s)
		while (K1 >= kappa){
			b2 <- b1
			K2 <- K1
			b1 <- (j/n + b2)/2
			K1 <- n*KSUT(j/n,b1,s)
		}
		while (K2 >= kappa && b2 - b1 > prec/n){
			bm <- (b1 + b2)/2
			Km <- n*KSUT(j/n,bm,s)
			if (Km < kappa){
				b1 <- bm
				K1 <- Km
			}else{
				b2 <- bm
				K2 <- Km
			}
		}
		b[j+1] <- b2
	}
	a <- 1 - b[(n+1):1]
	AB <- cbind('a'=a,'b'=b)
	dimnames(AB)[[1]] <- 0:n
	return(AB)
}

DW_band <- function(n,s=1,nu=1,kappa,prec=10^(-7))
	# For a given parameter s in (0,2] and given
	# critical value kappa, this procedure computes
	# (approximately) the coefficients of the DW-confidence
	# band resulting from the statistic
	#    T_{s,nu}(F) := sup_{x in R}
	#       (n K_s(Fhat(x), F(x)) - C_nu(Fhat(x),F(x))).
	# The inequality T_{s,nu}(F) <= kappa is equivalent to
	#    AB[i+1,1] <= F(x) <= AB[i+1,2]
	# for 0 <= i <= n and X_{n;i} <= x < X_{n:i+1}.
	# In other words, the random vector U with components
	#    U[i] = F(X_{n:i})
# satisfies
#    a[i] <= U[i] <= b[i]  for 1 <= i <= n,
# where a[i] := AB[i+1,1] and b[i] := AB[i,2].
{
	if (s <= 0){
		return("Sorry, the bands haven't been implemented for s <= 0.")
	}
	b <- rep(1,n+1)
	for (j in (n-1):0){
		b2 <- b[j+2]
		K2 <- n*KSUT(j/n,b2,s) - CNuUT(j/n,b2,nu)
		b1 <- (j/n + b2)/2
		K1 <- n*KSUT(j/n,b1,s) - CNuUT(j/n,b1,nu)
		while (K1 >= kappa){
			b2 <- b1
			K2 <- K1
			b1 <- (j/n + b2)/2
			K1 <- n*KSUT(j/n,b1,s) - CNuUT(j/n,b1,nu)
		}
		while (K2 >= kappa && b2 - b1 > prec/n){
			bm <- (b1 + b2)/2
			Km <- n*KSUT(j/n,bm,s) - CNuUT(j/n,bm,nu)
			if (Km < kappa){
				b1 <- bm
				K1 <- Km
			}else{
				b2 <- bm
				K2 <- Km
			}
		}
		b[j+1] <- b2
	}
	a <- 1 - b[(n+1):1]
	AB <- cbind('a'=a,'b'=b)
	dimnames(AB)[[1]] <- 0:n
	return(AB)
}

