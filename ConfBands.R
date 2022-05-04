# For a given test statistic T(X,F) and a critical
# value kappa, the following procedures compute the
# coefficients AB of the resulting confidence band.
# That is, if X contains the order statistics
#    X[1] <= X[2] <= ... <= X[n] ,
# any distribution function F such that T(X,F) <= kappa
# satisfies
#    AB[i+1,1] <= F(x) <= AB[i+1,2]
#    for i = 0,1,...,n and X[i] <= x < X[i+1],
# where X[0] := -Inf and X[n+1] := Inf.


# Kolmogorov-Smirnov confidence band ----

KS_band <- function(n,kappa)
{
	u <- (0:n)/n
	AB <- cbind('a'=pmax(u-kappa,0),'b'=pmin(u+kappa,1))
	dimnames(AB)[[1]] <- 0:n
	return(AB)
}

# Stepanova-Pavlenko confidence band ----

SP_band <- function(n,kappa,gamma=1)
{
	u <- (1:(n-1))/n
	v <- u*(1-u)
	qnv <- sqrt(v*log(log(gamma/v)))/sqrt(n)
	a <- c(0,pmax(u - kappa*qnv, 0), 0)
	a[n+1] <- a[n]
	b <- c(1,pmin(u + kappa*qnv, 1),1)
	b[1] <- b[2]
	AB <- cbind('a'=a,'b'=b)
	dimnames(AB)[[1]] <- 0:n
	return(AB)
}

# Berk-Jones-Owen confidence band ----

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

# Duembgen-Wellner confidence band ----

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

# Graphical displays ----

DisplayCBand <- function(AB,X=NULL,xlim=NULL,col='black',
						 main=NULL,
						 centered=FALSE,
						 reference=c('uniform','gaussian'),
						 AB2=NULL,AB3=NULL,
						 col2='blue',col3='red')
	# Display the confidence band AB for a distribution
	# function F for a sample X. If X==NULL, then the
	# prodedure generates a "sample" with components
	#    X[i] = i/(n+1)         if reference=="uniform",
	#    X[i] = qnorm(i/(n+1))  if reference=="uniform".
	# If centered==FALSE, the confidence band is displayed
	# together with the empirical distribution function of X.
	# If centered==TRUE, one sees the difference of the
	# boundaries and the empirical distribution function.
{
	n <- dim(AB)[1]-1
	if (!is.null(X)){
		X <- sort(X)
	}else{
		reference <- match.arg(reference,
							   c('uniform','gaussian'))
		X <- (1:n)/(n+1)
		if (reference=="gaussian"){
			X <- qnorm(X)
		}
	}
	if (is.null(xlim)){
		xlim <- c(2*X[1] - X[2], 2*X[n] - X[n-1])
	}
	XX <- c(min(xlim[1],2*X[1]-X[n]),
			rep(X,each=2),
			max(xlim[2],2*X[n]-X[1]))
	FhatXX <- rep((0:n)/n,each=2)
	if (centered){
		ABCXX <- cbind(rep(AB[,1],each=2)-FhatXX,
					   rep(AB[,2],each=2)-FhatXX)
		ylim=range(ABCXX)
		if (!is.null(AB2)){
			AB2CXX <- cbind(rep(AB2[,1],each=2)-FhatXX,
							rep(AB2[,2],each=2)-FhatXX)
			ylim=range(ylim,AB2CXX)
		}
		if (!is.null(AB3)){
			AB3CXX <- cbind(rep(AB3[,1],each=2)-FhatXX,
							rep(AB3[,2],each=2)-FhatXX)
			ylim=range(ylim,AB3CXX)
		}
		plot(range(XX),rep(0,2),type='l',col=col,
			 xlim=xlim,ylim=ylim,main=main,
			 xlab=expression(italic(x)),
			 ylab=expression(italic(F(x)-F[emp](x))))
		if (!is.null(AB3)){
			lines(XX,AB3CXX[,1],col=col3)
			lines(XX,AB3CXX[,2],col=col3)
		}
		if (!is.null(AB2)){
			lines(XX,AB2CXX[,1],col=col2)
			lines(XX,AB2CXX[,2],col=col2)
		}
		lines(XX,ABCXX[,1],col=col,lwd=2)
		lines(XX,ABCXX[,2],col=col,lwd=2)
	}else{
		plot(XX,FhatXX,type='l',col=col,
			 xlim=xlim,main=main,
			 xlab=expression(italic(x)),
			 ylab=expression(italic(F(x))))
		if (!is.null(AB3)){
			lines(XX,rep(AB3[,1],each=2),col=col3)
			lines(XX,rep(AB3[,2],each=2),col=col3)
		}
		if (!is.null(AB2)){
			lines(XX,rep(AB2[,1],each=2),col=col2)
			lines(XX,rep(AB2[,2],each=2),col=col2)
		}
		lines(XX,rep(AB[,1],each=2),col=col,lwd=2)
		lines(XX,rep(AB[,2],each=2),col=col,lwd=2)
	}
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
