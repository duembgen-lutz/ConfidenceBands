# For numerical accuracy (dealing with products and sums
# of very small positive numbers), our version of Noe's
# recursion uses log-probabilities internally.

LogSum <- function(la,lb)
# Reliable computation of log(exp(la) + exp(lb))
{
	return(pmax(la,lb) + log(1 + exp(-abs(la-lb))))
}
	
NoeRecursion <- function(alpha,beta,countdown=FALSE)
{
	if (sum(alpha >= beta) > 0){
		return(0)
	}
	n <- length(alpha)
	logfact <- c(0,cumsum(log(1:n)))
	gamma <- unique(sort(c(alpha,beta)))
	N <- length(gamma)-1
	logQ <- matrix(-Inf,n,N)
	if (countdown){
		print(paste(N,':',sep=''),quote=FALSE)
	}
	kmin <- 0
	ellmax <- 1
	while (ellmax < n && alpha[ellmax+1] <= gamma[1]){
		ellmax <- ellmax+1
	}
	logqf <- log(gamma[2] - gamma[1])
	logQ[1:ellmax,1] <- (1:ellmax)*logqf
	for (m in 2:N){
		if (countdown){
			print(paste(N+1-m,':',sep=''),quote=FALSE)
		}
		while (beta[kmin+1] < gamma[m+1]){
			kmin <- kmin+1
		}
		while (ellmax < n && alpha[ellmax+1] <= gamma[m]){
			ellmax <- ellmax+1
		}
		ellmin <- max(1,min(kmin,ellmax))
		logQ[ellmin:ellmax,m] <- logQ[ellmin:ellmax,m-1]
		logqf <- log(gamma[m+1] - gamma[m])
		if (kmin == 0){
			logQ[1:ellmax,m] <- LogSum(logQ[1:ellmax,m],
				(1:ellmax)*logqf)
		}
		k <- ellmin
		while (k < ellmax){
			tmpell <- (k+1):ellmax
			logQ[tmpell,m] <- LogSum(logQ[tmpell,m],
				(tmpell - k)*logqf +
					logfact[tmpell+1] - logfact[tmpell-k+1] -
					logfact[k+1] + logQ[k,m-1])
			k <- k+1
		}
	}
	return(exp(logQ[n,N]))
}