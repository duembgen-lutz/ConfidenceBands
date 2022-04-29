source("NoeRecursion.R")

n <- 30
k <- 7
a <- 0.4
b <- 0.6
alpha <- c(rep(0,k-1),rep(a,n-k+1))
beta <- c(rep(b,k),rep(1,n-k))
cbind(alpha,beta)

Q <- NoeRecursion(alpha,beta)
Q
pbeta(b,k,n+1-k) - pbeta(a,k,n+1-k)

n <- 25
c1 <- 0.4
c2 <- 0.5
j <- 10
alpha <- c(rep(0,j),rep(c2,n-j))
beta <- c(rep(c1,j),rep(1,n-j))
Q <- NoeRecursion(alpha,beta)
Q
choose(n,j)*c1^j * (1 - c2)^(n-j)

n <- 100
gamma <- 0.25
alpha <- gamma*(1:n)/n
beta <- rep(1,n)
Q <- NoeRecursion(alpha,beta)
Q

n <- 100
gamma <- 0.25
alpha <- rep(0,n)
beta <- 1 - gamma*(n:1)/n
Q <- NoeRecursion(alpha,beta)
Q

n <- 8000
kappa <- 1.3564/sqrt(n)
nkappa
alpha <- pmax((1:n)/n - kappa, 0)
beta <- pmin((0:(n-1))/n + kappa, 1)
Q <- NoeRecursion(alpha,beta,countdown=TRUE)
Q

mcsim <- 100000
Tv <- rep(0,mcsim)
for (s in 1:mcsim){
	U <- sort(runif(n))
	Tv[s] <- max(U - (0:(n-1))/n, (1:n)/n - U)
}
mean(Tv <= kappa)
binom.test(sum(Tv <= kappa),mcsim,Q)

