source("NoeRecursion.R")

# Description ----

# NoeRecursion(alpha,beta):
# For two vectors alpha and beta of real numbers
#   alpha[1] <= alpha[2] <= ... <= alpha[n],
#   beta[1]  <= beta[2]  <= ... <= beta[n]
# such that alpha[i] < beta[i], this function computes
# the probability Q that the order statistics
#   U[1] < U[2] < ... < U[n]
# of n independent random variables with uniform distribution
# on [0,1] satisfy
#    alpha[i] < U[i] <= beta[i]  for i=1,2,...,n .

# Test 1: a single order statistic ----

# If
#    alpha[i] = 1{i >= k}*a ,
#    beta[i]  = b + 1{i > k}*(1-b)
# for 0 <= a < b <= 1, then the probability Q equals the probability that
#    a < U[k] <= b ,
# and this is equal to the probability of (a,b] under the beta
# distribution with parameters k and k+1-k.

n <- 300
k <- 100
a <- 0.4
b <- 0.6
alpha <- c(rep(0,k-1),rep(a,n-k+1))
beta <- c(rep(b,k),rep(1,n-k))
cbind(alpha,beta)

Q <- NoeRecursion(alpha,beta)
Q
pbeta(b,k,n+1-k) - pbeta(a,k,n+1-k)

# Test 2: a multinomial probability ----

# If
#    alpha[i] = 1{i > j}*c2 ,
#    beta[i]  = c1 + 1{i > j}*(1 - c1)
# for real numbers 0 < c1 < c2 < 1, then the probability Q is the
# probability that in a standard uniform sample of size n,
#    j observations are in [0,c1]
# and
#    n-j observations are in (c2,1] ,
# and this probability equals
#    choose(n,j) * c1^j * (1 - c2)^(n-j) .

n <- 250
c1 <- 0.42
c2 <- 0.5
j <- 125
alpha <- c(rep(0,j),rep(c2,n-j))
beta <- c(rep(c1,j),rep(1,n-j))
Q <- NoeRecursion(alpha,beta)
Q
choose(n,j)*c1^j * (1 - c2)^(n-j)

# Test 3: Daniel's test ----

# For any number gamma in (0,1), the probability that
#   U[k] > gamma*k/n  for k=1,2,...,n
# is equal to 1 - gamma, see Daniel (19??).

n <- 500
gamma <- 0.25
alpha <- gamma*(1:n)/n
beta <- rep(1,n)
Q <- NoeRecursion(alpha,beta)
Q


# Test 4: Kolmogorov-Smirnov confidence bands ----

# True coverage probability of Kolmogorov-Smirnov confidence bands:

n <- 2000
kappa <- 1.3564/sqrt(n)
kappa
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

