source("ConfBands.R")

# First illustrations for s = 1 ----

par(cex=1.2,mai=c(0.5,0.5,0.05,0.05),mgp=c(2,1,0))

n <- 250

# Kolmogorov-Smirnov bands:

# Critical value (without factor sqrt(n)):
kappaKS <- 0.08520  # n=250
# Confidence band:
ABKS <- KS_band(n,kappa=kappaKS)

DisplayCBand(AB1=ABKS)
DisplayCBand(AB1=ABKS,reference="g")
DisplayCBand(AB1=ABKS,reference="g",centered=TRUE)

# Stepanova-Pavleko (2018) bands:

gamma_a <- 1
kappaSP_a <- 4.75708 # n=250,  gamma=1
ABSP_a <- SP_band(n,kappa=kappaSP_a,gamma=gamma_a)

# Modified version:
gamma_b <- 1.4579
kappaSP_b <- 4.06479 # n=250,  gamma=1.4579
ABSP_b <- SP_band(n,kappa=kappaSP_b,gamma=gamma_b)

DisplayCBand(AB1=ABSP_a)
DisplayCBand(AB1=ABSP_a,reference="g")
DisplayCBand(AB1=ABSP_a,reference="g",centered=TRUE)
DisplayCBand(AB1=ABSP_a,AB2=ABSP_b,AB3=ABKS,reference="g",
			 centered=TRUE)

# Berk-Jones-Owen bands:

kappaBJ <- 5.56533 # n=250, s=1.0
ABBJ <- BJ_band(n,s=1,kappa=kappaBJ)

DisplayCBand(AB1=ABBJ)
DisplayCBand(AB1=ABBJ,reference="g")
DisplayCBand(AB1=ABBJ,AB2=ABSP_a,AB3=ABKS,reference="g")
DisplayCBand(AB1=ABBJ,AB2=ABSP_a,AB3=ABKS,
			 reference="g",centered=TRUE)

# Duembgen-Wellner bands:

kappaDW <- 4.61583 # n=250, s=1.0
ABDW <- DW_band(n,s=1,kappa=kappaDW)

DisplayCBand(AB1=ABDW)
DisplayCBand(AB1=ABDW,reference="g")
DisplayCBand(AB1=ABDW,AB2=ABSP_a,
			 colbg=c('black','green'),
			 lwdbg=c(2,2),
			 colfg=c('black','black'),
			 lwdfg=c(2,2),
			 reference="g",
			 xlab='',ylab='')
DisplayCBand(AB1=ABDW,reference="g",centered=TRUE)
DisplayCBand2(AB1=ABDW,AB2=ABSP_a,AB3=ABKS,
			  colbg=c('black','green','orange'),
			  lwdbg=c(2,2,2),
			  colfg=c('black','black','darkgray'),
			  lwdfg=c(2,2,2),
			  reference="g",differences=FALSE,
			  xlab='',ylab='',main='')
legend("topright",c('KS','DW','BJO'),
	   lwd=c(2,2,2),
	   col=c("orange","black","green"))
DisplayCBand2(AB1=ABDW,AB2=ABBJ,AB3=ABKS,
			  colbg=c('black','green','orange'),
			  lwdbg=c(1,1,1),
			  lwdfg=c(2,2,2),
			  reference="g",differences=TRUE,
			  xlab='',ylab='')


# Impact of s ----

n <- 500
n <- 2000

kappaDW06 <- 5.12003 # n=500,s=0.6
kappaDW10 <- 4.61260 # n=500,s=1.0
kappaDW14 <- 5.48632 # n=500,s=1.4

ABDW06 <- DW_band(n,s=0.6,kappa=kappaDW06)
ABDW10 <- DW_band(n,s=1.0,kappa=kappaDW10)
ABDW14 <- DW_band(n,s=1.4,kappa=kappaDW14)
DisplayCBand2(AB1=ABDW10,AB2=ABDW14,AB3=ABDW06,
			  colbg=c('black','yellow','orange'),
			  lwdbg=c(1,1,1),
			  lwdfg=c(2,2,2),
			  reference="g",differences=FALSE,
			  xlab='',ylab='',main='')
DisplayCBand2(AB1=ABDW10,AB2=ABDW14,AB3=ABDW06,
			  colbg=c('black','yellow','orange'),
			  lwdbg=c(1,1,1),
			  lwdfg=c(2,2,2),
			  reference="g",differences=TRUE,
			  xlab='',ylab='',main='')

n <- 2000
sv <- c(0.6,0.8,1.0,1.2,1.4)
kappav <- c(5.00564,4.67497,4.61223,4.79777,5.38582)
colv <- c('forestgreen','blue','black','red','orange')
DisplayDWBandsS(n,sv[5:1],kappav[5:1],j0=3,
				col=colv[5:1],reference="g",
				xlab='',ylab='',main='')
DisplayDWBandsS(n,sv[5:1],kappav[5:1],j0=3,
				col=colv[5:1],reference="g",
				differences=TRUE,
				xlab='',ylab='',main='')
