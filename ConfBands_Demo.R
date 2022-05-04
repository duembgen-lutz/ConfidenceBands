source("ConfBands.R")

n <- 100
n <- 250

kappaKS <- 0.13403 # n=100
kappaKS <- 0.08520 # n=250
ABKS <- KS_band(n,kappa=kappaKS)

DisplayCBand(AB=ABKS)
DisplayCBand(AB=ABKS,reference="g",
			 main='KS band')
DisplayCBand(AB=ABKS,reference="g",centered=TRUE,
			 main='KS band (centered)')

kappaSP <- 4.76746 # n=100
kappaSP <- 4.75708 # n=250

ABSP <- SP_band(n,kappa=kappaSP)

DisplayCBand(AB=ABSP)
DisplayCBand(AB=ABSP,reference="g",
			 main='SP band')
DisplayCBand(AB=ABSP,reference="g",centered=TRUE,
			 main='SP band (centered)')
DisplayCBand(AB=ABSP,AB3=ABKS,reference="g",centered=TRUE,
			 main='SP (black) and KS (red) bands (centered)')


kappaBJ <- 5.37660 # n=100
kappaBJ <- 5.56533 # n=250,s=1.0

ABBJ <- BJ_band(n,s=1,kappa=kappaBJ)

DisplayCBand(AB=ABBJ)
DisplayCBand(AB=ABBJ,reference="g",
			 main='BJO band')
DisplayCBand(AB=ABBJ,reference="g",
			 AB2=ABSP,AB3=ABKS,
			 main='BJO (black), SP (blue) and KS (red) bands')
DisplayCBand(AB=ABBJ,reference="g",centered=TRUE,
			 main='BJO band (centered)')
DisplayCBand(AB=ABBJ,reference="g",centered=TRUE,
			 AB2=ABSP,AB3=ABKS,
			 main='BJO (black), SP (blue) and KS (red) bands (centered)')

kappaDW <- 4.62316 # n=100
kappaDW <- 4.61583 # n=250

ABDW <- DW_band(n,s=1,kappa=kappaDW)

DisplayCBand(AB=ABDW)
DisplayCBand(AB=ABDW,reference="g",
			 main='DW band')
DisplayCBand(AB=ABDW,reference="g",
			 AB2=ABBJ,
			 AB3=ABSP,
			 main='DW (black), BJO (blue) and SP (red) bands')
DisplayCBand(AB=ABDW,reference="g",centered=TRUE,
			 main='DW confidence band (centered)')
DisplayCBand(AB=ABDW,reference="g",centered=TRUE,
			 AB2=ABBJ,
			 AB3=ABSP,
			 main='DW (black), BJO (blue) and SP (red) confidence bands (centered)')

# Impact of s:

n <- 250
kappaBJ05 <- 7.61203 # n=250,s=0.5
kappaBJ10 <- 5.56533 # n=250,s=1.0
kappaBJ15 <- 8.08549 # n=250,s=1.5
ABBJ05 <- BJ_band(n,s=0.5,kappa=kappaBJ05)
ABBJ10 <- BJ_band(n,s=1.0,kappa=kappaBJ10)
ABBJ15 <- BJ_band(n,s=1.5,kappa=kappaBJ15)
DisplayCBand(AB=ABBJ10,reference="g",
			 AB2=ABBJ05,
			 AB3=ABBJ15,
			 main='BJO bands for s = 0.5 (blue), 1.0 (black), 1.5 (red)')
DisplayCBand(AB=ABBJ10,reference="g",centered=TRUE,
			 AB2=ABBJ05,
			 AB3=ABBJ15,
			 main='Centered BJO bands for s = 0.5 (blue), 1.0 (black), 1.5 (red)')

which(ABBJ05[,1] > ABBJ10[,1])
which(ABBJ15[,1] > ABBJ10[,1])


n <- 250
kappaDW05 <- 5.82036 # n=250,s=0.5
kappaDW10 <- 4.61583 # n=250,s=1.0
kappaDW15 <- 6.18853 # n=250,s=1.5

n <- 500
kappaDW05 <- 5.67719 # n=500,s=0.5
kappaDW10 <- 4.61260 # n=500,s=1.0
kappaDW15 <- 6.07882 # n=500,s=1.5

n <- 1000
kappaDW05 <- 5.56261 # n=1000,s=0.5
kappaDW10 <- 4.61152 # n=1000,s=1.0
kappaDW15 <- 5.99012 # n=1000,s=1.5

ABDW05 <- DW_band(n,s=0.5,kappa=kappaDW05)
ABDW10 <- DW_band(n,s=1.0,kappa=kappaDW10)
ABDW15 <- DW_band(n,s=1.5,kappa=kappaDW15)
DisplayCBand(AB=ABDW10,reference="g",
			 AB2=ABDW05,
			 AB3=ABDW15,
			 main='DW bands')
legend("topleft",c('s=1.5','s=1.0','s=0.5'),
	   lwd=c(1,2,1),
	   col=c("red","black","blue"))
DisplayCBand(AB=ABDW10,reference="g",centered=TRUE,
			 AB2=ABDW05,
			 AB3=ABDW15,
			 main='Centered DW bands')
legend("topleft",c('s=1.5','s=1.0','s=0.5'),
	   lwd=c(1,2,1),
	   col=c("red","black","blue"))

which(ABDW05[,1] > ABDW10[,1])
which(ABDW15[,1] > ABDW10[,1])


# Data example ----

X <- read.table("Galaxy-Dat.txt",header=FALSE)
X <- sort(as.vector(X$V1))
n <- length(X)
n

# Various 95%-confidence bands
kappaKS <- 0.14779
ABKS <- KS_band(n,kappa=kappaKS)
kappaSP <- 4.76530
ABSP <- SP_band(n,kappa=kappaSP)
kappaBJ <- 5.33160 # s=1.0
ABBJ <- BJ_band(n,s=1.0,kappaBJ)
kappaDW <- 4.62493 # s=1.0
ABDW <- DW_band(n,s=1.0,kappa=kappaDW)

par(cex=1,mai=c(0.9,0.9,0.5,0.05))
DisplayCBand(AB=ABKS,X=X,xlim=c(9,35),
			 main='KS band')
DisplayCBand(AB=ABSP,X=X,xlim=c(9,35),
			 main='SP band')
DisplayCBand(AB=ABBJ,X=X,xlim=c(9,35),
			 main='BJO band')
DisplayCBand(AB=ABDW,X=X,xlim=c(9,35),
			 main='DW band')

# Comparisons:
DisplayCBand(AB=ABDW,X=X,xlim=c(9,35),
			 AB2=ABSP,AB3=ABKS,
			 centered=TRUE,
			 main='DW (black), SP (blue) and KS (red) bands centered')
DisplayCBand(AB=ABDW,X=X,xlim=c(9,35),
			 AB2=ABBJ,AB3=ABKS,
			 centered=TRUE,
			 main='DW (black), BJO (blue) and KS (red) bands centered')


