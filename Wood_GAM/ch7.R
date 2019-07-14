#---------
# Brain data
library(mgcv)
library(MASS)
sm <- smoothCon(s(times,k=10),data=mcycle,knots=NULL)[[1]]
## use it to fit a regression spline model...
beta <- coef(lm(mcycle$accel~sm$X-1))
with(mcycle,plot(times,accel)) ## plot data
times <- seq(0,60,length=200) ## create prediction times
## Get matrix mapping beta to spline prediction at ’times’
Xp <- PredictMat(sm,data.frame(times=times))
lines(times,Xp%*%beta) ## add smooth to plot

#----------

library(gss)
data(wesdr)
k=10
b <- gam(ret ~ s(dur,k=k) + s(gly,k=k) + s(bmi,k=k) + ti(dur,gly,k=k) + ti(dur,bmi,k=k) + ti(gly,bmi,k=k), select=TRUE, data=wesdr, family=binomial(), method="ML")
b
mgcv::plot.gam(b)
mgcv::vis.gam(b, se=1, colors="bw")

#-------------
# Chicago data
library(gamair)
data(chicago)

#------------
# Try Poisson distribution
ap0 <- gam(death~s(time,bs="cr",k=200)+pm10median+so2median+
             o3median+tmpd,data=chicago,family=poisson)
gam.check(ap0)

par(mfrow=c(2,1))
plot(ap0,n=1000) # n increased to make plot smooth
plot(ap0,residuals=TRUE,n=1000)

# Try non-linear smoothing terms
ap1<-gam(death ~ s(time,bs="cr",k=200)+s(pm10median,bs="cr")+
           s(so2median,bs="cr")+s(o3median,bs="cr")+s(tmpd,bs="cr"),
         data=chicago,family=poisson)
plot(ap1,residuals=TRUE,n=1000)

#-----------
# Single index model


lagard <- function(x,n.lag=6) {
  n <- length(x); X <- matrix(NA,n,n.lag)
  for (i in 1:n.lag) X[i:n,i] <- x[i:n-i+1] 
  X
}
dat <- list(lag=matrix(0:5,nrow(chicago),6,byrow=TRUE),
            death=chicago$death,time=chicago$time)
dat$pm10 <- lagard(chicago$pm10median)
dat$tmp <- lagard(chicago$tmpd)
dat$o3 <- lagard(chicago$o3median)

si <- function(theta,dat,opt=TRUE) {
  ## Return ML if opt==TRUE or fitted gam otherwise.
  alpha <- c(1,theta) ## alpha defined via unconstrained theta
  kk <- sqrt(sum(alpha^2)); alpha <- alpha/kk  ## ||alpha||=1
  o3 <- dat$o3%*%alpha; tmp <- dat$tmp%*%alpha
  pm10 <- dat$pm10%*%alpha ## re-weight lagged covariates
  b<- bam(dat$death~s(dat$time,k=200,bs="cr")+s(pm10,bs="cr")+
            te(o3,tmp,k=8),family=poisson) ## fit model
  cat(".") ## give user something to watch
  if (opt) return(b$gcv.ubre) else {
    b$alpha <- alpha  ## add alpha to model object
    b$J <- outer(alpha,-theta/kk^2) ## get dalpha_i/dtheta_j
    for (j in 1:length(theta)) b$J[j+1,j] <- b$J[j+1,j] + 1/kk
    return(b)
  }
} ## si

## WARNING: the next line takes around half an hour to run

f1 <- optim(rep(1,5),si,method="BFGS",hessian=TRUE,dat=dat)

apsi <- si(f1$par,dat,opt=FALSE)
apsi$alpha

#--------------
# 7.4.2 distributed lag...

apl <- bam(death~s(time,bs="cr",k=200)+te(pm10,lag,k=c(10,5))+
             te(o3,tmp,lag,k=c(8,8,5)),family=poisson,data=dat)

#------
## 7.7.2 Cairo temperature
data(cairo)
ctamm <- gamm(temp~s(day.of.year,bs="cc",k=20)+s(time,bs="cr"),
              data=cairo,correlation=corAR1(form=~1|year))
summary(ctamm$gam)
intervals(ctamm$lme,which="var-cov")
ctamm$gam$sig2/ctamm$gam$sp
plot(ctamm$gam,scale=0,pages=1)

REML <- rho <- 0.6+0:20/100
for (i in 1:length(rho)) {
  ctbam <- bam(temp~s(day.of.year,bs="cc",k=20)+s(time,bs="cr"),
               data=cairo,rho=rho[i])
  REML[i] <- ctbam$gcv.ubre
}
rho[REML==min(REML)]
