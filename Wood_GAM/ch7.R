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


#------------

## 7.8.1 time dependent
## copy functions from ?cox.pht in mgcv...

app <- function(x,t,to) {
  ## wrapper to approx for calling from apply...
  y <- if (sum(!is.na(x))<1) rep(NA,length(to)) else
    approx(t,x,to,method="constant",rule=2)$y
  if (is.factor(x)) factor(levels(x)[y],levels=levels(x)) else y
} ## app

tdpois <- function(dat,event="z",et="futime",t="day",
                   status="status1",id="id") {
  ## dat is data frame. id is patient id; et is event time; t is
  ## observation time; status is 1 for death 0 otherwise;
  ## event is name for Poisson response.
  if (event %in% names(dat)) warning("event name in use")
  require(utils) ## for progress bar
  te <- sort(unique(dat[[et]][dat[[status]]==1])) ## event times
  sid <- unique(dat[[id]])
  prg <- txtProgressBar(min = 0, max = length(sid), initial = 0,
                        char = "=",width = NA, title="Progress", style = 3)
  ## create dataframe for poisson model data
  dat[[event]] <- 0; start <- 1
  dap <- dat[rep(1:length(sid),length(te)),]
  for (i in 1:length(sid)) { ## work through patients
    di <- dat[dat[[id]]==sid[i],] ## ith patient's data
    tr <- te[te <= di[[et]][1]] ## times required for this patient
    ## Now do the interpolation of covariates to event times...
    um <- data.frame(lapply(X=di,FUN=app,t=di[[t]],to=tr))
    ## Mark the actual event...
    if (um[[et]][1]==max(tr)&&um[[status]]==1) um[[event]][nrow(um)] <- 1 
    um[[et]] <- tr ## reset time to relevant event times
    dap[start:(start-1+nrow(um)),] <- um ## copy to dap
    start <- start + nrow(um)
    setTxtProgressBar(prg, i)
  }
  close(prg)
  dap[1:(start-1),]
} ## tdpois

## model fitting...
library(survival)
data(pbc)
pbcseq$status1 <- as.numeric(pbcseq$status==2) ## deaths
pb <- tdpois(pbcseq) ## conversion
pb$tf <- factor(pb$futime) ## add factor for event time

b <- bam(z ~ tf - 1  +  trt + s(sqrt(protime)) + s(platelet) + 
           s(age) + s(bili) + s(albumin) + s(sqrt(ast)),
         family=poisson,data=pb,discrete=TRUE,nthreads=2)

chaz <- tapply(fitted(b),pb$id,sum) ## cum. hazard by subject
d <- tapply(pb$z,pb$id,sum) ## censoring indicator
mrsd <- d - chaz ## Martingale residuals
drsd <- sign(mrsd)*sqrt(-2*(mrsd + d*log(chaz))) ## deviance

te <- sort(unique(pb$futime)) ## event times
di <- pbcseq[pbcseq$id==25,] ## data for subject 25
## interpolate to te using app from ?cox.pht...
pd <- data.frame(lapply(X=di,FUN=app,t=di$day,to=te)) 
pd$tf <- factor(te)
X <- predict(b,newdata=pd,type="lpmatrix")
eta <- drop(X%*%coef(b)); H <- cumsum(exp(eta))
J <- apply(exp(eta)*X,2,cumsum)
se <- diag(J%*%vcov(b)%*%t(J))^.5
par(mfrow=c(1,2))
plot(stepfun(te,c(1,exp(-H))),do.points=FALSE,ylim=c(0.7,1),
     ylab="S(t)",xlab="t (days)",main="",lwd=2)
lines(stepfun(te,c(1,exp(-H+se))),do.points=FALSE)
lines(stepfun(te,c(1,exp(-H-se))),do.points=FALSE)
rug(pbcseq$day[pbcseq$id==25]) ## measurement times

er <- pbcseq[pbcseq$id==25,]
plot(er$day,er$ast);lines(te,pd$ast)

#-------
# ch7 exercises
# 1
library(gamair)
library(mgcv)
data(hubble)
hub.mod <- lm(y ~ x - 1, data=hubble)
summary(hub.mod)
k=1
b = gam(y~ s(x), data = hubble)
gam.check(b)
plot(b,residuals=TRUE,n=1000, cex = 5)
qq.gam(b, cex = 5)
h0 <- gam(y~x,data=hubble)
b
AIC(h0, b)


gam.check(b) # oh dear
h2 <- gam(y~s(x),data=hubble,family=quasi(var=mu))
gam.check(h2) # not great, but better
h2
# The residual plots for h1 are problematic: there is a clear relationship between
# the mean and the variance. Perhaps a quasi-likelihood approach might solve
# this. m2 does have somewhat better residual plots, although they are still not
# perfect. All evidence for departure from Hubble’s law has now vanished.

# 2
library(MASS)
mcycle = mcycle
plot(mcycle$times, mcycle$accel)
mcycle0 = gam(accel~s(times, k = 30), data = mcycle, method = "GCV.Cp")
plot(mcycle0,  residuals = T, n = 1000, cex = 5, se = FALSE)
mcycle0

mcycle_lin <- lm(accel~ poly(times,11), data = mcycle)
plot(mcycle_lin, residuals = TRUE, n=1000, cex = 5)
termplot(mcycle_lin, partial.resid = T, rug = T)

mcycle0_unpenal = gam(accel~s(times), data = mcycle, method = "GCV.Cp")
plot(mcycle0_unpenal,  residuals = T, n = 1000, cex = 5, se = FALSE)
mcycle0_unpenal

mcycle0_unpenal_cubic = gam(accel~s(times, bs = "cc"), data = mcycle, method = "GCV.Cp")
plot(mcycle0_unpenal_cubic,  residuals = T, n = 1000, cex = 5, se = FALSE)
mcycle0_unpenal_cubic

length(mcycle$times)
mcycle$weights[1:133] = 0
mcycle$weights[1:20] = 2
mcycle0_unpenal_cubic_weight = gam(accel~s(times, bs = "cc", weight = mcycle$weights), data = mcycle, method = "GCV.Cp")
plot(mcycle0_unpenal_cubic_weight,  residuals = T, n = 1000, cex = 5, se = FALSE)
mcycle0_unpenal_cubic_weight
