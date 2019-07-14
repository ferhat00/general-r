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
ap0 <- gam(death~s(time,bs="cr",k=200)+pm10median+so2median+
             o3median+tmpd,data=chicago,family=poisson)
gam.check(ap0)

par(mfrow=c(2,1))
plot(ap0,n=1000) # n increased to make plot smooth
plot(ap0,residuals=TRUE,n=1000)

ap1<-gam(death ~ s(time,bs="cr",k=200)+s(pm10median,bs="cr")+
           s(so2median,bs="cr")+s(o3median,bs="cr")+s(tmpd,bs="cr"),
         data=chicago,family=poisson)
plot(ap1,residuals=TRUE,n=1000)
