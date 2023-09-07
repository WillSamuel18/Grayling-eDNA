rm(list=ls(all=TRUE))
library(unmarked)
setwd("C:/Users/npwil/OneDrive/Desktop/School/Grad School/Thesis/Data and Analysis/2022 Summer eDNA/Grayling-eDNA R")


###load data
lwd_pdata <- read.csv("nmix_data.csv")

####create unmarked data object for n-mixture model
y <- as.matrix(lwd_pdata[,2:3])

sitecovs <- data.frame(x1 = rep(1,29))

m <- 29
n <- 2
obscovs <- list(x2 = round(matrix(runif(m * n), m, n)))


time<-as.factor(rep(c(1,2),29))
time
blgr.obs<-data.frame(time)


lwd <- unmarkedFramePCount(y=y, siteCovs=sitecovs, 
                           obsCovs=blgr.obs)   
lwd
summary(lwd)

####fit models

lwd_nmix_pois <- pcount(~1 ~1, lwd, K=200) 
lwd_nmix_pois_t <- pcount(~time ~1, lwd, K=200) 

lwd_nmix_pois
lwd_nmix_pois_t

#####model selection
fmlist <- fitList(timeinvar = lwd_nmix_pois, timevar = lwd_nmix_pois_t)
modSel(fmlist)

#####detection probabilities by observer 1=Nate, 2=Brian plus 95% CI's
timeinvar <- backTransform(linearComb(lwd_nmix_pois_t, c(1, 1), type="det"))
timevar <- backTransform(linearComb(lwd_nmix_pois_t, c(1, 2), type="det"))

timeinvar_CI <- confint(backTransform(linearComb(lwd_nmix_pois_t, c(1, 1), type="det")))
timevar_CI <- confint(backTransform(linearComb(lwd_nmix_pois_t, c(1, 2), type="det")))

#####plot it
require(ggplot2)
df <- data.frame(obs = c("Nate","Brian"), p = c(0.751,0.84), upCI = c(0.79,0.88), lowCI = c(0.71,0.78))
ggplot(df, aes(x = obs, y = p)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upCI, ymin = lowCI))

