rm(list=ls(all=TRUE))
library(unmarked)
library(ggplot2)
setwd("C:/Users/npwil/OneDrive/Desktop/School/Grad School/Thesis/Data and Analysis/2022 Summer eDNA/Grayling-eDNA R")
occ2014 <- read.csv("edna_occ_2014.csv",header=T)


#<- read.csv(file.choose(),header=T)


attach(occ2014)
y <- occ2014[,2:4]
y

siteCovs <- data.frame(occ2014[,"cat"])
names(siteCovs)<-"cat"
area1 <- data.frame(cbind(a1 = occ2014$area, a2 = occ2014$area, a3 = occ2014$area))

obsCovs <- list(area=area1[,c("a1", "a2", "a3")])

ts <- unmarkedFrameOccu(y = y, siteCovs=siteCovs,obsCovs=obsCovs)

obsCovs(ts) <- scale(obsCovs(ts))
summary(ts)


ts1 <- occu(~1 ~1,ts)
ts1

backTransform(ts1, 'state')
backTransform(ts1, 'det')


##############################################for each category

ts2 <- occu(~1 ~cat,ts)
ts2

backTransform(linearComb(ts2, c(1,0,0),type = 'state'))
backTransform(linearComb(ts2, c(0,1,0),type = 'state'))
backTransform(linearComb(ts2, c(0,0,1),type = 'state'))
############################################################



ts3 <- occu(~area ~1,ts)
ts3



backTransform(linearComb(ts3, c(1,0),type = 'det'))


min(obsCovs(ts))



newData <- data.frame(area = rep(seq(from = min(obsCovs(ts)), to = max(obsCovs(ts)), length.out = 10)))

z <- predict(ts3,type='det',newdata=newData,appendData=T)
z
plot(z$area,z$Predicted)

z <- predict(ts3,type='det',newdata=newData,appendData=T)

z <- cbind(z,area3 = (z$area*sd(as.matrix(area1))+mean(as.matrix(area1)))*0.000001)

theme_set(theme_bw(base_size = 24)) 

det.plot <- ggplot(z, aes(area3, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size=2) +
  labs(x = "Area (square km)", y = "Detection probability")
det.plot
#Cat<-read.csv(file.choose(),header=T)

#attach(Cat)
#Cat
#ggplot(data=Cat, aes(x=Categories, y=Occupancy, legend=NULL )) + geom_bar(stat="identity", fill="grey", color="black")+ geom_errorbar(aes(ymin=Occupancy-SE, ymax=Occupancy+SE))
