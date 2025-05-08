library(riojaPlot)
library(readxl)
library(dplyr)
library(rioja)
library(analogue)
library(vegan)
library(ggpalaeo)
library(tidyverse)
library(mgcv)
setwd("C:\\Users\\Andrew\\Desktop\\Kathleen\\")
coreinput<-read.csv(file="spider.csv", row.names=1)

BCI<-coreinput

chron <- cbind(BCI[,1] , BCI[,2])
colnames(chron) <- c("Depth","Year")
chron=as.data.frame(chron)

BCI<-BCI[,-cbind(1:2)] #first two columns
core <- BCI / rowSums(BCI) * 100

# Fortin et al. 2015 used filtered >=2%-in-2-lakes spp criteria in their analyses
# this code is to retain only taxa that are >=2% in 2 intervals in the training set
N <- 2
M <- 2
i <- colSums(core >= M) >= N
spp_red <- core[, i, drop = FALSE]
# 61 chironomid taxa remain in training set
ncol(spp_red)

#this is for spider, check for others
fos<-spp_red[,-cbind(35:38)] #remove undifferentiated spider lake

nms <- colnames(fos)
nms3<-gsub("\\.", " ", nms)




depth<-as.numeric(chron$Depth)



#diversity

H<-diversity(core)
H<-as.data.frame(H)
J <- H/log(specnumber(core))
J<-as.data.frame(J)
#DCA and PCA
dca <- decorana(fos, iweigh=1)
sc1 <- scores(dca, display="sites", choices=1:2)
DC1<-as.data.frame(sc1)

H1.pca <- rda(fos)
sc.pca <-scores(H1.pca,display="sites", choices=1:1)
PC1<-as.data.frame(sc.pca)




isotopes_spider <-read.csv(file="spider_isotopes.csv", row.names=1)
#recon_left <-read.csv(file="reconstruction_left.csv", row.names=1)

iso1<-isotopes_spider[,-cbind(1, 2, 3, 5)] #first two columns
iso2<-isotopes_spider[,-cbind(1, 2, 4, 6)]#c13 and n15


diss <- dist(core/100)
clust <- chclust(diss)
bstick(clust, 5)


#gam1 <- mgcv::gam(chron$Year ~ s(re04$Comp02, k=5))
tiff("spider.tiff", width = 6, height = 4, units = 'in', res = 300)


rp1 <- riojaPlot(fos, chron, 
                 ymin = 1800, ymax=2023, yinterval = 20, 
                 yvar.name = "Year",
                 y.rev=FALSE,
                 yTop=0.6,
                 yBottom=0.20,
                 scale.percent=TRUE, plot.groups=FALSE, do.clust = TRUE,
                 plot.poly=FALSE, plot.line=FALSE, lwd.bar=3,
                 plot.clust=TRUE,  labels.italicise=TRUE,
                 plot.cumul=TRUE, cex.cumul=0.6, xlabels=nms3, labels.break.long=FALSE,
                 srt.xlabel=45, plot.bar=TRUE,
                 tcl=-0.1, cex.yaxis=0.5, cex.xlabel=0.4,  cex.ylabel=0.5, col.bar="black",
                 cex.xaxis=0.5, xRight = 0.7, plot.line=FALSE)

#par(fig=c(0.8,1,0,1), new=TRUE)


rp2 <- riojaPlot(iso1, isotopes_spider, yvar.name="Year",
                 riojaPlot=rp1,xGap = 0.01,
                 xRight=0.84, 
                 scale.minmax=FALSE, plot.bar=F, las.xaxis=2, cex.xlabel=0.5,
                 plot.line=T, plot.symb=TRUE, symb.cex=0.5)

rp3 <- riojaPlot(DC1, chron, yvar.name="Year",
                 riojaPlot=rp2,xGap = 0.01,
                 xRight=0.94, las.xaxis=2,
                 scale.minmax=FALSE, plot.bar=F, cex.xlabel=0.5,
                 plot.line=T, plot.symb=TRUE, symb.cex=0.5)

rp4 <- riojaPlot(H, chron, yvar.name="Year",
                 riojaPlot=rp3,xGap = 0.01,
                 xRight=0.98, las.xaxis=2,
                 scale.minmax=FALSE, plot.bar=F, cex.xlabel=0.5,
                 plot.line=T, plot.symb=TRUE, symb.cex=0.5)

rp3a <- addRPClustZone(rp3, clust, nZone=2, xRight=0.94, col="black")
rp3b <- addRPClustZone(rp4, clust, nZone=3, xRight=0.94, col="#595959", lty=2)


dev.off()



#
#
#


coreinput<-read.csv(file="settle.csv", row.names=1)

BCI<-coreinput

chron <- cbind(BCI[,1] , BCI[,2])
colnames(chron) <- c("Depth","Year")
chron=as.data.frame(chron)

BCI<-BCI[,-cbind(1:2)] #first two columns
core <- BCI / rowSums(BCI) * 100

# Fortin et al. 2015 used filtered >=2%-in-2-lakes spp criteria in their analyses
# this code is to retain only taxa that are >=2% in 2 intervals in the training set
N <- 2
M <- 2
i <- colSums(core >= M) >= N
spp_red <- core[, i, drop = FALSE]
# 61 chironomid taxa remain in training set
ncol(spp_red)

#this is for settle, check for others
fos<-spp_red[,-cbind(36:41)] #remove undifferentiated settle lake

nms <- colnames(fos)
nms3<-gsub("\\.", " ", nms)




depth<-as.numeric(chron$Depth)



#diversity

H<-diversity(core)
H<-as.data.frame(H)
J <- H/log(specnumber(core))
J<-as.data.frame(J)
#DCA and PCA
dca <- decorana(fos, iweigh=1)
sc1 <- scores(dca, display="sites", choices=1:2)
DC1<-as.data.frame(sc1)

H1.pca <- rda(fos)
sc.pca <-scores(H1.pca,display="sites", choices=1:1)
PC1<-as.data.frame(sc.pca)




isotopes_spider <-read.csv(file="settle_isotopes.csv", row.names=1)
#recon_left <-read.csv(file="reconstruction_left.csv", row.names=1)

iso1<-isotopes_spider[,-cbind(1, 2, 3, 5)]#  C and N
iso2<-isotopes_spider[,-cbind(1, 2, 4, 6)]#c13 and n15


diss <- dist(core/100)
clust <- chclust(diss)
bstick(clust, 5)



tiff("settle2.tiff", width = 6, height = 4, units = 'in', res = 300)


rp1 <- riojaPlot(fos, chron, 
                 ymin = 1780, ymax=2023, yinterval = 20, 
                 yvar.name = "Year",
                 y.rev=FALSE,
                 yTop=0.6,
                 yBottom=0.20,
                 scale.percent=TRUE, plot.groups=FALSE, do.clust = TRUE,
                 plot.poly=FALSE, plot.line=FALSE, lwd.bar=3,
                 plot.clust=TRUE,  labels.italicise=TRUE,
                 plot.cumul=TRUE, cex.cumul=0.6, xlabels=nms3, labels.break.long=FALSE,
                 srt.xlabel=45, plot.bar=TRUE,
                 tcl=-0.1, cex.yaxis=0.5, cex.xlabel=0.4,  cex.ylabel=0.5, col.bar="black",
                 cex.xaxis=0.5, xRight = 0.7, plot.line=FALSE)

#par(fig=c(0.8,1,0,1), new=TRUE)


rp2 <- riojaPlot(iso2, isotopes_spider, yvar.name="Year",
                 riojaPlot=rp1,xGap = 0.01,
                 xRight=0.84, 
                 scale.minmax=FALSE, plot.bar=F, las.xaxis=2, cex.xlabel=0.5,
                 plot.line=T, plot.symb=TRUE, symb.cex=0.5)

rp3 <- riojaPlot(DC1, chron, yvar.name="Year",
                 riojaPlot=rp2,xGap = 0.01,
                 xRight=0.94, las.xaxis=2,
                 scale.minmax=FALSE, plot.bar=F, cex.xlabel=0.5,
                 plot.line=T, plot.symb=TRUE, symb.cex=0.5)

rp4 <- riojaPlot(H, chron, yvar.name="Year",
                 riojaPlot=rp3,xGap = 0.01,
                 xRight=0.98, las.xaxis=2,
                 scale.minmax=FALSE, plot.bar=F, cex.xlabel=0.5,
                 plot.line=T, plot.symb=TRUE, symb.cex=0.5)

rp3a <- addRPClustZone(rp3, clust, nZone=2, xRight=0.94, col="black")
rp3b <- addRPClustZone(rp4, clust, nZone=3, xRight=0.94, col="#595959", lty=2)


dev.off()

#
#
#
#
#
#

coreinput<-read.csv(file="Chocolate.csv", row.names=1)

BCI<-coreinput

chron <- cbind(BCI[,1] , BCI[,2])
colnames(chron) <- c("Depth","Year")
chron=as.data.frame(chron)

BCI<-BCI[,-cbind(1:2)] #first two columns
core <- BCI / rowSums(BCI) * 100

# Fortin et al. 2015 used filtered >=2%-in-2-lakes spp criteria in their analyses
# this code is to retain only taxa that are >=2% in 2 intervals in the training set
N <- 2
M <- 2
i <- colSums(core >= M) >= N
spp_red <- core[, i, drop = FALSE]
# 61 chironomid taxa remain in training set
ncol(spp_red)

#this is for choc, check for others
fos<-spp_red[,-cbind(32:36)] #remove undifferentiated choc lake

nms <- colnames(fos)
nms3<-gsub("\\.", " ", nms)




depth<-as.numeric(chron$Depth)



#diversity

H<-diversity(core)
H<-as.data.frame(H)
J <- H/log(specnumber(core))
J<-as.data.frame(J)
#DCA and PCA
dca <- decorana(fos, iweigh=1)
sc1 <- scores(dca, display="sites", choices=1:2)
DC1<-as.data.frame(sc1)

H1.pca <- rda(fos)
sc.pca <-scores(H1.pca,display="sites", choices=1:1)
PC1<-as.data.frame(sc.pca)




loi <-read.csv(file="loi-chocolate.csv", row.names=1)
#recon_left <-read.csv(file="reconstruction_left.csv", row.names=1)

loi1<-loi[,-cbind(1, 2)]#  C and N
colnames(loi1) <- c("%Org","%Carb")


diss <- dist(core/100)
clust <- chclust(diss)
bstick(clust, 5)

#the broken stick looks funky because of the dating

tiff("choc.tiff", width = 6, height = 4, units = 'in', res = 300)


rp1 <- riojaPlot(fos, chron, 
                 ymin = 1780, ymax=2023, yinterval = 20, 
                 yvar.name = "Year",
                 y.rev=FALSE,
                 yTop=0.6,
                 yBottom=0.20,
                 scale.percent=TRUE, plot.groups=FALSE, do.clust = TRUE,
                 plot.poly=FALSE, plot.line=FALSE, lwd.bar=3,
                 plot.clust=TRUE,  labels.italicise=TRUE,
                 plot.cumul=TRUE, cex.cumul=0.6, xlabels=nms3, labels.break.long=FALSE,
                 srt.xlabel=45, plot.bar=TRUE,
                 tcl=-0.1, cex.yaxis=0.5, cex.xlabel=0.4,  cex.ylabel=0.5, col.bar="black",
                 cex.xaxis=0.5, xRight = 0.7, plot.line=FALSE)


rp2 <- riojaPlot(loi1, loi, yvar.name="Year",
                 riojaPlot=rp1,xGap = 0.01,
                 xRight=0.84, 
                 scale.minmax=FALSE, plot.bar=F, las.xaxis=2, cex.xlabel=0.5,
                 plot.line=T, plot.symb=TRUE, symb.cex=0.5)

rp3 <- riojaPlot(DC1, chron, yvar.name="Year",
                 riojaPlot=rp2,xGap = 0.01,
                 xRight=0.94, las.xaxis=2,
                 scale.minmax=FALSE, plot.bar=F, cex.xlabel=0.5,
                 plot.line=T, plot.symb=TRUE, symb.cex=0.5)

rp4 <- riojaPlot(H, chron, yvar.name="Year",
                 riojaPlot=rp3,xGap = 0.01,
                 xRight=0.98, las.xaxis=2,
                 scale.minmax=FALSE, plot.bar=F, cex.xlabel=0.5,
                 plot.line=T, plot.symb=TRUE, symb.cex=0.5)

rp3a <- addRPClustZone(rp3, clust, nZone=2, xRight=0.94, col="black")
rp3b <- addRPClustZone(rp4, clust, nZone=3, xRight=0.94, col="#595959", lty=2)


dev.off()

#
#
#
#
#
#

