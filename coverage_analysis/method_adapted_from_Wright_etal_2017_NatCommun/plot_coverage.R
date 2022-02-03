## This script creates moving average plots for coverage data

rm(list=ls())
ls() 

#Load packages
library(ggplot2)
library(RcppRoll)

#Define a moving average function
movingaverage <- function (x, window) {
  ma <- roll_mean(x, window, fill=NA)
}
#For all P. reticulata analyses and P. wingei analyses that used the X. maculatus geneome windowsizeCI = 20 and windowsizeMovAv = 20 (as example here)
#For P. wingei analyses that used the P. wingei de novo geneome, we used windowsizeCI = 50 and windowsizeMovAv = 60. 
windowsizeCI = 20
windowsizeMovAv = 20

#Load data
#sexchromosome is "NC_024342.1" for P. reticulata genomes or "8" for X. maculatus genome and P. wingie de novo genome. The example here is based on the P. reticulata genome.
sexchromosome = "NC_024342.1"
#The example script here is based on the Quare P. reticulata genome analysis. 
#If running the analyses based on the X. maculatus or P. reticulata reference genomes then the coverage_fold_change_inversion.txt file should be used and "WindowStart" should be replaced by "WinsowStart_inversion" throughout
datacov <- read.csv("coverage_fold_change.txt", row.names=1,header=T)
datacovsorted <- datacov[order(datacov$Chromosome,datacov$WindowStart),]

#Calculate confidence interval for moving average
datacovminusSexChr <- datacovsorted[datacovsorted$Chromosome!=sexchromosome,]
MFpermutescov <- replicate(1000,mean(sample(datacovminusSexChr$MFLogaverage,windowsizeCI,replace = FALSE)))
MFCI25cov <- quantile(MFpermutescov, c(.025, .5, .975))[[1]]
MFCI975cov <- quantile(MFpermutescov, c(.025, .5, .975))[[3]]

#plot sex chromosome coverage data
datacovSexChr <- datacovsorted[datacovsorted$Chromosome==sexchromosome,]
smoothlinefccov = movingaverage(datacovSexChr$MFLogaverage,windowsize)
paneldatacov <- as.data.frame(smoothlinefccov)
paneldatacov$startmb <- datacovSexChr$WindowStart/1000000
paneldatacov <- na.omit(paneldatacov)

plotfccov <- ggplot(datacovSexChr, aes(x= WindowStart/1000000, y= MFLogaverage)) +
	 geom_rect(xmax=60,xmin=-10,ymax=MFCI975cov,ymin=MFCI25cov,fill="grey", alpha=0.08)+
	 geom_rect(xmax=20,xmin=26,ymax=-0.4,ymin=0.4,fill="#d6cadd", alpha=0.005)+
	 geom_point(colour="#919191",fill="#919191", alpha=0.5, cex=0.4) +
     geom_line(data = paneldatacov, aes(x= startmb, y=smoothlinefccov), size=0.6) +
     coord_cartesian(ylim=c(-0.3,0.3)) +
     theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
     theme(
           text=element_text(size=12),
           plot.margin = unit(c(1,1,1,1),"lines"),
           axis.line.y = element_line(color="black", size = 0.3),
           axis.line.x = element_line(color="black", size = 0.3),
           axis.text.x = element_text(size=10),
           axis.text.y = element_text(size=10)
           ) +
       scale_y_continuous(labels=scales::number_format(accuracy=0.1), limits=c(-0.3,0.3),breaks=seq(-0.3,0.3,by=0.1)) + 
       scale_x_continuous(breaks=seq(0,35,5)) +
       xlab('Start position (Mb)') +
       ylab(expression('M:F log'[2]*' coverage'))
plotfccov
