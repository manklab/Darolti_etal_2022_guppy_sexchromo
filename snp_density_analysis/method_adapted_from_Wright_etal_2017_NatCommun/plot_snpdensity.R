rm(list=ls())
ls() 

library(ggplot2)
library(reshape)
library(gridExtra)
library(ggthemes)
library(extrafont)
loadfonts()
library(RcppRoll)
require(cowplot)

#snp data
Auto <- read.csv("snpdensity_autosomes.txt", header=T)
Sex <- read.csv("snpdensity_sexchromosomes.txt", header=T)

#sort data
AutoSorted <-Auto[order(Auto$LG,Auto$Start),]
SexSorted <-Sex[order(Sex$LG,Sex$Start),]

#significance tests
median(AutoSorted$MFLogaverage)
median(SexSorted$MFLogaverage)
wilcox.test(AutoSorted$MFLogaverage, SexSorted$MFLogaverage)

#plot
theme_set(theme_gray())

snpdensity_plot <- ggplot(AutoSorted, aes(x=MFLogaverage)) + 
     geom_density(colour="#000000",fill="#5F808D",size=0.7,alpha=0.7) +
     geom_density(data=SexSorted, aes(x=MFLogaverage),colour="#000000", fill="#E1B91A", alpha=0.7, size=0.7)+
     geom_vline(aes(xintercept=median(AutoSorted$MFLogaverage)),color="#5F808D",linetype="dashed",size=0.8) +
     geom_vline(aes(xintercept=median(SexSorted$MFLogaverage)),color="#E1B91A",linetype="dashed",size=0.8) +
     coord_cartesian(ylim = c(0,2),xlim=c(-1,1)) +
     theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
     theme(
     		axis.title.y=element_text(margin=margin(0,5,0,0)),
           text=element_text(size=12),
           plot.margin = unit(c(2.5,6,1,3),"lines"),
          axis.line.y = element_line(color="black", size = 0.5),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.text=element_text(size=10)
           ) +
       scale_y_continuous(expand = c(0,0)) + 
       xlab(expression('M:F log'[2]*' normalised SNP density')) +
       ylab("Density") +
       geom_text(x=0.5, y=1.5, label="***", size=5)
snpdensity_plot
