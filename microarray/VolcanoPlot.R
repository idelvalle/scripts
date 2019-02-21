
############# VOLCANO PLOT ADRENAL vs CONTROL ##############

rm(list = ls())

library(oligo)
library("ggplot2")
library("ggdendro")
library("RColorBrewer")
library(tidyr)
library(dplyr)
library(gplots)
library( "RColorBrewer" )
library("limma")
library(pheatmap)
library(Mfuzz)
library(extrafont)
library(sva)
library(pamr)
library(Biobase)

setwd("~/Desktop/Atlas.Array/analysis/all_up_to_10_weeks/paper")

##### The file adrenal-control-volcano.txt contains all (cut-off LogFC=0) 23307 genes after NA and dup removal ##########

volcano <- read.table(file="adrenal-control-volcano.txt", header=TRUE, sep ="\t", row.names =1)
volcano[1:3,]
nrow(volcano) # 23307 genes

pdf(file = "VolcanoPlot.Adrenal-Control.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
with(volcano, plot(logFC, -log10(P.Value), pch=20, xlab="Log2FC",  xlim=c(-4,9)))
with(subset(volcano, adj.P.Val<.05 & abs(logFC)>2), points(logFC, -log10(P.Value), pch=20, col="#D95F02"))
identify(volcano$logFC, -log10(volcano$P.Value), labels=volcano$SYMBOL, cex = 0.8)
dev.off()

#################### VOLCANO PLOT GGPLOT2 ###############

require(ggplot2)

##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off

volcano$threshold = as.factor(abs(volcano$logFC) > 2 & volcano$P.Value < .000005)
number<-subset(volcano,abs(volcano$logFC) > 2 & volcano$P.Value < .000005) 
nrow(number) # 285

##Construct the plot object
pdf(file = "VolcanoPlot.Adrenal-Control-ggplot.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
svg(file = "VolcanoPlot.Adrenal-Control-ggplot.svg", width = 7, height = 6)
ggplot(data=volcano, aes(x=logFC, y=-log10(P.Value), colour=threshold))+scale_color_manual(values=c("#000000", "#D95F02"))+ geom_point(alpha=0.4, size=1.75) + xlim(c(-4, 9)) + ylim(c(0, 50)) + xlab("log2 FC") + ylab("-log10 p-value")+
theme_bw() +theme(legend.position="none")
dev.off()

session <- sessionInfo()

save.image(file="VolcanoPlot.RData")
