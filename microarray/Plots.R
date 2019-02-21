
############ PLOT SRY AND SOX9 ###############

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
library(ggfortify)
library(plotly)

setwd("~/Desktop/Atlas.Array/analysis/all_up_to_10_weeks/paper/Time-course-graphs/")

##### The file COMPLETE-graphs.csv contains all batch-corrected 23307 genes after NA and dup removal ##########

genes <- read.csv(file="COMPLETE-graphs.csv", header=TRUE, row.names=1, sep =",")
genes <- t(genes)
genes[,1:3]

control <- genes[1:6,]
ovary <- genes[7:16,]
testis <- genes[17:36,]
adrenal <- genes[37:53,]

sry <- ggplot(testis, aes(x=STAGE, y=SRY)) + geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "SRY.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
sry + expand_limits(y=c(0,8.5), x=c(45.75,74))
sry + expand_limits(y=c(4,6), x=c(45.75,74))
sry + expand_limits(y=c(0,6), x=c(45.75,74))
dev.off()

lhcgr <- ggplot(testis, aes(x=STAGE, y=LHCGR))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "LHCGR.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
lhcgr + expand_limits(y=c(0,9), x=c(45.75,74))
dev.off()

star <- ggplot(testis, aes(x=STAGE, y=STAR))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "STAR.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
star + expand_limits(y=c(0,14), x=c(45.75,74))
dev.off()

cyp11a1 <- ggplot(testis, aes(x=STAGE, y=CYP11A1))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "CYP11A1.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
cyp11a1 + expand_limits(y=c(0,14), x=c(45.75,74))
dev.off()

hsd3b2 <- ggplot(testis, aes(x=STAGE, y=HSD3B2))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "HSD3B2.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
hsd3b2 + expand_limits(y=c(0,10), x=c(45.75,74))
dev.off()

cyp17a1 <- ggplot(testis, aes(x=STAGE, y=CYP17A1))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "CYP17A1.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
cyp17a1 + expand_limits(y=c(0,14), x=c(45.75,74))
dev.off()

hsd17b3 <- ggplot(testis, aes(x=STAGE, y=HSD17B3))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "HSD17B3.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
hsd17b3 + expand_limits(y=c(0,9), x=c(45.75,74))
dev.off()

sox9 <- ggplot(testis, aes(x=STAGE, y=SOX9))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "SOX9.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
sox9 + expand_limits(y=c(7,8.6), x=c(45.75,74))
dev.off()

#sox9 + expand_limits(y=c(5,8.6)) + ylab("SOX9 values") + xlab("Days")

cited1 <- ggplot(testis, aes(x=STAGE, y=CITED1))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "CITED1.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
cited1 + expand_limits(y=c(5,11), x=c(45.75,74))
dev.off()


aspn <- ggplot(testis, aes(x=STAGE, y=ASPN))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "ASPN.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
aspn + expand_limits(y=c(0,8), x=c(45.75,74))
dev.off()

slc52a3 <- ggplot(testis, aes(x=STAGE, y=SLC52A3))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "SLC52A3.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
slc52a3 + expand_limits(y=c(0,8), x=c(45.75,74))
dev.off()

rmdn2 <- ggplot(testis, aes(x=STAGE, y=RMDN2))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "RMDN2.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
rmdn2 + expand_limits(y=c(0,8), x=c(45.75,74))
dev.off()

map3k15 <- ggplot(testis, aes(x=STAGE, y=MAP3K15))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "MAP3K15.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
map3k15 + expand_limits(y=c(0,9), x=c(45.75,74))
dev.off()

foxo4 <- ggplot(testis, aes(x=STAGE, y=FOXO4))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "FOXO4.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
foxo4 + expand_limits(y=c(0,12), x=c(45.75,74)) 
dev.off()

gramd1b <- ggplot(testis, aes(x=STAGE, y=GRAMD1B))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "GRAMD1B.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
gramd1b + expand_limits(y=c(0,12), x=c(45.75,74))
dev.off()

znf280b <- ggplot(testis, aes(x=STAGE, y=ZNF280B))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "ZNF280B.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
znf280b + expand_limits(y=c(6.6,8.2), x=c(45.75,74))
dev.off()

txlngy <- ggplot(testis, aes(x=STAGE, y=TXLNGY))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "TXLNGY.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
txlngy + expand_limits(y=c(0,8.5), x=c(45.75,74))
dev.off()

ren <- ggplot(testis, aes(x=STAGE, y=REN))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "REN.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
ren + expand_limits(y=c(0,8.5), x=c(45.75,74))
dev.off()

ren <- ggplot(ovary, aes(x=STAGE, y=REN))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
ren + expand_limits(y=c(0,12), x=c(45.75,67))

rspo <- ggplot(ovary, aes(x=STAGE, y=RSPO1))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
rspo + expand_limits(y=c(2.5,10), x=c(45.75,67))

wnt4 <- ggplot(ovary, aes(x=STAGE, y=WNT4))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
wnt4 + expand_limits(y=c(5,9), x=c(45.75,67))

ctnnb1 <- ggplot(ovary, aes(x=STAGE, y=CTNNB1))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
ctnnb1 + expand_limits(y=c(9,13), x=c(45.75,67))

prps2 <- ggplot(testis, aes(x=STAGE, y=PRPS2))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
prps2 + expand_limits(y=c(7,9.5), x=c(45.75,74))

ank <- ggplot(testis, aes(x=STAGE, y=ANKRD18A))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
ank + expand_limits(y=c(5,8), x=c(45.75,74))

dhx37 <- ggplot(testis, aes(x=STAGE, y=DHX37))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()
pdf(file = "DHX37.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
dhx37 + expand_limits(y=c(6,8), x=c(45.75,74)) + labs(x= "dpc", y="Log2 normalised values", title="DHX37 Testis")
dev.off()

ggplot(control, aes(x=STAGE, y=DHX37)) + geom_violin()

aspn

session <- sessionInfo()

save.image("Graphs.RData")

ggplot(ovary,aes(x=STAGE, y=DCAF4L1))+geom_point()+geom_smooth(method="loess", color="red",level=0.95) + theme_bw()+expand_limits(y=c(0,7))
