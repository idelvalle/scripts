

############# CORRELATION PLOTS #########################

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

setwd("~/Desktop/Atlas.Array/analysis/all_up_to_10_weeks/paper")

##### The file COMPLETE.txt contains all batch-corrected 23307 genes after NA and dup removal ##########

complete <- read.table(file="COMPLETE.txt", header=TRUE, sep ="\t", row.names =1)
complete[1:3,]
nrow(complete) # 23307 genes


###### CORRELATION PLOT ALL DATA #######

pdf(file = "CorrelationPlot.Samples.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
pheatmap(cor(complete), color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), scale = "none", cluster_rows = FALSE,
         cluster_cols = FALSE,fontsize = 8,border_color="NA")
dev.off()

########### CORRELATION PLOTS BETWEEN SAMPLES ################

colnames(complete)

control <- complete[,c(1:6)]
control[1:3,]
length(colnames(control))

control$control.mean <- apply(control, 1 , mean)
control$control.log2 <- log2(control$control.mean)

ovary <- complete[,c(7:16)]
ovary[1:3,]
length(colnames(ovary))

ovary$ovary.mean <- apply(ovary, 1 , mean)
ovary$ovary.log2 <- log2(ovary$ovary.mean)

testis <- complete[,c(17:36)]
testis[1:3,]
length(colnames(testis))

testis$testis.mean <- apply(testis, 1 , mean)
testis$testis.log2 <- log2(testis$testis.mean)

adrenal <- complete[,c(37:53)]
adrenal[1:3,]
length(colnames(adrenal))

adrenal$adrenal.mean <- apply(adrenal, 1 , mean)
adrenal$adrenal.log2 <- log2(adrenal$adrenal.mean)

all.samples <- cbind.data.frame(control,ovary,testis,adrenal)

samples.log2FC <- data.frame(control$control.log2, ovary$ovary.log2, testis$testis.log2, adrenal$adrenal.log2, row.names = rownames(all.samples))
samples.log2FC[1:3,]
nrow(samples.log2FC) # 23307 genes
colnames(samples.log2FC) <- c("control", "ovary", "testis", "adrenal")

corr <- cor(samples.log2FC) # Pair-wise sample correlation higher than 96%
cor.test(samples.log2FC$control, samples.log2FC$adrenal)
cor.test(samples.log2FC$ovary, samples.log2FC$testis)
cor.test(samples.log2FC$control, samples.log2FC$ovary)
cor.test(samples.log2FC$control, samples.log2FC$testis)



################ CORRELATION PLOTS ################

pdf(file = "Adrenal-Control.corr.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
ggplot(all.samples, aes(x=control.log2, y=adrenal.log2)) + geom_point() + theme_bw() + labs(x="Control", y="Adrenal", title="Log2 Average Values")
dev.off()


p1 <- qplot(control.log2, adrenal.log2, data = all.samples, label = rownames(all.samples), xlab = "Control", ylab = "Adrenal", main = "Log2 Average Values")



pdf(file = "Testis-Control.corr.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
ggplot(all.samples, aes(x=control.log2, y=testis.log2)) + geom_point() + theme_bw() + labs(x="Control", y="Testis", title="Log2 Average Values")
dev.off()

p2 <- qplot(control.log2, testis.log2, data = all.samples, label = rownames(all.samples), xlab = "Control", ylab = "Testis", main = "Log2 Average Values")

pdf(file = "Ovary-Control.corr.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
ggplot(all.samples, aes(x=control.log2, y=ovary.log2)) + geom_point() + theme_bw() + labs(x="Control", y="Ovary", title="Log2 Average Values")
dev.off()

p3 <- qplot(control.log2, ovary.log2, data = all.samples, label = rownames(all.samples), xlab = "Control", ylab = "Ovary", main = "Log2 Average Values")

pdf(file = "Adrenal-Testis.corr.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
ggplot(all.samples, aes(x=adrenal.log2, y=testis.log2)) + geom_point() + theme_bw() + labs(x="Adrenal", y="Testis", title="Log2 Average Values")
dev.off()

p4 <- qplot(adrenal.log2, testis.log2, data = all.samples, label = rownames(all.samples), xlab = "Adrenal", ylab = "Testis", main = "Log2 Average Values")


pdf(file = "Adrenal-Ovary.corr.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
ggplot(all.samples, aes(x=adrenal.log2, y=ovary.log2)) + geom_point() + theme_bw() + labs(x="Adrenal", y="Ovary", title="Log2 Average Values")
dev.off()

p5 <- qplot(adrenal.log2, ovary.log2, data = all.samples, label = rownames(all.samples), xlab = "Adrenal", ylab = "Ovary", main = "Log2 Average Values")

pdf(file = "Ovary-Testis.corr.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
ggplot(all.samples, aes(x=ovary.log2, y=testis.log2)) + geom_point() + theme_bw() + labs(x="Ovary", y="Testis", title="Log2 Average Values")
dev.off()

p6 <- qplot(ovary.log2, testis.log2, data = all.samples, label = rownames(all.samples), xlab = "Ovary", ylab = "Testis", main = "Log2 Average Values")


ggplotly(p1)
ggplotly(p2)
ggplotly(p3)
ggplotly(p4)
ggplotly(p5)
ggplotly(p6)

session <- sessionInfo()

save.image(file="Correlation.RData")

