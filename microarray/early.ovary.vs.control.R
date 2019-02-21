
######## EARLY OVARY 6-7 WEEKS VS CONTROLS #######


rm(list=ls())

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
library(plotly)
library(extrafont)
library(Biobase)

setwd("~/Desktop/Atlas.Array/analysis/early.ovary.vs.control/paper")

data_matrix <- read.table(file = "expression.combat.all.txt", sep ="\t", header = TRUE, row.names =1)

data_matrix <- data_matrix[,c(18:22,23,27,33)]

colnames(data_matrix)
######################## DIFFERENTIAL EXPRESSION USING LIMMA ###############################

#limma

eset <-new("ExpressionSet", exprs= as.matrix(data_matrix))
annotation(eset) <- "pd.hugene.1.0.st.v1"
featureData(eset) <- getNetAffx(eset, type="transcript")
eset

targets <- c(rep("control",5), rep("ovary",2),"control")
targets
f <- factor(targets)
f
design <- model.matrix(~0+f) # Create Design Matrix
colnames(design) <- levels(f)
design
data.fit <- lmFit(eset,design) # Fit Linear Model
names(data.fit)
colnames(design)

#coef=1: ovary-control


contrast.matrix = makeContrasts(ovary-control,levels=design)
contrast.matrix
data.fit.con = contrasts.fit(data.fit,contrast.matrix) # Compute estimated coefficients and standard errors
data.fit.eb = eBayes(data.fit.con) # Bayes Statistics for Differential Expression
data.fit.eb


######### OVARY VS CONTROL ############

ovary_control <-topTable(data.fit.eb, coef=1, number=10000, adjust="fdr", lfc=1) # Extract genes from linear model fit
nrow(ovary_control) # 2069

ovary_control.final <- ovary_control[,c(1,8, 19:24)]
results.ovary_control <- separate(ovary_control.final, geneassignment, into = c("NM", "SYMBOL"), sep = " // ", extra = "drop") #Pay attention to SEPARATION FORMAT!!!
write.table(results.ovary_control, "early.ovary-control.complete.txt", row.names = TRUE, sep="\t")


################## HEATMAP TOP50 Early Testis VS CONTROL ################

complete <- read.table(file ="COMPLETE.txt", row.names = 1, header = TRUE, sep="\t")

complete.selected <- complete[,c(1:8)]

complete.selected[1:3,]

top.ovary <- results.ovary_control[!duplicated(results.ovary_control[,3]),] # Removing duplicated gene names
nrow(top.ovary) # 1549 genes

top.ovary <- top.ovary[complete.cases(top.ovary[,3]),] # Removing NA values
nrow(top.ovary) # 1548 genes


write.table(top.ovary, "early.ovary.vs.control.txt", row.names = TRUE, sep="\t")

top.ovary <- top.ovary[,-c(1:2)]
rownames(top.ovary) <- NULL
rownames(top.ovary) <- top.ovary[,1]
top.ovary <- top.ovary[,-1]

top.ovary <- top.ovary[order(-top.ovary$logFC),]
top.ovary[1:3,]

top50.ovary <- top.ovary[1:51,]
top50.ovary

top50.ovary_selected <- complete.selected[rownames(complete.selected) %in% rownames(top50.ovary),]

top50.ovary_selected[1:3,]


ovr.top50 <- top50.ovary_selected[rownames(top50.ovary),] #### IMPORTANT TO MAINTAIN ORDER AFTER MATCHING

ovr.top50 <- ovr.top50[-35,] ############### REMOVING ZNF676 #################

pdf(file = "Heatmap.Early.Ovary.vs.Control.Top50.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)

pheatmap(ovr.top50, color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), scale = "row", cluster_rows = FALSE,
         cluster_cols = FALSE,fontsize = 8,border_color="NA")

dev.off()

write.table(ovr.top50, file="Heatmap.Early.Ovary.vs.Control.Top50.txt", row.names=TRUE, sep="\t")


######################### DENDROGRAM AND PCA ANALYSIS ##############################################

#Dendrogram

distance <- dist(t(data_matrix), method = "maximum")
clusters <- hclust(distance, method = "ward.D2")
clusters$order



pdf(file = "dendrogram.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
lab <-c("Brain", "Heart", "Kidney","Liver", "Muscle", "Ov.CS18", "Ov.CS21", "Spine")
plot(clusters, labels=lab,hang = -1, cex = 0.8, main = "ward.D2")
dev.off()

#PCA

colnames(data_matrix)

color <- c(rep('#1B9E77',5),rep('#E7298A',2),'#1B9E77')
pca <- prcomp(t(data_matrix),scale.=TRUE, center=TRUE)

summary(pca) # First 2 component explain more than 51% of the variance
plot(pca, type="l")
loadings <- pca$rotation 
rownames(pca$x)

pdf(file = "PCA.expression.early.ovary.6-7weeks.vs.Controls.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
plot(pca$x[,1:2],col=color, pch = 16, cex = 1.5, ylab = "PC2 (25%)", xlab = "PC1 (34%)")
legend("topleft", legend = c("Control", " Early Ovary"), col = c('#1B9E77','#E7298A'), pch = 19, cex=.7)
dev.off()


########################### CORRELATION PLOTS ####################################

control <- complete.selected[, c(1:6)]

control$control.mean <- apply(control, 1 , mean)
control$control.log2 <- log2(control$control.mean)

ovary <- complete.selected[,c(7:8)]

ovary$ovary.mean <- apply(ovary, 1 , mean)
ovary$ovary.log2 <- log2(ovary$ovary.mean)

samples.log2 <- data.frame(control$control.log2,ovary$ovary.log2, row.names = rownames(control))
samples.log2[1:3,]
nrow(samples.log2) # 23307 genes
colnames(samples.log2) <- c("control", "ovary")

cor(samples.log2) # Pair-wise sample correlation higher than 96%

samples.log2["SRY",] # Control=2.1 & Testis=1.9


pdf(file = "corr.early.ovary.6-7weeks.vs.Controls.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)

qplot(control, ovary, data = samples.log2, 
      label = rownames(samples.log2), xlab = "Control", ylab = "Early Ovary", main = "Log2 Average Values")

dev.off()

p1 <- qplot(control, ovary, data = samples.log2, 
            label = rownames(samples.log2), xlab = "Control", ylab = "Early Ovary", main = "Log2 Average Values")

#ggplotly(p1)

session <- sessionInfo()

save.image(file = "Early.Ovary.6-7weeks.vs.Controls.RData")
