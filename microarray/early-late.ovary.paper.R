
######## #Early vs Late Ovary (6-7 weeks vs 8-9 weeks) 2 vs 5 samples #####


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

setwd("~/Desktop/Atlas.Array/analysis/early.ovary-late.ovary/paper/")

data_matrix <- read.table(file = "expression.combat.all.txt", sep ="\t", header = TRUE, row.names =1)

data_matrix <- data_matrix[,c(23,27,28:32)]

colnames(data_matrix)
######################## DIFFERENTIAL EXPRESSION USING LIMMA ###############################

#limma

eset <-new("ExpressionSet", exprs= as.matrix(data_matrix))
annotation(eset) <- "pd.hugene.1.0.st.v1"
featureData(eset) <- getNetAffx(eset, type="transcript")

targets <- c("early","early",rep("late",5))
targets
f <- factor(targets)
f
design <- model.matrix(~0+f) # Create Design Matrix
colnames(design) <- levels(f)
design
data.fit <- lmFit(eset,design) # Fit Linear Model
names(data.fit)
colnames(design)

#coef=1: late-early


contrast.matrix = makeContrasts(late-early,levels=design)
contrast.matrix
data.fit.con = contrasts.fit(data.fit,contrast.matrix) # Compute estimated coefficients and standard errors
data.fit.eb = eBayes(data.fit.con) # Bayes Statistics for Differential Expression
data.fit.eb


######### LATE TESTIS VS EARLY TESTIS ############

ovary <-topTable(data.fit.eb, coef=1, number=10000, adjust="fdr", lfc=1) # Extract genes from linear model fit
nrow(ovary) # 374

ovary.final <- ovary[,c(1,8, 19:24)]
results.ovary <- separate(ovary.final, geneassignment, into = c("NM", "SYMBOL"), sep = " // ", extra = "drop") #Pay attention to SEPARATION FORMAT!!!
write.table(results.ovary, "early-late.ovary.complete.txt", row.names = TRUE, sep="\t")


################## HEATMAP TOP50 Early Testis VS CONTROL ################

complete <- read.table(file ="COMPLETE.txt", row.names = 1, header = TRUE, sep="\t")

complete.selected <- complete[,c(7,8,12:16)]

complete.selected[1:3,]

top.ovary <- results.ovary[!duplicated(results.ovary[,3]),] # Removing duplicated gene names
nrow(top.ovary) # 317 genes

top.ovary <- top.ovary[complete.cases(top.ovary[,3]),] # Removing NA values
nrow(top.ovary) # 316 genes

write.table(top.ovary, "early-late.ovary.txt", row.names = TRUE, sep="\t")

top.ovary <- top.ovary[,-c(1:2)]
rownames(top.ovary) <- NULL
rownames(top.ovary) <- top.ovary[,1]
top.ovary <- top.ovary[,-1]

top.ovary <- top.ovary[order(-top.ovary$logFC),]
top.ovary[1:3,]

top50.ovary <- top.ovary[1:50,]
top50.ovary

top50.ovary_selected <- complete.selected[rownames(complete.selected) %in% rownames(top50.ovary),]

top50.ovary_selected[1:3,]


ovr.top50 <- top50.ovary_selected[rownames(top50.ovary),] #### IMPORTANT TO MAINTAIN ORDER AFTER MATCHING

colnames(ovr.top50)

pdf(file = "Heatmap.Early-Late.Ovary.Top50.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)

pheatmap(ovr.top50, color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), scale = "row", cluster_rows = FALSE,
         cluster_cols = FALSE,fontsize = 8,border_color="NA")

dev.off()

write.table(ovr.top50, file="Heatmap.Early-Late.Ovary.Top50.txt", row.names=TRUE, sep="\t")


######################### DENDROGRAM AND PCA ANALYSIS ##############################################

#Dendrogram

distance <- dist(t(data_matrix), method = "maximum")
clusters <- hclust(distance, method = "ward.D2")
clusters$order

colnames(data_matrix)

pdf(file = "dendrogram.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
plot(clusters, labels=c("Ovary CS18", "Ovary CS21","Ovary 8.5w.1","Ovary 8.5w.2","Ovary 9w.1","Ovary 9w.2","Ovary 9w.3"),hang = -1, cex = 0.8, main = "ward.D2")
dev.off()

#PCA

colnames(data_matrix)

color <- c(rep("#ef6eb0",2),rep('#b01463',5))
pca <- prcomp(t(data_matrix),scale.=TRUE, center=TRUE)

summary(pca) # First 2 component explain more than 50% of the variance
plot(pca, type="l")
loadings <- pca$rotation 
rownames(pca$x)

pdf(file = "PCA.expression.early-late.ovary.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
plot(pca$x[,1:2],col=color, pch = 16, cex = 1.5, ylab = "PC2 (20%)", xlab = "PC1 (30%)")
legend("topleft", legend = c("Early Ovary", "Late Ovary"), col = c('#ef6eb0','#b01463'), pch = 19, cex=.9)
dev.off()


########################### CORRELATION PLOTS ####################################

early <- complete.selected[, c(1:2)]

early$early.mean <- apply(early, 1 , mean)
early$early.log2 <- log2(early$early.mean)

late <- complete.selected[,c(3:7)]

late$late.mean <- apply(late, 1 , mean)
late$late.log2 <- log2(late$late.mean)

samples.log2 <- data.frame(early$early.log2,late$late.log2, row.names = rownames(early))
samples.log2[1:3,]
nrow(samples.log2) # 23307 genes
colnames(samples.log2) <- c("early", "late")

cor(samples.log2) # Pair-wise sample correlation higher than 98%

samples.log2["SRY",] # Control=1.9 & Testis=2


pdf(file = "corr.early.ovary-late.ovary.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)

qplot(early, late, data = samples.log2, 
      label = rownames(samples.log2), xlab = "Early Ovary", ylab = "Late Ovary", main = "Log2 Average Values")

dev.off()

p1 <- qplot(early, late, data = samples.log2, 
            label = rownames(samples.log2), xlab = "Early Ovary", ylab = "Late Ovary", main = "Log2 Average Values")

#ggplotly(p1)

session <- sessionInfo()

save.image(file = "Early-Late.Ovary.RData")
