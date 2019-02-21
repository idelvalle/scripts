## CLUSTERING OVARY ##

#All Data up to 10 weeks

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


setwd("~/Desktop/Atlas.Array/analysis/Clustering/")

ovary <- read.table(file="ovary.txt", sep="\t")
colnames(ovary) <- "genes"
nrow(ovary) # 247 genes

data_matrix <- read.table(file = "COMPLETE.txt", sep ="\t", header = TRUE, row.names =1)
head(data_matrix)
data_matrix <- data_matrix[,c(7:16)]
colnames(data_matrix)
ov <- data_matrix[rownames(data_matrix) %in% ovary$genes,]

pdf(file = "Heatmap.Ovary.Clustering.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
pheatmap(ov, scale="row",color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), cluster_cols = FALSE, cutree_rows = 12, show_rownames = TRUE, border_color="NA",fontsize_row = 1)
dev.off()

six <- ov[,1]
seven <- ov[,2]
seven.half <- apply(ov[,3:5], 1, mean)
eight <- apply (ov[,6:7],1,mean)
nine <- apply (ov[,8:10],1,mean)
ov.mean <- rbind(six,seven,seven.half,eight,nine) # Combines data by rows (colum of datasets must be the same)
ov.mean <- t(ov.mean)

k.max <- 15 # Maximal number of cluster

wss <- sapply(1:k.max, 
              function(k){kmeans(ov.mean, k, nstart=10 )$tot.withinss})
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
abline(v = 3, lty =2)

results.ovary <- ExpressionSet(assayData=ov.mean)
head(exprs(results.ovary))



#MFUZZ

results <- standardise(results.ovary)
m1 <- mestimate(results)
m1
tmp <- Dmin(results,m=2.28,crange=seq(2,20,2),repeats=3,visu=TRUE)

cl <- mfuzz(results,c=3,m=2.28)

attributes(cl)
clusters <- cl$cluster
class(clusters)
write.table(clusters, "3-clusters.ovary.txt", sep="\t")

pdf("mfuzz.3-clusters.ovary.plot2.pdf")
mfuzz.plot2(results, cl = cl, mfrow = c(3,3), time.labels=c(6,7,7.5,8,9),x11 = FALSE)
dev.off()

ovary.cluster <- read.csv(file="3-clusters.ovary.csv", header=TRUE)
colnames(ovary.cluster)
nrow(ovary.cluster) # 247 genes
cluster.1 <- ovary.cluster[ovary.cluster$cluster==1,]
cluster.2 <- ovary.cluster[ovary.cluster$cluster==2,]
cluster.3 <- ovary.cluster[ovary.cluster$cluster==3,]

data_matrix <- read.table(file = "COMPLETE.txt", sep ="\t", header = TRUE, row.names =1)
head(data_matrix)

clusters.ovary <- data_matrix[rownames(data_matrix) %in% ovary.cluster$symbol,]
rownames(ovary.cluster) <- ovary.cluster$symbol
cluster.heatmap <- clusters.ovary[rownames(ovary.cluster),]
head(cluster.heatmap)
pdf(file = "Heatmap.ovary.Clusters.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
pheatmap(cluster.heatmap, scale="row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, fontsize_row = 2,color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255),border_color="NA")
dev.off()

data_matrix <- data_matrix[,c(7:16)]
colnames(data_matrix)

clusters.ovary.1 <- data_matrix[rownames(data_matrix) %in% cluster.1$symbol,]
pdf(file = "Heatmap.Ovary.Cluster-1.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
pheatmap(clusters.ovary.1, scale="row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, fontsize_row = 4,color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255),border_color="NA")
dev.off()

# MFUZZ K-MEANS

km <- kmeans2(results, k = 3)
km
attributes(km)
cluster.km <- km$cluster

write.table(cluster.km, "km-3-clusters.ovary.txt", sep="\t")

kmeans2.plot(results,kl=km,mfrow=c(3,3))
dev.print(pdf, "km.adrenal.plot2.pdf")
dev.off()

### CLUSTERING TESTIS ####

testis <- read.table(file="testis.txt", sep="\t")
colnames(testis) <- "genes"
nrow(testis) # 157 genes

data_matrix <- read.table(file = "COMPLETE.txt", sep ="\t", header = TRUE, row.names =1)
head(data_matrix)

data_matrix <- data_matrix[,c(17:36)]
colnames(data_matrix)

tst <- data_matrix[rownames(data_matrix) %in% testis$genes,]

pdf(file = "Heatmap.Testis.Clustering.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
pheatmap(tst, scale="row", cluster_cols = FALSE, cutree_rows = 6, show_rownames = TRUE, fontsize_row = 1,color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255),border_color="NA")
dev.off()

six <- apply(tst[,1:4], 1, mean)
seven <- apply(tst[,5:8], 1, mean)
seven.half <- apply(tst[,9:12], 1, mean)
eight <- apply (tst[,13:16],1,mean)
nine <- apply (tst[,17:18],1,mean)
ten <- apply (tst[,19:20],1,mean)
tst.mean <- rbind(six,seven,seven.half,eight,nine,ten) # Combines data by rows (colum of datasets must be the same)
tst.mean <- t(tst.mean)

k.max <- 15 # Maximal number of cluster

wss <- sapply(1:k.max, 
              function(k){kmeans(tst.mean, k, nstart=10 )$tot.withinss})
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
abline(v = 3, lty =2)

results.testis <- ExpressionSet(assayData=tst.mean)
head(exprs(results.testis))



#MFUZZ

results <- standardise(results.testis)
m1 <- mestimate(results)
m1
tmp <- Dmin(results,m=2.045,crange=seq(2,20,2),repeats=3,visu=TRUE)

cl <- mfuzz(results,c=3,m=2.045)

attributes(cl)
clusters <- cl$cluster
class(clusters)
write.table(clusters, "3-clusters.testis.txt", sep="\t")

pdf("mfuzz.3-clusters.testis.plot2.pdf")
mfuzz.plot2(results, cl = cl, mfrow = c(3,3), time.labels=c(6,7,7.5,8,9,10),x11 = FALSE)
dev.off()

testis.cluster <- read.csv(file="3-clusters.testis.csv", header=TRUE)
colnames(testis.cluster)
nrow(testis.cluster) # 157 genes
cluster.1 <- testis.cluster[testis.cluster$cluster==1,]
cluster.2 <- testis.cluster[testis.cluster$cluster==2,]
cluster.3 <- testis.cluster[testis.cluster$cluster==3,]

data_matrix <- read.table(file = "COMPLETE.txt", sep ="\t", header = TRUE, row.names =1)
head(data_matrix)

clusters.testis <- data_matrix[rownames(data_matrix) %in% testis.cluster$symbol,]
rownames(testis.cluster) <- testis.cluster$symbol
cluster.heatmap <- clusters.testis[rownames(testis.cluster),]
head(cluster.heatmap)
pdf(file = "Heatmap.Testis.Clusters.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
pheatmap(cluster.heatmap, scale="row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, fontsize_row = 4,color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255),border_color="NA")
dev.off()

clusters.testis.1 <- data_matrix[rownames(data_matrix) %in% cluster.1$symbol,]
pdf(file = "Heatmap.Testis.Cluster-1.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
pheatmap(clusters.testis.1, scale="row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, fontsize_row = 4,color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255),border_color="NA")
dev.off()

# MFUZZ K-MEANS

km <- kmeans2(results, k = 3)
km
attributes(km)
cluster.km <- km$cluster

write.table(cluster.km, "km-3-clusters.testis.txt", sep="\t")

kmeans2.plot(results,kl=km,mfrow=c(3,3))
dev.print(pdf, "km-3-clusters.testis.pdf")
dev.off()


####### CLUSTERING OVARY-TESTIS ############

ovary <- read.table(file="ovary-testis.NODUP.txt", sep="\t")
nrow(ovary) # 273 genes

data_matrix <- read.table(file = "COMPLETE.txt", sep ="\t", header = TRUE, row.names =1)
head(data_matrix)
data_matrix <- data_matrix[,c(7:16)]
colnames(data_matrix)
ov <- data_matrix[rownames(data_matrix) %in% ovary$SYMBOL,]

pdf(file = "Heatmap.Ovary-Testis.Clustering.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
pheatmap(ov, scale="row",color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), cluster_cols = FALSE, show_rownames = TRUE, border_color="NA",fontsize_row = 1)
dev.off()

six <- ov[,1]
seven <- ov[,2]
seven.half <- apply(ov[,3:5], 1, mean)
eight <- apply (ov[,6:7],1,mean)
nine <- apply (ov[,8:10],1,mean)
ov.mean <- rbind(six,seven,seven.half,eight,nine) # Combines data by rows (colum of datasets must be the same)
ov.mean <- t(ov.mean)

k.max <- 15 # Maximal number of cluster

wss <- sapply(1:k.max, 
              function(k){kmeans(ov.mean, k, nstart=10 )$tot.withinss})
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
abline(v = 3, lty =2)

results.ovary <- ExpressionSet(assayData=ov.mean)
head(exprs(results.ovary))



#MFUZZ

results <- standardise(results.ovary)
m1 <- mestimate(results)
m1
tmp <- Dmin(results,m=2.26,crange=seq(2,20,2),repeats=3,visu=TRUE)

cl <- mfuzz(results,c=3,m=2.28)

attributes(cl)
clusters <- cl$cluster
class(clusters)
write.table(clusters, "3-clusters.ovary-testis.txt", sep="\t")

pdf("mfuzz.3-clusters.ovary-testis.plot2.pdf")
mfuzz.plot2(results, cl = cl, mfrow = c(3,3), time.labels=c(6,7,7.5,8,9),x11 = FALSE)
dev.off()

ovary.cluster <- read.csv(file="3-clusters.ovary-testis.csv", header=TRUE)
colnames(ovary.cluster)
nrow(ovary.cluster) # 273 genes
cluster.1 <- ovary.cluster[ovary.cluster$Cluster==1,]
cluster.2 <- ovary.cluster[ovary.cluster$Cluster==2,]
cluster.3 <- ovary.cluster[ovary.cluster$Cluster==3,]

data_matrix <- read.table(file = "COMPLETE.txt", sep ="\t", header = TRUE, row.names =1)
head(data_matrix)

clusters.ovary <- data_matrix[rownames(data_matrix) %in% ovary.cluster$Gene,]
rownames(ovary.cluster) <- ovary.cluster$Gene
cluster.heatmap <- clusters.ovary[rownames(ovary.cluster),]
head(cluster.heatmap)
pdf(file = "Heatmap.ovary-testis.Clusters.All.Samples.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
pheatmap(cluster.heatmap, scale="row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, fontsize_row = 2,color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255),border_color="NA")
dev.off()

data_matrix <- data_matrix[,c(7:16)]
colnames(data_matrix)

clusters.ovary.1 <- data_matrix[rownames(data_matrix) %in% cluster.1$Gene,]
pdf(file = "Heatmap.Ovary-Testis.Cluster-1.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
pheatmap(clusters.ovary.1, scale="row", cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, fontsize_row = 4,color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255),border_color="NA")
dev.off()

clusters.ovary.2 <- data_matrix[rownames(data_matrix) %in% cluster.2$Gene,]
pdf(file = "Heatmap.Ovary-Testis.Cluster-2.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
pheatmap(clusters.ovary.2, scale="row", cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, fontsize_row = 4,color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255),border_color="NA")
dev.off()

clusters.ovary.3 <- data_matrix[rownames(data_matrix) %in% cluster.3$Gene,]
pdf(file = "Heatmap.Ovary-Testis.Cluster-3.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
pheatmap(clusters.ovary.3, scale="row", cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, fontsize_row = 4,color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255),border_color="NA")
dev.off()

# MFUZZ K-MEANS

km <- kmeans2(results, k = 3)
km
attributes(km)
cluster.km <- km$cluster

write.table(cluster.km, "km-3-clusters.ovary-testis.txt", sep="\t")

kmeans2.plot(results,kl=km,mfrow=c(3,3))
dev.print(pdf, "km.ovary-testis.plot2.pdf")
dev.off()

