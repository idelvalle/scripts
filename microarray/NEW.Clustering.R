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


setwd("~/Desktop/Atlas.Array/analysis/Clustering/NEW")


### CLUSTERING TESTIS ####

data_matrix <- read.table(file = "COMPLETE.txt", sep ="\t", header = TRUE, row.names =1)
head(data_matrix)

data_matrix <- data_matrix[,c(17:36)]
colnames(data_matrix)

CS18 <- apply(data_matrix[,1:3], 1, mean)
CS19 <- data_matrix[,4]
CS21 <- apply(data_matrix[,5:8], 1, mean)
CS22 <- data_matrix[,9]
CS23 <- apply (data_matrix[,10:12],1,mean)
E8.5 <- apply (data_matrix[,13:16],1,mean)
E9 <- apply (data_matrix[,17:18],1,mean)
E10 <- apply (data_matrix[,19:20],1,mean)

data_matrix.mean <- rbind(CS18,CS19,CS21,CS22,CS23,E8.5,E9,E10) # Combines data by rows (colum of datasets must be the same)
data_matrix.mean <- t(data_matrix.mean)

k.max <- 15 # Maximal number of cluster

wss <- sapply(1:k.max, 
              function(k){kmeans(data_matrix.mean, k, nstart=10 )$tot.withinss})
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
abline(v = 3, lty =2)

results.testis <- ExpressionSet(assayData=data_matrix.mean)
head(exprs(results.testis))



#MFUZZ

results <- standardise(results.testis)
m1 <- mestimate(results)
m1
tmp <- Dmin(results,m=1.428,crange=seq(2,20,2),repeats=3,visu=TRUE)

cl <- mfuzz(results,c=4,m=1.428)

attributes(cl)
clusters <- cl$cluster
class(clusters)
write.table(clusters, "4-clusters.testis.txt", sep="\t")

pdf("mfuzz.4-clusters.testis.plot2.pdf")
mfuzz.plot2(results, cl = cl, mfrow = c(3,3), time.labels=c("CS18","CS19","CS21","CS22","CS23","E8.5","E9","E10"),x11 = FALSE)
dev.off()


# MFUZZ K-MEANS

km <- kmeans2(results, k = 4)
km
attributes(km)
cluster.km <- km$cluster

write.table(cluster.km, "km-4-clusters.testis.txt", sep="\t")

kmeans2.plot(results,kl=km,mfrow=c(3,3))
dev.print(pdf, "km-4-clusters.testis.pdf")
dev.off()

### CLUSTERING ADRENAL ####
rm(list = ls())
data_matrix <- read.table(file = "COMPLETE.txt", sep ="\t", header = TRUE, row.names =1)
head(data_matrix)

data_matrix <- data_matrix[,c(37:53)]
colnames(data_matrix)

CS17 <- data_matrix[,1]
CS18 <- apply(data_matrix[,2:3], 1, mean)
CS21 <- apply(data_matrix[,4:5], 1, mean)
CS23 <- apply (data_matrix[,6:8],1,mean)
E8.5 <- apply (data_matrix[,9:12],1,mean)
E9 <- apply (data_matrix[,13:16],1,mean)
E10 <- data_matrix[,17]

data_matrix.mean <- rbind(CS17,CS18,CS21,CS23,E8.5,E9,E10) # Combines data by rows (colum of datasets must be the same)
data_matrix.mean <- t(data_matrix.mean)

k.max <- 15 # Maximal number of cluster

wss <- sapply(1:k.max, 
              function(k){kmeans(data_matrix.mean, k, nstart=10 )$tot.withinss})
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
abline(v = 3, lty =2)

results.adrenal <- ExpressionSet(assayData=data_matrix.mean)
head(exprs(results.adrenal))



#MFUZZ

results <- standardise(results.adrenal)
m1 <- mestimate(results)
m1
tmp <- Dmin(results,m=1.539,crange=seq(2,20,2),repeats=3,visu=TRUE)

cl <- mfuzz(results,c=3,m=1.539)

attributes(cl)
clusters <- cl$cluster
class(clusters)
write.table(clusters, "3-clusters.adrenal.txt", sep="\t")

pdf("mfuzz.3-clusters.adrenal.plot2.pdf")
mfuzz.plot2(results, cl = cl, mfrow = c(3,3), time.labels=c("CS17","CS18","CS21","CS23","E8.5","E9","E10"),x11 = FALSE)
dev.off()


# MFUZZ K-MEANS

km <- kmeans2(results, k = 4)
km
attributes(km)
cluster.km <- km$cluster

write.table(cluster.km, "km-4-clusters.adrenal.txt", sep="\t")

kmeans2.plot(results,kl=km,mfrow=c(3,3))
dev.print(pdf, "km-4-clusters.adrenal.pdf")
dev.off()
