#All Data up to 10 weeks (53 samples)

#Excluding Adrenal_6w_988.CEL

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
library(plotly)
library(extrafont)
library(Biobase)
library(sva)

setwd("~/Desktop/Atlas.Array/analysis/all_up_to_10_weeks/paper")

############################## READING DATA ####################################

setwd("~/Desktop/Atlas.Array/analysis/all_up_to_10_weeks/paper")

celpath <- "~/Desktop/Atlas.Array/data/all_up_to_10_weeks"
list <- list.files(celpath,full.names=TRUE)
pd <-read.AnnotatedDataFrame(filename="phenodata.txt", sep="\t", path="~/Desktop/Atlas.Array/data/phenoData" ,header=TRUE)
data <- read.celfiles(list, phenoData = pd)
pData(data)

data # GeneFeatureSet

################ PLOTS BEFORE NORMALIZATION ######################################

pdf(file = "boxplot.before.normalization.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
boxplot(data, names = pd[[2]], las =2, cex.axis=.8, ylab = "Intensity", main = "Raw Data Expression", which= "all")
dev.off()

pdf(file = "histogram.before.normalization.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
hist(data, las =2, cex.axis=.8 , main = "Raw Data Histogram", which= "all")
dev.off()

############# RMA NORMALIZATION AND  PLOTS ########################

data_rma <- rma(data, target = "core")

pdf(file = "boxplot.after.normalization.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
boxplot(data_rma, names = pd[[2]], las =2, cex.axis=.8, ylab = "Intensity", main = "Raw Data Expression")
dev.off()

pdf(file = "histogram.after.normalization.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
hist(data_rma, las =2, cex.axis=.8 , main = "Raw Data Histogram")
dev.off()

############# ANNOTATION ########################

featureData(data_rma) <- getNetAffx(data_rma, type="transcript")

probes <- fData(data_rma)
probes[1:3,]
probes <- probes[,c(2,8)] # selecting annotation fields

edata <- exprs(data_rma)

write.table( edata, file = "raw.data.txt", row.names=TRUE, sep="\t") # saving raw data file

expression <- data.frame(probes, edata)

expression[1:3,]

expression.separate <- separate(expression, geneassignment, into = c("NM", "SYMBOL"), sep = " // ", extra = "drop") #Pay attention to SEPARATION FORMAT!!!

expression.separate[1:3,]
expression.separate <- expression.separate[,-c(1:2)]

write.table(expression.separate, file="non.corrected.txt", row.names=TRUE,sep="\t") # saving raw data annotated file

############# BATCH CORRECTION ACCORDING TO DATES (5) AND TISSUE (4) ########################

pheno = pData(data)

batch = pheno$Batch

covariate <- pheno$Target

modcombat = model.matrix(~covariate, data=pheno) # Apply ComBat function from SVA Package

combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)

combat_edata[1:2,]

write.table(combat_edata,file="expression.combat.all.txt",row.names=TRUE,sep="\t")

expression.combat <- data.frame(probes, combat_edata) 

expression.separate.combat <- separate(expression.combat, geneassignment, into = c("NM", "SYMBOL"), sep = " // ", extra = "drop") #Pay attention to SEPARATION FORMAT!!!

expression.separate.combat <- expression.separate.combat[,-c(1:2)]

write.table(expression.separate.combat, file="covariate.txt", row.names=TRUE,sep="\t") #saving batch-corrected values

######################### DENDROGRAM AND PCA ANALYSIS ##############################################

#Dendrogram

data_matrix <- combat_edata

distance <- dist(t(data_matrix), method = "maximum")
clusters <- hclust(distance, method = "ward.D2")
clusters$order



pdf(file = "dendrogram.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
svg(file = "dendrogram.svg", width = 7, height = 6, family = "sans")
plot(clusters, labels = pd[[2]], hang = -1, cex = 0.6, main = "ward.D2")
dev.off()

#PCA

colnames(data_matrix)

brewer.pal(n = 4, name = "Dark2") # Hexadecimal RGB Color specification

color <- c(rep('#D95F02',17),rep('#1B9E77',5),rep('#E7298A',10),'#1B9E77', rep('#7570B3',20))
pca <- prcomp(t(data_matrix),scale.=TRUE, center=TRUE)

summary(pca) # First 2 component explain more than 37% of the variance
plot(pca, type="l")
loadings <- pca$rotation 
rownames(pca$x)

pdf(file = "PCA.pdf", width = 7, height = 6, family = "sans", useDingbats = FALSE)
svg(file = "PCA.svg", width = 7, height = 6, family = "sans")
plot(pca$x[,1:2],col=color, pch = 16, cex = 1.5, ylab = "PC2 (16%)", xlab = "PC1 (22%)")
legend("topleft", legend = c("Control", "Ovary", "Testis", "Adrenal"), col = c('#1B9E77','#E7298A','#7570B3','#D95F02'), pch = 19, cex=.9)
dev.off()

a <- autoplot(pca, colour=color) + labs(x="PC1 (22%)", y="PC2 (16%)")

######################## DIFFERENTIAL EXPRESSION USING LIMMA ###############################

#limma

eset <-new("ExpressionSet", exprs= data_matrix)
featureData(eset) <- featureData(data_rma)
eset

targets <- readTargets("phenodata.txt", path ="~/Desktop/Atlas.Array/data/phenoData", row.names = 1, sep="\t")
targets
f <- paste(targets$Target,sep="")
f <- factor(f)
f
design <- model.matrix(~0+f) # Create Design Matrix
colnames(design) <- levels(f)
design
data.fit <- lmFit(eset,design) # Fit Linear Model
names(data.fit)
colnames(design)

#coef=1: adrenal-control
#coef=2: testis-control
#coef=3: ovary-control
#coef=4: adrenal-ovary
#coef=5: adrenal-testis
#coef=6: ovary-testis

contrast.matrix = makeContrasts(Adrenal-Control,Testis-Control,Ovary-Control,Adrenal-Ovary,Adrenal-Testis,Ovary-Testis,levels=design)
contrast.matrix
data.fit.con = contrasts.fit(data.fit,contrast.matrix) # Compute estimated coefficients and standard errors
data.fit.eb = eBayes(data.fit.con) # Bayes Statistics for Differential Expression
data.fit.eb

######### ADRENAL VS CONTROL ############

adrenal_control <-topTable(data.fit.eb, coef=1, number=10000, adjust="fdr", lfc=1) # Extract genes from linear model fit
nrow(adrenal_control) # 2118 genes (1403 previously)

adrenal_control.final <- adrenal_control[,c(1,8, 19:24)]
results.adrenal_control <- separate(adrenal_control.final, geneassignment, into = c("NM", "SYMBOL"), sep = " // ", extra = "drop") #Pay attention to SEPARATION FORMAT!!!
write.table(results.adrenal_control, "adrenal-control.txt", row.names = TRUE, sep="\t")

######### TESTIS VS CONTROL ############

testis_control <-topTable(data.fit.eb, coef=2, number=10000, adjust="fdr", lfc=1) # Extract genes from linear model fit
nrow(testis_control) # 1517 genes (1035 previously)

testis_control.final <- testis_control[,c(1,8, 19:24)]
results.testis_control <- separate(testis_control.final, geneassignment, into = c("NM", "SYMBOL"), sep = " // ", extra = "drop") #Pay attention to SEPARATION FORMAT!!!
write.table(results.testis_control, "testis-control.txt", row.names = TRUE, sep="\t")

################## HEATMAP TOP50 ADRENAL VS CONTROL ################

color.heatmap <- c(rep('#1B9E77',6),rep('#E7298A',10),rep('#7570B3',20), rep('#D95F02',17))

complete <- expression.separate.combat[complete.cases(expression.separate.combat$SYMBOL),]

complete <- complete[!duplicated(complete[,1]),] # 23307 genes after removing duplicates

rownames(complete) <- NULL

str(complete)

rownames(complete) <- complete[,1]

complete <- complete[,-1]

complete[1:3,]

colnames(complete) <- pd[[2]]

complete <- complete[,c(33,18,22,19:21,23,27,24:26,28:32,36:39,43:47,40:42,48:53,34,35,2:4,8,9,5:7,10:17,1)]

write.table(complete, "COMPLETE.txt", row.names = TRUE, sep="\t")


top.adrenal <- results.adrenal_control[!duplicated(results.adrenal_control[,3]),] # Removing duplicated gene names
nrow(top.adrenal) # 1924 genes

top.adrenal <- top.adrenal[complete.cases(top.adrenal[,3]),] # Removing NA values
nrow(top.adrenal) # 1923 genes (1281 previously)

#write.table(top.adrenal, "adrenal-control-volcano.txt", row.names = TRUE, sep="\t") Using logFC=0


top.adrenal <- top.adrenal[,-c(1:2)]
rownames(top.adrenal) <- NULL
rownames(top.adrenal) <- top.adrenal[,1]
top.adrenal <- top.adrenal[,-1]

top.adrenal <- top.adrenal[order(-top.adrenal$logFC),]
top.adrenal[1:3,]

top50.adrenal <- top.adrenal[1:50,]
top50.adrenal

top50.adrenal_selected <- complete[rownames(complete) %in% rownames(top50.adrenal),]

top50.adrenal_selected[1:3,]


adr.top50 <- top50.adrenal_selected[rownames(top50.adrenal),] #### IMPORTANT TO MAINTAIN ORDER AFTER MATCHING


pdf(file = "Heatmap.Adrenal.Top50.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
cairo_ps(file = "Heatmap.Adrenal.Top50.ps", family = "sans")
pheatmap(adr.top50, color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), scale = "row", cluster_rows = FALSE,
cluster_cols = FALSE,fontsize = 8,border_color="NA", width=7, height=6, filename = "test.pdf")

dev.off()

write.table(adr.top50, file="Heatmap.Adrenal.Top50.txt", row.names=TRUE, sep="\t")


######################### HEATMAP STEROIDOGENESIS ##############################

#Heatmap Steroidogenesis

steroidogenesis <- complete[c("LDLR","SCARB1","DHCR24","STAR","CYP11A1", "FDX1", "HSD3B2", "CYP17A1", "PAPSS2", "FOXO4", "MAP3K15", "INHA", "MGARP",
                              "RMDN2", "GRAMD1B", "SLC16A9", "SLC8B1"),]


pdf(file = "Heatmap-Steroidogenesis.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
svg(file = "Heatmap-Steroidogenesis.svg", width = 7, height = 6, family = "sans")
pheatmap(steroidogenesis, color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), scale = "row", cluster_rows = FALSE,
         cluster_cols = FALSE,fontsize = 8,border_color="NA")
dev.off()

write.table(steroidogenesis, file="Heatmap-Steroidogenesis.txt", row.names=TRUE, sep="\t")

dhx37 <- complete[c("DHX37"),]

svg(file = "DHX37.svg", width = 7, height = 6, family = "sans")
pheatmap(dhx37, color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), scale = "row", cluster_rows = FALSE,
         cluster_cols = FALSE,fontsize = 8,border_color="NA",cellheight =12)
dev.off()
######################### HEATMAP SRY  ##############################

#Heatmap Steroidogenesis

sex <- complete[c("SRY","TXLNGY","AMELY","KALP","ASMTL","ASMTL-AS1",
                  "BCORP1",
                  "BPY2",
                  "CD24",
                  "CD99",
                  "CDY2B",
                  "CSF2RA",
                  "CSPG4P1Y",
                  
                  "DDX3Y",
                  "DHRSX",
                  "DUX4L24",
                  "EIF1AY",
                  "ERVH-6",
                  "GOLGA6L9",
                  "HSFY2",
                  "IL3RA",
                  "IL9R",
                  "KDM5D",
                  "LOC105377221",
                  "NLGN4Y",
                  "P2RY8",
                  "PCDH11X",
                  "PLCXD1",
                  "PPP2R3B",
                  "PRKY",
                  "PRORY",
                  "PRY2",
                  "RBMY1A3P",
                  "RBMY1B",
                  "RBMY2BP",
                  "RFTN1",
                  "RPS4Y1",
                  "RPS4Y2",
                  "SHOX",
                  "SLC25A6",
                  "SPRY3",
                  "TBL1Y",
                  "TEKT4P2",
                  "TGIF2LX",
                  "TMSB4Y",
                  "TSPY2",
                  "TTTY10",
                  "TTTY11",
                  "TTTY12",
                  "TTTY13",
                  "TTTY14",
                  "TTTY1B",
                  "TTTY2",
                  "TTTY5",
                  "TTTY6",
                  "TTTY7",
                  "TTTY8",
                  "TTTY9B",
                  "USP9Y",
                  "UTY",
                  "VAMP7",
                  "VCX",
                  "XKRY",
                  "ZFY"),]
final.sex <- sex[complete.cases(sex$Spine),]
sex2 <- complete[c("SRY","TXLNGY"),]
pdf(file = "SRY-TXLNGY.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)
svg(file = "CHRY-heatmap.svg", width = 7, height = 6, family = "sans")
pheatmap(sex2, color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), scale = "row", cluster_rows = FALSE,
         cluster_cols = FALSE,fontsize = 12,cellheight=12,border_color="NA")
pheatmap(final.sex, color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), scale = "row", cluster_rows = TRUE,
         cluster_cols = FALSE,fontsize = 8,border_color="NA")
dev.off()

write.table(final.sex, file="CHRY.txt", row.names=TRUE, sep="\t")

####################### HEATMAP TRANSCRIPTION ##############################

#### I Removed ZNF676 (THERE ARE TWO PROBES WITH VERY DIFFERENT VALUES)

transcription<- complete[c("FOXO4","TBX3","RFX6","FOSL2","CCNE1","SCAP","NMI","TCEA3","RORA",
                           "LIN28B", "ARNTL", "CREM", "NR0B1", "NR5A1", "CITED1", "TSPYL2", "MAMLD1", "TCF21", "APOBEC3G", "NOTCH2",
                           "WT1", "BNC1", "NANOG", "POU5F1", "ZFP42", "GATA4", "BNC2", "TFAP2C", "LIN28A",
                           "PRDM1", "BRDT", "TFCP2L1", "LHX9", "LEF1", "ZNF208", "NR6A1", "CDCA7L", "EMX2"),]  

transcription[1:3,]

pdf(file = "Heatmap-Transcription.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)

pheatmap(transcription, color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), scale = "row", cluster_rows = FALSE,
         cluster_cols = FALSE,fontsize = 8,border_color="NA")
dev.off()

write.table(transcription, file="Heatmap-Transcription.txt", row.names=TRUE, sep="\t")


############################ HEATMAP SRY 16 FACTORS ###################################

sry<- complete[c("ASPN","GSTA1","CITED1","ANKRD18A","MOCOS","HIST1H2AA","MAPK4","NOSTRIN","G6PD",
                 "RASSF2", "SLC52A3", "KEL", "ZNF280B", "PRPS2", "SOX9", "INHBB"),]  

sry

pdf(file = "Heatmap-SRY.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)

pheatmap(sry, color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), scale = "row", cluster_rows = FALSE,
         cluster_cols = FALSE,fontsize = 8,border_color="NA", cellheight =12)
dev.off()

write.table(transcription, file="Heatmap-SRY.txt", row.names=TRUE, sep="\t")

############################ HEATMAP 45 FACTORS ###################################

factors <- complete[c("CYP17A1","CYP11A1","STAR","MGARP","GSTA1","MAP3K15","HPGD","C7","SERPINA5", "FDX1P1","SCARB1","FDX1","PAPSS2","DHCR24","MT2A","DLK1","INHA","MGST1","APOA1","FOXO4","APOC1","HSD3B2","GRAMD1B","RMDN2","MSMO1","AS3MT","SLC46A3","LDLR","SLC16A9","DHCR7","SLC8B1","GSTA3","POR","ACSS1","COL15A1","NPC1","AVPI1","PDK4","ABCA5","FRRS1","IDI1","VCAM1","FLJ38894","GPD1L","HMGCR"),] 

pdf(file = "45Factors.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)

pheatmap(factors, color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), scale = "row", cluster_rows = TRUE,
         cluster_cols = FALSE,fontsize = 8,border_color="NA")
dev.off()

write.table(factors, file="Heatmap-45Factors.txt", row.names=TRUE, sep="\t")
######### OVARY VS CONTROL ############

ovary_control <-topTable(data.fit.eb, coef=3, number=10000, adjust="fdr", lfc=1) # Extract genes from linear model fit
nrow(ovary_control) # 1665 genes (1214 previously)

ovary_control.final <- ovary_control[,c(1,8, 19:24)]
results.ovary_control <- separate(ovary_control.final, geneassignment, into = c("NM", "SYMBOL"), sep = " // ", extra = "drop") #Pay attention to SEPARATION FORMAT!!!
write.table(results.ovary_control, "ovary-control.txt", row.names = TRUE, sep="\t")


####################### HEATMAP TOP50 OVARY VS CONTROL ##################################

#### REMOVING ZNF676 ######

top.ovary <- results.ovary_control[!duplicated(results.ovary_control[,3]),] # Removing duplicated gene names
nrow(top.ovary) # 1365 genes

top.ovary <- top.ovary[complete.cases(top.ovary[,3]),] # Removing NA values
nrow(top.ovary) # 1364 genes (1025 previously)

top.ovary <- top.ovary[,-c(1:2)]
rownames(top.ovary) <- NULL
rownames(top.ovary) <- top.ovary[,1]
top.ovary <- top.ovary[,-1]

top.ovary <- top.ovary[order(-top.ovary$logFC),]
top.ovary[1:3,]

top50.ovary <- top.ovary[1:50,]
top50.ovary

top50.ovary_selected <- complete[rownames(complete) %in% rownames(top50.ovary),]

top50.ovary_selected[1:3,]

ovr.top50 <- top50.ovary_selected[rownames(top50.ovary),] #### IMPORTANT TO MAINTAIN ORDER AFTER MATCHING

ovr.top50 <- ovr.top50[-c(10),]


pdf(file = "Heatmap.Ovary.Top50.pdf", width = 7, height = 6, family = "Helvetica", useDingbats = FALSE)

pheatmap(ovr.top50, color = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255), scale = "row", cluster_rows = FALSE,
         cluster_cols = FALSE,fontsize = 8,border_color="NA")

dev.off()

write.table(ovr.top50, file="Heatmap.Ovary.Top50.txt", row.names=TRUE, sep="\t")

######### OVARY VS TESTIS ############

ovary_testis <-topTable(data.fit.eb, coef=6, number=10000, adjust="fdr", lfc=1) # Extract genes from linear model fit
nrow(ovary_testis) # 294 genes 

ovary_testis.final <- ovary_testis[,c(1,8, 19:24)]
results.ovary_testis <- separate(ovary_testis.final, geneassignment, into = c("NM", "SYMBOL"), sep = " // ", extra = "drop") #Pay attention to SEPARATION FORMAT!!!
write.table(results.ovary_testis, "ovary-testis.txt", row.names = TRUE, sep="\t")

top.ovary_testis <- results.ovary_testis[!duplicated(results.ovary_testis[,3]),] # Removing duplicated gene names
nrow(top.ovary_testis) # 274 genes

top.ovary_testis <- top.ovary_testis[complete.cases(top.ovary_testis[,3]),] # Removing NA values
nrow(top.ovary_testis)# 273 genes
write.table(top.ovary_testis, "ovary-testis.NODUP.txt", row.names = TRUE, sep="\t")

