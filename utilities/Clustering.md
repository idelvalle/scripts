EXPRESSION DATA MATRIX
======================

Gene expression data are usually presented in an expression matrix. Each column represents all the gene expression levels from a single experiment, and each row represents the expression of a gene across all experiments. Each element is a log ratio. The log ratio is defined as log2 (T/R), where T is the gene expression level in the testing sample, R is the gene expression level in the reference sample.

CLUSTERING
==========

Clustering is a data exploratory technique used for discovering groups or pattern in a dataset.  
There are two standard clustering strategies: partitioning methods and hierarchical clustering.

-	**K-means** clustering, in which, each cluster is represented by the center or means of the data points belonging to the cluster.  
-	**K-medoids** clustering or PAM (Partitioning Around Medoids), in which, each cluster is represented by one of the objects in the cluster.  

There is also a variant of *PAM* named *CLARA* (Clustering Large Applications) which is used for analyzing large data sets.

There are two main Algorithms for Clustering analysis:

Hierarchical Clustering
-----------------------

In hierarchical clustering, genes with similar expression patterns are grouped together and are connected by a series of branches (clustering tree or dendrogram). Experiments with similar expression profiles can also be grouped together using the same method.

-	**How to determine the similarity between two genes?**  

We calculate the distance between two expression vectors. A Gene Expression Vector consists of the expression of a gene over a set of experimental conditions. Distances are afterwards calculated using Euclidean or Manhattan methods.

-	**How to determine the similarity between clusters?**  

The method for determining cluster-to-cluster distance is called linkage method.

**1. Single Linkage**: Cluster-to-cluster distance is defined as the minimun distance between members of one cluster and members of the another cluster.  
**2. Complete Linkage**: Cluster-to-cluster distance is defined as the maximum distance between members of one cluster and members of the another cluster.  
**3. Average Linkage**: Cluster-to-cluster distance is defined as the average distance between members of one cluster and members of the another cluster.

There is no theoretical guideline for selecting the best linkage method. In practice, people almost always use the average linkage method.

K-Means/K-Medians Clustering
----------------------------

K-means clustering is the simplest and the most commonly used partitioning method for splitting a dataset into a set of k groups (i.e. clusters). It requires the analyst to specify the number of optimal clusters to be generated from the data.

In k-means clustering, each cluster is represented by its center (i.e, centroid) which corresponds to the mean of points assigned to the cluster.

A good clustering will generate clusters with a high **intra-class similarity** and a low **inter-class** similarity.

Hence, the basic idea behind K-means clustering consists of **defining clusters so that the total intra-cluster variation (known as total within-cluster variation) is minimized.**

Median is the middle number, i.e. the middle of the distribution. Median is more robust against outliers

For an odd number of numbers, the median is simply the middle number. For example, the median of 2, 4 and 7 is 4.

For an even number of numbers, the median is the average of the two middle numbers. Thus, the median of the numbers 2, 4, 7, 12 is (4+7)/2 = 5.5.

**1.** Specify number of clusters  
**2.** Randomly assign genes to clusters  
**3.** Calculate mean/median expression profile of each cluster  
**4.** Shuffle genes among clusters such that each gene is now in the cluster whose mean/median expression profile is the closest to that gene's expression profile  
**5.** Repeat steps 3 and 4 until genes cannot be shuffled around any more OR a user-specified number of iterations has been reached

**R Function for Cluster Analysis**

```R
kmeans(x, centers, iter.max = 10, nstart =1)
# x: Numeric matrix, data frame or vector  
# centers: Number of clusters  
# nstart: The number of random starting partitions when centers is a number.
```

**It’s strongly recommended to compute k-means clustering with a large value of nstart such as 25 or 50, in order to have a more stable result.**

CLUSTER ANALYSIS
================

A) Using MeV (Multi Experiment Viewer)
--------------------------------------

Follow MultiExperiment Viewer Quickstart Guide

B) ConsensusClusterPlus (BioConductor)
--------------------------------------

The input data format is a matrix where columns are samples (items), rows are features and cells are numerical values.  
I reduce the dataset to the top 5,000 most variable genes, measured by median absolute deviation
