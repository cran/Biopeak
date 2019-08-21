## ----load heat-----------------------------------------------------------
library(Biopeak)
# load the heat-shock data
exprmat <- heat
# Show first rows and columns of the expression matrix
head(exprmat)
# Transform data frame to numeric matrix
exprmat <- as.matrix(exprmat)

## ----peakdetection-------------------------------------------------------
# define condition series
series <- c(37,40,41,42,43)
# Find the expression limit for the 50 % highest expressing genes and selec them
exprlim <- quantile(exprmat, probs <- seq(0,1,0.01))[51]
exprmat <- exprmat[which(rowMeans(exprmat) > exprlim),]
# run the peak detection algorithm
peakdet <- peakDetection(exprmat, series, type ='rnaseq', actstrength = 1.5, prominence = 1.3,
                         bgcorr = F)

## ----output--------------------------------------------------------------
# number of genes detected by the peakDetection function
length(peakdet$peakgenes)
# filter for genes which peak at 43 °C
peakdet$peakgenes[which(peakdet$peakloc == 43)]
# return the peakheight and location of the gene CDK13
c(peakdet$peakheight[which(peakdet$peakgenes == 'CDK13')],
  peakdet$peakloc[which(peakdet$peakgenes == 'CDK13')])
# assess the main peak neighborhood for sustained activation
peakdet$neighbors[which(peakdet$peakgenes == 'CDK13'),]

## ----plotexpression,fig.cap="\\label{fig:plotexpression}Figure 1: Individual gene expression signal for CDK13. The dashed line marks the main peak location.", fig.align = 'center', fig.height = 4, fig.width = 5----
plotExpression(exprmat,'CDK13',series,peakdet)

## ----gettcormat,fig.cap="\\label{fig:getcormat}Figure 2: Bi-clustered heatmap of the gene correlation matrix", fig.align = 'center'----
dev.off()
corobjects <- getCormat(peakdet, exprmat, method = 'spearman')
# extract heatmap object
corheatmap <- corobjects$hm
# extract re-ordered correlation matrix
cormatrix <- corobjects$hm_cormat
# inspect the first 5 x 5 gene-wise correlation of the correlation matrix returned by plotCormat
cormatrix[1:5,1:5]

## ----findclusters, fig.height = 10, fig.width = 6,fig.cap="\\label{fig:findclusters}Figure 3: Clusters of biomarkers with similar temporal regulation based on three different algorithms (hclust, dbscan and kmeans).", fig.align = 'center'----
# display all plots in one graph
par(mfrow = c(3,1))
# cluster exploration using hierarchical clustering based on 3 clusters
clusters <-findClusters(peakdet, exprmat, method = 'hclust', clusters = 3)
# cluster exploration using dbscan (density-based) with an epsilon parameter of 0.02
clusters <- findClusters(peakdet, exprmat, method = 'dbscan', eps = 0.1)
# cluster exploration using kmeans with a maximum of 5 clusters to be assigned
clusters <- findClusters(peakdet, exprmat, maxclusters = 3, method = 'kmeans')

## ----plotheatmap,fig.cap="\\label{fig:plotheatmap}Figure 4: Heatmap of the second cluster returned by the findClusters function. The second cluster reflects primarily genes involved in the immediate LPS-induced response.",fig.height = 8, fig.align = 'center'----
heatmap <- plotHeatmap(peakdet, exprmat, clustermembers = clusters$clustermembers[[3]])

## ----saveOutput----------------------------------------------------------
saveOutput(peakdet,file.path(tempdir(),'heat_out.txt'))

