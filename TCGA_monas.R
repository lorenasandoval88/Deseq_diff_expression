#----packages-----------------------------------------
#BiocManager::install("TCGAutils")
# BiocManager::install("DESeq2")
#BiocManager::install('PCAtools')

library(TCGAutils)
library(curatedTCGAData)
library(MultiAssayExperiment)
library(RTCGAToolbox)
library(BiocFileCache)
library(rtracklayer)
library(R.utils)
library("DESeq2")
library(ggplot2)
library(dplyr)
library(biomaRt)
library(PCAtools)
library(DT)
library(IHW)
require(TCGAbiolinks)
library(SummarizedExperiment)
#-----------------------------------------------------

#set working directory
setwd("C:/Users/sandovall2/Downloads")

#The data can be loaded with the following command
tcga_data = readRDS(file = "tcga_dataAll.RDS")

# We will decide to keep genes that are expressed in at least 10 samples.
is_expressed <- assay(tcga_data) >= 5
sum(is_expressed[128,])

# rowSums will give the sum of each row in a matrix and return the result as a vector.
df <- data.frame(Expressed = rowSums(is_expressed))
ggplot(df, aes(x=Expressed)) + geom_bar()

keep <- rowSums(assay(tcga_data) >= 5) >= 5
table(keep)
tcga_data <- tcga_data[keep,]


# remove samples with unknown subtype (132 samples removed incl aa, asian, white, 1083 remain)
idx = which(!is.na(tcga_data@colData$paper_BRCA_Subtype_PAM50))
tcga_data = tcga_data[,idx]

# remove samples with unknown race (1 sample removed)
idx = which(!is.na(tcga_data@colData$race))
tcga_data = tcga_data[,idx]

## Subset 40 eu and 40 aa samples to reduce computational burden
idx_eu = sample(which(tcga_data@colData$race=="white"),40)
idx_aa = sample(which(tcga_data@colData$race=="black or african american"),40)

data = tcga_data[,c(idx_eu,idx_aa)]

# Pie Chart full data of all subtypes and races___fig1_pie

mytable <- table(tcga_data@colData$race)
sum(mytable)
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls,col=rainbow(length(lbls)),radius = 0.8,
    main="Pie Chart of TCGA-BRCA by race\n(1214 samples)")

# mytable <- table(tcga_data@colData$"paper_BRCA_Subtype_PAM50",useNA = "ifany")
# lbls <- paste(names(mytable), "\n", mytable, sep="")
# pie(mytable, labels = lbls,col=rainbow(length(lbls)),radius = 0.8,
#     main="Pie Chart of TCGA-BRCA by subtype\n(1214 samples)")

# mytable <- table(tcga_data@colData$"gender")
# lbls <- paste(names(mytable), "\n", mytable, sep="")
# pie(mytable, labels = lbls,col=rainbow(length(lbls)),radius = 0.8,
#     main="Pie Chart of TCGA-BRCA by gender\n(1214 samples)")

# # Pie Chart of basal, female by race (eu and aa)

idx = which(tcga_data@colData$paper_BRCA_Subtype_PAM50=="Basal"
            & tcga_data@colData$gender=="female")
basal_data = tcga_data[,idx]
mytable <- table(basal_data@colData$race)
sum(mytable)
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls,col=rainbow(length(lbls)),radius = 0.8,
    main="Pie Chart of TCGA-BRCA basal subtype by race\n (190 sample size)")
#--------------------------------------------------------------------

# Make subtype expression data with subtype as columnnames and genes as rownames
subtype_mat = assay(data)
rownames(subtype_mat) = rowData(data)$gene_id
colnames(subtype_mat) = colData(data)$paper_BRCA_Subtype_PAM50



#--------------------------------------------------------------------

# Make race expression data with race as columnnames and genes as rownames
race_mat = assay(data)
rownames(race_mat) = rowData(data)$gene_id
colnames(race_mat) = colData(data)$race


#--------------------------------------------------------------------
#PAM 50 list to find probes/rows in expression set
pam50_list = c("UBE2C","PTTG1","MYBL2","BIRC5","CCNB1","TYMS","MELK",
               "CEP55","KNTC2","UBE2T","RRM2","CDC6","ANLN","ORC6L",
               "KIF2C","EXO1","CDCA1","CENPF","CCNE1","MKI67","CDC20",
               "MMP11","GRB7","ERBB2","TMEM45B","BAG1","PGR","MAPT",
               "NAT1","GPR160","FOXA1","BLVRA","CXXC5",
               "ESR1","SLC39A6","KRT17","KRT5","SFRP1",
               "BCL2","KRT14","MLPH","MDM2","FGFR4","MYC",
               "MIA","FOXC1","ACTR3B","PHGDH","CDH3","EGFR")

# 47 genes found from PAM50 list
pam50 <- subtype_mat[ rownames(subtype_mat) %in% pam50_list, ]

#--------------------------------------------------------------------
# PCA based on pam50 and 40 samples
idxA = which(colnames(pam50)=="LumA")
idxBasal = which(colnames(pam50)=="Basal")
idxB = which(colnames(pam50)=="LumB")
idxHer2 = which(colnames(pam50)=="Her2")
idxNorm = which(colnames(pam50)=="Normal")

pcaData = pam50[,c(idxA[1:14],idxBasal[1:14],idxNorm,idxHer2,idxB)]

## scale
pcaData = scale(pcaData,center=T, scale = T)

# vsd <- vst(as.matrix(assay(round(pcaData))), nsub=nrow(pcaData))
# Z <- scale(vsd)
# pc <- prcomp(t(vsd))
pc <- prcomp(t(pcaData))

# load packages
#https://tavareshugo.github.io/data-carpentry-rnaseq/03_rnaseq_pca.html
library(tidyverse)

# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- pc$x
pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- pc$sdev^2
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

# Calc how much variation in the original data each
#principal component accounts for
pca.var = pc$sdev^2
pca.var.per = round(pca.var/sum(pca.var)*100,1)

# print the result
pc_scores
Subtype =factor(colnames(pcaData))
pc_scores %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2,
             color = Subtype)) +
  #geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA Graph")+
  geom_point()

#--------------------------------------------------------------------
# Heatmap based on pam50

#transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts,
#and which normalizes with respect to library size.
pam = rlog(round(pam50))
heatmap(pam)
#kmeans
cl <- kmeans(t(pam), 5)
plot(pam, col = cl$cluster,
     main= "K-means clustering with 5 clusters of sizes 191, 416, 467, 4, 5")
points(cl$centers, col = 1:5, pch = 8, cex=2)

#-------------------------------------------------------------
#PCA on AA vs EU basal subtypes

# Normalization 
subtypes_scaled<-scale(subtype_mat,T,T)

#The first five eigenvalues from the correlation matrix
eigen(cor(as.matrix(subtypes_scaled)))$values[1:30]

## cluster analysis on the component scores of the principal components
## The percentages of variance explained by the first two components canbe computed

sum(eigen(cor(as.matrix(subtypes_scaled)))$values[1:2])/1083*100
#the variance explained by the first 2 eigenvalues
# is 67%
#Hence, the data can allow for a reduction in dimensions from
#1083 to 2 while still capturing over 67% of the total variance.


#It can be checked that all the correlations between the patients are positive.
-eigen(cor(as.matrix(subtypes_scaled)))$vec[,1:2]
# There is positive and negative correlations between samples
#in eigenvector


#the genes with the largest expression values from the first principal
#component can be printed

pca <- princomp(as.matrix(subtypes_scaled), center = TRUE, cor=TRUE, scores=TRUE)
o <- order(pca$scores[,2])
top100genes = rowData(data)$gene_id[o[1:100]]
top100genes  = subset(subtypes_scaled, rownames(subtypes_scaled) %in% top100genes)
rownames(top100genes)


#BIPLOT
#Biplot. A useful manner to plot both genes (cases) and patients
#(variables) is the biplot, which is based on a two-dimensional approximation
#(reduction) of the data very similar to principal components analysis
#(PCA).
biplot(princomp(as.matrix(subtypes_scaled),
                cor=TRUE),pc.biplot=TRUE,
       cex=0.5,expand=0.7)
#no sig clusters in samples or genes found in biplot


#In order to select those that have an experimental
#effect, we use a two-sample t-test:
pt <- apply(top100genes, 1, function(x) t.test(x ~ as.factor(colnames(top100genes)))$p.value)
oo <- o[pt[o]<0.01]

#-------------------------------------------------------------
#PCA by race (180 aa and 746 eu samples)

# Normalization 
# center is TRUE then centering is done by subtracting the column means (omitting NAs) of x from their corresponding columns
# If scale is TRUE then scaling is done by dividing the (centered) columns of x by their standard deviations if center is TRUE, and the root mean square otherwise
race_scaled<-scale(race_mat,T,T)

#The first five eigenvalues from the correlation matrix
eigen(cor(as.matrix(race_scaled)))$values[1:30]

## cluster analysis on the component scores of the principal components
## The percentages of variance explained by the first two components can be computed

sum(eigen(cor(as.matrix(race_scaled)))$values[1:2])/926*100
#the variance explained by the first 2 eigenvalues
# is 67%
#Hence, the data can allow for a reduction in dimensions from
#1083 to 2 while still capturing over 67% of the total variance.


#It can be checked that all the correlations between the patients are positive.
-eigen(cor(as.matrix(race_scaled)))$vec[,1:2]
# There are positive and negative correlations between samples
#in eigenvector 2


#the genes with the largest expression values from the first principal
#component can be printed

pca_race <- princomp(as.matrix(race_scaled), center = TRUE, cor=TRUE, scores=TRUE)
o <- order(pca$scores[,2])
top950genes = race_scaled[o[1:950],]
rownames(top950genes)


#BIPLOT
#Biplot. A useful manner to plot both genes (cases) and patients
#(variables) is the biplot, which is based on a two-dimensional approximation
#(reduction) of the data very similar to principal components analysis
#(PCA).
colnames(top950genes)[colnames(top950genes) == 'white'] <- 'eu'
colnames(top950genes)[colnames(top950genes) == 'black or african american'] <- 'aa'
biplot(princomp(as.matrix(top950genes),
                cor=TRUE),pc.biplot=TRUE,
       cex=0.5,expand=0.7)
#no sig clusters in samples or genes found in biplot


#In order to select those that have an experimental (results in 299 genes with pval<0.001)
#effect, we use a two-sample t-test:
pt <- apply(top950genes, 1, function(x) t.test(x ~ as.factor(colnames(top950genes)))$p.value)
oo <- o[pt[o]<0.001 & !is.na(pt[o])]


#This yields 306 genes, of which the row numbers are selected in the vector
#oo. In order to identify genes in directions of large variation we use the scores
#on the first two principal components:
Z <- race_scaled
K <- eigen(cor(Z))
P <- Z %*% -K$vec[,1:2]
leu <- data.frame(P[oo,], row.names= oo)
plot(leu,xlim=c(-10,15), ylim=c(-10,10), pch=19, cex=1.2, xlab="Principal Component 1",
     ylab="Principal Component 2", col="darkgreen")
text(x = leu$X1, y=leu$X2, labels=rownames(leu), pos = 1, col="blue") #FIX o1
fac <- as.integer(oo %in% o1) + 2 * as.integer(oo %in% o2) + 3 * as.integer(oo %in% o3)
text(x = leu$X1, y=leu$X2, labels=fac, pos = 3, col="red")


#The scores on the first two principal components of the selected genes are
# stored in the data frame leu. From the plotted component scores in Figure
# 7.14, it seems that there are several sub-clusters of genes. The genes that
# belong to these clusters can be identified by hiearchical cluster analysis:

cl <- hclust(dist(leu,method="euclidian"),method="single")
plot(cl,
     lwd=3,
     col="blue",
     col.axis = "brown",
     ylab="Distance",
     xlab="Clustering of the expression of genes",
     hang=-1,
     main=NA,
     sub=NA,
     axes=FALSE)
axis(side = 2, at = seq(0, 5, 1), col = "brown",labels = TRUE, lwd = 4)


# From the resultant tree (dendrogram) in Figure 7.15 various clusters of genes
# are apparent that also appear in Figure 7.14.3 The ordered genes can be
# obtained from the object cl as follows:
a <- as.integer(rownames(leu)[cl$order])
for (i in 1:length(a)) { cat(a[i],rownames(race_scaled)[a[i]],"\n") }


#------------------------------------------------PCA
#https://www.kaggle.com/code/bharatrai0/m9-pca/notebook
#--------------------------------race pca
library (ISLR)
labs <- colnames(race_mat)
d <- t(race_mat) #when t() to samples as rows, get 0 var cols
dim(d)
d = d[ , which(apply(d, 2, var) != 0)]
dim(d)

pc <- prcomp(d , scale=TRUE)
Cols=function (vec ){cols=rainbow (length (unique (vec )))
return (cols[as.numeric (as.factor (vec))])}
plot(pc$x[,1:2], col = Cols(labs), pch =19, xlab ="Z1",ylab="Z2")
#--------------------------------Subtype pca
library (ISLR)
labs <- colnames(subtype_mat)
d <- t(subtype_mat) #when t() to samples as rows, get 0 var cols
dim(d)
d = d[ , which(apply(d, 2, var) != 0)]
dim(d)

pc <- prcomp(d , scale=TRUE, center=TRUE)
Cols=function (vec ){cols=rainbow (length (unique (vec )))
return (cols[as.numeric (as.factor (vec))])}
plot(pc$x[,1:2], col = Cols(labs), pch =19, xlab ="PC1",ylab="PC2")

plot(pc) #varaance explained by each PC vector
#--------------------------------------------------
#https://github.com/BarryDigby/TCGA_Biolinks/blob/master/TCGA_Biolinks.Rmd
## remove columns we dont need, keep counts
mrna_meta <- data$race
mrna_meta <- cbind(mrna_meta, data$paper_BRCA_Subtype_PAM50)
mrna_df <- assay(data)


## tidy matrix colnames 
delim_fn = function(x, n, i){
  do.call(c, lapply(x, function(X)
    paste(unlist(strsplit(X, "-"))[(n+1):(i)], collapse = "-")))
}
colnames(mrna_df) <- delim_fn(x = colnames(mrna_df), n = 0, i = 4)
mrna_meta <- as.data.frame(mrna_meta)
mrna_df <- as.data.frame(mrna_df)

colnames(mrna_meta) <- c("race","paper_BRCA_Subtype_PAM50")


## 3 gene names are duplicated 'PINX1' 'TMSB15B' 'MATR3'
## append numerics to names to make unique
rownames(mrna_df) <- make.names(rownames(mrna_df), unique = T)

## DESeq2 Analysis
mrna_dds <- DESeqDataSetFromMatrix(round(mrna_df), colData = mrna_meta, design = ~ paper_BRCA_Subtype_PAM50)
mrna_dds$race <- relevel(as.factor(mrna_dds$race), ref = "white")

mrna_dds <- DESeq(mrna_dds)
resultsNames(mrna_dds)

## DESeq2 results
# get plots next: https://www.youtube.com/watch?v=jLVI-gPEff8
mrna_res <- results(mrna_dds, filterFun = ihw, alpha = 0.05, name = "paper_BRCA_Subtype_PAM50_Normal_vs_Basal")
summary(mrna_res)
mrna_res_df <- as.data.frame(mrna_res)
mrna_upreg <- get_upregulated(mrna_res)
mrna_downreg <- get_downregulated(mrna_res)

## Write results for plots and analysis
#https://www.youtube.com/watch?v=jLVI-gPEff8
mrna_counts <- counts(mrna_dds, normalized = T)
write.table(mrna_counts, "~/Desktop/TCGA/mRNA/results/mRNA_norm.counts.txt", quote = F, sep = "\t")
