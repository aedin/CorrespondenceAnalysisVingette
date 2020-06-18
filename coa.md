---
title: "Relationship between glmPCA and Correspondence Analysis"
author: "Aedin Culhane"
date: "June 17, 2020"
output: html_document
bibliography: references.bib  
---

In the Townes et al., (2019), generalized principal component analysis (GLM-PCA), a generalization of PCA to exponential family likelihoods. In practice, Townes et al. report that both the multinomial and Dirichlet-multinomial are computationally intractable on single cell data and may be approximated by the Poisson and negative binomial likelihoods, and the paper focuses on Poisson glmPCA due to its performance. One disadvantage of GLM-PCA is it depends on an iterative algorithm to obtain estimates for the latent factors and is at least ten times slower than PCA. Therefore a fast approximation with deviance residuals or Pearson residuals are propsed. 

Single cell RNAseq data are counts and thus are naturally Poisson, and they are not normalized by log or scale transformation (as we recently described Hsu L & Culhane AC, 2020). PCA or decomposition of the Pearson residuals is a excellent method for exploratory analysis of count data and is also know as correspondence analysis. Because of its speed (which is comparable to PCA) and its performance it has been adopted by many fields. Correspondence analysis was originally proposed in 1935 (Hirschfeld, 1935), and was later developed by Benzécri as part of the French duality diagram framework for multivariate statistics (Benzécri, 1973; Holmes, 2008). His trainee Michael Greenacre later popularized its use with large, sparse count data in diverse settings and disciplines, including linguistics, business and marketing research, and archaeology (Greenacre, 1984, 2010). In numerical ecology, it is commonly applied to analyzing species abundance count matrices (M. Greenacre, 2013; Legendre et al., 1998). It has been successfully  applied to microarray transcriptomics data (Fellenberg et al., 2001). and a Bioconductor version was implemented in made4 Bioconductor package by Culhane et al., (2002) and the method was also used in joint and multi natrix factorization methods (Culhane et al.,  2003; Meng et al., 2014). Its first application on single cell RNA seq was in analysis of metagenomic scRNAseq microbiome census data data (McMurdie & Holmes, 2013). Correspondence analysis is available in the ade4 and vegan R packages and was listed in a review of single and multi table matrix factorization approaches (Meng C, Zeleznik OA. et al., 2016).

Correspondence analysis (COA) is considered a dual-scaling method, because both the rows and columns are scaled prior to singular value decomposition (svd). Lauren Hsu provided an excellent description of Correspondence analysis applied to one or an extension to joint decomposition of multi single cell RNA seq in her masters thesis (Biostatistics, Harvard School of Public Health, 2020) and has implemented an fast implementation of COA that supports sparse matrix operations and scalable svd (IRBLA) in her Bioconductor package Corral.   


Hafemeister and Satija in their characterization of the SCTransform method also propose
that Pearson residuals of a regularized negative binomial model (a generalized linear model with sequencing depth as a covariate) could be used to remove technical characteristics while preserving biological heterogeneity, with the residuals used directly as input for downstream analysis, in place of log-transformed counts. (Hafemeister & Satija, 2019) Neither Townes et al., (2019) or Hafemeister and Satija (2019) paper cites the extensive literature on correspondence analysis and its application to count data in other fields. However it is not the first time  correspondence analysis has been “rediscovered”, its facts its frequently rediscovered (Abdi & Valentin, 2007; Greenacre, 2010), and each rediscovery enforces the importance of it to the field 

## Dataset described in the glmPCA vignette

The function simData() created a simulated dataset that is provided in the glmPCA vignette and was orginally provided by [Jake Yeung](https://github.com/willtownes/scrna2019/issues/2). It has 4989 rows and 150 columns, 3 biological groups (clusters) of 50 cells each and 2 batches. 5000 rows are created but some are filtered to create a matrix of 4989 rows.




```r
simData<-function() {
    set.seed(202)
  ngenes <- 5000 #must be divisible by 10
  ngenes_informative<-ngenes*.1
  ncells <- 50 #number of cells per cluster, must be divisible by 2
  nclust<- 3
  # simulate two batches with different depths
  batch<-rep(1:2, each = nclust*ncells/2)
  ncounts <- rpois(ncells*nclust, lambda = 1000*batch)
  # generate profiles for 3 clusters
  profiles_informative <- replicate(nclust, exp(rnorm(ngenes_informative)))
  profiles_const<-matrix(ncol=nclust,rep(exp(rnorm(ngenes-ngenes_informative)),nclust))
  profiles <- rbind(profiles_informative,profiles_const)
  # generate cluster labels
  clust <- sample(rep(1:3, each = ncells))
  # generate single-cell transcriptomes 
  counts <- sapply(seq_along(clust), function(i){
    rmultinom(1, ncounts[i], prob = profiles[,clust[i]])
  })
  rownames(counts) <- paste("gene", seq(nrow(counts)), sep = "_")
  colnames(counts) <- paste("cell", seq(ncol(counts)), sep = "_")
  # clean up rows
  Y <- counts[rowSums(counts) > 0, ]
  sz<-colSums(Y)
  Ycpm<-1e6*t(t(Y)/sz)
  Yl2<-log2(1+Ycpm)
  z<-log10(sz)
  pz<-1-colMeans(Y>0)
  cm<-data.frame(total_counts=sz,zero_frac=pz,clust=factor(clust),batch=factor(batch))
  return(Y)
}
mat=simData()
dim(mat)
```

```
## [1] 4989  150
```



```r
mat=simData()
dim(mat)
```

```
## [1] 4989  150
```

#glmPCA  
glmPCA is iterative,so I set seed so its reproducible. I discuss about the iterative nature of glmPCA below


```r
library(glmpca)
set.seed(50)
system.time(res1<-glmpca(mat,L=2,fam="poi"))
```

```
##    user  system elapsed 
##   1.857   0.929   2.802
```

```r
dim(res1$factors)  # Column scores dim=c(150,   2)
```

```
## [1] 150   2
```

```r
dim(res1$loadings) #row scoress dim =c(4989 ,   2)
```

```
## [1] 4989    2
```

```r
plot(res1$factors$dim1, res1$factors$dim2)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)


#CORRAL
Correl from Lauren Hsu provides a fast implementation correspondence analysis in the function corral. Corral accepts input formats SingleCellExperiment (Bioconductor), SparseMatrix, matrix, data.frame etc


```r
library(corral)
system.time(res2<-corral(mat,ncomp=2))
```

```
##    user  system elapsed 
##   0.162   0.040   0.212
```

```r
res2
```

```
## corral output summary===========================================
##   Output "list" includes standard coordinates (SCu, SCv),
##   principal coordinates (PCu, PCv), & SVD output (u, d, v)
## Variance explained----------------------------------------------
##                           PC1  PC2
## percent.Var.explained    0.02 0.01
## cumulative.Var.explained 0.02 0.03
## 
## Dimensions of output elements-----------------------------------
##   Singular values (d) :: 2
##   Left singular vectors & coordinates (u, SCu, PCu) :: 4989 2
##   Right singular vectors & coordinates (v, SCv, PCv) :: 150 2
##   See corral help for details on each output element.
## ================================================================
```

```r
res2$PCu[1:2,]  # Row coordinates  dims =c(4989,10)
```

```
##           [,1]      [,2]
## [1,] 0.6098867 0.6836028
## [2,] 0.1179244 0.2704358
```

```r
res2$PCv[1:2,]  # Col coordinates  dims=c(150,10)
```

```
##            [,1]        [,2]
## [1,] -0.2770690 0.007619298
## [2,] -0.3045123 0.041767103
```

```r
plot(res2$PCv[,1], res2$PCv[,2])
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

Correpondence analysis is also implemented in ade4 (and made4 which is an extension of ade4). The implementation of Correpondence analysis in ade4 is fast, but is slightly slower than corral.  But this is a small dataset, far smaller than a typical scRNAseq data. 

## dudi.coa

```r
system.time(res3<-ade4::dudi.coa(mat,scannf = FALSE, n=2))
```

```
##    user  system elapsed 
##   0.419   0.076   0.506
```

```r
# user  system elapsed 
# 0.382   0.052   0.439

res3$li[1:2,]  #Row coordinates dims =c(4989,10)
```

```
##            Axis1      Axis2
## gene_1 0.6098867 -0.6836028
## gene_2 0.1179244 -0.2704358
```

```r
res3$co[1:2,]  #Col coordinates dims=c(150,10)
```

```
##             Comp1        Comp2
## cell_1 -0.2770690 -0.007619305
## cell_2 -0.3045123 -0.041767103
```

```r
res3$l1 [1:2,1:2]   #Row scores dims =c(4989,10)
```

```
##             RS1       RS2
## gene_1 2.714914 -3.110609
## gene_2 0.524941 -1.230568
```

```r
res3$c1 [1:2,1:2]   #Col scores dims=c(150,10)
```

```
##              CS1         CS2
## cell_1 -1.233374 -0.03467025
## cell_2 -1.355538 -0.19005352
```

```r
plot(res3$co[,1], res3$co[,2])
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

```r
summary(res3)
```

```
## Class: coa dudi
## Call: ade4::dudi.coa(df = mat, scannf = FALSE, nf = 2)
## 
## Total inertia: 3.364
## 
## Eigenvalues:
##     Ax1     Ax2     Ax3     Ax4     Ax5 
## 0.05046 0.04830 0.03039 0.03018 0.02958 
## 
## Projected inertia (%):
##     Ax1     Ax2     Ax3     Ax4     Ax5 
##  1.5003  1.4359  0.9035  0.8972  0.8794 
## 
## Cumulative projected inertia (%):
##     Ax1   Ax1:2   Ax1:3   Ax1:4   Ax1:5 
##   1.500   2.936   3.840   4.737   5.616 
## 
## (Only 5 dimensions (out of 149) are shown)
```

# Comparing the output of all 3 approaches


Correlation between glmPCA and corral. Note the PCs are flipped, PC1 in Corral shares r=0.99 correlation with PC2 of glmPCA(type = poi). 

```r
cor(res1$factors[,1], res2$PCv[,2])
```

```
## [1] 0.9912386
```

```r
cor(res1$factors[,2], res2$PCv[,1])
```

```
## [1] -0.9900288
```

```r
par(mfrow=c(1,2))
plot(-res1$factors[,1], res2$PCv[,2])
plot(res1$factors[,2], res2$PCv[,1])
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)


Correlation between correspondence analysis and corral is identical


```r
cor(res3$co[,1], res2$PCv[,1])
```

```
## [1] 1
```

```r
cor(res3$co[,2], res2$PCv[,2])
```

```
## [1] -1
```

```r
par(mfrow=c(1,2))
plot(res3$co[,1], res2$PCv[,1])
plot(res3$co[,2], res2$PCv[,2])
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)


So what does correspondence analysis do and how is it related to decomposotio of the Pearson Residuals. As mentioned above a critical step in decomposition is choosing the form of the data that makes more sense. SVD will detect linear vectors that can the most variance in the data.  So if the variance is skewed or varies with the mean, it needs normalization.  

I will describe it in two ways, you could consider the expected value of an element (Xij) in a matrix N to be the producet of its row (ri) and column (cj) weight, where the row weight is the rowsum/total sum and thus represents the contributes of that row to the total matrix.  Equally the column weight is the colsum/total sum, and is the contribution of that column to the total.  Sometimes the row and column weight are called the row and column mass respectively.  The expected weight of every element is the outer product of the row and column weights 


## expected matrix

```r
rw<-rowSums(mat)/sum(mat) # row weights
cw<-colSums(mat)/sum(mat)  # column weights

length(rw)
```

```
## [1] 4989
```

```r
length(cw)
```

```
## [1] 150
```

```r
expw<- outer(cw, rw)
dim(expw)
```

```
## [1]  150 4989
```

```r
expw[1:2,1:3]
```

```
##              gene_1       gene_2       gene_3
## cell_1 1.433571e-06 6.881143e-07 4.205143e-07
## cell_2 1.446955e-06 6.945386e-07 4.244403e-07
```

In COA, the data are treated like a contingency tables so the residuals are the difference between the observed data and the expected under the assumption that there is no relationship.The Pearson chi sq statistics is the observed-expected/ sqft (expected) and this matrix is subject to svd.

![COA]("./coa.png")   

A much more detailed explanation of correspondence analysis is provided at https://www.displayr.com/math-correspondence-analysis/


# Few Notes

One disadvantage of GLM-PCA is it uses on an iterative algorithm to obtain estimates for the latent factors and is at least ten times slower than PCA.  Becaause it is iterative, somes estimates are unstable, and thus in practice we run it multiple times to get a stable estimate.  In runs, it tends to flip the order of Prinipcal Components.

## glmPCA tend to flip PC order, relative to COA

```r
cor(res1$factors[,1], res2$PCv[,2])
```

```
## [1] 0.9912386
```

```r
cor(res1$factors[,2], res2$PCv[,1])
```

```
## [1] -0.9900288
```

```r
par(mfrow=c(1,2))
plot(-res1$factors[,1], res2$PCv[,2])
plot(res1$factors[,2], res2$PCv[,1])
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

## Variance of components in different runs of glmpca on the same dataset



```r
library(glmpca)
iter=20
tt<-lapply(1:iter, function(x) glmpca(mat, L=2, fam="poi",verbose=FALSE))
factors1<-sapply(1:iter, function(i) tt[[i]]$factors$dim1)  # first factor
factors2<-sapply(1:iter, function(i) tt[[i]]$factors$dim2)  # first factor

colnames(factors1) = paste0("PC1", 1:iter, "_")
colnames(factors2) = paste0("PC2", 1:iter, "_")

boxplot(factors1, ylab="factor PC1", xlab="iter")  # dim=c(150  , iter)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

```r
loadings1<-sapply(1:iter, function(i) tt[[i]]$loadings$dim1)  # first loadings
boxplot(loadings1,  ylab="loadings PC1", xlab="iter")# dim=c(4989  , iter)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-2.png)





# References

Abdi, H., & Valentin, D. (2007). Multiple Correspondence Analysis. In N. Salkind (Ed.), Encyclopedia of Measurement and Statistics (pp. 652–657). Sage Publications, Inc. https://doi.org/10.4135/9781412952644.n299

Benzécri, J.-P. (Ed.). (1973). L’analyse des données. 

Busold, C. H., Winter, S., Hauser, N., Bauer, A., Dippon, J., Hoheisel, J. D., & Fellenberg, K. (2005). Integration of GO annotations in Correspondence Analysis: Facilitating the interpretation of microarray data. Bioinformatics, 21(10), 2424–2429. https://doi.org/10.1093/bioinformatics/bti367

Culhane, A C, Perrière, G., & Higgins, D. G. (2003). Cross-platform comparison and visualisation of gene expression data using co-inertia analysis. BMC Bioinformatics, 4(1), 59. https://doi.org/10.1186/1471-2105- 4-59

Culhane, A. C., Perriere, G., Considine, E. C., Cotter, T. G., & Higgins, D. G. (2002). Between-group analysis of microarray data. Bioinformatics, 18(12), 1600–1608. https://doi.org/10.1093/bioinformatics/18.12.1600

Duò, A., Robinson, M. D., & Soneson, C. (2018). A systematic performance evaluation of clustering methods for single cell RNA-seq data. F1000Research, 7, 1141. https://doi.org/10.12688/f1000research.15666.2

Fellenberg, K., Hauser, N. C., Brors, B., Neutzner, A., Hoheisel, J. D., & Vingron, M. (2001). Correspondence analysis applied to microarray data. Proceedings of the National Academy of Sciences, 98(19), 10781–10786. https://doi.org/10.1073/pnas.181597298

Grantham, R., Gautier, C., & Gouy, M. (1980). Codon frequencies in 119 individual genes confirm corsistent choices of degenerate bases according to genome type. Nucleic Acids Research, 8(9), 1893–1912. https://doi.org/10.1093/nar/8.9.1893

Greenacre, M. (2013). The contributions of rare objects in correspondence analysis. Ecology, 94(1), 241–249. https://doi.org/10.1890/11-1730.1

Greenacre, M. J. (1984). Theory and applications of correspondence analysis. Academic Press.

Greenacre, M. J. (2010). Correspondence analysis: Correspondence analysis. Wiley Interdisciplinary Reviews: Computational Statistics, 2(5), 613–619. https://doi.org/10.1002/wics.114

Greenacre, M., & Hastie, T. (1987). The Geometric Interpretation of Correspondence Analysis. Journal of the American Statistical Association, 82(398), 437–447.https://doi.org/10.1080/01621459.1987.10478446

Hafemeister, C., & Satija, R. (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biology, 20(1), 296. https://doi.org/10.1186/s13059-019-1874-1

Hicks, S. C., Townes, F. W., Teng, M., & Irizarry, R. A. (2018). Missing data and technical variability in single-cell RNAsequencing experiments. Biostatistics, 19(4), 562–578. https://doi.org/10.1093/biostatistics/kxx053

Hirschfeld, H. O. (1935). A Connection between Correlation and Contingency. Mathematical Proceedings of the Cambridge Philosophical Society, 31(4), 520–524. https://doi.org/10.1017/S0305004100013517

Holmes, S. (2008). Multivariate data analysis: The French way. ArXiv:0805.2879 [Stat], 219–233. https://doi.org/10.1214/193940307000000455

Hsu L, Culhane A. Impact of Data Preprocessing on Integrative Matrix Factorization of Single Cell Data. Front Oncol. 2020;10:973. doi:10.3389/fonc.2020.00973.

Hubert, L., & Arabie, P. (1985). Comparing partitions. Journal of Classification, 2(1), 193–218. https://doi.org/10.1007/BF01908075

Legendre, P., Legendre, L., Legendre, L., & Legendre, L. (1998). Numerical ecology (2nd English ed). Elsevier.

Meng, C., Kuster, B., Culhane, A. C., & Gholami, A. (2014). A multivariate approach to the integration of multi-omics datasets. BMC Bioinformatics, 15(1), 162. https://doi.org/10.1186/1471-2105-15-162

Meng, C., Zeleznik, O. A., Thallinger, G. G., Kuster, B., Gholami, A. M., & Culhane, A. C. (2016). Dimension reduction techniques for the integrative analysis of multi-omics data. Briefings in Bioinformatics, 17(4), 628–641. https://doi.org/10.1093/bib/bbv108

Townes, F.W., Hicks, S.C., Aryee, M.J. et al. Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model. Genome Biol 20, 295 (2019). https://doi.org/10.1186/s13059-019-1861-6 

