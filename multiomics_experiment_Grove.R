#Write RScript_Multiomics.R that takes data x1, x2, and x3 by copying and pasting the code snippets from 
#https://compgenomr.github.io/book/use-case-multi-omics-data-from-colorectal-cancer.html and does three ways of JDR based on the code snippets from 
#https://compgenomr.github.io/book/matrix-factorization-methods-for-unsupervised-multi-omics-data-integration.html.

#copied from the github
devtools::install_github("compgenomr/compGenomRData")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c('qvalue','plot3D','ggplot2','pheatmap','cowplot',
                       'cluster', 'NbClust', 'fastICA', 'NMF','matrixStats',
                       'Rtsne', 'mosaic', 'knitr', 'genomation',
                       'ggbio', 'Gviz', 'DESeq2', 'RUVSeq',
                       'gProfileR', 'ggfortify', 'corrplot',
                       'gage', 'EDASeq', 'citr', 'formatR',
                       'svglite', 'Rqc', 'ShortRead', 'QuasR',
                       'methylKit','FactoMineR', 'iClusterPlus',
                       'enrichR','caret','xgboost','glmnet',
                       'DALEX','kernlab','pROC','nnet','RANN',
                       'ranger','GenomeInfoDb', 'GenomicRanges',
                       'GenomicAlignments', 'ComplexHeatmap', 'circlize', 
                       'rtracklayer', 'BSgenome.Hsapiens.UCSC.hg38',
                       'BSgenome.Hsapiens.UCSC.hg19','tidyr',
                       'AnnotationHub', 'GenomicFeatures', 'normr',
                       'MotifDb', 'TFBSTools', 'rGADEM', 'JASPAR2018'
))
# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_gex.csv", 
                       package="compGenomRData")
x1 <- read.csv(csvfile, row.names=1)
# Fix the gene names in the data frame
rownames(x1) <- sapply(strsplit(rownames(x1), "\\|"), function(x) x[1])
# Output a table
knitr::kable(head(t(head(x1))), caption="Example gene expression data (head)")

# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_muts.csv", 
                       package="compGenomRData")
x2 <- read.csv(csvfile, row.names=1)
# Set mutation data to be binary (so if a gene has more than 1 mutation,
# we only count one)
x2[x2>0]=1
# output a table
knitr::kable(head(t(head(x2))), caption="Example mutation data (head)")

# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_cnv.csv", 
                       package="compGenomRData")
x3 <- read.csv(csvfile, row.names=1)
# output a table
knitr::kable(head(t(head(x3))), 
             caption="Example copy number data for CRC samples")

# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_subtypes.csv",
                       package="compGenomRData")
covariates <- read.csv(csvfile, row.names=1)
# Fix the TCGA identifiers so they match up with the omics data
rownames(covariates) <- gsub(pattern = '-', replacement = '\\.',
                             rownames(covariates))
covariates <- covariates[colnames(x1),]
# create a dataframe which will be used to annotate later graphs
anno_col <- data.frame(cms=as.factor(covariates$cms_label))
rownames(anno_col) <- rownames(covariates)

pheatmap::pheatmap(x1,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Gene expression data")

pheatmap::pheatmap(x2,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Mutation data")

pheatmap::pheatmap(x3,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Copy number data")

#now we move on to the next github link to run the JDR
# run the MFA function from the FactoMineR package
r.mfa <- FactoMineR::MFA(
    t(rbind(x1,x2,x3)), # binding the omics types together
    c(dim(x1)[1], dim(x2)[1], dim(x3)[1]), # specifying the dimensions of each
    graph=FALSE)

# first, extract the H and W matrices from the MFA run result
mfa.h <- r.mfa$global.pca$ind$coord
mfa.w <- r.mfa$quanti.var$coord

# create a dataframe with the H matrix and the CMS label
mfa_df <- as.data.frame(mfa.h)
mfa_df$subtype <- factor(covariates[rownames(mfa_df),]$cms_label)

# create the plot
ggplot2::ggplot(mfa_df, ggplot2::aes(x=Dim.1, y=Dim.2, color=subtype)) +
    ggplot2::geom_point() + ggplot2::ggtitle("Scatter plot of MFA")

pheatmap::pheatmap(t(mfa.h)[1:2,], annotation_col = anno_col,
                   show_colnames = FALSE,
                   main="MFA for multi-omics integration")

# Feature-normalize the data
x1.featnorm <- x1 / rowSums(x1)
x2.featnorm <- x2 / rowSums(x2)
x3.featnorm <- x3 / rowSums(x3)

# Normalize by each omics type's frobenius norm
x1.featnorm.frobnorm <- x1.featnorm / norm(as.matrix(x1.featnorm), type="F")
x2.featnorm.frobnorm <- x2.featnorm / norm(as.matrix(x2.featnorm), type="F")
x3.featnorm.frobnorm <- x3.featnorm / norm(as.matrix(x3.featnorm), type="F")

#excise 11.6.1 MF methods
#In order to use NMF algorithms on data which can be negative, 
#we need to split each feature into two new features, one positive and one negative. 
#Implement the following function, and see that the included test does not fail: 
#[Difficulty: Intermediate/Advanced]
#gotten from HJ
split_neg_columns <- function(x,
                              pos_suffix = "_pos",
                              neg_suffix = "_neg",
                              drop_zero  = FALSE) {
    x <- as.matrix(x)
    if (!is.numeric(x)) stop("x must be numeric.")
    if (anyNA(x)) warning("x contains NA values; they are propagated to the output.")
    
    # Column names (create defaults if missing)
    cn <- colnames(x)
    if (is.null(cn)) cn <- paste0("V", seq_len(ncol(x)))
    
    pos <- pmax(x, 0)
    neg <- pmax(-x, 0)
    
    colnames(pos) <- paste0(cn, pos_suffix)
    colnames(neg) <- paste0(cn, neg_suffix)
    
    out <- cbind(pos, neg)
    
    # Interleave so each feature's pos/neg columns sit side by side
    ord <- as.vector(rbind(seq_len(ncol(x)), ncol(x) + seq_len(ncol(x))))
    out <- out[, ord, drop = FALSE]
    
    if (drop_zero) {
        keep <- colSums(out != 0, na.rm = TRUE) > 0
        out <- out[, keep, drop = FALSE]
    }
    
    rownames(out) <- rownames(x)
    out
}

# Split the features of the CNV matrix into two non-negative features each
x3.featnorm.frobnorm.nonneg <- t(split_neg_columns(t(x3.featnorm.frobnorm)))


# run the nmf function from the NMF package
#require(NMF)

r.nmf <- nmf(t(rbind(x1.featnorm.frobnorm,
                     x2.featnorm.frobnorm,
                     x3.featnorm.frobnorm.nonneg)),
             2,
             method='Frobenius')

# exctract the H and W matrices from the nmf run result
nmf.h <- NMF::basis(r.nmf)
nmf.w <- NMF::coef(r.nmf)
nmfw <- t(nmf.w)

# create a dataframe with the H matrix and the CMS label (subtype)
nmf_df <- as.data.frame(nmf.h)
colnames(nmf_df) <- c("dim1", "dim2")
nmf_df$subtype <- factor(covariates[rownames(nmf_df),]$cms_label)

# create the scatter plot
ggplot2::ggplot(nmf_df, ggplot2::aes(x=dim1, y=dim2, color=subtype)) +
    ggplot2::geom_point() +
    ggplot2::ggtitle("Scatter plot of 2-component NMF for multi-omics integration")


pheatmap::pheatmap(t(nmf_df[,1:2]),
                   annotation_col = anno_col,
                   show_colnames=FALSE,
                   main="Heatmap of 2-component NMF")

# run the iClusterPlus function
r.icluster <- iClusterPlus::iClusterPlus(
    t(x1), # Providing each omics type
    t(x2),
    t(x3),
    type=c("gaussian", "binomial", "multinomial"), # Providing the distributions
    K=2, # provide the number of factors to learn
    alpha=c(1,1,1), # as well as other model parameters
    lambda=c(.03,.03,.03))

# extract the H and W matrices from the run result
# here, we refer to H as z, to keep with iCluster terminology
icluster.z <- r.icluster$meanZ
rownames(icluster.z) <- rownames(covariates) # fix the row names
icluster.ws <- r.icluster$beta

# construct a dataframe with the H matrix (z) and the cancer subtypes
# for later plotting
icp_df <- as.data.frame(icluster.z)
colnames(icp_df) <- c("dim1", "dim2")
rownames(icp_df) <- colnames(x1)
icp_df$subtype <- factor(covariates[rownames(icp_df),]$cms_label)

# create the scatter plot
ggplot2::ggplot(icp_df, ggplot2::aes(x=dim1, y=dim2, color=subtype)) +
    ggplot2::geom_point() +
    ggplot2::ggtitle("Scatter plot of 2-component Icluster for multi-omics integration")


pheatmap::pheatmap(t(icp_df[,1:2]),
                   annotation_col = anno_col,
                   show_colnames=FALSE,
                   main="Heatmap of 2-component i cluster")


#calculate the euclidean distance between 1 and 3
#https://www.rdocumentation.org/packages/rsgeo/versions/0.1.7/topics/distance_euclidean_pairwise
#install.packages(
#    "rsgeo",
#    repos = c("https://josiahparry.r-universe.dev", "https://cloud.r-project.org")
#)
#library(rsgeo)

#distance_euclidean_matrix(icp_df$dim1, icp_df$dim2)

#https://r-spatial.github.io/sf/reference/geos_measures.html
#library(sf)#

#st_as_sf(icp_df, arg = subtypes, sf_column_name = )

#st_distance()


#https://stackoverflow.com/questions/12102899/quickest-way-to-calculate-the-euclidean-distance-between-two-points-in-n-dimensi
#library(ThreeWay)
#nrm2(icp_df)

library(Coreset)
library(vegan)
library(labdsv)

#test <- pivot_wider(icp_df, names_from = subtype, values_from = subtype)
#get the vegist
#vegdist(icp_df, method="euclidean") 
#meandist(icp_df)

#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/dist
#https://www.rdocumentation.org/packages/Coreset/versions/1.0.0/topics/MeanDist
#back to trying MeanDist
nmf_dim <- nmf_df %>% select("dim1", "dim2")
nmf_dist <- dist(nmf_dim, method = "euclidean", diag = FALSE, upper = FALSE)
nmf_mean <- MaxMean(nmf_dist)
MeanDist(nmf_dist, nmf_mean)

#back to trying MeanDist
mfa_dim <- mfa_df %>% select("Dim.1", "Dim.2")
mfa_dist <- dist(mfa_dim, method = "euclidean", diag = FALSE, upper = FALSE)
mfa_mean <- MaxMean(mfa_dist)
MeanDist(mfa_dist, mfa_mean)


#back to trying MeanDist
icp_dim <- icp_df %>% select("dim1", "dim2")
icp_dist <- dist(icp_dim, method = "euclidean", diag = FALSE, upper = FALSE)
icp_mean <- MaxMean(icp_dist)
MeanDist(icp_dist, icp_mean)


#running sources
#https://vegandevs.github.io/vegan/reference/vegdist.html
#https://cran.r-project.org/web/packages/vegan/vegan.pdf
#https://www.rdocumentation.org/packages/vegan/versions/2.7-6/topics/vegdist
#https://www.datacamp.com/datalab/preview?_tag=rdocs&rdocsPath=packages%2FThreeWay%2Fversions%2F1.0%2Ftopics%2Fnrm2&utm_source=r-docs&utm_medium=docs&utm_term=nrm2&utm_content=run_example_in_datalab
#https://www.rdocumentation.org/packages/ThreeWay/versions/1.1.4/topics/nrm2
#https://stackoverflow.com/questions/12102899/quickest-way-to-calculate-the-euclidean-distance-between-two-points-in-n-dimensi
#https://stackoverflow.com/questions/12102899/quickest-way-to-calculate-the-euclidean-distance-between-two-points-in-n-dimensi
#https://www.rdocumentation.org/packages/sf/versions/1.1-3/topics/sf
#https://github.com/r-spatial/sf/blob/main/R/geom-measures.R
#https://r-spatial.github.io/sf/reference/geos_measures.html
#https://stackoverflow.com/questions/71774492/how-to-calculate-euclidean-distance-from-sf-data-frame-points-in-r
#https://www.rdocumentation.org/packages/rsgeo/versions/0.1.7
#https://www.rdocumentation.org/packages/rsgeo/versions/0.1.7/topics/distance_euclidean_pairwise
#https://rdrr.io/cran/vegan/src/R/meandist.R
#https://www.rdocumentation.org/packages/Coreset/versions/1.0.0/topics/MeanDist








