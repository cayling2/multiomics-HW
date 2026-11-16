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
split_neg_columns <- function(x) {
    #input x is a dataframe
    #
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


