\name{merge}
\alias{merge}
\title{General method to merge different ExpressionSets}
\description{General method to merge different ExpressionSets by applying different techniques to remove inter-study bias.}

\usage{

merge(esets, method=NA);

}

\arguments{
\item{esets}{List of ExpressionSet objects.}
\item{method}{Merging method aimed at removing inter-study bias. Possible options are: BMC, COMBAT, DWD, GENENORM and XPN. If none are specified, the merging  More information about each method is given below in the details.}
}

\value{
  A (merged) ExpressionSet object.
}

\details{

Currently the following different merging techniques are provided:
\describe{
  \item{'BMC':}{In [1] they successfully applied a technique similar to z-score normalization for merging breast cancer datasets. They transformed the data by batch mean-centering, which means that the mean is subtracted.}
  \item{'COMBAT':}{Empirical Bayes [2] (also called EJLR or COMBAT) is a method that estimates the parameters of a model for mean and variance for each gene and then adjusts the genes in each batch to meet the assumed model. The parameters are estimated by pooling information from multiple genes in each batch.}
  \item{'DWD':}{In [3] they propose to use Distance Weighted Discrimination to find the hyperplane that separates the expression values of two studies. Assuming that the variation due to studies originating from different labs is bigger than any biological variation present in the data a translation in the direction of the normal vector of the hyperplane is used to remove bias.}
  \item{'GENENORM':}{One of the simplest mathematical transformations to make datasets more comparable is z-score normalization. In this method, for each gene expression value in each study separately all values are altered by subtracting the mean of the gene in that dataset divided by its standard deviation.}
  \item{NA:}{Combine esets without any additional transformation. Similar to 'combine' function.}
  \item{'XPN':}{The basic idea behind the cross-platform normalization [4] approach is to find blocks (clusters) of genes and samples in both studies that have similar expression characteristics. In XPN, a gene measurement can be considered as a scaled and shifted block mean.}
}

Note that after using any of those methods the resulting merged dataset only contains the common list of genes/probes between all studies.
}

\examples{

# retrieve two datasets:
library(inSilicoDb);
eset1 = getDataset("GSE19804", "GPL570", norm="FRMA", features="GENE");
eset2 = getDataset("GSE10072", "GPL96",  norm="FRMA", features="GENE");
esets = list(eset1,eset2);

# merge them using different methods:
library(inSilicoMerging);
eset_FRMA = merge(esets);
eset_BMC = merge(esets, method="BMC");

}

\references{
[1] A. Sims, \emph{et al.}, The removal of multiplicative, systematic bias allows integration of breast cancer gene expression datasets - improving meta-analysis and prediction of prognosis, \emph{BMC Medical Genomics}, vol. 1, no. 1, p. 42, 2008.

[2] C. Li and A. Rabinovic, Adjusting batch effects in microarray expression data using empirical bayes methods, \emph{Biostatistics}, vol. 8, no. 1, pp. 118-127, 2007.

[3] M. Benito, \emph{et al.}, Adjustment of systematic microarray data biases, \emph{Bioinformatics}, vol. 20, no. 1, pp. 105-114, 2004.

[4] A. A. Shabalin, \emph{et al.}, Merging two gene-expression studies via cross-platform normalization, \emph{Bioinformatics}, vol. 24, no. 9, pp. 1154-1160, 2008.
}
