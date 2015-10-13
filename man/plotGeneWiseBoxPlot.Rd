\name{plotGeneWiseBoxPlot}
\alias{plotGeneWiseBoxPlot}
\title{Create gene-wise boxplot from (merged) ExpressionSet}
\description{Gene-wise boxplots describe the gene-wise distribution of samples. Sample can be grouped together using the batchLabel parameter and can be coloroed usig the colLabel parameter for optimal visualization of the possible batch effects.}

\usage{
plotGeneWiseBoxPlot(eset, colLabel, batchLabel, gene=NULL, legend=TRUE, file=NULL, ...)
}

\arguments{
\item{eset}{ExpressionSet object.}
\item{colLabel}{colname in pData(eset) to retrieve information for the labeling of samples with a color. All samples with the same value in pData(eset)[,colLabel] will share the same color.}
\item{batchLabel}{colname in pData(eset) to retrieve information for the grouping of samples. All samples with the same value in pData(eset)[,batchLabel] and with the same color will be grouped together.}
\item{gene}{Gene for which the boxplot will be created. If not specified a random gene will be selected.}
\item{legend}{If TRUE a legend will be provided next to the gene-wise box plot.}
\item{file}{If defined, the resulting plot will be stored as a pdf file instead of shown interactively.}
\item{...}{Additional parameters for the 'plot' function (e.g. 'main').}
}

\examples{

# retrieve two datasets:
library(inSilicoDb);
InSilicoLogin("rpackage_tester@insilicodb.com", "5c4d0b231e5cba4a0bc54783b385cc9a");
eset1 = getDataset("GSE18842", "GPL570", norm="FRMA", features="gene");
eset2 = getDataset("GSE31547", "GPL96",  norm="FRMA",  features="gene");
esets = list(eset1,eset2);

# merge them using no additional merging technique and the 'COMBAT' method:
library(inSilicoMerging)
eset_FRMA = merge(esets);
eset_COMBAT = merge(esets, method="COMBAT");

# check available annotations:
colnames(pData(eset_FRMA))
table(pData(eset_FRMA)[,"Disease"]);
table(pData(eset_FRMA)[,"Study"]);

# Visual inspection of a random gene in the two merged datasets
gene = sample(rownames(exprs(eset_FRMA)), 1)
plotGeneWiseBoxPlot(eset_FRMA, colLabel="Disease", batchLabel="Study", gene=gene)
plotGeneWiseBoxPlot(eset_COMBAT, colLabel="Disease", batchLabel="Study", gene=gene)

}
