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
eset1 = getDataset("GSE19804", "GPL570", norm="FRMA", genes=TRUE);
eset2 = getDataset("GSE10072", "GPL96",  norm="FRMA", genes=TRUE);
esets = list(eset1,eset2);

# merge them using the 'NONE' and 'COMBAT' method:
library(inSilicoMerging)
eset_NONE = merge(esets, method="NONE");
eset_COMBAT = merge(esets, method="COMBAT");

# check available annotations:
colnames(pData(eset_NONE))
table(pData(eset_NONE)[,"Disease"]);
table(pData(eset_NONE)[,"Study"]);

# Visual inspection of a random gene in the two merged datasets
gene = sample(rownames(exprs(eset_NONE)), 1)
plotGeneWiseBoxPlot(eset_NONE, colLabel="Disease", batchLabel="Study", gene=gene)
plotGeneWiseBoxPlot(eset_COMBAT, colLabel="Disease", batchLabel="Study", gene=gene)

}
