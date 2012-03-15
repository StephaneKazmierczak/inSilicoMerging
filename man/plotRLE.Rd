\name{plotRLE}
\alias{plotRLE}
\title{Create RLE plot from (merged) ExpressionSet}
\description{Create relative log expression (RLE) plot from ExpressionSet. RLE plots were initially proposed to measure the overall quality of a dataset but can also be used to visualize the presence of unwanted batch effects in the data.}

\usage{
plotRLE(eset, colLabel, legend=TRUE, file=NULL, ...)
}

\arguments{
\item{eset}{ExpressionSet object.}
\item{colLabel}{colname in pData(eset) to retrieve information for the labeling of samples with a color. All samples with the same value in pData(eset)[,colLabel] will share the same color.}
\item{legend}{If TRUE a legend will be provided next to the RLE plot.}
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

# Visual inspection of the two merged datasets through an MDS plot
plotRLE(eset_NONE, colLabel="Disease")
plotRLE(eset_COMBAT, colLabel="Disease")

}
