# apply geneNormESet on list of eSets
#-------------------------------------------------------------------------------
mergeGENENORM = function(esets)
{
  esets = lapply(esets, geneNormESet);
  eset = mergeNONE(esets);  
}

#-------------------------------------------------------------------------------
geneNormESet = function(eset)
{
  exprs(eset) = apply(exprs(eset), 1, function(x){x-mean(x, na.rm=TRUE)});
  exprs(eset) = apply(exprs(eset), 1, function(x){x/sd(x, na.rm=TRUE)});
  return(eset);
}
