# apply geneBMCESet on list of eSets
#-------------------------------------------------------------------------------
mergeBMC = function(esets)
{
	esets = lapply(esets, geneBMCESet);
	
  eset = mergeNONE(esets);
  return(eset);
}

#-------------------------------------------------------------------------------
geneBMCESet = function(eset)
{
	# SM: Changed / in - since exprs already on log scale
	# This corresponds to shifting mean to 0
	exprs(eset) = exprs(eset) - rowMeans(exprs(eset));
	return(eset);
}
