library(inSilicoDb);

#-------------------------------------------------------------------------------
pp  = function(...) { paste(...,sep=""); }
err = function(...) { stop(...,call.=FALSE); }
msg = function(...) { message(paste("  INSILICOMERGING:",...)); } 
war = function(...) { msg(" ! WARNING ! ",...); }

#-------------------------------------------------------------------------------
MERGING_METHODS = c("BMC", "COMBAT", "DWD", "NONE", "GENENORM", "XPN");

#-------------------------------------------------------------------------------
isOdd = function(x) { as.logical(x%%2) };

# Combines recursively all elements of a list pairwise in the following manner:
# [ A ; B ; C ; D ; E ]
#   => [ E ; AB ; CD ]
#   => [ CD ; EAB ]
#   => [ CDEAB ]
#-------------------------------------------------------------------------------
combineByTwo = function(x, fun, rand=FALSE, ...)
{
  if(rand) { x = sample(x); }
  
  n = length(x);
  if(n==1) { return(x[[1]]); }
  
  res = NULL;
  if(isOdd(n)) { res = c(res,x[[n]]); }
  for(i in seq(1,n-1,by=2))
  {
    res = c(res,fun(x[[i]],x[[i+1]],...));
  }
  combineByTwo(res,fun);
}

#-------------------------------------------------------------------------------
commonGenes = function(lst)
{
  common_genes = identify_common_genes(lst);
  lst = lapply(lst,function(x){return(x[rownames(exprs(x))%in%common_genes,])});
  return(lst);
}

#-------------------------------------------------------------------------------
identify_common_genes = function(lst)
{
  temp = rownames(exprs(lst[[1]]))
  for(eset in lst)
  {
    temp = intersect(rownames(exprs(eset)),temp);
  }
  return(temp);
}

#-------------------------------------------------------------------------------
matchExprsPheno = function(eset) {
	#make sure the samples are ordered in the same way
	sampleNames = colnames(exprs(eset));  
	pData(eset) = pData(eset)[sampleNames, , drop = FALSE];
	return(eset);
}

#-------------------------------------------------------------------------------

