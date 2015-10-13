# apply xpn on list of eSets
#-------------------------------------------------------------------------------
mergeXPN = function(esets)
{
  aux = function(x,y)
  {
    res = xpn(x,y);
    return(mergeNONE(list(res[[1]], res[[2]])));
  }

  return(combineByTwo(esets,aux,rand=TRUE));
}
