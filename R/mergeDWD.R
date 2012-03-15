# apply dwd on list of eSets
#-------------------------------------------------------------------------------
mergeDWD = function(esets)
{
  aux = function(x,y)
  {
    res = suppressMessages(dwd(x,y));
    return(mergeNONE(list(res[[1]], res[[2]])));
  }

  return(combineByTwo(esets,aux,rand=TRUE));
}
