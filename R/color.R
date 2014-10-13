COLORS = c("red","green3","blue","cyan","magenta","yellow","gray","black",
           "orange","darkred","green","darkblue","darkcyan","darkmagenta",
           "darkorchid1","darkgoldenrod3","aquamarine","antiquewhite",
           "darkolivegreen3");

# create color map
#-------------------------------------------------------------------------------
makeColorMap = function(eset, label)
{
  colMap = list();
  vec = unique(as.vector(pData(eset)[,label]));
  for(i in 1:length(vec)) { colMap[[ vec[i] ]] = COLORS[i]; }
  return(colMap);
}  

# create color vector for plots
#-------------------------------------------------------------------------------
makeColorVec = function(eset, label, colMap)
{
  labels = as.vector(pData(eset)[,label]);
  return(as.vector(unlist(sapply(labels, function(x) { colMap[x]; }))));
}

