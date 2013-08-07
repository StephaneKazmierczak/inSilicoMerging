dwd = function(eset1, eset2)
{
  n1 = ncol(exprs(eset1));
  n2 = ncol(exprs(eset2));

  eset = mergeNONE(list(eset1, eset2));
  batchInfo = c(rep(1,ncol(eset1)),
                rep(2,ncol(eset2)));
  res = kdwd(x=t(exprs(eset)), y=as.factor(batchInfo), scaled=FALSE);
	
  w = res@w[[1]]

  eset1 = eset[,1:n1];
  eset2 = eset[,(n1+1):(n1+n2)];

  #  Project data
  vproj1 = t(exprs(eset1)) %*% w ;
  vproj2 = t(exprs(eset2)) %*% w ;
 
  meanproj1 = mean(vproj1) ;
  meanproj2 = mean(vproj2) ; 
 
  #  Subtract respective class means
  exprs(eset1) = apply(exprs(eset1), 2, function(x){x - meanproj1 * w});
  exprs(eset2) = apply(exprs(eset2), 2, function(x){x - meanproj2 * w});

  rownames(exprs(eset1)) = rownames(fData(eset1));
  rownames(exprs(eset2)) = rownames(fData(eset2));

  return(list(eset1,eset2));
}
