# Concatenate different esets. No data transformation is applied in this
# function.
#-------------------------------------------------------------------------------
mergeNONE = function(esets)
{
  eset1 = esets[[1]];
  annot1 = annotation(eset1)
  
  for(i in 2:length(esets))
  {
    eset2 = esets[[i]];
    d1 = exprs(eset1);
  	d2 = exprs(eset2);

    #-----------------------------------------------------
    # Rebuild fData
    #-----------------------------------------------------
    cg = sort(intersect(rownames(d1), rownames(d2))); 
    # If too few overlapping genes
    if(length(cg) < (min(dim(d1)[1],dim(d2)[1])/100))
    {
    	msg(" ! WARNING ! Number of common genes < 1%")
    }
    fData = fData(eset1)[cg,]
    
    #-----------------------------------------------------
    # Rebuild pData
    #-----------------------------------------------------
   
    # Create new pData: first common annotations, then unique ones
	  p1 = pData(eset1)
	  p2 = pData(eset2)
	  cp = sort(intersect(colnames(p1),colnames(p2)));
	  tp = sort(unique(union(colnames(p1),colnames(p2))))

	  sp1 = setdiff(colnames(p1),cp)
	  sp2 = setdiff(colnames(p2),cp)

	  pheno = matrix(NA,ncol=length(tp),nrow=nrow(p1)+nrow(p2))
	  rownames(pheno) = c(rownames(p1),rownames(p2))
	  colnames(pheno) = tp;
	
    if(length(cp)!=0)
	  {
		  pheno[1:nrow(p1),cp] = as.matrix(p1[,cp])
		  pheno[(nrow(p1)+1):(nrow(p1)+nrow(p2)),cp] = as.matrix(p2[,cp])
	  }
	
	  if(length(sp1)!=0)
	  {
		  pheno[1:nrow(p1),sp1] = as.matrix(p1[,sp1])
	  }
	  if(length(sp2)!=0)
	  {
		  pheno[(nrow(p1)+1):(nrow(p1)+nrow(p2)), sp2] = as.matrix(p2[,sp2])
	  }
	
	  pData = as.data.frame(pheno)

    #-----------------------------------------------------
    # Rebuild rest eset
    #-----------------------------------------------------
	
    # prevent to loose the column name in case there is only 1 sample in one exprs
    d1 = d1[cg, 1:ncol(d1), drop = FALSE];
    d2 = d2[cg, 1:ncol(d2), drop = FALSE];
    
    
    eset1 = new("ExpressionSet");    
    exprs(eset1) = cbind(d1, d2);
    pData(eset1) = pData;
    fData(eset1) = fData;

    annot1 = c(annot1, annotation(eset2));

  }
  
  #make sure the measurements matches
  measNames = sort(colnames(exprs(eset1)));  
  exprs(eset1) = exprs(eset1)[, measNames];
  pData(eset1) = pData(eset1)[measNames,];

  annotation(eset1) = unique(annot1);
  return(eset1);
}