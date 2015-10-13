# plot gene-wise box
#-------------------------------------------------------------------------------
plotGeneWiseBoxPlot = function(eset, colLabel, batchLabel, gene=NULL, 
                               legend=TRUE, file=NULL, ...)
{
  if(! is.null(file)) { pdf(file, width=12, height=7); }

  if(is.null(gene)) 
  { 
    s = sample(1:nrow(exprs(eset)), 1); 
    gene = rownames(exprs(eset))[s];
  }
  
  colMap = makeColorMap(eset, colLabel); 

  #--Add margin to the right for the legend 
  tmp = par()$mar;
  if(legend) { par(xpd=T, mar=par()$mar+c(0,0,0,8)); }
  
  #--Get values per class and batch for boxplot
  values = list();
  names = NULL;
  colVec = NULL;

  batches = sort(unique(pData(eset)[,batchLabel]));
  for(batch in batches)
  {  
    idx1 = pData(eset)[,batchLabel]==batch;
    classes = sort(unique(pData(eset)[idx1,colLabel]));
    for(class in classes)
    {
      idx2 = pData(eset)[,colLabel]==class;
      idx = idx1 & idx2;
      values[[length(values)+1]] = as.vector(exprs(eset)[gene,idx]);
      names = c(names,batch);
      colVec = c(colVec,colMap[[class]]);
    }
    # Mimick empty bar between batches
    values[[length(values)+1]] = NA
    names = c(names,"");
    colVec = c(colVec,"white");
  } 
 
  #--Calculate plot boundaries 
  min_x = 0;
  max_x = length(values);
  min_y = min(unlist(values), na.rm=TRUE);
  max_y = max(unlist(values), na.rm=TRUE);

  #--Create dummy empty plot and set axis
  plot(1,
       xlab="",
       ylab="",
       xaxt="n",
       yaxt="n",
       xlim=c(min_x,max_x),
       ylim=c(min_y,max_y));

  #--Add boxplots
  boxplot(values, 
          col=colVec,
          outline=FALSE, 
          add = TRUE,
          range=0,
          xaxt="n",
          panel.first={U = par("usr"); rect(U[1],U[3],U[2],U[4],
                       col="azure2",
                       border="black",
                       lwd=3);},
          ...,
          main=gene);

  #--Reset X-axis
  axis(1,at=min_x+1:max_x,labels=names, tck=0, las=2); 

  #--Add legend
  if(legend)
  {
    x = max_x + (max_x - min_x) * 0.1;
    y = max_y - (max_y - min_y) * 0.1; 

    legend(x,y,
           legend = names(colMap),
           pt.lwd=2,
           pch = 19,
           col = unlist(colMap),
           box.lwd=3,
           bg="azure2");

    #--Reset margin
    par(xpd=F, mar=tmp);
  }
  
  if(! is.null(file)) { dev.off(); }
}
