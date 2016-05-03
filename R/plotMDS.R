# plot colored MDS plot
#-------------------------------------------------------------------------------
plotMDS = function(eset, colLabel, symLabel, legend=TRUE, file=NULL, ...)
{
  if(! is.null(file)) { pdf(file, width=12, height=7); }
  
  mds = cmdscale(dist(t(exprs(eset))), eig=TRUE);
  
  colMap = makeColorMap(eset, colLabel); 
  colVec = makeColorVec(eset, colLabel, colMap);

  #--Add margin to the right for the legend 
  tmp = par()$mar;
  if(legend) { par(xpd=T, mar=par()$mar+c(0,0,0,8)); }

  range_x = range(mds$points[,1]);
  range_y = range(mds$points[,2]);

  plot(mds$points,
       col=colVec,
       pch=as.numeric(pData(eset)[,symLabel]),
       panel.first={ U = par("usr");
                     rect(U[1],U[3],U[2],U[4],
                          col="azure2",
                          border="black",
                          lwd=3)},
       lwd=2,
       xlab="",
       ylab="",
       xlim=range_x,
       ylim=range_y,
       ...);

  if(legend)
  {
    x = range_x[2] + (range_x[2]-range_x[1])*0.1;
    y = range_y[2] - (range_y[2]-range_y[1])*0.1;
  
    syms = unique(pData(eset)[,symLabel])
    legend(x,y,
           legend = syms,
           pt.lwd=2,
           pch = as.numeric(syms),
           box.lwd=3,
           bg="azure2");

    legend(x,y-(length(syms)*(range_y[2]-range_y[1])*0.1),
           legend = names(colMap),
           pt.lwd=2,
           pch=19,
           col = unlist(colMap),
           box.lwd=3,
           bg="azure2");
  
    #--Reset margin
    par(xpd=F, mar=tmp)
  }
  
  if(! is.null(file)) { dev.off(); }
}
