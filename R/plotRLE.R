# plot colored RLE plot
#-------------------------------------------------------------------------------
plotRLE = function(eset, colLabel, legend=TRUE, file=NULL, ...)
{
  if(! is.null(file)) { pdf(file, width=12, height=7); }
  
  colMap = makeColorMap(eset, colLabel); 
  colVec = makeColorVec(eset, colLabel, colMap);

  deviations = exprs(eset) - rowMedians(exprs(eset));

  #--Add margin to the right for the legend 
  tmp = par()$mar;
  if(legend) { par(xpd=T, mar=par()$mar+c(0,0,0,8)); }
 
  min_x = 1;
  max_x = ncol(eset);
  min_y = min(deviations);
  max_y = max(deviations);
  #min_y = -2;
  #max_y = 2;

  #--Create dummy empty plot and set axis
  plot(1,
       xlab="",
       ylab="",
       xaxt="n",
       xlim=c(1,ncol(eset)),
       ylim=c(min_y,max_y));

  #--Add layout
  U = par("usr");
  lines(c(U[1],U[2]),c(0,0),
        lwd=1,
        panel.first={rect(U[1],U[3],U[2],U[4],
                          col="azure2",
                          border="black",
                          lwd=3);});

  #--Add boxplots
  boxplot.matrix(deviations, 
                 col=colVec,
                 outline=FALSE, 
                 add = TRUE,
                 #names=rep("",ncol(eset)),
                 #range=1.5,
                 range=0,
                 xaxt="n",
                 ...);

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
