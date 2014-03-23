#-------------------------------------------------------------------------------
merge = function(esets, method=NULL, norm = "FRMA", features = "GENE")
{
  
  if( !is.null(method) && !is.element(method, MERGING_METHODS))
  {
    err("Merging Method: ",method," is not supported.");
  }
  
  esets = formatData(esets, norm, features);
  if( is.null(method)) {
    msg("Run with no additional merging technique...");
    return(mergeNONE(esets));
  }
  #if(method=="NONE")	 { msg("Run NONE..."); 	   return(mergeNONE(esets));	 }
  if(method=="BMC")      { msg("Run BMC...");      return(mergeBMC(esets));      }
  if(method=="COMBAT")   { msg("Run COMBAT...");   return(mergeCOMBAT(esets));   }
  if(method=="DWD")      { msg("Run DWD...");      return(mergeDWD(esets));      }
  if(method=="GENENORM") { msg("Run GENENORM..."); return(mergeGENENORM(esets)); }
  if(method=="XPN")      { msg("Run XPN...");      return(mergeXPN(esets));      }
}
