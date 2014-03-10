#-------------------------------------------------------------------------------
merge = function(esets, method=NULL, norm = "FRMA")
{
  
  if( !is.null(method) && !is.element(method, MERGING_METHODS))
  {
    err("Merging Method: ",method," is not supported.");
  }
  
  if( is.null(method) && (norm == "FRMA" || norm == "UPC" || norm == "SCAN")) {
    msg("Run with no additional merging technique...");
    esets = formatData(esets, method);
    return(mergeNONE(esets));
  }
  
  esets = formatData(esets, norm);
  #if(method=="NONE")	 { msg("Run NONE..."); 	   return(mergeNONE(esets));	 }
  if(method=="BMC")      { msg("Run BMC...");      return(mergeBMC(esets));      }
  if(method=="COMBAT")   { msg("Run COMBAT...");   return(mergeCOMBAT(esets));   }
  if(method=="DWD")      { msg("Run DWD...");      return(mergeDWD(esets));      }
  if(method=="GENENORM") { msg("Run GENENORM..."); return(mergeGENENORM(esets)); }
  if(method=="XPN")      { msg("Run XPN...");      return(mergeXPN(esets));      }
}
