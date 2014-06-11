#-------------------------------------------------------------------------------
merge = function(esets, method=NA)
{
  
  if( !is.na(method) && !is.element(method, MERGING_METHODS))
  {
    err("Merging Method: ",method," is not supported.");
  }
  
  if(is.na(method)) {
    msg("Run with no additional merging technique...");
    return(mergeNONE(esets));
  }
  if(method=="BMC")      { msg("Run BMC...");      return(mergeBMC(esets));      }
  if(method=="COMBAT")   { msg("Run COMBAT...");   return(mergeCOMBAT(esets));   }
  if(method=="DWD")      { msg("Run DWD...");      return(mergeDWD(esets));      }
  if(method=="GENENORM") { msg("Run GENENORM..."); return(mergeGENENORM(esets)); }
  if(method=="XPN")      { msg("Run XPN...");      return(mergeXPN(esets));      }
}
