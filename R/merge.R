#-------------------------------------------------------------------------------
merge = function(esets, method="NONE")
{
  if(! is.element(method, MERGING_METHODS))
  {
    err("Merging Method: ",method," is not supported.");
  }

  if(method=="BMC")      { msg("Run BMC...");      return(mergeBMC(esets));      }
  if(method=="COMBAT")   { msg("Run COMBAT...");   return(mergeCOMBAT(esets));   }
  if(method=="DWD")      { msg("Run DWD...");      return(mergeDWD(esets));      }
  if(method=="GENENORM") { msg("Run GENENORM..."); return(mergeGENENORM(esets)); }
  if(method=="NONE")     { msg("Run NONE...");     return(mergeNONE(esets));     }
  if(method=="XPN")      { msg("Run XPN...");      return(mergeXPN(esets));      }
}
