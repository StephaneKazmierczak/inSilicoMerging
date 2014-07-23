#-------------------------------------------------------------------------------
merge = function(esets, method=NA)
{
	
	if(!is.na(method) && !is.element(method, MERGING_METHODS))
	{
		err("Merging Method: ",method," is not supported.");
	}
	
	if(is.na(method)) {
		msg("Run with no additional merging technique...");
		return(matchExprsPheno(mergeNONE(esets)));
	}
	if(method=="BMC")      { msg("Run BMC...");      return(matchExprsPheno(mergeBMC(esets)));      }
	if(method=="COMBAT")   { msg("Run COMBAT...");   return(matchExprsPheno(mergeCOMBAT(esets)));   }
	if(method=="DWD")      { msg("Run DWD...");      return(matchExprsPheno(mergeDWD(esets)));      }
	if(method=="GENENORM") { msg("Run GENENORM..."); return(matchExprsPheno(mergeGENENORM(esets))); }
	if(method=="XPN")      { msg("Run XPN...");      return(matchExprsPheno(mergeXPN(esets)));      }
	
}
