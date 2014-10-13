#-------------------------------------------------------------------------------
merge = function(esets, method= "NONE")
{
	
	if(!is.element(method, MERGING_METHODS))
	{
		err("Merging Method: ",method," is not supported.");
	}
	
	
	if(method == "NONE") {
		msg("Run with no additional merging technique...");
		mergedEset = mergeNONE(esets);
	} 	else if(method=="BMC") 
	{ 
		msg("Run BMC...");
		mergedEset = mergeBMC(esets);      
	} else if(method=="COMBAT")   
	{ 
		msg("Run COMBAT...");
		mergedEset = mergeCOMBAT(esets);
	} else if(method=="GENENORM")
	{ 	
		msg("Run GENENORM..."); 
		mergedEset = mergeGENENORM(esets); 
	} else if(method=="XPN")
	{ 
		msg("Run XPN...");
		mergedEset = mergeXPN(esets);
	}
	
	return(matchExprsPheno(mergedEset));
	
}
