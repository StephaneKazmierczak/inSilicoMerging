
library (inSilicoDb)

# a very simple-minded test:  are the means of the null merge expression set different from the bmc merge?
test.compare_null_bmc_merge = function ()
{
  eset1 = getDataset("GSE19804", "GPL570", norm="FRMA", features = "gene"); 
  eset2 = getDataset("GSE10072", "GPL96", norm="FRMA", features = "gene"); 
  esets = list(eset1,eset2);
  eset_merged_null = merge(esets);

     ## disabled, 13 mar 2013, pshannon:
     ##    checkEquals (as.integer (nrow (eset_merged_null)), 12718);

  checkEquals (as.integer (nrow (eset_merged_null)), 12701); #12701
  checkEquals (as.integer (ncol (eset_merged_null)), 227);

  expectedColumnCount = as.integer (ncol (eset1)) + as.integer (ncol (eset2))
  actualColumnCount = as.integer (ncol (eset_merged_null))
  checkEquals (expectedColumnCount, actualColumnCount)

    # now merge using the BMC method
  eset_merged_bmc = merge (esets, method='BMC')
  
    # quick (dumb) sanity check on distribution of the values

  checkTrue (mean (exprs (eset_merged_null)) > 7)
  checkEqualsNumeric (mean (exprs (eset_merged_bmc)), 0)

} # test.compare_null_bmc_merge 


