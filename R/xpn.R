##--------------------------------------------------------------------
# xpn(eset1,eset2)
# accum_array
##--------------------------------------------------------------------
# Implemenation is based on XPN implementation in Matlab accompanying 
# the paper: A.A. Shabalin et al., Merging two gene-expression studies
# via cross-platform normalization.Bioinformatics, Vol. 24, no. 9, 
# 2008, pages 1154-1160
#
# Code based on matlab code taken from: https://genome.unc.edu/xpn/
##--------------------------------------------------------------------

repmat = function(a,n,m) {kronecker(matrix(1,n,m),a)}

#---------------------------------------------------------------------
xpn=function(eset1,eset2,geneClusters=25,sampleClusters=5,nRep=32)
{
  #Find common genes
  res = commonGenes(list(eset1,eset2));
  eset1 = res[[1]];
  eset2 = res[[2]];

  x1 = exprs(eset1);
  x2 = exprs(eset2);

	#Size check
	p=nrow(x1)
	n1=ncol(x1)
	n2=ncol(x2)
	if(nrow(x1)!=nrow(x2))
	{
		err("Different number of genes!")
	}

	#Check for NaNs
	if(any(is.nan(c(x1))) || any(is.nan(c(x2))))
	{
		err("Missing values in input data!")
	}
	
	#Check for geneClusters argument
	if(nargs()<3)
	{
		geneClusters=min(25,max(floor(p/50),1))
	}
	
	#Check for sampleClusters argument
	if(nargs()<4)
	{
		sampleClusters=min(5,floor(min(n1,n2)/4))
	}
	
	if(geneClusters <= 0)
	{
		err("geneClusters is too small!")
	}
	
	if(sampleClusters <= 0)
	{
		err("sampleClusters is too small!")
	}
	
	if(nRep <= 0)
	{
		err("nRep is too small!")
	}
	
	if(nRep < 15)
	{
		war("nRep is small, should be at least 15!")
	}
	
	#Check for constant rows
	constrows<-sd(t(cbind(x1,x2)))==0
	if(any(constrows))
	{
		if(all(constrows))
		{
			war("Constant input data! Will be removed.")
		}

		x1fix=x1[as.numeric(!constrows)*1:n1,]
		x2fix=x2[as.numeric(!constrows)*1:n2,]
	
    exprs(eset1) = x1fix;
    exprs(eset2) = x2fix;
	
		XPN_result=xpn(eset1,
                   eset2,
                   geneClusters,
                   sampleClusters,
                   nRep)
		y1=XPN_result[[1]]
		y2=XPN_result[[2]]

		z1=matrix(0,p,n1)
		z1[as.numeric(!constrows)*1:n1,]=y1
		z1[as.numeric(constrows),]=x1[as.numeric(constrows),]

		z2=matrix(0,p,n2)
		z2[as.numeric(!constrows)*1:n2,]=y2
		z2[as.numeric(constrows),]=x2[as.numeric(constrows),]	
		return(list(z1,z2))
	}
	
  #Check skewness
	z=c(x1)
	m=mean(z)
	s=sd(z)
	skewness1=mean((z-m)^3)/s^3
	
	z=c(x2)
	m=mean(z)
	s=sd(z)
	skewness2=mean((z-m)^3)/s^3
	
	if(abs(skewness1)>8 || abs(skewness2)>8)
	{
		war("Input data is excessively skewed")
		war("The data may require log-transformation")
	}

  #Transform platforms to have 0 mean
  x1med=matrix(0,p,1)
  x2med=matrix(0,p,1)
  for(i in 1:p)
  {
    x1med[i]=median(x1[i,])
    x2med[i]=median(x2[i,])
  }

  x1primeMed=x1-repmat(x1med,1,n1)
  x2primeMed=x2-repmat(x2med,1,n2)
  xxprimeMed=cbind(x1primeMed, x2primeMed) 

  TheSum1=matrix(0,nrow(x1),ncol(x1))
  TheSum2=matrix(0,nrow(x2),ncol(x2))
	
	#Begin Cycle
  for(cycle in 1:nRep)
  {
			IDR=kmeans(xxprimeMed,geneClusters,1000,1)[[1]]
			repeatje=TRUE
      counter = 0;  # Count the number of troublesome clusters
      maxCounter = 200;
			while(repeatje)
			{
        if(counter > maxCounter) 
        { err("Unable to find clusters. Number of samples too small..."); }
				IDC=kmeans(t(xxprimeMed),sampleClusters,1000,1)[[1]]
				IDC1=IDC[1:n1];
				IDC2=IDC[(n1+1):(n1+n2)];
				countL1 = matrix(table(IDC1),1); 
        countL2 = matrix(table(IDC2),1);
				if (length(countL1) < sampleClusters || length(countL2) < sampleClusters)
				{ 
          #msg("Troublesome cluster:",counter,"/",maxCounter); 
          counter = counter + 1; 
          repeatje=TRUE; 
        }
				else
				{ 
          if((cycle && 8)==0) { msg("Iteration:",cycle,"/",nRep); } 
          repeatje=FALSE; 
        }
			}
	
			#GLP means and GP variance G -gene, L -colCluster, P- study
			sumGL1=accum_array(x1,IDC1,sampleClusters)
			sumGL2=accum_array(x2,IDC2,sampleClusters)
	
			sumGL1sq=accum_array(x1^2,IDC1,sampleClusters)
			sumGL2sq=accum_array(x2^2,IDC2,sampleClusters)
			
			xbar1=sumGL1 / repmat(countL1,p,1)
			xbar2=sumGL2 / repmat(countL2,p,1)
	
			sigmaGL1sq=sumGL1sq / repmat(countL1,p,1) - xbar1^2
			sigmaGL2sq=sumGL2sq / repmat(countL2,p,1) - xbar2^2
	
			#Apply MLE to each gene cluster
			A1=matrix(0,geneClusters,sampleClusters)
			A2=matrix(0,geneClusters,sampleClusters)
			b1=matrix(0,p,1)
			b2=matrix(0,p,1)
			c1=matrix(0,p,1)
			c2=matrix(0,p,1)
			s1=matrix(0,p,1)
			s2=matrix(0,p,1)
			
			#print("Start core function XPN_MLE")
			for(i in 1:geneClusters)
			{
				xpn_result=XPN_MLE(xbar1[IDR==i,],sigmaGL1sq[IDR==i,],countL1)
				A1[i,]=xpn_result[[1]]
				b1[IDR==i] = xpn_result[[2]]
				c1[IDR==i] = xpn_result[[3]]
				s1[IDR==i] = xpn_result[[4]]
				
				xpn_result=XPN_MLE(xbar2[IDR==i,],sigmaGL2sq[IDR==i,],countL2)
				A2[i,]=xpn_result[[1]]; 
				b2[IDR==i] = xpn_result[[2]]; 
				c2[IDR==i] = xpn_result[[3]];
				s2[IDR==i] = xpn_result[[4]]; 
	
			}
	
			#Average the estimates accross platforms
			A = (A1 * repmat(countL1,geneClusters,1) + A2 * repmat(countL2,geneClusters,1))/repmat(countL1+countL2,geneClusters,1)
			b = (b1 %*% n1 + b2 %*% n2)/(n1+n2)
			c = (c1 %*% n1 + c2 %*% n2)/(n1+n2)
			s = (s1 %*% n1 + s2 %*% n2)/(n1+n2)
		
			#KsiP residual of the model - P - study
			Ksi1 = (x1 - A1[IDR,IDC1] * repmat(b1,1,n1) - repmat(c1,1,n1))/repmat(sqrt(s1),1,n1)
			Ksi2 = (x2 - A2[IDR,IDC2] * repmat(b2,1,n2) - repmat(c2,1,n2))/repmat(sqrt(s2),1,n2)
	
			#xPstar transformed data - P - study
			x1star = A[IDR,IDC1] * repmat(b,1,n1) + repmat(c,1,n1) + Ksi1 * repmat(sqrt(s),1,n1)
			x2star = A[IDR,IDC2] * repmat(b,1,n2) + repmat(c,1,n2) + Ksi2 * repmat(sqrt(s),1,n2)
		
			#Check whether NaNs are created by XPN code
			if(any(any(is.nan(x1star))) || any(any(is.nan(x2star))))
			{
				err("NaN values are produced by the XPN code!")
			}

      #Aggregate results of several tries in TheSumP
      TheSum1=TheSum1+x1star
      TheSum2=TheSum2+x2star

  }#end cycle 1:nRep
				
	#Take the average over several tries
  z1=TheSum1/nRep
  z2=TheSum2/nRep
	
  exprs(eset1) = z1;
  exprs(eset2) = z2;
  return(list(eset1,eset2))
}

#Internal MLE procedure for XPN
XPN_MLE=function(xbar, sigma2,nj)
{
	repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}
	#xbar matrix of averages of X within column clusters (given gene, platform)
	#sigma2 - variance of the elements within column clusters (given gene,platform)
	#nj - number of elements in the column clusters
	
	#Remark: the platform/gene cluster is fixed
	
	#n = number of samples in the platform
	n=sum(nj)
	
	#I - number of genes in the gene cluster
	#J - number of sample clusters

  ##--> UPDATE 30Nov2009 JONATAN

  if(!is.matrix(xbar))
  {
    I = 1   
    J = length(xbar)
  }
  else
  {
    I=nrow(xbar)
    J=ncol(xbar)
  }

	#I=nrow(xbar)
	#J=ncol(xbar)
  
  ##--> UPDATE 30Nov2009 JONATAN
	
	#A,b,c,s2-parameters of the model
	A=matrix(0,1,J)
	b=matrix(1,I,1)
	c=matrix(0,I,1)
	s2=matrix(1,I,1)

	#previous values
	old=cbind(A,t(b),t(c),t(s2))*0
	current=cbind(A,t(b),t(c),t(s2))

	while(sum((current-old)^2)>1e-16 * max(old))
	{
		old=current

		#iteratively update the values
		c = rowSums((xbar - repmat(b,1,J)*repmat(A,I,1)) * repmat(nj,I,1))/n

		#fix sign of b
		if(sum(b)<0) { b = -b; }

		A = colSums(repmat(b,1,J) * (xbar - repmat(c,1,J)) / repmat(s2,1,J),1)/sum(b^2/s2)
		
    #Enforce constraints on A
		A=A-mean(A)
		A = A * sqrt(J/sum(A^2))

		b = rowSums(repmat(t(A),I,1) * (xbar - repmat(c,1,J)) * repmat(nj,I,1),2)/sum(A^2*nj)
		s2=rowSums(((xbar - repmat(c,1,J) - repmat(t(A),I,1) * repmat(b,1,J))^2 + sigma2)*repmat(nj,I,1))/n
  
    ##--> UPDATE 30Nov2009 JONATAN

		#s2(s2==0) = realmin('double')
		#s2[s2==0]=1e-16
    s2[s2<1e-16] = 1e-16
    
    ##--> UPDATE 30Nov2009 JONATAN
		
		A=t(A)
		current=cbind(A,t(b),t(c),t(s2))
	}

	#Return something
	return(list(A,b,c,s2))
}

#-------------------------------------------------------------------------------

accum_array = function(my_data, clustering, nr_clusters)
{
  accumulated=matrix(0,nrow(my_data),nr_clusters)
  for(j in 1:nr_clusters)
  {
    accumulated[,j]=rowSums(my_data[,clustering==j,drop=F])
  }
  return(accumulated)
}
