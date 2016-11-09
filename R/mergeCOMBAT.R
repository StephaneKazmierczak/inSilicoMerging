# Apply Combat on list of eSets
#
# Code taken from: http://jlab.byu.edu/ComBat/Download_files/ComBat.R
# (slightly adapted to work with ExpressionSets)
#-------------------------------------------------------------------------------
mergeCOMBAT = function(esets)
{
	raw_merged = mergeNONE(esets)

  batchInfo = NULL;
  for(i in 1:length(esets))
  {
    batchInfo = c(batchInfo, rep(i,ncol(esets[[i]])));
  }
	
	saminfo = cbind(rownames(pData(raw_merged)),
                  rownames(pData(raw_merged)),
                  batchInfo)
	colnames(saminfo) = c("Array name", "Sample name", "Batch")
	
	dat = exprs(raw_merged)
	
	design <- design.mat(saminfo)
	
	batches <- list.batch(saminfo)
	n.batch <- length(batches)
	n.batches <- sapply(batches, length)
	n.array <- sum(n.batches)

	##Standardize Data across genes
	B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))
	grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
	var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
	
	stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
	if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}	
	s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))
	
	##Get regression batch effect parameters
	batch.design <- design[,1:n.batch]
	gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))
	delta.hat <- NULL
	for (i in batches)
	{
		delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
	}
	
	##Find Priors
	gamma.bar <- apply(gamma.hat, 1, mean)
	t2 <- apply(gamma.hat, 1, var)
	a.prior <- apply(delta.hat, 1, aprior)
	b.prior <- apply(delta.hat, 1, bprior)
	
	##Find EB batch adjustments
	gamma.star <- delta.star <- NULL
	for (i in 1:n.batch)
	{
			temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
	}
	
	### Normalize the Data ###
	bayesdata <- s.data
	j <- 1
	for (i in batches)
	{
		bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
		j <- j+1
	}
	
	bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
  
	print("Debug Combat")
	print(paste("dim bayesdata :",dim(bayesdata)))
	print(paste("dim exprs(eset) :",dim(exprs(raw_merged))))
	
	print(paste("colnames bayesdata :",colnames(bayesdata)))
	print(paste("colnames exprs(eset) :",colnames(exprs(raw_merged))))
	
	eset=raw_merged
	exprs(eset)=bayesdata
	return(eset)	
}

#Helper functions
#-------------------------------------------------------------------------------
build.design <- function(vec, des=NULL, start=2)
{
	tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
	for (i in 1:ncol(tmp)){tmp[,i] <- vec==levels(vec)[i+start-1]}
	cbind(des,tmp)
}

#-------------------------------------------------------------------------------
design.mat <- function(saminfo)
{
	tmp <- which(colnames(saminfo) == 'Batch')
	tmp1 <- as.factor(saminfo[,tmp])
	msg("  => Found",nlevels(tmp1),'batches');
	design <- build.design(tmp1,start=1)
	ncov <- ncol(as.matrix(saminfo[,-c(1:2,tmp)]))
	msg("  => Found",ncov,'covariate(s)');
	if(ncov>0){
		for (j in 1:ncov){
			tmp1 <- as.factor(as.matrix(saminfo[,-c(1:2,tmp)])[,j])
			design <- build.design(tmp1,des=design)
			}
		}
	design
}

#-------------------------------------------------------------------------------
list.batch <- function(saminfo)
{
	tmp1 <- as.factor(saminfo[,which(colnames(saminfo) == 'Batch')])
	batches <- NULL
	for (i in 1:nlevels(tmp1)){batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))}
	batches
}

#-------------------------------------------------------------------------------
aprior <- function(gamma.hat)
{
	m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2
}

#-------------------------------------------------------------------------------
it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001)
{
	n <- apply(!is.na(sdat),1,sum)
	g.old <- g.hat
	d.old <- d.hat
	change <- 1
	count <- 0
	while(change>conv){
		g.new <- postmean(g.hat,g.bar,n,d.old,t2)
		sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
		d.new <- postvar(sum2,n,a,b)
		change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
		g.old <- g.new
		d.old <- d.new
		count <- count+1
		}
	adjust <- rbind(g.new, d.new)
	rownames(adjust) <- c("g.star","d.star")
	adjust
}

#-------------------------------------------------------------------------------
bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}
postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}
