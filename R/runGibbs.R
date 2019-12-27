#' Gibbs Sampler Used for Beta Record Linkage 
#'
#' Run a Gibbs sampler to explore the posterior distribution of bipartite matchings 
#' that represent the linkage of the datafiles in beta record linkage.
#' 
#' @param cd a list with the same structure as the output of the function \code{\link{compData}}, containing: 
#' \describe{
#'   \item{\code{comparisons}}{
#' 						matrix with \code{n1*n2} rows, where the comparison pattern for record pair \eqn{(i,j)}
#'						appears in row \code{(j-1)*n1+i}, for \eqn{i} in \eqn{{1,\dots,n1}}, and \eqn{j} in \eqn{{1,\dots,n2}}.  A comparison field with \eqn{L+1} levels of disagreement, 
#'						is represented by \eqn{L+1} columns of TRUE/FALSE indicators.  Missing comparisons are coded as FALSE,
#'						which is justified under an assumption of ignorability of the missing comparisons, see Sadinle (2017).
#'						}
#'   \item{\code{n1,n2}}{the datafile sizes, \code{n1 = nrow(df1)} and \code{n2 = nrow(df2)}.}
#'   \item{\code{nDisagLevs}}{a vector containing the number of levels of
#'					 	disagreement per comparison field.}
#' }
#' @param nIter	number of iterations of Gibbs sampler.
#' @param a,b	hyper-parameters of the Dirichlet priors for the \eqn{m} and \eqn{u} parameters
#' 				in the model for the comparison data among matches and non-matches, respectively.
#' 				These can be vectors with as many 
#'				entries as disagreement levels among all comparison fields.  If specified as positive constants, they 
#'				get recycled to the required length.  If not specified, flat priors are taken.  
#' @param aBM,bBM	hyper-parameters of beta prior on bipartite matchings. Default is \code{aBM=bBM=1}.
#' @param seed	seed to be used for pseudo-random number generation.  By default it sets \code{seed=0}.
#' 
#' @return a list containing: 
#' \describe{
#'   \item{\code{Z}}{matrix with \code{n2} rows and \code{nIter} columns containing the chain of bipartite matchings.  
#'			A number smaller or equal to \code{n1} in row \code{j} indicates the record in datafile 1 to which record \code{j} in datafile 2 
#'			is linked at that iteration, otherwise \code{n1+j}.
#'			}
#'   \item{\code{m,u}}{chain of \eqn{m} and \eqn{u} parameters in the model for the comparison data among matches and non-matches, respectively.}
#' }
#' 
#' @references Mauricio Sadinle (2017). Bayesian Estimation of Bipartite Matchings for Record Linkage. \emph{Journal of the
#' American Statistical Association} 112(518), 600-612.
#' 
#' @examples
#' data(twoFiles)
#' 
#' myCompData <- compData(df1, df2)
#' 
#' chain <- runGibbs(myCompData)

runGibbs <- function(cd, nIter=1000, a=1, b=1, aBM=1, bBM=1, seed=0){
	
	if( !is.numeric(nIter) | (nIter<1) ) stop("nIter should be a positive integer")
	
	if(!is.numeric(a)) stop("'a' should be numeric")
	if(!is.numeric(b)) stop("'b' should be numeric")
	if(!all(a > 0)) stop("'a' should be positive")
	if(!all(b > 0)) stop("'b' should be positive")
	
	nBinAgrLevels <- dim(cd$comparisons)[2]
	
	if(length(a)==1) a <- rep(a,nBinAgrLevels)
	if(length(b)==1) b <- rep(b,nBinAgrLevels)
	
	if(length(a)!=nBinAgrLevels) 
		stop("length of 'a' should equal the total number of disagreement levels among all comparison fields")
	if(length(b)!=nBinAgrLevels) 
		stop("length of 'b' should equal the total number of disagreement levels among all comparison fields")
	
	if( !is.numeric(aBM) | !(aBM>0) ) stop("aBM should be a positive number")
	if( !is.numeric(bBM) | !(bBM>0) ) stop("bBM should be a positive number")
		
	set.seed(seed)
	
	n1 <- cd$n1
	n2 <- cd$n2
	nDisagLevs <- cd$nDisagLevs
	
	sumsBinAgrLevels <- colSums(cd$comparisons)
	
	# encoding binary indicators by decimal representation
	codeBinAgrLevels <- as.vector(cd$comparisons %*% (2^((nBinAgrLevels-1):0)))
	
	# a trick to transform the decimal rep into numbers between 1 and the number of 
	# 	unique binary comparison profiles
	codeBinAgrLevels <- as.numeric(as.factor(codeBinAgrLevels))
	
	# unique binary comparison profiles
	uniqueBinAgrLevels <- unique(cd$comparisons)
	codeUniqBinAgrLevels <- as.vector(uniqueBinAgrLevels %*% (2^((nBinAgrLevels-1):0)))
	# Same trick as before
	codeUniqBinAgrLevels <- as.numeric(as.factor(codeUniqBinAgrLevels))
	
	# we now want to relabel codeBinAgrLevels to indicate the 
	# 	corresponding rows of uniqueBinAgrLevels
	newCodes <- rep(NA,nrow(uniqueBinAgrLevels))
	newCodes[codeUniqBinAgrLevels] <- 1:nrow(uniqueBinAgrLevels)
	newCodeBinAgrLevels <- newCodes[codeBinAgrLevels]
	# now the numbers in newCodeBinAgrLevels correspond to the rows of uniqueBinAgrLevels
	
	#dyn.load(file.path(paste0("BRLGibbs", .Platform$dynlib.ext)))
	
	out <- .C("BRLGibbs", 
			nIter=as.integer(nIter), 	# number of iterations of Gibbs sampler
			n1=as.integer(n1), 			# size of datafile 1
			n2=as.integer(n2),			# size of datafile 2
			nFields=length(nDisagLevs),	# number of fields
			nDisagLevs=as.integer(nDisagLevs), # vector with number of agreement levels per field
			uCompPatts=as.integer(t(uniqueBinAgrLevels)), # unique comparison patterns in binary representation
			nCompPatts=nrow(uniqueBinAgrLevels), # number of unique comparison patterns
			codesCompPatts = as.integer(newCodeBinAgrLevels-1), # decimal codes of comp patterns for all pairs
			freqAgrLevs=as.integer(sumsBinAgrLevels), # frequency of agreement levels among all pairs
			a=as.double(a), 			# hyperparameters for Dirichlet prior for m parameters
			b=as.double(b), 			# hyperparameters for Dirichlet prior for u parameters
			aBM=as.double(aBM), 		# hyperparameter of beta prior on bipartite matchings
			bBM=as.double(bBM), 		# hyperparameter of beta prior on bipartite matchings
			Z=rep(0L,n2*nIter),   		# will contain the iterations of Z
			m=rep(0.0,nBinAgrLevels*nIter),  # will contain the iterations of m
			u=rep(0.0,nBinAgrLevels*nIter),	 # will contain the iterations of u
			PACKAGE="BRL")
	
	Z <- matrix(out$Z,n2,nIter)
	m <- matrix(out$m,nBinAgrLevels,nIter)
	u <- matrix(out$u,nBinAgrLevels,nIter)
	
	return(list(Z=Z+1,m=m,u=u))
}
