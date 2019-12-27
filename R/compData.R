#' Creation of Comparison Data 
#'
#' Create comparison vectors for all pairs of records coming from 
#' two datafiles to be linked.
#'
#' @param df1,df2 two data frames to be linked, containing the same number of columns,
#'				with the same variable names and classes, in the same order.  The columns represent 
#'				the comparison fields to be used for the linkage.  Currently, only columns of class
#'				\code{character} and \code{factor} are supported.  Without loss of generality, 
#' 				\code{df1} is assumed to have no less records than \code{df2}.
#' @param types	a vector of characters indicating the comparison types per comparison field.  The options
#'				are: \code{"Lev"} for comparisons based on the Levenshtein distance normalized to \eqn{[0,1]}, with \eqn{0}  
#'				indicating no disagreement and \eqn{1} indicating maximum disagreement; and 
#'				\code{"bin"} for binary comparisons (agreement/disagreement). 
#' 				The default is \code{"Lev"} for columns of class \code{character} and \code{"bin"} for columns of class \code{factor}.
#' @param breaks	for comparisons based on the normalized Levenshtein distance, a vector of length \eqn{L} of break 
#'				points for the interval \eqn{[0,1]} to obtain \eqn{L+1} levels of disagreement.
#'				The default is \code{breaks=c(.001,.25,.5)}.
#' @return a list containing: 
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
#' 
#' @references Mauricio Sadinle (2017). Bayesian Estimation of Bipartite Matchings for Record Linkage. \emph{Journal of the
#' American Statistical Association} 112(518), 600-612.
#' 
#' @examples
#' data(twoFiles)
#' 
#' myCompData <- compData(df1, df2, types=c("Lev","Lev","bin","bin"), breaks=c(.001,.25,.5))
#' 
#' ## same as 
#' myCompData <- compData(df1, df2)

compData <- function(df1, df2, types=NULL, breaks=c(.001,.25,.5)){
	
	# control the input
	
	if(!is.data.frame(df1)) stop("df1 is not a data frame")
	if(!is.data.frame(df2)) stop("df2 is not a data frame")
	
	if(ncol(df1)!=ncol(df2)) stop("df1 and df2 have different numbers of columns")
	
	if( any(names(df1) != names(df1)) ) stop("df1 and df2 have different column names")
	
	df1ColClass <- sapply(df1, class)
	df2ColClass <- sapply(df2, class)
	
	if( any(df1ColClass != df2ColClass) ) stop("df1 and df2 have different column classes")
	
	if( !(all(df1ColClass %in% c("character", "factor"))) )
		stop("Some columns in df1 are not of class character or factor")
		
	if( !(all(df2ColClass %in% c("character", "factor"))) )
		stop("Some columns in df2 are not of class character or factor")
	
	n1 <- nrow(df1)
	n2 <- nrow(df2)
	
	if(n2 > n1) stop("df2 has more records than df1")
	
	if(is.null(types)) types <- c("Lev","bin")[match(df1ColClass, c("character", "factor"))]
	
	F <- length(types)
	
	if(F != ncol(df1)) 
		stop("Vector length of specified comparison types does not match number of datafile columns")
		
	if( any( !(types %in% c("Lev","bin")) ) ) 
		stop("types should be a vector of 'l' or 'b' indicating Levenshtein-based or binary comparison")
	
	if( any(breaks < 0 | breaks > 1) ) stop("breaks should be a vector of numbers in [0,1]")
	
	n <- n1 + n2
	breaks <- c(-Inf,breaks,Inf)
	
	pairInds1 <- rep(1:n1, n2) 
	pairInds2 <- rep(1:n2, each=n1)
	
	compData <- list()
	
	for(fld in 1:F){
		if(types[fld]=="bin"){
		# computing agreement levels for binary comparisons
			same <- df1[pairInds1,fld] == df2[pairInds2,fld]
			AgrLev <- 1*same
			AgrLev[!same] <- 2
			compData[[fld]] <- as.factor(AgrLev)
		}else{
		# computing agreement levels for Levenshtein-based comparisons
			df1[,fld] <- as.character(df1[,fld])
			df2[,fld] <- as.character(df2[,fld])
			
			# adist in the utils package returns the matrix of Levenshtein distances
			lvd <- as.numeric(adist(df1[,fld], df2[,fld]))/
					pmax(nchar(df1[,fld])[pairInds1], nchar(df2[,fld])[pairInds2])
			
			AgrLev <- cut(lvd, breaks=breaks, labels=seq_len(length(breaks)-1))
			compData[[fld]] <- as.factor(AgrLev)
		}
	}
	
	nDisagLevs <- sapply(compData, FUN=nlevels)
		
	# next four lines compute the binary indicators from the agreement levels
	xfun <- function(f) paste("(compData[[",f,"]]==",seq_len(nDisagLevs[f]),")")
	expr1 <- paste(sapply(sapply(seq_len(F),xfun),FUN=paste,collapse=","),collapse=",")
	expr2 <- paste("cbind(",expr1,")")
	comparisons <- eval(parse(text=expr2))
	
	# replacing NAs by FALSE or zeroes is justified by ignorability and CI assumptions (see Sadinle 2017)
	comparisons[is.na(comparisons)] <- FALSE 

	return(list(comparisons=comparisons, n1=n1, n2=n2, nDisagLevs=nDisagLevs))
}
