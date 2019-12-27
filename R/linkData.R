#' Bayes Estimates of Bipartite Matchings
#' 
#' Bayes point estimates of bipartite matchings that can be obtained 
#' in closed form according to Theorems 1, 2 and 3 of Sadinle (2017).  
#' 
#' @param Zchain	matrix as the output \code{$Z} of the function \code{\link{runGibbs}}, with \code{n2} rows and \code{nIter} columns containing a chain 
#'		 		of draws from a posterior distribution on bipartite matchings.  Each column indicates the records in datafile 1 to which the records in datafile 2 are matched according to that draw.
#' 
#' @param n1		number of records in datafile 1.
#' 
#' @param lFNM		individual loss of a false non-match in the loss functions of Sadinle (2017), default \code{lFNM=1}.
#' 
#' @param lFM1		individual loss of a false match of type 1 in the loss functions of Sadinle (2017), default \code{lFM1=1}.
#' 
#' @param lFM2		individual loss of a false match of type 2 in the loss functions of Sadinle (2017), default \code{lFM2=2}.
#' 
#' @param lR		individual loss of 'rejecting' to make a decision in the loss functions of Sadinle (2017), default \code{lR=Inf}.
#' 
#' @details	Not all combinations of losses \code{lFNM, lFM1, lFM2, lR} 
#' 			are supported.  The losses have to be positive numbers and satisfy one of three conditions:
#' \enumerate{
#' 			\item conditions of Theorem 1 of Sadinle (2017):
#'			\code{(lR == Inf) & (lFNM <= lFM1) & (lFNM + lFM1 <= lFM2)}
#' 			\item conditions of Theorem 2 of Sadinle (2017):
#'			\code{((lFM2 >= lFM1) & (lFM1 >= 2*lR)) | ((lFM1 >= lFNM) & (lFM2 >= lFM1 + lFNM))}
#' 			\item conditions of Theorem 3 of Sadinle (2017):
#'			\code{(lFM2 >= lFM1) & (lFM1 >= 2*lR) & (lFNM >= 2*lR)}
#' }
#' If one of the last two conditions is satisfied, the point estimate might be partial, meaning that there
#' might be some records in datafile 2 for which the point estimate does not include a linkage decision.
#' For combinations of losses not supported here, the linear sum assignment problem outlined by Sadinle (2017)
#' needs to be solved.  
#' 
#' @return 	a vector containing the point estimate of the bipartite matching.  If \code{lR != Inf} the output might be a partial estimate.
#'			A number smaller or equal to \code{n1} in entry \code{j} indicates the record in datafile 1 to which record \code{j} in datafile 2 
#'			gets linked, a number \code{n1+j} indicates that record \code{j} does not get linked to any record in datafile 1, and the value \code{-1} 
#'			indicates a 'rejection' to link, meaning that the correct linkage decision is not clear.
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
#' 
#' ## discard first 100 iterations of Gibbs sampler
#' 
#' ## full estimate of bipartite matching (full linkage)
#' fullZhat <- linkData(chain$Z[,-c(1:100)], n1=nrow(df1), lFNM=1, lFM1=1, lFM2=2, lR=Inf)
#' 
#' ## partial estimate of bipartite matching (partial linkage)
#' partialZhat <- linkData(chain$Z[,-c(1:100)], n1=nrow(df1), lFNM=1, lFM1=1, lFM2=2, lR=.1)

linkData <- function(Zchain, n1, lFNM=1, lFM1=1, lFM2=2, lR=Inf){
	
	# control the input
	if(!is.matrix(Zchain)) stop("Zchain should be a matrix")
	n2 <- nrow(Zchain)
	# make sure the labels in Zchain are within the expected range
	if(max(Zchain) > n1 + n2) stop("Labels in Zchain exceed n1+n2")
	# - positive losses
	C0 <- (lFNM > 0) & (lFM1 > 0) & (lFM2 > 0) & (lR > 0)
	# - conditions of Theorem 1 of Sadinle (2017)
	C1 <- (lR == Inf) & (lFNM <= lFM1) & (lFNM + lFM1 <= lFM2)
	# - conditions of Theorem 2 of Sadinle (2017)
	C2 <- ((lFM2 >= lFM1) & (lFM1 >= 2*lR)) | ((lFM1 >= lFNM) & (lFM2 >= lFM1 + lFNM))
	# - conditions of Theorem 3 of Sadinle (2017)
	C3 <- (lFM2 >= lFM1) & (lFM1 >= 2*lR) & (lFNM >= 2*lR)
	# check we can handle the specified losses
	if(!C0) stop("Losses need to be positive")
	if(!any(c(C1,C2,C3))) stop("Invalid configuration of losses")
		
	# temporarily replace all nonlink labels by n1+1
	Zchain[Zchain > n1+1] <- n1+1
	tableLabels <- apply(Zchain, 1, tabulate, nbins=max(Zchain))
	tableLabels <- tableLabels/ncol(Zchain)
	probNoLink <- tableLabels[n1+1,]
	# find marginal best option for each record based only on probability
	maxProbOption <- apply(tableLabels, 2, which.max)
	maxProbOption[maxProbOption==n1+1] <- (n1+1:n2)[maxProbOption==n1+1]
	probMaxProbOption <- apply(tableLabels, 2, max)
	maxProbOptionIsLink <- maxProbOption <= n1
	
	if(C1){# if not using reject option and conditions of Theorem 1
	
		Zhat <- (n1+1):(n1+n2)
		tholdLink <- lFM1/(lFM1+lFNM) + 
			(lFM2-lFM1-lFNM)*(1 - probNoLink - probMaxProbOption)/(lFM1+lFNM)
		Zhat[maxProbOptionIsLink & (probMaxProbOption > tholdLink)] <- 
			maxProbOption[maxProbOptionIsLink & (probMaxProbOption > tholdLink)]
	
	}else{# if using reject option 
		if(C3){# if conditions of Theorem 3 are satisfied
		
			Zhat <- rep(-1,n2) # represents the reject option
			tholdLink <- 1 - lR/lFM1 + (lFM2-lFM1)*(1 - probNoLink - probMaxProbOption)/lFM1
			Zhat[maxProbOptionIsLink & (probMaxProbOption > tholdLink) ] <- 
				maxProbOption[maxProbOptionIsLink & (probMaxProbOption > tholdLink) ]
			noLinkDec <- probNoLink > 1-lR/lFNM
			Zhat[noLinkDec] <- ((n1+1):(n1+n2))[noLinkDec]
		
		}else{ # Theorem 2
			
			# compute equation (6) in Sadinle (2017)
			tableLabels[-n1-1,] <- t( lFM2*(t(1-tableLabels[-n1-1,])-tableLabels[n1+1,]) + 
										lFM1*tableLabels[n1+1,] )
			tableLabels[n1+1,] <- lFNM*(1-tableLabels[n1+1,])
			# find the options with the marginal minimal loss
			lossMinLossOption <- apply(tableLabels, 2, min)
			minLossOption <- apply(tableLabels, 2, which.min)
			noLinkDec <- minLossOption == n1+1
			minLossOption[noLinkDec] <- ((n1+1):(n1+n2))[noLinkDec]
			Zhat <- rep(-1,n2) # represents the reject option
			Zhat[lossMinLossOption < lR] <- minLossOption[lossMinLossOption < lR]
		
		}
	}
	
	return(Zhat)
}
