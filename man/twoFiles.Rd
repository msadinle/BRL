 \name{twoFiles}
 \docType{data}
 \alias{twoFiles}
 \alias{df1}
 \alias{df2}
 \alias{df2ID}
 \title{Two Datasets for Record Linkage}
 \description{
  Two data frames, \code{df1} and \code{df2}, containing 300 and 150 records of artificially created 
  individuals, where 50 of them are included in both datafiles.  In addition, the vector \code{df2ID} 
  contains one entry per record in \code{df2} indicating the true matching between the datafiles, codified as follows: 
  a number smaller or equal to \code{n1=300} in entry \code{j} 
  indicates the record in \code{df1} to which record \code{j} in \code{df2} 
	truly matches, and a number \code{n1+j} indicates that record \code{j} in \code{df2} does not match any record in \code{df1}.
 }
 \usage{data(twoFiles)}
 \source{Extracted from the datafiles used in the simulation studies of Sadinle (2017).  The datafiles were originally 
 generated using code provided by Peter Christen (\url{https://users.cecs.anu.edu.au/~Peter.Christen/}).}
 \examples{
data(twoFiles)

n1 <- nrow(df1)

## the true matches
cbind( df1[df2ID[df2ID<=n1],], df2[df2ID<=n1,] )

## alternatively
df1$ID <- 1:n1
df2$ID <- df2ID
merge(df1, df2, by="ID")

## all the records in a merged file
merge(df1, df2, by="ID", all=TRUE)

}
 \references{
 Mauricio Sadinle (2017). Bayesian Estimation of Bipartite Matchings for Record Linkage. \emph{Journal of the
American Statistical Association} 112(518), 600-612. [\href{https://doi.org/10.1080/01621459.2016.1148612}{Published}] [\href{https://arxiv.org/abs/1601.06630}{arXiv}]
 }