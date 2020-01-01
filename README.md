# BRL

An R package implementing the Beta Record Linkage methodology for probabilistic bipartite record linkage: the task of merging two duplicate-free datafiles that lack unique identifiers.  


## Details	

Beta Record Linkage (BRL, Sadinle, 2017) is a methodology for probabilistic bipartite record linkage, that is, the task of merging two duplicate-free datafiles that lack unique identifiers.  This is accomplished by using the common partially identifying information 
		of the entities contained in the datafiles.  The duplicate-free requirement means that we expect each entity to be represented maximum once in 
		each datafile.  This methodology should not be used with datafiles that contain duplicates nor should it be used for deduplicating a single datafile.  
		
The main function of the package is `BRL` which implements the three main steps of the BRL methodology.

The first step of BRL, accomplished by the function `compareRecords`, consists of constructing comparison vectors for each pair of records from the two datafiles. The current implementation allows binary comparisons (agree/disagree)	and comparisons based on the Levenshtein edit distance.  This can be easily extended to other comparison types, so a resourceful user should be able to construct an object that recreates the output of `compareRecords` for other types of comparisons (so long as they get transformed to levels of disagreement), and still be able to run the next step outside the function `BRL`.  Other types of comparisons will be implemented in the future.  

The second step of BRL, accomplished by the function `bipartiteGibbs`, consists of running a Gibbs sampler that explores the space of bipartite matchings 
 			representing the plausible ways of linking the datafiles.  This sampler is derived from a model for the comparison data and a _beta_ prior 
distribution on the space of bipartite matchings.  See Sadinle (2017) for details. 

The third step of BRL, accomplished by the function `linkRecords`, consists of deriving a point estimate of the bipartite matching 
		(which gives us the optimal way of linking the datafiles) 
			by minimizing the expected value of 
		a loss function that uses different penalties for different types of linkage errors.  The current implementation only supports the 
			Bayes point estimates of bipartite matchings that can be obtained in closed form according to Theorems 1, 2 and 3 of Sadinle (2017).
 			The losses have to be positive numbers and satisfy one of three conditions:

1. Conditions of Theorem 1 of Sadinle (2017):
`(lR == Inf) & (lFNM <= lFM1) & (lFNM + lFM1 <= lFM2)`
2. Conditions of Theorem 2 of Sadinle (2017):
`((lFM2 >= lFM1) & (lFM1 >= 2*lR)) | ((lFM1 >= lFNM) & (lFM2 >= lFM1 + lFNM))`
3. Conditions of Theorem 3 of Sadinle (2017):
`(lFM2 >= lFM1) & (lFM1 >= 2*lR) & (lFNM >= 2*lR)`


If one of the last two conditions is satisfied, the point estimate might be partial, meaning that there
			might be some records in datafile 2 for which the point estimate does not include a linkage decision.
 			For combinations of losses not supported here, the linear sum assignment problem outlined by Sadinle (2017)
			needs to be solved.


## Installation

The BRL R package is currently available only via GitHub.  To install it, use the following code:

```r
#install.packages("remotes") # run this once if you don't have this package already installed
library("remotes") # loads the 'remotes' package in your R session
install_github("msadinle/BRL") # installs BRL from GitHub
library("BRL") # loads the 'BRL' package in your R session
```

## Example

The BRL package comes with two toy datasets `df1` and `df2` along with an ID vector `df2ID`
for `df2` that indicates which entities of `df1` truly match the entities of `df2`.  These can be loaded  using `data(twoFiles)`.  Type `?twoFiles` in your R console for more information.  

The `BRL` function implements all the steps of the BRL methodology, and outputs a vector that indicates which rows of `df1` are linked to the rows of `df2`.
```r
data(twoFiles)
(Zhat <- BRL(df1, df2, flds=c("gname", "fname", "age", "occup"), 
 				types=c("lv","lv","bi","bi")))
n1 <- nrow(df1)
## the linked record pairs
cbind( df1[Zhat[Zhat<=n1],], df2[Zhat<=n1,] )
```
For more details type `?BRL` in your R console.


## Citation

If you use this methodology, please cite the following reference:

Mauricio Sadinle (2017). Bayesian Estimation of Bipartite Matchings for Record Linkage. _Journal of the American Statistical Association_ 112(518), 600-612. [[Published]](https://doi.org/10.1080/01621459.2016.1148612)
[[arXiv]](https://arxiv.org/abs/1601.06630)