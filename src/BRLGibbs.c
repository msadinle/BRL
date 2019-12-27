/*
 * 	- Change directory:
 *		cd C:/directory
 *	- Compile: 
 *		R CMD SHLIB BRLGibbs.c 
 */

#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <math.h>

void RandomPickNoSort(int *n, double *p, int *ans)
{
    double rU;
    int i;
    int nm1 = n[0] - 1;

    /* compute cumulative probabilities */
    for (i = 1 ; i < n[0]; i++)
	p[i] += p[i - 1];
	
	GetRNGstate();
    /* compute the sample */
	rU = runif(0,p[nm1]);
	for (i = 0; i < nm1; i++) {
	    if (rU <= p[i])
		break;
	}
	ans[0] = i;
    PutRNGstate();
}

/* BRLGibbs: full Gibbs sampler in C*/
void BRLGibbs(	int *nIter, int *n1, int *n2,
				int *nFields, int *nDisagLevs,
				int *uCompPatts, int *nCompPatts, 
				int *codesCompPatts, int *freqAgrLevs,
				double *a, double *b, 
				double *aBM, double *bBM,
				int *Z, double *m, double *u)
{   
	int iter,i,j,k,f,al,al1,co;
	int from[1];
	int pickedPosition[1];
	int nValidLabels[1];
	int nBinAgrLevels[1];
	int Znew[n2[0]];
	int freeLinks[n1[0]];
	int nFreeLinks[1];
	int freqCodesLinks[nCompPatts[0]];
	double uLambda[nCompPatts[0]];
	double smm[nFields[0]];
	double smu[nFields[0]];
	
	/* Initialization */
	
	nBinAgrLevels[0]=0;
	for(f=0; f<nFields[0]; f++){
		nBinAgrLevels[0] += nDisagLevs[f];
	}
	
	double logmMINUSlogu[nBinAgrLevels[0]];	
	double A[nBinAgrLevels[0]];
	double B[nBinAgrLevels[0]];
	double mnew[nBinAgrLevels[0]];
	double unew[nBinAgrLevels[0]];
	
	for(j=0; j<n2[0]; j++){	
		Znew[j] = n1[0] + j;
	}
	
	for(i=0; i<n1[0]; i++){
		freeLinks[i] = 1;
	}
	
	nFreeLinks[0]=n1[0];
	
	for(co=0; co<nCompPatts[0]; co++){
		freqCodesLinks[co]=0;
	}
	
	
	/* Gibbs iterations */
	
	for(iter=0; iter<nIter[0]; iter++){
	
		for(al=0; al<nBinAgrLevels[0]; al++){
			A[al] = 0;
			for(co=0; co<nCompPatts[0]; co++){
				if(uCompPatts[nBinAgrLevels[0]*co+al]) A[al] += freqCodesLinks[co];
			}
		}
		
		for(al=0; al<nBinAgrLevels[0]; al++){
			B[al] = freqAgrLevs[al] - A[al];
		}
		
		/* generate Dirichlet random vectors mnew and unew */
		for(al=0; al<nBinAgrLevels[0]; al++){
			mnew[al] = rgamma(A[al]+a[al],1);
			unew[al] = rgamma(B[al]+b[al],1);
		}
		
		al1=0;
		for(f=0; f<nFields[0]; f++){
			smm[f]=0;
			smu[f]=0;
			for(al=0; al<nDisagLevs[f]; al++){
				smm[f] = smm[f] + mnew[al1];
				smu[f] = smu[f] + unew[al1];
				al1 = al1+1;
			}	
		}
		
		al1=0;
		for(f=0; f<nFields[0]; f++){
			for(al=0; al<nDisagLevs[f]; al++){
				mnew[al1] = mnew[al1]/smm[f];
				unew[al1] = unew[al1]/smu[f];
				m[iter*nBinAgrLevels[0] + al1] = mnew[al1];
				u[iter*nBinAgrLevels[0] + al1] = unew[al1];
				al1 = al1+1;
			}	
		}
		
		/* Compute loglikelihood ratios for sampling links */
		for(al=0; al<nBinAgrLevels[0]; al++){
			logmMINUSlogu[al] = log(mnew[al]) - log(unew[al]);
		}
		
		for(co=0; co<nCompPatts[0]; co++){
			uLambda[co] = 0;
			for(al=0; al<nBinAgrLevels[0]; al++){
				if(uCompPatts[nBinAgrLevels[0]*co+al]) uLambda[co] = uLambda[co] + logmMINUSlogu[al];
			}
		}
		
		from[0]=0;
		
		for(j=0; j<n2[0]; j++){	
			
			if(Znew[j] < n1[0]){ /*release link occupied by j, if any*/
				freeLinks[Znew[j]] = 1;
				nFreeLinks[0] = nFreeLinks[0] + 1;
				freqCodesLinks[codesCompPatts[from[0]+Znew[j]]] -= 1; 
			}
			
			if(nFreeLinks[0]==0){/*not necessary, there are always free links if second file is smaller or equal size*/
				Znew[j]=n1[0]+j;
			}else{
				/*loglikelihood ratios for j and all the records from file 1 that are free to be linked*/
				double Lambda_i_free[nFreeLinks[0]+1]; 
				int which_freelinks[nFreeLinks[0]+1]; 
				
				k=1;
				for (i=0; i<n1[0]; i++){
					if(freeLinks[i]){
						Lambda_i_free[k] = exp(uLambda[codesCompPatts[from[0]+i]]);
						which_freelinks[k] = i;
						k++;
					}
				}
				
				/* *sample new link* */
				/*add the option of no-link*/
				Lambda_i_free[0] = 
					nFreeLinks[0]*(n2[0]-n1[0]+nFreeLinks[0]-1+bBM[0])/(n1[0]-nFreeLinks[0]+aBM[0]);
				nValidLabels[0] = nFreeLinks[0]+1;
				which_freelinks[0] = n1[0]+j;
				/*randomly pick one with prob Lambda_i_free*/
				RandomPickNoSort(nValidLabels, Lambda_i_free, pickedPosition);
				Znew[j] = which_freelinks[pickedPosition[0]];
				Z[iter*n2[0] + j] = Znew[j];
				
				if(Znew[j] < n1[0]){ /*update if j occupies any link*/
					freeLinks[Znew[j]] = 0;
					nFreeLinks[0] = nFreeLinks[0] - 1;
					freqCodesLinks[codesCompPatts[from[0]+Znew[j]]] += 1; 
				}
			}
			from[0] = from[0] + n1[0];
		}
	}
}