#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void BRLGibbs(int *, int *, int *,
				int *, int *,
				int *, int *, 
				int *, int *,
				double *, double *, 
				double *, double *,
				int *, double *, double *);

static const R_CMethodDef CEntries[] = {
    {"BRLGibbs", (DL_FUNC) &BRLGibbs, 16},
    {NULL, NULL, 0}
};

void R_init_BRL(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
