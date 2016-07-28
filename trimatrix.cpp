#include <stdlib.h>
#include "trimatrix.h"

/**************************************
 algorithm for solving system of linear equations with tridiagonal matrix
***************************************/
void solve_trimatrix_system (trimatrix_system tms)
{
	double *alf, *bet;
	int n;
	int i;

	n=tms.dim;
	alf = (double *) calloc(n+1, sizeof(double));
	bet = (double *) calloc(n+1, sizeof(double));

	alf[0]=bet[0]=0.0;
	for( i=0; i<=n-1; i++)
	{
		alf[i+1] = tms.up[i]/(-tms.diag[i] - tms.lo[i]*alf[i] );
		bet[i+1]=( tms.lo[i]*bet[i] - tms.rhs[i] )/(-tms.diag[i] - tms.lo[i]*alf[i] );
	}

	tms.c[n-1] = bet[n];
	for( i=n-2; i>=0; i--)
	{
		tms.c[i] = alf[i+1]*tms.c[i+1]+bet[i+1];
	}

	free(alf);
	free(bet);
}
