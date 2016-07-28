
/* data structure for a tri-diagonal matrix: */
typedef struct {
		double *lo;   /* lower diagonal           */
		double *diag; /* diagonal                 */
		double *up;   /* upper diagonal           */
		double *rhs;  /* right hand side          */
		double *c;    /* solution - concentration */
		int dim;      /* dimension of matrix =n   */
		} trimatrix_system;

//   system looks like:		/diag0 up0    0      0  ..   0    0    \  / c0   \   /  rhs0  \
//							| lo1  diag1  up1    0  ..   0    0     | | c1   |   |  rhs1  |
//							|  0    lo2  diag2  up2 ..   0    0     | | c2   |   |  rhs2  |
//							|            ....                       | |  ..  | = |   ..   |
//							|            ....                       | |  ..  |   |   ..   |
//							\  0     0     0  ..      lon-1 diagn-1 /  \cn-1/    \ rhsn-1/


extern trimatrix_system tms;

/* algorithm for solving system of linear equations with tridiagonal matrix */
void solve_trimatrix_system (trimatrix_system tms);