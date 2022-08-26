//---------------------------------------------------------------------------
#ifndef LinOptimizationH
#define LinOptimizationH


 int  LinNumericalSolver( int nvars,int bvars, double *f, int nrows,
							double *a, double *b,int  nrows_eq,
					double *a_eq, double *b_eq, double *lb,
							double *ub,int *ix, double *x, double &fval) ;

int QuantNotZeroElements(double *parr, const int lenarr) ;

 // ������ �����������     �������������, � ������� ������������ ����������
int  LinNumericalSolver( int nvars,int bvars, double *f, int nrows,
							double *a, double *b,int  nrows_eq,
					double *a_eq, double *b_eq, double *lb,
							double *ub,int *ix, double *x, double * dualVars, double &fval);

 #endif
