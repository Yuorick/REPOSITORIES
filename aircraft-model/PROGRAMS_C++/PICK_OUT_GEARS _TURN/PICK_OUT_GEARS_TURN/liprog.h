#ifndef LIPROG_H
#define LIPROG_H

int linprog11( int nvars,  double *f, int nrows, double *a,  double *b
             ,int  nrows_eq,  double *a_eq,  double *b_eq
              ,double *lb, double *ub
            , double *x, double &fval);

#endif // LIPROG_H
