#ifndef LINDIFFEQ_H
#define LINDIFFEQ_H


class QLinDiffEq
{
public:

    static void calcFundamentalMtrx(const double *arrA, const int dimA, const double VAlT, double *arrF);

    static void calcFundamentalMtrx_dim2(const double *arrA,  const double VAlT, double *arrF);
};

#endif // LINDIFFEQ_H
