#include "MinSquare.h"
#include "MatrixProccess.h"

QMinSquare::QMinSquare()
{

}
//----------------------------------
double QMinSquare::solvMinSquare( double *arrA,const int qRows , const int qCols , double *arrb, double *arrX )
{
    double *arrATr = new double [qCols * qRows];
    double *arrATrA = new double [qCols * qCols];
    double *arrATrAInv = new double [qCols * qCols];
    double *arrATrb = new double [qCols];
    MatrTransp(arrA, qRows , qCols, arrATr);
    MtrxMultMatrx( arrATr, qCols, qRows ,  arrA, qCols, arrATrA) ;
    MtrxMultMatrx( arrATr, qCols, qRows ,  arrb, 1, arrATrb) ;
    InverseMtrx( arrATrA, qCols,  arrATrAInv);
    MtrxMultMatrx(arrATrAInv, qCols, qCols ,   arrATrb, 1, arrX) ;
    MatrxMultScalar(arrX, qCols, 1, -1.,arrX);

   // double valNormb = NormVect(arrb, qRows);
    //double valReturn = valNormb * valNormb + ScalProduct(arrATrb , arrX, qCols) ;
    double *arrT0 = new double [qRows];
    double *arrT1 = new double [qRows];
    MtrxMultMatrx( arrA, qRows , qCols,   arrX, 1, arrT0) ;

    MtrxSumMatrx(arrT0, arrb,1, qRows, arrT1) ;
    double valReturn = NormVect(arrT1, qRows);
    valReturn = valReturn * valReturn;
    delete []arrATr;
    delete []arrATrA;
    delete []arrATrb;
    delete []arrATrAInv;
    delete []arrT0;
    delete []arrT1;
    return valReturn;
}


//----------------------------------
long double  QMinSquare::solvMinSquare( long double  *arrA,const int qRows , const int qCols , long double  *arrb, long double  *arrX )
{
    long double  *arrATr = new long double  [qCols * qRows];
    long double  *arrATrA = new long double  [qCols * qCols];
    long double  *arrATrAInv = new long double  [qCols * qCols];
    long double  *arrATrb = new long double  [qCols];
    MatrTransp(arrA, qRows , qCols, arrATr);
    MtrxMultMatrx( arrATr, qCols, qRows ,  arrA, qCols, arrATrA) ;
    MtrxMultMatrx( arrATr, qCols, qRows ,  arrb, 1, arrATrb) ;
    InverseMtrx( arrATrA, qCols,  arrATrAInv);
    MtrxMultMatrx(arrATrAInv, qCols, qCols ,   arrATrb, 1, arrX) ;
    MatrxMultScalar(arrX, qCols, 1, -1.,arrX);

   // long double  valNormb = NormVect(arrb, qRows);
    //long double  valReturn = valNormb * valNormb + ScalProduct(arrATrb , arrX, qCols) ;
    long double  *arrT0 = new long double  [qRows];
    long double  *arrT1 = new long double  [qRows];
    MtrxMultMatrx( arrA, qRows , qCols,   arrX, 1, arrT0) ;

    MtrxSumMatrx(arrT0, arrb,1, qRows, arrT1) ;
    long double  valReturn = NormVect(arrT1, qRows);
    valReturn = valReturn * valReturn;
    delete []arrATr;
    delete []arrATrA;
    delete []arrATrb;
    delete []arrATrAInv;
    delete []arrT0;
    delete []arrT1;
    return valReturn;
}
