#ifndef MINSQUARE_H
#define MINSQUARE_H


class QMinSquare
{
public:
    QMinSquare();
    static double solvMinSquare( double *arrA,const int qRows , const int qCols , double *arrb, double *arrX );

    static long double  solvMinSquare( long double  *arrA,const int qRows , const int qCols , long double  *arrb, long double  *arrX );
};

#endif // MINSQUARE_H
