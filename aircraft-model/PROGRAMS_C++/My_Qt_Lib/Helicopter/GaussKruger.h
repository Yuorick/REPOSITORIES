#ifndef GAUSSKRUGER_H
#define GAUSSKRUGER_H


void    mhGeo2Gauss_( double B, double L, double *X, double *Y );
void    mhGauss2Geo_( double X, double Y, double *B, double *L );
double  mhMod2Pi_( double a_fVal );

#endif // GAUSSKRUGER_H
