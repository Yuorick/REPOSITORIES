//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "DiagrSinX.h"
extern const double TET0707 = 1.39115573;

// нахождение обощенного угла при котором амплитуда диаграммы равна sqrt(2)/2
double findTet07()
{
 double tet0 = M_PI/2;

 int i =0;
 double a = sqrt(2.)/2.;
 for ( i =0; i < 10; i++)
 {
  double del = -(fncDiagrSinx_div_x(tet0) - a)/ fncDerivDiagrSinx_div_x (tet0);
  tet0 += del;
  if (fabs(del) < 0.0000001) break;

 }
 return tet0;
}

double fncDiagrSinx_div_x(double tet)
{
if (fabs(tet)< 0.0000001) return 1.;

return sin(tet)/tet;
}

double fncDerivDiagrSinx_div_x(double tet)
{
if (fabs(tet)< 0.0000001) return 0.;
return (tet * cos(tet) -sin(tet))/ tet/tet;
}

double fncDeriv2DiagrSinx_div_x(double tet)
{
if (fabs(tet)< 0.0000001) return -1./3.;
return (- sin(tet)/tet - 2. *cos(tet)/tet/tet + 2. * sin(tet)/tet/tet/tet);
}

double fncDeriv3DiagrSinx_div_x(double tet)
{
if (fabs(tet)< 0.0000001) return 0.;
return (- cos(tet)/tet + 3.*sin(tet)/tet/tet + 6.*cos(tet)/tet/tet/tet - 6.*sin(tet)/tet/tet/tet/tet);
}

// диаграмма sin(ax)/ax
double fncDiagrSimple(double valWidthDgr, double valTetRad)
{
double coeff =  TET0707 / valWidthDgr;
return fncDiagrSinx_div_x(valTetRad * coeff);
}
// производная диаграмма sin(ax)/ax
double fncDerivDiagrSimple(double valWidthDgr, double valTetRad)
{
double coeff =  TET0707 / valWidthDgr;
return coeff* fncDerivDiagrSinx_div_x(coeff *valTetRad);
}
// иторая производная диаграмма sin(ax)/ax
double fncDeriv2DiagrSimple(double valWidthDgr, double valTetRad)
{
double coeff =  TET0707 / valWidthDgr;
return coeff * coeff* fncDeriv2DiagrSinx_div_x(coeff *valTetRad);
}


#pragma package(smart_init)
