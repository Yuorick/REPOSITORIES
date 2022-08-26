//---------------------------------------------------------------------------


#pragma hdrstop

#include "Line.h"

#include "UrPointXY.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

//---------------------------------------------------------------------------
TLine::TLine():TURPolyLine(1,2)
{
}



TLine::TLine(TURPointXY * pPoints):TURPolyLine( pPoints,2)
{
}

TLine::TLine(const TURPointXY  pnt1, const TURPointXY  pnt2):TURPolyLine( pnt1, pnt2)
{
}

bool  TLine::IsVertical()
{
  return(fabs(Points[0].X -  Points[1].X )< 2.* DBL_MIN);
}

bool TLine::calcTang(double *pvalTang)
{
  if (IsVertical())
  {
	return false;
  }
   *pvalTang =  (Points[0].Y -  Points[1].Y )/(Points[0].X -  Points[1].X );
   return true;
}

bool  TLine::calcY(const double x, double *py)
{
  if (IsVertical())
  {
	return false;
  }
  double valTang ;
  calcTang(&valTang);
  *py =  valTang * (x -  Points[0].X) + Points[0].Y;
  return true;

}

bool  TLine::calcPointXY(const double x, TURPointXY *pPointXY)
{
   double y;
  if(!calcY( x, &y) )
  {
	  return false;
  }
  *pPointXY = TURPointXY(x,y);
   return true;

}

#pragma package(smart_init)
