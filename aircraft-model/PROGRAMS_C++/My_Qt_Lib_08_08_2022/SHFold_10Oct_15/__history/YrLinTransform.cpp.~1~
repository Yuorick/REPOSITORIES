//---------------------------------------------------------------------------


#pragma hdrstop

#include "YrLinTransform.h"
#include <vcl.h>
#include <math.h>

//---------------------------------------------------------------------------

#pragma package(smart_init)

const double  EPS_ABS = 0.000000001 ;
const double  PI = 3.1415926 ;
void LinTransf(const double x0, const double y0,
				const double alf0,
				const double x, const double y
				, double *xTr,  double *yTr
				)
// линейное преобразование координат
// x0, y0 - координаты точки (0,0) новой системы координат (— )
// alf0 - угол  между ос€ми 0X старой и новой — 
// исчисл€етс€ против часовой стрелки от 0X1 к 0X2
// x,y - координаты точки в старой — 
// xTr,  yTr  - координаты той же точки в новой — 
{
   double cosa = cos(alf0) ;
   double sina = sin(alf0) ;
	double dx0 = cosa * x0 + sina * y0  ;
   double dy0 = - sina * x0 + cosa * y0  ;
  /* *xTr = cosa * x + sina * y - x0 ;
   *yTr = - sina * x + cosa * y - y0 ; */
   double dx = cosa * x + sina * y - dx0 ;
   double dy = - sina * x + cosa * y - dy0 ;
  // cout<<"dx = "<< dx<< "dy = "<< dy<< endl ;
   *xTr = dx ;
   *yTr =dy ;
   //cout<<"cosa ="<< cosa<< "sina = "<< sina<< endl ;
  // cout<< "x0="<< x0<< "y0= " << y0<< endl ;
  //	 cout<< "*xTr = "<< *xTr<<" *yTr = "<<*yTr<< " x =  "<<x<< " y= "<< y<<" alf0= "<<alf0<<" x0 = "<< x0<< endl;system("pause") ;

	return ;
}
double AngleBetweenVectors(const double x1,const double y1,const double x2,const double y2)
// вычисление угла между 2-м€ векторами исход€щими из начала координат
// угол исчисл€етс€ против часовой стрелкипо направлению от 1-го ко 2-му
{
  if ( fabs( (x1 *x1 + y1*y1) * (x2 *x2 + y2*y2)) < EPS_ABS )
  {
	 ShowMessage(L"There are not possible to calculate angle between vectors\n in AngleBetweenVectors") ;
	 return 200;
  }
  double alf1 =  asin(y1/dist2(x1,y1,0,0)) ;
  if ( x1 < 0)  alf1 =  PI -alf1 ;
  double alf2 =  asin(y2/dist2(x2,y2,0,0)) ;
  if ( x2 < 0)  alf2 =  PI -alf2 ;

  return (alf2 -alf1) ;
}
double dist2(const double x1,const double y1,const double x2,const double y2)
{
	return sqrt( (x1-x2)* (x1-x2) + (y1-y2) *(y1-y2));
}
