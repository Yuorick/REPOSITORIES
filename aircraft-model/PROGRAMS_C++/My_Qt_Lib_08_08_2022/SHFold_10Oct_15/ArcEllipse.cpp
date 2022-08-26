//---------------------------------------------------------------------------




#include "ArcEllipse.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include <stdio.h>
 #include <math.h>
 #include "URPolyLineZ.h"
 #include "Equations.h"
 #include "URPolygon.h"
 #include "Comp.h"



//---------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
 TArcEllipse ::TArcEllipse()
{
	ma = 4.;
	mb = 3. ;
	m_x0 = -4.;
	m_x1 = 0.;
}

 // парам констр
TArcEllipse :: TArcEllipse( const double a, const double b)

 {
   if (b >= a)
   {
	//ShowMessage (L"Error Constructor TArcEllipse");
	return ;
   }
   ma = a;
   mb = b ;
   m_x0 = -ma  ;
   m_x1 = 0;


 }
 // парам констр
 TArcEllipse :: TArcEllipse( const double a, const double b, const double x0, const double x1)

 {
   if (b >= a)
   {
	//ShowMessage (L"Error Constructor TArcEllipse");
	return ;
   }
   ma = a;
   mb = b ;
   m_x0 = x0  ;
   m_x1 = x1;


 }
  // парам констр 3
 TArcEllipse :: TArcEllipse( const double a, const double b, const int isign)

 {
   if (b >= a)
   {
	//ShowMessage (L"Error Constructor TArcEllipse");
	return ;
   }
   ma = a;
   mb = b ;
   double f = fncFocus ();
   if (isign > 0)
   {
	 m_x0 = f ;
	 m_x1 = a ;
   }
   else
   {
	m_x0 = -a ;
	m_x1 = -f ;
   }

 }

 // парам констр  создания из сумы расстояний до фокусов  и фокуса
TArcEllipse :: TArcEllipse( const double l, const double f, double *p)

 {

   ma =  l / 2. ;
   double temp =  ma *  ma - f * f ;
   if (temp <0.)
   {
	//ShowMessage (L"Error Constructor TArcEllipse 1");
	return ;
   }
   mb = sqrt(temp) ;
   m_x0 = -ma ;
   m_x1 = 0;
 }


 // парам констр  создания из сумы расстояний до фокусов  и фокального расстояния
TArcEllipse :: TArcEllipse( const double l, const double f,  const int isign, double *p)

 {

   ma = l / 2. ;
   double temp =  ma *  ma - f * f ;
   if (temp <0.)
   {
	//ShowMessage (L"Error Constructor TArcEllipse 1");
	return ;
   }
   mb = sqrt(temp) ;
    if (isign > 0)
   {
	 m_x0 = f ;
	 m_x1 = ma ;
   }
   else
   {
	m_x0 = -ma ;
	m_x1 = -f ;
   }
 }



 // оператор присваивания
 TArcEllipse TArcEllipse::operator=(TArcEllipse  R)
 {
	ma  = R.ma ;
	mb  = R.mb;
	m_x0 = R.m_x0 ;
	m_x1 = R.m_x1 ;


	return *this ;
 }

 // конструктор копирования
 TArcEllipse::TArcEllipse (const TArcEllipse &R)
 {
   ma  = R.ma ;
   mb  = R.mb;
   m_x0 = R.m_x0 ;
   m_x1 = R.m_x1 ;

 }

double  TArcEllipse::fncFocus ()
{
	return sqrt( ma * ma - mb * mb);
}
// пересечение эллипса с прямой, заданнной 2-мя точками  pPontsInp
// возвращает к-во точек пересения
// если одна точка пересения (касание) то она хранится в  pPontsOut[0]
// если 2 точки, то они хранятся в  pPontsOut[0] и pPontsOut[1]
int  TArcEllipse::fncLineIntersectEllips (TURPointXY *pPontsInp,TURPointXY *pPontsOut)

{
 double arrV [2] = {0.} ;
  arrV [0] = pPontsInp[1].X -  pPontsInp[0].X  ;
  arrV [1] = pPontsInp[1].Y -  pPontsInp[0].Y  ;
  double x0 = pPontsInp[0].X ;
  double y0 = pPontsInp[0].Y ;
  double vala = arrV [1]* arrV [1]* ma * ma + arrV [0]* arrV [0]* mb * mb ;
  double valb = 2. *(x0 * arrV [0] * mb * mb + y0 * arrV [1] * ma * ma );
  double valc =  mb * mb * x0 * x0 +  ma * ma * y0 * y0 -  ma * ma *  mb * mb;
  // Решение квадраьного улавнения a*x*x + b*x + c =0
// возвращает
// 0 - 2 действительных некраьных корня
// 1 - 2 действительных кратных корня
// 2- имеется по крайней мере один нулевой корень a!= 0, c=0
// 3 - 2 комплексно сопряженных корня
// 4 -  1 действительный корень (a =0)
// 5  - несовместность a=b=0, c!=0
// 6 - тождество )a=b=c=0
	TComp cx1, cx2 ;
	double arrRoots[2] = {0.} ;
	int quantRoots = 0;
	int irez = SolvEq2( vala, valb, valc, cx1, cx2) ;
	switch(irez)
	{   case 0:
		 arrRoots[0] = cx1.m_Re ;
		 arrRoots[0] = cx2.m_Re ;
		 quantRoots = 2;
		 break;

		 case 1:
         arrRoots[0] = cx1.m_Re ;
		 quantRoots = 1;
		 break ;

		default:
		return 0 ;
	}
	for (int i = 0; i < quantRoots; i++)
	{
	  pPontsOut[i].X = x0 +  arrV [0] * arrRoots[i] ;
	  pPontsOut[i].Y = y0 +  arrV [1] * arrRoots[i] ;
	}
 return quantRoots ;
}

 // пересечение  отрезка и линии эллипса
// возвращает к-во точек пересечения
//arrPntRez [2] - массив точек пересечения

 int TArcEllipse::fncSegmentCutArcEllipse(TURPointXY pntSegm0, TURPointXY pntSegm1 // сегмент
							,TURPointXY *arrPntRez // точек ассив а пересечения
								)
{
  // перевод вершин сегмента  в СК в которойэллипс является
  // елдиничной окружностью
   TURPointXY arrSegmPnt[2];
   arrSegmPnt[0].X = pntSegm0.X / ma ;
   arrSegmPnt[0].Y = pntSegm0.Y / mb ;
   arrSegmPnt[1].X = pntSegm1.X / ma ;
   arrSegmPnt[1].Y = pntSegm1.Y / mb ;
   ///

 // пересечение  отрезка и линии окружности
// возвращает 0 если не пересекаются
// 10 если касаются  внутри отрезка
// 11 точка  pntSegm0 (=arrPntRez[0]) на окружности,  pntSegm1 снаружи, весь отрезок снаружи
// 12 точка pntSegm1 на окружности (=arrPntRez[0] ),  pntSegm0 снаружи, весь отрезок снаружи
// 13 точка pntSegm0 на окружности (=arrPntRez[0] ),  pntSegm1 внутри, весь отрезок внутри
// 14 точка pntSegm1 (=arrPntRez[0] )на окружности ,  pntSegm0 внутри, весь отрезок внутри
// 15 касание в точке  pntSegm1 (=arrPntRez[0] )
// 16 касание в точке  pntSegm0 (=arrPntRez[0] )
// 17 точка pntSegm0 внутри ,  точка pntSegm1 снаружи, 1 точка пересечения
// 18 точка pntSegm1внутри ,  точка pntSegm0 снаружи, 1 точка пересечения
// 20  точка pntSegm0 снаружи,  pntSegm1 снаружи, 2 точки пересечения
// 21 точка pntSegm0 на окружности (=arrPntRez[0] ),  точка pntSegm1 на окружности (=arrPntRez[1] ) 2 точки пересечения
// 23  точка pntSegm1 на окружности (=arrPntRez[1] ),  pntSegm0 снаружи, 2 точки пересечения
// 22  точка pntSegm0 на окружности (=arrPntRez[1] ),  pntSegm1 снаружи, 2 точки пересечения
// arrPntRez- масссив точе пересечения,
TURPointXY pntCrclCentre(0.,0.) ;
TURPointXY arrPntRez1 [2] ;
arrPntRez1 [0] =  pntCrclCentre ;
arrPntRez1 [1] =  pntCrclCentre ;
int irez = TURPolygon::fncSegmentCutCircleLine(  pntCrclCentre,1.// окружность
							, arrSegmPnt[0],  arrSegmPnt[1]
							,arrPntRez1 // точек ассив а пересечения
								) ;
   if(irez == 0) return 0;
   int irez1 = 0 ;
   for (int i = 0; i < 2; i++)
   {
	 arrPntRez[i].X =  arrPntRez1[i].X * ma ;
	 arrPntRez[i].Y =  arrPntRez1[i].Y * mb ;
   }

   if (irez < 20)
   {
	 if (!(Is_X_BelongeSegm( arrPntRez[0].X, m_x0, m_x1) && (arrPntRez[0].Y >= -TOLRNC) ) )
	  {
		 irez1 = 0 ;
		 return irez1;
	  }

	 return 1;
   }

   bool b0 =  (Is_X_BelongeSegm( arrPntRez[0].X, m_x0, m_x1) && (arrPntRez[0].Y >= -TOLRNC));
   bool b1 =  (Is_X_BelongeSegm( arrPntRez[1].X, m_x0, m_x1) && (arrPntRez[1].Y >= -TOLRNC));
   if ( b0 && (!b1 )) return 1;
   if ( b1 && (!b0) )
   {
	arrPntRez[0] = arrPntRez[1];
	return 1;
	}
	if ( (!b1) && (!b0) ) return 0;
	return 2;

}



double TArcEllipse::calcFncValue( const double valx)
{
	return mb * sqrt(ma* ma - valx*valx)/ ma ;
}
double TArcEllipse::calc_dY_po_dX( const double valx)
{
	return -mb * valx/ ma/sqrt(ma* ma - valx*valx) ;
}

bool TArcEllipse::Is_X_BelongeSegm( const double valx, const double vala,const double valb)
{
	bool ireturn = false;
	if ((valx <= (valb + TOLRNC) )&& ( (vala - TOLRNC) <= valx)) ireturn = true ;

	return ireturn ;
}
void TArcEllipse::ShowMe(wchar_t *FileName)
{
   int iN = 6000;
   TURPolyLine line(*this,iN);
   line.WriteSetSHPFiles(FileName,&line, 1 ) ;
}
void TArcEllipse::ShowFullGraph(wchar_t *FileName)
{
 TArcEllipse tempEll = *this;
 tempEll.m_x0 = -ma;
 tempEll.m_x1 =  ma;
 tempEll.ShowMe(FileName) ;

}



