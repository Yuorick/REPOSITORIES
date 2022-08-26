//---------------------------------------------------------------------------


#pragma hdrstop

#include "ArcParab.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include <stdio.h>
 #include <math.h>
 #include "URPolyLineZ.h"
 #include "Equations.h"
 #include "URPolygon.h"
#include <vcl.h>
#include "Comp.h"


//---------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
 TArcParab ::TArcParab()
{
	ma = 1.;
	ml = 3. ;
	m_x0 = -2. ;
	m_x1 =  0.;
}

 // парам констр
TArcParab :: TArcParab( const double a, const double l)

 {
   if (l <= 2.* a)
   {
	ShowMessage (L"Error Constructor TArcParab");
	return ;
   }
   ma = a;
   ml = l ;
   m_x0 = ma - ml ;
   m_x1 =  0.;
 }
 // парам констр 2
TArcParab :: TArcParab( const double a, const double l, const int isign)

 {
   if (l <= 2.* a)
   {
	ShowMessage (L"Error Constructor TArcParab");
	return ;
   }
   ma = a;
   ml = l ;
   if (isign > 0)
   {
	 m_x0 = 0;
	 m_x1 = -a;
   }
   else
   {
	m_x0 = -a;
	m_x1 = 0;
   }
 }





 // оператор присваивания
 TArcParab TArcParab::operator=(TArcParab  R)
 {
	ma  = R.ma ;
	ml  = R.ml;
	m_x0 = R.m_x0 ;
	m_x1 = R.m_x1 ;
	return *this ;
 }

 // конструктор копирования
 TArcParab::TArcParab (const TArcParab &R)
 {
	ma  = R.ma ;
	ml  = R.ml;
	m_x0 = R.m_x0 ;
	m_x1 = R.m_x1 ;
 }

void  TArcParab::clcCoeff ( double &vala, double &valb, double &valc)
{
	vala = - 1. / 2./ ml ;
	valb = ma / ml ;
	valc = ( ml * ml - ma * ma) / 2. / ml ;
	return ;
}

double  TArcParab::calcFncValue ( const double valx)
{
  double vala = 0., valb = 0., valc = 0. ;
  clcCoeff (vala,valb, valc) ;
  return  vala * valx* valx +  valb * valx + valc ;

}


// пересечение параболы  с прямой, заданнной 2-мя точками  pPontsInp0 ,pPontsInp1
// возвращает к-во точек пересения
// если одна точка пересения (касание) то она хранится в  pPontsOut[0]
// если 2 точки, то они хранятся в  pPontsOut[0] и pPontsOut[1]
int  TArcParab::fncLineIntersectParab (TURPointXY pPontsInp0,TURPointXY pPontsInp1,TURPointXY *pPontsOut)

{
   double vala = 0., valb = 0., valc = 0. ;
  clcCoeff (vala,valb, valc) ;
 double arrV [2] = {0.} ;
  arrV [0] = pPontsInp1.X -  pPontsInp0.X  ;
  arrV [1] = pPontsInp1.Y -  pPontsInp0.Y  ;
  double x0 = pPontsInp0.X ;
  double y0 = pPontsInp0.Y ;
  double vala1 = arrV [0]* arrV [0]* vala ;
  double valb1 = 2. * x0 * vala * arrV [0]   + valb * arrV [0] -  arrV [1] ;
  double valc1 = vala * x0 * x0 + valb * x0 + valc - y0;

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
	int irez = SolvEq2( vala1, valb1, valc1, cx1, cx2) ;
	switch(irez)
	{   case 0:
		case 2:
		 arrRoots[0] = cx1.m_Re ;
		 arrRoots[1] = cx2.m_Re ;
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

 // пересечение сегмента с дугой параболы
 // возвращает к-во точек пересечения
// arrPntRez- масссив точе пересечения,
 int TArcParab::fncSegmentCutArcParab( TURPointXY pntSegm0, TURPointXY pntSegm1 // сегмент
							,TURPointXY *arrPntRez // точек ассив а пересечения
								)
{
	TURPointXY pPontsOut [2] ;
	int irez = fncSegmentCutParabLine (pntSegm0,pntSegm1, arrPntRez);

	  // исключение точек не лежащих на дуге
	bool b0 = true, b1 = true ;
	switch (irez)
	{
	case 0:
	return 0;

	case 1:
	if (!(Is_X_BelongeSegm( arrPntRez[0].X, m_x0, m_x1) && (arrPntRez[0].Y >= -TOLRNC) ) )
	{
	 return 0;
	}
	return 1 ;

	case 2:
	b0 =  (Is_X_BelongeSegm( arrPntRez[0].X, m_x0, m_x1) && (arrPntRez[0].Y >= -TOLRNC));
	b1 =  (Is_X_BelongeSegm( arrPntRez[1].X, m_x0,
			m_x1) && (arrPntRez[1].Y >= -TOLRNC));
	if (b0 && (!b1))
	return 1;
	if ( b1 && (!b0) )
	{
	arrPntRez[0] = arrPntRez[1];
	return 1;
	}
	if ( (!b1) && (!b0) ) return 0;
	return 2;
	default:
	return -1;
        ;
	}

//  return irez ;
}

// пересечение сегмента arrPntSgm[2] c линией параболы
// возвращает к-во очек пересечения
// arrPntRez  - масссив точе пересечения,
int TArcParab::fncSegmentCutParabLine(TURPointXY pntSegm0, TURPointXY pntSegm1 // сегмент
							,TURPointXY *arrPntRez // точек ассив а пересечения
								)
{
	arrPntRez[0] = TURPointXY(0.,0.);
	arrPntRez[1] = TURPointXY(0.,0.);
	int irez = fncLineIntersectParab (pntSegm0, pntSegm1, arrPntRez);
	// исключение точек, не принадлежащих сегменту
	// делаем линейное преобразование координат - поворот оси X так, чтобы
	// ось X1 была бы параллельна линии сегимента
	double valAlf = 0.;
	 TURPolygon::Atan2_(pntSegm1.Y - pntSegm0.Y , pntSegm1.X - pntSegm0.X,  &valAlf ) ;
	 double arrMtrxPer[4] = {0.};
	 arrMtrxPer[0] = cos(valAlf ) ;
	 arrMtrxPer[1] = sin(valAlf ) ;
	 arrMtrxPer[2] = - arrMtrxPer[1];
	 arrMtrxPer[3] = arrMtrxPer[0] ;
	 TURPointXY pntSegmNew0,pntSegmNew1;
	 TURPointXY arrPntRezNew[2];
	 pntSegmNew0 = pntSegm0.fncLinTrasform(arrMtrxPer)  ;
	 pntSegmNew1 = pntSegm1.fncLinTrasform(arrMtrxPer)  ;
	 arrPntRezNew[0] = arrPntRez[0].fncLinTrasform(arrMtrxPer)  ;
	 arrPntRezNew[1] = arrPntRez[1].fncLinTrasform(arrMtrxPer)  ;
	bool b0 = true, b1 = true ;
	switch (irez)
	{
	 case 0:
	 return 0;

	 case 1:
	  if (!(Is_X_BelongeSegm( arrPntRezNew[0].X, pntSegmNew0.X, pntSegmNew1.X) && (arrPntRez[0].Y >= -TOLRNC) ) )
	  {
		 return 0;
	  }
	  return 1 ;

	  case 2:
	  b0 =  (Is_X_BelongeSegm( arrPntRezNew[0].X, pntSegmNew0.X, pntSegmNew1.X) && (arrPntRez[0].Y >= -TOLRNC));
	  b1 =  (Is_X_BelongeSegm( arrPntRezNew[1].X, pntSegmNew0.X, pntSegmNew1.X) && (arrPntRez[1].Y >= -TOLRNC));
	  if (b0 && (!b1))
	  return 1;
	  if ( b1 && (!b0) )
	  {
	  arrPntRez[0] = arrPntRez[1];
	  return 1;
	  }
	  if ( (!b1) && (!b0) ) return 0;
	  return 2;

	default:
	return -1;
        ;
	}
  //	return -1;
}



void TArcParab::ShowMe(wchar_t *FileName)
{
   int iN = 600;
   TURPolyLine line(*this,iN);
   line.WriteSetSHPFiles(FileName,&line, 1 ) ;

}

void TArcParab::ShowFullGraph(wchar_t *FileName)
{
 TArcParab tempPrb = *this;
 TURPointXY  arrPntSgm[2] ;
  arrPntSgm[0] =  TURPointXY(1.,0.);
  arrPntSgm[1] =  TURPointXY(-1.,0.);
  TURPointXY  arrPntRez[2] ;
 int irez = fncLineIntersectParab (arrPntSgm[0], arrPntSgm[1],arrPntRez);
  if (irez == 2)
  {
	tempPrb.m_x0 = min_(arrPntRez[0].X,arrPntRez[1].X) ;
	tempPrb.m_x1 = max_(arrPntRez[0].X,arrPntRez[1].X) ;
	tempPrb.ShowMe(FileName) ;
  }


}

 bool Is_X_BelongeSegm( const double valx, const double vala,const double valb)
{
	bool ireturn = false;
	if ((valx <= (valb + TOLRNC) )&& ( (vala - TOLRNC) <= valx)) ireturn = true ;

	return ireturn ;
}

double max_( const double x0, const double x1)
{
	return (x0 >= x1)? x0:x1;
}

double min_( const double x0, const double x1)
{
	return (x0 >= x1)? x1:x0;
}

#pragma package(smart_init)
