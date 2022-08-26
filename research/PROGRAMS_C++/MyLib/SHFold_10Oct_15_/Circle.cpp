//---------------------------------------------------------------------------


#pragma hdrstop
#include <vcl.h>


#include <math.h>
#include "Circle.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "URPolygon.h"
#include "MatrixProccess.h"
//---------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
 TCircle ::TCircle()
{
	mPntCentre = TURPointXY(0., 0.);
	mR = 3. ;

}

 // парам констр
TCircle :: TCircle( const TURPointXY PntCentre, const double R )

 {
   mPntCentre = PntCentre ;
   mR = R ;

 }





 // оператор присваивания
 TCircle TCircle::operator=(TCircle  R)
 {
	mPntCentre = R.mPntCentre ;
	mR = R.mR ;

	return *this ;
 }

 // конструктор копирования
 TCircle::TCircle (const TCircle &R)
 {
	mPntCentre = R.mPntCentre ;
	mR = R.mR ;

 }

void TCircle::ShowMe(wchar_t *FileName)
{
   int iN = 600;
   TURPolyLine line(*this,iN);
   line.WriteSetSHPFiles(FileName,&line, 1 ) ;
}

// пересечение параматрической линии и окружности
// возвращает 0 если не пересекаются
// 1 если касаются
// 2 если 2 точки пересечения
// arrPntRez- массив точек пересечения
int TCircle::fncParametricLineCutCircle(TURPointXY pointLine, double *arrVectLine // линия
							,TURPointXY *arrPntRez // точек ассив а пересечения
								)
{
	if (fabs(arrVectLine[0]) < TOLRNC)
	{
	  return fncParametricVerticalLineCutCircle( pointLine, arrVectLine // линия
							,arrPntRez // точек ассив а пересечения
								);
	}
	else
	{
		// сведение случая к случаю вертикальной прямой путем поворота СК на ПИ/2
	  double valFi = - M_PI / 2. + atan( arrVectLine[1]/ arrVectLine[0]) ;
	  double arrMtxPer[4] = {0.};
	  arrMtxPer[0] =  cos (valFi ) ;
	  arrMtxPer[1] =  sin (valFi ) ;
	  arrMtxPer[2] = -sin (valFi ) ;
	  arrMtxPer[3] =  cos (valFi ) ;

	  TURPointXY pointLine_new, pntCrclCentre_new, arrPntRez_new[2];
	  double arrVectLine_new[2] = {0.} ;
	  pointLine_new = pointLine.fncLinTrasform(arrMtxPer);
	  TCircle CircleTemp = *this;
	  CircleTemp.mPntCentre = mPntCentre.fncLinTrasform(arrMtxPer);
	  MtrxMultMatrx(arrMtxPer,2, 2, arrVectLine,1, arrVectLine_new) ;

	  int ireturn = CircleTemp.fncParametricVerticalLineCutCircle(   pointLine_new, arrVectLine_new // линия
							,arrPntRez_new // точек ассив а пересечения
								);
	 if (ireturn > 0)
	 {
	  arrMtxPer[1] = -arrMtxPer[1] ;
	  arrMtxPer[2] = -arrMtxPer[2] ;
	  arrPntRez[0] =  arrPntRez_new[0].fncLinTrasform(arrMtxPer);
	  arrPntRez[1] =  arrPntRez_new[1].fncLinTrasform(arrMtxPer);

	 }
	 return ireturn;


	}

}

// пересечение параматрической вертикальной линии и окружности
// возвращает 0 если не пересекаются
// 1 если касаются
// 2 если 2 точки пересечения
// arrPntRez- массив точек пересечения
int TCircle::fncParametricVerticalLineCutCircle(TURPointXY pointLine, double *arrVectLine // линия
							,TURPointXY *arrPntRez // точек массив а пересечения
							)
{
	if (fabs(arrVectLine[0]) > TOLRNC)
	{
		ShowMessage(L"Line is not vertical. TURPolygon::fncParametricVerticalLineCutCircle") ;
		return 0 ;
	}
	double valDiscr =  mR * mR - (pointLine.X - mPntCentre.X) * (pointLine.X - mPntCentre.X);
	if (fabs (valDiscr) < TOLRNC * TOLRNC)
	{
	  arrPntRez[0].X =  pointLine.X;
	  arrPntRez[0].Y =  mPntCentre.Y ;
	  return 1;
	}
	if (valDiscr > 0.)
	{
	  arrPntRez[0].X =  pointLine.X;
	  arrPntRez[0].Y =  mPntCentre.Y - sqrt( valDiscr);
	  arrPntRez[1].X =  pointLine.X;
	  arrPntRez[1].Y =  mPntCentre.Y + sqrt( valDiscr);
	  return 2;
	}
  return -2;
}

#pragma package(smart_init)
