//---------------------------------------------------------------------------


#pragma hdrstop
#include <stdio.h>
#include <vcl.h>
#include <math.h>
#include <dir.h>
#include "OpticParabAnt.h"
#include "RectWaveGuide.h"
#include "ArcParab.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#//include "URPolygon.h"
#include "AntFar.h"
#include "Equations.h"
#include "Comp.h"
#include "YrWriteShapeFile.h"

//#include "unit2.h"

//---------------------------------------------------------------------------

__fastcall TOpticParabAnt::~TOpticParabAnt()
{
	if(mpArrWaveGuide) //delete [] mpArrWaveGuide ;
	mpArrWaveGuide = NULL ;

}
//---------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
 TOpticParabAnt ::TOpticParabAnt()
{
	// фокальный параметр зеркала антенны (ml/2 = фокусное расстояние)
	 ml =0.;
	// диаметр зеркала антенны
	mParabDiam=0.;
	// к-во волноводов
	mNumWaveGuide=0.;
	// расстояние между центрами волноводов
	mWaveGuideDist=0.;
	//высота излучателя
	mWaveGuideHeight=0.;
	// массив волноводов
	mpArrWaveGuide = NULL;
}

 // парам констр
TOpticParabAnt :: TOpticParabAnt(const double l,const double ParabDiam,const int NumWaveGuide
  ,const double WaveGuideDist,const double WaveGuideHeight,const TRectWaveGuide WaveGuide)

 {
	// фокальный параметр зеркала антенны (ml/2 = фокусное расстояние)
	ml =l ;
	// диаметр зеркала антенны
	mParabDiam = ParabDiam;
	// к-во волноводов
	mNumWaveGuide = NumWaveGuide;
	// расстояние между центрами волноводов
	mWaveGuideDist = WaveGuideDist;
	//высота излучателя
	mWaveGuideHeight= WaveGuideHeight;
	// массив волноводов
	mpArrWaveGuide = NULL;

	mpArrWaveGuide = new TRectWaveGuide[mNumWaveGuide ];
	if(mpArrWaveGuide== NULL)
	{
	   ShowMessage(L"Not memory for mpArrWaveGuide") ;
	   Abort() ;
	}
	for (int i =0; i < mNumWaveGuide; i++)
	{
	mpArrWaveGuide[i] = WaveGuide;
	}

 }

  // парам констр
TOpticParabAnt :: TOpticParabAnt(const double l,const double ParabDiam,const int NumWaveGuide
  ,const double WaveGuideDist,const double WaveGuideHeight,const double b)

 {
	// фокальный параметр зеркала антенны (ml/2 = фокусное расстояние)
	ml =l;
	// диаметр зеркала антенны
	mParabDiam = ParabDiam;
	// к-во волноводов
	mNumWaveGuide = NumWaveGuide;
	// расстояние между центрами волноводов
	mWaveGuideDist = WaveGuideDist;
	//высота излучателя
	mWaveGuideHeight= WaveGuideHeight;
	// массив волноводов
	mpArrWaveGuide = NULL;

	mpArrWaveGuide = new TRectWaveGuide[mNumWaveGuide ];
	if(mpArrWaveGuide== NULL)
	{
	   ShowMessage(L"Not memory for mpArrWaveGuide") ;
	   Abort() ;
	}
	TRectWaveGuide WaveGuide(4., b);
	for (int i =0; i < mNumWaveGuide; i++)
	{
	mpArrWaveGuide[i] = WaveGuide;
	}
 }

 // оператор присваивания
 TOpticParabAnt TOpticParabAnt::operator=(TOpticParabAnt  R)
 {
	// фокальный параметр зеркала антенны (ml/2 = R.mфокусное расстояние)
	ml = R.ml;
	// диаметр зеркала антенны
	mParabDiam = R.mParabDiam;
	// к-во волноводов
	mNumWaveGuide = R.mNumWaveGuide;
	// расстояние между центрами волноводов
	mWaveGuideDist = R.mWaveGuideDist;
	//высота излучателя
	mWaveGuideHeight= R.mWaveGuideHeight;

	if (mpArrWaveGuide != NULL) mpArrWaveGuide= NULL;
	if(R.mpArrWaveGuide != NULL)
	{
	mpArrWaveGuide = new TRectWaveGuide[R.mNumWaveGuide];
	if(mpArrWaveGuide == NULL) 	   {
	ShowMessage(L"Not memory for mpArrWaveGuide") ;
	Abort() ;
	}
	memcpy( mpArrWaveGuide,R.mpArrWaveGuide, R.mNumWaveGuide  * sizeof(TRectWaveGuide));
	}
	return *this ;
 }

 // конструктор копирования
 TOpticParabAnt::TOpticParabAnt (const TOpticParabAnt &R)
 {
	// фокальный параметр зеркала антенны (ml/2 = R.mфокусное расстояние)
	ml = R.ml;
	// диаметр зеркала антенны
	mParabDiam = R.mParabDiam;
	// к-во волноводов
	mNumWaveGuide = R.mNumWaveGuide;
	// расстояние между центрами волноводов
	mWaveGuideDist = R.mWaveGuideDist;
	//высота излучателя
	mWaveGuideHeight= R.mWaveGuideHeight;

	if (mpArrWaveGuide != NULL) mpArrWaveGuide= NULL;
	if(R.mpArrWaveGuide != NULL)
	{
	mpArrWaveGuide = new TRectWaveGuide[R.mNumWaveGuide];
	if(mpArrWaveGuide == NULL) 	   {
	ShowMessage(L"Not memory for mpArrWaveGuide") ;
	Abort() ;
	}
	memcpy( mpArrWaveGuide,R.mpArrWaveGuide, R.mNumWaveGuide  * sizeof(TRectWaveGuide));
	}
}

bool TOpticParabAnt::createPictWithRay(wchar_t *wchFoldPict, const double VAlContrReflAng
  , const double VAlTargAng, const double VAlYf)
{
  double  valXf =0.;
  if (VAlTargAng > 0.)
  {
	 valXf =  - mParabDiam/2.* sin(VAlContrReflAng)
	   + ( mParabDiam/2.* cos(VAlContrReflAng) -VAlYf) * tan(VAlContrReflAng + VAlTargAng);
  }
  else
  {
   valXf =  mParabDiam/2.* sin(VAlContrReflAng)
       - (  VAlYf +  mParabDiam/2.* cos(VAlContrReflAng) ) * tan(VAlContrReflAng + VAlTargAng);
  }
  const double  VAlXf = valXf;
  //- mParabDiam/2.* sin(VAlContrReflAng)
	//   + ( mParabDiam/2.* cos(VAlContrReflAng) -VAlYf) * tan(VAlContrReflAng + VAlTargAng);
  // парабола
  TArcParab parab(0., ml);
  parab.m_x0 = - mParabDiam/2.;
  parab.m_x1 =  mParabDiam/2.;

	TURPolyLine plnParab(  parab,1000);


//  DrawShapeForm->DSH+= plnParab;


	TURPointXY pntSdvig(0.,0.);
	plnParab = plnParab.LinTransform(-M_PI/2. ,  pntSdvig,1.) ;
	 wchar_t FileName0[400] = {0};
	wcscpy( FileName0,  wchFoldPict);
	wcscat(FileName0, L"\\plnParab.shp");
	plnParab.WriteSetSHPFiles(FileName0,&plnParab, 1) ;
	///

	// фокус
	TURPointXY pntFocus(0.,0.);
	wcscpy( FileName0,  wchFoldPict);
	wcscat(FileName0, L"\\pntFocus.shp");
	pntFocus.WriteSetSHPFiles(FileName0,&pntFocus, 1) ;
	///

	// волноводы
	// INPUT:
// wchFoldName - путь к папке с файлами
// Emitter:
// ValLength0 - длина горизонтальной части эмиттера
// ValLength1 - длина наклонной части эмиттера
// ValFi - угол между горизонтальной и наклонной частями эмиттера
//  AM:
// QuantEmit - к-во эмиттеров
// DistEmit - расстояние мжде эмиттерами
// ValLengthBoardAM - расстояние между границами АМ
// SideWidthAM - толщина боковых частей полигона кожуха эмиттеров
// BackWidthAM - толщина задней части  полигона кожуха эмиттеров
// LenghPlnOutAM  - длина линии выхода с АМ
// FAR:
// QuantAM - к-во АМ
//
wchar_t wchFoldName[400] = {0};
wcscpy( wchFoldName,  wchFoldPict);
wcscat(wchFoldName, L"\\WaveGuids");
_wmkdir(wchFoldName);
double DistEmit = mWaveGuideHeight;
const double ValLengthBoardAM = (mWaveGuideDist - mWaveGuideHeight)/2.;
const double SideWidthAM = 0.01;
const double BackWidthAM = 0.01;
const double LenghPlnOutAM = 10.;
const double ValLength0 = 2. * 0.1;
const double ValFi =  M_PI/ 6.;
const double ValLength1 =  1./cos(ValFi) * 0.8 * 0.1;

/*
const double ValLengthBoardAM =  0.8;
	 const double SideWidthAM = 0.4;
	 const double BackWidthAM = 3.;
	 const double LenghPlnOutAM = 20.;
	 const double ValLength0 = 2.;
	 const double ValFi = M_PI/ 6.;
	 const double ValLength1 = 1./cos(ValFi) * 0.8;
	 */
 createFarPict(wchFoldName, 1, DistEmit
  , mNumWaveGuide,  ValLengthBoardAM
  , SideWidthAM,  BackWidthAM , LenghPlnOutAM
  , ValLength0 , ValLength1, ValFi) ;
  ///

  // линия волноводов
  double valLenWaveGuides =  mWaveGuideHeight + ((double)mNumWaveGuide -1.) * mWaveGuideDist;
   TURPointXY  pnt1(0., -valLenWaveGuides/2.);
   TURPointXY  pnt2 (0., valLenWaveGuides/2.);
	TURPolyLine plnWaveGuide(   pnt1,   pnt2) ;

//	DrawShapeForm->DSH+= plnWaveGuide;

	wcscpy( FileName0,  wchFoldName);
	wcscat(FileName0, L"\\plnWaveGuide.shp");
	plnWaveGuide.WriteSetSHPFiles(FileName0,&plnWaveGuide, 1) ;
	///

	// линия контррефлектора
	TURPointXY pntarr[4];
	pntarr[0] =  TURPointXY(0.,- mParabDiam/2.);
	pntarr[1] =  TURPointXY(0.,- valLenWaveGuides/2.);
	pntarr[2] =  TURPointXY(0., valLenWaveGuides/2.);
	pntarr[3] =  TURPointXY(0.,mParabDiam/2.);
	int iarrParts[2] = {0,2};
	TURPolyLine plnContrRf(  2,4,iarrParts	,pntarr)  ;

//	DrawShapeForm->DSH+= plnContrRf;

	plnContrRf = plnContrRf.LinTransform(VAlContrReflAng ,  pntSdvig,1.) ;
	wcscpy( FileName0,  wchFoldName);
	wcscat(FileName0, L"\\plnContrRf.shp");
	plnContrRf.WriteSetSHPFiles(FileName0,&plnContrRf, 1) ;
	///

	// фронт волны в момент 0
	double valXfr00 = 0.,valYfr00 = 0. , valXfr01 = 0.,valYfr01 = 0.;
	TURPointXY pntFront00(valXfr00 , valYfr00 );
	TURPointXY pntFront01(valXfr00 , valYfr00 );
	if (VAlTargAng >0.)
	{
	valXfr00 = -mParabDiam/2. * sin( VAlContrReflAng);
	valYfr00 = mParabDiam/2. * cos( VAlContrReflAng);
	pntFront00 = TURPointXY(valXfr00 , valYfr00 );
	valXfr01 =valXfr00 +  mParabDiam *  cos(VAlTargAng ) * sin(VAlTargAng +  VAlContrReflAng);
	valYfr01 =valYfr00 -  mParabDiam *  cos(VAlTargAng ) * cos(VAlTargAng +  VAlContrReflAng);
	pntFront01 = TURPointXY(valXfr01 , valYfr01 );
	}
	else
	{
	valXfr00 = mParabDiam/2. * sin( VAlContrReflAng);
	valYfr00 = -mParabDiam/2. * cos( VAlContrReflAng);
	pntFront00 = TURPointXY(valXfr00 , valYfr00 );
	valXfr01 =valXfr00 -  mParabDiam *  cos(VAlTargAng ) * sin(VAlTargAng +  VAlContrReflAng);
	valYfr01 =valYfr00 + mParabDiam *  cos(VAlTargAng ) * cos(VAlTargAng +  VAlContrReflAng);
	pntFront01 = TURPointXY(valXfr01 , valYfr01 );
	}

	TURPolyLine plnFront0(   pntFront00,   pntFront01) ;
	wcscpy( FileName0,  wchFoldName);
	wcscat(FileName0, L"\\plnFront0.shp");
	plnFront0.WriteSetSHPFiles(FileName0,&plnFront0, 1) ;

//	DrawShapeForm->DSH+= plnFront0;

	///

	// линия падающего луча
	 TURPointXY pntarrRay[5] ;  // массив точек хода луча
		 // начальная точка луча
	 //  double valYf =valYfr00 -  (VAlXf - valXfr00)/ tan( VAlTargAng +  VAlContrReflAng);
	//  double valYf = VAlYf;
		pntarrRay[0] =  TURPointXY(VAlXf, VAlYf);
		 // точка тотражение от контррефлектора
	 //  double valXOKRF = (VAlXf * tan(VAlTargAng +  VAlContrReflAng) - VAlYf)
	 //  / (tan(VAlTargAng +  VAlContrReflAng ) + 1./ tan(  VAlContrReflAng ));
	//  double valYOKRF =  -valXOKRF / tan(  VAlContrReflAng ) ;
	double valYOKRF =  (VAlYf - VAlXf * tan(VAlTargAng +  VAlContrReflAng ))
		 / (1. + tan( VAlContrReflAng) * tan(VAlTargAng +  VAlContrReflAng));
		 double valXOKRF = - valYOKRF * tan( VAlContrReflAng);
		pntarrRay[1] =  TURPointXY (valXOKRF,valYOKRF);  // точка отражениея от контррефлектора
		// отрезок нормали к контррефлектору
		 TURPointXY pntNorm1(valXOKRF + 10., valYOKRF + tan(VAlContrReflAng) * 10.);
		 TURPolyLine plnContrReflNorm(pntarrRay[1] , pntNorm1) ;
		 wcscpy( FileName0,  wchFoldName);
		wcscat(FileName0, L"\\ContrReflNorm.shp");
		plnContrReflNorm.WriteSetSHPFiles(FileName0,&plnContrReflNorm, 1) ;

//		DrawShapeForm->DSH+= plnContrReflNorm;

				///
	  // точка отражения от  параболы
	   // находится путем решения квадратного уравнения
		 // коэффиц квадратного уравнения a,b,c
		 double vala = tan( VAlContrReflAng - VAlTargAng);
		 double valb = 2. * ml ;
		 double valc = -2. * ml * valYOKRF  + 2. * ml * valXOKRF * tan( VAlContrReflAng - VAlTargAng)
			- ml* ml * tan( VAlContrReflAng - VAlTargAng);
		 TComp y1, y2 ;
		 int irez = SolvEq2(vala,valb,valc, y1, y2);
		 if ((irez == 3)||(irez == 5)||(irez == 6))
		 {
			//ShowMessage(L"EROOR");
			return false;
		 }
		 double valYOP = y1.m_Re;
		 if (fabs(valYOP) > mParabDiam/2. )
		 {
			valYOP = y2.m_Re;
			if (fabs(valYOP) > mParabDiam/2. )
			{
		   //	ShowMessage(L"ERROR1");
			return false;
			}
		 }
		 double valXOP = -valYOP * valYOP/ 2./ ml + ml /2.;
		 pntarrRay[2] =   TURPointXY (valXOP,valYOP);  // точка отражениея от параболы

		 /*// точка отражения от  параболы
	   // находится путем решения квадратного уравнения
		 // коэффиц квадратного уравнения a,b,c
		 double vala = tan( VAlContrReflAng - VAlTargAng);
		 double valb = 2. * ml ;
		 double valc = -2. * ml * valYOKRF  + 2. * ml * valXOKRF * tan( VAlContrReflAng - VAlTargAng)
			- ml* ml * tan( VAlContrReflAng - VAlTargAng);
		 TComp y1, y2 ;
		 int irez = SolvEq2(vala,valb,valc, y1, y2);
		 if (irez > 2)
		 {
		   //	ShowMessage(L"EROOR");
			return false;
		 }
		 double valYOP = y1.m_Re;
		 if (fabs(valYOP) > mParabDiam/2. )
		 {
			valYOP = y2.m_Re;
			if (fabs(valYOP) > mParabDiam/2. )
			{
		  //	ShowMessage(L"ERROR1");
			 return false;
			}
		 }
		 double valXOP = -valYOP * valYOP/ 2./ ml + ml /2.;
		 //( valYOP - valYOKRF )/ tan( VAlContrReflAng - VAlTargAng)  + valXOKRF;
		 pntarrRay[2] =   TURPointXY (valXOP,valYOP);  // точка отражениея от параболы
		 ///
		 */
		   // отрезок нормали к параболе
	   TURPointXY pntNorm2(valXOP - 2.* ml, valYOP  - valYOP * 2.);
		 TURPolyLine plnParabReflNorm(pntarrRay[2] , pntNorm2) ;
	   wcscpy( FileName0,  wchFoldName);
	   wcscat(FileName0, L"\\plnParabReflNorm.shp");
		 plnParabReflNorm.WriteSetSHPFiles(FileName0,&plnParabReflNorm, 1) ;

//		 DrawShapeForm->DSH+= plnParabReflNorm;


		///
		 // угол наклона нормали к параболе в точке касания
		 double valBet = atan(valYOP / ml);

		 // точка пересечения с фокальной плоскостью
		 double valXFoc = 0.;
		 double valYFoc = valYOP + tan(VAlContrReflAng - VAlTargAng - 2. * valBet)* valXOP;
		 pntarrRay[3] =  TURPointXY (valXFoc,valYFoc);  // точка пресмечения с фокальной плоскостью
		///

		TURPolyLine plnRayWay( pntarrRay	,4)  ;
		wcscpy( FileName0,  wchFoldName);
		wcscat(FileName0, L"\\plnRayWay.shp");
		plnRayWay.WriteSetSHPFiles(FileName0,&plnRayWay, 1) ;

//		DrawShapeForm->DSH+= plnRayWay;


	  return true;
}

bool  TOpticParabAnt::BuildRayWayPln(const double VAlContrReflAng
  , const double VAlTargAng, const double VAlYf, TURPolyLine &plnRez)
{

	double  valXf =0.;
  if (VAlTargAng > 0.)
  {
	 valXf =  - mParabDiam/2.* sin(VAlContrReflAng)
	   + ( mParabDiam/2.* cos(VAlContrReflAng) -VAlYf) * tan(VAlContrReflAng + VAlTargAng);
  }
  else
  {
   valXf =  mParabDiam/2.* sin(VAlContrReflAng)
       - (  VAlYf +  mParabDiam/2.* cos(VAlContrReflAng) ) * tan(VAlContrReflAng + VAlTargAng);
  }
  const double  VAlXf = valXf;

	// фронт волны в момент 0
	double valXfr00 = 0.,valYfr00 = 0. , valXfr01 = 0.,valYfr01 = 0.;
	TURPointXY pntFront00(valXfr00 , valYfr00 );
	TURPointXY pntFront01(valXfr00 , valYfr00 );
	if (VAlTargAng >0.)
	{
	valXfr00 = -mParabDiam/2. * sin( VAlContrReflAng);
	valYfr00 = mParabDiam/2. * cos( VAlContrReflAng);
	pntFront00 = TURPointXY(valXfr00 , valYfr00 );
	valXfr01 =valXfr00 +  mParabDiam *  cos(VAlTargAng ) * sin(VAlTargAng +  VAlContrReflAng);
	valYfr01 =valYfr00 -  mParabDiam *  cos(VAlTargAng ) * cos(VAlTargAng +  VAlContrReflAng);
	pntFront01 = TURPointXY(valXfr01 , valYfr01 );
	}
	else
	{
	valXfr00 = mParabDiam/2. * sin( VAlContrReflAng);
	valYfr00 = -mParabDiam/2. * cos( VAlContrReflAng);
	pntFront00 = TURPointXY(valXfr00 , valYfr00 );
	valXfr01 =valXfr00 -  mParabDiam *  cos(VAlTargAng ) * sin(VAlTargAng +  VAlContrReflAng);
	valYfr01 =valYfr00 + mParabDiam *  cos(VAlTargAng ) * cos(VAlTargAng +  VAlContrReflAng);
	pntFront01 = TURPointXY(valXfr01 , valYfr01 );
	}


	///

	// линия падающего луча
	 TURPointXY pntarrRay[5] ;  // массив точек хода луча
	   // начальная точка луча
	  pntarrRay[0] =  TURPointXY(VAlXf, VAlYf);

	   // точка отражение от контррефлектора

	double valYOKRF =  (VAlYf - VAlXf * tan(VAlTargAng +  VAlContrReflAng ))
	   / (1. + tan( VAlContrReflAng) * tan(VAlTargAng +  VAlContrReflAng));
	   double valXOKRF = - valYOKRF * tan( VAlContrReflAng);
	  pntarrRay[1] =  TURPointXY (valXOKRF,valYOKRF);  // точка отражениея от контррефлектора

	  // точка отражения от  параболы
	   // находится путем решения квадратного уравнения
		 // коэффиц квадратного уравнения a,b,c
		 double vala = tan( VAlContrReflAng - VAlTargAng);
		 double valb = 2. * ml ;
		 double valc = -2. * ml * valYOKRF  + 2. * ml * valXOKRF * tan( VAlContrReflAng - VAlTargAng)
			- ml* ml * tan( VAlContrReflAng - VAlTargAng);
		 TComp y1, y2 ;
		 int irez = SolvEq2(vala,valb,valc, y1, y2);
		 if ((irez == 3)||(irez == 5)||(irez == 6))
		 {
		   //	ShowMessage(L"EROOR");
			return false;
		 }
		 double valYOP = y1.m_Re;
		 if (fabs(valYOP) > mParabDiam/2. )
		 {
			valYOP = y2.m_Re;
			if (fabs(valYOP) > mParabDiam/2. )
			{
		  //	ShowMessage(L"ERROR1");
			 return false;
			}
		 }
		 double valXOP = -valYOP * valYOP/ 2./ ml + ml /2.;
		 //( valYOP - valYOKRF )/ tan( VAlContrReflAng - VAlTargAng)  + valXOKRF;
		 pntarrRay[2] =   TURPointXY (valXOP,valYOP);  // точка отражениея от параболы
		 ///

		 // угол наклона нормали к параболе в точке касания
	   double valBet = atan(valYOP / ml);

		 // точка пересечения с фокальной плоскостью
	   double valXFoc = 0.;
	   double valYFoc = valYOP + tan(VAlContrReflAng - VAlTargAng - 2. * valBet)* valXOP;
		 pntarrRay[3] =  TURPointXY (valXFoc,valYFoc);  // точка пресмечения с фокальной плоскостью
	  ///

	  TURPolyLine plnRayWay( pntarrRay	,4)  ;
	  plnRez = plnRayWay;
	  return true;
}


bool  TOpticParabAnt::createPhaseFrontPict(wchar_t *wchFoldPict, const double VAlContrReflAng
  , const double VAlTargAng , double valY0, double valY1)
{
   int iMaxRaysQuant = 10000;
   double valStep = mParabDiam / ((double)iMaxRaysQuant);
  // вычисление необходимого к-ва полилиний
  int quantRays = calcRaysQuant (VAlContrReflAng, VAlTargAng, valY0, valY1 , valStep ) ;
  ///
   if (quantRays == -1)
	{
	return false;
	}
  //
  TURPolyLine *plnarr = new   TURPolyLine [quantRays];
  createRaysArray( VAlContrReflAng , VAlTargAng
	,valY0, valY1, valStep, quantRays, plnarr) ;
	  // , -totalLengthWaveGuide()/2., totalLengthWaveGuide()/2., valStep, quantRays, plnarr) ;

	wchar_t FileName0[400] = {0};
	wcscpy( FileName0,  wchFoldPict);
	wcscat(FileName0, L"\\PhaseFrontLines.shp");
	TURPolyLine::WriteSetSHPFiles(FileName0,plnarr, quantRays) ;

	// вычисление кратчайшего пути
	double valMin = 100000000.;
	int iArg = -1;
	for (int i = 0; i < quantRays; i++)
	{
	  double dist =  plnarr[i].calcLeng();
	  if (dist <valMin )
	  {
	   valMin = dist;
	   iArg = i;
	  }
	}
	TURPointXY *pntarrFront = new TURPointXY[quantRays ];
    for (int i = 0; i < quantRays; i++)
	{
	   plnarr[i].calcWayPoint( valMin, pntarrFront[i]) ;
	}
	wcscpy( FileName0,  wchFoldPict);
	wcscat(FileName0, L"\\PhaseFrontPoints.shp");
	TURPointXY::WriteSetSHPFiles(FileName0,pntarrFront, quantRays) ;
  delete []plnarr ;
  delete []pntarrFront;


}

int   TOpticParabAnt::calcRaysQuant (const double VAlContrReflAng , const double VAlTargAng
  , double valY0, double valY1,double valStep)
{
	// нахождение полилинии фаронта в момент 0
	TURPolyLine plnFront0 = calcFront0 (VAlContrReflAng , VAlTargAng);
	///
	int iC = mParabDiam / valStep;
	int irez = -1;
	TURPolyLine plnRez;
	 for (int i = 0; i < iC; i++)
	 {
		// вычисление координаты  Y  начала i -го луча
		double valYf = plnFront0.Points[0].Y + ((double)i)/((double)iC)* (plnFront0.Points[1].Y -plnFront0.Points[0].Y);
		if (!BuildRayWayPln( VAlContrReflAng,  VAlTargAng,  valYf, plnRez))
		{
		 continue;
		}
		if ((plnRez.Points[3].Y >  valY1) ||(plnRez.Points[3].Y <  valY0 ))
		{
		   continue;
		}
		irez++;
	 }
	return irez;
}

void TOpticParabAnt::createRaysArray(const double VAlContrReflAng , const double VAlTargAng
  , double valY0, double valY1,double valStep, int quantRays, TURPolyLine *plnarr)
{
 	// нахождение полилинии фаронта в момент 0
	TURPolyLine plnFront0 = calcFront0 (VAlContrReflAng , VAlTargAng);
	///
	int iC = mParabDiam / valStep;

	TURPolyLine plnRez;
	int iCur = 0;
	 for (int i = 0; i < iC; i++)
	 {
		// вычисление координаты  Y  начала i -го луча
		double valYf = plnFront0.Points[0].Y + ((double)i)/((double)iC)* (plnFront0.Points[1].Y -plnFront0.Points[0].Y);
		if (!BuildRayWayPln( VAlContrReflAng,  VAlTargAng,  valYf, plnRez))
		{
		 continue;
		}
		if ((plnRez.Points[3].Y >  valY1) ||(plnRez.Points[3].Y <  valY0 ))
		{
		   continue;
		}
		plnarr [iCur] = plnRez;
		iCur++;
	 }
}

 // полилиния фронта волны в момент 0
TURPolyLine TOpticParabAnt::calcFront0 (const double VAlContrReflAng , const double VAlTargAng)
{
	// фронт волны в момент 0
	double valXfr00 = 0.,valYfr00 = 0. , valXfr01 = 0.,valYfr01 = 0.;
	TURPointXY pntFront00(valXfr00 , valYfr00 );
	TURPointXY pntFront01(valXfr00 , valYfr00 );
	if (VAlTargAng >0.)
	{
	valXfr00 = -mParabDiam/2. * sin( VAlContrReflAng);
	valYfr00 = mParabDiam/2. * cos( VAlContrReflAng);
	pntFront00 = TURPointXY(valXfr00 , valYfr00 );
	//valXfr01 =valXfr00 +  mParabDiam *  cos(VAlTargAng ) * sin(VAlTargAng +  VAlContrReflAng);
	//valYfr01 =valYfr00 -  mParabDiam *  cos(VAlTargAng ) * cos(VAlTargAng +  VAlContrReflAng);
	valXfr01 =valXfr00 +  mParabDiam  * sin(VAlTargAng +  VAlContrReflAng);
	valYfr01 =valYfr00 -  mParabDiam * cos(VAlTargAng +  VAlContrReflAng);
	pntFront01 = TURPointXY(valXfr01 , valYfr01 );
	}
	else
	{
	valXfr00 = mParabDiam/2. * sin( VAlContrReflAng);
	valYfr00 = -mParabDiam/2. * cos( VAlContrReflAng);
	pntFront00 = TURPointXY(valXfr00 , valYfr00 );
	//valXfr01 =valXfr00 -  mParabDiam *  cos(VAlTargAng ) * sin(VAlTargAng +  VAlContrReflAng);
   //	valYfr01 =valYfr00 + mParabDiam *  cos(VAlTargAng ) * cos(VAlTargAng +  VAlContrReflAng);
	valXfr01 =valXfr00  - mParabDiam  * sin(VAlTargAng +  VAlContrReflAng);
	valYfr01 =valYfr00  + mParabDiam  * cos(VAlTargAng +  VAlContrReflAng);
	pntFront01 = TURPointXY(valXfr01 , valYfr01 );
	}

	TURPolyLine plnFront0(   pntFront00,   pntFront01) ;
	return plnFront0;
}

double TOpticParabAnt::totalLengthWaveGuide()
{
	return ((double)mNumWaveGuide -1.) *mWaveGuideDist + mWaveGuideHeight;
}

// постронеие массива диаграмм
// ILenDiagrArr - к-во точек по каждой диаграмме
//
// OUTPUT:
//cmparrVeerDiagrs[ ILenDiagrArr * mNumWaveGuide] - массив диаграмм
// хранятся по очереди, начиная с ниждней(нижний волновод)
void  TOpticParabAnt::calcVeerDiagrams(wchar_t *wchFoldPict, const double  VAlLambda,  const double VAlContrReflAng
	 , const int ILenDiagrArr, TComp *cmparrVeerDiagrs )
{
  double valMaxTargAng = 0.135;//80./180 * M_PI;
  double valTargAngStep =  2. * valMaxTargAng /((double)ILenDiagrArr);
  TComp *pcmparrDiagrCur =  new TComp[mNumWaveGuide];
  double *parrTargAng = new double [mNumWaveGuide];
  for (int i =0; i < ILenDiagrArr; i++)
  {
	if (i == 1960)
	{
      int iii=0;
	}
	double valTargAngCur =  -valMaxTargAng + ((double)i) * valTargAngStep ;
	calcDiagrArrayFromDirect(VAlLambda, VAlContrReflAng,valTargAngCur, pcmparrDiagrCur );
	//calcDiagrArrayFromDirect(VAlLambda, VAlContrReflAng,0., pcmparrDiagrCur );
	for (int j =0; j < mNumWaveGuide; j++)
	{
	  cmparrVeerDiagrs[j * ILenDiagrArr + i] =   pcmparrDiagrCur[j];
	}
	parrTargAng[i] = valTargAngCur* 1000.;

  }
  wchar_t Fold[400] ={0.};
   wcscpy( Fold, wchFoldPict);
   wcscat(Fold, L"\\Diagrams");
	_wmkdir(Fold);


  for (int i=0; i < mNumWaveGuide; i++)
  {
	wchar_t FileName[400] ={0.};
	wcscpy( FileName, Fold);
	wcscat(FileName, L"\\ModulDiagrNo_");
	wchar_t string[10] = {0.};
	_itow(i,string, 10);
	wcscat(FileName, string);
	wcscat(FileName, L".shp");

	wchar_t FileName1[400] ={0.};
	wcscpy( FileName1, Fold);
	wcscat(FileName1, L"\\ArgDiagrNo_");
	wchar_t string1[10] = {0.};
	_itow(i,string1, 10);
	wcscat(FileName1, string1);
	wcscat(FileName1, L".shp");

	double *parrModul = new double [ILenDiagrArr];
	double *parrArg = new double [ILenDiagrArr];
	for (int j = 0;j < ILenDiagrArr; j++)
	{
	  parrModul[j] =  cmparrVeerDiagrs[i * ILenDiagrArr + j].modul();
	  if (parrModul[j] > 0.000001)
	  {
	   parrArg[j] =  cmparrVeerDiagrs[i * ILenDiagrArr + j].phase() ;
	  }
	  else
	  {
		  parrArg[j] = 0.;
	  }
	}
	double scalex =1.,  scaley =100.;
	TYrWriteShapeFile::CreateShpFile(FileName, parrModul, parrTargAng
	 ,ILenDiagrArr, scalex, scaley);
	  scaley =1.;
	 TYrWriteShapeFile::CreateShpFile(FileName1,parrArg, parrTargAng
	 ,ILenDiagrArr, scalex, scaley);

	delete parrModul;
	delete parrArg;

  }
  delete []pcmparrDiagrCur;
  delete parrTargAng ;

    wchar_t FileName1[400] ={0.};
  	wcscpy( FileName1, Fold);
	wcscat(FileName1, L"\\Axes.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(FileName1,-1000.,1000.
	 ,-100.,100.,5.) ;
}

void  TOpticParabAnt::calcDiagrArrayFromDirect(const double  VAlLambda
		 ,const double VAlContrReflAng, const double VAlTargAngCur, TComp *pcmparrDiagrCur )
{
	memset(pcmparrDiagrCur,0, mNumWaveGuide * sizeof(TComp));
	// нахождение полилинии фаронта в момент 0
	double valTargAngCur = VAlTargAngCur;
	TURPolyLine plnFront0 = calcFront0 (VAlContrReflAng , valTargAngCur);
	///
	//int iC = mParabDiam / valStep;
	int iC = 1000;

	TURPolyLine plnRez;
	// создание массива сетментов концов волноводов
	TURPointXY *ppntarrWaveGuideCentre =  new TURPointXY [ mNumWaveGuide];
	for (int i = 0; i < mNumWaveGuide; i++)
	{
	  ppntarrWaveGuideCentre [i] = createWaveGuideCentre(i );
	}
	 for (int i = 0; i < iC; i++)
	 {
		// вычисление координаты  Y  начала i -го луча
		double valYf = plnFront0.Points[0].Y + ((double)i)/((double)iC)* (plnFront0.Points[1].Y -plnFront0.Points[0].Y);
		if (!BuildRayWayPln( VAlContrReflAng,  VAlTargAngCur,  valYf, plnRez))
		{
		 continue;
		}

		for (int j =0; j < mNumWaveGuide; j++)
		{
		  plnRez.Points[3] = ppntarrWaveGuideCentre[j] ;
		  double valR = TURPointXY::dist(plnRez.Points[3],plnRez.Points[2]);
		  double valTetta = atan( (plnRez.Points[3].Y - plnRez.Points[2].Y)/ (plnRez.Points[3].X - plnRez.Points[2].X));
		  const double fi = plnRez.calcLeng()/ VAlLambda * 2.* M_PI;
		  pcmparrDiagrCur[j] += TComp(mpArrWaveGuide[j].fncDiagr(VAlLambda, valTetta)/valR* 50./valR* 50., 0.) * exp_(fi);

		}
	 }
	 delete []ppntarrWaveGuideCentre;
	 for (int j =0; j < mNumWaveGuide; j++)
	 {
	  pcmparrDiagrCur[j] = pcmparrDiagrCur[j]/TComp( ((double)iC), 0.);
     }
}

void  TOpticParabAnt::calcDiagrArrayFromDirect__(const double  VAlLambda
		 ,const double VAlContrReflAng, const double VAlTargAngCur, TComp *pcmparrDiagrCur )
{
	memset(pcmparrDiagrCur,0, mNumWaveGuide * sizeof(TComp));
	// нахождение полилинии фаронта в момент 0
	double valTargAngCur = VAlTargAngCur;
	TURPolyLine plnFront0 = calcFront0 (VAlContrReflAng , valTargAngCur);
	///
	//int iC = mParabDiam / valStep;
	int iC = 1000;

	TURPolyLine plnRez;
	// создание массива сетментов концов волноводов
	TURPolyLine *pplnarrWaveGuideSegm =  new TURPolyLine [ mNumWaveGuide];
	for (int i = 0; i < mNumWaveGuide; i++)
	{
	  pplnarrWaveGuideSegm [i] = createWaveGuideSegm(i );
	  pplnarrWaveGuideSegm [i].Points[0].X += 0.0001;
	  pplnarrWaveGuideSegm [i].Points[1].X += 0.001;

	}
	 for (int i = 0; i < iC; i++)
	 {
		// вычисление координаты  Y  начала i -го луча
		double valYf = plnFront0.Points[0].Y + ((double)i)/((double)iC)* (plnFront0.Points[1].Y -plnFront0.Points[0].Y);
		if (!BuildRayWayPln( VAlContrReflAng,  VAlTargAngCur,  valYf, plnRez))
		{
		 continue;
		}

		for (int j =0; j < mNumWaveGuide; j++)
		{

		  if( (plnRez.Points[3].Y - pplnarrWaveGuideSegm[j].Points[0].Y )
			* (plnRez.Points[3].Y - pplnarrWaveGuideSegm[j].Points[1].Y )
			> 0.)
			{
			  continue;
			}
		  double valTetta = atan( (plnRez.Points[3].Y - plnRez.Points[2].Y)/ (plnRez.Points[3].X - plnRez.Points[2].X));
		  const double fi = plnRez.calcLeng()/ VAlLambda * 2.* M_PI;
		  pcmparrDiagrCur[j] += TComp(mpArrWaveGuide[j].fncDiagr(VAlLambda, valTetta), 0.) * exp_(fi);
		  break;
		}
	 }
	 delete []pplnarrWaveGuideSegm;
	 for (int j =0; j < mNumWaveGuide; j++)
	 {
	  pcmparrDiagrCur[j] = pcmparrDiagrCur[j]/TComp( ((double)iC), 0.);
     }
}

// создание полилинии вероикального сегмента открытого конца волновода
//INPUT:
// INumWaveGuide - номер волновода, нумерация с 0 до (mNumWaveGuide-1)
// идет снизу вверх
// OUTPUT:
// возвращает полилинию сегмента
TURPolyLine TOpticParabAnt::createWaveGuideSegm(const int INumWaveGuide  )
{
   double valCentrePos0 = -((double)mNumWaveGuide -1.) * mWaveGuideDist/ 2.;
   double valButtomPos  = valCentrePos0 - mWaveGuideHeight/ 2.;
   double valUpperPos  = valCentrePos0 + mWaveGuideHeight/ 2.;
   TURPointXY p0(0., valButtomPos +  ((double)INumWaveGuide) * mWaveGuideDist);
   TURPointXY p1(0., valButtomPos +  ((double)INumWaveGuide + 1.) * mWaveGuideDist);
	 TURPolyLine plnReturn(p0,p1);
   return plnReturn;
}

// создание точки центра сегмента открытого конца волновода
//INPUT:
// INumWaveGuide - номер волновода, нумерация с 0 до (mNumWaveGuide-1)
// идет снизу вверх
// OUTPUT:
// возвращает точку центра сегмента
TURPointXY TOpticParabAnt::createWaveGuideCentre(const int INumWaveGuide  )
{
   double valCentrePos0 = -((double)mNumWaveGuide -1.) * mWaveGuideDist/ 2.;

   TURPointXY p0(0., valCentrePos0 +  ((double)INumWaveGuide) * mWaveGuideDist);

   return p0;
}


// постронеие массива диаграмм в соответствии с правилом КИрхгофа токовым методом
// ILenDiagrArr - к-во точек по каждой диаграмме
//
// OUTPUT:
//cmparrVeerDiagrs[ ILenDiagrArr * mNumWaveGuide] - массив диаграмм
// хранятся по очереди, начиная с ниждней(нижний волновод)
void  TOpticParabAnt::createVeerDiagrams(wchar_t *wchFoldPict, const double  VAlLambda,  const double VAlContrReflAng
	 , const int ILenDiagrArr, TComp *cmparrVeerDiagrs )
{
  double valMaxTargAng = 0.135;//80./180 * M_PI;
  double valTargAngStep =  2. * valMaxTargAng /((double)ILenDiagrArr);
  TComp *pcmparrDiagrCur =  new TComp[mNumWaveGuide];
  double *parrTargAng = new double [mNumWaveGuide];
  for (int i =0; i < ILenDiagrArr; i++)
  {
	if (i == 1960)
	{
      int iii=0;
	}
	double valTargAngCur =  -valMaxTargAng + ((double)i) * valTargAngStep ;
	createDiagrArrayFromDirect(VAlLambda, VAlContrReflAng,valTargAngCur, pcmparrDiagrCur );
	//calcDiagrArrayFromDirect(VAlLambda, VAlContrReflAng,0., pcmparrDiagrCur );
	for (int j =0; j < mNumWaveGuide; j++)
	{
	  cmparrVeerDiagrs[j * ILenDiagrArr + i] =   pcmparrDiagrCur[j];
	}
	parrTargAng[i] = valTargAngCur* 1000.;

  }
  wchar_t Fold[400] ={0.};
   wcscpy( Fold, wchFoldPict);
   wcscat(Fold, L"\\Diagrams");
	_wmkdir(Fold);


  for (int i=0; i < mNumWaveGuide; i++)
  {
	wchar_t FileName[400] ={0.};
	wcscpy( FileName, Fold);
	wcscat(FileName, L"\\ModulDiagrNo_");
	wchar_t string[10] = {0.};
	_itow(i,string, 10);
	wcscat(FileName, string);
	wcscat(FileName, L".shp");

	wchar_t FileName1[400] ={0.};
	wcscpy( FileName1, Fold);
	wcscat(FileName1, L"\\ArgDiagrNo_");
	wchar_t string1[10] = {0.};
	_itow(i,string1, 10);
	wcscat(FileName1, string1);
	wcscat(FileName1, L".shp");

	double *parrModul = new double [ILenDiagrArr];
	double *parrArg = new double [ILenDiagrArr];
	for (int j = 0;j < ILenDiagrArr; j++)
	{
	  parrModul[j] =  cmparrVeerDiagrs[i * ILenDiagrArr + j].modul();
	  if (parrModul[j] > 0.000001)
	  {
	   parrArg[j] =  cmparrVeerDiagrs[i * ILenDiagrArr + j].phase() ;
	  }
	  else
	  {
		  parrArg[j] = 0.;
	  }
	}
	double scalex =1.,  scaley = 0.01;
	TYrWriteShapeFile::CreateShpFile(FileName, parrModul, parrTargAng
	 ,ILenDiagrArr, scalex, scaley);
	  scaley =1.;
	 TYrWriteShapeFile::CreateShpFile(FileName1,parrArg, parrTargAng
	 ,ILenDiagrArr, scalex, scaley);

	delete parrModul;
	delete parrArg;

  }
  delete []pcmparrDiagrCur;
  delete parrTargAng ;

    wchar_t FileName1[400] ={0.};
  	wcscpy( FileName1, Fold);
	wcscat(FileName1, L"\\Axes.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(FileName1,-1000.,1000.
	 ,-100.,100.,5.) ;
}


void  TOpticParabAnt::createDiagrArrayFromDirect(const double  VAlLambda
		 ,const double VAlContrReflAlf, const double VAlTargTetta, TComp *pcmparrDiagrCur )
{
	memset(pcmparrDiagrCur,0, mNumWaveGuide * sizeof(TComp));

	const double  VAlA0 =1000.;// уровень амплитуды(коэффициент , который сокращается, но нужен для вычислений)

	int iC = 1200;
	double valYStep =  mParabDiam / ((double)iC);

	// создание массива сетментов концов волноводов
	TURPointXY *ppntarrWaveGuideCentre =  new TURPointXY [ mNumWaveGuide];
	for (int i = 0; i < mNumWaveGuide; i++)
	{
	  ppntarrWaveGuideCentre [i] = createWaveGuideCentre(i );
	}
	 for (int i = 0; i < iC; i++)
	 {
		 double valYTemp = -mParabDiam/2. + ((double)i)* valYStep ;
		 TURPointXY pntReflParab(- valYTemp* valYTemp/2./ ml + ml/2., valYTemp);
		 TComp cmpW = calcW( valYTemp, VAlLambda, VAlContrReflAlf, VAlTargTetta);
		for (int j =0; j < mNumWaveGuide; j++)
		{

		  double valR = TURPointXY::dist(pntReflParab,ppntarrWaveGuideCentre[j]);
		  double valFi = atan( (valYTemp - ppntarrWaveGuideCentre[j].Y)/pntReflParab.X);
		  double fi =  valR  * 2. * M_PI / VAlLambda;
		  pcmparrDiagrCur[j] +=  cmpW * (TComp(mpArrWaveGuide[j].fncDiagr(VAlLambda, valFi )/valR* 50.* valYStep, 0.)) * exp_(fi) ;

		}
	 }
	 delete []ppntarrWaveGuideCentre;

}

TComp TOpticParabAnt::calcW( double valYTemp, double VAlLambda, const double VAlContrReflAlf, const double VAlTargTetta)
{
	const double  VAlA0 =1000.;// уровень амплитуды(коэффициент , который сокращается, но нужен для вычислений)

	int iC = 1200 ;
	double valXStep =  mParabDiam / ((double)iC);
	TComp cmpRez (0.,0.);

	 TURPointXY pntReflParab(- valYTemp* valYTemp/2./ ml + ml/2., valYTemp);

	 for (int i = 0; i < iC; i++)
	 {
		 double valPosTemp = -mParabDiam/2. + ((double)i)* valXStep ;
		 TURPointXY pntTemp(valPosTemp * cos(M_PI/2. + VAlContrReflAlf), valPosTemp * sin(M_PI/2. + VAlContrReflAlf));
		  double valR = TURPointXY::dist(pntReflParab,pntTemp);
		  TURPointXY pntDiffer (pntReflParab.X - pntTemp.X, pntReflParab.Y - pntTemp.Y);
		  TURPointXY pntParabNorm(-ml, -valYTemp);
		  double valCosGam = TURPointXY::ScalMult( pntParabNorm,pntDiffer )/ pntParabNorm.Norm()/ pntDiffer.Norm();
		  cmpRez += exp_(-valPosTemp * sin(VAlTargTetta) * 2. * M_PI / VAlLambda)
			* exp_(valR * 2.* M_PI / VAlLambda) * TComp(valCosGam *VAlA0 * valXStep/ valR, 0.);

	 }
	 cmpRez = cmpRez * TComp(cos( VAlTargTetta), 0.);
	return cmpRez;
}

// постронеие массива диаграмм в соответствии с правилом КИрхгофа токовым методом
// ILenDiagrArr - к-во точек по каждой диаграмме
//
// OUTPUT:
//cmparrVeerDiagrs[ ILenDiagrArr * mNumWaveGuide] - массив диаграмм
// хранятся по очереди, начиная с ниждней(нижний волновод)
void  TOpticParabAnt::createVeerDiagramsNew(wchar_t *wchFoldPict, const double  VAlLambda,  const double VAlContrReflAng
	 , const double VAlAngDiagrStep1, const double VAlMaxTargTetta)// , TComp *cmparrVeerDiagrs )
{

  const int ILenDiagrArr =  2. * VAlMaxTargTetta / VAlAngDiagrStep1 +1;
  const double VAlAngDiagrStep = 2. * VAlMaxTargTetta / ((double)ILenDiagrArr);
  TComp *pcmparrDiagrCur =  new TComp[ILenDiagrArr];
  double *parrDerivMod = new double [ILenDiagrArr];
  double *parrDerivArg = new double [ILenDiagrArr];
  double *parrTargAng = new double [ILenDiagrArr];
  double *parrModul = new double [ILenDiagrArr];
  double *parrArg = new double [ILenDiagrArr];
  // массив распрелделения токов га входе в антенну на регулярной сетке по отрезку AB
  const int LEnInpCurrentArr = 1200;
  TComp *pcmparrFGr  = new  TComp[LEnInpCurrentArr];
  TComp *pcmparrFGr1 = new  TComp[LEnInpCurrentArr];
  ///

  // создание папки с файлами
  for (int i =0; i < ILenDiagrArr; i++)
  {
	double valTargAngCur =  -VAlMaxTargTetta + ((double)i) * VAlAngDiagrStep ;
	parrTargAng[i] = valTargAngCur* 1000.;
  }

  wchar_t Fold[400] ={0.};
   wcscpy( Fold, wchFoldPict);
   wcscat(Fold, L"\\Diagrams");
	_wmkdir(Fold);
   ///


  // массив весовой фугнкции антенны на входе в антенну на регулярной сетке по отрезку AB
  TComp *pcmparrQ = new  TComp[LEnInpCurrentArr];
  ///
  for (int n =0; n < mNumWaveGuide; n++)
  {
	 TURPointXY pntWaveGuideCentre = createWaveGuideCentre(n )  ;
	 createWeightFuncArrNew(VAlLambda,VAlContrReflAng, LEnInpCurrentArr, pntWaveGuideCentre, pcmparrQ );
	 for (int i =0; i < ILenDiagrArr; i++)
	 {
		const double valTargEps = 2. * VAlContrReflAng -VAlMaxTargTetta + ((double)i) * VAlAngDiagrStep ;
		createInpCurrentArray(VAlLambda,valTargEps, LEnInpCurrentArr, pcmparrFGr, pcmparrFGr1);
		//MtrxMultMatrx(pcmparrFGr ,1, LEnInpCurrentArr , pcmparrQ,1, &cmparrVeerDiagrs[n * ILenDiagrArr + i]) ;
		TComp cmpF , cmpDerivF;
		MtrxMultMatrx(pcmparrFGr ,1, LEnInpCurrentArr , pcmparrQ,1, &cmpF) ;
		pcmparrDiagrCur[ i] = cmpF;
		MtrxMultMatrx(pcmparrFGr1 ,1, LEnInpCurrentArr , pcmparrQ,1, &cmpDerivF) ;
		parrModul[i] =  cmpF.modul();
		if (parrModul[i] > 0.000001)
		{
		parrArg[i] =  cmpF.phase() ;
		}
		else
		{
		  parrArg[i] = 0.;
		}

		// вычисление производнорй модуля
		TComp cmpFSopr = cmpF.Sopr() ;
		TComp cmpTemp =  cmpFSopr * cmpDerivF ;
		parrDerivMod [i] = cmpTemp.m_Re /parrModul[i];
		parrDerivArg [i] = cmpTemp.m_Im /parrModul[i]/parrModul[i];
	 }


	wchar_t FileName[400] ={0.};
	wcscpy( FileName, Fold);
	wcscat(FileName, L"\\ModulDiagrNo_");
	wchar_t string[10] = {0.};
	_itow(n,string, 10);
	wcscat(FileName, string);
	wcscat(FileName, L".shp");

	wchar_t FileName1[400] ={0.};
	wcscpy( FileName1, Fold);
	wcscat(FileName1, L"\\ArgDiagrNo_");
	wchar_t string1[10] = {0.};
	_itow(n,string1, 10);
	wcscat(FileName1, string1);
	wcscat(FileName1, L".shp");

	wchar_t FileName2[400] ={0.};
	wcscpy( FileName2, Fold);
	wcscat(FileName2, L"\\DerivModDiagrNo_");
	wchar_t string2[10] = {0.};
	_itow(n,string2, 10);
	wcscat(FileName2, string2);
	wcscat(FileName2, L".shp");

	wchar_t FileName3[400] ={0.};
	wcscpy( FileName3, Fold);
	wcscat(FileName3, L"\\DerivArgDiagrNo_");
	wchar_t string3[10] = {0.};
	_itow(n,string3, 10);
	wcscat(FileName3, string3);
	wcscat(FileName3, L".shp");

	double scalex =1.,  scaley = 100.;//0.01;
	TYrWriteShapeFile::CreateShpFile(FileName, parrModul, parrTargAng
	 ,ILenDiagrArr, scalex, scaley);
	// scaley *= 0.001;
	TYrWriteShapeFile::CreateShpFile(FileName2, parrDerivMod, parrTargAng
	 ,ILenDiagrArr, scalex, scaley);
	  scaley =1.;
	 TYrWriteShapeFile::CreateShpFile(FileName1,parrArg, parrTargAng
	 ,ILenDiagrArr, scalex, scaley);
	 TYrWriteShapeFile::CreateShpFile(FileName3,parrDerivArg, parrTargAng
	 ,ILenDiagrArr, scalex, scaley);
  }


/*

  for (int i=0; i < mNumWaveGuide; i++)
  {
	wchar_t FileName[400] ={0.};
	wcscpy( FileName, Fold);
	wcscat(FileName, L"\\ModulDiagrNo_");
	wchar_t string[10] = {0.};
	_itow(i,string, 10);
	wcscat(FileName, string);
	wcscat(FileName, L".shp");

	wchar_t FileName1[400] ={0.};
	wcscpy( FileName1, Fold);
	wcscat(FileName1, L"\\ArgDiagrNo_");
	wchar_t string1[10] = {0.};
	_itow(i,string1, 10);
	wcscat(FileName1, string1);
	wcscat(FileName1, L".shp");

	double *parrModul = new double [ILenDiagrArr];
	double *parrArg = new double [ILenDiagrArr];
	for (int j = 0;j < ILenDiagrArr; j++)
	{
	  parrModul[j] =  cmparrVeerDiagrs[i * ILenDiagrArr + j].modul();
	  if (parrModul[j] > 0.000001)
	  {
	   parrArg[j] =  cmparrVeerDiagrs[i * ILenDiagrArr + j].phase() ;
	  }
	  else
	  {
		  parrArg[j] = 0.;
	  }
	}
	double scalex =1.,  scaley = 0.01;
	TYrWriteShapeFile::CreateShpFile(FileName, parrModul, parrTargAng
	 ,ILenDiagrArr, scalex, scaley);
	  scaley =1.;
	 TYrWriteShapeFile::CreateShpFile(FileName1,parrArg, parrTargAng
	 ,ILenDiagrArr, scalex, scaley);

	delete parrModul;
	delete parrArg;

  }
  */
  delete []pcmparrDiagrCur;
  delete parrTargAng ;

    wchar_t FileName1[400] ={0.};
  	wcscpy( FileName1, Fold);
	wcscat(FileName1, L"\\Axes.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(FileName1,-1000.,1000.
	 ,-200.,200.,5.) ;

	 delete []pcmparrFGr;
	 delete []pcmparrFGr1;
	 delete []pcmparrQ;
	 delete parrDerivMod;
	 delete parrDerivArg;
	 delete parrModul;
	 delete parrArg;
}

//  расчет массива весовой функции антенны на регулярной сетке на отрезке AB
// INPUT:
// VAlLambda - длина волны
// LEnInpCurrentArr - к-во узлов (отрезков разбиения) отрезка AB
// PNtWaveGuideCentre - точка-центр волновода
// OUTPUT:
//  pcmparrQ[LEnInpCurrentArr] - массив весовой функции
void TOpticParabAnt::createWeightFuncArr(const double  VAlLambda,  const double VAlContrReflAng,const int LEnInpCurrentArr
, const TURPointXY PNtWaveGuideCentre, TComp *pcmparrQ )
{
  double valStep = mParabDiam / ((double)LEnInpCurrentArr);
 // создание заготовок массивов
 double *parr_r = new double [LEnInpCurrentArr];
 double *parr_ro = new double [LEnInpCurrentArr * LEnInpCurrentArr];
 double *parr_R = new double [LEnInpCurrentArr];
 double *parr_cosBetta = new double [LEnInpCurrentArr];
 double *parr_SqrtDcosOmega= new double [LEnInpCurrentArr];
 double *parr_cosGamma = new double [LEnInpCurrentArr * LEnInpCurrentArr];
 ///

  // заполненеи массивов parr_cosGamma и parr_ro
  for (int i =0; i < LEnInpCurrentArr; i++) // цикл по отрезку CD (по х)
	 {
	   double valX =  calTay( LEnInpCurrentArr, i)  ;// расстояние X на контррефлектгоре
	   TURPointXY pntX( -valX* sin(VAlContrReflAng) ,valX* cos(VAlContrReflAng)); // точка на контррефлекторе

	   for (int j =0; j < LEnInpCurrentArr; j++) // цикл по параболе (по У)
	   {
		 double valY =  calTay( LEnInpCurrentArr, j)  ;// высота точки на параболе
		 TURPointXY pntY(- valY * valY/ ml/2. +ml/2. ,valY); // точка на параболе
		 parr_ro[ i *  LEnInpCurrentArr + j] = TURPointXY::dist(pntX, pntY) ;

		 TURPointXY pntDiffer(pntX.X - pntY.X, pntX.Y - pntY.Y);
		// TURPointXY pntN(cos(VAlContrReflAng), sin(VAlContrReflAng));
		 TURPointXY pntN(-ml, -valY);
		 parr_cosGamma[ i *  LEnInpCurrentArr + j] = TURPointXY::calcCosAng( pntN, pntDiffer) ;
	   }
	 }

  ///

  // заполненеи массивов parr_SqrtDcosOmega и parr_R

	   for (int j =0; j < LEnInpCurrentArr; j++) // цикл по параболе (по У)
	   {
		 double valY =  calTay( LEnInpCurrentArr, j)  ;// высота точки на параболе
		 TURPointXY pntY(- valY * valY/ ml/2. +ml/2. ,valY); // точка на параболе
		 parr_R[j] = TURPointXY::dist(PNtWaveGuideCentre, pntY) ;

		 TURPointXY pntDiffer(pntY.X - PNtWaveGuideCentre.X, pntY.Y - PNtWaveGuideCentre.Y);

		 double valCos = pntY.X /parr_R[j];
		 double valOmega = acos(valCos);
		 parr_SqrtDcosOmega [j] = valCos * sqrt(1. + valY * valY/ ml/ml)
			  * mpArrWaveGuide[0].fncDiagr(VAlLambda, valOmega) ;
	   }
  ///


 //
 for (int n=0; n < LEnInpCurrentArr; n++)// цикл по отрезку АВ (по тау)
 {
	 double valTay = calTay( LEnInpCurrentArr, n)  ;
   //	 TURPointXY pntTay( ml/ 2. ,valTay);// точка входа в которой считается весовая функция  ОШИБКА!!!
   TURPointXY pntTay( ml/ 2. - mParabDiam * mParabDiam / 8./ ml,valTay);// точка входа в которой считается весовая функция
	  // заполненеи массивов parr_cosBetta и parr_r
	  for (int i =0; i < LEnInpCurrentArr; i++) // цикл по отрезку CD (по х)
	 {
	   double valX =  calTay( LEnInpCurrentArr, i)  ;// расстояние X на контррефлектгоре
	   TURPointXY pntX( -valX* sin(VAlContrReflAng) ,valX* cos(VAlContrReflAng)); // точка на контррефлекторе
	   parr_r[i] =  TURPointXY::dist(pntX, pntTay) ;
	   TURPointXY pntDiffer(pntTay.X - pntX.X, pntTay.Y - pntX.Y);
	   TURPointXY pntN(cos(VAlContrReflAng), sin(VAlContrReflAng));
	   parr_cosBetta[ i ] = TURPointXY::calcCosAng( pntN, pntDiffer) ;
	 }
	 ///

	   pcmparrQ[n] = TComp(0.,0.);
	  ///

	 for (int i =0; i < LEnInpCurrentArr; i++) // цикл по отрезку CD (по х)
	 {

	   for (int j =0; j < LEnInpCurrentArr; j++) // цикл по параболе (по У)
	   {
	   double valt = parr_cosBetta[i] * parr_SqrtDcosOmega[j] * parr_cosGamma[i * LEnInpCurrentArr + j]
		 / (parr_r[i] * parr_ro[i * LEnInpCurrentArr + j] * parr_R[j]);
	   pcmparrQ[n] += exp_(2. * M_PI /VAlLambda * (parr_r[i] + parr_ro[i * LEnInpCurrentArr + j] + parr_R[j]) )
	  // * TComp(valt,0.) * TComp(valStep * valStep* 10000., 0.);
	   * TComp(valt,0.) * TComp(valStep * valStep* valStep, 0.);
	   }
	 }
 }

 delete parr_r ;
 delete parr_ro ;
 delete parr_R ;
 delete parr_cosBetta ;
 delete parr_SqrtDcosOmega;
 delete parr_cosGamma  ;
 // const double VAlSegmLength =1000.; // длина сегмента полулинии по которой идет интегрирование
 // const double STepIntegr = 0.05;

}


//  расчет массива весовой функции антенны на регулярной сетке на отрезке AB
// INPUT:
// VAlLambda - длина волны
// LEnInpCurrentArr - к-во узлов (отрезков разбиения) отрезка AB
// PNtWaveGuideCentre - точка-центр волновода
// OUTPUT:
//  pcmparrQ[LEnInpCurrentArr] - массив весовой функции
void TOpticParabAnt::createWeightFuncArrNew(const double  VAlLambda,  const double VAlContrReflAng,const int LEnInpCurrentArr
, const TURPointXY PNtWaveGuideCentre, TComp *pcmparrQ )
{
  double valStep = mParabDiam / ((double)LEnInpCurrentArr);
 // создание заготовок массивов
 double *parr_r = new double [LEnInpCurrentArr];
 double *parr_ro = new double [LEnInpCurrentArr * LEnInpCurrentArr];
 double *parr_R = new double [LEnInpCurrentArr];
 double *parr_cosBetta_cosDelta = new double [LEnInpCurrentArr];
 double *parr_SqrtDcosNu= new double [LEnInpCurrentArr];
 double *parr_cosGamma_cosAlfa = new double [LEnInpCurrentArr * LEnInpCurrentArr];
 ///

  // заполненеи массивов parr_cosGamma_cosAlfa и parr_ro
  for (int i =0; i < LEnInpCurrentArr; i++) // цикл по отрезку CD (по х)
	 {
	   double valX =  calTay( LEnInpCurrentArr, i)  ;// расстояние X на контррефлектгоре
	   TURPointXY pntX( -valX* sin(VAlContrReflAng) ,valX* cos(VAlContrReflAng)); // точка на контррефлекторе

	   for (int j =0; j < LEnInpCurrentArr; j++) // цикл по параболе (по У)
	   {
		 double valY =  calTay( LEnInpCurrentArr, j)  ;// высота точки на параболе
		 TURPointXY pntY(- valY * valY/ ml/2. +ml/2. ,valY); // точка на параболе
		 parr_ro[ i *  LEnInpCurrentArr + j] = TURPointXY::dist(pntX, pntY) ;

		 TURPointXY pntDiffer(pntX.X - pntY.X, pntX.Y - pntY.Y);

		 TURPointXY pntN(-ml, -valY);
		// double valCosGamma = TURPointXY::calcCosAng( pntN, pntDiffer) ;

		TURPointXY pntN1(cos(VAlContrReflAng), sin(VAlContrReflAng));
		double valCosAlfa = fabs(TURPointXY::calcCosAng( pntN1, pntDiffer)) ;
		parr_cosGamma_cosAlfa [ i *  LEnInpCurrentArr + j] = /*sqrt(valCosGamma) * */ valCosAlfa;
	   }
	 }

  ///

  // заполненеи массивов parr_SqrtDcosNu и parr_R

	   for (int j =0; j < LEnInpCurrentArr; j++) // цикл по параболе (по У)
	   {
		 double valY =  calTay( LEnInpCurrentArr, j)  ;// высота точки на параболе
		 TURPointXY pntY(- valY * valY/ ml/2. +ml/2. ,valY); // точка на параболе
		 parr_R[j] = TURPointXY::dist(PNtWaveGuideCentre, pntY) ;

		 TURPointXY pntDiffer(pntY.X - PNtWaveGuideCentre.X, pntY.Y - PNtWaveGuideCentre.Y);

		 double valCos = pntY.X /parr_R[j];
		 double valOmega = acos(valCos);
		 TURPointXY pntN(-ml, -valY);// нормаль к параболе
		 double valCosNu = TURPointXY::calcCosAng( pntN, pntDiffer)  ;
		 parr_SqrtDcosNu [j] = valCosNu * sqrt(1. + valY * valY/ ml/ml)
			  * mpArrWaveGuide[0].fncDiagr(VAlLambda, valOmega) ;
	   }
  ///


 //
 for (int n=0; n < LEnInpCurrentArr; n++)// цикл по отрезку АВ (по тау)
 {
	 double valTay = calTay( LEnInpCurrentArr, n)  ;
   //	 TURPointXY pntTay( ml/ 2. ,valTay);// точка входа в которой считается весовая функция  ОШИБКА!!!
   TURPointXY pntTay( ml/ 2. - mParabDiam * mParabDiam / 8./ ml,valTay);// точка входа в которой считается весовая функция
	  // заполненеи массивов parr_cosBetta_cosDelta и parr_r
	  for (int i =0; i < LEnInpCurrentArr; i++) // цикл по отрезку CD (по х)
	 {
	   double valX =  calTay( LEnInpCurrentArr, i)  ;// расстояние X на контррефлектгоре
	   TURPointXY pntX( -valX* sin(VAlContrReflAng) ,valX* cos(VAlContrReflAng)); // точка на контррефлекторе
	   parr_r[i] =  TURPointXY::dist(pntX, pntTay) ;
	   TURPointXY pntDiffer(pntTay.X - pntX.X, pntTay.Y - pntX.Y);
	   TURPointXY pntN(cos(VAlContrReflAng), sin(VAlContrReflAng));

	  // double valCosBetta = TURPointXY::calcCosAng( pntN, pntDiffer) ;
	   double valCosDelta = fabs(pntDiffer.X/ parr_r[i]);
	   parr_cosBetta_cosDelta[ i ] = /*sqrt(valCosBetta) **/  valCosDelta ;
	 }
	 ///

	   pcmparrQ[n] = TComp(0.,0.);
	  ///

	 for (int i =0; i < LEnInpCurrentArr; i++) // цикл по отрезку CD (по х)
	 {

	   for (int j =0; j < LEnInpCurrentArr; j++) // цикл по параболе (по У)
	   {
	   double valt = parr_cosBetta_cosDelta[i] * parr_SqrtDcosNu[j] * parr_cosGamma_cosAlfa[i * LEnInpCurrentArr + j]
		 / (parr_r[i] * parr_ro[i * LEnInpCurrentArr + j] * parr_R[j]);
	   pcmparrQ[n] += exp_(2. * M_PI /VAlLambda * (parr_r[i] + parr_ro[i * LEnInpCurrentArr + j] + parr_R[j]) )
	   * TComp(valt,0.) * TComp(valStep * valStep* valStep, 0.);
	   }
	 }
 }

 delete parr_r ;
 delete parr_ro ;
 delete parr_R ;
 delete parr_cosBetta_cosDelta ;
 delete parr_SqrtDcosNu;
 delete parr_cosGamma_cosAlfa  ;
 // const double VAlSegmLength =1000.; // длина сегмента полулинии по которой идет интегрирование
 // const double STepIntegr = 0.05;

}

void TOpticParabAnt::createInpCurrentArray(const double  VAlLambda,double valTargEps, const int LEnInpCurrentArr
  , TComp *pcmparrFGr, TComp *pcmparrFGr1)
{
 // const int ICircle =  VAlSegmLength/ STepIntegr;
  for (int i =0; i < LEnInpCurrentArr; i++)
  {
	 double valTay = calTay(LEnInpCurrentArr, i);
	// pcmparrFGr[i] = calcInpCurrentInPoint(VAlLambda,valTargEps, valTay , VAlSegmLength, STepIntegr);
	 pcmparrFGr[i] = TComp(cos(valTargEps), 0.) * exp_(-valTay * sin(valTargEps)/VAlLambda * 2. * M_PI);
	 TComp cmpTemp(- tan(valTargEps)/2., - 2. * M_PI * valTay * cos(valTargEps)/VAlLambda);
	 pcmparrFGr1[i] = cmpTemp * pcmparrFGr[i];
  }
  // для проверки
 /*  double *parrTau = new double [LEnInpCurrentArr];
	double *parrModul = new double [LEnInpCurrentArr];
	 double *parrArg = new double [LEnInpCurrentArr];

   for (int i =0; i < LEnInpCurrentArr; i++)
  {
	  parrTau[i]  = calTay(LEnInpCurrentArr, i);
	 parrModul[i] = pcmparrFGr[i].modul() ;
	 parrArg[i] =  pcmparrFGr[i].phase() ;
  }

	double scalex =1.,  scaley = 1.;
	TYrWriteShapeFile::CreateShpFile(L"E:\\PROJECTS_C++\\Optic_5P10-03\\STUDY\\Arg.shp",parrArg, parrTau
	 ,LEnInpCurrentArr, scalex, scaley);
	TYrWriteShapeFile::CreateShpFile(L"E:\\PROJECTS_C++\\Optic_5P10-03\\STUDY\\Modul.shp",parrModul, parrTau
	 ,LEnInpCurrentArr, scalex, scaley);
	 TYrWriteShapeFile::CreateShpArrowedAxes(L"E:\\PROJECTS_C++\\Optic_5P10-03\\STUDY\\Axes.shp",-1000.,1000.
	 ,-100.,100.,5.) ;

   delete parrTau;
   delete parrModul;
   delete parrArg ; */
   ///
}

// вычисление координаты тау (Y) узла регулярной сетки отрезка АВ входа антиенны с номером INum
// INPUT
// LEnInpCurrentArr - к-во узлов (отрезков разбиения)
double TOpticParabAnt::calTay(const int LEnInpCurrentArr, const int INum)
{
	return -mParabDiam/2. + mParabDiam/2./ ((double) LEnInpCurrentArr) + ((double)INum) *  mParabDiam/((double) LEnInpCurrentArr);
}

// вычисление поля в точке входа в антенну с координатой valTay
// INPUT:
// valTargEps - угл места цели относительно нормали к антенне
// VAlSegmLength - длина сегмента интегрирования
// STepIntegr - шагт интегрирования
TComp TOpticParabAnt::calcInpCurrentInPoint(const double  VAlLambda,double valTargEps, double valTay
   ,const double VAlSegmLength, const double STepIntegr )
{
   if (fabs(valTargEps) < 0.000001)
   {
	   return TComp(1.,0.);
   }

  const int ICircle =  VAlSegmLength / STepIntegr;
  TComp cmpRez(0.,0.);
  TURPointXY pntTau(0.,valTay);
  for (int i =0 ; i < ICircle; i++)
  {
	// вычисление координат центра отрезка разбиения с номером i на отрезке интегрирования
	double valt = STepIntegr/2. + ((double)i)* STepIntegr;
	TURPointXY pntCentre;
	if (valTargEps > 0.)
	{
	 pntCentre =  TURPointXY (valt * sin(valTargEps), mParabDiam/2. - valt * cos(valTargEps) );
	}
	else
	{
	  pntCentre =  TURPointXY (-valt * sin(valTargEps), -mParabDiam /2. + valt * cos(valTargEps) );
	}

	double val_l = TURPointXY::dist(pntTau,pntCentre);
	double valCosPsi = pntCentre.X / val_l;
	cmpRez += exp_(val_l * 2. * M_PI/ VAlLambda)* TComp(valCosPsi /val_l * STepIntegr, 0.);
  }
  return cmpRez;
}

// вычисление поля в плоской волны
// волна вертикальная
// поле вычисляется в точке (0,0)
// фронт волны расположен в точке valA0
// INPUT:
// [-VAlSegmLength ; VAlSegmLength ] -  сегмент интегрирования
// STepIntegr - шаг интегрирования
// VAlLambda - длина волны
// valA0 - положение фронта фолны по оси OX
TComp calcEMFieldPlaneWave_(const double  VAlLambda,double valA0
   ,const double VAlSegmLength, const double STepIntegr )
{
  const int ICircle =  2. * VAlSegmLength / STepIntegr;
  TComp cmpRez(0.,0.);
  double valt = 0.;
  for (int i =0 ; i < ICircle; i++)
  {
	// вычисление координат центра отрезка разбиения с номером i на отрезке интегрирования
	valt = -VAlSegmLength + ((double)i)* STepIntegr;
	TURPointXY pntCur(valA0,valt) ;


	double val_l = pntCur.Norm();
	double valCosPsi = valA0 / val_l;
   //	double psi = acos(valCosPsi) - M_PI/2.;
   //	valCosPsi = cos(psi);
	cmpRez += exp_(val_l * 2. * M_PI/ VAlLambda)* TComp(valCosPsi /val_l * STepIntegr, 0.);
   //	cmpRez += exp_(val_l * 2. * M_PI/ VAlLambda)* TComp(sqrt(valCosPsi) /val_l * STepIntegr, 0.);
	//cmpRez += exp_(val_l * 2. * M_PI/ VAlLambda)* TComp(1. /val_l * STepIntegr*1000., 0.);
  }
  cmpRez = cmpRez / TComp(0., VAlLambda);

  return cmpRez;
}

// вычисление поля в плоской волны  в произвольной точке на оси OY
// волна вертикальная
// поле вычисляется в точке (0,Tay)
// фронт волны проходит через точку (valWaveX0, 0) и имеет угол с осью OX AlfaWave
// угол нормали фронта волны равенн   valAlfaWave
// INPUT:
// valTay -
// [-VAlSegmLength ; VAlSegmLength ] -  сегмент по которому идет интегрирование  на плоской волне.
// за 0 принята точка наиболее близко расположенная к началу координат
// STepIntegr - шаг интегрирования
// VAlLambda - длина волны
// valWaveX0 - точка на оси OY через которую проходит фронт волны
// valAlfaWave -  угол наклона нормали фронта волны к оси OX

TComp calcEMFieldPlaneWave(const double  VAlLambda,double valTay, double valWaveX0,
   double valAlfaWave,const double VAlSegmLength, const double STepIntegr )
{
   // проекция начала координат на линию фронта волны
   TURPointXY pntProj0(valWaveX0*  cos (valAlfaWave)*  cos (valAlfaWave), valWaveX0 *  sin (valAlfaWave)*  cos (valAlfaWave));
   ///

   // нормаль к фронту
   TURPointXY pntNormFront( -cos( valAlfaWave), -sin( valAlfaWave));
   ///

   // точка пространства
   TURPointXY pntTay(0., valTay);
   ///


  const int ICircle =  2. * VAlSegmLength / STepIntegr;
  TComp cmpRez(0.,0.);
  double valt = 0.;
  for (int i =0 ; i < ICircle; i++)
  {
	// вычисление координат текущей точки интегрирования  с номером i на отрезке интегрирования
	valt = -VAlSegmLength + ((double)i)* STepIntegr; // расстояние от точки  pntProj0

	TURPointXY pntCur(pntProj0.X - valt * sin (valAlfaWave) ,pntProj0.Y + valt *cos (valAlfaWave)) ;
	TURPointXY pntDiffer (pntTay.X - pntCur.X, pntTay.Y - pntCur.Y);


	double val_r = pntDiffer.Norm();

	double valCosPsi = fabs(TURPointXY::calcCosAng( pntDiffer, pntNormFront) ) ;
	double cosAlf = fabs( pntDiffer.X / val_r);
  
	cmpRez += exp_(-val_r * 2. * M_PI/ VAlLambda)* TComp(valCosPsi /val_r * (STepIntegr) /** cosAlf*/, 0.);

  }

  double valA0 =1.;
  cmpRez = cmpRez / TComp(0., VAlLambda) * TComp (valA0, 0.);

  return cmpRez;
}
#pragma package(smart_init)
