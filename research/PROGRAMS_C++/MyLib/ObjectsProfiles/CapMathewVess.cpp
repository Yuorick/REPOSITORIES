//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include <string.h>
#include "CapMathewVess.h"

#include "URPolygon.h"
#include "URPointXY.h"
#include "URPolyLine.h"
#include "YrWriteShapeFile.h"
#include "UABPlg.h"

void drawVessMath( wchar_t *pwcharrFileOut,const double  valAngRotation, const TURPointXY pntSdvig,const double valRastigenie)
{

	TURPolygon PlgVessMath(NumPartsVessMath, NumPointsVessMath, PartsVessMath, PointsVessMath);


  TURPolygon PlgVessMath1 = PlgVessMath.LinTransform( valAngRotation,  pntSdvig  , valRastigenie );
  if (fabs(valAngRotation)> M_PI/ 2.)
  {
   PlgVessMath1 = PlgVessMath1.SimOtragenieTransform();
  }


   PlgVessMath1.WriteSetSHPFiles(pwcharrFileOut, &PlgVessMath1,1);

}

void drawCoord( wchar_t *pwcharrFileOut,const double  valAngRotation, const TURPointXY pntSdvig,const double valRastigenie)
{

 TURPolygon PlgCoord(NumPartsCoord, NumPointsCoord, PartsCoord, PointsCoord);


  TURPolygon PlgCoord1 = PlgCoord.LinTransform( valAngRotation,  pntSdvig  , valRastigenie );
  if (fabs(valAngRotation)> M_PI/ 2.)
  {
   PlgCoord1 = PlgCoord1.SimOtragenieTransform();
  }

   PlgCoord1.WriteSetSHPFiles(pwcharrFileOut, &PlgCoord1,1);

}

 void drawVesselMathewPicture( wchar_t *pwcharrPath,const double  valAngPsi,const double  valAngEps
	 , const double  valCentreVessX,const double valRastigenie)
{
	wchar_t pwcharrFileVess[300] ={0};
	wcscpy(pwcharrFileVess, pwcharrPath);
	wcscat(pwcharrFileVess, L"\\VessMatveev.shp");
	TURPolygon PlgVessMath(NumPartsVessMath, NumPointsVessMath, PartsVessMath, PointsVessMath);
	TURPointXY pntSdvig(valCentreVessX, 0.);
	TURPolygon PlgVessMath1 = PlgVessMath.LinTransform( valAngPsi,  pntSdvig  , valRastigenie );
	if (fabs(valAngPsi)> M_PI/ 2.)
  {
   PlgVessMath1 = PlgVessMath1.SimOtragenieTransform();
  }

	PlgVessMath1.WriteSetSHPFiles(pwcharrFileVess, &PlgVessMath1,1);

   wcscpy(pwcharrFileVess, pwcharrPath);
   wcscat(pwcharrFileVess, L"\\CoordMatveev.shp");
   TURPolygon PlgCoord(NumPartsCoord, NumPointsCoord, PartsCoord, PointsCoord);
   TURPointXY pntCentre(centreCoordX, centreCoordY);
   TURPointXY pntCentreNew =  pntCentre.LinTransform( valAngPsi,  pntSdvig  , valRastigenie );
  // TURPointXY pntSdvig(valCentreVessX, 0.);
   TURPolygon PlgCoord1 = PlgCoord.LinTransform( valAngPsi ,  pntSdvig  , valRastigenie );

   TURPointXY pntSdvig1( -pntCentreNew.X, -pntCentreNew.Y);
	TURPolygon PlgCoord2 = PlgCoord1.SdvigTransform( pntSdvig1 );

	TURPointXY pntSdvig00(0.,0.);
	TURPolygon PlgCoord3 = PlgCoord2.LinTransform(valAngEps ,  pntSdvig00,1. ) ;

	TURPolygon PlgCoord4 = PlgCoord3.SdvigTransform( pntCentreNew);
	if (fabs(valAngPsi)> M_PI/ 2.)
  {
   PlgCoord4 = PlgCoord4.SimOtragenieTransform();
  }

   PlgCoord4.WriteSetSHPFiles(pwcharrFileVess, &PlgCoord4,1);

   // построение диаграммы
	// диаграмма
   double valWidth = 12./180. * M_PI ;
   // Создание диаграммы с центром в точке  pointCentre
// радиусом  valR с шагом по углу iNUM
// valWidthMilrad -  ширина диаграммы в милирадианах
// valRotAng - угол поворота относительно горизонтали
// положит направление против час стрелки
const double valR = 1000.;
const double valRotAng = 0.;
const int iNUM = 3000;
TURPolygon Diagr = TURPolygon::fncCreateDiagr(pntCentreNew,  valR,valRotAng , valWidth, iNUM) ;
	wcscpy(pwcharrFileVess, pwcharrPath);
   wcscat(pwcharrFileVess, L"\\Diagr.shp");

	   Diagr.WriteSetSHPFiles(pwcharrFileVess, &Diagr ,1);

}
// перегруженная функция, дополнительно на выходе выдающая центр координатора
void drawVesselMathewPicture( wchar_t *pwcharrPath,const double  valAngPsi,const double  valAngEps
	 , const double  valCentreVessX,const double valRastigenie, TURPointXY &pntCentreNew)
{
	wchar_t pwcharrFileVess[300] ={0};
	wcscpy(pwcharrFileVess, pwcharrPath);
	wcscat(pwcharrFileVess, L"\\VessMatveev.shp");
	TURPolygon PlgVessMath(NumPartsVessMath, NumPointsVessMath, PartsVessMath, PointsVessMath);
	TURPointXY pntSdvig(valCentreVessX, 0.);
	TURPolygon PlgVessMath1 = PlgVessMath.LinTransform( valAngPsi,  pntSdvig  , valRastigenie );
	if (fabs(valAngPsi)> M_PI/ 2.)
  {
   PlgVessMath1 = PlgVessMath1.SimOtragenieTransform();
  }

	PlgVessMath1.WriteSetSHPFiles(pwcharrFileVess, &PlgVessMath1,1);

   wcscpy(pwcharrFileVess, pwcharrPath);
   wcscat(pwcharrFileVess, L"\\CoordMatveev.shp");
   TURPolygon PlgCoord(NumPartsCoord, NumPointsCoord, PartsCoord, PointsCoord);
   TURPointXY pntCentre(centreCoordX, centreCoordY);
   pntCentreNew =  pntCentre.LinTransform( valAngPsi,  pntSdvig  , valRastigenie );
  // TURPointXY pntSdvig(valCentreVessX, 0.);
   TURPolygon PlgCoord1 = PlgCoord.LinTransform( valAngPsi ,  pntSdvig  , valRastigenie );

   TURPointXY pntSdvig1( -pntCentreNew.X, -pntCentreNew.Y);
	TURPolygon PlgCoord2 = PlgCoord1.SdvigTransform( pntSdvig1 );

	TURPointXY pntSdvig00(0.,0.);
	TURPolygon PlgCoord3 = PlgCoord2.LinTransform(valAngEps ,  pntSdvig00,1. ) ;

	TURPolygon PlgCoord4 = PlgCoord3.SdvigTransform( pntCentreNew);
	if (fabs(valAngPsi)> M_PI/ 2.)
  {
   PlgCoord4 = PlgCoord4.SimOtragenieTransform();
  }

   PlgCoord4.WriteSetSHPFiles(pwcharrFileVess, &PlgCoord4,1);

	// Создание диаграммы с центром в точке  pointCentre
	// радиусом  valR с шагом по углу iNUM
	// valWidthMilrad -  ширина диаграммы в милирадианах
	// valRotAng - угол поворота относительно горизонтали
	// положит направление против час стрелки
	double valWidth = 4./180. * M_PI ;
	const double valR = 12000.;
	const double valRotAng = 0.;//-2./180. * M_PI ;
	const int iNUM = 3000;
	TURPolygon Diagr = TURPolygon::fncCreateDiagr(pntCentreNew,  valR,valRotAng , valWidth, iNUM) ;
	wcscpy(pwcharrFileVess, pwcharrPath);
	wcscat(pwcharrFileVess, L"\\Diagr.shp");
	Diagr.WriteSetSHPFiles(pwcharrFileVess, &Diagr ,1);

	// создание шейп файла с центральной осью диаграммы
	TURPointXY pntEndCentralAxe(pntCentreNew.X + valR * cos(valRotAng), pntCentreNew.Y + valR * sin(valRotAng));
	TURPolyLine plnCentral( pntCentreNew, pntEndCentralAxe) ;
	wcscpy(pwcharrFileVess, pwcharrPath);
	wcscat(pwcharrFileVess, L"\\CentralLine.shp");
	 plnCentral.WriteSetSHPFiles(pwcharrFileVess, & plnCentral ,1);
	 ///




	/*
	   const TURPointXY pointBeginX(-40., 0.);
  const TURPointXY pointEndX(40., 0.);
  const TURPointXY pointBeginY(0., -40.);
  const TURPointXY pointEndY (0., 40.);
  const double valLength = 5.;
  TURPolyLine AxesKGSK = TURPolyLine::fncCreateAxes( pointBeginX,  pointEndX
									   ,pointBeginY, pointEndY, valLength) ;
  wchar_t pwcharrKGSK[300] ={0};
  wcscpy(pwcharrKGSK, pwcharrPath);
  wcscat(pwcharrKGSK, L"\\KGSK.shp");
  AxesKGSK.WriteSetSHPFiles(pwcharrKGSK, &AxesKGSK ,1);


  const TURPointXY pntSdvig00(0.,0.);
  TURPolyLine AxesPSK = AxesKGSK.LinTransform(valAngPsi , pntSdvig00,1. ) ;
  wcscpy(pwcharrKGSK, pwcharrPath);
  wcscat(pwcharrKGSK, L"\\PSK.shp");
  AxesPSK.WriteSetSHPFiles(pwcharrKGSK, &AxesPSK ,1);

  // центр координатора
   TURPointXY pntCentre(centreCoordX, centreCoordY);
   TURPointXY pntCentreNew =  pntCentre.LinTransform( valAngPsi,  pntSdvig00  , 1. );
   ///

   // ПСК-АС
   TURPolyLine AxesPSK_AS = AxesPSK.SdvigTransform(pntCentreNew ) ;
   wcscpy(pwcharrKGSK, pwcharrPath);
   wcscat(pwcharrKGSK, L"\\PSK-AS.shp");
   AxesPSK_AS.WriteSetSHPFiles(pwcharrKGSK, &AxesPSK_AS ,1);
  ///

  ///АСК
 TURPolyLine AxesASK = AxesPSK.LinTransform(valAngEps , pntCentreNew,1. ) ;
   wcscpy(pwcharrKGSK, pwcharrPath);
   wcscat(pwcharrKGSK, L"\\ASK.shp");
  AxesASK.WriteSetSHPFiles(pwcharrKGSK, &AxesASK ,1);

  ///
   ///ЛСК
  TURPolyLine AxesLSK = AxesPSK.LinTransform(valAngEps + valTargU , pntCentreNew,1. ) ;
   wcscpy(pwcharrKGSK, pwcharrPath);
   wcscat(pwcharrKGSK, L"\\LSK.shp");
   AxesLSK.WriteSetSHPFiles(pwcharrKGSK, &AxesLSK ,1);
   ///

   // диаграмма
   double valFiDiagr = 12./180. * M_PI ;
   TURPolygon Diagr = TURPolygon::fncCreateSector(pntCentreNew, valTargR * 1.3,
					valAngPsi +valAngEps + valTargU -valFiDiagr/2.
	   ,  valAngPsi +valAngEps + valTargU +valFiDiagr/2., 500 );

	 wcscpy(pwcharrKGSK, pwcharrPath);
   wcscat(pwcharrKGSK, L"\\Diagr.shp");
	   Diagr.WriteSetSHPFiles(pwcharrKGSK, &Diagr ,1);

  // угл сектор
  wcscpy(pwcharrKGSK, pwcharrPath);
   wcscat(pwcharrKGSK, L"\\Sect1.shp");
 TYrWriteShapeFile::CreateAngleMarks(pwcharrKGSK, pntSdvig00
	 , AxesKGSK.Points[1], AxesPSK.Points[1],20.,22. ) ;
	 ///

 // угл сектор
  wcscpy(pwcharrKGSK, pwcharrPath);
   wcscat(pwcharrKGSK, L"\\Sect2.shp");
 TYrWriteShapeFile::CreateAngleMarks(pwcharrKGSK, pntCentreNew
	 , AxesPSK_AS.Points[1], AxesLSK.Points[1],10.,12. ) ;
	 ///
	*/

}


// перегруженная функция, дополнительно на выходе выдающая центр координатора  и конечную точку центральной оси диаграммы
void drawVesselMathewPicture( wchar_t *pwcharrPath,const double  valAngPsi,const double  valAngEps
	 , const double  valCentreVessX,const double valRastigenie, TURPointXY &pntCentreNew, TURPointXY &pntEngDiagrAxe)
{
	wchar_t pwcharrFileVess[300] ={0};
	wcscpy(pwcharrFileVess, pwcharrPath);
	wcscat(pwcharrFileVess, L"\\VessMatveev.shp");
	TURPolygon PlgVessMath(NumPartsVessMath, NumPointsVessMath, PartsVessMath, PointsVessMath);
	TURPointXY pntSdvig(valCentreVessX, 0.);
	TURPolygon PlgVessMath1 = PlgVessMath.LinTransform( valAngPsi,  pntSdvig  , valRastigenie );
	if (fabs(valAngPsi)> M_PI/ 2.)
  {
   PlgVessMath1 = PlgVessMath1.SimOtragenieTransform();
  }

	PlgVessMath1.WriteSetSHPFiles(pwcharrFileVess, &PlgVessMath1,1);

   wcscpy(pwcharrFileVess, pwcharrPath);
   wcscat(pwcharrFileVess, L"\\CoordMatveev.shp");
   TURPolygon PlgCoord(NumPartsCoord, NumPointsCoord, PartsCoord, PointsCoord);
   TURPointXY pntCentre(centreCoordX, centreCoordY);
   pntCentreNew =  pntCentre.LinTransform( valAngPsi,  pntSdvig  , valRastigenie );
  // TURPointXY pntSdvig(valCentreVessX, 0.);
   TURPolygon PlgCoord1 = PlgCoord.LinTransform( valAngPsi ,  pntSdvig  , valRastigenie );

   TURPointXY pntSdvig1( -pntCentreNew.X, -pntCentreNew.Y);
	TURPolygon PlgCoord2 = PlgCoord1.SdvigTransform( pntSdvig1 );

	TURPointXY pntSdvig00(0.,0.);
	TURPolygon PlgCoord3 = PlgCoord2.LinTransform(valAngEps ,  pntSdvig00,1. ) ;

	TURPolygon PlgCoord4 = PlgCoord3.SdvigTransform( pntCentreNew);
	if (fabs(valAngPsi)> M_PI/ 2.)
  {
   PlgCoord4 = PlgCoord4.SimOtragenieTransform();
  }

   PlgCoord4.WriteSetSHPFiles(pwcharrFileVess, &PlgCoord4,1);

	// Создание диаграммы с центром в точке  pointCentre
	// радиусом  valR с шагом по углу iNUM
	// valWidthMilrad -  ширина диаграммы в милирадианах
	// valRotAng - угол поворота относительно горизонтали
	// положит направление против час стрелки
	double valWidth = 4./180. * M_PI ;
	const double valR = 4000.;
	const double valRotAng = valAngEps;//-2./180. * M_PI ;
	const int iNUM = 3000;
	TURPolygon Diagr = TURPolygon::fncCreateDiagr(pntCentreNew,  valR,valRotAng , valWidth, iNUM) ;
	wcscpy(pwcharrFileVess, pwcharrPath);
	wcscat(pwcharrFileVess, L"\\Diagr.shp");
	Diagr.WriteSetSHPFiles(pwcharrFileVess, &Diagr ,1);

	// создание шейп файла с центральной осью диаграммы
	pntEngDiagrAxe = TURPointXY (pntCentreNew.X + valR * cos(valRotAng), pntCentreNew.Y + valR * sin(valRotAng));


}


// перегруженная функция, дополнительно на выходе выдающая центр координатора  и конечную точку центральной оси диаграммы
// плюс на входе предусмотрена дальность диаграммы и угол диаграммы
void drawVesselMathewPicture( wchar_t *pwcharrPath,const double  valAngPsi,const double  valAngEps
	 , const double  valCentreVessX,const double valRastigenie, double valR, double valWidth, TURPointXY &pntCentreNew, TURPointXY &pntEngDiagrAxe)
{
	wchar_t pwcharrFileVess[300] ={0};
	wcscpy(pwcharrFileVess, pwcharrPath);
	wcscat(pwcharrFileVess, L"\\VessMatveev.shp");
	TURPolygon PlgVessMath(NumPartsVessMath, NumPointsVessMath, PartsVessMath, PointsVessMath);
	TURPointXY pntSdvig(valCentreVessX, 0.);
	TURPolygon PlgVessMath1 = PlgVessMath.LinTransform( valAngPsi,  pntSdvig  , valRastigenie );
	if (fabs(valAngPsi)> M_PI/ 2.)
  {
   PlgVessMath1 = PlgVessMath1.SimOtragenieTransform();
  }

	PlgVessMath1.WriteSetSHPFiles(pwcharrFileVess, &PlgVessMath1,1);

   wcscpy(pwcharrFileVess, pwcharrPath);
   wcscat(pwcharrFileVess, L"\\CoordMatveev.shp");
   TURPolygon PlgCoord(NumPartsCoord, NumPointsCoord, PartsCoord, PointsCoord);
   TURPointXY pntCentre(centreCoordX, centreCoordY);
   pntCentreNew =  pntCentre.LinTransform( valAngPsi,  pntSdvig  , valRastigenie );
  // TURPointXY pntSdvig(valCentreVessX, 0.);
   TURPolygon PlgCoord1 = PlgCoord.LinTransform( valAngPsi ,  pntSdvig  , valRastigenie );

   TURPointXY pntSdvig1( -pntCentreNew.X, -pntCentreNew.Y);
	TURPolygon PlgCoord2 = PlgCoord1.SdvigTransform( pntSdvig1 );

	TURPointXY pntSdvig00(0.,0.);
	TURPolygon PlgCoord3 = PlgCoord2.LinTransform(valAngEps ,  pntSdvig00,1. ) ;

	TURPolygon PlgCoord4 = PlgCoord3.SdvigTransform( pntCentreNew);
	if (fabs(valAngPsi)> M_PI/ 2.)
  {
   PlgCoord4 = PlgCoord4.SimOtragenieTransform();
  }

   PlgCoord4.WriteSetSHPFiles(pwcharrFileVess, &PlgCoord4,1);

	// Создание диаграммы с центром в точке  pointCentre
	// радиусом  valR с шагом по углу iNUM
	// valWidthMilrad -  ширина диаграммы в милирадианах
	// valRotAng - угол поворота относительно горизонтали
	// положит направление против час стрелки

	const double valRotAng = valAngEps;//-2./180. * M_PI ;
	const int iNUM = 3000;
	TURPolygon Diagr = TURPolygon::fncCreateDiagr(pntCentreNew,  valR,valRotAng , valWidth, iNUM) ;
	wcscpy(pwcharrFileVess, pwcharrPath);
	wcscat(pwcharrFileVess, L"\\Diagr.shp");
	Diagr.WriteSetSHPFiles(pwcharrFileVess, &Diagr ,1);

	// создание шейп файла с центральной осью диаграммы
	pntEngDiagrAxe = TURPointXY (pntCentreNew.X + valR * cos(valRotAng), pntCentreNew.Y + valR * sin(valRotAng));


}


void drawVesselMathewPicture_With_SK( wchar_t *pwcharrPath,const double  valAngPsi,const double  valAngEps
	 ,const double  valTargU,const double  valTargR)
{
  drawVesselMathewPicture( pwcharrPath,  valAngPsi, valAngEps , 0.,1.) ;

  const TURPointXY pointBeginX(-40., 0.);
  const TURPointXY pointEndX(40., 0.);
  const TURPointXY pointBeginY(0., -40.);
  const TURPointXY pointEndY (0., 40.);
  const double valLength = 5.;
  TURPolyLine AxesKGSK = TURPolyLine::fncCreateAxes( pointBeginX,  pointEndX
									   ,pointBeginY, pointEndY, valLength) ;
  wchar_t pwcharrKGSK[300] ={0};
  wcscpy(pwcharrKGSK, pwcharrPath);
  wcscat(pwcharrKGSK, L"\\KGSK.shp");
  AxesKGSK.WriteSetSHPFiles(pwcharrKGSK, &AxesKGSK ,1);


  const TURPointXY pntSdvig00(0.,0.);
  TURPolyLine AxesPSK = AxesKGSK.LinTransform(valAngPsi , pntSdvig00,1. ) ;
  wcscpy(pwcharrKGSK, pwcharrPath);
  wcscat(pwcharrKGSK, L"\\PSK.shp");
  AxesPSK.WriteSetSHPFiles(pwcharrKGSK, &AxesPSK ,1);

  // центр координатора
   TURPointXY pntCentre(centreCoordX, centreCoordY);
   TURPointXY pntCentreNew =  pntCentre.LinTransform( valAngPsi,  pntSdvig00  , 1. );
   ///

   // ПСК-АС
   TURPolyLine AxesPSK_AS = AxesPSK.SdvigTransform(pntCentreNew ) ;
   wcscpy(pwcharrKGSK, pwcharrPath);
   wcscat(pwcharrKGSK, L"\\PSK-AS.shp");
   AxesPSK_AS.WriteSetSHPFiles(pwcharrKGSK, &AxesPSK_AS ,1);
  ///

  ///АСК
 TURPolyLine AxesASK = AxesPSK.LinTransform(valAngEps , pntCentreNew,1. ) ;
   wcscpy(pwcharrKGSK, pwcharrPath);
   wcscat(pwcharrKGSK, L"\\ASK.shp");
  AxesASK.WriteSetSHPFiles(pwcharrKGSK, &AxesASK ,1);

  ///
   ///ЛСК
  TURPolyLine AxesLSK = AxesPSK.LinTransform(valAngEps + valTargU , pntCentreNew,1. ) ;
   wcscpy(pwcharrKGSK, pwcharrPath);
   wcscat(pwcharrKGSK, L"\\LSK.shp");
   AxesLSK.WriteSetSHPFiles(pwcharrKGSK, &AxesLSK ,1);
   ///

   // диаграмма
   double valFiDiagr = 12./180. * M_PI ;
   TURPolygon Diagr = TURPolygon::fncCreateSector(pntCentreNew, valTargR * 1.3,
					valAngPsi +valAngEps + valTargU -valFiDiagr/2.
	   ,  valAngPsi +valAngEps + valTargU +valFiDiagr/2., 500 );

	 wcscpy(pwcharrKGSK, pwcharrPath);
   wcscat(pwcharrKGSK, L"\\Diagr.shp");
	   Diagr.WriteSetSHPFiles(pwcharrKGSK, &Diagr ,1);

  // угл сектор
  wcscpy(pwcharrKGSK, pwcharrPath);
   wcscat(pwcharrKGSK, L"\\Sect1.shp");
 TYrWriteShapeFile::CreateAngleMarks(pwcharrKGSK, pntSdvig00
	 , AxesKGSK.Points[1], AxesPSK.Points[1],20.,22. ) ;
	 ///

 // угл сектор
  wcscpy(pwcharrKGSK, pwcharrPath);
   wcscat(pwcharrKGSK, L"\\Sect2.shp");
 TYrWriteShapeFile::CreateAngleMarks(pwcharrKGSK, pntCentreNew
	 , AxesPSK_AS.Points[1], AxesLSK.Points[1],10.,12. ) ;
	 ///

	 // ВС цели
	 double arrS_PSK[2] ={0.},arrS_KGSK[2] ={0.};
	 arrS_PSK[0] = valTargR *  cos (valAngPsi +valAngEps + valTargU - 2.* M_PI/ 180.);
	 arrS_PSK[1] = valTargR *  sin (valAngPsi +valAngEps + valTargU- 2.* M_PI/ 180.);
	 arrS_KGSK[0] = arrS_PSK[0] + pntCentreNew.X;
	 arrS_KGSK[1] = arrS_PSK[1] + pntCentreNew.Y;
	 TURPointXY pntSdvigTarg(arrS_KGSK[0], arrS_KGSK[1]);
	 double tetta = atan2(arrS_KGSK[1], arrS_KGSK[0]);
	 wcscpy(pwcharrKGSK, pwcharrPath);
	 wcscat(pwcharrKGSK, L"\\Target.shp");
	 drawUAB( pwcharrKGSK,tetta + M_PI, pntSdvigTarg,2.)  ;
	 ///

	 // вектор параллакса
	 TURPolyLine Arrow0 = TURPolyLine::fncCreateArrow(pntSdvig00, pntCentreNew
					,4,5./180. * M_PI);
	 TURPolyLine Arrow1 = TURPolyLine::fncCreateArrow(pntCentreNew, pntSdvig00
					,4, 5./180. * M_PI);
	 wcscpy(pwcharrKGSK, pwcharrPath);
	 wcscat(pwcharrKGSK, L"\\VectParall0.shp");
	 Arrow0.WriteSetSHPFiles(pwcharrKGSK, &Arrow0 ,1);

	 wcscpy(pwcharrKGSK, pwcharrPath);
	 wcscat(pwcharrKGSK, L"\\VectParall1.shp");
	 Arrow1.WriteSetSHPFiles(pwcharrKGSK, &Arrow1 ,1);

	 // ГСК
  const TURPointXY pntSdvig111(-140.,0.);
  TURPolyLine AxesGSK = AxesKGSK.LinTransform(0. , pntSdvig111,1. ) ;
  wcscpy(pwcharrKGSK, pwcharrPath);
  wcscat(pwcharrKGSK, L"\\GSK.shp");
  AxesGSK.WriteSetSHPFiles(pwcharrKGSK, &AxesGSK ,1);

}
#pragma package(smart_init)
