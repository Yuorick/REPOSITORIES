//---------------------------------------------------------------------------


#pragma hdrstop
#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "URPolygon.h"


//---------------------------------------------------------------------------
// создание шейп файла графика функции из массива parrInf
// INPUT:
// wchFileName -  имя файла
// parrInf  - массив значений функции
// parrT - массив соответствующих времен
// lenarr - длина этих массивов
// scalex - масштаб по оси X
// scaley - масштаб по оси Y

void __fastcall TYrWriteShapeFile::CreateShpFile(wchar_t *wchFileName, double *parrInf, double *parrT
	 ,const int lenarr, double &scalex, double &scaley)
{

	TURPointXY *pPoints = new TURPointXY[lenarr ];
	for (int i = 0; i < lenarr ; i++)
	{
	  pPoints [ i ].Y = parrInf[i]*scaley ;
	  pPoints [ i ].X = parrT[ i]*scalex ;
	}
	//
	int numParts   = 1  ;
	int numPoints  = lenarr ;
	int *iarrParts = new int[ numParts] ;
	iarrParts[0] = 0 ;

	TURPolyLine Plyline( numParts, numPoints ,iarrParts, pPoints);
	TURPolyLine::WriteSetSHPFiles(wchFileName,&Plyline, 1) ;

	delete [] pPoints;
	delete  iarrParts;
}
double TYrWriteShapeFile::maxDoubleArr(double *parr, const int lenarr, int &irez)
{
  double rez = -DBL_MAX ;
  for (int i =0; i < lenarr; i++)
  {
	  if (parr[i] > rez)
	  {
		rez = parr[i];
		irez = i ;
	  }
  }
  return rez ;
}
double TYrWriteShapeFile::minDoubleArr(double *parr, const int lenarr, int &irez)
{
  double rez = DBL_MAX ;
  for (int i =0; i < lenarr; i++)
  {
	  if (parr[i] < rez)
	  {
		rez = parr[i];
		irez = i ;
	  }
  }
  return rez ;
}

void TYrWriteShapeFile::WriteOneReport(wchar_t *wcharrPath1  // путь к папке
								  ,double * parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,const int nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  ,const int nBuffRows //  - к-во точек
								  ,wchar_t *wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,const int lenName // максимальная длина имени переменной
								  ,const int numx  // номер переменной по оси X
								  ,const int numy  // номер переменной по оси Y
								  ,const double scalex  //  масштаб по оси X
								  ,const double scaley  // масштаб по оси Y
								   )
{
	wchar_t wcharrPath [300] = {0};
	wcscpy(wcharrPath, wcharrPath1);
	if (wcharrPath[wcslen(wcharrPath) -1] != L'\\')
	{
		wcharrPath[wcslen(wcharrPath)] =  L'\\';
	}
  wchar_t wcharrFileName[300] = {0} ;
	memset( wcharrFileName, 0, sizeof(wchar_t) * 300) ;
  // формирование имени файла(полного пути с названием)
  createFileName( wcharrPath  // путь к папке
				,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
				, lenName // максимальная длина имени переменной
				, numx  // номер переменной по оси X
				,numy  // номер переменной по оси Y
				,wcharrFileName  // w строка с именем
				);
  //создание массива точек для фозданпия полтлтнтт графика
  TURPointXY *pPoints = new TURPointXY[nBuffRows ];
	for (int i = 0; i < nBuffRows ; i++)

	{

	  pPoints [ i ].Y = parrBuff[ nBuffCols * i + numy ]*scaley ;
	  pPoints [ i ].X = parrBuff[ nBuffCols * i + numx ]*scalex ;
	}

   	//
	int numParts   = 1  ;
	int numPoints  = nBuffRows ;
	int *iarrParts = new int[ numParts] ;
	iarrParts[0] = 0 ;

	TURPolyLine Plyline( numParts, numPoints ,iarrParts, pPoints);
	TURPolyLine::WriteSetSHPFiles(wcharrFileName,&Plyline, 1) ;

  delete  [] pPoints ;
  delete  iarrParts;

}

void TYrWriteShapeFile::createFileName( wchar_t *wcharrPath  // путь к папке
										,wchar_t *wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
										,const int lenName // максимальная длина имени переменной
										,const int numx  // номер переменной по оси X
										,const int numy  // номер переменной по оси Y
										,wchar_t *wcharrFileName  // w строка с именем
										)
{
   wcscpy(wcharrFileName, wcharrPath);
   wcscat(wcharrFileName, &wcharrFileNames[lenName * numy]);
   wcscat(wcharrFileName, L"( ");
   wcscat(wcharrFileName, &wcharrFileNames[lenName * numx]);
   wcscat(wcharrFileName, L" )");
   wcscat(wcharrFileName, L".shp");
}

// создание shape файла с осями координат
void  TYrWriteShapeFile::CreateShpAxes(wchar_t *wchFileName,const double xmin,const double xmax
	 ,const double ymin,const double ymax)

{
	int numParts   =  2 ;
	int numPoints  =  4 ;
	int iarrParts[2] = {0,2};
	 TURPointXY arrPoints[4] ;
	 arrPoints[0] =  TURPointXY(xmin,0);
	 arrPoints[1] =  TURPointXY(xmax,0);
	 arrPoints[2] =  TURPointXY(0,ymin);
	 arrPoints[3] =  TURPointXY(0,ymax);

	TURPolyLine Plyline( numParts, numPoints ,iarrParts, arrPoints);
	TURPolyLine::WriteSetSHPFiles(wchFileName,&Plyline, 1) ;
}

// Создание  осями координат  со стрелками
void TYrWriteShapeFile::CreateShpArrowedAxes(wchar_t *wchFileName,const double xmin,const double xmax
	 ,const double ymin,const double ymax,const double valLength)
{
  const TURPointXY pointBeginX(xmin, 0.);
  const TURPointXY pointEndX (xmax, 0.);
  const TURPointXY pointBeginY(0., ymin);
  const TURPointXY pointEndY(0., ymax);
  TURPolyLine Axes = TURPolyLine::fncCreateAxes( pointBeginX, pointEndX
									   , pointBeginY,  pointEndY, valLength);
  TURPolyLine::WriteSetSHPFiles(wchFileName,&Axes, 1) ;
}

// Создание  осями координат  со стрелками  c началом координат в точке pnt00
void TYrWriteShapeFile::CreateShpArrowedAxes(wchar_t *wchFileName,const double xmin,const double xmax
	 ,const double ymin,const double ymax,const double valLength, const TURPointXY pnt00)
{
  const TURPointXY pointBeginX(xmin, 0.);
  const TURPointXY pointEndX (xmax, 0.);
  const TURPointXY pointBeginY(0., ymin);
  const TURPointXY pointEndY(0., ymax);
  TURPolyLine Axes1 = TURPolyLine::fncCreateAxes( pointBeginX, pointEndX
									   , pointBeginY,  pointEndY, valLength);
  TURPolyLine  Axes =  Axes1.SdvigTransform(pnt00 ) ;
  TURPolyLine::WriteSetSHPFiles(wchFileName,&Axes, 1) ;
}

// Создание  полилинии с секторами показывающими угол
//Pnt0 -  вершина угла
//Pnt0 - Pnt1 - строна угла
// Pnt0 - Pnt1  - строна угла
//  Dist0-  расстояние до первого сектора
//  Dist1 -  расстояние до второго сектора
void TYrWriteShapeFile::CreateAngleMarks(wchar_t *wchFileName, const TURPointXY Pnt0
	 , const TURPointXY Pnt1, const TURPointXY Pnt2,const double Dist0,const double Dist1 )
{
	double valFi0 = 0.;
	if (fabs(Pnt1.X - Pnt0.X) < 0.00001)
	{
	  valFi0 = Sign_( Pnt1.Y - Pnt0.Y)*M_PI/2.;
	}
	else
	{
	  valFi0 = atan2( Pnt1.Y - Pnt0.Y,  Pnt1.X - Pnt0.X);
	}
	///
	double valFi1 = 0.;
	if (fabs(Pnt2.X - Pnt0.X) < 0.00001)
	{
	  valFi1 = Sign_( Pnt2.Y - Pnt0.Y)*M_PI/2.;
	}
	else
	{
	  valFi1 = atan2( Pnt2.Y - Pnt0.Y,  Pnt2.X - Pnt0.X);
	}
	TURPolyLine Sect0 = TURPolyLine::fncCreateSector(Pnt0, Dist0,valFi0, valFi1,1500) ;
	TURPolyLine Sect1 = TURPolyLine::fncCreateSector(Pnt0, Dist1,valFi0, valFi1,1500) ;

	int iarrParts[2] = {0, 1500};
	const int iNumPoints = 3000;
	const int iNumParts = 2;
	TURPointXY arrPoints[3000];
	memcpy(arrPoints, Sect0.Points, 1500 * sizeof (TURPointXY));
	memcpy(&arrPoints[1500], Sect1.Points, 1500 * sizeof (TURPointXY));
	TURPolyLine plnRez( iNumParts, iNumPoints,iarrParts,arrPoints) ;

  TURPolyLine::WriteSetSHPFiles(wchFileName,&plnRez, 1) ;
}

void  TYrWriteShapeFile::ShowNormProbDistr(wchar_t *wchFileName,const double valA,const double valSigm2)
{
	const int N = 1000;
	double xmin =  valA - 5. * sqrt( valSigm2);
	double xmax =  valA + 5. * sqrt( valSigm2);
	double valDel = (xmax- xmin)/ ((double)N) ;
	int numParts = 1;
	int numPoints = N ;
	double xt = xmin;
	TURPolyLine Plyline( numParts, numPoints );
	for (int i = 0; i < N; i++)
	{
	  Plyline.Points[i].X = xt;
	  Plyline.Points[i].Y =  exp(- (xt -valA) * (xt - valA) / 2./valSigm2)/sqrt( 2. * M_PI *valSigm2) ;
	  xt += valDel ;
	}
		TURPolyLine::WriteSetSHPFiles(wchFileName,&Plyline, 1) ;
}

double  TYrWriteShapeFile::Sign_( const double x)
{
	if (x >0. )
	{
	  return 1. ;
	}
	else
	{
		if (x <0.)
		{
		 return -1.;
		}
	}
	return 0.;
}

// картиека с ФАР
void  TYrWriteShapeFile::PictFar(wchar_t *wchFoldName)
{
 // первые 4 точки - излучатели
 // arrPoints[4] - начало координат - фазовый центр

 TURPointXY arrPoints[5];
 arrPoints[0] = TURPointXY(0., -3./2.) ;
 arrPoints[1] = TURPointXY(0., -1./2.) ;
 arrPoints[2] = TURPointXY(0., 1./2.) ;
 arrPoints[3] = TURPointXY(0., 3./2) ;
 arrPoints[4] = TURPointXY(100.,0.);
 // угол цели
  double valAl = 10./ 180. * M_PI;
 // линия фазового фронта проходящая через фазовый центр
 TURPointXY pntPhaseCntr(0.,0.);
 double valxfront = 6.;
 TURPointXY pntFront0( pntPhaseCntr.X -valxfront*sin(valAl), pntPhaseCntr.Y +valxfront*cos(valAl) );
 TURPointXY pntFront1 (-pntFront0.X, -pntFront0.Y);
 TURPolyLine plnFront(pntFront0,pntFront1);
 // точки пересечения лучей с линией фазового фронта
 TURPointXY arrPointFront[5];
 for (int i = 0 ; i < 4; i++)
 {
  arrPointFront[i] = TURPointXY(pntPhaseCntr.X - arrPoints[i].Y * cos(valAl) * sin(valAl)
	  ,pntPhaseCntr.X + arrPoints[i].Y * cos(valAl) * cos(valAl));
 }
 arrPointFront[4] = pntPhaseCntr;

 // линии с сигналом от цели проходят через точки arrPoints
 TURPolyLine arrPlyline[5], arrSegmSdvig[5];

 double valx =10.;
 for (int i = 0; i < 5; i++)
 {
   TURPointXY  pnt2(valx,  valx* tan(valAl) + arrPoints[i].Y);
   arrPlyline[i] = TURPolyLine(arrPointFront[i],   pnt2) ;
   arrSegmSdvig[i] = TURPolyLine(arrPoints[i], arrPointFront[i]);
 }

 // нормаль к ФАР в точке фазового центра
 TURPointXY pntNorm(10.,0.);
 TURPolyLine plnNorm(pntNorm,pntPhaseCntr) ;
///

// массив антенных модулей
TURPolygon arrPlgAM[4];
for (int i =0; i < 4; i++)
{
 TURPointXY pntTopRight( 0.,arrPoints[i].Y + 0.5 );
 TURPointXY pntBottomLeft( -0.05,arrPoints[i].Y - 0.5 );
 arrPlgAM[i] = TURPolygon (pntTopRight, pntBottomLeft);
}

 // линии-нормали
 /*TURPolyLine arrNormPlyline[5];

 for (int i = 0; i < 4; i++)
 {
   double d = arrPoints[i+1].Y -arrPoints[i].Y;
   TURPointXY  pnt2(d * sin(valAl)* cos(valAl),  arrPoints[i].Y + d * cos(valAl)* cos(valAl));
   arrNormPlyline[i] = TURPolyLine(arrPoints[i],   pnt2) ;
 }
   double d = arrPoints[3].Y -arrPoints[1].Y;
   TURPointXY  pnt2(d * sin(valAl)* cos(valAl),  arrPoints[1].Y + d * cos(valAl)* cos(valAl));
   arrNormPlyline[4] = TURPolyLine(arrPoints[1],   pnt2) ; */
   // файлы
  wchar_t pwcharrFile[300] ={0};
  wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\Points.shp");
  TURPointXY::WriteSetSHPFiles(pwcharrFile,arrPoints, 4) ;

  wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\Arrays.shp");
  TURPolyLine::WriteSetSHPFiles(pwcharrFile,arrPlyline, 5) ;

	wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\plnNOrm.shp");
  TURPolyLine::WriteSetSHPFiles(pwcharrFile,&plnNorm, 1) ;

   TURPolyLine plnFar( arrPoints[0], arrPoints[4]);
	wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\Far.shp");
  TURPolygon::WriteSetSHPFiles(pwcharrFile,arrPlgAM, 4) ;

  // обозначение угла между нормалью и лучом
  wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\Sect0.shp");
  TURPointXY  pnt0(0.,0.);
  double Dist0 = 0.8, Dist1 =0.85;
  CreateAngleMarks(pwcharrFile, pnt0
	 , pntNorm, arrPlyline[4].Points[1],Dist0, Dist1 ) ;
	 ///

 // обозначение угла между фронтом и ФАР
  wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\Sect1.shp");


  CreateAngleMarks(pwcharrFile, pnt0
	 , arrPoints[0], pntFront1,Dist0, Dist1 ) ;
	 ///

	 // обозначение угла между фронтом и ФАР
  wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\Sect2.shp");


  CreateAngleMarks(pwcharrFile, pnt0
	 , arrPoints[3], pntFront0,Dist0, Dist1 ) ;
	 ///

  //TURPolyLine plnNorm( pnt0, arrPoints[5]);
   //	wcscpy(pwcharrFile, wchFoldName);
 // wcscat(pwcharrFile, L"\\plnNorm.shp");
 // TURPolyLine::WriteSetSHPFiles(pwcharrFile,& plnNorm, 1) ;

  // линия фронта волны
	wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\plnFront.shp");
  TURPolyLine::WriteSetSHPFiles(pwcharrFile,&plnFront, 1) ;
  // отрезки фазоовых сдвигов

	wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\sgmSdvig.shp");
  TURPolyLine::WriteSetSHPFiles(pwcharrFile,arrSegmSdvig, 4) ;
  // фазовый центр
  	wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\Centre.shp");
  TURPointXY::WriteSetSHPFiles(pwcharrFile,&pntPhaseCntr, 1) ;
  ///


}
#pragma package(smart_init)
