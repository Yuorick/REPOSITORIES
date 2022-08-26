//---------------------------------------------------------------------------


#pragma hdrstop
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "SeaTargPicturte.h"
#include "YrWriteShapeFile.h"
#include "URPolygon.h"
#include "URPointXY.h"
#include "MatrixProccess.h"
#include "URPolyLine.h"
#include "ArcEllipse.h"
#include "Target.h"
#include "Gauss.h"
#include "PlanePolygon.h"

// построение картинки для береговой цели
// wchOutFold - папака с шейп файлами
// arrCoMtrx_1_and_3_Groups [36] - коррел матрица рассеяния точки падения  1 и 3 группы
// arrCoMtrx_2_Group[36] - коррел матрица рассеяния точки падения  2 группа
// arrVectMiss_GSK [6] - средний векторт промаха
// Target - цель
//
//

void createNibourAppointmPointPictureForSeaTarg(wchar_t *wchOutFold
	,double *arrCoMtrx_1_and_3_Groups, double *arrCoMtrx_2_Group,double *arrVectMiss_GSK
	, const TTarget Target, const int QuantShells)
{
	const TURPointXY pntSdvig(0.,0.);
	double arrCorrMatrxVectMiss_GSK[36] = {0.};
	MtrxSumMatrx(arrCoMtrx_1_and_3_Groups, arrCoMtrx_2_Group, 6, 6, arrCorrMatrxVectMiss_GSK);
	wchar_t wchFoldName[300] ={0}, wchFileNamePolygonEll1[300] ={0}, wchFileNamePolygonEll2[300] ={0}
	, wchFileNamePolygonEll3[300] ={0}
	, wchFileArrow[300] ={0}, wchFilePoints[300] ={0}, wchTargShadows[300] = {0}
	,  wchFileNamePolygonEll_1_3_Groups[300] ={0}, wchFileNamePolygonEll_2_Group[300] ={0};

	wcscpy(  wchFoldName,  wchOutFold);
	wcscat(wchFoldName, L"\\");

	// 1. курсовой угол снаряда:
	const	double VAlPsi = atan2(arrVectMiss_GSK[4], arrVectMiss_GSK[3]);
	///
	// 2. Формирование матрицы эллипса рассеяния точки падения 1  и 3 группы в плоскости OXY ГСК

	TURPolygon plgonEll_1_3_Groups = createEllipsTochkiPadenia(arrCoMtrx_1_and_3_Groups, pntSdvig, 1.);
	wcscpy(  wchFileNamePolygonEll_1_3_Groups,  wchFoldName);
	wcscat( wchFileNamePolygonEll_1_3_Groups, L"plgonEll_1_3_Groups.shp");
	plgonEll_1_3_Groups.WriteSetSHPFiles( wchFileNamePolygonEll_1_3_Groups, &plgonEll_1_3_Groups,1);

	///

		// 2. Формирование матрицы эллипса рассеяния точки падения 2 группы в плоскости OXY ГСК
	TURPolygon plgonEll_2_Group = createEllipsTochkiPadenia(arrCoMtrx_2_Group, pntSdvig, 1.);
	wcscpy(  wchFileNamePolygonEll_2_Group,  wchFoldName);
	wcscat( wchFileNamePolygonEll_2_Group, L"plgonEll_2_Group.shp");
	plgonEll_2_Group.WriteSetSHPFiles( wchFileNamePolygonEll_2_Group, &plgonEll_2_Group,1);

	///

	// 2. Формирование матрицы эллипса рассеяния точки падения в плоскости OXY ГСК

	TURPolygon plgonEll1 = createEllipsTochkiPadenia(arrCorrMatrxVectMiss_GSK, pntSdvig, 1.);
	wcscpy(  wchFileNamePolygonEll1,  wchFoldName);
	wcscat( wchFileNamePolygonEll1, L"PolygonEll1.shp");
	plgonEll1.WriteSetSHPFiles( wchFileNamePolygonEll1, &plgonEll1,1);

	///

	// 3. формирование эллипса разбросов по уровню 2
	TURPolygon plgonEll2 = createEllipsTochkiPadenia(arrCorrMatrxVectMiss_GSK, pntSdvig, 2.);

	wcscpy(  wchFileNamePolygonEll2,  wchFoldName);
	wcscat( wchFileNamePolygonEll2, L"PolygonEll2.shp");

	plgonEll2.WriteSetSHPFiles( wchFileNamePolygonEll2, &plgonEll2,1);

	///

	// 4. формирование эллипса разбросов по уровню 3
	TURPolygon plgonEll3 = createEllipsTochkiPadenia(arrCorrMatrxVectMiss_GSK, pntSdvig, 3.);

	wcscpy(  wchFileNamePolygonEll3,  wchFoldName);
	wcscat( wchFileNamePolygonEll3, L"PolygonEll3.shp");

	plgonEll3.WriteSetSHPFiles( wchFileNamePolygonEll3, &plgonEll3,1);

	///

	// 5. проекция вектора скорости снаряда на гориз. плоскость
	double val_t  = -100.;
	TURPointXY pointBegin (val_t * cos(VAlPsi), val_t * sin(VAlPsi) );
	TURPointXY pointEnd(0.,0.);
	TURPolyLine plnArrow = TURPolyLine::fncCreateArrow( pointBegin,  pointEnd
	,10.,20. * M_PI / 180.);

	wcscpy( wchFileArrow,  wchFoldName);
	wcscat( wchFileArrow, L"ShellVelo.shp");
	plnArrow.WriteSetSHPFiles( wchFileArrow, &plnArrow,1);
	///

	// оси  координат
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-500., 500.
	,-500., 500.,30.) ;


	/////////////////////////////////////////////////////////////////////////////////

	// 7. Точки падения
	double arrElK[4] = {0.};
	arrElK[0] = arrCorrMatrxVectMiss_GSK[0] ;
	arrElK[1] = arrCorrMatrxVectMiss_GSK[1] ;
	arrElK[2] = arrElK[1] ;
	arrElK[3] = arrCorrMatrxVectMiss_GSK[7] ;

	//   arrF - матрица собственных векторов  коррел матрицы
	double arrF[4] = {0.} , arrMtrxLamb[4] = {0.};
	CalcProperVectors2(arrElK, arrF , arrMtrxLamb) ;
	TURPointXY *pPntArr = new TURPointXY[QuantShells];
	double  arrPos_GSK[2] = {0.};
	for (int i =0; i < QuantShells; i++)
	{
	getGaussVector(2, arrVectMiss_GSK, arrF , arrMtrxLamb, arrPos_GSK);
	pPntArr[i].X = arrPos_GSK[0];
	pPntArr[i].Y = arrPos_GSK[1];

	}

	wcscpy(  wchFilePoints,  wchFoldName);
	wcscat( wchFilePoints, L"Points.shp");
	pPntArr[0].WriteSetSHPFiles( wchFilePoints, pPntArr,QuantShells);
	delete []pPntArr;
	//////////////////////////////////////////////////////////////////////////////
 /*	 double arrTargV1 [3] = {0.}, arrVMiss[3] = {0.};
	 double valTargCourse = 180. * M_PI / 180.;
	 arrTargV1[0] = 10. * sin(valTargCourse);
	 arrTargV1 [1] = 10. * cos(valTargCourse);

	 double valTet = -45. * M_PI /180.;
	 double valBet = -90. * M_PI / 180.;
	 double valV0 = 300.;
	 arrVMiss[0] = valV0 * cos (valTet) * sin ( valBet);
	 arrVMiss[1] = valV0 * cos (valTet) * cos ( valBet);
	 arrVMiss[2] = valV0 * sin (valTet) ;

*/
	///////////////////////////////////////
	TURPolygon *pPlgArrShadows  = new  TURPolygon[Target.mLenArrPlanePolygon];
	for (int i = 0; i < Target.mLenArrPlanePolygon; i++)
	{
	TPlanePolygon PlanePolygon =  Target.mpArrPlanePolygon[i];
 	double arrTargV[3] = {0.};
	memcpy(arrTargV, &(Target.mTraject.marrVectSostGSK[3]),3 *sizeof(double));
	pPlgArrShadows[i] =  PlanePolygon.createShadowPlg_For_PlaneBody(arrTargV,  &arrVectMiss_GSK[3]);
	}
	wcscpy(  wchTargShadows,  wchFoldName);
	wcscat( wchTargShadows, L"TargShadows.shp");
	pPlgArrShadows[0].WriteSetSHPFiles( wchTargShadows, pPlgArrShadows,Target.mLenArrPlanePolygon);
	delete []pPlgArrShadows;


}


TURPolygon createEllipsTochkiPadenia(double *arrCorrMatrxVectMiss_GSK, const TURPointXY PNtSdvig, const double VAlLevel)
{
	double arrElK[4] = {0.};
	arrElK[0] = arrCorrMatrxVectMiss_GSK[0] ;
	arrElK[1] = arrCorrMatrxVectMiss_GSK[1] ;
	arrElK[2] = arrElK[1] ;
	arrElK[3] = arrCorrMatrxVectMiss_GSK[7] ;

	//   arrF - матрица собственных векторов  коррел матрицы
	double arrF[4] = {0.} , arrMtrxLamb[4] = {0.};
	CalcProperVectors2(arrElK, arrF , arrMtrxLamb) ;

	///

	// матрица линейного преобразования, формирующего вектор   разбросов координат точки падения
	// из двухмерного вектора некореллированных координат единичной дисперсии
	double arrLinTrasf [4] ={0.}, arrLambSq[4] ={0.};
	arrLambSq[0] =  sqrt(arrMtrxLamb[0]);
	arrLambSq[3] =  sqrt(arrMtrxLamb[3]);
	MtrxMultMatrx(arrF,2, 2, arrLambSq,2, arrLinTrasf) ;
	///

	// 2. формирование эллипса разбросов по уровню 1

	TURPolygon plgCircle1 = TURPolygon::fncCreateCircle(PNtSdvig,VAlLevel, 1001) ; // это единичный круг
	TURPolygon plgonEll1 = plgCircle1.fncLinTransform(arrLinTrasf );// это его линейное преобразование
	return plgonEll1;
}

	///

#pragma package(smart_init)
