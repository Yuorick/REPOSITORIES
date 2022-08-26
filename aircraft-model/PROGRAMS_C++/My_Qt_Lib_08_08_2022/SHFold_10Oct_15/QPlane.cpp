//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include  <string.h>
#include  <float.h>
#include  <math.h>
#include "QPlane.h"
#include "MatrixProccess.h"
#include <time.h>

//--------------------------------------------------------------------------------------
 QPlane ::QPlane()
{
	double arrT [9] = {0.};
	memcpy(marrS0, arrT, 3 * sizeof(double));
	memcpy(marrF, arrT, 9 * sizeof(double));
    marrF[0] = marrF[4] = marrF[8] = 1.;
	
}

 // оператор присваивания
 QPlane &QPlane::operator=(const QPlane  &R)
 {
	memcpy(marrS0, R.marrS0, 3 * sizeof(double));
	memcpy(marrF, R.marrF, 9 * sizeof(double));
	return *this ;
 }

 // конструктор копирования
 QPlane::QPlane (const QPlane &R)
 {
	memcpy(marrS0, R.marrS0, 3 * sizeof(double));
	memcpy(marrF, R.marrF, 9 * sizeof(double));
 }


 // парам констр
QPlane ::QPlane( double* arrS0, double* arrF)
{
	memcpy(marrS0, arrS0, 3 * sizeof(double));
	memcpy(marrF, arrF, 9 * sizeof(double));

}
//-------------------------------------------------------------
bool QPlane::findIntersectingPoint_with_Line(double *arrPosWorking, double *arrVeloWorking
                             , double *pXintersect, double *pYintersect)
{
    // 1. Нахождение точки пересения прямой и плоскости в исходной системе координат
    double arrN[3] = {0.} ; // выектор нормали к плоскости
    fillNormalVect(arrN);
    ///


    double arrPointINtersect_PrSK[3] = {0.}; // это будет вектор точки пересечения висходной прям СК
    double arrDelta [3] ={0.};
    MtrxMinusMatrx(arrPosWorking, marrS0,3, 1, arrDelta);
    double val_t0 = ScalProduct(arrN ,  arrDelta, 3) ;
    double val_t = ScalProduct(arrN , arrVeloWorking, 3) ;
   /* if (fabs(val_t0) < 0.000000001)
    {  // точка arrPosWorking лежит на плоскости
        memcpy(arrPointINtersect_PrSK, arrPosWorking, 3 * sizeof(double));

    }
    else
    {
        if ( fabs(val_t )/ Norm3(arrVeloWorking) < 0.000000001)
        { // точка на плоскости не лежит, прямая параллельна плоскости
        return false;
        }
        double val_par =  -val_t0 / val_t ;
        double arrt0[3] = {0.};
        MatrxMultScalar(arrVeloWorking, 1, 3, val_par,arrt0);
        MtrxSumMatrx(arrPosWorking, arrt0, 3, 1, arrPointINtersect_PrSK) ;
    }
    */
    if ((fabs(val_t0) > 0.000000001 > DBL_MIN)&&( fabs(val_t )/ Norm3(arrVeloWorking) < 0.000000001))
    {  // точка на плоскости не лежит, прямая параллельна плоскости
        return false;
    }
    double val_par =  -val_t0 / val_t ;
    double arrt0[3] = {0.};
    MatrxMultScalar(arrVeloWorking, 1, 3, val_par,arrt0);
    MtrxSumMatrx(arrPosWorking, arrt0, 3, 1, arrPointINtersect_PrSK) ;

    double arrPointINtersect_SKP [3] ={0.};
    transform_xyzSSK_to_xyzSKP(arrPointINtersect_PrSK, arrPointINtersect_SKP);
    *pXintersect = arrPointINtersect_SKP[0] ;
    *pYintersect = arrPointINtersect_SKP[1 ];
    if (fabs( arrPointINtersect_SKP[2]) > 0.00000001)
    {
        //ShowMessage(L" ERROR QPlane::findIntersectingPoint_with_Line");
    }
    return true;
}
//-------------------------------------------------------------
//задана пломкость(*this) и прямая - точка arrPosWorking и  скорость arrVeloWorking
// точка и скорость заданы в исходной прямоугольной системе координат
// требуется найти точку пересечения прямой и плоскости в связанной системе координат.
//INPUT:
//arrPosWorking[3] , arrVeloWorking[3]
// OUTPUT:
// palf - угол между вектором скорости и плоскостью
// pt - время
// *pXintersect, *pYintersect координаты точки пересечения в СКП
bool QPlane::findIntersectingPoint_with_Line(double *arrPosWorking, double *arrVeloWorking
                             , double *pXintersect, double *pYintersect, double *palf
                                             , double *pt)
{
    // 1. Нахождение точки пересения прямой и плоскости в исходной сиситеме координат
    double arrN[3] = {0.} ; // выектор нормали к плоскости
    fillNormalVect(arrN);
    ///
    double arrPointINtersect_PrSK[3] = {0.}; // это будет вектор точки пересечения в исходной прям СК
    double arrDelta [3] ={0.};
    MtrxMinusMatrx(arrPosWorking, marrS0,3, 1, arrDelta);
    double val_t0 = ScalProduct(arrN ,  arrDelta, 3) ;
    double val_t = ScalProduct(arrN , arrVeloWorking, 3) ;
    if (fabs(val_t0) < 0.00000001)
    {  // точка arrPosWorking лежит на плоскости
        memcpy(arrPointINtersect_PrSK, arrPosWorking, 3 * sizeof(double));
        *pt =0;
        double arrPos_SKP[3] ={0.};
        transform_xyzSSK_to_xyzSKP(arrPosWorking, arrPos_SKP);
        *pXintersect = arrPos_SKP[0] ;
        *pYintersect = arrPos_SKP[1 ];
        return true;
    }
    else
    {
        if ( fabs(val_t )/ Norm3(arrVeloWorking) < 0.0000001)
        { // точка на плоскости не лежит, прямая параллельна плоскости
        return false;
        }
        *pt =  -val_t0 / val_t ;
        double arrt0[3] = {0.};
        MatrxMultScalar(arrVeloWorking, 1, 3, *pt,arrt0);
        MtrxSumMatrx(arrPosWorking, arrt0, 3, 1, arrPointINtersect_PrSK) ;
    }

    double arrPointINtersect_SKP [3] ={0.};
    transform_xyzSSK_to_xyzSKP(arrPointINtersect_PrSK, arrPointINtersect_SKP);
    *pXintersect = arrPointINtersect_SKP[0] ;
    *pYintersect = arrPointINtersect_SKP[1 ];
    if (fabs( arrPointINtersect_SKP[2]) > 0.00000001)
    {
        //ShowMessage(L" ERROR QPlane::findIntersectingPoint_with_Line");
    }
    *palf = M_PI/2. - acos(val_t/Norm3(arrVeloWorking));
    return true;
}
//-------------------------------------------------------------
//задана пломкость(*this) и прямая - точка arrPosWorking и  скорость arrVeloWorking
// точка и скорость заданы в исходной прямоугольной системе координат
// требуется найти точку пересечения прямой и плоскости в связанной системе координат.
//INPUT:
//arrPosWorking[3] , arrVeloWorking[3]
// OUTPUT:
// palf - угол между вектором скорости и плоскостью
// pt - время
// arrPointINtersect_PrSK[3] -  координаты точки пересечения в исходной прямоуг системе координат
bool QPlane::findIntersectingPoint_with_Line_PrSK(double *arrPosWorking, double *arrVeloWorking
                             , double *arrPointINtersect_PrSK, double *palf
                                             , double *pt)
{
    memset(arrPointINtersect_PrSK, 0.,3 * sizeof(double));
    // 1. Нахождение точки пересения прямой и плоскости в исходной сиситеме координат
    double arrN[3] = {0.} ; // выектор нормали к плоскости
    fillNormalVect(arrN);
    ///

    double arrDelta [3] ={0.};
    MtrxMinusMatrx(arrPosWorking, marrS0,3, 1, arrDelta);
    double val_t0 = ScalProduct(arrN ,  arrDelta, 3) ;
    double val_t = ScalProduct(arrN , arrVeloWorking, 3) ;
    if (fabs(val_t0) < 0.00000001)
    {  // точка arrPosWorking лежит на плоскости
        memcpy(arrPointINtersect_PrSK, arrPosWorking, 3 * sizeof(double));
        *pt =0;
        double arrPos_SKP[3] ={0.};
        transform_xyzSSK_to_xyzSKP(arrPosWorking, arrPos_SKP);
        return true;
    }
    else
    {
        if ( fabs(val_t )/ Norm3(arrVeloWorking) < 0.0000001)
        { // точка на плоскости не лежит, прямая параллельна плоскости
        return false;
        }
        *pt =  -val_t0 / val_t ;
        double arrt0[3] = {0.};
        MatrxMultScalar(arrVeloWorking, 1, 3, *pt,arrt0);
        MtrxSumMatrx(arrPosWorking, arrt0, 3, 1, arrPointINtersect_PrSK) ;
    }
    *palf = M_PI/2. - acos(val_t/Norm3(arrVeloWorking));
    return true;
}
//-----------------------------------------------
// формирование вектора нормали - 3-го столбца матрицы  marrF
//
void QPlane::fillNormalVect(double *arrN)
{
 arrN [0]	=  marrF [ 2];
 arrN [1]	=  marrF [ 5];
 arrN [2]	=  marrF [ 8];
}

//----------------------------------------------
// преобразование вектора из исходной прямоуг сиситемы координат в сиситему координат плоскости
//INPUT:
// arrPoint_PrSK[3]
// OUTPUT:
// arrPoint_SKP[3]
void QPlane::transform_xyzSSK_to_xyzSKP(double *arrPoint_PrSK, double *arrPoint_SKP)
{
	double arrDelta [3] ={0.};
	MtrxMinusMatrx(arrPoint_PrSK, marrS0,3, 1, arrDelta);
	MtrxTranspMultMatrx(marrF,3, 3, arrDelta, 1, arrPoint_SKP) ;
}
//----------------------------------------------
// преобразование вектора из сиситемы координат плоскости в исходную прямоуг сиситему координат
//INPUT:
// arrPoint_PrSK[3]
// OUTPUT:
// arrPoint_SKP[3]
void QPlane::transform_xyzSKP_to_xyzSSK(double *arrPoint_SKP, double *arrPoint_PrSK)
{
	double arrDelta [3] ={0.};
	MtrxMultMatrx(marrF,3, 3, arrPoint_SKP, 1, arrDelta) ;
	MtrxSumMatrx(arrDelta, marrS0,3, 1, arrPoint_PrSK);
}
//----------------------------------------------------
int QPlane::createInputDataReport(wchar_t*FileName, const bool bHeader)
{
    int len = wcslen(FileName) ;

    if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'x') // проверка, что
     && (FileName[len - 3] == 't') ) )  // указанный файл имеет расширение .flt
    {

      return 1 ;
    }

    FILE *fw ;

    if ((fw = _wfopen(FileName,L"a"))== NULL)

    {

     return 1 ;
    }
if (bHeader)
{
   fprintf(fw,"  Дата и время формирования отчета\n");
   time_t t = time(NULL);
   struct tm* aTm = localtime(&t);
   fprintf(fw,"  Год = %04d\n",aTm->tm_year+1900);
   fprintf(fw,"  Mесяц = %02d\n",aTm->tm_mon+1);
   fprintf(fw,"  День = %02d\n",aTm->tm_mday);
   fprintf(fw,"  Время = %02d:%02d:%02d\n",aTm->tm_hour, aTm->tm_min, aTm->tm_sec);
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
}
fprintf(fw,"***************************************\n");
fprintf(fw,"****** ПЛОСКОСТЬ *********************************\n");
fprintf(fw,"  Плоскость задается вектором положения начала собственной\n");
fprintf(fw,"  системы координат S0 и матрицей перехода из собственной \n");
fprintf(fw,"  системы координат в исходную f \n");

fprintf(fw,"  S0 = { %8.2f,   %8.2f,  %8.2f}\n",marrS0[0],marrS0[1],marrS0[2]);
fprintf(fw,"  f =  { %6.4f    %6.4f,  %6.4f \n",marrF[0],marrF[1],marrF[2]);
fprintf(fw,"         %6.4f    %6.4f,  %6.4f \n",marrF[3],marrF[4],marrF[5]);
fprintf(fw,"         %6.4f    %6.4f,  %6.4f }\n ",marrF[6],marrF[7],marrF[8]);
fclose(fw);
}

