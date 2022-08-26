#include "LongPlane.h"

#include <stdio.h>
#include <stdlib.h>
#include  <string.h>
#include  <math.h>
#include "Plane.h"
#include "MatrixProccess.h"
//#include "UrPointXY.h"

//--------------------------------------------------------------------------------------
 TLongPlane ::TLongPlane()
{
    long double arrT [9] = {0.};
    memcpy(marrS0, arrT, 3 * sizeof(long double));
    memcpy(marrF, arrT, 9 * sizeof(long double));

}

 // оператор присваивания
 TLongPlane &TLongPlane::operator=(const TLongPlane & R)
 {
    memcpy(marrS0, R.marrS0, 3 * sizeof(long double));
    memcpy(marrF, R.marrF, 9 * sizeof(long double));
    return *this ;
 }

 // конструктор копирования
 TLongPlane::TLongPlane (const TLongPlane &R)
 {
    memcpy(marrS0, R.marrS0, 3 * sizeof(long double));
    memcpy(marrF, R.marrF, 9 * sizeof(long double));
 }


 // парам констр
TLongPlane ::TLongPlane(long  double* arrS0,long  double* arrF)
{
    memcpy(marrS0, arrS0, 3 * sizeof(long double));
    memcpy(marrF, arrF, 9 * sizeof(long double));

}

bool TLongPlane::findIntersectingPoint_with_Line(long double *arrPosWorking, long double *arrVeloWorking
                                                 , long double *pvalPointIntersectX, long double *pvalPointIntersectY)
{
    // 1. Нахождение точки пересения прямой и плоскости в исходной сиситеме координат
    long double arrN[3] = {0.} ; // выектор нормали к плоскости
    fillNormalVect(arrN);
    ///
    long double arrPointINtersect_PrSK[3] = {0.}; // это будет вектор точки пересечения висходной прям СК
    long double arrDelta [3] ={0.};
    MtrxMinusMatrx(arrPosWorking, marrS0,3, 1, arrDelta);
    double val_t0 = ScalProduct(arrN ,  arrDelta, 3) ;
    double val_t = ScalProduct(arrN , arrVeloWorking, 3) ;
    if (fabsl(val_t0) < 0.00000001)
    {  // точка arrPosWorking лежит на плоскости или вектор скорости параллелен плоскости и пересечения нет
        memcpy(arrPointINtersect_PrSK, arrPosWorking, 3 * sizeof(long double));
        return false;
    }
    else
    {
        if ( fabsl(val_t )/ Norm3(arrVeloWorking) < 0.0000001)
        { // точка на плоскости не лежит, прямая параллельна плоскости
        return false;
        }
        long double val_par =  -val_t0 / val_t ;
        long double arrt0[3] = {0.};
        MatrxMultScalar(arrVeloWorking, 1, 3, val_par,arrt0);
        MtrxSumMatrx(arrPosWorking, arrt0, 3, 1, arrPointINtersect_PrSK) ;
    }

    long double arrPointINtersect_SKP [3] ={0.};
    transform_xyzSSK_to_xyzSKP(arrPointINtersect_PrSK, arrPointINtersect_SKP);
    *pvalPointIntersectX = arrPointINtersect_SKP[0] ;
    *pvalPointIntersectY = arrPointINtersect_SKP[1 ];
    if (fabs( arrPointINtersect_SKP[2]) > 0.00000001)
    {
        //ShowMessage(L" ERROR TLongPlane::findIntersectingPoint_with_Line");
    }
    return true;
}
//-----------------------------------------------
// формирование вектора нормали - 3-го столбца матрицы  marrF
//
void TLongPlane::fillNormalVect(long double *arrN)
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
void TLongPlane::transform_xyzSSK_to_xyzSKP(long double *arrPoint_PrSK,long  double *arrPoint_SKP)
{
    long double arrDelta [3] ={0.};
    MtrxMinusMatrx(arrPoint_PrSK, marrS0,3, 1, arrDelta);
    MtrxTranspMultMatrx(marrF,3, 3, arrDelta, 1, arrPoint_SKP) ;
}
//----------------------------------------------
// преобразование вектора из сиситемы координат плоскости в исходную прямоуг сиситему координат
//INPUT:
// arrPoint_PrSK[3]
// OUTPUT:
// arrPoint_SKP[3]
void TLongPlane::transform_xyzSKP_to_xyzSSK(long double *arrPoint_SKP,long  double *arrPoint_PrSK)
{
    long double arrDelta [3] ={0.};
    MtrxMultMatrx(marrF,3, 3, arrPoint_SKP, 1, arrDelta) ;
    MtrxSumMatrx(arrDelta, marrS0,3, 1, arrPoint_PrSK);
}

