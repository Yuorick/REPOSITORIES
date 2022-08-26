#include "ContObj6.h"

#include <string.h>
#include <math.h>
#include <float.h>

#include "MatrixProccess.h"
#include "CalcCorMatrx.h"
#include "MinSquare.h"
#include "Gauss.h"

#define NUmAng  3
extern const double COeff = 5000.;

QContObj6::QContObj6()
{
    memset(marrEilerCntrKP, 0, 3 * sizeof(double));
    memset(marrDispEilerCntrKP, 0, 3 * sizeof(double));
    memset(marrXTrue_RLK_PSK, 0, QUANT_GDG_PRMS * sizeof(double));
    memset(marrXTrue_IU_PSK, 0, QUANT_GDG_PRMS * sizeof(double));
    memset(marrXZv_RLK_PSK, 0, QUANT_GDG_PRMS * sizeof(double));
    memset(marrXZv_IU_PSK, 0, QUANT_GDG_PRMS * sizeof(double));
    memset(marrS_KGSK, 0, 3 * sizeof(double));
    memset(marrCorMtrxRLK, 0, 9 * sizeof(double));
    memset(marrEilerCntrKP_Zv, 0, 3 * sizeof(double));
    memset(marrSZv_KGSK, 0, 3 * sizeof(double));

    mAngleSigIU = 0.0001;
    mDistSigIU =0.1;
    mVessSKZ = 0.;

    mCntrlObjSKZ= 0.;



}
// конструктор копирования
 QContObj6 ::QContObj6 (const QContObj6 &R)
 {
     memcpy(marrEilerCntrKP, R.marrEilerCntrKP, 3 * sizeof(double));

     memcpy(marrDispEilerCntrKP, R.marrDispEilerCntrKP, 3 * sizeof(double));

     memcpy(marrXTrue_RLK_PSK, R.marrXTrue_RLK_PSK, QUANT_GDG_PRMS * sizeof(double));

     memcpy(marrXTrue_IU_PSK, R.marrXTrue_IU_PSK, QUANT_GDG_PRMS * sizeof(double));

     memcpy(marrXZv_RLK_PSK, R.marrXZv_RLK_PSK, QUANT_GDG_PRMS * sizeof(double));

     memcpy(marrXZv_IU_PSK, R.marrXZv_IU_PSK, QUANT_GDG_PRMS * sizeof(double));

     memcpy(marrS_KGSK, R.marrS_KGSK, 3 * sizeof(double));

     memcpy(marrSZv_KGSK, R.marrSZv_KGSK, 3 * sizeof(double));

     memcpy(marrCorMtrxRLK, R.marrCorMtrxRLK, 9 * sizeof(double));

     memcpy(marrEilerCntrKP_Zv, R.marrEilerCntrKP_Zv, 3 * sizeof(double));


     mAngleSigIU = R.mAngleSigIU;
     mDistSigIU = R.mDistSigIU;
     mVessSKZ = R.mVessSKZ;
     mCntrlObjSKZ = R.mCntrlObjSKZ;
 }
 // оператор присваивания
 QContObj6 &QContObj6::operator=(const QContObj6  &R)
 {
     memcpy(marrEilerCntrKP, R.marrEilerCntrKP, 3 * sizeof(double));

     memcpy(marrDispEilerCntrKP, R.marrDispEilerCntrKP, 3 * sizeof(double));

     memcpy(marrXTrue_RLK_PSK, R.marrXTrue_RLK_PSK, QUANT_GDG_PRMS * sizeof(double));

     memcpy(marrXTrue_IU_PSK, R.marrXTrue_IU_PSK, QUANT_GDG_PRMS * sizeof(double));

     memcpy(marrXZv_RLK_PSK, R.marrXZv_RLK_PSK, QUANT_GDG_PRMS * sizeof(double));

     memcpy(marrXZv_IU_PSK, R.marrXZv_IU_PSK, QUANT_GDG_PRMS * sizeof(double));

     memcpy(marrS_KGSK, R.marrS_KGSK, 3 * sizeof(double));

     memcpy(marrSZv_KGSK, R.marrSZv_KGSK, 3 * sizeof(double));

     memcpy(marrCorMtrxRLK, R.marrCorMtrxRLK, 9 * sizeof(double));

    memcpy(marrEilerCntrKP_Zv, R.marrEilerCntrKP_Zv, 3 * sizeof(double));


     mAngleSigIU = R.mAngleSigIU;
     mDistSigIU = R.mDistSigIU;
     mVessSKZ = R.mVessSKZ;
     mCntrlObjSKZ = R.mCntrlObjSKZ;


   return *this ;
 }

 // парм консьтруктор 1
 QContObj6::QContObj6(const double *arrEilerCntrKP, const double *arrDispEilerCntrKP
            ,const double * arrXTrue_RLK_PSK,const double * arrXTrue_IU_PSK
            ,const double *arrXZv_RLK_PSK,const double *arrXZv_IU_PSK
            ,const double * arrS_KGSK,const double * arrCorMtrxRLK
            , const double AngleSigIU, const double DistSigIU
            ,const double VessSKZ,const double  CntrlObjSKZ )
 {
     memcpy(marrEilerCntrKP, arrEilerCntrKP, 3 * sizeof(double));
     memcpy(marrDispEilerCntrKP, arrDispEilerCntrKP, 3 * sizeof(double));
     memcpy(marrXTrue_RLK_PSK, arrXTrue_RLK_PSK, QUANT_GDG_PRMS * sizeof(double));
     memcpy(marrXTrue_IU_PSK, arrXTrue_IU_PSK, QUANT_GDG_PRMS * sizeof(double));
     memcpy(marrXZv_RLK_PSK, arrXZv_RLK_PSK, QUANT_GDG_PRMS * sizeof(double));
     memcpy(marrXZv_IU_PSK, arrXZv_IU_PSK, QUANT_GDG_PRMS * sizeof(double));
     memcpy(marrS_KGSK, arrS_KGSK, 3 * sizeof(double));
     memcpy(marrSZv_KGSK, arrS_KGSK, 3 * sizeof(double));
     memcpy(marrCorMtrxRLK, arrCorMtrxRLK, 9 * sizeof(double));
     memcpy(marrEilerCntrKP_Zv, marrEilerCntrKP, 3 * sizeof(double));


     mAngleSigIU = AngleSigIU;
     mDistSigIU = DistSigIU;
     mVessSKZ =VessSKZ;
     mCntrlObjSKZ = CntrlObjSKZ;

 }


//----------------------------------------------------
// имитация измерений РЛК и ИУ
void QContObj6::imitateMeasures(const bool BNoise,double *arrVZvRLK
                                  , double *arrVZvIU)
{

    // 1. истинное положение в АСфСК

    memcpy(marrSZv_KGSK, marrS_KGSK, 3 * sizeof(double));
    if ( BNoise)
    {
        double arrVessPos[3] = {0.};
        arrVessPos[0] = getGauss(0, mVessSKZ );
        arrVessPos[1] = getGauss(0, mVessSKZ );
    for (int i =0; i < 3; ++i)
    {

        marrSZv_KGSK[i] +=  getGauss(0, mCntrlObjSKZ) - arrVessPos[i];

    }
    }
    double arrVTrue_RLK[3] = {0.};
    recalcPositionFromKGSK_to_SphericalSK(marrEilerCntrKP,marrS_KGSK,marrXTrue_RLK_PSK,arrVTrue_RLK);

    // 2. имитация замера в АСфСК
    for (int i =0; i< 3; ++i)
    {
        if ( BNoise)
        {
         arrVZvRLK[i] =  arrVTrue_RLK[i] + getGauss(0, sqrt(marrCorMtrxRLK[3 *i +i]) );
        }
        else
        {
         arrVZvRLK[i] =  arrVTrue_RLK[i];
        }
    }
    ///

    // 1. истинное положение в ИУСфСК
    double arrVTrue_IU[3] = {0.};
    recalcPositionFromKGSK_to_SphericalSK(marrEilerCntrKP,marrS_KGSK, marrXTrue_IU_PSK,arrVTrue_IU);


    // 2. имитация замера в ИУСфСК
    arrVZvIU[0] =  arrVTrue_IU[0] ;
    arrVZvIU[1] =  arrVTrue_IU[2] ;
    if ( BNoise)
    {
        arrVZvIU[0] +=  getGauss(0, mAngleSigIU);
        arrVZvIU[1] +=  getGauss(0, mAngleSigIU);
    }
    // 3. имитация ошибок СИНС
    memcpy(marrEilerCntrKP_Zv,marrEilerCntrKP, 3 * sizeof(double));
    if(BNoise)
    {
    for (int i =0; i < 3; ++i)
    {
      marrEilerCntrKP_Zv[i] =  marrEilerCntrKP[i] + getGauss(0, sqrt(marrDispEilerCntrKP[i]));
    }
    }


 }



/*
//-----------------------------------------------------
//вычисление вектора положениея КО в сферической сиситеме координат гаджета
//INPUT:
//arrXGdg[5] - вектор позиционирования гаджета в ПСК
// OUTPUT:
//arrVGdg[3]- вектор сфкеер координат КО, Betta, R, Eps
void QContObj6::recalcPositionFromKGSK_to_SphericalSK(double *arrEilerCntrKP
         , double *arrS_KGSK,double *arrXGdg,double *arrVGdg)
{
// 1. истинное положение в ПСК
double matrPereh_PSK_V_KGSK[9] = {0.};
calcMatr_PSK_v_KGSK(arrEilerCntrKP,matrPereh_PSK_V_KGSK) ;
double arrSTrue_PSK[3] = {0.};
MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrS_KGSK,1, arrSTrue_PSK) ;
///
// 2. истинное положение в ПСК-УСТРОЙСТВО
double arrSTrue_PSK_Gdg[3] = {0.};
MtrxMinusMatrx(arrSTrue_PSK, arrXGdg,1, 3, arrSTrue_PSK_Gdg);
///

// 3. истинное положение  в АСПК
double matrPereh_ASK_V_PSK[9]= {0.};
calcMtrx3_ASPK_v_PSK(&(arrXGdg[3]),matrPereh_ASK_V_PSK) ;
double arrSTrue_ASPK[3] = {0.};
MtrxTranspMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrSTrue_PSK_Gdg,1, arrSTrue_ASPK) ;
///

// 4. положение в АСфСАК

recalcCoord_INTO_Spherical(arrSTrue_ASPK, arrVGdg[1], arrVGdg[0]
        , arrVGdg[2]);
///
}

*/
//-----------------------------------------------------
//вычисление вектора положениея КО в сферической сиситеме координат гаджета
//INPUT:
//arrXGdg[6] - вектор позиционирования гаджета в ПСК
// OUTPUT:
//arrVGdg[3]- вектор сфкеер координат КО, Betta, R, Eps
void QContObj6::recalcPositionFromKGSK_to_SphericalSK(double *arrEilerCntrKP0
         , double *arrS_KGSK0,double *arrXGdg0,double *arrVGdg0)
{
    long double arrEilerCntrKP[3] = {0.};
    long double arrS_KGSK[3] = {0.};
    long double arrXGdg[QUANT_GDG_PRMS] = {0.};
    long double arrVGdg[3] = {0.};
    setLongDblArray(arrEilerCntrKP0, 3,arrEilerCntrKP);
    setLongDblArray( arrS_KGSK0, 3,arrS_KGSK);
    setLongDblArray( arrXGdg0,QUANT_GDG_PRMS,arrXGdg);

// 1. истинное положение в ПСК
long double matrPereh_PSK_V_KGSK[9] = {0.};
calcMatr_PSK_v_KGSK(arrEilerCntrKP,matrPereh_PSK_V_KGSK) ;
long double arrSTrue_PSK[3] = {0.};
MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrS_KGSK,1, arrSTrue_PSK) ;
///
// 2. истинное положение в ПСК-УСТРОЙСТВО
long double arrSTrue_PSK_Gdg[3] = {0.};
MtrxMinusMatrx(arrSTrue_PSK, arrXGdg,1, 3, arrSTrue_PSK_Gdg);
///

// 3. истинное положение  в АСПК
long double matrPereh_ASK_V_PSK[9]= {0.};
calcMtrx3_ASPK_v_PSK(&(arrXGdg[3]),matrPereh_ASK_V_PSK) ;
long double arrSTrue_ASPK[3] = {0.};
MtrxTranspMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrSTrue_PSK_Gdg,1, arrSTrue_ASPK) ;
///

// 4. положение в АСфСАК

recalcCoord_INTO_Spherical(arrSTrue_ASPK, arrVGdg[1], arrVGdg[0]
        , arrVGdg[2]);

setDblArray(arrVGdg, 3, arrVGdg0);
///
}


//-----------------------------------------------------
//
// INPUT:
//BDalnomer -признак дальномера
//arrX_RLK- вектор параметров позиционирования РЛК в ПСК
//arrX_IU- вектор параметров позиционирования ИУ в ПСК
//arrVZv_RLK - вектор сферических координат точки в АСфСК
//OUTPUT:
//arrf[BDalnomer +2]- вектор сферических координат точки в ИУСфСК
void QContObj6::calc_f(const bool BDalnomer,double *arrX_RLK,double *arrX_IU
                         , double *arrVZv_RLK, double *arrf)
{

    ///
    // 1.  положение в АСПК
    double arrS_ASPK[3] = {0.};
    recalcSphericalCoord_INTO_Rectangular(arrVZv_RLK[1],arrVZv_RLK[0]
                                               ,arrVZv_RLK[2], arrS_ASPK);
    ///

    // 2.  положение  в ПСК-РЛК
    double matrPereh_ASK_V_PSK[9]= {0.};
    calcMatr_ASK_v_PSK(&(arrX_RLK[3]),matrPereh_ASK_V_PSK) ;
    double arrS_PSK_RLK[3] = {0.};
    MtrxMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrS_ASPK,1, arrS_PSK_RLK) ;
    ///

    // 3. положение в ПСК-ИУ
    double arrS_PSK_IU[3] = {0.},  arrPWave[3] ={0.};
    MtrxMinusMatrx(arrX_RLK, arrX_IU,1, 3, arrPWave);
    MtrxSumMatrx(arrS_PSK_RLK, arrPWave,1, 3, arrS_PSK_IU) ;
    ///

    // 4. положение в ИУСПК
    double matrPereh_ASK_V_PSK1[9]= {0.};
    calcMatr_ASK_v_PSK(&(arrX_IU[3]),matrPereh_ASK_V_PSK1) ;
    double arrS_IUSPK[3] = {0.};
    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK1,3, 3, arrS_PSK_IU,1, arrS_IUSPK) ;
    ///

    // 5.  положение в СфСК-ИУ
    double arrV_IU[3] = {0.};
    recalcCoord_INTO_Spherical(arrS_IUSPK, arrV_IU[1], arrV_IU[0]
            , arrV_IU[2]);
    ///

    // 6
    if (BDalnomer)
    {
        arrf[0] =  arrV_IU[0] ;
        arrf[2] =  arrV_IU[2] ;
        arrf[1] =  arrV_IU[1] ;


    }
    else
    {
        arrf[0] =  arrV_IU[0] ;
        arrf[1] =  arrV_IU[2] ;

    }


}
/*
//-------------------------------------------
// INPUT:
//BDalnomer -признак дальнгомера
//arrX[5]- вектор параметров позиционирования РЛК в ПСК-ИУ
//arrVZv_RLK - вектор сферических координат точки в АСфСК
//OUTPUT:
//arrViu[BDalnomer +2]- вектор сферических координат точки в ИУСфСК
void QContObj6::calc_Viu(double *arrXRLKCur,double *arrX_IU
                         , double *arrVZv_RLK, double *arrViu)
{

    ///
    // 1.  положение в АСПК
    double arrS_ASPK[3] = {0.};
    recalcSphericalCoord_INTO_Rectangular(arrVZv_RLK[1],arrVZv_RLK[0]
                                               ,arrVZv_RLK[2], arrS_ASPK);
    ///

    // 2.  положение  в ПСК-РЛК
    double matrPereh_ASK_V_PSK[9]= {0.};

    calcMtrx3_ASPK_v_PSK(&(arrXRLKCur[3]),matrPereh_ASK_V_PSK) ;
    double arrS_PSK_RLK[3] = {0.};
    MtrxMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrS_ASPK,1, arrS_PSK_RLK) ;
    ///

    // 3. положение в ПСК-ИУ
    double arrS_PSK_IU[3] = {0.},  arrPWave[3] ={0.};
    MtrxSumMatrx(arrS_PSK_RLK, arrXRLKCur,1, 3, arrS_PSK_IU) ;
    ///

    // 4. положение в ИУСПК
    double matrPereh_ASK_V_PSK1[9]= {0.};

    calcMtrx3_ASPK_v_PSK(&(arrX_IU[3]),matrPereh_ASK_V_PSK1) ;
    double arrS_IUSPK[3] = {0.};
    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK1,3, 3, arrS_PSK_IU,1, arrS_IUSPK) ;
    ///

    // 5.  положение в СфСК-ИУ
    double arrV_IU[3] = {0.};
    recalcCoord_INTO_Spherical(arrS_IUSPK, arrV_IU[1], arrV_IU[0]
            , arrV_IU[2]);
    ///

    // 6
    arrViu[0] =  arrV_IU[0] ;
    arrViu[1] =  arrV_IU[2] ;


}
*/
//-------------------------------------------
// INPUT:
//BDalnomer -признак дальнгомера
//arrX[5]- вектор параметров позиционирования РЛК в ПСК-ИУ
//arrVZv_RLK - вектор сферических координат точки в АСфСК
//OUTPUT:
//arrViu[BDalnomer +2]- вектор сферических координат точки в ИУСфСК
void QContObj6::calc_Viu(double *arrXRLKCur0,double *arrX_IU0
                         , double *arrVZv_RLK0, double *arrViu0)
{

    ///
    /// \brief setLongDblArray
    /// \param LEN

    int iiu = QUANT_GDG_PRMS;
    long double arrXRLKCur[QUANT_GDG_PRMS] = {0.};
    long double arrX_IU[QUANT_GDG_PRMS] = {0.};
    long double arrVZv_RLK[3] = {0.};
    long double arrViu[2] = {0.};

    setLongDblArray( arrXRLKCur0, QUANT_GDG_PRMS,arrXRLKCur);

    setLongDblArray( arrX_IU0, QUANT_GDG_PRMS,arrX_IU);
    setLongDblArray( arrVZv_RLK0, 3,arrVZv_RLK);

    // 1.  положение в АСПК
    long double arrS_ASPK[3] = {0.};
    recalcSphericalCoord_INTO_Rectangular(arrVZv_RLK[1],arrVZv_RLK[0]
                                               ,arrVZv_RLK[2], arrS_ASPK);
    ///

    // 2.  положение  в ПСК-РЛК
    long double matrPereh_ASK_V_PSK[9]= {0.};

    calcMtrx3_ASPK_v_PSK(&(arrXRLKCur[3]),matrPereh_ASK_V_PSK) ;
    long double arrS_PSK_RLK[3] = {0.};
    MtrxMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrS_ASPK,1, arrS_PSK_RLK) ;
    ///

    // 3. положение в ПСК-ИУ
    long double arrS_PSK_IU[3] = {0.},  arrPWave[3] ={0.};
    MtrxSumMatrx(arrS_PSK_RLK, arrXRLKCur,1, 3, arrS_PSK_IU) ;
    ///

    // 4. положение в ИУСПК
    long double matrPereh_ASK_V_PSK1[9]= {0.};

    calcMtrx3_ASPK_v_PSK(&(arrX_IU[3]),matrPereh_ASK_V_PSK1) ;
    long double arrS_IUSPK[3] = {0.};
    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK1,3, 3, arrS_PSK_IU,1, arrS_IUSPK) ;
    ///

    // 5.  положение в СфСК-ИУ
    long double arrV_IU[3] = {0.};
    recalcCoord_INTO_Spherical(arrS_IUSPK, arrV_IU[1], arrV_IU[0]
            , arrV_IU[2]);
    ///

    // 6
    arrViu[0] =  arrV_IU[0] ;
    arrViu[1] =  arrV_IU[2] ;

  setDblArray(arrViu, 2, arrViu0);
}
//-------------------------------------------

// вычисление разностной производной функции f по переменной
// с номером j
//функция f зависит от 7 переменных
//от переменных вектора позиционирования РЛК в ПСК - arrX_RLK
//углов Betta и Eps
//OUTPUT:
//arr_df_po_dXj[2 + BDalnomer] - вектор частных производных ыектор функции f
//по переменной с номером NUmj (нумерация с 0)
void QContObj6::calc_df_po_dXj(const bool BDalnomer,double *arrX_RLK,double *arrX_IU
                         , double *arrV_RLK,const int NUmj, double *arr_df_po_dXj)
{
    double arrX_RLK1[5] = {0.};
    double arrX_IU1[5] = {0.};
    memcpy(arrX_RLK1, arrX_RLK, 5 * sizeof(double));
    memcpy(arrX_IU1, arrX_IU, 5 * sizeof(double));
    double delta = 0.;
    if (NUmj < 3)
    {
     delta = 0.05;
     arrX_RLK1[NUmj]+= delta;
    }
    else
    {
       if (NUmj < 5)
       {
         delta = 0.0001;
         arrX_RLK1[NUmj]+= delta;
       }
       else
       {
         delta = 0.0001;
         arrX_IU1[NUmj -2]+= delta;
       }
    }

    double arrf[3] = {0.},arrf1[3] = {0.}, arrT[3] = {0.},arrDelta[3] = {0.};

    calc_f( BDalnomer,arrX_RLK,arrX_IU, arrV_RLK, arrf);
    calc_f( BDalnomer,arrX_RLK1,arrX_IU1, arrV_RLK, arrf1);
    MtrxMinusMatrx(arrf1, arrf,1, 3, arrT);
    MatrxMultScalar(arrT, 1, 3, 1./delta,arrDelta);

    for (int i =0; i < (2 +BDalnomer);++i)
    {
      arr_df_po_dXj[i] =  arrDelta[i] ;
    }


}
//------------------------------


//-----------------------------------------------------------------------
// вычисление оценки вектора позиционирования РЛК и ИУ в ПСК-ИУ по Варианту 1
// INPUT:
//pContObjectArr - массив КО
//LEnArr - кол-во КО
// BNoise - признак имитации случайных шумов в измерениях
//OUTPUT:
//arrXEst[5] - оценка вектора позиционирования РЛК в ИУ
//arrMtrxCor[25] - коррел матрица ошибок оценивания
double QContObj6::imitate_and_estimateX_Var1(QContObj6 *pContObjectArr,const int QuantObj
       , const bool BNoise, double *arrXEst, double *arrMtrxCor)
{
    // векторы измерений РЛК
 double *parrVZv_RLK = new double [QuantObj * 3 ];

 // векторы измерений ИУ
 double *parrVZv_IU = new double [QuantObj * 2 ];

 // векторы истинного положения цели в ИУ
 //double *parrVTrue_IU = new double [QuantObj * (BDalnomer +2) ];

 // имитация измерений РЛК и ИУ
  for(int i =0; i < QuantObj; ++i)
  {
    pContObjectArr[i].imitateMeasures(BNoise,&(parrVZv_RLK[i * 3])
            , &(parrVZv_IU[i * 2])) ;
  }

  double rez = estimateParams_Var_1(pContObjectArr,QuantObj
                             ,parrVZv_RLK,parrVZv_IU, arrXEst, arrMtrxCor);

 delete []parrVZv_RLK;
 delete []parrVZv_IU;

  return rez;
}

//----------------------------------------------------
 double QContObj6::estimateParams_Var_1(QContObj6 *pContObjectArr,const int QuantObj
             ,double *parrVZv_RLK,double *parrVZv_IU, double *arrXEst, double *arrMtrxCor)  //++
 {
      //0.
     double arrMtrx_dFGr_po_dX_Inv[QUANT_GDG_PRMS * QUANT_GDG_PRMS]= {0.};

     // 1. формирование вектора начальгного приближения
     int iiii = QUANT_GDG_PRMS;
     double arrXEstCur[QUANT_GDG_PRMS]= {0.};
     //combine(pContObjectArr[0].marrXZv_RLK_PSK,pContObjectArr[0].marrXZv_IU_PSK,arrXEstCur);
     memcpy(arrXEstCur, pContObjectArr[0].marrXZv_RLK_PSK,QUANT_GDG_PRMS * sizeof(double) );
     MtrxMinusMatrx(pContObjectArr[0].marrXZv_RLK_PSK, pContObjectArr[0].marrXZv_IU_PSK,3, 1, arrXEstCur);
     ///

     // 2. вычисление начальной невязки
     double valNeviaz0 = calcSumNeviazka_Var1(arrXEstCur,pContObjectArr
            , QuantObj, parrVZv_RLK,parrVZv_IU);
     ///

     // 3. итреац процесс
     double arrXCur[QUANT_GDG_PRMS]= {0.};
     double valNeviaz1=0.;
     double arrDel[QUANT_GDG_PRMS] = {0.};
     for (int i =0; i < 500; ++i)
     {

       doOneIteration_Var1_New(arrXEstCur,pContObjectArr, QuantObj
                               ,parrVZv_RLK,parrVZv_IU, arrDel, arrMtrx_dFGr_po_dX_Inv);
       double del = NormVect(arrDel, QUANT_GDG_PRMS);
         double coef = 1.;//0.25;
         MatrxMultScalar( arrDel, QUANT_GDG_PRMS, 1, coef, arrDel);
         MtrxMinusMatrx(arrXEstCur, arrDel,QUANT_GDG_PRMS, 1, arrXCur);
        valNeviaz1 = calcSumNeviazka_Var1(arrXCur,pContObjectArr, QuantObj
                    ,parrVZv_RLK,parrVZv_IU);

       if (valNeviaz1 > valNeviaz0)
       {
         int iou=0;
          break;
       }
       if ((del < 1.E-8)|| (valNeviaz1 < 1.E-10 ))
       {
           memcpy(arrXEstCur, arrXCur, QUANT_GDG_PRMS * sizeof(double));
           valNeviaz0 = valNeviaz1;
           break;
       }
       memcpy(arrXEstCur, arrXCur, QUANT_GDG_PRMS * sizeof(double));
       valNeviaz0 = valNeviaz1;
     }

     memcpy(arrXEst, arrXEstCur,QUANT_GDG_PRMS * sizeof(double));


     calcCorMtrx_Var1(pContObjectArr, QuantObj
                  ,parrVZv_RLK,parrVZv_IU, arrXEst,arrMtrx_dFGr_po_dX_Inv, arrMtrxCor);


     return valNeviaz0;
 }

 //-------------------------------------------------------------
bool QContObj6::calcCorMtrx_Var1(QContObj6 *pContObjectArr,const int QuantObj
             ,double *parrVZv_RLK,double *parrVZv_IU, double *arrXEst
             , double *arrMtrx_dFGr_po_dX_Inv, double *arrMtrxCor)
 {
    double arrSum1[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrSum2[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrQKQT_Cur[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrgTgKggT_Cur[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    for(int i =0; i < QuantObj; ++i)
    {
      pContObjectArr[i].calcMtrx_QKQT(&(parrVZv_RLK[3 * i]),arrXEst,arrQKQT_Cur);

      MtrxSumMatrx(arrSum1, arrQKQT_Cur,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrSum1) ;

      pContObjectArr[i].calcMtrx_gTKg(&(parrVZv_RLK[3 * i])
              ,&(parrVZv_IU[2 * i]),arrXEst,arrgTgKggT_Cur);
      MtrxSumMatrx(arrSum2, arrgTgKggT_Cur,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrSum2) ;

    }
    double arrK1[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrK2[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrMtrx_dFGr_po_dX_Inv1[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};

    memcpy(arrMtrx_dFGr_po_dX_Inv1, arrMtrx_dFGr_po_dX_Inv, QUANT_GDG_PRMS* QUANT_GDG_PRMS * sizeof(double));
    MtrxMultMatrx_MultMatrxTransp(arrMtrx_dFGr_po_dX_Inv,arrSum1
                                  ,arrMtrx_dFGr_po_dX_Inv1,QUANT_GDG_PRMS, arrK1);

    MtrxMultMatrx_MultMatrxTransp(arrMtrx_dFGr_po_dX_Inv,arrSum2
                                  ,arrMtrx_dFGr_po_dX_Inv1,QUANT_GDG_PRMS, arrK2);

    double arrMtrxCor0[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.};

    MtrxSumMatrx(arrK1, arrK2,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrMtrxCor0) ;

    // проверка
   // memcpy(arrMtrxCor0, arrK2, QUANT_GDG_PRMS* QUANT_GDG_PRMS * sizeof(double));
    ///

    double arrCoeff[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.}, arrCoeff1[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.};
    fillE(arrCoeff, QUANT_GDG_PRMS);
    fillE(arrCoeff1, QUANT_GDG_PRMS);
    for (int i =0; i < 3; ++i)
    {
        arrCoeff[ i *QUANT_GDG_PRMS + i] = COeff;
        arrCoeff1[i *QUANT_GDG_PRMS + i] = COeff;
    }
    MtrxMultMatrx_MultMatrx(arrCoeff,arrMtrxCor0,arrCoeff1,QUANT_GDG_PRMS, arrMtrxCor)  ;


}
//-----------------------------------------------------------------
void QContObj6::calcMtrx_QKQT(double *arrVZv_RLK, double *arrX,double *arrQKQT)
{
    // 1
    double arrg[2 * QUANT_GDG_PRMS] = {0.};
    calc_dViu_po_dX_Var1_(arrX,marrXZv_IU_PSK, arrVZv_RLK,arrg);
    /*for (int i =0; i < 2; ++i)
    {
        for (int j=0; j < 3; ++j)
        {
          arrg[i *QUANT_GDG_PRMS + j]  =  arrg[i *QUANT_GDG_PRMS + j]/COeff ;
        }
    }*/

    ///

    // 2 вырезание первых 3 столбцов
    double arrT0[6] = {0.};
    for (int i =0; i < 2; ++i)
    {
        for (int j=0; j < 3; ++j)
        {
          arrT0[i * 3 +j] =  arrg[i *QUANT_GDG_PRMS + j]/COeff ;
        }
    }

    // 3 ваычисление матрицы M
    double matrPereh_ASK_V_PSK_RLK[9] = {0.};
    calcMtrx3_ASPK_v_PSK(&(arrX[3]),matrPereh_ASK_V_PSK_RLK) ;
    ///

    // 4 вычисление dA_po_dV
    double arr_dA_po_dV[9]={0.};
    calc_dRectangularl_po_dSpherical(arrVZv_RLK, arr_dA_po_dV);
    ///

    double arrT1[6] = {0.},arrT2[6] = {0.}, arrQ[3 *QUANT_GDG_PRMS] = {0.};
    MtrxMultMatrx(arrT0,2, 3, matrPereh_ASK_V_PSK_RLK,3, arrT1) ;
    MtrxMultMatrx(arrT1,2, 3, arr_dA_po_dV,3, arrT2) ;
    MtrxTranspMultMatrx(arrg,2, QUANT_GDG_PRMS, arrT2,3, arrQ) ;


    double arrT3[3 *QUANT_GDG_PRMS] = {0.};
    MtrxMultMatrx(arrQ,QUANT_GDG_PRMS, 3, marrCorMtrxRLK,3, arrT3) ;
    MtrxMultMatrxTransp(arrT3,QUANT_GDG_PRMS, 3, arrQ,QUANT_GDG_PRMS, arrQKQT) ;

}
//-----------------------------------------------------------------
void QContObj6::calcMtrx_gTKg(double *arrVZv_RLK
              ,double *arrVZv_IU, double *arrX, double *arrgTgKg)
{
    // 1
    double arrg[2 * QUANT_GDG_PRMS] = {0.};
    calc_dViu_po_dX_Var1_(arrX,marrXZv_IU_PSK, arrVZv_RLK,arrg);
   /* for (int i =0; i < 2; ++i)
    {
        for (int j=0; j < 3; ++j)
        {
          arrg[i *QUANT_GDG_PRMS + j]  =  arrg[i *QUANT_GDG_PRMS + j]/COeff ;
        }
    }*/
    ///
    double arrT3[2 *QUANT_GDG_PRMS] = {0.}, arrK[4] = {0.};
    arrK[0] =arrK[3] = mAngleSigIU * mAngleSigIU;
    MtrxTranspMultMatrx(arrg,2, QUANT_GDG_PRMS, arrK,2, arrT3) ;
    MtrxMultMatrx(arrT3,QUANT_GDG_PRMS, 2, arrg,QUANT_GDG_PRMS, arrgTgKg) ;


}
 //---------------------------------------------------------

//-------------------------------------------------------------
bool QContObj6::calcCorMtrxAngs_Var1(QContObj6 *pContObjectArr,const int QuantObj
            ,double *parrVZv_RLK,double *parrVZv_IU, double *arrXEst
            , double *arrMtrx_dFGr_po_dX_Inv, double *arrMtrxCor)
{

   double arrSum1[NUmAng* NUmAng] = {0.};
   double arrSum2[NUmAng* NUmAng] = {0.};
   double arrQKQT_Cur[NUmAng* NUmAng] = {0.};
   double arrgTgKggT_Cur[NUmAng* NUmAng] = {0.};
   for(int i =0; i < QuantObj; ++i)
   {
     pContObjectArr[i].calcMtrxAngs_QKQT(&(parrVZv_RLK[3 * i]),arrXEst,arrQKQT_Cur);

     MtrxSumMatrx(arrSum1, arrQKQT_Cur,NUmAng, NUmAng, arrSum1) ;

     pContObjectArr[i].calcMtrxAngs_gTKg(&(parrVZv_RLK[3 * i])
             ,&(parrVZv_IU[2 * i]),arrXEst,arrgTgKggT_Cur);
     MtrxSumMatrx(arrSum2, arrgTgKggT_Cur,NUmAng, NUmAng, arrSum2) ;

   }
   double arrK1[NUmAng * NUmAng] = {0.};
   double arrK2[NUmAng * NUmAng] = {0.};
   double arrMtrx_dFGr_po_dX_Inv1[NUmAng * NUmAng] = {0.};

   memcpy(arrMtrx_dFGr_po_dX_Inv1, arrMtrx_dFGr_po_dX_Inv, NUmAng * NUmAng * sizeof(double));
   MtrxMultMatrx_MultMatrxTransp(arrMtrx_dFGr_po_dX_Inv,arrSum1
                                 ,arrMtrx_dFGr_po_dX_Inv1,NUmAng, arrK1);

   MtrxMultMatrx_MultMatrxTransp(arrMtrx_dFGr_po_dX_Inv,arrSum2
                                 ,arrMtrx_dFGr_po_dX_Inv1,NUmAng, arrK2);

   MtrxSumMatrx(arrK1, arrK2,NUmAng, NUmAng, arrMtrxCor) ;


}
//-----------------------------------------------------------------
void QContObj6::calcMtrxAngs_QKQT(double *arrVZv_RLK, double *arrX,double *arrQKQT)
{
   // 1
   double arrg[2 * NUmAng] = {0.},arrgBig[2 * QUANT_GDG_PRMS] = {0.};
   calc_dViu_po_dX_Var1_(arrX,marrXZv_IU_PSK, arrVZv_RLK,arrgBig);
   for (int i =0; i < 2; ++i)
   {
       for (int j=0; j < NUmAng; ++j)
       {
         arrg[i *NUmAng + j]  =  arrgBig[i *QUANT_GDG_PRMS + 3 + j] ;
       }
   }

   ///

   // 2 вырезание первых 3 столбцов
   double arrT0[6] = {0.};
   for (int i =0; i < 2; ++i)
   {
       for (int j=0; j < 3; ++j)
       {
         arrT0[i * 3 +j] =  arrgBig[i *QUANT_GDG_PRMS + j]/COeff ;
       }
   }

   // 3 ваычисление матрицы M
   double matrPereh_ASK_V_PSK_RLK[9] = {0.};
   calcMtrx3_ASPK_v_PSK(&(arrX[3]),matrPereh_ASK_V_PSK_RLK) ;
   ///

   // 4 вычисление dA_po_dV
   double arr_dA_po_dV[9]={0.};
   calc_dRectangularl_po_dSpherical(arrVZv_RLK, arr_dA_po_dV);
   ///

   double arrT1[6] = {0.},arrT2[6] = {0.}, arrQ[3 *NUmAng] = {0.};
   MtrxMultMatrx(arrT0,2, 3, matrPereh_ASK_V_PSK_RLK,3, arrT1) ;
   MtrxMultMatrx(arrT1,2, 3, arr_dA_po_dV,3, arrT2) ;
   MtrxTranspMultMatrx(arrg,2, NUmAng, arrT2,3, arrQ) ;


   double arrT3[3 *NUmAng] = {0.};
   MtrxMultMatrx(arrQ,NUmAng, 3, marrCorMtrxRLK,3, arrT3) ;
   MtrxMultMatrxTransp(arrT3,NUmAng, 3, arrQ,NUmAng, arrQKQT) ;

}
//-----------------------------------------------------------------
void QContObj6::calcMtrxAngs_gTKg(double *arrVZv_RLK
             ,double *arrVZv_IU, double *arrX, double *arrgTgKg)
{
    // 1
       double arrg[2 * NUmAng] = {0.},arrgBig[2 * QUANT_GDG_PRMS] = {0.};
       calc_dViu_po_dX_Var1_(arrX,marrXZv_IU_PSK, arrVZv_RLK,arrgBig);
       for (int i =0; i < 2; ++i)
       {
           for (int j=0; j < NUmAng; ++j)
           {
             arrg[i *NUmAng + j]  =  arrgBig[i *QUANT_GDG_PRMS + 3 + j] ;
           }
       }

   ///
   double arrT3[2 *NUmAng] = {0.}, arrK[4] = {0.};
   arrK[0] =arrK[3] = mAngleSigIU * mAngleSigIU;
   MtrxTranspMultMatrx(arrg,2, NUmAng, arrK,2, arrT3) ;
   MtrxMultMatrx(arrT3,NUmAng, 2, arrg,NUmAng, arrgTgKg) ;


}
//---------------------------------------------------------
double QContObj6::calcSumNeviazka_Var1(double *arrXCur,QContObj6 *pContObjectArr
             ,const int QuantObj ,double *parrVZv_RLK
             ,double *parrVZv_IU)
 {
    double valSum = 0.;
    for (int i =0; i < QuantObj; ++i)
    {
        double rez = pContObjectArr[i].calcNeviazka_Var1(arrXCur
         ,&(parrVZv_RLK[ 3 * i])
        ,&(parrVZv_IU[ 2 * i]));
        valSum += rez;
    }
     return valSum;
 }
//---------------------------------------------------------

//---------------------------------------------------------
bool QContObj6::doOneIteration_Var1_New(double *arrX0,QContObj6 *pContObjectArr
            ,const int QuantObj,double *parrVZv_RLK
            ,double *parrVZv_IU, double *arrDel, double *arrMtrx_dF_po_dX_Inv)
{

    double arrF [QUANT_GDG_PRMS] = {0.};
    double arrMtrx_dF_po_dX [QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.};
    calc_FGr_and_dFGr_po_dX_Var1(arrX0,pContObjectArr, QuantObj,parrVZv_RLK
                        ,parrVZv_IU, arrF, arrMtrx_dF_po_dX);



    bool brez =   InverseMtrx(arrMtrx_dF_po_dX, QUANT_GDG_PRMS, arrMtrx_dF_po_dX_Inv);
    // проверка
   // double aeeT[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.};
    // MtrxMultMatrx( arrMtrx_dF_po_dX,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrMtrx_dF_po_dX_Inv,QUANT_GDG_PRMS, aeeT) ;
    ///
    if (!brez)
    {
        return false;
    }

     double arrCoeff[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.}, arrT[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.};
     fillE(arrCoeff, QUANT_GDG_PRMS);
     for (int i =0; i < 3; ++i)
     {
         arrCoeff[i *QUANT_GDG_PRMS + i] = COeff;
     }

      MtrxMultMatrx( arrCoeff,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrMtrx_dF_po_dX_Inv,QUANT_GDG_PRMS, arrT) ;
     // memcpy( arrMtrx_dF_po_dX_Inv, arrT, QUANT_GDG_PRMS * QUANT_GDG_PRMS * sizeof(double));

     // MtrxMultMatrx(arrMtrx_dF_po_dX_Inv,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrF,1, arrDel) ;
       MtrxMultMatrx(arrT,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrF,1, arrDel) ;

    return true;

}

//-------------------------------------------
// вычисление разностной производной функции f по переменной
// с номером j
//функция f зависит от 7 переменных
//от переменных вектора позиционирования РЛК и ИУ  в ПСК - arrX

//OUTPUT:
//arr_df_po_dViu[2 + BDalnomer] - вектор частных производных ыектор функции Viu
//по переменной с номером NUmj (нумерация с 0)
void QContObj6::calc_dViu_po_dXj(double *arrX,double *arrX_IU
                         , double *arrV_RLK,const int NUmj, double *arr_dViu_po_dXj)
{
    double arrX1[QUANT_GDG_PRMS] = {0.};

    memcpy(arrX1, arrX, QUANT_GDG_PRMS * sizeof(double));

    double delta = 0.;
    if (NUmj < 3)
    {
     delta = 0.5;

    }
    else
    {
        delta = 0.001;

    }
    arrX1[NUmj]+= delta;
    double arrViu[2] = {0.},arrViu1[2] = {0.}, arrT[2] = {0.},arrDelta[2] = {0.};
    calc_Viu( arrX1,arrX_IU, arrV_RLK, arrViu1);
    calc_Viu(arrX,arrX_IU, arrV_RLK, arrViu);

    MtrxMinusMatrx(arrViu1, arrViu,1, 2, arrT);
    MatrxMultScalar(arrT, 1, 2, 1./delta,arrDelta);

    for (int i =0; i < 2;++i)
    {
      arr_dViu_po_dXj[i] =  arrDelta[i] ;
    }


}
//------------------------------
//-------------------------------------------
// вычисление разностной производной функции f по переменной
// с номером j
//функция f зависит от 7 переменных
//от переменных вектора позиционирования РЛК и ИУ  в ПСК - arrX

//OUTPUT:
//arr_df_po_dViu[2 + BDalnomer] - вектор частных производных ыектор функции Viu
//по переменной с номером NUmj (нумерация с 0)
void QContObj6::calc_dViu_po_dX_Razn(double *arrX,double *arrX_IU
                         , double *arrV_RLK, double *arr_dViu_po_dX)
{
  for(int j =0; j < QUANT_GDG_PRMS; ++j)
  {
      double arrT  [2] = {0.};

      calc_dViu_po_dXj(arrX,arrX_IU
                               , arrV_RLK,j, arrT ) ;
      changeCol(arr_dViu_po_dX, 2,  QUANT_GDG_PRMS ,j ,arrT);


  }

}
//------------------------------
double QContObj6::calcNeviazka_Var1(double *arrX
       ,double *arrVZv_RLK, double *parrVZv_IU)
{
    double arrViu[2] = {0.}, arrDel[2] = {0.};

    calc_Viu(arrX, marrXZv_IU_PSK,arrVZv_RLK, arrViu);
    MtrxMinusMatrx(arrViu, parrVZv_IU,1, 2, arrDel);

    double t = NormVect(arrDel, 2 );


    // проверка
    double arrViu1[2] = {0.};
    calc_Viu(marrXTrue_RLK_PSK, marrXTrue_IU_PSK,arrVZv_RLK, arrViu1);

    return t * t;
}

//-----------------------------------------------------------------------
// вычисление оценки вектора позиционирования по углам РЛК и ИУ в ПСК-ИУ по Варианту 1
// INPUT:
//pContObjectArr - массив КО
//LEnArr - кол-во КО
// BNoise - признак имитации случайных шумов в измерениях
//OUTPUT:
//arrXEst[4] - оценка вектора углов
//arrMtrxCor[16] - коррел матрица ошибок оценивания
double QContObj6::imitate_and_estimateAngles_Var1_(QContObj6 *pContObjectArr,const int QuantObj
       , const bool BNoise , double *arrXEst, double *arrMtrxCor)
{
    // векторы измерений РЛК
 double *parrVZv_RLK = new double [QuantObj * 3 ];

 // векторы измерений ИУ
 double *parrVZv_IU = new double [QuantObj * 2 ];

 // векторы истинного положения цели в ИУ
 //double *parrVTrue_IU = new double [QuantObj * (BDalnomer +2) ];

 // имитация измерений РЛК и ИУ
  for(int i =0; i < QuantObj; ++i)
  {
    pContObjectArr[i].imitateMeasures(BNoise,&(parrVZv_RLK[i * 3])
            , &(parrVZv_IU[i * 2])) ;
  }


  double rez = estimateAngles_Var_1(pContObjectArr,QuantObj
                             ,parrVZv_RLK,parrVZv_IU, arrXEst, arrMtrxCor);

 delete []parrVZv_RLK;
 delete []parrVZv_IU;

  return rez;
}

//----------------------------------------------------
 double QContObj6::estimateAngles_Var_1(QContObj6 *pContObjectArr,const int QuantObj
              ,double *parrVZv_RLK
              ,double *parrVZv_IU , double *arrXRLKEst, double *arrMtrxCor)
 {
     // 0
     double arrMtrx_dFGr_po_dX_Inv[NUmAng * NUmAng] = {0.};
     // 1. формирование вектора начальгного приближения
     double arrXRLKRelativeEst[QUANT_GDG_PRMS]= {0.};

     memcpy(arrXRLKRelativeEst, pContObjectArr[0].marrXZv_RLK_PSK, QUANT_GDG_PRMS * sizeof(double));

     MtrxMinusMatrx(pContObjectArr[0].marrXZv_RLK_PSK, pContObjectArr[0].marrXZv_IU_PSK,3, 1, arrXRLKRelativeEst);

     ///

     // 2. вычисление начальной невязки
     double valNeviaz0 = calcSumNeviazka_Var1(arrXRLKRelativeEst,pContObjectArr
            , QuantObj, parrVZv_RLK,parrVZv_IU);
     ///

     // 3. итреац процесс
     double arrXCur[QUANT_GDG_PRMS]= {0.};
     memcpy(arrXCur, arrXRLKRelativeEst, QUANT_GDG_PRMS * sizeof(double));
     double valNeviaz1=0.;
     double *arrDel = new double[NUmAng];
     for (int i =0; i < 500; ++i)
     {

       doOneIterationAngles_Var1(arrXRLKRelativeEst,pContObjectArr, QuantObj
                               ,parrVZv_RLK,parrVZv_IU, arrDel,arrMtrx_dFGr_po_dX_Inv);
         double coef = 1.;//0.25;
         MatrxMultScalar( arrDel, NUmAng, 1, coef, arrDel);
         MtrxMinusMatrx(&(arrXRLKRelativeEst[3]), arrDel,NUmAng, 1, &(arrXCur[3]));
        valNeviaz1 = calcSumNeviazka_Var1(arrXCur,pContObjectArr, QuantObj
                    ,parrVZv_RLK,parrVZv_IU);
       //double arrDel [ 5] ={0.};
      // MtrxMinusMatrx(arrXCur, arrXRLKRelativeEst,1, 5, arrDel);
       double del = NormVect(arrDel, NUmAng);
       if (valNeviaz1 > valNeviaz0)
       {

           break;
       }
       if ((del < 1.E-8))
       {
           memcpy(arrXRLKRelativeEst, arrXCur, QUANT_GDG_PRMS * sizeof(double));
           valNeviaz0 = valNeviaz1;
           break;
       }
       memcpy(arrXRLKRelativeEst, arrXCur, QUANT_GDG_PRMS * sizeof(double));
       valNeviaz0 = valNeviaz1;
     }

     delete []arrDel;
     memcpy(arrXRLKEst, arrXRLKRelativeEst,QUANT_GDG_PRMS * sizeof(double));
     arrXRLKEst[0] = arrXRLKRelativeEst [0] + pContObjectArr[0].marrXZv_IU_PSK[0];
     arrXRLKEst[1] = arrXRLKRelativeEst [1] + pContObjectArr[0].marrXZv_IU_PSK[1];
     arrXRLKEst[2] = arrXRLKRelativeEst [2] + pContObjectArr[0].marrXZv_IU_PSK[2];


     double arrMtrxCor0[NUmAng * NUmAng] = {0.};
      calcCorMtrxAngs_Var1(pContObjectArr, QuantObj
                  ,parrVZv_RLK,parrVZv_IU, arrXRLKRelativeEst,arrMtrx_dFGr_po_dX_Inv, arrMtrxCor0);
      memset(arrMtrxCor, 0,QUANT_GDG_PRMS *QUANT_GDG_PRMS * sizeof(double));
      for (int i =0; i < NUmAng; ++i )
      {
          for (int j =0; j < NUmAng; ++j)
          {
             arrMtrxCor[( 3 + i) *QUANT_GDG_PRMS + 3 + j]  = arrMtrxCor0[i *  NUmAng + j];
          }
      }


     return valNeviaz0;


 }

 //---------------------------------------------------------

 //------------------------------
//---------------------------------------------------------
 double QContObj6::doOneIterationAngles_Var1(double *arrX0,QContObj6 *pContObjectArr
             ,const int QuantObj,double *parrVZv_RLK
             ,double *parrVZv_IU,double *arrDel, double *arrMtrx_dFGr_po_dX_Inv)
 {
     // к-во оцениваемых углов
    // const int NUmAng = 4;
     double arrF[NUmAng] = {0.};
     double arrMtrx_dF_po_dX [NUmAng * NUmAng]= {0.};
     calc_FGrAngs_and_dFGr_po_dX_Var1(arrX0,pContObjectArr, QuantObj ,parrVZv_RLK
                         ,parrVZv_IU, arrF, arrMtrx_dF_po_dX);

     //double arrMtrx_dF_po_dX_Inv [NUmAng * NUmAng]= {0.};


     bool brez =   InverseMtrx(arrMtrx_dF_po_dX, NUmAng, arrMtrx_dFGr_po_dX_Inv);

     // проверка
   //  double aeeT[25] = {0.};
     // MtrxMultMatrx( arrMtrx_dF_po_dX,2, 2, arrMtrx_dF_po_dX_Inv,2, aeeT) ;
     ///

     if (!brez)
     {
         return false;
     }    

     MtrxMultMatrx( arrMtrx_dFGr_po_dX_Inv,NUmAng, NUmAng, arrF,1, arrDel) ;

     return true;
 }

 //------------------------------
 //---------------------------------------------
 void  QContObj6::calc_FGrAngs_and_dFGr_po_dX_Var1(double *arrX0,QContObj6 *pContObjectArr
                     ,const int QuantObj,double *parrVZv_RLK
                     ,double *parrVZv_IU,double *parrF,double *parrMtrx_dF_po_dX)
 {
     memset(parrF, 0, NUmAng * sizeof(double));
     memset(parrMtrx_dF_po_dX, 0, NUmAng *NUmAng * sizeof(double));

     for (int i = 0; i < QuantObj; ++i)
     {
             double parrf [NUmAng]= {0.};
             double parr_df_po_dX [NUmAng * NUmAng]= {0.};
             pContObjectArr[i].calc_fiAng_and_dfiAng_po_dx_Var1(arrX0, &(parrVZv_RLK[3 * i])
                     ,&(parrVZv_IU[2 * i]), parrf,parr_df_po_dX);
             MtrxSumMatrx(parrF, parrf,1, NUmAng, parrF) ;
             MtrxSumMatrx(parrMtrx_dF_po_dX, parr_df_po_dX,NUmAng, NUmAng, parrMtrx_dF_po_dX) ;

     }
 }
 //-------------------------------------------------
 //-------------------------------------------------
 void  QContObj6::calc_fiAng_and_dfiAng_po_dx_Var1(double *arrX0, double *arrVZv_RLK
           ,double *arrVZv_IU,double *arrfi,double *arrMtrx_dfi_po_dX)
 {
   //1 вычисление вектор-функции fi
     // 1.1 векторт невязок координат ИУ f-Viu
     double arrDelViu[2 ]= {0.};
     double arrf [2 ] = {0.};
     calc_Viu( arrX0,marrXZv_IU_PSK, arrVZv_RLK, arrf);
     MtrxMinusMatrx(arrf, arrVZv_IU,1, 2,  arrDelViu);


    ///

    //1.2 матрица частных проихводных Viu по X. обозначается g
     double arr_gBig[ 2* QUANT_GDG_PRMS ]= {0.};
     double arr_g[ 2 * NUmAng ] = {0.};
     calc_dViu_po_dX_Var1_( arrX0, marrXZv_IU_PSK,arrVZv_RLK, arr_gBig);
     for (int i =0; i < 2; ++i)
         for (int j =0; j < NUmAng; ++j)
     {
        arr_g[i *  NUmAng + j] = arr_gBig[i *QUANT_GDG_PRMS + 3 + j ];
     }
     ///


     ////////////////////////////////////
     // PROVERKA!!  /////////////////////
     double arr_gBig1[ 2 * QUANT_GDG_PRMS ]= {0.};
     calc_dViu_po_dX_Razn(arrX0, marrXZv_IU_PSK,arrVZv_RLK, arr_gBig1);
     double *arr_gv= new double[ 2 * NUmAng ];
     for (int i =0; i < 2; ++i)
         for (int j =0; j < NUmAng; ++j)
     {
        arr_gv[i *  NUmAng + j] = arr_gBig1[i *QUANT_GDG_PRMS + 3 + j ];
     }
    ///

     //1.3 вычисление arrf
     MtrxTranspMultMatrx( arr_g, 2 ,NUmAng, arrDelViu,1, arrfi) ;
     ///

 //2 Вычисление матрицы dfi_po_dx
     // 2.1
     double arr_g1[ NUmAng *2 ]= {0.};
     memcpy(arr_g1, arr_g,  NUmAng * 2 * sizeof(double));
     MtrxTranspMultMatrx( arr_g, 2,NUmAng
                          , arr_g1,NUmAng, arrMtrx_dfi_po_dX) ;
     ///
/*
   for(int j =0; j<NUmAng; ++j)
     {
         double arr_dg_po_dXj_Big = new double[QUANT_GDG_PRMS *(2 + BDalnomer) ] ;
          calcMatr_dg_po_dXj(arrX0,  arrVZv_RLK,j,  arr_dg_po_dXj_Big);

         double *arr_dg_po_dXj = new double[NUmAng *(2 + BDalnomer) ] ;
         for (int i =0; i < (2 + BDalnomer); ++i)
             for (int j =0; j < NUmAng; ++j)
         {
            arr_dg_po_dXj[i *  NUmAng + j] = arr_dg_po_dXj_Big[i *QUANT_GDG_PRMS + 3 + j ];
         }

         double *arrColumn= new double[NUmAng];

         MtrxTranspMultMatrx( arr_dg_po_dXj, (2 + BDalnomer),NUmAng ,arrDelViu,1, arrColumn) ;
         for(int i =0; i < NUmAng;++i)
         {
           arrMtrx_dfi_po_dX[i * NUmAng + j] += arrColumn[i];
         }

         delete []arr_dg_po_dXj;
         delete []arrColumn;
         delete []arr_dg_po_dXj_Big;
     }
     */



 }
 //-------------------------------------------------
 //Вычисление матрицы частных производных вектор-функции целеуказаний
 // по вектору параметров позиционирования РЛК
 //INPUT:
 //BDalnomer- признак дальномера
 //arrXCur[6] - вектор параметров позиционирования
 //arrVZv_RLK[3] - положение цели (КО) в АСфСК
 //OUTPUT:
//  arr_g[ (BDalnomer +2)*5]
 // вектор параллакса - относитльный в ПСК-ИУ
void  QContObj6::calc_dViu_po_dX_Var1_(double *arrXCur,double *arrX_IU
                           , double *arrVZv_RLK,double *arr_g)
{


 // 1. пересчет замера в АСПК
 double arrS_ASPKZv[3] = {0.};
 recalcSphericalCoord_INTO_Rectangular(arrVZv_RLK[1],arrVZv_RLK[0],arrVZv_RLK[2]
                                            , arrS_ASPKZv);
 ///
 // 2. пересчет замера в ПСК-РЛК

 double arrS_PSK_RLK_Zv[3] = {0.}, matrPereh_ASK_V_PSK_RLK[9] = {0.};
 calcMtrx3_ASPK_v_PSK(&(arrXCur[3]),matrPereh_ASK_V_PSK_RLK) ;
 MtrxMultMatrx(matrPereh_ASK_V_PSK_RLK,3, 3, arrS_ASPKZv,1, arrS_PSK_RLK_Zv) ;
 ///

 // 3. пересчет замера в ПСК-ИУ
 double arrS_PSK_IU_Zv[3] = {0.};
 MtrxSumMatrx(arrS_PSK_RLK_Zv, arrXCur,1, 3, arrS_PSK_IU_Zv) ;
 ///

 // 4. пересчет замера в ИУСПК
 double arrS_PSK_IUSPK_Zv[3] = {0.}, matrPereh_ASK_V_PSK_IU[9] = {0.};

 calcMtrx3_ASPK_v_PSK(&(arrX_IU[3]),matrPereh_ASK_V_PSK_IU) ;
 MtrxTranspMultMatrx(matrPereh_ASK_V_PSK_IU,3, 3,arrS_PSK_IU_Zv,1, arrS_PSK_IUSPK_Zv) ;
 ///

 // 5.вычисление dB_po_dS
 double arr_dB_po_dS[6] = {0.}, arr_dB_po_dS_MultC[6] = {0.};
 calc_dSpherical_po_Rectangular_AnglesOnly(arrS_PSK_IUSPK_Zv, arr_dB_po_dS);
 calc_dSpherical_po_Rectangular_AnglesOnly_MultCoeff(COeff, arrS_PSK_IUSPK_Zv, arr_dB_po_dS_MultC);

 ///

 // 6. вычисление первых трех столбцов матрицы С и внесение их на место
 double arrBlok1[9];
 MtrxMultMatrxTransp(arr_dB_po_dS_MultC,2, 3, matrPereh_ASK_V_PSK_IU,3, arrBlok1) ;
 for (int i =0; i <2; ++i)
     for (int j =0; j < 3; ++j)
     {
         arr_g[ i * QUANT_GDG_PRMS + j] = arrBlok1[i * 3 + j];
     }
 double arrTT0[2] = {0.}, arrTT1[2] = {0.}, arrTT2[2] = {0.};
 calc_dViu_po_dXj(arrXCur,arrX_IU, arrVZv_RLK,0, arrTT0);
 calc_dViu_po_dXj(arrXCur,arrX_IU, arrVZv_RLK,1, arrTT1);
 calc_dViu_po_dXj(arrXCur,arrX_IU, arrVZv_RLK,2, arrTT2);
 ///

 // 7.вычисление 4-го столбца матрицы С и внесение его на место
 double arr_dM_RLK_po_dBet[9] = {0.};

 calc_dMtrx3_ASPK_v_PSK_po_dBet(&(arrXCur[3]),arr_dM_RLK_po_dBet);

 double arrTemp0[9] = {0.}, arrTemp1[3] = {0.},  arrAlf4[3] = {0.};
 MtrxTranspMultMatrx(matrPereh_ASK_V_PSK_IU,3, 3, arr_dM_RLK_po_dBet,3, arrTemp0) ;
 MtrxMultMatrx( arrTemp0,3, 3, arrS_ASPKZv,1, arrTemp1) ;
 MtrxMultMatrx( arr_dB_po_dS, 2, 3, arrTemp1,1, arrAlf4) ;

 for (int i =0; i < 2; ++i)
 {
    arr_g[ i * QUANT_GDG_PRMS + 3] =arrAlf4[i];
 }

calc_dViu_po_dXj(arrXCur,arrX_IU, arrVZv_RLK,3, arrTT2);
 ///

 // 8.вычисление 5-го столбца матрицы С и внесение его на место
 double arr_dM_RLK_po_dEps[9] = {0.};

 calc_dMtrx3_ASPK_v_PSK_po_dEps(&(arrXCur[3]),arr_dM_RLK_po_dEps);

 double arrTemp2[9] = {0.}, arrTemp3[3] = {0.},  arrAlf5[3] = {0.};
 MtrxTranspMultMatrx(matrPereh_ASK_V_PSK_IU,3, 3, arr_dM_RLK_po_dEps,3, arrTemp2) ;
 MtrxMultMatrx( arrTemp2,3, 3, arrS_ASPKZv,1, arrTemp3) ;
 MtrxMultMatrx( arr_dB_po_dS, 2, 3, arrTemp3,1, arrAlf5) ;

 for (int i =0; i < 2; ++i)
 {
    arr_g[ i * QUANT_GDG_PRMS + 4] =arrAlf5[i];
 }
calc_dViu_po_dXj(arrXCur,arrX_IU, arrVZv_RLK,4, arrTT2);
 ///


 // 10.вычисление 6-го столбца матрицы С и внесение его на место
 double arr_dM_RLK_po_dAlf [9] = {0.};

 calc_dMtrx3_ASPK_v_PSK_po_dAlf(&(arrXCur[3]),arr_dM_RLK_po_dAlf);

 double arrTemp5[9] = {0.}, arrTemp6[3] = {0.},  arrAlf7[3] = {0.};
 MtrxTranspMultMatrx(matrPereh_ASK_V_PSK_IU,3, 3, arr_dM_RLK_po_dAlf,3, arrTemp5) ;
 MtrxMultMatrx( arrTemp5,3, 3, arrS_ASPKZv,1, arrTemp6) ;
 MtrxMultMatrx( arr_dB_po_dS, 2, 3, arrTemp6,1, arrAlf7) ;

 for (int i =0; i < 2; ++i)
 {
    arr_g[ i * QUANT_GDG_PRMS + 5] =arrAlf7[i];
 }
calc_dViu_po_dXj(arrXCur,arrX_IU, arrVZv_RLK,5, arrTT2);
 ///

}
//---------------------------------------------

//---------------------------------------------
void  QContObj6::calc_FGr_and_dFGr_po_dX_Var1(double *arrX0,QContObj6 *pContObjectArr
                    ,const int QuantObj,double *parrVZv_RLK
                    ,double *parrVZv_IU,double *parrF,double *parrMtrx_dF_po_dX)
{
    memset(parrF, 0, QUANT_GDG_PRMS * sizeof(double));
    memset(parrMtrx_dF_po_dX, 0, QUANT_GDG_PRMS *QUANT_GDG_PRMS * sizeof(double));

    for (int i = 0; i < QuantObj; ++i)
    {
        double parrf  [QUANT_GDG_PRMS] = {0.};
        double parr_df_po_dX [QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.};
            pContObjectArr[i].calc_fi_and_dfi_po_dx_Var1(arrX0, &(parrVZv_RLK[3 * i])
                    ,&(parrVZv_IU[2 * i]), parrf,parr_df_po_dX);
            MtrxSumMatrx(parrF, parrf,1, QUANT_GDG_PRMS, parrF) ;
            MtrxSumMatrx(parrMtrx_dF_po_dX, parr_df_po_dX,QUANT_GDG_PRMS, QUANT_GDG_PRMS, parrMtrx_dF_po_dX) ;

    }
}
//-------------------------------------------------
void  QContObj6::calc_fi_and_dfi_po_dx_Var1(double *arrX0, double *arrVZv_RLK
          ,double *arrVZv_IU,double *arrfi,double *arrMtrx_dfi_po_dX)
{
  //1 вычисление вектор-функции fi
    // 1.1 векторт невязок координат ИУ f-Viu
    double arrDelViu [2] = {0.};
    double arrf [2] = {0.};
    calc_Viu(arrX0,marrXZv_IU_PSK, arrVZv_RLK, arrf);
    MtrxMinusMatrx(arrf, arrVZv_IU,1, 2 ,  arrDelViu);

   ///

   //1.2 матрица частных проихводных Viu по X. обозначается g
    double arr_g[ 2 * QUANT_GDG_PRMS ] = {0.};
    calc_dViu_po_dX_Var1_( arrX0, marrXZv_IU_PSK,arrVZv_RLK, arr_g);
    ///

    // проверка
    // double arr_gRazn[ 2 * QUANT_GDG_PRMS ] = {0.};
   // calc_dViu_po_dX_Razn(arrX0, marrXZv_IU_PSK,arrVZv_RLK, arr_gRazn);

    //1.3 вычисление arrf
    MtrxTranspMultMatrx( arr_g, 2 ,QUANT_GDG_PRMS, arrDelViu,1, arrfi) ;

    ///

//2 Вычисление матрицы dfi_po_dx
    // 2.1
     double arr_g1[ 2 * QUANT_GDG_PRMS ] = {0.};
    memcpy(arr_g1, arr_g,  QUANT_GDG_PRMS *2 * sizeof( double));
    MtrxTranspMultMatrx( arr_g, 2,QUANT_GDG_PRMS
                         , arr_g1,QUANT_GDG_PRMS, arrMtrx_dfi_po_dX) ;
    ///
/*
  for(int j =0; j<QUANT_GDG_PRMS; ++j)
    {
         double *arr_dg_po_dXj = new  double[QUANT_GDG_PRMS *(2 + BDalnomer) ] ;
        calcMatr_dg_po_dXj(arrX0,  BDalnomer ,arrVZv_RLK,j,  arr_dg_po_dXj);
         double *arrColumn= new  double[QUANT_GDG_PRMS];

        MtrxTranspMultMatrx( arr_dg_po_dXj, (2 + BDalnomer),QUANT_GDG_PRMS ,arrDelViu,1, arrColumn) ;
        for(int i =0; i < QUANT_GDG_PRMS;++i)
        {
          arrMtrx_dfi_po_dX[i * QUANT_GDG_PRMS + j] += arrColumn[i];
        }

        delete []arr_dg_po_dXj;
        delete []arrColumn;
    }
*/

}
//------------------------------------------------------------
void QContObj6::calcMatr_dg_po_dXj(double *arrX0, double *arrVZv_RLK
                   ,const int NUmj,  double *arr_dg_po_dXj)
{


    double arrX1[QUANT_GDG_PRMS] = {0.};

    memcpy(arrX1, arrX0, QUANT_GDG_PRMS * sizeof(double));

    double delta = 0.;
    if (NUmj < 3)
    {
    delta = 0.05;
    }
    else
    {
    delta = 0.0001;
    }
    arrX1[NUmj]+= delta;
    double arr_g [QUANT_GDG_PRMS *2] = {0.};
    double arr_g1 [QUANT_GDG_PRMS *2] = {0.};
    double arrT [QUANT_GDG_PRMS *2] = {0.};

    calc_dViu_po_dX_Var1_( arrX0,marrXZv_IU_PSK, arrVZv_RLK, arr_g);
    calc_dViu_po_dX_Var1_( arrX1, marrXZv_IU_PSK,arrVZv_RLK, arr_g1);

    MtrxMinusMatrx(arr_g1, arr_g,QUANT_GDG_PRMS, 2, arrT);
    MatrxMultScalar(arrT, QUANT_GDG_PRMS, 2, 1./delta,arr_dg_po_dXj);


}

//---------------------------------------------------

//---------------------------------------------------
//вычисление координат цели, измеренных в АСфСК, в системе прямоуголтных координат ИУ
void QContObj6::calc_S_IUSPK(double *arrXRelative_RLK,double *arrX_IU
                         , double *arrVZv_RLK, double *arrS_IUSPK)
{

    ///
    // 1.  положение в АСПК
    double arrS_ASPK[3] = {0.};
    recalcSphericalCoord_INTO_Rectangular(arrVZv_RLK[1],arrVZv_RLK[0]
                                               ,arrVZv_RLK[2], arrS_ASPK);
    ///

    // 2.  положение  в ПСК-РЛК
    double matrPereh_ASK_V_PSK[9]= {0.};
    calcMatr_ASK_v_PSK(&(arrXRelative_RLK[3]),matrPereh_ASK_V_PSK) ;
    double arrS_PSK_RLK[3] = {0.};
    MtrxMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrS_ASPK,1, arrS_PSK_RLK) ;
    ///

    // 3. положение в ПСК-ИУ
    double arrS_PSK_IU[3] = {0.},  arrPWave[3] ={0.};
    MtrxSumMatrx(arrS_PSK_RLK, arrXRelative_RLK,1, 3, arrS_PSK_IU) ;
    ///

    // 4. положение в ИУСПК
    double matrPereh_ASK_V_PSK1[9]= {0.};
    calcMatr_ASK_v_PSK(&(arrX_IU[3]),matrPereh_ASK_V_PSK1) ;

    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK1,3, 3, arrS_PSK_IU,1, arrS_IUSPK) ;
    ///

}

//-------------------------------------------------------

void QContObj6::calc_RelativeEilers(double *arrEilers_RLKd,double *arrEilers_IUd
                         ,double *arrEilers_Relatived)
{
    long double arrEilers_RLK[3] = {0.},arrEilers_IU[3] ={0.} ;
    for (int i = 0; i < 3; ++i)
    {
      arrEilers_RLK[i] = (long double)arrEilers_RLKd[i];
      arrEilers_IU[i]  = (long double)arrEilers_IUd[i];
    }
    long double  arrM_RLK_V_PSK[9] = {0.}, arrM_IU_V_PSK[9] = {0}, arrM[9] = {0.};
    calcMtrx3_ASPK_v_PSK(arrEilers_RLK,arrM_RLK_V_PSK) ;
    calcMtrx3_ASPK_v_PSK(arrEilers_IU,arrM_IU_V_PSK) ;
    MtrxTranspMultMatrx(arrM_IU_V_PSK,3, 3, arrM_RLK_V_PSK,3, arrM) ;

    long double eps1 = asinl((long double)(arrM[7]));
    long double alf1 = asinl(-(long double)(arrM[6])/ cosl(eps1));
    long double bet1 = asinl((long double)(arrM[1])/ cosl(eps1));

    arrEilers_Relatived[0]= (double)bet1;
    arrEilers_Relatived[1]= (double)eps1;
    arrEilers_Relatived[2]= (double)alf1;


}
//-------------------------------------------------------

void QContObj6::calc_RelativePosition(double *arrX6_RLK,double *arrX6_IU
                                        , double *arrX6_Relative)
{
calc_RelativeEilers(&(arrX6_RLK[3]),&(arrX6_IU[3]),&(arrX6_Relative[3]));
calc_RelativeParalacs(arrX6_RLK,arrX6_IU, arrX6_Relative);

}
//-----------------------------------------
void QContObj6::calc_RelativeParalacs(double *arrX6_RLK,double *arrX6_IU
                                        , double *arrVectParalacs_Relative)
{
    double  arrM_IU_V_PSK[9] = {0};

    calcMatr_PSK_v_KGSK(&(arrX6_IU[3]),arrM_IU_V_PSK) ;

    double arrT[3] ={0.};
    MtrxMinusMatrx(arrX6_RLK, arrX6_IU,3, 1, arrT);

    MtrxTranspMultMatrx( arrM_IU_V_PSK,3, 3, arrT,1, arrVectParalacs_Relative) ;

}

//---------------------------------------------
void setDblArray(long double *arrInp, const int LEN, double *arrOut)
{
    for (int i =0 ; i < LEN; ++i)
    {
       arrOut[i] = (double)arrInp[i];
    }
}
void setLongDblArray( double *arrInp, const int LEN,long double *arrOut)
{
    for (int i =0 ; i < LEN; ++i)
    {
       arrOut[i] = (long double)arrInp[i];
    }
}

//---------------------------------------------------------
double QContObj6::calcSumNeviazka_Var2(double *arrXCur,QContObj6 *pContObjectArr
            ,const int QuantObj ,double *parrVZv_IU)
{
   double valSum = 0.;
   for (int i =0; i < QuantObj; ++i)
   {
       double rez = pContObjectArr[i].calcNeviazka_Var2(arrXCur
                                                    ,&(parrVZv_IU[ 2 * i])  );
       valSum += rez;
   }
    return valSum;
}
//---------------------------------------------------------

//------------------------------
double QContObj6::calcNeviazka_Var2(double *arrX,double *arrVZv_RLK)
{
   double arrV[3] = {0.};
 recalcPositionFromKGSK_to_SphericalSK(marrEilerCntrKP_Zv,marrSZv_KGSK,arrX,arrV);

 return (arrV[0]-arrVZv_RLK[0]) * (arrV[0]-arrVZv_RLK[0]) +(arrV[2]-arrVZv_RLK[1]) * (arrV[2]-arrVZv_RLK[1]);
}
//----------------------------------------------------------------
void QContObj6::calc_VGadget(double *arrEilerCntrKP,double *arrX0
                            , double *arrSZv_KGSK, double *arrf)
{
    double arrT[3] = {0.};
   recalcPositionFromKGSK_to_SphericalSK(arrEilerCntrKP,arrSZv_KGSK,arrX0,arrT);
   arrf[0] =  arrT[0];
   arrf[1] =  arrT[2];
}

//-----------------------------------------------------------------------
// вычисление оценки вектора позиционирования РЛК и ИУ в ПСК-ИУ по Варианту 1
// INPUT:
//pContObjectArr - массив КО
//LEnArr - кол-во КО
// BNoise - признак имитации случайных шумов в измерениях
//OUTPUT:
//arrXEst[5] - оценка вектора позиционирования РЛК в ИУ
//arrMtrxCor[25] - коррел матрица ошибок оценивания
double QContObj6::imitate_and_estimateAngs_Var2(QContObj6 *pContObjectArr,const int QuantObj
       , const bool BNoise, double *arrXEst, double *arrMtrxCor)
{
    // векторы измерений РЛК
 double *parrVZv_IU = new double [QuantObj * 2];

 // векторы измерений ИУ  они не нужны
 double arrVZv_RLK [3];

 // векторы истинного положения цели в ИУ
 //double *parrVTrue_IU = new double [QuantObj * (BDalnomer +2) ];

 // имитация измерений РЛК и ИУ
  for(int i =0; i < QuantObj; ++i)
  {
    pContObjectArr[i].imitateMeasures(BNoise,arrVZv_RLK
            ,&(parrVZv_IU[i * 2])) ;
  }

  double rez = estimateAngs_Var_2(pContObjectArr,QuantObj
                             ,parrVZv_IU,arrXEst, arrMtrxCor);

 delete []parrVZv_IU;


  return rez;
}


 double QContObj6::estimateAngs_Var_2(QContObj6 *pContObjectArr,const int QuantObj
             ,double *parrVZv_IU, double *arrXIUEst, double *arrMtrxCor)
 {
     // 0.
     double arrMtrx_dFGr_po_dX_Inv[NUmAng * NUmAng] = {0.};
     // 1. формирование вектора начальгного приближения
     double arrXEst[QUANT_GDG_PRMS]= {0.};

     memcpy(arrXEst, pContObjectArr[0].marrXZv_IU_PSK, QUANT_GDG_PRMS * sizeof(double));


     ///

     // 2. вычисление начальной невязки
    double valNeviaz0 = calcSumNeviazka_Var2(arrXEst,pContObjectArr, QuantObj
              ,parrVZv_IU);
     ///

     // 3. итреац процесс
    double arrXCur[QUANT_GDG_PRMS]= {0.};
    memcpy(arrXCur, arrXEst, QUANT_GDG_PRMS * sizeof(double));

    double valNeviaz1=0.;
     double arrDel [QUANT_GDG_PRMS]= {0.};
     for (int i =0; i < 500; ++i)
     {


      doOneIterationAngs_Var2(arrXEst,pContObjectArr, QuantObj
                               ,parrVZv_IU, arrDel,arrMtrx_dFGr_po_dX_Inv);
         double coef = 0.25;
         MatrxMultScalar( arrDel, QUANT_GDG_PRMS, 1, coef, arrDel);
         MtrxMinusMatrx(arrXEst, arrDel,QUANT_GDG_PRMS, 1, arrXCur);
      //  valNeviaz1 = calcSumNeviazka_Var1(arrXCur,pContObjectArr, QuantObj, BDalnomer
       //             ,parrVZv_RLK,parrVZv_IU);
       //double arrDel [ 5] ={0.};
      // MtrxMinusMatrx(arrXCur, arrXRLKRelativeEst,1, 5, arrDel);
       double del = NormVect(arrDel, QUANT_GDG_PRMS);
      // if (valNeviaz1 > valNeviaz0)
      // {

      //     break;
      // }
       if ((del < 1.E-8))
       {
           memcpy(arrXEst, arrXCur, QUANT_GDG_PRMS * sizeof(double));
          // valNeviaz0 = valNeviaz1;
           break;
       }
       memcpy(arrXEst, arrXCur, QUANT_GDG_PRMS * sizeof(double));
      // valNeviaz0 = valNeviaz1;
     }



     memcpy(arrXIUEst, arrXEst, QUANT_GDG_PRMS * sizeof(double));

    // calcCorMtrx_Var1(pContObjectArr, QuantObj, BDalnomer
       //           ,parrVZv_RLK,parrVZv_IU, arrXRLKRelativeEst, arrMtrxCor);
     double arrMtrxCor0[NUmAng * NUmAng] = {0.};
     calcCorMtrxAngs_Var2(pContObjectArr, QuantObj
                  , arrXEst, arrMtrx_dFGr_po_dX_Inv, arrMtrxCor0);

     memset(arrMtrxCor, 0,QUANT_GDG_PRMS *QUANT_GDG_PRMS * sizeof(double));
     for( int i=0; i< NUmAng; ++i)
     {
         for (int j =0; j < NUmAng; ++j)
         {
            arrMtrxCor[ ( 3 + i) * QUANT_GDG_PRMS + 3 + j] =arrMtrxCor0[ i * NUmAng + j];
         }
     }

     return 0.;//valNeviaz0;
 }
 //-------------------------------------------------------------

//-----------------------------------------------------------------

bool QContObj6::doOneIterationAngs_Var2(double *arrX0,QContObj6 *pContObjectArr
            ,const int QuantObj, double *parrVZv_IU
           , double *arrDel, double *arrMtrx_dFGr_po_dX_Inv)
{

    double arrF [NUmAng]= {0.};
    double  arrMtrx_dF_po_dX [NUmAng * NUmAng] = {0.};
    calc_FGrAngs_and_dFGrAngs_po_dX_Var2(arrX0,pContObjectArr, QuantObj, parrVZv_IU
                        , arrF, arrMtrx_dF_po_dX);


   // double valScal = 1000.;
   // MatrxMultScalar(arrMtrx_dF_po_dX, NUmAng, NUmAng, valScal,arrMtrx_dF_po_dX);
    bool brez =   InverseMtrx(arrMtrx_dF_po_dX, NUmAng, arrMtrx_dFGr_po_dX_Inv);
    // проверка
   // double aeeT[NUmAng * NUmAng] = {0.};
   //  MtrxMultMatrx( arrMtrx_dF_po_dX,NUmAng, NUmAng, arrMtrx_dF_po_dX_Inv,NUmAng, aeeT) ;
    ///
    if (!brez)
    {
        return false;
    }


  // double arrCoeff[NUmAng * NUmAng] = {0.}, arrT[NUmAng * NUmAng] = {0.};
   //fillE(arrCoeff, NUmAng);
   // MtrxMultMatrx( arrCoeff,NUmAng, NUmAng, arrMtrx_dF_po_dX_Inv,NUmAng, arrT) ;

    MtrxMultMatrx(arrMtrx_dFGr_po_dX_Inv,NUmAng, NUmAng, arrF,1, &(arrDel[3])) ;



    return true;

}

//---------------------------------------------
void  QContObj6::calc_FGrAngs_and_dFGrAngs_po_dX_Var2(double *arrX0,QContObj6 *pContObjectArr
                    ,const int QuantObj, double *parrVZv_RLK
                    ,double *parrF,double *parrMtrx_dF_po_dX)
{
    memset(parrF, 0, NUmAng * sizeof(double));
    memset(parrMtrx_dF_po_dX, 0, NUmAng *NUmAng * sizeof(double));

    for (int i = 0; i < QuantObj; ++i)
    {
        double parrf [NUmAng] = {0.};
        double parr_df_po_dX[NUmAng * NUmAng] = {0.};
            pContObjectArr[i].calc_fiAngs_and_dfiAngs_po_dx_Var2(arrX0, &(parrVZv_RLK[2 * i])
                    , parrf,parr_df_po_dX);
            MtrxSumMatrx(parrF, parrf,1, NUmAng, parrF) ;
            MtrxSumMatrx(parrMtrx_dF_po_dX, parr_df_po_dX,NUmAng, NUmAng, parrMtrx_dF_po_dX) ;

    }
}
//-------------------------------------------
//-------------------------------------------------
void  QContObj6::calc_fiAngs_and_dfiAngs_po_dx_Var2(double *arrX0, double *arrVZv_RLK
          ,double *arrfi,double *arrMtrx_dfi_po_dX)
{
  //1 вычисление вектор-функции fi
    // 1.1 векторт невязок координат ИУ f-Viu
    double arrDelV_Gdg[2] = {0.};
    double arrf [2 ] = {0.};
    calc_VGadget(marrEilerCntrKP_Zv, arrX0, marrSZv_KGSK, arrf);

    MtrxMinusMatrx(arrf, arrVZv_RLK,1, 2 ,  arrDelV_Gdg);

   ///

   //1.2 матрица частных проихводных Viu по X. обозначается g
    double arr_g[ 2 * NUmAng ] = {0.};

     double arr_gBig[ 2 * QUANT_GDG_PRMS ]= {0.};

     calc_dV_po_dX_Var2( arrX0, marrEilerCntrKP_Zv, marrSZv_KGSK, arr_gBig);
     for (int i = 0; i < 2; ++i)
     {
         for (int j =0; j < NUmAng; ++j)
         {
            arr_g[i * NUmAng +j ] =  arr_gBig[i * QUANT_GDG_PRMS + 3 + j];
         }
     }
    ///

    //1.3 вычисление arrf
    MtrxTranspMultMatrx( arr_g, 2,NUmAng, arrDelV_Gdg,1, arrfi) ;

    ///

//2 Вычисление матрицы dfi_po_dx
    // 2.1
    double arr_g1 [NUmAng *2 ] = {0.};
    memcpy(arr_g1, arr_g,  NUmAng *2 * sizeof( double));
    MtrxTranspMultMatrx( arr_g, 2,NUmAng
                         , arr_g1,NUmAng, arrMtrx_dfi_po_dX) ;

}



//-------------------------------------------------
//Вычисление матрицы частных производных вектор-функции целеуказаний
// по вектору параметров позиционирования РЛК
//INPUT:
//BDalnomer- признак дальномера
//arrXCur[7] - вектор параметров позиционирования
//arrVZv_RLK[3] - положение цели (КО) в АСфСК
//OUTPUT:
//  arr_g[ (BDalnomer +2)*5]
// вектор параллакса - относитльный в ПСК-ИУ
/*
void  QContObj6::calc_dVGdg_po_dX_Var2( double *arrX0, double *arrS_KGSK,double *arr_g)
{


// 1. пересчет замера в ПСК
   double arrS_PSK_Zv[3] = {0.};
   double matrPereh_PSK_V_KGSK[9] = {0.};
   calcMatr_PSK_v_KGSK(marrEilerCntrKP_Zv, matrPereh_PSK_V_KGSK) ;
   MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3,arrS_KGSK,1, arrS_PSK_Zv) ;
///

   //2. пересчет замера в ПСК-РЛК
     double arrS_PSK_Gdg_Zv[3] = {0.};
     MtrxMinusMatrx(arrS_PSK_Zv, arrX0,1, 3, arrS_PSK_Gdg_Zv);

// 3. пересчет замера в АСПК
  double arrS_ASPK_Zv[3] = {0.}, matrPereh_ASK_V_PSK_Gdg[9] = {0.};

calcMtrx3_ASPK_v_PSK(&(arrX0[3]),matrPereh_ASK_V_PSK_Gdg) ;
MtrxTranspMultMatrx(matrPereh_ASK_V_PSK_Gdg,3, 3, arrS_PSK_Gdg_Zv,1, arrS_ASPK_Zv) ;
///



// 5.вычисление dB_po_dS
double arr_dB_po_dS[9] = {0.};//, arr_dB_po_dS_MultC[9] = {0.};
 calc_dSpherical_po_Rectangular_AnglesOnly(arrS_ASPK_Zv, arr_dB_po_dS);
///




// 7.вычисление 1-го столбца матрицы С и внесение его на место
  double arr_dM_Gdg_po_dBet[9] = {0.};
calc_dMtrx3_ASPK_v_PSK_po_dBet(&(arrX0[3]),arr_dM_Gdg_po_dBet);

  double arrTemp0[3] = {0.},  arrAlf4[3] = {0.};
MtrxTranspMultMatrx(arr_dM_Gdg_po_dBet,3, 3, arrS_PSK_Gdg_Zv,1, arrTemp0) ;

MtrxMultMatrx( arr_dB_po_dS, 2, 3, arrTemp0,1, arrAlf4) ;

for (int i =0; i < 2; ++i)
{
   arr_g[ i * NUmAng ] =arrAlf4[i];
}


///
// 7.вычисление 2-го столбца матрицы С и внесение его на место
  double arr_dM_Gdg_po_dEps[9] = {0.};

calc_dMtrx3_ASPK_v_PSK_po_dEps(&(arrX0[3]),arr_dM_Gdg_po_dEps);
  double arrTemp1[3] = {0.},  arrAlf5[3] = {0.};
MtrxTranspMultMatrx(arr_dM_Gdg_po_dEps,3, 3, arrS_PSK_Gdg_Zv,1, arrTemp1) ;

MtrxMultMatrx( arr_dB_po_dS, 2, 3, arrTemp1,1, arrAlf5) ;

for (int i =0; i < 2; ++i)
{
   arr_g[ i * NUmAng + 1] =arrAlf5[i];
}

// 7.вычисление 3-го столбца матрицы С и внесение его на место
  double arr_dM_Gdg_po_dAlf[9] = {0.};

calc_dMtrx3_ASPK_v_PSK_po_dAlf(&(arrX0[3]),arr_dM_Gdg_po_dAlf);
  double arrTemp2[3] = {0.},  arrAlf6[3] = {0.};
MtrxTranspMultMatrx(arr_dM_Gdg_po_dAlf,3, 3, arrS_PSK_Gdg_Zv,1, arrTemp2) ;

MtrxMultMatrx( arr_dB_po_dS, 2, 3, arrTemp2,2, arrAlf6) ;

for (int i =0; i < 2; ++i)
{
   arr_g[ i * NUmAng + 2] =arrAlf6[i];
}

}
//------------------------------------------------------------

*/

//------------------------------------------------------
//-------------------------------------------------
//----------------------------------------------
//-------------------------------------------------
//-----------------------------------------------------------------------
// вычисление оценки вектора позиционирования РЛК и ИУ в ПСК-ИУ по Варианту 1
// INPUT:
//pContObjectArr - массив КО
//LEnArr - кол-во КО
// BNoise - признак имитации случайных шумов в измерениях
//OUTPUT:
//arrXEst[5] - оценка вектора позиционирования РЛК в ИУ
//arrMtrxCor[25] - коррел матрица ошибок оценивания
double QContObj6::imitate_and_estimateParams_Var2(QContObj6 *pContObjectArr,const int QuantObj
       , const bool BNoise, double *arrXEst, double *arrMtrxCor)
{
    // векторы измерений РЛК
 double *parrVZv_IU = new double [QuantObj * 2];

 // векторы измерений ИУ  они не нужны
 double arrVZv_RLK [3];

 // векторы истинного положения цели в ИУ
 //double *parrVTrue_IU = new double [QuantObj * (BDalnomer +2) ];

 // имитация измерений РЛК и ИУ
  for(int i =0; i < QuantObj; ++i)
  {
    pContObjectArr[i].imitateMeasures(BNoise,arrVZv_RLK
            ,&(parrVZv_IU[i * 2])) ;
  }

  double rez = estimateParams_Var_2(pContObjectArr,QuantObj
                             ,parrVZv_IU,arrXEst, arrMtrxCor);

 delete []parrVZv_IU;


  return rez;
}

//----------------------------------------------------
 double QContObj6::estimateParams_Var_2(QContObj6 *pContObjectArr,const int QuantObj
             ,double *parrVZv_IU, double *arrXIUEst, double *arrMtrxCor)
 {
     //0.
    double arrMtrx_dFGr_po_dX_Inv[QUANT_GDG_PRMS * QUANT_GDG_PRMS]= {0.};


     // 1. формирование вектора начальгного приближения
     double arrXEst[QUANT_GDG_PRMS]= {0.};

     memcpy(arrXEst, pContObjectArr[0].marrXZv_IU_PSK, QUANT_GDG_PRMS * sizeof(double));


     ///

     // 2. вычисление начальной невязки
    double valNeviaz0 = calcSumNeviazka_Var2(arrXEst,pContObjectArr, QuantObj,parrVZv_IU);
     ///

     // 3. итреац процесс
    double arrXCur[QUANT_GDG_PRMS]= {0.};
    memcpy(arrXCur, arrXEst, QUANT_GDG_PRMS * sizeof(double));

    double valNeviaz1=0.;
     double arrDel [QUANT_GDG_PRMS]= {0.};
     for (int i =0; i < 100; ++i)
     {

      doOneIteration_Var2(arrXEst,pContObjectArr, QuantObj
                               ,parrVZv_IU, arrDel,arrMtrx_dFGr_po_dX_Inv);
      double del = NormVect(arrDel, QUANT_GDG_PRMS);
         double coef = 0.25;
         MatrxMultScalar( arrDel, QUANT_GDG_PRMS, 1, coef, arrDel);
         MtrxMinusMatrx(arrXEst, arrDel,QUANT_GDG_PRMS, 1, arrXCur);


      // if (valNeviaz1 > valNeviaz0)
      // {

      //     break;
      // }
         memcpy(arrXEst, arrXCur, QUANT_GDG_PRMS * sizeof(double));
       if ((del < 1.E-8))
       {

           valNeviaz0 = valNeviaz1;
           break;
       }

       valNeviaz0 = valNeviaz1;
     }



     memcpy(arrXIUEst, arrXEst, QUANT_GDG_PRMS * sizeof(double));

     calcCorMtrx_Var2(pContObjectArr, QuantObj
                  , arrXIUEst,arrMtrx_dFGr_po_dX_Inv, arrMtrxCor);


     return valNeviaz0;
 }
 //-------------------------------------------------------------

 //-------------------------------------------------------------
bool QContObj6::calcCorMtrx_Var2(QContObj6 *pContObjectArr,const int QuantObj
             , double *arrXEst, double *arrMtrx_dFGr_po_dX_Inv, double *arrMtrxCor)
 {
    double arrSum1[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrSum2[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrSum3[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrJKJT_Cur[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrCKCT_Cur[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrg2TKg2_Cur[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    for(int i =0; i < QuantObj; ++i)
    {
        pContObjectArr[i].calcMtrx_JKJT(arrXEst,pContObjectArr[i].marrS_KGSK,arrJKJT_Cur);

        MtrxSumMatrx(arrSum1, arrJKJT_Cur,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrSum1) ;

        pContObjectArr[i].calcMtrx_CKCT(arrXEst, pContObjectArr[i].marrS_KGSK,arrCKCT_Cur);

        MtrxSumMatrx(arrSum2, arrCKCT_Cur,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrSum2) ;

        pContObjectArr[i].calcMtrx_g2TKg2(arrXEst, pContObjectArr[i].marrS_KGSK,arrg2TKg2_Cur);

        MtrxSumMatrx(arrSum3, arrg2TKg2_Cur,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrSum3) ;

    }
    double arrK1[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrK2[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrK3[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};
    double arrMtrx_dFGr_po_dX_Inv1[QUANT_GDG_PRMS* QUANT_GDG_PRMS] = {0.};

    memcpy(arrMtrx_dFGr_po_dX_Inv1, arrMtrx_dFGr_po_dX_Inv, QUANT_GDG_PRMS* QUANT_GDG_PRMS * sizeof(double));
    MtrxMultMatrx_MultMatrxTransp(arrMtrx_dFGr_po_dX_Inv,arrSum1
                                  ,arrMtrx_dFGr_po_dX_Inv1,QUANT_GDG_PRMS, arrK1);

    MtrxMultMatrx_MultMatrxTransp(arrMtrx_dFGr_po_dX_Inv,arrSum2
                                  ,arrMtrx_dFGr_po_dX_Inv1,QUANT_GDG_PRMS, arrK2);

    MtrxMultMatrx_MultMatrxTransp(arrMtrx_dFGr_po_dX_Inv,arrSum3
                                  ,arrMtrx_dFGr_po_dX_Inv1,QUANT_GDG_PRMS, arrK3);

    double arrMtrxCor0[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.},arrMtrxTemp0[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.};


    MtrxSumMatrx(arrK1, arrK2,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrMtrxTemp0) ;

   // MtrxSumMatrx(arrMtrxTemp0, arrK3,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrMtrxCor) ;


    MtrxSumMatrx(arrMtrxTemp0, arrK3,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrMtrxCor0) ;

    // проверка
    memcpy(arrMtrxCor0, arrK1, QUANT_GDG_PRMS * QUANT_GDG_PRMS * sizeof(double));
    ///

    double arrCoeff[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.}, arrCoeff1[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.};
    fillE(arrCoeff, QUANT_GDG_PRMS);
    fillE(arrCoeff1, QUANT_GDG_PRMS);
    for (int i =0; i < 3; ++i)
    {
        arrCoeff[ i *QUANT_GDG_PRMS + i] = COeff;
        arrCoeff1[i *QUANT_GDG_PRMS + i] = COeff;
    }
    MtrxMultMatrx_MultMatrx(arrCoeff,arrMtrxCor0,arrCoeff1,QUANT_GDG_PRMS, arrMtrxCor)  ;


}
//----------------------------------------------------
void QContObj6::calcMtrx_JKJT(double *arrX0, double *arrS_KGSK,double *arrJKJT_Cur)
{
    //0.
    double arr_g[2 * QUANT_GDG_PRMS ] = {0.};
    calc_dV_po_dX_Var2(arrX0,marrEilerCntrKP_Zv, arrS_KGSK,arr_g);


    // 1. пересчет замера в ПСК
       double arrS_PSK_Zv[3] = {0.};
       double matrPereh_PSK_V_KGSK[9] = {0.};
       calcMatr_PSK_v_KGSK(marrEilerCntrKP_Zv, matrPereh_PSK_V_KGSK) ;
       MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3,marrS_KGSK,1, arrS_PSK_Zv) ;
    ///

       //2. пересчет замера в ПСК-РЛК
         double arrS_PSK_Gdg_Zv[3] = {0.};
         MtrxMinusMatrx(arrS_PSK_Zv, arrX0,1, 3, arrS_PSK_Gdg_Zv);

    // 3. пересчет замера в АСПК
      double arrS_ASPK_Zv[3] = {0.}, matrPereh_ASK_V_PSK_Gdg[9] = {0.};
    calcMtrx3_ASPK_v_PSK(&(arrX0[3]),matrPereh_ASK_V_PSK_Gdg) ;
    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK_Gdg,3, 3, arrS_PSK_Gdg_Zv,1, arrS_ASPK_Zv) ;
    ///



    // 5.вычисление dB_po_dS
    double arr_dB_po_dS[9] = {0.};//, arr_dB_po_dS_MultC[9] = {0.};

    calc_dSpherical_po_Rectangular_AnglesOnly(arrS_ASPK_Zv, arr_dB_po_dS);
   // calc_dSpherical_po_Rectangular_AnglesOnly_MultCoeff(COeff, arrS_ASPK_Zv, arr_dB_po_dS_MultC);
    ///




    // 6. вычисление первых трех столбцов матрицы С и внесение их на место
    double arrT1[6] = {0.},arrT2[6] = {0.}, arrJ[QUANT_GDG_PRMS * 3] = {0.}
                          ,arrT3[QUANT_GDG_PRMS * 3] = {0.};
      //MtrxMultMatrxTransp(arr_dB_po_dS_MultC,(2 + BDalnomer), 3, matrPereh_ASK_V_PSK_Gdg,3, arrBlok1) ;
      MtrxMultMatrxTransp(arr_dB_po_dS,2, 3, matrPereh_ASK_V_PSK_Gdg,3, arrT1) ;
      MtrxMultMatrxTransp(arrT1,2, 3, matrPereh_PSK_V_KGSK,3, arrT2) ;
      MtrxTranspMultMatrx(arr_g,2, QUANT_GDG_PRMS, arrT2,3, arrJ) ;


      // 7.
      double arrK[9] = {0.};
      arrK[0] =arrK[3  + 1] =  mVessSKZ*mVessSKZ + mCntrlObjSKZ * mCntrlObjSKZ;
      arrK[8] =   mCntrlObjSKZ * mCntrlObjSKZ;

     MtrxMultMatrx(arrJ,QUANT_GDG_PRMS, 3, arrK,3, arrT3) ;
     MtrxMultMatrxTransp(arrT3,QUANT_GDG_PRMS, 3, arrJ,QUANT_GDG_PRMS, arrJKJT_Cur) ;

}
//-----------------------------------------------------------------
void QContObj6::calcMtrx_CKCT(double *arrX0,double *arrS_KGSK,double *arrCKCT_Cur)
{
    //0.
    double arr_g[2 * QUANT_GDG_PRMS ] = {0.};
    calc_dV_po_dX_Var2(arrX0, marrEilerCntrKP_Zv,arrS_KGSK,arr_g);

    // 1. пересчет замера в ПСК
       double arrS_PSK_Zv[3] = {0.};
       double matrPereh_PSK_V_KGSK[9] = {0.};
       calcMatr_PSK_v_KGSK(marrEilerCntrKP_Zv, matrPereh_PSK_V_KGSK) ;
       MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3,marrS_KGSK,1, arrS_PSK_Zv) ;
    ///

       //2. пересчет замера в ПСК-РЛК
         double arrS_PSK_Gdg_Zv[3] = {0.};
         MtrxMinusMatrx(arrS_PSK_Zv, arrX0,1, 3, arrS_PSK_Gdg_Zv);

    // 3. пересчет замера в АСПК
      double arrS_ASPK_Zv[3] = {0.}, matrPereh_ASK_V_PSK_Gdg[9] = {0.};
    calcMtrx3_ASPK_v_PSK(&(arrX0[3]),matrPereh_ASK_V_PSK_Gdg) ;
    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK_Gdg,3, 3, arrS_PSK_Gdg_Zv,1, arrS_ASPK_Zv) ;
    ///



    // 5.вычисление dB_po_dS
    double arr_dB_po_dS[9] = {0.};//, arr_dB_po_dS_MultC[9] = {0.};

    calc_dSpherical_po_Rectangular_AnglesOnly(arrS_ASPK_Zv, arr_dB_po_dS);
   //calc_dSpherical_po_Rectangular_AnglesOnly_MultCoeff(COeff, arrS_ASPK_Zv, arr_dB_po_dS_MultC);
    ///




    // 6. вычисление первых трех столбцов матрицы С и внесение их на место
    double arrT1[6] = {0.},arrT2[QUANT_GDG_PRMS * 3] = {0.}
                          ,arrT3[QUANT_GDG_PRMS * 3] = {0.};

      MtrxMultMatrxTransp(arr_dB_po_dS,2, 3, matrPereh_ASK_V_PSK_Gdg,3, arrT1) ;

      MtrxTranspMultMatrx(arr_g,2, QUANT_GDG_PRMS, arrT1,3, arrT2) ;

    // 7.
      double arrMiddleC[9] = {0.};
      // 8 вычисление dM/dQ*Skgsk и вставка на место
      double arr_dM_po_dQ[9] = {0.},  arrCol0[3] ={0.};
      calc_dM_psk_to_kgsk_po_dQ(marrEilerCntrKP_Zv,arr_dM_po_dQ);
      MtrxTranspMultMatrx(arr_dM_po_dQ,3, 3, arrS_KGSK,1, arrCol0) ;
      changeCol(arrMiddleC, 3,  3,0 ,arrCol0);

      // 9 вычисление dM/dPsi*Skgsk и вставка на место
      double arr_dM_po_dPsi[9] = {0.},  arrCol1[3] ={0.};
      calc_dM_psk_to_kgsk_po_dPsi(marrEilerCntrKP_Zv,arr_dM_po_dPsi);
      MtrxTranspMultMatrx(arr_dM_po_dPsi,3, 3, arrS_KGSK,1, arrCol1) ;
      changeCol(arrMiddleC, 3,  3,1 ,arrCol1);

      // 9 вычисление dM/dЗыш*Skgsk и вставка на место
      double arr_dM_po_dTet[9] = {0.},  arrCol2[3] ={0.};
      calc_dM_psk_to_kgsk_po_dTet(marrEilerCntrKP_Zv,arr_dM_po_dTet);
      MtrxTranspMultMatrx(arr_dM_po_dTet,3, 3, arrS_KGSK,1, arrCol2) ;
      changeCol(arrMiddleC, 3,  3,2 ,arrCol2);

      // 10 вычисление матрицы C
      double arrC[QUANT_GDG_PRMS * 3] = {0.};
      MtrxMultMatrx(arrT2,QUANT_GDG_PRMS, 3, arrMiddleC,3, arrC) ;

      // 11. ывчисление коррел матрицы
      double arrK[9] = {0.};
      for(int i = 0; i < 3; ++i)
      {
         arrK[3 * i + i] =  marrDispEilerCntrKP[i];
      }

      // 12 вычисление CKCT
      double arrT4[QUANT_GDG_PRMS * 3] = {0.};
      MtrxMultMatrx(arrC,QUANT_GDG_PRMS, 3, arrK,3, arrT4) ;

      MtrxMultMatrxTransp(arrT4,QUANT_GDG_PRMS, 3, arrC,QUANT_GDG_PRMS, arrCKCT_Cur) ;

}
//-----------------------------------------------------------------
void QContObj6::calcMtrx_g2TKg2(double *arrX0,double *arrS_KGSK,double *arrg2TKg2_Cur)
{
    //0.
    double arr_g[2 * QUANT_GDG_PRMS ] = {0.};
    calc_dV_po_dX_Var2(arrX0, marrEilerCntrKP_Zv,arrS_KGSK,arr_g);


    // 2. формирование матрицы K
    double arrK[4] = {0.};
    arrK[0] = arrK[3] = mAngleSigIU * mAngleSigIU;

    // 3
    double arrT0[QUANT_GDG_PRMS * 2] = {0.};
    MtrxTranspMultMatrx(arr_g,2, QUANT_GDG_PRMS, arrK,2, arrT0) ;
    MtrxMultMatrx( arrT0,QUANT_GDG_PRMS, 2, arr_g,QUANT_GDG_PRMS, arrg2TKg2_Cur) ;

}

//--------------------------------

bool QContObj6::calcCorMtrxAngs_Var2(QContObj6 *pContObjectArr,const int QuantObj
             , double *arrXEst, double *arrMtrx_dFGr_po_dX_Inv, double *arrMtrxCor)
 {
    double arrSum1[NUmAng* NUmAng] = {0.};
    double arrSum2[NUmAng* NUmAng] = {0.};
    double arrSum3[NUmAng* NUmAng] = {0.};
    double arrJKJT_Cur[NUmAng* NUmAng] = {0.};
    double arrCKCT_Cur[NUmAng* NUmAng] = {0.};
    double arrg2TKg2_Cur[NUmAng* NUmAng] = {0.};
    for(int i =0; i < QuantObj; ++i)
    {
        pContObjectArr[i].calcMtrxAngs_JKJT(arrXEst,pContObjectArr[i].marrS_KGSK,arrJKJT_Cur);

        MtrxSumMatrx(arrSum1, arrJKJT_Cur,NUmAng, NUmAng, arrSum1) ;

        pContObjectArr[i].calcMtrxAngs_CKCT(arrXEst, pContObjectArr[i].marrS_KGSK,arrCKCT_Cur);

        MtrxSumMatrx(arrSum2, arrCKCT_Cur,NUmAng, NUmAng, arrSum2) ;

        pContObjectArr[i].calcMtrxAngs_g2TKg2(arrXEst, pContObjectArr[i].marrS_KGSK,arrg2TKg2_Cur);

        MtrxSumMatrx(arrSum3, arrg2TKg2_Cur,NUmAng, NUmAng, arrSum3) ;

    }
    double arrK1[NUmAng* NUmAng] = {0.};
    double arrK2[NUmAng* NUmAng] = {0.};
    double arrK3[NUmAng* NUmAng] = {0.};
    double arrMtrx_dFGr_po_dX_Inv1[NUmAng* NUmAng] = {0.};

    memcpy(arrMtrx_dFGr_po_dX_Inv1, arrMtrx_dFGr_po_dX_Inv, NUmAng* NUmAng * sizeof(double));
    MtrxMultMatrx_MultMatrxTransp(arrMtrx_dFGr_po_dX_Inv,arrSum1
                                  ,arrMtrx_dFGr_po_dX_Inv1,NUmAng, arrK1);

    MtrxMultMatrx_MultMatrxTransp(arrMtrx_dFGr_po_dX_Inv,arrSum2
                                  ,arrMtrx_dFGr_po_dX_Inv1,NUmAng, arrK2);

    MtrxMultMatrx_MultMatrxTransp(arrMtrx_dFGr_po_dX_Inv,arrSum3
                                  ,arrMtrx_dFGr_po_dX_Inv1,NUmAng, arrK3);

    double arrMtrxTemp0[NUmAng * NUmAng] = {0.};


    MtrxSumMatrx(arrK1, arrK2,NUmAng, NUmAng, arrMtrxTemp0) ;

   // MtrxSumMatrx(arrMtrxTemp0, arrK3,NUmAng, NUmAng, arrMtrxCor) ;


    MtrxSumMatrx(arrMtrxTemp0, arrK3,NUmAng, NUmAng, arrMtrxCor) ;




}
//-----------------------------------------------------------------
//----------------------------------------------------
void QContObj6::calcMtrxAngs_JKJT(double *arrX0, double *arrS_KGSK,double *arrJKJT_Cur)
{
    //0.
    double arr_g[2 * NUmAng ] = {0.}, arr_gBig[2 * QUANT_GDG_PRMS ] = {0.};
    calc_dV_po_dX_Var2(arrX0, marrEilerCntrKP_Zv,arrS_KGSK,arr_gBig);
    for (int i =0; i < 2; ++i)
    {
        for (int j =0; j < NUmAng; ++j)
        {
           arr_g[i * NUmAng + j] =  arr_gBig[ i * QUANT_GDG_PRMS + 3 + j];
        }
    }


    // 1. пересчет замера в ПСК
       double arrS_PSK_Zv[3] = {0.};
       double matrPereh_PSK_V_KGSK[9] = {0.};
       calcMatr_PSK_v_KGSK(marrEilerCntrKP_Zv, matrPereh_PSK_V_KGSK) ;
       MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3,marrS_KGSK,1, arrS_PSK_Zv) ;
    ///

       //2. пересчет замера в ПСК-РЛК
         double arrS_PSK_Gdg_Zv[3] = {0.};
         MtrxMinusMatrx(arrS_PSK_Zv, arrX0,1, 3, arrS_PSK_Gdg_Zv);

    // 3. пересчет замера в АСПК
      double arrS_ASPK_Zv[3] = {0.}, matrPereh_ASK_V_PSK_Gdg[9] = {0.};
    calcMtrx3_ASPK_v_PSK(&(arrX0[3]),matrPereh_ASK_V_PSK_Gdg) ;
    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK_Gdg,3, 3, arrS_PSK_Gdg_Zv,1, arrS_ASPK_Zv) ;
    ///



    // 5.вычисление dB_po_dS
    double arr_dB_po_dS[9] = {0.};//, arr_dB_po_dS_MultC[9] = {0.};

    calc_dSpherical_po_Rectangular_AnglesOnly(arrS_ASPK_Zv, arr_dB_po_dS);
   // calc_dSpherical_po_Rectangular_AnglesOnly_MultCoeff(COeff, arrS_ASPK_Zv, arr_dB_po_dS_MultC);
    ///




    // 6. вычисление первых трех столбцов матрицы С и внесение их на место
    double arrT1[6] = {0.},arrT2[6] = {0.}, arrJ[NUmAng * 3] = {0.}
                          ,arrT3[NUmAng * 3] = {0.};
      //MtrxMultMatrxTransp(arr_dB_po_dS_MultC,(2 + BDalnomer), 3, matrPereh_ASK_V_PSK_Gdg,3, arrBlok1) ;
      MtrxMultMatrxTransp(arr_dB_po_dS,2, 3, matrPereh_ASK_V_PSK_Gdg,3, arrT1) ;
      MtrxMultMatrxTransp(arrT1,2, 3, matrPereh_PSK_V_KGSK,3, arrT2) ;
      MtrxTranspMultMatrx(arr_g,2, NUmAng, arrT2,3, arrJ) ;


      // 7.
      double arrK[9] = {0.};
      arrK[0] =arrK[3  + 1] =  mVessSKZ*mVessSKZ + mCntrlObjSKZ * mCntrlObjSKZ;
      arrK[8] =   mCntrlObjSKZ * mCntrlObjSKZ;

     MtrxMultMatrx(arrJ,NUmAng, 3, arrK,3, arrT3) ;
     MtrxMultMatrxTransp(arrT3,NUmAng, 3, arrJ,NUmAng, arrJKJT_Cur) ;

}
//-----------------------------------------------------------------

//-----------------------------------------------------------------
void QContObj6::calcMtrxAngs_CKCT(double *arrX0,double *arrS_KGSK,double *arrCKCT_Cur)
{
    //0.
    double arr_g[2 * NUmAng ] = {0.}, arr_gBig[2 * QUANT_GDG_PRMS ] = {0.};
    calc_dV_po_dX_Var2(arrX0,marrEilerCntrKP_Zv, arrS_KGSK,arr_gBig);
    for (int i =0; i < 2; ++i)
    {
        for (int j =0; j < NUmAng; ++j)
        {
           arr_g[i * NUmAng + j] =  arr_gBig[ i * QUANT_GDG_PRMS + 3 + j];
        }
    }

  /*  // 1. пересчет замера в ПСК
       double arrS_PSK_Zv[3] = {0.};
       double matrPereh_PSK_V_KGSK[9] = {0.};
       calcMatr_PSK_v_KGSK(marrEilerCntrKP_Zv, matrPereh_PSK_V_KGSK) ;
       MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3,marrS_KGSK,1, arrS_PSK_Zv) ;
    ///

       //2. пересчет замера в ПСК-РЛК
         double arrS_PSK_Gdg_Zv[3] = {0.};
         MtrxMinusMatrx(arrS_PSK_Zv, arrX0,1, 3, arrS_PSK_Gdg_Zv);

    // 3. пересчет замера в АСПК
      double arrS_ASPK_Zv[3] = {0.}, matrPereh_ASK_V_PSK_Gdg[9] = {0.};
    calcMtrx3_ASPK_v_PSK(&(arrX0[3]),matrPereh_ASK_V_PSK_Gdg) ;
    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK_Gdg,3, 3, arrS_PSK_Gdg_Zv,1, arrS_ASPK_Zv) ;
    ///



    // 5.вычисление dB_po_dS
    double arr_dB_po_dS[9] = {0.};//, arr_dB_po_dS_MultC[9] = {0.};

    calc_dSpherical_po_Rectangular_AnglesOnly(arrS_ASPK_Zv, arr_dB_po_dS);
   //calc_dSpherical_po_Rectangular_AnglesOnly_MultCoeff(COeff, arrS_ASPK_Zv, arr_dB_po_dS_MultC);
    ///




    // 6. вычисление первых трех столбцов матрицы С и внесение их на место
    double arrT1[6] = {0.},arrT2[NUmAng * 3] = {0.}
                          ,arrT3[NUmAng * 3] = {0.};

      MtrxMultMatrxTransp(arr_dB_po_dS,2, 3, matrPereh_ASK_V_PSK_Gdg,3, arrT1) ;

      MtrxTranspMultMatrx(arr_g,2, NUmAng, arrT1,3, arrT2) ;

    // 7.
      double arrMiddleC[9] = {0.};
      // 8 вычисление dM/dQ*Skgsk и вставка на место
      double arr_dM_po_dQ[9] = {0.},  arrCol0[3] ={0.};
      calc_dM_psk_to_kgsk_po_dQ(marrEilerCntrKP_Zv,arr_dM_po_dQ);
      MtrxTranspMultMatrx(arr_dM_po_dQ,3, 3, arrS_KGSK,1, arrCol0) ;
      changeCol(arrMiddleC, 3,  3,0 ,arrCol0);

      // 9 вычисление dM/dPsi*Skgsk и вставка на место
      double arr_dM_po_dPsi[9] = {0.},  arrCol1[3] ={0.};
      calc_dM_psk_to_kgsk_po_dPsi(marrEilerCntrKP_Zv,arr_dM_po_dPsi);
      MtrxTranspMultMatrx(arr_dM_po_dPsi,3, 3, arrS_KGSK,1, arrCol1) ;
      changeCol(arrMiddleC, 3,  3,1 ,arrCol1);

      // 9 вычисление dM/dTet*Skgsk и вставка на место
      double arr_dM_po_dTet[9] = {0.},  arrCol2[3] ={0.};
      calc_dM_psk_to_kgsk_po_dTet(marrEilerCntrKP_Zv,arr_dM_po_dTet);
      MtrxTranspMultMatrx(arr_dM_po_dTet,3, 3, arrS_KGSK,1, arrCol2) ;
      changeCol(arrMiddleC, 3,  3,2 ,arrCol2);

      // 10 вычисление матрицы C
      double arrC[NUmAng * 3] = {0.};
      MtrxMultMatrx(arrT2,NUmAng, 3, arrMiddleC,3, arrC) ;
*/
    double arr_dViu_po_dMu[6] = {0.};
    calc_dVGadget_po_dMu(arrX0, marrEilerCntrKP_Zv, arrS_KGSK, arr_dViu_po_dMu);

  //  calc_dVGadget_po_dMu_Razn(arrX0, marrEilerCntrKP_Zv
     //                                     , arrS_KGSK, arr_dViu_po_dMu);

     double arrC[NUmAng * 3] = {0.};
     MtrxTranspMultMatrx(arr_g,2, NUmAng, arr_dViu_po_dMu,3, arrC) ;
      // 11. ывчисление коррел матрицы
      double arrK[9] = {0.};
      for(int i = 0; i < 3; ++i)
      {
         arrK[3 * i + i] =  marrDispEilerCntrKP[i];
      }

      // 12 вычисление CKCT
      double arrT4[NUmAng * 3] = {0.};
      MtrxMultMatrx(arrC,NUmAng, 3, arrK,3, arrT4) ;

      MtrxMultMatrxTransp(arrT4,NUmAng, 3, arrC,NUmAng, arrCKCT_Cur) ;

      // проверка
    //  double arr_dViu_po_dMu1[6] = {0.};
     // calc_dVGadget_po_dMu_Razn(arrX0, marrEilerCntrKP_Zv
                                          //  , arrS_KGSK, arr_dViu_po_dMu1);
     // double arrCKCT_Cur1[9] = {0.};
     // MtrxTranspMultMatrx(arr_g,2, NUmAng, arr_dViu_po_dMu1,3, arrCKCT_Cur1) ;
     // int uyu = 0;

}
//-----------------------------------------------------------------
void QContObj6::calcMtrxAngs_g2TKg2(double *arrX0,double *arrS_KGSK,double *arrg2TKg2_Cur)
{
    //0.
    double arr_g[2 * NUmAng ] = {0.}, arr_gBig[2 * QUANT_GDG_PRMS ] = {0.};
    calc_dV_po_dX_Var2(arrX0, marrEilerCntrKP_Zv,arrS_KGSK,arr_gBig);
    for (int i =0; i < 2; ++i)
    {
        for (int j =0; j < NUmAng; ++j)
        {
           arr_g[i * NUmAng + j] =  arr_gBig[ i * QUANT_GDG_PRMS + 3 + j];
        }
    }


    // 2. формирование матрицы K
    double arrK[4] = {0.};
    arrK[0] = arrK[3] = mAngleSigIU * mAngleSigIU;

    // 3
    double arrT0[NUmAng * 2] = {0.};
    MtrxTranspMultMatrx(arr_g,2, NUmAng, arrK,2, arrT0) ;
    MtrxMultMatrx( arrT0,NUmAng, 2, arr_g,NUmAng, arrg2TKg2_Cur) ;

}

//--------------------------------

//-----------------------------------------------------------------

bool QContObj6::doOneIteration_Var2(double *arrX0,QContObj6 *pContObjectArr
            ,const int QuantObj,double *parrVZv_IU
           , double *arrDel, double *arrMtrx_dFGr_po_dX_Inv)
{

    double arrF [QUANT_GDG_PRMS]= {0.};
    double  arrMtrx_dFGr_po_dX [QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.};
    calc_FGr_and_dFGr_po_dX_Var2(arrX0,pContObjectArr, QuantObj, parrVZv_IU
                        , arrF, arrMtrx_dFGr_po_dX);


   // double valScal = 1000.;
   // MatrxMultScalar(arrMtrx_dF_po_dX, QUANT_GDG_PRMS, QUANT_GDG_PRMS, valScal,arrMtrx_dF_po_dX);
    bool brez =   InverseMtrx(arrMtrx_dFGr_po_dX, QUANT_GDG_PRMS, arrMtrx_dFGr_po_dX_Inv);
    // проверка
   // double aeeT[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.};
    // MtrxMultMatrx( arrMtrx_dF_po_dX,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrMtrx_dF_po_dX_Inv,QUANT_GDG_PRMS, aeeT) ;
    ///
    if (!brez)
    {
        return false;
    }


   double arrCoeff[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.}, arrT[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.};
   fillE(arrCoeff, QUANT_GDG_PRMS);
   for (int  i=0; i < 3; ++i)
   {
     arrCoeff[ QUANT_GDG_PRMS * i + i] = COeff;
   }
    MtrxMultMatrx( arrCoeff,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrMtrx_dFGr_po_dX_Inv,QUANT_GDG_PRMS, arrT) ;

    MtrxMultMatrx(arrT,QUANT_GDG_PRMS, QUANT_GDG_PRMS, arrF,1, arrDel) ;


   // memcpy(arrMtrx_dFGr_po_dX_Inv, arrT, QUANT_GDG_PRMS * QUANT_GDG_PRMS * sizeof(double));
    return true;

}

//---------------------------------------------
void  QContObj6::calc_FGr_and_dFGr_po_dX_Var2(double *arrX0,QContObj6 *pContObjectArr
                    ,const int QuantObj,double *parrVZv_RLK
                    ,double *parrF,double *parrMtrx_dF_po_dX)
{
    memset(parrF, 0, QUANT_GDG_PRMS * sizeof(double));
    memset(parrMtrx_dF_po_dX, 0, QUANT_GDG_PRMS *QUANT_GDG_PRMS * sizeof(double));

    for (int i = 0; i < QuantObj; ++i)
    {
        double parrf [QUANT_GDG_PRMS] = {0.};
        double parr_df_po_dX[QUANT_GDG_PRMS * QUANT_GDG_PRMS] = {0.};
            pContObjectArr[i].calc_fi_and_dfi_po_dx_Var2(arrX0, &(parrVZv_RLK[2 * i])
                    , parrf,parr_df_po_dX);
            MtrxSumMatrx(parrF, parrf,1, QUANT_GDG_PRMS, parrF) ;
            MtrxSumMatrx(parrMtrx_dF_po_dX, parr_df_po_dX,QUANT_GDG_PRMS, QUANT_GDG_PRMS, parrMtrx_dF_po_dX) ;

    }
}
//-------------------------------------------
//-------------------------------------------------
void  QContObj6::calc_fi_and_dfi_po_dx_Var2(double *arrX0, double *arrVZv_RLK
          ,double *arrfi,double *arrMtrx_dfi_po_dX)
{
  //1 вычисление вектор-функции fi
    // 1.1 векторт невязок координат ИУ f-Viu
    double arrDelV_Gdg  [2] = {0.};
    double arrf [2 ] = {0.};
    calc_VGadget(marrEilerCntrKP_Zv, arrX0, marrSZv_KGSK, arrf);

    MtrxMinusMatrx(arrf, arrVZv_RLK,1, 2,  arrDelV_Gdg);

   ///

   //1.2 матрица частных проихводных Viu по X. обозначается g
    double arr_g[ 2 * QUANT_GDG_PRMS ];
    calc_dV_po_dX_Var2( arrX0, marrEilerCntrKP_Zv,marrSZv_KGSK, arr_g);
    ///

    //1.3 вычисление arrf
    MtrxTranspMultMatrx( arr_g, 2,QUANT_GDG_PRMS, arrDelV_Gdg,1, arrfi) ;

    ///

//2 Вычисление матрицы dfi_po_dx
    // 2.1
    double arr_g1 [QUANT_GDG_PRMS *2 ] = {0.};
    memcpy(arr_g1, arr_g,  QUANT_GDG_PRMS *2 * sizeof( double));
    MtrxTranspMultMatrx( arr_g, 2,QUANT_GDG_PRMS
                         , arr_g1,QUANT_GDG_PRMS, arrMtrx_dfi_po_dX) ;
    ///
/*
  for(int j =0; j< QUANT_GDG_PRMS; ++j)
    {
         double *arr_dg_po_dXj = new  double[QUANT_GDG_PRMS *(2 + BDalnomer) ] ;
        calcMatr_dg_po_dXj_Var2(arrX0,  BDalnomer ,marrSZv_KGSK,j,  arr_dg_po_dXj);
         double *arrColumn= new  double[QUANT_GDG_PRMS];

        MtrxTranspMultMatrx( arr_dg_po_dXj, (2 + BDalnomer),QUANT_GDG_PRMS ,arrDelV_Gdg,1, arrColumn) ;
        for(int i =0; i < QUANT_GDG_PRMS;++i)
        {
          arrMtrx_dfi_po_dX[i * QUANT_GDG_PRMS + j] += arrColumn[i];
        }

        delete []arr_dg_po_dXj;
        delete []arrColumn;
    }
    */


}




//------------------------------------------------------------
//Вычисление матрицы частных производных вектор-функции целеуказаний
// по вектору параметров позиционирования РЛК
//INPUT:
//BDalnomer- признак дальномера
//arrXCur[7] - вектор параметров позиционирования
//arrVZv_RLK[3] - положение цели (КО) в АСфСК
//OUTPUT:
//  arr_g[ (BDalnomer +2)*5]
// вектор параллакса - относитльный в ПСК-ИУ

void  QContObj6::calc_dV_po_dX_Var2(double *arrX0, double *arrEilerCntrKP,double *arrS_KGSK,double *arr_g)
{
// 1. пересчет замера в ПСК
   double arrS_PSK_Zv[3] = {0.};
   double matrPereh_PSK_V_KGSK[9] = {0.};
   calcMatr_PSK_v_KGSK(arrEilerCntrKP, matrPereh_PSK_V_KGSK) ;
   MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3,arrS_KGSK,1, arrS_PSK_Zv) ;
///

   //2. пересчет замера в ПСК-РЛК
     double arrS_PSK_Gdg_Zv[3] = {0.};
     MtrxMinusMatrx(arrS_PSK_Zv, arrX0,1, 3, arrS_PSK_Gdg_Zv);

// 3. пересчет замера в АСПК
  double arrS_ASPK_Zv[3] = {0.}, matrPereh_ASK_V_PSK_Gdg[9] = {0.};
calcMtrx3_ASPK_v_PSK(&(arrX0[3]),matrPereh_ASK_V_PSK_Gdg) ;
MtrxTranspMultMatrx(matrPereh_ASK_V_PSK_Gdg,3, 3, arrS_PSK_Gdg_Zv,1, arrS_ASPK_Zv) ;
///



// 5.вычисление dB_po_dS
double arr_dB_po_dS[9] = {0.}, arr_dB_po_dS_MultC[9] = {0.};

calc_dSpherical_po_Rectangular_AnglesOnly(arrS_ASPK_Zv, arr_dB_po_dS);
calc_dSpherical_po_Rectangular_AnglesOnly_MultCoeff(COeff, arrS_ASPK_Zv, arr_dB_po_dS_MultC);
///




// 6. вычисление первых трех столбцов матрицы С и внесение их на место
  double arrBlok1[9];
  //MtrxMultMatrxTransp(arr_dB_po_dS_MultC,(2 + BDalnomer), 3, matrPereh_ASK_V_PSK_Gdg,3, arrBlok1) ;
  MtrxMultMatrxTransp(arr_dB_po_dS_MultC,2, 3, matrPereh_ASK_V_PSK_Gdg,3, arrBlok1) ;
for (int i =0; i < 2; ++i)
    for (int j =0; j < 3; ++j)
    {
        arr_g[ i * QUANT_GDG_PRMS + j] = -arrBlok1[i * 3 + j];
    }
///

// 7.вычисление 4-го столбца матрицы С и внесение его на место
  double arr_dM_Gdg_po_dBet[9] = {0.};
calc_dMtrx3_ASPK_v_PSK_po_dBet(&(arrX0[3]),arr_dM_Gdg_po_dBet);

 double arrTemp0[3] = {0.},  arrAlf4[3] = {0.};
MtrxTranspMultMatrx(arr_dM_Gdg_po_dBet,3, 3, arrS_PSK_Gdg_Zv,1, arrTemp0) ;

MtrxMultMatrx( arr_dB_po_dS, 2, 3, arrTemp0,1, arrAlf4) ;

for (int i =0; i < 2; ++i)
{
   arr_g[ i * QUANT_GDG_PRMS +3] =arrAlf4[i];
}


///
// 7.вычисление 5-го столбца матрицы С и внесение его на место
  double arr_dM_Gdg_po_dEps[9] = {0.};
calc_dMtrx3_ASPK_v_PSK_po_dEps(&(arrX0[3]),arr_dM_Gdg_po_dEps);
  double arrTemp1[3] = {0.},  arrAlf5[3] = {0.};
MtrxTranspMultMatrx(arr_dM_Gdg_po_dEps,3, 3, arrS_PSK_Gdg_Zv,1, arrTemp1) ;

MtrxMultMatrx( arr_dB_po_dS, 2, 3, arrTemp1,1, arrAlf5) ;

for (int i =0; i < 2; ++i)
{
   arr_g[ i * QUANT_GDG_PRMS + 4] =arrAlf5[i];
}


// 7.вычисление 6-го столбца матрицы С и внесение его на место
  double arr_dM_Gdg_po_dAlf[9] = {0.};
calc_dMtrx3_ASPK_v_PSK_po_dAlf(&(arrX0[3]),arr_dM_Gdg_po_dAlf);
  double arrTemp2[3] = {0.},  arrAlf6[3] = {0.};
MtrxTranspMultMatrx(arr_dM_Gdg_po_dAlf,3, 3, arrS_PSK_Gdg_Zv,1, arrTemp2) ;

MtrxMultMatrx( arr_dB_po_dS, 2, 3, arrTemp2,1, arrAlf6) ;

for (int i =0; i < 2; ++i)
{
   arr_g[ i * QUANT_GDG_PRMS + 5] =arrAlf6[i];
}

}
//------------------------------------------------------------

/*
void QContObj6::calcMatr_dg_po_dXj_Var2(double *arrX0, const bool BDalnomer ,double *arrSZv_KGSK
                   ,const int NUmj,  double *arr_dg_po_dXj)
{


    double arrX1[NUmAng] = {0.};

    memcpy(arrX1, arrX0, NUmAng * sizeof(double));

    double delta = 0.0001;

    arrX1[3 + NUmj]+= delta;
    double arr_g  [NUmAng *2 ] = {0.};
    double arr_g1 [NUmAng *2 ] = {0.};
    double arrT   [NUmAng *2 ] = {0.};

    calc_dVGdg_po_dX_Var2( arrX0,arrSZv_KGSK, arr_g);
    calc_dVGdg_po_dX_Var2(arrX1, arrSZv_KGSK, arr_g1);

    MtrxMinusMatrx(arr_g1, arr_g,NUmAng, 2, arrT);
    MatrxMultScalar(arrT, NUmAng, 2, 1./delta,arr_dg_po_dXj);


}
*/

//----------------------------------------------------
//----------------------------------------------------
//-----  ЗВЕЗДА  -----------------------------------------------
//----------------------------------------------------
//----------------------------------------------------


//-----------------------------------------------------------------------
// вычисление оценки вектора позиционирования РЛК и ИУ в ПСК-ИУ по Варианту 1
// INPUT:
//pContObjectArr - массив КО
//LEnArr - кол-во КО
// BNoise - признак имитации случайных шумов в измерениях
//OUTPUT:
//arrXEst[5] - оценка вектора позиционирования РЛК в ИУ
//arrMtrxCor[25] - коррел матрица ошибок оценивания
double QContObj6::imitate_and_estimateStarProcessing(QContObj6 *pContObjectArr,const int QuantObj
       , const bool BNoise, double *arrXEst, double *arrMtrxCor)
{
    // векторы измерений РЛК
 double *parrVZv_IU = new double [QuantObj * 2];

 // векторы измерений ИУ  они не нужны
 double arrVZv_RLK [3];

 // векторы истинного положения цели в ИУ
 //double *parrVTrue_IU = new double [QuantObj * (BDalnomer +2) ];

 // имитация измерений РЛК и ИУ
  for(int i =0; i < QuantObj; ++i)
  {
    pContObjectArr[i].imitateStarMeasure(BNoise,&(parrVZv_IU[i * 2])) ;
  }

  double rez = estimateAngs_Var_2(pContObjectArr,QuantObj
                             ,parrVZv_IU,arrXEst, arrMtrxCor);

 delete []parrVZv_IU;


  return rez;
}


//----------------------------------------------------
// имитация измерений РЛК и ИУ
void QContObj6::imitateStarMeasure(const bool BNoise, double *arrVZvIU)
{

    // 1. истинное положение в АСфСАК

    memcpy(marrSZv_KGSK, marrS_KGSK, 3 * sizeof(double));


    // 1. истинное положение в ИУСфСК
    double arrVTrue_IU[3] = {0.};
    recalcPositionFromKGSK_to_SphericalSK(marrEilerCntrKP,marrS_KGSK, marrXTrue_IU_PSK,arrVTrue_IU);

    // 2. имитация замера в ИУСфСК
    arrVZvIU[0] =  arrVTrue_IU[0] ;
    arrVZvIU[1] =  arrVTrue_IU[2] ;
    if ( BNoise)
    {
        arrVZvIU[0] +=  getGauss(0, mAngleSigIU);
        arrVZvIU[1] +=  getGauss(0, mAngleSigIU);
    }
    // 3. имитация ошибок СИНС
    memcpy(marrEilerCntrKP_Zv,marrEilerCntrKP, 3 * sizeof(double));
    if(BNoise)
    {
    for (int i =0; i < 3; ++i)
    {
      marrEilerCntrKP_Zv[i] =  marrEilerCntrKP[i] + getGauss(0, sqrt(marrDispEilerCntrKP[i]));
    }
    }

 }
//------------------------------------------------------
// раскидывания  вектора переменныъх задачи МНК в векторы позиционирования
// РЛК и ИУ
//INPUT:
// arrXInp[QUANT_GDG_PRMS]- вектор неизвестных параметров задачи МНК
//    последовательность переменных следующая
//     0. XRelative
//     1. YRelative
//     2. ZRelative
//     3. BetRelative
//     4. Eps RLK
//     5. Eps IU
//     6. Alf RLK
//     7. Alf IU
// arrXInp_IU[QUANT_GDG_PRMS] - вектор исходных параметров позиционирования ИУ
// OUTPUT:
// arrXOut_RLK[QUANT_GDG_PRMS]- вектор позиционирования РЛК обновленный
// arrXOut_IU[QUANT_GDG_PRMS] - вектор позиционирования ИУ обновленный
void QContObj6::recombine(double *arrXCombinedInp,double *arrXInp_IU, double *arrXOut_RLK,double *arrXOut_IU)
{
    memcpy(arrXOut_IU, arrXInp_IU, QUANT_GDG_PRMS* sizeof(double));
    for (int j =0; j < QUANT_GDG_PRMS; ++j )
    {
    if ( j <3)
    {
       arrXOut_RLK[j] =   arrXCombinedInp[j] + arrXInp_IU[j];

    }
    else
    {
        switch(j)
        {
        case 3: // Bet RLK
            arrXOut_RLK[3] = arrXCombinedInp[j];
            break;
        case 4: // Eps RLK
            arrXOut_RLK[4] = arrXCombinedInp[j];

            break;
        case 5: // Eps IU
            arrXOut_IU[4] = arrXCombinedInp[j];
            break;
        case 6: // Alf RLK
            arrXOut_RLK[5] = arrXCombinedInp[j];
            break;
        case 7: // Alf IU
            arrXOut_IU[5] = arrXCombinedInp[j];
            break;
        default:
            break;
        }

    }
    }
}


//------------------------------------------------------
// формирование вектора
// РЛК и ИУ
//INPUT:
// arrXInp_RLK[QUANT_GDG_PRMS] - вектор параметров позиционирования РЛК,
//        3 координаты вектора параалакса + 3 угла Эйлера
// arrXInp_IU[QUANT_GDG_PRMS] - вектор параметров позиционирования ИУ,
//        3 координаты вектора паралакса + 3 угла Эйлера
// OUTPUT:
// QUANT_GDG_PRMS[QUANT_GDG_PRMS] - вектор неизвестных параметров задачи МНК
// последовательность переменных следующая
// 0. XRelative
// 1. YRelative
// 2. ZRelative
// 3. BetRelative
// 4. Eps RLK
// 5. Eps IU
// 6. Alf RLK
// 7. Alf IU

void QContObj6::combine(double *arrXInp_RLK,double *arrXInp_IU,double *arrXOut)
{
    for (int j =0; j < QUANT_GDG_PRMS; ++j )
    {
    if ( j <3)
    {
        arrXOut[j] = arrXInp_RLK[j] - arrXInp_IU[j];

    }
    else
    {
        switch(j)
        {
        case 3:
            arrXOut[j] = arrXInp_RLK[3];
            break;
        case 4: // Eps RLK
            arrXOut[j] = arrXInp_RLK[4];

            break;
        case 5: // Eps IU
        arrXOut[j] = arrXInp_IU[4];
            break;
        case 6: // Alf RLK
            arrXOut[j] = arrXInp_RLK[5];
            break;
        case 7: // Alf IU
            arrXOut[j] = arrXInp_IU[5];
            break;
        default:
            break;
        }

    }
    }
}
//--------------------------------------
//-------------------------------------------
// вычисление разностной производной функции f по переменной
// с номером j
//функция f зависит от 7 переменных
//от переменных вектора позиционирования РЛК и ИУ  в ПСК - arrX

//OUTPUT:
//arr_df_po_dViu[2 + BDalnomer] - вектор частных производных ыектор функции Viu
//по переменной с номером NUmj (нумерация с 0)
void QContObj6::calc_dVGadget_po_dMu_Razn(double *arrX, double *arrMu
                         , double *arrS_KGSK, double *arr_dViu_po_dMu)
{
  for(int j =0; j < 3; ++j)
  {
      double arrT  [2] = {0.};

      calc_dVGadget_po_dMuj_Razn(arrX,arrMu, arrS_KGSK,j, arrT ) ;
      changeCol(arr_dViu_po_dMu, 2,  3 ,j ,arrT);


  }

}
//------------------------------
//-------------------------------------------
// вычисление разностной производной функции f по переменной
// с номером j
//функция f зависит от 7 переменных
//от переменных вектора позиционирования РЛК и ИУ  в ПСК - arrX

//OUTPUT:
//arr_df_po_dViu[2 + BDalnomer] - вектор частных производных ыектор функции Viu
//по переменной с номером NUmj (нумерация с 0)
void QContObj6::calc_dVGadget_po_dMuj_Razn(double *arrX, double *arrMu, double *arrS_KGSK,const int NUmj, double *arr_dViu_po_dMuj)
{
    double arrMu1[3] = {0.};

    memcpy(arrMu1, arrMu, 3 * sizeof(double));

    double delta = 0.0005;

    arrMu1[NUmj]+= delta;
    double arrViu[2] = {0.},arrViu1[2] = {0.}, arrT[2] = {0.},arrDelta[2] = {0.};
    calc_VGadget(arrMu, arrX, arrS_KGSK, arrViu);
    calc_VGadget(arrMu1, arrX, arrS_KGSK, arrViu1);


    MtrxMinusMatrx(arrViu1, arrViu,1, 2, arrT);
    MatrxMultScalar(arrT, 1, 2, 1./delta,arr_dViu_po_dMuj);


}
//------------------------------
//-----------------------------------------------------------------
void QContObj6::calc_dVGadget_po_dMu(double *arrX, double *arrMu
                                           , double *arrS_KGSK, double *arr_dViu_po_dMu)
{

    // 1. пересчет замера в ПСК
       double arrS_PSK_Zv[3] = {0.};
       double matrPereh_PSK_V_KGSK[9] = {0.};
       calcMatr_PSK_v_KGSK(arrMu, matrPereh_PSK_V_KGSK) ;
       MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3,arrS_KGSK,1, arrS_PSK_Zv) ;
    ///

       //2. пересчет замера в ПСК-РЛК
         double arrS_PSK_Gdg_Zv[3] = {0.};
         MtrxMinusMatrx(arrS_PSK_Zv, arrX,1, 3, arrS_PSK_Gdg_Zv);

    // 3. пересчет замера в АСПК
      double arrS_ASPK_Zv[3] = {0.}, matrPereh_ASK_V_PSK_Gdg[9] = {0.};
    calcMtrx3_ASPK_v_PSK(&(arrX[3]),matrPereh_ASK_V_PSK_Gdg) ;
    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK_Gdg,3, 3, arrS_PSK_Gdg_Zv,1, arrS_ASPK_Zv) ;
    ///



    // 5.вычисление dB_po_dS
    double arr_dB_po_dS[9] = {0.};//, arr_dB_po_dS_MultC[9] = {0.};

    calc_dSpherical_po_Rectangular_AnglesOnly(arrS_ASPK_Zv, arr_dB_po_dS);
    ///




    // 6. вычисление первых трех столбцов матрицы С и внесение их на место
    double arrT1[6] = {0.};

      MtrxMultMatrxTransp(arr_dB_po_dS,2, 3, matrPereh_ASK_V_PSK_Gdg,3, arrT1) ;



    // 7.
      double arrMiddleC[9] = {0.};
      // 8 вычисление dM/dQ*Skgsk и вставка на место
      double arr_dM_po_dQ[9] = {0.},  arrCol0[3] ={0.};
      calc_dM_psk_to_kgsk_po_dQ(arrMu,arr_dM_po_dQ);
      MtrxTranspMultMatrx(arr_dM_po_dQ,3, 3, arrS_KGSK,1, arrCol0) ;
      changeCol(arrMiddleC, 3,  3,0 ,arrCol0);

      // 9 вычисление dM/dPsi*Skgsk и вставка на место
      double arr_dM_po_dPsi[9] = {0.},  arrCol1[3] ={0.};
      calc_dM_psk_to_kgsk_po_dPsi(arrMu,arr_dM_po_dPsi);
      MtrxTranspMultMatrx(arr_dM_po_dPsi,3, 3, arrS_KGSK,1, arrCol1) ;
      changeCol(arrMiddleC, 3,  3,1 ,arrCol1);

      // 9 вычисление dM/dTet*Skgsk и вставка на место
      double arr_dM_po_dTet[9] = {0.},  arrCol2[3] ={0.};
      calc_dM_psk_to_kgsk_po_dTet(arrMu,arr_dM_po_dTet);
      MtrxTranspMultMatrx(arr_dM_po_dTet,3, 3, arrS_KGSK,1, arrCol2) ;
      changeCol(arrMiddleC, 3,  3,2 ,arrCol2);

      // 10 вычисление матрицы arr_dViu_po_dMu
      double arrC[NUmAng * 3] = {0.};
      MtrxMultMatrx(arrT1,NUmAng, 3, arrMiddleC,3, arr_dViu_po_dMu) ;


}

//-------------------------------------------------------
// INPUT:
// arrNu[3] - истинные углы РЛК
//arrMu[3]- истинные углы ИУ
// arrMuEst[3]- оценка углов ИУ
//OUTPUT:
//arrEilersRlkIdeal[3] - углы РЛК
void QContObj6::calc_IdealRelativeEilers_Var1(double *arrNu,double *arrMu
               ,double *arrMuEst ,double *arrEilersRlkIdeal)
{
    //1.
    double mtrxM_ot_Nu[9];
    calcMtrx3_ASPK_v_PSK(arrNu,mtrxM_ot_Nu);

    //2.
    double mtrxL_ot_Mu[9];
    calcMtrx3_ASPK_v_PSK(arrMu,mtrxL_ot_Mu);

    //3.
    double mtrxL_ot_EstMu[9];
    calcMtrx3_ASPK_v_PSK(arrMuEst,mtrxL_ot_EstMu);

    // 4.
    double mtrxM_ot_EstNu[9] = {0.}, arrTemp[9] = {0.};

    MtrxMultMatrxTransp(mtrxL_ot_EstMu,3, 3, mtrxL_ot_Mu,3, arrTemp) ;

    MtrxMultMatrx( arrTemp,3, 3, mtrxM_ot_Nu,3, mtrxM_ot_EstNu) ;

    // 5.
     double eps1 = asin(mtrxM_ot_EstNu[7]);
     double alf1 = asin(-mtrxM_ot_EstNu[6]/ cos(eps1));
     double bet1 = asin(mtrxM_ot_EstNu[1]/ cosl(eps1));

    arrEilersRlkIdeal[0]= (double)bet1;
    arrEilersRlkIdeal[1]= (double)eps1;
    arrEilersRlkIdeal[2]= (double)alf1;

}
//-------------------------------------------------------



