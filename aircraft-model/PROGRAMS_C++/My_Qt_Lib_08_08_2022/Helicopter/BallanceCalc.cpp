#include "BallanceCalc.h"
#include "Helic.h"
#include "Environment.h"
#include <string.h>

#include <math.h>
#include "MatrixProccess.h"
#include "PartHelicTraj.h"





TBallanceCalc::TBallanceCalc()
{

}
//----------------------------------------------------------------------------------
// нахождение балансировочных параметров
// для равномерного прямолинейного движения
// балнсировочные параметры представляют из себя вектор arrXRez
// arrXRez[0] - угол крена
// arrXRez[1] - угол тангажа
// arrXRez[2] - цикл шаг по тангажу Kappa
// arrXRez[3] - общий шаг Fi
// arrXRez[4] - цикл шаг по крену Etta
// arrXRez[5] - дифф шаг DelFi
// INPUT:
// Helic  -вертолет
// Environment - атмосфера
// VAlY - высота
// VAlVx - Скорость
// VAlPsi - угол рысканя
// arrX0[6] - гачальное значение вектора
//OUTPUT:
// arrXRez[6]  - искомый вектор
// возвращает true если решение найдено
// false в противном случае
bool TBallanceCalc::calcBallParamsForSteadyLineMoving(const THelic Helic,const TEnvironment Environment
                                                      ,const long double VAlY,const long double VAlVx,const long double VAlPsi
                                                      ,long double *arrX0, long double *arrXRez )
{
    long double arrXCur[6] = {0.}
                ,arrFGr[6] = {0.}  // правая часть сиситемы уравнений Fgr(x) = 0
                ,arrJacFgr[36] = {0.} // якобиан
                ,arrJacFGrInv[36] = {0.}
                ,arrDel[6] = {0.} ; // вектор невязки

    // вектор приращений аргументов для вычислкеения якобиана
    long double arr_dx[6] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};

    memcpy(arrXCur, arrX0 , 6 * sizeof(long double));
    for (int i = 0; i < 100; ++i)
    {
        calcFgr_and_JacFgr_ForSteadyLineMoving( Helic, Environment
                ,VAlY, VAlVx,  VAlPsi,arrXCur, arrFGr, arrJacFgr);
        if (!InverseMtrx6(arrJacFgr, arrJacFGrInv))
        {
            return false;
        }
        MtrxMultMatrx(arrJacFGrInv,6, 6, arrFGr,1, arrDel) ;
        MtrxMinusMatrx(arrXCur, arrDel,6, 1, arrXRez);
        long double valDel = NormVect(arrDel, 6);
        if(fabsl(valDel) < 0.00001)
        {
          return true;
        }

        memcpy(arrXCur, arrXRez , 6 * sizeof(long double));
    }

    return false;
}
//----------------------------------------------------------------------------------
//вычисление вектор функции правой части уравненипя балансировки
//и ее якобиана
// для равномерного прямолинейного движения
// балнсировочные параметры представляют из себя вектор arrXRez
// arrXCur[0] - угол крена
// arrXCur[1] - угол тангажа
// arrXCur[2] - цикл шаг по тангажу Kappa
// arrXCur[3] - общий шаг Fi
// arrXCur[4] - цикл шаг по крену Etta
// arrXCur[5] - дифф шаг DelFi
// INPUT:
// Helic  -вертолет
// Environment - атмосфера
// VAlY - высота
// VAlVx - Скорость
// VAlPsi - угол рыскания
// arrXCur[6]  - вектор переменных
//OUTPUT:
// arrFGr[6]  - вектор функция
// arrJacFgr [36] - ee якобиан
void TBallanceCalc::calcFgr_and_JacFgr_ForSteadyLineMoving(const THelic Helic,const TEnvironment Environment
        , long double VAlY, long double VAlVx, long double VAlPsi
        ,long double *arrX, long double *arrFGr, long double *arrJacFgr)
{
    long double arrJacFgrTransp[36] = {0.} // транспогнированный якобиан
                ,arrDel[6] = {0.}
                ,arrTemp[6] = {0.}             ; // вектор невязки

    // вектор приращений аргументов для вычислкеения якобиана
    long double arr_dx[6] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};

    long double arrFGrCur[6] = {0.}, arrXCur[6] = {0.};

    calcFgr_ForSteadyLineMoving( Helic, Environment
            ,  VAlY,  VAlVx, VAlPsi,arrX, arrFGr);

    for (int i = 0; i < 6; ++i)
    {
       memcpy(arrXCur, arrX , 6 * sizeof(long double));
       arrXCur[i] += arr_dx[i] ;
       calcFgr_ForSteadyLineMoving( Helic, Environment
               ,  VAlY,  VAlVx, VAlPsi,arrXCur, arrFGrCur);
       MtrxMinusMatrx(arrFGrCur,arrFGr,6, 1, arrTemp);
       MatrxDivideScalar(arrTemp, 1, 6, arr_dx[i], &(arrJacFgrTransp[ 6 * i]));
    }

    MatrTransp(arrJacFgrTransp, 6, 6, arrJacFgr);

}
//----------------------------------------------------------------------------------
//вычисление вектор функции правой части уравненипя балансировки
// для равномерного прямолинейного движения
// балнсировочные параметры представляют из себя вектор arrXRez
// arrXCur[0] - угол крена
// arrXCur[1] - угол тангажа
// arrXCur[2] - цикл шаг по тангажу Kappa
// arrXCur[3] - общий шаг Fi
// arrXCur[4] - цикл шаг по крену Etta
// arrXCur[5] - дифф шаг DelFi
// INPUT:
// Helic  -вертолет
// Environment - атмосфера
// VAlY - высота
// VAlVx - Скорость
// VAlPsi - угол рыскания
// arrXCur[6]  - вектор переменных
//OUTPUT:
// arrFGr[6]  - вектор функция
void TBallanceCalc::calcFgr_ForSteadyLineMoving(const THelic Helic,const TEnvironment Environment
        , long double VAlY, long double VAlVx, long double VAlPsi
        ,long double *arrXCur, long double *arrFGr)
{
  TPartHelicTraj PartHelicTraj( Helic, Environment, 0.,0.);
  memset(PartHelicTraj.marrPhaseVect, 0, 12 * sizeof(long double));
  PartHelicTraj.marrPhaseVect[1] = VAlY;
  PartHelicTraj.marrPhaseVect[3] = VAlVx;
  PartHelicTraj.marrPhaseVect[11]= VAlPsi;

  PartHelicTraj.marrPhaseVect[9]  = arrXCur[0];
  PartHelicTraj.marrPhaseVect[10] = arrXCur[1];
  PartHelicTraj.doAirDensity();

  long double    arrf[12] = {0.};
  long double arrSvSK_Force[3] = {0.}, arrSvSK_Moment[3] = {0.}
          ,  arrSvSK_ShaftForce0[3] = {0.}, arrSvSK_ShaftMoment0[3] = {0.}
          , arrSvSK_AirForce0 [3] = {0.}, arrSvSK_AirMoment0 [3] = {0.} ;


  PartHelicTraj.calcRightMemberOfDifEqSystem(arrXCur[3], arrXCur[2]
                                                  , arrXCur[4], arrXCur[5]
                                                  , arrf, arrSvSK_Force, arrSvSK_Moment
                                                 ,  arrSvSK_ShaftForce0, arrSvSK_ShaftMoment0
                                                 , arrSvSK_AirForce0, arrSvSK_AirMoment0);
  memcpy(arrFGr, &(arrf[3]), 6 * sizeof(long double));
  arrFGr[3] = arrFGr[3];
  arrFGr[4] = arrFGr[4];
  arrFGr[5] = arrFGr[5];

}
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// нахождение балансировочных параметров
// для виража
// балнсировочные параметры представляют из себя вектор arrXRez
// arrXRez[0] - угол крена
// arrXRez[1] - угол тангажа
// arrXRez[2] - цикл шаг по тангажу Kappa
// arrXRez[3] - общий шаг Fi
// arrXRez[4] - цикл шаг по крену Etta
// arrXRez[5] - дифф шаг DelFi
// INPUT:
// Helic  -вертолет
// Environment - атмосфера
// VAlY - высота
// VAlVx - Скорость
// VAlPsi - угол рыскания
// VAlRad - радиус
// arrX0[6] - начальное значение вектора
//OUTPUT:
// arrXRez[6]  - искомый вектор
// возвращает true если решение найдено
// false в противном случае
bool TBallanceCalc::calcBallParamsForTurn(const THelic Helic,const TEnvironment Environment
                                                      ,const long double VAlY,const long double VAlVx,const long double VAlPsi
                                                      ,const long double VAlRad,long double *arrX0, long double *arrXRez )
{
    long double arrXCur[6] = {0.}
                ,arrFGr[6] = {0.}  // правая часть сиситемы уравнений Fgr(x) = 0
                ,arrJacFgr[36] = {0.} // якобиан
                ,arrJacFGrInv[36] = {0.}
                ,arrDel[6] = {0.} ; // вектор невязки

    // вектор приращений аргументов для вычислкеения якобиана
    long double arr_dx[6] = {0.00001,0.00001,0.00001,0.00001,0.00001,0.00001};

    memcpy(arrXCur, arrX0 , 6 * sizeof(long double));
    for (int i = 0; i < 40; ++i)
    {
        calcFgr_and_JacFgr_ForTurn( Helic, Environment
                ,VAlY, VAlVx,  VAlPsi,VAlRad,arrXCur, arrFGr, arrJacFgr);
        if (!InverseMtrx6(arrJacFgr, arrJacFGrInv))
        {
            return false;
        }
        MtrxMultMatrx(arrJacFGrInv,6, 6, arrFGr,1, arrDel) ;
        MtrxMinusMatrx(arrXCur, arrDel,6, 1, arrXRez);
        long double valDel = NormVect(arrDel, 6);
        long double valDelF = NormVect(arrFGr, 6);

        if(fabsl(valDel) < 0.000001)
        {
            if(valDelF > 0.1)
            {
              return false;
            }
            else
            {
                // ДЛЯ ОТЛАДКИ ПОТОМ УБРАТЬ
                       //     calcFgr_and_JacFgr_ForTurn( Helic, Environment
                         //           ,VAlY, VAlVx,  VAlPsi,VAlRad,arrXCur, arrFGr, arrJacFgr);
                 ///

              return true;
            }

        }

        memcpy(arrXCur, arrXRez , 6 * sizeof(long double));
    }

    return false;
}
//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
//вычисление вектор функции правой части уравнения балансировки
//и ее якобиана
// для равномерного прямолинейного движения
// балансировочные параметры представляют из себя вектор arrXRez
// arrXCur[0] - угол крена
// arrXCur[1] - угол тангажа
// arrXCur[2] - цикл шаг по тангажу Kappa
// arrXCur[3] - общий шаг Fi
// arrXCur[4] - цикл шаг по крену Etta
// arrXCur[5] - дифф шаг DelFi
// INPUT:
// Helic  -вертолет
// Environment - атмосфера
// VAlY - высота
// VAlVx - Скорость
// VAlPsi - угол рыскания
// VAlRad - радиус
// arrXCur[6]  - вектор переменных
//OUTPUT:
// arrFGr[6]  - вектор функция
// arrJacFgr [36] - ee якобиан
void TBallanceCalc::calcFgr_and_JacFgr_ForTurn(const THelic Helic,const TEnvironment Environment
        , long double VAlY, long double VAlVx, long double VAlPsi, const long double VAlRad
        ,long double *arrX, long double *arrFGr, long double *arrJacFgr)
{
    long double arrJacFgrTransp[36] = {0.} // транспогнированный якобиан
                ,arrDel[6] = {0.}
                ,arrTemp[6] = {0.}             ; // вектор невязки

    // вектор приращений аргументов для вычислкеения якобиана
    long double arr_dx[6] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};

    long double arrFGrCur[6] = {0.}, arrXCur[6] = {0.};

    calcFgr_ForTurn( Helic, Environment
            ,  VAlY,  VAlVx, VAlPsi,VAlRad, arrX, arrFGr);

    for (int i = 0; i < 6; ++i)
    {
       memcpy(arrXCur, arrX , 6 * sizeof(long double));
       arrXCur[i] += arr_dx[i] ;
      // calcFgr_ForSteadyLineMoving( Helic, Environment
       //   ,  VAlY,  VAlVx, VAlPsi,arrXCur, arrFGrCur);
       calcFgr_ForTurn( Helic, Environment
                      ,  VAlY,  VAlVx, VAlPsi,VAlRad,arrXCur, arrFGrCur);
       MtrxMinusMatrx(arrFGrCur,arrFGr,6, 1, arrTemp);
       MatrxDivideScalar(arrTemp, 1, 6, arr_dx[i], &(arrJacFgrTransp[ 6 * i]));
    }

    MatrTransp(arrJacFgrTransp, 6, 6, arrJacFgr);

}
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//вычисление вектор функции правой части уравненипя балансировки
// для равномерного прямолинейного движения
// балнсировочные параметры представляют из себя вектор arrXRez
// arrXCur[0] - угол крена
// arrXCur[1] - угол тангажа
// arrXCur[2] - цикл шаг по тангажу Kappa
// arrXCur[3] - общий шаг Fi
// arrXCur[4] - цикл шаг по крену Etta
// arrXCur[5] - дифф шаг DelFi
// INPUT:
// Helic  -вертолет
// Environment - атмосфера
// VAlY - высота
// VAlVx - Скорость
// VAlPsi - угол рыскания
// VAlRad - радиус поворота
// arrXCur[6]  - вектор переменных
//OUTPUT:
// arrFGr[6]  - вектор функция
void TBallanceCalc::calcFgr_ForTurn(const THelic Helic,const TEnvironment Environment
                                    , long double VAlY, long double VAlVx, long double VAlPsi, const long double VAlRad
                                    ,long double *arrXCur, long double *arrFGr)
{
  TPartHelicTraj PartHelicTraj( Helic, Environment, 0.,0.);
  memset(PartHelicTraj.marrPhaseVect, 0, 12 * sizeof(long double));
  PartHelicTraj.marrPhaseVect[1] = VAlY;
  PartHelicTraj.marrPhaseVect[3] = VAlVx;
  PartHelicTraj.marrPhaseVect[11]= VAlPsi; // psi

  PartHelicTraj.marrPhaseVect[9]  = arrXCur[0];// Gamma
  PartHelicTraj.marrPhaseVect[10] = arrXCur[1]; // NU

  // вычисление вектора угловых скоростей вращения в СвСК

    long double val_dPsi_po_dt = VAlVx / VAlRad;
    PartHelicTraj.marrPhaseVect[6] = val_dPsi_po_dt * sinl(PartHelicTraj.marrPhaseVect[10]);
    PartHelicTraj.marrPhaseVect[7] = val_dPsi_po_dt * cosl(PartHelicTraj.marrPhaseVect[10]) * cosl(PartHelicTraj.marrPhaseVect[9]);
    PartHelicTraj.marrPhaseVect[8] = -val_dPsi_po_dt * cosl(PartHelicTraj.marrPhaseVect[10]) * sinl(PartHelicTraj.marrPhaseVect[9]);

  PartHelicTraj.doAirDensity();

  long double    arrf[12] = {0.};
  long double arrSvSK_Force[3] = {0.}, arrSvSK_Moment[3] = {0.}
          ,  arrSvSK_ShaftForce0[3] = {0.}, arrSvSK_ShaftMoment0[3] = {0.}
          , arrSvSK_AirForce0 [3] = {0.}, arrSvSK_AirMoment0 [3] = {0.} ;


  PartHelicTraj.calcRightMemberOfDifEqSystem(arrXCur[3], arrXCur[2]
                                                  , arrXCur[4], arrXCur[5]
                                                  , arrf, arrSvSK_Force, arrSvSK_Moment
                                                 ,  arrSvSK_ShaftForce0, arrSvSK_ShaftMoment0
                                                 , arrSvSK_AirForce0, arrSvSK_AirMoment0);
  memcpy(arrFGr, &(arrf[3]), 6 * sizeof(long double));
  arrFGr[2] -= VAlVx * VAlVx / VAlRad;

}
//----------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// нахождение балансировочных параметров
// для вращения
// балнсировочные параметры представляют из себя вектор arrXRez
// arrXRez[0] - угол крена
// arrXRez[1] - угол тангажа
// arrXRez[2] - цикл шаг по тангажу Kappa
// arrXRez[3] - общий шаг Fi
// arrXRez[4] - цикл шаг по крену Etta
// arrXRez[5] - дифф шаг DelFi
// INPUT:
// Helic  -вертолет
// Environment - атмосфера
// VAlY - высота
// VAlVx - Скорость
// VAlPsi - угол рыскания
// VAlRad - радиус
// arrX0[6] - начальное значение вектора
//OUTPUT:
// arrXRez[6]  - искомый вектор
// возвращает true если решение найдено
// false в противном случае
bool TBallanceCalc::calcBallParamsForRotation(const THelic Helic,const TEnvironment Environment
                                                      ,const long double VAlY,const long double VAl_dPsi_po_dt
                                                      ,long double *arrX0, long double *arrXRez )
{
    long double arrXCur[6] = {0.}
                ,arrFGr[6] = {0.}  // правая часть сиситемы уравнений Fgr(x) = 0
                ,arrJacFgr[36] = {0.} // якобиан
                ,arrJacFGrInv[36] = {0.}
                ,arrDel[6] = {0.} ; // вектор невязки

    // вектор приращений аргументов для вычислкеения якобиана
    long double arr_dx[6] = {0.00001,0.00001,0.00001,0.00001,0.00001,0.00001};

    memcpy(arrXCur, arrX0 , 6 * sizeof(long double));
    for (int i = 0; i < 40; ++i)
    {
        calcFgr_and_JacFgr_ForRotation( Helic, Environment
                ,VAlY,  VAl_dPsi_po_dt,arrXCur, arrFGr, arrJacFgr);
        if (!InverseMtrx6(arrJacFgr, arrJacFGrInv))
        {
            return false;
        }
        MtrxMultMatrx(arrJacFGrInv,6, 6, arrFGr,1, arrDel) ;
        MtrxMinusMatrx(arrXCur, arrDel,6, 1, arrXRez);
        long double valDel = NormVect(arrDel, 6);
        long double valDelF = NormVect(arrFGr, 6);

        if(fabsl(valDel) < 0.000001)
        {
            if(valDelF > 0.1)
            {
              return false;
            }
            else
            {
                // ДЛЯ ОТЛАДКИ ПОТОМ УБРАТЬ
                calcFgr_and_JacFgr_ForRotation( Helic, Environment
                        ,VAlY,  VAl_dPsi_po_dt,arrXCur, arrFGr, arrJacFgr);
                       //     calcFgr_and_JacFgr_ForTurn( Helic, Environment
                         //           ,VAlY, VAlVx,  VAlPsi,VAlRad,arrXCur, arrFGr, arrJacFgr);
                 ///

              return true;
            }

        }

        memcpy(arrXCur, arrXRez , 6 * sizeof(long double));
    }

    return false;
}
//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
//вычисление вектор функции правой части уравнения балансировки
//и ее якобиана
// для равномерного прямолинейного движения
// балансировочные параметры представляют из себя вектор arrXRez
// arrXCur[0] - угол крена
// arrXCur[1] - угол тангажа
// arrXCur[2] - цикл шаг по тангажу Kappa
// arrXCur[3] - общий шаг Fi
// arrXCur[4] - цикл шаг по крену Etta
// arrXCur[5] - дифф шаг DelFi
// INPUT:
// Helic  -вертолет
// Environment - атмосфера
// VAlY - высота
// VAlVx - Скорость
// VAlPsi - угол рыскания
// VAlRad - радиус
// arrXCur[6]  - вектор переменных
//OUTPUT:
// arrFGr[6]  - вектор функция
// arrJacFgr [36] - ee якобиан
void TBallanceCalc::calcFgr_and_JacFgr_ForRotation(const THelic Helic,const TEnvironment Environment
        , long double VAlY, const long double VAl_dPsi_po_dt
        ,long double *arrX, long double *arrFGr, long double *arrJacFgr)
{
    long double arrJacFgrTransp[36] = {0.} // транспогнированный якобиан
                ,arrDel[6] = {0.}
                ,arrTemp[6] = {0.}             ; // вектор невязки

    // вектор приращений аргументов для вычислкеения якобиана
    long double arr_dx[6] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};

    long double arrFGrCur[6] = {0.}, arrXCur[6] = {0.};

    calcFgr_ForRotation( Helic, Environment
            ,  VAlY,  VAl_dPsi_po_dt, arrX, arrFGr);

    for (int i = 0; i < 6; ++i)
    {
       memcpy(arrXCur, arrX , 6 * sizeof(long double));
       arrXCur[i] += arr_dx[i] ;
       calcFgr_ForRotation( Helic, Environment
               ,  VAlY,  VAl_dPsi_po_dt,arrXCur, arrFGrCur);
       MtrxMinusMatrx(arrFGrCur,arrFGr,6, 1, arrTemp);
       MatrxDivideScalar(arrTemp, 1, 6, arr_dx[i], &(arrJacFgrTransp[ 6 * i]));
    }

    MatrTransp(arrJacFgrTransp, 6, 6, arrJacFgr);

}
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//вычисление вектор функции правой части уравненипя балансировки
// для равномерного прямолинейного движения
// балнсировочные параметры представляют из себя вектор arrXRez
// arrXCur[0] - угол крена
// arrXCur[1] - угол тангажа
// arrXCur[2] - цикл шаг по тангажу Kappa
// arrXCur[3] - общий шаг Fi
// arrXCur[4] - цикл шаг по крену Etta
// arrXCur[5] - дифф шаг DelFi
// INPUT:
// Helic  -вертолет
// Environment - атмосфера
// VAlY - высота
// VAlVx - Скорость
// VAlPsi - угол рыскания
// VAlRad - радиус поворота
// arrXCur[6]  - вектор переменных
//OUTPUT:
// arrFGr[6]  - вектор функция
void TBallanceCalc::calcFgr_ForRotation(const THelic Helic,const TEnvironment Environment
                                    , long double VAlY,const long double VAl_dPsi_po_dt
                                    ,long double *arrXCur, long double *arrFGr)
{
  TPartHelicTraj PartHelicTraj( Helic, Environment, 0.,0.);
  memset(PartHelicTraj.marrPhaseVect, 0, 12 * sizeof(long double));
  PartHelicTraj.marrPhaseVect[1] = VAlY;
  PartHelicTraj.marrPhaseVect[3] = 0.;
  PartHelicTraj.marrPhaseVect[11]= 0.; // psi

  PartHelicTraj.marrPhaseVect[9]  = arrXCur[0];// Gamma
  PartHelicTraj.marrPhaseVect[10] = arrXCur[1]; // NU

  // вычисление вектора угловых скоростей вращения в СвСК


    PartHelicTraj.marrPhaseVect[6] = VAl_dPsi_po_dt * sinl(PartHelicTraj.marrPhaseVect[10]);
    PartHelicTraj.marrPhaseVect[7] = VAl_dPsi_po_dt * cosl(PartHelicTraj.marrPhaseVect[10]) * cosl(PartHelicTraj.marrPhaseVect[9]);
    PartHelicTraj.marrPhaseVect[8] = -VAl_dPsi_po_dt * cosl(PartHelicTraj.marrPhaseVect[10]) * sinl(PartHelicTraj.marrPhaseVect[9]);

  PartHelicTraj.doAirDensity();

  long double    arrf[12] = {0.};
  long double arrSvSK_Force[3] = {0.}, arrSvSK_Moment[3] = {0.}
          ,  arrSvSK_ShaftForce0[3] = {0.}, arrSvSK_ShaftMoment0[3] = {0.}
          , arrSvSK_AirForce0 [3] = {0.}, arrSvSK_AirMoment0 [3] = {0.} ;


  PartHelicTraj.calcRightMemberOfDifEqSystem(arrXCur[3], arrXCur[2]
                                                  , arrXCur[4], arrXCur[5]
                                                  , arrf, arrSvSK_Force, arrSvSK_Moment
                                                 ,  arrSvSK_ShaftForce0, arrSvSK_ShaftMoment0
                                                 , arrSvSK_AirForce0, arrSvSK_AirMoment0);
  memcpy(arrFGr, &(arrf[3]), 6 * sizeof(long double));


}
//----------------------------------------------------------------------------------

// задан угол тангажа при висении на нулевой скорости
// задана центровка и расстояние от точки приложения силы каппа до центра масс (отцентованного вертолета)
// точка приложения силы каппа - это средняя точка на НВ между верхним и нижним винтами
// тогда можно найти угол заклинения
// Это по картинке 7.28 РЛЭ
// INPUT:
//valNu - угол тангажа в рад
//valDist - расстояние
//valXT - центровка
long double TBallanceCalc::calcAlfaZakl(const long double valNu,const long double valXT,const long double valDist)
{
 /* long double a= valXT/ valDist;
  long double alfaZakl0 = 0.;
  for (int i = 0; i < 40;++i)
  {
    long double del = (tanl(valNu) * cosl(alfaZakl0)  - sinl(alfaZakl0) - a)
            /(-tanl(valNu) * sinl(alfaZakl0) - cosl(alfaZakl0));
    alfaZakl0 -= del;
    if (fabsl(del) <= 0.00001)
    {
        return alfaZakl0;
    }

  }
  return -9999999999999999.;*/
    return valNu- asinl(valXT * cosl(valNu)/valDist);
}

