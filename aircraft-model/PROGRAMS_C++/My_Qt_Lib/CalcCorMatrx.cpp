//---------------------------------------------------------------------------

#include "CalcCorMatrx.h"
#include <string.h>
#include <math.h>
#include "MatrixProccess.h"


//---------------------------------------------------------------------------
// поворот по часовой стрелке
// 1. заданы 2 правых системы координат СК1 и СК2.
// система координат СК2 получается из СК1 путем поворота на угол а
// по часовой стрелке относительно оси OX СК1
// тогда matrP - матрица перехода из СК2 в СК1,
// то есть, если S2 - вектор координат точки в СК2, а S1 вектор координат той же точки в СК2,
// то они связаны соотношением
// S1 = matrP * S2
// 2. заданы правая  система координат СК1.
// задан вектор S1 в СК1
// производится поворот СК1 относительно оси OX по часовой стрелке
// на угол а вместе с вектором S1
// обозначим через S2 координаты повернутого векьтора в СК1
// тогда, координаты векторов S1 и S2 связаны след соотношением
// S2 = matrP * S1
void calcMatrP1(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));
  matrP[0] = 1;
  matrP[4] = cos(a);
  matrP[5] = sin(a);
  matrP[7] = -sin(a);
  matrP[8] = cos(a);

}
//---------------------------------------------------------------------------
// поворот  против часовой стрелки
// 1. заданы 2 правых системы координат СК1 и СК2.
// система координат СК2 получается из СК1 путем поворота на угол а
// против часовой стрелки относительно оси OY СК1
// тогда matrP - матрица перехода из СК2 в СК1,
// то есть, если S2 - вектор координат точки в СК2, а S1 вектор координат той же точки в СК2,
// то они связаны соотношением
// S1 = matrP * S2
// 2. заданы правая  система координат СК1.
// задан вектор S1 в СК1
// производится поворот СК1 относительно оси OY против  часовой стрелки
// на угол а вместе с вектором S1
// обозначим через S2 координаты повернутого вектора в СК1
// тогда, координаты векторов S1 и S2 связаны след соотношением
// S2 = matrP * S1
void calcMatrP2(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));
  matrP[4] = 1;
  matrP[0] = cos(a);
  matrP[2] = sin(a);
  matrP[6] = -sin(a);
  matrP[8] = cos(a);

}

//---------------------------------------------------------------------------
// поворот  по  часовой стрелке относительно оси OY
// 1. заданы 2 правых системы координат СК1 и СК2.
// система координат СК2 получается из СК1 путем поворота на угол а
// по часовой стрелке относительно оси OY СК1
// тогда matrP - матрица перехода из СК2 в СК1,
// то есть, если S2 - вектор координат точки в СК2, а S1 вектор координат той же точки в СК2,
// то они связаны соотношением
// S1 = matrP * S2
// 2. заданы правая  система координат СК1.
// задан вектор S1 в СК1
// производится поворот СК1 относительно оси OY против  часовой стрелки
// на угол а вместе с вектором S1
// обозначим через S2 координаты повернутого вектора в СК1
// тогда, координаты векторов S1 и S2 связаны след соотношением
// S2 = matrP * S1
void calcMatrP2_AlogWatchArrow(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));
  matrP[4] = 1;
  matrP[0] = cos(a);
  matrP[2] = -sin(a);
  matrP[6] = sin(a);
  matrP[8] = cos(a);

}
//---------------------------------------------------------------------------
// поворот по часовой стрелке
// заданы 2 правых системы координат СК1 и СК2.
// система координат СК2 получается из СК1 путем поворота на угол а
// по часовой стрелке относительно оси OZ СК1
// тогда matrP - матрица перехода из СК2 в СК1,
// то есть, если S2 - вектор координат точки в СК2, а S1 вектор координат той же точки в СК2,
// то они связаны соотношением
// S1 = matrP * S2
// 2. заданы правая  система координат СК1.
// задан вектор S1 в СК1
// производится поворот СК1 относительно оси OZ по часовой стрелке
// на угол а вместе с вектором S1
// обозначим через S2 координаты повернутого векьтора в СК1
// тогда, координаты векторов S1 и S2 связаны след соотношением
// S2 = matrP * S1
void calcMatrP3(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));
  matrP[8] = 1;
  matrP[0] = cos(a);
  matrP[1] = sin(a);
  matrP[3] = -sin(a);
  matrP[4] = cos(a);
}

//---------------------------------------------------------------------------
void calcMatr_dP1_po_da(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));

  matrP[4] = -sin(a);
  matrP[5] = cos(a);
  matrP[7] = -cos(a);
  matrP[8] = -sin(a);

}
//---------------------------------------------------------------------------
void calcMatr_d2P1_po_da2(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));

  matrP[4] = -cos(a);
  matrP[5] = -sin(a);
  matrP[7] = sin(a);
  matrP[8] = -cos(a);

}
//---------------------------------------------------------------------------
void calcMatr_d2P1otMinus_a_po_da2(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));

  matrP[4] = -cos(a);
  matrP[5] = sin(a);
  matrP[7] = -sin(a);
  matrP[8] = -cos(a);


}
//---------------------------------------------------------------------------
void calcMatr_dP2_po_da(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));

  matrP[0] = -sin(a);
  matrP[2] = cos(a);
  matrP[6] = -cos(a);
  matrP[8] = -sin(a);

}

//---------------------------------------------------------------------------
void calcMatr_dP1_otMinus_a_po_da(const double a, double*matrP)
{
    memset(matrP,0,9*sizeof(double));

    matrP[4] = -sin(a);
    matrP[5] = -cos(a);
    matrP[7] = cos(a);
    matrP[8] = -sin(a);

}
//---------------------------------------------------------------------------
void calcMatr_dP3_po_da(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));

  matrP[0] = -sin(a);
  matrP[1] = cos(a);
  matrP[3] = -cos(a);
  matrP[4] = -sin(a);

}



//---------------------------------------------------------------------------
void calcMatr_d2P3_po_da2(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));

  matrP[0] = -cos(a);
  matrP[1] = -sin(a);
  matrP[3] = sin(a);
  matrP[4] = -cos(a);

}
//-------------------------------------


void calcMatr_ASK_v_KGSK(double *arrMu,double *matrPereh_ASK_V_KGSK)
{

     double  matrPereh_ASK_V_PSK[9] = {0.}, matrPereh_PSK_V_KGSK[9] ={0.} ;
     calcMatr_ASK_v_PSK( &arrMu[3],matrPereh_ASK_V_PSK) ;
     calcMatr_PSK_v_KGSK(arrMu,matrPereh_PSK_V_KGSK);
     MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, matrPereh_ASK_V_PSK,3, matrPereh_ASK_V_KGSK ) ;


}

void calcMatr_ASK_v_PSK(double *arrMu,double *matrPereh_ASK_V_PSK)
{

    const double valBet = arrMu[0] ;
    const double valEps = arrMu[1] ;
    double matrP3Bet[9]={0},matrP1Eps[9]={0};
    calcMatrP3(valBet, matrP3Bet) ;
    calcMatrP1(-valEps, matrP1Eps) ;
    MtrxMultMatrx( matrP3Bet,3, 3, matrP1Eps,3, matrPereh_ASK_V_PSK);

}

void calcMatr_PSK_v_KGSK(double *arrMu,double *matrPereh_PSK_V_KGSK)
{
    const double valQ   = arrMu[0] ;
    const double valPsi = arrMu[1] ;
    const double valTet = arrMu[2] ;

    double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

    calcMatrP3(valQ, matrP3Q) ;
    calcMatrP1(valPsi, matrP1Psi) ;
    calcMatrP2(valTet, matrP2Tet) ;
    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, matrPereh_PSK_V_KGSK);

}

void calcMatr_PSK_v_KGSK_Eilers(double *arrMu,double *matrPereh_PSK_V_KGSK)
{
    const double valQ   = arrMu[0] ;
    const double valPsi = arrMu[1] ;
    const double valTet = arrMu[2] ;

    double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};
  // поворот по часовой стрелке
    calcMatrP3(-valQ, matrP3Q) ;
    // поворот по часовой стрелке:
    calcMatrP1(-valPsi, matrP1Psi) ;
    // поворот по часовой стрелке:
    calcMatrP2(valTet, matrP2Tet) ;
    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, matrPereh_PSK_V_KGSK);

}



void calc_dF_po_dQ_sq(const double valBet,const double valEps,double *matrRez)
{
    double arr[3]={0},arr0[3]={0};
    arr[0] = cos(valBet) * cos ( valEps) ;
    arr[1] = -sin(valBet) * cos ( valEps) ;
    memcpy(arr0,arr,3*sizeof(double));
    MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);
}
void calc_dF_po_dPsi_sq(const double valBet,const double valEps,double *matrRez)
{
    double arr[3]={0},arr0[3]={0};
    arr[1] = sin(valEps) ;
    arr[2] = -cos(valBet) * cos ( valEps) ;
    memcpy(arr0,arr,3*sizeof(double));
    MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);

}
void calc_dF_po_dTet_sq(const double valBet,const double valEps,double *matrRez)
{
    double arr[3]={0},arr0[3]={0};
    arr[0] = sin(valEps) ;
    arr[2] = -sin(valBet) * cos ( valEps) ;
    memcpy(arr0,arr,3*sizeof(double));
    MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);

}
void calc_dF_po_dBet_sq(const double valBet,const double valEps,double *matrRez)
{
    double arr[3]={0},arr0[3]={0};
    arr[0] = cos(valBet) * cos ( valEps) ;
    arr[1] = -sin(valBet) * cos ( valEps) ;
    memcpy(arr0,arr,3*sizeof(double));
    MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);

}
void calc_dF_po_dEps_sq(const double valBet,const double valEps,double *matrRez)
{
    double arr[3]={0},arr0[3]={0};
    arr[0] =  -sin(valBet) *sin(valEps) ;
    arr[1] = -cos(valBet) * sin ( valEps) ;
    arr[2] =  cos ( valEps) ;
    memcpy(arr0,arr,3*sizeof(double));
    MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);

}
void createExtendMtrx(double *arrInp, double * arrOut)
{
  memset(arrOut, 0, 36 * sizeof(double));
 for (int i =0; i < 3; i++)
 for (int j =0; j < 3; j++)
 {
  arrOut[ i * 6 + j] =  arrInp [3 * i + j];
  arrOut[ (3 + i) * 6 + 3 + j] =  arrInp [3 * i + j];
 }
}
void calcMatrJ1(const double valBet,const double valEps,double *arrJ1)
{
 memset(arrJ1, 0, 9 * sizeof(double)) ;
 arrJ1[1] = cos(valEps);
 arrJ1[2] = -sin(valEps);
 arrJ1[3] = -cos(valEps);
 arrJ1[6] = sin(valEps);

}
void calcMatrJ2(const double valBet,const double valEps,double *arrJ2)
{
 memset(arrJ2, 0, 9 * sizeof(double)) ;
 arrJ2[1] = -sin(valBet)*sin(valEps);
 arrJ2[2] = -sin(valBet)*cos(valEps);
 arrJ2[3] = sin(valBet)*sin(valEps);
 arrJ2[5] = cos(valBet);
 arrJ2[6] = sin(valBet)*cos(valEps);
 arrJ2[7] = -cos(valBet);


}
void calcMatrJ3(const double valBet,const double valEps,double *arrJ3)
{
 memset(arrJ3, 0, 9 * sizeof(double)) ;
 arrJ3[1] = cos(valBet)*sin(valEps);
 arrJ3[2] = cos(valBet)*cos(valEps);
 arrJ3[3] = -cos(valBet)*sin(valEps);
 arrJ3[5] = sin(valBet);
 arrJ3[6] = -cos(valBet)*cos(valEps);
 arrJ3[7] = -sin(valBet);

}
void calcMatrJ4(const double valBet,const double valEps,double *arrJ4)
{
 calcMatrJ1( valBet, valEps,arrJ4) ;

}
void calcMatrJ5(const double valBet,const double valEps,double *arrJ5)
{
 memset(arrJ5, 0, 9 * sizeof(double)) ;

 arrJ5 [ 5] = -1 ;
 arrJ5 [ 7] = 1 ;

}
void calcMatrJ1W(const double valBet,const double valEps,double *arrJ1W)
{
   memset(arrJ1W, 0, 36 * sizeof(double)) ;
   double arr[9]={0} ;
   calcMatrJ1( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ1W);
}
void calcMatrJ2W(const double valBet,const double valEps,double *arrJ2W)
{
   memset(arrJ2W, 0, 36 * sizeof(double)) ;
   double arr[9]={0} ;
   calcMatrJ2( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ2W);
}
void calcMatrJ3W(const double valBet,const double valEps,double *arrJ3W)
{
   memset(arrJ3W, 0, 36 * sizeof(double)) ;
   double arr[9]={0} ;
   calcMatrJ3( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ3W);
}
void calcMatrJ4W(const double valBet,const double valEps,double *arrJ4W)
{
   memset(arrJ4W, 0, 36 * sizeof(double)) ;
   double arr[9]={0} ;
   calcMatrJ4( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ4W);
}
void calcMatrJ5W(const double valBet,const double valEps,double *arrJ5W)
{
   memset(arrJ5W, 0, 36 * sizeof(double)) ;
   double arr[9]={0} ;
   calcMatrJ5( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ5W);
}

//--------------------------------------------------------------------------------
// Вычисление корреляционнной матрицы ошибок , возникающих при переводе вектора из ПСК в КГСК.
// ОШибки возникают в следствие того, что палубные углы известны неточно, а с ошибками измерений.
// Задан вектор arrVS[3], измерения палубных углов и дисперсии ошибок измерения этих углов.
// Требуется определить коррел матрицу ошибок преобразованного вектора arrVS
// INPUT:
// VAlQ,  VAlPsi,  VAlTet  - измерения курсового угда, угла килевой и бортовой качек
// VAlDispQ, VAlDispPsi	 , VAlDispTet  - дисперсии ошибок измерений этих углов
// arrVS[3] - вектор в ПСК
// OUTPUT:
// arrVS_KGSK[3] - вектор в КГСК
// arrMtxCorr[9] - корреляционная матрица ошибок
void calcCorMatrx_PSK_KGSK(const double VAlQ, const double VAlPsi, const double VAlTet
     , const double VAlDispQ,const double VAlDispPsi	 , const double VAlDispTet
     , double  *arrVS_PSK,double  *arrVS_KGSK, double *arrMtxCorr)
{
    double arrP1[9] = {0.},arrP2[9] = {0.},arrP3[9] = {0.},
        arr_dP1[9] = {0.},arr_dP2[9] = {0.},arr_dP3[9] = {0.};
    calcMatrP1(VAlPsi, arrP1);
    calcMatr_dP1_po_da(VAlPsi, arr_dP1);

    calcMatrP2(VAlTet, arrP2);
    calcMatr_dP2_po_da(VAlTet, arr_dP2);

    calcMatrP3(VAlQ, arrP3);
    calcMatr_dP3_po_da(VAlQ, arr_dP3);

    double arrJ1[9] = {0.}, arrJ2[9] = {0.}, arrJ3[9] = {0.};
    MtrxMultMatrx_MultMatrx(arr_dP3, arrP1,arrP2,3, arrJ1)  ;
    MtrxMultMatrx_MultMatrx(arrP3, arr_dP1,arrP2,3, arrJ2)  ;
    MtrxMultMatrx_MultMatrx(arrP3, arrP1,arr_dP2,3, arrJ3)  ;
    ///
 memset(arrMtxCorr, 0, 9 * sizeof(double));
 double arrQ[9] = {0.}, arrPsi[9] = {0.}, arrTet[9] = {0.}, arrTemp0[9] = {0.}, arrVS_PSK_Copy[3] ={0.}, arrS_mult_ST[9] = {0.};
 memcpy(arrVS_PSK_Copy, arrVS_PSK, 3 * sizeof(double));
 MtrxMultMatrxTransp(arrVS_PSK,3, 1, arrVS_PSK_Copy,3, arrS_mult_ST) ;
 MtrxMultMatrx_MultMatrxTransp(arrJ1,arrS_mult_ST, arrJ1, 3, arrQ);
 MtrxMultMatrx_MultMatrxTransp(arrJ2,arrS_mult_ST, arrJ2, 3, arrPsi);
 MtrxMultMatrx_MultMatrxTransp(arrJ3,arrS_mult_ST, arrJ3, 3, arrTet);

 MatrxMultScalar(arrQ, 3, 3,VAlDispQ,arrQ) ;
 MatrxMultScalar(arrPsi, 3, 3,VAlDispPsi,arrPsi) ;
 MatrxMultScalar(arrTet, 3, 3,VAlDispTet,arrTet) ;

 MtrxSumMatrx(arrQ, arrPsi,3, 3, arrTemp0) ;
 MtrxSumMatrx(arrTemp0, arrTet,3, 3, arrMtxCorr) ;

 ///

MtrxMultMatrx_MultMatrx(arrP3, arrP1,arrP2,3, arrJ1)  ;
MtrxMultMatrx(arrJ1,3, 3, arrVS_PSK, 1, arrVS_KGSK) ;
}

// матрица перехода parrMtrxPer из исходной прямоугольной сиситемы координат
// в скоростную сиситему координат (CCК), задаваемую вектором arrV
// ось X  ССК направлена по вектору V, ось Y CCК праллельна плоскости OXY
// исходной прямоугольной сиситемы координат, а ось Z дополняет до правой тройки
// parrMtrxPer - это матрица состоящая из столбцов вектров единичных ортов
void calcMatrxPer_from_DecartPrSK_To_SSK(double *arrV, double * parrMtrxPer)
{
    double val_v0 = NormVect2(arrV);
    memset(parrMtrxPer, 0, 9 * sizeof(double));
    if (val_v0 < 0.0000000001)
  {
    parrMtrxPer[2] = 1.;
    parrMtrxPer[3] = 1.;
    parrMtrxPer[7] = 1.;
    return;
  }
  double arr_v[3] = {0.};
  memcpy(arr_v, arrV, 3 *  sizeof(double));
    NormalizeVect3(arr_v) ;
    double val_v = NormVect2(arr_v);
  parrMtrxPer[0] = arr_v [0];
  parrMtrxPer[1] = arr_v [1];
    parrMtrxPer[2] = arr_v[2];
    parrMtrxPer[3] = -arr_v [1] / val_v;
    parrMtrxPer[4] = arr_v [0] / val_v;
  parrMtrxPer[5] = 0.;
    parrMtrxPer[6] = -arr_v [2] * arr_v [0] / val_v;
    parrMtrxPer[7] = -arr_v [2] * arr_v [1] / val_v;
  parrMtrxPer[8] =  val_v ;
  return;
}


// матрица перехода parrMtrxPer из скоростной  сиситемы координат (CCК), задаваемой вектором arrV
// в исходную  прямоугольную  сиситему координат
// ось X  ССК направлена по вектору V, ось Y CCК праллельна плоскости OXY
// исходной прямоугольной сиситемы координат, а ось Z дополняет до правой тройки
// parrMtrxPer - это матрица состоящая из столбцов вектров единичных ортов
// CCК  в исходной прямоуг сиситемы координат
void calcMatrxPer_from_SSK_To_DecartPrSK(double *arrV, double * parrMtrxPer)
{
  double arrT[9] = {0.};
  calcMatrxPer_from_DecartPrSK_To_SSK(arrV, arrT) ;
  MatrTransp(arrT, 3, 3, parrMtrxPer);
  return;
}

// Пенресчет прямоугольных координат в сферические
//  arrInp[3] - прямоугольные координаты
// valBet, valR,valEps - сферические координаты
void recalcCoord_INTO_Spherical(double *arrInp, double &valR, double &valBet, double &valEps)
{
    valR = sqrt ( arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1] +arrInp[2] * arrInp[2] ) ;

    valEps  = asin(arrInp[2] / valR ) ;

    valBet  = asin(arrInp[0]/ sqrt (  arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1])) ;
   // valBet  = atan2(arrInp[0], arrInp[1]);
    if (arrInp[1] < 0.)
     {
        if (arrInp[0] < 0.)
        {
          valBet = -M_PI - valBet ;
        }
        else
        {
           valBet = M_PI - valBet ;
        }
     }


}
//-------------------------------------------
// Пенресчет сферических  координат в прямоугольные
//  arrS[3] - прямоугольные координаты
// valBet, valR,valEps - сферические координаты
void recalcSphericalCoord_INTO_Rectangular(const double valR,const  double valBet,const  double valEps
                                           , double *arrS)
{
   arrS[0] =  valR * cos(valEps)* sin(valBet);
   arrS[1] =  valR * cos(valEps)* cos(valBet);
   arrS[2] =  valR * sin(valEps);
}
//-------------------------------------------------------
//Вычисление матрицы частных производных векто-функции сферических кооординат
// по прямоугольным. порядок следования сферических координат -
// курс угол, дальность, угол места
void calc_dSpherical_po_Rectangular(double *arrS, double *arr_dB_po_dS)
{
    double r = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
    double R = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]+ arrS[2] * arrS[2]);

    arr_dB_po_dS[0] = arrS[1]/r/r;
    arr_dB_po_dS[1] = -arrS[0]/r/r;
    arr_dB_po_dS[2] = 0.;

    arr_dB_po_dS[3] = arrS[0]/R;
    arr_dB_po_dS[4] = arrS[1]/R;
    arr_dB_po_dS[5] = arrS[2]/R;

    arr_dB_po_dS[6] = -arrS[0]* arrS[2]/R/R/r;
    arr_dB_po_dS[7] = -arrS[1]* arrS[2]/R/R/r;
    arr_dB_po_dS[8] = r/R/R;
}


//Вычисление матрицы частных производных векто-функции сферических кооординат
// по прямоугольным. порядок следования сферических координат -
// курс угол, дальность, угол места
void calc_dSpherical_po_Rectangular_MultCoeff(double valCoeff, double *arrS, double *arr_dB_po_dS_MultC)
{
    double r = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
    double R = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]+ arrS[2] * arrS[2]);

    arr_dB_po_dS_MultC[0] = arrS[1]/r*valCoeff/r;
    arr_dB_po_dS_MultC[1] = -arrS[0]/r*valCoeff/r;
    arr_dB_po_dS_MultC[2] = 0.;

    arr_dB_po_dS_MultC[3] = arrS[0]/R *valCoeff;
    arr_dB_po_dS_MultC[4] = arrS[1]/R *valCoeff;
    arr_dB_po_dS_MultC[5] = arrS[2]/R *valCoeff;

    arr_dB_po_dS_MultC[6] = -arrS[0]/r*valCoeff/R * arrS[2]/R;
    arr_dB_po_dS_MultC[7] = -arrS[1]/r*valCoeff/R * arrS[2]/R;
    arr_dB_po_dS_MultC[8] = r/R* valCoeff/R;
}
//-------------------------------------------------------
//Вычисление матрицы частных производных векто-функции угловых сферических кооординат
// по прямоугольным. порядок следования сферических координат -
// курс угол,  угол места
void calc_dSpherical_po_Rectangular_AnglesOnly(double *arrS, double *arr_dB_po_dS)
{
    double r = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
    double R = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]+ arrS[2] * arrS[2]);

    arr_dB_po_dS[0] = arrS[1]/r/r;
    arr_dB_po_dS[1] = -arrS[0]/r/r;
    arr_dB_po_dS[2] = 0.;

    arr_dB_po_dS[3] = -arrS[0]* arrS[2]/R/R/r;
    arr_dB_po_dS[4] = -arrS[1]* arrS[2]/R/R/r;
    arr_dB_po_dS[5] = r/R/R;

}
//-------------------------------------------------------
//Вычисление матрицы частных производных векто-функции угловых сферических кооординат
// по прямоугольным. порядок следования сферических координат -
// курс угол,  угол места
void calc_dSpherical_po_Rectangular_AnglesOnly_MultCoeff(double valCoeff, double *arrS, double *arr_dB_po_dS_MultC)
{
    double r = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
    double R = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]+ arrS[2] * arrS[2]);

    arr_dB_po_dS_MultC[0] = arrS[1]/r*valCoeff/r;
    arr_dB_po_dS_MultC[1] = -arrS[0]/r*valCoeff/r;
    arr_dB_po_dS_MultC[2] = 0.;

    arr_dB_po_dS_MultC[3] = -arrS[0]/r*valCoeff/R * arrS[2]/R;
    arr_dB_po_dS_MultC[4] = -arrS[1]/r*valCoeff/R * arrS[2]/R;
    arr_dB_po_dS_MultC[5] = r/R* valCoeff/R;

}
//-------------------------------------------------------
//Вычисление матрицы  производных матрицы перехода их АСК в ПСК
// по углу Betta (курсовому)
void calc_dM_ask_to_psk_po_Betta(const double valBet,const double valEps,double *arr_dM_po_dBet)
{
    memset(arr_dM_po_dBet, 0, 9 * sizeof(double));
    arr_dM_po_dBet[0] = -sin(valBet);
    arr_dM_po_dBet[1] = cos(valBet) * cos (valEps);
    arr_dM_po_dBet[2] = -cos(valBet) * sin (valEps);

    arr_dM_po_dBet[3] = -cos(valBet);
    arr_dM_po_dBet[4] = -sin(valBet) * cos (valEps);
    arr_dM_po_dBet[5] = sin(valBet) * sin (valEps);
}
//-------------------------------------------------------
//Вычисление матрицы  производных матрицы перехода их АСК в ПСК
// по углу Eps (места)
void calc_dM_ask_to_psk_po_Eps(const double valBet,const double valEps,double *arr_dM_po_dEps)
{
    memset(arr_dM_po_dEps, 0, 9 * sizeof(double));

    arr_dM_po_dEps[1] = -sin(valBet) * sin (valEps);
    arr_dM_po_dEps[2] = -sin(valBet) * cos (valEps);

    arr_dM_po_dEps[4] = -cos(valBet) * sin (valEps);
    arr_dM_po_dEps[5] = -cos(valBet) * cos (valEps);

    arr_dM_po_dEps[7] =  cos (valEps);
    arr_dM_po_dEps[8] =  -sin (valEps);
}

//-------------------------------------------------------
//Вычисление матрицы частных производных вектор-функции прямоугольных кооординат
// по сферическим.
// порядок следования сферических координат -
//arrV[3]- курс угол, дальность, угол места
void calc_dRectangularl_po_dSpherical(const double *arrV, double *arr_dA_po_dV)
{
    memset(arr_dA_po_dV, 0, 9 * sizeof(double));
    double e = arrV[2];
    double q = arrV[0];
    double r = arrV[1];

    arr_dA_po_dV[0] = r * cos(q) * cos (e);
    arr_dA_po_dV[1] = sin(q) * cos (e);
    arr_dA_po_dV[2] = -r * sin(q) * sin (e);

    arr_dA_po_dV[3] = -r * sin(q) * cos (e);
    arr_dA_po_dV[4] = cos(q) * cos (e);
    arr_dA_po_dV[5] = -r * cos(q) * sin (e);
    arr_dA_po_dV[6] = 0;
    arr_dA_po_dV[7] = sin(e);
    arr_dA_po_dV[8] = r * cos(e);
}
//-------------------------------------------------
// вычисление частной призводной матрицы Якобиана
// преобразования прямоугольных координат в чферические по переменной
// X (первой переменной)
void calc_d2B_po_ds_po_dS0(double *arrS, double *arr_d2B_po_dS_po_dS0)
{

    double r = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
    double R = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]+ arrS[2] * arrS[2]);
    double x =arrS[0];
    double y =arrS[1];
    double z =arrS[2];
    arr_d2B_po_dS_po_dS0[0] = -2. * x * y / (r * r * r * r);
    arr_d2B_po_dS_po_dS0[1] = (y * y + z * z)/ (R * R * R);
    arr_d2B_po_dS_po_dS0[2] = - z * (R * R * r * r - 2. * x * x * r * r - R * R * x * x)
            /(R * R * R * R * r * r * r);

    arr_d2B_po_dS_po_dS0[3] = ( x * x - y * y)/( r * r * r* r);
    arr_d2B_po_dS_po_dS0[4] = - x * y /(R * R * R);
    arr_d2B_po_dS_po_dS0[5] =  x * y * z *(2.* r * r + R * R)
                /(R * R * R * R * r * r * r);

    arr_d2B_po_dS_po_dS0[6] = 0.;
    arr_d2B_po_dS_po_dS0[7] = - x * z / (R * R * R);
    arr_d2B_po_dS_po_dS0[8] = x * (R * R - 2. * r * r)/(R * R * R * R * r);

}

//-------------------------------------------------
// вычисление частной призводной матрицы Якобиана
// преобразования прямоугольных координат в чферические по переменной
// Y (второй переменной)
void calc_d2B_po_ds_po_dS1(double *arrS, double *arr_d2B_po_dS_po_dS1)
{

    double r = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
    double R = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]+ arrS[2] * arrS[2]);
    double x =arrS[0];
    double y =arrS[1];
    double z =arrS[2];
    arr_d2B_po_dS_po_dS1[0] = ( x * x - y * y)/( r * r * r* r); //
    arr_d2B_po_dS_po_dS1[1] = - x * y /(R * R * R); //
    arr_d2B_po_dS_po_dS1[2] = x * y * z *(2.* r * r + R * R)
            /(R * R * R * R * r * r * r); //

    arr_d2B_po_dS_po_dS1[3] = 2. * x * y / (r * r * r * r); //
    arr_d2B_po_dS_po_dS1[4] = (x * x + z * z)/ (R * R * R); //

    arr_d2B_po_dS_po_dS1[5] =   - z * (R * R * r * r - 2. * y * y * r * r - R * R * y * y )
            /(R * R * R * R * r * r * r);

    arr_d2B_po_dS_po_dS1[6] = 0.;//
    arr_d2B_po_dS_po_dS1[7] = - y * z / (R * R * R);
    arr_d2B_po_dS_po_dS1[8] = y * (R * R - 2. * r * r)/(R * R * R * R * r);

}

//-------------------------------------------------
// вычисление частной призводной матрицы Якобиана
// преобразования прямоугольных координат в чферические по переменной
// Z (третьей переменной)
void calc_d2B_po_ds_po_dS2(double *arrS, double *arr_d2B_po_dS_po_dS2)
{

    double r = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
    double R = sqrt(arrS[0] * arrS[0] + arrS[1] * arrS[1]+ arrS[2] * arrS[2]);
    double x =arrS[0];
    double y =arrS[1];
    double z =arrS[2];
    arr_d2B_po_dS_po_dS2[0] = 0.;//
    arr_d2B_po_dS_po_dS2[1] = - x * z /(R * R * R); //
    arr_d2B_po_dS_po_dS2[2] = - x / r * (x * x + y * y - z * z)
            /(R * R * R * R ); //

    arr_d2B_po_dS_po_dS2[3] = 0.;
    arr_d2B_po_dS_po_dS2[4] = - y * z /(R * R * R); //

    arr_d2B_po_dS_po_dS2[5] =   - y / r * (x * x + y * y - z * z)
            /(R * R * R * R ); //

    arr_d2B_po_dS_po_dS2[6] = 0.;//
    arr_d2B_po_dS_po_dS2[7] =  (x * x + y * y)/ (R * R * R); //
    arr_d2B_po_dS_po_dS2[8] = z * (R * R - 2. * r * r)/(R * R * R * R * r);

}
//---------------------------------

void calc_dM_psk_to_kgsk_po_dQ(double *arrMu,double *arr_dM_po_dBet)
{
        const double valBet  = arrMu[0] ;
        const double valEps = arrMu[1] ;
        const double valAlf = arrMu[2] ;
        double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

        calcMatr_dP3_po_da(valBet, matrP3Q);
        calcMatrP1(valEps, matrP1Psi) ;
        calcMatrP2(valAlf, matrP2Tet) ;
        MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
        MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr_dM_po_dBet);

}


void calc_dM_psk_to_kgsk_po_dPsi(double *arrMu,double *arr_dM_po_dEps)
{

    const double valBet  = arrMu[0] ;
    const double valEps = arrMu[1] ;
    const double valAlf = arrMu[2] ;
    double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

    calcMatrP3(valBet, matrP3Q);
    calcMatr_dP1_po_da(valEps, matrP1Psi) ;
    calcMatrP2(valAlf, matrP2Tet) ;
    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr_dM_po_dEps);
}

void calc_dM_psk_to_kgsk_po_dTet(double *arrMu,double *arr_dM_po_dAlf)
{
    const double valBet  = arrMu[0] ;
    const double valEps = arrMu[1] ;
    const double valAlf = arrMu[2] ;
    double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};


    calcMatrP3(valBet, matrP3Q);
    calcMatrP1(valEps, matrP1Psi) ;
    calcMatr_dP2_po_da(valAlf, matrP2Tet) ;

    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr_dM_po_dAlf);
}




void calc_dMtrx3_ASPK_v_PSK_po_dBet(double *arrMu,double *arr_dM_po_dBet)
{
        const double valBet  = arrMu[0] ;
        const double valEps = -arrMu[1] ;
        const double valAlf = arrMu[2] ;
        double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

        calcMatr_dP3_po_da(valBet, matrP3Q);
        calcMatrP1(valEps, matrP1Psi) ;
        calcMatrP2(valAlf, matrP2Tet) ;
        MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
        MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr_dM_po_dBet);

}


void calc_dMtrx3_ASPK_v_PSK_po_dEps(double *arrMu,double *arr_dM_po_dEps)
{

    const double valBet  = arrMu[0] ;
    const double valEps = arrMu[1] ;
    const double valAlf = arrMu[2] ;
    double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

    calcMatrP3(valBet, matrP3Q);
    calcMatr_dP1_otMinus_a_po_da(valEps, matrP1Psi) ;
    calcMatrP2(valAlf, matrP2Tet) ;
    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr_dM_po_dEps);
}

void calc_dMtrx3_ASPK_v_PSK_po_dAlf(double *arrMu,double *arr_dM_po_dAlf)
{
    const double valBet  = arrMu[0] ;
    const double valEps = -arrMu[1] ;
    const double valAlf = arrMu[2] ;
    double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

    calcMatrP3(valBet, matrP3Q);
    calcMatrP1(valEps, matrP1Psi) ;
    calcMatr_dP2_po_da(valAlf, matrP2Tet) ;

    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr_dM_po_dAlf);
}

//------------------------------------------------------------
// матрица перехода из АСПК в ПСК
//arrMu[3] - углы повороитов
// Betta, Eps, Alf
void calcMtrx3_ASPK_v_PSK(double *arrMu,double *matrPereh_PSK_V_KGSK)
{
    const double valQ   = arrMu[0] ;
    const double valPsi = -arrMu[1] ;// направление поворота по углу места !!!!!
    const double valTet = arrMu[2] ;

    double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

    calcMatrP3(valQ, matrP3Q) ;
    calcMatrP1(valPsi, matrP1Psi) ;
    calcMatrP2(valTet, matrP2Tet) ;
    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, matrPereh_PSK_V_KGSK);

}

//-----------------------------------------------------------------------------
//------------------------------------------------------------------------------


//---------------------------------------------------------------------------
void calcMatrP1(const long double  a, long double *matrP)
{
  memset(matrP,0,9*sizeof(long double ));
  matrP[0] = 1;
  matrP[4] = cosl(a);
  matrP[5] = sinl(a);
  matrP[7] = -sinl(a);
  matrP[8] = cosl(a);

}
void calcMatrP2(const long double  a, long double *matrP)
{
  memset(matrP,0,9*sizeof(long double ));
  matrP[4] = 1;
  matrP[0] = cosl(a);
  matrP[2] = sinl(a);
  matrP[6] = -sinl(a);
  matrP[8] = cosl(a);

}
void calcMatrP3(const long double  a, long double *matrP)
{
  memset(matrP,0,9*sizeof(long double ));
  matrP[8] = 1;
  matrP[0] = cosl(a);
  matrP[1] = sinl(a);
  matrP[3] = -sinl(a);
  matrP[4] = cosl(a);
}

//---------------------------------------------------------------------------
void calcMatr_dP1_po_da(const long double  a, long double *matrP)
{
  memset(matrP,0,9*sizeof(long double ));

  matrP[4] = -sinl(a);
  matrP[5] = cosl(a);
  matrP[7] = -cosl(a);
  matrP[8] = -sinl(a);

}
//---------------------------------------------------------------------------
void calcMatr_d2P1_po_da2(const long double  a, long double *matrP)
{
  memset(matrP,0,9*sizeof(long double ));

  matrP[4] = -cosl(a);
  matrP[5] = -sinl(a);
  matrP[7] = sinl(a);
  matrP[8] = -cosl(a);

}
//---------------------------------------------------------------------------
void calcMatr_d2P1otMinus_a_po_da2(const long double  a, long double *matrP)
{
  memset(matrP,0,9*sizeof(long double ));

  matrP[4] = -cosl(a);
  matrP[5] = sinl(a);
  matrP[7] = -sinl(a);
  matrP[8] = -cosl(a);


}
//---------------------------------------------------------------------------
void calcMatr_dP2_po_da(const long double  a, long double *matrP)
{
  memset(matrP,0,9*sizeof(long double ));

  matrP[0] = -sinl(a);
  matrP[2] = cosl(a);
  matrP[6] = -cosl(a);
  matrP[8] = -sinl(a);

}

//---------------------------------------------------------------------------
void calcMatr_dP1_otMinus_a_po_da(const long double  a, long double *matrP)
{
    memset(matrP,0,9*sizeof(long double ));

    matrP[4] = -sinl(a);
    matrP[5] = -cosl(a);
    matrP[7] = cosl(a);
    matrP[8] = -sinl(a);

}
//---------------------------------------------------------------------------
void calcMatr_dP3_po_da(const long double  a, long double *matrP)
{
  memset(matrP,0,9*sizeof(long double ));

  matrP[0] = -sinl(a);
  matrP[1] = cosl(a);
  matrP[3] = -cosl(a);
  matrP[4] = -sinl(a);

}



//---------------------------------------------------------------------------
void calcMatr_d2P3_po_da2(const long double  a, long double *matrP)
{
  memset(matrP,0,9*sizeof(long double ));

  matrP[0] = -cosl(a);
  matrP[1] = -sinl(a);
  matrP[3] = sinl(a);
  matrP[4] = -cosl(a);

}
//-------------------------------------


void calcMatr_ASK_v_KGSK(long double  *arrMu,long double  *matrPereh_ASK_V_KGSK)
{
    /*const long double  valQ   = arrMu[0] ;
    const long double  valPsi = arrMu[1] ;
    const long double  valTet = arrMu[2] ;
    const long double  valBet = arrMu[3] ;
    const long double  valEps = arrMu[4] ;
   long double  arr1[9]={0},arr0[9]={0},arr2[9]={0}
     ,matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0},matrP3Bet[9]={0},matrP1Eps[9]={0};

     calcMatrP3(valQ, matrP3Q) ;
     calcMatrP1(valPsi, matrP1Psi) ;
     calcMatrP2(valTet, matrP2Tet) ;
     calcMatrP3(valBet, matrP3Bet) ;
     calcMatrP1(-valEps, matrP1Eps) ;
     MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
     MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr1);
     MtrxMultMatrx( arr1,3, 3, matrP3Bet,3, arr2);
     MtrxMultMatrx( arr2,3, 3, matrP1Eps,3, matrPereh_ASK_V_KGSK); */
     long double   matrPereh_ASK_V_PSK[9] = {0.}, matrPereh_PSK_V_KGSK[9] ={0.} ;
     calcMatr_ASK_v_PSK( &arrMu[3],matrPereh_ASK_V_PSK) ;
     calcMatr_PSK_v_KGSK(arrMu,matrPereh_PSK_V_KGSK);
     MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, matrPereh_ASK_V_PSK,3, matrPereh_ASK_V_KGSK ) ;


}

void calcMatr_ASK_v_PSK(long double  *arrMu,long double  *matrPereh_ASK_V_PSK)
{

    const long double  valBet = arrMu[0] ;
    const long double  valEps = arrMu[1] ;
    long double  matrP3Bet[9]={0},matrP1Eps[9]={0};
    calcMatrP3(valBet, matrP3Bet) ;
    calcMatrP1(-valEps, matrP1Eps) ;
    MtrxMultMatrx( matrP3Bet,3, 3, matrP1Eps,3, matrPereh_ASK_V_PSK);

}

void calcMatr_PSK_v_KGSK(long double  *arrMu,long double  *matrPereh_PSK_V_KGSK)
{
    const long double  valQ   = arrMu[0] ;
    const long double  valPsi = arrMu[1] ;
    const long double  valTet = arrMu[2] ;

    long double  arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

    calcMatrP3(valQ, matrP3Q) ;
    calcMatrP1(valPsi, matrP1Psi) ;
    calcMatrP2(valTet, matrP2Tet) ;
    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, matrPereh_PSK_V_KGSK);

}


void calc_dF_po_dQ_sq(const long double  valBet,const long double  valEps,long double  *matrRez)
{
    long double  arr[3]={0},arr0[3]={0};
    arr[0] = cosl(valBet) * cosl ( valEps) ;
    arr[1] = -sinl(valBet) * cosl ( valEps) ;
    memcpy(arr0,arr,3*sizeof(long double ));
    MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);
}
void calc_dF_po_dPsi_sq(const long double  valBet,const long double  valEps,long double  *matrRez)
{
    long double  arr[3]={0},arr0[3]={0};
    arr[1] = sinl(valEps) ;
    arr[2] = -cosl(valBet) * cosl ( valEps) ;
    memcpy(arr0,arr,3*sizeof(long double ));
    MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);

}
void calc_dF_po_dTet_sq(const long double  valBet,const long double  valEps,long double  *matrRez)
{
    long double  arr[3]={0},arr0[3]={0};
    arr[0] = sinl(valEps) ;
    arr[2] = -sinl(valBet) * cosl ( valEps) ;
    memcpy(arr0,arr,3*sizeof(long double ));
    MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);

}
void calc_dF_po_dBet_sq(const long double  valBet,const long double  valEps,long double  *matrRez)
{
    long double  arr[3]={0},arr0[3]={0};
    arr[0] = cosl(valBet) * cosl ( valEps) ;
    arr[1] = -sinl(valBet) * cosl ( valEps) ;
    memcpy(arr0,arr,3*sizeof(long double ));
    MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);

}
void calc_dF_po_dEps_sq(const long double  valBet,const long double  valEps,long double  *matrRez)
{
    long double  arr[3]={0},arr0[3]={0};
    arr[0] =  -sinl(valBet) *sinl(valEps) ;
    arr[1] = -cosl(valBet) * sinl ( valEps) ;
    arr[2] =  cosl ( valEps) ;
    memcpy(arr0,arr,3*sizeof(long double ));
    MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);

}
void createExtendMtrx(long double  *arrInp, long double  * arrOut)
{
  memset(arrOut, 0, 36 * sizeof(long double ));
 for (int i =0; i < 3; i++)
 for (int j =0; j < 3; j++)
 {
  arrOut[ i * 6 + j] =  arrInp [3 * i + j];
  arrOut[ (3 + i) * 6 + 3 + j] =  arrInp [3 * i + j];
 }
}
void calcMatrJ1(const long double  valBet,const long double  valEps,long double  *arrJ1)
{
 memset(arrJ1, 0, 9 * sizeof(long double )) ;
 arrJ1[1] = cosl(valEps);
 arrJ1[2] = -sinl(valEps);
 arrJ1[3] = -cosl(valEps);
 arrJ1[6] = sinl(valEps);

}
void calcMatrJ2(const long double  valBet,const long double  valEps,long double  *arrJ2)
{
 memset(arrJ2, 0, 9 * sizeof(long double )) ;
 arrJ2[1] = -sinl(valBet)*sinl(valEps);
 arrJ2[2] = -sinl(valBet)*cosl(valEps);
 arrJ2[3] = sinl(valBet)*sinl(valEps);
 arrJ2[5] = cosl(valBet);
 arrJ2[6] = sinl(valBet)*cosl(valEps);
 arrJ2[7] = -cosl(valBet);


}
void calcMatrJ3(const long double  valBet,const long double  valEps,long double  *arrJ3)
{
 memset(arrJ3, 0, 9 * sizeof(long double )) ;
 arrJ3[1] = cosl(valBet)*sinl(valEps);
 arrJ3[2] = cosl(valBet)*cosl(valEps);
 arrJ3[3] = -cosl(valBet)*sinl(valEps);
 arrJ3[5] = sinl(valBet);
 arrJ3[6] = -cosl(valBet)*cosl(valEps);
 arrJ3[7] = -sinl(valBet);

}
void calcMatrJ4(const long double  valBet,const long double  valEps,long double  *arrJ4)
{
 calcMatrJ1( valBet, valEps,arrJ4) ;

}
void calcMatrJ5(const long double  valBet,const long double  valEps,long double  *arrJ5)
{
 memset(arrJ5, 0, 9 * sizeof(long double )) ;

 arrJ5 [ 5] = -1 ;
 arrJ5 [ 7] = 1 ;

}
void calcMatrJ1W(const long double  valBet,const long double  valEps,long double  *arrJ1W)
{
   memset(arrJ1W, 0, 36 * sizeof(long double )) ;
   long double  arr[9]={0} ;
   calcMatrJ1( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ1W);
}
void calcMatrJ2W(const long double  valBet,const long double  valEps,long double  *arrJ2W)
{
   memset(arrJ2W, 0, 36 * sizeof(long double )) ;
   long double  arr[9]={0} ;
   calcMatrJ2( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ2W);
}
void calcMatrJ3W(const long double  valBet,const long double  valEps,long double  *arrJ3W)
{
   memset(arrJ3W, 0, 36 * sizeof(long double )) ;
   long double  arr[9]={0} ;
   calcMatrJ3( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ3W);
}
void calcMatrJ4W(const long double  valBet,const long double  valEps,long double  *arrJ4W)
{
   memset(arrJ4W, 0, 36 * sizeof(long double )) ;
   long double  arr[9]={0} ;
   calcMatrJ4( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ4W);
}
void calcMatrJ5W(const long double  valBet,const long double  valEps,long double  *arrJ5W)
{
   memset(arrJ5W, 0, 36 * sizeof(long double )) ;
   long double  arr[9]={0} ;
   calcMatrJ5( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ5W);
}

//--------------------------------------------------------------------------------
// Вычисление корреляционнной матрицы ошибок , возникающих при переводе вектора из ПСК в КГСК.
// ОШибки возникают в следствие того, что палубные углы известны неточно, а с ошибками измерений.
// Задан вектор arrVS[3], измерения палубных углов и дисперсии ошибок измерения этих углов.
// Требуется определить коррел матрицу ошибок преобразованного вектора arrVS
// INPUT:
// VAlQ,  VAlPsi,  VAlTet  - измерения курсового угда, угла килевой и бортовой качек
// VAlDispQ, VAlDispPsi	 , VAlDispTet  - дисперсии ошибок измерений этих углов
// arrVS[3] - вектор в ПСК
// OUTPUT:
// arrVS_KGSK[3] - вектор в КГСК
// arrMtxCorr[9] - корреляционная матрица ошибок
void calcCorMatrx_PSK_KGSK(const long double  VAlQ, const long double  VAlPsi, const long double  VAlTet
     , const long double  VAlDispQ,const long double  VAlDispPsi	 , const long double  VAlDispTet
     , long double   *arrVS_PSK,long double   *arrVS_KGSK, long double  *arrMtxCorr)
{
    long double  arrP1[9] = {0.},arrP2[9] = {0.},arrP3[9] = {0.},
        arr_dP1[9] = {0.},arr_dP2[9] = {0.},arr_dP3[9] = {0.};
    calcMatrP1(VAlPsi, arrP1);
    calcMatr_dP1_po_da(VAlPsi, arr_dP1);

    calcMatrP2(VAlTet, arrP2);
    calcMatr_dP2_po_da(VAlTet, arr_dP2);

    calcMatrP3(VAlQ, arrP3);
    calcMatr_dP3_po_da(VAlQ, arr_dP3);

    long double  arrJ1[9] = {0.}, arrJ2[9] = {0.}, arrJ3[9] = {0.};
    MtrxMultMatrx_MultMatrx(arr_dP3, arrP1,arrP2,3, arrJ1)  ;
    MtrxMultMatrx_MultMatrx(arrP3, arr_dP1,arrP2,3, arrJ2)  ;
    MtrxMultMatrx_MultMatrx(arrP3, arrP1,arr_dP2,3, arrJ3)  ;
    ///
 memset(arrMtxCorr, 0, 9 * sizeof(long double ));
 long double  arrQ[9] = {0.}, arrPsi[9] = {0.}, arrTet[9] = {0.}, arrTemp0[9] = {0.}, arrVS_PSK_Copy[3] ={0.}, arrS_mult_ST[9] = {0.};
 memcpy(arrVS_PSK_Copy, arrVS_PSK, 3 * sizeof(long double ));
 MtrxMultMatrxTransp(arrVS_PSK,3, 1, arrVS_PSK_Copy,3, arrS_mult_ST) ;
 MtrxMultMatrx_MultMatrxTransp(arrJ1,arrS_mult_ST, arrJ1, 3, arrQ);
 MtrxMultMatrx_MultMatrxTransp(arrJ2,arrS_mult_ST, arrJ2, 3, arrPsi);
 MtrxMultMatrx_MultMatrxTransp(arrJ3,arrS_mult_ST, arrJ3, 3, arrTet);

 MatrxMultScalar(arrQ, 3, 3,VAlDispQ,arrQ) ;
 MatrxMultScalar(arrPsi, 3, 3,VAlDispPsi,arrPsi) ;
 MatrxMultScalar(arrTet, 3, 3,VAlDispTet,arrTet) ;

 MtrxSumMatrx(arrQ, arrPsi,3, 3, arrTemp0) ;
 MtrxSumMatrx(arrTemp0, arrTet,3, 3, arrMtxCorr) ;

 ///

MtrxMultMatrx_MultMatrx(arrP3, arrP1,arrP2,3, arrJ1)  ;
MtrxMultMatrx(arrJ1,3, 3, arrVS_PSK, 1, arrVS_KGSK) ;
}

// матрица перехода parrMtrxPer из исходной прямоугольной сиситемы координат
// в скоростную сиситему координат (CCК), задаваемую вектором arrV
// ось X  ССК направлена по вектору V, ось Y CCК праллельна плоскости OXY
// исходной прямоугольной сиситемы координат, а ось Z дополняет до правой тройки
// parrMtrxPer - это матрица состоящая из столбцов вектров единичных ортов
void calcMatrxPer_from_DecartPrSK_To_SSK(long double  *arrV, long double  * parrMtrxPer)
{
    long double  val_v0 = NormVect2(arrV);
    memset(parrMtrxPer, 0, 9 * sizeof(long double ));
    if (val_v0 < 0.0000000001)
  {
    parrMtrxPer[2] = 1.;
    parrMtrxPer[3] = 1.;
    parrMtrxPer[7] = 1.;
    return;
  }
  long double  arr_v[3] = {0.};
  memcpy(arr_v, arrV, 3 *  sizeof(long double ));
    NormalizeVect3(arr_v) ;
    long double  val_v = NormVect2(arr_v);
  parrMtrxPer[0] = arr_v [0];
  parrMtrxPer[1] = arr_v [1];
    parrMtrxPer[2] = arr_v[2];
    parrMtrxPer[3] = -arr_v [1] / val_v;
    parrMtrxPer[4] = arr_v [0] / val_v;
  parrMtrxPer[5] = 0.;
    parrMtrxPer[6] = -arr_v [2] * arr_v [0] / val_v;
    parrMtrxPer[7] = -arr_v [2] * arr_v [1] / val_v;
  parrMtrxPer[8] =  val_v ;
  return;
}


// матрица перехода parrMtrxPer из скоростной  сиситемы координат (CCК), задаваемой вектором arrV
// в исходную  прямоугольную  сиситему координат
// ось X  ССК направлена по вектору V, ось Y CCК праллельна плоскости OXY
// исходной прямоугольной сиситемы координат, а ось Z дополняет до правой тройки
// parrMtrxPer - это матрица состоящая из столбцов вектров единичных ортов
// CCК  в исходной прямоуг сиситемы координат
void calcMatrxPer_from_SSK_To_DecartPrSK(long double  *arrV, long double  * parrMtrxPer)
{
  long double  arrT[9] = {0.};
  calcMatrxPer_from_DecartPrSK_To_SSK(arrV, arrT) ;
  MatrTransp(arrT, 3, 3, parrMtrxPer);
  return;
}

// Пенресчет прямоугольных координат в сферические
//  arrInp[3] - прямоугольные координаты
// valBet, valR,valEps - сферические координаты
void recalcCoord_INTO_Spherical(long double  *arrInp, long double  &valR, long double  &valBet, long double  &valEps)
{
    valR = sqrtl ( arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1] +arrInp[2] * arrInp[2] ) ;

    valEps  = asinl(arrInp[2] / valR ) ;

    valBet  = asinl(arrInp[0]/ sqrtl (  arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1])) ;
    if (arrInp[1] < 0.)
     {
        if (arrInp[0] < 0.)
        {
          valBet = -M_PI - valBet ;
        }
        else
        {
           valBet = M_PI - valBet ;
        }
     }


}
//-------------------------------------------
// Пенресчет сферических  координат в прямоугольные
//  arrS[3] - прямоугольные координаты
// valBet, valR,valEps - сферические координаты
void recalcSphericalCoord_INTO_Rectangular(const long double  valR,const  long double  valBet,const  long double  valEps
                                           , long double  *arrS)
{
   arrS[0] =  valR * cosl(valEps)* sinl(valBet);
   arrS[1] =  valR * cosl(valEps)* cosl(valBet);
   arrS[2] =  valR * sinl(valEps);
}
//-------------------------------------------------------
//Вычисление матрицы частных производных векто-функции сферических кооординат
// по прямоугольным. порядок следования сферических координат -
// курс угол, дальность, угол места
void calc_dSpherical_po_Rectangular(long double  *arrS, long double  *arr_dB_po_dS)
{
    long double  r = sqrtl(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
    long double  R = sqrtl(arrS[0] * arrS[0] + arrS[1] * arrS[1]+ arrS[2] * arrS[2]);

    arr_dB_po_dS[0] = arrS[1]/r/r;
    arr_dB_po_dS[1] = -arrS[0]/r/r;
    arr_dB_po_dS[2] = 0.;

    arr_dB_po_dS[3] = arrS[0]/R;
    arr_dB_po_dS[4] = arrS[1]/R;
    arr_dB_po_dS[5] = arrS[2]/R;

    arr_dB_po_dS[6] = -arrS[0]* arrS[2]/R/R/r;
    arr_dB_po_dS[7] = -arrS[1]* arrS[2]/R/R/r;
    arr_dB_po_dS[8] = r/R/R;
}
//-------------------------------------------------------
//Вычисление матрицы частных производных векто-функции угловых сферических кооординат
// по прямоугольным. порядок следования сферических координат -
// курс угол,  угол места
void calc_dSpherical_po_Rectangular_AnglesOnly(long double  *arrS, long double  *arr_dB_po_dS)
{
    long double  r = sqrtl(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
    long double  R = sqrtl(arrS[0] * arrS[0] + arrS[1] * arrS[1]+ arrS[2] * arrS[2]);

    arr_dB_po_dS[0] = arrS[1]/r/r;
    arr_dB_po_dS[1] = -arrS[0]/r/r;
    arr_dB_po_dS[2] = 0.;

    arr_dB_po_dS[3] = -arrS[0]* arrS[2]/R/R/r;
    arr_dB_po_dS[4] = -arrS[1]* arrS[2]/R/R/r;
    arr_dB_po_dS[5] = r/R/R;

}
//-------------------------------------------------------
//Вычисление матрицы  производных матрицы перехода их АСК в ПСК
// по углу Betta (курсовому)
void calc_dM_ask_to_psk_po_Betta(const long double  valBet,const long double  valEps,long double  *arr_dM_po_dBet)
{
    memset(arr_dM_po_dBet, 0, 9 * sizeof(long double ));
    arr_dM_po_dBet[0] = -sinl(valBet);
    arr_dM_po_dBet[1] = cosl(valBet) * cosl (valEps);
    arr_dM_po_dBet[2] = -cosl(valBet) * sinl (valEps);

    arr_dM_po_dBet[3] = -cosl(valBet);
    arr_dM_po_dBet[4] = -sinl(valBet) * cosl (valEps);
    arr_dM_po_dBet[5] = sinl(valBet) * sinl (valEps);
}
//-------------------------------------------------------
//Вычисление матрицы  производных матрицы перехода их АСК в ПСК
// по углу Eps (места)
void calc_dM_ask_to_psk_po_Eps(const long double  valBet,const long double  valEps,long double  *arr_dM_po_dEps)
{
    memset(arr_dM_po_dEps, 0, 9 * sizeof(long double ));

    arr_dM_po_dEps[1] = -sinl(valBet) * sinl (valEps);
    arr_dM_po_dEps[2] = -sinl(valBet) * cosl (valEps);

    arr_dM_po_dEps[4] = -cosl(valBet) * sinl (valEps);
    arr_dM_po_dEps[5] = -cosl(valBet) * cosl (valEps);

    arr_dM_po_dEps[7] =  cosl (valEps);
    arr_dM_po_dEps[8] =  -sinl (valEps);
}

//-------------------------------------------------------
//Вычисление матрицы частных производных вектор-функции прямоугольных кооординат
// по сферическим.
// порядок следования сферических координат -
//arrV[3]- курс угол, дальность, угол места
void calc_dRectangularl_po_dSpherical(const long double  *arrV, long double  *arr_dA_po_dV)
{
    memset(arr_dA_po_dV, 0, 9 * sizeof(long double ));
    long double  e = arrV[2];
    long double  q = arrV[0];
    long double  r = arrV[1];

    arr_dA_po_dV[0] = r * cosl(q) * cosl (e);
    arr_dA_po_dV[1] = sinl(q) * cosl (e);
    arr_dA_po_dV[2] = -r * sinl(q) * sinl (e);

    arr_dA_po_dV[3] = -r * sinl(q) * cosl (e);
    arr_dA_po_dV[4] = cosl(q) * cosl (e);
    arr_dA_po_dV[5] = -r * cosl(q) * sinl (e);
    arr_dA_po_dV[6] = 0;
    arr_dA_po_dV[7] = sinl(e);
    arr_dA_po_dV[8] = r * cosl(e);
}
//-------------------------------------------------
// вычисление частной призводной матрицы Якобиана
// преобразования прямоугольных координат в чферические по переменной
// X (первой переменной)
void calc_d2B_po_ds_po_dS0(long double  *arrS, long double  *arr_d2B_po_dS_po_dS0)
{

    long double  r = sqrtl(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
    long double  R = sqrtl(arrS[0] * arrS[0] + arrS[1] * arrS[1]+ arrS[2] * arrS[2]);
    long double  x =arrS[0];
    long double  y =arrS[1];
    long double  z =arrS[2];
    arr_d2B_po_dS_po_dS0[0] = -2. * x * y / (r * r * r * r);
    arr_d2B_po_dS_po_dS0[1] = (y * y + z * z)/ (R * R * R);
    arr_d2B_po_dS_po_dS0[2] = - z * (R * R * r * r - 2. * x * x * r * r - R * R * x * x)
            /(R * R * R * R * r * r * r);

    arr_d2B_po_dS_po_dS0[3] = ( x * x - y * y)/( r * r * r* r);
    arr_d2B_po_dS_po_dS0[4] = - x * y /(R * R * R);
    arr_d2B_po_dS_po_dS0[5] =  x * y * z *(2.* r * r + R * R)
                /(R * R * R * R * r * r * r);

    arr_d2B_po_dS_po_dS0[6] = 0.;
    arr_d2B_po_dS_po_dS0[7] = - x * z / (R * R * R);
    arr_d2B_po_dS_po_dS0[8] = x * (R * R - 2. * r * r)/(R * R * R * R * r);

}

//-------------------------------------------------
// вычисление частной призводной матрицы Якобиана
// преобразования прямоугольных координат в чферические по переменной
// Y (второй переменной)
void calc_d2B_po_ds_po_dS1(long double  *arrS, long double  *arr_d2B_po_dS_po_dS1)
{

    long double  r = sqrtl(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
    long double  R = sqrtl(arrS[0] * arrS[0] + arrS[1] * arrS[1]+ arrS[2] * arrS[2]);
    long double  x =arrS[0];
    long double  y =arrS[1];
    long double  z =arrS[2];
    arr_d2B_po_dS_po_dS1[0] = ( x * x - y * y)/( r * r * r* r); //
    arr_d2B_po_dS_po_dS1[1] = - x * y /(R * R * R); //
    arr_d2B_po_dS_po_dS1[2] = x * y * z *(2.* r * r + R * R)
            /(R * R * R * R * r * r * r); //

    arr_d2B_po_dS_po_dS1[3] = 2. * x * y / (r * r * r * r); //
    arr_d2B_po_dS_po_dS1[4] = (x * x + z * z)/ (R * R * R); //

    arr_d2B_po_dS_po_dS1[5] =   - z * (R * R * r * r - 2. * y * y * r * r - R * R * y * y )
            /(R * R * R * R * r * r * r);

    arr_d2B_po_dS_po_dS1[6] = 0.;//
    arr_d2B_po_dS_po_dS1[7] = - y * z / (R * R * R);
    arr_d2B_po_dS_po_dS1[8] = y * (R * R - 2. * r * r)/(R * R * R * R * r);

}

//-------------------------------------------------
// вычисление частной призводной матрицы Якобиана
// преобразования прямоугольных координат в чферические по переменной
// Z (третьей переменной)
void calc_d2B_po_ds_po_dS2(long double  *arrS, long double  *arr_d2B_po_dS_po_dS2)
{

    long double  r = sqrtl(arrS[0] * arrS[0] + arrS[1] * arrS[1]);
    long double  R = sqrtl(arrS[0] * arrS[0] + arrS[1] * arrS[1]+ arrS[2] * arrS[2]);
    long double  x =arrS[0];
    long double  y =arrS[1];
    long double  z =arrS[2];
    arr_d2B_po_dS_po_dS2[0] = 0.;//
    arr_d2B_po_dS_po_dS2[1] = - x * z /(R * R * R); //
    arr_d2B_po_dS_po_dS2[2] = - x / r * (x * x + y * y - z * z)
            /(R * R * R * R ); //

    arr_d2B_po_dS_po_dS2[3] = 0.;
    arr_d2B_po_dS_po_dS2[4] = - y * z /(R * R * R); //

    arr_d2B_po_dS_po_dS2[5] =   - y / r * (x * x + y * y - z * z)
            /(R * R * R * R ); //

    arr_d2B_po_dS_po_dS2[6] = 0.;//
    arr_d2B_po_dS_po_dS2[7] =  (x * x + y * y)/ (R * R * R); //
    arr_d2B_po_dS_po_dS2[8] = z * (R * R - 2. * r * r)/(R * R * R * R * r);

}
//---------------------------------
void calc_dM_psk_to_kgsk_po_dBetta(long double  *arrMu,long double  *arr_dM_po_dBet)
{
        const long double  valBet  = arrMu[0] ;
        const long double  valEps = arrMu[1] ;
        const long double  valAlf = arrMu[2] ;
        long double  arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

        calcMatr_dP3_po_da(valBet, matrP3Q);
        calcMatrP1(valEps, matrP1Psi) ;
        calcMatrP2(valAlf, matrP2Tet) ;
        MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
        MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr_dM_po_dBet);

}

void calc_dM_psk_to_kgsk_po_dEps(long double  *arrMu,long double  *arr_dM_po_dEps)
{

    const long double  valBet  = arrMu[0] ;
    const long double  valEps = arrMu[1] ;
    const long double  valAlf = arrMu[2] ;
    long double  arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

    calcMatrP3(valBet, matrP3Q);
    calcMatr_dP1_po_da(valEps, matrP1Psi) ;
    calcMatrP2(valAlf, matrP2Tet) ;
    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr_dM_po_dEps);
}

void calc_dM_psk_to_kgsk_po_dAlf(long double  *arrMu,long double  *arr_dM_po_dAlf)
{
    const long double  valBet  = arrMu[0] ;
    const long double  valEps = arrMu[1] ;
    const long double  valAlf = arrMu[2] ;
    long double  arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

    calcMatr_dP3_po_da(valBet, matrP3Q);
    calcMatrP1(valEps, matrP1Psi) ;
    calcMatr_dP2_po_da(valAlf, matrP2Tet) ;

    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr_dM_po_dAlf);
}

//------------------------------------------------------------
// матрица перехода из АСПК в ПСК
//arrMu[3] - углы повороитов
// Betta, Eps, Alf
void calcMtrx3_ASPK_v_PSK(long double  *arrMu,long double  *matrPereh_PSK_V_KGSK)
{
    const long double  valQ   = arrMu[0] ;
    const long double  valPsi = -arrMu[1] ;// направление поворота по углу места !!!!!
    const long double  valTet = arrMu[2] ;

    long double  arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

    calcMatrP3(valQ, matrP3Q) ;
    calcMatrP1(valPsi, matrP1Psi) ;
    calcMatrP2(valTet, matrP2Tet) ;
    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, matrPereh_PSK_V_KGSK);

}




