#include "CoordSystTrsf.h"
#include <string.h>
#include <math.h>
#include "MatrixProccess.h"
#include <QPolygonF>

QCoordSystTrsf::QCoordSystTrsf()
{

}
//---------------------------------------------------------------------------

// заданы 2 правых системы координат СК1 и СК2.
// система координат СК2 получается из СК1 путем поворота на угол а
// против часовой стрелки относительно оси OX СК1
// тогда matrP - матрица перехода из СК2 в СК1,
// то есть, если S2 - вектор координат точки в СК2, а S1 вектор координат той же точки в СК2,
// то они связаны соотношением
// S1 = matrP * S2
void QCoordSystTrsf::calcMatrP1(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));
  matrP[0] = 1;
  matrP[4] = cos(a);
  matrP[5] = -sin(a);
  matrP[7] = sin(a);
  matrP[8] = cos(a);

}
//---------------------------------------------------------------------------

// заданы 2 правых системы координат СК1 и СК2.
// система координат СК2 получается из СК1 путем поворота на угол а
// против часовой стрелки относительно оси OY СК1
// тогда matrP - матрица перехода из СК2 в СК1,
// то есть, если S2 - вектор координат точки в СК2, а S1 вектор координат той же точки в СК2,
// то они связаны соотношением
// S1 = matrP * S2
void QCoordSystTrsf::calcMatrP2(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));
  matrP[4] = 1;
  matrP[0] = cos(a);
  matrP[2] = sin(a);
  matrP[6] = -sin(a);
  matrP[8] = cos(a);

}
//---------------------------------------------------------------------------

// заданы 2 правых системы координат СК1 и СК2.
// система координат СК2 получается из СК1 путем поворота на угол а
// против часовой стрелки относительно оси OZ СК1
// тогда matrP - матрица перехода из СК2 в СК1,
// то есть, если S2 - вектор координат точки в СК2, а S1 вектор координат той же точки в СК2,
// то они связаны соотношением
// S1 = matrP * S2
void QCoordSystTrsf::calcMatrP3(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));
  matrP[8] = 1;
  matrP[0] = cos(a);
  matrP[1] = -sin(a);
  matrP[3] = sin(a);
  matrP[4] = cos(a);
}

//---------------------------------------------------------------------------
void QCoordSystTrsf::calcMatr_dP1_po_da(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));

  matrP[4] = -sin(a);
  matrP[5] = -cos(a);
  matrP[7] = cos(a);
  matrP[8] = -sin(a);

}
//---------------------------------------------------------------------------

void QCoordSystTrsf::calcMatr_dP2_po_da(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));

  matrP[0] = -sin(a);
  matrP[2] = cos(a);
  matrP[6] = -cos(a);
  matrP[8] = -sin(a);

}

//---------------------------------------------------------------------------
void QCoordSystTrsf::calcMatr_dP3_po_da(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));

  matrP[0] = -sin(a);
  matrP[1] = -cos(a);
  matrP[3] = cos(a);
  matrP[4] = -sin(a);

}
//--------------------------------------------------------------
// задана система координат КГСК и ПСК
// ПСК получается из КГСК путем осуществления 3-х поворотов
// против часовой стрелки (левые повороты) на углы
// arrMu[0] - курс
// arrMu[1] - дифферент (кил качка)
// arrMu[2] - крен (бортовая качка)
// это по правилам теор механики!!!
// matrPereh_PSK_V_KGSK - матрица перехода из ПСК в КГСК
//иными словами, S1 - вектор координат точки в КГСК
//S2 - вектор координат точки в ПСК
// Тогда:
// S1 = matrPereh_PSK_V_KGSK * S2
void QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot(double *arrMu,double *matrPereh_PSK_V_KGSK)
{
    const double valQ   = arrMu[0] ;
    const double valPsi = arrMu[1] ;
    const double valTet = arrMu[2] ;

    double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};
  // поворот против часовой стрелки
    calcMatrP3(valQ, matrP3Q) ;
  // поворот против часовой стрелки
    calcMatrP1(valPsi, matrP1Psi) ;
   // поворот против часовой стрелки
    calcMatrP2(valTet, matrP2Tet) ;
    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, matrPereh_PSK_V_KGSK);

}

//--------------------------------------------------------------
// задана система координат КГСК и ПСК
// ПСК получается из КГСК путем осуществления 3-х поворотов
// по часовой стрелке (праввые повороты) на углы
// arrMu[0] - курс
// arrMu[1] - диффенрент (кил качка)
// arrMu[2] - крен (бортовая качка)
// так принято в морском деле !!!
// matrPereh_PSK_V_KGSK - матрица перехода из ПСК в КГСК
//иными словами, S1 - вектор координат точки в КГСК
//S2 - вектор координат точки в ПСК
// Тогда:
// S1 = matrPereh_PSK_V_KGSK * S2
void QCoordSystTrsf::calcMatr_PSK_v_KGSK_RightRot(double *arrMu,double *matrPereh_PSK_V_KGSK)
{
    double arrMu1[3] = {0.};
    MatrxMultScalar(arrMu, 1, 3, -1.,arrMu1);
    calcMatr_PSK_v_KGSK_LeftRot(arrMu1,matrPereh_PSK_V_KGSK);

}
// Пенресчет прямоугольных координат в сферические
//  arrInp[3] - прямоугольные координаты
// valBet, valR,valEps - сферические координаты
void QCoordSystTrsf::recalcCoord_INTO_Spherical(double *arrInp, double &valR, double &valBet, double &valEps)
{
    valR = sqrt ( arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1] +arrInp[2] * arrInp[2] ) ;

    valEps  = asin(arrInp[2] / valR ) ;

    valBet  = asin(arrInp[0]/ sqrt (  arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1])) ;
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
void QCoordSystTrsf::recalcSphericalCoord_INTO_Rectangular(const double valR,const  double valBet,const  double valEps
                                           , double *arrS)
{
   arrS[0] =  valR * cos(valEps)* sin(valBet);
   arrS[1] =  valR * cos(valEps)* cos(valBet);
   arrS[2] =  valR * sin(valEps);
}
//-------------------------------------------------------
//Вычисление матрицы частных производных вектор-функции сферических кооординат
// по прямоугольным. порядок следования сферических координат -
// курс угол, дальность, угол места
void QCoordSystTrsf::calc_dSpherical_po_Rectangular(double *arrS, double *arr_dB_po_dS)
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
void QCoordSystTrsf::calc_dSpherical_po_Rectangular_MultCoeff(double valCoeff, double *arrS, double *arr_dB_po_dS_MultC)
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
void QCoordSystTrsf::calc_dSpherical_po_Rectangular_AnglesOnly(double *arrS, double *arr_dB_po_dS)
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
//Вычисление матрицы частных производных вектор-функции угловых сферических кооординат
// по прямоугольным. порядок следования сферических координат -
// курс угол,  угол места
void QCoordSystTrsf::calc_dSpherical_po_Rectangular_AnglesOnly_MultCoeff(double valCoeff, double *arrS, double *arr_dB_po_dS_MultC)
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

//------------------------------------------------------------
// матрица перехода из АСПК в ПСК
//arrMu[3] - углы повороитов
// АСПК получается из ПСК путем омуществления 3-х поворотов:
// поворот относительно оси OZ ПСК на угол arrMu[0]по часовой стрелке
// поворот относительно оси OX ПСК на угол  arrMu[1] против часовой стрелки
// поворот относительно оси OY ПСК на угол arrMu[2]по часовой стрелке
void QCoordSystTrsf::calcMtrx3_ASPK_v_PSK(double *arrMu,double *matrPereh_PSK_V_KGSK)
{
    const double valQ   = -arrMu[0] ;
    const double valPsi = arrMu[1] ;// направление поворота по углу места !!!!!
    const double valTet = -arrMu[2] ;

    double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

    calcMatrP3(valQ, matrP3Q) ;
    calcMatrP1(valPsi, matrP1Psi) ;
    calcMatrP2(valTet, matrP2Tet) ;
    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, matrPereh_PSK_V_KGSK);

}
//-------------------------------------------------------
//Вычисление матрицы частных производных вектор-функции прямоугольных кооординат
// по сферическим.
// порядок следования сферических координат -
//arrV[3]- курс угол, дальность, угол места
void QCoordSystTrsf::calc_dRectangularl_po_dSpherical(const double *arrV, double *arr_dA_po_dV)
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

void QCoordSystTrsf::calc_dMLeft_psk_to_kgsk_po_dQ(double *arrMu,double *arr_dM_po_dBet)
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


void QCoordSystTrsf::calc_dMLeft_psk_to_kgsk_po_dPsi(double *arrMu,double *arr_dM_po_dEps)
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

void QCoordSystTrsf::calc_dMLeft_psk_to_kgsk_po_dTet(double *arrMu,double *arr_dM_po_dAlf)
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
//---------------------------------------------

void QCoordSystTrsf::calc_dMtrx3_ASPK_v_PSK_po_dQ(double *arrMu,double *arr_dM_po_dQ)
{
    const double valQ   = -arrMu[0] ;
    const double valPsi = arrMu[1] ;// направление поворота по углу места !!!!!
    const double valTet = -arrMu[2] ;

    double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

    calcMatr_dP3_po_da(valQ, matrP3Q);
    calcMatrP1(valPsi, matrP1Psi) ;
    calcMatrP2(valTet, matrP2Tet) ;
    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr_dM_po_dQ);
    MatrxMultScalar(arr_dM_po_dQ, 3, 3, -1.,arr_dM_po_dQ);

}
//---------------------------------------------
void QCoordSystTrsf::calc_dMtrx3_ASPK_v_PSK_po_dPsi(double *arrMu,double *arr_dM_po_dPsi)
{
    const double valQ   = -arrMu[0] ;
    const double valPsi = arrMu[1] ;// направление поворота по углу места !!!!!
    const double valTet = -arrMu[2] ;

    double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};



    calcMatrP3(valQ, matrP3Q);
    calcMatr_dP1_po_da(valPsi, matrP1Psi) ;
    calcMatrP2(valTet, matrP2Tet) ;
    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr_dM_po_dPsi);

}
//---------------------------------------------
void QCoordSystTrsf::calc_dMtrx3_ASPK_v_PSK_po_dTet(double *arrMu,double *arr_dM_po_dTet)
{
    const double valQ   = -arrMu[0] ;
    const double valPsi = arrMu[1] ;// направление поворота по углу места !!!!!
    const double valTet = -arrMu[2] ;

    double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

    calcMatrP3(valQ, matrP3Q);
    calcMatrP1(valPsi, matrP1Psi) ;
    calcMatr_dP2_po_da(valTet, matrP2Tet) ;

    MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
    MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr_dM_po_dTet);
    MatrxMultScalar(arr_dM_po_dTet, 3, 3, -1.,arr_dM_po_dTet);


}
//-----------------------------------
// матрица перехода parrMtrxPer из исходной прямоугольной сиситемы координат
// в скоростную сиситему координат (CCК), задаваемую вектором arrV
// ось X  ССК направлена по вектору V, ось Y CCК праллельна плоскости OXY
// исходной прямоугольной сиситемы координат, а ось Z дополняет до правой тройки
// parrMtrxPer - это матрица состоящая из столбцов вектров единичных ортов
void QCoordSystTrsf::calcMatrxPer_from_DecartPrSK_To_SSK(double *arrV, double * parrMtrxPer)
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

//------------------------------------------
// матрица перехода parrMtrxPer из скоростной  сиситемы координат (CCК), задаваемой вектором arrV
// в исходную  прямоугольную  сиситему координат
// ось X  ССК направлена по вектору V, ось Y CCК праллельна плоскости OXY
// исходной прямоугольной сиситемы координат, а ось Z дополняет до правой тройки
// parrMtrxPer - это матрица состоящая из столбцов вектров единичных ортов
// CCК  в исходной прямоуг сиситемы координат
void QCoordSystTrsf::calcMatrxPer_from_SSK_To_DecartPrSK(double *arrV, double * parrMtrxPer)
{
  double arrT[9] = {0.};
  calcMatrxPer_from_DecartPrSK_To_SSK(arrV, arrT) ;
  MatrTransp(arrT, 3, 3, parrMtrxPer);
  return;
}
//
//-----------------------------------
// матрица перехода parrMtrxPer из скоростной прямоугольной сиситемы координат
// в нормальную систему координат ГОСТ 20058-80
// Скоростная СК это подвижная СК ось OX которой совпадает с направлением скорости ЛА
// ось OY лежит в плоскости симметрии ЛА и направлена к верхней части ЛА
// psia - скоростной угол рыскания
// это угол между осью OXg ноормальной СК и проекцией скорости на горизонтальную
// плоскость OXgZg нормальной СК
// СКоростной угол рыскания следует считать положительным, когда ось OXg
//совмещается  с проекцией скорости на горизонт плоскость OXgZg
 // поворотом вокруг оси OYg против часовой стрелки
// nua - скоростной угол тангажа
// это угол между скоростной осью и горизоньтальной плоскостью
// OxgZg нормальной СК
// gamma_a - скоростной угол крена
// это угол между боковой осью OZa скоростнгой СК и осью OZg нормальной СК
void QCoordSystTrsf::calcMatrxPer_from_SkSK_To_NormSK(const double psia ,
    const double nua ,const double gamma_a  ,double* parrMtrxPer)
{
    parrMtrxPer[0] =  cos (psia) * cos (nua);

    parrMtrxPer[1] = -cos (psia) * sin (nua) * cos (gamma_a)
            + sin (psia) * sin(gamma_a);

    parrMtrxPer[2] = cos (psia) * sin (nua) * sin (gamma_a)
            + sin (psia) * cos(gamma_a);
    ///
    parrMtrxPer[3] =  sin (nua);

    parrMtrxPer[4] = cos (nua)* cos(gamma_a);

    parrMtrxPer[5] = - cos (nua)* sin(gamma_a);
    ///
    parrMtrxPer[6] = -sin (psia) * cos (nua);

    parrMtrxPer[7] = sin (psia) * sin (nua) * cos (gamma_a)
            + cos (psia) * sin(gamma_a);

    parrMtrxPer[8] = -sin (psia) * sin (nua) * sin (gamma_a)
            + cos (psia) * cos(gamma_a);
}
//-----------------------------------
// матрица перехода parrMtrxPer из связанной СК ЛА в скоргостную СК ГОСТ 20058-80
// в нормальную систему координат
// Скоростная СК это подвижная СК ось OX которой совпадает с направлением скорости ЛА
// ось OY лежит в плоскости симметрии ЛА и направлена к верхней части ЛА
// Связанная система координат это подвижная СК,
// осями которой являются продольная ось OX, нормальная ось OY
// и поперечная ось ОZ ЛА
// INPUT:
//alf - угол атаки, это угол между продольной осью ЛА и
// проекцией скорости ЛА на плоскость OXY связанной СК
// угол атаки следует считать положительным,
// если проекция скорости ЛА на нормальную ось отрицательна
//bet - угол скольжения
// это угол между направлением скорости ЛА и плоскостью OXY связанной СК
void QCoordSystTrsf::calcMatrxPer_from_SvSK_To_SkSK(const double alf ,
    const double bet ,double* parrMtrxPer)
{
    parrMtrxPer[0] = cos(alf) * cos (bet);
    parrMtrxPer[1] = - sin(alf) * cos (bet);
    parrMtrxPer[2] = sin (bet);

    parrMtrxPer[3] = sin (alf);
    parrMtrxPer[4] = cos(alf);
    parrMtrxPer[5] = 0.;

    parrMtrxPer[6] = -cos(alf) * sin (bet);
    parrMtrxPer[7] = sin(alf) * sin (bet);
    parrMtrxPer[8] = cos (bet);
}

//---------------------------------------
// матрица перехода из связанной  системы координат в нормальную (ГОСТ)
void QCoordSystTrsf::calcMatrxPer_from_SvSK_To_NormSK_(double *arrVeloNSK,const double gamma_a
                ,const double alf,const double bet ,double* parrMtrxPer)
{
// 1 матрица перехода из скоростной системы координат в нормальную
double arrMtrxPer_SkSK_to_NormSK[9] = {0.};
// psia - скоростной угол рыскания
//это угол между осью OXg ноормальной СК и проекцией скорости на горизонтальную
// плоскость OXgZg нормальной СК
//СКоростной угол рыскания следует считать положительным, когда ось OXg
//совмещается  с проекцией скорости на горизонт плоскость OXgZg
// поворотом вокруг оси OYg против часовой стрелке
double psia = -atan2(arrVeloNSK[2], arrVeloNSK[0]);
// nua - скоростной угол тангажа
// угол между скоростной осью и горизоньтальной плоскостью
// OxgZg нормальной СК
double nua = asin(arrVeloNSK[1]/ Norm3(arrVeloNSK));
calcMatrxPer_from_SkSK_To_NormSK( psia ,
nua ,gamma_a  ,arrMtrxPer_SkSK_to_NormSK);
///
//double  parrMtrxPer1[9] = {0.};
//calcMatrxPer_from_SSK_To_DecartPrSK(&(arrVectSostNSK[3]),  parrMtrxPer1);
// 2 матрица перехода из связанной СК в скоростную систему координат
double arrMtrxPer_SvSK_to_SkSK[9] = {0.};
QCoordSystTrsf::calcMatrxPer_from_SvSK_To_SkSK(alf ,
bet ,arrMtrxPer_SvSK_to_SkSK);
///

// 3 матрица перехода из связанной СК в нормальную систему координат

MtrxMultMatrx(arrMtrxPer_SkSK_to_NormSK, 3, 3,  arrMtrxPer_SvSK_to_SkSK,3, parrMtrxPer);
}
//---------------------------------------
// матрица перехода из нормальной  системы координат в связанную (ГОСТ)
void QCoordSystTrsf::calcMatrxPer_from_NormSK_To_SvSK_(double *arrVeloNSK,const double gamma_a
                ,const double alf,const double bet ,double* parrMtrxPer)
{
    double arrMtrxPer_SvSK_to_NormSK[9] = {0.};
    calcMatrxPer_from_SvSK_To_NormSK_(arrVeloNSK, gamma_a
                    , alf, bet , arrMtrxPer_SvSK_to_NormSK);

    MatrTransp( arrMtrxPer_SvSK_to_NormSK, 3, 3, parrMtrxPer);
}

//---------------------------------
QPolygonF QCoordSystTrsf::rotate(const QPolygonF plg0, const double a)
{
    QPolygonF plgRez = plg0;
    double matrP[4] = {0.};
    memset(matrP,0,4*sizeof(double));

    matrP[0] = cos(a);
    matrP[1] = -sin(a);
    matrP[2] = sin(a);
    matrP[3] = cos(a);
    for (int i =0 ; i < plg0.size(); ++i)
    {
        double arr[2] ={ plg0.at(i).x(), plg0.at(i).y()};
        double arr1[2] = {0.};
         MtrxMultMatrx(matrP,2, 2,  arr,1, arr1);
         plgRez.replace(i, QPointF(arr1[0], arr1[1]));
    }
    return plgRez;
}
//-----------------------------------------

// Вычисление производной матрицы перехода из ПСК в КГСК
//
//
void QCoordSystTrsf::calc_dM_PSK_to_KGSK_po_dt (double *arrEilerCntrKP0, double *arr_dEilers_po_dt0
         , double *arr_dM_PSK_to_KGSK_po_dt)
{
    double arrEilerCntrKP[3] = {0.}, arr_dEilers_po_dt[3] = {0.};
    MatrxMultScalar(arrEilerCntrKP0, 3, 1, -1.,arrEilerCntrKP);
    MatrxMultScalar(arr_dEilers_po_dt0, 3, 1, -1.,arr_dEilers_po_dt);
double arr_dM_po_dQ[9] = {0.},arr_dM_po_dPsi[9]= {0.}, arr_dM_po_dTet[9]= {0.};
QCoordSystTrsf::calc_dMLeft_psk_to_kgsk_po_dQ(arrEilerCntrKP,arr_dM_po_dQ);

QCoordSystTrsf::calc_dMLeft_psk_to_kgsk_po_dPsi(arrEilerCntrKP,arr_dM_po_dPsi);

QCoordSystTrsf::calc_dMLeft_psk_to_kgsk_po_dTet(arrEilerCntrKP,arr_dM_po_dTet);

MatrxMultScalar(arr_dM_po_dQ, 1, 9, arr_dEilers_po_dt[2],arr_dM_po_dQ);
MatrxMultScalar(arr_dM_po_dPsi, 1, 9, arr_dEilers_po_dt[0],arr_dM_po_dPsi);
MatrxMultScalar(arr_dM_po_dTet, 1, 9, arr_dEilers_po_dt[1],arr_dM_po_dTet);

double arrT0[9] = {0.},arrT1[9] = {0.},arrT2[3] = {0.};
MtrxSumMatrx(arr_dM_po_dQ, arr_dM_po_dPsi,1, 9, arrT0) ;
MtrxSumMatrx(arrT0, arr_dM_po_dTet,1, 9, arr_dM_PSK_to_KGSK_po_dt) ;
}



