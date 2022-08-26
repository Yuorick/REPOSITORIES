#include "Orientation.h"
#include "CalcCorMatrx.h"
#include "MatrixProccess.h"
#include <math.h>
#include <string.h>
#include <CavityAntenna.h>



extern const double arrBaseFreqData_ShenYan[] = {
 // f         U           V
    360,	-14.4675,	27.925
    ,390	,-12.9511	,18.2519
    ,420	,-11.4752	,9.5762
    ,450	,-10.1067	,1.5303
    ,480	,-8.7084	,-5.7743
    ,495	,-8.0777	,-9.4278
    ,510	,-7.4554	,-13.0074
    ,526	,-6.7941	,-16.7857
};

extern const double arrWorkFreq_ShenYan[]= {
    525.968
    ,519.872
    ,513.200
    ,507.056
    ,500.960
    ,494.768
    ,488.480
    ,482.288
    ,476.096
    ,470.384
    ,464.672
    ,458.864
    ,452.960
    ,446.960
    ,440.960
    ,435.488
    ,429.968
    ,424.976
    ,419.984
    ,414.656
    ,409.184
    ,404.096
    ,398.960
    ,393.968
    ,389.456
    ,384.176
    ,380.000
    ,375.392
    ,370.784
    ,366.368
    ,362.000

};


extern const int QuantBaseFreq_ShenYan = 8;

extern const int QuantWorkFreq_ShenYan = 31;

QOrientation::QOrientation()
{

}
//---------------------------------------
/*
// вычисление единичного вектора направления в АПСК
// INPUT:
//VAlV, VAlU - угловые координаты
//OUTPUT:
//arrUnit[3] - единичный вектор
void QOrientation::calcUnitVector(const double VAlV, const double VAlU, double *arrUnit)
{
    arrUnit[2] = sin(VAlU) ;
    arrUnit[0] = sin(VAlV) ;
    arrUnit[1] = 1. - arrUnit[0] *arrUnit[0] -arrUnit[2] *arrUnit[2] ;
    if (arrUnit[1] < 0.)
    {
      arrUnit[1] = 0.;
      return;
    }
    arrUnit[1] = sqrt(arrUnit[1]);

}*/

  //---------------------------------------
// вычисление единичного вектора направления в АПСК
// INPUT:
//VAlV, VAlU - угловые координаты
//OUTPUT:
//arrUnit[3] - единичный вектор
void QOrientation::calcUnitVector(const double VAlV, const double VAlU, double *arrUnit)
{
    arrUnit[0] = sin(VAlU) ;
    arrUnit[2] = sin(VAlV) ;
    arrUnit[1] = 1. - arrUnit[0] *arrUnit[0] -arrUnit[2] *arrUnit[2] ;
    if (arrUnit[1] < 0.)
    {
      arrUnit[1] = 0.;
      return;
    }
    arrUnit[1] = sqrt(arrUnit[1]);

}
//-----------------------------------------------------


//вычисление углов поворота антенны, устанавливающих единичный вектор arr_e0 по оси OY
// а единичный вектор arr_e1 в плоскость OYZ
// INPUT:
//arr_e0[3]- первый ед вектор
//arr_e1[1] - второй един вектор
//OUTPUT:
//arrMu[3] - вектор углов поворота
//arrMu[0] - поворот на курсовой угол относительно оси OZ по ч стр
// arrMu[1] - поворот на  угол места относительно оси OX против ч стр
// arrMu[2] - поворот на  угол места относительно оси OY по ч стр
void QOrientation::calcRotationAngs( double *arr_e0, double *arr_e1
                   , const int IPredeterminationType, const double* ARrPredeterminatedAngs
                   , double *arrMu, double *matrRez)
{
  double valR = 0., valBet  = 0., valEps  = 0.;
  recalcCoord_INTO_Spherical(arr_e0, valR, valBet, valEps);
  arrMu[0] = valBet;
  arrMu[1] = -valEps; // поворот относительно оси OX
  //по часовой стрелке, а здесь положительное
  //направление исчисляется против часовой стрелки.
  if (IPredeterminationType == 1)
  {
    arrMu[1] = ARrPredeterminatedAngs[0];
  }

  double matrP3Bet[9]={0},matrP1Eps[9]={0}, matrPereh[9] = {0.};
  calcMatrP3(-arrMu[0], matrP3Bet) ;
  calcMatrP1(-arrMu[1], matrP1Eps) ;
  MtrxMultMatrx(matrP1Eps,3, 3, matrP3Bet ,3, matrPereh);

  // для проверки
  double arrt5[3] = {0.};
  MtrxMultMatrx(matrPereh,3, 3, arr_e0 ,1, arrt5);
  ///

  double arr_e11[3] = {0.};
  MtrxMultMatrx(matrPereh,3, 3, arr_e1 ,1, arr_e11);
  double arr_e01[3] = {0.};
  MtrxMultMatrx(matrPereh,3, 3, arr_e0 ,1, arr_e01);

 // double valAlf = atan2(-(arr_e11[0]-arr_e01[0]), (arr_e11[2] -arr_e01[2]));
  double x11= arr_e01[0];
  double x12= arr_e01[1];
  double x13= arr_e01[2];

  double x21= arr_e11[0];
  double x22= arr_e11[1];
  double x23= arr_e11[2];
  //double valAlf = -atan2(x21*x12 - x11 * x22, x23 *x12 - x13 * x22);
  double valAlf = -atan((x21*x12 - x11 * x22)/( x23 *x12 - x13 * x22));
  arrMu[2] = valAlf;

  if (IPredeterminationType == 2)
  {
    arrMu[2] = ARrPredeterminatedAngs[1];
  }
  double matrP2Alf[9] = {0.};

  calcMatrP2_AlogWatchArrow(-arrMu[2], matrP2Alf);

  MtrxMultMatrx(matrP2Alf,3, 3, matrPereh ,3, matrRez);

  // проверк
 // double arrt0[3] = {0.}, arrt1[3] = {0.};
  // MtrxMultMatrx(matrRez,3, 3, arr_e0 ,1, arrt0);
  // MtrxMultMatrx(matrRez,3, 3, arr_e1 ,1, arrt1);


  arrMu[0] = -arrMu[0] ;
  arrMu[1] = arrMu[1] ;
  arrMu[2] = -arrMu[2] ;

  //recalcCoord_INTO_Spherical(arrt1, valR, valBet, valEps);
}
//----------------------------------------------------

void QOrientation::processingExperimentalData( QCavityAntenna cavityAntenna0, double *arrBaseData
           ,const int quantBaseFreq,double *arrWorkLamb, const int quantWorkFreq
       , const int IPredeterminationType, const double* ARrPredeterminatedAngs
      ,double *arrDataOut, double *arrOutMu0, QCavityAntenna *pCavityAntennaEst, double *arrNeviaz)
{   
    double val_aLin = 0., val_dLin = 0.;
    arrNeviaz[0] = estimate_aLin_and_dLin(arrBaseData, quantBaseFreq, &val_aLin, &val_dLin);

    double val_aDM= 0. ,val_dDM= 0., val_LDM= 0. ;
    int Ik = cavityAntenna0.calcDM_k(arrBaseData[0]);
    double arrDM_params[3] = {0.};
    cavityAntenna0.calc_DM_params_array(arrBaseData[0],arrDM_params);

    //arrNeviaz[1] = estimate_aDM_and_dDM_and_LDM(arrBaseData, quantBaseFreq,Ik, &val_aDM, &val_dDM,&val_LDM, arrDM_params);
    arrNeviaz[1] = estimate_aDM_and_dDM_and_LDM_( cavityAntenna0
             ,arrBaseData, quantBaseFreq
             , &val_aDM, &val_dDM, &val_LDM);

    *pCavityAntennaEst = QCavityAntenna(val_aLin, val_dLin,val_aDM,val_LDM,val_dDM);

    //-------------

    for (int i = 0; i < quantWorkFreq; ++i)
    {
        double *p = &arrDataOut[7 * i];
        double lamCur = arrWorkLamb[ i ];
         const double VAlV = pCavityAntennaEst->calcV(lamCur) ;
        // double vCur = fncV_params(lamCur, arrDM_params);
         const double VAlU = pCavityAntennaEst->calcU(lamCur) ;
         p[0] = lamCur;
         p[1] = VAlU;
         p[2] = VAlV;

    }
    int iNum0 = 0;
    const double VAlV0 = arrDataOut[iNum0 * 7 + 2] ;
    const double VAlU0 = arrDataOut[iNum0 * 7 + 1];
    double arr_e0[3] = {0.};
    QOrientation::calcUnitVector( VAlV0,  VAlU0, arr_e0);

    int iNum1 = quantWorkFreq- 1;
    const double VAlV1 = arrDataOut[iNum1 * 7 + 2] ;
    const double VAlU1 = arrDataOut[iNum1 * 7 + 1];
    double arr_e1[3] = {1.};
    QOrientation::calcUnitVector( VAlV1,  VAlU1, arr_e1);


    double  matrRez[9] = {0.};
    QOrientation::calcRotationAngs( arr_e0, arr_e1
                       ,  IPredeterminationType,  ARrPredeterminatedAngs
                       ,arrOutMu0, matrRez);
   // QOrientation::calcEilersAngs( arr_e0, arr_e1, arrOutMu0, matrRez);
    double arrOutMu1[3]= {0.}, matrRez1[9] = {0.};

        const double valQ   = arrOutMu0[0] ;
        const double valPsi = -arrOutMu0[1] ;
        const double valTet = -arrOutMu0[2] ;

        double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

        calcMatrP3(valQ, matrP3Q) ;
        calcMatrP1(valPsi, matrP1Psi) ;
        calcMatrP2(valTet, matrP2Tet) ;
        MtrxMultMatrx( matrP2Tet,3, 3, matrP1Psi,3, arr0);
        MtrxMultMatrx( arr0,3, 3, matrP3Q,3, matrRez1);






    //
    if (IPredeterminationType == 1)
    {
      arrOutMu0[1] =  ARrPredeterminatedAngs[0];
    }


    double arrt0[3] = {0.}, arrt1[3] = {0.};
     MtrxMultMatrx(matrRez,3, 3, arr_e0 ,1, arrt0);
     MtrxMultMatrx(matrRez,3, 3, arr_e1 ,1, arrt1);

    for (int i = 0; i < quantWorkFreq; ++i)
    {
        double *p = &arrDataOut[7 * i];


        double lamCur = arrWorkLamb[ i ];
         const double VAlV = pCavityAntennaEst->calcV(lamCur) ;
         const double VAlU = pCavityAntennaEst->calcU(lamCur) ;
         p[0] = lamCur;
         p[1] = VAlU;
         p[2] = VAlV;
         // единичный     вектор
        double arr_e[3] = {0.};
        QOrientation::calcUnitVector( VAlV,  VAlU, arr_e);
        double valR = 0., valBet  = 0., valEps  = 0.;
        recalcCoord_INTO_Spherical(arr_e, valR, valBet, valEps);
        p[3] = valEps;
        p[4] = valBet;

        // поврот единичного вектора
        double arr_e2[3] = {0.};
        MtrxMultMatrx(matrRez,3,3, arr_e,1, arr_e2)  ;
        recalcCoord_INTO_Spherical(arr_e2, valR, valBet, valEps);
        p[5] = valEps;
        p[6] = valBet;

    }





}
//-------------------------------------
double QOrientation::calcLamb_(const double f)
{
    return VELO_C/f;
}
//---------------------------------
// оценка параметра алин
double QOrientation::estimate_aLin_and_dLin(const double *ArrBaseData,const int quantBaseFreq
                      , double *pval_aLin, double *pval_dLin)
{
    // 1. вычисление начального приближения
    double arrMu[2] = {0.};// arrMu[0] - это aлин, arrMu[1] - это dлин
    calc_a_and_d_Lin_Prev(ArrBaseData[0], ArrBaseData[3 * (quantBaseFreq -1)]
            ,ArrBaseData[1], ArrBaseData[3 * (quantBaseFreq -1) + 1],arrMu );

    // 2.

    double valFGr_u0 = calcFGr_u(ArrBaseData, quantBaseFreq,arrMu);
    double valFGr_u =valFGr_u0;

    double eps = 1.E-8;
    double arrF[2] = {0.}, arrJacF[4] = {0.};
    for (int i = 0; i < 100; ++i)
    {
        valFGr_u = calcFGr_u(ArrBaseData, quantBaseFreq, arrMu);

        calcLin_dF_and_JacF(ArrBaseData, quantBaseFreq, arrMu, arrF, arrJacF);

        double arrJacF_Inv[4] = {0.};
        InverseMtrx2(arrJacF, arrJacF_Inv);

        double arrDel[2] = {0.};
        MtrxMultMatrx(arrJacF_Inv,2, 2,  arrF,1, arrDel) ;

        double del = NormVect2(arrDel);

        double coeff = 10.;
        for(int i = 0; i < 2; ++i)
        {
          arrMu[i] -=  arrDel[i] / coeff;
        }

        if (fabs (del) < eps)
        {

           break;
        }
    }


    *pval_aLin = arrMu[0];
    *pval_dLin = arrMu[1];
    return sqrt(valFGr_u / ((double)quantBaseFreq));
}

//------------------------------------------------------
double QOrientation::estimate_aDM_and_dDM_and_LDM(const double* ArrBaseData,const int quantBaseFreq
           ,const int Ik , double*val_aDM,double*val_dDM,double*val_LDM, double *arrDM_params)
{
    // 1. вычисление начального приближения
   // calc_DM_params_init(ArrBaseData[0], ArrBaseData[3 * quantBaseFreq/2], ArrBaseData[3 * (quantBaseFreq -1)]
         //   ,ArrBaseData[2], ArrBaseData[3 * quantBaseFreq/2 + 2], ArrBaseData[3 * (quantBaseFreq -1) + 2], arrDM_params );
    ///
    // 2.

    double valFGr_V0 = calcFGr_V(ArrBaseData, quantBaseFreq,arrDM_params);
    double valFGr_V =valFGr_V0;

    double eps = 1.E-6;
    double arrF[3] = {0.}, arrJacF[9] = {0.};
    double del = -1.;
    int i = 0;
    for (i = 0; i < 1000; ++i)
    {
        valFGr_V = calcFGr_V(ArrBaseData, quantBaseFreq, arrDM_params);

        calcDM_dF_and_JacF(ArrBaseData, quantBaseFreq, arrDM_params, arrF, arrJacF);

        double arrJacF_Inv[9] = {0.};
        InverseMtrx3(arrJacF, arrJacF_Inv);

        double arrDel[3] = {0.};
        MtrxMultMatrx(arrJacF_Inv,3, 3,  arrF,1, arrDel) ;

        del = Norm3(arrDel);

        double coeff = 10.;
        for(int i = 0; i < 3; ++i)
        {
          arrDM_params[i] -=  arrDel[i] / coeff;
        }

        if (fabs (del) < eps)
        {

           break;
        }
    }


   *val_aDM = 1. / 2. /  arrDM_params[1];
    *val_dDM = (double( Ik ))/ arrDM_params[2];
    *val_LDM = (*val_aDM) *  arrDM_params[0];
   return sqrt(valFGr_V / ((double)quantBaseFreq));

}


//------------------------------------------------------
double QOrientation::estimate_aDM_and_dDM_and_LDM_(const  QCavityAntenna cavityAntenna0
         ,const double* ArrBaseData,const int quantBaseFreq
         , double*val_aDM,double*val_dDM,double*val_LDM)
{
    *val_aDM  = cavityAntenna0.mDM_a;
    *val_LDM = cavityAntenna0.mDM_L;
    *val_dDM = cavityAntenna0.mDM_d;
    double valFGr0 = calNeviazV(cavityAntenna0, ArrBaseData, quantBaseFreq);

    QCavityAntenna cavityAntennaCur = cavityAntenna0;

    const double VAlInterva = 0.01;
    const double VAlIntervd = 0.01;
    const double VAlIntervL = 0.2;

    const double VAlStepL = 0.001;
    const double VAlStepa = 0.0001;
    const double VAlStepd = 0.0001;

    const int Na = 2. * VAlInterva/ VAlStepa;
    const int Nd = 2. * VAlIntervd/ VAlStepd;
    const int NL = 2. * VAlIntervL/ VAlStepL;

    const double VAl_a0 = ArrBaseData[0] /2. + 0.0001;
    const double VAl_d0 = cavityAntenna0.mDM_d - VAlIntervd;
    const double VAl_L0 = cavityAntenna0.mDM_L - VAlIntervL;
    for (int i = 0; i < Na; ++i)
    {
        cavityAntennaCur.mDM_a = VAl_a0 + ((double)i) * VAlStepa;
        for (int j =0; j < NL; ++j)
        {
           cavityAntennaCur.mDM_L = VAl_L0 + ((double)j) * VAlStepL;
           for (int k = 0; k < Nd; ++k)
           {
             cavityAntennaCur.mDM_d = VAl_d0 + ((double)k) * VAlStepd;
             double valFGr = calNeviazV(cavityAntennaCur, ArrBaseData, quantBaseFreq);
             if (valFGr < valFGr0)
             {
               valFGr0 = valFGr;
               *val_aDM = cavityAntennaCur.mDM_a;
               *val_LDM = cavityAntennaCur.mDM_L;
               *val_dDM = cavityAntennaCur.mDM_d;
             }

           }
        }
    }


   return sqrt(valFGr0/ ((double)quantBaseFreq));

}
//---------------------------------------------
double QOrientation::calNeviazV(  QCavityAntenna cavityAntenna
                                  ,const double* ArrBaseData,const int quantBaseFreq)
{

    double sum = 0.;
    for (int i = 0; i < quantBaseFreq; ++i)
    {
       double t = cavityAntenna.calcV(ArrBaseData[3 * i]) - ArrBaseData[3 * i +2] ;
       sum += t * t ;
    }
    return sum;



}
//---------------------------------------
void QOrientation::calcDM_dF_and_JacF(const double *ArrBaseData,const int quantBaseFreq
                ,const double * arrDM_params,double * arrF,double * arrJacF)

{
    memset(arrF, 0, 3 * sizeof(double));
    memset(arrJacF, 0, 9 * sizeof(double));
    double arrFi[3] = {0.}, arrJacFi[9] = {0.};
    for (int i = 0; i < quantBaseFreq; ++i)
    {
     calcDM_dFi_and_JacFi(ArrBaseData[3 * i], ArrBaseData[ 3 * i +2], arrDM_params,arrFi, arrJacFi );
     MtrxSumMatrx(arrF, arrFi,1, 3, arrF) ;
     MtrxSumMatrx(arrJacF, arrJacFi,1, 9, arrJacF) ;
    }
    return ;
}


//---------------------------------------
void QOrientation::calcDM_dFi_and_JacFi(const double lamb, const double V
                , const double* arrDM_params, double*arrdFi, double*arrJacFi )
{
    const double x0 = arrDM_params[0];
    const double x1 = arrDM_params[1];
    const double x2 = arrDM_params[2];

     double y = x0 * sqrt(1. - lamb* lamb * x1* x1) - x2 * lamb;
     if (y > 0.9999999)
     {
         y = 0.9999999;
     }
     if (y < -0.9999999)
     {
         y = -0.9999999;
     }

    const double z = asin(y) +V;

    const double a = sqrt(1. - lamb* lamb * x1* x1);

    const double b = sqrt(1. - y*y);

    double arr_dy_po_dx[3] = {0.};
    arr_dy_po_dx[0] = a;
    arr_dy_po_dx[1] = -x0 * x1 *lamb* lamb/a;
    arr_dy_po_dx[2]= -lamb;

    double arr_dz_po_dx[3] = {0.};
    MatrxMultScalar(arr_dy_po_dx, 1, 3, 1./ b, arr_dz_po_dx);

    MatrxMultScalar(arr_dy_po_dx, 1, 3, 2.* z/ b, arrdFi);

    double arrT0[3] = {0.};//d po dx0 {z/a}
    arrT0[0] = arr_dz_po_dx[0]/a + z * y/a/a/a *arr_dy_po_dx[0];
    arrT0[1] = arr_dz_po_dx[1]/a + z * y /a/a/a *arr_dy_po_dx[1];
    arrT0[2] = arr_dz_po_dx[2]/a + z * y /a/a/a *arr_dy_po_dx[2];

    arrJacFi[0] = 2. * b *arrT0[0];
    arrJacFi[1] = arrJacFi[3] = -2. * x1 * lamb* lamb/b
            *(arr_dy_po_dx[0] *x0 + z /a);

    double temp = 1./b/b/b;
    arrJacFi[4] = -x0 *lamb* lamb
         *(x1/b *arrT0[1] + z /a * temp);

    arrJacFi[2] = arrJacFi[6] = -2. *lamb *arrT0[0];
    arrJacFi[5] = arrJacFi[7] = -2. *lamb *arrT0[1];
    arrJacFi[8] = -2. *lamb *(arr_dz_po_dx[2]/a + z * y/a/a/a *arr_dy_po_dx[2]);

}

//---------------------------------------
double QOrientation::calcFGr_V(const double *ArrBaseFreqData,const int quantBaseFreq, const double *arrDM_params)
{

    double sum = 0.;
    for (int i = 0; i < quantBaseFreq; ++i)
    {

       sum += fncFi_DM(ArrBaseFreqData[i * 3], ArrBaseFreqData[i * 3 + 2], arrDM_params) ;
    }
    return sum;


}
//--------------------------------------
double QOrientation::fncFi_DM(const double lamb, const double VZv, const double *arrDM_params)
{
    double t1 = arrDM_params[0] *sqrt(1. - lamb * lamb *arrDM_params[1]*arrDM_params[1] ) - lamb* arrDM_params[2];
    if (t1 >0.9999999)
    {
        t1 = 0.9999999;
    }
    if (t1 < -0.9999999)
    {
        t1 = -0.9999999;
    }
    double t = asin(t1 )  + VZv;
    return t * t;
}

//--------------------------------------
double QOrientation::fncV_params(const double lamb, const double *arrDM_params)
{
    double t1 = arrDM_params[0] *sqrt(1. - lamb * lamb *arrDM_params[1]*arrDM_params[1] ) - lamb* arrDM_params[2];
    if (t1 >0.9999999)
    {
        t1 = 0.9999999;
    }
    if (t1 < -0.9999999)
    {
        t1 = -0.9999999;
    }
    double t = -asin(t1 ) ;
    return t ;
}
//--------------------------------------
//-------------------------------------------------------
bool QOrientation::calc_DM_params_init(const double lamb1, const double lamb2, const double lamb3
                                   ,const double v1, const double v2, const double v3,double*  arrx )
{
    bool breturn = false;
    // решение нелинейного уравнения относительно x[1]
    double x = 7.;
    double y1 = -sin(v1) *lamb2 + sin(v2) *lamb1;
    double y2 = -sin(v1) *lamb3 + sin(v3) *lamb1;
    double valF0 = calc_f(y1,y2, lamb1,  lamb2,  lamb3, x);
    for (int i =0; i < 2000;++i)
    {
        double valdF = calc_df( y1, y2,  lamb1
                               ,  lamb2,  lamb3, x);
        double valF = calc_f( y1, y2,  lamb1
                               ,  lamb2,  lamb3, x);
        double del = valF/valdF;
        x -= del/10.;

        if (fabs(del)< 1.E-8)
        {
            breturn = true;
            break;
        }
    }
    if(!breturn)
    {
        return false;
    }

    arrx [1] = x;
    arrx [0] = y1 / (lamb2 *calc_g(lamb1, x) - lamb1 *calc_g(lamb2, x) );
    arrx [2] = (arrx [0] *calc_g(lamb1, x) + v1)/ lamb1;
     double tem = y2 / (lamb3 *calc_g(lamb1, x) - lamb1 *calc_g(lamb3, x) );

  return true;
}
//-------------------------------------
double QOrientation::calc_f(const double y1,const double y2,const double  lamb1
           , const double lamb2, const double lamb3,const double x)
{
    return y1 * calcG(lamb1, lamb3, x) - y2 * calcG(lamb1, lamb2, x);

}
//-------------------------------------
double QOrientation::calc_df(const double y1,const double y2,const double  lamb1
           , const double lamb2, const double lamb3,const double x)
{
    return y1 * calcdG(lamb1, lamb3, x) - y2 * calcdG(lamb1, lamb2, x);

}
//-------------------------------------
double QOrientation::calc_d2f(const double y1,const double y2,const double  lamb1
           , const double lamb2, const double lamb3,const double x)
{
    return y1 * calcd2G(lamb1, lamb3, x) - y2 * calcd2G(lamb1, lamb2, x);

}
//----------------------------------------------
double QOrientation::calcG(const double  lam, const double mu, const double x)
{
  return mu * calc_g(lam, x) - lam * calc_g(mu, x);
}
//----------------------------------------------
double QOrientation::calcdG(const double  lam1, const double lam2, const double x)
{
  return lam2 * calc_dg(lam1, x) - lam1 * calc_dg(lam2, x);
}

//----------------------------------------------
double QOrientation::calcd2G(const double  lam, const double mu, const double x)
{
  return mu * calc_d2g(lam, x) - lam * calc_d2g(mu, x);
}
//-----------------------------
double QOrientation::calc_g(const double  lam, const double x)
{
    return sqrt(1. - x*x *lam*lam);
}
//-----------------------------
double QOrientation::calc_dg(const double  lam, const double x)
{
    return -x *lam*lam/sqrt(1. - x*x *lam*lam);
}
//-----------------------------
double QOrientation::calc_d2g(const double  lam, const double x)
{
    return -lam*lam/sqrt(1. - x*x *lam*lam)/(1. - x*x *lam*lam);
}
//--------------------------------------------------
void QOrientation::calc_a_and_d_Lin_Prev(const double lamb1, const double lamb2
        ,const double u1, const double u2, double* arrMu )
{
   double y =  lamb2 * sin (u1) - lamb1 * sin (u2);  

   double t0 = lamb1 * lamb2; // l1* l2

    double x = lamb1 * lamb1 - (lamb2*lamb2 - lamb1*lamb1 - y * y)
           * (lamb2*lamb2 - lamb1*lamb1 - y * y)/ 4./ y/y;

   double a = t0/ 2./ sqrt(x);

   double d = fnc_d(lamb1, u1, a);

   arrMu [0] = a;
   arrMu [1] = d;

}
//---------------------------------------
double QOrientation::fnc_d(const double lamb1,const double u1,const double a)
{
    return lamb1/ 2./ (sqrt( 1. - lamb1 * lamb1 /4./ a/a) - sin(u1));
}

//---------------------------------------
double QOrientation::calcFGr_u(const double *ArrBaseFreqData,const int quantBaseFreq, const double *arrMu)
{

    double sum = 0.;
    for (int i = 0; i < quantBaseFreq; ++i)
    {

       sum += fncFiMu(ArrBaseFreqData[i * 3], ArrBaseFreqData[i * 3 + 1], arrMu) ;
    }
    return sum;

    return 0.;
}
//---------------------------------------
void QOrientation::calcLin_dF_and_JacF(const double *ArrBaseData,const int quantBaseFreq
                ,const double * arrMu,double * arrF,double * arrJacF)

{
    memset(arrF, 0, 2 * sizeof(double));
    memset(arrJacF, 0, 4 * sizeof(double));
    double arrFi[2] = {0.}, arrJacFi[4] = {0.};
    for (int i = 0; i < quantBaseFreq; ++i)
    {
     calcLin_dFi_and_JacFi(ArrBaseData[3 * i], ArrBaseData[ 3 * i +1], arrMu,arrFi, arrJacFi );
     MtrxSumMatrx(arrF, arrFi,1, 2, arrF) ;
     MtrxSumMatrx(arrJacF, arrJacFi,1, 4, arrJacF) ;
    }
    return ;
}

//---------------------------------------
void QOrientation::calcLin_dFi_and_JacFi(const double lamb, const double u
                , const double* arrMu, double*arrdFi, double*arrJacFi )
{
    double a  = arrMu[0];
    double d  = arrMu[1];
    double y = 1. - lamb * lamb / a/a/4.;
    double x = sqrt(y) - lamb/ 2./ d;
    double z = asin(x) - u;

    double arr_dx_po_dMu[2] = {0.};
    arr_dx_po_dMu[0] = lamb * lamb / (4. * a * a * a * sqrt(y)); // dx_po_da
    arr_dx_po_dMu[1] =  lamb / (2. * d * d); // dx_po_da

    arrdFi[0] = 2. * z *arr_dx_po_dMu[0];
    arrdFi[1] = 2. * z *arr_dx_po_dMu[1];

    double arr_dz_po_dMu[2] = {0.};
    MatrxMultScalar(arr_dx_po_dMu, 1, 2, 1./ sqrt(1. - x * x),arr_dz_po_dMu);

    arrJacFi[0] = 1./2. * arr_dz_po_dMu[0] * lamb * lamb / ( a * a * a * sqrt(y))
            +1./2. *z * lamb * lamb*(-3./(a * a * a* a * sqrt(y)) + 1./2. * lamb * lamb/ (a *a * a * a* a * sqrt(y)*y) );

    arrJacFi[1] = lamb * arr_dz_po_dMu[0]/ d/ d;
    arrJacFi[2] = arrJacFi[1];

    arrJacFi[3] = lamb * arr_dz_po_dMu[1]/ d/ d - 2. *lamb  * z /(d * d * d);

}


//--------------------------------------
double QOrientation::fncFiMu(const double lamb, const double uZv, const double *arrMu)
{
    double t = asin(sqrt(1. - lamb * lamb /arrMu[0]/arrMu[0]/ 4. ) - lamb/ 2. /arrMu[1] ) - uZv;
    return t * t;
}
//--------------------------------------


