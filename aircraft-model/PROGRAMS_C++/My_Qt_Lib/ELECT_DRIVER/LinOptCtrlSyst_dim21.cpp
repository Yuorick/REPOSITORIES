#include "LinOptCtrlSyst_dim21.h"
#include "MatrixProccess.h"
#include  <string.h>
#include <math.h>
#include <float.h>
#include "Equations.h"
#include "Comp.h"

QLinOptCtrlSyst_dim21::QLinOptCtrlSyst_dim21():QLinOptCtrlSyst()
{
    memset(marr_ekt_b, 0., 4 * sizeof(double));
    memset(marr_t_ekt_b, 0., 4 * sizeof(double));
    memset(marr_b0Inv, 0., 4 * sizeof(double));
    mbInit = false;
    mbMultRoots = false;
    mbZeroRoot =   false;
    mk1 = 0.;
    mk2 = 0.;

}
// конструктор копирования
 QLinOptCtrlSyst_dim21 :: QLinOptCtrlSyst_dim21 (const  QLinOptCtrlSyst_dim21 &R):QLinOptCtrlSyst( R)
 {
     memcpy(marr_ekt_b, R.marr_ekt_b, 4 * sizeof(double));
     memcpy(marr_t_ekt_b,R.marr_t_ekt_b, 4 * sizeof(double));
     memcpy(marr_b0Inv, R.marr_b0Inv, 4 * sizeof(double));
     mbInit = R.mbInit;
     mbMultRoots = R.mbMultRoots;
     mbZeroRoot = R.mbZeroRoot;
     mk1 = R.mk1;
     mk2 = R.mk2;

 }

 // оператор присваивания
  QLinOptCtrlSyst_dim21  &QLinOptCtrlSyst_dim21::operator=( const QLinOptCtrlSyst_dim21  &R)
 {
      if(this == &R)
      {
          return *this;
      }
     QLinOptCtrlSyst:: operator= (R);
     memcpy(marr_ekt_b, R.marr_ekt_b, 4 * sizeof(double));
     memcpy(marr_t_ekt_b,R.marr_t_ekt_b, 4 * sizeof(double));
     memcpy(marr_b0Inv, R.marr_b0Inv, 4 * sizeof(double));
     mbInit = R.mbInit;
     mbMultRoots = R.mbMultRoots;
     mbZeroRoot = R.mbZeroRoot;
     mk1 = R.mk1;
     mk2 = R.mk2;
     return *this ;
 }

// параметрический конструктор
QLinOptCtrlSyst_dim21:: QLinOptCtrlSyst_dim21 (double *pA,
                                   double *pB,
                                   double *pC,
                                   const double T0,
                                   double *arrPhVect0,
                                   const double TCur,
                                   double *arrPhVect):QLinOptCtrlSyst( 2,1,pA, pB, pC,T0,arrPhVect0, TCur,arrPhVect)
{
    memset(marr_ekt_b, 0., 4 * sizeof(double));
    memset(marr_t_ekt_b, 0., 4 * sizeof(double));
    memset(marr_b0Inv, 0., 4 * sizeof(double));

    mbInit =calcMtrx_ForFundMtrxCalculation();
}
//--------------------------------------------------------
//-----------------------------------------
// вычисление фундамернтальной матрицы для лин диф уравнния
// 2-го порядка с постоянными коэффициентами
// INPUT:
//arrA[4] - матрица
//VAlT - время
// OUTPUT:
//arrF[4] - фундам матрица в момент T
//решается характ уравнение, опредляются корни k1, k2
// вычисляется матрица arrF0] = exp(k1*T) , arrF3] = exp(k2*T)
bool QLinOptCtrlSyst_dim21::calcMtrx_ForFundMtrxCalculation()
{

    const double a =1.;
    const double b = -(mpA[0] + mpA[3]);
    const double c = mpA[0] * mpA[3] - mpA[1] * mpA[2];
    TComp x1,x2;
    // возвращает
    // 0 - 2 действительных некраьных корня
    // 1 - 2 действительных кратных корня
    // 2- имеется по крайней мере один нулевой корень a!= 0, c=0
    // 3 - 2 комплексно сопряженных корня
    // 4 -  1 действительный корень (a =0)
    // 5  - несовместность a=b=0, c!=0
    // 6 - тождество )a=b=c=0
    int irez = SolvEq2( a, b, c,x1,x2);
    if ( 3 == irez)
    {
      return false;
    }
    if (2 == irez)
    {
        mbZeroRoot = true;

    }
    else
    {
      mbZeroRoot = false;
    }
    double arrT[4] ={0.};
        if ((fabs(mpA[1]) <0.00000001)&&(fabs(mpA[0] - mpA[3])< 0.0000001)) // кратные
        {
          marr_ekt_b[0] = marr_ekt_b[3] = 1.;
          marr_t_ekt_b[2] = 1.;
          marr_b0Inv[0] = marr_b0Inv[3] = 1.;
          mbMultRoots = true;
          mk1 = mpA[0];
          mk2 = mpA[0];

            return true;
        }

        if (fabs(mpA[1]) <0.00000001) // некратные
        {

            double temp2 = mpA[2]/ (mpA[0] - mpA[3]);
            marr_ekt_b[0] = 1;
            marr_ekt_b[1] = 0.;
            marr_ekt_b[2] = temp2;
            marr_ekt_b[3] = 1.;
            marr_b0Inv[0] = 1.;
            marr_b0Inv[1] = 0.;
            marr_b0Inv[2] = -temp2;
            marr_b0Inv[3] = 1.;
            mbMultRoots = false;
            mk1 = mpA[0];
            mk2 = mpA[3];
            return true;
        }

        double k =0., k1 = 0., k2 = 0.,temp =0., temp1 =0, temp2 =0., alf2 =0., bet2 =0., val_d =0.;
        switch(irez)
        {
        case 1:// кратные
             k = x1.m_Re;

             val_d = mpA[1] * mpA[2] - mpA[0]* mpA[3];
             temp1 = mpA[3]/ mpA[1]  + val_d/ mpA[1]/ k ;
            alf2 = -val_d/mpA[1]/k/k;
             bet2 = mpA[3]/ mpA[1] + val_d/ mpA[1]/ k ;


            marr_ekt_b[0] = 1.;
            marr_ekt_b[2] = temp1;
            marr_ekt_b[3] = alf2 ;
            marr_t_ekt_b[1]  = 1.;
            marr_t_ekt_b[3] = bet2;
            marr_b0Inv[0] =1. ;
            marr_b0Inv[1] = 0.;
            marr_b0Inv[2] = -temp1 /alf2;
            marr_b0Inv[3] = 1./ alf2;
            mbMultRoots = true;
            mk1 = k;
            mk2 = k;
            return true;
        break;

        case 0:  // не кратные
             k1 = x1.m_Re;
             k2 = x2.m_Re;


             val_d = mpA[1] * mpA[2] - mpA[0]* mpA[3];
             temp1 = mpA[3]/ mpA[1]   + val_d/ mpA[1]/ k1 ;
             temp2 = mpA[3]/ mpA[1]    + val_d/ mpA[1]/ k2 ;



            marr_ekt_b[0] = 1.;
            marr_ekt_b[1] = 1.;
            marr_ekt_b[2] = temp1 ;
            marr_ekt_b[3] = temp2 ;
            marr_b0Inv[0] =temp2 / (temp2 - temp1) ;
            marr_b0Inv[1] = -1./ (temp2 - temp1) ;
            marr_b0Inv[2] = -temp1 / (temp2 - temp1) ;
            marr_b0Inv[3] = 1./ (temp2 - temp1) ;
            mbMultRoots = false;
            mk1 = k1;
            mk2 = k2;
            return true;
            break;

        case 2:// один нулевой
            if ( (fabs(mpA[1])> 0.00000001)&&(fabs(mpA[0])< 0.00000001)) //один нулевой
            {   //a00 == 0
                k = x1.m_Re;

                 temp = mpA[1]/mpA[3];

                marr_ekt_b[0] = temp ;
                marr_ekt_b[1] = 1.  ;
                marr_ekt_b[2] = 1;
                marr_ekt_b[3] = 0.;
                marr_b0Inv[0] = 0. ;
                marr_b0Inv[1] = 1. ;
                marr_b0Inv[2] = 1. ;
                marr_b0Inv[3] = -temp ;
                mbMultRoots = false;
                mk1 = k;
                mk2 = 0.;
                return true;
            }
            else
            {
                double k = x1.m_Re;



                double temp0 = mpA[3]/ mpA[1] ;
                double temp2 = -mpA[2]/ mpA[3]  ;



                marr_ekt_b[0] = 1;
                marr_ekt_b[1] = 1.;
                marr_ekt_b[2] = temp0 ;
                marr_ekt_b[3] = temp2 ;
                marr_b0Inv[0] =temp2 / (temp2 - temp0) ;
                marr_b0Inv[1] = -1./ (temp2 - temp0) ;
                marr_b0Inv[2] = -temp0 / (temp2 - temp0) ;
                marr_b0Inv[3] = 1./ (temp2 - temp0) ;
                mbMultRoots = false;
                mk1 = k;
                mk2 = 0.;
                return true;
            }

        default:
        break;
        };

    return true;
}


//----------------------------------------------------
//решение задачи на быстродействие
//INPUT:
//VAlUqu0-
//VAlOmegaStat-
//OUTPUT:
//*pvalt0 - момент рервого переключения
//*pvalt1 - конечвное время
//*isignum0 - знак управления на первом участке
//
/*
bool QLinOptCtrlSyst_dim21::calcSwitchPointTimes_FixedSignum0_Old(double *arrxEnd,const double VAlUqu0
                                      ,double *pvalt0,double *pvalt1, const int isignum0)
{

    double arrt[2] = {0.001, 0.001}, arrFGr[2] = {0.}, arr_dFGr_po_dt[4] = {0.}
          , arr_dFGr_po_dtInv[4] = {0.}, arrtPrev[2] = {0.}, arrDelt[2] = {0.}, arrt1[2] = {0.};
    arrtPrev[0] = arrt[0];
    arrtPrev[1] = arrt[1];
    bool breturn = false;
    const double VAlEps = 0.001;
    double delta = -1., delta1 = -1.;
    double valCoeff = 0.4, arrDelt0[2] ={0.};
    for (int i =0; i < 600; ++i)
    {
      calcFgr_and_dFgr(VAlUqu0, isignum0,arrt,arrxEnd, arrFGr, arr_dFGr_po_dt) ;
      if(! InverseMtrx2(arr_dFGr_po_dt, arr_dFGr_po_dtInv))
      {
          return false;
      }
      MtrxMultMatrx(arr_dFGr_po_dtInv,2, 2, arrFGr,1, arrDelt);
      MatrxMultScalar(arrDelt, 1, 2, valCoeff ,arrDelt0);
      MtrxMinusMatrx(arrt, arrDelt0,1, 2, arrt1) ;

      memcpy(arrt, arrt1, 2 * sizeof(double));
      delta = NormVect2(arrDelt);
      delta1 = NormVect2(arrFGr);
      if((delta< VAlEps) && (delta1 < VAlEps))
      {
          *pvalt0 = arrt[0];
          *pvalt1 = arrt[0] +arrt[1];
          return true;
      }

    }
    return false;

}
*/
//-----------------------------------------------------

//вычисление матрицы перехода (фундаментальной)
//
//
//
void QLinOptCtrlSyst_dim21:: calcFundMtrx(const double VAlT,double *arrFundMtrx)
{
    double mtrxDiag0[4] = {0.},mtrxDiag1[4] = {0.};

    mtrxDiag0[0] = exp(VAlT * mk1);
    if(mbMultRoots)
    {
     mtrxDiag0[3] =  exp(VAlT * mk1);
     mtrxDiag1[0] = mtrxDiag1[3] =  VAlT *exp(VAlT * mk1);
    }
    else
    {
        mtrxDiag0[3] =  exp(VAlT * mk2);
    }

    double arrT0[4] = {0.}, arrT1[4] = {0.},arrT2[4] = {0.};
    MtrxMultMatrx(marr_ekt_b,2, 2, mtrxDiag0,2, arrT0) ;
    MtrxMultMatrx(marr_t_ekt_b,2, 2, mtrxDiag1,2, arrT1) ;
    MtrxSumMatrx(arrT0, arrT1,2, 2, arrT2) ;
    MtrxMultMatrx(arrT2,2, 2, marr_b0Inv,2, arrFundMtrx) ;
}
//---------------------------------------------------------------
void QLinOptCtrlSyst_dim21::calcFgr_and_dFgr(const double VAlUqu0, const int isignum0
                  , double *arrt, double *arrxEnd, double *arrFGr, double *arr_dFGr_po_dt)
{
    double arrvect0[2] = {0.},arrvect1[2] = {0.},arrT0[4] = {0.}
                                                            ,arrT1[4] = {0.},arrvect2[2] = {0.}
         ,arrvect3[2] = {0.},arrvect4[2] = {0.},arrvect5[2] = {0.},arrvect6[2] ={0.}, arrPost0[2] = {0.};
    // 1. F(t0)
    double arrMtrxFt0[4] = {0.};
    calcFundMtrx(arrt[0],arrMtrxFt0);
    ///

    // 5. F(t0)*x0
    MtrxMultMatrx(arrMtrxFt0,2, 2, marrPhVect,1, arrvect0) ;//  arrvect0 = F(t1)* x0

    // 2. Integr(F(tau)) limits 0 t0
    double arrMtrxDefInegral_0_t0[4] = {0.};
    calcDefiniteIntegralFundMtrx(0,arrt[0],arrMtrxDefInegral_0_t0);

    ///
    // 7. BU
    double arrvectBU[2] = {0.};
     MatrxMultScalar(mpB, 1, 2, VAlUqu0 * ((double)isignum0),arrvectBU) ;
    //

     // 5.
     MtrxMultMatrx(arrMtrxDefInegral_0_t0,2, 2, arrvectBU,1, arrvect1) ;//

     //
     MtrxSumMatrx(  arrvect0,  arrvect1,2, 1, arrPost0) ;// F(t0)x0+ Integr



    // 3. F(t1)
    double arrMtrxFt1[4] = {0.};
    calcFundMtrx(arrt[1],arrMtrxFt1);
    ///

    //
     MtrxMultMatrx(arrMtrxFt1,2, 2, arrPost0,1, arrvect3) ;//

     //

    // 4. Integr(F(tau)) limits 0 t1
    double arrMtrxDefInegral_0_t1[4] = {0.};
    calcDefiniteIntegralFundMtrx(0,arrt[1],arrMtrxDefInegral_0_t1);

    ///

    //
    MtrxMultMatrx(arrMtrxDefInegral_0_t1,2, 2, arrvectBU,1, arrvect4) ;//
    ///

    //
    MtrxMultMatrx(arrMtrxFt1,2, 2, arrPost0,1, arrvect5) ;

   MtrxMinusMatrx( arrvect5, arrvect4,2, 1, arrvect6) ;
   MtrxMinusMatrx( arrvect6, arrxEnd,2, 1, arrFGr) ;

  ///-----------------------
  ///
  ///
   double arr_dfGrTransp[4] = {0.};
   double arrvect7[2] = {0.},  arrvect8[2] = {0.},  arrvect9[2] = {0.},  arrvect10[2] = {0.}, arrvect11[2] = {0.};
   MtrxMultMatrx(mpA,2, 2,  arrvect0,1, arrvect7) ;
   MtrxMultMatrx(arrMtrxFt0,2, 2,  arrvectBU,1, arrvect8) ;
    MtrxSumMatrx( arrvect7, arrvect8,2, 1, arrvect10) ;
    MtrxMultMatrx(arrMtrxFt1,2, 2,  arrvect10,1, arr_dfGrTransp) ;

    MtrxMultMatrx(mpA,2, 2,  arrvect5,1, arrvect8) ;
    MtrxMultMatrx(arrMtrxFt1,2, 2,  arrvectBU,1, arrvect9) ;
    MtrxMinusMatrx( arrvect8, arrvect9,2, 1, &(arr_dfGrTransp[2])) ;

    MatrTransp(arr_dfGrTransp, 2, 2, arr_dFGr_po_dt);


}
//------------------------------------------------
bool QLinOptCtrlSyst_dim21::calcSwitchPointTimes(double *arrxEnd,const double VAlUqu0
                                      ,double *pvalt0,double *pvalt1, int *isignum0)
{
    *isignum0 = 1;

    if( calcSwitchPointTimes_FixedSignum0(arrxEnd, VAlUqu0
                                             ,pvalt0,pvalt1, *isignum0))
      {
          return true;
      }
      *isignum0 = -1;
      if( calcSwitchPointTimes_FixedSignum0(arrxEnd, VAlUqu0
                                             ,pvalt0,pvalt1, *isignum0))
      {
          return true;
      }
   return false;
}
//--------------------------------------
// вычисление неопределенного интеграла от матрицы перехода
void QLinOptCtrlSyst_dim21:: calcIntegralFundMtrx(const double VAlT,double *arrIntegralFundMtrx)
{
    double mtrxDiag0[4] = {0.},mtrxDiag1[4] = {0.};
    // корни некраьные, неравные нулю
    if ((!mbMultRoots)&&(!mbZeroRoot))
    {
      mtrxDiag0[0] = 1./mk1  *exp(VAlT * mk1);
      mtrxDiag0[3] = 1./mk2  *exp(VAlT * mk2);
    }
    ///

    //
    // корни некраьные, один равен нулю
    if ((!mbMultRoots)&&(mbZeroRoot))
    {
      mtrxDiag0[0] = 1./mk1  *exp(VAlT * mk1);
      mtrxDiag0[3] = VAlT;
    }
    ///

    // корнги кратные
    if (mbMultRoots)
    {
        mtrxDiag0[0] = 1./mk1  *exp(VAlT * mk1);
        mtrxDiag0[3] = 1./mk1  *exp(VAlT * mk2);
        mtrxDiag1[0] = 1./mk1  *(VAlT - 1./mk1) *exp(VAlT * mk1);
        mtrxDiag1[3] = mtrxDiag1[0];

    }

    double arrT0[4] = {0.}, arrT1[4] = {0.},arrT2[4] = {0.};
    MtrxMultMatrx(marr_ekt_b,2, 2, mtrxDiag0,2, arrT0) ;
    MtrxMultMatrx(marr_t_ekt_b,2, 2, mtrxDiag1,2, arrT1) ;
    MtrxSumMatrx(arrT0, arrT1,2, 2, arrT2) ;
    MtrxMultMatrx(arrT2,2, 2, marr_b0Inv,2, arrIntegralFundMtrx) ;
}
//------------------------------------------------------------

void QLinOptCtrlSyst_dim21:: calcDefiniteIntegralFundMtrx(const double VAla,const double VAlb,double *arrDefIntegralFundMtrx)
{
    double  arrmtrxT0[4] = {0.}, arrmtrxT1[4] = {0.};
    calcIntegralFundMtrx(VAlb,arrmtrxT0);
    calcIntegralFundMtrx(VAla,arrmtrxT1);
    MtrxMinusMatrx(arrmtrxT0, arrmtrxT1, 2, 2, arrDefIntegralFundMtrx);
}

//----------------------------------------------------------------------------
void QLinOptCtrlSyst_dim21:: calcSpecialIntegralFundMtrx(const double VAlx1,const double VAla,const double VAlb,double *arrDefIntegralFundMtrx)
{
    double  arrmtrxT0[4] = {0.}, arrmtrxT1[4] = {0.};
    calcDefiniteIntegralFundMtrx(VAla, VAlx1, arrmtrxT0);
    calcDefiniteIntegralFundMtrx(VAlb, VAlx1,arrmtrxT1);
    MtrxMinusMatrx(arrmtrxT0, arrmtrxT1, 2, 2, arrDefIntegralFundMtrx);
}
//--------------------------------------------------

void QLinOptCtrlSyst_dim21::drag( const double VAlIntegrStep, const double VAlU0
                 ,const double valt0,const double valt1)
{
    const long double VAlT = valt1;
    const long double VAlLongIntegrStep = VAlIntegrStep;

    int iQuantSteps = int(VAlT/VAlLongIntegrStep);
    double arrRightPart[2] = {0.};
    long double arrLongRightPart[2] = {0.}, arr_f [2] = {0.},arrt0[4] = {0.}, arrt1[4] = {0.};
    long double arrLongPhVect[2] = {0.};
    arrLongPhVect[0] = marrPhVect[0];
    arrLongPhVect[1] = marrPhVect[1];
    long double valLongU = -1.;

    for (int i = 0; i < (iQuantSteps -1); ++i )
    {

      if(mTCur < valt0)
      {
       valLongU = VAlU0;
      }
      else
      {
       valLongU = -VAlU0;
      }

      calcRightPart(valLongU,arrLongPhVect,arrLongRightPart);

      //
      MatrxMultScalar(arrLongRightPart, 1, 2, VAlLongIntegrStep,arrt0);

      MtrxSumMatrx(arrLongPhVect, arrt0,1, 2, arrt1) ;

      memcpy(arrLongPhVect, arrt1,2 * sizeof(long double));
      mTCur += VAlIntegrStep;

    }

    marrPhVect[0] = arrLongPhVect[0]  ;
    marrPhVect[1] = arrLongPhVect[1] ;
 }



//----------------------------------------
void QLinOptCtrlSyst_dim21::calcRightPart(const double VAlU, double *arrRightPart)
{
    double arrT0[2] = {0.}, arrT1[2] = {0.}, arrT2[2] = {0.};
    MtrxMultMatrx(mpA,2, 2, marrPhVect,1, arrT0) ;
    MtrxSumMatrx(arrT0, mpC,1, 2, arrT1) ;
    MatrxMultScalar(mpB, 1, 2, VAlU,arrT2) ;
    MtrxSumMatrx(arrT1, arrT2,1, 2, arrRightPart) ;
}
//------------------------------------------------
//----------------------------------------
void QLinOptCtrlSyst_dim21::calcRightPart(const long double VAlU
              , long double *arrPhVect,long  double *arrRightPart)
{
    long double arrT0[2] = {0.}, arrT1[2] = {0.}, arrT2[2] = {0.}, pA[4]= {0.}, pC[2]= {0.},pB[2] = {0.};
    pA[0] = mpA[0];
    pA[1] = mpA[1];
    pA[2] = mpA[2];
    pA[3] = mpA[3];
    pC[0] = mpC[0];
    pC[1] = mpC[1];
    pB[0] = mpB[0];
    pB[1] = mpB[1];

    MtrxMultMatrx(pA,2, 2, arrPhVect,1, arrT0) ;
    MtrxSumMatrx(arrT0, pC,1, 2, arrT1) ;
    MatrxMultScalar(pB, 1, 2, VAlU,arrT2) ;
    MtrxSumMatrx(arrT1, arrT2,1, 2, arrRightPart) ;
}
//------------------------------------------------


//передвигание класса с текущего сремени на  время VAlNextTime
//при помощи матрицы перехода и ее интеграла
//INPUT:
//VAlNextTime - мледующий момоент времени
// VAlU0 - величина постоянного управления
//
void QLinOptCtrlSyst_dim21::goToNextTime(const double VAlNextTime, const double VAlU0)
{
    const double VAlDeltaNextTime0 = VAlNextTime - mTCur;
    double arrvect0[2] = {0.},arrvect1[2] = {0.},arrT0[4] = {0.};
    // 1. F(t0)
    double arrMtrxFt0[4] = {0.};
    calcFundMtrx( VAlDeltaNextTime0,arrMtrxFt0);
    ///

    // 5. F(t0)*x0
    MtrxMultMatrx(arrMtrxFt0,2, 2, marrPhVect,1, arrvect0) ;//  arrvect0 = F(t1)* x0

    // 2. Integr(F(tau)) limits 0 t0
    double arrMtrxDefInegral_0_t0[4] = {0.};
    calcDefiniteIntegralFundMtrx(0, VAlDeltaNextTime0,arrMtrxDefInegral_0_t0);

    ///
    // 7. BU
    double arrvectBU[2] = {0.};
     MatrxMultScalar(mpB, 1, 2, VAlU0 ,arrvectBU) ;
    //

     // 5.
     MtrxMultMatrx(arrMtrxDefInegral_0_t0,2, 2, arrvectBU,1, arrvect1) ;//

     //
     MtrxSumMatrx(  arrvect0,  arrvect1,2, 1, marrPhVect) ;// F(t0)x0+ Integr
     mTCur = VAlNextTime;



}
//--------------------------------------------------------
//передвигание класса с текущего сремени на  время VAlNextTime
//при помощи матрицы перехода и ее интеграла
//INPUT:
//VAlNextTime - мледующий момоент времени
// VAlU0 - величина постоянного управления
//
void QLinOptCtrlSyst_dim21::goToNextTimeEiler(const double VAlNextTime, const double VAlU0
                                              , const double VAlIntegrStep)
{
    const long double VAlDeltaNextTime0 = VAlNextTime - mTCur;
    //const long double VAlT =VAlNextTime;
    const long double VAlLongIntegrStep = VAlIntegrStep;

    int iQuantSteps = int(VAlDeltaNextTime0/VAlLongIntegrStep);
    double arrRightPart[2] = {0.};
    long double arrLongRightPart[2] = {0.}, arr_f [2] = {0.},arrt0[4] = {0.}, arrt1[4] = {0.};
    long double arrLongPhVect[2] = {0.};
    arrLongPhVect[0] = marrPhVect[0];
    arrLongPhVect[1] = marrPhVect[1];
    long double valLongU = VAlU0;

    for (int i = 0; i < (iQuantSteps -1); ++i )
    {

      calcRightPart(valLongU,arrLongPhVect,arrLongRightPart);

      //
      MatrxMultScalar(arrLongRightPart, 1, 2, VAlLongIntegrStep,arrt0);

      MtrxSumMatrx(arrLongPhVect, arrt0,1, 2, arrt1) ;

      memcpy(arrLongPhVect, arrt1,2 * sizeof(long double));
      mTCur += VAlIntegrStep;

    }
    marrPhVect[0] = arrLongPhVect[0]  ;
    marrPhVect[1] = arrLongPhVect[1] ;

}

//------------------------------------------
//решение задачи на быстродействие
//INPUT:
//VAlUqu0-
//VAlOmegaStat-
//OUTPUT:
//*pvalt0 - момент рервого переключения
//*pvalt1 - конечвное время
//*isignum0 - знак управления на первом участке
//
bool QLinOptCtrlSyst_dim21::calcSwitchPointTimes_FixedSignum0(double *arrxEnd,const double VAlUqu0
                                      ,double *pvalt0,double *pvalt1, const int isignum0)
{

    double arrt[2] = {0.0001, 0.001}, arrFGr[2] = {0.}, arr_dFGr_po_dt[4] = {0.}
          , arr_dFGr_po_dtInv[4] = {0.}, arrtPrev[2] = {0.}, arrDelt[2] = {0.}, arrt1[2] = {0.};
    arrtPrev[0] = arrt[0];
    arrtPrev[1] = arrt[1];


    double delta = -1., delta1 = -1.;
    double valCoeff = 0.051, arrDelt0[2] ={0.};
    const double VAlEps = 0.001  * valCoeff;
    for (int i =0; i < 600; ++i)
    {
      calcFgr_and_dFgr_Full(VAlUqu0, isignum0,arrt,arrxEnd, arrFGr, arr_dFGr_po_dt) ;
      if(! InverseMtrx2(arr_dFGr_po_dt, arr_dFGr_po_dtInv))
      {
          return false;
      }
      MtrxMultMatrx(arr_dFGr_po_dtInv,2, 2, arrFGr,1, arrDelt);
      MatrxMultScalar(arrDelt, 1, 2, valCoeff ,arrDelt0);
      MtrxMinusMatrx(arrt, arrDelt0,1, 2, arrt1) ;

      memcpy(arrt, arrt1, 2 * sizeof(double));
      delta = NormVect2(arrDelt);
      delta1 = NormVect2(arrFGr);
      if((delta< VAlEps) && (delta1 < VAlEps))
      {
          *pvalt0 = arrt[0];
          *pvalt1 = arrt[0] +arrt[1];
          return true;
      }

    }
    return false;

}
//-----------------------------------------------------

//---------------------------------------------------------------
void QLinOptCtrlSyst_dim21::calcFgr_and_dFgr_Full(const double VAlUqu0, const int isignum0
                  , double *arrt, double *arrxEnd, double *arrFGr, double *arr_dFGr_po_dt)
{
    double arrvect0[2] = {0.},arrvect1[2] = {0.}
         ,arrvect3[2] = {0.},arrvect4[2] = {0.},arrvect5[2] = {0.}
        ,arrvect6[2] ={0.},arrvect14[2] ={0.},arrvect15[2] ={0.}, arrPost0[2] = {0.};
    // 1. F(t0)
    double arrMtrxFt0[4] = {0.};
    calcFundMtrx(arrt[0],arrMtrxFt0);
    ///

    // 5. F(t0)*x0
    MtrxMultMatrx(arrMtrxFt0,2, 2, marrPhVect,1, arrvect0) ;//  arrvect0 = F(t1)* x0

    // 2. Integr(F(tau)) limits 0 t0
    double arrMtrxDefInegral_0_t0[4] = {0.};
    calcDefiniteIntegralFundMtrx(0,arrt[0],arrMtrxDefInegral_0_t0);

    ///
    // 7. BU
    double arrvectBU[2] = {0.};
     MatrxMultScalar(mpB, 1, 2, VAlUqu0 * ((double)isignum0),arrvectBU) ;
    //

     // 5.
     MtrxMultMatrx(arrMtrxDefInegral_0_t0,2, 2, arrvectBU,1, arrvect1) ;//

     //
     MtrxSumMatrx(  arrvect0,  arrvect1,2, 1, arrPost0) ;// F(t0)x0+ Integr



    // 3. F(t1)
    double arrMtrxFt1[4] = {0.};
    calcFundMtrx(arrt[1],arrMtrxFt1);
    ///

    //
     MtrxMultMatrx(arrMtrxFt1,2, 2, arrPost0,1, arrvect3) ;//

     //

    // 4. Integr(F(tau)) limits 0 t1
    double arrMtrxDefInegral_0_t1[4] = {0.};
    calcDefiniteIntegralFundMtrx(0,arrt[1],arrMtrxDefInegral_0_t1);

    ///

    //
    MtrxMultMatrx(arrMtrxDefInegral_0_t1,2, 2, arrvectBU,1, arrvect4) ;//
    ///

    //
    MtrxMultMatrx(arrMtrxFt1,2, 2, arrPost0,1, arrvect5) ;

   MtrxMinusMatrx( arrvect5, arrvect4,2, 1, arrvect6) ;
   MtrxMinusMatrx( arrvect6, arrxEnd,2, 1, arrvect15) ;

   ///
   // 2. Integr(F(tau)) limits 0 t0
   double arrMtrxDefInegral_0_t0_plus_t1[4] = {0.};
   calcDefiniteIntegralFundMtrx(0,arrt[0] + arrt[1],arrMtrxDefInegral_0_t0_plus_t1);
   MtrxMultMatrx(arrMtrxDefInegral_0_t0_plus_t1,2, 2, mpC,1, arrvect14) ;
   MtrxSumMatrx( arrvect15, arrvect14,2, 1, arrFGr) ;


   ///
  ///-----------------------
  ///
  ///
   double arr_dfGrTransp[4] = {0.};
   double arrvect7[2] = {0.},  arrvect8[2] = {0.},  arrvect9[2] = {0.}
           ,  arrvect10[2] = {0.}, arrvect11[2] = {0.}, arrvect12[2] = {0.};
   MtrxMultMatrx(mpA,2, 2,  arrvect0,1, arrvect7) ;
   MtrxMultMatrx(arrMtrxFt0,2, 2,  arrvectBU,1, arrvect8) ;
    MtrxSumMatrx( arrvect7, arrvect8,2, 1, arrvect10) ;
    MtrxMultMatrx(arrMtrxFt1,2, 2,  arrvect10,1, arr_dfGrTransp) ;

    MtrxMultMatrx(mpA,2, 2,  arrvect5,1, arrvect8) ;
    MtrxMultMatrx(arrMtrxFt1,2, 2,  arrvectBU,1, arrvect9) ;

    MtrxMinusMatrx( arrvect8, arrvect9,2, 1, arrvect12) ;

    // 1. F(t0)
    double arrMtrxF_t0_plus_t1[4] = {0.};
    calcFundMtrx(arrt[0] + arrt[1],arrMtrxF_t0_plus_t1);
    MtrxMultMatrx(arrMtrxF_t0_plus_t1,2, 2,  mpC,1, arrvect11) ;
     MtrxSumMatrx( arrvect12, arrvect11,2, 1,  &(arr_dfGrTransp[2])) ;

    ///

    MatrTransp(arr_dfGrTransp, 2, 2, arr_dFGr_po_dt);


}
//------------------------------------------------
