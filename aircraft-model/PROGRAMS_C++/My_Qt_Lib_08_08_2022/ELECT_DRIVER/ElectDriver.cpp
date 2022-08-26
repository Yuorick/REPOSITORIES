#include <math.h>
#include  <string.h>
#include  <float.h>
#include <wchar.h>
#include <stdio.h>
#include <time.h>
#include "ElectDriver.h"
#include "MatrixProccess.h"
#include "Comp.h"
#include "ParamsID.h"


// коэффициент лобового сопротивления пластины
#define VAL_CX 0.3
// плотность атмосферы
#define RO_ATM 1.22
// к-во пар обмоток
extern const double CONST_ZP;
//extern const double CONST_GARMON_COEFF;

class QRezPointTraj;


//---------------------------------------------------------------------------

   QElectDriver:: QElectDriver()
{
    mElectMotor = QElectMotor();
   // mpLoad = QLoad();
    mpLoad = NULL;
    memset(marrSpreadParams, 0, QPARAMS * sizeof(double));
   //mMomOut = 0.;
    mInvSumJ0 = 0.;
}

//---------------------------------------------------------------------------

// конструктор копирования
     QElectDriver :: QElectDriver (const  QElectDriver &R)
 {
    mElectMotor = R.mElectMotor;
    mpLoad = R.mpLoad;
    memcpy(marrSpreadParams, R.marrSpreadParams, QPARAMS * sizeof(double));
   // mMomOut = R.mMomOut;
    mInvSumJ0 = R.mInvSumJ0;
 }

 // оператор присваивания
  QElectDriver  &QElectDriver::operator=( const QElectDriver  &R)
 {
      mElectMotor = R.mElectMotor;
      mpLoad = R.mpLoad;
      memcpy(marrSpreadParams, R.marrSpreadParams, QPARAMS* sizeof(double));
     // mMomOut = R.mMomOut;
      mInvSumJ0 = R.mInvSumJ0;

     return *this ;
 }


  // парам конструктор
     QElectDriver:: QElectDriver (const QElectMotor ElectMotor
         ,  QLoad*  Load, const double *arrSpreadParams)
 {
    mElectMotor = ElectMotor;
    mpLoad = Load;
    memcpy(marrSpreadParams, arrSpreadParams, 3 * sizeof(double));
   // mMomOut = MomOut;
    mInvSumJ0 = calcInvSumJ0();

 }

     // парам конструктор 2
   // QElectDriver:: QElectDriver (const QElectMotor ElectMotor ,  QBrls*  pBrls)
   // {
   //    mElectMotor = ElectMotor;
    //   mpLoad = pBrls;
    //}

    //
    void QElectDriver::calcMagnRightPart_Without_ResidualMom(double *arrPhVect, double *arrRightPart )
    {
    // const double  JPayLoad = mElectMotor.mJ0 + (*mpLoad).mJPayLoad;

     arrRightPart[0] = calcInductMom(arrPhVect[2])  * mInvSumJ0;
     arrRightPart[1] = - mElectMotor.mResist * mElectMotor.mInvL * arrPhVect[1]
             + mElectMotor.mZp *arrPhVect[0]*arrPhVect[2] ;
    arrRightPart[2] = - mElectMotor.mZp * arrPhVect[0]*arrPhVect[1]
            -  mElectMotor.mResist* mElectMotor.mInvL  *arrPhVect[2]
           -mElectMotor.mZp *mElectMotor.mPsi_f * mElectMotor.mInvL  * arrPhVect[0] ;
    arrRightPart[3] = arrPhVect[0];
    if (QVARS > 4)
    {
        // интеграл положения
     arrRightPart[4] = arrPhVect[3];
    }

    }

    //------------------------------------------------------
    // Вычисление индукционного крутящеог момента
    double  QElectDriver::calcInductMom(const double VAl_Iqu)
    {
      return 3./2. * mElectMotor.mZp* mElectMotor.mPsi_f* VAl_Iqu;
    }



    //-----------------------------------------------------------------

    double QElectDriver::max_(const double a,const double b)
    {
        return (a >b)?a:b;
    }
    double QElectDriver::sign_(const double a)
     {
         return (a >0)?1:0;
     }

    //---------------------------------------------------------------------
    void QElectDriver::calcFx( TEnvironment &Environment, double *arrPhVect, double *arrFx,const double VAlMomOut)
    {

      calcMagnRightPart_Without_ResidualMom(arrPhVect, arrFx );
      double valMom = calcMomResidual_plus_AirMom_plus_MomOut(Environment,arrPhVect[3],arrPhVect[0], VAlMomOut);

      //double valResistAirMom = mpLoad->calcAirResistMom(Environment,arrPhVect[3],arrPhVect[0]);

      //const double  JPayLoad = mElectMotor.mJ0 + mpLoad->mJPayLoad;
     // arrFx [0] = arrFx [0] + valResistAirMom/ JPayLoad ;

     // double valMomResidual = calcMomResidual(arrPhVect[3], arrPhVect[0]);

      //arrFx[0] =  arrFx[0]  + valMomResidual/JPayLoad ;
      arrFx[0] =  arrFx[0]  + valMom * mInvSumJ0 ;

    }

    //----------------------------------------------
    // вычисление остаточного момента с имитацией возмущений в правой части
   double QElectDriver::calcMomResidual(const double VAlTettaStat,const double VAlOmegaStat)
    {
        double valResidualMom = mElectMotor.calcSumMomNoise( VAlTettaStat, VAlOmegaStat) ;
        return valResidualMom;

    }

    //----------------------------------------------------------------------

  double QElectDriver::calcMomResidual_plus_AirMom_plus_MomOut(TEnvironment &Environment
                                 ,const double VAlTetta,const double VAlOmega,const double VAlMomOut)
{

  double valResistAirMom = mpLoad->calcAirResistMom(Environment,VAlTetta,VAlOmega);

  double valMomResidual = calcMomResidual(VAlTetta,VAlOmega);


  return valMomResidual + valResistAirMom + VAlMomOut;
 }

    // вычисление матрицы частных произыодных подсиситемы уравнен6ий
    void QElectDriver::fill_df_po_px_and_mtrxB_(const double *arrPhVect
                     ,double *arr_dF_po_dx,double * arr_dF_po_dW)
    {

        //const double JPayLoad = mElectMotor.mJ0 + mpLoad->mJPayLoad;
        memset(arr_dF_po_dW, 0, 2 * QVARS * sizeof(double));
        calcMtrxB(arr_dF_po_dW);

        memset(arr_dF_po_dx, 0, QVARS * QVARS * sizeof(double));

        TEnvironment Environment(0.,0.,0.);


        arr_dF_po_dx[0] = mpLoad->calc_dAirResistMom_po_dOmega( Environment, arrPhVect[3], arrPhVect[0]) * mInvSumJ0;
        arr_dF_po_dx[1]  = 0.;
        arr_dF_po_dx[2]  = 3./2. * mElectMotor.mZp * mElectMotor.mPsi_f *mInvSumJ0;

        arr_dF_po_dx[QVARS   ]  = mElectMotor.mZp *arrPhVect[2];
        arr_dF_po_dx[QVARS +1]  = - mElectMotor.mResist* mElectMotor.mInvL  ;
        arr_dF_po_dx[QVARS +2]  = mElectMotor.mZp *arrPhVect[0];

        arr_dF_po_dx[2 *QVARS   ]  = -mElectMotor.mZp *(arrPhVect[1] + mElectMotor.mPsi_f * mElectMotor.mInvL  );
        arr_dF_po_dx[2 *QVARS +1]  = -mElectMotor.mZp *arrPhVect[0];
        arr_dF_po_dx[2 *QVARS +2] = - mElectMotor.mResist * mElectMotor.mInvL ;

        arr_dF_po_dx[3 *QVARS] = 1.;

        if(QVARS >4)
        {
         arr_dF_po_dx[4 *QVARS +3] = 1.;
        }


    }
//--------------------------------------------

    // интегрирование сиситемы уравнений на интервале времени [0;VAlTime]
    // с шагом VAlIntegrStep0
    // arrPhVectBegin[QVARS] - начальный фазовый вектор
    // arrPhVectEnd[QVARS] - конечный фазовый вектор
    // arrU[2] - постоянные управления
void QElectDriver::dragPhVect(TEnvironment &Environment, double *arrPhVectBegin,double *arrU
                   ,const double VAlTime,const double VAlIntegrStep0,double *arrPhVectEnd,const double VAlMomOut)
{
    int iQuantSteps = VAlTime /VAlIntegrStep0;
    double arrPhVectCur[QVARS] = {0.};
    memcpy(arrPhVectCur, arrPhVectBegin, QVARS * sizeof(double));
    double valTCur = 0.;
    for (int i =0; i < iQuantSteps;++i)
    {
      if(valTCur > VAlTime)
      {
          break;
      }
      doEilerStep(Environment,arrPhVectCur, arrU,VAlIntegrStep0, VAlMomOut);

      valTCur += VAlIntegrStep0;

    }

    doEilerStep(Environment,arrPhVectCur, arrU,VAlTime - valTCur, VAlMomOut);
    memcpy(arrPhVectEnd, arrPhVectCur, QVARS * sizeof(double));

 }
//-------------------------------------
// шаг метода Эйлера
void QElectDriver::doEilerStep(TEnvironment &Environment,double *arrPhVectCur
 , double *arrU, const double VAlStep,const double VAlMomOut)
{
    double arrFx[QVARS] = {0.}, arrT0[QVARS] = {0.}, arrT1[QVARS] = {0.};
    calcFx(Environment, arrPhVectCur,arrFx,VAlMomOut);

    double arrB[2 * QVARS] = {0.}, arrBU[QVARS] = {0.}, arrRigthPart[QVARS] = {0.};
    calcMtrxB(arrB);
    MtrxMultMatrx(arrB, QVARS, 2,  arrU,1, arrBU);
    MtrxSumMatrx(arrFx, arrBU,1, QVARS, arrRigthPart) ;

    MatrxMultScalar(arrRigthPart, 1, QVARS, VAlStep,arrT0);
    MtrxSumMatrx(arrT0, arrPhVectCur,1, QVARS, arrT1) ;
    memcpy(arrPhVectCur, arrT1, QVARS * sizeof(double));

}

//-----------------------------------------------------
// вычисление матрицы частных производных правой части
// по параметрам LInv,R,Psi_f, J, Cx, Cv
// OUTPUT:
//arr_dF_po_dAlf[ QVARS *QPARAMS]
//
void QElectDriver::calc_dF_po_dAlf(double *arrEstPhVect,double *arr_dF_po_dAlf
            ,const double VAlMomOut)
{
 memset (arr_dF_po_dAlf, 0, QVARS *QPARAMS * sizeof(double));
 //const double  JPayLoad = mElectMotor.mJ0 + (*mpLoad).mJPayLoad;

 TEnvironment Environment(0.,0.,0.);
 double valResistAirMom = mpLoad->calcAirResistMom(Environment,arrEstPhVect[3],arrEstPhVect[0]);
 double valMomResidual = calcMomResidual(arrEstPhVect[3], arrEstPhVect[0]);
 double valInductMom = calcInductMom(arrEstPhVect[2]);
 double valMom = valInductMom + valResistAirMom +valMomResidual + VAlMomOut;

 //
 // проверка
// double arrFx[QVARS] = {0.};
//calcMagnRightPart_Without_ResidualMom(arrEstPhVect, arrFx );
// double valMomMagn = arrFx[0] * JPayLoad ;
// double valMom1 = calcMomResidual_plus_AirMom_plus_MomOut(Environment,arrEstPhVect[3],arrEstPhVect[0]);

// valMom1 =  valMom1 + valMomMagn ;

 //

 arr_dF_po_dAlf[2] =  3./2. * mElectMotor.mZp*  arrEstPhVect[2]  *mInvSumJ0;
 arr_dF_po_dAlf[3] =  valMom;
 arr_dF_po_dAlf[4] = (*mpLoad).calc_dMa_po_Cx(arrEstPhVect[0])*mInvSumJ0;
 arr_dF_po_dAlf[5] = (*mpLoad).calc_dMa_po_Cv(arrEstPhVect[0])* mInvSumJ0;
 arr_dF_po_dAlf[6] = mInvSumJ0;

 arr_dF_po_dAlf[QPARAMS] =  -mElectMotor.mResist * arrEstPhVect[1];
 arr_dF_po_dAlf[QPARAMS +1] = -arrEstPhVect[1]* mElectMotor.mInvL ;

 arr_dF_po_dAlf[2 * QPARAMS] = -(mElectMotor.mResist * arrEstPhVect[2]  + mElectMotor.mZp*  arrEstPhVect[0] *  mElectMotor.mPsi_f) ;
 arr_dF_po_dAlf[2 * QPARAMS + 1] =-arrEstPhVect[2]* mElectMotor.mInvL ;
 arr_dF_po_dAlf[2 * QPARAMS + 2] =-mElectMotor.mZp*  arrEstPhVect[0]* mElectMotor.mInvL ;

}

// вычисление вектора частных производных вектора BU по индуктивности
// arr_dBU_po_dIndL[4] - OUTPUT
void QElectDriver::calc_dBU_po_dIndL(double *arrU,double *arr_dBU_po_dIndL)
{

  memset(arr_dBU_po_dIndL, 0, QVARS  * sizeof(double));
  arr_dBU_po_dIndL[1] = arrU[0] ;
  arr_dBU_po_dIndL[2] = arrU[1] ;
}
//-----------------------------------------
//------------------------------------------------
void QElectDriver::calcMtrxB(double *arrB)
{
    memset(arrB, 0, 2 * QVARS * sizeof(double));
    arrB[2] = arrB[5] =  mElectMotor.mInvL ;
}

//----------------------------------------------
// получение параметрра с номером i
double  QElectDriver::getParam(const int i)
{
    switch(i)
    {
    case 0:
     return    mElectMotor.getLInv();
        break;
    case 1:
      return    mElectMotor.mResist;
        break;
    case 2:
      return    mElectMotor.mPsi_f;
        break;
    case 3:
       return  mInvSumJ0;
        break;
    case 4:
       return  mpLoad->mCx;
        break;
    case 5:
       return  mpLoad->mCv;
        break;
    case 6:
       //return  mMomOut;
        break;
    case 7:
       return  mElectMotor.getMaxHarmonAmp () ;
        break;
    case 8:
       return  mElectMotor.getHarmonPh () ;
        break;
    default: break;
    };
}
//----------------------------------------------
// установка параметрра с номером i
void  QElectDriver::setParam(const int i, const double VAlParam)
{
    switch(i)
    {
    case 0:
        mElectMotor.setLInv( VAlParam);
        break;
    case 1:
        mElectMotor.mResist = VAlParam;
        break;
    case 2:
        mElectMotor.mPsi_f = VAlParam;
        break;
    case 3:
       mInvSumJ0 = VAlParam ;
       mpLoad->mJPayLoad = 1./mInvSumJ0 - mElectMotor.getJ0 ();

        break;
    case 4:
       mpLoad->mCx = VAlParam;
        break;
    case 5:
        mpLoad->mCv = VAlParam;
        break;
    case 6:
       //mMomOut = VAlParam;
        break;
    case 7:
       mElectMotor.setMaxAmpGrmn ( VAlParam);

        break;
    case 8:
       mElectMotor.setHarmonPh( VAlParam);
        break;

    default: break;
    };
}

//-----------------------------------------------------------
void  QElectDriver::calcDynamicFgr_and_arrGrFgr(const double VAlIntegrStep0,QRezPointTraj *parrRezPointTraj
              ,const int LEnarrRezPointTraj,double  *pvallFgr, double *arrGrFgr,const double VAlMomOut)
{
    TEnvironment Environment(0.,0.,0.);
    double valTcur = -VAlIntegrStep0;
    double arrTruePhVect[QVARS] = {0.};

    memcpy(arrTruePhVect, parrRezPointTraj[0].marrX, QVARS * sizeof(double));

    double *arrMtrxY = new double[QVARS * (2 + QPARAMS)];
    memset(arrMtrxY, 0., QVARS * (2 + QPARAMS) * sizeof(double));
    double arr_dF_po_dx[QVARS * QVARS] = {0.},arr_dF_po_dW [2 * QVARS] = {0.}
     , arr_dBU_po_dIndL[QVARS] = {0.}, arr_dF_po_dAlf[QVARS *QPARAMS] = {0.}
     , arr_dF_po_dAlf_Full[QVARS *( 2 + QPARAMS)] = {0.};

    double arrT0[QVARS *( 2 + QPARAMS)] = {0.},arrT1[QVARS *( 2 + QPARAMS)] = {0.};

   // const double  JPayLoad = mElectMotor.mJ0 + (*mpLoad).mJPayLoad;

    memset(arrGrFgr, 0,( 2 + QPARAMS)* sizeof(double));

    *pvallFgr = 0.;


     for (int i = 0; i < (LEnarrRezPointTraj- 1) ; ++i )
     {
         // 1. вычисление матрицы частных производных  df_po_dx
         fill_df_po_px_and_mtrxB_(arrTruePhVect,arr_dF_po_dx,arr_dF_po_dW);
            // надо добавить слагаемые, вызванные гармоничевским моментом
         arr_dF_po_dx[0] += mElectMotor.calc_dHarmMom_po_dOm (arrTruePhVect[3]) * mInvSumJ0;
         arr_dF_po_dx[3] += mElectMotor.calc_dHarmMom_po_dTetta (arrTruePhVect[3], arrTruePhVect[0]) * mInvSumJ0;
         ///

         // 2. вычисление матрицы частных производных правой части по векторру параметров
         calc_dF_po_dAlf(arrTruePhVect,arr_dF_po_dAlf,0.);
         calc_dBU_po_dIndL(parrRezPointTraj[i].marrU,arr_dBU_po_dIndL);
         for (int j = 0; j < QVARS; ++j)
         {
           arr_dF_po_dAlf[j * QPARAMS]+= arr_dBU_po_dIndL[j];
         }
         memset(arr_dF_po_dAlf_Full, 0, QVARS *( 2 + QPARAMS) * sizeof(double));
         for (int j = 0; j < QVARS; ++j)
         {
             memcpy(&(arr_dF_po_dAlf_Full[ (2 + QPARAMS) * j]),&(arr_dF_po_dAlf[ QPARAMS * j]), QPARAMS * sizeof(double) );
         }
            // надо добавить слагаемые, вызванные гармоничевским моментом

        // arr_dF_po_dAlf_Full[3] += -mElectMotor.calc_HarmMom (arrTruePhVect[3], arrTruePhVect[0])/ JPayLoad/ JPayLoad;
         arr_dF_po_dAlf_Full[7] += mElectMotor.calc_dHarmMom_po_dAmp (arrTruePhVect[3], arrTruePhVect[0]) * mInvSumJ0;
         arr_dF_po_dAlf_Full[8] += mElectMotor.calc_dHarmMom_po_dPh  (arrTruePhVect[3], arrTruePhVect[0]) *mInvSumJ0;
       ///


        // 3. шаг интегрирования фазового вектора
         valTcur = valTcur + VAlIntegrStep0;
         double arrPhVectEnd[QVARS] = {0.};
         dragPhVect(Environment,arrTruePhVect, parrRezPointTraj[i].marrU
                    , VAlIntegrStep0, VAlIntegrStep0, arrPhVectEnd,VAlMomOut);
         memcpy(arrTruePhVect,arrPhVectEnd, QVARS * sizeof(double));
         ///

         // 4. шаг интегрирования матрицы разброса параметров
         MtrxMultMatrx( arr_dF_po_dx,QVARS  ,QVARS , arrMtrxY,(2 + QPARAMS), arrT0) ;
         MtrxSumMatrx(arrT0, arr_dF_po_dAlf_Full,QVARS, 2 + QPARAMS, arrT1) ;
         MatrxMultScalar(arrT1,QVARS, 2 + QPARAMS, VAlIntegrStep0, arrT0) ;
         MtrxSumMatrx(arrT0, arrMtrxY,QVARS, 2 + QPARAMS, arrT1) ;
         memcpy(arrMtrxY, arrT1, sizeof(double) * QVARS *( 2 + QPARAMS));
         ///


         // 5. вычисление ZCur = YTr * delX
         double arrDelXCur[QVARS] = {0.}, arrZCur[2 + QPARAMS] = {0.};
         MtrxMinusMatrx(arrTruePhVect, parrRezPointTraj[i+1].marrX,1, QVARS, arrDelXCur);
         MtrxTranspMultMatrx(arrMtrxY,QVARS , ( 2 + QPARAMS), arrDelXCur,1, arrZCur) ;
         *pvallFgr += NormSquareVect(arrDelXCur, QVARS);
         ///


         // 6. сложение
         MtrxSumMatrx(arrGrFgr, arrZCur,1, 2 + QPARAMS, arrT1) ;
         memcpy(arrGrFgr, arrT1, sizeof(double) * ( 2 + QPARAMS));


            }

      MatrxMultScalar(arrGrFgr ,1,  QPARAMS, 2., arrGrFgr) ;

}

//-----------------------------------------------------------
void  QElectDriver::calcDynamicFgr(const double VAlIntegrStep0,QRezPointTraj *parrRezPointTraj
              ,const int LEnarrRezPointTraj,double  *pvallFgr,const double VAlMomOut)
{
    TEnvironment Environment(0.,0.,0.);

    double arrTruePhVect[QVARS] = {0.};

    memcpy(arrTruePhVect, parrRezPointTraj[0].marrX, QVARS * sizeof(double));

    *pvallFgr = 0.;


     for (int i = 0; i < (LEnarrRezPointTraj- 1) ; ++i )
     {
         double arrPhVectEnd[QVARS] = {0.};
         dragPhVect(Environment,arrTruePhVect, parrRezPointTraj[i].marrU
                    , VAlIntegrStep0, VAlIntegrStep0, arrPhVectEnd, VAlMomOut);
         memcpy(arrTruePhVect,arrPhVectEnd, QVARS * sizeof(double));
         ///

         double arrDelXCur[QVARS] = {0.};
         MtrxMinusMatrx(arrTruePhVect, parrRezPointTraj[i+1].marrX,1, QVARS, arrDelXCur);

         *pvallFgr += NormSquareVect(arrDelXCur, QVARS);
      }
}
//-----------------------------------------------------------
void  QElectDriver::calcDynamicFgr_and_arrGrFgr_temp(const double VAlIntegrStep0,QRezPointTraj *parrRezPointTraj
              ,const int LEnarrRezPointTraj,double  *pvallFgr, double *arrGrFgr,const double VAlMomOut)
{
    TEnvironment Environment(0.,0.,0.);
    double valTcur = -VAlIntegrStep0;
    double arrTruePhVect[QVARS] = {0.};

    memcpy(arrTruePhVect, parrRezPointTraj[0].marrX, QVARS * sizeof(double));

    double *arrMtrxY = new double[QVARS * (2 + QPARAMS)];
    memset(arrMtrxY, 0., QVARS * (2 + QPARAMS) * sizeof(double));
    double arr_dF_po_dx[QVARS * QVARS] = {0.},arr_dF_po_dW [2 * QVARS] = {0.}
     , arr_dBU_po_dIndL[QVARS] = {0.}, arr_dF_po_dAlf[QVARS *QPARAMS] = {0.}
     , arr_dF_po_dAlf_Full[QVARS *( 2 + QPARAMS)] = {0.};

    double arrT0[QVARS *( 2 + QPARAMS)] = {0.},arrT1[QVARS *( 2 + QPARAMS)] = {0.};

   // const double  JPayLoad = mElectMotor.mJ0 + (*mpLoad).mJPayLoad;

    memset(arrGrFgr, 0,( 2 + QPARAMS)* sizeof(double));

    *pvallFgr = 0.;

    QElectDriver Driver_temp = *this;

    double del = 0.001;
    QLoad load = *mpLoad;
    Driver_temp.mpLoad = &load;
    double temp = Driver_temp.getParam(3);
    temp  += del;
    Driver_temp.setParam(3,temp);
    double arrTruePhVect_temp[QVARS] = {0.};
     for (int i = 0; i < (LEnarrRezPointTraj- 1) ; ++i )
     {
         // 1. вычисление матрицы частных производных  df_po_dx
         fill_df_po_px_and_mtrxB_(arrTruePhVect,arr_dF_po_dx,arr_dF_po_dW);
            // надо добавить слагаемые, вызванные гармоничевским моментом
         arr_dF_po_dx[0] += mElectMotor.calc_dHarmMom_po_dOm (arrTruePhVect[3]) * mInvSumJ0;
         arr_dF_po_dx[3] += mElectMotor.calc_dHarmMom_po_dTetta (arrTruePhVect[3], arrTruePhVect[0]) * mInvSumJ0;
         ///

         // 2. вычисление матрицы частных производных правой части по векторру параметров
         calc_dF_po_dAlf(arrTruePhVect,arr_dF_po_dAlf,0.);
         calc_dBU_po_dIndL(parrRezPointTraj[i].marrU,arr_dBU_po_dIndL);
         for (int j = 0; j < QVARS; ++j)
         {
           arr_dF_po_dAlf[j * QPARAMS]+= arr_dBU_po_dIndL[j];
         }
         memset(arr_dF_po_dAlf_Full, 0, QVARS *( 2 + QPARAMS) * sizeof(double));
         for (int j = 0; j < QVARS; ++j)
         {
             memcpy(&(arr_dF_po_dAlf_Full[ (2 + QPARAMS) * j]),&(arr_dF_po_dAlf[ QPARAMS * j]), QPARAMS * sizeof(double) );
         }
            // надо добавить слагаемые, вызванные гармоничевским моментом

        // arr_dF_po_dAlf_Full[3] += -mElectMotor.calc_HarmMom (arrTruePhVect[3], arrTruePhVect[0])/ JPayLoad/ JPayLoad;
         arr_dF_po_dAlf_Full[7] += mElectMotor.calc_dHarmMom_po_dAmp (arrTruePhVect[3], arrTruePhVect[0])* mInvSumJ0;
         arr_dF_po_dAlf_Full[8] += mElectMotor.calc_dHarmMom_po_dPh  (arrTruePhVect[3], arrTruePhVect[0])* mInvSumJ0;
       ///


        // 3. шаг интегрирования фазового вектора
         valTcur = valTcur + VAlIntegrStep0;
         double arrPhVectEnd[QVARS] = {0.};
         dragPhVect(Environment,arrTruePhVect, parrRezPointTraj[i].marrU
   , VAlIntegrStep0, VAlIntegrStep0, arrPhVectEnd, VAlMomOut);
         memcpy(arrTruePhVect,arrPhVectEnd, QVARS * sizeof(double));
         ///

         // 4. шаг интегрирования матрицы разброса параметров
         MtrxMultMatrx( arr_dF_po_dx,QVARS  ,QVARS , arrMtrxY,(2 + QPARAMS), arrT0) ;
         MtrxSumMatrx(arrT0, arr_dF_po_dAlf_Full,QVARS, 2 + QPARAMS, arrT1) ;
         MatrxMultScalar(arrT1,QVARS, 2 + QPARAMS, VAlIntegrStep0, arrT0) ;
         MtrxSumMatrx(arrT0, arrMtrxY,QVARS, 2 + QPARAMS, arrT1) ;
         memcpy(arrMtrxY, arrT1, sizeof(double) * QVARS *( 2 + QPARAMS));
         ///


         // 5. вычисление ZCur = YTr * delX
         double arrDelXCur[QVARS] = {0.}, arrZCur[2 + QPARAMS] = {0.};
         MtrxMinusMatrx(arrTruePhVect, parrRezPointTraj[i+1].marrX,1, QVARS, arrDelXCur);
         MtrxTranspMultMatrx(arrMtrxY,QVARS , ( 2 + QPARAMS), arrDelXCur,1, arrZCur) ;
         *pvallFgr += NormSquareVect(arrDelXCur, QVARS);
         ///


         // 6. сложение
         MtrxSumMatrx(arrGrFgr, arrZCur,1, 2 + QPARAMS, arrT1) ;
         memcpy(arrGrFgr, arrT1, sizeof(double) * ( 2 + QPARAMS));

         // интегрирование фаового вектора temp
         double arrPhVectEnd_temp[QVARS] = {0.};
         Driver_temp.dragPhVect(Environment,arrTruePhVect_temp, parrRezPointTraj[i].marrU
                                , VAlIntegrStep0, VAlIntegrStep0, arrPhVectEnd_temp, VAlMomOut);
         memcpy(arrTruePhVect_temp,arrPhVectEnd_temp, QVARS * sizeof(double));

         double arrDiff[ QVARS] = {0.};
         MtrxMinusMatrx(arrTruePhVect_temp, arrTruePhVect,QVARS , 1,arrDiff);
         MatrxMultScalar(arrDiff, QVARS , 1, 1./ del,arrDiff);

         double arrDeriv[QVARS ] = {0.};
         int j0 = 3;
         for(int j =0; j < QVARS; ++j)
         {
           arrDeriv[j ] = arrMtrxY[ j * ( 2 + QPARAMS) + j0];

         }

         if( (i%100)== 0)
         {
             int ttyu = 0;
         }
    }
      MatrxMultScalar(arrGrFgr ,1,  QPARAMS, 2., arrGrFgr) ;

}

//-----------------------------------------------------
double   QElectDriver::calcInvSumJ0()
{

  return 1./ (mElectMotor.getJ0 () + mpLoad->mJPayLoad);
}
//-------------------------------------------
double   QElectDriver::getInvSumJ0()
{
    return mInvSumJ0;
}
//-------------------------------------------------------------
double   QElectDriver::getJPayLoad()
{
    return mpLoad->mJPayLoad;
}
//-----------------------------------------------------
//,const double VAlMomOut)

double   QElectDriver::calcDynamicFgr_lambda(const double VAlIntegrStep,QRezPointTraj *parrRezPointTraj
                     ,const int LEnarrRezPointTraj, const double VAllamb, const int *IArrTargNums
                     , const int LEnTargNums, const double *ARrGrFgrCutting,const double VAlMomOut)
{

    QLoad LoadCur = *(mpLoad);
    QElectMotor MotorCur = mElectMotor;
    QElectDriver DriverCur (MotorCur,  &LoadCur, marrSpreadParams);
    for(int k =0; k < LEnTargNums; ++k)
    {
        double valParam = getParam(IArrTargNums[k]) - VAllamb * ARrGrFgrCutting[k];
       DriverCur.setParam(IArrTargNums[k], valParam);
    }

    double vallFgr = -1;
    DriverCur.calcDynamicFgr( VAlIntegrStep,parrRezPointTraj
            ,LEnarrRezPointTraj, &vallFgr, VAlMomOut );
    return vallFgr;

}

//-------------------------------------------------------

void QElectDriver::calcStationarySolution(const double *arrObjective
     , double*arrStatinaryPhVect,double *valUd, double*valUq, double *pvalDelFi,const double VAlMomOut)
{

    TEnvironment Environment0;
   // QSegment SegmentCur;
    //double arrGSKOrtZero[3] = {0.};
    // arrGSKOrtZero[0] = 1.;
    // TPlane PlaneCur;

    memset(arrStatinaryPhVect, 0, QVARS * sizeof(double));
    arrStatinaryPhVect [0] = arrObjective[1];
    arrStatinaryPhVect [1] = 0.;// Id


   // double temp = mpLoad->calcAirResistMom(Environment0, 0., arrStatinaryPhVect [0]);
    // Iq:
    double valMom = calcMomResidual_plus_AirMom_plus_MomOut(Environment0,arrObjective[0]
            , arrObjective[1], VAlMomOut);
   // double valMomResidual =calcMomResidual(arrObjective[0], arrObjective[1]);
    //double valCurrentIq = arrStatinaryPhVect [2] = -2./ 3. *
      // ( mpLoad->calcAirResistMom(Environment0, 0., arrStatinaryPhVect [0])
         //   + valMomResidual + mMomOut)
          //   / mElectMotor.mZp/mElectMotor.mPsi_f;
    double valCurrentIq = arrStatinaryPhVect [2] = -2./ 3. * valMom/ mElectMotor.mZp/mElectMotor.mPsi_f;
   arrStatinaryPhVect [3] = arrObjective[0];
   *valUd = - arrObjective[1] * valCurrentIq *mElectMotor.mZp / mElectMotor.mInvL;
   *valUq = mElectMotor.mZp * mElectMotor.mPsi_f * arrObjective[1]
           + mElectMotor.mResist * valCurrentIq;
   *pvalDelFi = atan2(-(*valUd),(*valUq));
}

//--------------------------------------------
int QElectDriver::createInputDataReport(wchar_t*FileName, const bool bHeader)

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
fprintf(fw," Электропривод \n");
fprintf(fw,"  обратная величина сум. момента инерции(кг*м*м) =  %6.4f\n",mInvSumJ0);
fclose(fw);
mElectMotor.createInputDataReport(FileName, false);
mpLoad->createInputDataReport(FileName, false);
createInputDataReport_Ctrl(FileName, false);
}
//-------------------------------------------------
int QElectDriver::createInputDataReport_Ctrl(wchar_t*FileName, const bool bHeader)
{

}

