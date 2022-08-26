#include "Ctrl.h"
#include "MatrixProccess.h"
#include <math.h>
#include <string.h>
#include <wchar.h>
#include <stdio.h>
#include <time.h>
#include "Equations.h"
#include "Comp.h"
#include "LinOptCtrlSyst_dim21.h"
#include "ElectDriver.h"
#include "LinDiffEq.h"
#include "CtrlFollow.h"

extern double aRRTruePhVect[];
//extern const bool BEZ_SHUMOV;
extern  bool BEZ_SHUMOV;
//extern const double constDbl_mh ;
extern  double constDbl_mh ;


class QMeasure;
class QDriverMeasure;
QCtrl::QCtrl():QElectDriver()
{    

    // начальное значение фазового вектора
     //memset (marrPhVect0, 0, 4 * sizeof(double));
    //
     mT0 = 0.;
    //
     mTCur = 0;
    //
    memset (mCmpArrLamb, 0, 4 * sizeof(TComp));
     //mStatSolutionParams = QStatSolutionParams();

     //marrObjective[0] = marrObjective[1] = 0.;

     mFiltr = QFiltr();
     //mEnvironment = TEnvironment();

}
//---------------------------------------------------------------------------



// конструктор копирования
 QCtrl :: QCtrl (const  QCtrl &R):QElectDriver( R)
 {
    mT0 = R.mT0;
    mTCur = R.mTCur;
   // memcpy (marrIntegrators, R.marrIntegrators, QVARS * sizeof(double));
    //mStatSolutionParams = R.mStatSolutionParams;
    //memcpy (marrObjective, R.marrObjective,2 * sizeof(double));
    mFiltr = R.mFiltr;
    memcpy (mCmpArrLamb, R.mCmpArrLamb, 4 * sizeof(TComp));

 }
 //-----------------------------------------------------------------------------------------------------

 // оператор присваивания
  QCtrl  &QCtrl::operator=( const QCtrl  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      QElectDriver:: operator= (R);
      mT0 = R.mT0;
      mTCur = R.mTCur;
     // memcpy (marrIntegrators, R.marrIntegrators, QVARS * sizeof(double));
      //mStatSolutionParams = R.mStatSolutionParams;
      //memcpy (marrObjective, R.marrObjective,2 * sizeof(double));
      mFiltr = R.mFiltr;
      memcpy (mCmpArrLamb, R.mCmpArrLamb, 4 * sizeof(TComp));

      return *this ;
 }

//---------------------------------------------
QCtrl:: QCtrl ( const QElectMotor ElectMotor ,  QLoad*  Load
                , const double *arrSpreadParams,const double VAlMomOut, const double VAlTettaBegin
     ,const double VAlOmegaBegin , const double T0, const double TCur, const double valh,TComp *pCmpArrLamb)
    :QElectDriver ( ElectMotor ,  Load, arrSpreadParams)
{

    mT0 = T0;
    mTCur = TCur;
    double  valUd=0.,  valUq=0.,  valDelFi=0., arrObjVect[2 ] ={0.};
    arrObjVect[0] = VAlTettaBegin;
    arrObjVect[1] = VAlOmegaBegin;
    double arrPhVectBegin[QVARS] = {0.};
    calcStationarySolution( arrObjVect, arrPhVectBegin,&valUd, &valUq, &valDelFi,VAlMomOut);


    mFiltr = QFiltr ( QVARS, 3,arrPhVectBegin, TCur, valh  );
    memcpy (mCmpArrLamb, pCmpArrLamb, 4 * sizeof(TComp));
   // MtrxMinusMatrx(arrPhVectBegin, mStatSolutionParams.marrStatPhVect,QVARS , 1, marrIntegrators);

}


//-----------------------------------------------------------------------------------------------------------
void QCtrl:: CalcEigenValues( const double VAlTettaStat, const double VAlOmegaStat
                              ,const double *arrC,const double VAlMomOut, TComp* cmpArrEigen)
{

}

//-------------------------------------------------------
/*
void QCtrl::calcStationarySolution(const double *arrObjective, double*arrStatinaryPhVect,double *valUd, double*valUq, double *pvalDelFi)
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
    double valMom = calcMomResidual_plus_AirMom_plus_MomOut(Environment0,arrObjective[0], arrObjective[1]);
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
}*/



//------------------------------------------------------

//матрица частных производных правой части по фазовым переменным
// INPUT:
//arrObjective[2] - угловое положение и скоррость цели

void QCtrl::calc_df_po_dx_plus_BC(const double *arrObjective, const double *arrC, double *arrFGr,const double VAlMomOut)
{
double arrStatPhVect[QVARS] = {0.}, arrUstat[2] = {0.};

double valIq, valUd=0.,  valUq=0.,  valDelFi=0.;

calcStationarySolution( arrObjective, arrStatPhVect,&valUd, &valUq, &valDelFi,VAlMomOut);

arrUstat[0] = valUd;
arrUstat[1] = valUq;
double arr_dF_po_dx[QVARS *QVARS] = {0.},arr_dF_po_dW[2 * QVARS] = {0.}
         , arrT[QVARS * QVARS] = {0.}, arrCC[2 *QVARS] ={0.};

fill_df_po_px_and_mtrxB_(arrStatPhVect
                       ,arr_dF_po_dx,arr_dF_po_dW);
memcpy(arrCC, arrC, 2 * QVARS * sizeof(double));
MtrxMultMatrx(arr_dF_po_dW,4, 2, arrCC, QVARS, arrT) ;
MtrxSumMatrx(arrT, arr_dF_po_dx, QVARS, QVARS, arrFGr) ;

}

//-------------------------------------
void QCtrl::calcGearMtrx(const double *arrObjective, TComp* CmpArrLamb,const double VAlMomOut, double* arrC)
{

}

void QCtrl::calcCurU(const double *arrObjective,const double VAlTObjective,const double VAlMomOut, double *arrU)
{

}


//------------------------------------------
void QCtrl::calcStationaryParams(const double *arrObjective,  TComp* CmpArrLamb
      , QStatSolutionParams* pQStatSolutionParams,const double VAlMomOut)
{
    calcGearMtrx(arrObjective, CmpArrLamb, VAlMomOut, (*pQStatSolutionParams).marrGears);
    double valDelFi = 0.;
    calcStationarySolution(arrObjective, (*pQStatSolutionParams).marrStatPhVect
    ,&((*pQStatSolutionParams).marrStatU[0]),&((*pQStatSolutionParams).marrStatU[1]), &valDelFi, VAlMomOut);
}

//----------------------------------------------
// обюработка измерения и вычисление управлений
// INPUT:
// MEasureCur, VAlTcur - измерение и время привязки измерения
//OUTPUT:
//arrU[2] - управления
void QCtrl::processMeasure_and_calcCurU(const double *arrObjective,const double VAlTObjective
       , QDriverMeasure &MEasureCur, double &VAlTcur, double *arrPrevU,const double VAlMomOut, double *arrU)
{
  double arrKww[QVARS * QVARS] = {0.};
  calc_Kww(VAlTcur - mFiltr.mTCur/*, arrPrevU*/,arrKww);
 // MatrxMultScalar(arrKww, QVARS, QVARS, 100.,arrKww);
  double arr_dF_po_dx[QVARS * QVARS] = {0.},arr_dF_po_dW[2 * QVARS] = {0.},arr_BU[QVARS] = {0.};
  fill_df_po_px_and_mtrxB_(mFiltr.marrCurEst,arr_dF_po_dx,arr_dF_po_dW);
  MtrxMultMatrx(arr_dF_po_dW,QVARS, 2,  arrPrevU,1, arr_BU);

  double arrRightPart[4] = {0.};
  calcMagnRightPart_Without_ResidualMom(mFiltr.marrCurEst, arrRightPart );

  ///
  mFiltr.processMeasure( MEasureCur, VAlTcur, arr_dF_po_dx,arr_BU,arrKww);
  correctMomOut(MEasureCur.marrYzv[2], MEasureCur.marrKYzv[8],arrRightPart[0], VAlTcur);
/*
  if(mFiltr.marrCurEst[3] >= 2. * M_PI)
  {
      mFiltr.marrCurEst[3] -= 2. * M_PI;
  }

  if(mFiltr.marrCurEst[3] < 0.)
  {
      mFiltr.marrCurEst[3] += 2. * M_PI;
  }
*/
  //if( BEZ_SHUMOV)
 // {
   // memcpy(mFiltr.marrCurEst,aRRTruePhVect, QVARS * sizeof(double)); // отладка)
 // }


  mTCur = VAlTcur;
  calcCurU(arrObjective,VAlTObjective,VAlMomOut,arrU);

}
//--------------------------------------
/*
// обюработка измерения и вычисление управлений
// INPUT:
// MEasureCur, VAlTcur - измерение и время привязки измерения
//OUTPUT:
//arrU[2] - управления
void QCtrl::processingDriverMeasure(QMeasure &MEasureCur, double &VAlTcur, double *arrPrevU)
{
  double arrKww[QVARS * QVARS] = {0.};
  calc_Kww(VAlTcur - mFiltr.mTCur,arrKww);
 // MatrxMultScalar(arrKww, QVARS, QVARS, 100.,arrKww);
  double arr_dF_po_dx[QVARS * QVARS] = {0.},arr_dF_po_dW[2 * QVARS] = {0.},arr_BU[QVARS] = {0.};
  fill_df_po_px_and_mtrxB_(mFiltr.marrCurEst,arr_dF_po_dx,arr_dF_po_dW);
  MtrxMultMatrx(arr_dF_po_dW,QVARS, 2,  arrPrevU,1, arr_BU);

  double arrRightPart[4] = {0.};
  calcMagnRightPart_Without_ResidualMom(mFiltr.marrCurEst, arrRightPart );

  ///
  mFiltr.processMeasure( MEasureCur, VAlTcur, arr_dF_po_dx,arr_BU,arrKww);

  double valEstMomOut = estimateMomOut(MEasureCur.marrYzv[2], MEasureCur.marrKYzv[8]
          ,arrRightPart[0], VAlTcur);

 // if( BEZ_SHUMOV)
  //{
  //  memcpy(mFiltr.marrCurEst,aRRTruePhVect, QVARS * sizeof(double)); // отладка)
  //}


  mTCur = VAlTcur;

}
*/
//--------------------------------------

// обюработка измерения и вычисление управлений
// INPUT:
// MEasureCur, VAlTcur - измерение и время привязки измерения
//OUTPUT:
//arrU[2] - управления
void QCtrl::processingDriverMeasure(QDriverMeasure &MEasureCur, double &VAlTcur, double *arrPrevU)
{
  double arrKww[QVARS * QVARS] = {0.};
  calc_Kww(VAlTcur - mFiltr.mTCur/*, arrPrevU*/,arrKww);
 // MatrxMultScalar(arrKww, QVARS, QVARS, 100.,arrKww);
  double arr_dF_po_dx[QVARS * QVARS] = {0.},arr_dF_po_dW[2 * QVARS] = {0.},arr_BU[QVARS] = {0.};
  fill_df_po_px_and_mtrxB_(mFiltr.marrCurEst,arr_dF_po_dx,arr_dF_po_dW);
  MtrxMultMatrx(arr_dF_po_dW,QVARS, 2,  arrPrevU,1, arr_BU);

  double arrRightPart[4] = {0.};
  calcMagnRightPart_Without_ResidualMom(mFiltr.marrCurEst, arrRightPart );

  ///
  mFiltr.processMeasure( MEasureCur, VAlTcur, arr_dF_po_dx,arr_BU,arrKww);
  correctMomOut(MEasureCur.marrYzv[2], MEasureCur.marrKYzv[8],arrRightPart[0], VAlTcur);
/*
  if(mFiltr.marrCurEst[3] >= 2. * M_PI)
  {
      mFiltr.marrCurEst[3] -= 2. * M_PI;
  }

  if(mFiltr.marrCurEst[3] < 0.)
  {
      mFiltr.marrCurEst[3] += 2. * M_PI;
  }
*/
 // if( BEZ_SHUMOV)
 // {
  //  memcpy(mFiltr.marrCurEst,aRRTruePhVect, QVARS * sizeof(double)); // отладка)
  //}
  mTCur = VAlTcur;

}
//--------------------------------------
// вычисление коррел матрицы ошибок правой части
void QCtrl::calc_Kww(const double VAlt/*, double *arrU*/, double *arrKww)
{

    memset(arrKww, 0, QVARS * QVARS * sizeof(double)) ;
    // вычисление лщрреляц матрицы случ процесса в правой части.
    // этот процесс постоянен на интервале [0; VAlt]
   // разброс по L
    double arrp[ QVARS * QVARS * QPARAMS] = {0.};

    double *parr_d2F_po_dx_po_dL =  arrp;
    double *parr_d2F_po_dx_po_dR = &arrp[QVARS * QVARS];
    double *parr_d2F_po_dx_po_dPsif = &arrp[ 2 *QVARS * QVARS];
    double *parr_d2F_po_dx_po_dJ = &arrp[ 3 *QVARS * QVARS];
    double *parr_d2F_po_dx_po_dCx = &arrp[ 4 *QVARS * QVARS];
    double *parr_d2F_po_dx_po_dCv = &arrp[ 5 *QVARS * QVARS];
        calc_d2F_po_dx_po_dL(parr_d2F_po_dx_po_dL);

        calc_d2F_po_dx_po_dR( parr_d2F_po_dx_po_dR  );

        calc_d2F_po_dx_po_dPsif(parr_d2F_po_dx_po_dPsif  );

        calc_d2F_po_dx_po_dJ( parr_d2F_po_dx_po_dJ  );

        calc_d2F_po_dx_po_dCx( parr_d2F_po_dx_po_dCx  );

        calc_d2F_po_dx_po_dCv( parr_d2F_po_dx_po_dCv  );





        double arrKSum[QVARS * QVARS] = {0.}, arrDisp[QPARAMS] = {0.};
        for (int i = 0; i < QPARAMS; ++i)
        {
          arrDisp[i] = marrSpreadParams[i] * marrSpreadParams[i] ;
        }
        for (int i =0; i < QPARAMS; ++i)
        {
            double parrRez[QVARS] = {0.}, parrRez0[QVARS *QVARS] = {0.}
           , parrRez1[QVARS *QVARS] = {0.}, parrRez2[QVARS *QVARS] = {0.};
          MtrxMultMatrx(&(arrp[ QVARS *QVARS * i]), QVARS, QVARS,  mFiltr.marrCurEst,1, parrRez) ;
          MtrxMultMatrxTransp(parrRez, QVARS, 1, parrRez, QVARS, parrRez0) ;
          MatrxMultScalar(parrRez0, QVARS, QVARS, arrDisp[i],parrRez1);
          MtrxSumMatrx(arrKSum, parrRez1, QVARS, QVARS, parrRez2) ;
          memcpy(arrKSum, parrRez2, QVARS * QVARS * sizeof(double));

        }
        const double  JInvPayLoad = getInvSumJ0();
        //arrKSum[0] += mElectMotor.getMomResidual()/sqrt(2.)/JPayLoad * mElectMotor.getMomResidual()/sqrt(2.)/JPayLoad;
        arrKSum[0] += mElectMotor.calcDisp_MomNoise (mFiltr.marrCurEst[0])* JInvPayLoad * JInvPayLoad;


       //

        arrKww[0] = arrKSum[0] * VAlt ;
        arrKww[3] = arrKww[QVARS * 3 ] = arrKSum[0] * VAlt * VAlt /2.;
        arrKww[QVARS * 3 + 3] = arrKSum[0] * VAlt * VAlt* VAlt /3.;
        if (QVARS > 4)
        {
         arrKww[4] = arrKww[QVARS * 4] =  arrKSum[0] * VAlt * VAlt* VAlt /6.;
         arrKww[QVARS * 3 + 4] = arrKww[QVARS * 4 + 3] =  arrKSum[0] * (VAlt * VAlt/2. + VAlt * VAlt *VAlt * VAlt/ 8.);
         arrKww[QVARS * 4 + 4] = arrKSum[0] * (VAlt +VAlt * VAlt* VAlt /3.+ VAlt * VAlt* VAlt* VAlt* VAlt /20.);
        }

        // выделение квадратной подматрицы, сщщтветствующей Id, Iqu
        double arrKdqu[4] = {0.};
        arrKdqu[0]=arrKSum[QVARS + 1] + (0.000001)* (0.000001);
        arrKdqu[1] = arrKdqu[2] = arrKSum[QVARS + 2];
        arrKdqu[3] = arrKSum[2 *QVARS + 2 ]+ (0.000001)* (0.000001);

        // интегрование L*arrKdqu*LT на интервале [0;VAlt]
        // L - фундаментальная матрица чичтемы второгот порядка
        // соответствующей Id и Iqu
        double arr_dF_po_dx[QVARS * QVARS] = {0.}, arr_dF_po_dW [2 * QVARS] = {0.};
        fill_df_po_px_and_mtrxB_(mFiltr.marrCurEst,arr_dF_po_dx,arr_dF_po_dW);
        double arrK_Integr_dqu[QVARS] = {0.}, arra[QVARS] = {0.};
        arra[0] = arr_dF_po_dx[QVARS +1];
        arra[1] = arr_dF_po_dx[QVARS +2];
        arra[2] = arr_dF_po_dx[2 *QVARS +1];
        arra[3] = arr_dF_po_dx[2 *QVARS +2];

        const int QUanSteps = 5;
        const double VAlStep = VAlt/ (double(QUanSteps));
        double arrFund[QVARS] = {0.}, parrOut[QVARS] = {0.}
             , arrFundCopy[QVARS] = {0.},parrOut0[QVARS] = {0.};//,arrFund1[4] = {0.};
        for (int i = 0; i < QUanSteps; ++i)
        {
         double valtcur = (double(i))*VAlStep;
         QLinDiffEq::calcFundamentalMtrx_dim2(arra,  valtcur, arrFund) ;
         memcpy(arrFundCopy, arrFund, 4 * sizeof(double));
        // QLinDiffEq::calcFundamentalMtrx(arra, 2, valtcur, arrFund1) ;
         MtrxMultMatrx_MultMatrxTransp(arrFund,arrKdqu,arrFundCopy,2,parrOut);
         MtrxSumMatrx(arrK_Integr_dqu, parrOut,2, 2, parrOut0) ;
         memcpy(arrK_Integr_dqu, parrOut0, 4 * sizeof(double));
         int iii0= 0;
        }
        MatrxMultScalar(arrK_Integr_dqu, 4, 4, VAlStep,arrK_Integr_dqu);

        arrKww[QVARS +1]  = arrK_Integr_dqu[0] ;
        arrKww[QVARS +2]  = arrKww[2 *QVARS +1] = arrK_Integr_dqu[1];
        arrKww[2 *QVARS +2] = arrK_Integr_dqu[3];


}
//---------------------------------------
void QCtrl::calc_d2F_po_dx_po_dL(double *arr_d2F_po_dx_po_dL)
{

memset(arr_d2F_po_dx_po_dL, 0, QVARS * QVARS * sizeof(double));
arr_d2F_po_dx_po_dL[QVARS + 1]  = -mElectMotor.mResist ;
arr_d2F_po_dx_po_dL[2 * QVARS]  = -mElectMotor.mPsi_f * mElectMotor.mZp;
arr_d2F_po_dx_po_dL[2 * QVARS + 2] = arr_d2F_po_dx_po_dL[QVARS + 1] ;
}
//---------------------------------------
void QCtrl::calc_d2F_po_dx_po_dR(double *arr_d2F_po_dx_po_dR)
{
    memset(arr_d2F_po_dx_po_dR, 0, QVARS * QVARS * sizeof(double));
    arr_d2F_po_dx_po_dR[QVARS + 1]  = arr_d2F_po_dx_po_dR [2 *QVARS + 2] = -mElectMotor.mInvL;
}
//---------------------------------------
void QCtrl::calc_d2F_po_dx_po_dPsif(double *arr_d2F_po_dx_po_dPsif)
{
  memset(arr_d2F_po_dx_po_dPsif, 0, QVARS * QVARS * sizeof(double));
  const double  JInvPayLoad = getInvSumJ0();//   mElectMotor.mJ0 + (*mpLoad).mJPayLoad;
  arr_d2F_po_dx_po_dPsif[2] = 3./ 2. * mElectMotor.mZp * JInvPayLoad;
  arr_d2F_po_dx_po_dPsif[2 * QVARS ] = - mElectMotor.mZp * mElectMotor.mInvL;
}
//---------------------------------------
void QCtrl::calc_d2F_po_dx_po_dJ(double *arr_d2F_po_dx_po_dJ)
{
  memset(arr_d2F_po_dx_po_dJ, 0, QVARS * QVARS * sizeof(double));
  // const double  JPayLoad = mElectMotor.mJ0 + (*mpLoad).mJPayLoad;
  TEnvironment Environment(0.,0.,0.);
  arr_d2F_po_dx_po_dJ[0] = mpLoad->calc_dAirResistMom_po_dOmega( Environment, mFiltr.marrCurEst[3], mFiltr.marrCurEst[0]);

}
//---------------------------------------
void QCtrl::calc_d2F_po_dx_po_dCx(double *arr_d2F_po_dx_po_dCx)
{
  memset(arr_d2F_po_dx_po_dCx, 0, QVARS * QVARS * sizeof(double));
  const double  JInvPayLoad =  getInvSumJ0();//  mElectMotor.mJ0 + (*mpLoad).mJPayLoad;
  arr_d2F_po_dx_po_dCx[0] = mpLoad->calc_d2Ma_po_dOm_po_Cx(mFiltr.marrCurEst[0]) * JInvPayLoad;
}
//---------------------------------------
void QCtrl::calc_d2F_po_dx_po_dCv(double *arr_d2F_po_dx_po_dCv)
{
  memset(arr_d2F_po_dx_po_dCv, 0, QVARS * QVARS * sizeof(double));
  const double  JInvPayLoad =  getInvSumJ0();// mElectMotor.mJ0 + (*mpLoad).mJPayLoad;
  arr_d2F_po_dx_po_dCv[0] = mpLoad->calc_d2Ma_po_dOm_po_Cv(mFiltr.marrCurEst[0])* JInvPayLoad;
}

//---------------------------------------
// вычисление остаточного момента с имитацией возмущений в правой части
double QCtrl::calcMomResidual(const double VAlTettaStat,const double VAlOmegaStat)
{
    return 0.;
}

//----------------------------------------------------------
//-----------------------------------------------------------
void  QCtrl::calcStatFgr_and_arrGrFgr( int quantExprs, QStatSolutionParams *parrStatSolutionParams
                                       , double *arrDynamicDelX,double  *pvallFgr
                                       , double *arrGrFgr,const double VAlMomOut)
{

    *pvallFgr = 0.;
    memset(arrGrFgr, 0, QPARAMS * sizeof(double));

    double *arr_dF_po_dx_plus_BC_delX = new double [quantExprs * QVARS ];
    memset(arr_dF_po_dx_plus_BC_delX, 0, quantExprs * QVARS  * sizeof(double));

    for(int i =0; i < quantExprs; ++i)
    {


        // вектор стационарного решения
        double arrXZv[QVARS] = {0.};
        MtrxSumMatrx(&(arrDynamicDelX[i * QVARS]), parrStatSolutionParams[i].marrStatPhVect,QVARS, 1, arrXZv) ;
        ///
        double arrFi[QVARS] = {0.},mtrx_dFi_po_dAlf[QVARS*QPARAMS] = {0.}, arrt0[QVARS] = {0.}, arrt1[QVARS] = {0.};
        calc_vectFi_and_mtrx_dFi_po_dAlf(parrStatSolutionParams[i], arrXZv,arrFi, mtrx_dFi_po_dAlf, VAlMomOut);
        *pvallFgr += NormSquareVect(arrFi, QVARS);
        //MtrxMultMatrx(arrB, QVARS, 2,  parrStatSolutionParams[i].marrGears, QVARS, arrT0);
        
        MtrxTranspMultMatrx(mtrx_dFi_po_dAlf,QVARS, QPARAMS, arrFi, 1, arrt0) ;
        MtrxSumMatrx(arrGrFgr, arrt0,1, QPARAMS, arrt1) ;
        memcpy(arrGrFgr, arrt1, QPARAMS * sizeof(double));


    }

}

//-----------------------------------------------------------

double   QCtrl::calcStatFgr_lambda(int quantExprs, QStatSolutionParams *parrStatSolutionParams
                                          , double *arrDynamicDelX,const double VAllamb, const int *IArrTargNums
                     , const int LEnTargNums, const double *ARrGrFgrCutting)
{

  /*  QLoad LoadCur = *(mpLoad);
    QElectMotor MotorCur = mElectMotor;
     double arrObjective[2] = {0.};
     QCtrlFollow ModelCtrlFollow (mElectMotor,mpLoad
                     ,marrSpreadParams, mMomOut, 0., 0., 0.,arrObjective,constDbl_mh);


     QCtrl *pCtrl = &ModelCtrlFollow;
    for(int k =0; k < LEnTargNums; ++k)
    {
       double valParam = getParam(IArrTargNums[k]) - VAllamb * ARrGrFgrCutting[k];
       (*pCtrl).setParam(IArrTargNums[k], valParam);
    }

    double vallFgr = -1;
    (*pCtrl).calcStatFgr( quantExprs, parrStatSolutionParams
                          , arrDynamicDelX,  &vallFgr );

    return vallFgr;
*/
    return 0;
}
//-----------------------------------------------------------
void  QCtrl::calcStatFgr( int quantExprs, QStatSolutionParams *parrStatSolutionParams
                          , double *arrDynamicDelX,double  *pvallFgr,const double VAlMomOut)
{

    *pvallFgr = 0.;


    double *arr_dF_po_dx_plus_BC_delX = new double [quantExprs * QVARS ];
    memset(arr_dF_po_dx_plus_BC_delX, 0, quantExprs * QVARS  * sizeof(double));

    for(int i =0; i < quantExprs; ++i)
    {
        // вектор стационарного решения
        double arrXZv[QVARS] = {0.};
        MtrxSumMatrx(&(arrDynamicDelX[i * QVARS]), parrStatSolutionParams[i].marrStatPhVect,QVARS, 1, arrXZv) ;
        ///
        double arrFi[QVARS] = {0.},mtrx_dFi_po_dAlf[QVARS*QPARAMS] = {0.}, arrt0[QVARS] = {0.}, arrt1[QVARS] = {0.};
        calc_vectFi_and_mtrx_dFi_po_dAlf(parrStatSolutionParams[i]
                                         , arrXZv,arrFi, mtrx_dFi_po_dAlf, VAlMomOut);
        *pvallFgr += NormSquareVect(arrFi, QVARS);
    }

}
//-----------------------------------------------------------

void  QCtrl::calc_vectFi_and_mtrx_dFi_po_dAlf( QStatSolutionParams StatSolutionParams
         , double *arrXZv, double *arrFi, double *mtrx_dFi_po_dAlf ,const double VAlMomOut)
{

    TEnvironment Environment(0.,0.,0.);
    memset(mtrx_dFi_po_dAlf, 0, QVARS *QPARAMS * sizeof(double));
    memset(arrFi, 0, QVARS  * sizeof(double));


    double arrB [2 * QVARS] = {0.};
    double *arr_dBU_po_dIndL = new double[QVARS];
    double *arrDelX = new double [QVARS];
    MtrxMinusMatrx(arrXZv, StatSolutionParams.marrStatPhVect,1, QVARS, arrDelX);

    calcMtrxB(arrB);
    double arrFx[QVARS] = {0.};
    calcFx( Environment, arrXZv, arrFx,VAlMomOut);
    double arrU[2] = {0.}, arrDelU[2] = {0.}, arrt0[QVARS] = {0.};
    MtrxMultMatrx(StatSolutionParams.marrGears,2, QVARS,   arrDelX, 1, arrDelU);
    MtrxSumMatrx(arrDelU, StatSolutionParams.marrStatU,1, 2, arrU) ;
    MtrxMultMatrx(arrB, QVARS, 2,  arrU, 1, arrt0);
    MtrxSumMatrx(arrFx, arrt0,1, QVARS, arrFi) ;

    calc_dF_po_dAlf(arrXZv, mtrx_dFi_po_dAlf, VAlMomOut);

    calc_dBU_po_dIndL(arrU,arr_dBU_po_dIndL);
    for(int j =0; j < QVARS; ++j)
    {
      mtrx_dFi_po_dAlf[ j * QPARAMS] +=   arr_dBU_po_dIndL[j];
    }


  delete []arrDelX;
  }

//--------------------------------------
//нахождение сатционарного решения при заданном векторе управлений
bool QCtrl::calcStableSolution(const double *arrObjective
      , const double valUd, const double valUq, double*arrStatinaryPhVect,const double VAlMomOut)
{
    bool breturn = false;
    // решение нелинейного уравнения относительно Omega
    double valOm = 0.;
    double valL = mElectMotor.getL();
    double a = mElectMotor.mResist * valUq/ valL;
    double b = -mElectMotor.mResist * mElectMotor.mZp * mElectMotor.mPsi_f/ valL - mElectMotor.mZp *valUd;
    double c = mElectMotor.mZp * mElectMotor.mZp * valL;
    double d = mElectMotor.mResist  * mElectMotor.mResist  / valL;
    TEnvironment Environment(0.,0.,0.);
    for (int i = 0; i < 100; ++i)
    {
       double valIqu = (a + b * valOm)/ (c * valOm * valOm + d);
       double valFi = 3./2. * mElectMotor.mZp * mElectMotor.mPsi_f * valIqu
               + mpLoad->calcAirResistMom( Environment,0., valOm) + VAlMomOut;
       double valDerivFi = 3./2. * mElectMotor.mZp * mElectMotor.mPsi_f
               *( b/ (c * valOm * valOm + d) - 2. * c * valOm *(a + b * valOm)/ (c * valOm * valOm + d)/ (c * valOm * valOm + d));


       double valdMa = mpLoad-> calc_dAirResistMom_po_dOmega( Environment,0., valOm);
        valDerivFi += valdMa;
       double delta = valFi/valDerivFi;

        valOm = valOm - 0.25 * delta;
       if (fabs(delta) < 0.0000001)
       {
           breturn = true;
           break;
       }

    }
    if (!breturn )
    {
        return false;
    }
    arrStatinaryPhVect[0] = valOm;
    arrStatinaryPhVect[2] = (a + b * valOm)/ (c * valOm * valOm + d);
    arrStatinaryPhVect[1] =  (valL * mElectMotor.mZp * valOm * arrStatinaryPhVect[2] + valUd) / mElectMotor.mResist ;
    arrStatinaryPhVect[3] = arrObjective[0];
}

//--------------------------------------------------
void QCtrl::correctMomOut(const double VAlTettaZv
                          ,const double VAlDispTetta    , const double VAlW, const double VAlTcur)
{

}
//-------------------------------
//
void QCtrl::doExtrapolation(const double VAlDeltaT,const double VAlIntegrStep0
             , double *arrUPrev,const double VAlMomOut)
{
    TEnvironment Environment0(0.,0.,0.);
    double arrEstPhVectExtrp[QVARS] = {0.};

    dragPhVect(Environment0,mFiltr.marrCurEst,arrUPrev
                       , VAlDeltaT, VAlIntegrStep0, arrEstPhVectExtrp, VAlMomOut);
    memcpy(mFiltr.marrCurEst,arrEstPhVectExtrp, QVARS * sizeof(double));

}

//----------------------------------------
double QCtrl::calcEstMomOut()
{
    return 0.;
}
//---------------------------------
int QCtrl::createInputDataReport_Ctrl(wchar_t*FileName, const bool bHeader)

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
fprintf(fw,"     Система управления приводом\n");
fprintf(fw,"  собственные числа характеристического уравнения:\n");
fprintf(fw,"  lambda[0] = %5.1f + i* %5.1f\n",mCmpArrLamb[0].m_Re, mCmpArrLamb[0].m_Im);
fprintf(fw,"  lambda[1] = %5.1f + i* %5.1f\n",mCmpArrLamb[1].m_Re, mCmpArrLamb[1].m_Im);
fprintf(fw,"  lambda[2] = %5.1f + i* %5.1f\n",mCmpArrLamb[2].m_Re, mCmpArrLamb[2].m_Im);
fprintf(fw,"  lambda[3] = %5.1f + i* %5.1f\n",mCmpArrLamb[3].m_Re, mCmpArrLamb[3].m_Im);
fprintf(fw,"***************************************\n");
fclose(fw);

}
//------------------------------------------------------------------------
int QCtrl::createInputDataReport_Ctrl_Inherited(wchar_t*FileName, const bool bHeader)
{

}





