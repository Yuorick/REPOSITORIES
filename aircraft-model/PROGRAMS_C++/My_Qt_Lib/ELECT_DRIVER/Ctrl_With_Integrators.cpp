#include "Ctrl_With_Integrators.h"
#include <string.h>
#include "MatrixProccess.h"


QCtrl_With_Integrators::QCtrl_With_Integrators():QCtrlFollow( )
{
    memset(marrIntegrators, 0, QVARS * sizeof(double));
    mTIntegrCur = 0.;
    mIntegratorTime = 0.;
}


// конструктор копирования
 QCtrl_With_Integrators::QCtrl_With_Integrators(const QCtrl_With_Integrators &R):QCtrlFollow( R)
 {

     memcpy(marrIntegrators, R.marrIntegrators, QVARS * sizeof(double));
     mTIntegrCur = R.mTIntegrCur;
     mIntegratorTime = R.mIntegratorTime;

 }

// оператор присваивания
 QCtrl_With_Integrators  &QCtrl_With_Integrators::operator=( const QCtrl_With_Integrators  &R)
{
     if(this == &R)
     {
         return *this;
     }
     memcpy(marrIntegrators, R.marrIntegrators, QVARS * sizeof(double));
     mTIntegrCur = R.mTIntegrCur;
     mIntegratorTime = R.mIntegratorTime;


    QCtrlFollow:: operator= (R);

    return *this ;
}
//-----------------------------------

// парам конструктор
 // параметрический конструктор
 QCtrl_With_Integrators::QCtrl_With_Integrators( const QElectMotor ElectMotor ,  QLoad*  Load
             , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
             ,const double VAlOmegaStatBegin , const double T0,const double TCur, const double valh
               , TComp *pCmpArrLamb,const double IntegratorTime )
     :QCtrlFollow (ElectMotor ,  Load
                    ,arrSpreadParams, MomOut,VAlTettaBegin
                    ,VAlOmegaStatBegin, T0 ,TCur, valh,pCmpArrLamb)
 {
  memset(marrIntegrators, 0, QVARS * sizeof(double));
  mTIntegrCur = TCur;
  mIntegratorTime = IntegratorTime;
 }

//------------------------------------------------------------------
 void QCtrl_With_Integrators::calcCurU(const double *arrObjective,const double VAlTObjective
          ,const double VAlMomOut, double *arrU)
 {

     // это виртуальная функция, позволяющая производить периодические вычисления параметров Mout и R
     // с целью их коррекции в процессе функционирования.
     // эти вычисления реализованы в наследуемом классе QCtrlFollow:QCtrlDnmkCorrect_M_R

     recalcIntegrators(arrObjective, VAlMomOut);
     fncRenewParams_M_R(arrObjective,VAlMomOut);
     ///


     double arrDelta[QVARS] = {0.}, arrdelU[2] = {0.};
     double arrTargPhVect[QVARS] ={0.};

     QStatSolutionParams StatSolutionParams;
     calcStationaryParams(arrObjective,  mCmpArrLamb, &StatSolutionParams, VAlMomOut);
     memcpy(arrTargPhVect, StatSolutionParams.marrStatPhVect, QVARS * sizeof(double));

     arrTargPhVect[3] = arrObjective[0];
     MtrxMinusMatrx(mFiltr.marrCurEst, arrTargPhVect,1,  QVARS, arrDelta);

     MtrxMultMatrx(StatSolutionParams.marrGears,2,  QVARS, arrDelta,1, arrdelU) ;
     MtrxSumMatrx(StatSolutionParams.marrStatU, arrdelU,1, 2, arrU) ;
     if(NormVect2(arrU) > mElectMotor.mUMax)
     {
         MatrxMultScalar(arrU, 1, 2, mElectMotor.mUMax/ NormVect2(arrU),arrU);
     }

     // это виртуальная функция, позволяющая производить периодическую коррекция
     // вектора управлениий arrU
     // с целью  коррекции ошибки по Tetta
     // эти вычисления реализованы в наследуемом классе QCtrlFollow:QCtrlDnmkCorrect_W
     double arr[2] = {0.};
     memcpy(arr, arrObjective, 2 * sizeof(double));
     addW(arr,VAlMomOut,arrU);
     ///
 }

 //---------------------------------------
 void QCtrl_With_Integrators::recalcIntegrators(const double *arrObjective,const double VAlMomOut)
 {
     // пересчет интеграторов
     double alf = 1. / mIntegratorTime;
     double arrTargPhVect[QVARS] = {0.};
     QStatSolutionParams StatSolutionParams;
     calcStationaryParams(arrObjective,  mCmpArrLamb, &StatSolutionParams,VAlMomOut);
     memcpy(arrTargPhVect, StatSolutionParams.marrStatPhVect, QVARS * sizeof(double));
     arrTargPhVect[3] = arrObjective[0] + (mTCur - mT0) * arrObjective[1];
     for (int j =0; j < QVARS ; ++j)
     {
        double valDelta = mFiltr.marrCurEst[j] - arrTargPhVect[j];
        marrIntegrators[j] = marrIntegrators[j]
     +   (- alf *  marrIntegrators[j]   + valDelta) *  (mTCur -mTIntegrCur) ;
    }

    mTIntegrCur = mTCur;


 }

 //коррекция вектора управлений
 void QCtrl_With_Integrators::addW(const double *arrObjective,const double VAlMomOut,double *arrU)
 {

 }

 void QCtrl_With_Integrators::fncRenewParams_M_R(const double *arrObjective,const double VAlMomOut)
 {

 }



