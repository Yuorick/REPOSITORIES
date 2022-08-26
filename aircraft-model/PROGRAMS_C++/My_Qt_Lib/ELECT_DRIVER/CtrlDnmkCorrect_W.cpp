#include "CtrlDnmkCorrect_W.h"
#include <string.h>
#include <math.h>
#include "MatrixProccess.h"

QCtrlDnmkCorrect_W::QCtrlDnmkCorrect_W():QCtrl_With_Integrators( )
{
    mTPeriodRenew_dU = 0.;
    mTLastRenew_dU = 0.;
    memset(marrW, 0, 2 * sizeof(double));
}
//-----------------------------------------
// конструктор копирования
 QCtrlDnmkCorrect_W::QCtrlDnmkCorrect_W(const QCtrlDnmkCorrect_W &R):QCtrl_With_Integrators( R)
 {
  mTPeriodRenew_dU = R.mTPeriodRenew_dU;
  mTLastRenew_dU = R.mTLastRenew_dU;
  memcpy(marrW, R.marrW, 2 * sizeof(double));
 }
// оператор присваивания
 QCtrlDnmkCorrect_W  &QCtrlDnmkCorrect_W::operator=( const QCtrlDnmkCorrect_W  &R)
{
     if(this == &R)
     {
         return *this;
     }

     mTPeriodRenew_dU = R.mTPeriodRenew_dU;
     mTLastRenew_dU = R.mTLastRenew_dU;
     memcpy(marrW, R.marrW, 2 * sizeof(double));
    QCtrl_With_Integrators:: operator= (R);

    return *this ;
}
//-----------------------------------

// парам конструктор
 // параметрический конструктор
QCtrlDnmkCorrect_W::QCtrlDnmkCorrect_W( const QElectMotor ElectMotor ,  QLoad*  Load
                      , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
                      ,const double VAlOmegaStatBegin , const double T0,const double TCur,const double valh,TComp *pCmpArrLamb
                      , const double IntegratorTime,const double TPeriodRenew_dU,const double TLastRenew_dU)
     :QCtrl_With_Integrators ( ElectMotor ,   Load,arrSpreadParams,  MomOut,  VAlTettaBegin
                   , VAlOmegaStatBegin ,T0,  TCur, valh,pCmpArrLamb, IntegratorTime)
 {
    mTPeriodRenew_dU = TPeriodRenew_dU;
    mTLastRenew_dU = TCur;
    memset(marrW, 0, 2 * sizeof(double));
 }

//------------------------------------------------------------------


//коррекция вектора управлений
void QCtrlDnmkCorrect_W::addW(const double *arrObjective,const double VAlMomOut,double *arrU)
{
    arrU[0] += marrW[0];
    arrU[1] += marrW[1];
    if(!(( mTCur - mTLastRenew_dU) >= mTPeriodRenew_dU))
    {
        return;
    }

    calcCorrectingDeltaW(arrObjective,VAlMomOut,marrW);

    mTLastRenew_dU= mTCur;

}

//---------------------------------
void QCtrlDnmkCorrect_W::calcCorrectingDeltaW(const double *arrObjective,const double VAlMomOut, double *arrW)
{

  double arr_dF_po_dx[QVARS *QVARS] = {0.},arr_dF_po_dW[2 * QVARS] = {0.}
           , arrT[QVARS * QVARS] = {0.}, arrCC[2 *QVARS] ={0.}, arrA[QVARS *QVARS] = {0.}
         , arrAInv[QVARS *QVARS] = {0.},arr_d[2 *QVARS] ={0.},arr_PhVectMiddle[QVARS] ={0.};
  QStatSolutionParams StatSolutionParams;
  calcStationaryParams(arrObjective, mCmpArrLamb, &StatSolutionParams,VAlMomOut);
  MtrxSumMatrx( marrIntegrators,  StatSolutionParams.marrStatPhVect,QVARS, 1, arr_PhVectMiddle) ;
  arr_PhVectMiddle[3] += arrObjective[1] * ( mTCur -  mT0);
   fill_df_po_px_and_mtrxB_(arr_PhVectMiddle
                         ,arr_dF_po_dx,arr_dF_po_dW);
  memcpy(arrCC,  StatSolutionParams.marrGears, 2 * QVARS * sizeof(double));
  MtrxMultMatrx(arr_dF_po_dW,QVARS, 2, arrCC, QVARS, arrT) ;
  MtrxSumMatrx(arrT, arr_dF_po_dx, QVARS, QVARS, arrA) ;
  bool brez = InverseMtrx(arrA, QVARS, arrAInv) ;


  MtrxMultMatrx(arrAInv, QVARS, QVARS, arr_dF_po_dW, 2, arr_d) ;
  double arra[4] ={0.}, arrb[2]  = {0.};
  arra[0] = -arr_d[0];
  arra[1] = -arr_d[1];
  arra[2] = -arr_d[6];
  arra[3] = -arr_d[7];

  //double arrTarg[2] = {0.};
 // arrTarg[0] =   marrObjective[0] +  marrObjective[1] * ( mTCur - mT0);
 // arrTarg[1] =   marrObjective[1];

  arrb[0] = - marrIntegrators[0];
  arrb[1] = - marrIntegrators[3] ;

  if(! SolvLinEq2(arra, arrb, arrW))
  {
    if (fabs(arra[3] > 0.00001))
    {
       arrW[0] = 0;
       arrW[1] = arrb[1] / arra[3];
    }
  }
  int uuu=0;


}
