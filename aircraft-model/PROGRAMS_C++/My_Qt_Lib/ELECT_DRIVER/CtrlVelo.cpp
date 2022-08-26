#include "CtrlVelo.h"
#include <math.h>
#include  <string.h>
#include "MatrixProccess.h"
#include "Equations.h"
#include "Comp.h"
#include "MeasurmentImitator.h"


QCtrlVelo::QCtrlVelo():QCtrl()
{

}

// конструктор копирования
 QCtrlVelo::QCtrlVelo(const QCtrlVelo &R):QCtrl( R)
 {

 }

 // оператор присваивания
  QCtrlVelo  &QCtrlVelo::operator=( const QCtrlVelo  &R)
 {
      if(this == &R)
      {
          return *this;
      }
     QCtrl:: operator= (R);

     return *this ;
 }
//-----------------------------------

// парам конструктор
  // параметрический конструктор
QCtrlVelo:: QCtrlVelo( const QElectMotor ElectMotor ,  QLoad*  Load
                       , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
                       ,const double VAlOmegaStatBegin , const double T0,const double TCur, const double valh,TComp *pCmpArrLamb)
      :QCtrl(  ElectMotor ,  Load, arrSpreadParams, MomOut,VAlTettaBegin
               , VAlOmegaStatBegin, T0 ,TCur, valh,pCmpArrLamb)
  {

  }

//------------------------------------------------------------------
//--------------------------------------------------------
void QCtrlVelo::CalcEigenValues( const double VAlTettaStat, const double VAlOmegaStat
                                 ,const double *arrC,const double VAlMomOut, TComp* cmparrEigenValues)
{

    double arrFGr[9] = {0.};
    calc_df_po_dx_plus_BC_Velo( VAlOmegaStat,arrC,VAlMomOut, arrFGr);
    double arrLamb[3] = {0.};
    switch(CalcProper_Numbers_R3(arrFGr, arrLamb))
    {
    case 1:
        cmparrEigenValues[0].m_Re =arrLamb[0];
        cmparrEigenValues[0].m_Im = 0.;
        cmparrEigenValues[1].m_Re =  arrLamb[1];
        cmparrEigenValues[1].m_Im =  arrLamb[2];
        cmparrEigenValues[2].m_Re =  arrLamb[1];
        cmparrEigenValues[2].m_Im = -arrLamb[2];
        break;

    default:
        cmparrEigenValues[0].m_Re =arrLamb[0];
        cmparrEigenValues[0].m_Im = 0.;
        cmparrEigenValues[1].m_Re =  arrLamb[1];
        cmparrEigenValues[1].m_Im = 0.;
        cmparrEigenValues[2].m_Re =  arrLamb[2];
        cmparrEigenValues[2].m_Im = -0.;
        break;
    }

  //  double x0 =  calcDet_A_minus_LambdaE_D3( arrFGr, -8.);
   // double x1 =  calcDet_A_minus_LambdaE_D3( arrFGr, -9.);
   // double x2 =  calcDet_A_minus_LambdaE_D3( arrFGr, -10.);

}




//---------------------------------------------------------------
// вычисление матрицы передаточных чисел
//INPUT:
//arrObjective[2] - вектор подожения и скорости в состоянии равновесия (нужна только скоросчть)
//CmpArrLamb [3] - корни характерист уравнения
// OUTPUT:
// arrC0[8] - матрица передаточных чисел

void QCtrlVelo::calcGearMtrx(const double *arrObjective,  TComp* CmpArrLamb,const double VAlMomOut, double* arrC0)
{
    double  valUd =0., valUq =0., valDelFi =0.;
    double arrStatinaryPhVect[QVARS] = {0.};
     calcStationarySolution( arrObjective, arrStatinaryPhVect, &valUd, &valUq, &valDelFi, VAlMomOut);

    double arrStatinaryU[2] = {0.};
    arrStatinaryU[0] =valUd;
    arrStatinaryU[1] =valUq;

    int quantControlledVars = 3;

     double *arr_dF_po_dx = new  double[quantControlledVars * quantControlledVars] ;
     double *arr_dF_po_dW = new  double[quantControlledVars * 2] ;
    memset(arr_dF_po_dx, 0, quantControlledVars * quantControlledVars * sizeof( double));
    memset(arr_dF_po_dW, 0, 2* quantControlledVars * sizeof( double));

   fill_df_po_px_and_mtrxB_Velocity(arrStatinaryPhVect,arr_dF_po_dx, arr_dF_po_dW);
   double arrC[6] = {0.};
   arrC[0] = -arr_dF_po_dx[3] /mElectMotor.mInvL;
   arrC[2] = -arr_dF_po_dx[5]/mElectMotor.mInvL;
   arrC[4] = - arr_dF_po_dx[7]/mElectMotor.mInvL;

   arrC[1] =  (CmpArrLamb[0].m_Re -arr_dF_po_dx[4])/mElectMotor.mInvL;
   double x3 = CmpArrLamb[1].m_Re + CmpArrLamb[2].m_Re -arr_dF_po_dx[0];
   arrC[5] =  (x3    -arr_dF_po_dx[8])/mElectMotor.mInvL;
   double x2 = (arr_dF_po_dx[0]* x3 - CmpArrLamb[1].m_Re *CmpArrLamb[2].m_Re + CmpArrLamb[1].m_Im * CmpArrLamb[2].m_Im)
                /arr_dF_po_dx[2];
   arrC[3] =  (x2 -arr_dF_po_dx[6])/mElectMotor.mInvL;

   delete []arr_dF_po_dx;
   delete []arr_dF_po_dW;

   memset (arrC0, 0, 2 * QVARS  * sizeof(double)) ;
   memcpy(arrC0, arrC, 3 * sizeof(double));
   memcpy(&(arrC0[QVARS]), &(arrC[3]), 3 * sizeof(double));



}


// вычисление матрицы частных произыодных подсиситемы уравнен6ий по скорости
void QCtrlVelo::fill_df_po_px_and_mtrxB_Velocity(const double *arrStatinaryPhVect
                 ,double *arr_dF_po_dx,double * arr_dF_po_dW)
{
    double arr_dF_po_dx_Big[QVARS * QVARS] = {0.}, arr_dF_po_dW_Big[2 * QVARS] = {0};
    fill_df_po_px_and_mtrxB_(arrStatinaryPhVect
                     ,arr_dF_po_dx_Big,arr_dF_po_dW_Big);
    memcpy(arr_dF_po_dW, arr_dF_po_dW_Big, 6 * sizeof(double));
    memcpy(&(arr_dF_po_dx[0]),&(arr_dF_po_dx_Big[0]), 3 * sizeof(double));
    memcpy(&(arr_dF_po_dx[3]),&(arr_dF_po_dx_Big[QVARS]), 3 * sizeof(double));
    memcpy(&(arr_dF_po_dx[6]),&(arr_dF_po_dx_Big[2 * QVARS]), 3 * sizeof(double));   

}
//---------------------------------------------------

//матрица частных производных правой части по фазовым переменным
// для первых 3-х уравнений _ угловая скорость, Id, Iq
void QCtrlVelo::calc_df_po_dx_plus_BC_Velo(const double VAlOmegaStat
                      ,const double *arrC,const double VAlMomOut, double *arrFGr)
{
double arrStatPhVect[QVARS] = {0.}, arrUstat[2] = {0.};

double valIq, valUd=0.,  valUq=0.,  valDelFi=0., arrObjVect[2 ] ={0.}, arrStatinaryPhVect[QVARS] = {0.};
arrObjVect[0] = 0.;
arrObjVect[1] = VAlOmegaStat;

calcStationarySolution( arrObjVect, arrStatPhVect,&valUd, &valUq, &valDelFi, VAlMomOut);

arrUstat[0] = valUd;
arrUstat[1] = valUq;
double arr_dF_po_dx[9] = {0.},arr_dF_po_dW[6] = {0.}, arrT[9] = {0.}, arrCCC[6] ={0.};

fill_df_po_px_and_mtrxB_Velocity(arrStatPhVect
                       ,arr_dF_po_dx,arr_dF_po_dW);


for (int i = 0; i < 2; ++i )
    for (int j =0; j < 3; j++)
{
  arrCCC[3 * i + j] = arrC[QVARS * i + j];
}

MtrxMultMatrx(arr_dF_po_dW,3, 2, arrCCC,3, arrT) ;
MtrxSumMatrx(arrT, arr_dF_po_dx, 3, 3, arrFGr) ;

}
/*
//----------------------------------------
//вычисление фазоового вектора и вектора напряжкений в состоянии равновесия
// INPUT:
//arrObjective[2] -вектор положения и скорости в состоянии равновесия
// OUTPUT:
//arrStatinaryPhVect[4] - фазовый вектор (угл скорость, Id, Iq, угл положение)
//valUd, valUq - напряжения
// *pvalDelFi - угол опрежения
//
void QCtrlVelo::calcStationarySolution(const double *arrObjective, double*arrStatinaryPhVect,double *valUd, double*valUq, double *pvalDelFi)
{

    TEnvironment Environment0;
    QSegment SegmentCur;
    double arrGSKOrtZero[3] = {0.};
     arrGSKOrtZero[0] = 1.;
     TPlane PlaneCur;
    arrStatinaryPhVect [0] = arrObjective[1];
    arrStatinaryPhVect [1] = 0.;// Id

    // Iq:
    double valCurrentIq = arrStatinaryPhVect [2] = 2./ 3. *
       (  mpLoad->calcAirResistMom(Environment0,PlaneCur,SegmentCur, 0., marrObjective[1])
            +  mElectMotor.mMomResidual * sign_(arrObjective[1]))
             /  mElectMotor.mZp/ mElectMotor.mPsi_f;
   arrStatinaryPhVect [3] = arrObjective[0];
   *valUd = - mElectMotor.mInductL* arrObjective[1] * valCurrentIq *  mElectMotor.mZp;
   *valUq =  mElectMotor.mZp *  mElectMotor.mPsi_f * arrObjective[1]
           +  mElectMotor.mResist * valCurrentIq;
   *pvalDelFi = atan2(-(*valUd),(*valUq));
}
*/
//----------------------------------------
/*
// вычисление матрицы частных произыодных подсиситемы уравнен6ий по скорости
void QCtrlVelo::fill_df_po_px_and_mtrxB_Velocity(const double *arrStatinaryPhVect
                 ,double *arr_dF_po_dx,double * arr_dF_po_dW)
{
  const double JPayLoad =  mElectMotor.mJ0 +  mpLoad->mJPayLoad;
  memset(arr_dF_po_dW, 0, 6 * sizeof(double));
  arr_dF_po_dW[2] = arr_dF_po_dW[5] = 1./ mElectMotor.mInductL;

  memset(arr_dF_po_dx, 0, 9 * sizeof(double));
  //arr_dF_po_dx[0] =  mpLoad->calc_dAirResistMom_po_dOmega(arrStatinaryPhVect[0])/ JPayLoad;
  TEnvironment Environment(0.,0.,0.);
  TPlane Plane;
  QSegment Sgm;
  arr_dF_po_dx[0] =  mpLoad
       ->calc_dAirResistMom_po_dOmega( Environment, Plane, Sgm, arrStatinaryPhVect[3], arrStatinaryPhVect[0])/ JPayLoad;

  arr_dF_po_dx[1] = 0.;
  arr_dF_po_dx[2] = 3./2. *  mElectMotor.mZp *  mElectMotor.mPsi_f / JPayLoad;
  arr_dF_po_dx[3] =  mElectMotor.mZp *arrStatinaryPhVect[2];
  arr_dF_po_dx[4] = -  mElectMotor.mResist / mElectMotor.mInductL ;
  arr_dF_po_dx[5] =  mElectMotor.mZp *arrStatinaryPhVect[0];
  arr_dF_po_dx[6] = - mElectMotor.mZp *(arrStatinaryPhVect[1] +  mElectMotor.mPsi_f/  mElectMotor.mInductL );
  arr_dF_po_dx[7] = - mElectMotor.mZp *arrStatinaryPhVect[0];
  arr_dF_po_dx[8] = -  mElectMotor.mResist /  mElectMotor.mInductL ;

}
*/

void QCtrlVelo::calcCurU(const double *arrObjective,const double VAlTObjective,const double VAlMomOut, double *arrU)
{
    double arrDelta[QVARS] = {0.}, arrdelU[2] = {0.};
    QStatSolutionParams StatSolutionParams;
    calcStationaryParams(arrObjective, mCmpArrLamb, &StatSolutionParams, VAlMomOut);
    MtrxMinusMatrx(mFiltr.marrCurEst, StatSolutionParams.marrStatPhVect,1, QVARS, arrDelta);
    MtrxMultMatrx(StatSolutionParams.marrGears,2, QVARS, arrDelta,1, arrdelU) ;
    MtrxSumMatrx(StatSolutionParams.marrStatU, arrdelU,1, 2, arrU) ;
    if(NormVect2(arrU) >  mElectMotor.mUMax)
    {
        MatrxMultScalar(arrU, 1, 2,  mElectMotor.mUMax/ NormVect2(arrU),arrU);
    }

}

//------------------------------------------
//void QCtrlVelo::calcStationaryParams(const double *arrObjective,  TComp* CmpArrLamb, QStatSolutionParams* pQStatSolutionParams)
//{

//}

