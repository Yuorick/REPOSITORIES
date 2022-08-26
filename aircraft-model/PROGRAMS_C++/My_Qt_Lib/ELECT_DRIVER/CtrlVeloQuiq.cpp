#include "CtrlVeloQuiq.h"
#include "LinOptCtrlSyst_dim21.h"
#include "MatrixProccess.h"
#include  <string.h>

QCtrlVeloQuiq::QCtrlVeloQuiq():QCtrlVelo()
{
  mMaxU = 0.;
  mVeloQuiqData = QVeloQuiqData();
}

// конструктор копирования
 QCtrlVeloQuiq::QCtrlVeloQuiq(const QCtrlVeloQuiq &R):QCtrlVelo( R)
 {
     mMaxU = R.mMaxU;
     mVeloQuiqData = R.mVeloQuiqData;
 }

 // оператор присваивания
  QCtrlVeloQuiq  &QCtrlVeloQuiq::operator=( const QCtrlVeloQuiq  &R)
 {
      if(this == &R)
      {
          return *this;
      }
     QCtrlVelo:: operator= (R);
     mMaxU = R.mMaxU;
     mVeloQuiqData = R.mVeloQuiqData;
     return *this ;
 }
//-----------------------------------


  // парам конструктор 2
     // параметрический конструктор
   QCtrlVeloQuiq:: QCtrlVeloQuiq( const QElectMotor ElectMotor ,  QLoad*  Load
              , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
                                  ,const double VAlOmegaStatBegin , const double T0,const double TCur, const double *arrObjective
                                  ,const double MaxU,const double valh,TComp *pCmpArrLamb)
                 :QCtrlVelo(  ElectMotor ,  Load,  arrSpreadParams, MomOut, VAlTettaBegin
                          , VAlOmegaStatBegin , T0, TCur,valh,pCmpArrLamb)

     {
     mMaxU = MaxU;
     memcpy(mVeloQuiqData.marrObjective, arrObjective, 2 * sizeof(double));
     calcSwithOverTimes(&mVeloQuiqData, MomOut);
     }

  //---------------------------------------
void QCtrlVeloQuiq::calcCurU(const double *arrObjective,const double VAlTObjective,const double VAlMomOut, double *arrU)

{
    double valU0 = 0.;
    if(mTCur < mVeloQuiqData.mEndT)
     {

     valU0 = -mMaxU * ((double)(mVeloQuiqData.misignum));
     if(mTCur < mVeloQuiqData.mOverSwitchT)
     {
       valU0 = mMaxU * ((double)(mVeloQuiqData.misignum));
     }
    arrU[1] = valU0;
    arrU[0] =- mFiltr.marrCurEst[0]* mFiltr.marrCurEst[2]
            * mElectMotor.mZp / mElectMotor.mInvL   ;


    }
     else
     {
        double arrDelta[QVARS] = {0.}, arrdelU[2] = {0.};
        QStatSolutionParams StatSolutionParams;
        calcStationaryParams(arrObjective,  mCmpArrLamb, &StatSolutionParams,VAlMomOut);
        MtrxMinusMatrx(mFiltr.marrCurEst, StatSolutionParams.marrStatPhVect,1, QVARS, arrDelta);
        MtrxMultMatrx(StatSolutionParams.marrGears,2, QVARS, arrDelta,1, arrdelU) ;
        MtrxSumMatrx(StatSolutionParams.marrStatU, arrdelU,1, 2, arrU) ;
        if(NormVect2(arrU) >  mElectMotor.mUMax)
        {
            MatrxMultScalar(arrU, 1, 2,  mElectMotor.mUMax/ NormVect2(arrU),arrU);
        }

     }
}
//---------------------------------------

// вычисление моментов переключения
// INPUT:
// VAlUqu0  - максимальная величина управления на переходном участке
// VAlOmegaStat - желаемая величина угловой скорости (на правлм конце)
// arrObjective[2]
//OUTPUT:
//*pvalt0 - момент переключения
//*pvalt1 - конечный момент времени
//*pisignum0 - знак управления на первом участке [0;*pvalt0 ]
bool QCtrlVeloQuiq:: calcSwithOverTimes(QVeloQuiqData *pVeloQuiqData,const double VAlMomOut)
{
double valIq = 0.,valUd = 0., valUqu = 0., valDelFi = 0.;
double arrStatPhVect[QVARS] = {0.}, arrUstat[2] = {0.};
calcStationarySolution(pVeloQuiqData->marrObjective, arrStatPhVect,arrUstat, &(arrUstat[1]), &valDelFi,VAlMomOut);

// формирование члена класса оптимальной сиситемы 2-го порядка
// формирование матрицы A
double arrAZv[4] = {0.};
const double JInvPayLoad = getInvSumJ0();

arrAZv[0] = -(*(mpLoad)).mCv * JInvPayLoad ;
arrAZv[1] = 3./2. * mElectMotor.mZp* mElectMotor.mPsi_f * JInvPayLoad ;
arrAZv[2] = -mElectMotor.mZp *mElectMotor.mPsi_f *mElectMotor.mInvL ;
arrAZv[3] = -mElectMotor.mResist*mElectMotor.mInvL ;

// формирование матрицы B
double arrBZv[2] = {0.};
 arrBZv[1] = mElectMotor.mInvL;

 // формирование матрицы C
 double arrCZv[2] = {0.};
 arrCZv[0] = 0.;// -mElectMotor.getMomResidual()/JPayLoad;
 ///


 double arrxBegin[2] = {0.},arrx[2] = {0.}, arrxEnd[2] ={0.};
 const double T0 = 0.;
  const double TCur = 0.;
  arrx [0] = arrxBegin[0] = mFiltr.marrCurEst[0];
  arrx [1] = arrxBegin[1] = mFiltr.marrCurEst[2];
 QLinOptCtrlSyst_dim21 OptCtrlSyst(arrAZv,arrBZv, arrCZv,T0, arrxBegin,TCur,arrx);
 ///



// конечное состояние фазового вектьора сиситемы
 arrxEnd[0] = arrStatPhVect[0];
 arrxEnd[1] = arrStatPhVect[2];
 ///

//поиск линейного оптимального быстродействия
mVeloQuiqData.mOverSwitchT = 2;
OptCtrlSyst.calcSwitchPointTimes(arrxEnd,mMaxU,&((*pVeloQuiqData).mOverSwitchT)
                                 , &((*pVeloQuiqData).mEndT), &((*pVeloQuiqData).misignum));
}



//-------------------------------------------------------------
//-------------------------------------------------------------
//-----    // CLASS QVeloQuiqData    ---------------------------
//-------------------------------------------------------------
//-------------------------------------------------------------
//-------------------------------------------------------------

QVeloQuiqData::QVeloQuiqData()
{
    mOverSwitchT =0.;
    mEndT = 0.;
    misignum = 0;
    memset(marrObjective, 0, 2 * sizeof(double));
}

QVeloQuiqData::QVeloQuiqData (const QVeloQuiqData &R)
{
    mOverSwitchT = R.mOverSwitchT;
    mEndT = R.mEndT;
    misignum = R.misignum;
    memcpy(marrObjective, R.marrObjective,2 * sizeof(double));
}

QVeloQuiqData  &QVeloQuiqData::operator=( const QVeloQuiqData  &R)
{
    mOverSwitchT = R.mOverSwitchT;
    mEndT = R.mEndT;
    misignum = R.misignum;
    memcpy(marrObjective, R.marrObjective,2 * sizeof(double));

return *this;
}

/*
QVeloQuiqData:: QVeloQuiqData(const double  OverSwitchT, const double  EndT,const int isignum)
{
  mOverSwitchT = OverSwitchT;
  mEndT = EndT;
  misignum = isignum;
}
*/
