#include "CtrlDnmkCorrect_M_R.h"
#include <string.h>
#include "MatrixProccess.h"

QCtrlDnmkCorrect_M_R::QCtrlDnmkCorrect_M_R():QCtrl_With_Integrators( )
{
   // memset(mCmpArrLamb, 0, 4 * sizeof(TComp));
    mTPeriodParamsRenew = 0.;
    mTLastRenew = 0.;
}


// конструктор копирования
 QCtrlDnmkCorrect_M_R::QCtrlDnmkCorrect_M_R(const QCtrlDnmkCorrect_M_R &R):QCtrl_With_Integrators( R)
 {

    // memcpy(mCmpArrLamb, R.mCmpArrLamb, 4 * sizeof(TComp));
     mTPeriodParamsRenew = R.mTPeriodParamsRenew;
     mTLastRenew = R.mTLastRenew;
 }
// оператор присваивания
 QCtrlDnmkCorrect_M_R  &QCtrlDnmkCorrect_M_R::operator=( const QCtrlDnmkCorrect_M_R  &R)
{
     if(this == &R)
     {
         return *this;
     }
    // memcpy(mCmpArrLamb, R.mCmpArrLamb, 4 * sizeof(TComp));
     mTPeriodParamsRenew = R.mTPeriodParamsRenew;
     mTLastRenew = R.mTLastRenew;

    QCtrl_With_Integrators::operator= (R);

    return *this ;
}
//-----------------------------------

// парам конструктор
 // параметрический конструктор
QCtrlDnmkCorrect_M_R::QCtrlDnmkCorrect_M_R(  const QElectMotor ElectMotor ,  QLoad*  Load
     , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
     ,const double VAlOmegaStatBegin , const double T0,const double TCur,const double valh,TComp *pCmpArrLamb
     , const double IntegratorTime, const double TPeriodParamsRenew
     , const double TLastRenew )
     :QCtrl_With_Integrators ( ElectMotor ,   Load,arrSpreadParams,  MomOut,  VAlTettaBegin
                   , VAlOmegaStatBegin ,T0,  TCur, valh,pCmpArrLamb, IntegratorTime)
 {

    mTPeriodParamsRenew = TPeriodParamsRenew;
    mTLastRenew = TLastRenew;

 }

//------------------------------------------------------------------
// коррекция параметров Mout и R по данным интеграторов
void QCtrlDnmkCorrect_M_R::fncRenewParams_M_R(const double *arrObjective,const double VAlMomOut)
{

    if(!(( mTCur - mTLastRenew) >= mTPeriodParamsRenew))
    {
        return;
    }
        const int QUantParams = 2;
        double arrParams[3] = {0.};
        calcNewParams(arrObjective, VAlMomOut,QUantParams, arrParams);
        setParam(6,arrParams[0]);
        setParam(1,arrParams[1]);
        setParam(2,arrParams[2]);
        QStatSolutionParams StatSolutionParams;


        calcStationaryParams(arrObjective,   mCmpArrLamb, &StatSolutionParams, VAlMomOut);
       // mStatSolutionParams = StatSolutionParams;
        mTLastRenew = mTCur;

}


void QCtrlDnmkCorrect_M_R::calcNewParams(const double *arrObjective,const double VAlMomOut
                                         ,const int QUantParams, double *arrParams)
{
    double arrMeanPhVect[QVARS] = {0.};
       createMeanPhVect(arrObjective, VAlMomOut,arrMeanPhVect);
        double  arrbWave[3] = {0.}, arrb[3] = {0.};

        TEnvironment Environment(0.,0.,0.) ;
        arrbWave[0] = mpLoad->calcAirResistMom(Environment, 0., arrMeanPhVect[0])
                * getInvSumJ0();
       arrbWave[1] =  mElectMotor.mZp * arrMeanPhVect[0] * arrMeanPhVect[2 ];
       arrbWave[2] = -mElectMotor.mZp * arrMeanPhVect[0] * arrMeanPhVect[1 ];

        double arrCC[ 2 * QVARS] = {0.},arrU[2] = {0.}
                    , arrT[2] = {0.}, arrT1[QVARS] = {0.}, arrB[QVARS *2] = {0.};
        QStatSolutionParams StatSolutionParams;
        calcStationaryParams(arrObjective,  mCmpArrLamb, &StatSolutionParams, VAlMomOut);
        memcpy(arrCC, StatSolutionParams.marrGears, 2 * QVARS * sizeof(double));
        MtrxMultMatrx(arrCC,2, QVARS, marrIntegrators, 1, arrT) ;
        MtrxSumMatrx(arrT, StatSolutionParams.marrStatU, 1, 2, arrU) ;
        calcMtrxB(arrB);
        MtrxMultMatrx(arrB,QVARS,2,  arrU, 1, arrT1) ;

        double valL = mElectMotor.getL();
//        /double valJ = 1. / getInvSumJ0();
       // double temp_LId = (1000000. * valL)/ (1000000. * arrMeanPhVect[1]);
       // double valZp = mElectMotor.mZp;
       // for (int i =0; i < 3; ++i)
       // {
        //  arrb [i] = arrbWave[i] +  arrT1[i] ;
        //}
        switch(QUantParams)
        {
        case 1:
            arrParams[1] = getParam(1);
            arrParams[2] = getParam(2);
            arrParams[0] = - 3./ 2. * mElectMotor.mZp *arrMeanPhVect[2] * getParam(2)
              - mpLoad->calcAirResistMom(Environment, 0., arrMeanPhVect[0]) ;

            break;
        case 2 :
            arrParams[2] = getParam(2);
            arrParams[0] =- 3./ 2. * mElectMotor.mZp *arrMeanPhVect[2] * getParam(2)
                    - mpLoad->calcAirResistMom(Environment, 0., arrMeanPhVect[0]) ;
            arrbWave[1] =  mElectMotor.mZp * arrMeanPhVect[0] * arrMeanPhVect[2 ];
            arrbWave[2] = -mElectMotor.mZp * arrMeanPhVect[0] * arrMeanPhVect[1 ]
                    - mElectMotor.mZp * arrMeanPhVect[0] * mElectMotor.mPsi_f/ valL;
            arrb [1] = arrbWave[1] +  arrT1[1] ;
            arrb [2] = arrbWave[2] +  arrT1[2] ;
            arrParams[1] =  valL *(arrMeanPhVect[1] *arrb[1] + arrMeanPhVect[2] *arrb[2] )
                    /(arrMeanPhVect[1] * arrMeanPhVect[1] + arrMeanPhVect[2] * arrMeanPhVect[2]);
            break;


        default:
            break;
        }

}

//---------------------------------------------

void QCtrlDnmkCorrect_M_R::createMeanPhVect(const double *arrObjective,const double VAlMomOut,double *arrMeanPhVect)
{
    QStatSolutionParams StatSolutionParams;
    calcStationaryParams(arrObjective,  mCmpArrLamb, &StatSolutionParams, VAlMomOut);
    for(int i = 0; i < QVARS ; ++i)
    {
      arrMeanPhVect[i] =  StatSolutionParams.marrStatPhVect[i] +  marrIntegrators[i];
    }
    arrMeanPhVect[3] += arrObjective[1] *(mTCur - mT0);// ??
}


