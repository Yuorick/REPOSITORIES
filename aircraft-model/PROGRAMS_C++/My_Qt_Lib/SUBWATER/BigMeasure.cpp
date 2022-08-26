#include "BigMeasure.h"
#include "string.h"
#include <math.h>
#include "MatrixProccess.h"
#include "PeaceVess.h"

//QBigMeasure::~QBigMeasure()
//{

//}

QBigMeasure::QBigMeasure()
{
    mTzaprZv = 0.;
    mTotvZv = 0.;
    mTobr = 0.;
    memset(marrSVessWaveZv, 0, 3 * sizeof(double));
    memset(marrMuWaveZv, 0, 3 * sizeof(double));
    memset(marrSVessZv, 0, 3 * sizeof(double));
    memset(marrMuZv, 0, 3 * sizeof(double));
    mDimMeas = 0;
    mqzv = 0.;
    mezv = 0.;
    // СКЗ измерения времени
     mSig_t= 0.;
    // СКЗ измерения КУ
     mSig_q= 0.;
    // СКЗ измерения УМ
     mSig_e= 0.;

}

// конструктор копирования
QBigMeasure :: QBigMeasure (const  QBigMeasure &R)
 {
    memcpy(marrSVessWaveZv, R.marrSVessWaveZv, 3 * sizeof(double));
    memcpy(marrMuWaveZv, R.marrMuWaveZv, 3 * sizeof(double));
    memcpy(marrSVessZv, R.marrSVessZv, 3 * sizeof(double));
    memcpy(marrMuZv, R.marrMuZv, 3 * sizeof(double));
    mTzaprZv = R.mTzaprZv;
    mTotvZv = R.mTotvZv;
    mTobr = R.mTobr;
    mDimMeas = R.mDimMeas;
    mqzv = R.mqzv;
    mezv = R.mezv;
    mSig_t  = R.mSig_t;
    mSig_q = R.mSig_q;
    mSig_e = R.mSig_e;

 }

 // оператор присваивания
QBigMeasure  &QBigMeasure::operator=( const QBigMeasure  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      memcpy(marrSVessWaveZv, R.marrSVessWaveZv, 3 * sizeof(double));
      memcpy(marrMuWaveZv, R.marrMuWaveZv, 3 * sizeof(double));
      memcpy(marrSVessZv, R.marrSVessZv, 3 * sizeof(double));
      memcpy(marrMuZv, R.marrMuZv, 3 * sizeof(double));
      mTzaprZv = R.mTzaprZv;
      mTotvZv = R.mTotvZv;
      mTobr = R.mTobr;
      mDimMeas = R.mDimMeas;
      mqzv = R.mqzv;
      mezv = R.mezv;
      mSig_t  = R.mSig_t;
      mSig_q = R.mSig_q;
      mSig_e = R.mSig_e;

     return *this ;
 }

  // парам конструктор 1
QBigMeasure:: QBigMeasure (const double  *arrSVessWaveZv,const double  *arrMuWaveZv,
                     const double  *arrSVessZv,const double  *arrMuZv,
                     const double  TzaprZv,const double  TotvZv,
                     const double  Tobr,const int  DimMeas
                     ,const double qzv, const double ezv
                     ,const double Sig_t,const double Sig_q,const double Sig_e)

 {
    memcpy(marrSVessWaveZv, arrSVessWaveZv, 3 * sizeof(double));
    memcpy(marrMuWaveZv, arrMuWaveZv, 3 * sizeof(double));
    memcpy(marrSVessZv, arrSVessZv, 3 * sizeof(double));
    memcpy(marrMuZv, arrMuZv, 3 * sizeof(double));
    mTzaprZv = TzaprZv;
    mTotvZv = TotvZv;
    mTobr = Tobr;
    mDimMeas = DimMeas;
    mqzv = qzv;
    mezv = ezv;
    mSig_t = Sig_t;
    mSig_q = Sig_q;
    mSig_e = Sig_e;
 }
//---------------------------------
bool QBigMeasure::isEqual(QBigMeasure &meas0, QBigMeasure &meas1, const double VAlTolerDist
                          , const double VAlTolerSinsAngs,  const double VAlTolerT, const double VAlTolerMeasAngs )
{
    if(fabs(meas0.mTzaprZv -meas1.mTzaprZv ) > VAlTolerT)
    {
        return false;
    }

    if(fabs(meas0.mTotvZv -meas1.mTotvZv ) > VAlTolerT)
    {
        return false;
    }

    if(fabs(meas0.mTobr -meas1.mTobr ) > VAlTolerMeasAngs)
    {
        return false;
    }

    if(fabs(meas0.mqzv -meas1.mqzv ) > VAlTolerMeasAngs)
    {
        return false;
    }

    if(fabs(meas0.mezv -meas1.mezv ) > VAlTolerMeasAngs)
    {
        return false;
    }


    if(fabs(meas0.mezv -meas1.mezv ) > VAlTolerMeasAngs)
    {
        return false;
    }


    if(fabs(Norm3_A_Minus_B( meas0.marrSVessWaveZv, meas1.marrSVessWaveZv) ) > VAlTolerDist)
    {
        return false;
    }

    if(fabs(Norm3_A_Minus_B( meas0.marrMuWaveZv, meas1.marrMuWaveZv) ) > VAlTolerSinsAngs)
    {
        return false;
    }


    if(fabs(Norm3_A_Minus_B( meas0.marrSVessZv, meas1.marrSVessZv) ) > VAlTolerDist)
    {
        return false;
    }

    if(fabs(Norm3_A_Minus_B( meas0.marrMuZv, meas1.marrMuZv) ) > VAlTolerSinsAngs)
    {
        return false;
    }


    if(fabs(meas0.mSig_t -meas1.mSig_t ) > VAlTolerT* 0.1)
    {
        return false;
    }
    if(fabs(meas0.mSig_q -meas1.mSig_q ) > VAlTolerMeasAngs * 0.1)
    {
        return false;
    }
    if(fabs(meas0.mSig_e -meas1.mSig_e ) > VAlTolerMeasAngs * 0.1)
    {
        return false;
    }
return true;

}
//------------------------
void QBigMeasure::fill_info_row( double *arrP_Gps, double *arrRow)
{
  arrRow[0] =mTzaprZv;

  double arrGps_KGSK[3] = {0.};
  QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( arrP_Gps, marrMuWaveZv, NULL,arrGps_KGSK,3 );
  double arrSVessWaveZv[3] = {0.};
  MtrxSumMatrx(marrSVessWaveZv, arrGps_KGSK,1, 3, arrSVessWaveZv);
  memcpy(&arrRow[1], arrSVessWaveZv, 3 * sizeof(double));

  arrRow[4] =  marrMuWaveZv[0] * 180./M_PI;
  arrRow[5] =  marrMuWaveZv[1] * 180./M_PI;
  arrRow[6] =  marrMuWaveZv[2] * 180./M_PI;


  arrRow[7] =  mTotvZv;

  double arrSVessZv[3] = {0.};
  QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( arrP_Gps, marrMuZv, NULL,arrGps_KGSK,3 );
  MtrxSumMatrx(marrSVessZv, arrGps_KGSK,1, 3, arrSVessZv);
  memcpy(&arrRow[8], arrSVessZv, 3 * sizeof(double));

  arrRow[11] =  marrMuZv[0] * 180./M_PI;
  arrRow[12] =  marrMuZv[1] * 180./M_PI;
  arrRow[13] =  marrMuZv[2] * 180./M_PI;


  arrRow[14] = mqzv* 180./M_PI;
  arrRow[15] = -mezv* 180./M_PI;


}

