#include "TrueMeasParams.h"
#include "string.h"

QTrueMeasParams::QTrueMeasParams()
{
    mTzapr = 0.;
    mTotv = 0.;
    memset(marrSVessWave, 0, 3 * sizeof(double));
    memset(marrMuWave, 0, 3 * sizeof(double));
    memset(marrMuWaveZv, 0, 3 * sizeof(double));
    memset(marrSVess, 0, 3 * sizeof(double));
    memset(marrMu, 0, 3 * sizeof(double));
    memset(marrMuZv, 0, 3 * sizeof(double));
    mq = 0.;
    me = 0.;
}

// конструктор копирования
QTrueMeasParams :: QTrueMeasParams (const  QTrueMeasParams &R)
 {
    memcpy(marrSVessWave, R.marrSVessWave, 3 * sizeof(double));
    memcpy(marrMuWave, R.marrMuWave, 3 * sizeof(double));
    memcpy(marrMuWaveZv, R.marrMuWaveZv, 3 * sizeof(double));
    memcpy(marrSVess, R.marrSVess, 3 * sizeof(double));
    memcpy(marrMu, R.marrMu, 3 * sizeof(double));
    memcpy(marrMuZv, R.marrMuZv, 3 * sizeof(double));

    mTzapr = R.mTzapr;
    mTotv = R.mTotv;
    mq = R.mq;
    me = R.me;

 }

 // оператор присваивания
QTrueMeasParams  &QTrueMeasParams::operator=( const QTrueMeasParams  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      memcpy(marrSVessWave, R.marrSVessWave, 3 * sizeof(double));
      memcpy(marrMuWave, R.marrMuWave, 3 * sizeof(double));
      memcpy(marrMuWaveZv, R.marrMuWaveZv, 3 * sizeof(double));
      memcpy(marrSVess, R.marrSVess, 3 * sizeof(double));
      memcpy(marrMu, R.marrMu, 3 * sizeof(double));
      memcpy(marrMuZv, R.marrMuZv, 3 * sizeof(double));

      mTzapr = R.mTzapr;
      mTotv = R.mTotv;
      mq = R.mq;
      me = R.me;

     return *this ;
 }

  // парам конструктор 1
QTrueMeasParams:: QTrueMeasParams (const double  *arrSVessWave
                                   ,const double  *arrMuWave,const double  *arrMuWaveZv,
                                     const double  *arrSVess
                                   ,const double  *arrMu,const double  *arrMuZv,
                                     const double  Tzapr,const double  Totv
                                     ,const double q, const double e)

 {
    memcpy(marrSVessWave, arrSVessWave, 3 * sizeof(double));
    memcpy(marrMuWave, arrMuWave, 3 * sizeof(double));
    memcpy(marrMuWaveZv, arrMuWaveZv, 3 * sizeof(double));
    memcpy(marrSVess, arrSVess, 3 * sizeof(double));
    memcpy(marrMu, arrMu, 3 * sizeof(double));
    memcpy(marrMuZv, arrMuZv, 3 * sizeof(double));
    mTzapr = Tzapr;
    mTotv = Totv;
    mq = q;
    me = e;
 }
