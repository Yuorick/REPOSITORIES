#include "Rls_Usbl2D.h"
#include "TrueMeasParams.h"
#include "Gauss.h"
#include "HidroRLS.h"

QRls_Usbl2D::QRls_Usbl2D():QHidroRLS()
{
  mSig_q = 0.01;

}

// конструктор копирования
QRls_Usbl2D :: QRls_Usbl2D (const  QRls_Usbl2D &R):QHidroRLS( R)
 {
  mSig_q = R.mSig_q;


 }

 // оператор присваивания
QRls_Usbl2D  &QRls_Usbl2D::operator=( const QRls_Usbl2D  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      QHidroRLS:: operator= (R);
      mSig_q = R.mSig_q;


     return *this ;
 }

  // парам конструктор 1
QRls_Usbl2D:: QRls_Usbl2D (const double SigR, const double DiagWidth
                           ,const double Sig_q)
  :QHidroRLS (SigR,  DiagWidth)
 {
  mSig_q = Sig_q;

 }
//------------------------------------------

// вычисление  измерения КУ
void QRls_Usbl2D::calc_qZv_and_eZv(const QTrueMeasParams trueMeasParams, double &qZv, double &eZv
                                   , double &Sig_q, double &Sig_e)
{
    qZv =  trueMeasParams.mq + getGauss(0, mSig_q);
    eZv = NODATA;
    Sig_q = mSig_q;
    Sig_e = NODATA;
}
