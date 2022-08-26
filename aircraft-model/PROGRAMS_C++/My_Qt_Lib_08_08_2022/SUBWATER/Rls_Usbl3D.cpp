#include "Rls_Usbl3D.h"
#include "TrueMeasParams.h"
#include "Gauss.h"

QRls_Usbl3D::QRls_Usbl3D():QHidroRLS()
{
  mSig_q = 0.01;
  mSig_e = 0.01;
}
//-------------------------------------------------------
// конструктор копирования
QRls_Usbl3D :: QRls_Usbl3D (const  QRls_Usbl3D &R):QHidroRLS( R)
 {
  mSig_q = R.mSig_q;
  mSig_e = R.mSig_e;
 }
//-------------------------------------------------------
 // оператор присваивания
QRls_Usbl3D  &QRls_Usbl3D::operator=( const QRls_Usbl3D  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      QHidroRLS:: operator= (R);
      mSig_q = R.mSig_q;
      mSig_e = R.mSig_e;
     return *this ;
 }
//--------------------------------------------------
  // парам конструктор 1
QRls_Usbl3D:: QRls_Usbl3D (const double SigR, const double DiagWidth
                           ,const double Sig_q, const double Sig_e)
  :QHidroRLS (SigR,  DiagWidth)
 {
  mSig_q = Sig_q;
  mSig_e = Sig_e;
 }
//------------------------------------------
// вычисление  измерения КУ
void QRls_Usbl3D::calc_qZv_and_eZv(const QTrueMeasParams trueMeasParams, double &qZv, double &eZv
                                   , double &Sig_q, double &Sig_e)
{
    // ТЕСТ
    qZv = trueMeasParams.mq +getGauss(0, mSig_q);
    eZv = trueMeasParams.me + getGauss(0, mSig_e);
    Sig_q = mSig_q;
    Sig_e = mSig_e;
}
