#include "Solver3D_Angs.h"
#include "Table_1D.h"
#include "BigMeasure.h"
//#include "DataExchange.h"
#include "MatrixProccess.h"
#include <math.h>

QSolver3D_Angs::QSolver3D_Angs():QSolver2D_Angs()
{
    mDimX = 3;
    mDimY = 2;
}
// конструктор копирования
QSolver3D_Angs :: QSolver3D_Angs (const  QSolver3D_Angs &R):QSolver2D_Angs( R)
 {

 }

 // оператор присваивания
QSolver3D_Angs  &QSolver3D_Angs::operator=( const QSolver3D_Angs  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      QSolver2D_Angs:: operator= (R);

     return *this ;
 }

// парам конструктор 1
QSolver3D_Angs::QSolver3D_Angs (const QBigMeasure *parrBigMeasures
                              , const int  QntMeas    ,const TTable_1D tblEstPrfl
                              , const double Deepth, const double Toler
                            ,const double *arrPAnt_PSK,const double *arrSBeacon_GSK)
:QSolver2D_Angs ( parrBigMeasures,  QntMeas    ,tblEstPrfl
            ,Deepth,Toler,arrPAnt_PSK,arrSBeacon_GSK)
{
  mDimX = 3;
  mDimY = 2;
}


//------------------------------------------
void  QSolver3D_Angs::collectGrlsMeasure(const QBigMeasure &Meas
                                             ,double *arrMeasure,double *arrD)
{

   arrMeasure[0] = Meas.mqzv;
   arrMeasure[1] = Meas.mezv;
   memset(arrD, 0, 4 * sizeof(double));
   arrD[0] = Meas.mSig_q * Meas.mSig_q;
   arrD[3] = Meas.mSig_e * Meas.mSig_e;

}
//------------------------------
bool QSolver3D_Angs::calc_arrNeviaz_and_dArrNeviaz_po_dX_3D( const int NUmBigMeasure
            ,const double e_est,double *arr_dV_po_dJ ,double *arrt7, double *arrNeviaz , double *dArrNeviaz_po_dX)
{
    QBigMeasure Meas = mVectBigMeasures.at(NUmBigMeasure);
    MtrxMultMatrx(&arr_dV_po_dJ[3],1, 3, arrt7,mDimX, &dArrNeviaz_po_dX[mDimX]) ;
    arrNeviaz [1] = e_est - Meas.mezv;
    if(fabs(arrNeviaz [1])> 1.)
    {
        int iii =0;
    }
}

//------------------------------
bool QSolver3D_Angs::calc_arrNeviaz_3D( const int NUmBigMeasure
            ,const double e_est, double *arrNeviaz )
{
    QBigMeasure Meas = mVectBigMeasures.at(NUmBigMeasure);
    arrNeviaz [1] = e_est - Meas.mezv;
    if(fabs(arrNeviaz [1])> 1.)
    {
        int iii =0;
    }
    return true;
}
