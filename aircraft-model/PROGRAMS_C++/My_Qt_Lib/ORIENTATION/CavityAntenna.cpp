#include "CavityAntenna.h"
#include <math.h>

QCavityAntenna::QCavityAntenna()
{
    // сечение канала линейки (большая стенка), [м]
     mLin_a  = 0.;
    // шаг между щелями, [м]
     mLin_d  = 0.;
    // сечение канала делителя мощности, [м]
     mDM_a  = 0.;
    // длина витка делителя мощности, [м]
     mDM_L  = 0.;
    // шаг выходов в делители мощности, [м]
     mDM_d  = 0.;
}

// конструктор копирования
 QCavityAntenna ::QCavityAntenna (const QCavityAntenna &R)
 {
    mLin_a = R.mLin_a;
    mLin_d = R.mLin_d;
    mDM_a = R.mDM_a;
    mDM_L = R.mDM_L;
    mDM_d = R.mDM_d;

 }


 // оператор присваивания
 QCavityAntenna &QCavityAntenna::operator=(const QCavityAntenna  &R)
 {
     mLin_a = R.mLin_a;
     mLin_d = R.mLin_d;
     mDM_a = R.mDM_a;
     mDM_L = R.mDM_L;
     mDM_d = R.mDM_d;

    return *this ;
 }

  // парам конструктор1
 QCavityAntenna::QCavityAntenna (const double Lin_a, const double Lin_d
                      ,  const double DM_a,  const double DM_L,  const double DM_d)
 {
    mLin_a = Lin_a;
    mLin_d = Lin_d;
    mDM_a = DM_a;
    mDM_L = DM_L;
    mDM_d = DM_d;
 }
 //------------------------------
 void QCavityAntenna::createInpDataArray (double *arrInpData)
{
    arrInpData[0] = mLin_a;
    arrInpData[1] = mLin_d;
    arrInpData[2] = mDM_a;
    arrInpData[3] = mDM_L;
    arrInpData[4] = mDM_d;
}
 //------------------------------
 double QCavityAntenna::calcLamb_w(const double VAlLamb)
 {
     return VAlLamb/ sqrt(1. - VAlLamb * VAlLamb/ mLin_a/ mLin_a/4.);
 }

 //------------------------------
 double QCavityAntenna::calcU(const double VAlLamb)
 {
     double valLamb_w = calcLamb_w( VAlLamb);
     double val_dfi = 2. * M_PI *mLin_d /valLamb_w - M_PI;
     double valKsi = val_dfi /(2. * M_PI *mLin_d /VAlLamb);
     //double temp = 1. - VAlLamb * VAlLamb/ mLin_a/ mLin_a/4.;
    // if (temp > 0.999999999)
     //{
     //  temp = 0.999999999;
     //}


     return asin(valKsi);
 }


 //------------------------------
 double QCavityAntenna::calcV(const double VAlLamb)
 {
     double lam_d = calcLambDM(VAlLamb);
     int k = mDM_L / lam_d;
     double del_fi = 2. * M_PI *(mDM_L / lam_d - ((double)k));// -M_PI;
     if (del_fi > M_PI)
     {
       del_fi -= 2. * M_PI;
     }
     double t3 = fmod(2. * M_PI *mDM_L/ lam_d, 2. * M_PI);// -M_PI;

     double valKsi = del_fi * VAlLamb/(2. * M_PI  *mDM_d);

     if (fabs(valKsi)> 0.999999999)
     {
       valKsi = (valKsi>0.)?0.99999999:-0.99999999;
     }

    // return acos(valKsi) - M_PI/2.;
     return -asin(valKsi);
 }
 //---------------------------
 double QCavityAntenna::calcLambDM(const double VAlLamb)
 {
     double t = 1. -VAlLamb* VAlLamb/ 4./mDM_a/ mDM_a;
     return VAlLamb/ sqrt(t);
 }

 //------------------------------------------------
 void QCavityAntenna::calc_DM_params_array(const double VAlLamb, double *arrParamsDM)
 {
     double lam_d = calcLambDM(VAlLamb);
     int k = mDM_L / lam_d;
     double del_fi = 2. * M_PI *(mDM_L / lam_d - ((double)k));// -M_PI;
     if (del_fi > M_PI)
     {
       del_fi -= 2. * M_PI;
       k++;
     }


     arrParamsDM[1] = 1./2./ mDM_a;
     arrParamsDM[0] = mDM_L/ mDM_d;
     arrParamsDM[2] = ((double)(  k))/ mDM_d;
 }
 //-------------------------------------
 int QCavityAntenna::calcDM_k(const double VAlLamb)
 {
     double lam_d = calcLambDM(VAlLamb);
     int k = mDM_L / lam_d;
     double del_fi = 2. * M_PI *(mDM_L / lam_d - ((double)k));// -M_PI;
     if (del_fi > M_PI)
     {

       k++;
     }
     return k;
 }

