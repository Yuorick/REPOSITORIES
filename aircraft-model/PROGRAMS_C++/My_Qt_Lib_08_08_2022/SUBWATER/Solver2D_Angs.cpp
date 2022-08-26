#include "Solver2D_Angs.h"
#include "BigMeasure.h"
#include <QDebug>
#include <math.h>
#include "MatrixProccess.h"
#include "CoordSystTrsf.h"
#include "SubWaterBeam.h"
#include "PeaceVess.h"
#include "LblSolver.h"

QSolver2D_Angs::QSolver2D_Angs():QPosSolver()
{
    mDimX = 3;
    mDimY =1;
    memset(marrPAnt_PSK, 0, 3 * sizeof(double));
    memset(marrSBeaconXY_GSK, 0, 2 * sizeof(double));

}
// конструктор копирования
  QSolver2D_Angs :: QSolver2D_Angs (const  QSolver2D_Angs &R):QPosSolver( R)
 {
      memcpy(marrPAnt_PSK, R.marrPAnt_PSK, 3 * sizeof(double));
      memcpy(marrSBeaconXY_GSK, R.marrSBeaconXY_GSK, 2 * sizeof(double));
      mvct_e_psk = R.mvct_e_psk;
 }

 // оператор присваивания
  QSolver2D_Angs  &QSolver2D_Angs::operator=( const QSolver2D_Angs  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      QPosSolver:: operator= (R);
      memcpy(marrPAnt_PSK, R.marrPAnt_PSK, 3 * sizeof(double));
      memcpy(marrSBeaconXY_GSK, R.marrSBeaconXY_GSK, 2 * sizeof(double));
      mvct_e_psk = R.mvct_e_psk;


     return *this ;
 }


  // парам конструктор 1
 QSolver2D_Angs:: QSolver2D_Angs (const QBigMeasure *parrBigMeasures
                                  , const int  QntMeas    ,const TTable_1D tblEstPrfl
                                  , const double Deepth
                                  , const double Toler,const double *arrPAnt_PSK
                                  ,const double *arrSBeaconXY_GSK)
  :QPosSolver (parrBigMeasures,  QntMeas    ,tblEstPrfl
                , Deepth, 3,  1,Toler)
 {
     mDimX = 3;
     mDimY = 1;
     memcpy(marrPAnt_PSK, arrPAnt_PSK, 3 * sizeof(double));
     memcpy(marrSBeaconXY_GSK, arrSBeaconXY_GSK, 2 * sizeof(double));
     fill_vct_e_PSK();
 }


//--------------------------------------------
void QSolver2D_Angs::fill_vct_e_PSK()
{
    double arr_e_psk[3] ={0.};
    mvct_e_psk.resize(3 * mVectBigMeasures.size());
    for (int i = 0; i < mVectBigMeasures.size(); ++i)
    {
       QBigMeasure Meas = mVectBigMeasures.at(i);
       calc_arr_e_PSK(Meas, arr_e_psk);
       mvct_e_psk.replace(3 * i    ,arr_e_psk[0] );
       mvct_e_psk.replace(3 * i + 1,arr_e_psk[1] );
       mvct_e_psk.replace(3 * i + 2,arr_e_psk[2] );
    }

}
//------------------------------------------
void  QSolver2D_Angs::collectGrlsMeasure(const QBigMeasure &Meas
                                             ,double *arrMeasure,double *arrD)
{

   arrMeasure[0] = Meas.mqzv;
   arrD[0] = Meas.mSig_t * Meas.mSig_t ;

}
//------------------------------

bool QSolver2D_Angs::calc_arrNeviaz_and_dArrNeviaz_po_dX(double *arrX, const int NUmBigMeasure
                                       , double *arrNeviaz , double *dArrNeviaz_po_dX)
{
    QBigMeasure Meas = mVectBigMeasures.at(NUmBigMeasure);
   // 1. вычисление массива Y
    // первые 3 координаты - вектор положения антенны
    // затем 3 координаты - вектор положения маяка
    // затем 3 угла ориентации антенны
    double arrY [9] = {0.};

    // персчет вектора положения антенны в ГСК

    // создание матрицы перехода из   КГСК в ПСК
     double arrEilers[3] = {0.},  arr_KGSK[3] = {0.};
    MatrxMultScalar(Meas.marrMuZv, 1, 3, -1.,arrEilers);
    double matrPereh_PSK_V_KGSK[9] = {0} ;
    QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilers, matrPereh_PSK_V_KGSK) ;
    // вектор положения в ПСК-центр тяжести
    // вычисление вектора положениея в ПСК

    MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, marrPAnt_PSK, 1, arr_KGSK) ;

    MtrxSumMatrx(arr_KGSK, Meas.marrSVessZv,1, 3, arrY) ;
    // !



    arrY [3] = marrSBeaconXY_GSK[0]; // X маяка
    arrY [4] = marrSBeaconXY_GSK[1]; // Y маяка
    arrY [5] = mDeepth; // Z маяка
    arrY[6] = arrX[0];// угол
    arrY[7] = arrX[1];// угол
    arrY[8] = arrX[2];// угол

    // 1!

    // 2. вычисление массива Z
    double arrZ [4] = {0.};

   calc_arrZ(arrY, arrZ);
  // const double VAlTetta = arrZ[6];


    // 2!

    // 3. Вычисление массива I
    double arrI[6] = {0.};
    calc_arrI(arrZ, arrI);

    // 4. вектор направления в ГСК
    double *parr_e_gsk =arrI;
    // !4

    // 5. вектор направления в ПСК
    double arr_e_psk[3]={0.};
    MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, parr_e_gsk, 1, arr_e_psk) ;
    // !5



    // 6. вычисление вектора единичного направления в АСПК  J
    double arrJ[3] = {0.};
    calcVectJ_(arrX,  arr_e_psk,arrJ);

    // 6 !

    //  ОТЛАДКА
   // calcVectJ(Meas.marrMuZv, arrI, arrJ1);
    // ОТЛАДКА !

  // вычисление вектора V
  double arrV[2] ={0.};
  //arrV[0 ] = asin(arrJ[0]/ sqrt(arrJ[0] * arrJ[0] + arrJ[1] * arrJ[1]));
  arrV[0 ] = QPeaceVess::calcCourseAngle(arrJ[0], arrJ[1]);
  arrV[1 ] = atan2(arrJ[2],sqrt(arrJ[0] * arrJ[0] + arrJ[1] * arrJ[1]));

// ОТЛАДКА
 // double arrGSK_XYZ[3] = {0.};
 // arrGSK_XYZ[0] = arrX[3];
 // arrGSK_XYZ[1] = arrX[4];
 // arrGSK_XYZ[2] = mDeepth;
 // double arrAntPosParams[6] = {0.};
 // arrAntPosParams[0] = arrX[0];
//  arrAntPosParams[1] = arrX[1];
 // arrAntPosParams[2] = arrX[2];
 // arrAntPosParams[3] = arrX[5];
 // arrAntPosParams[4] = arrX[6];
 // arrAntPosParams[5] = arrX[7];
//  double pval_q = 0., pval_e = 0.,  pval_t= 0.;
//  transf_GSK_XYZ_to_USBL3D(mtblEstPrfl, arrGSK_XYZ, Meas.marrSVessZv, Meas.marrMuZv
 //   ,arrAntPosParams,  &pval_q, &pval_e,  &pval_t);
//  int uiu = 0;
  // ОТЛАДКА !
  arrNeviaz[0] = calcCorseNeviaz(arrV[0 ], Meas.mqzv);
  if(fabs(arrNeviaz[0])> 1.)
  {
      int iii = 0;
  }


  // 5. вычисление матрицы dV_po_dJ[2*3]
  double arr_dV_po_dJ[6] = {0.};
  double val_r = sqrt(arrJ[0] * arrJ[0] + arrJ[1] * arrJ[1]);
  double val_R2 = (arrJ[0] * arrJ[0] + arrJ[1] * arrJ[1] + arrJ[2] * arrJ[2]);
    arr_dV_po_dJ[0] =  arrJ[1]/ (val_r * val_r);
    arr_dV_po_dJ[1] = -arrJ[0]/ (val_r * val_r);
    arr_dV_po_dJ[3] = -arrJ[0] * arrJ[2]/ (val_R2  * val_r);
    arr_dV_po_dJ[4] = -arrJ[1] * arrJ[2]/ (val_R2  * val_r);
    arr_dV_po_dJ[5] = val_r/ val_R2;
    // 5!

    // 6. Вычисление матрицы dJ_po_dX[3 *3]
    double arr_dJ_po_dX[3 *3] = {0.};

    calcMtrx_dJ_po_dX(arrX, arr_e_psk, arr_dJ_po_dX);



    // ОТЛАДКА
//double arr_dJ_po_dX1_T[6 *3] = {0.}, arrJ1[6] = {0.},arrX1[3] = {0.}
 //          ,arr_dJ_po_dX1[6 *3] = {0.},arrdel[3] = {0.}
   //     , arrd[6]= {0.001, 0.001, 0.001, 0.0001, 0.0001, 0.0001};
 // for (int i =0; i < 3;++i)
  // {
   //   memcpy(arrX1,arrX, 3 * sizeof(double));
    //   arrX1[i] += arrd[i];
   // calcVectJ_(arrX1,  arr_e_psk,arrJ1);

   //  MtrxMinusMatrx(arrJ1, arrJ,1, 3, arrdel);
    //   MatrxMultScalar(arrdel, 1, 3, 1./arrd[i],&arr_dJ_po_dX1_T[i * 3]);
   // }
  // MatrTransp( arr_dJ_po_dX1_T, 3, 3, arr_dJ_po_dX1);
  // double ii = 0;// !
    // ОТЛАДКА !

   MtrxMultMatrx(arr_dV_po_dJ,1, 3, arr_dJ_po_dX,mDimX, dArrNeviaz_po_dX) ;

   calc_arrNeviaz_and_dArrNeviaz_po_dX_3D(NUmBigMeasure, arrV[1]
                         ,arr_dV_po_dJ ,arr_dJ_po_dX, arrNeviaz , dArrNeviaz_po_dX);
    return true;
}
//--------------------------------------
void QSolver2D_Angs::calc_arr_e_PSK(QBigMeasure Meas, double *arr_e_psk)
{
    // 1. вычисление массива Y
     // первые 3 координаты - вектор положения антенны
     // затем 3 координаты - вектор положения маяка
     // затем 3 угла ориентации антенны
     double arrY [6] = {0.};

     // персчет вектора положения антенны в ГСК

     // создание матрицы перехода из   КГСК в ПСК
      double arrEilers[3] = {0.},  arr_KGSK[3] = {0.};
     MatrxMultScalar(Meas.marrMuZv, 1, 3, -1.,arrEilers);
     double matrPereh_PSK_V_KGSK[9] = {0} ;
     QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilers, matrPereh_PSK_V_KGSK) ;
     // вектор положения в ПСК-центр тяжести
     // вычисление вектора положениея в ПСК

     MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, marrPAnt_PSK, 1, arr_KGSK) ;

     MtrxSumMatrx(arr_KGSK, Meas.marrSVessZv,1, 3, arrY) ;
     // !



     arrY [3] = marrSBeaconXY_GSK[0]; // X маяка
     arrY [4] = marrSBeaconXY_GSK[1]; // Y маяка
     arrY [5] = mDeepth; // Z маяка


     // 1!

     // 2. вычисление массива Z
     double arrZ [4] = {0.};

    calc_arrZ(arrY, arrZ);
    const double VAlTetta = arrZ[3];


     // 2!

     // 3. Вычисление массива I
     double arrI[6] = {0.};
     calc_arrI(arrZ, arrI);

     // 4. вектор направления в ГСК
     double *parr_e_gsk =arrI;
     // !4

     // 5. вектор направления в ПСК

     MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, parr_e_gsk, 1, arr_e_psk) ;
     // !5


}


//------------------------------
// промежуточный вектор
//INPUT:
// arrY[9] -3 координаты антенны в ГСК, 3 координаты маяка в ГСК, 3 угла антенны
// OUTPUT:
//arrZ[4]
//arrZ[0] - горизонтальное расстояние от антенны до маяка

//arrZ[1] - расстояние от маяка до антенны по оси OX ГСК
//arrZ[2] - расстояние от маяка до антенны по оси OY ГСК
//arrZ[3] - угол скольжения
void QSolver2D_Angs::calc_arrZ(double *arrY, double *arrZ)
{
    arrZ[0]  = sqrt((arrY [3] - arrY [0]) * (arrY [3] - arrY [0]) + (arrY [4] - arrY [1]) *(arrY [4] - arrY [1]));
    arrZ[1]  = arrY [3] - arrY [0];
    arrZ[2]  = arrY [4] - arrY [1];
    calcTetta( -arrY[2],-arrY [5],arrZ[0], mtblEstPrfl, arrZ[3]);
}
//---------------------------------------

bool QSolver2D_Angs::calc_arrNeviaz_and_dArrNeviaz_po_dX_3D( const int NUmBigMeasure
,const double e_est,double *arr_dV_po_dJ ,double *arr_dJ_po_dX, double *arrNeviaz , double *dArrNeviaz_po_dX)
{

}
//----------------------------
void QSolver2D_Angs::calcVectJ_(double *arrX,  double *arr_e_psk,double *arrJ)
{

    double matrPereh_ASK_V_PSK[9]= {0.};
    QCoordSystTrsf::calcMtrx3_ASPK_v_PSK(arrX,matrPereh_ASK_V_PSK) ;

    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK,3, 3, arr_e_psk,1, arrJ) ;
}


//---------------------------------
void QSolver2D_Angs::calc_arrI(double *arrZ, double *arrI)
{
    memset(arrI, 0, 3 * sizeof(double));
    arrI[0] = arrZ[1]* cos(arrZ[3])/arrZ[0];
    arrI[1] = arrZ[2]* cos(arrZ[3])/arrZ[0];
    arrI[2] = -sin(arrZ[3]);

}
//-----------------------------
double QSolver2D_Angs::calcCorseNeviaz(const double qest, const double qzv)
{
    double nev = qest - qzv;
    double nev1 = nev + 2. * M_PI;
    if (fabs(nev) <= fabs(nev1))
    {
        return nev;
    }
    return nev1;

}
//------------------------

void QSolver2D_Angs::calcMtrx_dJ_po_dX(double *arrX, double *arr_e_psk
                                        , double *arr_dJ_po_dX)
{
    double arr_dM_po_dX0[9] = {0.}, arrt[3] = {0.};
   QCoordSystTrsf::calc_dMtrx3_ASPK_v_PSK_po_dQ(arrX,arr_dM_po_dX0);
   MtrxTranspMultMatrx(arr_dM_po_dX0,3, 3, arr_e_psk,1, arrt) ;
   for (int i =0;i< 3; ++i)
   {
      arr_dJ_po_dX[3 *i ] = arrt[i];
   }

   double arr_dM_po_dX1[9] = {0.};
   //QCoordSystTrsf::calc_dMLeft_psk_to_kgsk_po_dPsi(arrEilers1,arr_dM_po_dI5);
   QCoordSystTrsf::calc_dMtrx3_ASPK_v_PSK_po_dPsi(arrX,arr_dM_po_dX1);
   MtrxTranspMultMatrx(arr_dM_po_dX1,3, 3, arr_e_psk,1, arrt) ;
   for (int i =0;i< 3; ++i)
   {
      arr_dJ_po_dX[3 * i + 1 ] = arrt[i];
   }


   double arr_dM_po_dX2[9] = {0.};

    QCoordSystTrsf::calc_dMtrx3_ASPK_v_PSK_po_dTet(arrX,arr_dM_po_dX2);
   MtrxTranspMultMatrx(arr_dM_po_dX2,3, 3, arr_e_psk,1, arrt) ;
   for (int i =0;i< 3; ++i)
   {
      arr_dJ_po_dX[3 * i + 2 ] = arrt[i];
   }

}
//------------------------------

bool QSolver2D_Angs::calc_arrNeviaz(double *arrX, const int NUmBigMeasure
                                       , double *arrNeviaz)
{
    QBigMeasure Meas = mVectBigMeasures.at(NUmBigMeasure);
   // 1. вычисление массива Y
    // первые 3 координаты - вектор положения антенны
    // затем 3 координаты - вектор положения маяка
    // затем 3 угла ориентации антенны
    double arrY [9] = {0.};

    // персчет вектора положения антенны в ГСК

    // создание матрицы перехода из   КГСК в ПСК
     double arrEilers[3] = {0.},  arr_KGSK[3] = {0.};
    MatrxMultScalar(Meas.marrMuZv, 1, 3, -1.,arrEilers);
    double matrPereh_PSK_V_KGSK[9] = {0} ;
    QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilers, matrPereh_PSK_V_KGSK) ;
    // вектор положения в ПСК-центр тяжести
    // вычисление вектора положениея в ПСК

    MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, marrPAnt_PSK, 1, arr_KGSK) ;

    MtrxSumMatrx(arr_KGSK, Meas.marrSVessZv,1, 3, arrY) ;
    // !



    arrY [3] = marrSBeaconXY_GSK[0]; // X маяка
    arrY [4] = marrSBeaconXY_GSK[1]; // Y маяка
    arrY [5] = mDeepth; // Z маяка
    arrY[6] = arrX[0];// угол
    arrY[7] = arrX[1];// угол
    arrY[8] = arrX[2];// угол

    // 1!

    // 2. вычисление массива Z
    double arrZ [4] = {0.};

   calc_arrZ(arrY, arrZ);
  // const double VAlTetta = arrZ[6];


    // 2!

    // 3. Вычисление массива I
    double arrI[6] = {0.};
    calc_arrI(arrZ, arrI);

    // 4. вектор направления в ГСК
    double *parr_e_gsk =arrI;
    // !4

    // 5. вектор направления в ПСК
    double arr_e_psk[3]={0.};
    MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, parr_e_gsk, 1, arr_e_psk) ;
    // !5



    // 6. вычисление вектора единичного направления в АСПК  J
    double arrJ[3] = {0.};
    calcVectJ_(arrX,  arr_e_psk,arrJ);

    // 6 !

    //  ОТЛАДКА
   // calcVectJ(Meas.marrMuZv, arrI, arrJ1);
    // ОТЛАДКА !

  // вычисление вектора V
  double arrV[2] ={0.};
  //arrV[0 ] = asin(arrJ[0]/ sqrt(arrJ[0] * arrJ[0] + arrJ[1] * arrJ[1]));
  arrV[0 ] = QPeaceVess::calcCourseAngle(arrJ[0], arrJ[1]);
  arrV[1 ] = atan2(arrJ[2],sqrt(arrJ[0] * arrJ[0] + arrJ[1] * arrJ[1]));

// ОТЛАДКА
 // double arrGSK_XYZ[3] = {0.};
 // arrGSK_XYZ[0] = arrX[3];
 // arrGSK_XYZ[1] = arrX[4];
 // arrGSK_XYZ[2] = mDeepth;
 // double arrAntPosParams[6] = {0.};
 // arrAntPosParams[0] = arrX[0];
//  arrAntPosParams[1] = arrX[1];
 // arrAntPosParams[2] = arrX[2];
 // arrAntPosParams[3] = arrX[5];
 // arrAntPosParams[4] = arrX[6];
 // arrAntPosParams[5] = arrX[7];
//  double pval_q = 0., pval_e = 0.,  pval_t= 0.;
//  transf_GSK_XYZ_to_USBL3D(mtblEstPrfl, arrGSK_XYZ, Meas.marrSVessZv, Meas.marrMuZv
 //   ,arrAntPosParams,  &pval_q, &pval_e,  &pval_t);
//  int uiu = 0;
  // ОТЛАДКА !
  arrNeviaz[0] = calcCorseNeviaz(arrV[0 ], Meas.mqzv);
  if(fabs(arrNeviaz[0])> 1.)
  {
      int iii = 0;
  }


    // ОТЛАДКА
//double arr_dJ_po_dX1_T[6 *3] = {0.}, arrJ1[6] = {0.},arrX1[3] = {0.}
 //          ,arr_dJ_po_dX1[6 *3] = {0.},arrdel[3] = {0.}
   //     , arrd[6]= {0.001, 0.001, 0.001, 0.0001, 0.0001, 0.0001};
 // for (int i =0; i < 3;++i)
  // {
   //   memcpy(arrX1,arrX, 3 * sizeof(double));
    //   arrX1[i] += arrd[i];
   // calcVectJ_(arrX1,  arr_e_psk,arrJ1);

   //  MtrxMinusMatrx(arrJ1, arrJ,1, 3, arrdel);
    //   MatrxMultScalar(arrdel, 1, 3, 1./arrd[i],&arr_dJ_po_dX1_T[i * 3]);
   // }
  // MatrTransp( arr_dJ_po_dX1_T, 3, 3, arr_dJ_po_dX1);
  // double ii = 0;// !
    // ОТЛАДКА !

   calc_arrNeviaz_3D(NUmBigMeasure, arrV[1], arrNeviaz );
    return true;
}

//------------------------------
bool QSolver2D_Angs::calc_arrNeviaz_3D( const int NUmBigMeasure
            ,const double e_est, double *arrNeviaz )
{

    return true;
}
//--------------------------------------------------
void QSolver2D_Angs::calc_dAngs_po_dSgps(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv, QVector<double> *pvct_dAngs_po_dSgps)
{
    double arr_dAngs_po_dSgps[9]={0.};
    QBigMeasure *parrBigMeasures = new QBigMeasure[mVectBigMeasures.size()];
    for (int i = 0; i < mVectBigMeasures.size(); ++i)
    {
      parrBigMeasures[i] =   mVectBigMeasures.at(i);
    }
    QLblSolver LblSolver(parrBigMeasures, mVectBigMeasures.size()    ,mtblEstPrfl
                 , mDeepth , 0.);
    delete []parrBigMeasures;
    //

    // производная функции FGr_ang по вектру положения GPS

    double  *parrMtrx_dFGrA_po_dSgps = new double[mDimX *3];
    memset (parrMtrx_dFGrA_po_dSgps, 0, mDimX *3 * sizeof(double));
    double *parrMtrx_dAng0_po_dSgps = new double [mDimX *3];
    calc_dFGrA_po_dSgps(arrAngX0, parrMtrx_dFGrA_po_dSgps);
    //MtrxMultMatrx( arrMtrx_dFGr_po_dXang_Inv,mDimX, mDimX, parrMtrx_dFGrA_po_dSgps,3, parrMtrx_dAng0_po_dSgps) ;
    //MatrxMultScalar(parrMtrx_dAng0_po_dSgps, 3, mDimX, -1.,parrMtrx_dAng0_po_dSgps);
    // !

    // матрица Якоби функции FGr по вектору позиционирования X0 [5]
    double  *parrMtrx_dFGrA_po_dX0pos = new double[mDimX *5];
    calc_dFGrA_po_dX0pos(arrAngX0, parrMtrx_dFGrA_po_dX0pos);
    // !

    // матрица Якоби вектора позиционирования по вектору положения GPS
    double *arrX0pos = new double [LblSolver.mDimX];
    arrX0pos[0] = marrPAnt_PSK[0];
    arrX0pos[1] = marrPAnt_PSK[1];
    arrX0pos[2] = marrPAnt_PSK[2];
    arrX0pos[3] = marrSBeaconXY_GSK[0];
    arrX0pos[4] = marrSBeaconXY_GSK[1];

    double *arrFGrPos = new double [LblSolver.mDimX];

    double *arrMtrx_dFGrPos_po_dX_Inv = new double [LblSolver.mDimX * LblSolver.mDimX];

    double *arr_dXpos_po_dSgps = new double [LblSolver.mDimX * 3];

    LblSolver.calc_arrFGr_and_Mtrx_dFGr_po_dX_Inv(arrX0pos,arrFGrPos,arrMtrx_dFGrPos_po_dX_Inv);

    LblSolver.calc_dX_po_dSgps(arrX0pos,arrMtrx_dFGrPos_po_dX_Inv,arr_dXpos_po_dSgps);
    // !

    // вычисление результирующей матрицы Якоби
    //  parrMtrx_dFGrA_po_dX0pos * arr_dXpos_po_dSgps
    double *arrTemp = new double [mDimX * 3];
    MtrxMultMatrx( parrMtrx_dFGrA_po_dX0pos,mDimX, 5, arr_dXpos_po_dSgps,3, arrTemp) ;

    double *parrMtrx_dFGrA_po_dSgps_full = new double [mDimX * 3];
    MtrxSumMatrx(arrTemp, parrMtrx_dFGrA_po_dSgps,mDimX, 3, parrMtrx_dFGrA_po_dSgps_full) ;

    MtrxMultMatrx( arrMtrx_dFGr_po_dXang_Inv,mDimX, 3, parrMtrx_dFGrA_po_dSgps_full,3, arr_dAngs_po_dSgps) ;

    MatrxMultScalar(arr_dAngs_po_dSgps, 3, mDimX, -1.,arr_dAngs_po_dSgps);

    QPosSolver::createDblVect(arr_dAngs_po_dSgps,3 *mDimX ,*pvct_dAngs_po_dSgps);

    delete []parrMtrx_dFGrA_po_dSgps;
    delete []parrMtrx_dFGrA_po_dX0pos;
    delete []arrX0pos;
    delete []arrFGrPos;
    delete []arrMtrx_dFGrPos_po_dX_Inv;
    delete []arr_dXpos_po_dSgps;
    delete []parrMtrx_dAng0_po_dSgps;
    delete []arrTemp;
    delete []parrMtrx_dFGrA_po_dSgps_full;
}
//------------------------------------------------------
bool QSolver2D_Angs::calc_dFGrA_po_dX0pos(double *arrAngX0,double * parrMtrx_dFGrA_po_dX0pos)
{
    memset(parrMtrx_dFGrA_po_dX0pos, 0, mDimX *5 * sizeof(double));

    const int QntMeas = mVectBigMeasures.size();


    double *parr_dfA_po_dXpos = new double [QntMeas * mDimX* 5];

    bool *barRreturn = new bool[QntMeas];
    memset(barRreturn, false,QntMeas * sizeof(bool));
                    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
  #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
 { // OMP (начало блока, который выполняется в нескольких потоках
   #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


    for (int i = 0; i < QntMeas; ++i)
    {


        if (!calc_arrdfiA_po_dXpos(arrAngX0
                                            ,i,&parr_dfA_po_dXpos[i * mDimX* 5])
               )
        {

            barRreturn[i] =  true;
        }


    }

 } // OMP (начало блока, который выполняется в нескольких потоках !
    for (int i =0; i < QntMeas; ++i)
    {
        if(barRreturn[i])
        {
            delete []parr_dfA_po_dXpos;
            delete []barRreturn;
            return false;
        }
    }
    for (int i = 0; i < QntMeas; ++i)
    {

        MtrxSumMatrx(parrMtrx_dFGrA_po_dX0pos, &parr_dfA_po_dXpos[i * mDimX* 5],mDimX, 5, parrMtrx_dFGrA_po_dX0pos) ;

    }

    delete []parr_dfA_po_dXpos;
    delete []barRreturn;

     return true ;
}
//-------------------------------------------
bool QSolver2D_Angs::calc_dFGrA_po_dSgps(double *arrAngX0,double * parrMtrx_dFGrA_po_dSgps)
{
   memset(parrMtrx_dFGrA_po_dSgps, 0, mDimX *3 * sizeof(double));

   const int QntMeas = mVectBigMeasures.size();


   double *parr_dfA_po_dSgps = new double [QntMeas * mDimX* 3];

   bool *barRreturn = new bool[QntMeas];
   memset(barRreturn, false,QntMeas * sizeof(bool));
                   //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
  #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
   #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


   for (int i = 0; i < QntMeas; ++i)
   {


       if (!calc_arrfi_and_arrdfi_po_dSgps(arrAngX0,i,&parr_dfA_po_dSgps[i * mDimX* 3]))

       {

           barRreturn[i] =  true;
       }


   }


}
   // OMP (начало блока, который выполняется в нескольких потоках !
   for (int i =0; i < QntMeas; ++i)
   {
       if(barRreturn[i])
       {
           delete []parr_dfA_po_dSgps;
           delete []barRreturn;
           return false;
       }
   }
   for (int i = 0; i < QntMeas; ++i)
   {

       MtrxSumMatrx(parrMtrx_dFGrA_po_dSgps, &parr_dfA_po_dSgps[i * mDimX* 3],mDimX, 3, parrMtrx_dFGrA_po_dSgps) ;

   }

   delete []parr_dfA_po_dSgps;
   delete []barRreturn;

    return true ;
}
//-----------------------------------
bool  QSolver2D_Angs::calc_arrdfiA_po_dXpos(double *arrX
                       ,const int NUmBigMeasure,double *arrdfiA_po_dXpos)
{
   double *arrNeviaz= new double [mDimY];
   double *arrfi = new double [mDimX];

    if(!calc_arrfi(arrX, NUmBigMeasure, arrNeviaz,arrfi))
    {
        delete []arrNeviaz;
        delete []arrfi;
        return false;
    }
     // вектор приращений
    double arrPrirashen[] = {0.01,0.01,0.01,0.01,0.01};
    double *arrfiCur = new double [mDimX];
    double *arrRazn = new double [mDimX];
    double *arrTemp = new double [5 * mDimX];
    double *arrNeviazCur = new double [mDimY];
    for (int i =0; i < 5; ++i)
    {

        if (i<3)
        {
        marrPAnt_PSK[i] += arrPrirashen[i];
        }
        else
        {
            marrSBeaconXY_GSK[i-3]  += arrPrirashen[i];
        }

      if (!calc_arrfi(arrX, NUmBigMeasure, arrNeviazCur,arrfiCur))
      {
          if (i<3)
          {
          marrPAnt_PSK[i] -= arrPrirashen[i];
          }
          else
          {
              marrSBeaconXY_GSK[i-3]  -= arrPrirashen[i];
          }
          delete []arrNeviaz;
          delete []arrfi;
          delete []arrfiCur;
          delete []arrRazn;
          delete [] arrTemp;
          delete []arrNeviazCur;
          return false;
      }
     MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
     MatrxMultScalar(arrRazn, 1, mDimX, 1./arrPrirashen[i],&arrTemp[mDimX * i]);
     if (i<3)
     {
     marrPAnt_PSK[i] -= arrPrirashen[i];
     }
     else
     {
         marrSBeaconXY_GSK[i-3]  -= arrPrirashen[i];
     }
    }

    MatrTransp( arrTemp, 5, mDimX, arrdfiA_po_dXpos);

 delete []arrfiCur;
 delete []arrRazn;
 delete [] arrTemp;
 delete []arrNeviazCur;
 delete []arrNeviaz;
 delete []arrfi;

 return true;

}
//----------------------------------------------
void QSolver2D_Angs::calc_dAngs_po_dHBeacon(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv
                                            ,QVector<double> &vct_dAngs_po_dH)
//,double *arr_dAngs_po_dH)
{
    double arr_dAngs_po_dH [3] ={0.};
    QBigMeasure *parrBigMeasures = new QBigMeasure[mVectBigMeasures.size()];
    for (int i = 0; i < mVectBigMeasures.size(); ++i)
    {
      parrBigMeasures[i] =   mVectBigMeasures.at(i);
    }
    QLblSolver LblSolver(parrBigMeasures, mVectBigMeasures.size()    ,mtblEstPrfl
                 , mDeepth , 0.);
    delete []parrBigMeasures;
    //

    // произвлдная функции FGr_ang по H

    double  *parrMtrx_dFGrA_po_dH = new double[mDimX ];
    memset (parrMtrx_dFGrA_po_dH, 0, mDimX  * sizeof(double));
    double *parrMtrx_dAng0_po_dH = new double [mDimX ];
    calc_dFGrA_po_dH(arrAngX0, parrMtrx_dFGrA_po_dH);

    // !

    // матрица Якоби функции FGr по вектору позиционирования X0 [5]
    double  *parrMtrx_dFGrA_po_dX0pos = new double[mDimX *5];
    calc_dFGrA_po_dX0pos(arrAngX0, parrMtrx_dFGrA_po_dX0pos);
    // !

    // матрица Якоби вектора позицимонирования по H
    double *arrX0pos = new double [LblSolver.mDimX];
    arrX0pos[0] = marrPAnt_PSK[0];
    arrX0pos[1] = marrPAnt_PSK[1];
    arrX0pos[2] = marrPAnt_PSK[2];
    arrX0pos[3] = marrSBeaconXY_GSK[0];
    arrX0pos[4] = marrSBeaconXY_GSK[1];

    double *arrFGrPos = new double [LblSolver.mDimX];

    double *arrMtrx_dFGrPos_po_dX_Inv = new double [LblSolver.mDimX * LblSolver.mDimX];

    double *arr_dXpos_po_dH = new double [LblSolver.mDimX ];

    LblSolver.calc_arrFGr_and_Mtrx_dFGr_po_dX_Inv(arrX0pos,arrFGrPos,arrMtrx_dFGrPos_po_dX_Inv);

    LblSolver.calc_dX_po_dHBeacon(arrX0pos,arrMtrx_dFGrPos_po_dX_Inv,arr_dXpos_po_dH);
    // !

    // вычисление результирующей матрицы Якоби
    //  parrMtrx_dFGrA_po_dX0pos * arr_dXpos_po_dSgps
    double *arrTemp = new double [mDimX];
    MtrxMultMatrx( parrMtrx_dFGrA_po_dX0pos,mDimX, 5, arr_dXpos_po_dH,1, arrTemp) ;

    double *parrMtrx_dFGrA_po_dH_full = new double [mDimX ];
    MtrxSumMatrx(arrTemp, parrMtrx_dFGrA_po_dH,mDimX, 1, parrMtrx_dFGrA_po_dH_full) ;

    MtrxMultMatrx( arrMtrx_dFGr_po_dXang_Inv,mDimX, mDimX, parrMtrx_dFGrA_po_dH_full,1, arr_dAngs_po_dH) ;

    MatrxMultScalar(arr_dAngs_po_dH, mDimX, mDimX, -1.,arr_dAngs_po_dH);

    QPosSolver::createDblVect(arr_dAngs_po_dH, 3, vct_dAngs_po_dH);

    delete []parrMtrx_dFGrA_po_dH;
    delete []parrMtrx_dFGrA_po_dX0pos;
    delete []arrX0pos;
    delete []arrFGrPos;
    delete []arrMtrx_dFGrPos_po_dX_Inv;
    delete []arr_dXpos_po_dH;
    delete []parrMtrx_dAng0_po_dH;
    delete []arrTemp;
    delete []parrMtrx_dFGrA_po_dH_full;
}
//------------------------------------
//-------------------------------------------
bool QSolver2D_Angs::calc_dFGrA_po_dH(double *arrAngX0,double * parrMtrx_dFGrA_po_dH)
{
   memset(parrMtrx_dFGrA_po_dH, 0, mDimX  * sizeof(double));

   const int QntMeas = mVectBigMeasures.size();


   double *parr_dfA_po_dH = new double [QntMeas * mDimX];

   bool *barRreturn = new bool[QntMeas];
   memset(barRreturn, false,QntMeas * sizeof(bool));
                   //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
// ! #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
// !  #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


   for (int i = 0; i < QntMeas; ++i)
   {


       if (!calc_arrdfiA_po_dH(arrAngX0
                                           ,i,&parr_dfA_po_dH[i * mDimX])
              )
       {

           barRreturn[i] =  true;
       }


   }


} // OMP (начало блока, который выполняется в нескольких потоках !
   for (int i =0; i < QntMeas; ++i)
   {
       if(barRreturn[i])
       {
           delete []parr_dfA_po_dH;
           delete []barRreturn;
           return false;
       }
   }
   for (int i = 0; i < QntMeas; ++i)
   {

       MtrxSumMatrx(parrMtrx_dFGrA_po_dH, &parr_dfA_po_dH[i * mDimX],mDimX, 1, parrMtrx_dFGrA_po_dH) ;

   }

   delete []parr_dfA_po_dH;
   delete []barRreturn;

    return true ;
}
//-------------------------------------------------
bool  QSolver2D_Angs::calc_arrdfiA_po_dH(double *arrX
                       ,const int NUmBigMeasure,double *arrdfi_po_dH)
{
   double *arrNeviaz= new double [mDimY];
   double *arrfi = new double [mDimX];

    if(!calc_arrfi(arrX, NUmBigMeasure, arrNeviaz,arrfi))
    {
        delete []arrNeviaz;
        delete []arrfi;
        return false;
    }
     // вектор приращений
    double valPrirashen = 0.01;

    double *arrfiCur = new double [mDimX];
    double *arrRazn = new double [mDimX];
    double *arrTemp = new double [ mDimX];
    double *arrNeviazCur = new double [mDimY];


      mDeepth += valPrirashen;
        if (!calc_arrfi(arrX, NUmBigMeasure, arrNeviazCur,arrfiCur))
      {
          delete []arrNeviaz;
          delete []arrfi;
          delete []arrfiCur;
          delete []arrRazn;
          delete [] arrTemp;
          delete []arrNeviazCur;

          mDeepth -= valPrirashen;
          return false;
      }
     MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
     MatrxMultScalar(arrRazn, 1, mDimX, 1./valPrirashen,arrdfi_po_dH);
     mDeepth -= valPrirashen;


 delete []arrfiCur;
 delete []arrRazn;
 delete [] arrTemp;
 delete []arrNeviazCur;
 delete []arrNeviaz;
 delete []arrfi;
 return true;
}
//------------------------------------
// Вычислениекореляционной матрицы ошибок рассеяний оценки вектора X
// вызванных ошибками измерения GPS c сигмой 1 м
//
//
bool  QSolver2D_Angs:: calc_K_qe(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv,QVector<double> &vct_K_dqe)
{
    double arr_K_dqe[9] = {0.};
    memset(arr_K_dqe, 0, mDimX *mDimX * sizeof(double));

    const int QntMeas = mVectBigMeasures.size();


    double *parr_df_po_dqe= new double [QntMeas * mDimX* mDimY];

    bool *barRreturn = new bool[QntMeas];
    memset(barRreturn, false,QntMeas * sizeof(bool));
                    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
 #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
 { // OMP (начало блока, который выполняется в нескольких потоках
 #pragma omp for // OMP (директива для распределения итераций цикла между потоками)


    for (int i = 0; i < QntMeas; ++i)
    {


        if (!calc_arrfi_and_arrdfi_po_dqe(arrAngX0
                                            ,i,&parr_df_po_dqe[i * mDimX* mDimY])
               )
        {

             barRreturn[i] =  true;
        }

    }


 }
// OMP (начало блока, который выполняется в нескольких потоках !

    for (int i =0; i < QntMeas; ++i)
    {
        if(barRreturn[i])
        {

            delete []parr_df_po_dqe;
            delete []barRreturn;
            return false;
        }
    }


    double *arrTemp0 = new double [mDimX *mDimX];
    double *arrTemp1 = new double [mDimX *mDimX];
    memset(arrTemp1, 0, sizeof(double) *mDimX *mDimX);

    double *arrMeasure = new double [mDimY];
    double *arrD = new double [mDimY *mDimY];
    for (int i = 0; i < QntMeas; ++i)
    {
        QBigMeasure Meas = mVectBigMeasures.at(i);
        collectGrlsMeasure(Meas ,arrMeasure,arrD);

        double *arrT  = new double [mDimX *mDimY];
        MtrxMultMatrx(&parr_df_po_dqe[i * mDimX* mDimY],mDimX , mDimY, arrD,mDimY, arrT ) ;
        MtrxMultMatrxTransp(arrT,mDimX , mDimY, &parr_df_po_dqe[i * mDimX* mDimY],mDimX, arrTemp0) ;
        delete []arrT;

        MtrxSumMatrx(arrTemp1, arrTemp0,mDimX, mDimX, arrTemp1) ;
    }
    delete []arrMeasure;
    delete []parr_df_po_dqe;
    delete []barRreturn;
    delete []arrD;


    calcF_D_FTransp(arrMtrx_dFGr_po_dXang_Inv,arrTemp1,mDimX, arr_K_dqe);

    QPosSolver::createDblVect(arr_K_dqe, 9, vct_K_dqe);
    delete []arrTemp0;
    delete []arrTemp1;

     return true ;
}
//----------------------------

bool  QSolver2D_Angs::calc_arrfi_and_arrdfi_po_dqe(double *arrAngX0
                                            ,const int NUmBigMeasure,double *arr_df_po_dqe)
{
    double *arrNeviaz= new double [mDimY];
    double *arrfi = new double [mDimX];

     if(!calc_arrfi(arrAngX0, NUmBigMeasure, arrNeviaz,arrfi))
     {
         delete []arrNeviaz;
         delete []arrfi;
         return false;
     }
      // вектор приращений
     double valPrirashen = 0.0001;

     double *arrfiCur = new double [mDimX];
     double *arrRazn = new double [mDimX];
     double *arrTemp = new double [mDimY * mDimX];
     double *arrNeviazCur = new double [mDimY];

     for (int i =0; i < mDimY;++i)
     {
         QBigMeasure meas = mVectBigMeasures.at(i);
         if (i ==0)
         {
             meas.mqzv +=valPrirashen;
         }
         else
         {
            meas.mezv +=valPrirashen;
         }
         mVectBigMeasures.replace(NUmBigMeasure, meas);
         if (!calc_arrfi(arrAngX0, NUmBigMeasure, arrNeviazCur,arrfiCur))
       {
           delete []arrNeviaz;
           delete []arrfi;
           delete []arrfiCur;
           delete []arrRazn;
           delete [] arrTemp;
           delete []arrNeviazCur;

             if (i ==0)
             {
                 meas.mqzv -=valPrirashen;
             }
             else
             {
                meas.mezv -=valPrirashen;
             }
             mVectBigMeasures.replace(NUmBigMeasure, meas);
           return false;
       }
      MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
      MatrxMultScalar(arrRazn, 1, mDimX, 1./valPrirashen,&arrTemp[i * mDimX]);
      if (i ==0)
      {
          meas.mqzv -=valPrirashen;
      }
      else
      {
         meas.mezv -=valPrirashen;
      }
      mVectBigMeasures.replace(NUmBigMeasure, meas);

     }

     MatrTransp( arrTemp, mDimY, mDimX, arr_df_po_dqe);

  delete []arrfiCur;
  delete []arrRazn;
  delete [] arrTemp;
  delete []arrNeviazCur;
  delete []arrNeviaz;
  delete []arrfi;
  return true;
}
//--------------------------------------
bool  QSolver2D_Angs:: calc_dAngs_po_dAngSins(double* arrAngX0,double* arrMtrx_dFGr_po_dXang_Inv
                                              ,QVector<double> &vct_dXa_po_AngSins)
{
    double arr_dXa_po_AngSins[9] = {0.};

    QBigMeasure *parrBigMeasures = new QBigMeasure[mVectBigMeasures.size()];
    for (int i = 0; i < mVectBigMeasures.size(); ++i)
    {
      parrBigMeasures[i] =   mVectBigMeasures.at(i);
    }
    QLblSolver LblSolver(parrBigMeasures, mVectBigMeasures.size()    ,mtblEstPrfl
                 , mDeepth , 0.);
    delete []parrBigMeasures;
    //

    // производная функции FGr_ang по вектору сдвига ИС

    double  *parrMtrx_dFGrA_po_dSins = new double[mDimX *3];
    memset (parrMtrx_dFGrA_po_dSins, 0, mDimX *3 * sizeof(double));
    double *parrMtrx_dAng0_po_dSins = new double [mDimX *3];
   // calc_dFGrA_po_dSins(arrAngX0, parrMtrx_dFGrA_po_dSins);



    calc_dFGr_po_dSins_(arrAngX0, parrMtrx_dFGrA_po_dSins);
    //MtrxMultMatrx( arrMtrx_dFGr_po_dXang_Inv,mDimX, mDimX, parrMtrx_dFGrA_po_dSgps,3, parrMtrx_dAng0_po_dSgps) ;
    //MatrxMultScalar(parrMtrx_dAng0_po_dSgps, 3, mDimX, -1.,parrMtrx_dAng0_po_dSgps);
    // !

    // матрица Якоби функции FGr по вектору позиционирования X0 [5]
    double  *parrMtrx_dFGrA_po_dX0pos = new double[mDimX *5];
    calc_dFGrA_po_dX0pos(arrAngX0, parrMtrx_dFGrA_po_dX0pos);


    //double arrMtrx_dFusbl_po_dXlbl[15] = {0.};
   // calc_dFusbl_po_dXlbl(arrAngX0,arrMtrx_dFusbl_po_dXlbl);
    // !

    // матрица Якоби вектора позицимонирования по вектору положения GPS
    double *arrX0pos = new double [LblSolver.mDimX];
    arrX0pos[0] = marrPAnt_PSK[0];
    arrX0pos[1] = marrPAnt_PSK[1];
    arrX0pos[2] = marrPAnt_PSK[2];
    arrX0pos[3] = marrSBeaconXY_GSK[0];
    arrX0pos[4] = marrSBeaconXY_GSK[1];

    double *arrFGrPos = new double [LblSolver.mDimX];

    double *arrMtrx_dFGrPos_po_dX_Inv = new double [LblSolver.mDimX * LblSolver.mDimX];

    double *arr_dXpos_po_dSins = new double [LblSolver.mDimX * 3];

    LblSolver.calc_arrFGr_and_Mtrx_dFGr_po_dX_Inv(arrX0pos,arrFGrPos,arrMtrx_dFGrPos_po_dX_Inv);

    LblSolver.calc_dX_po_dAngSins(arrX0pos,arrMtrx_dFGrPos_po_dX_Inv,arr_dXpos_po_dSins);
    // !

    // вычисление результирующей матрицы Якоби

    double *arrTemp = new double [mDimX * 3];
    double *arrTemp1 = new double [mDimX * 3];
    double *arrTemp2 = new double [mDimX * 3];
    MtrxMultMatrx( parrMtrx_dFGrA_po_dX0pos,mDimX, 5, arr_dXpos_po_dSins,3, arrTemp) ;
    MtrxMultMatrx( arrMtrx_dFGr_po_dXang_Inv,mDimX, mDimX, arrTemp,3, arrTemp1) ;
    MatrxMultScalar(arrTemp1, 3, mDimX, -1.,arrTemp1);

    MtrxMultMatrx( arrMtrx_dFGr_po_dXang_Inv,mDimX, mDimX, parrMtrx_dFGrA_po_dSins,3, arrTemp2) ;
    MatrxMultScalar(arrTemp2, 3, mDimX, -1.,arrTemp2);

    MtrxSumMatrx(arrTemp1, arrTemp2,mDimX, 3, arr_dXa_po_AngSins) ;

    QPosSolver::createDblVect(arr_dXa_po_AngSins, 9,vct_dXa_po_AngSins);

    delete []parrMtrx_dFGrA_po_dSins;
    delete []parrMtrx_dFGrA_po_dX0pos;
    delete []arrX0pos;
    delete []arrFGrPos;
    delete []arrMtrx_dFGrPos_po_dX_Inv;
    delete []arr_dXpos_po_dSins;
    delete []parrMtrx_dAng0_po_dSins;
    delete []arrTemp;


    delete []arrTemp1;
    delete []arrTemp2;
}
//----------------------------------------
bool  QSolver2D_Angs:: calc_dAngs_po_dProfile(double* arrAngX0,double* arrMtrx_dFGr_po_dXang_Inv
                                              ,QVector<double> &vct_dAngs_po_dProfile)
{
    QBigMeasure *parrBigMeasures = new QBigMeasure[mVectBigMeasures.size()];
    for (int i = 0; i < mVectBigMeasures.size(); ++i)
    {
      parrBigMeasures[i] =   mVectBigMeasures.at(i);
    }
    QLblSolver LblSolver(parrBigMeasures, mVectBigMeasures.size()    ,mtblEstPrfl
                 , mDeepth , 0.);
    delete []parrBigMeasures;
    //

    // производная функции FGr_ang по вектору систематической ошибки профиля

    double  *parrMtrx_dFGrA_po_dProfile = new double[mDimX];
    memset (parrMtrx_dFGrA_po_dProfile, 0, mDimX * sizeof(double));
    double *parrMtrx_dAng0_po_dProfile = new double [mDimX];
    calc_dFGr_po_dProfile(arrAngX0, parrMtrx_dFGrA_po_dProfile);
    //MtrxMultMatrx( arrMtrx_dFGr_po_dXang_Inv,mDimX, mDimX, parrMtrx_dFGrA_po_dSgps,3, parrMtrx_dAng0_po_dSgps) ;
    //MatrxMultScalar(parrMtrx_dAng0_po_dSgps, 3, mDimX, -1.,parrMtrx_dAng0_po_dSgps);
    // !

    // матрица Якоби функции FGr по вектору позиционирования X0 [5]
    double  *parrMtrx_dFGrA_po_dX0pos = new double[mDimX *5];
    calc_dFGrA_po_dX0pos(arrAngX0, parrMtrx_dFGrA_po_dX0pos);
    // !

    // матрица Якоби вектора позицимонирования по вектору систематической ошибки профиля
    double *arrX0pos = new double [LblSolver.mDimX];
    arrX0pos[0] = marrPAnt_PSK[0];
    arrX0pos[1] = marrPAnt_PSK[1];
    arrX0pos[2] = marrPAnt_PSK[2];
    arrX0pos[3] = marrSBeaconXY_GSK[0];
    arrX0pos[4] = marrSBeaconXY_GSK[1];

    double *arrFGrPos = new double [LblSolver.mDimX];

    double *arrMtrx_dFGrPos_po_dX_Inv = new double [LblSolver.mDimX * LblSolver.mDimX];

    double *arr_dXpos_po_dProfile = new double [LblSolver.mDimX ];

    LblSolver.calc_arrFGr_and_Mtrx_dFGr_po_dX_Inv(arrX0pos,arrFGrPos,arrMtrx_dFGrPos_po_dX_Inv);

    LblSolver.calc_dX_po_dProfile(arrX0pos,arrMtrx_dFGrPos_po_dX_Inv,arr_dXpos_po_dProfile);
    // !

    // вычисление результирующей матрицы Якоби
    //  parrMtrx_dFGrA_po_dX0pos * arr_dXpos_po_dSgps
    double *arrTemp = new double [mDimX ];
    MtrxMultMatrx( parrMtrx_dFGrA_po_dX0pos,mDimX, 5, arr_dXpos_po_dProfile,1, arrTemp) ;

    double *parrMtrx_dFGrA_po_dProfile_full = new double [mDimX];
    MtrxSumMatrx(arrTemp, parrMtrx_dFGrA_po_dProfile,mDimX, 1, parrMtrx_dFGrA_po_dProfile_full) ;

    double arr_dAngs_po_dProfile[3] = {0.};
    MtrxMultMatrx( arrMtrx_dFGr_po_dXang_Inv,mDimX, 1, parrMtrx_dFGrA_po_dProfile_full,1, arr_dAngs_po_dProfile) ;

    MatrxMultScalar(arr_dAngs_po_dProfile, 1, mDimX, -1.,arr_dAngs_po_dProfile);

    QPosSolver::createDblVect(arr_dAngs_po_dProfile, 3,vct_dAngs_po_dProfile);

    delete []parrMtrx_dFGrA_po_dProfile;
    delete []parrMtrx_dFGrA_po_dX0pos;
    delete []arrX0pos;
    delete []arrFGrPos;
    delete []arrMtrx_dFGrPos_po_dX_Inv;
    delete []arr_dXpos_po_dProfile;
    delete []parrMtrx_dAng0_po_dProfile;
    delete []arrTemp;
    delete []parrMtrx_dFGrA_po_dProfile_full;
}

//-------------------------------------------------
bool  QSolver2D_Angs::calc_arrdfiA_po_dProfile(double *arrX
                       ,const int NUmBigMeasure,double *arrdfi_po_dProfile)
{
   double *arrNeviaz= new double [mDimY];
   double *arrfi = new double [mDimX];

    if(!calc_arrfi(arrX, NUmBigMeasure, arrNeviaz,arrfi))
    {
        delete []arrNeviaz;
        delete []arrfi;
        return false;
    }
     // вектор приращений
    double valPrirashen= 0.1;
    QSolver2D_Angs Solver2D_AngsCur = *this;
    for (int i = 0; i <mtblEstPrfl.mNumCols;++i )
    {
       Solver2D_AngsCur.mtblEstPrfl.mparrVal[i] +=  valPrirashen;
    }
    double *arrNeviazCur = new double [mDimY];
    double *arrfiCur = new double [mDimX];
    if(!Solver2D_AngsCur.calc_arrfi(arrX, NUmBigMeasure, arrNeviazCur,arrfiCur))
    {
        delete []arrNeviaz;
        delete []arrfi;
        delete []arrNeviazCur;
        delete []arrfiCur;
        return false;
    }



    double *arrRazn = new double [mDimX];
    double *arrTemp = new double [ mDimX];


    MtrxMinusMatrx(arrfiCur, arrfi,1, mDimX, arrRazn);
    MatrxMultScalar(arrRazn, 1, mDimX, 1./valPrirashen,arrTemp);


 delete []arrfiCur;
 delete []arrRazn;
 delete [] arrTemp;
 delete []arrNeviazCur;
 delete []arrNeviaz;
 delete []arrfi;
 return true;

}


//----------------------------------------------------
bool  QSolver2D_Angs::calc_AngKgps(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv
                                ,double *arrMtrx_dFGrLBL_po_dX_Inv,QVector<double> &vct__Kgps_per_1m)

{
    // 1. Вычисление dF_po_dP
    const int dimP = 5;
    double *arr_dF_po_dP = new double[mDimX *dimP];
    calc_dFusbl_po_dXlbl(arrAngX0,arr_dF_po_dP);
    //1 !

    // 2. Вычисление матриц dP_po_dai
    QLblSolver LBLsolv(mVectBigMeasures
               ,mtblEstPrfl, mDeepth, 0.1 );
    double *parr_dfiLBL_po_dai = new double [mVectBigMeasures.size() * dimP * 3];
    double *arrXlbl = new double [LBLsolv.mDimX];

    memcpy(arrXlbl, marrPAnt_PSK,3 * sizeof(double));
    memcpy(&arrXlbl[3], marrSBeaconXY_GSK,2 * sizeof(double));

    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
 #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
 #pragma omp for // OMP (директива для распределения итераций цикла между потоками)
        for (int i =0; i < mVectBigMeasures.size();++i)
        {
        LBLsolv.calc_darrfi_po_dSvess(arrXlbl, i
                                               ,  0.1, &parr_dfiLBL_po_dai[i * dimP * 3]);
        }
}
    //2!

// 3. Вычисление матриц dfi_po_dai
    double *parr_dfi_po_dai = new double [mVectBigMeasures.size() * mDimX *  3];
    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
 #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
 #pragma omp for // OMP (директива для распределения итераций цикла между потоками)
    for (int i =0; i < mVectBigMeasures.size();++i)
    {
    calc_darrfi_po_dSvess(arrAngX0, i
                                           ,  0.1, &parr_dfi_po_dai[i * mDimX * 3]);
    }
 }
    //3!

    //4. основной цикл

    double *arr_Sum = new double [ mDimX * mDimX];
    memset(arr_Sum, 0, sizeof(double) * mDimX * mDimX);

    double  arrGpsDisp[9] = {1.,0.,0.
                              ,0.,1.,0.
                              ,0.,0.,1.};
    for (int i = 0; i < mVectBigMeasures.size();++i)
    {
        double *parrT0 = new double [dimP * 3];
        double *parrT1 = new double [mDimX * 3];
        double *parrT2 = new double [mDimX * mDimX];
        MtrxMultMatrx( arrMtrx_dFGrLBL_po_dX_Inv,dimP, dimP, &parr_dfiLBL_po_dai[i * dimP * 3],3, parrT0);
        MtrxMultMatrx( arr_dF_po_dP,mDimX, dimP, parrT0,3, parrT1);
        MtrxSumMatrx(parrT1, &parr_dfi_po_dai[i * mDimX * 3],mDimX, 3, parrT1) ;
        calcF_Mult_FTransp(parrT1,mDimX, 3, parrT2);
        //calcF_D_FTransp_(parrT1,mDimX,3,arrGpsDisp, parrT2);
        MtrxSumMatrx(arr_Sum, parrT2,mDimX, mDimX, arr_Sum) ;
        delete []parrT0;
        delete []parrT1;
        delete []parrT2;
    }

    double arr_Kgps_per_1m[9] = {0.};
    calcF_D_FTransp_(arrMtrx_dFGr_po_dXang_Inv,mDimX,mDimX,arr_Sum,arr_Kgps_per_1m);

    QPosSolver::createDblVect(arr_Kgps_per_1m, 9, vct__Kgps_per_1m);


    delete []arr_dF_po_dP;

    delete []parr_dfiLBL_po_dai;
    delete []arrXlbl;
    delete []parr_dfi_po_dai;
    delete []arr_Sum;
}
//-----------------------------------------------------------

bool  QSolver2D_Angs::calc_AngKsinsQ(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv
                              ,double *arrMtrx_dFGrLBL_po_dX_Inv ,QVector<double> &vct_AngK_SinsQ)

{
    // 1. Вычисление dF_po_dP
        const int dimP = 5;
        double *arr_dF_po_dP = new double[mDimX *dimP];
        calc_dFusbl_po_dXlbl(arrAngX0,arr_dF_po_dP);
    //1 !

    // 2. Вычисление матриц dP_po_dQi
    QLblSolver LBLsolv(mVectBigMeasures
               ,mtblEstPrfl, mDeepth, 0.1 );
    double *parr_dfiLBL_po_dQi = new double [mVectBigMeasures.size() * dimP];
    double *arrXlbl = new double [LBLsolv.mDimX];


    memcpy(arrXlbl, marrPAnt_PSK,3 * sizeof(double));
    memcpy(&arrXlbl[3], marrSBeaconXY_GSK,2 * sizeof(double));
    double del = 0.001;
    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
 #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
 #pragma omp for // OMP (директива для распределения итераций цикла между потоками)
        for (int i =0; i < mVectBigMeasures.size();++i)
        {

        LBLsolv.calc_arrfi_and_arrdfi_po_dQ(arrXlbl, i, del,&parr_dfiLBL_po_dQi[i * dimP]);

        }
}
    //2!

// 3. Вычисление матриц dfi_po_dai
    double *parr_dfi_po_dQi = new double [mVectBigMeasures.size() * mDimX];

    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
  #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
 #pragma omp for // OMP (директива для распределения итераций цикла между потоками)
    for (int i =0; i < mVectBigMeasures.size();++i)
    {

        calc_arrfi_and_arrdfi_po_dQ(arrAngX0, i, del,&parr_dfi_po_dQi[i * mDimX ]);

    }
 }
    //3!

    //4. основной цикл
    double *arr_Sum = new double [ mDimX * mDimX];
    memset(arr_Sum, 0, sizeof(double) *  mDimX * mDimX);
    for (int i = 0; i < mVectBigMeasures.size();++i)
    {
        double *parrT0 = new double [dimP ];
        double *parrT1 = new double [mDimX ];
        double *parrT2 = new double [mDimX * mDimX];
        MtrxMultMatrx( arrMtrx_dFGrLBL_po_dX_Inv,dimP, dimP, &parr_dfiLBL_po_dQi[i * dimP],1, parrT0);
        MtrxMultMatrx( arr_dF_po_dP,mDimX, dimP, parrT0,1, parrT1);
        MtrxSumMatrx(parrT1, &parr_dfi_po_dQi[i * mDimX ],mDimX, 1, parrT1) ;
        double disp[1] = {1.};
        calcF_D_FTransp_(parrT1,mDimX,1,disp, parrT2);
        MtrxSumMatrx(arr_Sum, parrT2,mDimX, mDimX, arr_Sum) ;
        delete []parrT0;
        delete []parrT1;
        delete []parrT2;
    }

    double arr_AngK_SinsQ[9] = {0.};
    calcF_D_FTransp_(arrMtrx_dFGr_po_dXang_Inv,mDimX,mDimX,arr_Sum, arr_AngK_SinsQ);

    QPosSolver::createDblVect(arr_AngK_SinsQ, 9, vct_AngK_SinsQ);


    delete []arr_dF_po_dP;
    delete []parr_dfiLBL_po_dQi;
    delete []arrXlbl;


     delete []parr_dfi_po_dQi;
    delete []arr_Sum;
}
//------------------------------------------------------------
bool  QSolver2D_Angs::calc_dFusbl_po_dXlbl(double *arrAngX0,double *arrMtrx_dFusbl_po_dXlbl)
{
    // 1. Вычисление dF_po_dP
    const int dimP = 5;

    double *arr_dF_po_dP_T = new double[dimP * mDimX];
    double *arr_F0 = new double[mDimX ];
    double *arr_FCur= new double[mDimX];
    double * parrMtrx_dF_po_dX = new double [mDimX*dimP];
    double valSumVeviaz2 = 0.;
    calc_FGr_and_dFGr_po_dX(arrAngX0
                       ,arr_F0,parrMtrx_dF_po_dX, valSumVeviaz2);
    double del = 0.1;
    for (int i=0; i <dimP; ++i)
    {
        if (i < 3)
        {
        marrPAnt_PSK[i] +=del;
        }
        else
        {
          marrSBeaconXY_GSK[i -3] +=del;
        }
        calc_FGr_and_dFGr_po_dX(arrAngX0
                           ,arr_FCur,parrMtrx_dF_po_dX, valSumVeviaz2);
        MtrxMinusMatrx(arr_FCur, arr_F0,1, mDimX, arr_FCur);
        MatrxMultScalar(arr_FCur, 1, mDimX, 1./del,&arr_dF_po_dP_T[i *mDimX]);

        if (i < 3)
        {
        marrPAnt_PSK[i] -=del;
        }
        else
        {
          marrSBeaconXY_GSK[i -3] -=del;
        }
    }
    MatrTransp( arr_dF_po_dP_T, dimP, mDimX, arrMtrx_dFusbl_po_dXlbl);
    //1 !

    delete []arr_dF_po_dP_T;
    delete []arr_F0;
    delete []arr_FCur;
    delete []parrMtrx_dF_po_dX;

}
//-----------------------------------------
// Вычисление кореляционной матрицы ошибок рассеяний оценки вектора X
// вызванных ошибками измерения GPS c сигмой 1 м
bool  QSolver2D_Angs::calc_AngKsins_Psi_Tetta(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv
                                 ,double *arrMtrx_dFGrLBL_po_dX_Inv,QVector<double> &vct_AngKgps_Psi_Tetta_per_1_rad)
                                              //double *arr_AngKgps_Psi_Tetta_per_1_rad)

{

    // 1. Вычисление dF_po_dP
        const int dimP = 5;
        double *arr_dF_po_dP = new double[mDimX *dimP];
        calc_dFusbl_po_dXlbl(arrAngX0,arr_dF_po_dP);
    //1 !

    // 2. Вычисление матриц dP_po_dPsiTet
    QLblSolver LBLsolv(mVectBigMeasures
               ,mtblEstPrfl, mDeepth, 0.1 );
    double *parr_dfiLBL_po_dPsiTeti = new double [mVectBigMeasures.size() * dimP *2];
    double *arrXlbl = new double [LBLsolv.mDimX];


    memcpy(arrXlbl, marrPAnt_PSK,3 * sizeof(double));
    memcpy(&arrXlbl[3], marrSBeaconXY_GSK,2 * sizeof(double));
    double arrPrirashen[2] = {0.001,0.001};
    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
 #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
 #pragma omp for // OMP (директива для распределения итераций цикла между потоками)
        for (int i =0; i < mVectBigMeasures.size();++i)
        {
        LBLsolv.calc_arrfi_and_arrdfi_po_dPsiTetta(arrXlbl, i,arrPrirashen
                                                   ,&parr_dfiLBL_po_dPsiTeti[i * dimP *2]);

        }
}
    //2!

// 3. Вычисление матриц dfi_po_dai
    double *parr_dfi_po_dPsiTeti = new double [mVectBigMeasures.size() * mDimX *2];

    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
 #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
 #pragma omp for // OMP (директива для распределения итераций цикла между потоками)
    for (int i =0; i < mVectBigMeasures.size();++i)
    {
        calc_arrfi_and_arrdfi_po_dPsiTetta(arrAngX0, i,arrPrirashen
                                                           ,&parr_dfi_po_dPsiTeti[i * mDimX *2]);

    }
 }
    //3!

    //4. основной цикл

    double *arr_Sum = new double [ mDimX * mDimX];
    memset(arr_Sum, 0, sizeof(double) *  mDimX * mDimX);
    for (int i = 0; i < mVectBigMeasures.size();++i)
    {
        double *parrT0 = new double [dimP *2 ];
        double *parrT1 = new double [mDimX *2];
        double *parrT2 = new double [mDimX * mDimX];
        MtrxMultMatrx( arrMtrx_dFGrLBL_po_dX_Inv,dimP, dimP, &parr_dfiLBL_po_dPsiTeti[i * dimP *2],2, parrT0);
        MtrxMultMatrx( arr_dF_po_dP,mDimX, dimP, parrT0,2, parrT1);
        MtrxSumMatrx(parrT1, &parr_dfi_po_dPsiTeti[i * mDimX *2],mDimX, 2, parrT1) ;
        double disp[4] = {1., 0., 0.,1.};
        calcF_D_FTransp_(parrT1,mDimX,2,disp, parrT2);
        MtrxSumMatrx(arr_Sum, parrT2,mDimX, mDimX, arr_Sum) ;
        delete []parrT0;
        delete []parrT1;
        delete []parrT2;
    }

    double arr_AngKgps_Psi_Tetta_per_1_rad[9] = {0.};
    calcF_D_FTransp_(arrMtrx_dFGr_po_dXang_Inv,mDimX,mDimX,arr_Sum, arr_AngKgps_Psi_Tetta_per_1_rad);

    QPosSolver::createDblVect(arr_AngKgps_Psi_Tetta_per_1_rad,9,vct_AngKgps_Psi_Tetta_per_1_rad);

    delete []arr_dF_po_dP;
    delete []parr_dfiLBL_po_dPsiTeti;
    delete []arrXlbl;   
    delete []parr_dfi_po_dPsiTeti;
    delete []arr_Sum;
    return true;
}
//--------------------------------------------------------
bool  QSolver2D_Angs::calc_AngKt(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv
                                 ,double *arrMtrx_dFGrLBL_po_dX_Inv ,QVector<double> &vct_AngKt_per_1sec)
{
    qDebug() << "IN: calc_AngKt";
     // 1. Вычисление dF_po_dP (по вектору позиционирования)
     const int dimP = 5;
     double *arr_dF_po_dP = new double[mDimX *dimP];
     memset(arr_dF_po_dP, 0, mDimX *dimP * sizeof(double));
    calc_dFusbl_po_dXlbl(arrAngX0,arr_dF_po_dP);
    //1 !

    // 2. Вычисление матриц dP_po_dt
     QBigMeasure *parrBigMeasures = new QBigMeasure[mVectBigMeasures.size()];
     for (int i = 0; i < mVectBigMeasures.size(); ++i)
     {
       parrBigMeasures[i] =   mVectBigMeasures.at(i);
     }
     QLblSolver LBLsolv(parrBigMeasures, mVectBigMeasures.size()    ,mtblEstPrfl
                  , mDeepth , 0.);
     delete []parrBigMeasures;
     //

    qDebug() << "calc_AngKt: LBLsolv CREATED";
    double *parr_dfiLBL_po_dt = new double [mVectBigMeasures.size() * dimP ];
    memset(parr_dfiLBL_po_dt, 0, mVectBigMeasures.size() * dimP * sizeof(double));
    double *arrXlbl = new double [LBLsolv.mDimX];

    memcpy(arrXlbl, marrPAnt_PSK,3 * sizeof(double));
    memcpy(&arrXlbl[3], marrSBeaconXY_GSK,2 * sizeof(double));
    double valPrirashen = -0.0001;
    const int Ncircle = mVectBigMeasures.size();
    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
 #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
 #pragma omp for // OMP (директива для распределения итераций цикла между потоками)
        for (int i =0; i < Ncircle;++i)
        {
            double *arrfi_ = new double [LBLsolv.mDimX];

      LBLsolv.calc_arrfi_and_arrdfi_po_dt(arrXlbl, i, valPrirashen,arrfi_,&parr_dfiLBL_po_dt[i * dimP]);
      delete []arrfi_;

       qDebug() << "LBLsolv.calc_arrfi_and_arrdfi_po_dt" << i;
        }
}
    //2!

// 3. Вычисление матриц dfi_po_dt
    double *parr_dfi_po_dt= new double [Ncircle * mDimX ];
    memset(parr_dfi_po_dt, 0, Ncircle * mDimX * sizeof(double));

    //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
#pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
{ // OMP (начало блока, который выполняется в нескольких потоках
 #pragma omp for // OMP (директива для распределения итераций цикла между потоками)
    for (int i =0; i < Ncircle;++i)
    {
        double *arrfi = new double [mDimX];
       calc_arrfi_and_arrdfi_po_dt(arrAngX0, i,valPrirashen
                                                         ,arrfi,&parr_dfi_po_dt[i * mDimX]);
       delete []arrfi;
      // qDebug() << "calc_arrfi_and_arrdfi_po_dt" << i;
    }
 }
    //3!

    //4. основной цикл
    qDebug() << "In Main Circle calc_AngKt";
    double *arr_Sum = new double [ mDimX * mDimX];
    memset(arr_Sum, 0, sizeof(double) *  mDimX * mDimX);
    for (int i = 0; i < Ncircle;++i)
    {
        double *parrT0 = new double [dimP  ];
        double *parrT1 = new double [mDimX ];
        double *parrT11 = new double [mDimX ];
        double *parrT2 = new double [mDimX * mDimX];
        MtrxMultMatrx( arrMtrx_dFGrLBL_po_dX_Inv,dimP, dimP, &parr_dfiLBL_po_dt[i * dimP ],1, parrT0);
        MtrxMultMatrx( arr_dF_po_dP,mDimX, dimP, parrT0,1, parrT1);
        MtrxSumMatrx(parrT1, &parr_dfi_po_dt[i * mDimX],mDimX, 1, parrT1) ;
        memcpy(parrT11, parrT1, sizeof(double) * mDimX);
        qDebug() << "1";
        MtrxMultMatrxTransp(parrT1,mDimX, 1, parrT11,mDimX, parrT2) ;

        MtrxSumMatrx(arr_Sum, parrT2,mDimX, mDimX, arr_Sum) ;
        qDebug() << "2";
        delete []parrT0;
        delete []parrT1;
        delete []parrT2;
        delete []parrT11;
    }
 qDebug() << "OUT FROM CIRCLE calc_AngKt";
    double arr_AngKt_per_1sec[9] = {0.};
    calcF_D_FTransp_(arrMtrx_dFGr_po_dXang_Inv,mDimX,mDimX,arr_Sum, arr_AngKt_per_1sec);
    QPosSolver::createDblVect(arr_AngKt_per_1sec,9,vct_AngKt_per_1sec);

    delete []arr_dF_po_dP;
    delete []parr_dfiLBL_po_dt;
    delete []arrXlbl;    
    delete []parr_dfi_po_dt;
    delete []arr_Sum;
  return true;
}



