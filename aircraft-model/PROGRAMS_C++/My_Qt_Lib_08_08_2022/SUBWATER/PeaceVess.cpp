#include "PeaceVess.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <dir.h>
#include <wchar.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "MatrixProccess.h"


#include "Gauss.h"


#include "CoordSystTrsf.h"
#include "BigMeasure.h"
#include "Table_1D.h"
#include "SubWaterBeam.h"
#include "TrueMeasParams.h"



extern const int QUantCols_VessTrajReport =  14;




//---------------------------------------------------------------------------

QPeaceVess::QPeaceVess()
{
   mPlatform = QPlatform();

   mSins = QPeaceSins () ;

   mGps = QGps();

  //	mTraceFlt = TTraceFlt() ;
	// константы
	mWidth = 0; // ширина(м)
	mLength = 0; // длина(м)

	mMaxQ = 3./180.*M_PI; // максимальный угол отклонения от генерального курса(амплитуда угла рыскания)
	mT_Q = 18.; // период рыскания
	mMaxPsi = 3./180.*M_PI;// максимальный угол килевой качки(амплитуда)
	mT_Psi = 12; // период килевой качки
	mMaxTet = 12./180.*M_PI; //максимальный угол боротовой качки(амплитуда)
    //mT_Tet = 6; // период бортовой качки
    mT_Tet = 8.; // период бортовой качки
	mMaxVert = 1. ;

	// парамеитры движения
	mQ0 =0. ; // генеральный курс
	mVVess = 0.514 * 20 ;// скорость корабля своего 20 узлов
	mTVess = 0.; // время привязки траекторной информации
	memset(marrVectSost, 0 , 9 * sizeof(double)); // вектор траектории( положения, скорости и ускорения в ГСК )
    //memset(marrVectSost0, 0 , 9 * sizeof(double));
    //memset(marrEstVectSost, 0 , 9 * sizeof(double)); // вектор траектории( положения, скорости и ускорения в ГСК )
	marrVectSost[4] = mVVess ;
	mQ = 0. ;// угол курса
	mPsi = 0. ; // угол килевой качки
	mTet = 0. ; // угол бортовой качки
	mVQ = 0. ; //скорость изменения угла курса
	mVPsi = 0. ; // скорость изменения угла килевой качки
	mVTet = 0.; // вкорость изменения угла бортовой качуки

    memset(marrDelt, 0, 4 * sizeof(double)) ;


}

//---------------------------------------------------------------------------
// конструктор копирования
 QPeaceVess ::QPeaceVess (const QPeaceVess &R)
 {	
    mPlatform = R.mPlatform;
    mWidth = R.mWidth; // ширина(м)
	mLength = R.mLength; // длина(м)

	mMaxQ = R.mMaxQ; // максимальный угол отклонения от генерального курса(амплитуда угла рыскания)
	mT_Q = R.mT_Q; // период рыскания
	mMaxPsi = R.mMaxPsi;// максимальный угол килевой качки(амплитуда)
	mT_Psi = R.mT_Psi; // период килевой качки
	mMaxTet = R.mMaxTet; //максимальный угол боротовой качки(амплитуда)
	mT_Tet = R.mT_Tet; // период бортовой качки
	mMaxVert = R.mMaxVert ;

	// парамеитры движения
	mQ0 =R.mQ0 ; // генеральный курс
	mVVess = R.mVVess ;
	mTVess = R.mTVess; // время привязки траекторной информации
    memcpy(marrVectSost,R.marrVectSost , 9 * sizeof(double)); // вектор  положения в ГСК )
    //memcpy(marrVectSost0,R.marrVectSost0 , 9 * sizeof(double)); // вектор  положения в ГСК )
//	memcpy(marrEstVectSost,R.marrEstVectSost , 9 * sizeof(double)); // вектор  положения в ГСК )

	mQ = R.mQ ;
	mPsi = R.mPsi ; // угол килевой качки
	mTet = R.mTet ; // угол бортовой качки
	mVQ = R.mVQ ; //скорость изменения угла курса
	mVPsi = R.mVPsi ; // скорость изменения угла килевой качки
	mVTet = R.mVTet; // вкорость изменения угла бортовой качуки



	memcpy(marrDelt, R.marrDelt, 4 * sizeof(double)) ;
    mSins = R.mSins;
    mGps =R.mGps;

 }
 // оператор присваивания
 QPeaceVess &QPeaceVess::operator=(const QPeaceVess  &R)
 {
     mPlatform = R.mPlatform;

     mWidth = R.mWidth; // ширина(м)
     mLength = R.mLength; // длина(м)

     mMaxQ = R.mMaxQ; // максимальный угол отклонения от генерального курса(амплитуда угла рыскания)
     mT_Q = R.mT_Q; // период рыскания
     mMaxPsi = R.mMaxPsi;// максимальный угол килевой качки(амплитуда)
     mT_Psi = R.mT_Psi; // период килевой качки
     mMaxTet = R.mMaxTet; //максимальный угол боротовой качки(амплитуда)
     mT_Tet = R.mT_Tet; // период бортовой качки
     mMaxVert = R.mMaxVert ;

     // парамеитры движения
     mQ0 =R.mQ0 ; // генеральный курс
     mVVess = R.mVVess ;
     mTVess = R.mTVess; // время привязки траекторной информации
     memcpy(marrVectSost,R.marrVectSost , 9 * sizeof(double)); // вектор  положения в ГСК ) 
    // memcpy(marrVectSost0,R.marrVectSost0 , 9 * sizeof(double)); // вектор  положения в ГСК )

     mQ = R.mQ ;
     mPsi = R.mPsi ; // угол килевой качки
     mTet = R.mTet ; // угол бортовой качки
     mVQ = R.mVQ ; //скорость изменения угла курса
     mVPsi = R.mVPsi ; // скорость изменения угла килевой качки
     mVTet = R.mVTet; // вкорость изменения угла бортовой качуки



     memcpy(marrDelt, R.marrDelt, 4 * sizeof(double)) ;

     mSins = R.mSins;
     mGps =R.mGps;

	return *this ;
 }


 // парам конструктор 1
QPeaceVess::QPeaceVess (const TEnvironment Environment,const double Width,const double Length
                  ,const  double MaxQ ,const  double T_Q
                  ,const double MaxPsi,const double T_Psi ,const  double MaxTet
                  ,const double T_Tet,const double MaxVert, const double Q0,const double VVess
                  , const QPeaceSins Sins
                  ,const  QPlatform Platform, const double TVess, const QGps Gps)
{
  mTVess = TVess;
  const	double valCoeffEnv = ((double) (Environment.mBallWave))/9.;
  mWidth  =Width ;
  mLength =Length ;


  mMaxQ =MaxQ * valCoeffEnv ;
  mT_Q  =T_Q;
  mMaxPsi =MaxPsi * valCoeffEnv ;
  mT_Psi = T_Psi ;
  mMaxTet =MaxTet * valCoeffEnv ;
  mT_Tet =T_Tet ;
  mMaxVert = MaxVert * valCoeffEnv ;



  mQ0 =Q0 ;
  mVVess = VVess ;

  //for (int i = 0 ; i < 4; i++)
 // {
  // marrDelt [i] = -M_PI/2. + getRand01() * 2. * M_PI;
  //}
  marrDelt [0] = -M_PI/4. + getRand01() * 2. * M_PI;
  marrDelt [1] = -M_PI/4. + getRand01() * 2. * M_PI;
  marrDelt [2] = -M_PI/4. + getRand01() * 2. * M_PI;
  marrDelt [3] = -M_PI/2. + getRand01() * 2. * M_PI;

  memset(marrVectSost,0, 9 * sizeof(double)) ;
  mQ = mQ0 +  mMaxQ * cos(marrDelt[0]);
  mVQ =   2. * M_PI/mT_Q * mMaxQ * sin(marrDelt[0]);

  mPsi = mMaxPsi * cos(marrDelt[1]);
  mVPsi =  2. * M_PI/mT_Psi * mMaxPsi * sin( marrDelt[1]);

  mTet =  mMaxTet * cos( marrDelt[2]);
  mVTet =  2. * M_PI/mT_Tet * mMaxTet * sin( marrDelt[2]);


  // вектор состояния корабля в ГСК истинный
   marrVectSost[0] = 0;
   marrVectSost[1] = 0;
   marrVectSost[2] = mMaxVert * sin( marrDelt[3]) ;
   marrVectSost[3] = mVVess * sin(mQ) ;
   marrVectSost[4] = mVVess * cos(mQ) ;
   marrVectSost[5] = 2. * M_PI/ mT_Psi * mMaxVert * cos( marrDelt[3]) ;
   mSins = Sins;
   tuneCurrentSinsMeasures();
   mPlatform = Platform;
   mGps = Gps;

}
// парам конструктор 1
QPeaceVess::QPeaceVess (const TEnvironment Environment,const double Width,const double Length
                 ,const  double MaxQ ,const  double T_Q
                 ,const double MaxPsi,const double T_Psi ,const  double MaxTet
                 ,const double T_Tet,const double MaxVert, const double Q0,const double VVess
                 , const QPeaceSins Sins,const  QPlatform Platform, const double TVess
                 , const QGps Gps, const double *arrS_Vess0)
{
 mTVess = TVess;
 const	double valCoeffEnv = ((double) (Environment.mBallWave))/9.;
 mWidth  =Width ;
 mLength =Length ;


 mMaxQ =MaxQ * valCoeffEnv ;
 mT_Q  =T_Q;
 mMaxPsi =MaxPsi * valCoeffEnv ;
 mT_Psi = T_Psi ;
 mMaxTet =MaxTet * valCoeffEnv ;
 mT_Tet =T_Tet ;
 mMaxVert = MaxVert * valCoeffEnv ;



 mQ0 =Q0 ;
 mVVess = VVess ;

 //for (int i = 0 ; i < 4; i++)
// {
 // marrDelt [i] = -M_PI/2. + getRand01() * 2. * M_PI;
 //}
 marrDelt [0] = -M_PI/4. + getRand01() * 2. * M_PI;
 marrDelt [1] = -M_PI/4. + getRand01() * 2. * M_PI;
 marrDelt [2] = -M_PI/4. + getRand01() * 2. * M_PI;
 marrDelt [3] = -M_PI/2. + getRand01() * 2. * M_PI;

 memset(marrVectSost,0, 9 * sizeof(double)) ;
 mQ = mQ0 +  mMaxQ * cos(marrDelt[0]);
 mVQ =   2. * M_PI/mT_Q * mMaxQ * sin(marrDelt[0]);

 mPsi = mMaxPsi * cos(marrDelt[1]);
 mVPsi =  2. * M_PI/mT_Psi * mMaxPsi * sin( marrDelt[1]);

 mTet =  mMaxTet * cos( marrDelt[2]);
 mVTet =  2. * M_PI/mT_Tet * mMaxTet * sin( marrDelt[2]);


 // вектор состояния корабля в ГСК истинный
  marrVectSost[0] = arrS_Vess0[0];
  marrVectSost[1] = arrS_Vess0[1];
  marrVectSost[2] = mMaxVert * sin( marrDelt[3]) ;
  marrVectSost[3] = mVVess * sin(mQ) ;
  marrVectSost[4] = mVVess * cos(mQ) ;
  marrVectSost[5] = 2. * M_PI/ mT_Psi * mMaxVert * cos( marrDelt[3]) ;
  //memcpy(marrVectSost0, marrVectSost, 9 * sizeof(double));
  mSins = Sins;
  tuneCurrentSinsMeasures();
  mPlatform = Platform;
  mGps = Gps;

}
//-----------------------------
// начальное заполнение mSins в соответствии с параметрами траектории судна
void QPeaceVess::tuneCurrentSinsMeasures()
{
mSins.fillValues_Delts_and_Ests
 (mTVess,  mVVess,	 mQ, mPsi,	 mTet, mVQ
 ,mVPsi,  mVTet, marrVectSost[2], marrVectSost[5]
 , mT_Q,  mT_Psi, mT_Tet, marrDelt) ;

}



// линейная экстраполяция вектора состояния на время  VAlTExtr впероед
 void QPeaceVess::extrapolateTrueVS_GSK(const double VAlTExtr, double *arrVessExtrapVS_GSK)
 {
     memcpy(arrVessExtrapVS_GSK, marrVectSost, 9 * sizeof(double));
     double arrT[3] = {0.};
     MatrxMultScalar(&marrVectSost[3], 1, 3, VAlTExtr, arrT);
     MtrxSumMatrx(arrT, marrVectSost,3, 1, arrVessExtrapVS_GSK) ;
 }

//------------------------------------------
 bool QPeaceVess::recalcVess(const double valT)
 {
     const double h = valT - mTVess ;
     if (h < 0)
     {
        return false;
     }
     calcCentreDeckAngles(valT) ;
    // пересчет вектора состояния корабля в ГСК истинного
     marrVectSost[0] += marrVectSost[3]*h;
     marrVectSost[1] += marrVectSost[4]*h;
     marrVectSost[2] = marrVectSost[2];//mMaxVert * sin(2. * M_PI / mT_Psi*valT + marrDelt[3]) ;
     marrVectSost[3] = mVVess * sin(mQ0) ;
     marrVectSost[4] = mVVess * cos(mQ0) ;
     marrVectSost[5] = 0.;//2. * M_PI/ mT_Psi * mMaxVert * cos(2. * M_PI/ mT_Psi*valT + marrDelt[3]) ;
     mTVess = valT ;

 }
 // вычисление углов ориентации палубы в центре качания в момент valT
 void QPeaceVess::calcCentreDeckAngles(const double valT)
 {
   mQ = mQ0 +  mMaxQ * cos(2. * M_PI/mT_Q * valT -  marrDelt[0]);
   mVQ =  -2. * M_PI/mT_Q * mMaxQ * sin(2. * M_PI/mT_Q * valT -  marrDelt[0]);

   mPsi = mMaxPsi * cos(2. * M_PI/mT_Psi * valT -  marrDelt[1]); // килевая
   mVPsi =  -2. * M_PI/mT_Psi * mMaxPsi * sin(2. * M_PI/mT_Psi * valT -  marrDelt[1]);

   mTet =  mMaxTet * cos(2. * M_PI/mT_Tet * valT -  marrDelt[2]); // бортовая
   mVTet =  -2. * M_PI/mT_Tet * mMaxTet * sin(2. * M_PI/mT_Tet * valT -  marrDelt[2]);

 }
 // вычисление углов ориентации палубы в точке в момент valT
 void QPeaceVess::calcDeckAngles(const double valT, double *arrPointPositionPSK
    ,double *pvalQ, double *pvalVQ,double *pvalPsi,double *pvalVPsi
    ,double *pvalTet ,double *pvalVTet)
 {
    calcCentreDeckAngles(valT);
     *pvalQ = mQ ;
     *pvalVQ =mVQ ;
     //

     *pvalPsi = mPsi ;
     *pvalVPsi = mVPsi ;

     *pvalTet = mTet ;
     *pvalVTet = mVTet ;

 }
 //------------------------------------------------------------------
 // передвижение корабля на время  valT
  void QPeaceVess::Move(const double valT)
 {



   recalcVess(valT ) ;
   mSins.recalcPeaceSins_v0(valT
                              ,mVVess //VVess
                              ,mQ //	const double Q
                              ,mPsi// const double Psi
                              ,mTet// Tet
                              ,mVQ //  / VQ
                              , mVPsi //   VPsi
                              ,mVTet//  VTet
                              ,marrVectSost[2] //    H
                              ,marrVectSost[5] // VH
                              ,mT_Q// c  T_Q
                              ,mT_Psi //   T_Psi
                              ,mT_Tet //   T_Tet
                              ,marrDelt// double *arrDelt
                                                     ) ;
 }


 //---------------------------------------------------
 // передвижение корабля на время  valT
  /*
  void QPeaceVess::calcPosition(TEnvironment Environment,const double valT, double *arrEilers)
 {
      /*
      // время прохода одной стороны квадрата
      double valTSide = mVVess/VALLenghtQdrtSide;

      int numside = (valT )/valTSide;
      double tcur = valT - ((double) numside)*valTSide;

      int n0 = numside -(numside/ 4) * 4;
      mQ0 = M_PI/2. + ((double)n0) *M_PI/2.;
      double arrVectSost0[6] ={0.};
      double valSignx = ((n0== 0)||(n0== 2))?-1.: 1.;
      double valSigny = ((n0== 0)||(n0== 1))? 1.:-1.;
      arrVectSost0[0] = marrVectSost0[0] + valSignx *  VALLenghtQdrtSide/ 2.;
      arrVectSost0[1] = marrVectSost0[1] + valSigny *  VALLenghtQdrtSide/ 2.;
      arrVectSost0[2] = marrVectSost0[2];
      arrVectSost0[3] = mVVess * sin(mQ0) ;
      arrVectSost0[4] = mVVess * cos(mQ0) ;
      arrVectSost0[5] = 0.;

      marrVectSost[0] = arrVectSost0[0] + arrVectSost0[3] * tcur;
      marrVectSost[1] = arrVectSost0[1] + arrVectSost0[4] * tcur;
      marrVectSost[2] = arrVectSost0[2];
      marrVectSost[3] = arrVectSost0[3];
      marrVectSost[4] = arrVectSost0[4];
      marrVectSost[5] = arrVectSost0[5];

   mSins.recalcPeaceSins_v0(valT
                              ,mVVess //VVess
                              ,mQ //	const double Q
                              ,mPsi// const double Psi
                              ,mTet// Tet
                              ,mVQ //  / VQ
                              , mVPsi //   VPsi
                              ,mVTet//  VTet
                              ,marrVectSost[2] //    H
                              ,marrVectSost[5] // VH
                              ,mT_Q// c  T_Q
                              ,mT_Psi //   T_Psi
                              ,mT_Tet //   T_Tet
                              ,marrDelt// double *arrDelt
                                                     ) ;



    integrateEilersVect (arrEilers);

 }
 */
  //-------------------------------------------------
  //--------------------------------------

 // перечсчет вектора соcтояния цели  из KGSK в ПСК c центром в центре гаджета кажущуюся
 // если  lenarrKGSK == 6 , то пересчитывается положение и скорость
 // если   lenarrKGSK == 3   , то пересчитывается только положение
 // на вход подается вектор состоящий из 3 или 6  координат.
 // первые 3 координаты представляют из себя  положение точки в КГС
 // последние 3 координаты представляют из себя скорость точки в КГСК
 // INPUT:
 //arrKGSK[lenarrKGSK] - вектор положния (положения и скорости) цели
 //arrGdgPosParams[3]  - вектор положение гаджета в ПСК
 // arrEilerCntrKP[3] -вектора углов ориентации корабля (Q, Tetta, Psi)
 // arr_dEilers0_po_dt[3] - вектор  скоростей углов (Psi, Tet, Q)
 // OUTPUT:
 // arrPSK- вектор положения (или положениея и скорости) цели в ПСК -гаджет
 // вычисления производятся в соответствие с правилами теор механики
 void  QPeaceVess::RecalcVect_KGSK_INTO_GdgPSK( double *arrKGSK
                    , double *arrEilers0, double *arr_dEilers0_po_dt,  double *arrGdgPosParams
                     ,double *arrPSK,const int lenarrKGSK )
{

    // персчет вектора положения в АСК
    // создание матрицы перехода из   КГСК в ПСК
    double arrEilers[3] = {0.}, arr_dEilers_po_dt[3] = {0.};
    MatrxMultScalar(arrEilers0, 1, 3, -1.,arrEilers);
    double matrPereh_PSK_V_KGSK[9] = {0} ;
    QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilers, matrPereh_PSK_V_KGSK) ;
    // вектор положения в ПСК-центр тяжести
    // вычисление вектора положениея в ПСК
    double arrPosPSK_CT[3] = {0.};
    MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrKGSK, 1, arrPosPSK_CT) ;
    // создание вектора положения в ПСК гаджета
    MtrxMinusMatrx(arrPosPSK_CT, arrGdgPosParams,1, 3, arrPSK);
    ///

    // пересчет вектора скорости
    if (lenarrKGSK == 6)
    {
        MatrxMultScalar(arr_dEilers0_po_dt, 1, 3, -1.,arr_dEilers_po_dt);
        double arrOm[3] = {0.};
        //calcVectOm(arrEilers, arr_dEilers_po_dt, arrOm);
        MtrxTranspMultMatrx( matrPereh_PSK_V_KGSK,3, 3, &arrKGSK[3],1, &arrPSK[3]) ;


        double arrV_Per[3];
        OuterProduct(arrOm,arrPosPSK_CT ,  arrV_Per);
        MtrxMinusMatrx(&arrPSK[3], arrV_Per,1, 3, &arrPSK[3]);
        ///

    }
    return ;

}
 //---------------------------------------
 // 26.10.2021
 // пересчет вектора соcтояния цели  из KGSK в ПСК_CT
 // если  lenarrKGSK == 6 , то пересчитывается положение и скорость
 // если   lenarrKGSK == 3   , то пересчитывается только положение
 // на вход подается вектор состоящий из 3 или 6  координат.
 // первые 3 координаты представляют из себя  положение точки в КГС
 // последние 3 координаты представляют из себя скорость точки в КГСК
 // INPUT:
 //arrKGSK[lenarrKGSK] - вектор положния (положения и скорости) цели
 //arrGdgPosParams[3]  - вектор положение гаджета в ПСК
 // arrEilerCntrKP[3] -вектора углов ориентации корабля (Q, Tetta, Psi)
 // arrOmegaPSK[3] - вектор угловых скоростей (Psi, Tet, Q)
 // OUTPUT:
 // arrPSK- вектор положения (или положениея и скорости) цели в ПСК_CT
 // вычисления производятся в соответствие с правилами теор механики
 void  QPeaceVess::RecalcVect_KGSK_INTO_PSK_CT( double *arrKGSK, double *arrEilers0
             , double *arr_dEilers0_po_dt,double *arrPSK_CT,const int lenarrKGSK )
{


     // персчет вектора положения в АСК
     // создание матрицы перехода из   КГСК в ПСК
     double arrEilers[3] = {0.}, arr_dEilers_po_dt[3] = {0.};
     MatrxMultScalar(arrEilers0, 1, 3, -1.,arrEilers);
     double matrPereh_PSK_V_KGSK[9] = {0} ;
     QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilers, matrPereh_PSK_V_KGSK) ;
     // вектор положения в ПСК-центр тяжести
     // вычисление вектора положениея в ПСК
    // double arrPosPSK_CT[3] = {0.};
     MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrKGSK, 1, arrPSK_CT) ;
     // создание вектора положения в ПСК гаджета
    // MtrxMinusMatrx(arrPosPSK_CT, arrGdgPosParams,1, 3, arrPSK);
     ///

     // пересчет вектора скорости
     if (lenarrKGSK == 6)
     {
         MatrxMultScalar(arr_dEilers0_po_dt, 1, 3, -1.,arr_dEilers_po_dt);
         double arrOm[3] = {0.};
         //calcVectOm(arrEilers, arr_dEilers_po_dt, arrOm);
         MtrxTranspMultMatrx( matrPereh_PSK_V_KGSK,3, 3, &arrKGSK[3],1, &arrPSK_CT[3]) ;


         double arrV_Per[3];
         OuterProduct(arrOm,arrPSK_CT ,  arrV_Per);
         MtrxMinusMatrx(&arrPSK_CT[3], arrV_Per,1, 3, &arrPSK_CT[3]);
         ///

     }
     return ;

}
 //---------------------------------------
 // 26.10.2021
 // пересчет вектора соcтояния цели  из KGSK в ПСК_CT
 // если  lenarrKGSK == 6 , то пересчитывается положение и скорость
 // если   lenarrKGSK == 3   , то пересчитывается только положение
 // на вход подается вектор состоящий из 3 или 6  координат.
 // первые 3 координаты представляют из себя  положение точки в КГС
 // последние 3 координаты представляют из себя скорость точки в КГСК
 // INPUT:
 //arrKGSK[lenarrKGSK] - вектор положния (положения и скорости) цели
 //arrGdgPosParams[3]  - вектор положение гаджета в ПСК
 // arrEilerCntrKP[3] -вектора углов ориентации корабля (Q, Tetta, Psi)
 // arrOmegaPSK[3] - вектор угловых скоростей (Psi, Tet, Q)
 // OUTPUT:
 // arrPSK- вектор положения (или положениея и скорости) цели в ПСК_CT
 // вычисления производятся в соответствие с правилами теор механики
 void  QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( double *arr_PSK_CT, double *arrEilers0
             , double *arr_dEilers0_po_dt,double *arr_KGSK,const int lenarrKGSK )
{

    // персчет вектора положения в АСК
    // создание матрицы перехода из   КГСК в ПСК
    double arrEilers[3] = {0.}, arr_dEilers_po_dt[3] = {0.};
    MatrxMultScalar(arrEilers0, 1, 3, -1.,arrEilers);
    double matrPereh_PSK_V_KGSK[9] = {0} ;
    QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilers, matrPereh_PSK_V_KGSK) ;
    // вектор положения в ПСК-центр тяжести
    // вычисление вектора положениея в ПСК

    MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arr_PSK_CT, 1, arr_KGSK) ;

    ///

    // пересчет вектора скорости
    if (lenarrKGSK == 6)
    {
        MatrxMultScalar(arr_dEilers0_po_dt, 1, 3, -1.,arr_dEilers_po_dt);
        double arrOm[3] = {0.};
      //  calcVectOm(arrEilers, arr_dEilers_po_dt, arrOm);



        double arrV_Per[3] = {0.}, arr_temp[3] = {0.};
        OuterProduct(arrOm,arr_PSK_CT ,  arrV_Per);
        MtrxSumMatrx(&arr_PSK_CT[3], arrV_Per,1, 3, &arr_temp[3]);
        MtrxMultMatrx( matrPereh_PSK_V_KGSK,3, 3, &arr_temp[3],1, &arr_KGSK[3]) ;
        ///

    }
    return ;  

}

 //-----------------------------------
 //---------------------------------------
 // 26.10.2021
 // пересчет вектора соcтояния цели  из KGSK в ПСК_CT
 // если  lenarrKGSK == 6 , то пересчитывается положение и скорость
 // если   lenarrKGSK == 3   , то пересчитывается только положение
 // на вход подается вектор состоящий из 3 или 6  координат.
 // первые 3 координаты представляют из себя  положение точки в КГС
 // последние 3 координаты представляют из себя скорость точки в КГСК
 // INPUT:
 //arrKGSK[lenarrKGSK] - вектор положния (положения и скорости) цели
 //arrGdgPosParams[3]  - вектор положение гаджета в ПСК
 // arrEilerCntrKP[3] -вектора углов ориентации корабля (Q, Tetta, Psi)
 // arrOmegaPSK[3] - вектор угловых скоростей (Psi, Tet, Q)
 // OUTPUT:
 // arrPSK- вектор положения (или положениея и скорости) цели в ПСК_CT
 // вычисления производятся в соответствие с правилами теор механики
 void  QPeaceVess::RecalcVect_PSK_CT_INTO_GSK( double *arr_PSK_CT, double *arrEilers0
             , double *arr_dEilers0_po_dt,double *arr_GSK,const int lenarrKGSK )
{

    // персчет вектора положения в АСК
    // создание матрицы перехода из   КГСК в ПСК
     double arrEilers[3] = {0.}, arr_dEilers_po_dt[3] = {0.}, arr_KGSK[3] = {0.};
    MatrxMultScalar(arrEilers0, 1, 3, -1.,arrEilers);
    double matrPereh_PSK_V_KGSK[9] = {0} ;
    QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilers, matrPereh_PSK_V_KGSK) ;
    // вектор положения в ПСК-центр тяжести
    // вычисление вектора положениея в ПСК

    MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arr_PSK_CT, 1, arr_KGSK) ;

    MtrxSumMatrx(arr_KGSK, marrVectSost,1, 3, arr_GSK) ;

    ///

    // пересчет вектора скорости
    if (lenarrKGSK == 6)
    {
        MatrxMultScalar(arr_dEilers0_po_dt, 1, 3, -1.,arr_dEilers_po_dt);
        double arrOm[3] = {0.}, arrV_KGSK[3] = {0.};
      //  calcVectOm(arrEilers, arr_dEilers_po_dt, arrOm);



        double arrV_Per[3] = {0.}, arr_temp[3] = {0.};
        OuterProduct(arrOm,arr_PSK_CT ,  arrV_Per);
        MtrxSumMatrx(&arr_PSK_CT[3], arrV_Per,1, 3, &arr_temp[3]);
        MtrxMultMatrx( matrPereh_PSK_V_KGSK,3, 3, &arr_temp[3],1, &arrV_KGSK[3]) ;
        MtrxSumMatrx(arrV_KGSK, &(marrVectSost[3]),1, 3, &(arr_GSK[3])) ;
        ///

    }
    return ;

}
 // 26.10.2021
 // перечсчет вектора соcтояния цели  из ПSK в ПСК c центром в центре гаджета кажущуюся
 // если  lenarrKGSK == 6 , то пересчитывается положение и скорость
 // если   lenarrKGSK == 3   , то пересчитывается только положение
 // на вход подается вектор состоящий из 3 или 6  координат.
 // первые 3 координаты представляют из себя  положение точки в ПГС
 // последние 3 координаты представляют из себя скорость точки в ПСК
 // INPUT:
 //arrPSK[lenarrKGSK] - вектор положния (положения и скорости) цели
 //arrGdgPosParams[3]  - вектор положение гаджета в ПСК
 // arrEilerCntrKP[3] -вектора углов ориентации корабля (Q, Tetta, Psi)
 // arrOmegaPSK[3] - вектор угловых скоростей (Psi, Tet, Q)
 // OUTPUT:
 // arrGdgPSK- вектор положения (или положениея и скорости) цели в ПСК -гаджет
 // вычисления производятся в соответствие с правилами теор механики
 void  QPeaceVess::RecalcVect_PSK_INTO_GdgPSK( double *arrPSK_CT, double *arrGdgPosParams
                     ,double *arrGdgPSK,const int lenarrKGSK )
{
    // создание вектора положения в ПСК гаджета
    MtrxMinusMatrx(arrPSK_CT, arrGdgPosParams,1, 3, arrGdgPSK);
    ///

    // пересчет вектора скорости
    if (lenarrKGSK == 6)
    {
        memcpy(&arrGdgPSK[3], &arrPSK_CT[3], 3 *sizeof(double));
    }
    return ;

}

 //------------------------------------------------------
  // перечсчет вектора соcтояния из KGSK в ПСК c центром в центре гаджета кажущуюся
  // если  lenarrKGSK == 6 , то пересчитывается положение и скорость
  // если   lenarrKGSK == 3   , то пересчитывается только положение
  // на вход подается вектор состоящий из 3 или 6  координат.
  // первые 3 координаты представляют из себя  положение точки в КГС
  // последние 3 координаты представляют из себя скорость точки в КГСК
  // INPUT:
  //arrKGSK[lenarrKGSK] - вектор положния (положения и скорости) цели
  //arrGdgPosParams[3]  - вектор положение гаджета в ПСК
  // arrEilerCntrKP[3] -вектора углов ориентации корабля (Q, Tetta, Psi)
  // arrOmegaPSK[3] - вектор угловых скоростей (Psi, Tet, Q)
  // OUTPUT:
  // arrPSK- вектор положения (или положениея и скорости) цели в ПСК -гаджет
  // вычисления производятся при помощи непосредственного дифференцирования
void  QPeaceVess::RecalcVect_KGSK_INTO_GdgPSK_differentiation( double *arrKGSK
                     , double *arrEilerCntrKP0, double *arr_dEilers_po_dt0,  double *arrGdgPosParams
                      ,double *arrPSK,const int lenarrKGSK )
{
    // персчет вектора положения в АСК
    // создание матрицы перехода из   КГСК в ПСК
    double arrEilerCntrKP[3] = {0.}, arr_dEilers_po_dt[3] = {0.};
    MatrxMultScalar(arrEilerCntrKP0, 3, 1, -1,arrEilerCntrKP);

    double matrPereh_PSK_V_KGSK[9] = {0} ;
    QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilerCntrKP, matrPereh_PSK_V_KGSK) ;
    // вектор положения в ПСК-центр тяжести
    // вычисление вектора положениея в ПСК
    double arrPosPSK_CT[3] = {0.};
    MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrKGSK, 1, arrPosPSK_CT) ;
    // создание вектора положения в ПСК гаджета
    MtrxMinusMatrx(arrPosPSK_CT, arrGdgPosParams,1, 3, arrPSK);
    ///

    // пересчет вектора скорости
    if (lenarrKGSK == 6)
    {
        MatrxMultScalar(arr_dEilers_po_dt0, 3, 1, -1,arr_dEilers_po_dt);
        MtrxTranspMultMatrx( matrPereh_PSK_V_KGSK,3, 3, &arrKGSK[3],1, &arrPSK[3]) ;
        double arr_dM_po_dQ[9] = {0.},arr_dM_po_dPsi[9]= {0.}, arr_dM_po_dTet[9]= {0.};
        QCoordSystTrsf::calc_dMLeft_psk_to_kgsk_po_dQ(arrEilerCntrKP,arr_dM_po_dQ);

        QCoordSystTrsf::calc_dMLeft_psk_to_kgsk_po_dPsi(arrEilerCntrKP,arr_dM_po_dPsi);

        QCoordSystTrsf::calc_dMLeft_psk_to_kgsk_po_dTet(arrEilerCntrKP,arr_dM_po_dTet);

        MatrxMultScalar(arr_dM_po_dQ, 1, 9, arr_dEilers_po_dt[2],arr_dM_po_dQ);
        MatrxMultScalar(arr_dM_po_dPsi, 1, 9, arr_dEilers_po_dt[0],arr_dM_po_dPsi);
        MatrxMultScalar(arr_dM_po_dTet, 1, 9, arr_dEilers_po_dt[1],arr_dM_po_dTet);

        double arrT0[9] = {0.},arrT1[9] = {0.},arrT2[3] = {0.};
        MtrxSumMatrx(arr_dM_po_dQ, arr_dM_po_dPsi,1, 9, arrT0) ;
        MtrxSumMatrx(arrT0, arr_dM_po_dTet,1, 9, arrT1) ;

        MtrxTranspMultMatrx(arrT1,3, 3,  arrKGSK,1, arrT2);
        MtrxSumMatrx(&arrPSK[3], arrT2,1, 3, &arrPSK[3]) ;
        ///
    }
    return ;

}
 //-------------------------------------------------
// пересчет вектора соcтояния из KGSK в собственную систему
// если  lenarrKGSK == 6 , то пересчитывается положение и скорость
// если   lenarrKGSK == 3   , то пересчитывается только положение
// на вход подается вектор состоящий из 3 или 6  координат.
// первые 3 координаты представляют из себя  положение точки в КГС
// последние 3 координаты представляют из себя скорость точки в КГСК
// INPUT:
//arrKGSK[lenarrKGSK] - вектор положния (положения и скорости) цели
//arrGdgPosParams[6]  - вектор положение гаджета в ПСК, вектор углов ориентации гаджета в ПСК
// arrEilerCntrKP[3] -вектора углов ориентации корабля (Q, Tetta, Psi)
// arrOmegaPSK0[3] - вектор угловых скоростей (Psi, Tet, Q)
// OUTPUT:
// arrSobSK- вектор положения (или положениея и скорости) цели в ПСК -гаджет
void QPeaceVess::recalcVect_KGSK_INTO_GdgSobSK( double *arrKGSK
           , double *arrEilerCntrKP, double *arrVectOmegaPSK0
           ,  double *arrGdgPosParams,double *arrSobSK,const int lenarrKGSK )
{

// пересчет вектора arrKGSK в палубную (подвижную) систему координат, привязанную к центру гаджета
    double arrPSK[6] = {0.};
    RecalcVect_KGSK_INTO_GdgPSK(arrKGSK,arrEilerCntrKP,arrVectOmegaPSK0, arrGdgPosParams
                            ,arrPSK,lenarrKGSK );
    ///

// пересчет вектора arrPSK из ПСК в собственную СК гаджета

    double matrPereh_ASK_V_PSK[9]= {0.};
    QCoordSystTrsf::calcMtrx3_ASPK_v_PSK(&(arrGdgPosParams[3]),matrPereh_ASK_V_PSK) ;

    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrPSK,1, arrSobSK) ;
    if (lenarrKGSK == 6)
    {
       MtrxTranspMultMatrx(matrPereh_ASK_V_PSK,3, 3, &arrPSK[3],1, &arrSobSK[3]) ;
    }
    return;
}

//------------------------------------------------
// 26.10.2021
// пересчет вектора соcтояния из PSK_CT в собственную систему
// если  lenarrKGSK == 6 , то пересчитывается положение и скорость
// если   lenarrKGSK == 3   , то пересчитывается только положение
// на вход подается вектор состоящий из 3 или 6  координат.
// первые 3 координаты представляют из себя  положение точки в КГС
// последние 3 координаты представляют из себя скорость точки в КГСК
// INPUT:
//arrKGSK[lenarrKGSK] - вектор положния (положения и скорости) цели
//arrGdgPosParams[6]  - вектор положение гаджета в ПСК, вектор углов ориентации гаджета в ПСК

// OUTPUT:
// arrSobSK- вектор положения (или положениея и скорости) цели в ПСК -гаджет
void QPeaceVess::recalcVect_PSK_INTO_GdgSobSK( double *arrPSK_CT
           ,  double *arrGdgPosParams,double *arrSobSK,const int lenarrKGSK )
{

// пересчет вектора arrKGSK в палубную (подвижную) систему координат, привязанную к центру гаджета
    double arrPSK[6] = {0.};
    RecalcVect_PSK_INTO_GdgPSK(arrPSK_CT, arrGdgPosParams
                            ,arrPSK,lenarrKGSK );
    ///

// пересчет вектора arrPSK из ПСК в собственную СК гаджета

    double matrPereh_ASK_V_PSK[9]= {0.};
    QCoordSystTrsf::calcMtrx3_ASPK_v_PSK(&(arrGdgPosParams[3]),matrPereh_ASK_V_PSK) ;

    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrPSK,1, arrSobSK) ;
    if (lenarrKGSK == 6)
    {
       MtrxTranspMultMatrx(matrPereh_ASK_V_PSK,3, 3, &arrPSK[3],1, &arrSobSK[3]) ;
    }
    return;
}

//------------------------------
//
// возвращает true, если сигнал проходит
bool QPeaceVess::createTrueMeasureParams( TTable_1D &tblPrfl,const double TObrabotki
                                , double *arrTrueBeaconPos,QTrueMeasParams &trueMeasParams  )
{
  // 1.вычисление момента прихода ответного сигнала
    double valOtvetT = -1.;
  if (!calcOtvetT(tblPrfl, TObrabotki, arrTrueBeaconPos, valOtvetT))
  {

     return false;
  }
  // 1!


  // 2. формирование истинного измерения по запросному сигналу
  integrateEilersVect (trueMeasParams.marrMuWave);
  MtrxSumMatrx(trueMeasParams.marrMuWave, mSins.marrSinsSyst,1, 3, trueMeasParams.marrMuWave) ;
  trueMeasParams.mTzapr = mTVess;
  memcpy(trueMeasParams.marrSVessWave,marrVectSost, 3 * sizeof(double) );
  mSins.getEstArrEilers(mTVess,trueMeasParams.marrMuWaveZv);
  // 2!

  // 3. передвижение корабля
  Move(valOtvetT);
  // 3!


  // 4. вычисление вектора истинного положения антенны на момент
    // прихода ответного сигнала
    double arrY[6] = {0.}, arrEilers0[3] = {0.};
    integrateEilersVect (arrEilers0);
    RecalcVect_PSK_CT_INTO_GSK( mPlatform.marrPosParams, arrEilers0, NULL,arrY,3 );
    //4!

    // 5. вычисление угла скольжения valTetta и
    // курсового угла направления на маяк
    double valTetta =0., val_q_GSK = 0., val_e_GSK = 0.;
    calcTetta_InpGSK(tblPrfl, arrY,arrTrueBeaconPos,valTetta);
    val_q_GSK =calcCourseAngle(arrTrueBeaconPos[0] - arrY[0], arrTrueBeaconPos[1] - arrY[1]);
    val_e_GSK =  - valTetta;
    // 5!

    // 6. формирование единичного вектора направления
    // принимаемомго сигнала в ГСК
    double arr_gsk[3] = {0.} ;

    arr_gsk[0] = cos(val_e_GSK) * sin(val_q_GSK);
    arr_gsk[1] = cos(val_e_GSK) * cos(val_q_GSK);
    arr_gsk[2] = sin(val_e_GSK) ;
    // 6!

    // 7. пересчет вектора arr_gsk в АСПК
     double arr_aspk[3] = {0.}, arrPosParamsTemp[6] ={0.};
    memcpy(&arrPosParamsTemp[3], &(mPlatform.marrPosParams[3]), 3 * sizeof(double));
    recalcVect_KGSK_INTO_GdgSobSK( arr_gsk
               , arrEilers0, NULL
               ,  arrPosParamsTemp, arr_aspk, 3 );

    // 7!

    // 8. вычисление истинных угловых измерений
    trueMeasParams.me = asin(arr_aspk[2]);

    //trueMeasParams.mq = atan2(arr_aspk[0], arr_aspk[1]);
    trueMeasParams.mq = calcCourseAngle(arr_aspk[0], arr_aspk[1]);
    // 8!

    // 9.завершение заполнения структуры trueMeasParams
    integrateEilersVect (trueMeasParams.marrMu);
    MtrxSumMatrx(trueMeasParams.marrMu, mSins.marrSinsSyst,1, 3, trueMeasParams.marrMu) ;
    trueMeasParams.mTotv = mTVess;
    memcpy(trueMeasParams.marrSVess,marrVectSost, 3 * sizeof(double) );
    mSins.getEstArrEilers(mTVess,trueMeasParams.marrMuZv);

    // ОТЛАДКА
    double pval_q =0., pval_e=0.,  pval_t=0.;
    transf_GSK_XYZ_to_USBL3D(tblPrfl, arrTrueBeaconPos, trueMeasParams.marrSVess, trueMeasParams.marrMuZv
      ,mPlatform.marrPosParams,  &pval_q, &pval_e,  &pval_t);
    trueMeasParams.me =pval_e;
    trueMeasParams.mq = pval_q;
    // !ОТЛАДКА
    return true;
}
//------------------------------
//
// возвращает true, если сигнал проходит
bool QPeaceVess::createTrueMeasureParams( TTable_1D &tblPrfl,const double TObrabotki
         , double *arrTrueBeaconPos, double *arrTrueBeaconVelo,QTrueMeasParams &trueMeasParams  )
{
  // 1.вычисление момента прихода ответного сигнала
    double valOtvetT = -1.;

    double val_tOtv = -1., val_tZapr = -1., arrBeaconPosOut[3] = {0.}, arrBeaconPosInp[3] = {0.};
  if (!calcOtvetT_(tblPrfl, TObrabotki, arrTrueBeaconPos, arrTrueBeaconVelo,valOtvetT
                   ,val_tOtv, val_tZapr,arrBeaconPosInp, arrBeaconPosOut))
  {

     return false;
  }
  // 1!


  // 2. формирование истинного измерения по запросному сигналу
  integrateEilersVect (trueMeasParams.marrMuWave);
  trueMeasParams.mTzapr = mTVess;
  memcpy(trueMeasParams.marrSVessWave,marrVectSost, 3 * sizeof(double) );
  mSins.getEstArrEilers(mTVess,trueMeasParams.marrMuWaveZv);
  // 2!

  // 3. передвижение корабля
  Move(valOtvetT);
  // 3!


  // 4. вычисление вектора истинного положения антенны на момент
    // прихода ответного сигнала
    double arrY[6] = {0.}, arrEilers0[3] = {0.};
    integrateEilersVect (arrEilers0);
    RecalcVect_PSK_CT_INTO_GSK( mPlatform.marrPosParams, arrEilers0, NULL,arrY,3 );
    //4!

    // 5. вычисление угла скольжения valTetta и
    // курсового угла направления на маяк
    double valTetta =0., val_q_GSK = 0., val_e_GSK = 0.;
    calcTetta_InpGSK(tblPrfl, arrY,arrBeaconPosOut,valTetta);
    val_q_GSK =calcCourseAngle(arrBeaconPosOut[0] - arrY[0], arrBeaconPosOut[1] - arrY[1]);
    val_e_GSK =  - valTetta;
    // 5!

    // 6. формирование единичного вектора направления
    // принимаемомго сигнала в ГСК
    double arr_gsk[3] = {0.} ;

    arr_gsk[0] = cos(val_e_GSK) * sin(val_q_GSK);
    arr_gsk[1] = cos(val_e_GSK) * cos(val_q_GSK);
    arr_gsk[2] = sin(val_e_GSK) ;
    // 6!

    // 7. пересчет вектора arr_gsk в АСПК
     double arr_aspk[3] = {0.}, arrPosParamsTemp[6] ={0.};
    memcpy(&arrPosParamsTemp[3], &(mPlatform.marrPosParams[3]), 3 * sizeof(double));
    recalcVect_KGSK_INTO_GdgSobSK( arr_gsk
               , arrEilers0, NULL
               ,  arrPosParamsTemp, arr_aspk, 3 );
    // 7!

    // 8. вычисление истинных угловых измерений
    trueMeasParams.me = asin(arr_aspk[2]);

    //trueMeasParams.mq = atan2(arr_aspk[0], arr_aspk[1]);
    trueMeasParams.mq = calcCourseAngle(arr_aspk[0], arr_aspk[1]);
    // 8!

    // 9.завершение заполнения структуры trueMeasParams
    integrateEilersVect (trueMeasParams.marrMu);
    trueMeasParams.mTotv = mTVess;
    memcpy(trueMeasParams.marrSVess,marrVectSost, 3 * sizeof(double) );
    mSins.getEstArrEilers(mTVess,trueMeasParams.marrMuZv);

    return true;
}
//-----------------------------------------
//вычисление момента времени поступления ответного сигнала маяка
//arrTrueBeaconPos [3] - истинные координаты маяка в ГСК
// возвращает true, если сигнал проходит и возвращается
bool QPeaceVess::calcOtvetT( TTable_1D &tblPrfl,const double TObrabotki
                             ,double *arrTrueBeaconPos, double &valOtvetT)
{
   double val_tZapr = -1.,val_dtZapr_po_dT = -1.;
   // вычисление времени прохождения запросного сигнала
   if(!calc_tBeam_and_dtBeam_po_dT(  tblPrfl,arrTrueBeaconPos
                                    , val_tZapr, val_dtZapr_po_dT))
   {
       return false;
   }
   //  1!
   valOtvetT =  mTVess + 2. * val_tZapr +TObrabotki;
   bool breturn = false;
   QPeaceVess vessCur;
   for (int i =0; i< 20; ++i)
   {
       double fi =0., dfi =0.;
       // вычисление времени прохождения ответного сигнала
       //и его производной по моменту времени получения
       double val_tOtv = 0., val_dtOtv= 0.;

       vessCur = *this;
       vessCur.Move(valOtvetT);
       double val_tBeam = -1., val_dtBeam_po_dT = -1.;
       if(!vessCur.calc_tBeam_and_dtBeam_po_dT( tblPrfl,arrTrueBeaconPos
                                            , val_tOtv, val_dtOtv))
       {
           return false;
       }


       fi = valOtvetT - mTVess - val_tOtv -val_tZapr  - TObrabotki;
       dfi =1. -val_dtOtv;
       double dt = -fi/dfi;
       valOtvetT += dt;
       if (fabs(dt) < 1.E-9)
       {
           breturn  = true;
           break;
       }
    }
   return breturn;
}


//-----------------------------------------
//вычисление момента времени поступления ответного сигнала маяка
// INPUT:
//arrTrueBeaconPos [3] - истинные координаты маяка в ГСК
//arrBeaconVelo[3] - истинная корость маяка
// TObrabotki - время обработки принятого сигнала на маяке
//OUTPUT:
// valOtvetT - моментвремени приема ответного сигнала кораблем
// val_tOtv - время прохождения ответного сигнала
// val_tZapr - варемя прохождения запросного сигнала
// arrBeaconPosOtv[3] -положение маяка на мосент передачи ответногго сигнала
// возвращает true, если сигнал проходит и возвращается
bool QPeaceVess::calcOtvetT_( TTable_1D &tblPrfl,const double TObrabotki
      ,double *arrBeaconPos,double *arrBeaconVelo, double &valOtvetT
      , double &val_tOtv, double &val_tZapr, double *arrBeaconPosInp, double *arrBeaconPosOut)
{
   double val_dtZapr_po_dT = -1.;
   // вычисление времени прохождения запросного сигнала

   // 1. вычисление вектора истинного положения антенны на момент
     // передачи запросного сигнала
     double arrSAnt_GSK[6] = {0.}, arrEilers0[3] = {0.};
     integrateEilersVect (arrEilers0);
     RecalcVect_PSK_CT_INTO_GSK( mPlatform.marrPosParams, arrEilers0, NULL,arrSAnt_GSK,3 );

     if(!calc_tZapr(tblPrfl, arrSAnt_GSK
                     ,arrBeaconPos,arrBeaconVelo,  val_tZapr))
     {
         return false;
     }
     //1!

   valOtvetT =  mTVess + 2. * val_tZapr +TObrabotki;
   bool breturn = false;
   QPeaceVess vessCur;

   double  val_dtOtv= 0.;

   for (int i =0; i< 20; ++i)
   {
       memcpy(arrBeaconPosOut, arrBeaconPos, 3 * sizeof(double));
       memcpy(arrBeaconPosInp, arrBeaconPos, 3 * sizeof(double));
       if (arrBeaconVelo != NULL)
       {
           arrBeaconPosOut[0] +=  (val_tZapr +TObrabotki)*arrBeaconVelo[0];
           arrBeaconPosOut[1] +=  (val_tZapr +TObrabotki)*arrBeaconVelo[1];
           arrBeaconPosOut[2] +=  (val_tZapr +TObrabotki)*arrBeaconVelo[2];

           arrBeaconPosInp[0] +=  (val_tZapr )*arrBeaconVelo[0];
           arrBeaconPosInp[1] +=  (val_tZapr )*arrBeaconVelo[1];
           arrBeaconPosInp[2] +=  (val_tZapr )*arrBeaconVelo[2];
       }
       double fi =0., dfi =0.;
       // вычисление времени прохождения ответного сигнала
       //и его производной по моменту времени получения

       vessCur = *this;
       vessCur.Move(valOtvetT);
       double val_tBeam = -1., val_dtBeam_po_dT = -1.;
       if(!vessCur.calc_tBeam_and_dtBeam_po_dT( tblPrfl,arrBeaconPosOut
                                            , val_tOtv, val_dtOtv))
       {
           return false;
       }


       fi = valOtvetT - mTVess - val_tOtv -val_tZapr  - TObrabotki;
       dfi =1. -val_dtOtv;
       double dt = -fi/dfi;
       valOtvetT += dt;
       if (fabs(dt) < 1.E-7)
       {
           breturn  = true;
           break;
       }
    }
   return breturn;
}


//-----------------------------

//вычисление времени прохождения запросного сигнала до маяка
//возвращает false, если сигнал не проходит
bool QPeaceVess::calc_tBeam_and_dtBeam_po_dT( TTable_1D &tblPrfl,const double *arrTrueBeaconPos
                     , double &val_tBeam, double &val_dtBeam_po_dT)
{
  // 1. вычисление вектора истинного положения антенны на момент
    // передачи запросного сигнала
    double arrY[6] = {0.}, arrEilers0[3] = {0.};
    integrateEilersVect (arrEilers0);
    RecalcVect_PSK_CT_INTO_GSK( mPlatform.marrPosParams, arrEilers0, NULL,arrY,3 );


    if (arrY[2] >= -2.)
    {
        return false;//антенна над водой
    }
    //1!

    // 2.
    double val_t = 0., arrdt_po_dYAnt[3] = {0.};
    if(!calc_t_and_dt_po_dYAnt_InpGSK(tblPrfl, arrY
                 ,arrTrueBeaconPos,  val_tBeam,  arrdt_po_dYAnt))
    {
        return false;
    }
    // 2!

    // 3. вычисление dYAnt_po_dT

    double arr_dYAnt_po_dT[3]={0.};
    calc_dYAnt_po_dT(arr_dYAnt_po_dT);
    // 3!

    // 4. вычисление val_dtBeam_po_dT
    val_dtBeam_po_dT = ScalProduct(arrdt_po_dYAnt , arr_dYAnt_po_dT, 3) ;
    return true;

}
//-----------------------------

//вычисление вектора dYAnt(T)/dT производной по времени от вектора YAnt(T)
// YAnt(T) - вектор положения центра антенны в ГСК
// OUTPUT:
// arr_dYAnt_po_dT[3]
void QPeaceVess::calc_dYAnt_po_dT(double *arr_dYAnt_po_dT)
{
 // 1. вычисление вектора скорости судна
    memcpy(arr_dYAnt_po_dT, &(marrVectSost[3]),3 * sizeof(double));
    // 1!

    // 2. вычисление производной от матрицы перехода ПСК-КГСК ао времени
    double arr_dMtrxPtr_PSK_to_KGSK_po_dT[9]= {0.}, arrMu[3] = {0.}, arr_dMu[3] = {0.};
    integrateEilersVect (arrMu);

    integrate_arr_dEilers_po_dt(arr_dMu);

    QCoordSystTrsf::calc_dM_PSK_to_KGSK_po_dt (arrMu, arr_dMu, arr_dMtrxPtr_PSK_to_KGSK_po_dT);
    double arrt[3] = {0.};
    MtrxMultMatrx( arr_dMtrxPtr_PSK_to_KGSK_po_dT,3, 3, mPlatform.marrPosParams,1, arrt) ;
    MtrxSumMatrx(arr_dYAnt_po_dT, arrt, 3, 1, arr_dYAnt_po_dT) ;
}
//-------------------------
double QPeaceVess::calcCourseAngle(const double x, const double y)
{
    if(fabs(x) > 1.E6 *fabs(y) )
    {
       if (y<0.)
       {
           return -M_PI/2.;
       }
       return M_PI/2.;
    }
    double temp = atan2(x,y);
   // double temp1 = asin(x/sqrt(x * x + y * y));
    if (temp < 0.)
    {
        return temp + 2. *M_PI;
    }
    return temp;

}

//-----------------------------------------------------
//вычисление вектора положениея объекта в сферической сиситеме координат гаджета
//INPUT:
//arrXGdg[6] - вектор позиционирования гаджета в ПСК
// arrEilerCntrKP [3]- вектор углов Эйлера
//arrKGSK[3]-вектор положениея объекта в КГСК
//arrGdgPosParams[6] - вектор параметров позиционирования гаджета
// OUTPUT:
//arrVGdg[3]- вектор сфкеер координат Оb, Betta, R, Eps
void QPeaceVess::recalcPositionFromKGSK_to_GdgSphericalSK(double *arrEilerCntrKP
         , double *arrKGSK,double *arrGdgPosParams,double *arrVGdg)
{
    double arrSobSK[3] ={0.};
    recalcVect_KGSK_INTO_GdgSobSK( arrKGSK
               , arrEilerCntrKP, NULL,  arrGdgPosParams,arrSobSK,3 );
// 4. положение в АСфСАК

 QCoordSystTrsf::recalcCoord_INTO_Spherical(arrSobSK, arrVGdg[1], arrVGdg[0]
        , arrVGdg[2]);

}

//-----------------------------------------------------
//вычисление вектора положениея объекта в сферической сиситеме координат гаджета
//INPUT:
//arrXGdg[6] - вектор позиционирования гаджета в ПСК
//arrSVess_gsk[3] - вектор положения корабля в ГСК
// arrEilerCntrKP [3]- вектор углов Эйлера
//arrSObj_gsk[3]-вектор положениея объекта в ГСК
//arrGdgPosParams[6] - вектор параметров позиционирования гаджета
// OUTPUT:
//arrVGdg[3]- вектор сфкеер координат Оb, Betta, R, Eps
void QPeaceVess::recalcPositionFromGSK_to_GdgSphericalSK(double *arrSVess_gsk,double *arrEilerCntrKP
         , double *arrSObj_gsk,double *arrGdgPosParams,double *arrVGdg)
{
    // вектор положения точки в КГСК
    double arrKGSK[3] = {0.};
    MtrxMinusMatrx(arrSObj_gsk, arrSVess_gsk,1, 3, arrKGSK);
    //
    recalcPositionFromKGSK_to_GdgSphericalSK(arrEilerCntrKP
             , arrKGSK,arrGdgPosParams,arrVGdg);

}

//---------------------------------------
//вычисление вектора положениея объекта в КГСК
//INPUT:
// arrV[3]- вектор сфкеер координат объекта, Betta, R, Eps в сист координат гаджета
// arrEilerCntrKP [3]- вектор углов Эйлера
//arrGdgPosParams[6] - вектор параметров позиционирования гаджета в ПСК
// OUTPUT:
//arrKGSK[3]-вектор положениея объекта в КГСК
void QPeaceVess::recalcPositionFromGdgSphericalSK_to_KGSK(double *arrV,double *arrEilerCntrKP
                                                        ,double *arrGdgPosParams, double *arrKGSK)
{

    ///
    // 1.  положение в АСПК
    double arrS_ASPK[3] = {0.};
    QCoordSystTrsf::recalcSphericalCoord_INTO_Rectangular(arrV[1],arrV[0]
                                               ,arrV[2], arrS_ASPK);
    ///

    // 2.  положение  в ПСК-гаджет
    double matrPereh_ASK_V_PSK[9]= {0.};
    QCoordSystTrsf::calcMtrx3_ASPK_v_PSK(&arrGdgPosParams[3],matrPereh_ASK_V_PSK) ;
    double arrS_PSK_RLK[3] = {0.};
    MtrxMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrS_ASPK,1, arrS_PSK_RLK) ;
    ///

    // 3. положение в ПСК-ЦТ
    double arrS_PSK_CT[3] = {0.};

    MtrxSumMatrx(arrS_PSK_RLK, arrGdgPosParams,1, 3, arrS_PSK_CT) ;
    ///

    // 4. положение в КГСК
    double matrPereh_PSK_V_KGSK[9]= {0.};
    QCoordSystTrsf::calcMatr_PSK_v_KGSK_RightRot(arrEilerCntrKP,matrPereh_PSK_V_KGSK) ;
    MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrS_PSK_CT,1, arrKGSK) ;
    ///
}
//-------------------------------------------------------





// формирование вектора углов Эйлера в ПСК
void QPeaceVess::integrateEilersVect (double *arrMu)
{
    arrMu[0] = mQ ;
    arrMu[1] = mPsi ;
    arrMu[2] = mTet ;
}
//--------------------------------------------------------------


//------------------------------------------------
///void QPeaceVess::initTraj(const double *arrAprioriBeaconPos,const double  ZonaR
 //                         ,const enumTypeOfVessTraj TypeOfVessTraj,const double  T0)


//---------------------------------------------

int QPeaceVess::createInputDataReport(wchar_t*FileName, const bool bHeader)

{
    int len = wcslen(FileName) ;

    if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'x') // проверка, что
     && (FileName[len - 3] == 't') ) )  // указанный файл имеет расширение .txt
    {
      return 1 ;
    }

    FILE *fw ;

    if ((fw = _wfopen(FileName,L"a"))== NULL)

    {
     return 1 ;
    }
if (bHeader)
{
   fprintf(fw,"  Дата и время формирования отчета\n");
   time_t t = time(NULL);
   struct tm* aTm = localtime(&t);
   fprintf(fw,"  Год = %04d\n",aTm->tm_year+1900);
   fprintf(fw,"  Mесяц = %02d\n",aTm->tm_mon+1);
   fprintf(fw,"  День = %02d\n",aTm->tm_mday);
   fprintf(fw,"  Время = %02d:%02d:%02d\n",aTm->tm_hour, aTm->tm_min, aTm->tm_sec);
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
}
fprintf(fw,"***************************************\n");
fprintf(fw,"      Корабль свой\n");

fprintf(fw,"      Геометрические размеры:\n");
fprintf(fw,"  ширина(м) = %3.1f\n",mWidth);
fprintf(fw,"  длина(м) = %3.1f\n\n",mLength);

fprintf(fw,"      Хар-ки качек, рыскания:\n");
fprintf(fw,"  амплитуда угла рыскания (рад) = %6.5f\n",mMaxQ);
fprintf(fw,"  период рыскания (с) = %4.1f\n",mT_Q);
fprintf(fw,"  амплитуда угла килевой качки (рад) = %6.5f\n",mMaxPsi);
fprintf(fw,"  период килевой качки (с) = %4.1f\n",mT_Psi);
fprintf(fw,"  амплитуда угла бортовой качки (рад) = %6.5f\n",mMaxTet);
fprintf(fw,"  период бортовой качки (с) = %4.1f\n",mT_Tet);
fprintf(fw,"  амплитуда вертикальной качки (м) = %2.2f\n\n",mMaxVert);

fprintf(fw,"      Параметры движения:\n");
fprintf(fw,"  генеральный курс (рад) = %6.4f\n",mQ0);
fprintf(fw,"  скорость корабля своего (м/с) = %5.2f\n\n",mVVess);

fprintf(fw,"      Параметры качек:\n");
fprintf(fw,"  нач. фазы колебаний палубных углов(рад) { %6.4f; %6.4f; %6.4f }\n",marrDelt[0], marrDelt[1], marrDelt[2]);
fprintf(fw,"  нач. фаза колебаний центра тяжести(рад) = %6.4f\n",marrDelt[4]);


fprintf(fw,"      Параметры деформаций корпуса:\n");



 fclose(fw);
 mSins.createInputDataReport(FileName, false);
 //if(mpPltfOec)
// {
//    mpPltfOec->createInputDataReport(FileName, false);
 //}

}

// Аккумуляция информации о траектории корабля
// для последующего отображениия в SHP файлах
// VAlBeginT, VAlEndT - время начала и конца движения
// VAlStepT - шаг по времени
// parrData
//
void QPeaceVess::createTrajectReport(const double VAlBeginT,const double VAlEndT
                   ,const double VAlStepT,double* parrData, int *piQuantRows )
{
    int quantCircle = int(VAlEndT/VAlStepT) -1;
    int quantRows = 0;
    for (int i =0; i < quantCircle; ++i )
    {

        double valT = ((double)i) *VAlStepT;
        recalcVess( valT);
        if ((valT >= VAlBeginT)&& (valT <= VAlEndT))
        {

           int num0 =  quantRows * QUantCols_VessTrajReport;

           parrData [num0] =     mTVess ;

           parrData [num0 +1]  =  marrVectSost [0];  // arrShipGSK_X

           parrData [num0 +2]  =  marrVectSost [1] ; // arrShipGSK_Y

           parrData [num0 +3]  =  marrVectSost [2] ; // arrShipGSK_Z

           parrData [num0 +4]  = marrVectSost [3] ;   // arrShipGSK_VX

           parrData [num0 +5]  = marrVectSost [4];  // arrShipGSK_VY

           parrData [num0 +6]  = marrVectSost [5]  ; // arrShipGSK_VZ

            parrData [num0 + 7]  =  mVVess;

           parrData [num0 + 8]  =  mQ; //

           parrData [num0 + 9]  =  mPsi;

           parrData [num0 + 10]  = mTet;

           parrData [num0 + 11]  = mVQ;

           parrData [num0 + 12]  = mVPsi;

           parrData [num0 + 13]  = mVTet;


            quantRows++ ;
        }
    }
    *piQuantRows = quantRows;
//
}

//---------------------------------------
// формирование вектора угловой скорости в подвижной СК
void  QPeaceVess::integrate_arr_dEilers_po_dt(double *arr_dMu)
{
    arr_dMu[0] = mVPsi;
    arr_dMu[1] = mVTet;
    arr_dMu[2] = mVQ;
}

//---------------------------------------------
 // Сивухин Д.В. Общий курс физики, т 1, стр 339
    // вычисление вектора ускорения инерции цели в ПСК
    //OUTPUT
    // arrInertAccel[3]

    // Aабс = Аотн +Апер+ Акор
    // Акор = 2[Om;Vотн]
    // Апер = dV0/dt +[Om;[Om;r]] +[dOm/dt;r]
    // надо найти Аотн
    // Известно, что:
    // 1. Aабс =0
    // 2. dV0/dt = 0
    // dOm/dt не известно, положим dOm/dt = 0
    // тогда Аотн =  -Апер- Акор= -2[Om;Vотн]-[Om;[Om;r]]
    // необходимо принять во внимание, что в морском деле углы исчисляются по часовой стрелке
    // то есть, положительное напроавление отсчета углов ведется по
    // часовой стрелке, вместо того, чтобы вестись против.
    // Скорости измерения углов (угловые скорости) по этой же причине
    // также имеют противоположный знак.
    // По этой причине, результирующая формула для инерциального ускорения будет иметь вид
    // Аотн =  -Апер- Акор= 2[Om;Vотн]-[Om;[Om;r]] = [Om;2*Vотн-[Om;r]]
    // Om - вектор угловой скорости подвижной системы координат относительно неподвижной
    // Vотн - вектор скорости тела в подвижной системе координат
    // r - радиус вектор положения цели в подвижной системе координат
    // ВЫЧИСЛЕНИЕ ДИСПЕРСИИ РАЗБРОСА ОЦЕНКИ ВЕКТОРА ИНЕРЦИАЛЬНОГО УСКОРЕНИЯ
    // del = dV0/dt + [dOm/dt;r]
    // -> |del| <= |dV0/dt| +|[dOm/dt;r]|<= 0.2 + |dOm/dt| * |r|
    // dOm/dt [i] = Ai *(2*M_PI/Ti) * (2*M_PI/Ti) * (2*M_PI/Ti)
/*
void QPeaceVess::calcInertAccel_PSK( double *arrS_PSK, double *arrEilers0
                      ,double *arr_dEilers_po_dt0
                      ,double *arrInertAccel)
{
    double arrOm0[3] = {0.};
    calcVectOm(arrEilers0, arr_dEilers_po_dt0, arrOm0);


    double arrEilers[3] = {0.};
    MatrxMultScalar(arrEilers0, 1, 3, -1.,arrEilers);
    double matrPereh_PSK_V_KGSK[9] = {0} ;
    QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilers, matrPereh_PSK_V_KGSK) ;


    double arrOm[3] = {0.};
    MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrOm0,1, arrOm) ;
       // double arr_dEilers_po_dt0[3] = {0.};
      //  memcpy(arr_dEilers_po_dt0, arr_dEilers_po_dt, 3 * sizeof(double));
   // MatrxMultScalar(arr_dEilers_po_dt0, 1, 3, -1., arr_dEilers_po_dt0);

    double arrT0[3] = {0.}, arrTCoriol[3] = {0.},arrT2[3] = {0.},arrT3[3] = {0.}
             ,arrT4[3] = {0.},arrT5[3] = {0.}, arrPeren[3] = {0.};
    // кориолисово ускорение Акор = 2[Om;Vотн]
    MatrxMultScalar(&(arrS_PSK[3]), 3, 1, 2.,arrT0);
   // MatrxMultScalar(&(arrS_PSK[3]), 3, 1, 1.,arrT0);
    OuterProduct(arrOm , arrT0, arrTCoriol) ;
    ///
    // Апер = dV0/dt +[Om;[Om;r]] +[dOm/dt;r]
    // dV0/dt = 0
    OuterProduct(arrOm , arrS_PSK, arrT2) ;
    OuterProduct(arrOm, arrT2, arrPeren) ;
    ///
   MatrxMultScalar(arrTCoriol, 1, 3, -1., arrTCoriol);
    MtrxMinusMatrx(arrTCoriol, arrPeren,1, 3, arrInertAccel);

    double arr_dOm[3] = {0.};
    double psi = arrEilers0[0];// Q
    double nu = arrEilers0[1]; // Psi
    double gam = arrEilers0[2]; // Tet
    double psit = arr_dEilers_po_dt0[2];
    double nut = arr_dEilers_po_dt0[0];
    double gamt = arr_dEilers_po_dt0[1];

    arr_dOm[1] =  psit *nut *cos(nu);
    arr_dOm[2] = nut * gamt * cos(gam) - psit *(gamt * sin(gam) * cos(nu) + nut * cos(gam) * sin(nu));
    arr_dOm[0] = -nut * gamt * sin(gam) - psit * gamt * cos(gam)*cos(nu) + psit * nut * sin(gam)*sin(nu);

    OuterProduct(arr_dOm , arrS_PSK, arrT3) ;

    MtrxMinusMatrx(arrInertAccel, arrT3,1, 3, arrInertAccel);

}
*/



    /*
// экстраполяция вектора состояния объекта класса QPeaceVess на время  valTExtr
// INPUT:
// valTExtr - момент времени, на короторый требуется произвести экстраполяцию
// OUTPUT:
// arrVSVessExtrap - экстраполированный вектор состояния
void QPeaceVess::VSProlong(const double valTExtr, double *arrVSVessExtrap)
{

  // экстраполяция вектора состояния корабля на момент  valT

  double h = valTExtr - mTVess ;
  memcpy ( arrVSVessExtrap, marrEstVectSost, 6 * sizeof(double)) ;
  const double valModVessV = sqrt(marrEstVectSost[3] * marrEstVectSost[3] + marrEstVectSost[4] * marrEstVectSost[4]) ;
  // оценка курса корабля
   double valEst_Q = mSins.mEstQ  ;
   double valEst_dQ_Po_dt = mSins.mEstVQ ;

  const double valQExtr = valEst_Q + valEst_dQ_Po_dt * h ;
  arrVSVessExtrap [2] +=  arrVSVessExtrap [5] * h ;

  if (fabs(valEst_dQ_Po_dt) < 1E-15)
  {
  arrVSVessExtrap [0] +=  valModVessV *  sin (valEst_Q ) * h ;
  arrVSVessExtrap [1] +=  valModVessV *  cos (valEst_Q ) * h ;
  }
  else
  {
  arrVSVessExtrap [0] -=  valModVessV * ( cos (valQExtr ) - cos (valEst_Q) )/ valEst_dQ_Po_dt ; //      ( sin (valEst_Q) *  h +  valEst_dQ_Po_dt * cos (valEst_Q) *  h *h / 2.);
  arrVSVessExtrap [1] +=  valModVessV * ( sin (valQExtr ) - sin (valEst_Q) )/ valEst_dQ_Po_dt ; // ( cos (valEst_Q) *  h  -  valEst_dQ_Po_dt * sin (valEst_Q) *  h *h / 2.);
  arrVSVessExtrap [3] =  valModVessV * sin (valEst_Q + valEst_dQ_Po_dt * h) ;
  arrVSVessExtrap [4] =  valModVessV * cos (valEst_Q + valEst_dQ_Po_dt * h) ;
  }
}


// пересчет замера   из кажущейся АСфСК в кажущуюся ГСК  с корреляц мвтрицей
//  OUTPUT:
// arrZam_GSK - вектор замера в ГСК - 3-х мерный массив
// arrKZam_GSK - коррел матрица замера, 9-ти мерный массив
//
void  QPeaceVess::GetZamer_IN_GSK (TZamer InpASKZamer,  TZamer *pOutGSKZamer )

{
    GetZamer_IN_KGSK (InpASKZamer, pOutGSKZamer);
  for (int i = 0; i < 3; i++) (*pOutGSKZamer).marrMeas[i] +=  marrEstVectSost [i] ;



}

// расчет результирующей корреляционной матрицы ошибок ищзмерения замера, включающней ошибки
 // флуктуацинного происх ,  ошибки привода и ошибки определония качек в ГСК  и прогиба корпуса корабля
// все ошибки складываютя под корнем
// INPUT:
//InpZamer - входной замер в АСК
// InpZamer.marrMeas [0] - угол V
// InpZamer.marrMeas [1] - дальность R
// InpZamer.marrMeas[2] - угол U
 void  QPeaceVess::calcSummarizedCorMtrx_ErrMes_In_GSK (TZamer InpASKZamer,double *arrCorrMtrx_GSK )
{
   double ar_dF_po_dBet_sq [9],ar_dF_po_dEps_sq [9], matrPereh_ASK_V_KGSK[9],arrMu[5] ;
       double  arrDT0[9] = {0}, arrDT1[9] = {0}, arrDT2[9] = {0}
                ,arrDT3[9] = {0},arrDT4[9] = {0},arrDT5[9] = {0},arrDT51[9] = {0};
    arrMu[0] = mSins.mEstQ;
    arrMu[1] = mSins.mEstPsi ;
    arrMu[2] = mSins.mEstTet ;
    arrMu[3] = mDriver.mEstBet ;
    arrMu[4] = mDriver.mEstEps ;

    double valDispV = InpASKZamer.marrCorr[0];
    double valDispR = InpASKZamer.marrCorr[4];
    double valDispU = InpASKZamer.marrCorr[8];
         calc_dF_po_dBet_sq(mDriver.mEstBet,  mDriver.mEstEps,ar_dF_po_dBet_sq ) ;
       calc_dF_po_dEps_sq(mDriver.mEstBet,  mDriver.mEstEps,ar_dF_po_dEps_sq ) ;
       calcMatr_ASK_v_KGSK(arrMu,matrPereh_ASK_V_KGSK) ;
         arrDT0[0] = valDispV * InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1] ;
         arrDT0[4] = valDispR ;
         arrDT0[8] = valDispU * InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1] ;
       MtrxMultMatrx(matrPereh_ASK_V_KGSK,3, 3, arrDT0,3, arrDT1) ;
       MtrxMultMatrxTransp(arrDT1,3, 3, matrPereh_ASK_V_KGSK, 3, arrDT2) ;

         double temp0 =  InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1]
         * mDriver.mSigBet* mDriver.mSigBet ;
       MatrxMultScalar(ar_dF_po_dBet_sq , 3, 3, temp0,arrDT3);
       MtrxSumMatrx(arrDT2, arrDT3,3,3, arrDT4) ;

         double temp1 =  InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1]
       * mDriver.mSigEps* mDriver.mSigEps ;
       MatrxMultScalar(ar_dF_po_dEps_sq , 3, 3, temp1,arrDT5);
       MtrxSumMatrx(arrDT4, arrDT5,3,3, arrDT51) ;

       // расчет матриц, вызванных ошибками определения качек
       double ar_dF_po_dTet_sq [9] = {0},ar_dF_po_dQ_sq [9] = {0},ar_dF_po_dPsi_sq [9] = {0}, arrDT6[9] = {0}
          , arrDT7[9]  = {0}, arrDT8[9] = {0}, arrDT9[9] = {0}, arrDT10[9]  = {0};
       calc_dF_po_dQ_sq  (mDriver.mEstBet,  mDriver.mEstEps, ar_dF_po_dQ_sq) ;
       calc_dF_po_dPsi_sq(mDriver.mEstBet,  mDriver.mEstEps, ar_dF_po_dPsi_sq) ;
       calc_dF_po_dTet_sq(mDriver.mEstBet,  mDriver.mEstEps, ar_dF_po_dTet_sq) ;

         double temp2 = InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1]
         * mSins.mSig_Q * mSins.mSig_Q  ;
         double temp5 = calcAmpAftFlexure(marrParral[1]);
       double temp3 = InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1]
         *( mSins.mSig_Psi* mSins.mSig_Psi + temp5 * temp5/2.);
         double temp6 = calcAmpBoardFlexure(marrParral[1]);
       double temp4 = InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1]
         * (mSins.mSig_Tet* mSins.mSig_Tet + temp6 * temp6/2.);

       MatrxMultScalar(ar_dF_po_dQ_sq   , 3, 3, temp2,arrDT6);
       MatrxMultScalar(ar_dF_po_dPsi_sq , 3, 3, temp3,arrDT7);
       MatrxMultScalar(ar_dF_po_dTet_sq , 3, 3, temp4,arrDT8);

       MtrxSumMatrx(arrDT51, arrDT6,3,3, arrDT9) ;
       MtrxSumMatrx(arrDT9, arrDT7,3,3, arrDT10) ;
       MtrxSumMatrx(arrDT10, arrDT8,3,3, arrCorrMtrx_GSK) ;
 }


 // имитация замера в АСК c учетом одного антипода
void QPeaceVess::createMeasure(const double VAlTrue_R,const double VAlTrue_V, const double VAlTrue_U
    , const double VAlT, const double VAlTargDesEps, const double VAlTargDesBet
     , TEtalonSign EtalonSign , const double VAlTargEPR
    , const double VAlPowerPrd, const double VAlKYPrd, TZamer *pOutASKZamer )
{
    // вычисление высоты антенны над морской поверхностью в заданный момент времени
    double arrFarKGSK[3] ={0.};
    RecalcVect_PSK_CT_True_INTO_KGSK_True (marrParral
    ,VAlTargDesEps,  VAlTargDesBet,arrFarKGSK,3) ;

    ///
    double arrMeas[3] = {0.};
    arrMeas[0] = VAlTrue_V;
    arrMeas[1] = VAlTrue_R ;
    arrMeas[2] = VAlTrue_U;
    double arrCorr [9] = {0.};
    TZamer InpASKZamer( arrMeas,arrCorr, 0.) ;
    TZamer OutKGSKZamer;
    GetZamer_IN_KGSK( InpASKZamer, &OutKGSKZamer);
     ///

     double valAntpPhaze = 0.;
     double valSigmaEps = -1., valSigmaBet = -1.0;
        mFar_2D.calc_SKZ_LAT(VAlTrue_R ,OutKGSKZamer.marrMeas[2]
        ,arrFarKGSK[2],  VAlTargDesEps, VAlTargEPR
        ,EtalonSign ,VAlPowerPrd,  VAlKYPrd,  valAntpPhaze
        ,&valSigmaEps ,&valSigmaBet);


     (*pOutASKZamer).marrMeas[0] = VAlTrue_V + getGauss(0., valSigmaBet );
     (*pOutASKZamer).marrMeas[1] = VAlTrue_R + getGauss(0., 10. );
     (*pOutASKZamer).marrMeas[2] = VAlTrue_U + getGauss(0., valSigmaEps );
     memset((*pOutASKZamer).marrCorr, 0, 9 * sizeof(double));
     (*pOutASKZamer).marrCorr[0] = valSigmaBet * valSigmaBet;
     (*pOutASKZamer).marrCorr[4] = 100.;
     (*pOutASKZamer).marrCorr[8] = valSigmaEps * valSigmaEps;
     (*pOutASKZamer).mT = VAlT;
     ///

//	mFar_2D.createMeasure_Targ_Plus_Measure_Antp(const double VAlTrue_R,const double VAlTrue_V, const double VAlTrue_U
 //	, const double VAlT,
}


void QPeaceVess::createMeasure_ForSingleTarg(const double VAlTrue_R,const double VAlTrue_V, const double VAlTrue_U
    , const double VAlT, const double VAlTargDesEps, const double VAlTargDesBet
     , TEtalonSign EtalonSign , const double VAlTargEPR
    , const double VAlPowerPrd, const double VAlKYPrd, TZamer *pOutASKZamer )
{
    // вычисление высоты антенны над морской поверхностью в заданный момент времени
    double arrFarKGSK[3] ={0.};
    RecalcVect_PSK_CT_True_INTO_KGSK_True (marrParral
    ,VAlTargDesEps,  VAlTargDesBet,arrFarKGSK,3) ;

    double arrMeas[3] = {0.};
    arrMeas[0] = VAlTrue_V;
    arrMeas[1] = VAlTrue_R ;
    arrMeas[2] = VAlTrue_U;
    double arrCorr [9] = {0.};
    TZamer InpASKZamer( arrMeas,arrCorr, 0.) ;
    TZamer OutKGSKZamer;
    GetZamer_IN_KGSK( InpASKZamer, &OutKGSKZamer);
     ///
      double valSigmaEps = -1., valSigmaBet = -1.0;
    mFar_2D.calc_SKZ_SingleTarg(VAlTrue_R, VAlTargDesEps,  VAlTargEPR
  ,  EtalonSign ,VAlPowerPrd,  VAlKYPrd ,&valSigmaEps ,&valSigmaBet);

     (*pOutASKZamer).marrMeas[0] = VAlTrue_V + getGauss(0., valSigmaBet );
     (*pOutASKZamer).marrMeas[1] = VAlTrue_R + getGauss(0., 10. );
     (*pOutASKZamer).marrMeas[2] = VAlTrue_U + getGauss(0., valSigmaEps );
     memset((*pOutASKZamer).marrCorr, 0, 9 * sizeof(double));
     (*pOutASKZamer).marrCorr[0] = valSigmaBet * valSigmaBet;
     (*pOutASKZamer).marrCorr[4] = 100.;
     (*pOutASKZamer).marrCorr[8] = valSigmaEps * valSigmaEps;
     (*pOutASKZamer).mT = VAlT;
     ///

//	mFar_2D.createMeasure_Targ_Plus_Measure_Antp(const double VAlTrue_R,const double VAlTrue_V, const double VAlTrue_U
 //	, const double VAlT,
}


//определение положенпия АУ в текущий момоент в КГСК
// ЭТО ЗАГЛУШКА!!!!!
void  __fastcall QPeaceVess::calcAY_Position(double *arr_AY_Position_KGSK)
{
    memset(arr_AY_Position_KGSK, 0, sizeof(double) * 3);
}


 //--------------------------------------------------------------
 // вычисление дисперсий ишибок  линейной экстраполяцмии палубных углов СИНС
 // INPUT:
 // VAlTExtrap - время экстраполяции
 // OUTPUT:
 // pvalDispDEltaExtrIS_Psi, *pvalDispDEltaExtrIS_Tet , *pvalDispDEltaExtrIS_Q - дисперсии
 void QPeaceVess::calcDispExtrapDeckAngles_SINS(const double VAlTExtrap, double *pvalDispDEltaExtrIS_Psi
     , double *pvalDispDEltaExtrIS_Tet, double *pvalDispDEltaExtrIS_Q)
{
    *pvalDispDEltaExtrIS_Psi = mSins.mSig_Psi * mSins.mSig_Psi +  VAlTExtrap *VAlTExtrap * mSins.mSig_dPsidt * mSins.mSig_dPsidt
        + VAlTExtrap * VAlTExtrap *VAlTExtrap * VAlTExtrap/4. * mMaxPsi* mMaxPsi * ( 2. * M_PI / mT_Psi)* ( 2. * M_PI / mT_Psi)
        * ( 2. * M_PI / mT_Psi)* ( 2. * M_PI / mT_Psi)/ 2.;

    *pvalDispDEltaExtrIS_Tet = mSins.mSig_Tet * mSins.mSig_Tet +  VAlTExtrap *VAlTExtrap * mSins.mSig_dTetdt * mSins.mSig_dTetdt
        + VAlTExtrap * VAlTExtrap *VAlTExtrap * VAlTExtrap/4. * mMaxTet* mMaxTet * ( 2. * M_PI / mT_Tet)* ( 2. * M_PI / mT_Tet)
        * ( 2. * M_PI / mT_Tet)* ( 2. * M_PI / mT_Tet)/ 2.;

  *pvalDispDEltaExtrIS_Q = mSins.mSig_Q * mSins.mSig_Q +  VAlTExtrap * VAlTExtrap * mSins.mSig_dQdt * mSins.mSig_dQdt
        + VAlTExtrap * VAlTExtrap *VAlTExtrap * VAlTExtrap/4. * mMaxQ* mMaxQ * ( 2. * M_PI / mT_Q)* ( 2. * M_PI / mT_Q)
        * ( 2. * M_PI / mT_Q)* ( 2. * M_PI / mT_Q)/ 2.;
}


 //Вычисление истинных углов наклона палубы в точке установки РЛС
 // OUTPUT:
 // arrMu [5] - массив углов
 //Q, Psi, Tet, Bet, Eps
 //
  void QPeaceVess::calcVectTrueDeckAngles_For_Far(double *arrMu)
  {
      double  valVQ =0., valVPsi =0., valVTet = 0.;
     calcDeckAngles(mTVess, marrParral,&arrMu[0], &valVQ,&arrMu[1], &valVPsi,&arrMu[2] ,&valVTet   );
     arrMu[3] = mDriver.mRealBet ;
     arrMu[4] = mDriver.mRealEps ;
  }
*/
//---------------------------------------
//пересчет измерения из сферической системы координат гаджета в КГСК
//INPUT:
// *pZamerGdgSh-  измерение в сист координат гаджета
// arrEilerCntrKP [3]- вектор углов Эйлера
//arrGdgPosParams[6] - вектор параметров позиционирования гаджета в ПСК
// OUTPUT:
//*pZamerKGSK-измерение объекта в КГСК
/*
void QPeaceVess::recalcZamerFrom_GdgSphericalSK_to_KGSK(QZamer *pZamerGdgSh,double *arrEilerCntrKP
                                                        ,double *arrGdgPosParams, QZamer**ppZamerKGSK)
{
    if (pZamerGdgSh == NULL)
    {
      *ppZamerKGSK = NULL;
      return;
    }
    ///
    // 1.  положение в АСПК
    double arrS_ASPK[3] = {0.};
    QCoordSystTrsf::recalcSphericalCoord_INTO_Rectangular((*pZamerGdgSh).marrMeas[1],(*pZamerGdgSh).marrMeas[0]
                                               ,(*pZamerGdgSh).marrMeas[2], arrS_ASPK);
    ///

    // 2.  положение  в ПСК-гаджет
    double matrPereh_ASK_V_PSK[9]= {0.};
    QCoordSystTrsf::calcMtrx3_ASPK_v_PSK(&arrGdgPosParams[3],matrPereh_ASK_V_PSK) ;
    double arrS_PSK_RLK[3] = {0.};
    MtrxMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrS_ASPK,1, arrS_PSK_RLK) ;
    ///

    // 3. положение в ПСК-ЦТ
    double arrS_PSK_CT[3] = {0.};

    MtrxSumMatrx(arrS_PSK_RLK, arrGdgPosParams,1, 3, arrS_PSK_CT) ;
    ///

    // 4. положение в КГСК
    double matrPereh_PSK_V_KGSK[9]= {0.};
    QCoordSystTrsf::calcMatr_PSK_v_KGSK_RightRot(arrEilerCntrKP,matrPereh_PSK_V_KGSK) ;
    MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrS_PSK_CT,1, (**ppZamerKGSK).marrMeas) ;
    ///
    //5. коррел матрица
    double parrRez[9] = {0.}, parrRez0[9] = {0.}, parrRez00[9] = {0.};
    MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3,  matrPereh_ASK_V_PSK,3, parrRez);

            double arr_dA_po_dV[9] = {0.};
            QCoordSystTrsf::calc_dRectangularl_po_dSpherical((*pZamerGdgSh).marrMeas, arr_dA_po_dV);
   MtrxMultMatrx(parrRez,3, 3,  arr_dA_po_dV,3, parrRez0);
   memcpy(parrRez00, parrRez0, 9 * sizeof(double));

   //MtrxMultMatrx_MultMatrxTransp(parrRez0,(*pZamerGdgSh).marrCorr,parrRez00,3, (**ppZamerKGSK).marrCorr)  ;
  // (**ppZamerKGSK).mT = (*pZamerGdgSh).mT;

}

*/
/*
//---------------------------------------
//пересчет измерения из сферической системы координат гаджета в ПСК
//INPUT:
// *pZamerGdgSh-  измерение в сист координат гаджета
// arrEilerCntrKP [3]- вектор углов Эйлера
//arrGdgPosParams[6] - вектор параметров позиционирования гаджета в ПСК
// OUTPUT:
//*pZamerPSK-измерение объекта в КГСК
void QPeaceVess::recalcZamerFrom_GdgSphericalSK_to_PSK(QZamer *pZamerGdgSh,double *arrEilerCntrKP
                                                        ,double *arrGdgPosParams, QZamer** ppZamerPSK)
{
    if (pZamerGdgSh == NULL)
    {
      *ppZamerPSK = NULL;
      return;
    }
    ///
    // 1.  положение в АСПК
    double arrS_ASPK[3] = {0.};
    QCoordSystTrsf::recalcSphericalCoord_INTO_Rectangular((*pZamerGdgSh).marrMeas[1],(*pZamerGdgSh).marrMeas[0]
                                               ,(*pZamerGdgSh).marrMeas[2], arrS_ASPK);
    ///

    // 2.  положение  в ПСК-гаджет
    double matrPereh_ASK_V_PSK[9]= {0.};
    QCoordSystTrsf::calcMtrx3_ASPK_v_PSK(&arrGdgPosParams[3],matrPereh_ASK_V_PSK) ;
    double arrS_PSK_RLK[3] = {0.};
    MtrxMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrS_ASPK,1, arrS_PSK_RLK) ;
    ///

    // 3. положение в ПСК-ЦТ
    MtrxSumMatrx(arrS_PSK_RLK, arrGdgPosParams,1, 3, (**ppZamerPSK).marrMeas) ;
    ///


    //5. коррел матрица
    double  parrRez0[9] = {0.}, parrRez00[9] = {0.};

    double arr_dA_po_dV[9] = {0.};
    QCoordSystTrsf::calc_dRectangularl_po_dSpherical((*pZamerGdgSh).marrMeas, arr_dA_po_dV);
   MtrxMultMatrx(matrPereh_ASK_V_PSK,3, 3,  arr_dA_po_dV,3, parrRez0);
   memcpy(parrRez00, parrRez0, 9 * sizeof(double));

   MtrxMultMatrx_MultMatrxTransp(parrRez0,(*pZamerGdgSh).marrCorr,parrRez00,3, (**ppZamerPSK).marrCorr)  ;
   (**ppZamerPSK).mT = (*pZamerGdgSh).mT;

}
*/
//------------------------------------------



