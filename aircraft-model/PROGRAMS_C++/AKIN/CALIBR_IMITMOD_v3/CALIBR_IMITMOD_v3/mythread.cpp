#include "mythread.h"
#include <QDebug>
#include "LblSolver.h"
#include "Solver2D_Angs.h"
#include "Solver3D_Angs.h"
#include "BigMeasure.h"
#include "MatrixProccess.h"
#include<math.h>

MyThread::MyThread(QObject *parent) :
    QThread(parent),
  userBreakRequested(false)
 // ,maxStep(1),
 // errorThreshold(0.1)
{    
}


void MyThread::run()
{
    qDebug() << "MyThread::run(), in";

    int result = MYTHREAD_RESULT_MAX_STEP; // результат по умолчанию
    int icur = 0;

    int  numPart = 1;



        // !следующую строку следует заменить на реальные вычисления!
        //QThread::msleep(100); // этот вызов здесь нужен для имитации длительности вычислений
      //****************************
        //const int QuantMeas = ;
         QBigMeasure *parrMeas = (QBigMeasure *)malloc((mvctMeas.size()) * sizeof(QBigMeasure) );
         for (int i = 0; i < mvctMeas.size(); ++i )
         {
            parrMeas[i] =  mvctMeas.at(i);
         }
      QLblSolver LblSolver(parrMeas, mvctMeas.size()    ,mtblEstPrfl
                                ,mHBeaconEst, 0.015);



        double arrX0[8] ={0.},arrXRez [8] = {0.}, valNeviaz0 = -1.;
        for (int i = 0; i < 8; ++i)
        {
           arrX0[i] =  mvctX0.at(i);
        }
        QPosSolver *pPosSolver = &LblSolver;
        double sigmaTreshold = 2.5;
        const int NUmVar = 0;
        QVector<int> vctGoodZameresNum;
        pPosSolver->repaireZamerArray_(arrX0, sigmaTreshold
                                       ,NUmVar, &vctGoodZameresNum);

        mvctMeas = pPosSolver->mVectBigMeasures;      

        double arrMtrx_dFGr_po_dX_Inv[25] = {0.};
        double valSumVeviaz2__ = 0.;
        double  parrF [5] = {0.};

        //  ОТЛАДКА
      //  double arrMtrx_dFGr_po_dX_Inv_[25] = {0.},arrt[25] = {0.};
     //  pPosSolver->calc_FGr_and_dFGr_po_dX(arrX0
                //         ,parrF,arrt, valSumVeviaz2__);
      // InverseMtrx( arrt, 5, arrMtrx_dFGr_po_dX_Inv_) ;
      // memcpy(arrXRez, arrX0, 8 * sizeof(double));

    // MYTHREAD_RESULT THREAD_RESULT = MYTHREAD_RESULT_OK;
       // ОТЛАДКА !

   MYTHREAD_RESULT THREAD_RESULT = estimateParams(pPosSolver,  arrX0
        ,1,mMaxQuantIter_LBL, mLBLcoeff,arrXRez, valNeviaz0, arrMtrx_dFGr_po_dX_Inv);
        qDebug() << "MyThread::run(), LBL done";
      switch (THREAD_RESULT)
      {

      case MYTHREAD_RESULT_MAX_STEP:   // вычисления провалены (достигли максимально допустимого кол-ва итераций)
      case MYTHREAD_RESULT_USER_CANCEL: // вычисления провалены (прервано пользователем)
      case MYTHREAD_RESULT_ERROR:
          delete []parrMeas;
          emit calcfinished(THREAD_RESULT);
          return;
          break;

      default:
          break;
      }
      QVector <double>vct_arr_dX_po_dSgps(15);
      QVector <double>vct_arr_dX_po_dHBeacon(5);
      QVector <double>vct_arr_Kgps(25);
      QVector <double>vct_arr_Kt(25);
      QVector <double>vct_arr_dX_po_AngSins(15);
      QVector <double>vct_arr_KsinsQ_per_1_rad(25);
      QVector <double>vct_arr_Kgps_Psi_Tetta_per_1_rad(25);
      QVector <double>vct_arr_dX_po_dProfile(5);
      QVector <double>vct_arr_KProfile_per_1_m(25);


      // передаются векторы без содержания.
      // это надо,чтобы вывесить текст с информацией о том,что идет расчет точности
      emit LBLSucceeded(              true,
                                      vct_arr_dX_po_dSgps
                                   , vct_arr_dX_po_dHBeacon
                                   , vct_arr_Kgps
                                   , vct_arr_Kt
                                   , vct_arr_dX_po_AngSins
                                   , vct_arr_KsinsQ_per_1_rad
                                   , vct_arr_Kgps_Psi_Tetta_per_1_rad
                                   , vct_arr_dX_po_dProfile
                                   , vct_arr_KProfile_per_1_m
                                   );

      const double VAlMaxHeightTreshold = -2.;
  //  pPosSolver->cleanZamersInAccordanceWithAntHeight(arrXRez, VAlMaxHeightTreshold);
   //   double arrFGr[5] = {0.};
   // double arrMtrx_dFGr_po_dX_Inv_[100] = {0.};
    // LblSolver.calc_arrFGr_and_Mtrx_dFGr_po_dX_Inv(arrXRez,arrFGr,arrMtrx_dFGr_po_dX_Inv_);


      double arr_dX_po_dSgps[15]= {0.};
    LblSolver.calc_dX_po_dSgps(arrXRez,arrMtrx_dFGr_po_dX_Inv,arr_dX_po_dSgps);
        qDebug() << "1.LBL calc_dX_po_dSgps done" ;

      double arr_dX_po_dHBeacon[5] = {0.};
  LblSolver.calc_dX_po_dHBeacon(arrXRez,arrMtrx_dFGr_po_dX_Inv,arr_dX_po_dHBeacon);
    qDebug() << "2.LBL calc_dX_po_dHBeacon" ;

      double arr_Kgps[25] = {0.};
    LblSolver.calc_Kgps(arrXRez,arrMtrx_dFGr_po_dX_Inv,arr_Kgps);
     qDebug() << "3.LBL calc_Kgps" ;

      double arr_Kt[25] = {0.};
    LblSolver.calc_Kt(arrXRez,arrMtrx_dFGr_po_dX_Inv,arr_Kt); //+
    qDebug() << "4.LBL calc_Kt" ;

      double arr_dX_po_AngSins[5*3] = {0.};
   LblSolver.calc_dX_po_dAngSins(arrXRez,arrMtrx_dFGr_po_dX_Inv,arr_dX_po_AngSins);//+
      qDebug() << "5.LBL calc_dX_po_dAngSins" ;

      double arr_KsinsQ_per_1_rad[5*5] = {0.};
   LblSolver.calc_KsinsQ(arrXRez,arrMtrx_dFGr_po_dX_Inv,arr_KsinsQ_per_1_rad);//+
      qDebug() << "6.LBL calc_KsinsQ" ;

      double arr_Kgps_Psi_Tetta_per_1_rad[25] ={0.};
   LblSolver.calc_Ksins_Psi_Tetta(arrXRez,arrMtrx_dFGr_po_dX_Inv,arr_Kgps_Psi_Tetta_per_1_rad);//+
      qDebug() << "7.LBL calc_Ksins_Psi_Tetta" ;

      double arr_dX_po_dProfile[5] = {0.};
  LblSolver.calc_dX_po_dProfile(arrXRez,arrMtrx_dFGr_po_dX_Inv,arr_dX_po_dProfile);//+
     qDebug() << "8. LBL calc_dX_po_dProfile" ;


      double arr_KProfile_per_1_m[25] ={0.};
     // LblSolver.calc_KProfile(arrXRez,arrMtrx_dFGr_po_dX_Inv,arr_KProfile_per_1_m);

      QPosSolver::createDblVect(arr_dX_po_dSgps,15,vct_arr_dX_po_dSgps);//+

      QPosSolver::createDblVect(arr_dX_po_dHBeacon,5,vct_arr_dX_po_dHBeacon);//+

      QPosSolver::createDblVect(arr_Kgps,25,vct_arr_Kgps);//+

      QPosSolver::createDblVect(arr_Kt,25,vct_arr_Kt);//+

      QPosSolver::createDblVect(arr_dX_po_AngSins,15,vct_arr_dX_po_AngSins);//+

      QPosSolver::createDblVect(arr_KsinsQ_per_1_rad,25,vct_arr_KsinsQ_per_1_rad);//+

      QPosSolver::createDblVect(arr_Kgps_Psi_Tetta_per_1_rad,25,vct_arr_Kgps_Psi_Tetta_per_1_rad);//+

      QPosSolver::createDblVect(arr_dX_po_dProfile,5,vct_arr_dX_po_dProfile); //+

      QPosSolver::createDblVect(arr_KProfile_per_1_m,25,vct_arr_KProfile_per_1_m);


     // передача корреляционных матриц
      emit LBLSucceeded(              false,
                                     vct_arr_dX_po_dSgps
                                  , vct_arr_dX_po_dHBeacon
                                  , vct_arr_Kgps
                                  , vct_arr_Kt
                                  , vct_arr_dX_po_AngSins
                                  , vct_arr_KsinsQ_per_1_rad
                                  , vct_arr_Kgps_Psi_Tetta_per_1_rad
                                  , vct_arr_dX_po_dProfile
                                  , vct_arr_KProfile_per_1_m
                                  );

      if (mTypeOfSolverTask == LBL)
      {
          delete []parrMeas;
          emit calcfinished(THREAD_RESULT);
          return;
      }


//***************************************
      numPart = 2;

    qDebug() << "MyThread::run(), out";

    //*********************************
    double arrSBeacon_GSK[3] = {0.};
    arrSBeacon_GSK[0] = arrXRez[3];
    arrSBeacon_GSK[1] = arrXRez[4];
    arrSBeacon_GSK[2] = mHBeaconEst;

    QSolver3D_Angs Solver3D_Angs;
    QSolver2D_Angs Solver2D_Angs;
    QSolver2D_Angs *pSolver2D_Angs;


    switch(mTypeOfSolverTask )
        {

        case USBL_2D:
            Solver2D_Angs = QSolver2D_Angs(parrMeas, mvctMeas.size(),mtblEstPrfl
                 ,mHBeaconEst,0.0009,arrXRez,arrSBeacon_GSK);

            pSolver2D_Angs = &Solver2D_Angs;
            pPosSolver = &Solver2D_Angs;

            break;

        case 2:
        Solver3D_Angs = QSolver3D_Angs(parrMeas, mvctMeas.size(),mtblEstPrfl
             ,mHBeaconEst,0.0009,arrXRez,arrSBeacon_GSK);

        pSolver2D_Angs = &Solver3D_Angs;
        pPosSolver = &Solver3D_Angs;

            break;

        default:

            break;
        }
    delete []parrMeas;


    mvctMeas = pPosSolver->mVectBigMeasures;

       double arrMtrx_dFGrA_po_dX_Inv[9] = {0.};
        if (mALGOR_TYPE == PEREBOR)
        {
            THREAD_RESULT = estimateAngs(pSolver2D_Angs,  &arrX0[5]
            , numPart, mMaxQuantIter_USBL,&arrXRez[5]
            , valNeviaz0, arrMtrx_dFGrA_po_dX_Inv);
        }
        if (mALGOR_TYPE == NUTON)
        {

           THREAD_RESULT = estimateParams(pSolver2D_Angs,  &arrX0[5]
           ,numPart ,mMaxQuantIter_USBL, mUSBLcoeff, &arrXRez[5]
            , valNeviaz0, arrMtrx_dFGrA_po_dX_Inv);
        }


        QVector <double>vct_arr_dAngs_po_dSgps( pSolver2D_Angs->mDimX *3);
        QVector <double>vct_arr_dAngs_po_dH( pSolver2D_Angs->mDimX);
        QVector <double>vct_arr_K_dqe((pSolver2D_Angs->mDimX) *(pSolver2D_Angs->mDimX));
        QVector <double>vct_arr_dAngs_po_dAngSins((pSolver2D_Angs->mDimX) );

        QVector <double>vct_arr_dAngs_po_dProfile((pSolver2D_Angs->mDimX) );

        QVector <double>vct_arr_AngKgps((pSolver2D_Angs->mDimX) * (pSolver2D_Angs->mDimX) );

        QVector <double>vct_arr_AngKsinsQ((pSolver2D_Angs->mDimX) * (pSolver2D_Angs->mDimX) );

        QVector <double>vct_arr_AngKgps_Psi_Tetta_per_1_rad((pSolver2D_Angs->mDimX) * (pSolver2D_Angs->mDimX) );

        QVector <double>vct_arr_arr_AngKt_per_1sec((pSolver2D_Angs->mDimX) * (pSolver2D_Angs->mDimX) );



     emit USBLSucceeded(  true
                           ,vct_arr_dAngs_po_dSgps
                           ,vct_arr_dAngs_po_dH
                           ,vct_arr_K_dqe
                           ,vct_arr_dAngs_po_dAngSins
                           ,vct_arr_dAngs_po_dProfile
                           ,vct_arr_AngKgps
                           ,vct_arr_AngKsinsQ
                           ,vct_arr_AngKgps_Psi_Tetta_per_1_rad
                           ,vct_arr_arr_AngKt_per_1sec
                           );
        // 1
    QVector <double>vct_dAngs_po_dSgps(9);
    vct_dAngs_po_dSgps.fill(0.);
 pSolver2D_Angs->calc_dAngs_po_dSgps(&arrXRez[5],arrMtrx_dFGrA_po_dX_Inv, arr_dX_po_dSgps,&vct_dAngs_po_dSgps);//
   qDebug() << "1. calc_dAngs_po_dSgps done" << vct_dAngs_po_dSgps.at(0);
       // 2
    QVector <double>vct_dAngs_po_dH(3);
    vct_dAngs_po_dH.fill(0.);
 pSolver2D_Angs->calc_dAngs_po_dHBeacon(&arrXRez[5],arrMtrx_dFGrA_po_dX_Inv, arr_dX_po_dHBeacon,vct_dAngs_po_dH);//
   qDebug() << "2. calc_dAngs_po_dHBeacon done";
    // 3
    QVector <double>vct_K_dqe(9) ;
    vct_K_dqe.fill(0.);
  pSolver2D_Angs->calc_K_qe(&arrXRez[5],arrMtrx_dFGrA_po_dX_Inv, vct_K_dqe);//- !!
    qDebug() << "3. calc_K_qe done";
    // 4
    QVector <double>vct_dAngs_po_dAngSins(9) ;//
    vct_dAngs_po_dAngSins.fill(0.);
  pSolver2D_Angs->calc_dAngs_po_dAngSins(&arrXRez[5],arrMtrx_dFGrA_po_dX_Inv, arr_dX_po_AngSins, vct_dAngs_po_dAngSins);//+
    qDebug() << "4. calc_dAngs_po_dAngSins done";
    // 5
    QVector <double>vct_dAngs_po_dProfile(3);
    vct_dAngs_po_dProfile.fill(0.);
  pSolver2D_Angs->calc_dAngs_po_dProfile(&arrXRez[5],arrMtrx_dFGrA_po_dX_Inv, arr_dX_po_dProfile,vct_dAngs_po_dProfile);//
    qDebug() << "5. calc_dAngs_po_dProfile done";
    // 6
    QVector <double>vct_AngKgps(9);//
    vct_AngKgps.fill(0.);
   pSolver2D_Angs->calc_AngKgps(&arrXRez[5],arrMtrx_dFGrA_po_dX_Inv,arrMtrx_dFGr_po_dX_Inv, vct_AngKgps);//
    qDebug() << "6. calc_AngKgps done";

// 7
    QVector <double>vct_AngKsinsQ_per_1_rad(9,0.);
  pSolver2D_Angs->calc_AngKsinsQ(&arrXRez[5],arrMtrx_dFGrA_po_dX_Inv,arrMtrx_dFGr_po_dX_Inv, vct_AngKsinsQ_per_1_rad);//++

  qDebug() << "7. calc_AngKsinsQ";
// 8
    QVector <double>vct_AngKsins_Psi_Tetta_per_1_rad(9,0.);
 pSolver2D_Angs->calc_AngKsins_Psi_Tetta(&arrXRez[5],arrMtrx_dFGrA_po_dX_Inv,arrMtrx_dFGr_po_dX_Inv, vct_AngKsins_Psi_Tetta_per_1_rad);//
   qDebug() << "8. calc_AngKsins_Psi_Tetta";

    // 9
    QVector <double>vct_AngKt_per_1sec(9,0.);
  pSolver2D_Angs->calc_AngKt(&arrXRez[5],arrMtrx_dFGrA_po_dX_Inv ,arrMtrx_dFGr_po_dX_Inv,  vct_AngKt_per_1sec);
    qDebug() << "9. calc_AngKt";
    emit USBLSucceeded(  false
                          ,vct_dAngs_po_dSgps
                          ,vct_dAngs_po_dH
                          ,vct_K_dqe
                          ,vct_dAngs_po_dAngSins
                          ,vct_dAngs_po_dProfile
                          ,vct_AngKgps
                          ,vct_AngKsinsQ_per_1_rad
                          ,vct_AngKsins_Psi_Tetta_per_1_rad
                          ,vct_AngKt_per_1sec
                          );


    // отправляем в интерфейс результат выполнения задачи
    emit calcfinished(THREAD_RESULT);
}


//--------------------------------
void MyThread::setParams(QVector <QBigMeasure> vctMeas, TYPE_OF_SOLVER_TASK TypeOfSolverTask
                        ,QVector <double> vctX0,TTable_1D tblEstPrfl, double *arrGPSAprioriPosParams
                       ,double HBeaconEst,int quantIter_LBL, int quantIter_USBL
                         ,double Gape,double AngStep, ALGOR_TYPE nALGOR_TYPE,double LBLcoeff,double USBLcoeff)
{
   mvctMeas = vctMeas;
   mTypeOfSolverTask = TypeOfSolverTask;
   mvctX0 = vctX0;
   mtblEstPrfl = tblEstPrfl;
   memcpy(marrGPSAprioriPosParams, arrGPSAprioriPosParams, 3 * sizeof(double));
   mHBeaconEst = HBeaconEst;

   mMaxQuantIter_LBL = quantIter_LBL;
   mMaxQuantIter_USBL = quantIter_USBL;
   mGape =Gape;
   mAngStep = AngStep;
   mALGOR_TYPE = nALGOR_TYPE;
   mLBLcoeff = LBLcoeff;
   mUSBLcoeff = USBLcoeff;
}
//--------------------------------

// вычисление неизвестного вектора позиционирования методом наименьших квадратов
//arrX[0],arrX[1],arrX[2] - вектор параллакса антенны ПСК
//arrX[3],arrX[4],arrX[5] - вектор координат маяка в ГСК
//arrX[6],arrX[7],arrX[8] - вектор углов ориентации антенны в ПСК
// INPUT:
//arrX0[mDimX] - начальное приближение
// OUTPUT:
//arrXEst[mDimX] - решение
// valNeviaz0 - невязка
MYTHREAD_RESULT MyThread::estimateParams(QPosSolver *pPosSolver,  double *arrX0
     , int numPart,const int MaxQuantIter, double coeff,double *arrXRez
     , double &valNeviaz0, double *arrMtrx_dFGr_po_dX_Inv)
{

    int mDimX = pPosSolver->mDimX;
    int mDimY = pPosSolver->mDimY;
    //0.

   memset(arrMtrx_dFGr_po_dX_Inv, 0, mDimX * mDimX * sizeof(double));

   double *arrMean = new double [mDimY];
   double *arrDisp = new double [mDimY];
   double *arrNeviazSquare = new double [mDimY];


   // 1. формирование вектора начальгного приближения
   memcpy(arrXRez, arrX0, mDimX * sizeof(double));

   ///

   // 2. вычисление начальной невязки
   int quantGoodZamers0 = 0;
   valNeviaz0 = pPosSolver->calcArrNeviaz_Mean_and_Disp(arrXRez, quantGoodZamers0
                                                , arrMean, arrDisp, arrNeviazSquare);
   QVector<double> vctMean(mDimY);
   QVector<double> vctDisp(mDimY);
   QVector<double> vctNeviazSquare(mDimY);
   for (int i =0; i < mDimY;++i)
   {
     vctMean.replace(i,arrMean[i]) ;
     vctDisp.replace(i,arrDisp[i]) ;
     vctNeviazSquare.replace(i,arrNeviazSquare[i]) ;
   }
   ///

   // 3. итреац процесс
  double *arrXCur = new double[mDimX];
  memcpy(arrXCur, arrXRez, mDimX * sizeof(double));

  double valNeviaz1=0.;
  double *arrDel = new double[mDimX];

  memset(arrDel, 0, mDimX* sizeof(double));

  bool breturn = false;

   for (int i =0; i < MaxQuantIter; ++i)
   {
       if(!(pPosSolver->doOneIteration(arrXRez, arrDel,arrMtrx_dFGr_po_dX_Inv)))
       {
       delete []arrXCur;
       delete []arrDel;
       return MYTHREAD_RESULT_ERROR;

       }

       double del1 = NormVect(arrDel, mDimX);

       MatrxMultScalar( arrDel, mDimX, 1, coeff, arrDel);

       MtrxMinusMatrx(arrXRez, arrDel,mDimX, 1, arrXCur);

       memcpy(arrXRez, arrXCur, mDimX * sizeof(double));


      pPosSolver->calcArrNeviaz_Mean_and_Disp(arrXRez, quantGoodZamers0
                                                   , arrMean, arrDisp, arrNeviazSquare);

      for (int j = 0; j < mDimX; ++j)
      {
          mvctX0.replace((numPart -1) *5 +j, arrXRez[j]);

      }
      for (int j =0; j < mDimY;++j)
      {
        vctMean.replace(j,arrMean[j]) ;
        vctDisp.replace(j,arrDisp[j]) ;
        vctNeviazSquare.replace(j,arrNeviazSquare[j]) ;
      }
       if (del1 < pPosSolver->mToler)
       {
           delete []arrXCur;
           delete []arrDel;
           delete[]arrMean ;
           delete[]arrDisp ;
           delete[]arrNeviazSquare ;

        emit progress(i , vctNeviazSquare,vctMean, vctDisp,numPart, true,mvctX0);

       return MYTHREAD_RESULT_OK;

       }

      if(userBreakRequested)
      {
          delete []arrXCur;
          delete []arrDel;
          delete[]arrMean ;
          delete[]arrDisp ;
          delete[]arrNeviazSquare ;
          return MYTHREAD_RESULT_USER_CANCEL;
      }

       emit progress(i , vctNeviazSquare,vctMean, vctDisp,numPart, true,mvctX0);
   }



   delete []arrXCur;
   delete []arrDel;
   delete[]arrMean ;
   delete[]arrDisp ;
   delete[]arrNeviazSquare ;

    return MYTHREAD_RESULT_MAX_STEP;
}

//--------------------------------------
// вычисление неизвестного вектора позиционирования методом наименьших квадратов
//arrX[0],arrX[1],arrX[2] - вектор параллакса антенны ПСК
//arrX[3],arrX[4],arrX[5] - вектор координат маяка в ГСК
//arrX[6],arrX[7],arrX[8] - вектор углов ориентации антенны в ПСК
// INPUT:
//arrX0[mDimX] - начальное приближение
// OUTPUT:
//arrXEst[mDimX] - решение
// valNeviaz0 - невязка
MYTHREAD_RESULT MyThread::estimateAngs(QSolver2D_Angs *pPosSolver,  double *arrX0
     , int numPart,const int MaxQuantIter,double *arrXRez
     , double &valNeviaz0, double *arrMtrx_dFGr_po_dX_Inv)
{


    int mDimX = pPosSolver->mDimX;
    int mDimY = pPosSolver->mDimY;

    double *arrMean = new double [mDimY];
    double *arrDisp = new double [mDimY];
    double *arrNeviazSquare = new double [mDimY];
    //0.

   memset(arrMtrx_dFGr_po_dX_Inv, 0, mDimX * mDimX * sizeof(double));


   // 1. формирование вектора начальгного приближения
   memcpy(arrXRez, arrX0, mDimX * sizeof(double));

   ///

   // 2. вычисление начальной невязки
   int quantGoodZamers0 = 0;
   valNeviaz0 = pPosSolver->calcArrNeviaz_Mean_and_Disp(arrXRez, quantGoodZamers0
                                                        , arrMean, arrDisp, arrNeviazSquare);
   QVector<double> vctMean(mDimY);
   QVector<double> vctDisp(mDimY);
   QVector<double> vctNeviazSquare(mDimY);
   for (int i =0; i < mDimY;++i)
   {
     vctMean.replace(i,arrMean[i]) ;
     vctDisp.replace(i,arrDisp[i]) ;
     vctNeviazSquare.replace(i,arrNeviazSquare[i]) ;
   }
   ///
   // перебор 3-ной цикл
   double arg_step = mAngStep/ 180. * M_PI;
   double valGape = mGape/180. * M_PI;
   double in_v_storonu = valGape/arg_step;
   int in = 2* in_v_storonu +1;
   double bet_min = arrX0[0] - ((double)in_v_storonu) * arg_step;
   double eps_min = arrX0[1] - ((double)in_v_storonu) * arg_step;
   double alf_min = arrX0[2] - ((double)in_v_storonu) * arg_step;

   int nSt = in * in * in;

   double valNevCur = 0.;

   int iChetchik = 0;
   for (int i =0; i < in; ++i)
       for (int j =0; j < in; ++j)
           for (int k =0; k < in; ++k)
   {
      double arrXcur[3] = {0.};
      arrXcur[0] = bet_min + ((double)i) * arg_step;
      arrXcur[1] = eps_min + ((double)j) * arg_step;
      arrXcur[2] = alf_min + ((double)k) * arg_step;

      valNevCur = pPosSolver->calcArrNeviaz_Mean_and_Disp(arrXRez, quantGoodZamers0
                                                          , arrMean, arrDisp, arrNeviazSquare);
      if (valNevCur < valNeviaz0)
      {
         valNeviaz0 =  valNevCur;
         memcpy(arrXRez, arrXcur,3 * sizeof(double));
         for (int i = 0; i < mDimX; ++i)
         {
             mvctX0.replace((numPart -1) *5 +i, arrXRez[i]);

         }
         for (int i =0; i < mDimY;++i)
         {
           vctMean.replace(i,arrMean[i]) ;
           vctDisp.replace(i,arrDisp[i]) ;
           vctNeviazSquare.replace(i,arrNeviazSquare[i]) ;
         }

      }
      ++iChetchik;
      if ((iChetchik % 50) == 0)
      {
         emit progress(iChetchik , vctNeviazSquare,vctMean, vctDisp,numPart, true,mvctX0);

          if(userBreakRequested)
          {
              delete[]arrMean ;
              delete[]arrDisp ;
              delete[]arrNeviazSquare ;

              return MYTHREAD_RESULT_USER_CANCEL;
          }
      }

   }


   delete[]arrMean ;
   delete[]arrDisp ;
   delete[]arrNeviazSquare ;

emit progress(nSt, vctNeviazSquare,vctMean, vctDisp,numPart, true,mvctX0);
return MYTHREAD_RESULT_OK;

}



