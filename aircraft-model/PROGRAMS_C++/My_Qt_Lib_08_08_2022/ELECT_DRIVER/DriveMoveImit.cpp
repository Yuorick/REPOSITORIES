#include "DriveMoveImit.h"
#include <string.h>
#include <math.h>
#include "Ctrl.h"
#include "MatrixProccess.h"
#include "ElectMotor.h"

#include "LinDiffEq.h"

extern const int QUantColsCSVReport = 14;
extern double valOmReal =0.;

extern  bool BEZ_SHUMOV;

extern double aRRTruePhVect[QVARS] = {0.}; // отладка
QDriveMoveImit::QDriveMoveImit()
{
 mpCtrl = NULL;
 mRealDriver = QElectDriver();
 mMeasurmentImitator = QMeasurmentImitator();
 mEnvironment = TEnvironment();
 mpTargTraj = NULL;
 mMomOutReal = 0.;
 mMomOutModel = 0.;

}
//-------------------------------------
// конструктор копирования
 QDriveMoveImit :: QDriveMoveImit (const  QDriveMoveImit &R)
 {
         mpCtrl  = R.mpCtrl ;
         mRealDriver = R.mRealDriver;
         mMeasurmentImitator = R.mMeasurmentImitator;
         mEnvironment = R.mEnvironment;
         mpTargTraj = R.mpTargTraj;
         mMomOutReal = R. mMomOutReal;
         mMomOutModel = R.mMomOutModel;

 }
//-------------------------------------------------------------------
 // оператор присваивания
  QDriveMoveImit  &QDriveMoveImit::operator=( const QDriveMoveImit &R)
 {
      mpCtrl  = R.mpCtrl ;
      mRealDriver = R.mRealDriver;
      mMeasurmentImitator = R.mMeasurmentImitator;
      mEnvironment = R.mEnvironment;
      mpTargTraj = R.mpTargTraj;
      mMomOutReal = R. mMomOutReal;
      mMomOutModel = R.mMomOutModel;

     return *this ;
 }
  //-------------------------------------------------------------------
    // парам конструктор 1
   QDriveMoveImit:: QDriveMoveImit ( QCtrl *pCtrl, const QElectDriver RealDriver,
     const QMeasurmentImitator MeasurmentImitator,   QAbstractTraj  *pTargTraj
     ,const double MomOutReal,const double MomOutModel)
   {
      mpCtrl = pCtrl;
      mRealDriver = RealDriver;
      mMeasurmentImitator = MeasurmentImitator;
      mEnvironment = TEnvironment();
      mpTargTraj = pTargTraj;
      mMomOutReal = MomOutReal;
      mMomOutModel = MomOutModel;

   }

   //-------------------------------------------------------------------
     // парам конструктор 2
    QDriveMoveImit:: QDriveMoveImit ( QCtrl *pCtrl, const QElectDriver RealDriver,
      const QMeasurmentImitator MeasurmentImitator, const TEnvironment Environment
       ,QAbstractTraj *pTargTraj,const double MomOutReal,const double MomOutModel)
    {
       mpCtrl = pCtrl;
       mRealDriver = RealDriver;
       mMeasurmentImitator = MeasurmentImitator;
       mEnvironment = Environment;
       mpTargTraj = pTargTraj;
       mMomOutReal = MomOutReal;
       mMomOutModel = MomOutModel;


    }

 //-----------------------------------------------------------------------------------------------------
 // передвижение класса в заданную точку
 // INPUT:
 // VAlOmegaStat- заданная скорость
 // VAlIntegrStep - шаг интегрирования
 // VAlU0 - максимальная величина модуля U_qu
 // OUTPUT:
 // arrPhVectEnd[QVARS] - конечный фазовый вектор
 // содержит 7 столбцов - время, фазовый вектор, 2 управления, производную угловой скорости.
 void QDriveMoveImit:: move(double *arrQObjective0,const double VAlTBegin,const double VAlMovingT,const double VAlIntegrStep0
                  ,double *arrBuff, int *piQuantRows)
 {
     double arrQObjectiveCur[2] = {0.};
     arrQObjectiveCur[1] = arrQObjective0[1];
 int iQuantSteps = int((VAlMovingT - VAlTBegin)/VAlIntegrStep0);

 double   arrU[2] = {0.}, arrPrevU[2] = {0.};
 mpCtrl->mT0 = VAlTBegin;
 mpCtrl->mTCur = VAlTBegin;

// mpCtrl->mFiltr.init(valTcur, 1./20000.,mpCtrl->marrPhVect0);
 // истинный фазовый вектор системы текущий
 double arrTruePhVect[QVARS] = {0.};
 memcpy(arrTruePhVect, mpCtrl->mFiltr.marrCurEst, QVARS * sizeof(double));
 //arrTruePhVect[3] += mpCtrl->marrObjective[0];

 if(BEZ_SHUMOV)
 {
  memcpy(aRRTruePhVect, arrTruePhVect,  QVARS * sizeof(double));
 }
 // интеграторы dx = -Alf*x +u
 // Alf = 1 /T
// double val_Tet_Itegr   = mpCtrl->mFiltr.marrCurEst[3];
 //double val_Om_Itegr   = mpCtrl->mFiltr.marrCurEst[0];
 //double val_Id_Itegr   = mpCtrl->mFiltr.marrCurEst[1];
 //double val_Iqu_Itegr   = mpCtrl->mFiltr.marrCurEst[2];
 //double alf = 1. / VAlIntegratorTime;
 //double valTIntegrCur = mpCtrl->mTCur;

 //
 double arrKy[9] = {0.};
 QDriverMeasure measureCur;
 mMeasurmentImitator.createMeasure(arrTruePhVect, &measureCur);

mpTargTraj->recaclDecartCoord(mpCtrl->mTCur);
// нахождение углового положения цели и угловой скорости по КУ относительно привода
//double arrObjective[2] = {0.}; // угл положегние и скорость
//mpTargTraj->calcQ_And_Omega_z(arrObjective);
///
arrQObjectiveCur[0] = arrQObjective0[0] + arrQObjective0[1] *(mpCtrl->mTCur -mpCtrl->mT0);
//mpCtrl->processMeasure_and_calcCurU(arrObjective, mpCtrl->mTCur,measureCur,mpCtrl->mTCur,arrPrevU,arrU);
 mpCtrl->processMeasure_and_calcCurU(arrQObjectiveCur, mpCtrl->mTCur, measureCur
                                     ,mpCtrl->mTCur,arrPrevU,mMomOutModel,arrU);
 //
 if(arrBuff != NULL)
  {
      double *p = arrBuff;
      p[0] = mpCtrl->mTCur;
      memcpy(&p[1],arrTruePhVect, 4 * sizeof(double));
      memcpy(&p[5],arrU, 2 * sizeof(double));
      memcpy(&p[8],mpCtrl->mFiltr.marrCurEst, QVARS * sizeof(double));
      memcpy(&p[12],arrQObjectiveCur, 2 * sizeof(double));
      (*piQuantRows)++;
  }

TEnvironment Environment0(0.,0.,0.);

 for (int i = 1; i < iQuantSteps ; ++i )
 {

     mpCtrl->mTCur = mpCtrl->mTCur + VAlIntegrStep0;
     double arrPhVectEnd[QVARS] = {0.};
     mRealDriver.dragPhVect(mEnvironment,arrTruePhVect, arrU
                            , VAlIntegrStep0, VAlIntegrStep0, arrPhVectEnd,mMomOutReal);
     memcpy(arrTruePhVect,arrPhVectEnd, QVARS * sizeof(double));    

     double arrEstPhVectExtrp[QVARS] = {0.};

     mpCtrl->dragPhVect(Environment0,mpCtrl->mFiltr.marrCurEst,arrU
                        , VAlIntegrStep0, VAlIntegrStep0, arrEstPhVectExtrp,mMomOutModel);
     memcpy(mpCtrl->mFiltr.marrCurEst,arrEstPhVectExtrp, QVARS * sizeof(double));
     if(mpCtrl->mTCur >= (mpCtrl->mFiltr.mTCur + mpCtrl->mFiltr.mh))
     {

         if (mpCtrl->mTCur >= 20)
         {
             int uuui = 0;
         }
   mRealDriver.dragPhVect(mEnvironment,arrTruePhVect, arrU, mpCtrl->mTCur - (mpCtrl->mFiltr.mTCur+ mpCtrl->mFiltr.mh)
       , VAlIntegrStep0, arrPhVectEnd, mMomOutReal);
         memcpy(arrTruePhVect,arrPhVectEnd, QVARS * sizeof(double));
  mpCtrl->dragPhVect(Environment0,mpCtrl->mFiltr.marrCurEst,arrU, mpCtrl->mTCur - (mpCtrl->mFiltr.mTCur+ mpCtrl->mFiltr.mh)
       , VAlIntegrStep0, arrEstPhVectExtrp, mMomOutModel);
         memcpy(mpCtrl->mFiltr.marrCurEst,arrEstPhVectExtrp, QVARS * sizeof(double));
         mMeasurmentImitator.createMeasure(arrTruePhVect, &measureCur);


         arrPrevU[0]  = arrU[0];
         arrPrevU[1]  = arrU[1];
         //
         if(BEZ_SHUMOV)
         {
           memcpy(aRRTruePhVect, arrTruePhVect, QVARS * sizeof(double));// отладка
         }


         // нахождение углового положения цели и угловой скорости по КУ относительно привода
         mpTargTraj->recaclDecartCoord(mpCtrl->mTCur);
        // double arrObjective[2] = {0.}; // угл положегние и скорость
        // mpTargTraj->calcQ_And_Omega_z(arrObjective);

         ///
         arrQObjectiveCur[0] = arrQObjective0[0] + arrQObjective0[1] *(mpCtrl->mTCur -mpCtrl->mT0);
      mpCtrl->processMeasure_and_calcCurU(arrQObjectiveCur, mpCtrl->mTCur
        , measureCur,mpCtrl->mTCur,arrPrevU,mMomOutModel,arrU);

         if(arrBuff != NULL)
         {
             double arrFx[QVARS];
             mRealDriver.calcFx( mEnvironment, arrTruePhVect, arrFx,mMomOutReal);
             double *p =&(arrBuff[((*piQuantRows)) * QUantColsCSVReport]);
             p[0] = mpCtrl->mTCur;
             memcpy(&p[1],arrTruePhVect, 4 * sizeof(double));
             memcpy(&p[5],arrU, 2 * sizeof(double));
             p[7] = arrFx[0];
             memcpy(&p[8],mpCtrl->mFiltr.marrCurEst, QVARS * sizeof(double));
             memcpy(&p[12],arrQObjectiveCur, 2 * sizeof(double));

             (*piQuantRows)++;
         }

         // пересчет интеграторов
      /*   double arrTargPhVect[QVARS] = {0.};
         memcpy(arrTargPhVect, mpCtrl->mStatSolutionParams.marrStatPhVect, QVARS * sizeof(double));
         arrTargPhVect[3] = mpCtrl->marrObjective[0] + (mpCtrl->mTCur - mpCtrl->mT0) * mpCtrl->marrObjective[1];
         for (int j =0; j < QVARS ; ++j)
         {
            double valDelta = mpCtrl->mFiltr.marrCurEst[j] - arrTargPhVect[j];
            mpCtrl->marrIntegrators[j] = mpCtrl->marrIntegrators[j]
         +   (- alf *  mpCtrl->marrIntegrators[j]   + valDelta) *  (mpCtrl->mTCur -valTIntegrCur) ;
        }

         valTIntegrCur = mpCtrl->mTCur;
        */

     }

 }

 }

 //---------------------------------------------
 // анализ разбросов траектории на последнем интервале времени
  // [mMovingT - TPeriod; mMovingT]
 //INPUT:
 //VAl_TPeriod - ПЕРИОД гармонической составляющей остат момента
 // arrBuff[quantDoneSteps * QUantColsCSVReport] - буфер с выходной и нформацией имитац модели

  //OUTPUT:
 // arrSyst[2] - вектор систематической ошибки по угловому положению и угл скоргости
 // arrDisp[2] - вектор диспесий M(x-xsyst)*(x-xsyst)

 void QDriveMoveImit::StatisticProcess(double *arrBuff, const int QUantRows,const double VAl_TPeriod, double *arrSyst, double *arrDisp)
 {

     // вычисление среднего вектора
     double arrE[2] = {0.};
     memset (arrSyst, 0,2 * sizeof(double));
     memset (arrDisp, 0,2 * sizeof(double));
     double arrT0[2] = {0.}, arrT1[2] = {0.}, arrT2[2] = {0.};
     int quantPoints = 0;
     const double VAlTEnd = arrBuff[QUantRows * QUantColsCSVReport -  QUantColsCSVReport];
     for (int i = 0; i < QUantRows; ++i)
     {
         double *p = &(arrBuff[QUantRows * QUantColsCSVReport - (i + 1) * QUantColsCSVReport]);

       if (p[0] < (VAlTEnd - VAl_TPeriod) )
       {
           break;
       }

       arrT0[0] = p[4]- p[12] ; // разн угл положение
       arrT0[1] = p[1]- p[13] ;  // разн угл скорость




      // arrT0[3] -= (mpCtrl->marrObjective[0] + p[0] * mpCtrl->marrObjective[1]);
       MtrxSumMatrx(arrE, arrT0, 1, 2, arrT1) ;

       memcpy(arrE, arrT1,2 * sizeof(double));

       for (int j =0; j < 2; ++j)
       {
         arrT2[j] =  arrT0[j] * arrT0[j];

       }

       MtrxSumMatrx(arrDisp,arrT2, 1, 2, arrT0) ;
       memcpy(arrDisp, arrT0,2 * sizeof(double));
       quantPoints++;

     }
     ///
     MatrxMultScalar(arrE, 1, 2, 1./ (double(quantPoints)),arrE);
     MatrxMultScalar(arrDisp, 1, 2, 1./ (double(quantPoints)),arrDisp);
     for (int j =0; j < 2 ; ++j)
     {
       arrDisp[j] = arrDisp[j] - arrE[j] * arrE[j];
     }


     arrSyst[0] = arrE[0];
     arrSyst[1] = arrE[1];
    // MtrxMinusMatrx(arrE, mpCtrl->mStatSolutionParams.marrStatPhVect,1, QVARS-1, arrSyst);
     //arrSyst[3] = arrE[3] ;

 // второй способ
    memset (arrDisp, 0,2 * sizeof(double));
     for (int i = 0; i < QUantRows; ++i)
     {
         double *p = &(arrBuff[QUantRows * QUantColsCSVReport - (i + 1) * QUantColsCSVReport]);

       if (p[0] < (VAlTEnd - VAl_TPeriod) )
       {
           break;
       }

       arrT0[0] = p[4]- p[12] ; // разн угл положение
       arrT0[1] = p[1]- p[13] ;  // разн угл скорость
       MtrxMinusMatrx(arrT0, arrE,1, 2, arrT1) ;


       for (int j =0; j < 2; ++j)
       {
         arrT2[j] =  arrT1[j] * arrT1[j];

       }

       MtrxSumMatrx(arrDisp,arrT2, 1, 2, arrT0) ;
       memcpy(arrDisp, arrT0,2 * sizeof(double));


     }
 MatrxMultScalar(arrDisp, 1, 2, 1./ (double(quantPoints)),arrDisp);


 }
 //----------------------------------------------------
 // вычисление разбросов траектории (точности) аналитическим методом
 // INPUT:
 // arrW[2] - вектор добавочных управлений
 // OUTPUT:
 // arr_dx_po_dAlf[QVARS * QPARAMS] - матрица частных производных фазового вектора по вектору параметров
 // arrAHarmSquare[QVARS ] - векторов квадратов модулей амплитуд гармоник фазового вектора
 // arrSystErr[QVARS ] - вектор сиситематических ошибок фазового вектора
 // arrSumSig[QVARS ] - результирующий вектор СКЗ отклонений фазового вектора
  void QDriveMoveImit::calcAnaliticSumTrajScatters(double *arrObjective, double *arrW,double *arr_dx_po_dAlf, double *arrAHarm
                                           , double *arrSystErr,double *arrRandomSig, double *arrSumSig)
  {
   calcSystTrajScatters(arrObjective,arrW,arr_dx_po_dAlf,  arrSystErr) ;
    calcHarmAmpVect(arrObjective,arrW, arrAHarm );

    double arrFluktK[QVARS * QVARS ] ={0.};
    calcFluctK(arrW,  arrFluktK);
    for (int i = 0; i < QVARS; ++i)
    {
      arrRandomSig[i] =  sqrt(arrAHarm[i] * arrAHarm[i]/2. + arrFluktK[ i *QVARS + i] );
      arrSumSig[i] = sqrt(arrSystErr[i] *arrSystErr[i] + arrAHarm[i] * arrAHarm[i]/2. + arrFluktK[ i *QVARS + i] );
    }
  }

  //---------------------------------
    // вычисление систематических разбросов траектории (точности) аналитическим методом
    //INPUT:
    // arrW[2] - вектор добавочных управлений
    // OUTPUT:
    // arr_dx_po_dAlf[QVARS * QPARAMS] - матрица частных производных фазового вектора по вектору параметров
    // arrSystErr[QVARS ] - вектор сиситематических ошибок фазового вектора

    void QDriveMoveImit::calcSystTrajScatters(double *arrObjective,double *arrW,double *arr_dx_po_dAlf,  double *arrSystErr)
    {
        // вычисление нового стайионарного решения
        double arrStablePhVect[QVARS ] ={0.};
        QStatSolutionParams StatSolutionParams;
        mpCtrl->calcStationaryParams(arrObjective, mpCtrl->mCmpArrLamb, &StatSolutionParams,mMomOutModel);
        const double valUd = StatSolutionParams.marrStatU[0] +arrW[0];
        const double valUq = StatSolutionParams.marrStatU[1] +arrW[1];
        mpCtrl->calcStableSolution(arrObjective, valUd,  valUq, arrStablePhVect,mMomOutModel);

        double arr_dF_po_dx[QVARS *QVARS] = {0.},arr_dF_po_dW[2 * QVARS] = {0.}
                 , arrT[QVARS * QVARS] = {0.}, arrCC[2 *QVARS] ={0.}, arrA[QVARS *QVARS] = {0.}, arrAInv[QVARS *QVARS] = {0.};

        mpCtrl->fill_df_po_px_and_mtrxB_(arrStablePhVect
                               ,arr_dF_po_dx,arr_dF_po_dW);
        memcpy(arrCC, StatSolutionParams.marrGears, 2 * QVARS * sizeof(double));
        MtrxMultMatrx(arr_dF_po_dW,QVARS, 2, arrCC, QVARS, arrT) ;
        MtrxSumMatrx(arrT, arr_dF_po_dx, QVARS, QVARS, arrA) ;

        bool brez = InverseMtrx(arrA, QVARS, arrAInv) ;


        double arr_dF_po_dAlf[QVARS * QPARAMS] = {0.};
        mpCtrl->calc_dF_po_dAlf(arrStablePhVect, arr_dF_po_dAlf,mMomOutModel);

        MtrxMultMatrx(arrAInv,QVARS, QVARS, arr_dF_po_dAlf,QPARAMS, arr_dx_po_dAlf);
        MatrxMultScalar(arr_dx_po_dAlf, QVARS, QPARAMS, -1., arr_dx_po_dAlf);

        double arrParamsDiffer[QPARAMS] = {0.};
        createVectParamsDiffers(arrParamsDiffer);
        arrParamsDiffer[6] -= mpCtrl->mElectMotor.getMomDryFriction ();
                ///
         MtrxMultMatrx(arr_dx_po_dAlf,QVARS, QPARAMS, arrParamsDiffer,1, arrSystErr);



    }


    //---------------------------------
      // вычисление амплитуд разбросов траектории
    //, вызванных гармонической составляющей остатосчного момента
    //(точности) аналитическим методом
    //INPUT:
    // arrW[2] - вектор добавочных управлений
      // OUTPUT:
      // arrAHarmSquare[QVARS] - вектор  амлитуд
void QDriveMoveImit::calcHarmAmpVect(double *arrObjective,double *arrW, double *arrAHarm )
{
        // вычисление нового стационарного решения
          double arrStablePhVect[QVARS ] ={0.};
          QStatSolutionParams StatSolutionParams;
          mpCtrl->calcStationaryParams(arrObjective, mpCtrl->mCmpArrLamb, &StatSolutionParams,mMomOutModel);
          const double valUd = StatSolutionParams.marrStatU[0] +arrW[0];
          const double valUq = StatSolutionParams.marrStatU[1] +arrW[1];
          mpCtrl->calcStableSolution(arrObjective, valUd,  valUq, arrStablePhVect,mMomOutModel);

          double arr_dF_po_dx[QVARS *QVARS] = {0.},arr_dF_po_dW[2 * QVARS] = {0.}
                   , arrT[QVARS * QVARS] = {0.}, arrCC[2 *QVARS] ={0.}, arrA[QVARS *QVARS] = {0.}, arrAInv[QVARS *QVARS] = {0.};

          mpCtrl->fill_df_po_px_and_mtrxB_(arrStablePhVect
                                 ,arr_dF_po_dx,arr_dF_po_dW);
          memcpy(arrCC, StatSolutionParams.marrGears, 2 * QVARS * sizeof(double));
          MtrxMultMatrx(arr_dF_po_dW,QVARS, 2, arrCC, QVARS, arrT) ;
          MtrxSumMatrx(arrT, arr_dF_po_dx, QVARS, QVARS, arrA) ;
          TComp cmparrA[QVARS * QVARS], cmparrAInv[QVARS * QVARS];
          for (int i =0; i < QVARS * QVARS; ++i)
          {
           cmparrA[i] = TComp(arrA[i], 0.);
          }
          for (int i =0; i < QVARS ; ++i)
          {
             cmparrA[i * QVARS + i].m_Im +=
                     StatSolutionParams.marrStatPhVect[0] * mpCtrl->mElectMotor.mZp ;
          }


          bool brez = InverseMtrx(cmparrA, QVARS, cmparrAInv) ;

          TComp parrRez[QVARS * QVARS];
          double temp = (mpCtrl->mElectMotor.calc_AmpHarmMom (0., arrStablePhVect[0]))*mpCtrl->getInvSumJ0();
          TComp valScal(temp ,0.);
          MatrxMultScalar(cmparrAInv, QVARS, QVARS, valScal,parrRez);

          for (int i =0; i < QVARS; ++i)
          {
             arrAHarm[i] = parrRez[i * QVARS].modul();
          }

 }

  //-----------------------------
  void QDriveMoveImit::createVectParamsDiffers(double *arrParamsDiffer)
  {

      // mInvL - marrSpreadParams[0]
      // mResist  - marrSpreadParams[1]
      // mPsi_f   - marrSpreadParams[2]
      // InvJ        - marrSpreadParams[3]
      // Cx       -marrSpreadParams[4]
      // Cv       - marrSpreadParams[5]
      // MomOut       - marrSpreadParams[6]
      for(int i =0; i < QPARAMS; ++i)
      {
         double valReal =  mRealDriver.getParam(i);
         double valModel=  mpCtrl->getParam(i);
         arrParamsDiffer[i] = valReal - valModel;

      }
  }
//--------------------------------------------------------------
void QDriveMoveImit::calcFluctK(double *arrW,  double *arrFluktK)
  {
    /*  // вычисление нового стайионарного решения
      double arrStablePhVect[QVARS ] ={0.};
      const double valUd = mpCtrl->mStatSolutionParams.marrStatU[0] +arrW[0];
      const double valUq = mpCtrl->mStatSolutionParams.marrStatU[1] +arrW[1];
      mpCtrl->calcStableSolution( valUd,  valUq, arrStablePhVect);

      double arr_dF_po_dx[QVARS *QVARS] = {0.},arr_dF_po_dW[2 * QVARS] = {0.}
               , arrT[QVARS * QVARS] = {0.}, arrCC[2 *QVARS] ={0.}, arrA[QVARS *QVARS] = {0.}, arrAInv[QVARS *QVARS] = {0.};

      mpCtrl->fill_df_po_px_and_mtrxB_(arrStablePhVect
                             ,arr_dF_po_dx,arr_dF_po_dW);
      memcpy(arrCC, mpCtrl->mStatSolutionParams.marrGears, 2 * QVARS * sizeof(double));
      MtrxMultMatrx(arr_dF_po_dW,QVARS, 2, arrCC, QVARS, arrT) ;
      MtrxSumMatrx(arrT, arr_dF_po_dx, QVARS, QVARS, arrA) ;

      //вычисление интеграла
      memset(arrFluktK, 0,QVARS * QVARS * sizeof(double) );
      double arrKw[QVARS * QVARS] ={0.},arrT0[QVARS * QVARS] ={0.}
       ,arrT1[QVARS * QVARS] ={0.},arrT2[QVARS * QVARS] ={0.},arrKwDelT[QVARS * QVARS] ={0.};

      double valJInvSum = mpCtrl->getInvSumJ0();
      arrKw[0] = mpCtrl->mElectMotor.calcDisp_FluctMomNoise(arrStablePhVect[0])*valJInvSum*valJInvSum;
      double valTIntegr = 5.;
      double valStepIntegr = 1./5000.;
      int nst = valTIntegr/valStepIntegr;

      // матрица перехода

      double  arr_L [QVARS * QVARS] = {0.};
     MatrxMultScalar(arrA, QVARS, QVARS, valStepIntegr,arr_L);
     for(int i =0; i < QVARS; ++i)
     {
       arr_L[ QVARS * i + i] += 1.;
     }
     ///

     memcpy(arrT0, arr_L, QVARS * QVARS * sizeof(double) );
     MatrxMultScalar(arrKw, QVARS, QVARS, valStepIntegr,arrKwDelT);
      for (int i =0; i< nst; ++i)
      {
        MtrxMultMatrx_MultMatrxTransp(arr_L, arrFluktK, arrT0,QVARS, arrT1);
        MtrxSumMatrx(arrFluktK, arrKwDelT,QVARS, QVARS, arrT0) ;
        memcpy(arrFluktK, arrT0, QVARS * QVARS * sizeof(double) );

      }
*/
  }


