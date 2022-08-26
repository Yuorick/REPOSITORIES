#include "ParamsID.h"
#include <math.h>
#include <float.h>
#include "StatSolutionParams.h"
#include "DriveMoveImit.h"
#include "Comp.h"
#include "CtrlFollow.h"
#include "Ctrl.h"
#include "ElectMotor.h"
#include "Load.h"
#include  <string.h>
#include "MatrixProccess.h"
#include "MinSquare.h"
#include "CtrlTrad.h"

//extern const double constDbl_mh;
extern double constDbl_mh;
extern const int QUantColsCSVReport;

QParamsID::QParamsID()
{
    mDriverModel = QElectDriver  ();
    mDriverReal = QElectDriver  ();
    mMeasurmentImitator = QMeasurmentImitator();
    memset(marrSpreadParams, 0, QPARAMS* sizeof(double));
    // массив номеров идентифицируемых параметров
    memset(miarrTargNums, 0, (QPARAMS + 2) * sizeof(int));
    mlenTargNums = 0;

}

// конструктор копирования
 QParamsID ::QParamsID (const QParamsID &R)
 {
     mDriverModel = R.mDriverModel;
     mDriverReal = R.mDriverReal;
     mMeasurmentImitator = R.mMeasurmentImitator;
     memcpy(marrSpreadParams, R.marrSpreadParams, QPARAMS * sizeof(double));
     memcpy(miarrTargNums, R.miarrTargNums, (QPARAMS + 2) * sizeof(int));
     mlenTargNums = R. mlenTargNums;
 }
 // оператор присваивания
 QParamsID &QParamsID::operator=(const QParamsID  &R)
 {
     mDriverModel = R.mDriverModel;
     mDriverReal = R.mDriverReal;
     mMeasurmentImitator = R.mMeasurmentImitator;
     memcpy(marrSpreadParams, R.marrSpreadParams, QPARAMS * sizeof(double));
     memcpy(miarrTargNums, R.miarrTargNums, (QPARAMS + 2) * sizeof(int));
     mlenTargNums = R.mlenTargNums;

   return *this ;
 }
// парам констр
QParamsID :: QParamsID(  const QElectDriver DriverModel,const QElectDriver DriverReal
              , const QMeasurmentImitator MeasurmentImitator, const double *arrSpreadParams
              , const int *iarrTargNums , const double lenTargNums )
{
    mDriverModel = DriverModel;
    mDriverReal = DriverReal;
    mMeasurmentImitator = MeasurmentImitator;
    memcpy(marrSpreadParams, arrSpreadParams, QPARAMS * sizeof(double));
    memcpy(miarrTargNums, iarrTargNums, (QPARAMS + 2) * sizeof(int));
    mlenTargNums = lenTargNums;
}

//моделирование идентификации параметров привода методом стационарных решений
// DriverModel - это представление о приводе, заложенное в модель
// на основании моельных данных рассчитываются коэффициенты передачи (управления)
// DriverReal - это реальный привод, с параметрами, которые могут отличаться от заложенных в модели
// ОПИСАНИЕ МЕТОДА
// моделируются набор экспериментов, по выводу сиситемы на стационарный режим.
// описание одного эксперимента:
//  1.заданы начальные условия сиситемы
//  2.заданы целевые условия задачи слежения - Tetta0, Omega0
//  3. задан набор корней характеристического уравнения
//  4. производится моелирование выывода реального привода на стационарный режим
//  5.данные экспееримента апоминаются
//После того как проиведено моеделирование реультатов экспериментов,
// производится оцуенка реальных параметров привода методом наименьших квадратов
// INPUT:
// DriverModel - модельный привод
// DriverReal - реальный привод
// MeasurmentImitator - модель измерительной сиситемы
// arrSpreadParams - вектор СКЗ разбросов параметров привода (в расчета не участвует)
// deltaTetta - шаг по Tetta0
// quantTetta - к-во шагов по Tetta0
// deltaOmega - шаг по Omega0
// quantOm - к-во шагов по Omega0
// valMovingT - врмя вывода на стационарный режим
 // cmpArrLamb[QVARS * QUantLambSets] - массив наборов корней характеристич уравнения
 // QUantLambSets - к-во наборов корней
// замечание - общее к-во экспериметов равно quantTetta * quantOm *QUantLambSets
 // iarrTargNums[lenTargNums] - номера параметров сиситемы, по которым производится идентификация
 // lenTargNums - длина массива номеров параметров сиситемы, по которым производится идентификация
// замечание - вектор  параметров сиситемы ореленв классе QElectDriver
// в него входят, в частности дуктивность и соппроивление, с номерами 0 и 1 соответственно
// но производить идентфикацию этих параметов нет нобходимости,
// т.к. они могум быть измерены с очен выскй тчностью мульимтом ил RCL-метром
// OUTPUT:
// возвращеет оценочный привод
QElectDriver QParamsID::static_identify_by_StatMeth(const QElectDriver DriverModel,const QElectDriver DriverReal
    , const QMeasurmentImitator MeasurmentImitator, const double *arrSpreadParams
    ,const double DEltaTetta, const int QUantTetta,const double VAlOmega0,const double DEltaOmega, const int QUantOm
    ,const double  VAlMovingT,TComp *cmpArrLamb, const int QUantLambSets,int *iarrTargNums
    ,const int LEnTargNums, const double VAlIntegrStep )
{
    QStatSolutionParams *parrStatSolutionParams = new QStatSolutionParams[QUantLambSets * QUantTetta * QUantOm];
    double *arrDynamicDelX = new double[QUantLambSets * QUantTetta * QUantOm * QVARS];

    QParamsID ParamsID( DriverModel, DriverReal
                         , MeasurmentImitator,arrSpreadParams, iarrTargNums,LEnTargNums);

    ParamsID.imitateStatExperimentsRez(DEltaTetta, QUantTetta, VAlOmega0,DEltaOmega, QUantOm
                                       , VAlMovingT,cmpArrLamb, QUantLambSets, parrStatSolutionParams
                                       , arrDynamicDelX, VAlIntegrStep );
    double arrDElAlf[QPARAMS] = {0.};
    int quantExprs = QUantLambSets * QUantTetta * QUantOm;
    //ParamsID.calcDelAlf_StatMeth(quantExprs, parrStatSolutionParams, arrDynamicDelX,arrDElAlf);
    QElectMotor MotorRez;
    QLoad LoadRez;
    double valMomOut = -1.;


    double valNeviaz = ParamsID.identify_by_StatMeth(quantExprs, parrStatSolutionParams
             , arrDynamicDelX, &MotorRez, &LoadRez, &valMomOut);

    delete []parrStatSolutionParams;
    delete []arrDynamicDelX;
 QElectDriver DriverReturn = DriverModel;
 for (int i =0; i < LEnTargNums; ++i)
 {
     double valParam = DriverReturn.getParam(iarrTargNums[i]) + arrDElAlf[i];
     DriverReturn.setParam(iarrTargNums[i], valParam);
 }
  return DriverReturn ;
}

// имтация результатов экспериментов
// вывода сиситемы на стационарный режим

void QParamsID::imitateStatExperimentsRez(const double DEltaTetta, const int QUantTetta,const double VAlOmega0,const double DEltaOmega
    , const int QUantOm,const double  VAlMovingT, TComp *cmpArrLamb, const int QUantLambSets
      , QStatSolutionParams *parrStatSolutionParams,double *arrDynamicDelX, const double VAlIntegrStep )
{
     double arrObjective[2] = {0.};

     //начальнвая  угловая скорость вращения
     double valOmegaBegin = 0.;
     //  начальное угловое пложение
     double valTettaBegin =0.;




    for (int iz = 0; iz <QUantLambSets; iz++)
    {

        for (int itet = 0; itet< QUantTetta; ++itet)
        {
           arrObjective[0] = (((double)itet) + 1.) * DEltaTetta;
           for (int iOm = 0; iOm < QUantOm; iOm++ )
           {
             arrObjective[1] = VAlOmega0 + ((double)iOm)  * DEltaOmega;


             QCtrlFollow ModelCtrlFollow (mDriverModel.mElectMotor,mDriverModel.mpLoad
                             ,marrSpreadParams, mDriverModel.mMomOut, valTettaBegin, valOmegaBegin, 0.,arrObjective,constDbl_mh);
             ModelCtrlFollow.calcStationaryParams(arrObjective,&(cmpArrLamb[ QVARS * iz])
               , &(parrStatSolutionParams[iz * QUantTetta * QUantOm + itet * QUantOm + iOm]));
             ModelCtrlFollow.mStatSolutionParams = parrStatSolutionParams[iz * QUantTetta * QUantOm + itet * QUantOm + iOm];

             QCtrl *pCtrl = &ModelCtrlFollow;
             QDriveMoveImit DriveMoveImit (pCtrl,mDriverReal,mMeasurmentImitator );

            double  arrDElX[QVARS] = {0.};
            const double VAlIntegratorTime = 1.;
            int iQuantRows = -1;
            DriveMoveImit.move(0., VAlMovingT, VAlIntegrStep
                             ,VAlIntegratorTime ,NULL, &iQuantRows);
            // double *parrDynamicPhVect = DriveMoveImit.mpCtrl->marrIntegrators;
            // MtrxMinusMatrx(parrDynamicPhVect
            //  , parrStatSolutionParams[iz * QUantTetta * QUantOm + itet * QUantOm + iOm].marrStatPhVect,QVARS, 1, arrDElX);
            // arrDElX [3] -= arrObjective[1] * VAlMovingT;
             memcpy(&(arrDynamicDelX[QVARS * (iz * QUantTetta * QUantOm + itet * QUantOm + iOm )]), DriveMoveImit.mpCtrl->marrIntegrators, QVARS * sizeof(double));
           }
        }
    }
}
//------------------------------------------------------------
double QParamsID::calcDelAlf_StatMeth(int quantExprs, QStatSolutionParams *parrStatSolutionParams
        , double *arrDynamicDelX,  double *arrDElAlf)
{


    double *arrCutting = new double [QPARAMS *mlenTargNums];
    memset(arrCutting, 0,QPARAMS * mlenTargNums * sizeof(double));
    for (int i =0; i < mlenTargNums; ++i)
    {
      arrCutting[miarrTargNums[i] * mlenTargNums + i] = 1.;
    }

    ///

double *arrA = new double [quantExprs * QVARS * mlenTargNums];

memset(arrA, 0, quantExprs * QVARS * mlenTargNums * sizeof(double));
double *arr_dF_po_dx_plus_BC_delX = new double [quantExprs * QVARS ];
memset(arr_dF_po_dx_plus_BC_delX, 0, quantExprs * QVARS  * sizeof(double));

double *arr_dF_po_dAlf = new double[QVARS *QPARAMS];
double *arr_dF_po_dx = new double[QVARS *QVARS];
double *arrB = new double[2 * QVARS];
double *arrT0 = new double[QVARS * QVARS];
double *arrT1 = new double[QVARS * QVARS];
//double *arrDelX = new double[QVARS];
double *arr_dBU_po_dIndL = new double[QVARS];
double *arr_dF_po_dAlf_Targ = new double [QVARS *mlenTargNums];
mDriverModel.calcMtrxB(arrB);

for(int i =0; i < quantExprs; ++i)
{
    double arrObjective[2] = {0.};
    QCtrlFollow ModelCtrlFollow (mDriverModel.mElectMotor,mDriverModel.mpLoad
                    ,marrSpreadParams, mDriverModel.mMomOut, 0., 0., 0.,arrObjective,constDbl_mh);


    QCtrl *pCtrl = &ModelCtrlFollow;
    //заполнение b

  (*pCtrl).fill_df_po_px_and_mtrxB_(parrStatSolutionParams[i].marrStatPhVect
                          ,arr_dF_po_dx,arrB);
  MtrxMultMatrx(arrB, QVARS, 2,  parrStatSolutionParams[i].marrGears, QVARS, arrT0);
  MtrxSumMatrx(arr_dF_po_dx, arrT0,QVARS, QVARS, arrT1) ;


  // провеерка

// TComp cmparrRoots[4];
// solvCharactEq_4_New(arrT1, cmparrRoots);
 //int yyui=0;

  double *parrDelX = &(arrDynamicDelX[i * QVARS]);

  MtrxMultMatrx(arrT1, QVARS, QVARS,  parrDelX, 1, &(arr_dF_po_dx_plus_BC_delX[i * QVARS]));
  ///

  //заполнение a
  (*pCtrl).calc_dF_po_dAlf(parrStatSolutionParams[i].marrStatPhVect, arr_dF_po_dAlf);
  double arrU[2] = {0.}, arrDelU[2] = {0.};
  MtrxMultMatrx(parrStatSolutionParams[i].marrGears,2, QVARS,   parrDelX, 1, arrDelU);
  MtrxSumMatrx(arrDelU, parrStatSolutionParams[i].marrStatU,1, 2, arrU) ;
  (*pCtrl).calc_dBU_po_dIndL(arrU,arr_dBU_po_dIndL);
  for(int j =0; j < QVARS; ++j)
  {
    arr_dF_po_dAlf[ j * QPARAMS] +=   arr_dBU_po_dIndL[j];
  }
  MtrxMultMatrx(arr_dF_po_dAlf, QVARS, QPARAMS,  arrCutting, mlenTargNums, &(arrA[i * QVARS *mlenTargNums]));


}
double valFGr0 = NormVect(arr_dF_po_dx_plus_BC_delX, quantExprs * QVARS);
valFGr0 = valFGr0 * valFGr0;
double valReturn = QMinSquare::solvMinSquare(arrA,quantExprs *QVARS , mlenTargNums, arr_dF_po_dx_plus_BC_delX, arrDElAlf );

delete []arrA;
delete []arr_dF_po_dx_plus_BC_delX ;
delete []arr_dF_po_dAlf;
delete []arr_dF_po_dx;
delete []arrB ;
delete []arrT0 ;
delete []arrT1 ;

delete []arr_dBU_po_dIndL;
delete []arr_dF_po_dAlf_Targ;
delete []arrCutting;
return valReturn;
}

//--------------------------------------------------
//-------------------------------------------------------

double QParamsID::identify_by_StatMeth(int quantExprs, QStatSolutionParams *parrStatSolutionParams
                                       , double *arrStatDelX, QElectMotor *pMotorRez, QLoad *pLoadRez, double *pMomOut)

{
    double *arrCutting = new double [( QPARAMS) *mlenTargNums];
    memset(arrCutting, 0,( QPARAMS) * mlenTargNums * sizeof(double));
    for (int i =0; i < mlenTargNums; ++i)
    {
      arrCutting[miarrTargNums[i] * mlenTargNums + i] = 1.;
    }


   // QElectDriver DriverCur (MotorCur,  &LoadCur, mDriverModel.marrSpreadParams, mDriverModel.mMomOut);
    double arrGrFgr  [ QPARAMS] = {0.};


    double valFgr0 = DBL_MAX,valFgr = DBL_MAX, val_dFgr_po_dLamb = 0.;

    double *arrGrFgrCutting = new double [mlenTargNums];
    double *arrZCutting = new double [mlenTargNums];





    double *arrGrFgrCutting1 = new double [mlenTargNums];
    double *arrZCutting1 = new double [mlenTargNums];

    const double EPS = 0.0001;
    double arrObjective[2] = {0.};
    QLoad LoadCur = *(mDriverModel.mpLoad);
    QElectMotor MotorCur = mDriverModel.mElectMotor;
    QCtrlFollow ModelCtrlFollow0 (MotorCur,&LoadCur
                    ,marrSpreadParams, mDriverModel.mMomOut, 0., 0., 0.,arrObjective,constDbl_mh);

    QCtrl *pCtrl0 = &ModelCtrlFollow0;
    int i =0;
    for ( i =0; i < 1000; ++i)
    {
        bool bend = false;
      //  DriverCur.calcDynamicFgr_and_arrGrFgr( VAlIntegrStep,parrRezPointTraj
        //        ,quantDoneSteps, &valFgr , arrGrFgr);
        QCtrlFollow ModelCtrlFollow ((*pCtrl0).mElectMotor,mDriverModel.mpLoad
                        ,marrSpreadParams, (*pCtrl0).mMomOut, 0., 0., 0.,arrObjective,constDbl_mh);

        QCtrl *pCtrl = &ModelCtrlFollow;

        (*pCtrl).calcStatFgr_and_arrGrFgr( quantExprs, parrStatSolutionParams
                                            ,arrStatDelX, &valFgr , arrGrFgr);


      MtrxTranspMultMatrx(arrGrFgr, QPARAMS, 1, arrCutting,mlenTargNums , arrGrFgrCutting) ;
      if(NormSquareVect(arrGrFgrCutting, mlenTargNums) <= 0.000000000000001)
      {
          break;
      }
      MtrxTranspMultMatrx(arrGrFgr, QPARAMS, 1, arrCutting,mlenTargNums , arrZCutting) ;
      val_dFgr_po_dLamb = ScalProduct(arrGrFgrCutting , arrZCutting, mlenTargNums) ;

      // |arrGrFgrCutting| = 1
      NormalizeVect(arrGrFgrCutting, mlenTargNums);
      ///

     // поиск минимума по направлению
      double valFgrLocalMin = valFgr;

      QLoad load_temp = (*pCtrl0->mpLoad);
      QElectDriver Driver_temp((*pCtrl0).mElectMotor,  &load_temp, mDriverModel.marrSpreadParams, (*pCtrl0).mMomOut);
      QCtrlFollow ModelCtrlFollow_temp (mDriverModel.mElectMotor,mDriverModel.mpLoad
                      ,marrSpreadParams, mDriverModel.mMomOut, 0., 0., 0.,arrObjective,constDbl_mh);

      QCtrl *pCtrl_temp = &ModelCtrlFollow_temp;
      ///////////////////////////////////////////////////////////////
      //////////  поиск интервала минимизации /////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////
       double coeff = 2.;
      double lamb0 =  0.;
      double lamb1 =  0.0001/ coeff;
      double lamb2 =  0.0001/ coeff;
      double valFgrTemp0 = valFgr;
      double valFgrTemp1 = valFgr;
      double valFgrTemp2 = 0.;
      for (int j =0; j < 100; ++j)
      {
          lamb2 = lamb1 * coeff;
          if(lamb2 >= 10000.)
          {
             bend = true;
             break;
          }

          for(int k =0; k < mlenTargNums; ++k)
          {
              double valParam = (*pCtrl_temp).getParam(miarrTargNums[k]) - lamb2 * arrGrFgrCutting[k];
             (*pCtrl_temp).setParam(miarrTargNums[k], valParam);
          }




              (*pCtrl_temp).calcStatFgr( quantExprs, parrStatSolutionParams, arrStatDelX,&valFgrTemp2 );

          if (valFgrTemp2 > valFgrTemp1)
          {
            break;
          }

          lamb0 = lamb1;
          lamb1 = lamb2;
          valFgrTemp0 = valFgrTemp1;
          valFgrTemp1 = valFgrTemp2;
      }
   /// интервал минимизации -  [lamb0; lamb2]
         if (bend)
         {
           break;
         }
        valFgrLocalMin = valFgrTemp0;

      ////////////////////////////////////////////////////////////////////////
      /// \brief lamb
      ///
        double valGold = (1. + sqrt(5))/2.;
        double a = lamb0;
        double b = lamb2;
        double x1 = b - (b-a)/ valGold;
        double x2 = a + (b-a)/ valGold;
        double y1 = (*pCtrl_temp).calcStatFgr_lambda(quantExprs, parrStatSolutionParams
                                                      , arrStatDelX,x1, miarrTargNums
                                                , mlenTargNums, arrGrFgrCutting);
        double y2 = (*pCtrl_temp).calcStatFgr_lambda(quantExprs, parrStatSolutionParams
                                                      , arrStatDelX,x2, miarrTargNums
                                                , mlenTargNums, arrGrFgrCutting);

      double lamb_opt = -1.;
      const int NUmi = 1500;
      for (int j =0; j < NUmi; ++j)
      {
          lamb_opt = (x1 + x2)/ 2.;
          if (i%50 ==0)
          {
             y1 = (*pCtrl_temp).calcStatFgr_lambda(quantExprs, parrStatSolutionParams
                                                            , arrStatDelX,x1, miarrTargNums
                                                      , mlenTargNums, arrGrFgrCutting);
             y2 = (*pCtrl_temp).calcStatFgr_lambda(quantExprs, parrStatSolutionParams
                                                            , arrStatDelX,x2, miarrTargNums
                                                      , mlenTargNums, arrGrFgrCutting);

          }
          if ((fabs(y1 -y2)<= EPS)||(fabs(x1 -  x2)< 0.00000001))
          {

            break;
          }
          if (y1 >= y2)
          {
            a = x1;
            x1 = x2;
            y1 = y2;
            x2 = a + (b-a)/ valGold;
            y2 = (*pCtrl_temp).calcStatFgr_lambda(quantExprs, parrStatSolutionParams
                                                           , arrStatDelX,x2, miarrTargNums
                                                    , mlenTargNums, arrGrFgrCutting);
          }
          else
          {
            b= x2;
            x2 = x1;
            y2 = y1;
            x1 = b - (b-a)/ valGold;
            y1 = (*pCtrl_temp).calcStatFgr_lambda(quantExprs, parrStatSolutionParams
                                                           , arrStatDelX,x1, miarrTargNums
                                                     , mlenTargNums, arrGrFgrCutting);
          }
      }

       valFgr = (*pCtrl).calcStatFgr_lambda(quantExprs, parrStatSolutionParams
                                               , arrStatDelX,lamb_opt, miarrTargNums
                                         , mlenTargNums, arrGrFgrCutting);

       if(valFgr > valFgr0)
       {
           break;
       }

      for(int k =0; k < mlenTargNums; ++k)
      {
         double valParam = (*pCtrl0).getParam(miarrTargNums[k]) - lamb_opt * arrGrFgrCutting[k];
         (*pCtrl0).setParam(miarrTargNums[k], valParam);
      }

      if(fabs(valFgr - valFgr0) <= EPS)
      {
          break;
      }

      valFgr0 = valFgr ;
    }


    delete []arrCutting;

    delete []arrGrFgrCutting;
    delete []arrZCutting ;


    delete []arrGrFgrCutting1;
    delete []arrZCutting1 ;
    *pMotorRez = (*pCtrl0).mElectMotor;
    *pLoadRez = *(pCtrl0->mpLoad);

    *pMomOut =  (*pCtrl0).mMomOut;
    return sqrt(valFgr0/ (double(quantExprs)));

}
//моделирование идентификации параметров привода методом динамическим
// DriverModel - это представление о приводе, заложенное в модель
// на основании моельных данных рассчитываются коэффициенты передачи (управления)
// DriverReal - это реальный привод, с параметрами, которые могут отличаться от заложенных в модели
// ОПИСАНИЕ МЕТОДА
// моделируются набор экспериментов, по выводу сиситемы на стационарный режим.
// описание одного эксперимента:
//  1.заданы начальные условия сиситемы
//  2.заданы целевые условия задачи слежения - Tetta0, Omega0
//  3. задан набор корней характеристического уравнения
//  4. производится моелирование выывода реального привода на стационарный режим
//  5.данные экспееримента апоминаются
//После того как проиведено моеделирование реультатов экспериментов,
// производится оцуенка реальных параметров привода методом наименьших квадратов
// INPUT:
// DriverModel - модельный привод
// DriverReal - реальный привод
// MeasurmentImitator - модель измерительной сиситемы
// arrSpreadParams - вектор СКЗ разбросов параметров привода (в расчета не участвует)
// deltaTetta - шаг по Tetta0
// quantTetta - к-во шагов по Tetta0
// deltaOmega - шаг по Omega0
// quantOm - к-во шагов по Omega0
// valMovingT - врмя вывода на стационарный режим
 // cmpArrLamb[QVARS * QUantLambSets] - массив наборов корней характеристич уравнения
 // QUantLambSets - к-во наборов корней
// замечание - общее к-во экспериметов равно quantTetta * quantOm *QUantLambSets
 // iarrTargNums[lenTargNums] - номера параметров сиситемы, по которым производится идентификация
 // lenTargNums - длина массива номеров параметров сиситемы, по которым производится идентификация
// замечание - вектор  параметров сиситемы ореленв классе QElectDriver
// в него входят, в частности дуктивность и соппроивление, с номерами 0 и 1 соответственно
// но производить идентфикацию этих параметов нет нобходимости,
// т.к. они могум быть измерены с очен выскй тчностью мульимтом ил RCL-метром
// OUTPUT:
// возвращеет оценочный привод

QElectDriver QParamsID::identify_by_DynamicMeth(const QElectDriver DriverModel,const QElectDriver DriverReal
        , const QMeasurmentImitator MeasurmentImitator, const double *arrSpreadParams
        ,const double  VAlU0,const double   VAlT1,const double VAlMovingT,const double TrgOm, const double TetRotor0
        , const int QuantIsp, int *iarrTargNums,const int LEnTargNums, const double VAlIntegrStep,double *arrBuff
        , int *pquantDoneSteps)

{
    QParamsID ParamsID( DriverModel, DriverReal
                         , MeasurmentImitator,arrSpreadParams,iarrTargNums,LEnTargNums);


    ParamsID.imitateDynamicExperimentsRez( VAlU0,VAlT1
         , VAlMovingT,TrgOm,TetRotor0, QuantIsp, VAlIntegrStep ,arrBuff, pquantDoneSteps);


    return DriverModel;
}

//--------------------------------------
void QParamsID::imitateDynamicExperimentsRez(const double  VAlU0, const double VAlT1
        ,const double  VAlMovingT,const double TrgOm, const double TetRotor0
        , const int QUantIsp,const double  VAlIntegrStep ,double *arrBuff,int *quantDoneSteps)
{


    double arrObjective[2] = {0.};

    arrObjective[1] = TrgOm;

   /* QCtrlTrad ModelCtrlTrad (mDriverModel.mElectMotor,mDriverModel.mpLoad
                    ,marrSpreadParams, mDriverModel.mMomOut, valTettaBegin, valOmegaBegin, 0.
                     ,arrObjective,constDbl_mh, VAlU0, VAlT1, VAlMovingT);


    QCtrl *pCtrl = &ModelCtrlTrad;
    QDriveMoveImit DriveMoveImit (pCtrl,mDriverReal,mMeasurmentImitator );

   // DriveMoveImit.move(0.,valMovingT, VAlIntegrStep, &(arrDynamicPhVect[QVARS * (iz * quantTetta * quantOm + itet * quantOm + iOm )]));
   double arrDynamicPhVect[QVARS] = {0.}, arrDElX[QVARS] = {0.};
*/

   const int QUantRows = int(VAlMovingT/VAlIntegrStep ) + 1;

   double *arrBuffCur = new double[QUantColsCSVReport * QUantRows];
   memset(arrBuffCur, 0, sizeof(double) * QUantColsCSVReport *(QUantRows ));
   memset(arrBuff, 0, sizeof(double) * QUantColsCSVReport *(QUantRows ));


   for (int i =0; i < QUantIsp; ++i)
   {

       //  начальное угловое пложение
       double valOmegaBegin =0.;

       QCtrlTrad ModelCtrlTrad (mDriverModel.mElectMotor,mDriverModel.mpLoad
                       ,marrSpreadParams, mDriverModel.mMomOut, TetRotor0, valOmegaBegin, 0.
                        ,arrObjective,constDbl_mh, VAlU0, VAlT1, VAlMovingT);
      memset(arrBuffCur, 0, sizeof(double) * QUantColsCSVReport *(QUantRows ));
       QCtrl *pCtrl = &ModelCtrlTrad;
       QDriveMoveImit DriveMoveImit (pCtrl,mDriverReal,mMeasurmentImitator );

       const double VAlIntegratorTime = 1.;

       DriveMoveImit.move(0., VAlMovingT, VAlIntegrStep
                        ,VAlIntegratorTime ,arrBuffCur, quantDoneSteps);
       MatrxMultScalar(arrBuffCur, QUantColsCSVReport, QUantRows , 1./((double)QUantIsp),arrBuffCur);
       MtrxSumMatrx(arrBuff, arrBuffCur,QUantColsCSVReport, QUantRows, arrBuff) ;
   }

   delete []arrBuffCur;

}

QElectDriver QParamsID::identify_by_DynamicMeth_( const double VAlIntegrStep,QRezPointTraj *parrRezPointTraj
        ,const int quantDoneSteps)

{
    double *arrCutting = new double [(2 + QPARAMS) *mlenTargNums];
    memset(arrCutting, 0,(2 + QPARAMS) * mlenTargNums * sizeof(double));
    for (int i =0; i < mlenTargNums; ++i)
    {
      arrCutting[miarrTargNums[i] * mlenTargNums + i] = 1.;
    }
    QElectDriver DriverCur = mDriverModel;
    double arrGrFgr  [2 + QPARAMS] = {0.};
   // double *arrZ = new double [2 + QPARAMS];
   // double *arrZCur = new double [2 + QPARAMS];

    double valFgr0 = DBL_MAX,valFgr = DBL_MAX, val_dFgr_po_dLamb = 0.;

    double *arrGrFgrCutting = new double [mlenTargNums];
    double *arrZCutting = new double [mlenTargNums];

    double arrGrFgr1 [2 + QPARAMS]= {0.};
   // double *arrZ1 = new double [2 + QPARAMS];
   // double *arrZCur1 = new double [2 + QPARAMS];

    double valFgr1 = 0., val_dFgr_po_dLamb1 = 0.;

    double *arrGrFgrCutting1 = new double [mlenTargNums];
    double *arrZCutting1 = new double [mlenTargNums];

    const double EPS = 0.0001;
    for (int i =0; i < 100; ++i)
    {
        DriverCur.calcDynamicFgr_and_arrGrFgr( VAlIntegrStep,parrRezPointTraj
                ,quantDoneSteps, &valFgr , arrGrFgr);

      //  DriverCur.calcDynamicFgr_and_arrGrFgr_temp( VAlIntegrStep,parrRezPointTraj
        //        ,quantDoneSteps, &valFgr , arrGrFgr);
      if(valFgr > valFgr0)
      {
          break;
      }
      if(fabs(valFgr - valFgr0) <= EPS)
      {
          break;
      }
      MtrxTranspMultMatrx(arrGrFgr,2 + QPARAMS, 1, arrCutting,mlenTargNums , arrGrFgrCutting) ;
      if(NormSquareVect(arrGrFgrCutting, mlenTargNums) <= EPS)
      {
          break;
      }
      MtrxTranspMultMatrx(arrGrFgr,2 + QPARAMS, 1, arrCutting,mlenTargNums , arrZCutting) ;
      val_dFgr_po_dLamb = ScalProduct(arrGrFgrCutting , arrZCutting, mlenTargNums) ;


      QElectDriver DriverCur1 = DriverCur;
      double deltaLamb = 0.0000001;
      for(int j = 0; j < mlenTargNums; ++j)
      {
          double valParam = DriverCur1.getParam(miarrTargNums[i]) + deltaLamb *arrGrFgr[miarrTargNums[i]];
         DriverCur1.setParam(miarrTargNums[i], valParam);
      }
      DriverCur1.calcDynamicFgr_and_arrGrFgr( VAlIntegrStep,parrRezPointTraj
              ,quantDoneSteps, &valFgr1 , arrGrFgr1);
      MtrxTranspMultMatrx(arrGrFgr1,2 + QPARAMS, 1, arrCutting,mlenTargNums , arrGrFgrCutting1) ;
      MtrxTranspMultMatrx(arrGrFgr1,2 + QPARAMS, 1, arrCutting,mlenTargNums , arrZCutting1) ;
      val_dFgr_po_dLamb1 = ScalProduct(arrGrFgrCutting1 , arrZCutting1, mlenTargNums) ;

    double val_d2Ggr_po_d2Lamb = (val_dFgr_po_dLamb1 - val_dFgr_po_dLamb)/deltaLamb;
     double valLambStep = -val_dFgr_po_dLamb/val_d2Ggr_po_d2Lamb;
     // double valLambStep = - valFgr/val_dFgr_po_dLamb;
      for(int j = 0; j < mlenTargNums; ++j)
      {
          double valParam = DriverCur.getParam(miarrTargNums[i]) + valLambStep*arrGrFgr[miarrTargNums[i]];
         DriverCur.setParam(miarrTargNums[i], valParam);
      }
      valFgr0 = valFgr ;
    }


    delete []arrCutting;

    delete []arrGrFgrCutting;
    delete []arrZCutting ;


    delete []arrGrFgrCutting1;
    delete []arrZCutting1 ;
    return DriverCur;
}


void QParamsID::identify_by_DynamicMeth_1( const double VAlIntegrStep,QRezPointTraj *parrRezPointTraj
        ,const int quantDoneSteps, QElectMotor *pMotorRez, QLoad *pLoadRez, double *pMomOut)

{
    double *arrCutting = new double [(2 + QPARAMS) *mlenTargNums];
    memset(arrCutting, 0,(2 + QPARAMS) * mlenTargNums * sizeof(double));
    for (int i =0; i < mlenTargNums; ++i)
    {
      arrCutting[miarrTargNums[i] * mlenTargNums + i] = 1.;
    }

    QLoad LoadCur = *(mDriverModel.mpLoad);
    QElectMotor MotorCur = mDriverModel.mElectMotor;
    QElectDriver DriverCur (MotorCur,  &LoadCur, mDriverModel.marrSpreadParams, mDriverModel.mMomOut);
   // QElectDriver DriverCur = mDriverModel;
    double arrGrFgr  [2 + QPARAMS] = {0.};
   // double *arrZ = new double [2 + QPARAMS];
   // double *arrZCur = new double [2 + QPARAMS];

    double valFgr0 = DBL_MAX,valFgr = DBL_MAX, val_dFgr_po_dLamb = 0.;

    double *arrGrFgrCutting = new double [mlenTargNums];
    double *arrZCutting = new double [mlenTargNums];

    double arrGrFgr1 [2 + QPARAMS]= {0.};
   // double *arrZ1 = new double [2 + QPARAMS];
   // double *arrZCur1 = new double [2 + QPARAMS];

    double valFgr1 = 0., val_dFgr_po_dLamb1 = 0.;

    double *arrGrFgrCutting1 = new double [mlenTargNums];
    double *arrZCutting1 = new double [mlenTargNums];

    const double EPS = 0.0001;

    for (int i =0; i < 100; ++i)
    {
        DriverCur.calcDynamicFgr_and_arrGrFgr( VAlIntegrStep,parrRezPointTraj
                ,quantDoneSteps, &valFgr , arrGrFgr);

      //  DriverCur.calcDynamicFgr_and_arrGrFgr_temp( VAlIntegrStep,parrRezPointTraj
        //        ,quantDoneSteps, &valFgr , arrGrFgr);

      MtrxTranspMultMatrx(arrGrFgr,2 + QPARAMS, 1, arrCutting,mlenTargNums , arrGrFgrCutting) ;
      if(NormSquareVect(arrGrFgrCutting, mlenTargNums) <= EPS)
      {
          break;
      }
      MtrxTranspMultMatrx(arrGrFgr,2 + QPARAMS, 1, arrCutting,mlenTargNums , arrZCutting) ;
      val_dFgr_po_dLamb = ScalProduct(arrGrFgrCutting , arrZCutting, mlenTargNums) ;

      // |arrGrFgrCutting| = 1
      NormalizeVect(arrGrFgrCutting, mlenTargNums);
      ///

     // поиск минимума по направлению
      double valFgrLocalMin = valFgr;

      QLoad load_temp = *(DriverCur.mpLoad);
      QElectDriver Driver_temp(DriverCur.mElectMotor,  &load_temp, mDriverModel.marrSpreadParams, DriverCur.mMomOut);


      double lamb = -1.;
      for (int j =0; j < 1000; ++j)
      {
          lamb = (double(j) +1) * 0.00001;
          // Driver_temp = QElectDriver (MotorCur,  &load_temp, mDriverModel.marrSpreadParams, mDriverModel.mMomOut);

         // QLoad load = *(DriverCur.mpLoad);
         // Driver_temp.mpLoad = &load;
          for(int k =0; k < mlenTargNums; ++k)
          {
              double valParam = DriverCur.getParam(miarrTargNums[k]) - lamb * arrGrFgrCutting[k];
             Driver_temp.setParam(miarrTargNums[k], valParam);
          }
           double valFgr1 = -1;
          Driver_temp.calcDynamicFgr( VAlIntegrStep,parrRezPointTraj
                  ,quantDoneSteps, &valFgr1 );



          if(valFgr1 > valFgrLocalMin)
          {
            valFgrLocalMin = valFgr1;
            break;
          }
         valFgrLocalMin = valFgr1 ;
      }
       valFgr =valFgrLocalMin;


      for(int k =0; k < mlenTargNums; ++k)
      {
         double valParam = DriverCur.getParam(miarrTargNums[k]) - lamb * arrGrFgrCutting[k];
         DriverCur.setParam(miarrTargNums[k], valParam);
      }
      if(valFgr > valFgr0)
      {
          break;
      }
      if(fabs(valFgr - valFgr0) <= EPS)
      {
          break;
      }

      valFgr0 = valFgr ;
    }


    delete []arrCutting;

    delete []arrGrFgrCutting;
    delete []arrZCutting ;


    delete []arrGrFgrCutting1;
    delete []arrZCutting1 ;
    *pMotorRez = DriverCur.mElectMotor;
    *pLoadRez = *(DriverCur.mpLoad);
    *pMomOut =  DriverCur.mMomOut;

}
//-------------------------------------------------------

double QParamsID::identify_by_DynamicMeth_2( const double VAlIntegrStep,QRezPointTraj *parrRezPointTraj
        ,const int quantDoneSteps, QElectMotor *pMotorRez, QLoad *pLoadRez, double *pMomOut)

{
    double *arrCutting = new double [(2 + QPARAMS) *mlenTargNums];
    memset(arrCutting, 0,(2 + QPARAMS) * mlenTargNums * sizeof(double));
    for (int i =0; i < mlenTargNums; ++i)
    {
      arrCutting[miarrTargNums[i] * mlenTargNums + i] = 1.;
    }

    QLoad LoadCur = *(mDriverModel.mpLoad);
    QElectMotor MotorCur = mDriverModel.mElectMotor;
    QElectDriver DriverCur (MotorCur,  &LoadCur, mDriverModel.marrSpreadParams, mDriverModel.mMomOut);
   // QElectDriver DriverCur = mDriverModel;
    double arrGrFgr  [2 + QPARAMS] = {0.};
   // double *arrZ = new double [2 + QPARAMS];
   // double *arrZCur = new double [2 + QPARAMS];

    double valFgr0 = DBL_MAX,valFgr = DBL_MAX, val_dFgr_po_dLamb = 0.;

    double *arrGrFgrCutting = new double [mlenTargNums];
    double *arrZCutting = new double [mlenTargNums];

    double arrGrFgr1 [2 + QPARAMS]= {0.};
   // double *arrZ1 = new double [2 + QPARAMS];
   // double *arrZCur1 = new double [2 + QPARAMS];

    double valFgr1 = 0., val_dFgr_po_dLamb1 = 0.;

    double *arrGrFgrCutting1 = new double [mlenTargNums];
    double *arrZCutting1 = new double [mlenTargNums];

    const double EPS = 0.0001;

    for (int i =0; i < 1000; ++i)
    {
        bool bend = false;
      //  DriverCur.calcDynamicFgr_and_arrGrFgr( VAlIntegrStep,parrRezPointTraj
        //        ,quantDoneSteps, &valFgr , arrGrFgr);


        DriverCur.calcDynamicFgr_and_arrGrFgr( VAlIntegrStep,parrRezPointTraj
                ,quantDoneSteps, &valFgr , arrGrFgr);

      //  DriverCur.calcDynamicFgr_and_arrGrFgr_temp( VAlIntegrStep,parrRezPointTraj
        //        ,quantDoneSteps, &valFgr , arrGrFgr);

      MtrxTranspMultMatrx(arrGrFgr,2 + QPARAMS, 1, arrCutting,mlenTargNums , arrGrFgrCutting) ;
      if(NormSquareVect(arrGrFgrCutting, mlenTargNums) <= EPS * EPS)
      {
          break;
      }
      MtrxTranspMultMatrx(arrGrFgr,2 + QPARAMS, 1, arrCutting,mlenTargNums , arrZCutting) ;
      val_dFgr_po_dLamb = ScalProduct(arrGrFgrCutting , arrZCutting, mlenTargNums) ;

      // |arrGrFgrCutting| = 1
      NormalizeVect(arrGrFgrCutting, mlenTargNums);
      ///

     // поиск минимума по направлению
      double valFgrLocalMin = valFgr;

      QLoad load_temp = *(DriverCur.mpLoad);
      QElectDriver Driver_temp(DriverCur.mElectMotor,  &load_temp, mDriverModel.marrSpreadParams, DriverCur.mMomOut);
      ///////////////////////////////////////////////////////////////
      //////////  поиск интервала минимизации /////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////
       double coeff = 2.;
      double lamb0 =  0.;
      double lamb1 =  0.1/ coeff;
      double lamb2 =  0.1/ coeff;
      double valFgrTemp0 = valFgr;
      double valFgrTemp1 = valFgr;
      double valFgrTemp2 = 0.;
      for (int j =0; j < 100; ++j)
      {
          lamb2 = lamb1 * coeff;
          if(lamb2 >= 10000.)
          {
             bend = true;
             break;
          }
          // Driver_temp = QElectDriver (MotorCur,  &load_temp, mDriverModel.marrSpreadParams, mDriverModel.mMomOut);

         // QLoad load = *(DriverCur.mpLoad);
         // Driver_temp.mpLoad = &load;
          for(int k =0; k < mlenTargNums; ++k)
          {
              double valParam = DriverCur.getParam(miarrTargNums[k]) - lamb2 * arrGrFgrCutting[k];
             Driver_temp.setParam(miarrTargNums[k], valParam);
          }

          Driver_temp.calcDynamicFgr( VAlIntegrStep,parrRezPointTraj
                  ,quantDoneSteps, &valFgrTemp2 );
          if (valFgrTemp2 > valFgrTemp1)
          {
            break;
          }

          lamb0 = lamb1;
          lamb1 = lamb2;
          valFgrTemp0 = valFgrTemp1;
          valFgrTemp1 = valFgrTemp2;
      }
   /// интервал минимизации -  [lamb0; lamb2]
         if (bend)
         {
           break;
         }
        valFgrLocalMin = valFgrTemp0;
         // if(valFgr1 > valFgrLocalMin)
        //  {
          //  valFgrLocalMin = valFgr1;
           // break;
         // }
         //valFgrLocalMin = valFgr1 ;



      ////////////////////////////////////////////////////////////////////////
      /// \brief lamb
      ///


        double valGold = (1. + sqrt(5))/2.;
        double a = lamb0;
        double b = lamb2;
        double x1 = b - (b-a)/ valGold;
        double x2 = a + (b-a)/ valGold;
        double y1 = DriverCur.calcDynamicFgr_lambda( VAlIntegrStep,parrRezPointTraj
                                          , quantDoneSteps, x1, miarrTargNums
                                          , mlenTargNums, arrGrFgrCutting);
        double y2 = DriverCur.calcDynamicFgr_lambda( VAlIntegrStep,parrRezPointTraj
                                          , quantDoneSteps, x2, miarrTargNums
                                          , mlenTargNums, arrGrFgrCutting);

      double lamb_opt = -1.;
      const int NUmi = 500;
      for (int j =0; j < NUmi; ++j)
      {
          lamb_opt = (x1 + x2)/ 2.;
          if (i%50 ==0)
          {
              y1 = DriverCur.calcDynamicFgr_lambda( VAlIntegrStep,parrRezPointTraj
                                                        , quantDoneSteps, x1, miarrTargNums
                                                        , mlenTargNums, arrGrFgrCutting);
              y2 = DriverCur.calcDynamicFgr_lambda( VAlIntegrStep,parrRezPointTraj
                                                        , quantDoneSteps, x2, miarrTargNums
                                                        , mlenTargNums, arrGrFgrCutting);

          }
          if ((fabs(y1 -y2)<= EPS)||(fabs(x1 -  x2)< 0.00000001))
          {

            break;
          }
          if (y1 >= y2)
          {
            a = x1;
            x1 = x2;
            y1 = y2;
            x2 = a + (b-a)/ valGold;
            y2 = DriverCur.calcDynamicFgr_lambda( VAlIntegrStep,parrRezPointTraj
                                                      , quantDoneSteps, x2, miarrTargNums
                                                      , mlenTargNums, arrGrFgrCutting);


          }
          else
          {
            b= x2;
            x2 = x1;
            y2 = y1;
            x1 = b - (b-a)/ valGold;
            y1 = DriverCur.calcDynamicFgr_lambda( VAlIntegrStep,parrRezPointTraj
                                                      , quantDoneSteps,x1, miarrTargNums
                                                      , mlenTargNums, arrGrFgrCutting);


          }
      }
         /* lamb_opt = lamb0 + (lamb2 - lamb0) / (double (NUmi)) *(double(j) +1.) ;

          for(int k =0; k < mlenTargNums; ++k)
          {
              double valParam = DriverCur.getParam(miarrTargNums[k]) - lamb_opt * arrGrFgrCutting[k];
             Driver_temp.setParam(miarrTargNums[k], valParam);
          }
           double valFgr1 = -1;
          Driver_temp.calcDynamicFgr( VAlIntegrStep,parrRezPointTraj
                  ,quantDoneSteps, &valFgr1 );



          if(valFgr1 > valFgrLocalMin)
          {
            valFgrLocalMin = valFgr1;
            break;
          }
         valFgrLocalMin = valFgr1 ;
      }
      */
     //  valFgr = valFgrLocalMin;
       valFgr = DriverCur.calcDynamicFgr_lambda( VAlIntegrStep,parrRezPointTraj
                                                 , quantDoneSteps, lamb_opt, miarrTargNums
                                                 , mlenTargNums, arrGrFgrCutting);
       if(valFgr > valFgr0)
       {
           break;
       }

      for(int k =0; k < mlenTargNums; ++k)
      {
         double valParam = DriverCur.getParam(miarrTargNums[k]) - lamb_opt * arrGrFgrCutting[k];
         DriverCur.setParam(miarrTargNums[k], valParam);
      }

      if(fabs(valFgr - valFgr0) <= EPS)
      {
          break;
      }

      valFgr0 = valFgr ;
    }


    delete []arrCutting;

    delete []arrGrFgrCutting;
    delete []arrZCutting ;


    delete []arrGrFgrCutting1;
    delete []arrZCutting1 ;
    *pMotorRez = DriverCur.mElectMotor;
    *pLoadRez = *(DriverCur.mpLoad);
    *pMomOut =  DriverCur.mMomOut;
    return sqrt(valFgr0/ (double(quantDoneSteps)));

}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

QRezPointTraj::QRezPointTraj()
{
    mt =0.; // время привязки точки траектории
    memset(marrX, 0, QVARS * sizeof(double));
    memset(marrU, 0, 2 * sizeof(double));

}

// конструктор копирования
 QRezPointTraj ::QRezPointTraj (const QRezPointTraj &R)
 {
     mt = R.mt;
     memcpy(marrX, R.marrX, QVARS * sizeof(double));
     memcpy(marrU, R.marrU, 2 * sizeof(double));
 }
 // оператор присваивания
 QRezPointTraj &QRezPointTraj::operator=(const QRezPointTraj  &R)
 {
     mt = R.mt;
     memcpy(marrX, R.marrX, QVARS * sizeof(double));
     memcpy(marrU, R.marrU, 2 * sizeof(double));
   return *this ;
 }
// парам констр
QRezPointTraj :: QRezPointTraj(  const double t, const double *arrX ,  const double *arrU)
{
    mt = t;
    memcpy(marrX, arrX,QVARS * sizeof(double));
    memcpy(marrU, arrU, 2 * sizeof(double));
}
