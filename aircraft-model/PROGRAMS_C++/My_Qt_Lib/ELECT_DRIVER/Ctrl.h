#ifndef CTRL_H
#define CTRL_H

#include "ElectDriver.h"
#include "Segment.h"
#include "Environment.h"
#include "Plane.h"
#include "Comp.h"
#include "StatSolutionParams.h"
#include "MeasurmentImitator.h"
#include "Filtr.h"




// класс описывает динамику изменения фазового вектора
// описывающего привод
// marrPhVect[0] - угловая скорость ротора
// marrPhVect[1] - ток статора во вращ СК Id
// marrPhVect[2] - ток статора во вращ СК Iq
// marrPhVect[3] - угол проворота ротора Tetta




class QElectDriver;
class QSegment;
class TEnvironment;
class TPlane;
class QVeloCtrl;
class QCtrl;
class QStatSolutionParams;
class QMeasurmentImitator;
class QFiltr;
class TComp;
class QCtrl:public QElectDriver
{
public:


    // время привязки начального вектора состояния
    double mT0;
    // текущее время
    double mTCur;

  //  QStatSolutionParams mStatSolutionParams;
    TComp mCmpArrLamb[4];

    //цель управления - уговое положения и и угловая скорость
    // 1.если управление производиться по угловой скорости, то
    // величина marrObjective[0] не имеет значения, а marrObjective[1]
    // равна желаемой угловой скорости вращения
    //2. если управление производиться по угловому положению, то
    // угловая скорость: marrObjective[1] = 0
    // а marrObjective[0] равно желаемому угловому положению
    // 3. если управления производиться в интересах слежения за
    // объектом, двигающимся равномерно и прямолиненйно, то
    // marrObjective[1] - угловая скорость ,
    // а marrObjective[0] -  угловое положения объекта в момент времени mT0
    //
  //  double marrObjective[2];

    // имитатор измерений фазового вектора


    // структура с параметрами фильтра
    QFiltr mFiltr;

      QCtrl();
    // конструктор копирования
      QCtrl (const  QCtrl &R);

     // оператор присваивания
      QCtrl  &operator=( const QCtrl  &R);



      //
      QCtrl ( const QElectMotor ElectMotor ,  QLoad*  Load
                    , const double *arrSpreadParams,const double VAlMomOut, const double VAlTettaBegin
         ,const double VAlOmegaBegin , const double T0, const double TCur, const double valh,TComp *pCmpArrLamb);
      //-------------------------------------------------------------------------
      //-------------------------------------------------------------------------
      //-------------------------------------------------------------------------


     bool calcSwithOverTimes(const double VAlUqu0
                       ,const double VAlOmegaStat, double *pvalt0,double *pvalt1, int *pisignum0);



     inline double sign_(const double x)
     {
         return (x >=0.) ?1.:-1.;
     }

     void calcGearMtrx_forPosCtrl(const double VAlTettaStat, const double VAlLamb0,const double VAlLamb1
                                                ,const TComp VAlLamb2,const TComp VAlLamb3, double *arrC0 );


     virtual void  CalcEigenValues( const double VAlTettaStat, const double VAlOmegaStat
                                       ,const double *arrC,const double VAlMomOut, TComp* cmpArrEigen);

     virtual void calcGearMtrx(const double *arrObjective, TComp* CmpArrLamb,const double VAlMomOut, double* arrC);

     virtual void calcCurU(const double *arrObjective,const double VAlTObjective,const double VAlMomOut, double *arrU);



     //void calcStationarySolution(const double *arrObjective, double*arrStatinaryPhVect,double *valUd, double*valUq, double *pvalDelFi);

     void calcStationaryParams(const double *arrObjective,  TComp* CmpArrLamb
              , QStatSolutionParams* pQStatSolutionParams,const double VAlMomOut);

     void  calc_df_po_dx_plus_BC(const double *arrObjective,const double *arrC, double *arrFGr,const double VAlMomOut);

     void processMeasure_and_calcCurU(const double *arrObjective,const double VAlTObjective
        ,QDriverMeasure &MEasureCur, double &VAlTcur, double *arrPrevU,const double VAlMomOut, double *arrU);

     void calc_Kww(const double VAlt/*,double *arrU*/, double *arrKww);

     void calc_d2F_po_dx_po_dL(double *arr_d2F_po_dx_po_dL);

     void   calc_d2F_po_dx_po_dR(double *arr_d2F_po_dx_po_dR  );

     //---------------------------------------
     void   calc_d2F_po_dx_po_dPsif(double *arr_d2F_po_dx_po_dPsif  );

     //---------------------------------------
     void   calc_d2F_po_dx_po_dJ(double *arr_d2F_po_dx_po_dJ  );

     //---------------------------------------
     void   calc_d2F_po_dx_po_dCx(double *arr_d2F_po_dx_po_dCx  );

     //---------------------------------------
     void   calc_d2F_po_dx_po_dCv(double *arr_d2F_po_dx_po_dCv  );

     virtual double calcMomResidual(const double VAlTettaStat,const double VAlOmegaStat);

     void  calcStatFgr_and_arrGrFgr( int quantExprs, QStatSolutionParams *parrStatSolutionParams
     , double *arrDynamicDelX,double  *pvallFgr, double *arrGrFgr,const double VAlMomOut);

     double   calcStatFgr_lambda(int quantExprs, QStatSolutionParams *parrStatSolutionParams
                                               , double *arrDynamicDelX,const double VAllamb, const int *IArrTargNums
                          , const int LEnTargNums, const double *ARrGrFgrCutting);

     void  calcStatFgr( int quantExprs, QStatSolutionParams *parrStatSolutionParams
                        , double *arrDynamicDelX,double  *pvallFgr,const double VAlMomOut);

     void  calc_vectFi_and_mtrx_dFi_po_dAlf( QStatSolutionParams StatSolutionParams
              , double *arrXZv, double *arrFi, double *mtrx_dFi_po_dAlf ,const double VAlMomOut);

     bool calcStableSolution(const double *arrObjective, const double valUd
    , const double valUq, double*arrStatinaryPhVect,const double VAlMomOut);

     virtual void correctMomOut(const double VAlTettaZv
                                ,const double VAlDispTetta    , const double VAlW, const double VAlTcur);

     void doExtrapolation(const double VAlT,const double VAlIntegrStep0, double *arrUPrev,const double VAlMomOut);

     //void processingDriverMeasure(QMeasure &MEasureCur, double &VAlTcur, double *arrPrevU);

     void processingDriverMeasure(QDriverMeasure &MEasureCur, double &VAlTcur, double *arrPrevU);

     virtual double calcEstMomOut();

     virtual  int createInputDataReport_Ctrl(wchar_t*FileName, const bool bHeader);

     virtual int createInputDataReport_Ctrl_Inherited(wchar_t*FileName, const bool bHeader);

};



#endif // DRIVETRAJ_H
