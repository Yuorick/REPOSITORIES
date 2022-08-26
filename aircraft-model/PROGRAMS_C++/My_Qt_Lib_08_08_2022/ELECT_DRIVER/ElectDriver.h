#ifndef ELECTDRIVER_H
#define ELECTDRIVER_H
#include "ElectMotor.h"
#include "Load.h"
#include "Brls.h"


// к-во переменных в фазовом векторе
#define QVARS 4

#define QPARAMS 7


class QElectMotor;
class QLoad;
class TComp;
class QBrls;
class QRezPointTraj ;


// marrSpreadParams[7] - вектор разбросов параметров привода
// перенумеруем переметры элдвигателя alf[7] в следующем порядке:
// mInvL - marrSpreadParams[0]
// mResist  - marrSpreadParams[1]
// mPsi_f   - marrSpreadParams[2]
// InvJ        - marrSpreadParams[3]
// Cx       -marrSpreadParams[4]
// Cv       - marrSpreadParams[5]
// MomOut       - marrSpreadParams[6]
// эти параметры известны с точностью d = delta
//

class QElectDriver
{
     private:

    double  calcInductMom(const double VAl_Iqu);
    double mInvSumJ0;


    ///
public:
    QElectMotor mElectMotor;
    QLoad* mpLoad;
    // дополнительный внешний момент
   // double mMomOut;


    QElectDriver();

    // конструктор копирования
      QElectDriver (const  QElectDriver &R);

     // оператор присваивания
      QElectDriver  &operator=( const QElectDriver  &R);

      // парам конструктор
      QElectDriver  (const QElectMotor ElectMotor
                     ,  QLoad*  Load, const double *arrSpreadParams);


      //вектор разбросов параметров электродвигателя - L, R, Psif
      double marrSpreadParams[QPARAMS];


      void calcMagnRightPart_Without_ResidualMom(double *arrPhVect, double *arrRightPart );

      double max_(const double a,const double b);

      double sign_(const double a);


       void calcFx( TEnvironment &Environment, double *arrPhVect, double *arrFx,const double VAlMomOut);

       virtual  double calcMomResidual(const double VAlTettaStat,const double VAlOmegaStat);

       void fill_df_po_px_and_mtrxB_(const double *arrStatinaryPhVect
                        ,double *arr_dF_po_dx,double * arr_dF_po_dW);
       void dragPhVect(TEnvironment &Environment, double *arrPhVectBegin, double *arrU, const double VAlTime
                            , const double VAlIntegrStep0, double *arrPhVectEnd,const double VAlMomOut) ;

       void doEilerStep(TEnvironment &Environment,double *arrPhVectCur
        , double *arrU, const double VAlStep,const double VAlMomOut);

       void calcMtrxB(double *arrB);

       void calc_dF_po_dAlf(double *arrEstPhVect,double *arr_dF_po_dAlf,const double VAlMomOut);

       void calc_dBU_po_dIndL(double *arrU,double *arr_dBU_po_dIndL);

       double  getParam(const int i);

       void  setParam(const int i, const double VAlParam);

       double calcMomResidual_plus_AirMom_plus_MomOut(TEnvironment &Environment
               ,const double VAlTetta,const double VAlOmega,const double VAlMomOut);

       void  calcDynamicFgr_and_arrGrFgr(const double VAlIntegrStep,QRezPointTraj *parrRezPointTraj
            ,const int LEnarrRezPointTraj,double  *pvallFgr, double *arrGrFgr,const double VAlMomOut);

       void  calcDynamicFgr_and_arrGrFgr_temp(const double VAlIntegrStep0,QRezPointTraj *parrRezPointTraj
                     ,const int LEnarrRezPointTraj,double  *pvallFgr, double *arrGrFgr,const double VAlMomOut);

       void  calcDynamicFgr(const double VAlIntegrStep0,QRezPointTraj *parrRezPointTraj
                     ,const int LEnarrRezPointTraj,double  *pvallFgr,const double VAlMomOut);

       double   calcInvSumJ0();

       double   getInvSumJ0();

       double   getJPayLoad();

       double   calcDynamicFgr_lambda(const double VAlIntegrStep,QRezPointTraj *parrRezPointTraj
                            ,const int LEnarrRezPointTraj, const double VAllamb, const int *IArrTargNums
                            , const int LEnTargNums, const double *ARrGrFgrCutting,const double VAlMomOut);

       void calcStationarySolution(const double *arrObjective
                       , double*arrStatinaryPhVect,double *valUd, double*valUq, double *pvalDelFi,const double VAlMomOut);

       int createInputDataReport(wchar_t*FileName, const bool bHeader);

       virtual int createInputDataReport_Ctrl(wchar_t*FileName, const bool bHeader);

};

#endif // ELECTDRIVER_H
