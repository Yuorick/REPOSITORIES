#ifndef PARTHELICTRAJ_H
#define PARTHELICTRAJ_H


#include "Helic.h"
#include "Environment.h"
//#include "FlyTask.h"

enum enumTypeOfManeuvre  {enLine, enTurn, enSpotHovering, enRotating};
// к-во переменных фазового вектора вертолета
#define  QUantCurNZSKVarsVS 13
// класс описывает операцию полета вертолета в соответствии с полетным заданием

class THelic;
class TEnvironment;
class TComp;



class TPartHelicTraj
{
public:
    //  вертолет
    THelic mHelic;
    // атмосфера
    TEnvironment  mEnvironment;



    // фазовый вектор вертолета
    // marrPhaseVect[0]- marrPhaseVect[2]  - вектор центра масс в ТНЗСК (текущей НЗСК)
    // marrPhaseVect[3]- marrPhaseVect[5]  - вектор скорости центра масс в ТНЗСК
    // marrPhaseVect[6]- marrPhaseVect[8]  - вектор угловой скорости в СвСК
    // marrPhaseVect[9]- marrPhaseVect[11] - углы крена, тангажа и рыскания соответственно
    // marrPhaseVect[12] - плотность атмосферы

    long double marrPhaseVect[QUantCurNZSKVarsVS];
    // время привязки траектории
    long double mTimeCur;

    // время начала траектории
    long double mTBegin;



    // вектор равнодействующей силы в СвСК
    long double marrSvSK_Force[3];




  TPartHelicTraj() ;
// Конструктор копирования
  TPartHelicTraj (const TPartHelicTraj &R2) ;

 // оператор присваивания
 TPartHelicTraj  operator=(TPartHelicTraj  R2) ;

  // парам констр
 TPartHelicTraj( const THelic Helic,const TEnvironment  Environment
                                 , const long double VAlTimeCur, const long double VAlTBegin);

 TPartHelicTraj( const THelic Helic,const TEnvironment  Environment , const long double VAlY
                                 , const long  double VAlVx, const long double VAlTimeCur, const long double VAlTBegin);

 void move(const double VAlFlyTime, const double VAlStIntegr, double *arrBuff, int *pquantSteps);



 void calcRightMemberOfDifEqSystem( const long double   valFi, const  long double   valKappa
                                    , const  long double   valEtta, const  long double   valDelFi
                                    , long double   * arrf, long double   *arrSvSK_Force, long double   *arrSvSK_Moment
                                   , long double   * arrSvSK_ShaftForce0, long double   *arrSvSK_ShaftMoment0
                                   , long double   *arrSvSK_AirForce0, long double   *arrSvSK_AirMoment0);


 static void calcRightPartForEilers(const long double gamma, const  long double nu
                                    , const  long  double psi, long  double *arrOmega, long  double *arrRightPart);

 //void calcCtrlParamsForVerticalTakeoff_And_TurnUp_New (const double VAlY,const double VAl_VeloY
   //              ,  double *arrCoefDifEq, double *arrEilers, double *arrCtrlParams);


 void calcUaSvSK_Case_NZSK(long double *arrUaSvSK);

 void calcCntrRulsPos_New(long double *arrCntrImpacts,long  double *pvalFi
                          ,long  double *pvalKappa ,long  double *pvalEtta,long  double *pvalDelFi);

 void calcMtrxTransf_SvSK_To_CurNZSK (long double* arrMtrxTransf);

 void calcMtrxTransf_CurNZSK_To_SvSK (long double* arrMtrxTransf);



// void calcGraph_FGr_ForSteadyState_Move(const long double VAlVx, const long double  VAlY
                       //            ,const long double   VAlRangNu, const int QUantPoints, double *arrBuff1);

 //long double fncFGr_HorizBalance(const long double VAlNu
         //                        ,long  double *pvalFRotX,long  double *pvalFRotY
          //                       ,long  double *arrSvSK_AirForce,long  double * arrSvSK_AirMoment);

// long double fncDerivF_HorizBalance(const long double VAlNu, const long  double VAl_dx);

// bool calcBalanceParams_ForSteadyState_Move(const long double VAlVx, const long double  VAlY
       //                                     , const long double VAlNu0 ,long double *arrCntrRulsPos
       //                                     ,long double *arrCntrImpacts,long double *pvalNu);

 void HorizontalSteadyMove(const long double VAlFlyTime, const long double VAlStIntegr
                          , double *arrBuff, int *pquantSteps);

 int makeStep(const long double valFi, const long  double valKappa
                           , const long  double valEtta, const long  double valDelFi
                           , const long   double dt, long  double *arrSvSK_Force, long double *arrSvSK_Moment
                          ,long  double * arrSvSK_ShaftForce0,long  double *arrSvSK_ShaftMoment0
                          ,long  double *arrSvSK_AirForce0,long  double *arrSvSK_AirMoment0);



 void create_df_po_dS(long  double *arrCntrRulsPos, long  double *arr_df_po_dS);

 void create_df_po_dW(long  double *arrCntrRulsPos, long  double *arr_df_po_dW);

 void calc_df_po_dSj(int j,long double del,long double *arrCntrRulsPos, long double *arr_df_po_dSj );

 void calc_df_po_dWj(int j,long double del,long double *arrCntrRulsPos,long double *arr_df_po_dWj );

 void HorizontalSteadyMove(const long double VAlY,const long double   VAlVx
                           ,long double *arrC, const int LEnarrC, const long double VAlFlyTime
                          , const long double VAlStIntegr, double *arrBuff, int *pquantSteps);

  void doAirDensity();

  bool replaceYourSelf(const long double VAlTend, const long double VAlStIntegr);

  void makeSingleStep( const long double VAlStIntegr);

  long double calcMach( long double valTay,long   double valVVozd);

  static long double SIGNUM(long double x);

//  void calc_dFa_and_dMa_po_dxj(const long double  Val_dxj, const int j,long double  *arrFa
   //               ,long double  *arrMa,long double  * arr_dFa_po_dxj,long double  *arr_dMa_po_dxj);

  static  bool  IsStability_(const long double z1 ,const long double  z2
                                     ,const long double val_A,long double*arrA,const  int dimArrA);

  static  bool IsStability(const long double z1 ,const long double  z2
                                   ,const long double val_A,long double*arrA,const  int dimArrA);

  static bool  IsRootsSuit_dim5(long double*arrA,const  int dimA, const  double VAlMinX
                         ,const  double VAlMaxX, const  double VAlMaxTang, TComp *cmparrRoots);

  //void fill_df_po_px_and_df_po_dW_LittleTask(const long double VAlY ,   const long double VAlVx
   //                                                         , long double *arr_dF_po_dx, long double *arr_dF_po_dW);
 // void fill_df_po_px_and_df_po_dW_BigTask(const long double VAlY ,   const long double VAlVx
      //                                                   , long double *arr_dF_po_dx, long double *arr_dF_po_dW);

  bool calcMtrxTransf_SvSK_To_SkorSK (long double* arrMtrxTr);

  bool calcMtrxTransf_SkorSK_To_SvSK (long double* arrMtrxTr);

  long double calcPeregrNorm();

  static void calcMtrxTransf_SvSK_To_CurNZSK_ (const long double gamma
                                              ,const long double nu,const long double psi, long double* arrMtrxTransf);

  bool TryHorizontalTurn(const long double VAlR, const long double VAlPsi, const long double VAlFlyTime, const long double VAlStIntegr
                           , double *arrBuff, int *pquantSteps);

  void fill_df_po_px_and_df_po_dW_Turn(const long double VAlRad, long double *arrW
                                                         , long double *arr_dF_po_dx, long double *arr_dF_po_dW);

   static void fillUpOmega_StableTurn(const long double VAlV,const long double VAlRad,const long double VAlGamma
                               ,const long double VAlNu, long double *arrOmega);

   void fill_df_po_px_and_df_po_dW_Turn_Little(const long double VAlRad, long double *arrW
                                                          , long double *arr_dF_po_dx, long double *arr_dF_po_dW);

   virtual void calc_W(long double *arrW);

   void Moving( const long double VAlFlyTime, const long double VAlStIntegr
                ,const long double VAlModelStep, double *arrBuff, int *pquantSteps);

   void fill_df_po_px_and_df_po_dW_TangCanale(long double *arr_dF_po_dx, long double *arr_dF_po_dW);

   virtual void get_arrSteadyW(long double *arrW);

   void fill_df_po_px_and_df_po_dW( long double *arr_dF_po_dx, long double *arr_dF_po_dW);

   virtual void get_arrayOfControlledVars(int *piNum, int *iarr);

   virtual void get_arrayOfControlledVarsTang(int *piNum, int *iarr);

   virtual void get_QuantOfControlledVarsTang(int *piNumr);

   virtual void get_QuantOfControlledVars(int *piNumr);

   virtual void findingCircleTang(long double valz1 ,  long double valz2, long double val_A, long double *arr_dF_po_dx, long double *arr_dF_po_dW
                                  ,double *arrDataBuf ,const int maxQuant, int *pquantRows);

   void selectGearsTang(long double valz1 ,  long double valz2, long double val_A
                        ,double *arrDataBuf ,const int maxQuant, int *quantRows);

   static void insertLocalGearMtrx_In_TotalGearMtrx(const int QuantControlledVars, int *iarrNumsControlledVars
                                                                        ,const int QuantControls, int *iarrNumsControls
                                                                        ,long double *arrLocalGearMtrx, double *arrRezTotalGearMtrx);

   static void extractLocalGearMtrx_From_TotalGearMtrx(double *arrInpTotalGearMtrx,const int QuantControlledVars, int *iarrNumsControlledVars
                                                                       ,const int QuantControls, int *iarrNumsControls
                                                                       ,long double *arrOutLocalGearMtrx );

   void selectGears(long double valz1 ,  long double valz2, long double val_A
                                   , double *arrC, double *arrDataBuf ,const int maxQuant, int *quantRows);

   virtual void findingCircle(long double valz1 ,  long double valz2, long double val_A
                                              , long double *arrC, long double *arr_dF_po_dx, long double *arr_dF_po_dW
                                             ,double *arrDataBuf ,const int maxQuant, int *pquantRows);

   virtual  bool IsEndOfPart();

   static  double  fncMod2Pi__( double a_fVal );

   static double  fnc_Minus_PI_Plus_PI( double a_fVal );


};


#endif // PARTHELICTRAJ_H
