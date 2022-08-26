#ifndef HELICTRAJ_H
#define HELICTRAJ_H



#include <QVector>
#include "Helic.h"
#include "Environment.h"
#include "FlyTask.h"
#include "PartHelicTraj.h"
#include "TurnMove.h"
#include "LineMove.h"
#include "Hover.h"
#include "Rotating.h"


// к-во переменных фазового вектора вертолета в НЗСК
//#define  QUantCurNZSKVarsVS 6

// максимально возможное к-во участков траектории, которое можно сеебе предствавить на текущий момент
#define  MAX_POSSIBLE_TRAJPART_QUANT 100


#define  QUANT_COLS_OUT_VECT  30

// класс описывает операцию полета вертолета в соответствии с полетным заданием

class THelic;
class TEnvironment;
class TFlyTask;
class TComp ;
class TPartHelicTraj;
class TTurnMove;
class TLineMove;
class THover;
class TRotating;

class THelicTraj
{
public:
    //  1.вертолет
    THelic mHelic;
    // 2.атмосфера
    TEnvironment  mEnvironment;
    // 3.полетное задание
    TFlyTask mFlyTask;

    // 4.вектор координат нулевой (первой по счету) точки полетного задания в ГСК
    long double marr_GSKxyz_Pnt0[3];
    ///

    // 5, 6 долгота и широта нулевой (первой по счету) точки полетного задания в ГСК
    long double mLongitude; // долгота
    long double mLatitude; // широта


    // 7. фазовый вектор вертолета в НЗСК (положение и скорость)
    // marrPhaseVect[0]- marrPhaseVect[2]  - вектор центра масс в НЗСК
    // marrPhaseVect[3]- marrPhaseVect[5]  - вектор скорости центра масс в НЗСК


    long double marrPhaseVect[QUantCurNZSKVarsVS];
    // 8. время привязки траектории
    long double mTimeCur;

  // 9., 10, 11 таблица магнитного склонения
    QVector<double> mqvectTabMagDev;
    int mQuantRowsTabMagDev;
    int mQuantColsTabMagDev;
    ///



  // 12 номер текущей части траектории
   int mNumPart;

   // 13 массив номеров матриц передаточных чисел для участков траектории
  // int miarr[MAX_POSSIBLE_TRAJPART_QUANT];

   // 14 признак иницуиализации текущей части траектории
   bool mbPartInit;


  // 15 матрица перехода из текущей НЗСК в НЗСК
  long double marrTransf_CurNZSK_to_NZSK[9];
  // 16 матрица перехода из  НЗСК в текущую НЗСК
  long double marrTransf_NZSK_to_CurNZSK[9];
  // 17 вектор координат начала НЗСК
  long double marrO_CurNZSK[3];

  // 18 участок траектории
  TPartHelicTraj *mpPartTraj;
  TTurnMove mTurnMove;
  TLineMove mLineMove;
  THover mHover;
  TRotating mRotating;


  //19  шаг интегрирования
  long double mStepIntegr;



  THelicTraj() ;
// Конструктор копирования
  THelicTraj (const THelicTraj &R2) ;

 // оператор присваивания
 THelicTraj  operator=(THelicTraj  R2) ;

  // парам констр
 THelicTraj( const THelic Helic,const TEnvironment  Environment
             , const TFlyTask FlyTask, const long double Longitude
              ,const long double Latitude,const long double StepIntegr);

 //THelicTraj(const int ITypeof_Environment, const int ITYpeof_FLyTask
  //                      ,const bool bTypeofHelic, const bool bTypeofRotor) ;

 THelicTraj(const TEnvironment Environment , const int ITYpeof_FLyTask
                        ,const bool bTypeofHelic, const bool bTypeofRotor);

 void move(const double VAlFlyTime, const double VAlStIntegr, double *arrBuff, int *pquantSteps);

 void fillMagDevTab();

 double calcMagDevValue(const double VAlLongtitude, const double VAlLatitude);

 static double BilinearValue (double* arrSh,const double cellsize
                                          ,const double xt,const double yt);



//void assignNumberOfSetGears(const TFlyTask FlyTask, int *iarrNums);

void create_transfMtrx_NZSK_to_CurNZSK(const int INumPart);

void create_transfMtrx_CurNZSK_to_NZSK(const int INumPart);

void create_vectO(const int INumPart);

void transf_xyzNZSK_to_xyzCurNZSK(long double *arrNZSKInp, long double *arrCurNZSKOut);

void initPartTraj(const int NUmPart);

void imitateMoving_and_collectData( const long double VAlFlyTime ,const long double  VAl_OutPut_dt
                                    , double *arrBuff, const int QUantColsBuff
                                    , int *pquantSteps, double *arrOutSideInfo);

 bool goAhead_in_accordance_with_time(const long double VAlTend);

 void fillOutPutBuffRow(double *arrBuffRow);

 static long double SIGNUM(long double x);

 void transform_VS_from_CurNZSK_to_NZSK();

 void transf_xyzCurNZSK_to_xyzNZSK(long double *arrCurNZSKInp, long double *arrNZSKOut);

 void fillOutPutVect(double *arrOutputVect);

 static void transf_xyzNZSK_to_xyzGSK(long double *arrNZSKInp, long double *arrGSKOut);

 static void transf_xyzGSK_to_xyzNZSK(long double *arrGSKInp, long double *arrNZSKOut);

 static void transf_xyzNZSK_to_xyzGSK( double *arrNZSKInp,  double *arrGSKOut);

 static void transf_xyzGSK_to_xyzNZSK( double *arrGSKInp,  double *arrNZSKOut);

 long double recalcPsi_NZSK_to_CurNZSK(long double valPsi_NZSK);

 long double recalcPsi_CurNZSK_to_NZSK(long double valPsi_CurNZSK);

};


#endif // HELICTRAJ_H
