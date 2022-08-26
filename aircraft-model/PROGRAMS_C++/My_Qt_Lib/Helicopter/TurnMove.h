#ifndef TURNMOVE_H
#define TURNMOVE_H

#include "PartHelicTraj.h"
class TPartHelicTraj;

class TTurnMove : public TPartHelicTraj
{
public:

    long double mY;
    long  double mVx;
    long double mAngYaw;
    long double mRadius;

    // координаты точки центра вращения в ТНЗСК
    long double marrRotCentre_CurNZSK[3];

    // вектор номеров управляемых переменных
    int mQuantControlledVars;
    int miarrNumsControlledVars[QUantCurNZSKVarsVS -1];


    // вектор номеров управляемых для задачи управления по тангажу
    int mQuantControlledVarsTang;
    int miarrNumsControlledVarsTang[QUantCurNZSKVarsVS -1];

    // матрица коэффициентов передачи (48)
    long double marrC0[(QUantCurNZSKVarsVS -1) * 4];

    // установившееся решение
    // вектор установившегося решения
    long double marrSteadySolution0[QUantCurNZSKVarsVS -1];

    // установившиеся положения рулей - вектор положения управляющих рулей (общ шаг, каппа, этта, педали)
    long double marrSteadyW[4];

    // конечный  угол курса в ТНЗСК
    long double mCourseEnd;

    // текущий интегрированный курс за время движения (без пересчета от -пи до пи)
    // нужен для движения по окружности и спирали
    long double mCurrentCourseIntegrated;
    TTurnMove();
    TTurnMove  operator=(TTurnMove  R2) ;

    TTurnMove (const TTurnMove &R);
    // парам констр виража
    TTurnMove( const THelic Helic,const TEnvironment  Environment , const long double VAlY
               , const long  double VAlVx, const long double VAlAngYaw
     , const long double VAlRadius,long double *arrRotCentre_CurNZSK
     ,const long double  TBegin,const long double  CourseEnd);



    virtual void calc_W(long double *arrW);

   // static double  fncMod2Pi__( double a_fVal );

    static void  calcMtrxPer_CurNZSK_to_RotatedSK (const long double val_Ang_OX_Rot, long double *arrMtrxPer_CurNZSK_to_RotatedSK);

    virtual void get_arrSteadyW(long double *arrW);

    virtual void get_arrayOfControlledVars(int *piNum, int *iarr);

    virtual void get_arrayOfControlledVarsTang(int *piNum, int *iarr);

    virtual void get_QuantOfControlledVarsTang(int *piNumr);

    virtual void get_QuantOfControlledVars(int *piNumr);

    virtual void findingCircleTang(long double valz1 ,  long double valz2, long double val_A, long double *arr_dF_po_dx, long double *arr_dF_po_dW
                                        ,double *arrDataBuf ,const int maxQuant, int *pquantRows);

    virtual void findingCircle(long double valz1 ,  long double valz2, long double val_A
                                          , long double *arrC, long double *arr_dF_po_dx, long double *arr_dF_po_dW
                                         ,double *arrDataBuf ,const int maxQuant, int *pquantRows);



  //  static double  fnc_Minus_PI_Plus_PI( double a_fVal );

    virtual  bool IsEndOfPart();

    void calcCurrentIntegtatedCourse();

};

#endif // TURNMOVE_H
