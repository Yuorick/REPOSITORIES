#ifndef HOVER_H
#define HOVER_H
//#include "LineMove.h"
#include "PartHelicTraj.h"

#include "PartHelicTraj.h"
class TPartHelicTraj;

class THover : public TPartHelicTraj
{
public:

    // высота
    long double mY;

    long double mTHovering;


   // long double mPsi0;




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



    THover();
    THover  operator=(THover  R2) ;

    THover (const THover &R);
    // парам констр виража
    THover( const THelic Helic,const TEnvironment  Environment , const long double VAlY
             , const long double VAlTHovering ,const long double  TBegin);



    virtual void calc_W(long double *arrW);

   // static double  fncMod2Pi__( double a_fVal );

   // static void  calcMtrxPer_CurNZSK_to_RotatedSK (const long double val_Ang_OX_Rot, long double *arrMtrxPer_CurNZSK_to_RotatedSK);

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

/*
class TLineMove;
class TPartHelicTraj;
class THover : public TLineMove
{
public:
    long double mTHovering;
   // long double mPsi0;

    THover();
    THover( const THelic Helic,const TEnvironment  Environment , const long double VAlY
                                    , const long  double VAlVx, const long double VAlAngYaw
                    , const long double VAlTHovering, const long double VAlTBegin);

    THover  operator=(THover  R2) ;
    THover (const THover &R);
    virtual  bool IsEndOfPart();
};
*/
#endif // HOVER_H
