#ifndef LINEMOVE_H
#define LINEMOVE_H
#include "PartHelicTraj.h"
class TPartHelicTraj;

class TLineMove : public TPartHelicTraj
{
public:

    long double mY;
    long  double mVx;
    long double mAngYaw;
    // матрица коэффициентов передачи (44)
    long double marrC0[(QUantCurNZSKVarsVS -1) * 4];

    // установившееся решение
    // вектор установившегося решения
    long double marrSteadySolution0[QUantCurNZSKVarsVS -1];

    // установившиеся положения рулей - вектор положения управляющих рулей ( каппа, общ шаг,этта, педали)
    long double marrSteadyW[4];

    // вектор номеров управляемых переменных
    int mQuantControlledVars;
    int miarrNumsControlledVars[QUantCurNZSKVarsVS -1];


    // вектор номеров управляемых для задачи управления по тангажу
    int mQuantControlledVarsTang;
    int miarrNumsControlledVarsTang[QUantCurNZSKVarsVS -1];
    TLineMove();
    TLineMove  operator=(TLineMove  R2) ;

    TLineMove (const TLineMove &R);
    // парам констр равном прямолин движения
    TLineMove( const THelic Helic,const TEnvironment  Environment , const long double VAlY
                                    , const long  double VAlVx, const long double VAlAngYaw,const long double  TBegin);

    virtual void calc_W(long double *arrW);

    virtual void get_arrayOfControlledVars(int *piNum, int *iarr);

    virtual void get_arrayOfControlledVarsTang(int *piNum, int *iarr);

    virtual  void get_arrSteadyW(long double *arrW);

    virtual void get_QuantOfControlledVarsTang(int *piNumr);

    virtual void get_QuantOfControlledVars(int *piNum);

    virtual void findingCircleTang(long double valz1 ,  long double valz2, long double val_A, long double *arr_dF_po_dx, long double *arr_dF_po_dW
                                   ,double *arrDataBuf ,const int maxQuant, int *pquantRows);

    virtual void findingCircle(long double valz1 ,  long double valz2, long double val_A
                                          , long double *arrC, long double *arr_dF_po_dx, long double *arr_dF_po_dW
                                         ,double *arrDataBuf ,const int maxQuant, int *pquantRows);

    virtual  bool IsEndOfPart();
};

#endif // LINEMOVE_H
