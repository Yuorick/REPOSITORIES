#ifndef ROTATING_H
#define ROTATING_H


#include "PartHelicTraj.h"
class TPartHelicTraj;

class TRotating : public TPartHelicTraj
{
public:

    // высота
    long double mY;

    // угловая скорость вращения
    long double m_dPsi_po_dt;

    // время вращения
    long double mTRotation;

    // начальный угол рыскания в НЗСК
    long double mPsi0;




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



    TRotating();
    TRotating  operator=(TRotating  R2) ;

    TRotating (const TRotating &R);
    // парам констр виража
  //  TRotating( const THelic Helic,const TEnvironment  Environment , const long double VAlY
        //                ,const long double  VAl_dPsi_po_dt , const long double VAlTRotation
         //               ,const long double  TBegin);

    TRotating( const THelic Helic,const TEnvironment  Environment , const long double VAlY
                        ,const long double  VAl_dPsi_po_dt , const long double VAlTRotation
                        ,const long double  TBegin,const long double  Psi0);



    virtual void calc_W(long double *arrW);



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


#endif // ROTATING_H
