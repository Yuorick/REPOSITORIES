#include "liprog.h"
#include <iostream>
#include <engine.h> // MatLab include  !!!!!!



//#pragma link "mylibeng.lib"
//#pragma link "mylibmx.lib"


// nvars -число переменных
// f - вектор целевой функции
// nrows - число строк в матрице неравенств
// nrows_eq - число строк в матрице равенств
// a_eq - матрица равенст
// b_eq - вектор правой части равенст
// lb - вектор нижних границ переменных
// ub - вектор верхних границ переменных
// x - [output] вектор значений переменных
// fval - [output] оптимальное значение целевой функции

 // Пшеница озимая
//  Ячмень
//  Пшеница яровая
//  Рожь
//  Картофель
//  Многолетние травы
//  Однолетние травы
//  Кукуруза на силос
//  Кукуруза на зерно
//  Лен-долгунец
//  Сахарная свекла
//  Подсолнечник
//  Соя
//  Чистый пар

// Функция для MatLab
//#ifdef MATLAB
int linprog11( int nvars,  double *f, int nrows, double *a,  double *b
             ,int  nrows_eq,  double *a_eq,  double *b_eq
              ,double *lb, double *ub
            , double *x, double &fval)
{


    Engine *peng = NULL ;
    peng = engOpen(NULL);
    if(!peng)
    {
       //	AfxMessageBox(_T("Can't start Matlab engine ERROR!!!"),MB_OK) ;
     //	AfxAbort() ;
        return 1;
    }
    mxArray *mxf = mxCreateDoubleMatrix(nvars, 1, mxREAL);
    memcpy(mxGetPr(mxf), f, sizeof(double)*nvars);

    //s_22 = pBuf ;
   //  AfxMessageBox(s_22) ;

    engPutVariable(peng, "f", mxf);
    //s_22 = pBuf ;
 //    AfxMessageBox(s_22) ;

    mxArray *mxAeq = mxCreateDoubleMatrix(nvars, nrows_eq, mxREAL);
    memcpy(mxGetPr(mxAeq), a_eq, sizeof(double)*nrows_eq*nvars);
    engPutVariable(peng, "Aeq", mxAeq);
    engEvalString(peng, "Aeq = Aeq';");

    mxArray *mxbeq = mxCreateDoubleMatrix(nrows_eq, 1, mxREAL);
    memcpy(mxGetPr(mxbeq), b_eq, sizeof(double)*nrows_eq);
    engPutVariable(peng, "beq", mxbeq);

    mxArray *mxA = mxCreateDoubleMatrix(nvars, nrows, mxREAL);
    memcpy(mxGetPr(mxA), a, sizeof(double)*nrows*nvars);
    engPutVariable(peng, "A", mxA);
    engEvalString(peng, "A = A';");

    mxArray *mxb = mxCreateDoubleMatrix(nrows, 1, mxREAL);
    memcpy(mxGetPr(mxb), b, sizeof(double)*nrows);
    engPutVariable(peng, "b", mxb);


    mxArray *mxlb = mxCreateDoubleMatrix(nvars, 1, mxREAL);
    memcpy(mxGetPr(mxlb), lb, sizeof(double)*nvars);
    engPutVariable(peng, "lb", mxlb);

    mxArray *mxub = mxCreateDoubleMatrix(nvars, 1, mxREAL);
    memcpy(mxGetPr(mxub), ub, sizeof(double)*nvars);
    engPutVariable(peng, "ub", mxub);



    engEvalString(peng, "[x, fval] = linprog(f, A, b, Aeq, beq, lb,ub)");

    mxArray *mxFval = engGetVariable(peng, "fval");
    if ( mxFval == NULL)
    {
     //  AfxMessageBox(_T("ERROR mxFval == NULL !!!!"),MB_OK);
    }
    mxArray *mxX = engGetVariable(peng, "x");
    fval = *mxGetPr(mxFval);


    memcpy(x, mxGetPr(mxX), sizeof(double)*nvars);

    mxDestroyArray(mxf);
    mxDestroyArray(mxA);
    mxDestroyArray(mxb);

    mxDestroyArray(mxAeq);
    mxDestroyArray(mxbeq);
    mxDestroyArray(mxlb);
    mxDestroyArray(mxub);

    mxDestroyArray(mxFval);
    mxDestroyArray(mxX);
    engClose(peng);
    return 0;
}


