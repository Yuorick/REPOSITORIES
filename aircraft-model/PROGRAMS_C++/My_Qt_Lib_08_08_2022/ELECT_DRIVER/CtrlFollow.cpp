#include "CtrlFollow.h"
#include "CtrlPos.h"
#include "Comp.h"

#include <math.h>
#include  <string.h>
#include <wchar.h>
#include <stdio.h>
#include <time.h>
#include "MatrixProccess.h"
#include "Equations.h"
#include "Comp.h"

QCtrlFollow::QCtrlFollow()
{

}

// конструктор копирования
 QCtrlFollow::QCtrlFollow(const QCtrlFollow &R):QCtrl( R)
 {

 }

 // оператор присваивания
  QCtrlFollow  &QCtrlFollow::operator=( const QCtrlFollow  &R)
 {
      if(this == &R)
      {
          return *this;
      }
     QCtrl:: operator= (R);

     return *this ;
 }
//-----------------------------------

// парам конструктор
  // параметрический конструктор
QCtrlFollow:: QCtrlFollow( const QElectMotor ElectMotor ,  QLoad*  Load
                           , const double *arrSpreadParams,const double VAlMomOut, const double VAlTettaBegin
                           ,const double VAlOmegaStatBegin , const double T0, const double TCur
                           , const double valh,TComp *pCmpArrLamb)
      :QCtrl (ElectMotor ,  Load, arrSpreadParams, VAlMomOut, VAlTettaBegin
               ,VAlOmegaStatBegin ,T0,  TCur, valh,pCmpArrLamb)
  {

  }

//------------------------------------------------------------------
//--------------------------------------------------------
void QCtrlFollow::CalcEigenValues( const double VAlTettaStat, const double VAlOmegaStat
                                 ,const double *arrC,const double VAlMomOut, TComp* cmparrEigenValues)
{


    double arrFGr[16] = {0.};

    calc_df_po_dx_plus_BC_Follow(VAlTettaStat, VAlOmegaStat, arrC,VAlMomOut, arrFGr);

    solvCharactEq_4(arrFGr, cmparrEigenValues);

}

//---------------------------------------------------------------
// вычисление матрицы передаточных чисел
//INPUT:
//arrObjective[2] - вектор подожения и скорости в состоянии равновесия (нужна только скоросчть)
//CmpArrLamb [4] - корни характерист уравнения? gthdst
// OUTPUT:
// arrC0[8] - матрица передаточных чисел
void QCtrlFollow::calcGearMtrx(const double *arrObjective, TComp* CmpArrLamb,const double VAlMomOut, double* arrC)
{

    //1. нахождение стационарного решения по регулировке положения

    double valUd=0.,  valUq=0.,  valDelFi=0., arrStatPhVect[QVARS ] = {0.};
    calcStationarySolution( arrObjective, arrStatPhVect,&valUd, &valUq, &valDelFi, VAlMomOut);

    double arr_dF_po_dx [QVARS  * QVARS ] = {0.} ;
    double arr_dF_po_dW [QVARS  * 2] = {0.};

   fill_df_po_px_and_mtrxB_(arrStatPhVect,arr_dF_po_dx, arr_dF_po_dW);

   double arrTemp[12]= {0.}, arrTemp0[9]= {0.}, arr_dF_po_dx_Pos[16] = {0.};
   for (int i = 0; i < 4; ++i)
       for (int j = 0; j < 4; ++j)
       {
         arr_dF_po_dx_Pos[ 4 * i + j] =  arr_dF_po_dx [ QVARS * i + j];
       }
   // вычисление d
   double arrd[8] = {0.};
  excludeRow(arr_dF_po_dx_Pos,4, 4, 1, arrTemp);
  excludeCol( arrTemp,3, 4, 1, arrTemp0);

  calcSubMtrx(arrTemp0, &(CmpArrLamb[1]), &(arrd[4]), &(arrd[6]), &(arrd[7]));
  // проверка
  double arrLamb[6];
  arrTemp0[3] += arrd[4];
 arrTemp0[4] += arrd[6];
 arrTemp0[5] += arrd[7];
  CalcProper_Numbers_R3(arrTemp0, arrLamb);

  ///

    arrd[0] =  -arr_dF_po_dx_Pos[4];//d1//
    arrd[1] =  CmpArrLamb[0].m_Re - arr_dF_po_dx_Pos[5];//d2
    arrd[2] = -arr_dF_po_dx_Pos[6];//d3//
    arrd[3] =  0.;//d4//
    arrd[5] = - arr_dF_po_dx_Pos[9];//d6

    ///
    /// \brief arrC
    memset(arrC, 0, 2 * QVARS * sizeof(double));
    double arrC0[8];
    MatrxMultScalar(arrd, 1, 8, 1./mElectMotor.mInvL,arrC0);
    memcpy(arrC, arrC0, 4 * sizeof(double));
    memcpy(&(arrC[QVARS]), &(arrC0[4]), 4 * sizeof(double));


}


void QCtrlFollow::calcCurU(const double *arrObjective,const double VAlTObjective,const double VAlMomOut, double *arrU)
{
    double arrDelta[QVARS] = {0.}, arrdelU[2] = {0.};
    double arrTargPhVect[QVARS] ={0.};
    QStatSolutionParams StatSolutionParams;
    calcStationaryParams(arrObjective, mCmpArrLamb, &StatSolutionParams,VAlMomOut);
    memcpy(arrTargPhVect, StatSolutionParams.marrStatPhVect, QVARS * sizeof(double));

    arrTargPhVect[3] = arrObjective[0];
    MtrxMinusMatrx(mFiltr.marrCurEst, arrTargPhVect,1,  QVARS, arrDelta);

    MtrxMultMatrx(StatSolutionParams.marrGears,2,  QVARS, arrDelta,1, arrdelU) ;
    MtrxSumMatrx(StatSolutionParams.marrStatU, arrdelU,1, 2, arrU) ;
    if(NormVect2(arrU) > mElectMotor.mUMax)
    {
        MatrxMultScalar(arrU, 1, 2, mElectMotor.mUMax/ NormVect2(arrU),arrU);
    }
}


//матрица частных производных правой части по фазовым переменным
// для первых 3-х уравнений _ угловая скорость, Id, Iq
void QCtrlFollow::calc_df_po_dx_plus_BC_Follow(const double VAlTettaStat, const double VAlOmegaStat
                      ,const double *arrC,const double VAlMomOut, double *arrFGr)
{
double arrStatPhVect[QVARS] = {0.}, arrUstat[2] = {0.};

double valUd=0.,  valUq=0.,  valDelFi=0., arrObjVect[2 ] ={0.};
arrObjVect[0] = VAlTettaStat;
arrObjVect[1] = VAlOmegaStat;

calcStationarySolution( arrObjVect, arrStatPhVect,&valUd, &valUq, &valDelFi, VAlMomOut);

arrUstat[0] = valUd;
arrUstat[1] = valUq;
double arr_dF_po_dx[16] = {0.},arr_dF_po_dW[8] = {0.}, arrT[16] = {0.}, arrCCC[8] ={0.};

fill_df_po_px_and_mtrxB_Follow(arrStatPhVect
                       ,arr_dF_po_dx,arr_dF_po_dW);

for (int i = 0; i < 2; ++i )
    for (int j =0; j < 4; j++)
{
  arrCCC[4 * i + j] = arrC[QVARS * i + j];
}

MtrxMultMatrx(arr_dF_po_dW,4, 2, arrCCC,4, arrT) ;
MtrxSumMatrx(arrT, arr_dF_po_dx, 4, 4, arrFGr) ;

}
// вычисление матрицы частных произыодных подсиситемы уравнен6ий по скорости
void QCtrlFollow::fill_df_po_px_and_mtrxB_Follow(const double *arrStatinaryPhVect
                 ,double *arr_dF_po_dx,double * arr_dF_po_dW)
{
    double arr_dF_po_dx_Big[QVARS * QVARS] = {0.}, arr_dF_po_dW_Big[2 * QVARS] = {0};
    fill_df_po_px_and_mtrxB_(arrStatinaryPhVect
                     ,arr_dF_po_dx_Big,arr_dF_po_dW_Big);
    memcpy(arr_dF_po_dW, arr_dF_po_dW_Big, 8 * sizeof(double));
    memcpy(&(arr_dF_po_dx[0]),&(arr_dF_po_dx_Big[0]), 4 * sizeof(double));
    memcpy(&(arr_dF_po_dx[4]),&(arr_dF_po_dx_Big[QVARS]), 4 * sizeof(double));
    memcpy(&(arr_dF_po_dx[8]),&(arr_dF_po_dx_Big[2 * QVARS]), 4 * sizeof(double));
    memcpy(&(arr_dF_po_dx[12]),&(arr_dF_po_dx_Big[3 * QVARS]), 4 * sizeof(double));
}
//---------------------------------------------------
//------------------------------------------
// задана матрица cпециального вида
//                A = {a11 a12 a13
//                     a21 a22 a23
//                       1 0   0 }
//задан массив корней CmpArrLamb[3] характ уравнения  матрицы
//                B = {a11     a12    a13
//                     a21+r1  a22+r2 a23+r2
//                     1       0      0 }
// r1, r2, r3 - неизвкстные числа
// требуется найти такие числа r1, r2, r3, чтобы характеристическо еуравнение матрицы B
// имело заданные коррни CmpArrLamb[3].
// CmpArrLamb[0] - действительное, CmpArrLamb[3].m_Im = 0
void QCtrlFollow::calcSubMtrx(const double *arrA, TComp* CmpArrLamb, double* pvalr1, double* pvalr2, double* pvalr3)
{
   TComp cmp0(0.,0.);
   cmp0 += CmpArrLamb[0];
   cmp0 += CmpArrLamb[1];
   cmp0 += CmpArrLamb[2];

   TComp cmp1(1.,0.);
   cmp1 *= CmpArrLamb[0];
   cmp1 *= CmpArrLamb[1];
   cmp1 *= CmpArrLamb[2];

   TComp cmpT12 = CmpArrLamb[0] *CmpArrLamb[1];
   TComp cmpT13 = CmpArrLamb[0] *CmpArrLamb[2];
   TComp cmpT23 = CmpArrLamb[1] *CmpArrLamb[2];

   TComp cmp2(0.,0.);
   cmp2 += cmpT12;
   cmp2 += cmpT13;
   cmp2 += cmpT23;
   double y = cmp0.m_Re - arrA[0];
   double x = (arrA[2] + arrA[0]* y - cmp2.m_Re)/ arrA[1];
   double z =  (cmp1.m_Re + arrA[2] * y )/ arrA[1];
   *pvalr1 = x - arrA[3];
   *pvalr2 = y - arrA[4];
   *pvalr3 = z - arrA[5];
}
//--------------------------------------------------
void QCtrlFollow::correctMomOut(const double VAlTettaZv
                          ,const double VAlDispTetta    , const double VAlW, const double VAlTcur)
{

}
//----------------------------------------
double QCtrlFollow::estimateMomOut(const double VAlTettaZv
     ,const double VAlDispTetta    , const double VAlW, const double VAlTcur)
{
return 0.;
}

//----------------------------------------
double QCtrlFollow::calcEstMomOut()
{
    return 0.;
}
//-----------------------------
int QCtrlFollow::createInputDataReport_Ctrl_Inherited(wchar_t*FileName, const bool bHeader)
{
    int len = wcslen(FileName) ;

    if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'x') // проверка, что
     && (FileName[len - 3] == 't') ) )  // указанный файл имеет расширение .flt
    {

      return 1 ;
    }

    FILE *fw ;

    if ((fw = _wfopen(FileName,L"a"))== NULL)

    {

     return 1 ;
    }
if (bHeader)
{
   fprintf(fw,"  Дата и время формирования отчета\n");
   time_t t = time(NULL);
   struct tm* aTm = localtime(&t);
   fprintf(fw,"  Год = %04d\n",aTm->tm_year+1900);
   fprintf(fw,"  Mесяц = %02d\n",aTm->tm_mon+1);
   fprintf(fw,"  День = %02d\n",aTm->tm_mday);
   fprintf(fw,"  Время = %02d:%02d:%02d\n",aTm->tm_hour, aTm->tm_min, aTm->tm_sec);
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
}
fprintf(fw,"***************************************\n");
fprintf(fw,"Тип системы управления - следящая\n");
fclose(fw);

createInputDataReport_Ctrl_Inherited(FileName, false);
return 0;
}
//-------------------------------------------------------
int QCtrlFollow::createInputDataReport_CtrlFollow_Inherited(wchar_t*FileName, const bool bHeader)
{
return 0;
}





