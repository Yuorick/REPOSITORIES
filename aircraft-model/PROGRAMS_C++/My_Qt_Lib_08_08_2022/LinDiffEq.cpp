#include "LinDiffEq.h"
#include "MatrixProccess.h"
#include "Comp.h"
#include "Equations.h"

#include <string.h>
#include <math.h>#include <math.h>
#include <float.h>


// вычисление фундамернтальной матрицы для лин диф уравнния
//  с постоянными коэффициентами
// INPUT:
//arrA[dimA*dimA] - матрица
//VAlT - время
// OUTPUT:
//arrF[dimA*dimA] - фундам матрица в момент T
void QLinDiffEq::calcFundamentalMtrx(const double *arrA, const int dimA, const double VAlT, double *arrF)
{
    double *parra = new double [dimA*dimA];
    double *parrT_n = new double [dimA*dimA];
    double *parrCur = new double [dimA*dimA];
    double *parrT_Sum = new double [dimA*dimA];

    double *parrT = new double [dimA*dimA];
    double *parrT1 = new double [dimA*dimA];
    memcpy(parra, arrA, dimA*dimA * sizeof(double));
    MatrxMultScalar(parra, dimA, dimA, VAlT,parrT);
    memcpy(parra, parrT, dimA*dimA * sizeof(double));

    fillE (parrCur,dimA);
    fillE (parrT_Sum,dimA);
    //memcpy(parrT_Sum,parrCur, dimA*dimA * sizeof(double));
  for (int i =0; i < 100; ++i)
  {
   MtrxMultMatrx(parrCur,dimA, dimA, parra,dimA, parrT1) ;
   MatrxMultScalar(parrT1, dimA, dimA, 1./((double)(i + 1.)),parrT_n);
   if(NormSquareVect(parrT_n, dimA *dimA)< DBL_MIN * 1000.)
   {
       break;
   }
   memcpy(parrCur, parrT_n, dimA*dimA * sizeof(double));
   MtrxSumMatrx(parrCur, parrT_Sum,dimA, dimA, parrT ) ;
   memcpy(parrT_Sum,parrT , dimA*dimA * sizeof(double));
  }
   memcpy(arrF,parrT_Sum, dimA*dimA * sizeof(double));
  delete []parra ;
  delete [] parrT_n ;
  delete [] parrCur ;
  delete [] parrT_Sum ;
  delete [] parrT ;
  delete [] parrT1 ;
}

//-----------------------------------------
// вычисление фундамернтальной матрицы для лин диф уравнния
// 2-го порядка с постоянными коэффициентами
// INPUT:
//arrA[4] - матрица
//VAlT - время
// OUTPUT:
//arrF[4] - фундам матрица в момент T
//решается характ уравнение, опредляются корни k1, k2
// вычисляется матрица arrF0] = exp(k1*T) , arrF3] = exp(k2*T)
void QLinDiffEq::calcFundamentalMtrx_dim2(const double *arrA,  const double VAlT, double *arrF)
{
    double arrT[4] ={0.};
    if ((fabs(arrA[1]) <0.00000001)&&(fabs(arrA[0] - arrA[3])< 0.0000001))
    {   double temp = exp(VAlT * arrA[0]) ;
        arrF[0] =arrF[3] = temp;
        arrF[2] = VAlT * arrA[2]*temp;
        return;
    }

    if (fabs(arrA[1]) <0.00000001)
    {
        double temp0 = exp(VAlT * arrA[0]) ;
        double temp1 = exp(VAlT * arrA[3]) ;
        double temp2 = arrA[2]/ (arrA[0] - arrA[3]);
        double arrB[4] = {0.},  arrB0Inv[4] = {0.};
        arrB[0] = temp0;
        arrB[1] = 0.;
        arrB[2] = temp0 *temp2;
        arrB[3] = temp1;
        arrB0Inv[0] = 1.;
        arrB0Inv[1] = 0.;
        arrB0Inv[2] = -temp2;
        arrB0Inv[3] = 1.;
        MtrxMultMatrx(arrB,2, 2, arrB0Inv,2, arrF) ;
        return;
    }

    const double a =1.;
    const double b = -(arrA[0] + arrA[3]);
    const double c =arrA[0] * arrA[3] - arrA[1] * arrA[2];
    TComp x1,x2;
    // возвращает
    // 0 - 2 действительных некраьных корня
    // 1 - 2 действительных кратных корня
    // 2- имеется по крайней мере один нулевой корень a!= 0, c=0
    // 3 - 2 комплексно сопряженных корня
    // 4 -  1 действительный корень (a =0)
    // 5  - несовместность a=b=0, c!=0
    // 6 - тождество )a=b=c=0
    int irez = SolvEq2( a, b, c,x1,x2);

    if (1 == irez) // кратные
    {
        double k = x1.m_Re;
        double val_exptk = exp(VAlT * k) ;
        double val_d = arrA[1] * arrA[2] - arrA[0]* arrA[3];
        double temp1 = arrA[3]/ arrA[1]  + val_d/ arrA[1]/ k ;
       double alf2 = -val_d/arrA[1]/k/k;
  double bet2 = arrA[3]/ arrA[1] + val_d/ arrA[1]/ k ;

        double arrB[4] = {0.},  arrB0Inv[4] = {0.};
        arrB[0] = val_exptk;
        arrB[1] = val_exptk * VAlT;
        arrB[2] = temp1 *val_exptk;
        arrB[3] = (alf2 + bet2 * VAlT) * val_exptk;
        arrB0Inv[0] =1. ;
        arrB0Inv[1] = 0.;
        arrB0Inv[2] = -temp1 /alf2;
        arrB0Inv[3] = 1./ alf2;
        MtrxMultMatrx(arrB,2, 2, arrB0Inv,2, arrF) ;
        return;
    }

    if (0 == irez)// не кратные
    {
        double k1 = x1.m_Re;
        double k2 = x2.m_Re;
        double val_exptk1 = exp(VAlT * k1) ;
        double val_exptk2 = exp(VAlT * k2) ;
        double val_d = arrA[1] * arrA[2] - arrA[0]* arrA[3];
        double temp1 = arrA[3]/ arrA[1]   + val_d/ arrA[1]/ k1 ;
        double temp2 = arrA[3]/ arrA[1]    + val_d/ arrA[1]/ k2 ;


        double arrB[4] = {0.}, arrB0Inv[4] = {0.};
        arrB[0] = val_exptk1;
        arrB[1] = val_exptk2;
        arrB[2] = temp1 *val_exptk1;
        arrB[3] = temp2 *val_exptk2;
        arrB0Inv[0] =temp2 / (temp2 - temp1) ;
        arrB0Inv[1] = -1./ (temp2 - temp1) ;
        arrB0Inv[2] = -temp1 / (temp2 - temp1) ;
        arrB0Inv[3] = 1./ (temp2 - temp1) ;
        MtrxMultMatrx(arrB,2, 2,arrB0Inv ,2, arrF) ;
        return;
    }

/*
    if ( (2 == irez)&&(fabs(arrA[1])< 0.00000001)&&(fabs(arrA[0])> 0.00000001)) //один нулевой
    {
        double k = x1.m_Re;
        double val_exptk = exp(VAlT * k) ;
        double temp0 = arrA[1]/ arrA[0] ;
        double arrB[4] = {0.}, arrB0Inv[4] = {0.};
        arrB[0] = val_exptk;
        arrB[1] = val_exptk;
        arrB[2] = temp0 *val_exptk;
        arrB[3] = temp0 *val_exptk +1;
        arrB0Inv[0] = 2. ;
        arrB0Inv[1] = -1 ;
        arrB0Inv[2] = -1;
        arrB0Inv[3] = 1. ;
        MtrxMultMatrx(arrB,2, 2,arrB0Inv ,2, arrF) ;
        return;
    }
*/
    if ( (2 == irez)&&(fabs(arrA[1])> 0.00000001)&&(fabs(arrA[0])< 0.00000001)) //один нулевой
    {
        double k = x1.m_Re;
        double val_exptk = exp(VAlT * k) ;
        double temp = arrA[1]/arrA[3];
        double arrB[4] = {0.}, arrB0Inv[4] = {0.};
        arrB[0] = temp *val_exptk;
        arrB[1] = temp *val_exptk +1;
        arrB[2] = val_exptk;
        arrB[3] = val_exptk ;
        arrB0Inv[0] = -1. ;
        arrB0Inv[1] = temp +1 ;
        arrB0Inv[2] =1. ;
        arrB0Inv[3] = -temp ;
        MtrxMultMatrx(arrB,2, 2,arrB0Inv ,2, arrF) ;
        return;
    }
    if (2 == irez) //один нулевой
    {

        double k = x1.m_Re;
        double val_exptk = exp(VAlT * k) ;


        double temp0 = arrA[3]/ arrA[1] ;
        double temp2 = -arrA[2]/ arrA[3]  ;


        double arrB[4] = {0.}, arrB0Inv[4] = {0.};
        arrB[0] = val_exptk;
        arrB[1] = 1.;
        arrB[2] = temp0 *val_exptk;
        arrB[3] = temp2 ;
        arrB0Inv[0] =temp2 / (temp2 - temp0) ;
        arrB0Inv[1] = -1./ (temp2 - temp0) ;
        arrB0Inv[2] = -temp0 / (temp2 - temp0) ;
        arrB0Inv[3] = 1./ (temp2 - temp0) ;
        MtrxMultMatrx(arrB,2, 2,arrB0Inv ,2, arrF) ;
        return;
    }

    if (3 == irez)
    {
        double k = x1.m_Re;
        double om = fabs(x1.m_Im);
         double val_d = arrA[1] * arrA[2] - arrA[0]* arrA[3];
        double val_exptk = exp(VAlT * k) ;
        double t1 = k * k + om * om;
        double c1 = arrA[3]/ arrA[1]   + val_d *k / arrA[1]/t1;
         double c2 = val_d/ arrA[1] *om/t1;

         double c3 = -om * val_d / arrA[1] /t1;
          double c4 = arrA[3]/ arrA[1]   +  val_d *k / arrA[1]/t1;


        double arrB[4] = {0.},  arrB0Inv[4] = {0.};
        arrB[0] =  cos (om *VAlT )*val_exptk;
        arrB[1] =  sin (om *VAlT )*val_exptk;
        arrB[2] = (c1 *cos (om *VAlT )+ c2* sin (om *VAlT ))*val_exptk;
        arrB[3] = (c3 *cos (om *VAlT )+ c4* sin (om *VAlT ))*val_exptk;
        arrB0Inv[0] = 1. ;
        arrB0Inv[1] = 0. ;
        arrB0Inv[2] = -c1 / c3;
        arrB0Inv[3] = 1./ c3 ;
        MtrxMultMatrx(arrB,2, 2,arrB0Inv ,2, arrF) ;
        return;

    }

}


