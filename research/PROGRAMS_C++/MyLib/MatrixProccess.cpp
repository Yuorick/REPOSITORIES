//---------------------------------------------------------------------------


#pragma hdrstop

#include "MatrixProccess.h"
#include  <string.h>
#include <math.h>
#include <float.h>
#include "Equations.h"
//#include <vcl.h>
#include "Comp.h"
 const double  EPS =  0.00000001 ;

//---------------------------------------------------------------------------
// ��������� ������� �� �������
void MtrxMultMatrx(double *parrA,int nRowsA, int nColsA, double * parrB,int nColsB, double *parrRez)
{
  memset (parrRez,0,sizeof(double) * nRowsA *  nColsB) ;
 double *pRez = parrRez ;
  double *pA = parrA ;
  double *pB = parrB ;

  for (int i = 0; i < nRowsA; i++)
  {
	double *pAt = pA;
   double *pBt =pB ;
  for (int j = 0; j < nColsB; j++)
  {
	pAt = pA ;
	pBt = parrB + j ;

	for (int k = 0; k < nColsA ; k++)
	{
	  *pRez += (*pAt) * (*pBt);
	  pAt ++;
	  pBt += nColsB ;
	}
	pRez++;
  }
  pA += nColsA ;
  }

}

// ��������� ������� �� ����������������� �������
void MtrxMultMatrxTransp(double *parrA,int nRowsA, int nColsA, double * parrB,int nRowsB, double *parrRez)
{
/*  memset (parrRez,0,sizeof(double) * nRowsA *  nRowsB) ;
 for (int i = 0; i < nRowsA; i++)
  {
	 for (int j = 0; j < nRowsB; j++)
	 {
		int ktemp = i* nRowsB + j ;
		for (int k = 0; k < nColsA ; k++)
		{
		  parrRez [ktemp] += parrA[i * nColsA + k]*  parrB [ j *  nColsA + k ] ;
		}
	 }
  } */

   memset (parrRez,0,sizeof(double) * nRowsA *  nRowsB) ;
  double *pRez = parrRez ;
  double *pA = parrA ;
  double *pB = parrB ;

  for (int i = 0; i < nRowsA; i++)
  {
	double *pAt = pA;
   double *pBt =pB ;
  for (int j = 0; j < nRowsB; j++)
  {
	pAt = pA ;


	for (int k = 0; k < nColsA ; k++)
	{
	  *pRez += (*pAt) * (*pBt);
	  pAt ++;
	  pBt ++;
	}
	pRez++;
  }
  pA += nColsA ;
  }
}
// ��������� ����������������� ������� ��  �������
void MtrxTranspMultMatrx(double *parrA,int nRowsA, int nColsA, double * parrB,int nColsB, double *parrRez)
{
  memset (parrRez,0,sizeof(double) * nColsA *  nColsB) ;
 for (int i = 0; i < nColsA ; i++)
  {
	 for (int j = 0; j < nColsB; j++)
	 {
		int ktemp = i* nColsB + j ;
		for (int k = 0; k < nRowsA ; k++)
		{
		  parrRez [ktemp] += parrA[k * nColsA + i]*  parrB [ k *  nColsB + j ] ;
		}
	 }
  }
  /*    memset (parrRez,0,sizeof(double) * nColsA *  nColsB) ;

  double *pRez = parrRez ;
  double *pA = parrA ;
  double *pB = parrB ;

  for (int i = 0; i < nColsA; i++)
  {
	//double *pAt = parrA + i;

  for (int j = 0; j < nColsB; j++)
  {
	double *pBt =parrB + j ;
	double *pAt = parrA + i;

	for (int k = 0; k < nRowsA ; k++)
	{
	  *pRez += (*pAt) * (*pBt);
	  pAt += nColsA ;
	  pBt += nColsB ;
	}
	pRez++;
  }

  }*/
}
void MtrxSumMatrx(double *pA, double * pB,int nRows, int nCols, double *pRez)
{
	double *parrA = pA;
	double * parrB = pB ;
	double *parrRez = pRez;
	for (int i = 0; i < nRows * nCols; i++)
	{
	  *parrRez = (*parrA)+ (* parrB);
	  parrRez++;
	  parrA++ ;
	  parrB++ ;

	}
}
void MtrxMinusMatrx(double *pA, double * pB,int nRows, int nCols, double *pRez)
{
	double *parrA = pA;
	double * parrB = pB ;
	double *parrRez = pRez;

	for (int i = 0; i < nRows * nCols; i++)
	{
	  *parrRez = (*parrA)- (* parrB);
	  parrRez++;
	  parrA++ ;
	  parrB++ ;

	}
}

void MatrxMultScalar(double *parrA, int nRows, int nCols, double valScal,double *parrRez)
{

   for (int i = 0; i < nRows * nCols; i++)
	{
	  parrRez [i] =  parrA [i] * valScal ;
	}
}
void MatrxDivideScalar(double *parrA, int nRows, int nCols, double valScal,double *parrRez)
{
   for (int i = 0; i < nRows * nCols; i++)
	{
	  *parrRez = (*parrA)/ valScal ;
	  parrRez++;
	  parrA++ ;


	}
}

//--------------------------------------------------------------------------------------
void MatrTransp(double *parrA, int nRows, int nCols, double *parrRez)
{
	double *pA = parrA ;
	for (int i = 0; i < nRows; i++)
	{

	double *pRezt =parrRez + i ;
		for (int j = 0; j < nCols; j++)
		{
		*pRezt = *pA;
		pA ++;
		pRezt += nRows ;
		}

	}
}
//--------------------------------------------------------------------------------------
//arrOut= arrInp0 * arrInp1 * arrInp2
// ������� ����������
void MtrxMultMatrx_MultMatrx(double *parrInp0,double *parrInp1,double *parrInp2,int nDim, double *parrOut)
{
	double *arrT  = new double [nDim * nDim];
	MtrxMultMatrx(parrInp0,nDim, nDim, parrInp1,nDim, arrT ) ;
	MtrxMultMatrx(arrT,nDim, nDim, parrInp2,nDim, parrOut ) ;
	delete arrT;
}

//--------------------------------------------------------------------------------------
//arrOut= arrInp0 * arrInp1 * arrInp2Transp
// ������� ����������
void MtrxMultMatrx_MultMatrxTransp(double *parrInp0,double *parrInp1,double *parrInp2,int nDim, double *parrOut)
{
	double *arrT  = new double [nDim * nDim];
	MtrxMultMatrx(parrInp0,nDim, nDim, parrInp1,nDim, arrT ) ;
	MtrxMultMatrxTransp(arrT,nDim, nDim, parrInp2,nDim, parrOut ) ;
	delete arrT;
}


//--------------------------------------------------------------------------------------

//arrOut= parrFInp * parrDInp * parrFInpTransp
// ������� ����������
void calcF_D_FTransp(double *parrFInp,double *parrDInp,int nDim, double *parrOut)
{
	double *arrT  = new double [nDim * nDim];
	MtrxMultMatrx(parrFInp,nDim, nDim, parrDInp,nDim, arrT ) ;
	MtrxMultMatrxTransp(arrT,nDim, nDim, parrFInp,nDim, parrOut ) ;
	delete arrT;
}


//--------------------------------------------------------------------------------------


// ��������� ������������
void OuterProduct(double *pVect0 , double *pVect1, double *pVectRez)
{
 pVectRez [0] = pVect0 [1] * pVect1 [2] -  pVect0 [2] * pVect1 [1] ;
 pVectRez [1] = pVect0 [2] * pVect1 [0] -  pVect0 [0] * pVect1 [2] ;
 pVectRez [2] = pVect0 [0] * pVect1 [1] -  pVect0 [1] * pVect1 [0] ;
}

// ���������  ������������
double ScalProduct(double *pVect0 , double *pVect1, const int len)
{
 double prod = 0. ;
 for (int i = 0; i < len; i++)
 {
  prod += pVect0[i] * pVect1[i] ;
 }

 return prod ;
}

double Norm3( double *arrA)
{
	return sqrt( arrA[0] * arrA[0] + arrA[1] * arrA[1] + arrA[2] * arrA[2]) ;
}

// ��������� ������������ 2-� ������ �������� - ������ �������
// >0 ���� ���������� ������� �� ������� ������� �� ������� �����
// ������ ������� �������  (���� ������� ������)
double OuterProduct_2(double *pVect0 , double *pVect1)
{
 return pVect0 [0] * pVect1[1] - pVect0 [1] * pVect1[0];
}

// ������������ �������� ���������� ������� ���� double
// INPUT:
// parrDoubleInp - ������ ���� double
// lenArrayInp - ����� ������� lenArrayInp
// OUTPUT:
// pNumArgMin - ����� ��������� �������, �� ������� ����������� �������
// ������� ���������� ����������� ��������
double MinDoubleArray(double *parrDoubleInp, const int lenArrayInp, int *pNumArgMin)
{
	*pNumArgMin = -1;
	double valRez = DBL_MAX ;
	for (int i = 0; i < lenArrayInp; i++)
	{
	if (parrDoubleInp [i] < valRez)
	{
	 *pNumArgMin = i ;
	 valRez = parrDoubleInp [i] ;
	}
	}
	return valRez ;
}

// ������������ �������� ���������� ������� ���� double
// INPUT:
// parrDoubleInp - ������ ���� double
// lenArrayInp - ����� ������� lenArrayInp
// OUTPUT:
// pNumArgax - ����� ��������� �������, �� ������� ����������� �������
// ������� ���������� ������������ ��������
double MaxDoubleArray(double *parrDoubleInp, const int lenArrayInp, int *pNumArgMax)
{
  *pNumArgMax = -1;
	double valRez = -DBL_MAX ;
  for (int i = 0; i < lenArrayInp; i++)
  {
	if (parrDoubleInp [i] > valRez)
	{
	 *pNumArgMax = i ;
	 valRez = parrDoubleInp [i] ;
	}
  }
  return valRez ;
}

//---------------------------------------------------------------------------
// ������������ �������� ���������� ������� ���� int
// INPUT:
// parrDoubleInp - ������ ���� double
// lenArrayInp - ����� ������� lenArrayInp
// OUTPUT:
// pNumArgax - ����� ��������� �������, �� ������� ����������� �������
// ������� ���������� ������������ ��������
int MaxIntArray(int *parrDoubleInp, const int lenArrayInp, int *pNumArgMax)
{
  *pNumArgMax = -1;
	int valRez = -1000000 ;
  for (int i = 0; i < lenArrayInp; i++)
  {
	if (parrDoubleInp [i] > valRez)
	{
	 *pNumArgMax = i ;
	 valRez = parrDoubleInp [i] ;
	}
  }
  return valRez ;
}

// ���������� �������� �������  3x3
// ���� ������������  arrInp != ����, �� True
// � ��������� ������ false
 bool   InverseMtrx3(double *arrInp, double *arrOut)
{
  // ���������� ������������
  double det = arrInp[0] * ( arrInp[4] * arrInp[8] - arrInp[5] * arrInp[7] )
			-  arrInp[1] * ( arrInp[3] * arrInp[8] - arrInp[5] * arrInp[6] )
			+  arrInp[2] * ( arrInp[3] * arrInp[7] - arrInp[4] * arrInp[6] )  ;
  ///
  if (fabs (det) < 0.0000000000000001) return false ;

  double arr0[9]= {0};

  arr0[0]=  ( arrInp[4] * arrInp[8] - arrInp[5] * arrInp[7] ) ;
  arr0[3]= -( arrInp[3] * arrInp[8] - arrInp[5] * arrInp[6] ) ;
  arr0[6]=  ( arrInp[3] * arrInp[7] - arrInp[4] * arrInp[6] ) ;
  arr0[1]= -( arrInp[1] * arrInp[8] - arrInp[2] * arrInp[7] ) ;
  arr0[4]=  ( arrInp[0] * arrInp[8] - arrInp[2] * arrInp[6] ) ;
  arr0[7]= -( arrInp[0] * arrInp[7] - arrInp[1] * arrInp[6] ) ;
  arr0[2]=  ( arrInp[1] * arrInp[5] - arrInp[2] * arrInp[4] ) ;
  arr0[5]= -( arrInp[0] * arrInp[5] - arrInp[2] * arrInp[3] ) ;
  arr0[8]=  ( arrInp[0] * arrInp[4] - arrInp[1] * arrInp[3] ) ;

  MatrxMultScalar(arr0, 3, 3, 1./ det, arrOut);



  return true ;
}
 // ���������� ������������
double   calcDet3( double *arrInp)
{
 return arrInp[0] * ( arrInp[4] * arrInp[8] - arrInp[5] * arrInp[7] )
			-  arrInp[1] * ( arrInp[3] * arrInp[8] - arrInp[5] * arrInp[6] )
			+  arrInp[2] * ( arrInp[3] * arrInp[7] - arrInp[4] * arrInp[6] )  ;
}

// ������� ������ �������� ��������� 6-�� �������� ��������� ����������
bool   SolvFrobenius6(double *arrInpMtrxA00,double *arrInpMtrxA01
	   ,double *arrInpMtrxA10,double *arrInpMtrxA11, double *arrInpVectB,double *arrOutX)
{
   long	double arrMtrxAInv[36] = {0.} ;
   long		double arrT00[9] = {0.},arrT01[9] = {0.},arrT10[9] = {0.}
	   ,arrK[9] = {0.},arrKInv[9] = {0.}, arrA00Inv[9] ={0.};
   long		double arrT0[9] = {0.},arrT1[9] = {0.},arrT2[9] = {0.},arrT3[9] = {0.},arrT4[9] = {0.},arrT5[9] = {0.};
   long	double  arrMtrxA00[9] ={0.},arrMtrxA01[9] ={0.},arrMtrxA10[9] ={0.},arrMtrxA11[9] ={0.}, arrVectB [6] = {0.};
   for (int i = 0; i < 9; i++)
   {
	 arrMtrxA00[i] = arrInpMtrxA00[i] ;
	 arrMtrxA10[i] = arrInpMtrxA10[i] ;
	 arrMtrxA01[i] = arrInpMtrxA01[i] ;
	 arrMtrxA11[i] = arrInpMtrxA11[i] ;
   }
   for (int i = 0; i < 6; i++)
   {
	 arrVectB[i] = arrInpVectB [i] ;
   }

	if(!InverseMtrx3(arrMtrxA00, arrA00Inv)) return false ;

	// ���������� ������ �
	MtrxMultMatrx(arrMtrxA10,3, 3, arrA00Inv, 3, arrT0) ;
	MtrxMultMatrx(arrT0 ,3, 3, arrMtrxA01, 3, arrT1) ;
	MtrxMinusMatrx(arrMtrxA11, arrT1, 3, 3, arrK);
	if(!InverseMtrx3(arrK, arrKInv)) return false ;
	///

	// ����������  arrT00
	 MtrxMultMatrx(arrA00Inv,3, 3, arrMtrxA01, 3, arrT0) ; //arrT0 =  arrA00Inv * arrMtrxA01
	 MtrxMultMatrx(arrT0,3, 3, arrKInv, 3, arrT1) ;       // arrT1 = arrA00Inv * arrMtrxA01 * arrKInv
	 MtrxMultMatrx(arrT1,3, 3, arrMtrxA10, 3, arrT2) ;    //arrT2 = arrA00Inv * arrMtrxA01 * arrKInv *  arrMtrxA10
	 MtrxMultMatrx(arrT2,3, 3, arrA00Inv, 3, arrT3) ;  //arrT3 = arrA00Inv * arrMtrxA01 * arrKInv *  arrMtrxA10 * arrA00Inv
	 MtrxSumMatrx(arrA00Inv, arrT3,3, 3, arrT00);
	 ///

	 // ���������� arrT01
	 MatrxMultScalar(arrT1, 3, 3, -1.,arrT01);
	 ///

	 // arrT10
	 MtrxMultMatrx(arrKInv,3, 3, arrMtrxA10, 3, arrT4) ;
	 MtrxMultMatrx(arrT4,3, 3, arrA00Inv, 3, arrT5) ;
	 MatrxMultScalar(arrT5, 3, 3, -1.,arrT10);
	 ///

	 // ������������ �������� �������
	 for (int i = 0; i < 3; i++)
	 for (int j = 0; j < 3; j++)
	 {
	   arrMtrxAInv [ 6 * i + j]           = arrT00  [ 3 * i + j] ;
	   arrMtrxAInv [ 6 * i + 3 + j]       = arrT01  [ 3 * i + j] ;
	   arrMtrxAInv [ 6 * (i + 3) + j]     = arrT10  [ 3 * i + j] ;
	   arrMtrxAInv [ 6 * (i + 3) + 3 + j] = arrKInv [ 3 * i + j] ;
	 }
	 ///

	 // ������� ������� ��������� ���������
	 long double arrX[6] = {0.} ;
	  MtrxMultMatrx(arrMtrxAInv, 6, 6, arrVectB, 1, arrX) ;
	  for (int i = 0; i < 6; i++)
	  {
		arrOutX [i] = arrX [i] ;
	  }


	return true ;
}

// ���������� YT * D * Y ��� ������� 3�3
double  calcYT_D_Y(double *arrY, double *arrD )
{
  double arrT0[3] = {0.} ;
  double valReturn = 0 ;
  MtrxMultMatrx(arrY,1, 3, arrD,3, arrT0) ;
  MtrxMultMatrx(arrT0,1, 3, arrY,1, &valReturn) ;
  return valReturn ;
}

// ���������� YT * D * Y ��� ������� IDim�IDim
double  calcYT_D_Y(double *arrY, double *arrD, const int IDim )
{
	double *arrT0 = new double [IDim];
  double valReturn = 0 ;
	MtrxMultMatrx(arrY,1, IDim, arrD,IDim, arrT0) ;
	MtrxMultMatrx(arrT0,1, IDim, arrY,1, &valReturn) ;
	delete arrT0 ;
  return valReturn ;
}

 // ����� ����������� ��������  � ����������� ����� ������������ ������������ ������� arrKInp 3� 3
//  arrV - ������� ����������� ��������
//  arrLamb - ��������������� ������� ����������� �����  3 � 3
//  Returns:
// 31 - ���� arrKInp = 3, �-�� ������ ����� = 1
// 32 - ���� arrKInp = 3, �-�� ������ ����� = 2
// 33 - ���� arrKInp = 3, �-�� ������ ����� = 3
// 22 - ���� arrKInp = 2, �-�� ������ ����� = 2
// 23 - ���� arrKInp = 2, �-�� ������ ����� = 3
// 12 - ���� arrKInp = 1, �-�� ������ ����� = 2
// 0  - ���� arrKInp = 0, �-�� ������ ����� = 1
// -1 - ������� arrKInp �� ������������ ����������
// OutPut:
// arrMatrProperVect - ������� ����������� ��������
// arrMatrProperNumb - ������������� ������� ��������������� ������ �����
// ����������� ����� ������������� � ������� �������� - �� ����
// arrMatrProperNumb[0] >= arrMatrProperNumb[4] >= arrMatrProperNumb[8]
int   CalcProperVectors_And_Numbers_R3(double * arrKInp,double *arrMatrProperVect , double *arrMatrProperNumb)
{
	 // 1. ������� ������������������� ��-� ��� �������  arrKInp

	 double arrLamb[3] ={0.} ; // ������ ������ ���������
	/* double a = - ( arrKInp [0] + arrKInp [4]  + arrKInp [8] ) ;
	 double b =  arrKInp [0] * arrKInp [4] +  arrKInp [0] * arrKInp [8] + arrKInp [4] * arrKInp [8]
			   - arrKInp [5] * arrKInp [5] - arrKInp [1] * arrKInp [1] - arrKInp [2] * arrKInp [2] ;
	 double c = -( arrKInp [0] * arrKInp [4]  * arrKInp [8] + 2 *  arrKInp [5] * arrKInp [1] * arrKInp [2]
			 - arrKInp [0] * arrKInp [5] * arrKInp [5] - arrKInp [8] * arrKInp [1] * arrKInp [1] - arrKInp [4] * arrKInp [2] * arrKInp [2]) ;
	 */
	  double a = - ( arrKInp [0] + arrKInp [4]  + arrKInp [8] ) ;
     double b =  arrKInp [0] * arrKInp [4] +  arrKInp [0] * arrKInp [8] + arrKInp [4] * arrKInp [8]
               - arrKInp [5] * arrKInp [7] - arrKInp [3] * arrKInp [1] - arrKInp [2] * arrKInp [6] ;
     double c = arrKInp [0] * arrKInp [5] * arrKInp [8]
             + arrKInp [4] * arrKInp [2] * arrKInp [6]
             +arrKInp [4] * arrKInp [3] * arrKInp [1]
              - arrKInp [0] * arrKInp [4]  * arrKInp [8]
              -arrKInp [2] * arrKInp [3]  * arrKInp [7]
             -arrKInp [3] * arrKInp [5]  * arrKInp [7];
	 int irez = Cubic(arrLamb, a, b, c) ;
	  ///


	 if (irez == 0) return -1 ; // ������, �������� arrKInp �� ������������-�����������

	 if( ( arrLamb[0] < -EPS) || ( arrLamb[1] < -EPS)|| ( arrLamb[2] < -EPS) ) return -1 ; // ������, �������� arrKInp �� ������������-�������������
	 int quantNull = 0 ; // �-�� ������� ������ �����
	  for (int i =0; i < 3; i++)
	  {
		  double dd = fabs(arrLamb[i]) ;
		  if( (dd - EPS) < 0.)
		  {
		   quantNull++;
          }

	  }
	 ///

	//qsort(arrLamb, 3, sizeof(double),&(fCompare)); // ���������� � ������� �����������
	 int iRang = 0 ;
	 int iReturn = -100;
	 for (int i =0; i < 3; i++)
	 {
	   if (arrLamb [i] > EPS) iRang++ ;
	   else arrLamb [i] = 0 ;

	 }

	 double mtrA[9] = {0}, det0,det1,det2, arrfTr[9] = {0} ;
	 int num1 = -1, num2 = -1 ;
	switch (irez)
	{
	case 3:
		for (int i = 0; i < 3; i++)
		{
			memcpy(mtrA, arrKInp, 9 * sizeof(double)) ;
			for (int j = 0; j < 3; j++) mtrA [ 3 * j + j ] -=  arrLamb[ i ];
			det0 =  mtrA[0] * mtrA[4 ] - mtrA[1] * mtrA[1 ] ;
			if (fabs(det0) < EPS)
			{
				det1 = mtrA[1] * mtrA[5 ] - mtrA[2] * mtrA[4 ] ;
				if (fabs(det1) < EPS)
				{
				det2               =  mtrA[0] * mtrA[8 ] - mtrA[2] * mtrA[2 ] ;
				arrfTr[3 * i  + 1] = 1 ;
				arrfTr[3 * i  ]    = ( - mtrA[1] * mtrA[8 ] + mtrA[2] * mtrA[5 ]) / det2;
				arrfTr[3 * i  + 2] = ( - mtrA[0] * mtrA[5 ] + mtrA[1] * mtrA[2 ]) / det2;
				}
				else
				{
				arrfTr[3 * i     ] = 1 ;
				arrfTr[3 * i + 1 ] = ( - mtrA[0] * mtrA[5 ] + mtrA[2] * mtrA[1 ]) / det1;
				arrfTr[3 * i + 2 ] = -det0 / det1;
				}

			}
			else
			{
			arrfTr[3 * i  + 2 ] = 1 ;
			arrfTr[3 * i   ]    = ( -mtrA[2] * mtrA[4 ] + mtrA[5] * mtrA[1 ]) / det0; ;
			arrfTr[3 * i + 1  ] = ( -mtrA[0] * mtrA[5 ] + mtrA[2] * mtrA[1 ]) / det0; ;

			}
		}
		if (iRang == 2)
		{
		  iReturn = 23 ;
		}
		else
		{
		  iReturn = 33 ;
		}

	break;
	case 2:
		memcpy(mtrA, arrKInp, 9 * sizeof(double)) ;
		for (int j = 0; j < 3; j++) mtrA [ 3 * j + j ] -=  arrLamb[ 0];
		det0 =  mtrA[0] * mtrA[4 ] - mtrA[1] * mtrA[1 ] ;
		//������ ������ ������� ��� ���������� ������ �����
		if (fabs(det0) < EPS)
		{
			det1 = mtrA[1] * mtrA[5 ] - mtrA[2] * mtrA[4 ] ;
			if (fabs(det1) < EPS)
			{
			det2               =  mtrA[0] * mtrA[8 ] - mtrA[2] * mtrA[2 ] ;
			arrfTr[ 1] = 1 ;
			arrfTr[0  ]    = ( - mtrA[1] * mtrA[8 ] + mtrA[2] * mtrA[5 ]) / det2;
			arrfTr[ 2] = ( - mtrA[0] * mtrA[5 ] + mtrA[1] * mtrA[2 ]) / det2;
			}
			else
			{
			arrfTr[0    ] = 1 ;
			arrfTr[ 1 ] = ( - mtrA[0] * mtrA[5 ] + mtrA[2] * mtrA[1 ]) / det1;
			arrfTr[ 2 ] = -det0 / det1;
			}

		}
		else
		{
		arrfTr[ 2 ] = 1 ;
		arrfTr[0   ]    = ( -mtrA[2] * mtrA[4 ] + mtrA[5] * mtrA[1 ]) / det0; ;
		arrfTr[ 1  ] = ( -mtrA[0] * mtrA[5 ] + mtrA[2] * mtrA[1 ]) / det0; ;

		}

		///

		// ���������� �������������� �������
		  // ����� ������� 2-� ��������� ����������� ������  arrfTr

		  for (int i = 0; i < 3; i++)
		  {
			if (fabs(arrfTr[i]) > EPS )
			{
			  if (num1 < 0)
			  {
				num1 = i ;
			  }
			  else
			  {
				  num2 = i ;
				  break ;
              }
			}
		  }
		///

		if (num2 == -1)
		{ // � ������� arrfTr ������� ���� ��������� ����������
		  num2 = (num1 + 1)% 3 ;
		  num1 = (num1 + 2)% 3 ;
		  arrfTr[3 + num1] = 1 ;
		  arrfTr[6 + num2] = 1 ;

		}
		else
		{
		  arrfTr[3 + num1] =  arrfTr[ num1] ;
		  arrfTr[3 + num2] = -arrfTr[num2]  ;
			OuterProduct(arrfTr , &arrfTr [3] , &arrfTr [6]) ;

        }

		 switch(iRang )
		 {
			 case 1:
			 iReturn = 12 ;
			 break;
			 case 2:
			 iReturn = 22 ;
			 break;
			 case 3:
			 iReturn = 32 ;
			 break;
			 case 0:
			 default:
			 return -5 ;

		 }
    break ;
	case 1:
		switch (quantNull)
		{
		case 0:
		  arrfTr[0    ] = 1 ;
		  arrfTr[3 + 1] = 1 ;
		  arrfTr[6 + 2] = 1 ;
		break ;
		case 3:
	 //	ShowMessage(L"Error -3") ;
		return -3;

		case 2:
		case 1:
		default:
	   //	ShowMessage(L"Error -1") ;
		return -1;
		;
		}
		switch(iRang )
		 {
			 case 3:
			 iReturn = 31 ;
			 break;
			 case 0:
			 iReturn = 0 ;
			 break;
			 case 1:
			 case 2:
			 default:
			 return -6 ;

		 }

	break ;
	default:
  //	ShowMessage(L"Error -2") ;
	 return -2;
	;
	 }

	// ���������� �� ��������
	bool bPr = true;
	while(bPr)
	{
		bPr = false ;
		for (int i = 0; i < 2; i++)
		{
		  if (arrLamb [ i + 1] > arrLamb [ i ])
		  {
			swap(&arrLamb [ i + 1], & arrLamb [ i ] ) ;
			swap_vect(&arrfTr [ 3 * (i + 1)], & arrfTr [3 *  i ], 3 ) ;
			bPr = true ;
			break ;
		  }
		}
	}
///
	 // �������������, ���������������
	  memset(arrMatrProperNumb , 0 , 9 * sizeof(double)) ;
	 for (int i = 0; i < 3; i++)
	 {
	  NormalizeVect3( &arrfTr [ i * 3]) ;
	  arrMatrProperNumb [ 3 * i + i ] =  arrLamb [ i ];
	 }
	 MatrTransp(arrfTr, 3, 3, arrMatrProperVect ) ;


	return iReturn ;
}

void   swap(double *a0, double *a1)
{
	double b = *a1;
	*a1 = *a0 ;
	*a0 = b ;
}

void   swap_vect(double *pa0, double *pa1, const int len)
{
	for (int i = 0; i < len; i++) swap( &pa0[i],&pa1[i]);


}

double   NormVect3(double *p)
{
	return sqrt(p[0] * p [0] + p[1] * p [1] + p[2] * p [2] ) ;
}

double   NormVect2(double *p)
{
	return sqrt(p[0] * p [0] + p[1] * p [1]  ) ;
}

void   NormalizeVect3(double *p)
{
	double valNorm =  NormVect3(p) ;
	for (int i = 0; i < 3; i++)
	{
	  p[i] = p[i] / valNorm ;
	}
}

// ����� ����������� ��������  � ����������� ����� ������������ ������������ ������� arrKInp
//  arrV - ������� ����������� ��������
//  arrLamb - ��������������� ������� ����������� �����
int  CalcProperVectors2(double * arrKInp,double *arrV , double *arrLamb)
{

	memset(arrLamb,0,4 * sizeof(double)) ;
	double a  = 1;
	double b = -(arrKInp[0] + arrKInp [3]);
	double c = arrKInp[0] * arrKInp [3] - arrKInp[1] * arrKInp [1] ; ;
	TComp x0,x1;
	int i0 = SolvEq2( a, b, c,x0,x1) ;
	double arrx[2];
	arrx[0] = x0.m_Re;
	arrx[1] = x1.m_Re;

	// ����������� ������
	try
	{
	for (int i =0; i < 2; i++)
	{
	  if (fabs (arrx[i]) < 0.0001 )
	  {
		int j = (i + 1)%2 ;
		double valLambPos = arrx[j] ;
		if (valLambPos < 0) return -1 ;


		if ((fabs( arrKInp[0] - valLambPos) > 0.000001) || (fabs( arrKInp[1] ) > 0.000001))
		{
		  double arre[2] = {0.} ;
		  double d = sqrt( arrKInp[1]* arrKInp[1] + (arrKInp[0] - valLambPos) * (arrKInp[0] - valLambPos)) ;
		  arre[0] =  -arrKInp[1] /d;
		  arre[1] =  (arrKInp[0] - valLambPos)/d  ;
		  arrV [0] = arre[0] ;
		  arrV [2] = arre[1] ;
		  arrV [1] = -arre[1] ;
		  arrV [3] =  arre[0] ;
		  arrLamb [0] = valLambPos ;

		  return 0 ;
		}
		else
		{
		   return -1 ;

		}
	  }
	}
	}
	catch(...)
	{
   //		ShowMessage( L"TFragmTraj::CalcProperVectors2 Error No 0");
		return -1;
	}
  ///

   try
   {
	switch(i0)
	{
		case 0:

		for (int i = 0; i < 2 ; i++)
		{
		  if  ( arrx[i] < 0)
		  {
			  int j = (i + 1)%2 ;
			  double arrProp[2] ={0},arrPropTemp[4] = {0} ;
			  if ( fabs (arrKInp[0] - arrx[j])> EPS)
			  {
				arrProp[0] = - arrKInp[1]/ (arrKInp[0] - arrx[j]);
				arrProp[1] = 1 ;
			  }
			  else
			  {
				arrProp[0] = 1 ;
				arrProp[1] =  -  (arrKInp[0] - arrx[j])/arrKInp[1];
              }
			  double d = sqrt(arrProp[0] * arrProp[0] + arrProp[1] * arrProp[1]) ;
			  for (int ii = 0; ii < 2; ii++) arrProp[ii] = arrProp[ii] / d ;


			  MtrxMultMatrxTransp(arrProp,2, 2, arrProp,2, arrPropTemp) ;
			  MatrxMultScalar(arrPropTemp, 2 , 2, arrx[j],arrKInp);
			  return 1 ;
		  }
		}
		for (int i = 0; i < 2; i++)
		{
		 arrLamb[i * 2 + i] = arrx[i];
		 if (fabs(arrKInp[0] - arrx[i] ) > 0.00000000001 )
		 {
		   arrV [ i ] = - arrKInp[1] / (arrKInp[0] - arrx[i] );
		   arrV [ 2 + i ] = 1 ;
		   double r = sqrt (arrV [ i ] *  arrV [ i ] + arrV [ 2 + i ] * arrV [ 2 + i ] );
		   arrV [ i ] = arrV [ i ] /r ;
		   arrV [ 2 + i ] = arrV [ 2 + i ] /r ;

		 }
		 else
		 {
		   arrV [ i ]     = 1 ;
		   arrV [ 2 + i ] = 0 ;
         }
		}
		break ;
		case 1:
		arrV [0] = 1 ;
		arrV [1] = 0 ;
		arrV [2] = 0 ;
		arrV [3] = 1 ;
		arrLamb [0] = arrKInp[0] ;
		arrLamb [1] = 0 ;
		arrLamb [2] = 0 ;
		arrLamb [3] = arrKInp[0] ;
		break ;

		default:

		return 1;

	}
	}
	catch(...)
	{
   //		ShowMessage( L"TFragmTraj::CalcProperVectors2 Error No 2");
		return -1;

    }
	return 0 ;
}
// ������� ������� �������� ��������� ������� �������  A * x = B
// Returns:
// true , ���� det(A) != 0
// false - � ��������� ������
bool   SolvLinEq2(double *arrA, double *arrB,double *arrX)
{
   double det = arrA [ 0] * arrA [ 3] - arrA [ 1] * arrA [ 2] ;
   if (fabs(det) < EPS) return false ;
   arrX [ 0] =   ( arrB [ 0] *  arrA [ 3] - arrA [ 1] * arrB [ 1] ) / det ;
   arrX [ 1] = ( arrA [ 0] * arrB [ 1]  - arrA [ 2] * arrB [ 0] ) / det ;

  return true;
}

// ��������� ������� 2-�� �������
// Returns:
// true , ���� det(A) != 0
// false - � ��������� ������
bool   InverseMtrx2(double *arrA, double *arrOut)
{
   double det = arrA [ 0] * arrA [ 3] - arrA [ 1] * arrA [ 2] ;
   if (fabs(det) < EPS) return false ;
   arrOut [ 0] =  arrA [ 3] / det ;
   arrOut [ 1] = -arrA [ 1] / det ;
   arrOut [ 2] = -arrA [ 2] / det ;
   arrOut [ 3] =  arrA [ 0] / det ;
  return true;
}

double NormVect(double *arr, const int lenarr)
{
  double sum = 0;
  for (int i = 0; i < lenarr; i++)
  {
	sum += arr[i] * arr[i] ;
  }
  return sqrt(sum) ;
}



// ������������ ��������������� ������ � �������� � ������  iRowElem, iColElem
// parrMtrxInp - �������� �������
// ndim - �����������  �������� �������
// iRowElem , jColElem  - ����� ��������
// OUTPUT:
// parrMinor - �������������� ����� �����������  [ndim -1]x [ndim -1]
void buildAddMinor(double *parrMtrxInp, const int ndim
	, const int iRowElem,const int jColElem, double *parrMinor)
{
	int iMinor = 0 ;
	int jMinor = 0 ;
	for (int i = 0; i < ndim; i++)
	{
	 if (i == iRowElem) continue;

	 for (int j = 0; j < ndim; j++)
	 {
	   if (j == jColElem) continue ;

		parrMinor [ (ndim-1) * iMinor + jMinor] =  parrMtrxInp [ ndim * i + j] ;
	   jMinor++;
	 }
	 iMinor++;
	 jMinor = 0;
	}
}

// ������������� �������� 4�4
double   calcDet4( double *arrInp)
{
	double det = 0.;
	double parrMinor [9] = {0.} ;
	double valIndex = 1.;
	for (int j = 0; j < 4; j++)
	{
	  buildAddMinor(arrInp, 4, 0, j, parrMinor) ;
	  det += valIndex * arrInp[j] * calcDet3( parrMinor) ;
	  valIndex = -valIndex;
	}
	return det;
}

//  ��������� ������ 4�4
//
bool   InverseMtrx4(double *arrInp, double *arrOut)
{
	double det =   calcDet4( arrInp);
	if (fabs(det) < 1.E-15) return false ;
	double parrMinor [9] = {0.} ;
	double valIndexI = 1.;
	double valIndexJ = 1.;
	for (int i = 0; i < 4; i++)
	{
	for (int j = 0; j < 4; j++)
	{
	 buildAddMinor(arrInp, 4, j,i, parrMinor) ;
	 arrOut[ 4 * i + j] =  valIndexI * valIndexJ * calcDet3( parrMinor)/ det;
	 valIndexJ = -valIndexJ;
	}
	valIndexI = -valIndexI;
	valIndexJ = 1.;
	}
	return true;
}

// ������������� �������� 5x5
double   calcDet5( double *arrInp)
{
	double det = 0.;
	double parrMinor [16] = {0.} ;
	double valIndex = 1.;
	for (int j = 0; j < 5; j++)
	{
	  buildAddMinor(arrInp, 5, 0, j, parrMinor) ;
	  det += valIndex * arrInp[j] * calcDet4( parrMinor) ;
	  valIndex = -valIndex;
	}
	return det;
}
//  ��������� ������ 5�5
//
bool   InverseMtrx5(double *arrInp, double *arrOut)
{
	double det =   calcDet5( arrInp);
	if (fabs(det) < 1.E-15) return false ;
	double parrMinor [16] = {0.} ;
	double valIndexI = 1.;
	double valIndexJ = 1.;
	for (int i = 0; i < 5; i++)
	{
	for (int j = 0; j < 5; j++)
	{
	 buildAddMinor(arrInp, 5, j,i, parrMinor) ;
	 arrOut[ 5 * i + j] =  valIndexI * valIndexJ * calcDet4( parrMinor)/ det;
	 valIndexJ = -valIndexJ;
	}
	valIndexI = -valIndexI;
	valIndexJ = 1.;
	}
	return true;
}

// ������������� �������� 6x6
double   calcDet6( double *arrInp)
{
	double det = 0.;
	double parrMinor [25] = {0.} ;
	double valIndex = 1.;
	for (int j = 0; j < 6; j++)
	{
	  buildAddMinor(arrInp, 6, 0, j, parrMinor) ;
	  det += valIndex * arrInp[j] * calcDet5( parrMinor) ;
	  valIndex = -valIndex;
	}
	return det;

}

// ������������� �������� 7x7
double   calcDet7( double *arrInp)
{
	double det = 0.;
	double parrMinor [36] = {0.} ;
	double valIndex = 1.;
	for (int j = 0; j < 7; j++)
	{
	  buildAddMinor(arrInp, 7, 0, j, parrMinor) ;
	  det += valIndex * arrInp[j] * calcDet6( parrMinor) ;
	  valIndex = -valIndex;
	}
	return det;

}

// ������������� �������� 8x8
double   calcDet8( double *arrInp)
{
	double det = 0.;
	double parrMinor [49] = {0.} ;
	double valIndex = 1.;
	for (int j = 0; j < 8; j++)
	{
	  buildAddMinor(arrInp, 8, 0, j, parrMinor) ;
	  det += valIndex * arrInp[j] * calcDet7( parrMinor) ;
	  valIndex = -valIndex;
	}
	return det;

}

//  ��������� ������ 6x6
//
bool   InverseMtrx6(double *arrInp, double *arrOut)
{
	double det =   calcDet6( arrInp);
	if (fabs(det) < 1.E-15) return false ;
	double parrMinor [25] = {0.} ;
	double valIndexI = 1.;
	double valIndexJ = 1.;
	for (int i = 0; i < 6; i++)
	{
	for (int j = 0; j < 6; j++)
	{
	 buildAddMinor(arrInp, 6, j,i, parrMinor) ;
	 arrOut[ 6 * i + j] =  valIndexI * valIndexJ * calcDet5( parrMinor)/ det;
	 valIndexJ = -valIndexJ;
	}
	valIndexI = -valIndexI;
	valIndexJ = 1.;
	}
	return true;
}

//------------------------------------------------------------------------
// ���������� ������������ �������� ������� �������
// parrDoubleMtrx - ������ ��������� �������� �������
// nRows , nCols - �-�� ����� � �������� �������
// NumCol - ����� ������� � �������� ������� �������
// OUTPUT:
// *pNumArgMin - ����� ������ � ������������ ���������
// ���������� �������� ������������� ��������
//
//
 double MinColDoubleMtrx(double *parrDoubleMtrx, const int nRows, const int nCols
	 , const int NumCol, int *pNumArgMin)
{
   *pNumArgMin = -1;
	double valRez = DBL_MAX ;
	for (int i = 0; i < nRows; i++)
	{
	if (parrDoubleMtrx [i * nCols + NumCol] < valRez)
	{
	 *pNumArgMin = i ;
	 valRez = parrDoubleMtrx [i * nCols + NumCol] ;
	}
	}
	return valRez ;
}
//---------------------------------------------------------------------
// ��������� ��������� ����������� ������� �� ��������
void flipArray(double *parr, const int len)
{
	for (int i = 0; i < len/ 2; i++)
	{
		 swap(&parr[i], &parr[len - 1 -i]) ;

	}
}
// ��������� ��������� ���������� ������� �� �������� (�� �������)
void flipTwoDimArray(double *parr, const int lenrows, const int lencols)
{
	for (int i = 0; i < lenrows/ 2; i++)
	{
		 swap_vect(&parr[ i *lencols] , &parr[(lenrows - 1 - i) *lencols] , lencols) ;

	}
}



////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

//  ��������� ������ 7x7
//
bool   InverseMtrx7(long double  *arrInp, long double  *arrOut)
{
	long double  det =   calcDet7( arrInp);
	if (fabs(det) < 1.E-15) return false ;
	long double  parrMinor [36] = {0.} ;
	long double  valIndexI = 1.;
	long double  valIndexJ = 1.;
	for (int i = 0; i < 7; i++)
	{
	for (int j = 0; j < 7; j++)
	{
	 buildAddMinor(arrInp, 7, j,i, parrMinor) ;
	 arrOut[ 7 * i + j] =  valIndexI * valIndexJ * calcDet6( parrMinor)/ det;
	 valIndexJ = -valIndexJ;
	}
	valIndexI = -valIndexI;
	valIndexJ = 1.;
	}
	return true;
}

//  ��������� ������ 8x8
//
bool   InverseMtrx8(long double  *arrInp, long double  *arrOut)
{
	long double  det =   calcDet8( arrInp);
	if (fabs(det) < 1.E-15) return false ;
	long double  parrMinor [49] = {0.} ;
	long double  valIndexI = 1.;
	long double  valIndexJ = 1.;
	for (int i = 0; i < 8; i++)
	{
	for (int j = 0; j < 8; j++)
	{
	 buildAddMinor(arrInp, 8, j,i, parrMinor) ;
	 arrOut[ 8 * i + j] =  valIndexI * valIndexJ * calcDet7( parrMinor)/ det;
	 valIndexJ = -valIndexJ;
	}
	valIndexI = -valIndexI;
	valIndexJ = 1.;
	}
	return true;
}


// ���������� ��������� ������� ������������ �����������
//
bool   InverseMtrx( long double  *arrInp, const int nDimen,  long double  *arrOut)
{
  switch (nDimen)
  {
  case 1:
  if (fabs(*arrInp) > 0.000000001)
  {
	*arrOut = 1./ (*arrInp);
	return true;
  }
  else
  {
	  return false;
  }
  case 2:
  return InverseMtrx2(arrInp, arrOut);

  case 3:
  return InverseMtrx3(arrInp, arrOut);


  case 4:
  return InverseMtrx4(arrInp, arrOut);


  case 5:
  return InverseMtrx5(arrInp, arrOut);


  case 6:
  return InverseMtrx6(arrInp, arrOut);

  case 7:
  return InverseMtrx7(arrInp, arrOut);

  case 8:
  return InverseMtrx8(arrInp, arrOut);

  default: return false ;
      ;
  }
 // return false ;
}


bool GaussMeth(double *parrA,int nDimA,  double * parrB, double *parrRez)
 {
	 double *arrAUn  = new double [ nDimA * (nDimA +1)];
	 int *iarrNums  = new int [ nDimA ];
	 for (int i =0; i < nDimA; i++)
	 {
	   memcpy( &arrAUn[ i *  (nDimA +1)], & parrA[i * nDimA ], sizeof(double) * nDimA);
	   arrAUn[ i *  (nDimA +1) + nDimA] =  parrB[i];
	   iarrNums [i] = i;
	 }
	 for (int i = 0; i < nDimA; i++)
	 {
		int numVedStroka = fncFindNumVedushiaStroka(arrAUn,  nDimA, i);
		if (numVedStroka == -1)
		{
			delete arrAUn;
			delete iarrNums;
			return false;
		}
		Swap (iarrNums,1, numVedStroka, i);
		Swap (arrAUn,nDimA +1, numVedStroka, i);
		for (int j = i + 1; j < nDimA; j++)
		{
		  if (fabs(arrAUn [j *(nDimA + 1) +i]) < 10. * DBL_MIN)
		  {
			continue;
		  }
		  if (!fncPodstanovka( arrAUn, nDimA,i, j) )
		  {
			delete arrAUn;
			delete iarrNums;
			return false;
          }
		}
	 }

	 double *arrRezTemp = new double [nDimA];
	 memset(arrRezTemp, 0, nDimA * sizeof(double));
	 arrRezTemp[nDimA -1]= arrAUn[ (nDimA-1) *  (nDimA +1) + nDimA]/  arrAUn[ (nDimA-1) *  (nDimA +1) + nDimA -1];
	 for (int i = nDimA -2; i >= 0; i--)
	 {
	   double sum = arrAUn[ i *  (nDimA +1) + nDimA];
	   for (int j =i+1; j < (nDimA -1); j++)
	   {
		 sum -=  arrRezTemp[j] * arrAUn[ i *  (nDimA +1) + nDimA -j];
	   }
	   arrRezTemp[i] = sum/  arrAUn[ i *  (nDimA +1) + i];
	 }
	 for (int i =0; i < nDimA; i++)
	 {
	   parrRez [i] =  arrRezTemp[iarrNums[i]];
	 }
	 delete arrAUn;
	 delete iarrNums;
	 delete arrRezTemp;
	 return true;
 }

bool fncPodstanovka( double *arrAUn, int nDimA, int i, int j)
{

	double temp = -  arrAUn [i *(nDimA + 1) +i] /arrAUn [j *(nDimA + 1) +i];
	arrAUn [j *(nDimA + 1) +i] = 0.;
	for (int k =i + 1; k < (nDimA + 1); k++)
	{
	  arrAUn [j *(nDimA + 1) +k] = arrAUn [i *(nDimA + 1) +k] + temp * arrAUn [j *(nDimA + 1) +k] ;
	}
	return true;
}

int fncFindNumVedushiaStroka(double *arrAUn, int  nDimA, int i0)
{
	int irez = -1;
	for (int i = i0; i < nDimA; i++)
	{
	  if (fabs(arrAUn [i *( nDimA +1) + i0]) > 10. * DBL_MIN )
	  {
		return i;
	  }
	}
	return irez ;
}

void Swap (int *iarrNums, int iNumCols, int j,int  i)
{
	int *iarrTemp  = new int [iNumCols];
	memcpy(iarrTemp ,&iarrNums[j * iNumCols], iNumCols * sizeof(int));
	memcpy(&iarrNums[j * iNumCols] ,&iarrNums[i * iNumCols], iNumCols * sizeof(int));
	memcpy(&iarrNums[i * iNumCols] ,iarrTemp, iNumCols * sizeof(int));

	delete iarrTemp ;
}

void Swap (double *parrA, int iNumCols, int j,int  i)
{
	double *parrTemp  = new double [iNumCols];
	memcpy(parrTemp ,&parrA[j * iNumCols], iNumCols * sizeof(double));
	memcpy(&parrA[j * iNumCols] ,&parrA[i * iNumCols], iNumCols * sizeof(double));
	memcpy(&parrA[i * iNumCols] ,parrTemp, iNumCols * sizeof(double));

	delete parrTemp ;
}

//-------------------------------------------------------------------------------------
// ������������ ���������� ������� �����������  IdimE
void formMatrxE(const int IdimE, double *arrOut)
{
 memset(arrOut, 0, IdimE * IdimE * sizeof(double));
 for (int i = 0; i < IdimE; i++)
 {
	 arrOut [ i * IdimE + i] = 1.;
 }
}

//---------------------------------------------------------------------------------------------
// �������� �������������� ������
// ������ �������   arrVectMissV[3] �  arrTargV[3]
// 1.������� �����������
// 2. ������������ ��������������� ������� arrTargV ������������  arrVectMissV,
// �� ����  arrTargV1 = arrTargV - <arrTargV , arrVectMissV> arrVectMissV
// 3. ������ ������ ������ ��� ��������� ������������
//
//
void createOrthogBasis_dim3 (double *arrVectMissV, double *arrTargV, double *arrBasis)
{
	memset (arrBasis, 0, 9 * sizeof(double));
	double arrx1[3] = {0.}, arrx21[3] = {0.}, arrx2[3] = {0.}, arrx3[3] = {0.}, arrt[3] ={0.};
	memcpy(arrx1, arrVectMissV, 3 * sizeof(double));
	memcpy(arrx2, arrTargV, 3 * sizeof(double));
	NormalizeVect3(arrx1) ;
	NormalizeVect3(arrx2) ;
	double valScal = ScalProduct(arrx1 , arrx2, 3) ;
	if (fabs(valScal) > 0.9999999999)  // arrVectMissV �  arrTargV �����������
	{
		if (fabs (arrx1[2]) > 0.999999999) // arrVectMissV �  arrTargV ����������� � �����������
		{
			arrBasis[1] = 1.;
			arrBasis[5] = 1.;
			arrBasis[6] = 1.;
			return;
		}
		arrx21[0] =  -arrx1[1];
		arrx21[1] =   arrx1[0];
	}
	else
	{
	MatrxMultScalar(arrx1, 1, 3,  valScal,arrt);
	MtrxMinusMatrx(arrx2, arrt,1, 3, arrx21);
	}
	NormalizeVect3(arrx21) ;
	OuterProduct(arrx21 , arrx1, arrx3)  ;

	double arrTemp[9] = {0.};
	memcpy(arrTemp, arrx1, 3 * sizeof(double));
	memcpy(&arrTemp[3], arrx21, 3 * sizeof(double));
	memcpy(&arrTemp[6], arrx3, 3 * sizeof(double));
	MatrTransp(arrTemp, 3, 3, arrBasis);
}

//---------------------------------------------------------------------------------------------
// ���������� ���� ����� ���������
// ������ �������   arr1 �  arr2

double  calcAngBetweenVect (double *arr1, double *arr2,const  int len)
{
	return acos(ScalProduct(arr1, arr2,  len)/ NormVect(arr1,  len)/ NormVect(arr2,  len) );
}



void   NormalizeVect(double *p, const int lenp)
{
double norm = NormVect(p, lenp);
for (int i =0; i < lenp; i++)
{
  p[i] =  p[i]/norm ;
}

}


double _sqrt_(const double a)
{
  return (a <0.)?0.:sqrt(a);
}
 //**********************************************************************************************************************
 //**********************************************************************************************************************
 //************   ��� Long double       **********************************************************************************************************
 //**********************************************************************************************************************
 //**********************************************************************************************************************
 //**********************************************************************************************************************
//
//---------------------------------------------------------------------------
// ��������� ������� �� �������
void MtrxMultMatrx(long double *parrA,int nRowsA, int nColsA, long double * parrB,int nColsB, long double *parrRez)
{
  memset (parrRez,0,sizeof(long double) * nRowsA *  nColsB) ;
 long double *pRez = parrRez ;
  long double *pA = parrA ;
  long double *pB = parrB ;

  for (int i = 0; i < nRowsA; i++)
  {
	long double *pAt = pA;
   long double *pBt =pB ;
  for (int j = 0; j < nColsB; j++)
  {
	pAt = pA ;
	pBt = parrB + j ;

	for (int k = 0; k < nColsA ; k++)
	{
	  *pRez += (*pAt) * (*pBt);
	  pAt ++;
	  pBt += nColsB ;
	}
	pRez++;
  }
  pA += nColsA ;
  }
/*  for (int i = 0; i < nRowsA; i++)
  {
	 for (int j = 0; j < nColsB; j++)
	 {
		int ktemp = i* nColsB + j ;
		for (int k = 0; k < nColsA ; k++)
		{
		  parrRez [ktemp] += parrA[i * nColsA + k]*  parrB [ k *  nColsB + j ] ;
		}
	 }
  }*/
}

// ��������� ������� �� ����������������� �������
void MtrxMultMatrxTransp(long double *parrA,int nRowsA, int nColsA, long double * parrB,int nRowsB, long double *parrRez)
{
/*  memset (parrRez,0,sizeof(long double) * nRowsA *  nRowsB) ;
 for (int i = 0; i < nRowsA; i++)
  {
	 for (int j = 0; j < nRowsB; j++)
	 {
		int ktemp = i* nRowsB + j ;
		for (int k = 0; k < nColsA ; k++)
		{
		  parrRez [ktemp] += parrA[i * nColsA + k]*  parrB [ j *  nColsA + k ] ;
		}
	 }
  } */

   memset (parrRez,0,sizeof(long double) * nRowsA *  nRowsB) ;
  long double *pRez = parrRez ;
  long double *pA = parrA ;
  long double *pB = parrB ;

  for (int i = 0; i < nRowsA; i++)
  {
	long double *pAt = pA;
   long double *pBt =pB ;
  for (int j = 0; j < nRowsB; j++)
  {
	pAt = pA ;


	for (int k = 0; k < nColsA ; k++)
	{
	  *pRez += (*pAt) * (*pBt);
	  pAt ++;
	  pBt ++;
	}
	pRez++;
  }
  pA += nColsA ;
  }
}
// ��������� ����������������� ������� ��  �������
void MtrxTranspMultMatrx(long double *parrA,int nRowsA, int nColsA, long double * parrB,int nColsB, long double *parrRez)
{
  memset (parrRez,0,sizeof(long double) * nColsA *  nColsB) ;
 for (int i = 0; i < nColsA ; i++)
  {
	 for (int j = 0; j < nColsB; j++)
	 {
		int ktemp = i* nColsB + j ;
		for (int k = 0; k < nRowsA ; k++)
		{
		  parrRez [ktemp] += parrA[k * nColsA + i]*  parrB [ k *  nColsB + j ] ;
		}
	 }
  }
  /*    memset (parrRez,0,sizeof(long double) * nColsA *  nColsB) ;

  long double *pRez = parrRez ;
  long double *pA = parrA ;
  long double *pB = parrB ;

  for (int i = 0; i < nColsA; i++)
  {
	//long double *pAt = parrA + i;

  for (int j = 0; j < nColsB; j++)
  {
	long double *pBt =parrB + j ;
	long double *pAt = parrA + i;

	for (int k = 0; k < nRowsA ; k++)
	{
	  *pRez += (*pAt) * (*pBt);
	  pAt += nColsA ;
	  pBt += nColsB ;
	}
	pRez++;
  }

  }*/
}
void MtrxSumMatrx(long double *pA, long double * pB,int nRows, int nCols, long double *pRez)
{
	long double *parrA = pA;
	long double * parrB = pB ;
	long double *parrRez = pRez;
	for (int i = 0; i < nRows * nCols; i++)
	{
	  *parrRez = (*parrA)+ (* parrB);
	  parrRez++;
	  parrA++ ;
	  parrB++ ;

	}
}
void MtrxMinusMatrx(long double *pA, long double * pB,int nRows, int nCols, long double *pRez)
{
	long double *parrA = pA;
	long double * parrB = pB ;
	long double *parrRez = pRez;

	for (int i = 0; i < nRows * nCols; i++)
	{
	  *parrRez = (*parrA)- (* parrB);
	  parrRez++;
	  parrA++ ;
	  parrB++ ;

	}
}

void MatrxMultScalar(long double *parrA, int nRows, int nCols, long double valScal,long double *parrRez)
{
   for (int i = 0; i < nRows * nCols; i++)
	{
	  *parrRez = (*parrA)* valScal ;
	  parrRez++;
	  parrA++ ;


	}
}
void MatrxDivideScalar(long double *parrA, int nRows, int nCols, long double valScal,long double *parrRez)
{
   for (int i = 0; i < nRows * nCols; i++)
	{
	  *parrRez = (*parrA)/ valScal ;
	  parrRez++;
	  parrA++ ;


	}
}
void MatrTransp(long double *parrA, int nRows, int nCols, long double *parrRez)
{
	long double *pA = parrA ;
	for (int i = 0; i < nRows; i++)
	{

	long double *pRezt =parrRez + i ;
		for (int j = 0; j < nCols; j++)
		{
		*pRezt = *pA;
		pA ++;
		pRezt += nRows ;
		}

	}
}

// ��������� ������������
void OuterProduct(long double *pVect0 , long double *pVect1, long double *pVectRez)
{
 pVectRez [0] = pVect0 [1] * pVect1 [2] -  pVect0 [2] * pVect1 [1] ;
 pVectRez [1] = pVect0 [2] * pVect1 [0] -  pVect0 [0] * pVect1 [2] ;
 pVectRez [2] = pVect0 [0] * pVect1 [1] -  pVect0 [1] * pVect1 [0] ;
}

// ���������  ������������
long double ScalProduct(long double *pVect0 , long double *pVect1, const int len)
{
 long double prod = 0. ;
 for (int i = 0; i < len; i++)
 {
  prod += pVect0[i] * pVect1[i] ;
 }

 return prod ;
}

long double Norm3( long double *arrA)
{
	return sqrt( arrA[0] * arrA[0] + arrA[1] * arrA[1] + arrA[2] * arrA[2]) ;
}

// ��������� ������������ 2-� ������ �������� - ������ �������
// >0 ���� ���������� ������� �� ������� ������� �� ������� �����
// ������ ������� �������  (���� ������� ������)
long double OuterProduct_2(long double *pVect0 , long double *pVect1)
{
 return pVect0 [0] * pVect1[1] - pVect0 [1] * pVect1[0];
}

// ������������ �������� ���������� ������� ���� long double
// INPUT:
// parrlong doubleInp - ������ ���� long double
// lenArrayInp - ����� ������� lenArrayInp
// OUTPUT:
// pNumArgMin - ����� ��������� �������, �� ������� ����������� �������
// ������� ���������� ����������� ��������
long double MinDoubleArray(long double *parrDblInp, const int lenArrayInp, int *pNumArgMin)
{
  *pNumArgMin = -1;
  long double valRez = DBL_MAX ;
  for (int i = 0; i < lenArrayInp; i++)
  {
	if (parrDblInp [i] < valRez)
	{
	 *pNumArgMin = i ;
	 valRez = parrDblInp [i] ;
	}
  }
  return valRez ;
}

// ������������ �������� ���������� ������� ���� long double
// INPUT:
// parrlong doubleInp - ������ ���� long double
// lenArrayInp - ����� ������� lenArrayInp
// OUTPUT:
// pNumArgax - ����� ��������� �������, �� ������� ����������� �������
// ������� ���������� ������������ ��������
long double MaxDoubleArray(long double *parrDblInp, const int lenArrayInp, int *pNumArgMax)
{
  *pNumArgMax = -1;
  long double valRez = DBL_MIN ;
  for (int i = 0; i < lenArrayInp; i++)
  {
	if (parrDblInp [i] > valRez)
	{
	 *pNumArgMax = i ;
	 valRez =parrDblInp [i] ;
	}
  }
  return valRez ;
}






// ���������� �������� �������  3x3
// ���� ������������  arrInp != ����, �� True
// � ��������� ������ false
 bool   InverseMtrx3(long double *arrInp, long double *arrOut)
{
  // ���������� ������������
  long double det = arrInp[0] * ( arrInp[4] * arrInp[8] - arrInp[5] * arrInp[7] )
			-  arrInp[1] * ( arrInp[3] * arrInp[8] - arrInp[5] * arrInp[6] )
			+  arrInp[2] * ( arrInp[3] * arrInp[7] - arrInp[4] * arrInp[6] )  ;
  ///
  if (fabsl (det) < 1E-10) return false ;

  long double arr0[9]= {0};

  arr0[0]=  ( arrInp[4] * arrInp[8] - arrInp[5] * arrInp[7] ) ;
  arr0[3]= -( arrInp[3] * arrInp[8] - arrInp[5] * arrInp[6] ) ;
  arr0[6]=  ( arrInp[3] * arrInp[7] - arrInp[4] * arrInp[6] ) ;
  arr0[1]= -( arrInp[1] * arrInp[8] - arrInp[2] * arrInp[7] ) ;
  arr0[4]=  ( arrInp[0] * arrInp[8] - arrInp[2] * arrInp[6] ) ;
  arr0[7]= -( arrInp[0] * arrInp[7] - arrInp[1] * arrInp[6] ) ;
  arr0[2]=  ( arrInp[1] * arrInp[5] - arrInp[2] * arrInp[4] ) ;
  arr0[5]= -( arrInp[0] * arrInp[5] - arrInp[2] * arrInp[3] ) ;
  arr0[8]=  ( arrInp[0] * arrInp[4] - arrInp[1] * arrInp[3] ) ;

  MatrxMultScalar(arr0, 3, 3, 1./ det, arrOut);



  return true ;
}
 // ���������� ������������
long double   calcDet3( long double *arrInp)
{
 return arrInp[0] * ( arrInp[4] * arrInp[8] - arrInp[5] * arrInp[7] )
			-  arrInp[1] * ( arrInp[3] * arrInp[8] - arrInp[5] * arrInp[6] )
			+  arrInp[2] * ( arrInp[3] * arrInp[7] - arrInp[4] * arrInp[6] )  ;
}









// ���������� YT * D * Y ��� ������� 3�3
long double  calcYT_D_Y(long double *arrY, long double *arrD )
{
  long double arrT0[3] = {0.} ;
  long double valReturn = 0 ;
  MtrxMultMatrx(arrY,1, 3, arrD,3, arrT0) ;
  MtrxMultMatrx(arrT0,1, 3, arrY,1, &valReturn) ;
  return valReturn ;
}

 // ����� ����������� ��������  � ����������� ����� ������������ ������������ ������� arrKInp 3� 3
//  arrV - ������� ����������� ��������
//  arrLamb - ��������������� ������� ����������� �����  3 � 3
//  Returns:
// 31 - ���� arrKInp = 3, �-�� ������ ����� = 1
// 32 - ���� arrKInp = 3, �-�� ������ ����� = 2
// 33 - ���� arrKInp = 3, �-�� ������ ����� = 3
// 22 - ���� arrKInp = 2, �-�� ������ ����� = 2
// 23 - ���� arrKInp = 2, �-�� ������ ����� = 3
// 12 - ���� arrKInp = 1, �-�� ������ ����� = 2
// 0  - ���� arrKInp = 0, �-�� ������ ����� = 1
// -1 - ������� arrKInp �� ������������ ����������
// OutPut:
// arrMatrProperVect - ������� ����������� ��������
// arrMatrProperNumb - ������������� ������� ��������������� ������ �����
// ����������� ����� ������������� � ������� �������� - �� ����
// arrMatrProperNumb[0] >= arrMatrProperNumb[4] >= arrMatrProperNumb[8]
int   CalcProperVectors_And_Numbers_R3(long double * arrKInp,long double *arrMatrProperVect , long double *arrMatrProperNumb)
{
	 // 1. ������� ������������������� ��-� ��� �������  arrKInp

	 long double arrLamb[3] ={0.} ; // ������ ������ ���������
	/* long double a = - ( arrKInp [0] + arrKInp [4]  + arrKInp [8] ) ;
	 long double b =  arrKInp [0] * arrKInp [4] +  arrKInp [0] * arrKInp [8] + arrKInp [4] * arrKInp [8]
			   - arrKInp [5] * arrKInp [5] - arrKInp [1] * arrKInp [1] - arrKInp [2] * arrKInp [2] ;
	 long double c = -( arrKInp [0] * arrKInp [4]  * arrKInp [8] + 2 *  arrKInp [5] * arrKInp [1] * arrKInp [2]
			 - arrKInp [0] * arrKInp [5] * arrKInp [5] - arrKInp [8] * arrKInp [1] * arrKInp [1] - arrKInp [4] * arrKInp [2] * arrKInp [2]) ;
	  */
	 long   double a = - ( arrKInp [0] + arrKInp [4]  + arrKInp [8] ) ;
	 long double b =  arrKInp [0] * arrKInp [4] +  arrKInp [0] * arrKInp [8] + arrKInp [4] * arrKInp [8]
               - arrKInp [5] * arrKInp [7] - arrKInp [3] * arrKInp [1] - arrKInp [2] * arrKInp [6] ;
     long double c = arrKInp [0] * arrKInp [5] * arrKInp [8]
             + arrKInp [4] * arrKInp [2] * arrKInp [6]
             +arrKInp [4] * arrKInp [3] * arrKInp [1]
              - arrKInp [0] * arrKInp [4]  * arrKInp [8]
              -arrKInp [2] * arrKInp [3]  * arrKInp [7]
			 -arrKInp [3] * arrKInp [5]  * arrKInp [7];
	 int irez = Cubic(arrLamb, a, b, c) ;
	  ///


	 if (irez == 0) return -1 ; // ������, �������� arrKInp �� ������������-�����������

	 if( ( arrLamb[0] < -EPS) || ( arrLamb[1] < -EPS)|| ( arrLamb[2] < -EPS) ) return -1 ; // ������, �������� arrKInp �� ������������-�������������
	 int quantNull = 0 ; // �-�� ������� ������ �����
	  for (int i =0; i < 3; i++)
	  {
		  long double dd = fabsl(arrLamb[i]) ;
		  if( (dd - EPS) < 0.)
		  {
		   quantNull++;
          }

	  }
	 ///

	//qsort(arrLamb, 3, sizeof(long double),&(fCompare)); // ���������� � ������� �����������
	 int iRang = 0 ;
	 int iReturn = -100;
	 for (int i =0; i < 3; i++)
	 {
	   if (arrLamb [i] > EPS) iRang++ ;
	   else arrLamb [i] = 0 ;

	 }

	 long double mtrA[9] = {0}, det0,det1,det2, arrfTr[9] = {0} ;
	 int num1 = -1, num2 = -1 ;
	switch (irez)
	{
	case 3:
		for (int i = 0; i < 3; i++)
		{
			memcpy(mtrA, arrKInp, 9 * sizeof(long double)) ;
			for (int j = 0; j < 3; j++) mtrA [ 3 * j + j ] -=  arrLamb[ i ];
			det0 =  mtrA[0] * mtrA[4 ] - mtrA[1] * mtrA[1 ] ;
			if (fabsl(det0) < EPS)
			{
				det1 = mtrA[1] * mtrA[5 ] - mtrA[2] * mtrA[4 ] ;
				if (fabsl(det1) < EPS)
				{
				det2               =  mtrA[0] * mtrA[8 ] - mtrA[2] * mtrA[2 ] ;
				arrfTr[3 * i  + 1] = 1 ;
				arrfTr[3 * i  ]    = ( - mtrA[1] * mtrA[8 ] + mtrA[2] * mtrA[5 ]) / det2;
				arrfTr[3 * i  + 2] = ( - mtrA[0] * mtrA[5 ] + mtrA[1] * mtrA[2 ]) / det2;
				}
				else
				{
				arrfTr[3 * i     ] = 1 ;
				arrfTr[3 * i + 1 ] = ( - mtrA[0] * mtrA[5 ] + mtrA[2] * mtrA[1 ]) / det1;
				arrfTr[3 * i + 2 ] = -det0 / det1;
				}

			}
			else
			{
			arrfTr[3 * i  + 2 ] = 1 ;
			arrfTr[3 * i   ]    = ( -mtrA[2] * mtrA[4 ] + mtrA[5] * mtrA[1 ]) / det0; ;
			arrfTr[3 * i + 1  ] = ( -mtrA[0] * mtrA[5 ] + mtrA[2] * mtrA[1 ]) / det0; ;

			}
		}
		if (iRang == 2)
		{
		  iReturn = 23 ;
		}
		else
		{
		  iReturn = 33 ;
		}

	break;
	case 2:
		memcpy(mtrA, arrKInp, 9 * sizeof(long double)) ;
		for (int j = 0; j < 3; j++) mtrA [ 3 * j + j ] -=  arrLamb[ 0];
		det0 =  mtrA[0] * mtrA[4 ] - mtrA[1] * mtrA[1 ] ;
		//������ ������ ������� ��� ���������� ������ �����
		if (fabs(det0) < EPS)
		{
			det1 = mtrA[1] * mtrA[5 ] - mtrA[2] * mtrA[4 ] ;
			if (fabsl(det1) < EPS)
			{
			det2               =  mtrA[0] * mtrA[8 ] - mtrA[2] * mtrA[2 ] ;
			arrfTr[ 1] = 1 ;
			arrfTr[0  ]    = ( - mtrA[1] * mtrA[8 ] + mtrA[2] * mtrA[5 ]) / det2;
			arrfTr[ 2] = ( - mtrA[0] * mtrA[5 ] + mtrA[1] * mtrA[2 ]) / det2;
			}
			else
			{
			arrfTr[0    ] = 1 ;
			arrfTr[ 1 ] = ( - mtrA[0] * mtrA[5 ] + mtrA[2] * mtrA[1 ]) / det1;
			arrfTr[ 2 ] = -det0 / det1;
			}

		}
		else
		{
		arrfTr[ 2 ] = 1 ;
		arrfTr[0   ]    = ( -mtrA[2] * mtrA[4 ] + mtrA[5] * mtrA[1 ]) / det0; ;
		arrfTr[ 1  ] = ( -mtrA[0] * mtrA[5 ] + mtrA[2] * mtrA[1 ]) / det0; ;

		}

		///

		// ���������� �������������� �������
		  // ����� ������� 2-� ��������� ����������� ������  arrfTr

		  for (int i = 0; i < 3; i++)
		  {
			if (fabsl(arrfTr[i]) > EPS )
			{
			  if (num1 < 0)
			  {
				num1 = i ;
			  }
			  else
			  {
				  num2 = i ;
				  break ;
              }
			}
		  }
		///

		if (num2 == -1)
		{ // � ������� arrfTr ������� ���� ��������� ����������
		  num2 = (num1 + 1)% 3 ;
		  num1 = (num1 + 2)% 3 ;
		  arrfTr[3 + num1] = 1 ;
		  arrfTr[6 + num2] = 1 ;

		}
		else
		{
		  arrfTr[3 + num1] =  arrfTr[ num1] ;
		  arrfTr[3 + num2] = -arrfTr[num2]  ;
			OuterProduct(arrfTr , &arrfTr [3] , &arrfTr [6]) ;

        }

		 switch(iRang )
		 {
			 case 1:
			 iReturn = 12 ;
			 break;
			 case 2:
			 iReturn = 22 ;
			 break;
			 case 3:
			 iReturn = 32 ;
			 break;
			 case 0:
			 default:
			 return -5 ;

		 }
    break ;
	case 1:
		switch (quantNull)
		{
		case 0:
		  arrfTr[0    ] = 1 ;
		  arrfTr[3 + 1] = 1 ;
		  arrfTr[6 + 2] = 1 ;
		break ;
		case 3:
	//	ShowMessage(L"Error -3") ;
		return -3;

		case 2:
		case 1:
		default:
	 //	ShowMessage(L"Error -1") ;
		return -1;
		;
		}
		switch(iRang )
		 {
			 case 3:
			 iReturn = 31 ;
			 break;
			 case 0:
			 iReturn = 0 ;
			 break;
			 case 1:
			 case 2:
			 default:
			 return -6 ;

		 }

	break ;
	default:
  //	ShowMessage(L"Error -2") ;
	 return -2;
	;
	 }

	// ���������� �� ��������
	bool bPr = true;
	while(bPr)
	{
		bPr = false ;
		for (int i = 0; i < 2; i++)
		{
		  if (arrLamb [ i + 1] > arrLamb [ i ])
		  {
			swap(&arrLamb [ i + 1], & arrLamb [ i ] ) ;
			swap_vect(&arrfTr [ 3 * (i + 1)], & arrfTr [3 *  i ], 3 ) ;
			bPr = true ;
			break ;
		  }
		}
	}
///
	 // �������������, ���������������
	  memset(arrMatrProperNumb , 0 , 9 * sizeof(long double)) ;
	 for (int i = 0; i < 3; i++)
	 {
	  NormalizeVect3( &arrfTr [ i * 3]) ;
	  arrMatrProperNumb [ 3 * i + i ] =  arrLamb [ i ];
	 }
	 MatrTransp(arrfTr, 3, 3, arrMatrProperVect ) ;


	return iReturn ;
}

void   swap(long double *a0, long double *a1)
{
	long double b = *a1;
	*a1 = *a0 ;
	*a0 = b ;
}

void   swap_vect(long double *pa0, long double *pa1, const int len)
{
	for (int i = 0; i < len; i++) swap( &pa0[i],&pa1[i]);


}

long double   NormVect3(long double *p)
{
	return sqrtl(p[0] * p [0] + p[1] * p [1] + p[2] * p [2] ) ;
}

long double   NormVect2(long double *p)
{
	return sqrtl(p[0] * p [0] + p[1] * p [1]  ) ;
}

void   NormalizeVect3(long double *p)
{
	long double valNorm =  NormVect3(p) ;
	for (int i = 0; i < 3; i++)
	{
	  p[i] = p[i] / valNorm ;
	}
}

// ����� ����������� ��������  � ����������� ����� ������������ ������������ ������� arrKInp
//  arrV - ������� ����������� ��������
//  arrLamb - ��������������� ������� ����������� �����
int  CalcProperVectors2(long double * arrKInp,long double *arrV , long double *arrLamb)
{

	memset(arrLamb,0,4 * sizeof(long double)) ;
	long double a  = 1;
	long double b = -(arrKInp[0] + arrKInp [3]);
	long double c = arrKInp[0] * arrKInp [3] - arrKInp[1] * arrKInp [1] ; ;
	TCompLong x0,x1;
	int i0 = SolvEq2( a, b, c,x0,x1) ;
	long double arrx[2];
	arrx[0] = x0.m_Re;
	arrx[1] = x1.m_Re;

	// ����������� ������
	try
	{
	for (int i =0; i < 2; i++)
	{
	  if (fabsl (arrx[i]) < 0.0001 )
	  {
		int j = (i + 1)%2 ;
		long double valLambPos = arrx[j] ;
		if (valLambPos < 0) return -1 ;


		if ((fabsl( arrKInp[0] - valLambPos) > 0.000001) || (fabsl( arrKInp[1] ) > 0.000001))
		{
		  long double arre[2] = {0.} ;
		  long double d = sqrtl( arrKInp[1]* arrKInp[1] + (arrKInp[0] - valLambPos) * (arrKInp[0] - valLambPos)) ;
		  arre[0] =  -arrKInp[1] /d;
		  arre[1] =  (arrKInp[0] - valLambPos)/d  ;
		  arrV [0] = arre[0] ;
		  arrV [2] = arre[1] ;
		  arrV [1] = -arre[1] ;
		  arrV [3] =  arre[0] ;
		  arrLamb [0] = valLambPos ;

		  return 0 ;
		}
		else
		{
		   return -1 ;

		}
	  }
	}
	}
	catch(...)
	{
	 //	ShowMessage( L"TFragmTraj::CalcProperVectors2 Error No 0");
		return -1;
	}
  ///

   try
   {
	switch(i0)
	{
		case 0:

		for (int i = 0; i < 2 ; i++)
		{
		  if  ( arrx[i] < 0)
		  {
			  int j = (i + 1)%2 ;
			  long double arrProp[2] ={0},arrPropTemp[4] = {0} ;
			  if ( fabsl (arrKInp[0] - arrx[j])> EPS)
			  {
				arrProp[0] = - arrKInp[1]/ (arrKInp[0] - arrx[j]);
				arrProp[1] = 1 ;
			  }
			  else
			  {
				arrProp[0] = 1 ;
				arrProp[1] =  -  (arrKInp[0] - arrx[j])/arrKInp[1];
              }
			  long double d = sqrtl(arrProp[0] * arrProp[0] + arrProp[1] * arrProp[1]) ;
			  for (int ii = 0; ii < 2; ii++) arrProp[ii] = arrProp[ii] / d ;


			  MtrxMultMatrxTransp(arrProp,2, 2, arrProp,2, arrPropTemp) ;
			  MatrxMultScalar(arrPropTemp, 2 , 2, arrx[j],arrKInp);
			  return 1 ;
		  }
		}
		for (int i = 0; i < 2; i++)
		{
		 arrLamb[i * 2 + i] = arrx[i];
		 if (fabsl(arrKInp[0] - arrx[i] ) > 0.00000000001 )
		 {
		   arrV [ i ] = - arrKInp[1] / (arrKInp[0] - arrx[i] );
		   arrV [ 2 + i ] = 1 ;
		   long double r = sqrtl (arrV [ i ] *  arrV [ i ] + arrV [ 2 + i ] * arrV [ 2 + i ] );
		   arrV [ i ] = arrV [ i ] /r ;
		   arrV [ 2 + i ] = arrV [ 2 + i ] /r ;

		 }
		 else
		 {
		   arrV [ i ]     = 1 ;
		   arrV [ 2 + i ] = 0 ;
         }
		}
		break ;
		case 1:
		arrV [0] = 1 ;
		arrV [1] = 0 ;
		arrV [2] = 0 ;
		arrV [3] = 1 ;
		arrLamb [0] = arrKInp[0] ;
		arrLamb [1] = 0 ;
		arrLamb [2] = 0 ;
		arrLamb [3] = arrKInp[0] ;
		break ;

		default:

		return 1;

	}
	}
	catch(...)
	{
	//	ShowMessage( L"TFragmTraj::CalcProperVectors2 Error No 2");
		return -1;

    }
	return 0 ;
}
// ������� ������� �������� ��������� ������� �������  A * x = B
// Returns:
// true , ���� det(A) != 0
// false - � ��������� ������
bool   SolvLinEq2(long double *arrA, long double *arrB,long double *arrX)
{
   long double det = arrA [ 0] * arrA [ 3] - arrA [ 1] * arrA [ 2] ;
   if (fabsl(det) < EPS) return false ;
   arrX [ 0] =   ( arrB [ 0] *  arrA [ 3] - arrA [ 1] * arrB [ 1] ) / det ;
   arrX [ 1] = ( arrA [ 0] * arrB [ 1]  - arrA [ 2] * arrB [ 0] ) / det ;

  return true;
}

long double NormVect(long double *arr, const int lenarr)
{
  long double sum = 0;
  for (int i = 0; i < lenarr; i++)
  {
	sum += arr[i] * arr[i] ;
  }
  return sqrtl(sum) ;
}


// ������������ ��������������� ������ � �������� � ������  iRowElem, iColElem
// parrMtrxInp - �������� �������
// ndim - �����������  �������� �������
// iRowElem , jColElem  - ����� ��������
// OUTPUT:
// parrMinor - �������������� ����� �����������  [ndim -1]x [ndim -1]
void buildAddMinor(long double *parrMtrxInp, const int ndim
	, const int iRowElem,const int jColElem, long double *parrMinor)
{
	int iMinor = 0 ;
	int jMinor = 0 ;
	for (int i = 0; i < ndim; i++)
	{
	 if (i == iRowElem) continue;

	 for (int j = 0; j < ndim; j++)
	 {
	   if (j == jColElem) continue ;

		parrMinor [ (ndim-1) * iMinor + jMinor] =  parrMtrxInp [ ndim * i + j] ;
	   jMinor++;
	 }
	 iMinor++;
	 jMinor = 0;
	}
}

// ������������� �������� 4�4
long double   calcDet4( long double *arrInp)
{
	long double det = 0.;
	long double parrMinor [9] = {0.} ;
	long double valIndex = 1.;
	for (int j = 0; j < 4; j++)
	{
	  buildAddMinor(arrInp, 4, 0, j, parrMinor) ;
	  det += valIndex * arrInp[j] * calcDet3( parrMinor) ;
	  valIndex = -valIndex;
	}
	return det;
}

//  ��������� ������ 4�4
//
bool   InverseMtrx4(long double *arrInp, long double *arrOut)
{
	long double det =   calcDet4( arrInp);
	if (fabs(det) < 1.E-15) return false ;
	long double parrMinor [9] = {0.} ;
	long double valIndexI = 1.;
	long double valIndexJ = 1.;
	for (int i = 0; i < 4; i++)
	{
	for (int j = 0; j < 4; j++)
	{
	 buildAddMinor(arrInp, 4, j,i, parrMinor) ;
	 arrOut[ 4 * i + j] =  valIndexI * valIndexJ * calcDet3( parrMinor)/ det;
	 valIndexJ = -valIndexJ;
	}
	valIndexI = -valIndexI;
	valIndexJ = 1.;
	}
	return true;
}

// ������������� �������� 5x5
long double   calcDet5( long double *arrInp)
{
	long double det = 0.;
	long double parrMinor [16] = {0.} ;
	long double valIndex = 1.;
	for (int j = 0; j < 5; j++)
	{
	  buildAddMinor(arrInp, 5, 0, j, parrMinor) ;
	  det += valIndex * arrInp[j] * calcDet4( parrMinor) ;
	  valIndex = -valIndex;
	}
	return det;
}
//  ��������� ������ 5�5
//
bool   InverseMtrx5(long double *arrInp, long double *arrOut)
{
	long double det =   calcDet5( arrInp);
	if (fabs(det) < 1.E-15) return false ;
	long double parrMinor [16] = {0.} ;
	long double valIndexI = 1.;
	long double valIndexJ = 1.;
	for (int i = 0; i < 5; i++)
	{
	for (int j = 0; j < 5; j++)
	{
	 buildAddMinor(arrInp, 5, j,i, parrMinor) ;
	 arrOut[ 5 * i + j] =  valIndexI * valIndexJ * calcDet4( parrMinor)/ det;
	 valIndexJ = -valIndexJ;
	}
	valIndexI = -valIndexI;
	valIndexJ = 1.;
	}
	return true;
}

// ��������� ������� 2-�� �������
// Returns:
// true , ���� det(A) != 0
// false - � ��������� ������
bool   InverseMtrx2(long double *arrA,long  double *arrOut)
{
	long  double det = arrA [ 0] * arrA [ 3] - arrA [ 1] * arrA [ 2] ;
   if (fabs(det) < EPS) return false ;
   arrOut [ 0] =  arrA [ 3] / det ;
   arrOut [ 1] = -arrA [ 1] / det ;
   arrOut [ 2] = -arrA [ 2] / det ;
   arrOut [ 3] =  arrA [ 0] / det ;
  return true;
}


// ������������� �������� 6x6
long double    calcDet6( long double  *arrInp)
{
	long double  det = 0.;
	long double  parrMinor [25] = {0.} ;
	long double  valIndex = 1.;
	for (int j = 0; j < 6; j++)
	{
	  buildAddMinor(arrInp, 6, 0, j, parrMinor) ;
	  det += valIndex * arrInp[j] * calcDet5( parrMinor) ;
	  valIndex = -valIndex;
	}
	return det;

}
// ������������� �������� 7x7
long double    calcDet7( long double  *arrInp)
{
	long double  det = 0.;
	long double  parrMinor [36] = {0.} ;
	long double  valIndex = 1.;
	for (int j = 0; j < 7; j++)
	{
	  buildAddMinor(arrInp, 7, 0, j, parrMinor) ;
	  det += valIndex * arrInp[j] * calcDet6( parrMinor) ;
	  valIndex = -valIndex;
	}
	return det;

}

//  ��������� ������ 6x6
//
bool   InverseMtrx6(long double  *arrInp, long double  *arrOut)
{
	long double  det =   calcDet6( arrInp);
	if (fabs(det) < 1.E-15) return false ;
	long double  parrMinor [25] = {0.} ;
	long double  valIndexI = 1.;
	long double  valIndexJ = 1.;
	for (int i = 0; i < 6; i++)
	{
	for (int j = 0; j < 6; j++)
	{
	 buildAddMinor(arrInp, 6, j,i, parrMinor) ;
	 arrOut[ 6 * i + j] =  valIndexI * valIndexJ * calcDet5( parrMinor)/ det;
	 valIndexJ = -valIndexJ;
	}
	valIndexI = -valIndexI;
	valIndexJ = 1.;
	}
	return true;
}


/// ������������� �������� 8x8
long double    calcDet8( long double  *arrInp)
{
	long double  det = 0.;
	long double  parrMinor [49] = {0.} ;
	long double  valIndex = 1.;
	for (int j = 0; j < 8; j++)
	{
	  buildAddMinor(arrInp, 8, 0, j, parrMinor) ;
	  det += valIndex * arrInp[j] * calcDet7( parrMinor) ;
	  valIndex = -valIndex;
	}
	return det;

}




//  ��������� ������ 8x8
//
bool   InverseMtrx8(double *arrInp, double *arrOut)
{
	double det =   calcDet8( arrInp);
	if (fabs(det) < 1.E-15) return false ;
	double parrMinor [49] = {0.} ;
	double valIndexI = 1.;
	double valIndexJ = 1.;
	for (int i = 0; i < 8; i++)
	{
	for (int j = 0; j < 8; j++)
	{
	 buildAddMinor(arrInp, 8, j,i, parrMinor) ;
	 arrOut[ 8 * i + j] =  valIndexI * valIndexJ * calcDet7( parrMinor)/ det;
	 valIndexJ = -valIndexJ;
	}
	valIndexI = -valIndexI;
	valIndexJ = 1.;
	}
	return true;
}


// ���������� ��������� ������� ������������ �����������
//
bool   InverseMtrx( double *arrInp, const int nDimen,  double *arrOut)
{
  switch (nDimen)
  {
  case 1:
  if (fabs(*arrInp) > 0.000000001)
  {
	*arrOut = 1./ (*arrInp);
	return true;
  }
  else
  {
	  return false;
  }
  case 2:
  return InverseMtrx2(arrInp, arrOut);

  case 3:
  return InverseMtrx3(arrInp, arrOut);


  case 4:
  return InverseMtrx4(arrInp, arrOut);


  case 5:
  return InverseMtrx5(arrInp, arrOut);


  case 6:
  return InverseMtrx6(arrInp, arrOut);

  case 7:
  return InverseMtrx7(arrInp, arrOut);

  case 8:
  return InverseMtrx8(arrInp, arrOut);

  default: return false ;
      ;
  }
 // return false ;
}


bool GaussMeth(long double *parrA,int nDimA,  long double * parrB, long double *parrRez)
 {
	 long double *arrAUn  = new long double [ nDimA * (nDimA +1)];
	 int *iarrNums  = new int [ nDimA ];
	 for (int i =0; i < nDimA; i++)
	 {
	   memcpy( &arrAUn[ i *  (nDimA +1)], & parrA[i * nDimA ], sizeof(long double) * nDimA);
	   arrAUn[ i *  (nDimA +1) + nDimA] =  parrB[i];
	   iarrNums [i] = i;
	 }
	 for (int i = 0; i < nDimA; i++)
	 {
		int numVedStroka = fncFindNumVedushiaStroka(arrAUn,  nDimA, i);
		if (numVedStroka == -1)
		{
			delete arrAUn;
			delete iarrNums;
			return false;
		}
		Swap (iarrNums,1, numVedStroka, i);
		Swap (arrAUn,nDimA +1, numVedStroka, i);
		for (int j = i + 1; j < nDimA; j++)
		{
		  if (fabs(arrAUn [j *(nDimA + 1) +i]) < 10. * DBL_MIN)
		  {
			continue;
		  }
		  if (!fncPodstanovka( arrAUn, nDimA,i, j) )
		  {
			delete arrAUn;
			delete iarrNums;
			return false;
		  }
		}
	 }

	 long double *arrRezTemp = new long double [nDimA];
	 memset(arrRezTemp, 0, nDimA * sizeof(long double));
	 arrRezTemp[nDimA -1]= arrAUn[ (nDimA-1) *  (nDimA +1) + nDimA]/  arrAUn[ (nDimA-1) *  (nDimA +1) + nDimA -1];
	 for (int i = nDimA -2; i >= 0; i--)
	 {
	   long double sum = arrAUn[ i *  (nDimA +1) + nDimA];
	   for (int j =i+1; j < nDimA ; j++)
	   {
		 sum -=  arrRezTemp[j] * arrAUn[ i *  (nDimA +1) + j];
	   }
	   arrRezTemp[i] = sum/  arrAUn[ i *  (nDimA +1) + i];
	 }
	 for (int i =0; i < nDimA; i++)
	 {
	   parrRez [i] =  arrRezTemp[iarrNums[i]];
	 }
	 delete arrAUn;
	 delete iarrNums;
	 delete arrRezTemp;
	 return true;
 }

bool fncPodstanovka( long double *arrAUn, int nDimA, int i, int j)
{

	long double temp = -  arrAUn [i *(nDimA + 1) +i] /arrAUn [j *(nDimA + 1) +i];
	arrAUn [j *(nDimA + 1) +i] = 0.;
	for (int k =i + 1; k < (nDimA + 1); k++)
	{
	  arrAUn [j *(nDimA + 1) +k] = arrAUn [i *(nDimA + 1) +k] + temp * arrAUn [j *(nDimA + 1) +k] ;
	}
	return true;
}

int fncFindNumVedushiaStroka(long double *arrAUn, int  nDimA, int i0)
{
	int irez = -1;
	for (int i = i0; i < nDimA; i++)
	{
	  if (fabs(arrAUn [i *( nDimA +1) + i0]) > 10. * DBL_MIN )
	  {
		return i;
	  }
	}
	return irez ;
}



void Swap (long double *parrA, int iNumCols, int j,int  i)
{
	long double *parrTemp  = new long double [iNumCols];
	memcpy(parrTemp ,&parrA[j * iNumCols], iNumCols * sizeof(long double));
	memcpy(&parrA[j * iNumCols] ,&parrA[i * iNumCols], iNumCols * sizeof(long double));
	memcpy(&parrA[i * iNumCols] ,parrTemp, iNumCols * sizeof(long double));

	delete parrTemp ;
}

#pragma package(smart_init)
