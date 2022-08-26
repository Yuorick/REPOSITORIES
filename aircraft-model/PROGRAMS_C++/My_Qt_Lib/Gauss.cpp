
#include "Gauss.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <time.h>
//#define DEBUG
const bool BDEbug = true;
//const bool BDEbug = false;
extern bool BEZ_SHUMOV;

 double  getGauss(const double   a, const  double   sig )

{
     if(BEZ_SHUMOV)
     {
         return a;
     }
     else
     {
        double   val = 0;
        int num = 12;

        for (int i = 0; i < num; i++)
        {
        val += ( Rand_())/(( double  )RAND_MAX) - 0.5;

        }


        val = val * sig + a ;

        return val ;
     }
}

// получение случайного числа с равномерным на отрезке [0;1]  распределением
 double   getRand01( )
{

  return (Rand_())/ (( double  )RAND_MAX);
}
//
 double  Rand_()
{
 	#ifndef DEBUG
 //	static int rand_init = 0;
 //	if( 0  == rand_init )
//	{
//		rand_init = time( NULL );
 //		srand((unsigned) time(NULL));
//	}

	#endif

  static  int rand_init = 0;
	if (!BDEbug )
	{
		if( 0  == rand_init )
	{
		rand_init = time( NULL );
		srand((unsigned) time(NULL));
	}

	}

 // return (double)rand();
 // double breturn =  (double)rand(); // для отладки
    if(BEZ_SHUMOV)
    {
        return 0.;
    }
    else
    {
      return (double)rand();
    }

}

// рекуррентный пересчет среднего и дисперсии
// для случайного процесса
//   INPUT:
// valX
// N
// valAver
// valAverSquare
// OUTPUT:
// valAver
// valAverSquare
void  CalcStatParams(const double valX, const int N
	, double &valAver,double  &valAverSquare)
{
   if (N == 1)
   {
	valAver = valX;
	valAverSquare = valX * valX;
	return;
   }
   valAver = valAver * (N - 1.)/ N + valX / N ;
    valAverSquare = valAverSquare * (N - 1.)/ N +  valX * valX / N ;
}

//----------------------------------------------------------------------------
 void __fastcall getGaussVector_dim3(double *arrMean, double *arrF, double *arrMtrxLamb, double *arrOut)
{
	double arrKsi[3] = {0.},  arrTemp[3]  ={0.};
	arrKsi[0] = getGauss(0., 1. )* sqrt(arrMtrxLamb[0]);;
	arrKsi[1] = getGauss(0., 1. )* sqrt(arrMtrxLamb[4]);
	arrKsi[2] = getGauss(0., 1. )* sqrt(arrMtrxLamb[8]);

	MtrxMultMatrx_(arrF,3, 3, arrKsi,1, arrTemp) ;
	MtrxSumMatrx_(arrMean, arrTemp,3, 1, arrOut) ;
}

//----------------------------------------------------------------------------
 void __fastcall getGaussVector(int iDim, double *arrMean, double *arrF, double *arrMtrxLamb, double *arrOut)
{
	 double *arrKsi = new double [iDim];
	 memset( arrKsi, 0,  iDim * sizeof(double));
	 double *arrTemp = new double [iDim];
	 memset( arrTemp, 0,  iDim * sizeof(double));
	 for (int i = 0; i < iDim; i++)
	 {
	   arrKsi[i] = getGauss(0., 1. )* sqrt(arrMtrxLamb[i + i *iDim ]);
	 }
	MtrxMultMatrx_(arrF,iDim, iDim, arrKsi,1, arrTemp) ;
	MtrxSumMatrx_(arrMean, arrTemp,iDim, 1, arrOut) ;
	delete arrTemp;
	delete arrKsi;
}

// умножение матрицы на матрицу
void __fastcall MtrxMultMatrx_(double *parrA,int nRowsA, int nColsA, double * parrB,int nColsB, double *parrRez)
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

void __fastcall MtrxSumMatrx_(double *pA, double * pB,int nRows, int nCols, double *pRez)
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

// Умножение транспонированной матрицы на  матрицу
void __fastcall MtrxTranspMultMatrx_(double *parrA,int nRowsA, int nColsA, double * parrB,int nColsB, double *parrRez)
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

}

