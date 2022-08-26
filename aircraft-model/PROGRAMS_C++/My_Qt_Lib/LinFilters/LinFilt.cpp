//---------------------------------------------------------------------------



#include <stdio.h>
#include <string.h>
#include <math.h>

#include "LinFilt.h"
#include "Measure.h"
#include "Comp.h"
#include "MatrixProccess.h"

//---------------------------------------------------------------------------
TLinFilt::TLinFilt()
{
    // размерность вектора состо¤ни¤ объекта
    mDimX = 2;
    // размерность вектора наблюдений
    mDimY = 1;
    // ƒинамическа¤ информаци¤
//  1.¬рем¤ прив¤зки оценок вектора состо¤ни¤
     mTf = 0.;
//  2.ќценка вектора состо¤ни¤   на момент mTf
    mparrEstX = new double [mDimX];

 // 3.коррел¤ционна¤ матрица ошибок оценивани¤
    mparrMtrxK = new double [mDimX * mDimX];

}

 TLinFilt::~TLinFilt()
{
        if(mparrEstX) delete []mparrEstX ;
    mparrEstX = NULL ;
        if(mparrMtrxK) delete []mparrMtrxK ;
    mparrMtrxK = NULL ;
}

TLinFilt::TLinFilt(const int DimX, const int DimY, const double ValT,   double *parrEstX, double *parrMtrxK)
{
    // размерность вектора состо¤ни¤ объекта
    mDimX = DimX;

    // размерность вектора наблюдений
    mDimY = DimY;
    mTf = ValT;


  //	2.ќценка вектора состо¤ни¤   на момент mTf
    mparrEstX = new double [mDimX];
    memcpy(mparrEstX, parrEstX, mDimX * sizeof(double));

 // 3.коррел¤ционна¤ матрица ошибок оценивани¤
    mparrMtrxK = new double [mDimX * mDimX];
    memcpy(mparrMtrxK, parrMtrxK, mDimX *mDimX * sizeof(double));

}


int TLinFilt::fncStep( TMeasure Measure, const double Val_Tau2)
{
    const double Val_h = Measure.mTYZv - mTf;
    if (Val_h <= 0.)
    {
    return -1;
    }

    double *pMtrxA = new double [ mDimX * mDimX];
    double *pMtrxCorrU = new double [ mDimX * mDimX];
    double *pMtrxC = new double [ mDimY * mDimX];

    createMtrxA(Measure ,Val_h, pMtrxA);
    createMtrxCorrU(Measure ,Val_h, Val_Tau2, pMtrxCorrU);
    createMtrxC(Measure ,Val_h, pMtrxC);  // mDimY строк и mDimX стобцов

    // экстрапол¤ци¤
    double *parrXExtr = new double [ mDimX];
    double *parrKExtr = new double [ mDimX * mDimX];
    double *parrKTemp = new double [ mDimX * mDimX];
    double *parrKTemp0 = new double [ mDimX * mDimX];
    MtrxMultMatrx(pMtrxA, mDimX, mDimX, mparrEstX, 1, parrXExtr) ; // parrXExtr [mDimX]
    MtrxMultMatrx(pMtrxA, mDimX, mDimX,mparrMtrxK, mDimX, parrKTemp) ; // parrKTemp [mDimX * mDimX]
    MtrxMultMatrxTransp(parrKTemp, mDimX, mDimX, pMtrxA,mDimX, parrKTemp0 ) ; // parrKTemp0 [mDimX * mDimX]
    MtrxSumMatrx(parrKTemp0, pMtrxCorrU, mDimX, mDimX, parrKExtr) ;
    ///
    // учет замера
    //  расчет коэфф усилени¤
    double *parrP = new double [mDimX * mDimY];
    double *parrT0 =   new double [mDimX * mDimY];
    double *parrT1 =   new double [mDimY * mDimY];
    double *parrT2 =   new double [mDimY * mDimY];
    double *parrT2Inv =   new double [mDimY * mDimY];
    MtrxMultMatrxTransp(parrKExtr, mDimX, mDimX, pMtrxC,mDimY, parrT0 ) ;//parrT0[ mDimX *mDimY ] = parrKExtr * (pMtrxC)T
    MtrxMultMatrx(pMtrxC, mDimY, mDimX, parrT0, mDimY, parrT1) ;// parrT1[ mDimY *mDimY ]=pMtrxC *parrKExtr * (pMtrxC)T
    MtrxSumMatrx(parrT1, Measure.mparrBMO_K, mDimY, mDimY, parrT2) ;//parrT2 [ mDimY *mDimY ] =  pMtrxC *parrKExtr * (pMtrxC)T +N

    if(!InverseMtrx( parrT2, mDimY, parrT2Inv)) //parrT2Inv[ mDimY *mDimY ] =  (parrT2)Inv
    {
    return -2;
    }
    MtrxMultMatrx(parrT0, mDimX, mDimY, parrT2Inv,mDimY, parrP ) ;
    /// пересчет коррел матрицы
    double *parrT3 =   new double [mDimX * mDimX];
    double *parrT4 =   new double [mDimX * mDimX];

    MtrxMultMatrx(parrP, mDimX, mDimY, pMtrxC,mDimX, parrT3 ) ;//parrT3 [mDimX * mDimX] =  parrP * pMtrxC
    MtrxMultMatrx(parrT3, mDimX, mDimX, parrKExtr,mDimX, parrT4) ;//parrT4 [mDimX * mDimX] = parrP * pMtrxC *parrKExtr
    MtrxMinusMatrx(parrKExtr, parrT4, mDimX, mDimX,  mparrMtrxK);//mparrMtrxK [mDimX * mDimX] = parrKExtr - parrP * pMtrxC *parrKExtr

    ///

    // учет замера
    double *parrT6 =   new double [mDimY ];
    double *parrT7 =   new double [mDimY ];
    double *parrT8 =   new double [mDimX ];

    MtrxMultMatrx(pMtrxC, mDimY, mDimX, parrXExtr,1, parrT6) ;// parrT6 [mDimY ] = pMtrxC * parrXExtr
    MtrxMinusMatrx(Measure.mparrYZv  , parrT6, mDimY, 1, parrT7); //parrT7 = mparrYZv - pMtrxC * parrXExtr
    MtrxMultMatrx(parrP, mDimX, mDimY, parrT7,1, parrT8) ; // parrT8 = parrP *( mparrYZv - pMtrxC * parrXExtr)
    MtrxSumMatrx(parrXExtr, parrT8 ,mDimX, 1, mparrEstX) ; // mparrEstX = parrXExtr +  parrP *( mparrYZv - pMtrxC * parrXExtr)
    mTf = Measure.mTYZv;

    delete pMtrxA ;
    delete pMtrxCorrU ;
    delete pMtrxC ;

    delete parrP ;
    delete parrT0 ;
    delete parrT1 ;
    delete parrT2 ;
    delete parrT2Inv ;
    delete parrT3 ;
    delete parrT4 ;

    delete parrT6 ;
    delete parrT7;
    delete parrT8;
    delete parrKTemp0;
    return 0;
}

void TLinFilt::fncStabSolutionDynamic(const double Val_h, const double Val_Tau2
	 ,const TMeasure Measure)
{
	TMeasure MeasureTemp = Measure;
	MeasureTemp.mTYZv =  Val_h;

	fncStep(MeasureTemp, Val_Tau2);
	MeasureTemp.mTYZv +=  Val_h;
	double *arrKTemp = new double[mDimX * mDimX];
	memcpy(arrKTemp, mparrMtrxK, mDimX * mDimX * sizeof(double));
	for (int i = 0; i < 1000000; i++)
	{

	fncStep( MeasureTemp, Val_Tau2) ;
	bool bEnd = true;;
	 for (int j = 0; j< mDimX; j++)
	 {
	   if  ((fabs(mparrMtrxK[j * mDimX + j] - arrKTemp [j * mDimX + j])) > 0.00005 * mparrMtrxK[j * mDimX + j] )
	   {
		   bEnd = false;
		   break;
       }

	 }


	memcpy(arrKTemp,mparrMtrxK,sizeof(double) * 4) ;
	MeasureTemp.mTYZv +=  Val_h;
	if (bEnd)
	{
	  break;
	}

	}
	delete arrKTemp;
}



