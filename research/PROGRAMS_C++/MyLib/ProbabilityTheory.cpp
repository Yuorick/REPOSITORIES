//---------------------------------------------------------------------------


#pragma hdrstop

#include "ProbabilityTheory.h"

#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "MatrixProccess.h"


//---------------------------------------------------------------------------
// ¬ычисление веро€тности того, что число успехов   из последовательности
// QUantIsp независимых испытаний биномиального распределени€ случ величины
// будет меньше VAlNesessarySuccess
//INPUT:
// QUantIsp  - обющее к-во испытаний
// VAlP - веро€тность успеха в одном испытании
// VAlNesessarySuccess - необходимое число успехов
double __fastcall TProbabilityTheory::funcDistributionOfProbabilitiesForBinom (const int QUantIsp,
   const double VAlP, const double VAlNesessarySuccess )
{
double valSig = sqrt(((double)QUantIsp) * VAlP * (1. - VAlP));
double valMean = ((double)QUantIsp) * VAlP /valSig;
double valx =  VAlNesessarySuccess / valSig - valMean;
double val_return =  fncGaussDistrib01_Dim1 ( valx);

return val_return;
}

// функци€ распределени€ одномерного гауссовго распределени€  c нулевым средним и един дисперсией
// INPUT:
// VAlx  - аргумент
// ¬ычисл€ет веро€тность тогоЅ что сл величина будет меньше  VAlx
double __fastcall TProbabilityTheory::fncGaussDistrib01_Dim1 (const double VAlx)
{
	double valGaussIntegral = calcGaussIntegral( VAlx );
	double valreturn = 0.5 +  valGaussIntegral;
 //	if (VAlx < 0.)
 //	{
 //	 valreturn = 0.5 -  valGaussIntegral;
 //	}
	return valreturn ;
}


// вычисление интеграла от функции exp(-x*x/2.)/sqrt(2.* M_PI)
// на отрезке [0; VAlx]
double __fastcall TProbabilityTheory::calcGaussIntegral(const double VAlx )
{
	if (fabs(VAlx) > 3.)
	{
	 return 0.5 * SIGNUM(VAlx);
	}
	double valStep = 0.0001;
	int iC = VAlx / valStep;
	double valSum = 0.;
	for (int i =0; i < iC; i++)
	{
	  double temp = ((double)i) * valStep;
	  valSum += exp(-temp * temp /2.);
	}
	 valSum *=  valStep / sqrt(2. * M_PI);
	return valSum;
}


// ¬ычисление веро€тности к-ва успехов биномиального распределени€
// INPUT:
//VAl_q - веро€тность успеха обдного испытани€
// QUantM - число испытаний
// QUantN - число успехов
// возвращает веро€тность осуществлени€ QUantN  успешных испытаний
double __fastcall TProbabilityTheory::calcBinomSuccessProbability(const double VAl_q, const int QUantM, const int QUantN )
{
	const long double VAl_ql = VAl_q;
	const long double VAl_M =  (long double) QUantM;
	const long double VAl_N =  (long double) QUantN;
	long double valSumLn = VAl_N *  logl(VAl_ql) + (VAl_M -VAl_N) * logl( 1. - VAl_ql);
	for (int i = 0; i < QUantN; i++)
	{
	 valSumLn += logl(VAl_M - (long double)i) - logl((long double)(i + 1));
	}
	long double temp = expl(valSumLn);
	double temp1 = (double) temp;
	return temp1;

}

//------------------------------------------------------------------
// ¬ычисдление интьеграла от нормальной плотности одномерного распределени€
// VAla - среднее(мат ожидание)
// VAlSig - — «
// VAlX0, VAlX0 - нижний и верхний пределы интегрировани€
// VAlStepIntegr - шаг интегрировани€
//
//
double TProbabilityTheory::calcIntegralNormalDensity(const double VAla, const double VAlSig
	, const double VAlX0, const double VAlX1,  const double VAlStepIntegr)
{
 double iC = (VAlX1 - VAlX0) / VAlStepIntegr;
 double sum = 0.;
 for (int i =0 ; i < iC; i++)
 {
	double temp = (VAlX0 + ((double)i) * VAlStepIntegr -VAla) /VAlSig;
	sum +=  exp(- temp * temp / 2.);
 }
 sum = sum * VAlStepIntegr / sqrt(2. * M_PI) / VAlSig;
 return sum ;
 }

//----------------------------------------------------------------------------------
// вычисление двойного иньтеграла от 2-х мероной плотности нормальеого распределени€
// по кругу с центрлм в точке  (VAlX, VAlY) радиусом  VAlKillingRange
// VAlStepIntegr - шаг сетки интегрроовани€
// arrElK - коррел€ционна€ матрица
double  TProbabilityTheory::calcIntegralNormalDensity_2D(const double VAlX, const double VAlY
			,double * arrElK, const double  VAlKillingRange,const double VAlStepIntegr)
{
 int iC = 2. * VAlKillingRange /VAlStepIntegr + 1.;
 double valh =  2. * VAlKillingRange / ((double)iC );
 double xll = VAlX - VAlKillingRange;
 double yll = VAlY - VAlKillingRange;
 double arr[2] ={0.};
 double  arrElKInv[4] = {0.};
 InverseMtrx2(arrElK, arrElKInv);
 double sum = 0.;
 for (int i = 0; i < iC; i++)
 {

 arr[0] = xll + ((double)i) * valh;
 double valDx = VAlX - arr[0];
	for (int j = 0; j < iC; j++)
	{
		arr[1] = yll + ((double)j) * valh;
		double valDy = VAlY - arr[1];
		if (sqrt(valDx * valDx +valDy * valDy)  > VAlKillingRange)
		{
			continue;
		}
		sum+= exp(-calcYT_D_Y(arr, arrElKInv, 2 )/2.);
  }
 }
 double det = arrElK[0] * arrElK[3] - arrElK[1]* arrElK[2];
 double xreturn = sum/ 2. / M_PI / sqrt(det) * valh * valh;
return xreturn;
}

// дисперси€ равномерного распределени€ на отрезке VAlL
double  TProbabilityTheory::calcDispRavnomern(const double VAlL)
{
	return  VAlL * VAlL / 12.;
}
#pragma package(smart_init)
