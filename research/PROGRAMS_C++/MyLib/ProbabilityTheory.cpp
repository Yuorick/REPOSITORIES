//---------------------------------------------------------------------------


#pragma hdrstop

#include "ProbabilityTheory.h"

#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "MatrixProccess.h"


//---------------------------------------------------------------------------
// ���������� ����������� ����, ��� ����� �������   �� ������������������
// QUantIsp ����������� ��������� ������������� ������������� ���� ��������
// ����� ������ VAlNesessarySuccess
//INPUT:
// QUantIsp  - ������ �-�� ���������
// VAlP - ����������� ������ � ����� ���������
// VAlNesessarySuccess - ����������� ����� �������
double __fastcall TProbabilityTheory::funcDistributionOfProbabilitiesForBinom (const int QUantIsp,
   const double VAlP, const double VAlNesessarySuccess )
{
double valSig = sqrt(((double)QUantIsp) * VAlP * (1. - VAlP));
double valMean = ((double)QUantIsp) * VAlP /valSig;
double valx =  VAlNesessarySuccess / valSig - valMean;
double val_return =  fncGaussDistrib01_Dim1 ( valx);

return val_return;
}

// ������� ������������� ����������� ��������� �������������  c ������� ������� � ���� ����������
// INPUT:
// VAlx  - ��������
// ��������� ����������� ����� ��� �� �������� ����� ������  VAlx
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


// ���������� ��������� �� ������� exp(-x*x/2.)/sqrt(2.* M_PI)
// �� ������� [0; VAlx]
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


// ���������� ����������� �-�� ������� ������������� �������������
// INPUT:
//VAl_q - ����������� ������ ������� ���������
// QUantM - ����� ���������
// QUantN - ����� �������
// ���������� ����������� ������������� QUantN  �������� ���������
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
// ����������� ���������� �� ���������� ��������� ����������� �������������
// VAla - �������(��� ��������)
// VAlSig - ���
// VAlX0, VAlX0 - ������ � ������� ������� ��������������
// VAlStepIntegr - ��� ��������������
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
// ���������� �������� ���������� �� 2-� ������� ��������� ����������� �������������
// �� ����� � ������� � �����  (VAlX, VAlY) ��������  VAlKillingRange
// VAlStepIntegr - ��� ����� ��������������
// arrElK - �������������� �������
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

// ��������� ������������ ������������� �� ������� VAlL
double  TProbabilityTheory::calcDispRavnomern(const double VAlL)
{
	return  VAlL * VAlL / 12.;
}
#pragma package(smart_init)
