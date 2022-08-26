//---------------------------------------------------------------------------


#pragma hdrstop
 #include <vcl.h>
#include <stdio.h>
#include "Measure.h"


//---------------------------------------------------------------------------
TMeasure::TMeasure()
{
		// ����������� ������� ��������� �������
	 mDimY =1;


//  1.����� �������� ������
	mTYZv =0.;
//  2.������ ������
	mparrYZv = new double[mDimY];

 // 3.�������������� ������� ���
	mparrBMO_K = new double[mDimY * mDimY];
	// 3.�������������� ������� ���
	mparrMMO_K = new double[mDimY * mDimY];



}

 TMeasure::~TMeasure()
{
	if(mparrYZv) delete mparrYZv ;
	mparrYZv = NULL ;
	if(mparrBMO_K) delete mparrBMO_K ;
	mparrBMO_K= NULL ;
	if(mparrMMO_K) delete mparrMMO_K ;
	mparrMMO_K= NULL ;
}


// ����������� �����������
 TMeasure ::TMeasure (const TMeasure &R)
 {
	// ����������� ������� ��������� �������
	 mDimY = R.mDimY;
	// ����������� ������� ���� � ������ ����� ��������� �������� �������
	mTYZv = R.mTYZv;


//  2.������ ������� ���������   �� ������ mTf
	if (mparrYZv  != NULL) mparrYZv  = NULL;
	if(R.mparrYZv  != NULL)
	{
		mparrYZv  = new double[R.mDimY  ];
		if(mparrYZv  == NULL)
		{
		ShowMessage(L"Not memory for mparrEstX") ;
		Abort() ;
		}
	memcpy( mparrYZv , R.mparrYZv , R.mDimY  * sizeof(double));
	}

 // 3.�������������� ������� ������ ����������

	if (mparrBMO_K != NULL) mparrBMO_K = NULL;
	if(R.mparrBMO_K != NULL)
	{
		mparrBMO_K = new double[R. mDimY * R. mDimY];
		if(mparrBMO_K == NULL)
		{
		ShowMessage(L"Not memory for mparrMtrxK") ;
		Abort() ;
		}
	memcpy( mparrBMO_K, R.mparrBMO_K, R. mDimY * R. mDimY  * sizeof(double));
	}
 // 3.�������������� ������� ������ ����������

	if (mparrMMO_K != NULL) mparrMMO_K = NULL;
	if(R.mparrMMO_K != NULL)
	{
		mparrMMO_K = new double[R. mDimY * R. mDimY];
		if(mparrMMO_K == NULL)
		{
		ShowMessage(L"Not memory for mparrMtrxK") ;
		Abort() ;
		}
	memcpy( mparrMMO_K, R.mparrMMO_K, R. mDimY * R. mDimY  * sizeof(double));
	}

 }

 // �������� ������������
 TMeasure TMeasure::operator=(TMeasure  R)
 {
	// ����������� ������� ��������� �������
	 mDimY = R.mDimY;
	// ����������� ������� ���� � ������ ����� ��������� �������� �������
	mTYZv = R.mTYZv;


//  2.������ ������� ���������   �� ������ mTf
	if (mparrYZv  != NULL) mparrYZv  = NULL;
	if(R.mparrYZv  != NULL)
	{
		mparrYZv  = new double[R.mDimY];
		if(mparrYZv  == NULL)
		{
		ShowMessage(L"Not memory for mparrEstX") ;
		Abort() ;
		}
	memcpy( mparrYZv , R.mparrYZv , R.mDimY  * sizeof(double));
	}

 // 3.�������������� ������� ������ ����������

	if (mparrBMO_K != NULL) mparrBMO_K = NULL;
	if(R.mparrBMO_K != NULL)
	{
		mparrBMO_K = new double[R. mDimY * R. mDimY];
		if(mparrBMO_K == NULL)
		{
		ShowMessage(L"Not memory for mparrMtrxK") ;
		Abort() ;
		}
	memcpy( mparrBMO_K, R.mparrBMO_K, R. mDimY * R. mDimY  * sizeof(double));
	}
 // 3.�������������� ������� ������ ����������

	if (mparrMMO_K != NULL) mparrMMO_K = NULL;
	if(R.mparrMMO_K != NULL)
	{
		mparrMMO_K = new double[R. mDimY * R. mDimY];
		if(mparrMMO_K == NULL)
		{
		ShowMessage(L"Not memory for mparrMtrxK") ;
		Abort() ;
		}
	memcpy( mparrMMO_K, R.mparrMMO_K, R. mDimY * R. mDimY  * sizeof(double));
	}


	return *this ;
 }

 TMeasure::TMeasure(const int DimY, const double ValDispBMO )
 {
  mDimY = DimY ;
  mTYZv =0.;
  mparrYZv  = new double[mDimY];
  memset(mparrYZv, 0,  DimY * sizeof(double));

  mparrBMO_K = new double[ mDimY * DimY];
  memset(mparrBMO_K, 0,  mDimY * DimY * sizeof(double));
  for (int i =0; i < mDimY; i++)
  {
   mparrBMO_K[ i * mDimY + i] =  ValDispBMO;
  }

  mparrMMO_K = new double[ mDimY *  mDimY];
  memset(mparrMMO_K, 0,  mDimY * DimY * sizeof(double));

 }

  TMeasure::TMeasure(const int DimY, const double ValDispBMO , const double ValDispMMO )
 {
  mDimY = DimY ;
  mTYZv =0.;
  mparrYZv  = new double[mDimY];
  memset(mparrYZv, 0,  DimY * sizeof(double));

  mparrBMO_K = new double[ mDimY * DimY];
  memset(mparrBMO_K, 0,  mDimY * DimY * sizeof(double));

  mparrMMO_K = new double[ mDimY *  mDimY];
  memset(mparrMMO_K, 0,  mDimY * DimY * sizeof(double));

  for (int i =0; i < mDimY; i++)
  {
   mparrBMO_K[ i * mDimY + i] =  ValDispBMO;
   mparrMMO_K[ i * mDimY + i] =  ValDispMMO;
  }



 }

#pragma package(smart_init)
