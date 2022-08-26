//---------------------------------------------------------------------------


//---------------------------------------------------------------------------


#pragma hdrstop

#include "Bomb.h"
#include <string.h>
#include <math.h>
#include "Atmosphere.h"

 // ������� �������� ������������� ��� � = 150/ 331
  const long double Cx00 =  0.2;
   // ����� ��� ���  mCx00
  const long double Ma�h00 =150./ 331.;
 // ������� �������� ������������� ��� � = 1.
  const long double Cx01 =  0.4;
   // ����� ��� ���  mCx00
  const long double  Ma�h01  = 1.2;

  // ������� ��������� ���� �������������� �������� (������� + �������)
  const long double CyFuz =  0.035;
  // ����������� ��������� ���� C�_���� = 0,8 * ( (D1/D)* (D1/D) -1) *2 * cos(tet1)  ** cos(tet1)
  // D = 0.45, D1 = 0.65,  tet1 = 30 ����. ������� � ������� ������� ������
  const long double CyKorma =  1.303;
   // ��������� �������� Cy����� ��� sqrt(1 -M*�) = 0.
  const long double Cy10 = 0.027;
 // ������� �������� ������������� ��� sqrt(1 -M*�) = Lamb_k
  const long double Cy11 =  0.018;



//---------------------------------------------------------------------------
TBomb::TBomb()
{
 // ������� ����� ��������
	 mMidSq = M_PI * 0.58 * 0.58 / 4.;//M_PI * 0.45 * 0.45 / 4.;
// ������� ������� ����������� �������� ( M_PI * D * L)
	 mSurfFuzSq = M_PI * sqrt( 4. * mMidSq / M_PI)* 4.63 ;
// �����
	mMass = 2200.;
// ���� ���� �����:
	mMaxAt = 20. * M_PI / 180. ;
 // ����� ��������� ����� L*L/S= 1.5*1.5/ (1.5*0.75)
	mLamb_k = 2.25/4.5;
 // �������� �������� S = 3. * 1.5
	mSq_k = 4.5;  //
 // ������� ������� ���
   //	mbHomeHead =true ;
 // ���
	 mHomingHead = THomingHead();
 // ���
	 mEbomb = TEbomb();
 //  ���� ������� ������� ������������ �� ���� �����
	 mTStabSyst = 0.5;

}

//---------------------------------------------------------------------------


// ����������� �����������
 TBomb ::TBomb (const TBomb &R)
 {
	 mMidSq = R.mMidSq ;
	 mSurfFuzSq = R.mSurfFuzSq;
	 mMass = R.mMass;
	 mMaxAt = R.mMaxAt ;
	 mSq_k = R.mSq_k ;
	 mLamb_k = R.mLamb_k ;
	 mHomingHead = R.mHomingHead ;
	 mEbomb = R.mEbomb ;
	 mTStabSyst = R.mTStabSyst ;

 }

 // �������� ������������
 TBomb TBomb::operator=(TBomb  R)
 {
	 mMidSq = R.mMidSq ;
	 mSurfFuzSq = R.mSurfFuzSq;
	 mMass = R.mMass;
	 mMaxAt = R.mMaxAt ;
	 mSq_k = R.mSq_k ;
	 mLamb_k = R.mLamb_k ;
	 mHomingHead = R.mHomingHead ;
	 mEbomb = R.mEbomb ;
	 mTStabSyst = R.mTStabSyst ;

	return *this ;
 }
// ���������������� ����������� 1
 TBomb::TBomb(
	 long double MidSq // ������� ����� ��������
	,long double SurfFuzSq // ������� ������� ����������� ��������
	,long double Mass// �����
	,long double MaxAt // ���� ���� �����
	,long double Lamb_k // ���������� �������
	,long double Sq_k // ������� �������, � ��
	,bool bHomeHead // ������� ������� ���
	,THomingHead HomingHead // ���
	,TEbomb Ebomb// ���
	)
{
	mMidSq = MidSq;
	mSurfFuzSq = SurfFuzSq ;
	mMass = Mass ;
	mMaxAt = MaxAt ;
	mLamb_k = Lamb_k ;
	mSq_k  = Sq_k  ;
	mHomingHead = HomingHead ;
	mEbomb  = Ebomb ;

}
// ���������������� ����������� 2
 TBomb::TBomb(THomingHead HomingHead )
{
	*this =  TBomb();
	mHomingHead = HomingHead ;

}

// ���������������� ����������� 3
 TBomb::TBomb(THomingHead HomingHead,TEbomb Ebomb )
{
	*this =  TBomb();
	mHomingHead = HomingHead ;
	mEbomb = Ebomb ;
}

// ���������������� ����������� 4
 TBomb::TBomb(THomingHead HomingHead,TEbomb Ebomb, const long double TStabSyst )
{
	*this =  TBomb();
	mHomingHead = HomingHead ;
	mEbomb = Ebomb ;
	mTStabSyst = TStabSyst;
}



// ������ ����� �������� ������
 long double TBomb::fncCx( long double valMach)
{
	if (valMach >= Ma�h01 )
	 return  Cx01;

	if (valMach >= Ma�h00 ) return  Cx00 + (Cx01 - Cx00)/(Ma�h01 - Ma�h00) * (valMach - Ma�h00);
	return  Cx00;
}


// ������ ����������� ����� �������� ������  �� ����� ����
 long double TBomb::fnc_dCx_po_dMach( long double valMach)
{
	if (valMach >= Ma�h01 ) return  0.;

	 return  2. *(Cx01 - Cx00)/(Ma�h01 - Ma�h00) ;

}

// ������ ����� ��������� ����
 long double TBomb::fncCy( long double valMach)
{
 /*  try
   {
	long double vala = 0.;
	if (valMach >= 1. )
	{
	  vala =  sqrt ( valMach * valMach - 1. ) ;
	}
	else
	{
	 vala = sqrt ( 1. -  valMach * valMach) ;
    }

	return 1.42 *((CyFuz * mSurfFuzSq/ mMidSq + CyKorma +  mSq_k/ mMidSq * mLamb_k * (Cy10 + (Cy11 - Cy10) * vala)));
   //	return 1. *((CyFuz * mSurfFuzSq/ mMidSq + CyKorma +  mSq_k/ mMidSq * mLamb_k * (Cy10 + (Cy11 - Cy10) * vala)));
	}
	catch(...)
	{
		int iii = 0;
	} */
	return 4.5;

}

// ������ ����������� ����� ��������� ����  �� ����� ����
 long double TBomb::fnc_dCy_po_dMach( long double valMach)
{
 /* try
  {
	if (valMach >= 1. ) return  0.;
	double vala = sqrt ( 1. -  valMach * valMach) ;

	return 1.42*( -mSq_k/ mMidSq * mLamb_k * (Cy11 - Cy10) *valMach /sqrt ( 1. -  valMach * valMach)) ;
   //return 1.*( -mSq_k/ mMidSq * mLamb_k * (Cy11 - Cy10) *valMach /sqrt ( 1. -  valMach * valMach)) ;
  }
  catch(...)
  {
	  int ii = 0;
  }  */
  return 0.;

}


//---------------------------------------------------------------------------
#pragma package(smart_init)
