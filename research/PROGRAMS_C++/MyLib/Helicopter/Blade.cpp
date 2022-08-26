//---------------------------------------------------------------------------


#pragma hdrstop

#include "Blade.h"


__fastcall  TBlade::TBlade()
{
  // ������ �����
	 mBladeR  = 0.;
	// ���������� �� ��� �������� ��  ������ ��
	 mRadHorizHsarnir =0.;
	// ������ ������������� �������  ������� ����� � ������ ��������
	 mPofile_d0 =0.;
	// ������ ������������� �������  ������� ������ �� ������ ��������
	 mPofile_d1=0.;
	// ����� �������
	 mBladeM=0.;
	// ����� �������
	 mBlade_b=0.;

}
// ����������� �����������
__fastcall  TBlade::TBlade (const TBlade &R)
 {
	// ������ �����
	 mBladeR  = R.mBladeR ;
	// ���������� �� ��� �������� ��  ������ ��
	 mRadHorizHsarnir = R.mRadHorizHsarnir;
	// ������ ������������� �������  ������� ����� � ������ ��������
	 mPofile_d0 = R.mPofile_d0;
	// ������ ������������� �������  ������� ������ �� ������ ��������
	 mPofile_d1 = R.mPofile_d1;
	// ����� �������
	 mBladeM = R.mBladeM;
	// ����� �������
	 mBlade_b = R.mBlade_b;
   }
 // �������� ������������
  TBlade TBlade::operator=(TBlade  R)
 {
	// ������ �����
	 mBladeR  = R.mBladeR ;
	// ���������� �� ��� �������� ��  ������ ��
	 mRadHorizHsarnir = R.mRadHorizHsarnir;
	// ������ ������������� �������  ������� ����� � ������ ��������
	 mPofile_d0 = R.mPofile_d0;
	// ������ ������������� �������  ������� ������ �� ������ ��������
	 mPofile_d1 = R.mPofile_d1;
	// ����� �������
	 mBladeM = R.mBladeM;
	// ����� �������
	 mBlade_b = R.mBlade_b;
	 return *this ;
 }

 // ����� ������ 1
 __fastcall TBlade::TBlade(const double  BladeR,const double  RadHorizHsarnir
   ,const double  Pofile_d0,const double  Pofile_d1,const double  BladeM  ,const double  Blade_b)
 {
   mBladeR = BladeR;
   mRadHorizHsarnir = RadHorizHsarnir;
   mPofile_d0 = Pofile_d0 ;
   mPofile_d1 = Pofile_d1;
   mBladeM = BladeM;
   mBlade_b = Blade_b;
 }
//-------------------------------------------------------------------------------
// ���������� ���������� �� ����� ���������� ���������������� ����
// ������������� �� ������ ������  ��� ��
double TBlade::calc_rQ()
{
   double x0 =  mRadHorizHsarnir;
   double x1 = mBladeR;
   double valK = (mPofile_d1  - mPofile_d0)/(x1 - x0);
   double temp =  mPofile_d0 - valK * x0;
   double temp2 = (x1 * x1 - x0 * x0) / 2.;
   double temp3 = (x1 * x1* x1 - x0 * x0* x0) / 3.;
   double temp4 = (x1 * x1* x1 * x1 - x0 * x0 * x0 * x0) / 4.;
   double val_rQ = (temp * temp3 + valK * temp4 )/ ( temp * temp2 + valK * temp3);
   return val_rQ;
}

//-------------------------------------------------------------
// ���������� ���������� �� ������ ���� �������  �� ��������� ����� ������� � ������ ��������
double TBlade::calc_X_BladeCentreMass()
{
	const double VAl_L = mBladeR - mRadHorizHsarnir;// VAl_L - ����� �������
	return VAl_L * (mPofile_d0 + 2. * mPofile_d1)/ (mPofile_d0 +  mPofile_d1) /3.;
}

//-------------------------------------------------------------
// ���������� ������������ ������� ������������� ������ ��������
double TBlade::calc_X_StatMoment_Sg()
{
	double valCentreMass = calc_X_BladeCentreMass();
	return   mBladeM * ( valCentreMass + mRadHorizHsarnir);
}

//-------------------------------------------------------------
// ���������� ������� �������  ������������ ������ �������
double TBlade::calcInertiaMoment0()
{
	const double VAlL = mBladeR - mRadHorizHsarnir;// VAl_L - ����� �������
	return   mBladeM * VAlL * VAlL *( mPofile_d0 * mPofile_d0 + mPofile_d0 * mPofile_d1 * 4.
	 + mPofile_d1* mPofile_d1)
	   / (mPofile_d0 + mPofile_d1) / (mPofile_d0 + mPofile_d1)/18.;
}

//-------------------------------------------------------------
// ���������� ������� �������  ������������ ������ ������ ��������������� �������
double TBlade::calcInertiaMomentHorSharnir()
{
	const double VAl_L = mBladeR - mRadHorizHsarnir;// VAl_L - ����� �������
	return   mBladeM* VAl_L * VAl_L *( mPofile_d0  +  3. *mPofile_d1)
	   / (mPofile_d0 + mPofile_d1) /6.;
}

//-------------------------------------------------------------
// ���������� ������� �������  ������������ ������ �������� �������
double  TBlade::calcInertiaMoment()
{
	double valXCentre = calc_X_BladeCentreMass() ;
	return   calcInertiaMoment0() + (valXCentre  +  mRadHorizHsarnir)
	 * (valXCentre  +  mRadHorizHsarnir) * mBladeM;
}


#pragma package(smart_init)