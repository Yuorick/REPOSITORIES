//---------------------------------------------------------------------------


#pragma hdrstop

#include "TransmitAnt.h"


__fastcall  TTransmitAnt::TTransmitAnt()
{
  // �������� �� ��������
	 mPowerPrd = 4000.;
	// �� �� ��������
	 mKYPrd = 3000.;
}
// ����������� �����������
__fastcall  TTransmitAnt::TTransmitAnt (const TTransmitAnt &R2)
 {
	 mPowerPrd = R2.mPowerPrd ;
	 mKYPrd  = R2.mKYPrd;
 }
 // �������� ������������
  TTransmitAnt &TTransmitAnt::operator=(const TTransmitAnt  &R2)
 {
	 mPowerPrd = R2.mPowerPrd ;
	 mKYPrd  = R2.mKYPrd;
	 return *this ;
 }

 // ����� ������ 1
 __fastcall TTransmitAnt::TTransmitAnt(const double PowerPrd,const double KYPrd)
 {
	 mPowerPrd = PowerPrd;
	 mKYPrd = KYPrd;
 }

#pragma package(smart_init)
