//---------------------------------------------------------------------------


#pragma hdrstop

#include "Glonass.h"

TGlonass ::TGlonass()
{
	//  ��� �� ���������
 mSKOPos = 3. ;
 // ��������� ���������
 mSKOVeloc = 0.4;

}

 // ����� ������
TGlonass :: TGlonass(long const double SKOPos,long const double SKOVeloc )

 {
   mSKOPos = SKOPos;
   mSKOVeloc = SKOVeloc;

 }

 // �������� ������������
 TGlonass TGlonass::operator=(TGlonass  R)
 {
   mSKOPos = R.mSKOPos ;
   mSKOVeloc = R.mSKOVeloc;
	return *this ;
 }

 // ����������� �����������
 TGlonass::TGlonass (const TGlonass &R)
 {
   mSKOPos = R.mSKOPos ;
   mSKOVeloc = R.mSKOVeloc;
 }

#pragma package(smart_init)
