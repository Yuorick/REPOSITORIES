//---------------------------------------------------------------------------


#pragma hdrstop

#include "TargData.h"

TTargData::TTargData()
{
	  mLDanger = 30.;// ����� ������� ����� ����
	  mLTarg = 200.;

}

//---------------------------------------------------------------------------


// ����������� �����������
 TTargData ::TTargData (const TTargData &R)
 {
	mLDanger  =  R.mLDanger;
	mLTarg = R.mLTarg;

 }
 // �������� ������������
 TTargData TTargData::operator=(TTargData  R)
 {
	mLDanger  =  R.mLDanger;
	mLTarg = R.mLTarg;


	return *this ;
 }

  // ����� �����������
 TTargData::TTargData (const double LDanger, const double LTarg ) {
	mLDanger = LDanger;
	mLTarg = LTarg;

}

#pragma package(smart_init)
