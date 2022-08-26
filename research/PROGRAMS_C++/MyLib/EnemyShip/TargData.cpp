//---------------------------------------------------------------------------


#pragma hdrstop

#include "TargData.h"

TTargData::TTargData()
{
	  mLDanger = 30.;// длина опасной части цели
	  mLTarg = 200.;

}

//---------------------------------------------------------------------------


// конструктор копирования
 TTargData ::TTargData (const TTargData &R)
 {
	mLDanger  =  R.mLDanger;
	mLTarg = R.mLTarg;

 }
 // оператор присваивания
 TTargData TTargData::operator=(TTargData  R)
 {
	mLDanger  =  R.mLDanger;
	mLTarg = R.mLTarg;


	return *this ;
 }

  // парам конструктор
 TTargData::TTargData (const double LDanger, const double LTarg ) {
	mLDanger = LDanger;
	mLTarg = LTarg;

}

#pragma package(smart_init)
