//---------------------------------------------------------------------------


#pragma hdrstop

#include "Glonass.h"

TGlonass ::TGlonass()
{
	//  СКО по положению
 mSKOPos = 3. ;
 // дальность диаграммы
 mSKOVeloc = 0.4;

}

 // парам констр
TGlonass :: TGlonass(long const double SKOPos,long const double SKOVeloc )

 {
   mSKOPos = SKOPos;
   mSKOVeloc = SKOVeloc;

 }

 // оператор присваивания
 TGlonass TGlonass::operator=(TGlonass  R)
 {
   mSKOPos = R.mSKOPos ;
   mSKOVeloc = R.mSKOVeloc;
	return *this ;
 }

 // конструктор копирования
 TGlonass::TGlonass (const TGlonass &R)
 {
   mSKOPos = R.mSKOPos ;
   mSKOVeloc = R.mSKOVeloc;
 }

#pragma package(smart_init)
