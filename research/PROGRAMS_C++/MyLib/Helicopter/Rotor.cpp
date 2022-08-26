//---------------------------------------------------------------------------


#pragma hdrstop

#include "Rotor.h"


__fastcall  TRotor::TRotor()
{
// лопасть винта
	 mBlade = TBlade();
	// к-во лопастей
	 mQuantBlades = 0;
	 //
	 mBasePLane  = TPlane ();
}
// Конструктор копирования
__fastcall  TRotor::TRotor (const TRotor &R)
 {
	 mBlade = R.mBlade;
	 mQuantBlades   = R.mQuantBlades ;
	 mBasePLane  = R.mBasePLane ;

 }
 // оператор присваивания
  TRotor TRotor::operator=(TRotor  R)
 {
	 mBlade = R.mBlade;
	 mQuantBlades   = R.mQuantBlades ;
	 mBasePLane  = R.mBasePLane ;
	 return *this ;
 }

 // парам констр 1
 __fastcall TRotor::TRotor(const TBlade Blade, const double QuantBlades, const TPlane BasePLane)
 {
	 mBlade = Blade;
	 mQuantBlades = QuantBlades;
	 mBasePLane = BasePLane ;
 }
 //-------------------------------------------------------
#pragma package(smart_init)
