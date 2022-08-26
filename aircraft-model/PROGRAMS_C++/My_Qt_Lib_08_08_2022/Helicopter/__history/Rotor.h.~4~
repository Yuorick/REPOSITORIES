//---------------------------------------------------------------------------

#ifndef RotorH
#define RotorH
#include "Plane.h"
#include "Blade.h"

// класс описывает винт вертолета
class TBlade;
class TPlane;
class TRotor
{
public:
	// лопасть винта
	TBlade mBlade;
	// к-во лопастей
	int mQuantBlades;
	// базовачя система координат привязанная к оси вращения вала  винта
	// ось Z направлена по оси вращения вала
	 TPlane mBasePLane;




 __fastcall  TRotor() ;
// Конструктор копирования
__fastcall  TRotor (const TRotor &R2) ;

 // оператор присваивания
 TRotor   operator=(TRotor  R2) ;

  // парам констр
   __fastcall TRotor::TRotor(const TBlade Blade, const double QuantBlades, const TPlane BasePLane);

};
#endif
