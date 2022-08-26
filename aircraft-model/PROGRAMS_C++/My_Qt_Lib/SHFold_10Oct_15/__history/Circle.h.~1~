//---------------------------------------------------------------------------

#ifndef CircleH
#define CircleH
#include "UrPointXY.h"
class TURPointXY ;
class TCircle
{
public:

 TURPointXY mPntCentre; // центр
 double mR ;  //  радиус


   TCircle();
 // конструктор копирования
 TCircle (const TCircle &R) ;
 // оператор присваивания
 TCircle operator=(TCircle  R);
 // параметр констр
  TCircle(const TURPointXY PntCentre, const double R );
	// парам констр 2

  void ShowMe(wchar_t *FileName);
  int fncParametricLineCutCircle( TURPointXY pointLine, double *arrVectLine // линия
							,TURPointXY *arrPntRez // точек ассив а пересечения
								);
  int fncParametricVerticalLineCutCircle( TURPointXY pointLine, double *arrVectLine // линия
							,TURPointXY *arrPntRez // точек ассив а пересечения
							) ;



};
#endif
