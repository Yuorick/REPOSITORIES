//---------------------------------------------------------------------------

#ifndef SectorH
#define SectorH
#include "UrPointXY.h"
class TURPointXY ;
class TSector
{
public:

 TURPointXY mPntCentre; // центр
 double mR ;  //  радиус
 double mFi0 ; // больший угол  mFi0 > mFi1 , направление по часовой стрелке
 double mFi1 ; //  уголменьший


   TSector();
 // конструктор копирования
 TSector (const TSector &R) ;
 // оператор присваивания
 TSector &operator=(const TSector  &R);
 // параметр констр
  TSector(const TURPointXY PntCentre, const double R,  const double Fi0,  const double Fi1 );
 	// парам констр 2

  void ShowMe(wchar_t *FileName);


};

#endif
