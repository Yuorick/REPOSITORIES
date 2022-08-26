//---------------------------------------------------------------------------

#ifndef GlonassH
#define GlonassH
class TGlonass
{
public:
 //  CКО ПО ПОЛОЖЕНИЮ
 long double mSKOPos;
 // CКО ПО Cкорости
 long double mSKOVeloc;




   TGlonass();
 // конструктор копирования
 TGlonass (const TGlonass &R) ;
 // оператор присваивания
 TGlonass operator=(TGlonass  R);
 // параметр констр
 TGlonass(long const double SKOPos,long const double SKOVeloc );

 };
#endif
