//---------------------------------------------------------------------------

#ifndef RectRectWaveGuideH
#define RectRectWaveGuideH
// открытый конец волновода прямоугольного сечения
class TRectWaveGuide
{
public:


 double ma ;  // большая сторона прямоугольника
 double mb ; // меньшая сторона прямоугольника
 double mFi1 ; //  уголменьший


   TRectWaveGuide();
 // конструктор копирования
 TRectWaveGuide (const TRectWaveGuide &R) ;
 // оператор присваивания
 TRectWaveGuide operator=(TRectWaveGuide  R);
 // параметр констр
  TRectWaveGuide(  const double a,  const double b );

 // void ShowMe(wchar_t *FileName);
 double fncDiagr(const double VAlLambda, const double VAlTetta);

};
#endif
