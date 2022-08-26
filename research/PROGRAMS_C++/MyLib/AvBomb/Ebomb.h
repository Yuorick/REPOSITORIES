//---------------------------------------------------------------------------

#ifndef EbombH
#define EbombH
class TEbomb
{
public:
 //  угол конуса поражения
 long double mFi;
 // дальность конуса поражения
 long double mR ;
 // поток мощности

 long double mW;



   TEbomb();
 // конструктор копирования
 TEbomb (const TEbomb &R) ;
 // оператор присваивания
 TEbomb operator=(TEbomb  R);
 // параметр констр
  TEbomb(const long  double Fi,const long  double R,const long  double W );
 long double TEbomb::fncCalcPower ();

static void createGraphR_of_Fi (wchar_t *pwchOutFile, double valW, double valP);

 };
#endif
