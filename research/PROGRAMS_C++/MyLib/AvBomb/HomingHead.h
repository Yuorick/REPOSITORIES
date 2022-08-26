//---------------------------------------------------------------------------

#ifndef HomingHeadH
#define HomingHeadH
// класс описываеть самонаводящийся носитель на последне участке траектории
class THomingHead
{
public:
 //  угол диаграммы
 long double mFi;
 // дальность диаграммы
 long double mR ;
 // сигма по дальности, м
 long double mSigFl_R ;
 // сигма по углу, мрад
 long double mSigFl_U ;
 // интервал врекмени между замерами
 long double m_h ;

 THomingHead();
 // конструктор копирования
 THomingHead (const THomingHead &R) ;
 // оператор присваивания
 THomingHead operator=(THomingHead  R);
 // параметр констр
 THomingHead( const long  double Fi, const long  double R
   ,  const long  double SigFl_R, long  double SigFl_U , long  double h   );
 };
#endif
