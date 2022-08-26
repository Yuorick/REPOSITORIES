// 3 - Vx//---------------------------------------------------------------------------

#ifndef BladeH
#define BladeH
// класс описывает лопасть винта вертолета
class TBlade
{
public:

	// радиус ометаемой площади
    long double  mBladeR ;
	// расстояние от оси вращения до  втулки ГШ
    long double  mRadHorizHsarnir;
	// высота вертикального сечения  лопасти ближе к центру вращения
    long double  mPofile_d0;
	// высота вертикального сечения  лопасти дальше от центра вращения
    long double  mPofile_d1;
	// масса лопасти
    long double  mBladeM;

	// хорда лопасти
    long double  mBlade_b;

    // коэффициент подъемной силы
    long double  mCyalfa;

    // коэффициент лобового сопротивления
    long double  mCX0;

    // статический момент
    long double  mSg;

    // момент инерции
    long double  mIg;



      TBlade() ;
	// Конструктор копирования
      TBlade (const TBlade &R2) ;

	// оператор присваивания
	TBlade   operator=(TBlade  R2) ;

	// парам констр
     TBlade(const long double   BladeR,const long double   RadHorizHsarnir
            ,const long double   Pofile_d0,const long double   Pofile_d1,const long double   BladeM
                         ,const long double   Blade_b, const long double  CX0);

   long double  calc_rQ();

   long double  calc_X_BladeCentreMass()  ;

   long double  calc_X_StatMoment_Sg() ;

   long double   calcInertiaMoment0();

   long double  calcInertiaMomentHorSharnir() ;

   long double   calcInertiaMoment() ;
};
#endif
