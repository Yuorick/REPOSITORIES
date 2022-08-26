//---------------------------------------------------------------------------

#ifndef BombH
#define BombH
//---------------------------------------------------------------------------

#include "HomingHead.h"
#include "Ebomb.h"

class THomingHead;
class TEbomb ;


// Класс описывает характеристики УАБ
// УАБ - тело вращения с расширяющейся кормовой частью (Лебедев Ченрнобровкин, стр 152, п 9)
//  и прямоунгольными крыльями
class TBomb
{
public:
 // 1.площадь иделя фюзеляжа
	long double mMidSq;
// 2.площадь боковой поверзности фюзеляжа
	long double mSurfFuzSq ;

// 3.масса
	long double mMass ;
// 4.макс угол атаки:
	long double mMaxAt ;
 // 5.уцдлинение крыльев
	long double mLamb_k ;
 // 6.площадб крыльев, м кв
	long double mSq_k ;
 // 7. ГСН
	THomingHead mHomingHead;
 // 8.ЭМБ
	TEbomb mEbomb;
 // 9. пост времени системы стабилизации по углу атаки
	long double mTStabSyst;


	// конструктор по умолчанию
	TBomb () ;
	// конструктор копирования
	TBomb  (const TBomb  &R) ;

	// оператор присваивания
	TBomb  operator=(TBomb   R2) ;
   // параметрическийц конструктор
	TBomb(
	 long double MidSq // площадь иделя фюзеляжа
	,long double SurfFuzSq // площадь боковой поверзности фюзеляжа
	,long double Mass// масса
	,long double MaxAt // макс угол атаки
	,long double Lamb_k // уцдлинение крыльев
	,long double Sq_k // площадб крыльев, м кв
	,bool bHomeHead // признак наличия ГСН
	,THomingHead HomingHead // ГСН
	,TEbomb Ebomb// ЭМБ
	);
	// параметрическийц конструктор 2
   TBomb(THomingHead HomingHead );
	// параметрическийц конструктор 3
   TBomb(THomingHead HomingHead,TEbomb Ebomb ) ;

   // параметрическийц конструктор 4
   TBomb(THomingHead HomingHead,TEbomb Ebomb, const long double TStabSyst );



	long double fncCx( long double valMach);
	long double fncCy( long double valMach);
	long double fnc_dCx_po_dMach( long double valMach);
	long double fnc_dCy_po_dMach( long double valMach);

}  ;
#endif
