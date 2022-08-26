//---------------------------------------------------------------------------

#ifndef TargDataH
#define TargDataH
class TTargData
{
public:
	 double mLDanger ;// длина опасной части цели
	 double mLTarg  ; // длина цели

	// конструктор по умолчанию
	TTargData () ;
	// конструктор копирования
	TTargData  (const TTargData  &R) ;
	// оператор присваивания
	TTargData  operator=(TTargData   R2) ;

	// парам конструктор
	 TTargData (const double LDanger, const double LTarg)  ;

}  ;
#endif
