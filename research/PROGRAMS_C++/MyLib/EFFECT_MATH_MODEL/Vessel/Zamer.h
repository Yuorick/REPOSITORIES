//---------------------------------------------------------------------------

#ifndef ZamerH
#define ZamerH
// одиночный сигнал
class TZamer
{
public:

	 // массив измерений
	double marrMeas[3];
	// корреляционная матрица ошибок
	double marrCorr[9];
	// время замера
	double mT;





 __fastcall  TZamer() ;
// Конструктор копирования
__fastcall  TZamer (const TZamer &R2) ;
 // парам констр
 __fastcall TZamer ( double *arrMeas,double *arrCorr, const double VAlT) ;

 // оператор присваивания
 TZamer   &operator=(const TZamer  &R2) ;




};
#endif
