//---------------------------------------------------------------------------

#ifndef LineMoveH
#define LineMoveH
class TLineMove
{
public:
// Врея привязки информации
	double mT0;
//  вектор состояния в момент mT
	double marrVS[9];


	 __fastcall ~TLineMove() ;
	// конструктор по умолчанию
	TLineMove () ;
	// конструктор копирования
	TLineMove  (const TLineMove  &R) ;
	// оператор присваивания
	TLineMove  operator=(TLineMove   R2) ;

	// парам конструктор1
	TLineMove (const double T,  double *arrVS) ;

	void ExtrapolateVS (const double tExtr, double *arrVSTargExtr_GSK);



}  ;
#endif
