//---------------------------------------------------------------------------

#ifndef TrajectH
#define TrajectH
#include "InitTargData.h"

//#define LEN_TARG_VS  9
class TInitTargData;
class TTraject
{
public:
	double mTCur; // текущее время
	double mTBegin ; // время начала траектории

	double marrSigW[3]; // скз шума объекта
	double marrVectSostGSK_Begin[9]; // вектор состояния в начальный момент - нач условия
	double marrVectSostGSK[9];     // вектор состояния на момент mTCur
	double marrVectDeviations[9];   // вектор возмущений на мометн mTCur
 //	THomingHead3D mHomingHead3D ;


	// конструктор по умолчанию
	TTraject () ;
	// конструктор копирования
	TTraject  (const TTraject  &R) ;
	// оператор присваивания
	TTraject  &operator=(const TTraject   &R2) ;

	// парам конструктор
	TTraject (const double TCur,  const double 	SigW
		,double *arrVectDeviations,double *arrVectSostGSK,double *arrVectSostGSK_Begin);

	TTraject (const double TCur,  double 	*arrSigW,	TInitTargData InitData );

	//
   TTraject (const double TCur,  const double 	SigW,TInitTargData InitData );

   TTraject (const double TBegin,  const double SigW ,double *arrVectSostGSK_Begin );

	int Sign(const double a);



	void recalcTrajPoint(const double tNext) ;

   void extrapolateTargVS(const double VAlTExtr, double *arrTargExtrapVS);

   void extrapolateTarg_BeginVS(const double VAlTExtr, double *arrTargExtrapVS);




}  ;
#endif
