//---------------------------------------------------------------------------

#ifndef MeasureH
#define MeasureH
class TMeasure
{
public:
	// размерность вектора состояния объекта
	int mDimY;


//  1.Время привязки замера
	double mTYZv;
//  2.вектор замера
	double *mparrYZv;

 // 3.корреляционная матрица БМО
	double *mparrBMO_K;
	// 3.корреляционная матрица БМО
	double *mparrMMO_K;





   ~TMeasure() ;

	TMeasure () ;

	// конструктор копирования
	TMeasure  (const TMeasure  &R) ;
	// оператор присваивания
	TMeasure  operator=(TMeasure   R2) ;

	TMeasure(const int DimY, const double ValDispBMO );

	TMeasure(const int DimY, const double ValDispBMO , const double ValDispMMO )  ;


};
#endif
