//---------------------------------------------------------------------------

#ifndef Table_2DH
#define Table_2DH

// класс описывает 2-х мерную таблицу
// таблица представляет из себя массив (вектор) одномерных таблиц в точках
// значений аргумента mparrArg
// значения функции хранятся в массиве  mparrVal
// значения аргументов в массиве mparrArg
// массив mparrArg представляет из себя монотонную последовательность аргументьов
// mNumCols - число точек в котороых задана функция
class TTable_1D;
//
class TTable_2D
{
public:

	int mNumCols ;//
	double *mparrArg; // вектор значений столбцов
	TTable_1D *mpArrTable_1D; // вектор одномерных таблиц-столбцов
	double Box[2];


	__fastcall ~TTable_2D() ;

	TTable_2D();
 // конструктор копирования
 TTable_2D (const TTable_2D &R) ;
 // оператор присваивания
 TTable_2D operator=(TTable_2D  R);
// парам констр
TTable_2D( double *parrArg, TTable_1D *pArrTable_1D, const int NumCols);
TTable_2D( double *parrArgTab1, const int NumColsTab1, double *parrArgTab2, const int NumColsTab2, double *parrVal );

double calcValue(const double VAlRowArg, const double VAlColArg);

void calcBoundBox();

};
#endif
