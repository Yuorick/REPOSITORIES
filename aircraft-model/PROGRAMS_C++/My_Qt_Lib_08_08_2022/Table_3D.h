//---------------------------------------------------------------------------

#ifndef Table_3DH
#define Table_3DH
// класс описывает 2-х мерную таблицу
// таблица представляет из себя массив (вектор) 2-х мерных таблиц в точках
// значений аргумента mparrArg
// значения функции хранятся в массиве  mparrVal
// значения аргументов в массиве mparrArg
// массив mparrArg представляет из себя монотонную последовательность аргументьов
// mNumCols - число точек в котороых задана функция
class TTable_2D;
//
class TTable_3D
{
public:

	int mNumCols ;//
	double *mparrArg; // вектор значений столбцов
	TTable_2D *mpArrTable_2D; // вектор одномерных таблиц-столбцов
	double Box[2];


	__fastcall ~TTable_3D() ;

	TTable_3D();
 // конструктор копирования
 TTable_3D (const TTable_3D &R) ;
 // оператор присваивания
 TTable_3D operator=(TTable_3D  R);
// парам констр
TTable_3D( double *parrArg, TTable_2D *pArrTable_2D, const int NumCols);
//
TTable_3D( double *parrArgTab1, const int NumColsTab1, double *parrArgTab2, const int NumColsTab2
 ,double *parrArgTab3, const int NumColsTab3, double *parrVal );

double calcValue(const double VAlTab3Arg,const double VAlRowArg, const double VAlColArg);

void calcBoundBox();

};
#endif
