//---------------------------------------------------------------------------

#ifndef Table_1DH
#define Table_1DH
// класс описывает однмерную таблицу
// таблица представляет из себя массив (вектор) значений некоторой функции в точках
// значения функции хранятся в массиве  mparrVal
// значения аргументов в массиве mparrArg
// mNumCols - число ьточек в котороых задана функция

//
class TTable_1D
{
public:

	int mNumCols ;// к-во столбцов
	double *mparrArg; // вектор значений столбцов
	double *mparrVal; // вектор аргументов столбцов
    double Box[4];



 ~TTable_1D() ;

 TTable_1D();
 // конструктор копирования
 TTable_1D (const TTable_1D &R) ;
 // оператор присваивания
 TTable_1D &operator=(const TTable_1D  R);
// парам констр1
TTable_1D( double *parrArg, double *parrVal, const int NumCols);
// парам констр2
TTable_1D( double *parrBuff, const int NumCols);

double calcValue(const double VAlArg);

void calcBoundBox();

double calcLinearValueApprox(const double x);

int getSegmentNum(const double x);

double calc_d_po_dx(const double x);

void multValue(const double coeff);

void multRandValue(const double sig);

TTable_1D makeMiddleProfile();

double calcMiddleValue();

void addValue(const double val);



};
#endif
