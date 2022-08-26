//File filter.h
#ifndef filterH
#define filterH
extern int workFilter;
//Файл описания данных для алгоритма фильтрации
typedef struct Filter_Input
{
 double tk;//Время 
 double Xkg,Ykg,Hkg;//Первичные оценки координат цели в КГСК-ЦТ
 double sbet,seps,sD;//Дисперсии первичных измерений
 double Vcx,Vcy,Vch;//Вектор скорости своего корабля
}Filter_Input;//Структура входных данных в алгоритм фильтрации

typedef struct Filter_Output
{
 double X, Y, H;//Вторичные оценки координат цели в КГСК-ЦТ
 double Vx,Vy,Vh;//Вектор оценок скорости цели
}Filter_Output;//Структура выходных данных из алгоритма фильтрации

//Описание функции, реализующей алгоритм фильтрации
void filtr_ukf(Filter_Input *TFilter_Input, Filter_Output *TFilter_Output, bool bSliga, const double VAlSigW, const double VAlSigMMO);
void  fncExtrapolateCorMtrx (const double  p11h, const double p12h, const double p22h
		,const double  t,const double  VAlDispW, double *sp11h, double *sp12h, double *sp22h );

void  calcCoefUsilenia (const double  sp11h, const double sp12h
		,const double  VAlSigBMO, const double  VAlSigMMO, double *pUs0, double *pUs1, double *pvalP, double *pvalAlf);
void calcCorMtrxCompleted(const double VAlKExtr00, const double VAlKExtr01, const double VAlKExtr11
	,const double VAlP,const double VAlAlf,const double VAlP1
	, double* pvalK00, double* pvalK01, double* pvalK11);
#endif


 //End filter.h
