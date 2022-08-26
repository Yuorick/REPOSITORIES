//---------------------------------------------------------------------------

#ifndef ManulEffectImproved_v8H
#define ManulEffectImproved_v8H
//---------------------------------------------------------------------------

#include "TargBearing0.h"
#include "EtalonSign.h"
#include "Far_2D.h"
#include "Traject.h"
#include "Target.h"
//Файл данных для расчета эффективности
//File MODEL_RASSCHET_EFF.h

 //******************************************************************************
 // Изменение № 10
 //******************************************************************************
 #define k 9 //Число интегрируемых уравнений +1
//******************************************************************************
// Конец изменения № 10
//******************************************************************************
   #include <math.h>


 class TFar_2D ;
 class TEtalonSign;
 class TTargBearing0;
//enum enumTargetType ;//{NOMOVING, ABOVEWATER, PLANE};//Типы целей


 //int //Kr,//Количество реализаций
  //	 Kcn;//Количество снарядов

//===============================================================================
//Константы для получения случайных чисел


//Функция расчета координат и скоростей своего корабля
 void raschet_coord_swoego_corablja (double ti, const double VAlSigSins);
  void fnc_Deriv_f(const double VAlCalibro, double *y,double *yp);
 //  void raschet_gauss_raspredelenie ();
 void fncPreviousArrangments(TTarget Targ0
   , const double VAlCalibro, const double VAlDn,const double VAlDk
	, const int n, const double dti,
	double &t00, double &tp00,  double &tk, double &tpkk);

 	void ResetFun();

 //Функция расчета координат и скоростей цели
 void raschet_coord_zeli (enumTargetType TargType,  const double WTargSkz, TTargBearing0  TargBear0, const double ti,  const double dti
	,double &targX, double &targY, double &targH
	,double &targVX, double &targVY, double &targVH);
 // расчет начальных координат цели
void raschet_nach_coord_zeli (enumTargetType TargType,TTargBearing0  TargBear0
	,double &targX0, double &targY0, double &targH0);

void 	raschet_nach_velo_zeli (TTargBearing0 TargBear0, double &valTargGSK_VX0
	,  double &valTargGSK_VY0,  double &valTargGSK_VH0);
void raschet_coord_zeli_GSK_ideal (enumTargetType TargType,  TTargBearing0  TargBear0, const double VAlT
	,double &valTargGSK_X, double &valTargGSK_Y, double &valTargGSK_H
	,double &valTargGSK_VX, double &valTargGSK_VY, double &valTargGSK_VH);

double fncLinExtrapolation(const double VAlF0, const double VAlF1
  ,const double VAlX,const double VAlDelX);
 //===============================================================================
 //===============================================================================
 int calcEffect(wchar_t *wchOutFold, TFar_2D Far_2D,  const double HAntenna, TTarget Targ , TInitTargData InitTargData
	 ,const double PowerPrd, const double KYPrd, TEtalonSign EtalonSign, const double VAlDT_filtr
		 , const int QuantIspit, const int QuantShells,const double VAlCalibro
	, const double VAlDn0,  const double VAlDk, const double VAlSigSins, const bool VAlBSkaliga
		, double *pvalProb, double *pvalProb_Gladk, double *pvalDistBeginSopr, int *pQuantShots);

 double MAX__(double a , double b)
 {
 return (a > b)? a:b;
 }

 double MIN__(double a , double b)
 {
 return (a < b)? a:b;
 }

 int MIN__(int a , int b)
 {
 return (a < b)? a:b;
 }
 //Функция расчета числа, равномерно распределенного на [0,2pi]
// void raschet_randu_raspredelenie ();
 //Функция расчета наблюдаемых координат цели
 void raschet_coord_zeli_nabl(double *arrVSTargKGSK_True, double *arrVSTargKGSK_Zv );

typedef struct INTEGRIR_Input
{
	double tp0,tp;//Начало и конец интегрирования
	double h0,b;//Шаг интегрирования со знаком b
//******************************************************************************
 // Изменение № 11
 //******************************************************************************
  	double y[9];//Начальные условия на систему интегрируемых уравнений
 //******************************************************************************
 // Конец изменения № 11
 //******************************************************************************
}INTEGRIR_Input;//Структура входных данных в блок интегрирования




 typedef struct RZW_Input
 {
	 double Xzel,Yzel,Hzel;//Координаты цели
	 double Vxzel,Vyzel,Vhzel;//Составляющие скорости цели
	 double fi,q;//Начальные углы наведения снаряда
 }RZW_Input;//Структура входных данных в RZW

 typedef struct RZW_Output
 {
	 double fi,q;//Углы наведения снаряда
	 double tp;//Полетное время
//******************************************************************************
 // Изменение № 16
 //******************************************************************************
	 double y[9];//Координаты (y[1],y[2],y[3]), x,h,y
				 //составляющие скорости (y[4],y[5],y[6]),
	             //плотность (y[7]) воздуха
 //******************************************************************************
 // Конец изменения № 16
 //******************************************************************************
	 double xy,yy,hy;//Упрежденные координаты цели
	 double dy;   // дальность до цели в точке встречи в КГСК
 }RZW_Output;//Структура выходных данных из RZW и любого блока из RZW



//--------------------------------------------------------------


 double  RZW__(const double VAlCalibro,double *arrVS_Vessel_GSK00,  double *arrVSTarg0, RZW_Output *TRZW_Output,
	double& Fi0, // начальный угол места точки бросания
	double& Qu0, // начальный азимут точки бросания
	double kMuldFi_2, // множитель приращения угла места для 2-го этапа алгоритма
	double kMuldQu_2 // множитель приращения азимута для 2-го этапа алгоритма
	 //,	double Vb0 // нач.скорость снаряда
		);

double AimFun_RZW__(	const double VAlCalibro
										, double *arrVS_Vessel_GSK00
													, double *arrVS_KGSK_Targ0
													,double Fi // угол места
													,double Qu // азимут
													,double *arrVS_Shell_KGSK// вектор состояния снаоряда в точке встречи
													,double *arrVS_Shell_SSK// вектор состояния снаоряда в точке встречи
													, double *arrVS_Targ_KGSK// вектор состояния цели  в точке встречи
													,double* tp_manul // подлётное время
													,double val_dtInt
												);


 void From_xyhKGSK_To_xyhSSK(double cosQu, double sinQu, double x,
	double y, double h, double& x1, double& y1, double& h1);

 void From_xyhSSK_To_xyhKGSK(double cosQu, double sinQu, double x,
	double y, double h, double& x1, double& y1, double& h1);

void InitFi0Qu0(double *arrVSTarg_KGSK0, double& Fi0, double& Qu0, double Vb0);





 void From_xyhKGSK_To_xyhSSK(double cosQu, double sinQu, double *arrKGSKInp, double *arrOut);
 void From_xyhSSK_To_xyhKGSK(double cosQu, double sinQu, double *arrSSKInp, double *arrOut) ;
int SZn,SZnn,SZnn1;



double fncSignum(double x);

void fncMove_130_Cal_Shell_TO_ZeroAlt_AND_ShowGraphs(double *arrStrSK_VS
   ,const double VAlCoefCx, const double VAlCoefCy, const double VAlCoefCz
   ,const double VAlVessV, const double VAlVessQ, const  double StepInt
		, wchar_t *wcharrPath,  double &valDHoriz );

void calcEilerStep(double *arrStrSK_VS/*,const double VAlVessV, const double VAlVessQ*/, const  double valStepInt);



/*
//Тип А1, новая, 300 м/c
 double Kr_Glad_1[5]={0,0.24,0.47,0.73,1};
 double Mas_Dal_1[5]={429.4,576.47,629.4,652.9,703};
*/
 /*
 //Тип B, новая,300м/c
 double Kr_Glad_1[13]={0,0.14,0.33,0.66,0.88,0.91,0.93,0.94,0.95,0.96,0.975,0.99,1};
 double Mas_Dal_1[13]={1980,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000};
 */
/*
 //Тип А1, новая, 700 м/c
 double Kr_Glad_1[5]={0,0.25,0.5,0.743,1};
 double Mas_Dal_1[5]={700.0,900.0,1000.0,1100.0,1200.0};
*/
/*
 //Тип А1, новая, 1000 м/c
 double Kr_Glad_1[5]={0,0.25,0.5,0.743,1};
 double Mas_Dal_1[5]={490.0,630.0,700.0,770.0,840.0};
*/

/*
 //Тип B, новая,  700м/c
 double Kr_Glad_2[6]={0,0.108,0.215,0.554,0.877,1};
 double Mas_Dal_2[6]={4500.0,5000.0,5500.0,6000.0,6500.0,7000.0};
*/
/*
 //Новый условный закон поражения, 100 мм, тип А1+B
 double X_yzp[9]={0.0,0.5,0.75,1.0,1.5,2.0,2.5,3.0,3.5};
 double Y_yzp[9]={1.0,1.0,0.96,0.904,0.704,0.496,0.296,0.144,0.064};
*/
/*
 //Новый условный закон поражения, 100 мм, тип А1
 double X_yzp[6]={0.0,0.5,1.0,1.5,2.0,2.5};
 double Y_yzp[6]={1.0,1.0,0.85,0.45,0.1,0.0};
*/
/*
 //Новый условный закон поражения, 100 мм, тип В
 double X_yzp[8]={0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5};
 double Y_yzp[8]={1.0,0.6,0.3,0.45,0.43,0.3,0.15,0.06};
*/



/*
 //Тип B, старая,1000 м/c
 double Kr_Glad_1[16]={0,0.04,0.1,0.14,0.25,0.26,0.44,0.52,0.62,0.68,0.78,0.8,0.844,0.94,0.96,1};
 double Mas_Dal_1[16]={0,210,322,420,560,700,875,1050,1225,1400,1750,1862,2100,2800,3500,4200};
*/

double calcHittingProbabylity(TTarget Targ , const double VAlCalibro, double valProm);
//End MODEL_RASSCHET_EFF.h
#endif
