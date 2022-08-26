//---------------------------------------------------------------------------


#pragma hdrstop
#include <vcl.h>
#include "StabSyst2.h"
#include  <string.h>
#include "Equations.h"
#include "math.h"
#include "MatrixProccess.h"
#include <float.h>
#include "URPolyLine.h"
#include "UrPointXY.h"



extern double EPS = 0.0001;
extern double PI = 3.1415926;
const int INET1 = 100;
const int INET2 = 100;

__fastcall TStabSyst2::~TStabSyst2()
{

}
//---------------------------------------------------------------------------


  TStabSyst2 ::TStabSyst2()
{
	mh = 0.05;
	msigm_w = 30 ;   //  скз w
	 msigm_mmo = 10 ; //  скз ммо
	  msigm_bmo = 10 ;
	double arrT0[4] = {1,0,0,1};
	arrT0[1] = mh ;
	memcpy(marrA, arrT0, 4 * sizeof(double)) ;    // матрица перехода за один такт;
	 marrB [0] = mh * mh/2 ;    // а атрицм помехи в канале объекта
	 marrB [1] = mh ;
	 marrC[0] = 1 ;
	 marrC[1] = 0;
	 msigm_fl = 0.;// флукт добавкка

}

// конструктор копирования
 TStabSyst2 ::TStabSyst2 (const TStabSyst2 &R)
 {
   memcpy(marrA,R.marrA,4 * sizeof(double));
   memcpy(marrB,R.marrB,2 * sizeof(double));
   memcpy(marrC,R.marrC,2 * sizeof(double));
   mh = R.mh ;
   msigm_w = R.msigm_w;
   msigm_mmo = R.msigm_mmo;
   msigm_bmo = R.msigm_bmo ;
   msigm_fl = R.msigm_fl;
 }
 // оператор присваивания
 TStabSyst2 TStabSyst2::operator=(TStabSyst2  R)
 {
   memcpy(marrA,R.marrA,4 * sizeof(double));
   memcpy(marrB,R.marrB,2 * sizeof(double));
   memcpy(marrC,R.marrC,2 * sizeof(double));
   mh = R.mh ;
   msigm_w = R.msigm_w;
   msigm_mmo = R.msigm_mmo;
   msigm_bmo = R.msigm_bmo ;
   msigm_fl = R.msigm_fl;
   //
   return *this ;
 }
// парам констр
TStabSyst2 :: TStabSyst2(  const double h,const double sigm_w, const double sigm_mmo
	  ,  const double sigm_bmo )

 {
   mh = h;
   msigm_w = sigm_w;
   msigm_mmo = sigm_mmo;
   msigm_bmo = sigm_bmo ;
  double arrT0[4] = {1,0,0,1};
	arrT0[1] = mh ;
	memcpy(marrA, arrT0, 4 * sizeof(double)) ;    // матрица перехода за один такт;
	 marrB [0] = mh * mh/2 ;    // а атрицм помехи в канале объекта
	 marrB [1] = mh ;
	 marrC[0] = 1 ;
	 marrC[1] = 0;
 }
 // парам констр
TStabSyst2 :: TStabSyst2(  const double h,const double sigm_w, const double sigm_mmo
	  ,  const double sigm_bmo, const double sigm_fl)

 {
   mh = h;
   msigm_w = sigm_w;
   msigm_mmo = sigm_mmo;
   msigm_bmo = sigm_bmo ;
   msigm_fl = sigm_fl ;
  double arrT0[4] = {1,0,0,1};
	arrT0[1] = mh ;
	memcpy(marrA, arrT0, 4 * sizeof(double)) ;    // матрица перехода за один такт;
	 marrB [0] = mh * mh/2 ;    // а атрицм помехи в канале объекта
	 marrB [1] = mh ;
	 marrC[0] = 1 ;
	 marrC[1] = 0;
 }
// аналитический расчет фильтра Калмана в установиви режиме
// Input:
// sigW, sigB - скз шума обекта и шума измерителя
// Output:
// arrK - коррелл матрица
// arrP - коефф усиления
// rt1, rt2 - корни характерист уравнения
// valT1, valT2 - память фильтра по положению и по скорости
// если система неустойчивая, то  valT1 = valT2 =  -1
void TStabSyst2 ::StabSolutionKalm( double *arrK,double * arrP,TComp &rt1,TComp &rt2,
	 double &valT1, double &valT2)
{
	double c =  msigm_w * mh * mh / msigm_bmo ;
	double x = SolvCharactEqKlm(c);
	arrK[0] =  msigm_bmo * msigm_bmo * (x * x -1 )/ x / x;
	arrK[1] =arrK[2] =  msigm_bmo * msigm_w * mh / x ;
	arrK[3] = msigm_w * msigm_bmo * ( x * x -1 ) /x - msigm_w * msigm_w * mh * mh /2 ;
	arrP [0] =  (x * x -1 )/ x / x;
	arrP [1] =  msigm_w  * mh/ msigm_bmo / x ;
	// расчет корней характ уравнения
	double a1 =1;
	double b1 = 2 - arrP [0] - mh * arrP [1] ;
	double c1 = 1- arrP [0] ;
	SolvEq2( a1, b1, c1,rt1,rt2 ) ;
	double vmod1 = rt1.modul();
	double vmod2 = rt2.modul();
	if ((vmod1 > 1) || (vmod2 > 1))
	{
	  valT1 = -1 ;
	  valT2 = -1 ;
	  return ;
	}
	// расчет памяти фильтра
	double arrL[4]={0},arrV[4] = {0},arrLamb[4]={0} ;
	arrL[0] =  1- arrP [0] ;
	arrL[1] =  mh * (1- arrP [0] );
	arrL[2] =  - arrP [1] ;
	arrL[3] =  1 - mh *  arrP [1] ;
	CalcProperVectors(arrL,arrV , arrLamb);
	valT1 = 0;
	valT2 = 0;
		double vmod1_n =  1;
		double vmod2_n =  1;

		double arr_q[2] = {0} ;
		MtrxTranspMultMatrx(arrV,2, 2, arrP,1, arr_q) ;  // arr_q = Ф(0)T * P
		double val1 = 0.05 * (fabs(arrV[0]  * arr_q[0] ) + fabs(arrV[1] *  arr_q[1]) ) ;
		double val2 = 0.05 * (fabs(arrV[2]  * arr_q[0] ) + fabs(arrV[3] *  arr_q[1]) ) ;
   // расчет памяти по положению
	while (1)
	{
		valT1 += mh;

		 vmod1_n *=vmod1;
		 vmod2_n *=vmod2;

		if ((vmod1_n * fabs(arrV[0] * arr_q[0]) + vmod2_n * fabs(arrV[1] * arr_q[1])  ) <=  val1) break ;

	}
	 // расчет памяти по скорости
	  vmod1_n =  1;
	  vmod2_n =  1;
	while (1)
	{
		valT2 += mh;

		 vmod1_n *=vmod1;
		 vmod2_n *=vmod2;
		if ((vmod1_n * fabs(arrV[2] * arr_q[0]) + vmod2_n * fabs(arrV[3] * arr_q[1])  ) <=  val2) break ;

	}

}
double  TStabSyst2 ::SolvCharactEqKlm(const double c)
{
	if (c >= 4)
	{
	   ShowMessage(L"C >= 4") ;
	   return - 1001;
	}
	return 1 + c /4 + sqrt(c /2 + c * c /16) ;
}

// поиск собственных векторов  и собственных чисел положительно поределенной матрица arrKInp
//  arrV - матрица собственных векторов
//  arrLamb - диагагнональная матрица собственных чисел
int  TStabSyst2 ::CalcProperVectors(double * arrKInp,double *arrV , double *arrLamb)
{
	
	memset(arrLamb,0,4 * sizeof(double)) ;
	double a  = 1;
	double b = -(arrKInp[0] + arrKInp [3]);
	double c = arrKInp[0] * arrKInp [3] - arrKInp[1] * arrKInp [1] ; ;
	TComp x0,x1;
	int i0 = SolvEq2( a, b, c,x0,x1) ;
	double arrx[2];
	arrx[0] = x0.m_Re;
	arrx[1] = x1.m_Re;

	switch(i0)
	{
		case 0:
		
		for (int i = 0; i < 2 ; i++)
		{
		  if  ( arrx[i] < 0)
		  {
			  int j = (i + 1)%2 ;
			  double arrProp[2] ={0},arrPropTemp[4] = {0} ;
			  if ( fabs (arrKInp[0] - arrx[j])> EPS)
			  {
				arrProp[0] = - arrKInp[1]/ (arrKInp[0] - arrx[j]);
				arrProp[1] = 1 ;
			  }
			  else
			  {
				arrProp[0] = 1 ;
				arrProp[1] =  -  (arrKInp[0] - arrx[j])/arrKInp[1];
              }
			  double d = sqrt(arrProp[0] * arrProp[0] + arrProp[1] * arrProp[1]) ;
			  for (int ii = 0; ii < 2; ii++) arrProp[ii] = arrProp[ii] / d ;
                  

			  MtrxMultMatrxTransp(arrProp,2, 2, arrProp,2, arrPropTemp) ;
			  MatrxMultScalar(arrPropTemp, 2 , 2, arrx[j],arrKInp);
			  return 1 ;
		  }
		}
		for (int i = 0; i < 2; i++)
		{
		 arrLamb[i * 2 + i] = arrx[i];
		 if (fabs(arrKInp[0] - arrx[i] ) > 0.00000000001 )
		 {
		   arrV [ i ] = - arrKInp[1] / (arrKInp[0] - arrx[i] );
		   arrV [ 2 + i ] = 1 ;
		   double r = sqrt (arrV [ i ] *  arrV [ i ] + arrV [ 2 + i ] * arrV [ 2 + i ] );
		   arrV [ i ] = arrV [ i ] /r ;
		   arrV [ 2 + i ] = arrV [ 2 + i ] /r ;

		 }
		 else
		 {
		   arrV [ i ]     = 1 ;
		   arrV [ 2 + i ] = 0 ;
         }
		}
		break ;
		case 1:
		arrV [0] = 1 ;
		arrV [1] = 0 ;
		arrV [2] = 0 ;
		arrV [3] = 1 ;
		arrLamb [0] = arrKInp[0] ;
		arrLamb [1] = 0 ;
		arrLamb [2] = 0 ;
		arrLamb [3] = arrKInp[0] ;
		break ;

		default:

		return 1;
		break ;
	}
	return 0 ;
}
// вывод параметров фильтра калмана в установившийся режим динамическим способом
// матирица B = {h*h/2 h}, вариант ДИСКРЕТНОГО БЕЛОГО ШУМА В КАНАЛЕ ОБЪЕКТА
void TStabSyst2 ::FltKlm( double *arrK,double * arrP)
{

  double arrKtemp[4] ,arrKtemp1[4];
  memcpy(arrKtemp,arrK,4 * sizeof(double)) ;


   for (int i = 0; i < 1000000; i++)
  {

	 OneStepFltrKlm(  arrKtemp,arrK, arrP) ;

	if ( (fabs(arrK[0] - arrKtemp [0]) < 0.001 ) && (fabs(arrK[3] - arrKtemp [3]) < 0.001 ))
	{
	  break;
	}

	 memcpy(arrKtemp,arrK,sizeof(double) * 4) ;

  }
}
// миатирица B = {h*h/2 h}, вариант ДИСКРЕТНОГО БЕЛОГО ШУМА В КАНАЛЕ ОБЪЕКТА
// расчет корреляционной матрицы и коэфф усиления фильтра Калмана для системы
// 2-го пеорядка на одном шаге
// Input:
// arrKInp - К(n)
// Output:
// arrKOut - K(n+1)
// arrP - P(n)
void TStabSyst2 ::OneStepFltrKlm( double *arrKInp,double *arrKOut,double * arrP)
{
  double arrKextr[4] = {0} , arrKt1[4] = {0},arrF[4] = {0} ;

  double arrKtemp[4] ,arrKtemp1[4];
  memcpy(arrKtemp ,arrKInp, 4 * sizeof(double)) ;
  arrF[0] = msigm_w * msigm_w * mh * mh * mh * mh/4 ;
  arrF [ 1] =  msigm_w * msigm_w * mh * mh * mh /2 ;
  arrF [ 2] =arrF [ 1] ;
  arrF [ 3] = msigm_w * msigm_w * mh * mh ;

	MtrxMultMatrx(marrA ,2, 2, arrKInp,2, arrKt1) ;
	MtrxMultMatrxTransp(arrKt1,2, 2,marrA ,2, arrKtemp1) ;
	MtrxSumMatrx(arrKtemp1, arrF,2,2, arrKextr) ;
	//
	arrP [0] = arrKextr[ 0 ]/( msigm_bmo * msigm_bmo + arrKextr[ 0 ]) ;
	arrP [1] = arrKextr[ 1 ]/( msigm_bmo * msigm_bmo + arrKextr[ 0 ]) ;
	arrKOut [0] =  arrKextr[ 0 ] - arrKextr[ 0 ] * arrKextr[ 0 ] / ( msigm_bmo * msigm_bmo + arrKextr[ 0 ]) ;
	arrKOut [1] =  arrKextr[ 1 ] - arrKextr[ 0 ] * arrKextr[ 1 ] / ( msigm_bmo * msigm_bmo + arrKextr[ 0 ]) ;
	arrKOut [2] = arrKOut [1] ;
	arrKOut [3] =  arrKextr[ 3 ] - arrKextr[ 1 ] * arrKextr[ 1 ] / ( msigm_bmo * msigm_bmo + arrKextr[ 0 ]) ;

}
// нахождение гарантированной точнолсти фильтрации для  фильтра Калмана
// НПРЕРЫВНЫЙ БЕЛЫЙ ШУМ В КАНАЛЕ НАБЛЮДЕНИЯ. ВАРИАНТ КНЯЗЕВА.
//Input:
// sigm_w, sigm_mmo, sigm_bmo   - реально действующие СКЗ
// arrK - нач корреляц матрица для фильтра
//
// Output :
// arrPStab - коефф усилнгия в установившемся режиме
// psi1,psi2 - углы формирования наихудших возмущений
// arrKTol - корреляционная функция(функция максимума)
// pparrK - двойной указатель на массив корреляционой матрицы по времени
// при наихудшем формированиия возмущений
// сначала идет массив K[0], затем К[1], затем K[3]
// plenparrarrK - общая длина массива parrK
//  Если pparrK == NULL, то заполнение массива *pparrK не производится
 bool TStabSyst2 ::CalcFuncMaximumKlm_Real( const double sigm_w,const double sigm_mmo
	  ,const double sigm_bmo ,double *arrK, double *arrPStab
	   , double &psi1, double &psi2,double *arrKTol , double **pparrK,int * plenparrarrK)
{
   memcpy(arrKTol,arrK,4 * sizeof(double)) ;
   const int iNet = 100;
   double arrK_t[4] ={0},arrKTol_t[4] ={0},arrKTolRez[4] = {0} ;
   double psi1_t = 0, psi2_t = 0 ;
   for (int i = 0; i < iNet; i++)
   {

   for (int j = 0; j < iNet; j++)
   {

	 memcpy(arrK_t,arrK,4 * sizeof(double)) ;
	 memcpy(arrKTol_t,arrKTol,4 * sizeof(double)) ;
	 psi1_t = 2 * PI * i/ iNet ;
	 psi2_t = 2 * PI * j/ iNet ;
	 CalcAccuracyKlm_Real( sigm_w,sigm_mmo, sigm_bmo, arrK_t, arrPStab
	   , psi1_t,psi2_t, arrKTol_t , NULL,NULL);
	 if (arrKTol_t[3] >  arrKTolRez[3])
	 {
	   memcpy(arrKTolRez, arrKTol_t,4 * sizeof(double)) ;
	   psi1 = psi1_t;
	   psi2 = psi2_t;
	 }
   }
   }
   memcpy(arrK_t,arrK,4 * sizeof(double)) ;
   memcpy(arrKTol_t,arrKTol,4 * sizeof(double)) ;
   CalcAccuracyKlm_Real( sigm_w,sigm_mmo, sigm_bmo, arrK_t, arrPStab
	   , psi1,psi2, arrKTol_t , pparrK, plenparrarrK);
   memcpy(arrK, arrK_t,4 * sizeof(double)) ;
   memcpy(arrKTol, arrKTolRez,4 * sizeof(double)) ;

  return true;
}
 // НЕПРЕРЫВНЫЙ БЕЛЫЙ ШУМ В КАНАЛЕ НАБЛЮДЕНИЯ. ВАРИАНТ КНЯЗЕВА.
// расчет корреляционной матрицы ошибок для фильтра Kaлмана настроенного на белые шумы
// , когда возмущения формируются с углами psi1 и psi2
// arrK0 - на входе априорная корреляц матрица для фильтра ,
// на выходе корреляционная матрица ошибок фильтрации в установившемся режиме с фильтра
// &plenparrarrK - начальная величина зарезервированной памяти под  *pparrK
// arrPStab - коэффициенты усиления фиольтра в установивишемся режиме
// psi1 и psi2 - углы формирования процессов
// arrKTol - на входе априорная корреляционная матрицв, на выходе
// корреляционная матрица ошибок фильтрации в установившемс ярежиме
//иначе говоря, это расчет функции риска - т к фильтр фиксирован и возмущения фиксированы
// pparrK - двойной указатель на массив корреляционой матрицы по времени
// сначала идет массив K[0], затем К[1], затем K[3]
// plenparrarrK - общая длина массива parrK
//  Если pparrK == NULL, то заполнение массива *pparrK не производится
bool TStabSyst2 ::CalcAccuracyKlm_Real(const double sigm_w,const double sigm_mmo,const double sigm_bmo
	   ,double *arrK, double *arrPStab
	   ,const double psi1,const double psi2,double *arrKTol , double **pparrK,int * plenparrK)
{
	int lenMemAllocated = 0 ;
   if (pparrK != NULL) 	 lenMemAllocated = *plenparrK;

   double arrK_t[4] = {0} ;
   double arrKTol_t[4] = {0} ;
   double arrS_Flt[2] = {0}, arrQ_Flt[2] = {0} ;
   double arrS_Tol[2] ={0}, arrQ_Tol[2] ={0} ;
   int i = 0;
   for ( i = 0; i < 10000; i++)
   {

		OneStepFltrKlm_Real( arrK,arrK_t, arrPStab) ;
	 // пересчет корреляц матрицы точности

	 if ( !OneStepRecalcTol( psi1,psi2,sigm_w, sigm_mmo, sigm_bmo, arrPStab,arrKTol,arrKTol_t
		,arrS_Tol,arrQ_Tol) )
		{
		 memset( arrKTol,0,4 * sizeof(double)) ;
		 return false;
		}
	// проверка окончания процесса по невязке
	if ((fabs( arrK[0] - arrK_t[0]) < EPS) && (fabs( arrK[3] - arrK_t[3]) < EPS)
	  &&(fabs( arrKTol[0] - arrKTol_t[0]) < EPS/10) && (fabs( arrKTol[3] - arrKTol_t[3]) < EPS/10) )
	{
	  memcpy( arrK,arrK_t,4 * sizeof(double)) ;
	  memcpy( arrKTol,arrKTol_t,4 * sizeof(double)) ;
	  break ;
	}
   memcpy( arrK,arrK_t,4 * sizeof(double)) ;
   memcpy( arrKTol,arrKTol_t,4 * sizeof(double)) ;
	if (pparrK != NULL)
	{




			if (i * 3 > lenMemAllocated -3)
			{
			lenMemAllocated += 1000;
			*pparrK =  (double *)realloc(*pparrK,lenMemAllocated * sizeof(double)) ;
			}

			(*pparrK)[i * 3     ] = arrKTol[0] ;
			(*pparrK)[i * 3 + 1 ] = arrKTol[1] ;
			(*pparrK)[i * 3 + 2 ] = arrKTol[3] ;

	}

   }
	if (pparrK != NULL)
	{
	 *pparrK =  (double *)realloc(*pparrK,i * 3 * sizeof(double)) ;
	*plenparrK = i * 3 ;
	double *parrTemp = new double [ *plenparrK];
	MatrTransp(*pparrK, i, 3, parrTemp);
	memcpy(*pparrK, parrTemp, i * 3 * sizeof(double)) ;
	delete parrTemp ;
	}


   return true ;
}
//  пересччет корреляционной матирицы ошибок фильтрации при формировании возмущений с углами psi1, psi2
// на шаге фильтрации с коеффициентами усиления arrP
 //Задана:
 // 1. корреляционная матрица ошибок фильтрации на n- ом шаге K(n)
 // 2. коэффициенты усиления насчитанные на какм-то фильтре P
 // sigm_w, sigm_mmo, sigm_bmo, h, - параметры системы
 // psi1 -угол формирования возмущения w - между w и остатком фильтрации по положению
 // psi2 -угол формирования возмущения ММО - между ММО и остатком фильтрации по положению
 // Функция выполняет следующие дейчствия:
 // 1.находит arrS [2] -вектор разложения w по остаткам фильтрации, соответствующий углу psi1
 // 2. производит экстраполяцию корреляционной матрицы  K(n+1|n)=(A + BST)K(n)(A + BST)T
 // 3. находит arrQ[2] - вектортразложения ММО по остаткам фильтрации, соответствующий углу psi2
 // 4. персчитывает корреляционную матрицу:
 //     K(n+1)= (E - PCT -PqT)K(n+1|n)(E - PCT -PqT)T + sigm_bmo*sigm_bmo*P*PT
 // INPUT:
// arrKInp - корреляционная матр K(n)
// Input : psi1,psi2,sigm_w, sigm_mmo, sigm_bmo, h, arrP,arrKInp
//OutPut:
// arrKInp - корр матр K(n+1)
// arrS,arrQ - векторы разложений возмущений

bool TStabSyst2::OneStepRecalcTol( const double psi1,const double psi2,const double sigm_w
   ,const double  sigm_mmo,const double sigm_bmo,  double *arrP
		,double *arrKInp,double *arrKOut,double *arrS,double *arrQ)
{
  double arrKExtr[4] = {0} ;
  double arrE[4] = {1,0,0,1} ;

   double arrT1[4]={0},arrT2[4] ={0},arrT0[4] = {0},arrT3[4]={0} ;
   double arrT4[4]={0},arrT5[4] ={0},arrT6[4] = {0}
	   ,arrT7[4]={0},arrT8[4]={0}, arrT9[4]={0},arrT10[4] ={0},arrT11[4] = {0};
   // Экстраполяция
	if ( !CalcCoeff( psi1,sigm_w,arrKInp,arrS)) return false;

	MtrxMultMatrxTransp(marrB,2, 1, arrS, 2, arrT1) ;  // arrT1 = B * ST
	MtrxSumMatrx(marrA, arrT1,2, 2, arrT2); // arrT2 = A + B*(ST)
	// персчет K(n+1|n)
	MtrxMultMatrx(arrT2,2, 2, arrKInp,2, arrT3) ;// arrT3 = (A + B*(ST)) * K
	MtrxMultMatrxTransp(arrT3,2, 2, arrT2,2, arrKExtr) ;// arrKExtr = (A + B*(ST)) * K *  (A + B*(ST))T
	if ( (arrKExtr[0]< EPS) || (arrKExtr[3]< EPS))
	{
	  // SmhowMessage(L"(arrKExtr[0]< EPS) || (arrKExtr[3]< EPS)" );
	  return false ;
	}

		// Определение ММО
	if( sigm_mmo > 0.000000001)
	{
		if ( !CalcCoeff(psi2,sigm_mmo,arrKExtr,arrQ)) return false;
	}
	else memset(arrQ , 0, 2 * sizeof(double)) ;
	   ///

		MtrxMultMatrxTransp (arrP,2, 1, marrC,2, arrT4) ;// arrT4  = P * CT
		MtrxMinusMatrx(arrE, arrT4,2, 2, arrT5); //arrT5 = E - P*CT
		MtrxMultMatrxTransp(arrP,2, 1, arrQ,2, arrT6) ; // arrT6 = P * QT
		MtrxMinusMatrx(arrT5 , arrT6,2,2, arrT7);      // arrT7 = E - P*CT - P * QT
		MtrxMultMatrxTransp(arrP,2, 1, arrP,2, arrT8) ; // arrT8 = P * PT
		// персчет матрицы персчет K(n+1)
		MtrxMultMatrx(arrT7,2, 2, arrKExtr,2, arrT9) ;// arrT9 = (E - P*CT - P * QT) * K (n+1|n)
		MtrxMultMatrxTransp(arrT9,2, 2, arrT7,2, arrT10) ;// arrT10 = (E - P*CT - P * QT) * K (n+1|n) *(E - P*CT - P * QT)T
		MatrxMultScalar(arrT8 , 2,2, sigm_bmo * sigm_bmo,arrT11); // arrT11 = sigm_bmo * sigm_bmo m* P * PT
		MtrxSumMatrx(arrT10, arrT11,2,2, arrKOut) ;

	return true ;

}
 // Расчет разложения возмущения при задангном угле между возмущением и
 // остатком фильторации по положению
// INPUT :
// fi - угол между возмущением и остатком по положению
// sigm - СКЗ возмужщения
// arrK - корреляционная матрица  остатков фильтрации по положению  и по скорости
// OUTPUT:
// arrS [2] - вектор разложения возмущения по остаткам фильтпрации
bool  TStabSyst2 ::CalcCoeff(double fi,double sigm,double *arrK,double *arrS)
{
  double cs = cos(fi) ;
  double sn = sin(fi) ;
  double sq0 = sqrt(arrK[0]) ;
  if((arrK[0] < EPS) || (arrK[3] < EPS))  return false ;

  double sq2 = arrK[3]  - arrK[1] * arrK[1]/arrK[0];
  if( sq2 < EPS)
  {
	if(fabs(fi) < EPS)
	{
	  arrS[ 0 ] =   sigm/sq0;
	  arrS[ 1 ] = 0;
	  return true;
	}
	if( fabs(fi - PI) < EPS)
	{
	  arrS[ 0 ] =   -sigm/sq0;
	  arrS[ 1 ] = 0;
	  return true;
    }
	return false;
  }
  else
  {
	double sq = (double)sqrt( sq2) ;
	arrS[ 0 ] =  (cs -  arrK[1]*sn/ sq/sq0) *  sigm/sq0 ;
	arrS[ 1 ] =  sn* sigm /sq;
	return true ;
  }
}


// аналитический расчет фильтра Калмана сопровождения в установиви режиме
// Input:

// pparrWeight - двойной указатель на массив весовой йункции
// lenparrWeight - начальная длина массива весов функции
// Output:
// arrK - коррелл матрица
// arrP - коефф усиления
// rt1, rt2 - корни характерист уравнения
// valT1, valT2 - память фильтра по положению и по скорости
// если система неустойчивая, то  valT1 = valT2 =  -1
// pparrWeight0 - двойной указатель на массив весовой  функции  по подложениею
// lenparrWeight - окончательная  длина массива весов функций  по подложениею и скорости
// pparrWeight1 - двойной указатель на массив весовой  функции  по скорости
//  pparrOut -двоЗйной указатель на массивы P1,P2,K[0], K[1],K[3]
// lenparrWeight - общая длина массива  (*pparrOut)
// lenparrWeight  = 5 * lenparrWeight
void TStabSyst2 ::StabSolutionKalm_Real( double *arrK,double * arrP,TComp &rt1,TComp &rt2,
	 double &valT1, double &valT2,double **pparrOut, int *plenparrOut, double **pparrWeight0
	 , double **pparrWeight1,int * plenparrWeight)
{
	double mu =  msigm_w* mh * sqrt(mh) / msigm_bmo;
	double x = SolvCharactEqKlm_Real(mu);
	arrK[0] =  msigm_bmo* msigm_bmo* (x * x -1 )/ x / x;
	arrK[1] =  msigm_bmo* msigm_w* sqrt(mh )/ x ;
	arrK[2] =arrK[1] ;
	arrK[3] = msigm_w* msigm_bmo* ( x * x -1 ) /x / sqrt(mh ) - msigm_w* msigm_w* mh  /2 ;
	arrP [0] =  (x * x -1 )/ x / x;
	arrP [1] =  msigm_w * sqrt( mh ) / msigm_bmo/ x ;
	// расчет корней характ уравнения
	double a1 =1;
	double b1 = 2 - arrP [0] - mh * arrP [1] ;
	double c1 = 1- arrP [0] ;
	SolvEq2( a1, b1, c1,rt1,rt2 ) ;
	double vmod1 = rt1.modul();
	double vmod2 = rt2.modul();
	if ((vmod1 > 1) || (vmod2 > 1))
	{
	  valT1 = -1 ;
	  valT2 = -1 ;
	  return ;
	}
	// расчет памяти фильтра
	double arrL[4]={0},arrLamb[4]={0}, arrV[4] ={0} ;
	arrL[0] =  1- arrP [0] ;
	arrL[1] =  mh * (1- arrP [0] );
	arrL[2] =  - arrP [1] ;
	arrL[3] =  1 - mh *  arrP [1] ;
	CalcProperVectors(arrL,arrV , arrLamb);
	valT1 = 0;
	valT2 = 0;
		double vmod1_n =  1;
		double vmod2_n =  1;

		double arr_q[2] = {0} ;
		MtrxTranspMultMatrx(arrV,2, 2, arrP,1, arr_q) ;  // arr_q = Ф(0)T * P
		double val1 = 0.05 * (fabs(arrV[0]  * arr_q[0] ) + fabs(arrV[1] *  arr_q[1]) ) ;
		double val2 = 0.05 * (fabs(arrV[2]  * arr_q[0] ) + fabs(arrV[3] *  arr_q[1]) ) ;
   // расчет памяти по положению
	while (1)
	{
		valT1 += mh;

		 vmod1_n *=vmod1;
		 vmod2_n *=vmod2;

		if ((vmod1_n * fabs(arrV[0] * arr_q[0]) + vmod2_n * fabs(arrV[1] * arr_q[1])  ) <=  val1) break ;

	}
	 // расчет памяти по скорости
	  vmod1_n =  1;
	  vmod2_n =  1;
	while (1)
	{
		valT2 += mh;

		 vmod1_n *=vmod1;
		 vmod2_n *=vmod2;
		if ((vmod1_n * fabs(arrV[2] * arr_q[0]) + vmod2_n * fabs(arrV[3] * arr_q[1])  ) <=  val2) break ;

	}
   if  (pparrWeight0 != NULL )
   {
	CalcWeightArray( arrP, pparrWeight0
	 ,pparrWeight1, plenparrWeight);
	 *plenparrOut = 5 *  (*plenparrWeight) ;
	 *pparrOut = (double*)realloc(*pparrOut, 5 * (*plenparrWeight) * sizeof(double)) ;

	  double arrKtemp[4] = {100000,50000,50000,100000},arrKtemp1[4] = {0};



   for (int i = 0; i < (*plenparrWeight); i++)
  {

	 OneStepFltrKlm_Real(  arrKtemp,arrKtemp1, arrP) ;
	 (*pparrOut)[                       i ] = arrP[0] ;
	 (*pparrOut)[(*plenparrWeight)    + i ] = arrP[1] ;
	 (*pparrOut)[(*plenparrWeight)*2  + i ] = arrKtemp1[0] ;
	 (*pparrOut)[(*plenparrWeight)*3  + i ] = arrKtemp1[1] ;
	 (*pparrOut)[(*plenparrWeight)*4  + i ] = arrKtemp1[3] ;


	 memcpy(arrKtemp,arrKtemp1,sizeof(double) * 4) ;

  }

   }

}
//решение нелинейного уравнения для реального фильтра калмана
double  TStabSyst2 ::SolvCharactEqKlm_Real(const double mu)
{
	double t0 = mu/2 + sqrt(mu * mu /12 + 4) ;
	return t0/2 + sqrt(t0*t0/4 -1) ;
}

// вывод параметров реального фильтра калмана в установившийся режим динамическим способом
void TStabSyst2 ::FltKlm_Real(   double *arrK,double * arrP)
{

  double arrKtemp[4] ,arrKtemp1[4];
  memcpy(arrKtemp,arrK,4 * sizeof(double)) ;


   for (int i = 0; i < 10000000; i++)
  {

	 OneStepFltrKlm_Real(  arrKtemp,arrK, arrP) ;

	if ( (fabs(arrK[0] - arrKtemp [0]) < 0.00001 ) && (fabs(arrK[3] - arrKtemp [3]) < 0.0000001 ))
	{
	  break;
	}

	 memcpy(arrKtemp,arrK,sizeof(double) * 4) ;

  }
}
// расчет корреляционной матрицы и коэфф усиления реального фильтра Калмана для системы
// 2-го пеорядка на одном шаге
// Input:
// arrKInp - К(n)
// Output:
// arrKOut - K(n+1)
// arrP - P(n)
void TStabSyst2 ::OneStepFltrKlm_Real(double *arrKInp,double *arrKOut,double * arrP)
{
  double arrKextr[4] = {0} , arrKt1[4] = {0},arrF[4] = {0} ;

  double arrKtemp[4] ,arrKtemp1[4];
  memcpy(arrKtemp ,arrKInp, 4 * sizeof(double)) ;



  arrF[0] = msigm_w* msigm_w* mh * mh * mh /3 ;
  arrF [ 1] =  msigm_w* msigm_w* mh * mh  /2 ;
  arrF [ 2] =arrF [ 1] ;
  arrF [ 3] = msigm_w* msigm_w* mh  ;

	MtrxMultMatrx(marrA ,2, 2, arrKInp,2, arrKt1) ;
	MtrxMultMatrxTransp(arrKt1,2, 2,marrA ,2, arrKtemp1) ;
	MtrxSumMatrx(arrKtemp1, arrF,2,2, arrKextr) ;
	arrP [0] = arrKextr[ 0 ]/( msigm_bmo * msigm_bmo + arrKextr[ 0 ]) ;
	arrP [1] = arrKextr[ 1 ]/( msigm_bmo * msigm_bmo + arrKextr[ 0 ]) ;
   //	MtrxMultMatrxTransp(arrP,2, 1, arrC,2, arrKt1) ;
   //	MtrxMultMatrx(arrKt1,2, 2, arrKextr,2, arrKtemp1) ;
   //	MtrxMinusMatrx(arrKextr, arrKtemp1,2, 2, arrK);
	arrKOut [0] =  arrKextr[ 0 ] - arrKextr[ 0 ] * arrKextr[ 0 ] / ( msigm_bmo * msigm_bmo + arrKextr[ 0 ]) ;
	arrKOut [1] =  arrKextr[ 1 ] - arrKextr[ 0 ] * arrKextr[ 1 ] / ( msigm_bmo * msigm_bmo + arrKextr[ 0 ]) ;
	arrKOut [2] = arrKOut [1] ;
	arrKOut [3] =  arrKextr[ 3 ] - arrKextr[ 1 ] * arrKextr[ 1 ] / ( msigm_bmo * msigm_bmo + arrKextr[ 0 ]) ;

}

// расчет весовой фуекции фильтра с коеф усиления arrP
// pparrWeight0 - двойной указатель на массив весовой йункции по положениею
// lenparrWeight0 - окончат7 ельная  длина массива весов функции  по положению
// pparrWeight1 - двойной указатель на массив весовой йункции по скорости
// lenparrWeight1  - окончат7 ельная  длина массива весов функции  по скорости
void TStabSyst2 ::CalcWeightArray(double * arrP, double **pparrWeight0
	 ,double **pparrWeight1, int *plenparrWeight)
{
     double arrL[4] ={0} ;
	arrL[0] =  1- arrP [0] ;
	arrL[1] =  mh * (1- arrP [0] );
	arrL[2] =  - arrP [1] ;
	arrL[3] =  1 - mh *  arrP [1] ;
		// рассчет весовой функцмии фильтра
	double arrLN [4] = {1,0,0,1} ;//L  в степени N
	double arrT[4] ={0},arrWT[2]={0} ;
	int it = 0;
	double valW0 = fabs(arrP[0]), valW1 = fabs(arrP[1]);
	int lenMemory = *plenparrWeight;
	while(1)
	{
	  it++;
	  if (it > lenMemory )
	  {
		*pparrWeight0 = (double *)realloc(*pparrWeight0 , (lenMemory + 1000) * sizeof(double));
		*pparrWeight1 = (double *)realloc(*pparrWeight1 , (lenMemory + 1000) * sizeof(double));
		lenMemory  += 1000;
	  }
	  MtrxMultMatrx(arrL,2,2, arrLN,2, arrT) ;
	  memcpy(arrLN,arrT, 4 * sizeof(double)) ;// arrLN =  arrL ^ N
	  MtrxMultMatrx(arrLN,2,2, arrP,1,arrWT) ;
	  (*pparrWeight0)[ it-1 ] = arrWT [0] ;
	  (*pparrWeight1)[ it -1] = arrWT [1] ;
	  if( (fabs(arrWT [0]) < 0.01 * fabs(valW0))&& (fabs(arrWT [1]) < 0.01 * fabs(valW1))) break;
	  //valW0 = arrWT [0] ;
	  //valW1 = arrWT [1] ;

	}
		*pparrWeight0 = (double *)realloc(*pparrWeight0 ,it * sizeof(double));
		*pparrWeight1 = (double *)realloc(*pparrWeight1 , it * sizeof(double));
		*plenparrWeight = it ;
}

void TStabSyst2 ::CreateShpFile(wchar_t *wchFileName, double *parrInf
	 ,const int lenarr, const int quantParts,const double h,double &scalex, double &scaley)
{
	int irez =0 ;
	double maxarr = maxDoubleArr(parrInf,  lenarr, irez);
	const int quanCols = lenarr / quantParts ;
	if (scaley < 0)
	{
	scalex = 200/h/ quanCols;
	scaley = 100/ maxarr ;
	}
	else
	{    if (scalex < 0)
		 {
           scalex =  maxarr *3/h/ quanCols;
		 }

	}
	TURPointXY *pPoints = new TURPointXY[lenarr ];
    int iii = 0;
	for (int i = 0; i < lenarr ; i++)

	{
	  int j = i % quanCols    ;
	  pPoints [ i ].Y = parrInf[i]*scaley ;
	  pPoints [ i ].X = j * h * scalex;
	}


	//
	int numParts   = quantParts  ;
	int numPoints  = lenarr ;
	int *iarrParts = new int[ numParts] ;
	for (int i = 0; i < quantParts; i++) iarrParts[i] = i * quanCols ;

	TURPolyLine Plyline( numParts, numPoints ,iarrParts, pPoints);
	TURPolyLine::WriteSetSHPFiles(wchFileName,&Plyline, 1) ;
//Plyline.WriteToASCII(L"D:\\PROJECTS_C++\\Pict\\One\\K22.txt");

	delete [] pPoints;
	delete  iarrParts;
}
// создание shape файла с осями координат
void TStabSyst2 ::CreateShpAxes(wchar_t *wchFileName,const double xmin,const double xmax
	 ,const double ymin,const double ymax)

{
	int numParts   =  2 ;
	int numPoints  =  4 ;
	int iarrParts[2] = {0,2};
	 TURPointXY arrPoints[4] ;
	 arrPoints[0] =  TURPointXY(xmin,0);
	 arrPoints[1] =  TURPointXY(xmax,0);
	 arrPoints[2] =  TURPointXY(0,ymin);
	 arrPoints[3] =  TURPointXY(0,ymax);

	TURPolyLine Plyline( numParts, numPoints ,iarrParts, arrPoints);
	TURPolyLine::WriteSetSHPFiles(wchFileName,&Plyline, 1) ;
}
double TStabSyst2 ::maxDoubleArr(double *parr, const int lenarr, int &irez)
{
  double rez = -DBL_MAX ;
  for (int i =0; i < lenarr; i++)
  {
	  if (parr[i] > rez)
	  {
		rez = parr[i];
		irez = i ;
	  }
  }
  return rez ;
}
//******************************************************************************************************************//
//******************************************************************************************************************//
//******************************************************************************************************************//
//************ ФУНКЦИИ МИНИМАКСНОЙ ФИЛЬТРАЦИИ        ******************************************************************************************************//
//******************************************************************************************************************//
//******************************************************************************************************************//
//******************************************************************************************************************//

// Реализация фильтра Голубева в чистом виде, в соответсвии с формулами Брисенко
// если pparrOut == NULL, то расчитываются тоько arrKStab и arrPStab
// если pparrOut != NULL, то расчитываются pparrOut - массив коэффициентов
// усиления и корреляционных матриц,  в котором информация хранится в следующем порядке:
// P1,P2,K11,K12,K22
// длина каждогго массива *plenparrWeight
// длина массива (*pparrOut ) равна *plenparrOut = 5 *(* plenparrWeight)
 // *pparrWeight0 и *pparrWeight1 - массивы с весовыми функциями длины *plenparrWeight
 // Если требуется сформипровать массивы *pparrOut, *pparrWeight0 и *pparrWeight1
 // необходимо вне программы создать двойные указатели на массивы длиной *plenparrOut и  *plenparrWeight
// в массиве pparrOut сохранены эти массивы один за другим.
// Input:
// arrKStab - начальная коррел матрица
//Output :
// arrKStab-  коррел матрица  в установивишеимся режиме
// arrPStab - коефииц усиления в установившемся режиме

 void  TStabSyst2 ::CalcCorrMatrGolubev(double  *arrKStab,double *arrPStab
		,double **pparrOut, int *plenparrOut, double **pparrWeight0
		,double **pparrWeight1,int * plenparrWeight)
{
  double  arrK1[4] = {1000000,250000,250000,250000};

 // memcpy( arrK1,arrKStab, 4 * sizeof(double)) ;
  memcpy(arrKStab,arrK1, 4 * sizeof(double)) ;
  int i ;
   if  (pparrWeight0 == NULL )
   {
	for ( i = 0; i < 10000000; i++)
	{

	  OneStepGolubev(arrKStab,arrPStab);
	  if ( (fabs(arrKStab[0] - arrK1[0]) < EPS) && (fabs(arrKStab[3] - arrK1[3]) < EPS) )  break;
	   memcpy( arrK1,arrKStab, 4 * sizeof(double)) ;
	}
   }
   else
   {

     for ( i = 0; i < 10000000; i++)
	{

	  OneStepGolubev(arrKStab,arrPStab);
	  if ( (fabs(arrKStab[0] - arrK1[0]) < EPS) && (fabs(arrKStab[3] - arrK1[3]) < EPS) )  break;
	   memcpy( arrK1,arrKStab, 4 * sizeof(double)) ;
	}

	CalcWeightArray( arrPStab, pparrWeight0
	 ,pparrWeight1, plenparrWeight);
	 *plenparrOut = 5 *  (*plenparrWeight) ;
	 *pparrOut = (double*)realloc(*pparrOut, 5 * (*plenparrWeight) * sizeof(double)) ;

	  double arrKtemp[4] = {100000,50000,50000,100000},arrKtemp1[4] = {0};



   for (int i = 0; i < (*plenparrWeight); i++)
  {

	 OneStepGolubev(arrKtemp,arrPStab);//(  arrKtemp,arrKtemp1, arrP) ;
	 (*pparrOut)[                       i ] = arrPStab[0] ;
	 (*pparrOut)[(*plenparrWeight)    + i ] = arrPStab[1] ;
	 (*pparrOut)[(*plenparrWeight)*2  + i ] = arrKtemp[0] ;
	 (*pparrOut)[(*plenparrWeight)*3  + i ] = arrKtemp[1] ;
	 (*pparrOut)[(*plenparrWeight)*4  + i ] = arrKtemp[3] ;

  }


   }
}
// шаг пересчета фильтра Голубева в чистом виде  , в соответсвии с формулами Борисенко
void 	TStabSyst2 ::OneStepGolubev(double *arrK,double *arrP)
{
	// Экстраполяция
	double arrKExtr[4] ;
	double sq = sqrt(arrK [3] - arrK [1] * arrK [1] /arrK [0]) ;
	arrKExtr[0] =  arrK [0] + 2 * mh * arrK[1] + mh * mh * arrK[3] +  msigm_w * msigm_w * mh* mh * mh * mh/4 + msigm_w * mh * mh * mh * sq;
	arrKExtr[1] =  arrK [1] + mh *arrK [3] + 3 * msigm_w * mh * mh * sq/2 + msigm_w * msigm_w * mh * mh * mh /2;
	arrKExtr[2] = arrKExtr[1] ;
	arrKExtr[3] =  arrK [3] + 2 * msigm_w * mh * sq + msigm_w * msigm_w * mh * mh ;
	// Выбор ММО
	double alf =  msigm_mmo/ sqrt( arrKExtr[0]) ;
	if (alf > 1) alf = 1 ;
	double valF = ( 1- alf)/( ( 1- alf)* ( 1- alf) * arrKExtr[0] + msigm_bmo * msigm_bmo) ;

	arrP[0] = arrKExtr[0]* valF ;
	arrP[1] =  valF * arrKExtr[1] ;
	if (arrP[0] > 1)
	{
	  arrP[0] = 1 ;
	  arrP[1] = arrKExtr[1] / arrKExtr[0];
	  alf = ( msigm_mmo * msigm_mmo + msigm_bmo * msigm_bmo) / arrKExtr[0] ;
	}
	double valP = 1 - (1- alf) * arrP[0] ;
	arrK [ 0] = valP * arrKExtr[0] ;
	arrK [ 1] = valP * arrKExtr[1] ;
	arrK [ 2] = arrK [ 1] ;
	arrK [ 3] = arrKExtr[3] -  ( 1 - alf) * arrP[1] * arrKExtr[1] ;
}

void TStabSyst2::CalcCorrMatrGolubev_flAdd(double  *arrKStab,double *arrPStab)
{
  double  arrK1[4] = {1000000,250000,250000,250000};

 // memcpy( arrK1,arrKStab, 4 * sizeof(double)) ;
  memcpy(arrKStab,arrK1, 4 * sizeof(double)) ;
  int i ;


	for ( i = 0; i < 10000000; i++)
	{

	  OneStepGolubev_flAdd(arrKStab,arrPStab);
	  if ( (fabs(arrKStab[0] - arrK1[0]) < EPS) && (fabs(arrKStab[3] - arrK1[3]) < EPS) )  break;
	   memcpy( arrK1,arrKStab, 4 * sizeof(double)) ;
	}



}

// шаг пересчета фильтра Голубева c добавкой в уравнении движения в пр части по положнению
void 	TStabSyst2 ::OneStepGolubev_flAdd(double *arrK,double *arrP)
{
	// Экстраполяция
	double arrKExtr[4] ;
	double sq = sqrt(arrK [3] - arrK [1] * arrK [1] /arrK [0]) ;
	arrKExtr[0] =  arrK [0] + 2 * mh * arrK[1] + mh * mh * arrK[3] +  msigm_w * msigm_w * mh* mh * mh * mh/4 + msigm_w * mh * mh * mh * sq;
	arrKExtr[0] = arrKExtr[0] + msigm_fl*msigm_fl;
	arrKExtr[1] =  arrK [1] + mh *arrK [3] + 3 * msigm_w * mh * mh * sq/2 + msigm_w * msigm_w * mh * mh * mh /2;
	arrKExtr[2] = arrKExtr[1] ;
	arrKExtr[3] =  arrK [3] + 2 * msigm_w * mh * sq + msigm_w * msigm_w * mh * mh ;
	// Выбор ММО
	double alf =  msigm_mmo/ sqrt( arrKExtr[0]) ;
	if (alf > 1) alf = 1 ;
	double valF = ( 1- alf)/( ( 1- alf)* ( 1- alf) * arrKExtr[0] + msigm_bmo * msigm_bmo) ;

	arrP[0] = arrKExtr[0]* valF ;
	arrP[1] =  valF * arrKExtr[1] ;
	if (arrP[0] > 1)
	{
	  arrP[0] = 1 ;
	  arrP[1] = arrKExtr[1] / arrKExtr[0];
	  alf = ( msigm_mmo * msigm_mmo + msigm_bmo * msigm_bmo) / arrKExtr[0] ;
	}
	double valP = 1 - (1- alf) * arrP[0] ;
	arrK [ 0] = valP * arrKExtr[0] ;
	arrK [ 1] = valP * arrKExtr[1] ;
	arrK [ 2] = arrK [ 1] ;
	arrK [ 3] = arrKExtr[3] -  ( 1 - alf) * arrP[1] * arrKExtr[1] ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Вывод корреляционной матрицы и креффициентов усиления на стационарный режим
// для фильтра с углами формиорвания fi1 и fi2
// Input:
//  fi1  - угол между w и остатком x1
//  fi2  - угол между ммо и остатком x1
//  arrK  - начальная корреляц матрица
//  Output:
//  arrK  - корреляц матрица в установившемся режиме
//  arrP - коефф усиления   в установившемся режиме
//  arrS - вектор формирования w ( w = ST * Xостаток)   в установившемся режиме
//  arrQ - вектор формирования ммо ( ммо = ST * Xостаток)   в установившемся режиме
// Возмущения заданы формированиями. Под эти ворзмущения ищется
// оптимаольный фильтр и выводится в установивишмийся режим
bool TStabSyst2::CalcFuncMinimum(double fi1,double fi2,double *arrK, double *arrP,double *arrS,double *arrQ)
{

  double arrKtemp[4] ;
  memcpy(arrKtemp,arrK,4 * sizeof(double)) ;



  for (int i = 0; i < 10000; i++)
  {

	 if (!OneStepAlfaBettaFlt( fi1, fi2,arrK,arrKtemp,arrP,arrS,arrQ))
	 {
		return false;
	 }

	if( (fabs(arrKtemp[0] -arrK[0])< EPS) && (fabs(arrKtemp[3] -arrK[3])< EPS) )
	{
	  memcpy(arrK,arrKtemp,4 * sizeof(double)) ;
	  break;
	}
	memcpy(arrK,arrKtemp,4 * sizeof(double)) ;


  }
  memcpy(arrK,arrKtemp,4 * sizeof(double)) ;
  if (  ((fabs(arrP[0])< 0.00000001) ||(fabs(arrP[1])< 0.00000001))) return false ;
  TComp rt1, rt2 ;
  if (!FindRoots( arrP, arrS,arrQ, rt1,rt2))
   return false ;  ;

}
// персчет корреляционной матрицы  фильтра с углами формирования fi1 и fi2
// fi1, fi2 - углы формирований возмущений w и ммо
bool TStabSyst2::OneStepAlfaBettaFlt(double fi1,double fi2,double *arrKInp,double *arrKOut, double *arrP
		,double *arrS,double *arrQ)
{
  // формирование матриц динам системы
  double arrL[4] = {1,0,0,1};
  arrL[1] = mh ;
  double arrB[2] ;
  arrB[0] = mh*mh/2;
  arrB[1] = mh ;
  double arrC[2] = {1,0};
  //
  double arrKtemp[4] ;
  memcpy (arrKtemp, arrKInp, 4 * sizeof(double)) ;

  double arrTemp[4]= {0} ;
  double arrTemp1[4]= {0} ;
  double arrTemp2[4]= {0} ;
  double arrKExtr[4],arrR[2],valN, arr_C_Plus_Q[2];
 // Экстраполяция
	CalcCoeff(fi1,msigm_w,arrKtemp,arrS);

	MtrxMultMatrxTransp(arrB,2, 1, arrS, 2, arrTemp) ;  // arrTemp = B * ST
	MtrxSumMatrx(arrL, arrTemp,2, 2, arrTemp1); // arrTemp1 = L + B*(ST)
	MtrxMultMatrx(arrTemp1,2, 2, arrKtemp,2, arrTemp) ;// arrTemp = (L + B*(ST)) * K
	MtrxMultMatrxTransp(arrTemp,2, 2, arrTemp1,2, arrKExtr) ;// arrKExtr = (L + B*(ST)) * K *  (L + B*(ST))T
	arrKExtr[1] = (arrKExtr[1] + arrKExtr[2])/2 ;
	arrKExtr[2] = arrKExtr[1] ;
	// Определение ММО
	CalcCoeff(fi2,msigm_mmo,arrKExtr,arrQ);
	if ( arrQ [ 0] < -1)
	{
	   arrQ [ 0] = - 1;
	   arrQ [ 1] = 0 ;
	}
	MtrxSumMatrx(arrC, arrQ,2, 1, arr_C_Plus_Q);  //arr_C_Plus_Q =  C + Q
	MtrxTranspMultMatrx(arr_C_Plus_Q,2, 1, arrKExtr,2, arrTemp2) ; // arrTemp2 = (C + Q)T * KExtr
	MtrxMultMatrx(arrTemp2,1, 2, arr_C_Plus_Q,1, &valN) ;//valN =  (C + Q)T * KExtr * (C+Q)
	valN +=  msigm_bmo * msigm_bmo ;
	MtrxMultMatrx(arrKExtr,2, 2,arr_C_Plus_Q,1, arrR) ; // R = KExtr * ( C+Q)
	MatrxDivideScalar(arrR, 2, 1, valN,arrP);  // P = R/N
	if(arrP[0] >1 )
	{
	   arrP[0] = 1 ;
	   arrP[1] = arrKExtr[1] /arrKExtr[0];
	   arrQ[0] = -( msigm_mmo * msigm_mmo + msigm_bmo * msigm_bmo)/ arrKExtr[0] ;
	   arrQ[1] = 0 ;
	   arr_C_Plus_Q[0] =  1 +  arrQ[0] ;
	   arr_C_Plus_Q[1] = 0 ;

	}
	MtrxMultMatrxTransp(arrP,2, 1, arr_C_Plus_Q,2, arrTemp1) ;  // P * (C + Q)T
	MtrxMultMatrx(arrTemp1,2, 2, arrKExtr,2, arrTemp) ;
	MtrxMinusMatrx(arrKExtr, arrTemp,2, 2, arrTemp1);
	arrTemp1[1] = (arrTemp1[1] + arrTemp1[2])/2;
	arrTemp1[2] = arrTemp1[1] ;
	if ( fabs(arrTemp1[0] *arrTemp1[1] - arrTemp1[3] * arrTemp1[1]) <  EPS)
	{
		return false ;
		int jjj = 0;
    }
	if ( (fabs(arrTemp1[1]/arrTemp1[0] *arrTemp1[1]/arrTemp1[3]) >  0.99)
	 || ( arrTemp1[0] < EPS) || ( arrTemp1[3] < EPS))
	 // if ( (fabs(arrTemp1[0] *arrTemp1[1] - arrTemp1[3] * arrTemp1[1]) <  EPS)
	  //	 || ( arrTemp1[0] < EPS) || ( arrTemp1[3] < EPS))

	{
	memset(arrKOut,0,4*sizeof(double)) ;
	return false ;

	}

	memcpy(arrKOut,arrTemp1,4 * sizeof(double)) ;
	return true ;

}

// нахождение гарантированной точнолсти фильтрации для заданного фильтра, являющегося
//  ответом на возмущения с углами формирований fi1, fi2
//Input:
// arrK - нач корреляц матрица для фильтра
// arrKTol - нач корреляц матрица для для расчета функции
//  fi1,fi2 - углы формирования возмущений  на которые настроен фильтр
// Output :
// arrPStab - коефф усилнгия в установившемся режиме
// psi1,psi2 - углы формирования наихудших возмущений
// arrKTol - корреляционная функция(функция максимума)
 bool TStabSyst2::CalcFuncMaximum(double *arrK, double *arrPStab,const double fi1,const double fi2
	   , double &psi1, double &psi2,double *arrKTol )
{
  // const int iNet1 = 50;
  // const int iNet2 = 1;
   double arrK_t[4] ={0},arrKTol_t[4] ={0},arrKTolRez[4] = {0} ;
   double psi1_t = 0, psi2_t = 0 ;
   for (int i = 0; i < INET1; i++)
   {

   for (int j = 0; j < INET2; j++)
   {

	 memcpy(arrK_t,arrK,4 * sizeof(double)) ;
	 memcpy(arrKTol_t,arrKTol,4 * sizeof(double)) ;
	 psi1_t = 2 * PI * i/ INET1 ;
	 psi2_t = 2 * PI * j/ INET1 ;
	 CalcAccuracy(arrK_t, arrPStab, fi1, fi2
	   , psi1_t, psi2_t,arrKTol_t );
	 if (arrKTol_t[3] >  arrKTolRez[3])
	 {
	   memcpy(arrKTolRez, arrKTol_t,4 * sizeof(double)) ;
	   psi1 = psi1_t;
	   psi2 = psi2_t;
	 }
   }
   }

   memcpy(arrK, arrK_t,4 * sizeof(double)) ;
   memcpy(arrKTol, arrKTolRez,4 * sizeof(double)) ;
   return true;
}
//Отыскание минимаксного фильтра путем поиска минимума функции максимума
// Output:
// arrPStab  -  коефф усиления в установившемся режиме
// fi1,fi2  - углы фрормирования возмущений на корорые настроен фильтр
// psi1,psi2 -  углы формирования наихудших для оптимаольного фильтра возмущений
// - должны совпадать с  fi1,fi2
// arrKMinMax -  коррел матрица в установившемся режиме
//
 bool TStabSyst2 ::FindMinMax(double *arrPStab, double &fi1, double & fi2
	   , double &psi1, double &psi2,double *arrKMinMax )
{
   double arrK[4] = {10000,0,0,10000} ;// нач значений коррел матрицы
  // const int iNet1 = 50;
 // const int iNet2 = 1;  // для системы без ММО
   double fi1_t = 0, fi2_t = 0 ;
   double psi1_t = 0, psi2_t = 0 ;
   arrKMinMax[0] = 10000000;
   arrKMinMax[1] = 0;
   arrKMinMax[2] = 0;
   arrKMinMax[3] = 10000000;
   double arrK_t[4] ={0};
   double arrKMinMax_t[4] = {0} ;
   double arrP_t[2] = {0} ;
   memcpy( arrK_t,arrK, 4 * sizeof(double)) ;

   for (int i = 0; i < INET1; i++)
   {

   for (int j = 0; j < INET2; j++)
   {
	fi1_t = 2 * PI * i/ INET1 ;
	fi2_t = 2 * PI * j/ INET2 ;

	memcpy(arrK_t, arrK,4 * sizeof(double)) ;
	memcpy(arrKMinMax_t, arrK,4 * sizeof(double)) ;
	CalcFuncMaximum( arrK_t, arrP_t, fi1_t, fi2_t
	   , psi1_t, psi2_t,arrKMinMax_t ) ;
	if (arrKMinMax_t [3] < arrKMinMax[3])
	{
	  memcpy(arrKMinMax, arrKMinMax_t, 4 * sizeof(double)) ;
	  memcpy( arrPStab, arrP_t, 2 * sizeof(double)) ;
	  psi1 = psi1_t;
	  psi2 = psi2_t;
	  fi1 = fi1_t ;
	  fi2 = fi2_t ;
	}

   }
   }

  return true ;
}

bool TStabSyst2::FindRoots( double *arrP
		   , double *arrS,double *arrQ, TComp &rt1,TComp &rt2)
{
  double arrE[4] = {1,0,0,1} ;
  double arrA [4] = {1,0,0,1};
  arrA [1] = mh ;
  double arrB[2];
   arrB[0] = mh*mh/2 ;
   arrB[1] = mh;
   double arrC[2] = {1,0} ;
   double arrT1[4]={0},arrT2[4] ={0},arrT0[4] = {0},arrT3[4]={0} ;
   double arrT4[4]={0},arrT5[4] ={0},arrT6[4] = {0},arrT7[4]={0} ;
   MtrxMultMatrxTransp(arrP,2, 1, arrC,2, arrT0) ;// arrT0  = P * CT
   MtrxMinusMatrx(arrE, arrT0,2, 2, arrT1); //arrT1 = E - P*CT
   MtrxMultMatrxTransp(arrP,2, 1, arrQ,2, arrT2) ; // arrT2 = P * QT
   MtrxMinusMatrx(arrT1 , arrT2,2,2, arrT3);      // arrT3 = E - P*CT - P * QT
   //
   MtrxMultMatrxTransp(arrB,2, 1, arrS,2, arrT4) ;// arrT4  = B * ST
   MtrxSumMatrx(arrA, arrT4,2, 2, arrT5); //arrT3 = A  + B * ST
   //
   MtrxMultMatrx(arrT3,2, 2, arrT5,2, arrT6) ; // arrT5 = ( E - PCT - P * QT ) *  (A  + B * ST)
   //


   double a =1;
   double b = -(arrT6 [0] + arrT6 [3]);
   double c =  arrT6 [0] * arrT6 [3] - arrT6 [1] * arrT6 [2] ;
   SolvEq2( a, b, c,rt1,rt2);
   if ( (rt1.modul()> 0.9999999999999) || (rt1.modul() > 0.9999999999999))return false;

   return true ;
}
// поиск максимина функции выигрыша (риска)
// на сетке перебираются углы формирования возмущений
// Для каждой пары возмущений ищется оптимальный фильтр
// Выбирается такакя пара формирования возмущений, котороая дает максиму
// дисперсии ошибки фильтрации скорости в сктановившемся режиме
void TStabSyst2::FindMaxMin ( double *arrP
	   ,double &fiRez1,double &fiRez2,double *arrS,double *arrQ)

{
	double arrK0[4] = {10000,0,0,10000} ; // Нач значений коррел матрицы
	double arrK[4]={0} ;

	const int Nit = 100; // количество разбиений по каждому из углов формирования

	double arrSt[2] ={0}, arrQt[2] ={0} ; //векторы формирования возмущений
	int iRez,jRez;
	double arrKrez[4] = {0},arrPrez[2] = {0},arrPt[2] = {0};

	 double fi1,fi2 ;
	 int i0,j0;
	for (int i = 0; i < Nit; i++)
	for (int j = 0; j < Nit; j++)
	{
		fi1 = 2 * PI* i /  Nit;
		fi2 = 2 * PI* j /  Nit;

		memcpy( arrK,arrK0,4 * sizeof(double)) ;
		if (!CalcFuncMinimum(fi1,fi2,arrK,arrPt,arrSt,arrQt)) continue ;
		if ( arrK[3] > arrKrez[3] )
		{
		memcpy( arrKrez,arrK,4 * sizeof(double)) ;
		memcpy( arrP,arrPt,2 * sizeof(double)) ;
		memcpy( arrS,arrSt,2 * sizeof(double)) ;
		memcpy( arrQ,arrQt,2 * sizeof(double)) ;
		fiRez1 = fi1 ;
		fiRez2 = fi2 ;
		i0 = i;
		j0= j;

		}
	}
   memcpy( arrK0,arrKrez,4 * sizeof(double)) ;

}

// расчет корреляционной матрицы ошибок для фильтра с углами формирования
// возмущений fi1 и fi2, когда возмущения формируются с углами psi1 и psi2
// arrK0 - на входе априорная корреляц матрица для фильтра ,
// на выходе корреляционная матрица ошибок фильтрации в установившемся режиме с фильтра
// arrPStab - коэффициенты усиления фиольтра в установивишемся режиме
// fi1 и fi2 - углы формирования процессов на которые настроен фильтр
// psi1 и psi2 - углы формирования процессов
// arrKTol - на входе априорная корреляционная матрицв, на выходе
// корреляционная матрица ошибок фильтрации в установившемс ярежиме
//иначе говоря, это расчет функции риска - т к фильтр фиксирован и возмущения фиксированы

bool TStabSyst2::CalcAccuracy(double *arrK, double *arrPStab,const double fi1,const double fi2
	   ,const double psi1,const double psi2,double *arrKTol )
{
   double arrK_t[4] = {0} ;
   double arrKTol_t[4] = {0} ;
   double arrS_Flt[2] = {0}, arrQ_Flt[2] = {0} ;
   double arrS_Tol[2] ={0}, arrQ_Tol[2] ={0} ;
   int i = 0;
   for ( i = 0; i < 10000; i++)
   {

	 // персчет фильтра

	 OneStepAlfaBettaFlt( fi1, fi2, arrK,arrK_t, arrPStab
		,arrS_Flt,arrQ_Flt) ;
	 // пересчет корреляц матрицы точности
		 if ( ! OneStepRecalcTol(  psi1,psi2, msigm_w,msigm_mmo,msigm_bmo
		 ,arrPStab,arrKTol,arrKTol_t,arrS_Tol,arrQ_Tol) )
		{
		 memset( arrKTol,0,4 * sizeof(double)) ;
		 return false;
		}
	// проверка окончания процесса по невязке
	if ((fabs( arrK[0] - arrK_t[0]) < EPS) && (fabs( arrK[3] - arrK_t[3]) < EPS)
	  &&(fabs( arrKTol[0] - arrKTol_t[0]) < EPS/10) && (fabs( arrKTol[3] - arrKTol_t[3]) < EPS/10) )
	{
	  memcpy( arrK,arrK_t,4 * sizeof(double)) ;
	  memcpy( arrKTol,arrKTol_t,4 * sizeof(double)) ;
	  break ;
	}
   memcpy( arrK,arrK_t,4 * sizeof(double)) ;
   memcpy( arrKTol,arrKTol_t,4 * sizeof(double)) ;
   }
   return true ;
}
// ММО НЕТ!!!
// нахождение гарантированной точнолсти фильтрации для  фильтра Калмана
// НА ПОДМНОЖЕСТВЕ ВОЗМУЩЕНИЙ сформиорванный из дискретного белого шума
// и динамической системы 1 го порядка из остатков фильтрации
// То есть: u = w + ksi
// Mu*u = SigU * SigU
// w(n+1 ) = mu * w (n) + s1 * x1 + s2 * x2;
//Input:
// sigm_u,  sigm_bmo   - реально действующие СКЗ
// Output :
// arrPStab - коефф усилнгия в установившемся режиме
// psi- угол формирования возмущения
// (1 - valLamb *valLamb) - CСКЗ ksi
//  valMu - параметр динамической системы формирующей w
// arrKTol - корреляционная функция(функция максимума)
// pparrK - двойной указатель на массив корреляционой матрицы по времени
// при наихудшем формированиия возмущений
// сначала идет массив K[0], затем К[1], затем K[3]
// plenparrarrK - общая длина массива parrK
//  Если pparrK == NULL, то заполнение массива *pparrK не производится
 bool TStabSyst2 ::F_MaxForKlmFilt( const double sigm_u
	  ,const double sigm_bmo ,double *arrK, double *arrPStab
	   , double *valPsi,double *valLamb,double * valMu,
	   double *arrKTol , double **pparrK,int * plenparrarrK)
{
   arrK [0] = 100000;
   arrK [1] = 0;
   arrK [2] = 0;
   arrK [3] = 100000;
   memcpy(arrKTol,arrK,4 * sizeof(double)) ;
   const int NPsi = 50;// шаг по psi
   const int NLam = 40;// шаг по Lambda
   const int NMu = 40;// шаг по mu
   double arrK_t[4] ={0},arrKTol_t[4] ={0},arrKTolRez[4] = {0} ;
   double psi1_t = 0, psi2_t = 0 ;
   double LambCur = 0 ;
   double MuCur = 0 ;
   double  PsiCur = 0 ;
   // Расчет когда все возмущение объекта белый шум
   Klm_Real_WhiteNoiseInput(sigm_u, sigm_bmo,arrK, arrKTol_t ) ;
   ///


	  LambCur = 1;
	  double sig_ksi = sqrt(1 - LambCur * LambCur)* sigm_u ;

   for (int i2 = 0; i2 < NMu-1; i2++)
   {
	   MuCur = -1.0 + 2.0 * (i2 + 1)/NMu;
	   for (int i3 = 0; i3 < NPsi; i3++)
	   {
		   PsiCur =2 * PI * i3/ NPsi ;

	 memcpy(arrK_t,arrK,4 * sizeof(double)) ;
	// memcpy(arrKTol_t,arrKTol,4 * sizeof(double)) ;

	 CalcToleranceKalmFilt( sigm_u, sigm_bmo, arrK_t, arrPStab
	   ,LambCur,MuCur, PsiCur, arrKTol_t , NULL,NULL);
	 if (arrKTol_t[3] >  arrKTolRez[3])
	 {
	   memcpy(arrKTolRez, arrKTol_t,4 * sizeof(double)) ;
	   *valPsi = PsiCur;
	   *valLamb = LambCur;
	   *valMu =  MuCur ;
	 }
   }
   }
  //
   memcpy(arrK_t,arrK,4 * sizeof(double)) ;
   memcpy(arrKTol_t,arrKTol,4 * sizeof(double)) ;
   CalcToleranceKalmFilt( sigm_u, sigm_bmo, arrK_t, arrPStab
	   ,LambCur,MuCur, PsiCur, arrKTol_t , pparrK, plenparrarrK);
   memcpy(arrK, arrK_t,4 * sizeof(double)) ;
   memcpy(arrKTol, arrKTolRez,4 * sizeof(double)) ;

  return true;

}
 bool TStabSyst2 ::CalcToleranceKalmFilt( const double sigm_u, const double sigm_bmo  ,double *arrK, double *arrPStab
		 ,const double LambCur,const double MuCur, const double PsiCur
	  ,double*  arrKTol , double **pparrK,int *plenparrK)

{


  // присваивание начальных значений матрице
  memcpy( arrKTol,arrK,4 * sizeof(double)) ;
  double arrP[2] ;
  // расширеннляц матроиуаая корред
  double arrKExtend [9] = {0} ;
  arrKExtend [0] = arrK [0] ;
  arrKExtend [1] = arrK [1] ;
  arrKExtend [2] = 0;
  arrKExtend [3] = arrK [1] ;
  arrKExtend [4] = arrK [3] ;
  arrKExtend [5] =  0;
  arrKExtend [6] =  0;
  arrKExtend [7] =  0;
  arrKExtend [8] =  LambCur * LambCur* sigm_u * sigm_u ;
  double arrKExtend_t[9] ={0} ;
  memcpy(arrKExtend_t,arrKExtend, 9* sizeof(double)) ;

	int lenMemAllocated = 0 ;
   if (pparrK != NULL) 	 lenMemAllocated = *plenparrK;

   double arrK_t[4] = {0} ;
   double arrKTol_t[9] = {0} ;

   double arrS_Tol[2] ={0} ;
   int i = 0;
   for ( i = 0; i < 10000; i++)
   {

		OneStepFltrKlm_Real( arrK,arrK_t, arrP) ;
	 // пересчет корреляц матрицы точности

	 if ( !CorrMtrxExtendedRecalc( sigm_u,sigm_bmo  ,arrP
		 , LambCur, MuCur,  PsiCur ,arrKExtend_t , arrS_Tol) )
		{
		 memset( arrKTol,0,4 * sizeof(double)) ;
		 return false;
		}
	// проверка окончания процесса по невязке
	if ((fabs( arrK[0] - arrK_t[0]) < EPS) && (fabs( arrK[3] - arrK_t[3]) < EPS)
	  &&(fabs( arrKExtend[0] - arrKExtend_t[0]) < EPS) && (fabs( arrKExtend[4] - arrKExtend_t[4]) < EPS) )
	{
	  memcpy( arrK,arrK_t,4 * sizeof(double)) ;
	  memcpy( arrKExtend,arrKExtend_t,9 * sizeof(double)) ;
	  break ;
	}
   memcpy( arrK,arrK_t,4 * sizeof(double)) ;
   memcpy( arrKExtend,arrKExtend_t,9 * sizeof(double)) ;
	if (pparrK != NULL)
	{




			if (i * 3 > lenMemAllocated -3)
			{
			lenMemAllocated += 1000;
			*pparrK =  (double *)realloc(*pparrK,lenMemAllocated * sizeof(double)) ;
			}

			(*pparrK)[i * 3     ] = arrKExtend[0] ;
			(*pparrK)[i * 3 + 1 ] = arrKExtend[1] ;
			(*pparrK)[i * 3 + 2 ] = arrKExtend[4] ;

	}

   }

   arrKTol[0] = arrKExtend_t [0];
   arrKTol[1] = arrKExtend_t [1];
   arrKTol[2] = arrKExtend_t [1];
   arrKTol[3] = arrKExtend_t [4];
	if (pparrK != NULL)
	{
	 *pparrK =  (double *)realloc(*pparrK,i * 3 * sizeof(double)) ;
	*plenparrK = i * 3 ;
	double *parrTemp = new double [ *plenparrK];
	MatrTransp(*pparrK, i, 3, parrTemp);
	memcpy(*pparrK, parrTemp, i * 3 * sizeof(double)) ;
	delete parrTemp ;
	}


   return true ;


}

 bool TStabSyst2 ::CorrMtrxExtendedRecalc(const double sigm_u,const double sigm_bmo , double *arrP
		 , const double LambCur, const double MuCur,  const double FiCur ,double *arrKExtend ,double * arrS)
{

  const double sigm_w = LambCur * sigm_u ;
  const double sigm_ksi = sqrt(1 - LambCur * LambCur ) * sigm_u ;
  // 1-ый спец случай
  if (fabs(MuCur)< EPS)
  {
	double arrKInp[4] ;
	arrKInp[0] = arrKExtend[0] ;
	arrKInp[1] = arrKExtend[1] ;
	arrKInp[2] = arrKExtend[1] ;
	arrKInp[3] = arrKExtend[4] ;
	 double arrKExtr[4] = {0} ;
	double arrE[4] = {1,0,0,1} ;

   double arrT1[4]={0},arrT2[4] ={0},arrT0[4] = {0},arrT3[4]={0} ;
   double arrT4[4]={0},arrT5[4] ={0},arrT6[4] = {0}
	   ,arrT7[4]={0},arrT8[4]={0}, arrT9[4]={0},arrT10[4] ={0},arrT11[4] = {0}
	   ,arrT12[4] = {0},arrT13[4] = {0},arrT14[4] = {0},arrT15[4] = {0},arrT16[4] = {0};
   // Экстраполяция

	MtrxMultMatrxTransp(marrB,2, 1, arrS, 2, arrT1) ;  // arrT1 = B * ST
	MtrxSumMatrx(marrA, arrT1,2, 2, arrT2); // arrT2 = A + B*(ST)
	// персчет K(n+1|n)
	MtrxMultMatrx(arrT2,2, 2, arrKInp,2, arrT3) ;// arrT3 = (A + B*(ST)) * K
	MtrxMultMatrxTransp(arrT3,2, 2, arrT2,2, arrKExtr) ;// arrKExtr = (A + B*(ST)) * K *  (A + B*(ST))T
	if ( (arrKExtr[0]< EPS) || (arrKExtr[3]< EPS))
	{
	  // SmhowMessage(L"(arrKExtr[0]< EPS) || (arrKExtr[3]< EPS)" );
	  return false ;
	}
		MtrxMultMatrxTransp (arrP,2, 1, marrC,2, arrT4) ;// arrT4  = P * CT
		MtrxMinusMatrx(arrE, arrT4,2, 2, arrT5); //arrT5 = E - P*CT

		MtrxMultMatrxTransp(arrP,2, 1, arrP,2, arrT8) ; // arrT8 = P * PT
		// персчет матрицы персчет K(n+1)
		MtrxMultMatrx(arrT5,2, 2, arrKExtr,2, arrT9) ;// arrT9 = (E - P*CT ) * K (n+1|n)
		MtrxMultMatrxTransp(arrT9,2, 2, arrT5,2, arrT10) ;// arrT10 = (E - P*CT) * K (n+1|n) *(E - P*CT )T
		MatrxMultScalar(arrT8 , 2,2, sigm_bmo * sigm_bmo,arrT11); // arrT11 = sigm_bmo * sigm_bmo m* P * PT
		MtrxSumMatrx(arrT10, arrT11,2,2, arrT12) ; // arrT12 =  (E - P*CT) * K (n+1|n) *(E - P*CT )T + sigm_bmo * sigm_bmo m* P * PT
		MtrxMultMatrx(arrT5,2, 2, marrB,2, arrT13) ;  // arrT13 =  (E - P*CT) * B
		MtrxMultMatrxTransp(arrT13,2, 1, arrT13,2, arrT14) ; // arrT14 =  (E - P*CT) * B * BT * (E - P*CT)T
		MatrxMultScalar(arrT14 , 2,2, sigm_ksi * sigm_ksi,arrT15);// arrT15 = sigm_ksi * sigm_ksi* (E - P*CT) * B * BT * (E - P*CT)T
		MtrxSumMatrx(arrT12, arrT15,2,2, arrT16) ;
		arrKExtend[0] = arrT16[0];
		arrKExtend[1] = arrT16[1];
		arrKExtend[3] = arrT16[2];
		arrKExtend[4] = arrT16[3];
		arrKExtend[2] = 0 ;
		arrKExtend[5] = 0;
		arrKExtend[6] = 0 ;
		arrKExtend[7] = 0;
		arrKExtend[8] = sigm_w * sigm_w ;

	return true ;
  }
  /*****************************************************************************************/
  /*****************************************************************************************/
  /*****************************************************************************************/
  // помсчк разложение
 double arrS1[4] = {0},arrS2[4] = {0};
 int quantRoots = 0;
 double a1,b1,c1,a,b,c,d,c3,c2,c4;
 TComp ct0,ct1;
 TComp cx0,cx1,cx2,cx3 ;
 int iEq,iEq1,iEq2,iEq3 ;
 double arrRoots[4] = {0};
 if (fabs (FiCur) < EPS)
 {
   quantRoots = 1 ;
   arrS1[0] = sqrt ( (1 - MuCur * MuCur )/ arrKExtend [ 0]) * sigm_w * FSIGN(cos( FiCur)) ;
   arrS2[0] =  0;

 }
 else
 {
		if (fabs(arrKExtend[2])< EPS)
		{   // K13 == 0

			if (fabs(arrKExtend[5])< EPS)
			{ // K13== 0   K23 == 0
			double arr[4] ={0};
			arr[0] = arrKExtend[0];
			arr[1] =  arrKExtend[1];
			arr[2] = arrKExtend[1];
			arr[3] = arrKExtend[4];
			if ( !CalcCoeff( FiCur,sqrt(1 - MuCur * MuCur) *sigm_w,arr,arrS)) return false;
			quantRoots = 1 ;
			arrS1[0] =  arrS[0] ;
			arrS2[0] =  arrS[1] ;
			}
			else
			{ // K13== 0   K23 != 0
			double valA = sqrt(1 - MuCur * MuCur ) * sigm_w ;
			a = 1;
			d = arrKExtend [ 0] * arrKExtend [ 4] - arrKExtend[1]* arrKExtend [1];
			b = 2 * MuCur *  arrKExtend [5] * arrKExtend [0] * sin(FiCur) * sin(FiCur)/d;
			c = - arrKExtend [ 0]* valA * valA * sin(FiCur) * sin(FiCur)/d;
           // Решение квадраьного улавнения a*x*x + b*x + c =0
// возвращает
// 0 - 2 действительных некраьных корня
// 1 - 2 действительных кратных корня
// 2- имеется по крайней мере один нулевой корень a!= 0, c=0
// 3 - 2 комплексно сопряженных корня
// 4 -  1 действительный корень (a =0)
// 5  - несовместность a=b=0, c!=0
// 6 - тождество )a=b=c=0
			int iEq = SolvEq2( a, b, c,cx0,cx1);
			quantRoots = 0 ;
			switch(iEq)
			{
				case 2:
				case 0:
				arrRoots[0] = cx0.m_Re;
				arrRoots[1] = cx1.m_Re;

				 for (int j  = 0; j < 2; j++)
				 {
					a1 = 1;
					b1 = 2 * arrRoots[j] * arrKExtend[1]/ arrKExtend[0]  ;
				   double c1 = (arrRoots[j] * arrRoots[j] * arrKExtend[4] +
					  2 * arrRoots[j] * MuCur * arrKExtend[5] + (MuCur *MuCur -1) *sigm_w*sigm_w )/ arrKExtend[0]  ;
				   iEq1 = SolvEq2( a1, b1, c1,cx2,cx3);
				   switch(iEq1)
				   {
					   case 2:
					   case 0:
					   arrS1[quantRoots ]     =  cx2.m_Re ;
					   arrS2[quantRoots ]     =  arrRoots[j] ;
					   arrS1[quantRoots  + 1] =  cx3.m_Re ;
					   arrS2[quantRoots  + 1] =  arrRoots[j] ;
					   quantRoots += 2 ;
					   break;
					   case 1:
					   arrS1[quantRoots ]     =  cx2.m_Re ;
					   arrS2[quantRoots ]     =  arrRoots[j] ;
					   quantRoots++ ;
					   break;
					   default:
					   break;
				   }
				 }
				break;
				case 1:
				 a1 = 1;
				 b1 = 2 * cx0.m_Re * arrKExtend[1]/ arrKExtend[0]  ;
				 c1 = (cx0.m_Re* cx0.m_Re * arrKExtend[4] +
				2 * cx0.m_Re * MuCur * arrKExtend[5] + (MuCur *MuCur -1) *sigm_w*sigm_w )/ arrKExtend[0]  ;
				 iEq3 = SolvEq2( a1, b1, c1,cx2,cx3);
				switch (iEq3)
				{
				 case 2:
				 case 0:
				 arrS2[0 ] = cx0.m_Re ;
				 arrS2[1 ] = cx0.m_Re ;
				 arrS1[0 ] = cx2.m_Re ;
				 arrS1[1 ] = cx3.m_Re ;
				 quantRoots = 2 ;
				 break;
				 case 1:
				 arrS2[0 ] = cx0.m_Re ;
				 arrS1[0 ] = cx2.m_Re ;
				 quantRoots = 1 ;
				 break;

				default:
				break    ;
				}
				break;
				default: break ;
            }
			}



	  }
	  else
	  { // K13 != 0  общий случай
	  // рапсчет коэффиц  c4, c3, c2 возвратного уравнения 4-ой степени
	   a = sqrt(1 - MuCur * MuCur ) * sigm_w ;
	   d = sqrt( (arrKExtend [ 0] * arrKExtend [ 4] - arrKExtend[1]* arrKExtend [1] )/ arrKExtend [ 0]) ;
	   c4 =  arrKExtend [ 0] * a * a * a * a/MuCur /MuCur /4 ;
	   c3 = (arrKExtend [ 0] * arrKExtend [ 5] - arrKExtend[1]* arrKExtend [2] )
		  * a * a * a* sin( FiCur) / MuCur/d ;
	   c2 = (arrKExtend [ 0]* arrKExtend [ 5] *arrKExtend [ 5]* sin( FiCur)* sin( FiCur)/d/d
	  - arrKExtend [ 0]* a * a/MuCur / MuCur / 2 - 2 * arrKExtend [ 1]* arrKExtend [ 5] *arrKExtend [ 2]* sin( FiCur)* sin( FiCur)/d/d
	  + arrKExtend [ 4]* arrKExtend [ 2] *arrKExtend [ 2]* sin( FiCur)* sin( FiCur)/d/d - arrKExtend [ 2] *arrKExtend [ 2]) * a * a ;

	  if (fabs(c2) > DBL_MAX/100)
	  {
         int id = 100000;
	  }
	  // решение квадратного уравнения относительно t = x - 1/x

	  double arrX[4] = {0};
	  int quantX = 0;
	    iEq2 = SolvEq2( c4, c3, c2 + 2 * c4,ct0,ct1);
	   switch(iEq2)
	   {
		   case 2:
		   case 0:
		   quantX = 4 ;
		   arrX[0] = (ct0.m_Re  + sqrt ( ct0.m_Re * ct0.m_Re + 4 ))/2 ;
		   arrX[1] = (ct0.m_Re  - sqrt ( ct0.m_Re * ct0.m_Re + 4 ))/2 ;
		   arrX[2] = (ct1.m_Re  + sqrt ( ct1.m_Re * ct1.m_Re + 4 ))/2 ;
		   arrX[3] = (ct1.m_Re  - sqrt ( ct1.m_Re * ct1.m_Re + 4 ))/2 ;
		   break ;
		   case 1:
		   quantX = 2 ;
		   arrX[0] = (ct0.m_Re  + sqrt ( ct0.m_Re * ct0.m_Re + 4 ))/2 ;
		   arrX[1] = (ct0.m_Re  - sqrt ( ct0.m_Re * ct0.m_Re + 4 ))/2 ;
		   break;
		   default:
		   break;

	   }
	   for (int j =0;j < quantX; j++)
	   {
		 arrS2[ j ] = arrX[j] * a * sin(FiCur)/d ;
		 arrS1[ j ] = (a * a /MuCur / 2 - d * d * arrS2[ j ] * arrS2[ j ] /sin(FiCur)/sin(FiCur)/MuCur / 2
		   - arrKExtend [ 5] * arrS2[ j ]) /arrKExtend [ 2] ;
	   }
	   quantRoots = quantX ;
	}

}
  	//отсев построниих корней
	//arrS2Good , arrS1Good - хороштие корни, quantGoodRoots  - их количествовщгиду
	 double arrS1Good[4] ={0},arrS2Good [4] ={0};
	 int quantGoodRoots = 0 ;
	 for (int j =0; j < quantRoots; j++)
	 {
	   if (IsRootGood( sigm_w,
		   LambCur,  MuCur,   FiCur ,arrKExtend,arrS1[j],arrS2[j] ))
	   {
		arrS1Good[quantGoodRoots] = arrS1[j] ;
		arrS2Good[quantGoodRoots] = arrS2[j] ;
		quantGoodRoots++ ;
	   }
	 }

//***************************************************************************
//*** В масиивах arrS1Good и arrS2Good длиной quantGoodRoots сохранены все возможные пары S1,S2************************************************************************
//****Для каждоой пары S1,S2 вычисляем корреляц матирицу и выбираем с наибольшим коэфф K22 ( или K11)***********************************************************************
 if (quantGoodRoots == 0)
 {
   //	ShowMessage(L"quantGoodRoots == 0") ;
	quantGoodRoots = 1;
	arrS1Good[0] = 0;
	arrS2Good[0] = 0;
 }
 double arrKt[9] ={0},arrKRez[9] = {0} ;
 for (int i = 0; i < quantGoodRoots; i++)
 {
	 RecalcExtendedMtrx(arrS1Good[i],arrS2Good[i], sigm_u, sigm_bmo , arrP
		 , LambCur,  MuCur,  FiCur ,arrKExtend ,arrKt) ;

		 if (arrKt[4] > arrKRez[4])
		 {
		   arrS[0] = arrS1Good[i]  ;
		   arrS[1] = arrS2Good[i]  ;
		   memcpy(arrKRez,arrKt,9 * sizeof(double)) ;
		 }
 }


 memcpy(arrKExtend,arrKRez,9 * sizeof(double)) ;
   arrKExtend[8] = sigm_w * sigm_w ;
  return true ;
 }
void TStabSyst2 ::RecalcExtendedMtrx(const double valS1,const double valS2,const double sigm_u,const double sigm_bmo , double *arrP
		 , const double LambCur, const double MuCur,  const double FiCur ,double *arrKExtend ,double *arrKExtendOut)
{
  double arrL[9] = {0},arrD[3] = {0},arrP3[3] = {0}
	  ,arrT0[4] = {0},arrT1[4] = {0},arrT2[4] = {0},arrT3[4] = {0},arrT4[4] = {0},arrT5[4] = {0},arrT6[4] = {0},arrT7[4] = {0};
  double arrE[4] = {1,0,0,1} ;
  double arrValS[2];
  arrValS[0] = valS1;
  arrValS[1] = valS2;
  MtrxMultMatrxTransp(arrP,2, 1, marrC,2, arrT0) ; //arrT0 = P * CT
  MtrxMinusMatrx(arrE, arrT0,2, 2, arrT1);   //arrT1 = E -P * CT
  MtrxMultMatrxTransp(marrB,2, 1,arrValS,2, arrT2) ; // arrT2 = B * ST
  MtrxSumMatrx(arrT2, marrA,2,2, arrT3) ; // arrT3 = A +  B * ST
  MtrxMultMatrx(arrT1,2,2, arrT3, 2, arrT4) ; //arrT4 = (E -P * CT) * (A +  B * ST)

  MtrxMultMatrx(arrT1,2,2, marrB, 1, arrT5) ; //arrT5 = (E -P * CT) * B
  MatrxMultScalar(arrT5, 2,1, MuCur,arrT6);  // arrT6 = MuCur * (E -P * CT) * B
  // формирование L
  arrL[0] = arrT4 [0];
  arrL[1] = arrT4 [1];
  arrL[3] = arrT4 [2];
  arrL[4] = arrT4 [3];
  arrL[2] = arrT6 [0] ;
  arrL[5] = arrT6 [1] ;
  arrL[6] = valS1;
  arrL[7] = valS2 ;
  arrL[8] = MuCur ;
  ///

  arrD[0] = arrT5 [0];
  arrD[1] = arrT5 [1];
  arrP3[0] = arrP [0] ;
  arrP3[1] = arrP [1] ;

  // пересчет
  double arr0[9] = {0} ,arr1[9] = {0},arr2[9] = {0},arr3[9] = {0},arr4[9] = {0}
		,arr5[9] = {0},arr6[9] = {0},arr7[9] = {0},arr8[9] = {0} ;
  MtrxMultMatrx(arrL, 3, 3, arrKExtend, 3, arr0) ; // arr0 = L * K
  MtrxMultMatrxTransp(arr0,3, 3,arrL,3, arr1) ; //arr1 =  L * K * LT
  MtrxMultMatrxTransp(arrD,3, 1,arrD,3, arr2) ; // arr2 = D * DT
  MatrxMultScalar(arr2, 3,3, (1 - LambCur * LambCur) *sigm_u * sigm_u ,arr3); // arr3 = (1 - LambCur * LambCur)*sigm_u * sigm_u * D * DT
  MtrxMultMatrxTransp(arrP3,3, 1,arrP3,3, arr4) ; // arr4 = P * PT
  MatrxMultScalar(arr4, 3,3, sigm_bmo * sigm_bmo ,arr5); // arr5 = sigm_bmo * sigm_bmo * P * PT
  MtrxSumMatrx(arr1, arr3,3,3, arr6) ; //arr6 =
  MtrxSumMatrx(arr6, arr5,3,3, arrKExtendOut) ;

}
bool TStabSyst2 ::IsRootGood(const double sigm_w
		 , const double LambCur, const double MuCur,  const double FiCur ,double *arrKExtend,const double s0, const double s1 )
{
	double f1 = s0 * s0 * arrKExtend[0] + s1 * s1* arrKExtend[4]
	  + 2 * s0 * s1 * arrKExtend[1] + 2 * s0 * MuCur * arrKExtend[2] + 2 * s1 * MuCur * arrKExtend[5]
	  + ( MuCur * MuCur - 1)* sigm_w * sigm_w;
	  double s3 = s0;
	  if (fabs(s1)> fabs( s0)) s3 = s1;
	   double z0 = s0/s3;
	   double z1 = s1/s3;

	double f2 = (z0 * arrKExtend[0] + z1 * arrKExtend[1]) / sqrt ( arrKExtend[0])
				  / sqrt( z0  * arrKExtend[0]* z0 +  2 * z0*z1 * arrKExtend[1] + z1* arrKExtend[4]*z1 ) - cos (FiCur);
   int ii = 0;
	if ((fabs (f1) < 0.001) && (fabs (f2) < 0.001)) { return true;}
	else { return false ; }
  	return true;

}
int TStabSyst2 ::FSIGN(const double x)
{
	if (fabs(x) < EPS) return 0;


	return  (x > 0)?1:-1;
}
 /////////////////////////////////////////////////////////////




 // расчет корреляционной матрицы ошибок для фильтра Kaлмана настроенного на белые шумы
 // объекта сформированный непрерывным ьелым шумом, когда на вход подаются дискретные белые шумы
// INPUT:
// sigm_wksi - дисперсия дискретного белолго шума на входе на объекте
// sigm_bmo - дисперсия дискретного белолго шума на входе  в измерениях
// OUTPUT:
//arrKRez - - коррел матрица реальнызх ошибок
bool TStabSyst2 ::Klm_Real_WhiteNoiseInput(const double sigm_wksi,const double sigm_bmo
	  ,double *arrK ,double *arrKRez )
{

   double arrP[2] = {0} ;
   double arrK1[4]= {1000000,0,0,100000};
   memcpy(arrK,arrK1,4 * sizeof(double)) ;
   double arrK_t [4] ={0} ;
   double arrKRez1[4]={0} ;
   memcpy(arrKRez,arrK,4 * sizeof(double)) ;

	double arrT0[4] = {0},arrT1[4] = {0},arrT2[4] = {0},arrT3[4] = {0},arrT4[4] = {0},arrT5[4] = {0},arrT6[4] = {0},arrT7[4] = {0};
  double arrE[4] = {1,0,0,1} ;


   int i = 0;
   for ( i = 0; i < 10000; i++)
   {

		OneStepFltrKlm_Real( arrK,arrK_t, arrP) ;
	 // пересчет корреляц матрицы точности

  MtrxMultMatrxTransp(arrP,2, 1, marrC,2, arrT0) ; //arrT0 = P * CT
  MtrxMinusMatrx(arrE, arrT0,2, 2, arrT1);   //arrT1 = E -P * CT


  MtrxMultMatrx(arrT1,2,2, marrA, 2, arrT4) ; //arrT4 = (E -P * CT) * A

  MtrxMultMatrx(arrT1,2,2, marrB, 1, arrT5) ; //arrT5 = (E -P * CT) * B
 // MatrxMultScalar(arrT5, 2,1, MuCur,arrT6);  // arrT6 = MuCur * (E -P * CT) * B


  double arr0[9] = {0} ,arr1[9] = {0},arr2[9] = {0},arr3[9] = {0},arr4[9] = {0}
		,arr5[9] = {0},arr6[9] = {0},arr7[9] = {0},arr8[9] = {0} ;
  MtrxMultMatrx(arrT4, 2, 2, arrKRez, 2, arr0) ; // arr0 = (E -P * CT) * A * K
	MtrxMultMatrxTransp(arr0,2, 2,arrT4,2, arr1) ; //arr1 =  (E -P * CT) * A *K* AT* (E -P * CT)T

	MtrxMultMatrxTransp(arrT5,2, 1, arrT5,2, arr2) ; //arr2 = (E -P * CT) * B*BT* (E -P * CT)T
	MatrxMultScalar(arr2, 2,2,sigm_wksi*sigm_wksi,arrT6);  // arrT6 = sigm_wksi*sigm_wksi *(E -P * CT) * B*BT* (E -P * CT)T

	MtrxMultMatrxTransp(arrP,2, 1,arrP,2, arr4) ; // arr4 = P * PT
	MatrxMultScalar(arr4, 2,2,sigm_bmo*sigm_bmo,arrT5);  // arrT5 = sigm_bmo*sigm_bmo *P * PT


	MtrxSumMatrx(arr1, arrT6,2,2, arr6) ; //arr6 =  (E -P * CT) * A *K* AT* (E -P * CT)T +  sigm_wksi*sigm_wksi *(E -P * CT) * B*BT* (E -P * CT)T
		MtrxSumMatrx(arr6, arrT5,2,2, arrKRez1) ; //arr6 =  (E -P * CT




	// проверка окончания процесса по невязке
	if ((fabs( arrK[0] - arrK_t[0]) < EPS) && (fabs( arrK[3] - arrK_t[3]) < EPS)
	  &&(fabs( arrKRez1[0] - arrKRez[0]) < EPS) && (fabs( arrKRez1[3] - arrKRez[3]) < EPS) )
	{
	  memcpy( arrK,arrK_t,4 * sizeof(double)) ;
	  memcpy( arrKRez,arrKRez1,4 * sizeof(double)) ;
	  break ;
	}
   memcpy( arrK,arrK_t,4 * sizeof(double)) ;
   memcpy( arrKRez,arrKRez1,4 * sizeof(double)) ;

   }

   return true ;
}

// Пгоиск установившегося решения  в случае когда корреляционная матрица ошибок оценивания
// вырождена
bool TStabSyst2 ::CaseOfDegeneracy(double *arrKRez ,int *iarrParamsRez )
{
  bool breturn = false ;
  int iarrParTemp[4] = {0};
 double arrParTemp[4] = {0} ;
 memset(arrKRez,0,sizeof(double)* 4);
 double arrKTemp[4] ={0} ;
 double eps = 0.001;
 arrParTemp[0] = mh ;
 arrParTemp[1] = msigm_w ;
 arrParTemp[2] = msigm_mmo ;
 arrParTemp[3] = msigm_bmo ;
 for (int iCase =0; iCase < 2; iCase++)
 {
   for (int iD1 =0; iD1 < 2; iD1++)
   {
	   for (int iD2 =0; iD2 < 2; iD2++)
	   {
		 for (int iNu = 0; iNu < 2; iNu++)
		 {


		 iarrParTemp[0] = iCase;
		 iarrParTemp[1] = iD1 ;
		 iarrParTemp[2] = iD2 ;
		 iarrParTemp[2] = iNu ;
		 for (int i = 0; i < 4; i++) if (iarrParTemp[i] == 0) iarrParTemp[i] = -1 ;
		 double x = 1000000;
			 if(TangMethod(x,func, pf,iarrParTemp,arrParTemp,eps))
			 {
				CalcMatrK(x,iarrParTemp,arrKTemp) ;

				if (arrKTemp[3] > arrKRez[3])
				{
				memcpy( arrKRez, arrKTemp,4 * sizeof(double)) ;
				memcpy( iarrParamsRez,iarrParTemp,4 * sizeof(int)) ;
				breturn = true;
				}

			 }
		 }
	  }
   }
 }


 return  breturn;
}

void TStabSyst2 ::CalcMatrK(double x,int *iarrPar,double *arrKTemp)
{
 int iMu =  iarrPar[0];
 int iD1 =  iarrPar[1];
 int iD2 =  iarrPar[2];
 int iNu =  iarrPar[3];
 double e2 = x;
 double e1 = iD1*iMu *e2*e2 / msigm_w  - mh * e2/2 ;

	double val = msigm_bmo / sqrt(msigm_bmo * msigm_bmo + (e1 + iarrPar[2]* msigm_mmo)*(e1 + iarrPar[2]* msigm_mmo));
	double k1 =  iNu * val * e1 ;
	double k2 =  iNu * val * e2 ;
	k1 = iD1 * mh * mh* msigm_w/2 + iMu *(e1 - mh * e2);
	k2 = iMu * e2 - iD1 * mh * msigm_w;
	arrKTemp[0] = k1 * k1 ;
	arrKTemp[1] = k1 * k2 ;
	arrKTemp[2] = k1 * k2 ;
	arrKTemp[3] = k2 * k2 ;
}
bool TStabSyst2 ::TangMethod(double &x,double (*f)( double x, int *iarrPar,double *arrPar)
   ,double (*pf)( double x,int *iarrPar,double *arrPar),int *iarrPar
	,double *arrPar,const double eps)
{

 double xt = x ;
 x = x + 10 * eps;
 double df;
 while (fabs(x- xt) > eps)
 {
   x = xt;
  if( fabs((df = (*pf)( xt,iarrPar,arrPar))) < 0.000001) return false;
  xt = xt - (*f)(  xt, iarrPar,arrPar)/df;

 }

return true ;
}

double TStabSyst2 ::func(double x,int *iarrPar,double *arrPar)
{
 double h = arrPar[0]  ;
 double sigm_w =arrPar[1];
 double sigm_mmo=arrPar[2];
 double sigm_bmo =arrPar[3];
 int iMu =  iarrPar[0];
 int iD1 =  iarrPar[1];
 int iD2 =  iarrPar[2];
 int iNu =  iarrPar[3];
 double val0 = (iD1 * iMu * x * x/sigm_w - h * x /2 + iD2 * sigm_mmo );
 double val1 =sqrt(sigm_bmo * sigm_bmo + val0 * val0 ) ;
 return x * iMu -iNu * sigm_bmo *x / val1 - iD1 * sigm_w *h ;

}
double TStabSyst2 ::pf(double x,int *iarrPar,double *arrPar)
{
 double h = arrPar[0]  ;
 double sigm_w =arrPar[1];
 double sigm_mmo=arrPar[2];
 double sigm_bmo =arrPar[3];
 int iMu =  iarrPar[0];
 int iD1 =  iarrPar[1];
 int iD2 =  iarrPar[2];
 int iNu =  iarrPar[3];
 double val0 = (iD1 * iMu * x * x/sigm_w - h * x /2 + iD2 * sigm_mmo );
 double val1 =sqrt(sigm_bmo * sigm_bmo + val0 * val0 ) ;
return iMu - iNu * sigm_bmo * (1/val1 - x *val0 *(2 *iD1 * iMu /sigm_w - h /2  )/val1/val1/val1 )  ;
}

#pragma package(smart_init)
