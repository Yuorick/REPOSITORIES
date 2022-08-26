
#include <math.h>
#include "Comp.h"
#include  <string.h>

//---------------------------------------------------------------------------
//***********************************************************************************************
//******** Класс копмпл чисел для DOUBLE *************************************************
//***********************************************************************************************

TComp::TComp():
 m_Re(0)
,m_Im(0)
{

}
// конструктор копирования
 TComp ::TComp (const TComp &R)
 {
   m_Re = R.m_Re ;
   m_Im = R.m_Im;
 }
 // оператор присваивания
 TComp &TComp::operator=(const TComp  &R)
 {
   m_Re = R.m_Re ;
   m_Im = R.m_Im;
   return *this ;
 }
// парам констр
TComp :: TComp(  const double x,const double y)

 {
   m_Re = x ;
   m_Im = y;
}


TComp TComp:: operator +(const TComp &cmp0)
{
  TComp cmpTemp;
  cmpTemp.m_Re = cmp0.m_Re + m_Re;
  cmpTemp.m_Im = cmp0.m_Im + m_Im;
  return cmpTemp;
}
TComp TComp:: operator -(const TComp &cmp0)
{
  TComp cmpTemp;
  cmpTemp.m_Re = -cmp0.m_Re + m_Re;
  cmpTemp.m_Im = -cmp0.m_Im + m_Im;
  return cmpTemp;
}
 TComp  TComp:: operator *(const TComp &cmp0)
{
  TComp cmpTemp;
  cmpTemp.m_Re = cmp0.m_Re * m_Re - m_Im * cmp0.m_Im;
  cmpTemp.m_Im = m_Im * cmp0.m_Re + m_Re * cmp0.m_Im;
  return cmpTemp;
}

TComp TComp:: operator /(const TComp &cmp0)
{
  TComp cmpTemp;
  double temp = cmp0.m_Re* cmp0.m_Re + cmp0.m_Im * cmp0.m_Im;
  cmpTemp.m_Re = (m_Re * cmp0.m_Re  + m_Im * cmp0.m_Im)/ temp  ;
  cmpTemp.m_Im = (m_Im * cmp0.m_Re - m_Re * cmp0.m_Im )/ temp  ;
  return cmpTemp;
}
const TComp &TComp:: operator +=(const TComp &cmp0)
{
  TComp cmpTemp;
  m_Re += cmp0.m_Re;
  m_Im += cmp0.m_Im ;
  return *this;
}
const TComp &TComp:: operator *=(const TComp &cmp0)
{
  TComp cmpTemp;
  cmpTemp = (*this) * cmp0;
  *this =  cmpTemp;
  return *this;
}

const TComp &TComp:: operator -=(const TComp &cmp0)
{
  TComp cmpTemp;
  m_Re -= cmp0.m_Re;
  m_Im -= cmp0.m_Im ;
  return *this;
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

TComp TComp::scalProd(TComp *cmparrA, TComp *cmparrB, const int LEn)
{
	 TComp cmpRez(0.,0.);
	for (int i = 0; i < LEn; i++)
	{
	  cmpRez += cmparrA[i] * cmparrB[i].Sopr();
	}
	return cmpRez;
}
double TComp::modul()
{
	return sqrt(m_Re* m_Re + m_Im * m_Im);
}

double TComp::phase()
{
   if (fabs(m_Re) < 0.0000000001)
   {
	 if (m_Im > 0)
	 {
	   return M_PI / 2.;
	 }
	 return -M_PI / 2.;
   }
  
	return atan2(m_Im, m_Re);
}

TComp exp_(const TComp cmp)
{
	double r = exp(cmp.m_Re);
	return TComp( r * cos(cmp.m_Im), r * sin(cmp.m_Im) );
}

TComp exp_(const double fi)
{

	return TComp(cos(fi), sin (fi) );
}
TComp TComp::Sqrt_()
{
	double ph =  phase();
	ph = ph/2.;
	double mod = sqrt(modul());
	return TComp( mod * cos(ph), mod * sin(ph));
}

// корень 3 степени   из комплексного числа
// OUTPUT:
// cmparrRoots[3] - массив корней
void TComp::root3( TComp *cmparrRoots)
{
	double ph =  phase() / 3.;
	double arrPh[3] = {0.};
	arrPh[0] = ph ;
	arrPh[1] = ph + 2. * M_PI / 3.;
	arrPh[2] = ph + 4. * M_PI / 3.;
	double mod =  pow(modul(), 1. / 3.);
	for (int i = 0; i < 3; i++)
	{
	 cmparrRoots[i] = TComp( mod * cos(arrPh[i]), mod * sin(arrPh[i]));
	}

}

TComp TComp::Ln()
{
	TComp cmpRez(0.,0.);
	cmpRez.m_Im =  phase();
	cmpRez.m_Re = log(modul());

	return cmpRez;
}
// комплексно сопряженное
TComp TComp::Sopr()
{
	TComp cmpRez(0.,0.);
	cmpRez.m_Im =  -m_Im;
	cmpRez.m_Re = m_Re;

	return cmpRez;
}

// эрмитово сопряженная матрица
void TComp::HermiteSoprMatr(TComp*parrA, int nRows, int nCols, TComp*parrRez)
{

	for (int i = 0; i < nRows; i++)
	for (int j = 0; j < nCols; j++)
	{
	 parrRez[ j * nRows + i] =  parrA[i * nCols + j].Sopr() ;
	}

}

// обоащение матрицы 2-го порядка
// Returns:
// true , если det(A) != 0
// false - в противном случае
bool   TComp::InverseCmpMtrx2( TComp*arrA, TComp*arrOut)
{
   TComp det = arrA [ 0] * arrA [ 3] - arrA [ 1] * arrA [ 2] ;
   TComp cmpMin1(-1.,0.);
   if (det.modul() < 0.000000000001) return false ;
   arrOut [ 0] =  arrA [ 3] / det ;
   arrOut [ 1] =  arrA [ 1] / det * cmpMin1;
   arrOut [ 2] =  arrA [ 2] / det * cmpMin1;
   arrOut [ 3] =  arrA [ 0] / det ;
  return true;
}

//---------------------------------------------------------------------------


// вычисление угла между двумя комплексными числами
// ЕСли одно из чисел == 0 возвращает false
// в пртиивном случае true
// INPUT:
// cmpInp1, cmpInp2 - комплексные числа
// OUTPUT:
// *pvalOut - угол
//
bool TComp::angBetveenComps(TComp cmpInp1, TComp cmpInp2, double *pvalOut)
{
  if (( cmpInp1.modul()< 0.000000001) || ( cmpInp2.modul()< 0.000000001))
  {
	return false;
  }
  TComp cmpt = cmpInp1 * (cmpInp2.Sopr());
  *pvalOut = acos(cmpt.m_Re / cmpInp1.modul()/ cmpInp2.modul());
  return true;
}


//
//---------------------------------------------------------------------------
// умножение матрицы на матрицу
void MtrxMultMatrx(TComp*parrA,int nRowsA, int nColsA, TComp* parrB,int nColsB, TComp*parrRez)
{
  memset (parrRez,0,sizeof(TComp) * nRowsA *  nColsB) ;
  TComp*pRez = parrRez ;
  TComp*pA = parrA ;
  TComp*pB = parrB ;

  for (int i = 0; i < nRowsA; i++)
  {
	TComp*pAt = pA;
   TComp*pBt =pB ;
  for (int j = 0; j < nColsB; j++)
  {
	pAt = pA ;
	pBt = parrB + j ;

	for (int k = 0; k < nColsA ; k++)
	{
	  *pRez += (*pAt) * (*pBt);
	  pAt ++;
	  pBt += nColsB ;
	}
	pRez++;
  }
  pA += nColsA ;
  }

}


// Умножение матрицы на транспонированную матрицу
void MtrxMultMatrxTransp(TComp*parrA,int nRowsA, int nColsA, TComp* parrB,int nRowsB, TComp*parrRez)
{
  memset (parrRez,0,sizeof(TComp) * nRowsA *  nRowsB) ;
  TComp*pRez = parrRez ;
  TComp*pA = parrA ;
  TComp*pB = parrB ;

  for (int i = 0; i < nRowsA; i++)
  {
	TComp*pAt = pA;
   TComp*pBt =pB ;
  for (int j = 0; j < nRowsB; j++)
  {
	pAt = pA ;


	for (int k = 0; k < nColsA ; k++)
	{
	  *pRez += (*pAt) * (*pBt);
	  pAt ++;
	  pBt ++;
	}
	pRez++;
  }
  pA += nColsA ;
  }
}
// Умножение транспонированной матрицы на  матрицу
void MtrxTranspMultMatrx(TComp*parrA,int nRowsA, int nColsA, TComp* parrB,int nColsB, TComp*parrRez)
{
  memset (parrRez,0,sizeof(TComp) * nColsA *  nColsB) ;
 for (int i = 0; i < nColsA ; i++)
  {
	 for (int j = 0; j < nColsB; j++)
	 {
		int ktemp = i* nColsB + j ;
		for (int k = 0; k < nRowsA ; k++)
		{
		  parrRez [ktemp] += parrA[k * nColsA + i]*  parrB [ k *  nColsB + j ] ;
		}
	 }
  }

}
void MtrxSumMatrx(TComp*pA, TComp* pB,int nRows, int nCols, TComp*pRez)
{
	TComp*parrA = pA;
	TComp* parrB = pB ;
	TComp*parrRez = pRez;
	for (int i = 0; i < nRows * nCols; i++)
	{
	  *parrRez = (*parrA)+ (* parrB);
	  parrRez++;
	  parrA++ ;
	  parrB++ ;

	}
}
void MtrxMinusMatrx(TComp*pA, TComp* pB,int nRows, int nCols, TComp*pRez)
{
	TComp*parrA = pA;
	TComp* parrB = pB ;
	TComp*parrRez = pRez;

	for (int i = 0; i < nRows * nCols; i++)
	{
	  *parrRez = (*parrA)- (* parrB);
	  parrRez++;
	  parrA++ ;
	  parrB++ ;

	}
}

void MatrxMultScalar(TComp*parrA, int nRows, int nCols, TComp valScal,TComp*parrRez)
{

   for (int i = 0; i < nRows * nCols; i++)
	{
	  parrRez [i] =  parrA [i] * valScal ;
	}
}
void MatrxDivideScalar(TComp*parrA, int nRows, int nCols, TComp valScal,TComp*parrRez)
{
   for (int i = 0; i < nRows * nCols; i++)
	{
	  *parrRez = (*parrA)/ valScal ;
	  parrRez++;
	  parrA++ ;


	}
}
void MatrTransp(TComp*parrA, int nRows, int nCols, TComp*parrRez)
{
	TComp*pA = parrA ;
	for (int i = 0; i < nRows; i++)
	{

	TComp*pRezt =parrRez + i ;
		for (int j = 0; j < nCols; j++)
		{
		*pRezt = *pA;
		pA ++;
		pRezt += nRows ;
		}

	}
}



void MatrRe(TComp*parrA, int nRows, int nCols, double *parrRez)
{
	for (int i = 0; i < nRows *nCols ; i++)
	{
	 parrRez [i] =  parrA[i].m_Re ;
	}
}

void MatrIm(TComp*parrA, int nRows, int nCols, double *parrRez)
{
	for (int i = 0; i < nRows *nCols ; i++)
	{
	 parrRez [i] =  parrA[i].m_Im ;
	}
}

void   swap(TComp *a0, TComp *a1)
{
	TComp  b = *a1;
	*a1 = *a0 ;
	*a0 = b ;
}


 // нахождение обратной матрицы  3x3
// если определитель  arrInp != нулю, то True
// в противном случае false
 bool   InverseMtrx3(TComp *arrInp, TComp *arrOut)
{
  // вычисление определителя
  TComp det = arrInp[0] * ( arrInp[4] * arrInp[8] - arrInp[5] * arrInp[7] )
			-  arrInp[1] * ( arrInp[3] * arrInp[8] - arrInp[5] * arrInp[6] )
			+  arrInp[2] * ( arrInp[3] * arrInp[7] - arrInp[4] * arrInp[6] )  ;
  ///
  if (fabs (det.modul()) < 0.0000000000000001) return false ;

  TComp arr0[9], cmpMin1(-1.,0.);

  arrOut[0]=  ( ( arrInp[4] * arrInp[8]) - (arrInp[5] * arrInp[7] ))/ det ;
  arrOut[3]= cmpMin1* (( arrInp[3] * arrInp[8]) - ( arrInp[5] * arrInp[6] ))/ det ;
  arrOut[6]=  (( arrInp[3] * arrInp[7]) - (arrInp[4] * arrInp[6] ))/ det ;
  arrOut[1]=  cmpMin1* (( arrInp[1] * arrInp[8]) -( arrInp[2] * arrInp[7] )) / det;
  arrOut[4]=  (( arrInp[0] * arrInp[8]) - (arrInp[2] * arrInp[6] ))/ det ;
  arrOut[7]=  cmpMin1* (( arrInp[0] * arrInp[7]) -( arrInp[1] * arrInp[6] ))/ det ;
  arrOut[2]=  (( arrInp[1] * arrInp[5]) - (arrInp[2] * arrInp[4] ))/ det ;
  arrOut[5]=  cmpMin1* (( arrInp[0] * arrInp[5]) - (arrInp[2] * arrInp[3] )) / det;
  arrOut[8]=  (( arrInp[0] * arrInp[4]) - (arrInp[1] * arrInp[3] ))/ det ;

  return true ;
}

// формирование дополнительного минора к элементу с нмером  iRowElem, iColElem
// parrMtrxInp - исходная матрица
// ndim - размерность  исходная матрица
// iRowElem , jColElem  - номер элемента
// OUTPUT:
// parrMinor - дополнительный минор размерности  [ndim -1]x [ndim -1]
void buildAddMinor(TComp *parrMtrxInp, const int ndim
	, const int iRowElem,const int jColElem, TComp *parrMinor)
{
	int iMinor = 0 ;
	int jMinor = 0 ;
	for (int i = 0; i < ndim; i++)
	{
	 if (i == iRowElem) continue;

	 for (int j = 0; j < ndim; j++)
	 {
	   if (j == jColElem) continue ;

		parrMinor [ (ndim-1) * iMinor + jMinor] =  parrMtrxInp [ ndim * i + j] ;
	   jMinor++;
	 }
	 iMinor++;
	 jMinor = 0;
	}
}

// определеитель матнрицы 4х4
TComp   calcDet4( TComp *arrInp)
{
	TComp det(0.,0.);
	TComp parrMinor [9]  ;
	memset(parrMinor, 0, 9 * sizeof(TComp));
	TComp valIndex ( 1.,0.), cmpMin1(-1.,0.);
	for (int j = 0; j < 4; j++)
	{
	  buildAddMinor(arrInp, 4, 0, j, parrMinor) ;
	  det += valIndex * arrInp[j] * calcDet3( parrMinor) ;
	  valIndex = cmpMin1 *valIndex;
	}
	return det;
}

//  обращение атрицы 4х4
//
bool   InverseMtrx4(TComp *arrInp, TComp *arrOut)
{
	TComp det =   calcDet4( arrInp);
	if (fabs(det.modul()) < 1.E-15) return false ;
	TComp parrMinor [9]  ;
	memset(parrMinor, 0, 9 * sizeof(TComp));
	TComp valIndexI ( 1., 0.);
	TComp valIndexJ  ( 1., 0.);
	TComp cmpMin1(-1.,0.);
	for (int i = 0; i < 4; i++)
	{
	for (int j = 0; j < 4; j++)
	{
	 buildAddMinor(arrInp, 4, j,i, parrMinor) ;
	 arrOut[ 4 * i + j] =  valIndexI * valIndexJ * calcDet3( parrMinor)/ det;
	 valIndexJ = cmpMin1 * valIndexJ;
	}
	valIndexI = cmpMin1 *valIndexI;
	valIndexJ = TComp(1., 0.);
	}
	return true;
}

// вычисление определителя
TComp   calcDet3( TComp *arrInp)
{
 return arrInp[0] * ( arrInp[4] * arrInp[8] - arrInp[5] * arrInp[7] )
			-  arrInp[1] * ( arrInp[3] * arrInp[8] - arrInp[5] * arrInp[6] )
			+  arrInp[2] * ( arrInp[3] * arrInp[7] - arrInp[4] * arrInp[6] )  ;
}


// создание матрицы ГАнкеля
// INPUT:
//cmparrS - массив компл чисел
// lenS  - длина массива
//OUTPUT:
// cmparrHankel [lenS/2 x lenS/2 ] - матирца Ганкеля
 void createHankel( int lenS, TComp *cmparrS, TComp *cmparrHankel )
{
 for (int i =0; i < lenS/2; i++)
 {
   memcpy(&cmparrHankel[ i *lenS/2], & cmparrS[i],  lenS/2 * sizeof(TComp));
 }
}

// эрмитово сопряженная матрица
void HermiteSoprMatr(TComp*parrA, int nRows, int nCols, TComp*parrRez)
{

	for (int i = 0; i < nRows; i++)
	for (int j = 0; j < nCols; j++)
	{
	 parrRez[ j * nRows + i] =  parrA[i * nCols + j].Sopr() ;
	}

}

// обоащение матрицы 2-го порядка
// Returns:
// true , если det(A) != 0
// false - в противном случае
bool   InverseMtrx2( TComp*arrA, TComp*arrOut)
{
   TComp det = arrA [ 0] * arrA [ 3] - arrA [ 1] * arrA [ 2] ;
   TComp cmpMin1(-1.,0.);
   if (det.modul() < 0.000000000001) return false ;
   arrOut [ 0] =  arrA [ 3] / det ;
   arrOut [ 1] =  arrA [ 1] / det * cmpMin1;
   arrOut [ 2] =  arrA [ 2] / det * cmpMin1;
   arrOut [ 3] =  arrA [ 0] / det ;
  return true;
}

// дро-линейная функциябно
// y = (a*Z + b)/(c*Z + d)
// INPUT:
// cmpZ
// OUTPUT:
// *pcmpRez
// возвращает false если знаменатель равен ную
// true , если не равен
bool fncLinFrac(TComp cmpa, TComp cmpb, TComp cmpc, TComp cmpd, TComp cmpZ, TComp *pcmpRez)
{
	TComp cmpt1 = cmpc* cmpZ + cmpd;
	if (fabs(cmpt1.modul()) < 0.00000000001)
	{
	 return false;
	}
	TComp cmpt0 = cmpa* cmpZ + cmpb;
	*pcmpRez = cmpt0/ cmpt1;
	return true;
}

// Нахождение параметров окружности описанной вокруг 3 точек комплексной плоскости
//INPUT:
// cmp0, cmp1, cmp2  - точки компл плоскости
// OUTPUT:
// pcmpCentre - центр
// pvalRadius - радиус
// возвращает true если точки не лежат на одной прямой
// false, если лежат
bool  findCircleParams (TComp cmp0,TComp  cmp1,TComp  cmp2
	, TComp  *pcmpCentre, double *pvalRadius)
{
 TComp cmpT0 =  cmp1 - cmp0;
 TComp cmpT1 =  cmp2 - cmp0;
 double arrA[4] ={0.};
 arrA[0] = cmpT0.m_Re;
 arrA[1] = cmpT0.m_Im;
 arrA[2] = cmpT1.m_Re;
 arrA[3] = cmpT1.m_Im;

 double arrB[2] = {0.};
 arrB[0] = ( cmp1.modul()* cmp1.modul() - cmp0.modul()* cmp0.modul())/ 2.;
 arrB[1] = ( cmp2.modul()* cmp2.modul() - cmp0.modul()* cmp0.modul())/ 2.;

 double arrX[2] = {0.};
 if( !SolvLinEq2__(arrA, arrB, arrX))
 {
	 return false;
 }
 (*pcmpCentre).m_Re = arrX[0] ;
 (*pcmpCentre).m_Im = arrX[1] ;
  TComp cmpT3 =  (*pcmpCentre) - cmp0;
  *pvalRadius =  cmpT3.modul()   ;
  return true;
}

// решение системы линейных уравнений второго порядка  A * x = B
// Returns:
// true , если det(A) != 0
// false - в противном случае
bool   SolvLinEq2__(double *arrA, double *arrB,double *arrX)
{
   double det = arrA [ 0] * arrA [ 3] - arrA [ 1] * arrA [ 2] ;
   if (fabs(det) < 0.000000001) return false ;
   arrX [ 0] =   ( arrB [ 0] *  arrA [ 3] - arrA [ 1] * arrB [ 1] ) / det ;
   arrX [ 1] = ( arrA [ 0] * arrB [ 1]  - arrA [ 2] * arrB [ 0] ) / det ;

  return true;
}

 // преобразование комплексного числа заданного в тригоном форме в обычную
 TComp transfTrigonForm(double valMod, double valArg)
 {
	 return TComp(valMod * cos(valArg),valMod * sin(valArg));
 }

// скалярное прозведение 2-х векторов
TComp scalProd(TComp *cmparr0,TComp *cmparr1, int lenarr)
{
	TComp cmpRez(0.,0.);
	for (int i = 0; i < lenarr; i++)
	{
	 cmpRez = cmpRez + cmparr0[i] * cmparr1[i].Sopr();
	}
	return cmpRez;
}

// норма вектора
TComp vectNorm(TComp *cmparr0, int lenarr)
{
	double  valRez= 0.;
	for (int i = 0; i < lenarr; i++)
	{
	 valRez +=  cmparr0[i].m_Re * cmparr0[i].m_Re + cmparr0[i].m_Im * cmparr0[i].m_Im;
	}
	return TComp(sqrt(valRez), 0.);
}

// норма вектора
double vectNorm_(TComp *cmparr0, int lenarr)
{
	double  valRez= 0.;
	for (int i = 0; i < lenarr; i++)
	{
	 valRez +=  cmparr0[i].m_Re * cmparr0[i].m_Re + cmparr0[i].m_Im * cmparr0[i].m_Im;
	}
	return sqrt(valRez);
}

// проектирование 3-х мерного вектора на плоскость двух векторов
// INPUT:
// cmparrBasis0[3], cmparrBasis1[3] - базисные векторы, на которые производится проектирование
// cmparrS[3] - вектор, который проектируется
// OUTPUT:
// pcmpAlf0, pcmpAlf1 - разложение проекции вектора cmparrInp по базисным векторам
// pcmparrPerp[3] - перпендикуляр
// возвращает True, если базтсные векторы линейно независимы и false в протиивном случае
bool TComp::projectVect_3Dim(TComp *cmparrBasis0, TComp *cmparrBasis1, TComp *cmparrS
  ,TComp *pcmpAlf0, TComp *pcmpAlf1, TComp *pcmparrPerp)
{
 // формирование мартрицы Грама
  TComp cmparrGram[4];
  cmparrGram[0] =  scalProd(cmparrBasis0, cmparrBasis0, 3);
  cmparrGram[1] =  scalProd(cmparrBasis1, cmparrBasis0, 3);
  cmparrGram[2] =  scalProd(cmparrBasis0, cmparrBasis1, 3);
  cmparrGram[3] =  scalProd(cmparrBasis1, cmparrBasis1, 3);
  ///

  // формирование правой части
  TComp cmparrB[2];
  cmparrB[0] =  scalProd(cmparrS, cmparrBasis0, 3);
  cmparrB[1] =  scalProd(cmparrS, cmparrBasis1, 3);
  ///

  // проверка на вырожденность
  TComp cmparrGramInv[4];
  if( !InverseMtrx2( cmparrGram, cmparrGramInv))
  {
	  *pcmpAlf0 = cmparrB[0]/  vectNorm(cmparrBasis0, 3) ;
	  *pcmpAlf1 = TComp(0.,0.);
	  TComp cmpTemp[3];
	  MatrxMultScalar(cmparrBasis0, 3, 1, *pcmpAlf0,cmpTemp);
	  MtrxMinusMatrx(cmparrS, cmpTemp,1, 3, pcmparrPerp);
	  return false;
  }
  else
  {
	  TComp cmparrAlf[2];
	  MtrxMultMatrx(cmparrGramInv,2, 2, cmparrB,1, cmparrAlf);
	  *pcmpAlf0 =  cmparrAlf[0];
	  *pcmpAlf1 =  cmparrAlf[1];

	  TComp cmpTemp0[3], cmpTemp1[3], cmpTemp2[3];
	  MatrxMultScalar(cmparrBasis0, 3, 1, *pcmpAlf0,cmpTemp0);
	  MatrxMultScalar(cmparrBasis1, 3, 1, *pcmpAlf1,cmpTemp1);
	  MtrxSumMatrx(cmpTemp0, cmpTemp1,1, 3, cmpTemp2);
	  MtrxMinusMatrx(cmparrS, cmpTemp2,1, 3, pcmparrPerp);

	  return true;
  }

}

// вычисление обратнолй матрицы произвольной размерности
//
bool   InverseMtrx( TComp *arrInp, const int nDimen,  TComp *arrOut)
{
  switch (nDimen)
  {
  case 1:
  if (fabs(arrInp[0].modul()) > 0.000000001)
  {
    arrOut[0] = TComp(1.,0.)/ arrInp[0];
    return true;
  }
  else
  {
      return false;
  }
  case 2:
  return InverseMtrx2(arrInp, arrOut);

  case 3:
  return InverseMtrx3(arrInp, arrOut);


  case 4:
  return InverseMtrx4(arrInp, arrOut);

  case 5:
  case 6:
  case 7:
  case 8:

  default: return false ;
      ;
  }
 // return false ;
}















//***********************************************************************************************
//********  Класс копмпл чисел для LONG DOUBLE ***************************************************************************************
//***********************************************************************************************
//---------------------------------------------------------------------------
TCompLong::TCompLong():
 m_Re(0)
,m_Im(0)
{

}
// конструктор копирования
 TCompLong ::TCompLong (const TCompLong &R)
 {
   m_Re = R.m_Re ;
   m_Im = R.m_Im;
 }
 // оператор присваивания
 TCompLong TCompLong::operator=(TCompLong  R)
 {
   m_Re = R.m_Re ;
   m_Im = R.m_Im;
   return *this ;
 }
// парам констр
TCompLong :: TCompLong(  const double x,const double y)

 {
   m_Re = x ;
   m_Im = y;
}
double TCompLong::modul()
{
	return sqrt(m_Re* m_Re + m_Im * m_Im);
}


