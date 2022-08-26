﻿
#include "Equations.h"
#include "math.h"
#include <float.h>
//#include <vcl.h>
#include "Comp.h"
#include  <string.h>
const double EPS = 0.00000001;
//extern double PI;


//***********************************************************************************************
//********* ФУНКЦИИ ДЛЯ DOUBLE**************************************************************************************
//***********************************************************************************************

// Решение квадраьного улавнения a*x*x + b*x + c =0
// возвращает
// 0 - 2 действительных некраьных корня
// 1 - 2 действительных кратных корня
// 2- имеется по крайней мере один нулевой корень a!= 0, c=0
// 3 - 2 комплексно сопряженных корня
// 4 -  1 действительный корень (a =0)
// 5  - несовместность a=b=0, c!=0
// 6 - тождество )a=b=c=0
int SolvEq2(const double a,const double b,const double c,TComp &x1,TComp &x2)
{
	x1.m_Re = 0;
	x1.m_Im = 0;
	x2.m_Re = 0;
	x2.m_Im = 0;
	if (fabs(a) < 0.000000000001)
	{
	  if (fabs(b) < 0.000000000001)
	  {
		if (fabs(c) < 0.000000000001) return 6; //тождество a=b=c=0

		return 5 ;  // несовместность a=b=0, c!=0
	  }
	  x1.m_Re = -c/b ;  // 1 действительный корень (a =0)
	  return 4 ;
	}
    if (fabs(c) < DBL_MIN *10.)
	{
	  x1.m_Re = - b/a;
	 return 2;       // имеется по крайней мере один нулевой корень a!= 0, c=0
	}

	if(fabs(c) > DBL_MAX/10)
	{
		////String(L"fabs(c) > DBL_MAX/10");
		return 10;
    }
	double d2 = b * b - 4 * c * a ;
	if (d2 > DBL_MIN)
	{
	 double d = sqrt(d2);
	 x1.m_Re = (-b + d)/ 2/a ;
	 x2.m_Re = (-b - d)/ 2/a ;
	 return 0;

	}
	if (d2 < - DBL_MIN)
	{
	  x1.m_Re = x2.m_Re = -b/2/a ;
      x1.m_Im = sqrt(-d2)/2./a;
      x2.m_Im = -x1.m_Im;
	  return 3;
	}
	x1.m_Re = x2.m_Re =  -b/2/a ;//sqrt(d2);
	return 1;

}

// Решение квадраьного улавнения a*x*x + b*x + c =0
// возвращает
// 0 - 2 некраьных корня
// 1 - 2  кратных корня
// 2- имеется по крайней мере один нулевой корень a!= 0, c=0

// 4 -  1  корень (a =0)
// 5  - несовместность a=b=0, c!=0
// 6 - тождество )a=b=c=0
int SolvEq2( TComp a, TComp b,  TComp c,TComp &x1,TComp &x2)
{
	x1.m_Re = 0;
	x1.m_Im = 0;
	x2.m_Re = 0;
	x2.m_Im = 0;
	TComp cmpMinus1(-1., 0.);
	if (a.modul() < 0.000000000001)
	{
	  if (b.modul() < 0.000000000001)
	  {
		if (c.modul() < 0.000000000001) return 6; //тождество a=b=c=0

		return 5 ;  // несовместность a=b=0, c!=0
	  }
	  x1 = cmpMinus1* c/b ;  // 1  корень (a =0)
	  return 4 ;
	}
	if (c.modul() < 0.000000000001)
	{
	  x1 = cmpMinus1* b/a;
	  x2 = TComp(0.,0.);
	 return 2;       // имеется по крайней мере один нулевой корень a!= 0, c=0
	}

	if(c.modul() > DBL_MAX/10)
	{
		////String(L"c.modul() > DBL_MAX/10");
		return 10;
	}
	TComp cmp4(4., 0.);
	TComp cmp2(2., 0.);
	TComp d2 = b * b - cmp4 * c * a ;
	if (d2.modul() > DBL_MIN)
	{
	 TComp d = d2.Sqrt_();
	 x1 = ( d -b)/ cmp2/a ;
	 x2 = cmpMinus1 * (b +d)/ cmp2/a ;
	 return 0;

	}
	if (d2.modul()  < - DBL_MIN)
	{
	  x1 = cmpMinus1 * b/cmp2/a ;
	  x2 = x1 ;
	  return 1;
	}

  return 7;
}

/*  http://algolist.manual.ru/maths/findroot/cubic.php
Cubic equation solution. Real coefficients case.
  x3+a*x2+b*x+c=0.
   int Cubic(double *x,double a,double b,double c);
   Parameters:
   x - solution array (size 3). On output:
       3 real roots -> then x is filled with them;
       1 real + 2 complex -> x[0] is real, x[1] is real part of
                             complex roots, x[2] - non-negative
                             imaginary part.
   a, b, c - coefficients, as described
   Returns:
			3 - 3 real different roots;
			2 - 2 real root
			1 - 1 real root ;
			0 - 1 real root + 2 complex;


*/

// поиск собственных чисел  матрицы arrKInp 3х 3
// arrLamb[3]  массив корней уравнения new 30.01.2020
int   CalcProper_Numbers_R3(double * arrKInp, double *arrLamb)
{
	 // Решение характеристического ур-я для матрицы  arrKInp


	 double a = - ( arrKInp [0] + arrKInp [4]  + arrKInp [8] ) ;
	 double b =  arrKInp [0] * arrKInp [4] +  arrKInp [0] * arrKInp [8] + arrKInp [4] * arrKInp [8]
               - arrKInp [5] * arrKInp [7] - arrKInp [3] * arrKInp [1] - arrKInp [2] * arrKInp [6] ;
     double c = arrKInp [0] * arrKInp [5] * arrKInp [7]
             + arrKInp [1] * arrKInp [3] * arrKInp [8]
             +arrKInp [2] * arrKInp [4] * arrKInp [6]
              - arrKInp [0] * arrKInp [4]  * arrKInp [8]
              -arrKInp [1] * arrKInp [6]  * arrKInp [5]
             -arrKInp [2] * arrKInp [3]  * arrKInp [7];



	return Cubic(arrLamb, a, b, c) ;
}
/*
// поиск собственных чисел  матрицы arrKInp 3х 3
// arrLamb[3]  массив корней уравнения  old
int   CalcProper_Numbers_R3(double * arrKInp, double *arrLamb)
{
     // Решение характеристического ур-я для матрицы  arrKInp


     double a = - ( arrKInp [0] + arrKInp [4]  + arrKInp [8] ) ;
     double b =  arrKInp [0] * arrKInp [4] +  arrKInp [0] * arrKInp [8] + arrKInp [4] * arrKInp [8]
               - arrKInp [5] * arrKInp [7] - arrKInp [3] * arrKInp [1] - arrKInp [2] * arrKInp [6] ;
     double c = arrKInp [0] * arrKInp [5] * arrKInp [8]
             + arrKInp [4] * arrKInp [2] * arrKInp [6]
             +arrKInp [4] * arrKInp [3] * arrKInp [1]
              - arrKInp [0] * arrKInp [4]  * arrKInp [8]
              -arrKInp [2] * arrKInp [3]  * arrKInp [7]
             -arrKInp [3] * arrKInp [5]  * arrKInp [7];



    return Cubic(arrLamb, a, b, c) ;
}*/


int Cubic(double *x,double a,double b,double c)
{
  if (fabs(c) < EPS)
  {
	c = 0 ;
	if (fabs(b) < EPS)
	{
		if (fabs( a )< EPS)
		{
		  // все корни равны 0
		  x [0] = 0. ;
		  x [1] = 0  ;
		  x [2] = 0 ;
          return 2 ;
		}
	// решение линейного уравнения x * x * x + a * x * x = 0 , плюс 2 нулевых корня
	   x [0] = 0. ;
	   x [1] = 0  ;
	   x [2] = -a ;
	   return 2 ;

	}
     // решение квадр уравнения x * x * x + a * x * x  + b * x = 0, плюс нулевой корень
	TComp x1,x2 ;
	switch(SolvEq2(1,a,b,x1,x2) )
	{
	   case 0:   // 0 - 2 действительных некраьных корня
	   x [0] = 0. ;
	   x [1] = x1.m_Re  ;
	   x [2] = x2.m_Re ;
	   return 3 ;


	   case 1:  // 1 - 2 действительных кратных корня
	   x [0] = 0. ;
	   x [1] = x1.m_Re  ;
	   x [2] = x [1] ;
	   return 2 ;


	   case 3:  // 3 - 2 комплексно сопряженных корня
	   x [0] = 0. ;
	   x [1] = x1.m_Re  ;
	   x [2] = x1.m_Im ;
	   return 0 ;

	   default:
	   break;
    }


  }
  if (fabs(c) < EPS) c = 0 ;
  if (fabs(a) < EPS) a = 0 ;

 /* double q = ( a * a - 3 * b ) / 9. ;
  double r = ( 2 * a * a *a - 9 * a * b + 27 * c )/54.;
  double s = q * q * q - r * r ;  */
  double q = a/9. * a - b  / 3. ;
  double r =  a/ 27. * a *a -  a /6. * b +  c /2.;
  double s =  q - pow(fabs(r), 1/3.) * pow(fabs(r), 1/3.) ;

  if( fabs (s ) < EPS)
  {  // уравнение вырождено
	 if (fabs (sqrt(q)< EPS ))
	 {    // уравнение вырождено и имеет 1 корень кратности 3
	   x[0] = - 2. * Sign(r) * sqrt (q) - a/ 3. ;
	   x[1] = x[0];
	   x[2] = x[0] ;
       return 2 ;
	 }
	 else
	 {   // уравнение вырождено и имеет второй корень кратности 2
	  x[0] = - 2. * Sign(r) * sqrt (q) - a/ 3. ;
	  x[1] = Sign(r) * sqrt (q) - a/ 3. ;
	  x[2] = x[1] ;
	  return 2 ;
	 }

  }
  else
  {    if (s > 0)
	   {   // 3 действительных корня
		 double fi;// = acos(r / sqrt( q * q * q) )/ 3. ;
		 double temp = r / sqrt( q)/ sqrt( q) / sqrt(q);
		 if (temp > 0.999999999999999)
		 {
			fi = 0 ;
		 }
		 else
		 {
			 if (temp < -0.999999999999999)
			 {
			   fi = M_PI/3.;
			 }
			 else
			 {
			   fi = acos(temp )/ 3. ;
             }
         }
		 //double fi = acos(Sign( q ) * r / sqrt( fabs(q))/ sqrt( fabs(q)) / sqrt( fabs(q)) )/ 3. ;
		 x[0] = - 2. * sqrt( q) * cos (fi)  - a/ 3. ;
		 x[1] = - 2. * sqrt( q) * cos (fi + 2 * M_PI/3)  - a/ 3. ;
		 x[2] = - 2. * sqrt( q) * cos (fi - 2 * M_PI/3)  - a/ 3. ;
		return 3 ;

	   }
	   else
	   { //случай 2-х комплексных корней
		 if (fabs (q ) < EPS)
		 {
			x[0] = - pow(c - a * a * a / 27, 1./ 3.)  - a/ 3. ;
			x[1] =  - ( a + x[0])/2 ;
			x[2] = sqrt( fabs(( a - 3 * x [0]) * ( a + x [0]) - 4 * b)) / 2.;
		  return 1 ;
		 }
		 else
		 {
		   if (q > 0 )
		   {   // ?               ath.h?Еm СДЕЛАНО. В <
			 double t = fabs(r)/sqrt(q * q * q) ;
			 double fi = log( t + sqrt( t * t -1)) /3. ;
			 x[0] = - 2. * Sign(r) * sqrt(q) * cosh(fi) - a / 3. ;
			 x[1] =   Sign(r) * sqrt(q) * cosh(fi) - a / 3. ;
			 x[2] = sqrt(q * 3) * sinh(fi);

			   return 1 ;
		   }
		   else
		   {
			 double t = fabs(r)/sqrt(fabs(q * q * q))  ;
			 double fi = log( t + sqrt( t * t +1 )) /3. ;
			 x[0] = - 2. * Sign(r) * sqrt(fabs(q)) * sinh(fi) - a / 3. ;
			 x[1] =   Sign(r) * sqrt(fabs(q)) * sinh(fi) - a / 3. ;
			 x[2] = sqrt( fabs(q)  * 3) * cosh(fi);

			  return 1 ;
		   }
		  }

	   }

   }
}


// решение кубического уравнения с комплексными коэффициентами
//  a*z^3 + b * z^2 + c* z + d =0
// Input:
//cmpa, cmpb, cmpc, cmpd - коэффициенты
// OUTPUT:
// cmparrZ[3] - массив корней
// RETURN:
// 0 - 2 некраьных корня (cmpa =0)
// 1 - 2  кратных корня  (cmpa =0)
// 2- имеется по крайней мере один нулевой корень (cmpa =0 cmpb!= 0, cmpd =0
// 4 -  1  корень (cmpa =0, cmpb= 0)
// 5  - несовместность cmpa =0, cmpb= 0,cmpc =0, cmpd!=0
// 6 - тождество cmpa =0, cmpb= 0,cmpc =0, cmpd =0
// 7 -  cmpa !=0, cmpb= 0,cmpc =0, cmpd =0 3 нулевых корня
// 8 -  cmpa !=0, cmpb != 0,cmpc =0, cmpd =0 2 нулевых корня и один не нулевой
// 9 - имеется один нулевой корень
// 10 - имеется два нулевых кореня
// 11 - имеется 3 корня
int SolvEq3(TComp *cmparrZ,TComp cmpa,TComp cmpb,TComp cmpc,TComp cmpd)
{
   if (cmpa.modul()< EPS)
   {
		// Решение квадраьного улавнения a*x*x + b*x + c =0
		// возвращает
		// 0 - 2 некраьных корня
		// 1 - 2  кратных корня
		// 2- имеется по крайней мере один нулевой корень a!= 0, c=0

		// 4 -  1  корень (a =0)
		// 5  - несовместность a=b=0, c!=0
		// 6 - тождество )a=b=c=0
	  return   SolvEq2(  cmpb,  cmpc,   cmpd, cmparrZ[0],cmparrZ[1]);
   }
   TComp a = cmpb/ cmpa;
   TComp b = cmpc/ cmpa;
   TComp c = cmpd/ cmpa;
  if (c.modul() < EPS)
  {

	if (b.modul() < EPS)
	{
		if (a.modul()< EPS)
		{
		  // все корни равны 0
		  cmparrZ[0] = TComp(0.,0.);
		  cmparrZ[1] = TComp(0.,0.);
		  cmparrZ[2] = TComp(0.,0.);
		  return 7 ;
		}
	// решение линейного уравнения x * x * x + a * x * x = 0 , плюс 2 нулевых корня
	   cmparrZ[0] = TComp(0.,0.);
	   cmparrZ[1] = TComp(0.,0.);
	   cmparrZ[2] = a * TComp(-1.,0.);
	   return 8 ;

	}
	 // решение квадр уравнения x * x  + a * x  + b  = 0, плюс нулевой корень
	TComp x1,x2 ;
	cmparrZ[0] = TComp(0.,0.);
	TComp cmp1(1.,0.);

		// 0 - 2 некраьных корня
		// 1 - 2  кратных корня
		// 2- имеется по крайней мере один нулевой корень a!= 0, c=0

		// 4 -  1  корень (a =0) не может быть
		// 5  - несовместность a=b=0, c!=0  не может быть
		// 6 - тождество )a=b=c=0    не может быть
	switch(SolvEq2(cmp1,a,b,cmparrZ[1],cmparrZ[2]) )
	{
	   case 0:   // 0 - 2  некраьных корня
	   case 1:  // 1 - 2 х кратных корня
	   return 9 ;


	   case 2:  // 3 - 2 комплексно сопряженных корня

	   return 10 ;

	   default:
	   return -1;
	 //  break;
    }


  }
 if (c.modul() < EPS) c = TComp(0.,0.) ;
 if (a.modul() < EPS) a = TComp(0.,0.) ;



  TComp p =  b  - a * a / TComp(3.,0.) ;
  TComp  q =  c - a * b / TComp(3.,0.) + TComp(2.,0.) * a *a*a/TComp(27.,0.);

  TComp alfCub(0.,0.), betCub(0.,0.);
  TComp discr = q * q / TComp(4.,0.) + p * p * p / TComp(27.,0.);
  TComp root = discr.Sqrt_() ;
  alfCub = q / TComp(-2.,0.) + root ;
  betCub = q / TComp(-2.,0.) - root ;

  TComp cmparrAlf[3],cmparrBet[3];
  alfCub.root3(cmparrAlf);
  betCub.root3(cmparrBet);

  // выбор корня бетта
  TComp cmpBet1;
  for (int i = 0; i < 3; i++)
  {
	TComp cmpT = cmparrAlf[0] * cmparrBet[i] + p / TComp(3.,0.);
    if (cmpT.modul() < EPS * 1000.)
	{
	 cmpBet1 = cmparrBet[i];
	 break;
	}
  }
  TComp cmpEps (cos (2. * M_PI/ 3.), sin (2. * M_PI/ 3.));
   cmparrBet[0] = cmpBet1;
   cmparrAlf[1] = cmparrAlf[0] * cmpEps;
   cmparrAlf[2] = cmparrAlf[1] * cmpEps;
   cmparrBet[2] = cmparrBet[0] * cmpEps;
   cmparrBet[1] = cmparrBet[2] * cmpEps;

   // корни наполного кубического уравнения
   TComp cmparrX[3];
   for (int i = 0; i < 3; i++)
   {
	 cmparrX[i] = cmparrAlf[i] + cmparrBet[i];
   }
  ///

  // корни уравнения
  for (int i = 0; i < 3; i++)
   {
	 cmparrZ[i] = cmparrX[i] - a / TComp(3.,0.);
   }
  return 11;

}

double Sign(double a)
{
	if (fabs (a)< EPS) return 0;

	return (a > 0)?1:-1;
}

int fCompare (const void *e1, const void *e2)
{

	double *p1 =  (double*)e1 ;
	double *p2 =  (double*) e2 ;
	double e = (*p1) - (*p2) ;
	if (fabs (e) < EPS) return 0 ;
	 return ( e > 0 )?1 : -1;

}


//решение уравнения 5-ой степени с действит коэффициентами
// arra[5] *X^5 + arra[4] *X^4+...+ arra[0] = 0
// VAlMin , VAlMax - границы в которых лежит дейстьвительный корень
int SolvEq5(TComp *cmparrRoots, double *arra, const double VAlMin, const double VAlMax)
{
    int ireturn = -1;
     // нахождение дейсьтвительного корня
    double step =  0.001;
    int iNc = (VAlMax - VAlMin) / step;
    double x0 = VAlMin;
    double f0 = polinom(x0, arra, 5);
    double f1 = -1.;
    double x =VAlMax;
    int i =0;
    for ( i =0; i < iNc; ++i)
    {
         x = VAlMin + ((double)i +1.) * step;
         f1 =  polinom(x, arra, 5);
         if (f0 * f1 >0.)
         {
          x0 = x;
          f0 = f1;
          continue;
         }
         else
         {
           ireturn = 0;
           break;
         }

    }
    if ( -1 == ireturn)
    {
        return ireturn;
    }

    // метод ньютона
   /* double valF = polinom(x0, arra, 5);
    int j = 0 ;
    ireturn = -2;
    for ( j = 0; j < 1000; ++j)
    {
        double valDeriv = polinomDeriv(x0, arra, 5);
        double del = valF/ valDeriv;
        x0 = x0 - del;
        valF = polinom(x0, arra, 5);
        if (fabs(valF )< 0.000001)
        {
            ireturn = 0;
            break;
        }
    }*/
    int j = 0 ;
    ireturn = -2;
    f1 = f1 / 100000.;
    f0 = f0 / 100000.;
    for ( j = 0; j < 1000; ++j)
    {
        double valk = (f1 - f0)/ (x - x0);
        double x2 = -f0 / valk + x0;
        double f2 = polinom(x2, arra, 5)/ 100000.;
        if (fabs(f2 )< 0.000001)
        {
            ireturn = 0;
            break;
        }
        if (f2 * f0 > 0.)
        {
            x0 = x2;
            f0 = f2;
        }
        else
        {
            x = x2;
            f1 = f2;
        }


        }
    ///

    if ( -2 == ireturn)
    {
        return ireturn;
    }
   double t1 =  polinom(x0, arra, 5);

    //
    double arra1[6] = {0.};
    memcpy(arra1, arra, 6 * sizeof(double));
    if (fabs(arra1[5]) > 0.00000001)
    {
        for (int j = 0; j < 5; ++j)
        {
          arra1[j] /= arra1[5];
        }
    }
    ///

    // нахождение коэффициентов уравнения 4-ой степени
    double arrb0[5] = {0.};
     arrb0[3] = arra1[4] + x0;
     arrb0[2] = arra1[3] + x0 * arrb0[3];
     arrb0[1] = arra1[2] + x0 * arrb0[2];
     arrb0[0] = arra1[1] + x0 * arrb0[1];
     arrb0[4] = 1.;
    SolvEq4(cmparrRoots,  arrb0);
    cmparrRoots[4].m_Re = x0;
    cmparrRoots[4].m_Im = 0.;
    // проверка
    TComp cmparra[5];
    for(int i =0; i < 5; ++i)
    {
      cmparra[i] = TComp(arrb0[i], 0.);
    }
    for(int i =0; i < 4; ++i)
    {
      TComp cmprez = polinom(cmparrRoots[i], cmparra, 4);
      int iuu=0;
   }
    return ireturn;
}


//
// arra[len +1]
double polinom(const double x, double *arra, const int len)
{
    double sum = arra[0];
    double temp = 1.;
    for (int i =0; i < len; ++i)
    {
      temp = temp *x;
      sum += arra[i + 1] * temp;
    }
    return sum;
}


// arra[len +1]
TComp polinom(const TComp x, TComp *arra, const int len)
{
    TComp sum = arra[0];
    TComp temp (1., 0.);
    for (int i =0; i < len; ++i)
    {
      temp = temp *x;
      sum =sum + arra[i + 1] * temp;
    }
    return sum;
}

// arra[len +1]
TComp polinom(const TComp x, double *arra0, const int len)
{
    TComp *arra = new TComp[len +1];
    for (int i = 0; i < (len +1); ++i)
    {
      arra [i] =  TComp(arra0[i], 0.);
    }
    TComp sum = arra[0];
    TComp temp (1., 0.);
    for (int i =0; i < len; ++i)
    {
      temp = temp *x;
      sum =sum + arra[i + 1] * temp;
    }
    delete []arra;
    return sum;
}

double polinomDeriv(const double x, double *arra, const int len)
{
    double *arrb = new double[len];
    for (int i = 0; i < len; ++i)
    {
      arrb[ len -1 -i] = ((double)(len  -i)) * arra [ len  -i] ;
    }
    double xreturn = polinom( x, arrb,  len-1);
    delete []arrb;
    return xreturn;
}


//решение уравнения 4-ой степени с действит коэффициентами
//  arra[4] *X^4+...+ arra[0] = 0
int SolvEq4(TComp *cmparrRoots, double *arra)
{
   if(fabs(arra[4]) < 0.000000001)
   {
       return -1;
   }
    double a =arra[3] /arra[4] ;
    double b = arra[2] /arra[4] ;
    double c = arra[1]/arra[4] ;
    double d = arra[0] /arra[4] ;


    TComp cmpa(1.,0.);
            TComp cmpb(-b,0.);
            TComp cmpc((a * c - 4.*d),0.);
            TComp cmpd(-a * a * d + 4. * b * d - c * c,0.);
    TComp cmparrRoots3[3];

    SolvEq3(cmparrRoots3,cmpa, cmpb, cmpc, cmpd);

    // проверка
    TComp cmparrCoef[4] ;
    cmparrCoef[0] =cmpd;
    cmparrCoef[1] =cmpc;
   cmparrCoef[2] =cmpb;
    cmparrCoef[3] =cmpa;

    for(int i =0; i < 3; ++i)
    {
      TComp cmprez = polinom(cmparrRoots3[i], cmparrCoef, 3);
      int iii =0;
    }
    //////////////////////////////////////////////////



   // нахождение действительного корня
    double y = cmparrRoots3[0].m_Re;

    for (int n = 0; n < 3; ++n)
    {
        if(fabs(cmparrRoots3[n].m_Im) < 0.000000000001)
        {
            y = cmparrRoots3[n].m_Re;
            break;
        }
        int n0 = n;
        int n1 = (n+1) % 3;
        TComp cmpTemp = cmparrRoots3[n1].Sopr();
        TComp cmpTemp1 = cmparrRoots3[n0] - cmpTemp;
        if(cmpTemp1.modul() < 0.00001)
        {
            y = cmparrRoots3[(n+2) % 3].m_Re;
            break;
        }
    }


    double t = a * a / 4. - b +y;
    double u = sqrt(fabs(t));

   // double t1 = ((a / 2. * y - c) > 0.)?1.:-1.;
    t = y * y/4. -d;
    double v = sqrt(fabs(t));
  //  v = v * t1;

    // первое квадратное уравнение
    double a0 =1.;
    double b0 = a/ 2. - u;
    double c0 = y/2. - v;
    SolvEq2(a0, b0, c0,cmparrRoots[0],cmparrRoots[1]);
    // проверка
    TComp cmparra[5];
   for (int i = 0; i < 5; ++i)
    {
      cmparra[i] = TComp(arra[i], 0.);
    }
   TComp z0 = polinom(cmparrRoots[0], cmparra, 4);
   if (z0.modul() > 0.0001)
   {
       v = -v;
       a0 =1.;
       b0 = a/ 2. - u;
       c0 = y/2. - v;
       SolvEq2(a0, b0, c0,cmparrRoots[0],cmparrRoots[1]);
   }

     // второе  квадратное уравнение
     b0 = a/ 2. + u;
     c0 = y/2. + v;
    SolvEq2(a0, b0, c0,cmparrRoots[2],cmparrRoots[3]);


    return 0.;
}


//решение уравнения 4-ой степени с действит коэффициентами
//  arra[4] *X^4+...+ arra[0] = 0
int SolvEq4_(TComp *cmparrRoots, double *arra)
{
   if(fabs(arra[4]) < 0.000000001)
   {
       return -1;
   }
    double b =arra[3] /arra[4] ;
    double c = arra[2] /arra[4] ;
    double d = arra[1]/arra[4] ;
    double e = arra[0] /arra[4] ;


    TComp cmpa(8.,0.);
            TComp cmpb(-4. * c,0.);
            TComp cmpc((2. * b * d - 8.*e),0.);
            TComp cmpd(e * (4. * c - b*b) - d * d,0.);
    TComp cmparrRoots3[3];
    SolvEq3(cmparrRoots3,cmpa, cmpb, cmpc, cmpd);
   // нахождение действительного корня
    double y = cmparrRoots3[0].m_Re;

    for (int n = 0; n < 3; ++n)
    {
        if(fabs(cmparrRoots3[n].m_Im) < 0.000000000001)
        {
            y = cmparrRoots3[n].m_Re;
            break;
        }
        int n0 = n;
        int n1 = (n+1) % 3;
        TComp cmpTemp = cmparrRoots3[n1].Sopr();
        TComp cmpTemp1 = cmparrRoots3[n0] - cmpTemp;
        if(cmpTemp1.modul() < 0.00001)
        {
            y = cmparrRoots3[(n+2) % 3].m_Re;
            break;
        }
    }

    double valA = sqrt(8. * y + b * b - 4. * c);
    // первое квадратное уравнение
    double a0 =1.;
    double b0 = (b  + valA)/2.;
    double c0 =y + ( b * y -d)/ valA;
    SolvEq2(a0, b0, c0,cmparrRoots[0],cmparrRoots[1]);
    // второе  квадратное уравнение
    valA = -valA;
    b0 = (b  + valA)/2.;
    c0 =y + ( b * y -d)/ valA;
    SolvEq2(a0, b0, c0,cmparrRoots[2],cmparrRoots[3]);

   /* double t = a * a / 4. - b +y;
    double u = sqrt(fabs(t));

   // double t1 = ((a / 2. * y - c) > 0.)?1.:-1.;
    t = y * y/4. -d;
    double v = sqrt(fabs(t));
  //  v = v * t1;

    // первое квадратное уравнение
    double a0 =1.;
    double b0 = a/ 2. - u;
    double c0 = y/2. - v;
    SolvEq2(a0, b0, c0,cmparrRoots[0],cmparrRoots[1]);
    // проверка
    TComp cmparra[5];
   for (int i = 0; i < 5; ++i)
    {
      cmparra[i] = TComp(arra[i], 0.);
    }
   TComp z0 = polinom(cmparrRoots[0], cmparra, 4);
   if (z0.modul() > 0.0001)
   {
       v = -v;
       a0 =1.;
       b0 = a/ 2. - u;
       c0 = y/2. - v;
       SolvEq2(a0, b0, c0,cmparrRoots[0],cmparrRoots[1]);
   }

     // второе  квадратное уравнение
     b0 = a/ 2. + u;
     c0 = y/2. + v;
    SolvEq2(a0, b0, c0,cmparrRoots[2],cmparrRoots[3]);
*/

    return 0.;
}


//-------------------------------------------------------
//вычисление коэффициентов полинома 4-ой степени x^4 + b * x^3 + c * x^2 + d * x + e =0
// по его корням, по теореме Виета
void calcPolinom_Coef_4Degree_Vieta( TComp*CmpArrLamb, double *pb, double *pc, double *pd, double *pe)
{
 // вычисленпие коэффициентов полинома по теореме Виета
 TComp cmpb = ((CmpArrLamb[0] + CmpArrLamb[1]) + CmpArrLamb[2]) + CmpArrLamb[3];
 *pb = -cmpb.m_Re;

 TComp cmpc = (((((CmpArrLamb[0] * CmpArrLamb[1]) + (CmpArrLamb[0] * CmpArrLamb[2])) + (CmpArrLamb[0] * CmpArrLamb[3]))
           +(CmpArrLamb[1] * CmpArrLamb[2])) + (CmpArrLamb[1] * CmpArrLamb[3])) + (CmpArrLamb[2] * CmpArrLamb[3]);
 *pc = cmpc.m_Re;

 TComp cmpd = ((((CmpArrLamb[0] * CmpArrLamb[1])* CmpArrLamb[2]) + ((CmpArrLamb[0] * CmpArrLamb[1])* CmpArrLamb[3]))
         +((CmpArrLamb[0] * CmpArrLamb[2])* CmpArrLamb[3])) +((CmpArrLamb[1] * CmpArrLamb[2])* CmpArrLamb[3]);
 *pd = -cmpd.m_Re;

 TComp cmpe = ((CmpArrLamb[0] * CmpArrLamb[1])* CmpArrLamb[2])* CmpArrLamb[3];
 *pe = cmpe.m_Re;
}



//***********************************************************************************************
//********* ФУНКЦИИ ДЛЯ LONG DOUBLE**************************************************************************************
//***********************************************************************************************
 // Решение квадраьного улавнения a*x*x + b*x + c =0
// возвращает
// 0 - 2 действительных некраьных корня
// 1 - 2 действительных кратных корня
// 2- имеется по крайней мере один нулевой корень a!= 0, c=0
// 3 - 2 комплексно сопряженных корня
// 4 -  1 действительный корень (a =0)
// 5  - несовместность a=b=0, c!=0
// 6 - тождество )a=b=c=0
int SolvEq2(const long double a,const long double b,const long double c,TCompLong &x1,TCompLong &x2)
{
	x1.m_Re = 0;
	x1.m_Im = 0;
	x2.m_Re = 0;
	x2.m_Im = 0;
	if (fabsl(a) < 0.000000000001)
	{
	  if (fabsl(b) < 0.000000000001)
	  {
		if (fabsl(c) < 0.000000000001) return 6; //тождество a=b=c=0

		return 5 ;  // несовместность a=b=0, c!=0
	  }
	  x1.m_Re = -c/b ;  // 1 действительный корень (a =0)
	  return 4 ;
	}
	if (fabsl(c) < 0.000000000001)
	{
	  x1.m_Re = - b/a;
	 return 2;       // имеется по крайней мере один нулевой корень a!= 0, c=0
	}

	if(fabsl(c) > DBL_MAX/10)
	{
		////String(L"fabsl(c) > DBL_MAX/10");
		return 10;
    }
	long double d2 = b * b - 4 * c * a ;
	if (d2 > DBL_MIN)
	{
	 long double d = sqrtl(d2);
	 x1.m_Re = (-b + d)/ 2/a ;
	 x2.m_Re = (-b - d)/ 2/a ;
	 return 0;

	}
	if (d2 < - DBL_MIN)
	{
	  x1.m_Re = x2.m_Re = -b/2/a ;
	  x1.m_Im = sqrtl(-d2);
	  x2.m_Im = -x1.m_Im;
	  return 3;
	}
	x1.m_Re = x2.m_Re =  -b/2/a ;//sqrt(d2);
	return 1;

}



/*  http://algolist.manual.ru/maths/findroot/cubic.php
Cubic equation solution. Real coefficients case.
  x3+a*x2+b*x+c=0.
   int Cubic(long double *x,long double a,long double b,long double c);
   Parameters:
   x - solution array (size 3). On output:
       3 real roots -> then x is filled with them;
       1 real + 2 complex -> x[0] is real, x[1] is real part of
                             complex roots, x[2] - non-negative
                             imaginary part.
   a, b, c - coefficients, as described
   Returns:
			3 - 3 real different roots;
			2 - 2 real root
			1 - 1 real root ;
			0 - 1 real root + 2 complex;


*/

// поиск собственных чисел  матрицы arrKInp 3х 3


// arrLamb[3]  массив корней уравнения
int   CalcProper_Numbers_R3(long double * arrKInp, long double *arrLamb)
{
	 // Решение характеристического ур-я для матрицы  arrKInp


	 long double a = - ( arrKInp [0] + arrKInp [4]  + arrKInp [8] ) ;
	 long double b =  arrKInp [0] * arrKInp [4] +  arrKInp [0] * arrKInp [8] + arrKInp [4] * arrKInp [8]
			   - arrKInp [5] * arrKInp [5] - arrKInp [1] * arrKInp [1] - arrKInp [2] * arrKInp [2] ;
	 long double c = -( arrKInp [0] * arrKInp [4]  * arrKInp [8] + 2 *  arrKInp [5] * arrKInp [1] * arrKInp [2]
			 - arrKInp [0] * arrKInp [5] * arrKInp [5] - arrKInp [8] * arrKInp [1] * arrKInp [1] - arrKInp [4] * arrKInp [2] * arrKInp [2]) ;


	return Cubic(arrLamb, a, b, c) ;
}


int Cubic(long double *x,long double a,long double b,long double c)
{
  if (fabsl(c) < EPS)
  {
	c = 0 ;
	if (fabsl(b) < EPS)
	{
		if (fabsl( a )< EPS)
		{
		  // все корни равны 0
		  x [0] = 0. ;
		  x [1] = 0  ;
		  x [2] = 0 ;
		  return 1 ;
		}
	// решение линейного уравнения x * x * x + a * x * x = 0 , плюс 2 нулевых корня
	   x [0] = 0. ;
	   x [1] = 0  ;
	   x [2] = -a ;
	   return 2 ;

	}
     // решение квадр уравнения x * x * x + a * x * x  + b * x = 0, плюс нулевой корень
	TCompLong x1,x2 ;
	switch(SolvEq2(1,a,b,x1,x2) )
	{
	   case 0:   // 0 - 2 действительных некраьных корня
	   x [0] = 0. ;
	   x [1] = x1.m_Re  ;
	   x [2] = x2.m_Re ;
	   return 3 ;

	   case 1:  // 1 - 2 действительных кратных корня
	   x [0] = 0. ;
	   x [1] = x1.m_Re  ;
	   x [2] = x [1] ;
	   return 2 ;


	   case 3:  // 3 - 2 комплексно сопряженных корня
	   x [0] = 0. ;
	   x [1] = x1.m_Re  ;
	   x [2] = x1.m_Im ;
	   return 0 ;

	   default:
	   break;
    }


  }
  if (fabsl(c) < EPS) c = 0 ;
  if (fabsl(a) < EPS) a = 0 ;

 /* long double q = ( a * a - 3 * b ) / 9. ;
  long double r = ( 2 * a * a *a - 9 * a * b + 27 * c )/54.;
  long double s = q * q * q - r * r ;  */
  long double q = a/9. * a - b  / 3. ;
  long double r =  a/ 27. * a *a -  a /6. * b +  c /2.;
  long double s =  q - powl(fabsl(r), 1/3.) * powl(fabsl(r), 1/3.) ;

  if( fabsl (s ) < EPS)
  {  // уравнение вырождено
	 if (fabs (sqrtl(q)< EPS ))
	 {    // уравнение вырождено и имеет 1 корень кратности 3
	   x[0] = - 2. * Sign(r) * sqrtl (q) - a/ 3. ;
	   x[1] = x[0];
	   x[2] = x[0] ;
	   return 1 ;
	 }
	 else
	 {   // уравнение вырождено и имеет второй корень кратности 2
	  x[0] = - 2. * Sign(r) * sqrtl (q) - a/ 3. ;
	  x[1] = Sign(r) * sqrtl (q) - a/ 3. ;
	  x[2] = x[1] ;
	  return 2 ;
	 }

  }
  else
  {    if (s > 0)
	   {   // 3 действительных корня
		 long double fi;// = acos(r / sqrt( q * q * q) )/ 3. ;
		 long double temp = r / sqrtl( q)/ sqrtl( q) / sqrtl(q);
		 if (temp > 0.999999999999999)
		 {
			fi = 0 ;
		 }
		 else
		 {
			 if (temp < -0.999999999999999)
			 {
			   fi = M_PI/3.;
			 }
			 else
			 {
			   fi = acosl(temp )/ 3. ;
             }
         }
		 //long double fi = acos(Sign( q ) * r / sqrt( fabs(q))/ sqrt( fabs(q)) / sqrt( fabs(q)) )/ 3. ;
		 x[0] = - 2. * sqrtl( q) * cosl (fi)  - a/ 3. ;
		 x[1] = - 2. * sqrtl( q) * cosl (fi + 2 * M_PI/3)  - a/ 3. ;
		 x[2] = - 2. * sqrtl( q) * cosl (fi - 2 * M_PI/3)  - a/ 3. ;
		return 3 ;

	   }
	   else
	   { //случай 2-х комплексных корней
		 if (fabsl (q ) < EPS)
		 {
			x[0] = - powl(c - a * a * a / 27, 1./ 3.)  - a/ 3. ;
			x[1] =  - ( a + x[0])/2 ;
			x[2] = sqrtl( fabsl(( a - 3 * x [0]) * ( a + x [0]) - 4 * b)) / 2.;
		  return 1 ;
		 }
		 else
		 {
		   if (q > 0 )
		   {   // ?               ath.h?Еm СДЕЛАНО. В <
			 long double t = fabsl(r)/sqrtl(q * q * q) ;
			 long double fi = logl( t + sqrtl( t * t -1)) /3. ;
			 x[0] = - 2. * Sign(r) * sqrtl(q) * coshl(fi) - a / 3. ;
			 x[1] =   Sign(r) * sqrtl(q) * coshl(fi) - a / 3. ;
			 x[2] = sqrtl(q * 3) * sinhl(fi);

			   return 1 ;
		   }
		   else
		   {
			 long double t = fabsl(r)/sqrtl(fabsl(q * q * q))  ;
			 long double fi = logl( t + sqrtl( t * t +1 )) /3. ;
			 x[0] = - 2. * Sign(r) * sqrtl(fabsl(q)) * sinhl(fi) - a / 3. ;
			 x[1] =   Sign(r) * sqrtl(fabs(q)) * sinhl(fi) - a / 3. ;
			 x[2] = sqrtl( fabsl(q)  * 3) * coshl(fi);

			  return 1 ;
		   }
		  }

	   }

   }
}

long double Sign(long double a)
{
	if (fabsl (a)< EPS) return 0;

	return (a > 0)?1:-1;
}






//***********************************************************************************************
//********* ФУНКЦИИ ДЛЯ float**************************************************************************************
//***********************************************************************************************

// Решение квадраьного улавнения a*x*x + b*x + c =0
// возвращает
// 0 - 2 действительных некраьных корня
// 1 - 2 действительных кратных корня
// 2- имеется по крайней мере один нулевой корень a!= 0, c=0
// 3 - 2 комплексно сопряженных корня
// 4 -  1 действительный корень (a =0)
// 5  - несовместность a=b=0, c!=0
// 6 - тождество )a=b=c=0
int SolvEq2(const float a,const float b,const float c,TComp &x1,TComp &x2)
{
    x1.m_Re = 0;
    x1.m_Im = 0;
    x2.m_Re = 0;
    x2.m_Im = 0;
    if (fabs(a) < 0.000000000001)
    {
      if (fabs(b) < 0.000000000001)
      {
        if (fabs(c) < 0.000000000001) return 6; //тождество a=b=c=0

        return 5 ;  // несовместность a=b=0, c!=0
      }
      x1.m_Re = -c/b ;  // 1 действительный корень (a =0)
      return 4 ;
    }
    if (fabs(c) < 0.000000000001)
    {
      x1.m_Re = - b/a;
     return 2;       // имеется по крайней мере один нулевой корень a!= 0, c=0
    }

    if(fabs(c) > DBL_MAX/10)
    {
        ////String(L"fabs(c) > DBL_MAX/10");
        return 10;
    }
    float d2 = b * b - 4 * c * a ;
    if (d2 > DBL_MIN)
    {
     float d = sqrt(d2);
     x1.m_Re = (-b + d)/ 2/a ;
     x2.m_Re = (-b - d)/ 2/a ;
     return 0;

    }
    if (d2 < - DBL_MIN)
    {
      x1.m_Re = x2.m_Re = -b/2/a ;
      x1.m_Im = sqrt(-d2);
      x2.m_Im = -x1.m_Im;
      return 3;
    }
    x1.m_Re = x2.m_Re =  -b/2/a ;//sqrt(d2);
    return 1;

}

/*  http://algolist.manual.ru/maths/findroot/cubic.php
Cubic equation solution. Real coefficients case.
  x3+a*x2+b*x+c=0.
   int Cubic(float *x,float a,float b,float c);
   Parameters:
   x - solution array (size 3). On output:
       3 real roots -> then x is filled with them;
       1 real + 2 complex -> x[0] is real, x[1] is real part of
                             complex roots, x[2] - non-negative
                             imaginary part.
   a, b, c - coefficients, as described
   Returns:
            3 - 3 real different roots;
            2 - 2 real root
            1 - 1 real root ;
            0 - 1 real root + 2 complex;


*/

// поиск собственных чисел  матрицы arrKInp 3х 3


// arrLamb[3]  массив корней уравнения
int   CalcProper_Numbers_R3(float * arrKInp, float *arrLamb)
{
     // Решение характеристического ур-я для матрицы  arrKInp


     float a = - ( arrKInp [0] + arrKInp [4]  + arrKInp [8] ) ;
     float b =  arrKInp [0] * arrKInp [4] +  arrKInp [0] * arrKInp [8] + arrKInp [4] * arrKInp [8]
               - arrKInp [5] * arrKInp [5] - arrKInp [1] * arrKInp [1] - arrKInp [2] * arrKInp [2] ;
     float c = -( arrKInp [0] * arrKInp [4]  * arrKInp [8] + 2 *  arrKInp [5] * arrKInp [1] * arrKInp [2]
             - arrKInp [0] * arrKInp [5] * arrKInp [5] - arrKInp [8] * arrKInp [1] * arrKInp [1] - arrKInp [4] * arrKInp [2] * arrKInp [2]) ;


    return Cubic(arrLamb, a, b, c) ;
}


int Cubic(float *x,float a,float b,float c)
{
  if (fabs(c) < EPS)
  {
    c = 0 ;
    if (fabs(b) < EPS)
    {
        if (fabs( a )< EPS)
        {
          // все корни равны 0
          x [0] = 0. ;
          x [1] = 0  ;
          x [2] = 0 ;
          return 1 ;
        }
    // решение линейного уравнения x * x * x + a * x * x = 0 , плюс 2 нулевых корня
       x [0] = 0. ;
       x [1] = 0  ;
       x [2] = -a ;
       return 2 ;

    }
     // решение квадр уравнения x * x * x + a * x * x  + b * x = 0, плюс нулевой корень
    TComp x1,x2 ;
    switch(SolvEq2(1,a,b,x1,x2) )
    {
       case 0:   // 0 - 2 действительных некраьных корня
       x [0] = 0. ;
       x [1] = x1.m_Re  ;
       x [2] = x2.m_Re ;
       return 3 ;


       case 1:  // 1 - 2 действительных кратных корня
       x [0] = 0. ;
       x [1] = x1.m_Re  ;
       x [2] = x [1] ;
       return 2 ;


       case 3:  // 3 - 2 комплексно сопряженных корня
       x [0] = 0. ;
       x [1] = x1.m_Re  ;
       x [2] = x1.m_Im ;
       return 0 ;

       default:
       break;
    }


  }
  if (fabs(c) < EPS) c = 0 ;
  if (fabs(a) < EPS) a = 0 ;

 /* float q = ( a * a - 3 * b ) / 9. ;
  float r = ( 2 * a * a *a - 9 * a * b + 27 * c )/54.;
  float s = q * q * q - r * r ;  */
  float q = a/9. * a - b  / 3. ;
  float r =  a/ 27. * a *a -  a /6. * b +  c /2.;
  float s =  q - pow(fabs(r), 1/3.) * pow(fabs(r), 1/3.) ;

  if( fabs (s ) < EPS)
  {  // уравнение вырождено
     if (fabs (sqrt(q)< EPS ))
     {    // уравнение вырождено и имеет 1 корень кратности 3
       x[0] = - 2. * Sign(r) * sqrt (q) - a/ 3. ;
       x[1] = x[0];
       x[2] = x[0] ;
       return 1 ;
     }
     else
     {   // уравнение вырождено и имеет второй корень кратности 2
      x[0] = - 2. * Sign(r) * sqrt (q) - a/ 3. ;
      x[1] = Sign(r) * sqrt (q) - a/ 3. ;
      x[2] = x[1] ;
      return 2 ;
     }

  }
  else
  {    if (s > 0)
       {   // 3 действительных корня
         float fi;// = acos(r / sqrt( q * q * q) )/ 3. ;
         float temp = r / sqrt( q)/ sqrt( q) / sqrt(q);
         if (temp > 0.999999999999999)
         {
            fi = 0 ;
         }
         else
         {
             if (temp < -0.999999999999999)
             {
               fi = M_PI/3.;
             }
             else
             {
               fi = acos(temp )/ 3. ;
             }
         }
         //float fi = acos(Sign( q ) * r / sqrt( fabs(q))/ sqrt( fabs(q)) / sqrt( fabs(q)) )/ 3. ;
         x[0] = - 2. * sqrt( q) * cos (fi)  - a/ 3. ;
         x[1] = - 2. * sqrt( q) * cos (fi + 2 * M_PI/3)  - a/ 3. ;
         x[2] = - 2. * sqrt( q) * cos (fi - 2 * M_PI/3)  - a/ 3. ;
        return 3 ;

       }
       else
       { //случай 2-х комплексных корней
         if (fabs (q ) < EPS)
         {
            x[0] = - pow(c - a * a * a / 27, 1./ 3.)  - a/ 3. ;
            x[1] =  - ( a + x[0])/2 ;
            x[2] = sqrt( fabs(( a - 3 * x [0]) * ( a + x [0]) - 4 * b)) / 2.;
          return 1 ;
         }
         else
         {
           if (q > 0 )
           {   // ?               ath.h?Еm СДЕЛАНО. В <
             float t = fabs(r)/sqrt(q * q * q) ;
             float fi = log( t + sqrt( t * t -1)) /3. ;
             x[0] = - 2. * Sign(r) * sqrt(q) * cosh(fi) - a / 3. ;
             x[1] =   Sign(r) * sqrt(q) * cosh(fi) - a / 3. ;
             x[2] = sqrt(q * 3) * sinh(fi);

               return 1 ;
           }
           else
           {
             float t = fabs(r)/sqrt(fabs(q * q * q))  ;
             float fi = log( t + sqrt( t * t +1 )) /3. ;
             x[0] = - 2. * Sign(r) * sqrt(fabs(q)) * sinh(fi) - a / 3. ;
             x[1] =   Sign(r) * sqrt(fabs(q)) * sinh(fi) - a / 3. ;
             x[2] = sqrt( fabs(q)  * 3) * cosh(fi);

              return 1 ;
           }
          }

       }

   }
}


float Sign(float a)
{
    if (fabs (a)< EPS) return 0;

    return (a > 0)?1:-1;
}



//#pragma package(smart_init)
