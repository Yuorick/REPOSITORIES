 #pragma hdrstop
#include <math.h>
#include "filter.h"
#include "StabSyst2.h"

//Параметры фильтра
const double q1x = 0.0;
const double q1y = 0.0;
const double q1h = 0.0;
const double q2x = 0.1;//0.01;//0.1;
const double q2y = 5.0;//1.0;//5.0;
const double q2h = 0.1;//0.01;//0.1;

// Минимальная высота
const double minH = -2.0;

// Ноль
const double ZERO = 1E-10;

//const double  R3 = 6371000.0;//Радиус Земли

// Признак начала работы фильтра
 int workFilter = 0;
 /* static double tops;
	static double k1x, k2x, k1y, k2y, k1h, k2h;

	// корреляционные матрицы
	static double p11x, p12x, p22x, p11y, p12y, p22y, p11h, p12h,
					p22h;

  static double x_g,   y_g,   h_g,   // coordinates
		x_vg,  y_vg,  h_vg,  // absolute velocity
		x_rvg, y_rvg, h_rvg; // relative velocity
  */
// INPUT:
// bSliga - если true, то работает старый вариант, если false то минимакс
// VAlSigW - скз шума цели (м*м/с/с/с/с)
// VAlSigMMO - скз ошибки оизмерения углов качки(рад*рад)
void filtr_ukf(Filter_Input *TFilter_Input, Filter_Output *TFilter_Output, bool bSliga, const double VAlDispW, const double VAlDispMMO)
{


   double x_s, y_s, h_s,//
	   x_vs, y_vs, h_vs,
	   x_ps,y_ps,h_ps;
	double se,ce,sb,cb;

  double x_kg,  y_kg,  h_kg, d_n;

  double Xkg,Ykg,Hkg,sbet,seps,sD,tk;

  static double tops;
	static double k1x, k2x, k1y, k2y, k1h, k2h;

	// корреляционные матрицы
	static double p11x, p12x, p22x, p11y, p12y, p22y, p11h, p12h,
					p22h;

  static double x_g,   y_g,   h_g,   // coordinates
		x_vg,  y_vg,  h_vg,  // absolute velocity
		x_rvg, y_rvg, h_rvg; // relative velocity

 double ys2 = (TFilter_Input->Xkg) * (TFilter_Input->Xkg) + (TFilter_Input->Ykg) * (TFilter_Input->Ykg)  + (TFilter_Input->Hkg) * (TFilter_Input->Hkg) ;
 double valPx =0.,valPy =0., valPh =0., valAlfx =0., valAlfy =0., valAlfh =0.;
	if(!workFilter)
	{
	  //Входные данные в алгоритм фильтрации
	  	Xkg=TFilter_Input->Xkg;
			Ykg=TFilter_Input->Ykg;
	    Hkg=TFilter_Input->Hkg;
	    sbet=TFilter_Input->sbet;
	    seps=TFilter_Input->seps;
	    sD=TFilter_Input->sD;
	    tk=TFilter_Input->tk;

		//Присвоение начальных условий на фильтр Калмана
		y_g  = Ykg;
		x_g  = Xkg;
		h_g  = Hkg;
		x_vg = h_vg = 0.0;
		y_vg = 0.0;


		 //	double ys2 = Xkg*Xkg+Ykg*Ykg+Hkg*Hkg;
			 if (bSliga)
			 {
				p11x = ys2*sbet;
				p11y = sD;
				p11h = ys2*seps;
				p12x = 0.0;
				p12y = 0.0;
				p12h = 0.0;
				p22h = 100.0;
				p22x = 1000.0;
				p22y = 10000.0;
			 }
			 else
			 {
				p11x = ys2* (sbet + VAlDispMMO);
				p11y = sD;
				p11h = ys2* (seps+ VAlDispMMO);
				p22h = 250000.;
				p22x = 250000.;
				p22y = 250000.;
				p12x = -0.9 * sqrt(p11x * p22x);
				p12y = -0.9 * sqrt(p11y * p22y);
				p12h = -0.9 * sqrt(p11h * p22h);
       }

		workFilter = 1;
	}
	else
	{
		double sp11x, sp11y, sp11h, sp12x, sp12y, sp12h, sp22x,
		sp22y, sp22h;
		{
			// Экстраполяция оценок координат и элементов матриц ковариаций
			//на время t
			double t = (TFilter_Input->tk)-tops, t2=t*t;

			x_g += t*x_rvg;
			y_g += t*y_rvg;
			h_g += t*h_rvg;
			if (bSliga)
			{
			if(h_g<minH)
			h_g=minH;

			sp11x = p11x+t2*p22x+t*(2.0*p12x+q1x);
			sp11y = p11y+t2*p22y+t*(2.0*p12y+q1y);
			sp11h = p11h+t2*p22h+t*(2.0*p12h+q1h);

			sp12x = p12x+t*p22x;
			sp12y = p12y+t*p22y;
			sp12h = p12h+t*p22h;

			sp22x = p22x+t*q2x;
			sp22y = p22y+t*q2y;
			sp22h = p22h+t*q2h;
			}
			else
			{

				 fncExtrapolateCorMtrx ( p11x, p12x, p22x, t, VAlDispW, &sp11x, &sp12x, &sp22x );
				 fncExtrapolateCorMtrx ( p11y, p12y, p22y, t, VAlDispW, &sp11y, &sp12y, &sp22y );
				 fncExtrapolateCorMtrx ( p11h, p12h, p22h, t, VAlDispW, &sp11h, &sp12h, &sp22h );
			}



		}

		//Вычисление коэффициентов усиления
			{
			double od;

			Xkg=TFilter_Input->Xkg;
			Ykg=TFilter_Input->Ykg;
			Hkg=TFilter_Input->Hkg;
			sbet=TFilter_Input->sbet;
			seps=TFilter_Input->seps;
			sD=TFilter_Input->sD;
			tk=TFilter_Input->tk;

			double ysk2 = Xkg*Xkg+Ykg*Ykg+Hkg*Hkg;
			if (bSliga)
			{
				od = sp11h+ysk2*seps;
				if( od>ZERO )
				{
				k1h   = sp11h/od;
				k2h   = sp12h/od;
				}

				od = sp11y+sD;
				if( od>ZERO )
				{
				k1y   = sp11y/od;
				k2y   = sp12y/od;
				}

				od = sp11x+ysk2*sbet;
				if( od>ZERO )
				{
				k1x   = sp11x/od;
				k2x   = sp12x/od;
				}
			}
			else
			{
			calcCoefUsilenia (sp11x, sp12x,sqrt(ysk2*sbet), sqrt(ysk2 * VAlDispMMO), &k1x, &k2x, &valPx, &valAlfx) ;
			calcCoefUsilenia (sp11y, sp12y, sD , 0., &k1y, &k2y, &valPy, &valAlfy) ;
			calcCoefUsilenia (sp11h, sp12h,sqrt(ysk2 * seps), sqrt(ysk2 * VAlDispMMO), &k1h, &k2h, &valPh, &valAlfh) ;
			}
		}

		// Вычисление корреляционных матриц на момент привязки последнего замера
		{

			if (bSliga)
			{
			double mk1;
			mk1   = 1.0-k1h;
			p11h = sp11h*mk1;
			p12h = sp12h*mk1;

			mk1   = 1.0-k1y;
			p11y = sp11y*mk1;
			p12y = sp12y*mk1;

			mk1   = 1.0-k1x;
			p11x = sp11x*mk1;
			p12x = sp12x*mk1;

			p22x = sp22x-k2x*sp12x;
			p22y = sp22y-k2y*sp12y;
			p22h = sp22h-k2h*sp12h;
			}
			else
			{
			 calcCorMtrxCompleted(sp11x, sp12x, sp22x,valPx, valAlfx, k2x, &p11x, &p12x, &p22x) ;
			 calcCorMtrxCompleted(sp11y, sp12y, sp22y,valPy, valAlfy, k2y, &p11y, &p12y, &p22y) ;
			 calcCorMtrxCompleted(sp11h, sp12h, sp22h,valPh, valAlfh, k2h, &p11h, &p12h, &p22h) ;
			}

		}


		x_kg=Xkg;
		y_kg=Ykg;
		h_kg=Hkg;
		h_kg = h_kg/*-(pow(x_kg,2)+pow(y_kg,2))/(2.0*R3)*/;

		{//Вычисление направляющих векторов целевой системы координат(ЦСК)
			double x_n, y_n, h_n, h_pr;
			h_pr = h_g/*-(pow(x_g,2)+pow(y_g,2))/(2.0*R3)*/;
			// d_n = sqrt(pow(x_g,2)+pow(y_g,2)+pow(h_pr,2));
			d_n = sqrt(x_g * x_g + y_g * y_g + h_pr * h_pr);
			x_n = x_g/d_n;
			y_n = y_g/d_n;
			h_n = h_pr/d_n;
			se = h_n;
			ce = ( fabs(se) > 0.999999999999)?0: sqrt(1.0-se*se);
			sb = x_n/ce;
			cb = y_n/ce;
		}

		//Преобразование первичных оценок в ЦСК
		x_ps = x_kg*cb - y_kg*sb;
		y_ps = x_kg*(sb*ce) + y_kg*(cb*ce) + h_kg*se;
		h_ps = -x_kg*(sb*se) - y_kg*(cb*se) + h_kg*ce;

		//Преобразование экстраполированных оценок координат и скоростей в ЦСК
		x_s = 0.0;
		y_s = d_n;
		h_s = 0.0;

		x_vs = x_vg*cb - y_vg*sb;
		y_vs = x_vg*(sb*ce) + y_vg*(cb*ce) + h_vg*se;
		h_vs = -x_vg*(sb*se) - y_vg*(cb*se) + h_vg*ce;

		//Вычисление оценок координат и скоростей цели
		double ra;
		ra = h_ps - h_s;
		h_s  += k1h*ra;
		h_vs += k2h*ra;

		ra = y_ps - y_s;
		y_s += k1y*ra;
		y_vs += k2y*ra;

		ra = x_ps - x_s;
		x_s += k1x*ra;
		x_vs += k2x*ra;

		//Преобразование оценок координат и скоростей цели в КГСК-ЦТ
		x_g = x_s*cb + y_s*(sb*ce) - h_s*(sb*se);
		y_g = -x_s*sb + y_s*(cb*ce) - h_s*(cb*se);
		h_g = y_s*se + h_s*ce;

		x_vg = x_vs*cb + y_vs*(sb*ce) - h_vs*(sb*se);
		y_vg = -x_vs*sb + y_vs*(cb*ce) - h_vs*(cb*se);
		h_vg = y_vs*se + h_vs*ce;

  }

  //Запоминаем предыдущее время съема данных
  tops=tk;

  //Вычисление относительных скоростей цели
  x_rvg = x_vg-(TFilter_Input->Vcx);
  y_rvg = y_vg-(TFilter_Input->Vcy);
  h_rvg = h_vg-(TFilter_Input->Vch);

  if(h_g<minH)
	h_g=minH;

  //Пересылка в структуру выходных данных
  TFilter_Output->X = x_g;
  TFilter_Output->Y = y_g;
  TFilter_Output->H = h_g;

  TFilter_Output->Vx = x_rvg;
  TFilter_Output->Vy = y_rvg;
  TFilter_Output->Vh = h_rvg;

}

void  fncExtrapolateCorMtrx (const double  p11h, const double p12h, const double p22h
		,const double  t,const double  VAlDispW, double *sp11h, double *sp12h, double *sp22h )
{
				TStabSyst2 StabSyst2(t, sqrt(VAlDispW), 0.,0. ) ;
				double arrK[4] ={0.}, arrKExtr[4] ={0.};
				arrK[0] = p11h;
				arrK[1] = p12h;
				arrK[2] = p12h;
				arrK[3] = p22h;
				StabSyst2.OneStepGolubevExtrapolation(arrK, t, arrKExtr) ;
				*sp11h = arrKExtr[0];
				*sp12h = arrKExtr[1];
				*sp22h = arrKExtr[3];
 }

 void  calcCoefUsilenia (const double  sp11h, const double sp12h
		,const double  VAlSigBMO, const double  VAlSigMMO, double *pUs0, double *pUs1, double *pvalP, double *pvalAlf)
{

 // Выбор ММО
	double alf =  VAlSigMMO/ sqrt( sp11h) ;
	if (alf > 1.) alf = 1. ;
	double valF = ( 1. - alf)/( ( 1.- alf)* ( 1.- alf) * sp11h + VAlSigBMO * VAlSigBMO) ;

	*pUs0 = sp11h* valF ;
	*pUs1 =  valF * sp12h ;
	if ((*pUs0) > 1. )
	{
		*pUs0 = 1. ;
		*pUs1 = sp12h / sp11h;
		alf = (VAlSigMMO * VAlSigMMO + VAlSigBMO * VAlSigBMO) / sp11h ;
	}
	*pvalP = 1. - (1.- alf) * (*pUs0);
	*pvalAlf = alf;
}

void calcCorMtrxCompleted(const double VAlKExtr00, const double VAlKExtr01, const double VAlKExtr11
	,const double VAlP,const double VAlAlf,const double VAlP1
	, double* pvalK00, double* pvalK01, double* pvalK11)
{
 *pvalK00 = VAlP * VAlKExtr00 ;
 *pvalK01 = VAlP * VAlKExtr01 ;
 *pvalK11 = VAlKExtr11 -  ( 1. - VAlAlf) * VAlP1 * VAlKExtr01 ;
}

 #pragma package(smart_init)



