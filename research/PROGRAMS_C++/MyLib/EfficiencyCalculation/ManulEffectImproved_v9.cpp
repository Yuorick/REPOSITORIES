//---------------------------------------------------------------------------

#pragma hdrstop
#include <vcl.h>
#include "ManulEffectImproved_v9.h"
//#include "Far_2D.h"
 #include "Far.h"
 #include "Comp.h"
 #include "URPointZ.h"
 #include "MatrixProccess.h"


//---------------------------------------------------------------------------

//#pragma hdrstop

//#include <tchar.h>

//---------------------------------------------------------------------------

//#pragma argsused
 #include <math.h>

 #include <exception>

 #include "filter.h"
 #include "Gauss.h"
 #include "YrWriteShapeFile.h"
 #include "Target.h"
 #include "InitTargData.h"
 #include "Diagrams.h"
 #include "Equations.h"

extern int workFilter;


 extern const double DTi = 0.01; // Временной интервал   ????
 //const double DT_filtr = 1.0/8.0; //Шаг фильтрации Калмана
 const double DT_rzv = 1.0/10.0; // Шаг РЗВ
 const double tau_eks = 0.0;//1.0/64.0; // Время экстраполяции для РЗВ
 const double perk = 3000./M_PI; // Коэффициент перевода угловых величин из рад в т.д
 extern const double tb1   = 0.1;
 const double tb2   = 0.1;
 const double K01 =0.01;// 0.01;//Коэффициент динамической ошибки качки (0.003)

 //Координаты АС в ПСК (константы для данного проекта корабля)
 const double Xac = 0;
 const double Yac = 0;
 const double Hac = 0;

// Координаты АУ в ПСК
const double Xay = 0;
const double Yay = 0;
const double Hay = 0;

//------Вероятность по дальностям-------
double masppp[100];
double obrmasppp[100];
double masVpo[100];
//------Вероятность по дальностям-------
/*
//76 mm
double sigmaA1x = 1.0;// %
double sigmaA1y = 0.4; // %
double sigmaA1z = 0.0; // %
double sigmavbn = 0.5; //%
*/
/*
//100 mm
double sigmaA1x = 0.38;//0.6; // %
double sigmaA1y = 0.38; // %
double sigmaA1z = 0.0; // %
double sigmavbn = 0.33; //%
*/

//130 mm

double sigmaA1x = 1.0; // %
double sigmaA1y = 0.4; // %
double sigmaA1z = 0.0; // %
double sigmavbn = 0.4; //%

/*
double sigmaA1x = 0.0; // %
double sigmaA1y = 0.0; // %
double sigmaA1z = 0.0; // %
double sigmavbn = 0.0; //%
*/
double Osch1;

double Ax,Ay,Az;

double dob = 0.0001;

//using std::cout; using std::endl;

//ofstream ofs;
 double	ti,//Текущее время
		dti, dttb1, dttb2,//временные интервалы
		qbc,//Курс свой
		qbc_,//производная qbc по ti
		sqbc,cqbc,//Синус и косинус курса своего
		shipVX,shipVY,shipVH,//Составляющие скорости своего корабля
		shipX,shipY,shipH,//Координаты своего корабля
		shipX0,shipY0,shipH0,
		psi,//Килевая качка
		psi_,//производная psi по ti
		teta,//Бортовая качка
		teta_,//производная teta по ti
		qssk,//Измеренный курс свой
		psisk,//Измеренная килевая качка
		psisk_AY,//Измеренная килевая качка
		tetask,//Измеренная бортовая качка
		Hcsk,//Измеренная высота своего корабля
		Vcxsk,Vcysk,Vchsk;//Измеренные составляющие скорости своего корабля
 unsigned long ixf1,ixf2,ixf3;//Используются при вычислении фаз
 double fa1,fa2,fa3;//Фазы рысканья,бортовой и килевой качек

 const double M_2PI = 6.28318530717958647693;//2*Пи
 const double GRAD_TO_RAD = 0.01745329251994329577;//Коэффициент перевода
														 //угловых величин из градусов в радианы
 const double RAD_TO_GRAD = 57.2957795130823208768;//Коэффициент перевода
														 //угловых величин из радианов в градусы
 const double ZERO = 1.0E-10;

 double delta;//Порог
 const unsigned long iy1  = 65539L;
 const double        iyf1 = 0.4656613e-9;


 unsigned long ixf;//Нечетное число, содержащее не более 9 значащих цифр
 double fa;//Искомое число
	 double  fii,fi0i,Qz,temi,tem0i,lambi;

	//Переменные, использующиеся для вычисления гауссовых случайных чисел
   //	unsigned long ixm[17];

	double slm[20];//[17];

  //	double xkg,ykg,hkg;
 double SIGMAd, SIGMAq, SIGMAe;

// double Xzy,Yzy,Hzy;//Наблюдаемые координаты цели
//Глобальные переменные для решения задачи встречи
 double y[k],yp[k];
 double y_rzw1[k],yp_rzw1[k];

 double bl,//Широта Ленинграда
				om,//
				g0,//Ускорение свободного падения
				r3,//Радиус Земли
				a,//Нормальная скорость звука
				t0,//Начальная виртуальная температура
				tkn,//Градиент начальной виртуальной температуры
				r;//Газовая постоянная
 double e1,//Точность интегрирования
				dd0,//Точность для сведения баланса по дальности, (м)
				df0,dq0,//Точности для сведения баланса по углам, (рад)
				h01;//Шаг интегрирования

//Параметры снаряда
 double //dk,//Диаметр миделя
		vb,//Начальная скорость
				qb;//Вес

 //******************************************************************************
 // Изменение № 12
 //******************************************************************************
 double ll,lm,lmz,lmg,ibx,hg0;
 double cbz,cbz1;
 double w0;
 //******************************************************************************
 // Конец изменения № 12
 //******************************************************************************

 double cbxm[31],cbym[31],fim[31];
 double aks[26]; //16.12.2010

	//******************************************************************************
 // Изменение № 13
 //******************************************************************************
	double cbzm[31];
 //******************************************************************************
 // Конец изменения № 13

	double coc,cl,sl,omx,omy,omz,
			qs1,qc1,sf,cf,sf1,cf1,fi,q,qs,qc, //cf1,sf1 - синусы и косинусы угн и увн
		qs1_rzw1,qc1_rzw1,sf_rzw1,cf_rzw1,sf1_rzw1,
		cf1_rzw1,fi_rzw1,q_rzw1,qs_rzw1,qc_rzw1,
				wet,aw,saw,caw,salf,calf,
				bky7,vbb,tp,dd1;
	double cbx,cby;
	double Xayk_,Yayk_,Hayk_; //20101115
	double Yzn,Xzn,Hzn
				 ,Xayk,Yayk,Hayk
				 ,Vyz,Vxz,Vhz;

 // int isch;//Счетчик тактов обращения к RZW
	int nd;//Размерность матрицы коэффициентов аппроксимации для деривации //16.12.2010
	int k1;//Число элементов в массивах коэффициентов согласования по координатам
	int k51;//Используется (как число элементов в массиве) в расчете коэффициента аэродинамического аксиального демпфирующего
								//момента (mmx) в зависимости  от числа Маха (va)
	//Составляющие скорости ветра:wx-продольный ветер, wz-поперечный ветер
	double wx,wy,wz,
			 wx_,wy_,wz_;//добавки к ним от производных координат АУ
	//Ошибки по продольному и поперечному ветру
	double sigmawx,sigmawz;




 double *Xz, *Yz, *Hz, *VXz, *VYz,*VHz, *tm;
 double *Xz1, *Yz1, *Hz1; //координаты цели без учета параваксов
 double t_rzw1,t00_rzw1;
// double  *VXzg, *VYzg,*VHzg;

// положение центра корабля (привязка КГСК)
// double *parrShip_X, *parrShip_Y,*parrShip_H,*parrShip_VX_True, *parrShip_VY_True, *parrShip_VH_True
//,*parrShip_VX_Zv, *parrShip_VY_Zv, *parrShip_VH_Zv;

	 //Кривые Гладковского
 double Kr_Glad[16];
 double Mas_Dal[16];

 //Функция расчета гауссового случайного числа с математическим ожиданием,
//равным 0 и дисперсией, равной 1

// unsigned long ix;//Нечетное число, содержащее не более 9 значащих цифр
 //double slu;//Искомое число


//===============================================================================
//Параметры своего корабля
 double	Vc,//Модуль скорости своего корабля
		Qco,//Генеральный курс
		Qcm,//Амплитуда рысканья
		PSIm,//Амплитуда килевой качки
		TETAm,//Амплитуда бортовой качки
		Tp,//Период рысканья
		Tb,//Период бортовой качки
		Tk,//Период килевой качки
		Hm,//Амплитуда вертикальной качки
		K0,//Коэффициент ошибки измерителя качки
		K11;// 1.+K0




// ФУНКЦИЯ РАСЧЕТА ЭФФЕКТИВНОСТИ
// INPUT:
// VAlDT_filtr - шаг фильтрации
//  QuantIspit  - к-во сипытаний
//  QuantShells - к-во снарядов
//  VAlCalibro  - калибр
//  VAlDn0 - дальность начала стрельбы
//  VAlDk  - дальность конца стрельбы
//  HAntenna  - высота антенны
//  Far_2D   - антенна
//  EtalonSign - эталонный сигнал, применяемый для расчета амплитуды сигнала  :
		 // амлитуда
	 // EtalonSign.mEtalonAmp;
	 //дальность
	 // EtalonSign.mEtalonDist;
	 //ЭПР
	 // EtalonSign.mEtalonEPR;
	 //СКО внутр шума суммарной диаграммы 5П10
	 // EtalonSign.mNoiseSKZ_5P10;
	 // СКЗ разброса коэффиц усиления суммарной диаграммы 5П10
	 // EtalonSign.mEtalonSigAmplFact_5P10;

 // TargBear - данные пеленга цели
		//  угол пеленга цели в ГСК
	// mBearing ;
	// угол курса цели в ГСК
	// mTargCourse ;
	// угол между ветором скорости цели и напрвлением в зенит
// mTargHorizAng ;
	//- скорость
	// mV;
	// дальность,
	// mR ;
	// высота
	// mH ;
	//  сигма w
	// mMy ;

	// ЭПР
	// mTargEPR;
	// тип цели
	// mTargType;
	// VAlSigSin - ошибка измерения углов качки и рыскания (СИНС)


//  OUTPUT:
// *pvalProb-    вероятность поражения
// *pvalDistBeginSopr - дальность начала разрешения цели
int calcEffect(wchar_t *wchOutFold, TFar_2D Far_2D,  const double HAntenna, TTarget Targ , TInitTargData InitTargData
	 ,const double PowerPrd, const double KYPrd, TEtalonSign EtalonSign, const double VAlDT_filtr
		 , const int QuantIspit, const int QuantShells,const double VAlCalibro
	, const double VAlDn0,  const double VAlDk, const double VAlSigSins, const bool VAlBSkaliga
		, double *pvalProb, double *pvalProb_Gladk, double *pvalDistBeginSopr, int *pQuantShots)
{


	ResetFun();


	*pvalDistBeginSopr = InitTargData.mR;
	TInitTargData InitTargData0 = InitTargData;

	// линейная ФАР в вертикальном направлении
	bool bLat = false;

	double valSigE = -1.,valSigQ  = -1.0;


		double valAntpPhaze =0;   // !!!!!!!!!!!!!!!!!!!!! ИСПРАВИТЬ!!!
		///

	*pvalDistBeginSopr = Far_2D.calc_TwoTargsZahvatDist(InitTargData0.mH, Targ.mTargEPR
  , PowerPrd,  KYPrd, EtalonSign, HAntenna
  ,&bLat  ,&valSigE ,&valSigQ  ) ;



	TFar Far0(Far_2D, true);
	TFar Far(Far0, 4) ;
	const double VAlAMDiagrWidth = Far.mFaceta.findDiagrWidth(); // ширина дигнраммы напрвл АМ

	double valDnCur = 1000000000.;

	if ((*pvalDistBeginSopr) < 0.)
	{
		*pvalDistBeginSopr = 0;
		*pvalProb = 0.;
		*pvalProb_Gladk= 0.0;
		return -1;
	}
	 InitTargData0.mR =  MIN__(InitTargData0.mR , *pvalDistBeginSopr);
	if(VAlDn0> (InitTargData0.mR -  InitTargData0.mV))
	{
		valDnCur=  InitTargData0.mR -  InitTargData0.mV;
		if(valDnCur< VAlDk)
		{
			*pvalProb = 0.;
			*pvalProb_Gladk= 0.0;
			 return -2;
		}
	}
	const double VAlDn = MIN__(valDnCur , VAlDn0)  ;

		//Установка параметров своего корабля
			double targVX0,targVY0,targVH0,shipVX0,shipVY0,shipVH0;
	Vc= 10.;//Модуль скорости своего корабля
	Qco= 0.0*GRAD_TO_RAD;//Генеральный курс
	Qcm= 1.0*GRAD_TO_RAD;//Амплитуда рысканья
	PSIm=	5.0*GRAD_TO_RAD;//Амплитуда килевой качки
	TETAm= 12.0*GRAD_TO_RAD;//Амплитуда бортовой качки

	Tp= 5.0;			//Период рысканья
	if(Tp>=ZERO) 	Tp=	M_2PI/Tp;
	else 					Tp= 0.0;
	Tb= 5.0;			//Период бортовой качки
	if(Tb>=ZERO)  Tb=	M_2PI/Tb;
	else 					Tb= 0.0;
	Tk= 3.0;			//Период килевой качки
	if(Tk>=ZERO)  Tk=	M_2PI/Tk;
	else 					Tk= 0.0;

	Hm=	0.0;//1.0;//Амплитуда вертикальной качки
	K0=	0.0004;//0.00015;//0.001;//Коэффициент ошибки измерителя качки
	K11= 1.0 + K01;

	//Начальные значения координат и скоростей своего корабля
	shipX0= shipY0= shipH0=0.0;
	shipVX0= Vc*sin(Qco);
	shipVY0= Vc*cos(Qco);
	shipVH0= 0.0;



	// пересчет начала сопровождение
	 // надо найти
	 double arrVOtn[3] = {0.};
		arrVOtn[0]= Targ.mTraject.marrVectSostGSK[3] - shipVX0;
		arrVOtn[1]= Targ.mTraject.marrVectSostGSK[4] - shipVY0;
		arrVOtn[2]= Targ.mTraject.marrVectSostGSK[5];
		double vala = Norm3 (arrVOtn) * Norm3 (arrVOtn);
		double valb = 2. * ScalProduct(Targ.mTraject.marrVectSostGSK , arrVOtn, 3) ;
		double valc =   Norm3 (Targ.mTraject.marrVectSostGSK) * Norm3 (Targ.mTraject.marrVectSostGSK) - valDnCur* valDnCur;
		TComp cmpx1, cmpx2;
		int irez00 =  SolvEq2(vala,valb,valc, cmpx1, cmpx2);
		if (irez00 >2)
		{
			 ShowMessage(L"Error1");
			 *pvalProb = 0.;
			 *pvalProb_Gladk= 0.0;
			 return -1;
		}
		double valTBegin = MIN__(cmpx1.m_Re, cmpx2.m_Re);
		TTarget  Targ0 = Targ;
		double arrPositionOtn[3] ={0.};
		arrPositionOtn[0] =  Targ.mTraject.marrVectSostGSK[0] + arrVOtn[0]* valTBegin ;
		arrPositionOtn[1] =  Targ.mTraject.marrVectSostGSK[1] + arrVOtn[1]* valTBegin ;
		arrPositionOtn[2] =  Targ.mTraject.marrVectSostGSK[2] + arrVOtn[2]* valTBegin ;
		InitTargData0.mBearing =  atan2( arrPositionOtn[0], arrPositionOtn[1]);


		Targ0.mTraject = TTraject (0.,  Targ.mTraject.mSigW, InitTargData0 );

	// расчет СКЗ ошибок измерения углов цели на эталонной дальности и угле места и курсовом угле цели 5 мрд
	const double VAlSigEtalonRSM_Eps = Far_2D.calcSigEtalonRSM_Eps(EtalonSign, 0.005);
	const double VAlSigEtalonRSM_Bet = Far_2D.calcSigEtalonRSM_Bet(EtalonSign, 0.005);
	const double VAl_NWaveEtalon = Far_2D.calcNWaveEtalon ( EtalonSign);
	/// конец расчета



	RZW_Input TRZW_Input;
	RZW_Output TRZW_Output;


	Filter_Input TFilter_Input;
	Filter_Output TFilter_Output;

	unsigned long i,n;
	int it0,itk,iisch,ksn,itt,iitt,nn,ii,/*Ns,*/nn1;
	double	x_kg,y_kg,h_kg,Dwstr0,
					Dwstr,t00,tk,
					vx_kg,vy_kg,vh_kg;
	double 	*Xf,*Yf,*Hf,*VXf,*VYf,*VHf,
					*tmf,*massq,*mascq,*masfi,*masx,*masy,
					*mashh,*mastp,*mastr,*masDy;


	double fii0, fi0i0, Qz0, temi0, tem0i0, lambi0;

	double Vpo,Vpc,Ppp,Vpg,Ppg;
	double Wes;

	double Dndelta, Dkdelta;

	double ttekk,ts;
	double tpst1;
	double vb0,Ssigmawx,Ssigmawz,sigmap,sigmavb, Ssigmaq, Ssigmafi;
	double sqssk,cqssk,spsisk,cpsisk,stetask,ctetask;
	double Xay_ksk,Yay_ksk;//,Xayk,Yayk,Hayk;
	double sqs,cqs,spsi,cpsi,steta,cteta;
	double /* tpRZW,*/ tp00, tpkk;

	int prisnak;
	double sqrzv1,cqrzv1, firzv1,tprzv1,dyrzv1;

	//Установка исходных данных
	//VAlDn = 10000.0;//Начальная дальность стрельбы
	//Dk = 6000.0;//Конечная дальность стрельбы
	delta = 50.0;//Порог
	//Kr = 500;//Количество реализаций

	//Ошибка, учитывающая деформацию корпуса корабля
	Osch1=0.00041;//(2.7/60.0)*GRAD_TO_RAD; //АП
	double Osch2=0.00041;//(18.5/60.0)*GRAD_TO_RAD; //АУ

	//Параметры снаряда
	// dk=0.13;//0.13;//0.0762;//0.1;//Диаметр миделя
	if(fabs(VAlCalibro-0.1)<= 0.000001)
	{
		//Kcn = 45;//174;//Количество снарядов
		ts= 0.75;//1.0;//0.8;//Темп стрельбы
		k1= 31;//Число элементов в массивах коэффициентов согласования по координатам
		k51= 31;
		lm= 0.522; // Характерная длина снаряда, м
		lmz= 0.21871;// расст.от основания головной части снаряда до ц.масс <м>
		lmg= 0.12229; // Длина головной части снаряда, м
		hg0= lmz+0.57*lmg-0.16* VAlCalibro;
		ibx= 0.02163;// аксиальный момент инерции снаряда,кг*м*м
		w0= 1843.1;  // начальная аксиальная угловая скорость в рад/с:
	}

	if(fabs(VAlCalibro-0.13)<= 0.000001)
	{
		// Kcn = 180;//Количество снарядов
		ts= 2.0;//0.67;//1.0;//2.0;//0.67;//2.0;//0.7;//Темп стрельбы
		k1= 31;//Число элементов в массивах коэффициентов согласования по координатам
		k51= 31;
		lm= 0.6833; // Характерная длина снаряда, м
		lmz= 0.3295;// расст.от основания головной части снаряда до ц.масс <м>
		lmg= 0.1263; // Длина головной части снаряда, м
		hg0= lmz+0.57*lmg-0.16 * VAlCalibro;
		ibx= 0.07896;// аксиальный момент инерции снаряда,кг*м*м
		w0= 1643.3;  // начальная аксиальная угловая скорость в рад/с:
	}

	if( fabs(VAlCalibro - 0.0762) <= 0.000001)
	{
		// Kcn = 152;//Количество снарядов
		ts= 0.5;	//Темп стрельбы
		k1= 16;	//Число элементов в массивах коэффициентов согласования по координатам
		k51= 16;
	}

	double SIGMAq0, SIGMAe0;
	// Установка шумов наблюдения
	SIGMAd= 10.0;//5.0; //5.0;

	// Установка ошибок рассеивания снаряда
	Ssigmaq= 1.0;//2.0*1.48;//1.0;//1.5;//0.6; //1.5; // в миллирадианах
	Ssigmafi= 1.0;//2.0*1.48;//1.0;//1.5; //1.5; // в миллирадианах

	//Установка метеоошибок
	Ssigmawx= 2.3;//в метрах
	Ssigmawz= 2.3;//в метрах
	sigmap= 0.9; //10;//в процентах
	sigmavb= 0.1;//0.5; //10;//в процентах

	//Присвоение начальных значений переменным, использующимся при вычислении
	//случайных фаз рысканья и качек своего корабля
	ixf1 = 31745; ixf2 = 24357; ixf3 = 12349;


	//Установка параметров цели
	double Par=0.0;
	const double VAlE0=asin(InitTargData0.mH / InitTargData0.mR);


	//Установка параметров для решения задачи встречи
	//Параметры снаряда
	if( fabs(VAlCalibro - 0.1) <= 0.000001)
	{
		vb0= vb= 880.0;	//Начальная скорость
		qb= 15.6;//Вес
		for(i=0;i<k1;i++)
		{
			cbxm[i]=cbxm100[i];
			cbym[i]=cbym100[i];
			cbzm[i]=cbzm100[i];
			fim[i]=fim100[i];
		}
	}
	if( fabs(VAlCalibro - 0.13) <= 0.000001)
	{
		vb0=vb=850.;//Начальная скорость
		qb=33.4;//Вес
		for(i=0;i<k1;i++)
		{
			cbxm[i]=cbxm130[i];
			cbym[i]=cbym130[i];
			cbzm[i]=cbzm130[i];
			fim[i]=fim130[i];
		}
	}
	if( fabs(VAlCalibro - 0.0762) <0.000001)
	{
		vb0=vb=980.;//Начальная скорость
		qb=5.9;//Вес
		for(i=0;i<k1;i++)
		{
			cbxm[i]=cbxm76[i];
			cbym[i]=cbym76[i];
			fim[i]=fim76[i];
		}
	}

	bky7= 1.0;//Плотность воздуха
	coc= 0.474 * VAlCalibro * VAlCalibro/qb;	// Баллистический коэффициент

	bl=M_PI/3.;//Широта Ленинграда
	om=0.7292/10000.;
	g0=9.81;//Ускорение свободного падения
	r3=6371.*1000.;//Радиус Земли
	a=340.8;//Нормальная скорость звука
	t0=288.9;//Начальная виртуальная температура
	tkn=0.006328;//Градиент начальной виртуальной температуры
	r=29.27;//Газовая постоянная

	e1=0.000001;//Точность интегрирования
	dd0=0.5;//Точность для сведения баланса по дальности, (м)
	df0=dq0=0.001;//Тoчность для сведения баланса по углам, (рад)
	h01=0.0005;//Шаг интегрирования
	dd1=0;

	cl=cos(bl);sl=sin(bl);//Синус и косинус широты Ленинграда
	fi= VAlE0 +0.1;//Начальный угол места наведения снаряда
	q = InitTargData.mBearing;//Начальный угол наведения снаряда в горизонтальной плоскости

	sf=sin(fi); cf=cos(fi);//Синус и косинус угла места наведения снаряда
	qs=sin(q); qc=cos(q);//Синус и косинус угла наведения снвряда в горизонт. плоскости
	qs1=qs; qc1=qc;

	wet=0.0;//Модуль скорости ветра
	aw=0.0;//Азимут ветра
	caw=cos(aw); saw=sin(aw);
 	vbb=vb;

	omy= 2.0*om*sl;

	wy=0.0;

	if( fabs(InitTargData0.mV) <= 0.0000001) n=5000;
	else //n=floor(D0/(V0*DTi));
	n=floor((InitTargData0.mR-VAlDk)/(InitTargData0.mV  * DTi));

	nn=floor((InitTargData0.mR - VAlDk)/(InitTargData0.mV  * VAlDT_filtr));
//  printf("\n n=%d nn=%d",n,nn);
	nn1=floor((InitTargData0.mR - VAlDk)/(InitTargData0.mV  * DT_rzv));

	// НАхождение длины отрезка траектории цели от момента обнаружения до конца обстрела
	int n_new = 0, nn_new = 0, nn1_new = 0;
	//double arrX0[6] ={0.};
	//TargBear0.raschet_nach_coord (arrX0) ;
	TURPointZ PNtPos(Targ0.mTraject.marrVectSostGSK[0], Targ0.mTraject.marrVectSostGSK[1], Targ0.mTraject.marrVectSostGSK[2]);
	//TURPointZ PNtVelo (arrX0[3] - shipVX0, arrX0[4] - shipVY0, arrX0[5] - shipVH0);
	TURPointZ PNtVelo (Targ0.mTraject.marrVectSostGSK[3], Targ0.mTraject.marrVectSostGSK[4], Targ0.mTraject.marrVectSostGSK[5] );
	 double pT [2] = {0.};
	 TURPointZ pPnt[2];
	int irez = TURPointZ::findIntertsectLine_And_Sphere( PNtPos, PNtVelo,  VAlDk
	, pT, pPnt);
	if ((irez != 2)|| (pT[0] <= 0.) )
	{
	 ShowMessage(L"Ошибка в исходных дагнных пеленга цели");
	 return -10;
	}

	const double VAlLengTrack = TURPointZ::Dist(PNtPos ,pPnt[0]) ;
	if( fabs(InitTargData0.mV) <= 0.0000001)
	{
	 n_new=5000;
	}
	else
	{
	n_new = VAlLengTrack/(TURPointZ::Norm(PNtVelo)  * DTi);
	}

	nn_new= VAlLengTrack/(TURPointZ::Norm(PNtVelo)  * VAlDT_filtr);
	nn1_new= VAlLengTrack/(TURPointZ::Norm(PNtVelo)  * DT_rzv);

	///


	//Выделение памяти для массивов
	SZn  = n  +1;
	SZnn = nn +1;
	SZnn1= nn1+1;

	Xz = new double [n+2];  // +
	Yz = new double [n+2]; // +
	Hz = new double [n+2];  // +
	Xz1 = new double [n+2]; // +
	Yz1 = new double [n+2];  // +
	Hz1 = new double [n+2];  // +
	VXz = new double [n+2]; // +
	VYz = new double [n+2]; // +
	VHz = new double [n+2]; // +
 //	VXzg = new double [n+2];
 //	VYzg = new double [n+2];
 //	VHzg = new double [n+2];
	tm = new double [n+2];  // +

	Xf = new double [nn+1]; // +
	Yf = new double [nn+1];  // +
	Hf = new double [nn+1];  // +
	VXf = new double [nn+1]; // +
	VYf = new double [nn+1]; // +
	VHf = new double [nn+1]; // +
	tmf = new double [nn+1];//Время фильтрацц   // +

	massq = new double [nn1+1]; //sin ПУГН   // +
	mascq = new double [nn1+1]; // cos ПУГН  // +
	masfi = new double [nn1+1]; //ПУВН     // +
	mastp = new double [nn1+1];//Полетное время  // +
	mastr = new double [nn1+1];//Время решения ЗВ  // +
	masDy = new double [nn1+1];//Упрежденная дальность  // +

	masx = new double [nn1+1];//для экстраполяции углов   // +
	masy = new double [nn1+1];//для экстраполяции углов   // +
	mashh = new double [nn1+1];//для экстраполяции углов // +


	double *parrTrueShipGSK_X = new double [n+2];
	double *parrTrueShipGSK_Y = new double [n+2];
	double *parrTrueShipGSK_Z = new double [n+2];
	iisch=0;
	ti=0.0;
	it0=itk=0;
	tp=0;

	Dndelta = VAlDn + delta;
	Dkdelta = VAlDk - delta;

	dti=DTi;
	dttb1 = dti/tb1;
	dttb2 = dti/tb2;

	sigmawx=sigmawz=0.0;

	Dwstr0 = InitTargData0.mR;

	//Начальные значения координат и скоростей своего корабля
	shipX=shipX0; shipY=shipY0; shipH=shipH0;
	shipVX=shipVX0; shipVY=shipVY0;shipVH=shipVH0;




  //    подгот цикл
	fncPreviousArrangments( Targ0,  VAlCalibro,VAlDn, VAlDk, n,  dti, t00, tp00,  tk, tpkk) ;
	ttekk= dti*(double)n;
	*pQuantShots = MIN__(int((tk - t00 + 0.0000001)/ ts) , QuantShells);
	if ((*pQuantShots) <= 0)
	{
		*pQuantShots = -1;
		*pvalProb = 0.;
		*pvalProb_Gladk= 0.0;
		 return -3;
	}
	///


	if(fabs(InitTargData0.mV)<= 0.00000001)
	{ t00=0.2; tk=ttekk=3.2;
	}

	t00_rzw1= t00;

	Ppp= 0.0;
	Ppg = 0.;

//------Вероятность по дальностям-------
	for(i=0;i<100;i++) {masppp[i] = 0.0; obrmasppp[i] = 0.0;}
//------Вероятность по дальностям-------
			int quantBuffCols1 = 4;
			int quantPoints = (*pQuantShots) ;
			double * parrBuff1 = (double *) malloc(quantBuffCols1  * quantPoints *  sizeof(double));
			memset(parrBuff1, 0,quantBuffCols1  * sizeof(double));
			wchar_t *wcharrFileNames1 = (wchar_t *) malloc(quantBuffCols1 * 30 *sizeof(wchar_t));
			memset(wcharrFileNames1, 0,quantBuffCols1 * 30 *sizeof(wchar_t));
			wcscpy(&wcharrFileNames1[0],L"NumShot");
			wcscpy(&wcharrFileNames1[30],L"DistTochkiVstr");
			wcscpy(&wcharrFileNames1[30 * 2],L"Probab_V_Tochke");
			wcscpy(&wcharrFileNames1[30 * 3],L"Promach");

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// ЦИКЛ ПО КОЛИЧЕСТВУ РЕАЛИЗАЦИЙ  ////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(int ikr = 1; ikr <= QuantIspit; ikr++)
	{
		TTarget TargCur = Targ0;
		workFilter = 0;
		fi= VAlE0 +0.1;//Начальный угол места наведения снаряда
		q= InitTargData0.mBearing;//Начальный угол наведения снаряда в горизонтальной плоскости

		sf= sin(fi); cf= cos(fi);	//Синус и косинус угла места наведения снаряда
		qs= sin(q);  qc= cos(q);	//Синус и косинус угла наведения снвряда в горизонт. плоскости
		qs1=qs; qc1=qc;

		dd1= 0;
		tp= 0;
		tpst1= 0;

		//Начальные значения координат и скоростей своего корабля
		shipX=  shipX0;  shipY = shipY0;  shipH = shipH0;
		shipVX= shipVX0; shipVY= shipVY0; shipVH= shipVH0;

		slm[6] = getGauss(0., 1.);
		slm[7] = getGauss(0., 1. );
		slm[8] = getGauss(0., 1. );
		slm[9] = getGauss(0., 1. );


		sigmawx= Ssigmawx*slm[6];
		sigmawz= Ssigmawz*slm[7];
		bky7= 1.0+sigmap*slm[8]/100.0;
		vb=   vb0+ sigmavb*vb0*slm[9]/100.0;

		ti=-dti; itt=-1; iitt=-1;

		//Получение равномерно распределенных на [0,2pi] фаз своего корабля
		//для курса своего

		fa1 = getRand01() * M_PI * 2.;
		fa2 = getRand01() * M_PI * 2.;
		fa3 = getRand01() * M_PI * 2.;



		Vpo=0;
		ksn=0;
		Vpg=0;
		for(i=0;i<100;i++) masVpo[i] = 0.0; // вероятность по дальностям

		vbb= vb;

		ii= -1;
		for(i=0;ti<=ttekk;i++)
		{
			ti=ti+dti;
			raschet_coord_swoego_corablja(ti, VAlSigSins); // ф-я расчета места и скорости своего корабля

			TargCur.recalcTrajPoint(ti) ;
			 x_kg  = TargCur.mTraject.marrVectSostGSK[0] - shipX;
			 y_kg  = TargCur.mTraject.marrVectSostGSK[1] -shipY;
			 h_kg  = TargCur.mTraject.marrVectSostGSK[2]  -shipH;
			 vx_kg = TargCur.mTraject.marrVectSostGSK[3] -shipVX;
			 vy_kg = TargCur.mTraject.marrVectSostGSK[4] -shipVY;
			 vh_kg = TargCur.mTraject.marrVectSostGSK[5]  -shipVH;

			//Вычисление координат артустановки в КГСК
			sqs= sin(qbc); spsi= sin(psi); steta= sin(teta);
			cqs= cos(qbc); cpsi= cos(psi); cteta= cos(teta);

			Xay_ksk=  Xay*cteta+Hay*steta;
			Yay_ksk= -Xay*spsi*steta+Yay*cpsi+Hay*spsi*cteta;

			Hayk= -Xay*cpsi*steta-Yay*spsi+Hay*cpsi*cteta;
			Xayk=  Xay_ksk*cqs+Yay_ksk*sqs;
			Yayk= -Xay_ksk*sqs+Yay_ksk*cqs;

			ii++;
			if(ii> SZn)
			{
			 break;   // ГАРАНТИЯ ПОПАДАНИЯ В ГРАНИЦЫ МАССИВОВ
			}
			Xz[ii] = x_kg;
			Yz[ii] = y_kg;
			Hz[ii] = h_kg;
			VXz[ii]= vx_kg;
			VYz[ii]= vy_kg;
			VHz[ii]= vh_kg;
			Xz1[ii]= x_kg;
			Yz1[ii]= y_kg;
			Hz1[ii]= h_kg; // место цели в КГСК
		//	VXzg[ii]= TargCur.mTraject.marrVectSostGSK[3];
		//	VYzg[ii]= TargCur.mTraject.marrVectSostGSK[4];
		//	VHzg[ii]= TargCur.mTraject.marrVectSostGSK[5];  // скорость цели в ГСК
			tm[ii]= ti;     // время привязки

			// истинные параметры движения корабля
			 parrTrueShipGSK_X[ii]  = shipX;
			 parrTrueShipGSK_Y[ii]  = shipY;
			 parrTrueShipGSK_Z[ii]  = shipH;
			// parrShip_VX_True[ii]  = shipVX;
			// parrShip_VY_True[ii]  = shipVY;
			// parrShip_VH_True[ii]  = shipVH;
			// измеренная скорость корабля
			// parrShip_VX_Zv[ii]  = Vcxsk ;
			// parrShip_VY_Zv[ii]  = Vcysk ;
			// parrShip_VH_Zv[ii]  = Vchsk ;
		}//конец цикла по i

		//Начальные значения координат и скоростей своего корабля
		shipX = shipX0;  shipY = shipY0;  shipH = shipH0;
		shipVX= shipVX0; shipVY= shipVY0; shipVH= shipVH0;

		ti= -VAlDT_filtr;

		// графики изменения СКЗ ошибок измерения углов от дальности
			int quantBuffCols = 4;
			double * parrBuff = (double *) malloc(quantBuffCols  * nn *  sizeof(double));
			memset(parrBuff, 0,quantBuffCols  * sizeof(double));
			wchar_t *wcharrFileNames = (wchar_t *) malloc(quantBuffCols * 30 *sizeof(wchar_t));
			memset(wcharrFileNames, 0,quantBuffCols * 30 *sizeof(wchar_t));
			wcscpy(&wcharrFileNames[0],L"D");
			wcscpy(&wcharrFileNames[30],L"SIGMAe");
			wcscpy(&wcharrFileNames[30 * 2],L"SIGMAq");
			wcscpy(&wcharrFileNames[30 * 3],L"brez0");

		///


		for(i=0;i<nn;i++)		//Начало цикла по времени для фильтра Калмана
		{
			if ( 145 == i)
			{
			 int iiii =0;
			}
			ti+= VAlDT_filtr;
			//Обращение к функции расчета координат и скоростей своего корабля
		 //	raschet_coord_swoego_corablja(ti, VAlSigSins);

			int NOM1,NOM2;
			NOM1=(int)(ti/dti+dob); NOM2=NOM1+1;

			if(NOM2> SZn)
			{
			 break;   	// ГАРАНТИЯ ПОПАДАНИЯ В ГРАНИЦЫ МАССИВОВ
			}
			double arrVSTargKGSK_True[3] ={0.};
		 /*	arrVSTargKGSK_True[0]= (Xz1[NOM2]-Xz1[NOM1])*ti/dti+(Xz1[NOM1]*tm[NOM2]-Xz1[NOM2]*tm[NOM1])/dti;
			arrVSTargKGSK_True[1]= (Yz1[NOM2]-Yz1[NOM1])*ti/dti+(Yz1[NOM1]*tm[NOM2]-Yz1[NOM2]*tm[NOM1])/dti;
			arrVSTargKGSK_True[2]= (Hz1[NOM2]-Hz1[NOM1])*ti/dti+(Hz1[NOM1]*tm[NOM2]-Hz1[NOM2]*tm[NOM1])/dti;*/


			double valDelT = ti - ( double(NOM1))* dti;
			arrVSTargKGSK_True[0]= fncLinExtrapolation(Xz[NOM1], Xz[NOM1 + 1],dti,valDelT) ;
			arrVSTargKGSK_True[1]= fncLinExtrapolation(Yz[NOM1], Yz[NOM1 + 1],dti,valDelT) ;
			arrVSTargKGSK_True[2]= fncLinExtrapolation(Hz[NOM1], Hz[NOM1 + 1],dti,valDelT) ;;

			// fncLinExtrapolation(Xz1[NOM1, Xz1[NOM1 + 1],dti,valDelT) ;

			//Обращение к функции расчета наблюдаемых координат, получены в КГСК
			double arrVSTargKGSK_Zv[3] = {0.};// наблюдаемые координаты цели в КГСК

			double valDCur =  Norm3(arrVSTargKGSK_True);

		const double VAlAntNormalAng = 0.;


		bool brez0 =  Far_2D.calc_SKZ_LAT(valDCur ,arrVSTargKGSK_True[2]
		,HAntenna,  VAlAntNormalAng, Targ.mTargEPR
		,EtalonSign ,PowerPrd, KYPrd,  valAntpPhaze
		,&SIGMAe ,&SIGMAq);
		 //	SIGMAe= 0.0005 ;
		 //	SIGMAq= 0.0005 ;
		 raschet_coord_zeli_nabl(arrVSTargKGSK_True, arrVSTargKGSK_Zv);


		//
			 parrBuff[ i * quantBuffCols    ] = valDCur;
			 parrBuff[ i * quantBuffCols + 1] = SIGMAe;
			 parrBuff[ i * quantBuffCols + 2] = SIGMAq;
			 parrBuff[ i * quantBuffCols + 3] = ((double)brez0 );

		///
			double valSigEps = MAX__(SIGMAe , 0.00075);
			double valSigBet = MAX__(SIGMAq , 0.00075);
		 //	double valSigEps = 0.0005;
		 //	double valSigBet = 0.0005 ;

			TFilter_Input.tk=ti;//Время
			TFilter_Input.Xkg = arrVSTargKGSK_Zv[0];
			TFilter_Input.Ykg = arrVSTargKGSK_Zv[1];
			TFilter_Input.Hkg = arrVSTargKGSK_Zv[2];
			TFilter_Input.sbet= valSigBet * valSigBet;//disp;//0.25e-06;
			TFilter_Input.seps=  valSigEps * valSigEps; //disp;//0.25e-06;
			TFilter_Input.sD=100;
			TFilter_Input.Vcx= Vcxsk;
			TFilter_Input.Vcy= Vcysk;
			TFilter_Input.Vch= Vchsk;

		 //	filtr_ukf(&TFilter_Input,&TFilter_Output);

			filtr_ukf(&TFilter_Input,&TFilter_Output, VAlBSkaliga, Targ0.mTraject.mSigW * Targ0.mTraject.mSigW, VAlSigSins*VAlSigSins);
			//Оценки координат экстраполируем на tau_eks
			Xf[i] = TFilter_Output.X + TFilter_Output.Vx*tau_eks;
			Yf[i] = TFilter_Output.Y + TFilter_Output.Vy*tau_eks;
			Hf[i] = TFilter_Output.H + TFilter_Output.Vh*tau_eks;
			VXf[i]= TFilter_Output.Vx;
			VYf[i]= TFilter_Output.Vy;
			VHf[i]= TFilter_Output.Vh;
			tmf[i]= ti;
	 	}	//конец цикла по i
/*
		// Построение графиков для фиольтра
		double scalex =1., scaley =1000000.;
		if(wchOutFold[0] != 0)
		{
		double arrscaley[4] = {0.};
		 arrscaley[0] = 1.;
		 arrscaley[1] =1000000.;
		 arrscaley[2] = 1000000.;
		 arrscaley[3] = 100.;
		wchar_t wchFoldName[400] = {0};
		wcscpy(wchFoldName , wchOutFold) ;
		wcscat(wchFoldName, L"\\");

		for (int iq = 1; iq < quantBuffCols; iq++)
	 {
		TYrWriteShapeFile::WriteOneReport(wchFoldName  // путь к папке
									, parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
									,quantBuffCols// - к-во переменных о корорых накоплена информация в буфере
									,nn//  - к-во точек
									,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30 // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,iq   // номер переменной по оси Y
									, arrscaley[0]  //  масштаб по оси X
									, arrscaley[iq]  // масштаб по оси Y
								   ) ;
	 }

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,0., 100000.
	 ,0.,50000., 5.) ;
	}
	*/
	 free (wcharrFileNames);
	 free(parrBuff);
///
		//	TYrWriteShapeFile::CreateShpFile(L"E:\\30-08-17\\Temp\\Yf.shp", Yf, tmf,nn, scalex, scaley);
		//	TYrWriteShapeFile::CreateShpFile(L"E:\\30-08-17\\Temp\\VYf.shp", VYf, tmf,nn, scalex, scaley);
		//	TYrWriteShapeFile::CreateShpAxes(L"E:\\30-08-17\\Temp\\Axes.shp",-1000.,1000.
		//	 ,-10000.,10000.);


		double Xfi,Yfi,Hfi,VXfi,VYfi,VHfi;
		//Начальные значения координат и скоростей своего корабля
		shipX = shipX0;  shipY = shipY0;  shipH = shipH0;
		shipVX= shipVX0; shipVY= shipVY0; shipVH= shipVH0;

		ti=-DT_rzv;
		// ОБЪЕКТ !!!!  РАБОТЫ
		double valFi0,valQu0, valMulFi, valMulQu,valPrevFi0,valPrevQu0;
		int ichetchik1 = 0;
		// for(i=0;i<nn1-10;i++)//Начало цикла по времени для РЗВ
//-------------------------------------------------------------------------
// --RZW--RZW--RZW--RZW--RZW--RZW--RZW--RZW--RZW--RZW--RZW--RZW-RZW-RZW--RZW-
//-------------------------------------------------------------------------
		for(i=0;i<nn1;i++)	//Начало цикла по времени для РЗВ        t00-ts
		{
			ti= ti+DT_rzv;
		 //	if (ti <  (t00- 2. * ts - DT_rzv))
		//	{
		//		continue;
		//	}
			if (i == 72)
			{
				int iii = 0;
			}
			if(i>0)
			{ valFi0= valPrevFi0;
				valQu0= valPrevQu0;
				valFi0 = valPrevFi0 -0.03;
				valQu0 = valPrevQu0;
				valMulFi = 0.95;
				valMulQu = 1.0;
			}

			ichetchik1= i;


			int NOM1;
			NOM1= (int)(ti/VAlDT_filtr + dob);

			 double valDelT =  ti - VAlDT_filtr * ((double)NOM1);
			if((NOM1 + 1)> SZnn)
			{
			 break;   // ГАРАНТИЯ ПОПАДАНИЯ В ГРАНИЦЫ МАССИВОВ
			}



			Xfi= fncLinExtrapolation(Xf[NOM1], Xf[NOM1 + 1],VAlDT_filtr,  valDelT );
			Yfi= fncLinExtrapolation(Yf[NOM1], Yf[NOM1 + 1],VAlDT_filtr,  valDelT );
			Hfi= fncLinExtrapolation(Hf[NOM1], Hf[NOM1 + 1],VAlDT_filtr,  valDelT );
			VXfi= fncLinExtrapolation(VXf[NOM1], VXf[NOM1 + 1],VAlDT_filtr,  valDelT );
			VYfi= fncLinExtrapolation(VYf[NOM1], VYf[NOM1 + 1],VAlDT_filtr,  valDelT );
			VHfi= fncLinExtrapolation(VHf[NOM1], VHf[NOM1 + 1],VAlDT_filtr,  valDelT );;
			//Обращение к функции расчета координат и скоростей своего корабля
			raschet_coord_swoego_corablja(ti, VAlSigSins);
			//Вычисление координат артустановки в КГСК
			sqssk  = sin(qssk);   cqssk  = cos(qssk);
			spsisk = sin(psisk);  cpsisk = cos(psisk);
			stetask= sin(tetask); ctetask= cos(tetask);

			Xay_ksk=  Xay*ctetask+Hay*stetask;
			Yay_ksk= -Xay*spsisk*stetask+Yay*cpsisk+Hay*spsisk*ctetask;
			Hayk= -Xay*cpsisk*stetask-Yay*spsisk+Hay*cpsisk*ctetask;
			Xayk=  Xay_ksk*cqssk+Yay_ksk*sqssk;
			Yayk= -Xay_ksk*sqssk+Yay_ksk*cqssk;



			sqbc= sin(qbc);
			cqbc= cos(qbc);

			prisnak= 0;

			sigmawx= Ssigmawx*prisnak;
			sigmawz= Ssigmawz*prisnak;
			bky7 = 1.0+sigmap*prisnak/100.0;
			vb = vb0+ sigmavb*vb0*prisnak/100.0;
			wy= shipVH;

			tp=0.0;//tpst1;
			vbb=vb0;
		 double arrVSTarg_KGSK0[6] ={0.};
		 arrVSTarg_KGSK0[0] = Xfi-Xayk;
		 arrVSTarg_KGSK0[1] = Yfi-Yayk;
		 arrVSTarg_KGSK0[2] = Hfi-Hayk;
		 arrVSTarg_KGSK0[3] = VXfi ;
		 arrVSTarg_KGSK0[4] = VYfi ;
		 arrVSTarg_KGSK0[5] = VHfi;
			if(i==0)
			{	InitFi0Qu0(arrVSTarg_KGSK0,valFi0,valQu0, vb0);  /// !!!!!! БЫЛА ОШИБКА!!!! 20.10.2017
				valFi0= 0.05;
				valMulFi= 0.7;
				valMulQu= 0.9;
			}

      // формирование вектора состояния корабля в ГСК
	double arrVS_Vessel_GSK00[6] ={0.};
	arrVS_Vessel_GSK00[0] = shipX ;
	arrVS_Vessel_GSK00[1] = shipY ;
	arrVS_Vessel_GSK00[2] = shipH ;
	arrVS_Vessel_GSK00[3] = Vcxsk ;
	arrVS_Vessel_GSK00[4] = Vcysk  ;
	arrVS_Vessel_GSK00[5] = Vchsk ;

	///

			double valRez =  RZW__( VAlCalibro
                  , arrVS_Vessel_GSK00
									,arrVSTarg_KGSK0
									,&TRZW_Output
									,valFi0   // начальный угол места точки бросания
									,valQu0   // начальный азимут точки бросания
									// множитель приращения угла места для 2-го этапа алгоритма
									,valMulFi
									// множитель приращения азимута для 2-го этапа алгоритма
									,valMulQu
												);



			valPrevFi0= valFi0;
			valPrevQu0= valQu0;

			tpst1=TRZW_Output.tp;

			massq[i] = sin(TRZW_Output.q);
			mascq[i] = cos(TRZW_Output.q);
			masfi[i] = TRZW_Output.fi;
			mastp[i] = TRZW_Output.tp;
			mastr[i] = ti;
			masDy[i] = TRZW_Output.dy;
			masx[i] = cos(TRZW_Output.fi)*sin(TRZW_Output.q);
			masy[i] = cos(TRZW_Output.fi)*cos(TRZW_Output.q);
			mashh[i] = sin(TRZW_Output.fi);

		}	//конец цикла по i

		//Начальные значения координат и скоростей своего корабля
		shipX=shipX0; shipY=shipY0; shipH=shipH0;
		shipVX=shipVX0; shipVY=shipVY0;shipVH=shipVH0;

		double xrzv1,yrzv1,hrzv1;




		ti=t00-ts;
		double kof1,kof2,kof3,kof4;

//-------------------------------------------------------------------------
// ---------------------------------RZW1-----------------------------------
//-------------------------------------------------------------------------
		int ichetchik = 0 ;

		for(int ii0 =0; ii0 < (*pQuantShots) ; ii0++)          //  цикл по к-ву заданных снарядов
		{

			ti= ti+ts;

			double mmy,mmy1,mmy2;
			int mmyn;

			mmy1= ti/DT_rzv;
			mmyn= (int)(mmy1+dob);

			if(mmyn> SZnn1)
			{
			 break;		// ГАРАНТИЯ ПОПАДАНИЯ В ГРАНИЦЫ МАССИВОВ
			}

			mmy = mmy1-mmyn;
			mmy2= mmy*mmy;

			int NOM1,NOM2;
     /*
			if(mmyn==0) // парамтры берутся  по предыдущему значению (оно единственное, с нулевым номером)
			{
				sqrzv1 = massq[0];  // синус курс угла  наведения
				cqrzv1 = mascq[0];  // косинус курс угла  наведения
				firzv1 = masfi[0];  // уголт места  наведения
				tprzv1 = mastp[0];  // полетное время
				dyrzv1 = masDy[0];  // дальность точки встречи
				xrzv1 = masx[0];    // cos(TRZW_Output.fi)*sin(TRZW_Output.q)
				yrzv1 = masy[0];    // cos(TRZW_Output.fi)*cos(TRZW_Output.q)
				hrzv1 = mashh[0];   // sin(TRZW_Output.fi)
			}

			if(mmyn==1)
			{
				sqrzv1 = massq[0]+mmy*(2.0*(massq[1]-massq[0]));
				cqrzv1 = mascq[0]+mmy*(2.0*(mascq[1]-mascq[0]));
				firzv1 = masfi[0]+mmy*(2.0*(masfi[1]-masfi[0]));
				tprzv1 = mastp[0]+mmy*(2.0*(mastp[1]-mastp[0]));
				dyrzv1 = masDy[0]+mmy*(2.0*(masDy[1]-masDy[0]));
				xrzv1 = masx[0]+mmy*(2.0*(masx[1]-masx[0]));
				yrzv1 = masy[0]+mmy*(2.0*(masy[1]-masy[0]));
				hrzv1 = mashh[0]+mmy*(2.0*(mashh[1]-mashh[0]));
			}

			if(mmyn==2)
			{
				kof1 = 2.5*mmy+0.5*mmy2;
				kof2 = 2-4.0*mmy-mmy2;
				kof3 = 1.0-1.5*mmy-0.5*mmy2;
				sqrzv1 = massq[2]*kof1+massq[1]*kof2-massq[0]*kof3;
				cqrzv1 = mascq[2]*kof1+mascq[1]*kof2-mascq[0]*kof3;
				firzv1 = masfi[2]*kof1+masfi[1]*kof2-masfi[0]*kof3;
				tprzv1 = mastp[2]*kof1+mastp[1]*kof2-mastp[0]*kof3;
				dyrzv1 = masDy[2]*kof1+masDy[1]*kof2-masDy[0]*kof3;

				xrzv1 = masx[2]*kof1+masx[1]*kof2-masx[0]*kof3;
				yrzv1 = masy[2]*kof1+masy[1]*kof2-masy[0]*kof3;
				hrzv1 = mashh[2]*kof1+mashh[1]*kof2-mashh[0]*kof3;
			}

			if(mmyn>2)
			{
				kof1 = 2.5*mmy+0.5*mmy2;
				kof2 = 3.0-5.0*mmy-mmy2;
				kof3 = 3.0-3.5*mmy-0.5*mmy2;
				kof4 = 1.0-mmy;
				sqrzv1 = massq[mmyn]*kof1+massq[mmyn-1]*kof2-massq[mmyn-2]*kof3+massq[mmyn-3]*kof4;
				cqrzv1 = mascq[mmyn]*kof1+mascq[mmyn-1]*kof2-mascq[mmyn-2]*kof3+mascq[mmyn-3]*kof4;
				firzv1 = masfi[mmyn]*kof1+masfi[mmyn-1]*kof2-masfi[mmyn-2]*kof3+masfi[mmyn-3]*kof4;
				tprzv1 = mastp[mmyn]*kof1+mastp[mmyn-1]*kof2-mastp[mmyn-2]*kof3+mastp[mmyn-3]*kof4;
				dyrzv1 = masDy[mmyn]*kof1+masDy[mmyn-1]*kof2-masDy[mmyn-2]*kof3+masDy[mmyn-3]*kof4;

				xrzv1 = masx[mmyn]*kof1+masx[mmyn-1]*kof2-masx[mmyn-2]*kof3+masx[mmyn-3]*kof4;
				yrzv1 = masy[mmyn]*kof1+masy[mmyn-1]*kof2-masy[mmyn-2]*kof3+masy[mmyn-3]*kof4;
				hrzv1 = mashh[mmyn]*kof1+mashh[mmyn-1]*kof2-mashh[mmyn-2]*kof3+mashh[mmyn-3]*kof4;
			}
			*/

				sqrzv1 = fncLinExtrapolation(massq[mmyn-1], massq[mmyn],DT_rzv, mmy) ;
				cqrzv1 = fncLinExtrapolation(mascq[mmyn-1], mascq[mmyn],DT_rzv, mmy) ;//mascq[0]+mmy*(2.0*(mascq[1]-mascq[0]));
				firzv1 = fncLinExtrapolation(masfi[mmyn-1], masfi[mmyn],DT_rzv, mmy) ;//masfi[0]+mmy*(2.0*(masfi[1]-masfi[0]));
				tprzv1 = fncLinExtrapolation(mastp[mmyn-1], mastp[mmyn],DT_rzv, mmy) ;//mastp[0]+mmy*(2.0*(mastp[1]-mastp[0]));
				dyrzv1 = fncLinExtrapolation(masDy[mmyn-1], masDy[mmyn],DT_rzv, mmy) ;//masDy[0]+mmy*(2.0*(masDy[1]-masDy[0]));

				xrzv1 = fncLinExtrapolation(masx[mmyn-1], masx[mmyn],DT_rzv, mmy) ;//masx[0]+mmy*(2.0*(masx[1]-masx[0]));
				yrzv1 = fncLinExtrapolation(masy[mmyn-1], masy[mmyn],DT_rzv, mmy) ;//masy[0]+mmy*(2.0*(masy[1]-masy[0]));
				hrzv1 = fncLinExtrapolation(mashh[mmyn-1], mashh[mmyn],DT_rzv, mmy) ;//mashh[0]+mmy*(2.0*(mashh[1]-mashh[0]));


			//Экстраполяция полетного времени по линейке без склейки
			if(mmyn==0)
			{
				tprzv1 = mastp[0];
			}
			else
			{
				double dtmy = DT_rzv*mmy;
				double sktp = (mastp[mmyn]-mastp[mmyn-1])/DT_rzv;//dtmy; В.И. 30.08.2011
				tprzv1 = mastp[mmyn]+sktp*dtmy;
			}

			double dalg = sqrt(xrzv1*xrzv1+yrzv1*yrzv1);
			if (dalg < 0.0000001)
			{
				int iii =0.;
			}
			qs1= xrzv1/dalg;
			qc1= yrzv1/dalg;
			sf1= hrzv1/sqrt(dalg*dalg+hrzv1*hrzv1);
			cf1= sqrt(1.0-sf1*sf1);

			// Обращение к функции расчета координат и скоростей своего корабля
			raschet_coord_swoego_corablja (ti, VAlSigSins);

			//Вычисление координат артустановки в КГСК
			sqssk = sin(qssk);
			cqssk = cos(qssk);
			spsisk = sin(psisk);
			cpsisk = cos(psisk);
			stetask = sin(tetask);
			ctetask = cos(tetask);

			Xay_ksk = Xay*ctetask+Hay*stetask;
			Yay_ksk = -Xay*spsisk*stetask+Yay*cpsisk+Hay*spsisk*ctetask;
			Hayk = -Xay*cpsisk*stetask-Yay*spsisk+Hay*cpsisk*ctetask;
			Xayk = Xay_ksk*cqssk+Yay_ksk*sqssk;
			Yayk = -Xay_ksk*sqssk+Yay_ksk*cqssk;

			sqbc = sin(qbc);
			cqbc = cos(qbc);


			//Пересчет УВН и УГН с учетом ошибок рысканья и качек
			double sfn,cfn,sqn,cqn;

			//Синусы и косинусы точных курса своего и качек
			sqs  = sin(qbc);  cqs  = cos(qbc);
			spsi = sin(psi);  cpsi = cos(psi);
			steta= sin(teta); cteta= cos(teta);

			double sfn1,cfn1,sqn1,cqn1;
			//Второй способ (поворот координат )
			double Xay_k,Yay_k,Hay_k,Xay_kk,Yay_kk,Xay_p,Yay_p,Hay_p;

			// вектор направления АУ в КГСК
			Xay_k= cf1*qs1;
			Yay_k= cf1*qc1;
			Hay_k= sf1;

			//Преобразование вектора направления АУ в ПСК по измеренным качкам
			Xay_kk = Xay_k*cqssk-Yay_k*sqssk;
			Yay_kk = Xay_k*sqssk+Yay_k*cqssk;
			Xay_p = Xay_kk*ctetask-Yay_kk*stetask*spsisk-Hay_k*stetask*cpsisk;
			Yay_p = Yay_kk*cpsisk-Hay_k*spsisk;
			Hay_p = Xay_kk*stetask+Yay_kk*ctetask*spsisk+Hay_k*ctetask*cpsisk;

			//Обратные преобразования вектора направления АУ в КГСК по  точным качкам
			Xay_kk = Xay_p*cteta+Hay_p*steta;
			Yay_kk = -Xay_p*spsi*steta+Yay_p*cpsi+Hay_p*spsi*cteta;
			Hay_k = -Xay_p*steta*cpsi-Yay_p*spsi+Hay_p*cteta*cpsi;
			Xay_k = Xay_kk*cqs+Yay_kk*sqs;
			Yay_k = -Xay_kk*sqs+Yay_kk*cqs;

			//Вычисление координат артустановки в КГСК по точным качкам
			sqs = sin(qbc);
			cqs = cos(qbc);
			spsi = sin(psi);
			cpsi = cos(psi);
			steta = sin(teta);
			cteta = cos(teta);

			Xay_ksk = Xay*cteta+Hay*steta;
			Yay_ksk = -Xay*spsi*steta+Yay*cpsi+Hay*spsi*cteta;
			Hayk = -Xay*cpsi*steta-Yay*spsi+Hay*cpsi*cteta;
			Xayk = Xay_ksk*cqs+Yay_ksk*sqs;
			Yayk = -Xay_ksk*sqs+Yay_ksk*cqs;

			//Вычисление УВН и УГН с учетом ошибок рысканья и качек
			sfn1= Hay_k;
			cfn1= sqrt(1.0-sfn1*sfn1);

			sqn1=Xay_k/cfn1;
			if(sqn1<-1.0)	sqn1=-1.0; if(sqn1>1.0)	sqn1=1.0;

			cqn1= Yay_k/cfn1;
			if(cqn1<-1.0)	cqn1= -1.0; if(cqn1>1.0)	cqn1=  1.0;

			sfn= sfn1; cfn= cfn1;
			sqn= sqn1; cqn= cqn1;

			//Расчет ошибок рассеивания
			
		slm[10] = getGauss(0., 1. );
		slm[11] = getGauss(0., 1. );


			double sigmaq,sigmafi;
			sigmaq = Ssigmaq *0.001 * slm[10];
			sigmafi= Ssigmafi*0.001 * slm[11];

			//Пересчет УВН и УГН с учетом ошибок рассеивания снаряда
			double fin,qn;
			fin= asin(sfn); qn= acos(cqn); if(sqn<0.0) qn = M_2PI-qn;
			fin= fin+sigmafi; qn= qn+sigmaq;

			//Сведение баланса по дальности
			sf_rzw1= sf1_rzw1= sin(fin);
			cf_rzw1= cf1_rzw1= cos(fin);

			qs_rzw1= qs1_rzw1= sin(qn);
			qc_rzw1= qc1_rzw1= cos(qn);
			if(qs_rzw1< -1.0) qs_rzw1= -1.0; if(qs_rzw1> 1.0)	qs_rzw1= 1.0;
			if(qc_rzw1< -1.0) qc_rzw1= -1.0; if(qc_rzw1> 1.0)	qc_rzw1= 1.0;


			prisnak = 1;
			sigmawx= Ssigmawx*slm[6]*prisnak;
			sigmawz= Ssigmawz*slm[7]*prisnak;
			bky7= 1.0+sigmap*slm[8]*prisnak/100.0;
			vb= vb0+ sigmavb*slm[9]*vb0*prisnak/100.0;
			wy= Vchsk;


			//Расчет случайных чисел для учета ошибок табличного рассеивания

		slm[12] = getGauss(0., 1. );
		slm[13] = getGauss(0., 1. );
		slm[14] = getGauss(0., 1. );
		slm[15] = getGauss(0., 1. );


			vb= vb+sigmavbn*slm[12]*vb0/100.0;
			Ax= 1.0 + sigmaA1x*slm[13]/100.0;
			Ay= 1.0 + sigmaA1y*slm[14]/100.0;
			Az= 1.0 + sigmaA1z*slm[15]/100.0;

			tp= tprzv1;
			t_rzw1= ti;

			//При последнем сведении баланса по дальности используем данные
			// по чистой цели, а не с фильтра
			double Vyc,Vxc,Vhc,Voty,Votx,Voth;
			double Yc,Xc,Hc,r1,r2,r3,zzz,Prom1,Prom,Vot,cProm1;

			NOM1= (int)((t_rzw1+dob)/dti);
			NOM2= NOM1+1;
			if(NOM2> SZn) Vpc= 0.0;      // ГАРАНТИЯ ПОПАДАНИЯ В ГРАНИЦЫ МАССИВОВ
			else
			{


			 /*	Yzn=(Yz[NOM2]-Yz[NOM1])*(t_rzw1)/dti+(Yz[NOM1]*tm[NOM2]-Yz[NOM2]*tm[NOM1])/dti;
				Xzn=(Xz[NOM2]-Xz[NOM1])*(t_rzw1)/dti+(Xz[NOM1]*tm[NOM2]-Xz[NOM2]*tm[NOM1])/dti;
				Hzn=(Hz[NOM2]-Hz[NOM1])*(t_rzw1)/dti+(Hz[NOM1]*tm[NOM2]-Hz[NOM2]*tm[NOM1])/dti;

				Vyz=(VYz[NOM2]-VYz[NOM1])*t_rzw1/dti+(VYz[NOM1]*tm[NOM2]-VYz[NOM2]*tm[NOM1])/dti;
				Vxz=(VXz[NOM2]-VXz[NOM1])*t_rzw1/dti+(VXz[NOM1]*tm[NOM2]-VXz[NOM2]*tm[NOM1])/dti;
				Vhz=(VHz[NOM2]-VHz[NOM1])*t_rzw1/dti+(VHz[NOM1]*tm[NOM2]-VHz[NOM2]*tm[NOM1])/dti;
				*/
				double valDelT =  t_rzw1 - tm[NOM1];
				Yzn = fncLinExtrapolation(Yz[NOM1], Yz[NOM1 + 1],dti,valDelT) ;
				Xzn = fncLinExtrapolation(Xz[NOM1], Xz[NOM1 + 1],dti,valDelT) ;
				Hzn = fncLinExtrapolation(Hz[NOM1], Hz[NOM1 + 1],dti,valDelT) ;
				Vyz = fncLinExtrapolation(VYz[NOM1], VYz[NOM1 + 1],dti,valDelT) ;
				Vxz = fncLinExtrapolation(VXz[NOM1], VXz[NOM1 + 1],dti,valDelT) ;;
				Vhz = fncLinExtrapolation(VHz[NOM1], VHz[NOM1 + 1],dti,valDelT) ;;



				TRZW_Input.Xzel = Xzn;//-Xayk;
				TRZW_Input.Yzel = Yzn;//-Yayk;
				TRZW_Input.Hzel = Hzn;//-Hayk;
				TRZW_Input.Vxzel = Vxz ;
				TRZW_Input.Vyzel = Vyz ;
				TRZW_Input.Vhzel = Vhz ;


				if((t_rzw1/DTi+1)> SZn)
				{
				 break;    // ГАРАНТИЯ ПОПАДАНИЯ В ГРАНИЦЫ МАССИВОВ
				}



 double arrVS_KGSK_Targ0[6] ={0.}, arrVS_Shell_GSK[6] ={0.}, arrVS_Shell_SSK[9] ={0.}, arrVS_Targ_GSK[6] = {0.};

 double val_dtInt = 0.0005;
 arrVS_KGSK_Targ0[0] = Xzn;
 arrVS_KGSK_Targ0[1] = Yzn;
 arrVS_KGSK_Targ0[2] = Hzn;
 arrVS_KGSK_Targ0[3] = Vxz ;
 arrVS_KGSK_Targ0[4] = Vyz ;
 arrVS_KGSK_Targ0[5] = Vhz ;
	double valFlightTime = -1.;

	// формирование вектора состояния корабля в ГСК
	double arrVS_Vessel_GSK00[6] ={0.};
	arrVS_Vessel_GSK00[0] = shipX ;
	arrVS_Vessel_GSK00[1] = shipY ;
	arrVS_Vessel_GSK00[2] = shipH ;
	arrVS_Vessel_GSK00[3] = shipVX ;
	arrVS_Vessel_GSK00[4] = shipVY ;
	arrVS_Vessel_GSK00[5] = shipVH ;
	///


Prom = AimFun_RZW__(	 VAlCalibro
													, arrVS_Vessel_GSK00
													, arrVS_KGSK_Targ0
													,fin // угол места
													,qn // азимут
													,arrVS_Shell_GSK// вектор состояния снаоряда в точке встречи
													,arrVS_Shell_SSK// вектор состояния снаоряда в точке встречи
													, arrVS_Targ_GSK// вектор состояния снаоряда в точке встречи
													,&valFlightTime // подлётное время
													, val_dtInt);


	 double arrrDelt[6] ={0.};

	 MtrxMinusMatrx(arrVS_Targ_GSK, arrVS_Shell_GSK,6, 1, arrrDelt);
	 double valScal =ScalProduct(arrrDelt , &arrrDelt[3], 3);
	 double Promach = sqrt(Norm3( arrrDelt) * Norm3( arrrDelt) - valScal * valScal/ (Norm3(&arrrDelt[3]) * Norm3(&arrrDelt[3])));
	 int iiiij = 0;


			 double valAppointmentPoint = valFlightTime +  ti;
			 int numCur = valAppointmentPoint/dti;
			 double valDElT1 =  valAppointmentPoint - ((double)numCur) * dti;

			 double arrTrueShipGSK[3] ={0.}, arrTemp[3] ={0.};

			 arrTrueShipGSK[0] =  fncLinExtrapolation(parrTrueShipGSK_X [numCur], parrTrueShipGSK_X [numCur +1],dti,valDElT1) ;
			 arrTrueShipGSK[1] =  fncLinExtrapolation(parrTrueShipGSK_Y [numCur], parrTrueShipGSK_Y [numCur +1],dti,valDElT1) ;
			 arrTrueShipGSK[2] =  fncLinExtrapolation(parrTrueShipGSK_Z [numCur], parrTrueShipGSK_Z [numCur +1],dti,valDElT1) ;

			 MtrxMinusMatrx(arrTrueShipGSK, arrVS_Shell_GSK,3, 1, arrTemp);
			 const double VAlDistAppPoint = Norm3(arrTemp);

				Vpc = calcHittingProbabylity(Targ0, VAlCalibro, Prom);

				Vpc=Vpc*0.9; // с учётом в-сти срабатывания взрывателя
				// Умножение вероятности поражения цели одним снарядом на кривую Гладковского

				Wes = Targ0.mplnGlagkovsky.LinearValueApprox(VAlDistAppPoint);

				//-----------------------------------------===================
			//Вычисление вероятности поражения цели очередью снарядов
			Vpg= Vpg+(1-Vpo)*Vpc*Wes;
			Vpo= Vpo+(1-Vpo)*Vpc;

			masVpo[ii0]= Vpo;//Vpg;//Vpo;
//-----------------------------------------===================
					
			double *pCurBuff1 =   &parrBuff1 [ ii0 * quantBuffCols1];
			pCurBuff1[0] = (double) ii0;
			pCurBuff1[1] +=  VAlDistAppPoint;
			pCurBuff1[2] +=  Vpc;
			pCurBuff1[3] +=  Prom;

      

		}

	
		} // конец цикла по количеству заданных снарядов


		for(i=0;i<(*pQuantShots) ;i++) masppp[i] += masVpo[i];
		Ppp= Ppp+Vpo;
		Ppg= Ppg+Vpg;

		}//Конец цикла по реализациям

	//------Вероятность по дальностям-------
	for(i=0;i<(*pQuantShots) ;i++)  masppp[i]/= QuantIspit;
	//------Вероятность по дальностям-------

	//Вычисление вероятности поражения цели
	*pvalProb =Ppp/QuantIspit;
	*pvalProb_Gladk =Ppg/QuantIspit;

 // построение графиков для вероятностей

	for (int ii = 0; ii <(*pQuantShots) ; ii++)
		{
			for (int jj = 1; jj < quantBuffCols1; jj++)
			{
				 parrBuff1[ quantBuffCols1 * ii + jj] = parrBuff1[ quantBuffCols1 * ii + jj] / QuantIspit;
			}
		}

			// Построение графиков
 /*		double scalex =1., scaley =1000000.;
		if(wchOutFold[0] != 0)
		{
		double arrscaley[4] = {0.};
		 arrscaley[0] = 100.;
		 arrscaley[1] =1.;
		 arrscaley[2] = 100.;
		 arrscaley[3] = 100.;
		wchar_t wchFoldName[400] = {0};
		wcscpy(wchFoldName , wchOutFold) ;
		wcscat(wchFoldName, L"\\");

		for (int iq = 1; iq < quantBuffCols1 ; iq++)
	 {
		TYrWriteShapeFile::WriteOneReport(wchFoldName  // путь к папке
									, parrBuff1 // массив с информацией - матрица nBuffRows x nBuffCols
									,quantBuffCols1// - к-во переменных о корорых накоплена информация в буфере
									,(*pQuantShots)//  - к-во точек
									,wcharrFileNames1 //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30 // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,iq   // номер переменной по оси Y
									, arrscaley[0]  //  масштаб по оси X
									, arrscaley[iq]  // масштаб по оси Y
								   ) ;
	 }

	wchar_t wchAxesFileName1[300] ={0};
	wcscpy(  wchAxesFileName1,  wchFoldName);
	wcscat(wchAxesFileName1, L"AxesArr1.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName1,0., 100000.
	 ,0.,50000., 5.) ;
	}
    */
	free(parrBuff1);
	free(wcharrFileNames1);


	delete  Xz;
	delete  Yz;
	delete  Hz;
	delete  VXz;
	delete  VYz;
	delete  VHz;
	delete  tm;

	delete  Xf;
	delete  Yf;
	delete  Hf;
	delete  VXf;
	delete  VYf;
	delete  VHf;

	delete  tmf;
	delete  massq;
	delete  mascq;
	delete  masfi;
	delete  mastp;
	delete  masDy;
	delete  masx;
	delete  masy;
	delete  mashh;

	delete  	Xz1;
	delete  Yz1 ;
	delete  Hz1 ;
	delete  mastr;
	delete parrTrueShipGSK_X ;
	delete parrTrueShipGSK_Y ;
	delete parrTrueShipGSK_Z ;


	return(0);
}
//---------------------------------------------------------------------------------------------------
void raschet_coord_zeli (enumTargetType TargType,  const double WTargSkz, TTargBearing0  TargBear0, const double ti,  const double dti
	,double &targX, double &targY, double &targH
	,double &targVX, double &targVY, double &targVH)
{
	double vclamb =0.;
	switch( TargType )
	{
		case NOMOVING://неподвижная цель
			{
			targVX = targVY = targVH = 0.0;
			break;
			}
	  case ABOVEWATER://морская цель
		{
		  Qz     = TargBear0.mBearing + Qcm*cos(Tp*ti); // 2
		  targVX = TargBear0.mV *sin(Qz);
		  targVY = TargBear0.mV *cos(Qz);
          targVH = 0.0;
			break;
		 }

	  case PLANE://самолет
	    {

		// Вычисление курса самолета
		Qz += dttb2*(fii-Qz);

		 // Вычисление угла наклона траектории
	    lambi += dttb2*(temi-lambi);

		// Вычисление направления продольной оси цели
	   //	ix = ixm[1];
		 //	raschet_gauss_raspredelenie ();
	   //	slm[1] = slu; ixm[1]=ix;
		slm[1] = getGauss(0., 1. );

		fii   += dttb1*(fi0i-fii);//+ WTargSkz * slm[1];

	   //	ix = ixm[2];
		 //	raschet_gauss_raspredelenie ();
	   //	slm[2] = slu; ixm[2]=ix;
		slm[2] = getGauss(0., 1. );

			temi  += dttb1*(tem0i-temi);//+ WTargSkz * slm[2];

	    vclamb = TargBear0.mV  * cos(lambi);
	    targVX = vclamb*sin(Qz);
	    targVY = vclamb*cos(Qz);
	    targVH = TargBear0.mV  *sin(lambi);
						break;
            }
     }
	//Вычисление координат цели
	targX += dti*targVX;
	targY += dti*targVY;
	targH += dti*targVH;

	if (ti==0)
	{
		//targX = targX0;
	   //	targY = targY0;
		//targH = targH0;
		raschet_nach_coord_zeli (TargType,TargBear0,targX,targY, targH) ;
	}
}


//--------------------------------------------------------------------------------------------------
//******************************************************************************
 // Изменения № 6
 //******************************************************************************

// Подпрограмма вычисления правых частей уравнений баллистики
void fnc_Deriv_f(const double VAlCalibro, double *y,double *yp)
{
 // Описание локальных переменных, используемых в функции f
   int i,i1;
   double z2,z4,z5,z6,z7,z8,g,tan,z29,vt,vtt,va,va1,
	 cx1,cx2,cx,e,r3z2,tau1y,z4w,z5w,z6w,vakw,vaky,
	 spsi,stet,ctet,vbw,vnp,cyk,cpsi,kwsk,mn1,mn2,
	 fzk,fzm,qsbm,knmm,mmx0,mmx,va2,cocvt;
	{
   z2=y[2];
   z4=y[4];
   z5=y[5];
   z6=y[6];
   z7=y[7];
   z8=y[8];
   g=9.78049-3.086e-06*z2+0.0517*sl*sl;

   // Вычисление нормальной виртуальной температуры воздуха и ее
   //   градиента изменения в зависимости от высоты
	 if(z2<=9300.) {tan=t0-tkn*z2;
		  tau1y=-tkn;
		 }
   else if(z2<=12000.){z29=z2-9300.;
		       tan=230.-tkn*z29+1.172e-06*z29*z29;
					 tau1y=-tkn+0.2344e-05*z29;
		      }
	 else {tan=221.5;
	       tau1y=0.;
				}

////
	// Вычисление числа Маха (va)
   z4w=z4+wx; z5w=z5+wy; z6w=z6+wz;
   vt=sqrt(z4w*z4w+z5w*z5w+z6w*z6w);
   vtt=vt*sqrt(t0/tan);
   va=vtt/a;
   vakw=va*va; vaky=vakw*va;
	//



  //Расчет коэффициента сопротивления воздуха (cx) в зависимости от (va)
   if(va<=my3[1]) cx=alf0[1];
   else{
	if(va>=my3[7]) cx=alf0[8];
	else
	    {for(i=1;i<=6;i++)
	      {i1=i+1;
	       if((va>=my3[i])&&(va<my3[i1])) {
		  cx=alf0[i1]+alf1[i1]*va+alf2[i1]*vakw+alf3[i1]*vaky;
		  goto ff;                    }
	      }
	    }
	 } ff:;


 // Вычисление деривационной функции (fzm) в зависимости от (va)
   if(va<=ky3[1]) knmm=kn0[1];
   else {
	 if(va>=ky3[6]) knmm=kn0[7];
	 else {
	       for(i=1;i<=5;i++)
		{i1=i+1;
		 if((va>=ky3[i])&&(va<ky3[i1]))
		   {
		    knmm=kn0[i1]+kn1[i1]*va+kn2[i1]*vakw;
		    goto ff1;
		   }
		 }
	       }
	 } ff1:;
	 fzm=lm*lm*knmm/(VAlCalibro *hg0);

	 //	fzk=cbz*fzm; // Коэффициент деривации
		fzk=cbz1*fzm; // Коэффициент деривации
   // Расчет коэффициента аэродинамического аксиального демпфирующего
   //   момента (mmx) в зависимости  от числа Маха (va)
   if(va<=mash[1]) mmx0=mmash[1];
   else{
	if(va>=mash[k51-1]) mmx0=mmash[k51-1];
	else {
	      for(i=1;i<=k51-2;i++)
	       {va1=mash[i]; i1=i+1; va2=mash[i1];
		cx1=mmash[i];
		cx2=mmash[i1];
		if((va>=va1)&&(va<va2))
		  {
		   mmx0=cx1+(cx2-cx1)*(va-va1)/(va2-va1);
		   goto ff2;
		  }
		}
	      }
	} ff2:;


   mmx = mmx0 * VAlCalibro/(2.* vb);

   // Вычисление значений вспомогательных переменных, необходимых для
   //  расчета производных
   cocvt=coc*z7*vt;
   qsbm=cocvt*vt;
   vnp=sqrt(z4*z4+z6*z6);
   vbw=sqrt(z4*z4+z5*z5+z6*z6);
   spsi=z6/vnp; cpsi=sqrt(1.-spsi*spsi);
   stet=z5/vbw; ctet=vnp/vbw;
   cyk=-cbx*wx*stet/vt;
   r3z2=r3+z2;
   double delteta;
   delteta=-omz*cpsi+omx*spsi;
   kwsk=(-g*ctet-cyk*qsbm)/vbw+delteta+vnp/r3z2;
   mn2=mmx*qsbm*qb*lm/ibx;
   mn1=fzk*ibx*z8/(lm*qb);
   e=cocvt*cx;

   yp[1]=z4; //Производная горизонтальной дальности
	 yp[2]=z5;// Производная высоты
   yp[3]=z6;// Производная бокового отклонения

   // Производные составляющих скорости
   yp[4]=-cbx*e*z4w-z4*z5/r3z2-omy*z6+omz*z5+spsi*mn1*kwsk;
   yp[5]=-cby*e*z5+(z4*z4+z6*z6)/r3z2-g-omz*z4+omx*z6;
	 yp[6]=-cbz*e*z6w-z6*z5/r3z2-omx*z5+omy*z4-cpsi*mn1*kwsk;

   yp[7]=-z7*yp[2]/tan*(1./r+tau1y);// Производная относительной
				     // плотности  воздуха
	 yp[8]=-mn2*z8;// Производная аксиальной угловой скорости
 }
}


//-------------------------------------------------------------------------------------------------
//******************************************************************************
 // Конец изменения № 6
 //******************************************************************************

//-------------------------------------------------------------------------------------------------
 /*
 void raschet_randu_raspredelenie ()
//Функция расчета числа, равномерно распределенного на [0,2pi]
{
	  double yfl;
      ixf *= iy1;
	  yfl=ixf*iyf1;
      yfl -= floorl(yfl);
      fa=M_2PI*yfl;
}
 */
//--------------------------------------------------------------------------------------------------
void raschet_coord_swoego_corablja (double ti, const double VAlSigSins)
//Функция расчета координат и скоростей своего корабля
{
	double prh;

	//Вычисление текущего значения курса своего
	qbc  = Qco+Qcm*cos(Tp*ti-fa1);

	sqbc = sin(qbc);
	cqbc = cos(qbc);

	//Вычисление вектора скорости
	shipVX = Vc*sqbc;
	shipVY = Vc*cqbc;
	prh   = Tk*ti-fa2;
	shipVH = Hm*Tk*cos(prh);

	//Вычисление координат в ГСК

	shipX += dti*shipVX;

	shipY += dti*shipVY;

	shipH  = Hm*sin(prh);

	//Вычисление углов килевой и бортовой качек

	psi   = PSIm*cos(prh);
	teta  = TETAm*cos(Tb*ti-fa3);

//Определение измеренных рысканья, качек и высоты


		qssk = qbc+0.00041*sin(Tp*ti-fa1);
		psisk = psi + Osch1*sin(prh);
		tetask = teta+0.00041*sin(Tb*ti-fa3);
	//double temp = sqrt(Osch1 * Osch1 + VAlSigSins * VAlSigSins);
 //	qssk = qbc + temp* sin(Tp*ti-fa1);
 //	psisk = psi + temp * sin(prh);
 //	tetask = teta + temp *sin(Tb*ti-fa3);
//    Hcsk = 1.1*Hm*sin(prh);
	Hcsk = 0.1*sin(prh);
	shipH  = shipH+Hcsk;

	//Вычисление измеренных составляющих скорости
	Vcxsk = K11*Vc*sin(qssk);
	Vcysk = K11*Vc*cos(qssk);
//	Vchsk = K11*Tk*Hm*cos(prh);
	Vchsk = K01*Tk*cos(prh);

	if (ti==0)
	{
		shipX = shipX0;
		shipY = shipY0;
		shipH = shipH0;
	}
}



void raschet_coord_zeli_GSK_ideal (enumTargetType TargType, TTargBearing0  TargBear0, const double VAlT
	,double &valTargGSK_X, double &valTargGSK_Y, double &valTargGSK_H
	,double &valTargGSK_VX, double &valTargGSK_VY, double &valTargGSK_VH)
//Функция расчета координат и скоростей целей
{
	double vclamb;
	double valTargGSK_X0 = 0., valTargGSK_Y0 =0., valTargGSK_H0 =0.
	,valTargGSK_VX0 = 0., valTargGSK_VY0 =0., valTargGSK_VH0 =0.;
	raschet_nach_coord_zeli (TargType,TargBear0,valTargGSK_X0, valTargGSK_Y0, valTargGSK_H0);
	raschet_nach_velo_zeli (TargBear0,valTargGSK_VX0, valTargGSK_VY0, valTargGSK_VH0);
	if (TargType ==  NOMOVING)
	{
	  valTargGSK_VX = valTargGSK_VY = valTargGSK_VH = 0.0;
	}
	valTargGSK_VX = valTargGSK_VX0;
	valTargGSK_VY = valTargGSK_VY0;
	valTargGSK_VH = valTargGSK_VH0;

	valTargGSK_X =  valTargGSK_X0 +  VAlT * valTargGSK_VX0;
	valTargGSK_Y =  valTargGSK_Y0 +  VAlT * valTargGSK_VY0;
	valTargGSK_H =  valTargGSK_H0 +  VAlT * valTargGSK_VH0;

}
//-------------------------------------------------------------------------
// расчет начальных координат цели
void raschet_nach_coord_zeli (enumTargetType TargType,TTargBearing0  TargBear0
	,double &targX0, double &targY0, double &targH0)
{
	   const double VAlE0=asin(TargBear0.mH / TargBear0.mR);
  // Вычисление начальных координат цели
  targX0 = TargBear0.mR *cos(VAlE0)*sin(TargBear0.mBearing);
  targY0 = (TargBear0.mR) *cos(VAlE0)*cos(TargBear0.mBearing);
  if(TargType == ABOVEWATER) targH0=5.;
  else targH0 = TargBear0.mR *sin(VAlE0);


}

//----------------------------------------------------------------------------------------------------

void 	raschet_nach_velo_zeli (TTargBearing0 TargBear0, double &valTargGSK_VX0
	,  double &valTargGSK_VY0,  double &valTargGSK_VH0)
{
  valTargGSK_VX0 =  TargBear0.mV * cos(TargBear0.mTargHorizAng) * sin (TargBear0.mTargCourse);
  valTargGSK_VY0 =  TargBear0.mV * cos(TargBear0.mTargHorizAng) * cos (TargBear0.mTargCourse);
  valTargGSK_VH0 =  TargBear0.mV * sin(TargBear0.mTargHorizAng) ;
}
//Функция, вычисляющая наблюдаемые координаты цели
void raschet_coord_zeli_nabl(double *arrVSTargKGSK_True, double *arrVSTargKGSK_Zv )
{
	double prh,sezy,cezy,sqzy,cqzy,xzsk,yzsk,hzsk;
	double Xac_ksk,Yac_ksk,Xac_kgsk,Yac_kgsk,Hac_kgsk,
		   xkg_ac,ykg_ac,hkg_ac;
	double sqsk,cqsk,spsisk,cpsisk,stetask,ctetask,
		   sqs,cqs,spsi,cpsi,steta,cteta;

		sqsk = sin(qssk);
		cqsk = cos(qssk);
		spsisk = sin(psisk);
		cpsisk = cos(psisk);
		stetask = sin(tetask);
		ctetask = cos(tetask);

		sqs = sin(qbc);
        cqs = cos(qbc);
        spsi = sin(psi);
        cpsi = cos(psi);
        steta = sin(teta);
        cteta = cos(teta);



		slm[3] = getGauss(0., 1. );

		slm[4] = getGauss(0., 1. );

		slm[5] = getGauss(0., 1. );




	//Преобразование координат цели из КГСК в ПСК по точным качкам
		double xk,yk,xkgsk,ykgsk,hkgsk,
		       xap,yap,hap;
		xk = arrVSTargKGSK_True [0] *cqs - arrVSTargKGSK_True[1] *sqs;
		yk = arrVSTargKGSK_True [0] *sqs + arrVSTargKGSK_True[1] *cqs;

		xap = xk*cteta-yk*(spsi*steta) - arrVSTargKGSK_True[2] *(cpsi*steta);
		yap = yk*cpsi - arrVSTargKGSK_True[2] * spsi;
		hap = xk*steta+yk*(spsi*cteta) + arrVSTargKGSK_True[2] *(cpsi*cteta);


	//Преобразование координат цели из ПСК в КГСК по измеренным качкам
		xk = xap*ctetask+hap*stetask;
		yk = -xap*spsisk*stetask+yap*cpsisk+hap*spsisk*ctetask;
		hkgsk = -xap*cpsisk*stetask-yap*spsisk+hap*cpsisk*ctetask;
		xkgsk = xk*cqsk+yk*sqsk;
		ykgsk = -xk*sqsk+yk*cqsk;

		//Преобразование координат АС из ПСК в КГСК-ЦТ
		Xac_ksk = Xac*ctetask+Hac*stetask;
		Yac_ksk = -Xac*spsisk*stetask+Yac*cpsisk+Hac*spsisk*ctetask;
		Hac_kgsk = -Xac*cpsisk*stetask-Yac*spsisk+Hac*cpsisk*ctetask;
		Xac_kgsk = Xac_ksk*cqsk+Yac_ksk*sqsk;
		Yac_kgsk = -Xac_ksk*sqsk+Yac_ksk*cqsk;

		//Вычисление координат цели в КГСК-АС !!!!!!!!!!!!!!!!!
		xkg_ac = xkgsk - Xac_kgsk;
		ykg_ac = ykgsk - Yac_kgsk;
	    hkg_ac = hkgsk - Hac_kgsk;

		//Вычисление направляющих векторов целевой системы координат
		//hkg = hkg - (xkg*xkg+ykg*ykg)/R32;
		prh=sqrt(xkg_ac*xkg_ac+ykg_ac*ykg_ac+hkg_ac*hkg_ac);

		if (prh < 0.001)
		{
		 int iii = 0;
		 return;
		}
		sezy = hkg_ac/prh;
		cezy = sqrt(1.-sezy*sezy);
		sqzy = xkg_ac/(prh*cezy);
		cqzy = ykg_ac/(prh*cezy);

		//Вычисление наблюдаемых координат цели в ЦСК
		xzsk = prh*SIGMAq*slm[5];
		yzsk = prh+SIGMAd*slm[3];
		hzsk = prh*SIGMAe*slm[4];

		//Преобразование наблюдаемых координат цели из ЦСК в КГСК
		arrVSTargKGSK_Zv[0] = xzsk*cqzy + yzsk*(sqzy*cezy) - hzsk*(sqzy*sezy) + Xac_kgsk;
		arrVSTargKGSK_Zv[1] = -xzsk*sqzy + yzsk*(cqzy*cezy) - hzsk*(cqzy*sezy) + Yac_kgsk;
		arrVSTargKGSK_Zv[2] = yzsk*sezy + hzsk*cezy + Hac_kgsk;


}

//-------------------------------------------------------------------------------------------------
double fncLinExtrapolation(const double VAlF0, const double VAlF1
	,const double VAlStepX,const double VAlDelX)
{
	return VAlF0 + (VAlF1 -VAlF0)/ VAlStepX *  VAlDelX;
}

//--------------------------------------------------------------------------------------------------
//OUTPUT:
//t00 - момент первого выстрела
// tp00 - полетное время первого выстрела
// tk - момент последнего выстрела
// tpkk - полетное время последнего выстрела
//
//
//
//
//
//
void fncPreviousArrangments(TTarget Targ0
   , const double VAlCalibro, const double VAlDn,const double VAlDk
	, const int n, const double dti,
	double &t00, double &tp00,  double &tk, double &tpkk)
{
		TTarget Targ =  Targ0;

		double Dwstr0 =  Norm3(Targ0.mTraject.marrVectSostGSK_Begin);//   InitTargData0.mR;

		int it0 =0, itk = 0;
		double  val_ti= 0.0;
		double valFi0, valQu0, valPrevFi0, valPrevQu0;
		double valMulFi;
		double valMulQu;
		RZW_Input  TRZW_Input;
		RZW_Output TRZW_Output;



 for( int i=0; i < n;i++)
 {
		if( i == (n-1))
		{
      int iiii = 0;
    }
		if(i>0)
		{
			valFi0 = valPrevFi0 -0.03;
			valQu0 = valPrevQu0;
			valMulFi = 0.5;
			valMulQu = 1.0;
		}


		val_ti= dti * ((double) (i+1));

		//Обращение к функции расчета координат и скоростей своего корабля
		raschet_coord_swoego_corablja (val_ti, 0.);

		//Обращение к функции расчета координат и скоростей целей

		  Targ.recalcTrajPoint(val_ti);
		//Координаты цели в КГСК
		double x_kg  = Targ.mTraject.marrVectSostGSK[0] - shipX;
		double y_kg  = Targ.mTraject.marrVectSostGSK[1] -shipY;
		double h_kg  = Targ.mTraject.marrVectSostGSK[2]  -shipH;
		double vx_kg = Targ.mTraject.marrVectSostGSK[3] -shipVX;
		double vy_kg = Targ.mTraject.marrVectSostGSK[4] -shipVY;
		double vh_kg = Targ.mTraject.marrVectSostGSK[5]  -shipVH;


		if((i==1)||((i%10)==0))
		{

			//Вычисление координат артустановки в КГСК
			double sqs = sin(qbc);
			double cqs = cos(qbc);
			double spsi = sin(psi);
			double cpsi = cos(psi);
			double steta = sin(teta);
			double cteta = cos(teta);

			double Xay_ksk = Xay*cteta+Hay*steta;
			double Yay_ksk = -Xay*spsi*steta+Yay*cpsi+Hay*spsi*cteta;
			double Hayk = -Xay*cpsi*steta-Yay*spsi+Hay*cpsi*cteta;
			double Xayk = Xay_ksk*cqs+Yay_ksk*sqs;
			double Yayk = -Xay_ksk*sqs+Yay_ksk*cqs;

		 double arrVSTarg_KGSK0[6] ={0.};
		 arrVSTarg_KGSK0[0] = x_kg - Xayk;
		 arrVSTarg_KGSK0[1] = y_kg - Yayk;
		 arrVSTarg_KGSK0[2] = h_kg - Hayk;
		 arrVSTarg_KGSK0[3] = vx_kg;
		 arrVSTarg_KGSK0[4] = vy_kg;
		 arrVSTarg_KGSK0[5] = vh_kg;
			if(itk==1) return;

			if(i==0)
			{
				InitFi0Qu0(arrVSTarg_KGSK0,valFi0,valQu0, vb);  // БЫЛА ОШ(ИБКА!!!!! 28ю102017
				valFi0= 0.05;
				valMulFi= 0.5;
				valMulQu= 0.9;
			}

	double arrVSTarg0[6] = {0.}; // вектор состояния цели в в КГСК  на момент выстрела
			arrVSTarg0[0] = x_kg - Xayk;
			arrVSTarg0[1] = y_kg - Yayk;
			arrVSTarg0[2] = h_kg - Hayk;
			arrVSTarg0[3] = vx_kg;
			arrVSTarg0[4] = vy_kg;
			arrVSTarg0[5] = vh_kg;

	// формирование вектора состояния корабля в ГСК
	double arrVS_Vessel_GSK00[6] ={0.};
	arrVS_Vessel_GSK00[0] = shipX ;
	arrVS_Vessel_GSK00[1] = shipY ;
	arrVS_Vessel_GSK00[2] = shipH ;
	arrVS_Vessel_GSK00[3] = shipVX ;
	arrVS_Vessel_GSK00[4] = shipVY ;
	arrVS_Vessel_GSK00[5] = shipVH ;
	///

	 double valRez =  RZW__( VAlCalibro
                  ,arrVS_Vessel_GSK00
									,arrVSTarg0
									,&TRZW_Output
									,valFi0   // начальный угол места точки бросания
									,valQu0   // начальный азимут точки бросания
									// множитель приращения угла места для 2-го этапа алгоритма
									,valMulFi
									// множитель приращения азимута для 2-го этапа алгоритма
									,valMulQu
	                     	);


			valPrevFi0= valFi0;
			valPrevQu0= valQu0;


			double  Dwstr = TRZW_Output.dy;

			double Vd = (x_kg*vx_kg+y_kg*vy_kg+h_kg*vh_kg)/sqrt(x_kg*x_kg+y_kg*y_kg+h_kg*h_kg);
			if(it0==0)
			{
				if(Dwstr < VAlDn)
				{
				t00=val_ti;
				it0=1;
				tp00 = TRZW_Output.tp;

				}
			}

			if(itk==0)
			{
				if(Dwstr <= VAlDk)
				{
					tk = val_ti;
					itk = 1;
					tpkk = TRZW_Output.tp;
					return;
				}
				else
				{
					if(it0!=0)
					{
						if((Dwstr >= Dwstr0)||(Vd > -20.0)||(TRZW_Output.fi>(GRAD_TO_RAD*85)))
						{
						tk=val_ti;
						itk=1;
						tpkk = TRZW_Output.tp;
						return;
						}
					}
				}
			}
			Dwstr0=Dwstr;
		}

 }
 int iii=0;
 }

//--------------------------------------------------------------------------

// --------------------------------------------------------------------------------------------------

// --------------------------------------------------------------------------
double AimFun_RZW__(	const double VAlCalibro
										, double *arrVS_Vessel_GSK0
													, double *arrVS_KGSK_Targ00 // вектор состояния в КГСК
													,double Fi // угол места
													,double Qu // азимут
													,double *arrVS_Shell_KGSK// вектор состояния снаоряда в точке встречи
													,double *arrVS_Shell_SSK// вектор состояния снаоряда в точке встречи
													, double *arrVS_Targ_KGSK// вектор состояния цели  в точке встречи
													,double* tp_manul // подлётное время
													,double val_dtInt
												)
{
	int i;

	double vxc, vyc, vhc; // компоненты скорости цели в КГСК
	double xc_prev, yc_prev, hc_prev; // предыдущие координаты  цели в КГСК
	double xsn_prev, ysn_prev, hsn_prev; // предыдущие координаты снаряда в КГСК
	double T = 0.0; //тек.полетное время снаряда

	double y[10], yp[10], y_prev[10];
	double dR = 0.5; // расстояние между снарядом и целью
	double dR_2X = dR * dR; // расстояние в квадрате

	double F0prev = 10000000.0; // невязка- кв-т расст. между снарядом и целью на момент T
	double F0; // невязка- кв-т расст. между снарядом и целью на момент T+dt

	double cosQu, sinQu;

	// Вычисление коэффициентов согласования
	cbx = cbxm[k1 - 1];
	cby = cbym[k1 - 1];
	cbz1 = cbzm[k1 - 1];
	double fm1, fk1, fk2, fi2;

	fi2= Fi * 180.0 / M_PI;
	for(int m1=1;m1<k1;m1++)
	{ // fi2= asin(sin(Fi))*180.0/M_PI;
		if(fi2<=fim[m1])
		{ if(m1==1)
			{ cbx= cbxm[1]; cby= cbym[1]; cbz1= cbzm[1]; goto met;}
			else
			{ fm1= fim[m1 - 1];
				fk1= cbxm[m1 - 1];
				fk2= cbxm[m1];
				cbx= fk1 + (fk2 - fk1) * (fi2 - fm1) / (fim[m1] - fm1);
				fk1= cbym[m1 - 1];
				fk2= cbym[m1];
				cby= fk1 + (fk2 - fk1) * (fi2 - fm1) / (fim[m1] - fm1);
				fk1= cbzm[m1 - 1];
				fk2= cbzm[m1];
				cbz1 = fk1 + (fk2 - fk1) * (fi2 - fm1) / (fim[m1] - fm1);
				goto met;
			}
		}
	}
	met: ;
	cbz= cbx;
    // cbz1 = cbz;

	cosQu= cos(Qu);  sinQu= sin(Qu);

	omx=  2.0*om * cl * cosQu;
	omz= -2.0*om * cl * sinQu;
	omy=  2.0*om * sl;

	wx=  wet*(caw*calf+saw*salf)+sigmawx;//+Vc*(cqbc*calf+sqbc*salf);
	wz=  wet*(saw*calf-caw*salf)+sigmawz;//+Vc*(sqbc*calf-cqbc*salf);

	// пересчет скорости корабля своего в ССК
	double arrVessVelocity_SSK[3] ={0.};
	From_xyhSSK_To_xyhKGSK(cosQu, sinQu, &arrVS_Vessel_GSK0[3],arrVessVelocity_SSK);
	///

	// добавление к отнгосительной скорости цели( в КГСК) скорости корабля
	double arrVS_KGSK_Targ0[6] = {0.}, arrT0[3] ={0.};
	memcpy (arrVS_KGSK_Targ0, arrVS_KGSK_Targ00, 6 * sizeof(double));
	MtrxSumMatrx(&arrVS_KGSK_Targ0[3], arrVessVelocity_SSK,1, 3, arrT0) ;
	memcpy (&arrVS_KGSK_Targ0[3], arrT0, 3* sizeof(double));
	///

	// начальные зачения параметров СДУ  (системы диф.уравнений)
	y[1]= y[2] = y[3] = 0.0; // нач.точка траектории
	y[4]= vb * cos(Fi) + arrVessVelocity_SSK[0]; // нач.значение скорости снаряда вдоль оси X
	y[5]= vb * sin(Fi) + arrVessVelocity_SSK[1]; // нач.значение вертикальной скорости снаряда
	y[6]= arrVessVelocity_SSK[2]; // нач.значение боковой скорости снаряда
	y[7]= bky7; // плотность воздуха в нач. точке траектории
	y[8]= w0; // нач. скорость вращения снаряда

	F0prev = Norm3( arrVS_KGSK_Targ0);
	F0prev *= 2.;

	// while(F0>=dR*dR)
	bool brez = false;
	int j;

	double val_dtInt_2X = val_dtInt * val_dtInt;
	double a, b, c;
	double arrVS_Shell_KGSK_Cur[6] ={0.}, arrVS_Shell_KGSK_Prev[6] ={0.}, arrVS_TargExtr_Prev[6] ={0.};
  double arrPosTarg_KGSK_ExtrCur[6] ={0.}, arrTemp[6] ={0.}, arrTemp1[6] ={0.};
	memcpy(arrPosTarg_KGSK_ExtrCur, arrVS_KGSK_Targ0, 6 * sizeof(double));



	for(j=0;j<1000000000;j++)
	{
		for(i=1;i<=8;i++) y_prev[i]= y[i];

		fnc_Deriv_f(VAlCalibro, y, yp);
		T+= val_dtInt;

		// for(i=1;i<=8;i++) y[i]+= dt*yp[i];  // следуещее значение параметров СДУ
		for(i=1;i<=3;i++) y[i]+= val_dtInt*yp[i]+yp[i+3]*val_dtInt_2X/2.0;
		for(i=4;i<=8;i++) y[i]+= val_dtInt*yp[i];

		// Обратное преобразование координат снаряда из ССК в КГСК

		From_xyhSSK_To_xyhKGSK(cosQu, sinQu, &y[1],arrVS_Shell_KGSK_Cur);


		MatrxMultScalar(&arrVS_KGSK_Targ0[3], 3, 1, T, arrTemp);
		MtrxSumMatrx(arrTemp, arrVS_KGSK_Targ0,3, 1, arrPosTarg_KGSK_ExtrCur) ;


		MtrxMinusMatrx(arrPosTarg_KGSK_ExtrCur, arrVS_Shell_KGSK_Cur, 3, 1, arrTemp1);

		 F0 = Norm3( arrTemp1) ;
		if(F0>F0prev)
		{
			brez= true;
			break;

		}


	//	if(0.001> fabs(F0-F0prev))
	 //	{
	 //		brez= true;
	 //		break;

	// }
		 F0prev = F0;
		}

 From_xyhSSK_To_xyhKGSK(cosQu, sinQu, &y[4],&arrVS_Shell_KGSK_Cur[3]);
 memcpy(arrVS_Shell_KGSK, arrVS_Shell_KGSK_Cur, 6 * sizeof(double));
	*tp_manul= T;
	memcpy(arrVS_Shell_SSK, y, 9 * sizeof(double));
	memcpy(arrVS_Targ_KGSK, arrPosTarg_KGSK_ExtrCur, 6 * sizeof(double));
	return(F0);
}

// --------------------------------------------------------------------------------------------------
void From_xyhKGSK_To_xyhSSK(double cosQu, double sinQu, double *arrKGSKInp, double *arrOut)
{
	arrOut[1] = arrKGSKInp[2];
	arrOut[0] = sinQu * arrKGSKInp[0] + cosQu * arrKGSKInp[1];
	arrOut[2] = cosQu * arrKGSKInp[0] - sinQu * arrKGSKInp[1];
}

// ------------------------------------------------------------------------------
 void From_xyhSSK_To_xyhKGSK(double cosQu, double sinQu, double *arrSSKInp, double *arrOut)
{
	arrOut[2] = arrSSKInp[1];
	arrOut[0] = sinQu * arrSSKInp[0] + cosQu * arrSSKInp[2];
	arrOut[1] = cosQu * arrSSKInp[0] - sinQu * arrSSKInp[2];
}
 void From_xyhKGSK_To_xyhSSK(double cosQu, double sinQu, double x,
	double y, double h, double& x1, double& y1, double& h1)
{
	x1 = x * sinQu + y * cosQu;
	h1 = x * cosQu - y * sinQu;
	y1 = h;
}

// ------------------------------------------------------------------------------
 void From_xyhSSK_To_xyhKGSK(double cosQu, double sinQu, double x,
	double y, double h, double& x1, double& y1, double& h1)
{
	x1 = x * sinQu + h * cosQu;
	y1 = x * cosQu - h * sinQu;
	h1 = y;
}




// ------------------------------------------------------------------------------
double  RZW__(const double VAlCalibro,double *arrVS_Vessel_GSK00,  double *arrVSTarg0, RZW_Output *TRZW_Output,
	double& Fi0, // начальный угол места точки бросания
	double& Qu0, // начальный азимут точки бросания
	double kMuldFi_2, // множитель приращения угла места для 2-го этапа алгоритма
	double kMuldQu_2 // множитель приращения азимута для 2-го этапа алгоритма
	 //,	double Vb0 // нач.скорость снаряда
		)
{
	FILE* pf;

	int sch=0;

	int i;
	int step; // текущее число вызовов целевой ф-ии

	double dR = 0.5; // расстояние достаточное для поражения цели
	double dFi, dQu; // приращения по углам наведения (рад.)
	double d1, d2, dFi1;

	double Xc, Yc, Hc, Vcx, Vcy, Vch; // коор-ты снаряда и компоненты скорости цели в КГСК !!! arrTargVS_KGSK  [6]
	double Xsn, Ysn, Hsn; // координаты снаряда в КГСК
	double F0, F0prev; // текущее и предыдущее значение целевой ф-ии

	double arrVS_Shell_KGSK[6] = {0.}, arrVS_Targ_KGSK[6] = {0.}, arrVS_Shell_SSK[9] = {0.};
	//double arrTargExtr_KGSK[3] ={0.};
	double tp_manul; // подлётное время

 //	double x1, y1, h1, x2, y2, h2; // сделать векторы в ССК !!!

	pf = fopen("report RZW_MANUL.txt", "w+t");

	step = 0;
	F0 = 0;
	F0prev = 1000000000.0;
	double val_dtInt = 0.01;
	bool bShag = false;
	double cosQu, sinQu;
	double Fi0prev, Qu0prev;

	while(1) // (F0prev>=F0)
	{
		sch++;

		F0= AimFun_RZW__(VAlCalibro,arrVS_Vessel_GSK00, arrVSTarg0, Fi0, Qu0,
													arrVS_Shell_KGSK, arrVS_Shell_SSK, arrVS_Targ_KGSK
													, &tp_manul, val_dtInt );

		Fi0prev= Fi0;
		Qu0prev= Qu0;

		step++;

		if(sch>30) break;

		if ((F0< dR) && bShag)
		{
		 break;
		}

		if((F0>F0prev)||(fabs(F0-F0prev)<0.0001))
		{
			Fi0= Fi0prev;
			Qu0= Qu0prev;

													F0= AimFun_RZW__(VAlCalibro, arrVS_Vessel_GSK00, arrVSTarg0, Fi0, Qu0,
													arrVS_Shell_KGSK,arrVS_Shell_SSK,arrVS_Targ_KGSK
													, &tp_manul, val_dtInt );
			if(bShag)
			{
			 break;
			}
			else
			{
			bShag= true;
			val_dtInt= 0.0005;
			F0prev = F0;
			continue;
			}
		}

		cosQu= cos(Qu0); sinQu= sin(Qu0);

		double valTargEps = asin(arrVS_Targ_KGSK[2]/ Norm3( arrVS_Targ_KGSK));
		double valShellEps = asin(arrVS_Shell_KGSK[2] / Norm3( arrVS_Shell_KGSK));

		// dFi= acos((x1*x2+z1*z2)/(sqrt(x1*x1+z1*z1)*sqrt(x2*x2+z2*z2)));
		 dFi=	valTargEps -valShellEps;

		 Fi0= Fi0+dFi*kMuldFi_2;
		 dQu= atan2(arrVS_Targ_KGSK[0],arrVS_Targ_KGSK[1])-atan2(arrVS_Shell_KGSK[0],arrVS_Shell_KGSK[1]);
		 Qu0= Qu0 + dQu*kMuldQu_2;


		F0prev = F0;
	}

  if(Qu0<0.0) Qu0= 2.0*M_PI-Qu0;

	TRZW_Output->fi= Fi0;
	TRZW_Output->q = Qu0;
	TRZW_Output->xy= arrVS_Targ_KGSK[0];
	TRZW_Output->yy= arrVS_Targ_KGSK[1];
	TRZW_Output->hy= arrVS_Targ_KGSK[2];
	TRZW_Output->dy= Norm3( arrVS_Targ_KGSK);//sqrt(xy_manul*xy_manul+yy_manul*yy_manul+hy_manul*hy_manul);
	TRZW_Output->tp = tp_manul;
	memcpy(TRZW_Output->y, arrVS_Shell_SSK, 9 * sizeof(double));


	fclose(pf);
	return F0;
}

// ---------------------------------------------------------------------------
void InitFi0Qu0(double *arrVSTarg_KGSK0, double& Fi0, double& Qu0, double Vb0)
 {
	double Xc, Yc, Hc, Vcx, Vcy, Vch;
	double A, B, C, D;
	double S1, S2, T0;
	double Vsnx, Vsny, Vsnh, Vsn;

	// Вычисление начальных значений углов наведения Fi0, Qu0
	Xc =   arrVSTarg_KGSK0[0];
	Yc =   arrVSTarg_KGSK0[1];
	Hc =   arrVSTarg_KGSK0[2];
	Vcx =  arrVSTarg_KGSK0[3];
	Vcy =  arrVSTarg_KGSK0[4];
	Vch =  arrVSTarg_KGSK0[5];

	A = Xc * Xc + Yc * Yc + Hc * Hc;
	B = 2 * (Xc * Vcx + Yc * Vcy + Hc * Vch);
	C = Vcx * Vcx + Vcy * Vcy + Vch * Vch - Vb0 * Vb0;
	A /= 10000.0;
	B /= 10000.0;
	C /= 10000.0;

	D = sqrt(B * B - 4 * A * C);

	S1 = (-B - D) / (2.0 * A);
	S2 = (-B + D) / (2.0 * A);

	T0 = 1.0 / S2;

	Vsnx = Xc * S2 + Vcx;
	Vsny = Yc * S2 + Vcy;
	Vsnh = Hc * S2 + Vch;
	Vsn = sqrt(Vsnx * Vsnx + Vsny * Vsny + Vsnh * Vsnh);

	// Начальные значения углов бросания Fi0 Qu0 для решения задачи встречи
	Xc = Xc + Vcx * T0;
	Yc = Yc + Vcy * T0;
	Hc = Hc + Vch * T0;
	D = sqrt(Xc * Xc + Yc * Yc + Hc * Hc);
	Fi0 = asin(Hc / D);
	D = sqrt(Xc * Xc + Yc * Yc);
	if (Xc >= 0)
		Qu0 = acos(Yc / D);
	else
		Qu0 = 2.0 * M_PI - acos(Yc / D);
}


// -------------------------------------------------------------------------------------------------
double fncSignum(double x)
{
	return (x >0.)?1.:-1.;
}
//---------------------------------------------------------------------------
void ResetFun()
{

//	memset(masppp,0,sizeof(masppp));
 //	memset(obrmasppp,0,sizeof(obrmasppp));
 //	memset(masVpo,0,sizeof(masVpo));


	sigmaA1x = 1.0; // %
	sigmaA1y = 0.4; // %
	sigmaA1z = 0.0; // %
	sigmavbn = 0.4; //%

	Osch1= 0.0;
	Ax= Ay= Az= 0.0;
	dob = 0.0001;

	ti= dti=  dttb1=  dttb2= qbc= qbc_= sqbc= cqbc=
	shipVX= shipVY= shipVH= shipX= shipY= shipH= shipX0= shipY0= shipH0=
	psi= psi_= teta= teta_= qssk= psisk= psisk_AY= tetask=
	Hcsk=
	Vcxsk= Vcysk= Vchsk= 0.0;

   //	ixf1= ixf2= ixf3= 0;//Используются при вычислении фаз
	fa1,fa2,fa3= 0.0;//Фазы рысканья,бортовой и килевой качек
	delta= 0.0;
   //	ixf= 0;
	fa= 0.0;
	fii= fi0i= Qz= temi= tem0i= lambi= 0.0;

	//memset(ixm,0,sizeof(ixm));
 //	memset(slm,0,20 * sizeof(double));

	//xkg= ykg= hkg= 0.0;
	SIGMAd= SIGMAq= SIGMAe= 0.0;

  //	Xzy= Yzy= Hzy= 0.0;
 //	memset(y,0,sizeof(y));
//	memset(yp,0,sizeof(yp));
//	memset(y_rzw1,0,sizeof(y_rzw1));
 //	memset(yp_rzw1,0,sizeof(yp_rzw1));


	bl= om= g0= r3= t0= tkn= r= e1= dd0=	df0= dq0= h01= 0.0;
	ll= lm= lmz= lmg= ibx= hg0= cbz= w0= 0.0;

	vb= qb= 0.0;
//	memset(cbxm,0,sizeof(cbxm));
//	memset(cbym,0,sizeof(cbym));
 //	memset(cbzm,0,sizeof(cbzm));
 //	memset(fim,0,sizeof(fim));
 //	memset(aks,0,sizeof(aks));

	coc= cl= sl= omx= omy= omz=
	qs1= qc1= sf= cf= sf1= cf1= fi= q= qs= qc=
	qs1_rzw1= qc1_rzw1= sf_rzw1= cf_rzw1= sf1_rzw1=
	cf1_rzw1= fi_rzw1= q_rzw1= qs_rzw1= qc_rzw1=
	wet= aw= saw= caw= salf= calf=
	bky7= vbb= tp= dd1= 0.0;
	cbx= cby= Xayk_= Yayk_= Hayk_= Yzn= Xzn= Hzn= Xayk= Yayk= Hayk=
	Vyz= Vxz= Vhz= 0.0;

//	int nd;//Размерность матрицы коэффициентов аппроксимации для деривации //16.12.2010
//	int k1;//Число элементов в массивах коэффициентов согласования по координатам
//	int k51;//Используется (как число элементов в массиве) в расчете коэффициента аэродинамического аксиального демпфирующего

	wx= wy= wz= wx_= wy_= wz_= sigmawx= sigmawz= 0.0;

// double *Xz, *Yz, *Hz, *VXz, *VYz,*VHz, *tm;
// double *Xz1, *Yz1, *Hz1; //координаты цели без учета параваксов
// double  *VXzg, *VYzg,*VHzg;

	t_rzw1= t00_rzw1=0;
 //	memset(Kr_Glad,0,sizeof(Kr_Glad));
 //	memset(Mas_Dal,0,sizeof(Mas_Dal));

   //	unsigned long ix;//Нечетное число, содержащее не более 9 значащих цифр
   //	slu= 0.0;
	Vc= Qco= Qcm= PSIm= TETAm=	Tp= Tb=	Tk= Hm= K0=	K11=0.0;

 //	Vpo=Vpc=Ppp=Vpg=Ppg=Wes=0.0;


}


//---------------------------------------------------------------------------------------------------


// интгрирование уравнений до момента  падения
// OUTPUT:
// valDHoriz - горизонтальтеая дальн точки падения
void fncMove_130_Cal_Shell_TO_ZeroAlt_AND_ShowGraphs(double *arrStrSK_VS
   ,const double VAlCoefCx, const double VAlCoefCy, const double VAlCoefCz
   ,const double VAlVessV, const double VAlVessQ, const  double StepInt
		, wchar_t *wcharrPath,  double &valDHoriz )
{

  const int QUANT_COLS = 9 , QUANT_POINTS_MAX = 3000;
  const double DEL_T = 0.1;
	double *parrBuff;
	double *pscaleY;
	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// максимальная длина имени переменной
   if (wcharrPath)
   {

	parrBuff = new double [QUANT_COLS * QUANT_POINTS_MAX] ;
	memset (parrBuff, 0, QUANT_COLS * QUANT_POINTS_MAX * sizeof(double)) ;
	pwcharrFileNames = new wchar_t [ QUANT_COLS * lenName] ;
	memset (pwcharrFileNames, 0, QUANT_COLS * lenName* sizeof(wchar_t)) ;

	wcscpy( &pwcharrFileNames[ 0 * 30], L"t");
	wcscpy( &pwcharrFileNames[ 1 * 30], L"X");
	wcscpy( &pwcharrFileNames[ 2* 30],  L"Y");
	wcscpy( &pwcharrFileNames[ 3 * 30], L"Z");
	wcscpy( &pwcharrFileNames[ 4 * 30], L"Vx");
	wcscpy( &pwcharrFileNames[ 5 * 30], L"Vy");
	wcscpy( &pwcharrFileNames[ 6 * 30], L"Vz");
	wcscpy( &pwcharrFileNames[ 7 * 30], L"P");
	wcscpy( &pwcharrFileNames[ 8 * 30], L"Omega");

	pscaleY = new double  [QUANT_COLS] ;
	pscaleY[1] = 0.1;
	pscaleY[2] = 0.1;
	pscaleY[3] = 0.1;
	pscaleY[4] = 1000.;
	pscaleY[5] = 1.;
	pscaleY[6] = 1.;
	pscaleY[7] = 1.;
	pscaleY[8] = 1000.;

	double iNumPoints = 0;
	 }

  int iCirc = 1000. / StepInt ;
 //  double valTTemp = mTCur + (( double)iCirc) * StepInt ;
  int i = 0 ;
  int iNupPointsOut = 0;
  double valTOut = 0. ;

  for ( i = 0; i < iCirc; i++)
  {
	 // calcEilerStep( VAlCoefCx,  VAlCoefCy,  VAlCoefCz,VAlVessV,  VAlVessQ,StepInt, arrStrSK_VS);
	  calcEilerStep(arrStrSK_VS/*,const double VAlVessV, const double VAlVessQ*/,  StepInt);
	  valTOut +=  StepInt;
	  if (arrStrSK_VS [2] < 0.)
	  {
		 break ;
	  }

	  if ( iNupPointsOut < QUANT_POINTS_MAX )
	  {

	 // valTOut = (double)valTCur ;

	   if (wcharrPath)
	   {
	   double *p = &parrBuff[iNupPointsOut *  QUANT_COLS] ;
	  // p[0] =  (double)mTCur ;
	   p[0] =  (double)valTOut ;
	   p[1] = (double)arrStrSK_VS [1];
	   p[2] = (double)arrStrSK_VS [2];
	   p[3] = (double)arrStrSK_VS [3];
	   p[4] = (double)arrStrSK_VS [4];
	   p[5] = (double)arrStrSK_VS [5];
	   p[6] = (double)arrStrSK_VS [6];
	   p[7] = (double)arrStrSK_VS [7];
	   p[8] = (double)arrStrSK_VS [8];
		}

	  iNupPointsOut++;

	  }



  }
   valDHoriz = sqrt(arrStrSK_VS [1]*arrStrSK_VS [1]+ arrStrSK_VS [3]*arrStrSK_VS [3]) ;

	 if (wcharrPath)
	 {
	 for (int j = 1; j < QUANT_COLS -1; j++)
	 {


	TYrWriteShapeFile::WriteOneReport(wcharrPath  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS // - к-во переменных о корорых накоплена информация в буфере
								  ,iNupPointsOut //  - к-во точек
								  ,pwcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,j  // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,pscaleY[j]  // масштаб по оси Y
								   );
	 }

	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS // - к-во переменных о корорых накоплена информация в буфере
								  ,iNupPointsOut //  - к-во точек
								  ,pwcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,1  // номер переменной по оси X
								  ,2  // номер переменной по оси Y
								  ,1.//  масштаб по оси X
								  ,1.  // масштаб по оси Y
								   );
   TYrWriteShapeFile::WriteOneReport(wcharrPath  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS // - к-во переменных о корорых накоплена информация в буфере
								  ,iNupPointsOut //  - к-во точек
								  ,pwcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,1  // номер переменной по оси X
								  ,3  // номер переменной по оси Y
								  ,1.//  масштаб по оси X
								  ,1.  // масштаб по оси Y
								   );
  wchar_t wchFileName4[300] = {0} ;
  wcscpy(wchFileName4, wcharrPath );
  wcscat(wchFileName4, L"\\Axes.shp");
  TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;

  delete parrBuff ;
	delete pwcharrFileNames ;
	delete pscaleY ;

   }

  // fncCalcMtrxPartialDeriv( VAlVessV,  VAlVessQ,  StepInt ) ;
  // fncCalcVectPartialDeriv_CoeffForm( marrDelta[8] ) ;
  // fncCalcVectPartialDeriv_Mass(  VAlVessV,  VAlVessQ,  StepInt, marrDelta[9] * mShellBody.mMass ) ;

}

// шаг етода эйлера на врея   mTCur +valStepInt
 void calcEilerStep(double *arrStrSK_VS/*,const double VAlVessV, const double VAlVessQ*/, const  double valStepInt)
 {

	double arrF[9] ={0.} ;
   // шаг интегрирования фазового вектора
   // f(130., arrStrSK_VS, mAlfDir, VAlVessV,  VAlVessQ, arrF) ;
	 fnc_Deriv_f(0.13, arrStrSK_VS, arrF) ;

	double arrT0[9] = {0.},arrT1[9] = {0.};
   MatrxMultScalar(arrF, 9, 1, valStepInt,arrT0); // f * dt
   MtrxSumMatrx(arrT0, arrStrSK_VS,9, 1, arrT1) ;
   memcpy(arrStrSK_VS, arrT1, 9 * sizeof( double)) ;   ///

 // mTCur += valStepInt ;
 }

double calcHittingProbabylity(TTarget Targ , const double VAlCalibro, double valProm)
{
	 TURPolyLine plnUZP;
	if (fabs (VAlCalibro - 0.13)< 0.001)
	{
		if (Targ.menumTargetType == PLANE)
		{
		//Вычисление So по заданному закону поражения (по R и h)
		double eps,R,h,x,Fx,x1,x2,Fx1,So;
		eps = 0.001;
		R = (12.0+1.2);//(9.0 + 1.2);///sqrt(2.5);//(9.0+1.2)/sqrt(2.5);//3.0;//12.0 + 1.2;
		h = 0.26;//0.26;//0.33;//0.21;//0.21;//0.16;

		x = 18.0;
		Fx = (R*R * h) / (1.0 - exp(-R*R/x)) - x;
		for (int  i=1; (fabs(Fx) > eps); i++)
		{
		x1 = Fx + x;
		Fx1 = (R*R * h) / (1.0 - exp(-R*R/x1)) - x1;
		x2 = Fx1 + x1;
		x = (x * x2 - x1*x1) / (x2 - 2 * x1 +x);
		Fx = (R*R * h) / (1.0 - exp(-R*R/x)) - x;
		}
		So=sqrt(0.5*x);
		if(valProm >= R) return 0.;
		else return exp(-valProm  * valProm/(2*So*So));


		}

		if ((Targ.menumTargetType == GARPUN_V300) || (Targ.menumTargetType == GARPUN_V700))
		{
			plnUZP = TURPolyLine( X_yzp, Y_yzp,7);
			return plnUZP.LinearValueApprox(valProm);
		}

	

		return 0.;



	}
}

#pragma package(smart_init)
