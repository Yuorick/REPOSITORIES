#include "HelicConstants.h"
#include <math.h>

// матрица моментов инерции для ьлчной модели
extern  long double ARrINertMtrx[9]=
                             /* { 4049.33, -3744.85,     0.
                            ,-3744.85, 65864.19,     0.
                            ,    0.  ,       0., 66309.42};*/
        //как у вована
                                { 19000, 0.,     0.
                                ,0., 28000,     0.
                                ,    0.  ,       0., 46000};

extern  long double ARrINertMtrxPrived[9]=

        //как у вована
                                { 1.9000, 0.,     0.
                                ,0., 2.8000,     0.
                                ,    0.  ,       0., 4.6000};

// матрица моментов инерции упрощенная
extern  long double ARrINertMtrx_[9]= { 4049.33, -3744.85,     0.
                                        ,-3744.85, 65864.19,     0.
                                        ,    0.  ,       0., 66309.42};;
// ПАРАМЕТРЫ ПЛАНЕРА
extern  long double constArrHelicPlanersData_[8 * 13] = {
   //                                 S	     Ox 	Oy  	Oz  	Px 	 Py	     Pz	   Cx0	   Cy0	  AlfCrit Alx	   Aly	   Alz
/* Горизонт. проекция фюзеляжа*/  19.26,  -0.528, 0.000, 0.000,  0.000, 0.000,  0.000, 0.450, 0.500, 10.000,  0.000,  0.000,  0.000,
/* Передняя проекция фюзеляжа */   4.66,   2.850, 0.283, 0.000,  0.000, 0.000,  0.000, 0.640, 0.000, 10.000,  0.000,  0.000, -90.000,
/* Боковая проекция фюзеляжа  */  18.96,  -1.171, 0.000, 0.000,  0.000, 0.000,  0.000, 0.640, 0.000, 10.000, 90.000,  0.000,   0.000,
/*  Правый стабилизатор       */   1.22,  -6.090, 0.661, 0.568,  0.795, 0.000,  0.384, 0.640, 1.500, 10.000,  0.000,  0.000,   0.000,
/* Левый стабилизатор         */   1.22,  -6.090, 0.661,-0.568,  0.795, 0.000,  0.384, 0.640, 1.500, 10.000,  0.000,  0.000,   0.000,
/* Правое крыло               */   2.48,   0.120, 0.290, 1.508, -0.588, 0.000,  0.732, 0.640, 1.500, 15.000,  0.000,  0.000,   6.000,
/* Левое крыло                */   2.48,   0.120, 0.290, -1.508, 0.000, 0.000, -0.732, 0.640, 1.500, 15.000,  0.000,  0.000,   6.000,
/*Руль направления            */   1.76,  -7.685, 0.424, 0.000,  1.170, 0.000, -0.296, 0.640, 1.500, 30.000,-90.000, -10.840,  0.000
};

// полная модель
extern  long double constArrHelicPlanersData[8 * 13] = {
   //                                 S	     Ox 	Oy  	Oz  	Px 	 Py	     Pz	   Cx0	   Cy0	  AlfCrit Alx	   Aly	   Alz
/* Горизонт. проекция фюзеляжа*/  19.26,   0.  ,   0.000,  0.000,  0.000, 0.000,  0.000, 0.640,   1., 10.000,    0.000,  0.000,  0.000,
/* Передняя проекция фюзеляжа */   4.66,   /*3.2873*/ 0.,   /*0.083*/0.   ,  0.000,  0.000, 0.000,  0.000, /*0.640*/0.35, 0.000, 10.000,  0.000,  0.000, -90.000,
/* Боковая проекция фюзеляжа  */  18.96,  -0.7337, 0.000,  0.000,  0.000, 0.000,  0.000, 0.64,   0.000, 10.000, 90.000,  0.000,   0.000,
/*  Правый стабилизатор       */   1.22,  -5.6527, 0.461,  0.568,  0.795, 0.000,  0.384, 0.640, 1.5, 10.000,  0.000,  0.000,   0.000,
/* Левый стабилизатор         */   1.22,  -5.6527, 0.461, -0.568,  0.795, 0.000,  0.384, 0.640, 1.5, 10.000,  0.000,  0.000,   0.000,
/* Правое крыло               */   2.48,   1.5,    0.2   , 1.508, -0.588, 0.000,  0.732, 0.640, 1.5,   15.000,  0.000,  0.000,   6.000,
/* Левое крыло                */   2.48,   1.5,    0.2  , -1.508,  0.000, 0.000, -0.732, 0.640, 1.5 ,  15.000,  0.000,  0.000,   6.000,
/*Руль направления            */   1.76,  -8.122,  0.224,  0.000,  1.170, 0.000, -0.296, 0.640, 1.500, 30.000,-90.000, -10.840,  0.000
};


extern const int NUmRowsArrHelicPlanersData = 8;
extern const int NUmColsArrHelicPlanersData = 13;
///
/// \brief The HelicConstants class
///
///
// ПАРАМЕТРЫ лопасти винта
extern const long double constBladeR = 7.2; // радиус ометаемой винтом площади
extern const long double constRadHorizHsarnir = 0.283; //расстояние от центра гориз шарнира до оси вращения вала винта marrDblSpinBoxBlade[0]
extern const long double constPofile_d0 = 0.142; // высота профиля у оси вала винта marrDblSpinBoxBlade[1]
extern const long double constPofile_d1 = 0.067; // высота профиля на конце marrDblSpinBoxBlade[2]
extern const long double constBlade_b = 0.46; // хорда     marrDblSpinBoxBlade[3]
extern const long double constBladeM = 26.; // масса лопасти
extern const long double constBladeCX0 = 1.;
///


// параметры  винтов  !!!!!!!!!!!!!!!!!!!!!
extern const int constQuantBlades = 3;
// максимально допустимый общий шаг НВ
extern const long double constFiMax = 15./ 180. * M_PI;

// угловая скорость вращениея винта
extern const long double constRotorOmega = 4.15 ; // уцгдловая скорость вращения винта в об/с
// базовачя система координат привязанная к оси вращения вала  винта
// ось Z направлена по оси вращения вала

// координаты точки основания оси вала винта в СвСК
extern const long double constShaftAxeX_ = 0.; // ЭТО ДЛЯ ОТЛАДКИ!!!!
extern const long double constShaftAxeY_ = 1.186;

//extern const long double constShaftAxeX =  -0.37 +0.4373 ;
//extern const long double constShaftAxeY = 1.186-0.2;//= 1.35;
//Это по измерениям у Володко для К26 максимум 0,28 м !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!.

///
// угол заклинения, рад., это угол между осью OZ СвСК вертолета и осью НВ
extern const long double constZaclinAng_ = 0.; // ЭТО упрощенная модель ДЛЯ ОТЛАДКИ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111111111111111111111111
//extern const long double constZaclinAng =  0.068215;//0.055;//!!0.05905; //  // ЭТО полная модель
// расстояние от основания вала винта до центра вращения верхнего винта
extern const long double constShaftUpperL_ = 2.347;
//extern const long double constShaftUpperL = 1.822;

// расстояние от основания вала винта до центра вращения нижнего винта
extern const long double constShaftLowerL_ = 0.646;
//extern const long double constShaftLowerL = 0.563;
///
/// \brief The HelicConstants class
///

// константы для инициализации вертолета по-другому, ось вала винта проходит через ЦМ и
// задается смщение по оси OX
   // угол заклинения
extern const  long double constAlfaZaklNew = 3.3 * M_PI/180.;

  // рассточние от ЦМ до втулки нижненего винта
extern   const  long double constXdistLow = 1.549;

  // рассточние от между втулками верхнего и нижнего винтов
extern   const  long double constDeltaXUp_Low = 1.3;
 ///




// масса вертолета
extern const long double constHelicMass = 8500.;

// массив наборов передаточных чисел
extern  long double constArrGearSets[] = {
    // 0.
    // набор передаточных чисел для Y=0 V=0
    // файл D:\REPOSITORIES\aircraft-model\PROGRAMS_C++\STABILITY\PeredChislaBolshZadachi\Chisla_Y=0_V=0_c=7_a=7_L=-9,9.csv
    // N0 500
 0.000000,0.000000,	0.000000,	-0.005100,	0.000000,	0.000000,	0.000000,	0.000000,	0.456000,	0.000000,	0.261000,	0.000000
,0.000000,-0.015100,	0.000000,	0.000000,	-0.036000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000
,0.000000,0.000000,	-0.009560,	0.000000,	0.000000,	-0.050986,	-0.030000,	0.000000,	0.000000,	-0.500000,	0.000000,	0.000000
,0.000000,0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	-0.040000,	0.000000,	0.000000,	0.000000,	-0.040000

    // 1.
    // набор передаточных чисел для Y=20 V=10
    // файл
   // D:\REPOSITORIES\aircraft-model\PROGRAMS_C++\STABILITY\PeredChislaBolshZadachi\Bol_Y=20_V=10\Bol_Vtor_Y=20_V=10_a=2,9_C=4,35_L=0,121_OTOBRANNIE.csv
   // набор чисел №1
  ,0.000000 ,0.000000	,0.000000	,-0.008833	,0.000000	,0.000000	,0.000000	,0.000000	,0.238000	,0.000000	,0.406500	,0.000000
  ,0.000000 , -0.005520	,0.000000	,0.000000	,-0.020280	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
  ,0.000000 ,0.000000	,-0.007011	,0.000000	,0.000000	,-0.050986	,-0.100000	,0.000000	,0.000000	,-0.350000	,0.000000	,0.000000
  ,0.000000  ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.050000	,0.000000	,0.000000	,0.000000	,-0.051000

    // 2.
    // подъем с (20,10) на (300,40)
   // файл
   // D:\REPOSITORIES\aircraft-model\PROGRAMS_C++\STABILITY\PeredChislaBolshZadachi\Bol_Y=300_V=40\Bol_Perv_c=2,69_a=4,05_L=0,1026_Otobrannie.csv
   // набор №2
    ,0.000000,0.000000	,0.000000	,-0.008833	,0.000000	,0.000000	,0.000000	,0.000000	,0.238000	,0.000000	,0.406500	,0.000000
    ,0.000000,-0.003680	,0.000000	,0.000000	,-0.016900	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000,0.000000	,-0.009518	,0.000000	,0.000000	,-0.050765	,-0.125000	,0.000000	,0.000000	,-0.300000	,0.000000	,0.000000
    ,0.000000,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.050000	,0.000000	,0.000000	,0.000000	,-0.051000





    // 3.
    //поворот на 90 град (300,40)
   // файл
   // D:\REPOSITORIES\aircraft-model\PROGRAMS_C++\STABILITY\PeredChislaBolshZadachi\Bol_Y=300_V=40\Bol_Perv_c=2,69_a=4,05_L=0,1026_Otobrannie.csv
   // набор №20
    //z0 = 400 Tперех = 27с


   /* ,0.000000,0.000000	,0.000000	,-0.008833	,0.000000	,0.000000	,0.000000	,0.000000	,0.238000	,0.000000	,0.406500	,0.000000
    ,0.000000,-0.003680	,0.000000	,0.000000	,-0.016900	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000,0.000000	,-0.038074	,0.000000	,0.000000	,-0.101531	,-0.150000	,0.000000	,0.000000	,-0.500000	,0.000000	,0.000000
    ,0.000000,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.050000	,0.000000	,0.000000	,0.000000	,-0.051000*/


    /* ,0.000000      , 0.000000	,0.000000	,-0.008833	,0.000000	,0.000000	,0.000000	,0.000000	,0.238000	,0.000000	,0.406500	,0.000000
     ,0.000000       ,-0.003680	,0.000000	,0.000000	,-0.016900	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
     ,0.000000       ,0.000000	,-0.006980	,0.000000	,0.000000	,-0.050765	,-0.125000	,0.000000	,0.000000	,-0.300000	,0.000000	,0.000000
     ,0.000000       ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.050000	,0.000000	,0.000000	,0.000000	,-0.051000*/
   /* ,0.000000	,0.000000	,0.000000	,-0.023000	,0.000000	,0.000000	,0.000000	,0.000000	,0.750000	,0.000000	,0.750000	,0.000000
    ,0.000000	,-0.020500	,0.000000	,0.000000	,-0.035500	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000	,0.000000	,-0.121000	,0.000000	,0.000000	,-0.241000	,-0.200000	,0.000000	,0.000000	,-1.210000	,0.000000	,0.000000
    ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.080100	,0.000000	,0.000000	,0.000000	,-0.016000*/


//D:\REPOSITORIES\aircraft-model\ФАЙЛЫ_С_ПЕРЕД_ЧИСЛАМИ\Line_Y=300_V=40\Bol_v2-1  номер 1
    ,0.000000	,0.000000	,0.000000	,-0.023000	,0.000000	,0.000000	,0.000000	,0.000000	,0.750000	,0.000000	,0.750000	,0.000000
    ,0.000000	,-0.020500	,0.000000	,0.000000	,-0.035500	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000	,0.000000	,-0.001000	,0.000000	,0.000000	,-0.011000	,-0.050000	,0.000000	,0.000000	,-0.070000	,0.000000	,0.000000
    ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.020100	,0.000000	,0.000000	,0.000000	,-0.001000






    // 4.
   // поворот на 90 град (300,50)
   // файл
   // D:\REPOSITORIES\aircraft-model\PROGRAMS_C++\STABILITY\PeredChislaBolshZadachi\Bol_Y=300_V=50\Bol_Perv_c=3_a=4,5_L=0,106_Otobrannie.csv
   //набор №80 (нумерация с 1)
    ,0.000000,0.000000	,0.000000	,-0.008833	,0.000000	,0.000000	,0.000000	,0.000000	,0.238000	,0.000000	,0.406500	,0.000000
    ,0.000000,-0.003680	,0.000000	,0.000000	,-0.016900	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000,0.000000	,-0.009518	,0.000000	,0.000000	,-0.050765	,-0.175000	,0.000000	,0.000000	,-0.700000	,0.000000	,0.000000
    ,0.000000,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.050000	,0.000000	,0.000000	,0.000000	,-0.051000





    // 5.
   // снижение с высоты 300 со скоростью 40 на высоту 100 со скоростью 40
   // файл
   // D:\REPOSITORIES\aircraft-model\PROGRAMS_C++\STABILITY\PeredChislaBolshZadachi\Bol_Y=100_V=40\2.csv
   //набор №1
    ,0.000000,0.000000,	0.000000,	-0.008833,	0.000000,	0.000000,	0.000000,	0.000000,	0.238000,	0.000000,	0.406500,	0.000000
    ,0.000000,-0.003680,	0.000000,	0.000000	,-0.016900	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000,0.000000	,-0.009518	,0.000000	,0.000000	,-0.050765	,-0.125000	,0.000000	,0.000000	,-0.300000	,0.000000	,0.000000
    ,0.000000,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.050000	,0.000000	,0.000000	,0.000000	,-0.051000

    // 6.
   // снижение с высоты 100 со скоростью 40 на высоту 20 со скоростью 10
   // файл
   // D:\REPOSITORIES\aircraft-model\PROGRAMS_C++\STABILITY\PeredChislaBolshZadachi\Bol_Y=20_V=10\Bol_Vtor_Y=20_V=10_a=2,9_C=4,35_L=0,121_OTOBRANNIE.csv
   //набор №1
    ,0.000000,0.000000	,0.000000	,-0.008833	,0.000000	,0.000000	,0.000000	,0.000000	,0.238000	,0.000000	,0.406500	,0.000000
    ,0.000000,-0.005520	,0.000000	,0.000000	,-0.020280	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000,0.000000	,-0.007011	,0.000000	,0.000000	,-0.050986	,-0.100000	,0.000000	,0.000000	,-0.350000	,0.000000	,0.000000
    ,0.000000,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.050000	,0.000000	,0.000000	,0.000000	,-0.051000

   // 7.
    // вираж на высоте 300 м со скоростью 40 м/с R =675
    ,-0.010500	,0.000000	,0.000000	,-0.043000	,0.000000	,0.000000	,0.000000	,0.000000	,0.550000	,0.000000	,0.750000	,0.000000
    ,0.000000	,-0.010500	,0.000000	,0.000000	,-0.020500	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000	,0.000000	,-0.061000	,0.000000	,0.000000	,-0.121000	,-0.150000	,0.000000	,0.000000	,-0.660000	,0.000000	,0.000000
    ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.090100	,0.000000	,0.000000	,0.000000	,-0.026000

    // 8. вираж на высоте 300 м со скоростью 50 м/с R = 820
    // D:\REPOSITORIES\aircraft-model\ФАЙЛЫ_С_ПЕРЕД_ЧИСЛАМИ\Turn_Y=300_V=50_R=820\Bol_v2 nomer 43
   /* ,-0.003000	,0.000000	,0.000000	,-0.023000	,0.000000	,0.000000	,0.000000	,0.000000	,0.650000	,0.000000	,0.650000	,0.000000
    ,0.000000	,-0.005500	,0.000000	,0.000000	,-0.020500	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000	,0.000000	,-0.061000	,0.000000	,0.000000	,-0.121000	,-0.175000	,0.000000	,0.000000	,-0.810000	,0.000000	,0.000000
    ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.090100	,0.000000	,0.000000	,0.000000	,-0.046000
    */
    ,-0.003000	,0.000000	,0.000000	,-0.023000	,0.000000	,0.000000	,0.000000	,0.000000	,0.650000	,0.000000	,0.650000	,0.000000
    ,0.000000	,-0.005500	,0.000000	,0.000000	,-0.020500	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000	,0.000000	,-0.061000	,0.000000	,0.000000	,-0.121000	,-0.175000	,0.000000	,0.000000	,-0.810000	,0.000000	,0.000000
    ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.110100	,0.000000	,0.000000	,0.000000	,-0.091000

    // 9. Висение на высотах 0 - 300 м
    ,0.000000	,0.000000	,0.000000	,-0.005500	,0.000000	,0.000000	,0.000000	,0.000000	,0.450000	,0.000000	,0.250000	,0.000000
    ,0.000000	,-0.000500	,0.000000	,0.000000	,-0.005500	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000	,0.000000	,-0.001000	,0.000000	,0.000000	,-0.011000	,-0.087500	,0.000000	,0.000000	,-0.170000	,0.000000	,0.000000
    ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.100100	,0.000000	,0.000000	,0.000000	,-0.011000

    // 10. разворот на висении
    // D:\REPOSITORIES\aircraft-model\ФАЙЛЫ_С_ПЕРЕД_ЧИСЛАМИ\Hover_Y=300\Bol_v1 №4
    ,0.000000	,0.000000	,0.000000	,-0.023000	,0.000000	,0.000000	,0.000000	,0.000000	,0.750000	,0.000000	,0.750000	,0.000000
    ,0.000000	,-0.010500	,0.000000	,0.000000	,-0.025500	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000	,0.000000	,-0.011000	,0.000000	,0.000000	,-0.021000	,-0.080000	,0.000000	,0.000000	,-0.190000	,0.000000	,0.000000
    ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.110100	,0.000000	,0.000000	,0.000000	,-0.041000


    // 11. вращение на высоте 300 м
    //
    ,-0.050500	,0.000000	,0.000000	,-0.093000	,0.000000	,0.000000	,0.000000	,0.000000	,0.450000	,0.000000	,1.050000	,0.000000
    ,0.000000	,-0.020500	,0.000000	,0.000000	,-0.060500	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000	,0.000000	,-0.661000	,0.000000	,0.000000	,-0.601000	,-0.550000	,0.000000	,0.000000	,-1.910000	,0.000000	,0.000000
    ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.120100	,0.000000	,0.000000	,0.000000	,-0.301000

// 12ю висение по новому (12 параметров под контолролем) на высоте 300 м
 // D:\REPOSITORIES\aircraft-model\ФАЙЛЫ_С_ПЕРЕД_ЧИСЛАМИ\Hover_Y=300\Bol_v1.csv №9
  //  ,-0.015500	,0.000000	,0.000000	,-0.088000	,0.000000	,0.000000	,0.000000	,0.000000	,0.750000	,0.000000	,0.950000	,0.000000
  //  ,0.000000	,-0.020500	,0.000000	,0.000000	,-0.040500	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
   // ,0.000000	,0.000000	,-0.001000	,0.000000	,0.000000	,-0.007000	,-0.285000	,0.000000	,0.000000	,-0.660000	,0.000000	,0.000000
   // ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.031000	,0.000000	,0.000000	,0.000000	,-0.036000


    // D:\REPOSITORIES\aircraft-model\ФАЙЛЫ_С_ПЕРЕД_ЧИСЛАМИ\Hover_Y=300\Bol_v1.csv № 190
    ,-0.015500	,0.000000	,0.000000	,-0.088000	,0.000000	,0.000000	,0.000000	,0.000000	,0.750000	,0.000000	,0.950000	,0.000000
    ,0.000000	,-0.020500	,0.000000	,0.000000	,-0.040500	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000	,0.000000	,-0.001000	,0.000000	,0.000000	,-0.007000	,-0.285000	,0.000000	,0.000000	,-0.660000	,0.000000	,0.000000
    ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.086000	,0.000000	,0.000000	,0.000000	,-0.036000


    // 13. висение по новому (12 параметров под контолролем) на высоте 0 м
     // D:\REPOSITORIES\aircraft-model\ФАЙЛЫ_С_ПЕРЕД_ЧИСЛАМИ\Hover_Y=0\Bol_v2.csv №10
    ,-0.010500	,0.000000	,0.000000	,-0.048000	,0.000000	,0.000000	,0.000000	,0.000000	,0.650000	,0.000000	,0.850000	,0.000000
    ,0.000000	,-0.030500	,0.000000	,0.000000	,-0.060500	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
    ,0.000000	,0.000000	,-0.037000	,0.000000	,0.000000	,-0.067000	,-0.135000	,0.000000	,0.000000	,-0.460000	,0.000000	,0.000000
    ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.086000	,0.000000	,0.000000	,0.000000	,-0.026000

    // 14. висение по новому (12 параметров под контолролем) на высоте 100 м
     // file:///D:/REPOSITORIES/aircraft-model/ФАЙЛЫ_С_ПЕРЕД_ЧИСЛАМИ/Hover_Y=100/Bol_v2.csv\Bol_v2.csv №7
    ,-0.010500	,0.000000	,.000000	,-0.048000	,0.000000	,0.000000	,0.000000	,0.000000	,0.650000	,0.000000	,0.850000	,0.000000
    ,0.000000	,-0.050500	,0.000000	,0.000000	,-0.070500	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
   ,0.000000	,0.000000	,-0.025000	,0.000000	,0.000000	,-0.055000	,-0.135000	,0.000000	,0.000000	,-0.435000	,0.000000	,0.000000
    ,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,-0.086000	,0.000000	,0.000000	,0.000000	,-0.026000

};

// параметры виража для упрощенной аэродинамики и винта
extern const long double consrArrTurnParams_SimpleHalic_SimpleRotor[] = {
    // Y    V       R    UgolSkolg
    300.,    40.,   675.,  0.
   ,300.,    50.,   820.,  0.
};

//HelicConstants::HelicConstants()
//{

//}
