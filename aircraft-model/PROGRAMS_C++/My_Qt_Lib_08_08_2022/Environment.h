//---------------------------------------------------------------------------

#ifndef EnvironmentH
#define EnvironmentH
#define ATM_PN0  1000.  //нормальное наземное давление воздуха
#define ATM_RoN0  1.2054  // нормальная назеная плотность воздуха
#define ATM_AN0  340.8   //нормальная назеная скорость звука в воздухе
#define ATM_TAYN0  289   //нормальная виртуальная температура на 0 высоте
#define ATM_R_UNIVER  287.05287   //удельная газовая постоянная воздуха
#define R_ZEMLI 6371000.
#define G_ZEMLI 9.80665
//---------------------------------------------------------------------------
// 9 балов волнения соответствует 20 м/с скорости ветра

class TEnvironment
{
public:
// Волнение на море по шкале Significance Wave Height (SWH). от 0 до 9 баллов
     int  mBallWave;
//  Сила ветра м/с
	 double mWind_V  ;
         // направление откуда дует ветер, отсчитывается от направления на Север по часовой стрелке
	 double mWind_Alf ;
     // скорость вертикального ветра
	 double mWind_VertV;
     // темпрература воздуха у поверхности земли по С
     double mAirT0;


	// конструктор по умолчанию
	   TEnvironment () ;
	// конструктор копирования
	   TEnvironment  (const TEnvironment  &R) ;
	// оператор присваивания
      TEnvironment  &operator=(const TEnvironment   &R2) ;

	 // парам конструктор
     TEnvironment (const int BallWave, const double Wind_V);

     TEnvironment (const double Wind_V, const double Wind_Alf, const double  Wind_VertV);

     TEnvironment (const double Wind_V, const double Wind_Alf
                                 , const double  Wind_VertV, const double  AirT0 );

     double calcAirRelativeDensity(const double  HeliH);

     double calcAirDensity(const double  HeliH);

     long double calcAirDensity(const long double  HeliH);

     void createVectWindV(long double *arrWindV);

     void createVectWindV( double *arrWindV);

     int createInputDataReport(wchar_t*FileName, const bool bHeader);
}  ;
#endif
