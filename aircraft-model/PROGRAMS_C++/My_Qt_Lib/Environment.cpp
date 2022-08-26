
#include <math.h>
#include "Environment.h"
#include <wchar.h>
#include <stdio.h>
#include <time.h>

TEnvironment::TEnvironment()
{
	// Волнение на море по шкале Significance Wave Height (SWH). от 0 до 9 баллов
	 mBallWave = 0.;
//  Скорость ветра  горизонтальная  м/с
	 mWind_V = 0. ;

//  Скорость ветра  вертикальная  м/с
	 mWind_VertV = 0. ;
	 // направление откуда дует веьтер
	 mWind_Alf = 0. ;

     mAirT0 = 18.;

}

//---------------------------------------------------------------------------


// конструктор копирования
    TEnvironment ::TEnvironment (const TEnvironment &R)
 {
   mBallWave = R.mBallWave;
   mWind_V = R.mWind_V ;
   mWind_Alf = R.mWind_Alf;
   mWind_VertV = R.mWind_VertV;
   mAirT0 = R.mAirT0;
 }
 // оператор присваивания
 TEnvironment &TEnvironment::operator=(const TEnvironment  &R)
 {
   mBallWave = R.mBallWave;
   mWind_V = R.mWind_V ;
   mWind_Alf = R.mWind_Alf;
   mWind_VertV = R.mWind_VertV;
   mAirT0 = R.mAirT0;
	return *this ;
 }


  // парам конструктор
    TEnvironment::TEnvironment (const int BallWave, const double Wind_V)
 {
	mBallWave  = BallWave ;
	mWind_V  = Wind_V  ;
	mWind_VertV = 0.;
    mAirT0 = 18.;

 }

    // парам конструктор
      TEnvironment::TEnvironment (const double Wind_V, const double Wind_Alf, const double  Wind_VertV)
   {
      mWind_V  = Wind_V  ;
      mWind_Alf = Wind_Alf;
      mBallWave  = int (Wind_V / 20. * 9. +0.0001)  ;
      mWind_VertV = Wind_VertV;
      mAirT0 = 18.;

   }

//-------------------------------------------------------------------
// парам конструктор
TEnvironment::TEnvironment (const double Wind_V, const double Wind_Alf
                            , const double  Wind_VertV, const double  AirT0 )
{
    mWind_V  = Wind_V  ;
    mWind_Alf = Wind_Alf;
    mBallWave  = int (Wind_V / 20. * 9. +0.0001)  ;
    mWind_VertV = Wind_VertV;
    mAirT0 = AirT0 ;

}

//-------------------------------------------------------------

    // HeliH -  высота
double TEnvironment::calcAirRelativeDensity(const double  HeliH)
{
    double valTKel0 = 273.15 + mAirT0;
    double valTKel = valTKel0  - 0.0065 * HeliH;

    //double valPAm = 101325. * exp(log(valTKel/ valTKel0) * 9.81* 0.028964/ 8.31/ 0.0065);
    //double valRo =   valPAm * 0.028964/8.31/  valTKel;
    return pow(valTKel/valTKel0, 4.395);
}
//-------------------------------------------------------------

// HeliH -  высота
double TEnvironment::calcAirDensity(const double  HeliH)//, double  temp )(const double TKel0,const double HeliH)
{
    return ATM_RoN0 *calcAirRelativeDensity(HeliH);

   //return ATM_RoN0 *calcAirRelativeDensity(0.);
}

// HeliH -  высота
long double TEnvironment::calcAirDensity(const long double  HeliH)//, double  temp )(const double TKel0,const double HeliH)
{
    long double valTKel0 = 273.15 + mAirT0;
    long double valTKel = valTKel0  - 0.0065 * HeliH;


    long double RelDens =  powl(valTKel/valTKel0, 4.395);
    return ATM_RoN0 *RelDens ;


}

// вектор скорости ветра в ГСК
void TEnvironment::createVectWindV(long double *arrWindV)
{
   arrWindV[0]= -mWind_V * sin(mWind_Alf);
   arrWindV[1]= -mWind_V * cos(mWind_Alf);
   arrWindV[2] = -mWind_VertV;
}

// вектор скорости ветра в ГСК
void TEnvironment::createVectWindV( double *arrWindV)
{
   arrWindV[0]= -mWind_V * sin(mWind_Alf);
   arrWindV[1]= -mWind_V * cos(mWind_Alf);
   arrWindV[2] = -mWind_VertV;
}
//---------------------------------
int TEnvironment::createInputDataReport(wchar_t*FileName, const bool bHeader)

{
    int len = wcslen(FileName) ;

    if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'x') // проверка, что
     && (FileName[len - 3] == 't') ) )  // указанный файл имеет расширение .flt
    {

      return 1 ;
    }

    FILE *fw ;

    if ((fw = _wfopen(FileName,L"a"))== NULL)

    {

     return 1 ;
    }
if (bHeader)
{
   fprintf(fw,"  Дата и время формирования отчета\n");
   time_t t = time(NULL);
   struct tm* aTm = localtime(&t);
   fprintf(fw,"  Год = %04d\n",aTm->tm_year+1900);
   fprintf(fw,"  Mесяц = %02d\n",aTm->tm_mon+1);
   fprintf(fw,"  День = %02d\n",aTm->tm_mday);
   fprintf(fw,"  Время = %02d:%02d:%02d\n",aTm->tm_hour, aTm->tm_min, aTm->tm_sec);
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
}
fprintf(fw,"***************************************\n");
fprintf(fw,"      Параметры   атмосферы\n");

fprintf(fw,"  Сила ветра (м/с) = %5.2f\n",mWind_V);
fprintf(fw,"  направление откуда дует ветер, отсчитывается  \n");
fprintf(fw,"  от направления на Север по часовой стрелке (град) = %5.2f\n",mWind_Alf *180./M_PI);
fprintf(fw,"  скорость вертикального ветра (м/с) = %5.2f\n",mWind_VertV);
fprintf(fw,"  темпрература воздуха у поверхности земли по С (град) = %4.1f\n",mAirT0);
fclose(fw);
}


