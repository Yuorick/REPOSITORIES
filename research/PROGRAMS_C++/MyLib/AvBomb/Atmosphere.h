//---------------------------------------------------------------------------

#ifndef AtmosphereH
#define AtmosphereH
#define ATM_PN0  1000.  //нормальное наземное давление воздуха
#define ATM_RoN0  1.2054  // нормальная назеная плотность воздуха
#define ATM_AN0  340.8   //нормальная назеная скорость звука в воздухе
#define ATM_TAYN0  289   //нормальная виртуальная температура на 0 высоте
#define ATM_R_UNIVER  287.05287   //удельная газовая постоянная воздуха
#define R_ZEMLI 6371000.
#define G_ZEMLI 9.80665
void fncCalcNormTemperature(const long double valy,long double &valTay, long double &valGradTay)  ;
long double fncCalcOtnositP(const long double H) ;


#endif
