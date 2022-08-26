//---------------------------------------------------------------------------

#ifndef AtmosphereH
#define AtmosphereH
#define ATM_PN0  1000.  //нормальное наземное давление воздуха
#define ATM_RoN0  1.2054  // нормальная назеная плотность воздуха
#define ATM_AN0  340.8   //нормальная назеная скорость звука в воздухе
#define ATM_TAYN0  289   //нормальная виртуальная температура на 0 высоте
#define ATM_R_UNIVER  287.05287   //удельная газовая постоянная воздуха
void fncCalcNormTemperature(const  double valy, double &valTay,  double &valGradTay)  ;

#endif
