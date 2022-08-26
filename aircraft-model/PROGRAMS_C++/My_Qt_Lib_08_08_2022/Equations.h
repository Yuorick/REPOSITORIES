//---------------------------------------------------------------------------

#ifndef EquationsH
#define EquationsH
class TCompLong;
class TComp;

 int SolvEq2(const double a,const double b,const double c,TComp &x1,TComp &x2);
 int SolvEq2( TComp cmpa, TComp cmpb,  TComp cmpc,TComp &x1,TComp &x2);
 int Cubic(double *x,double a,double b,double c) ;
 double Sign(double a);
 int fCompare (const void *e1, const void *e2) ;
 int   CalcProper_Numbers_R3(double * arrKInp, double *arrLamb);
 int SolvEq3(TComp *cmparrRoots,TComp a,TComp b,TComp c,TComp d);
 double polinom(const double x, double *arra, const int len);
 double polinomDeriv(const double x, double *arra, const int len);
 int SolvEq4(TComp *cmparrRoots, double *arra);
 int SolvEq4_(TComp *cmparrRoots, double *arra);
 TComp polinom(const TComp x, TComp *arra, const int len);
 TComp polinom(const TComp x, double *arra0, const int len);

 int SolvEq2(const long double a,const long double b,const long double c,TCompLong &x1,TCompLong &x2);
 int Cubic(long double *x,long double a,long double b,long double c) ;
 long double Sign(long double a);
 int   CalcProper_Numbers_R3(long double * arrKInp, long double *arrLamb);
 int SolvEq5(TComp *cmparrRoots, double *arra, const double VAlMin, const double VAlMax);
 void calcPolinom_Coef_4Degree_Vieta( TComp*CmpArrLamb, double *pb, double *pc, double *pd, double *pe);



 int SolvEq2(const float a,const float b,const float c,TComp &x1,TComp &x2);

 int Cubic(float *x,float a,float b,float c) ;
 float Sign(float a);

 int   CalcProper_Numbers_R3(float * arrKInp, float *arrLamb);



#endif
