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



 int SolvEq2(const long double a,const long double b,const long double c,TCompLong &x1,TCompLong &x2);
 int Cubic(long double *x,long double a,long double b,long double c) ;
 long double Sign(long double a);
 int   CalcProper_Numbers_R3(long double * arrKInp, long double *arrLamb);






#endif
