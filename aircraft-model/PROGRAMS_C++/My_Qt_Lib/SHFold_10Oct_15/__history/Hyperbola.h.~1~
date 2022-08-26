//---------------------------------------------------------------------------

#ifndef HyperbolaH
#define HyperbolaH
// уравнение в канонической форме:
// x*x/a/a - y*y/b/b =0
//
class TURPointXY ;
class THyperbola
{
public:

 double ma; // большая полуось
 double mb ;  //  малая полуось


   THyperbola();
 // конструктор копирования
 THyperbola (const THyperbola &R) ;
 // оператор присваивания
 THyperbola operator=(THyperbola  R);
 // параметр констр
 THyperbola ( const double a, const double b);



 double  fncFocus ();



 void ShowMe(wchar_t *FileName, const int iNumPoints0, const double VAlDiap);
 TURPointXY  ProjectPointOnHyperbola (TURPointXY  pntInp);
 double fncTemp(double q, double xZv, double yZv, TURPointXY  pntInp);
 int findTangencyPoints(TURPointXY  pntInp, TURPointXY  &pntOut0, TURPointXY  &pntOut1);

};

bool chekRoot (TURPointXY  pntInp0, TURPointXY pntOut);
#endif
