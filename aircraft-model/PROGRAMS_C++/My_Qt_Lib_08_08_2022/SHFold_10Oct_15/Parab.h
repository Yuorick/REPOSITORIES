//---------------------------------------------------------------------------

#ifndef ParabH
#define ParabH
class TURPointXY ;
class TParab
{
public:

 double ma; // точка на оси X
 double ml ;  //  сумма расстояний



   TParab();
 // конструктор копирования
 TParab (const TParab &R) ;
 // оператор присваивания
 TParab operator=(TParab  R);
 // параметр констр
 TParab ( const double a, const double b);

 void  clcCoeff ( double &vala, double &valb, double &valc) ;
 double  TParab::clcFncVal ( const double valx) ;

 int  fncLineIntersectParab (TURPointXY *pPontsInp,TURPointXY *pPontsOut) ;
 int fncSegmentCutParab( TURPointXY pntSegm0,  TURPointXY pntSegm1
							,TURPointXY *arrPntRez // точек ассив а пересечения
								) ;



};
#endif
