//---------------------------------------------------------------------------

#ifndef ArcEllipseH
#define ArcEllipseH
class TURPointXY ;
class TArcEllipse
{
public:

 double ma; // большая полуось
 double mb ;  //  малая полуось
 double m_x0 ;
 double m_x1 ;


   TArcEllipse();
 // конструктор копирования
 TArcEllipse (const TArcEllipse &R) ;
 // оператор присваивания
 TArcEllipse operator=(TArcEllipse  R);
 // параметр констр
 TArcEllipse ( const double a, const double b);
 // параметр констр
 TArcEllipse( const double l, const double f, double *p);
 // парам констр
 TArcEllipse( const double a, const double b, const double x0, const double x1);
  // парам констр 3
 TArcEllipse( const double a, const double b, const int isign);
  // парам констр 4  создания из сумы расстояний до фокусов  и фокального расстояния
  TArcEllipse( const double l, const double f,  const int isign, double *p) ;



 double  fncFocus ();



 int  fncLineIntersectEllips (TURPointXY *pPontsInp,TURPointXY *pPontsOut) ;

 double calcFncValue( const double valx);
 double calc_dY_po_dX( const double valx);
 void ShowMe(wchar_t *FileName);
 void ShowFullGraph(wchar_t *FileName);
 int fncSegmentCutArcEllipse(TURPointXY pntSegm0, TURPointXY pntSegm1 // сегмент
							,TURPointXY *arrPntRez1 // точек ассив а пересечения
								);

};

bool Is_X_BelongeSegm( const double valx, const double vala,const double valb);
#endif
