//---------------------------------------------------------------------------

#ifndef ArcParabH
#define ArcParabH

// Парабола это геометр место точек равноудалённых от данной прямой (называемой директрисой параболы)
//  и данной точки (называемой фокусом параболы).
// Наряду с эллипсом и гиперболой, парабола является коническим сечением.
// Она может быть определена как коническое сечение с единичным эксцентриситетом.
///////////////////////////////////////////////////////////////////////////////////
// парабола описываемая этим классом
// ml - фокальный параметр
// ml/2 - фокусное расстояние
// ma - координата вершины параболы по оси OX
// парабола типа X*X = -2PY
// Уравнение параболы
// Y=-(X-ma)*(X-ma)/(2*ml)+ ml/2
// то есть фокус расположен на оси OX
// координаты вершины (ma, ml/2)
// ветви направлены вниз
class TURPointXY ;
class TArcParab
{
public:

 double ma; // точка на оси X
 double ml ;  //  сумма расстояний
 double m_x0 ;
 double m_x1 ;


   TArcParab();
 // конструктор копирования
 TArcParab (const TArcParab &R) ;
 // оператор присваивания
 TArcParab &operator=(const TArcParab  &R);
 // параметр констр
 TArcParab ( const double a, const double b);
 	// парам констр 2
TArcParab( const double a, const double l, const int isign) ;


 void  clcCoeff ( double &vala, double &valb, double &valc) ;
 double  calcFncValue ( const double valx) ;

 int  fncLineIntersectParab (TURPointXY pPontsInp0,TURPointXY pPontsInp1,TURPointXY *pPontsOut);

 int fncSegmentCutArcParab( TURPointXY pntSegm0, TURPointXY pntSegm1 // сегмент
							,TURPointXY *arrPntRez // точек ассив а пересечения
								);
  void ShowMe(wchar_t *FileName);
  int fncSegmentCutParabLine(TURPointXY pntSegm0, TURPointXY pntSegm1 // сегмент
							,TURPointXY *arrPntRez // точек ассив а пересечения
								) ;
  void ShowFullGraph(wchar_t *FileName);



};

 bool Is_X_BelongeSegm( const double valx, const double vala,const double valb);
 double max_( const double x0, const double x1);
 double min_( const double x0, const double x1);


#endif
