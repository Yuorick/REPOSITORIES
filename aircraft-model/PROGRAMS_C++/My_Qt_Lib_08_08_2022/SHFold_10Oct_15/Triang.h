//---------------------------------------------------------------------------

#ifndef TriangH
#define TriangH


#endif
 #include "UrPointXY.h"
 #include <math.h>
class TURPolygon;
class TURPointXY;
class TTriang
{
public:
	 TURPointXY m_pVert[4];
	 TTriang() ;
      // конструктор копирования
	  TTriang  (const TTriang  &R) ;
	  // оператор присваивания
      TTriang  &operator=(const TTriang   &R2) ;
	// парам констр
	TTriang( TURPointXY  *arrPoints);
	// парам констр
 //	TTriang( TURPointXY  p0,TURPointXY  p1,TURPointXY  p2);
	// парам констр
	TTriang( TURPointXY  p0,TURPointXY  p1,TURPointXY  p2);
	// площадь треугольника
	static  double calcVectS(TURPointXY P1,TURPointXY P2) ;
	double calcSq();

	static double max_d(const double x0, const double x1 ) ;
    static double min_d(const double x0, const double x1 );

	static void SubTwoTrsToTrs( TTriang *ptrT0, TTriang *ptrT1, TTriang *pTr,int &quantTr ) ;

	static void SubTwoTrsToRings( TTriang *ptrT0, TTriang *ptrT1, TURPolygon *urplgArr,int &quantPlg ) ;

	static bool IntersectTwoSegments(
					const TURPointXY p00,const  TURPointXY p01 //1-ый сегмент
					,const TURPointXY p10,const TURPointXY p11  // 2-ой сегмент
					, TURPointXY *p0  // точка пересечения

					, int *ipCaseType //
					);
   int PntInTriangle(  TURPointXY p0);

 static	void  InsertFalse(bool *barrOut,const int i0 ,const bool b);

 static	void InsertVert(TURPolygon *plg, const int i0 ,const  TURPointXY p)  ;

 static void InsertNum(int *iarrNums,const int i,const int j);

 static int urSign(double a);

 static int PntInSegment(double a,double x);

 static int IntersectTwoSegmInLine(double a0,double b0,double b1);

 static double dist(TURPointXY p0 ,TURPointXY p1) ;
 static bool IsTrianglesIntersectsLine(TTriang *ptrT0,const int i0,TTriang *ptrT1,const int i1) ;
 void CutTriangle(TURPointXY p0 ,TURPointXY p1, int *iarrNumEdge,TURPolygon *purplgRings );
 static void SubTwoEnclosedTrToRings( TTriang *trT0, TTriang *trT1, TURPolygon *urplgArr ) ;
 void ChangeDir();
static  void FinishDisassembleTriangleArr(TTriang **ppTr,int *lenTrArr) ;
static int TypeOfTrianglesInersection( TTriang *pTr0, TTriang *pTr1
		,TURPointXY *urpntP, int *i0,int * i1) ;
static int Signum(const double a)
{
	if (a > 0) return 1;
	else
	if (a < 0) return -1;
	else return 0;
}

static double dist2( TURPointXY *p0, TURPointXY *p1)
{
 return (((*p0).X - (*p1).X) * ((*p0).X - (*p1).X) + ((*p0).Y - (*p1).Y) * ((*p0).Y - (*p1).Y) );
}
static double dist(TURPointXY *p0,TURPointXY *p1)
{
 return sqrt(dist2(p0,p1));
}


};
