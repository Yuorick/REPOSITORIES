#ifndef DrawShapeH
#define DrawShapeH

//---------------------------------------------------------------------------
#include <vcl.h>


#include "URPointXY.h"
#include "URPolyLine.h"
#include "URPolygon.h"


#define MAX_NPTS 5000
#define MAX_NPLL 100
#define MAX_NPLG 100

class TDrawShape
{ public:

	TURPointXY  Pts[MAX_NPTS];    // массив рисуемых точек
	TURPolyLine Pll[MAX_NPLL];    // массив рисуемых полиний
	TURPolygon  Plg[MAX_NPLG];    // массив рисуемых полигонов
	int nPts;         // к-во точек
	int nPll;         // к-во полиний
	int nPlg;         // к-во полигонов

	double Box[4];      // ограничивающий фрейм-бокс  {xmin,ymin,xmax,ymax}


	bool bPts[MAX_NPTS];    // признаки присутствия на холсте
	bool bPll[MAX_NPLL];    // признаки присутствия на холсте
	bool bPlg[MAX_NPLG];    // признаки присутствия на холсте

	TColor ColPts[MAX_NPTS];   	// цвета точек
	TColor ColPll[MAX_NPLL];   	// цвета полилиний
	TColor ColPlg[MAX_NPLG];   	// цвета полигонов

	int RPts[MAX_NPTS];        	// толщины точек
	int RPll[MAX_NPLL];        	// толщины полиний
	int RPlg[MAX_NPLG];        	// толщины полигонов


	// члены визуализации
	TImage* Holst; 			// орган управления TImage
	int Wframe;         // отступ рамки

	double dkScale;     // коэфф.масштабирования
	double dkx,dky;     // доп.коэфф.масштабир-я отдельно по X и Y



  bool IS_FRAME_NEED;
	bool IS_AXES_NEED;

	TDrawShape(void);
	TDrawShape(	TImage* Im);
	//TDrawShape(TDrawShape &dsh);


	void CalcFrameBox();
	void DrawPoints();
	void DrawPolylines();
	void DrawPolygons();
	void DrawObjects();


	void Save(AnsiString fn);     // сохранить объект в виде файла draw
	void Open(AnsiString fn);     // считать объект из файла draw


void __fastcall TDrawShape::DrawTimeFun
														(	double *pF,   //  массив значений ф-ии времени - ось ординат
															double *pT,   //  массив замеров времени - ось абсцисс
															int N,   			//  размер массивов
															double &dT,   //  масштабирование по абсциссе
															double &dF    //  масштабирование по ординате
														);
double __fastcall TDrawShape::	maxDoubleArr
																(	double *parr,
																	int N,
																	int &irez
																);
double __fastcall TDrawShape::	minDoubleArr
																(	double *parr,
																	int N,
																	int &irez
																);

void __fastcall TDrawShape::	WriteOneReport
															(	double *mtxData, 	// матрица данных nRows x nCols
																const int nCols, 	// - к-во функций
																const int nRows, 	// - к-во значений абсцисс
																wchar_t *Names, 		// массив с именаими переменных - матрица nCols X lenName
																int lenName, // макс.длина имени переменной
																int nX,  	// № столбца переменной, явл-ся абсциссой
																int nY,  	// № столбца переменной, явл-ся ординатой
																double dX,  	// масштаб по абсциссе
																double dY  	// масштаб по ординате
															);

void  TDrawShape::	DrawAxes(	double xmin,
															double xmax,
															double ymin,
															double ymax
														);

void TDrawShape::		DrawAxes(	double xmin,
															double xmax,
															double ymin,
															double ymax,
															double Len
														);


void TDrawShape::		DrawAxes(	double xmin,
															double xmax,
															double ymin,
															double ymax,
															double Len,
															TURPointXY XY0    // начало координат
														);
void TDrawShape::	 CreateAngleMarks
									 ( TURPointXY P0,  // вершина угла
										 TURPointXY P1,  // 1-я точка
										 TURPointXY P2,  // 2-я точка
										 double D0,      // расстояние до 1-го сектора
										 double D1       // расстояние до 2-го сектора
									 );
double TDrawShape:: SGN(double x);

void  TDrawShape:: 	ShowNormProbDistr( double A,double Sig2);

void  TDrawShape:: 	PictFar();


TDrawShape operator = (TDrawShape dsh);
TDrawShape operator = (TDrawShape& dsh);

TDrawShape operator + (TURPolyLine pll);
TDrawShape operator += (TURPointXY pt);
TDrawShape operator += (TURPolyLine pll);
TDrawShape operator += (TURPolygon plg);
void TDrawShape::	ClearPts(void);
void TDrawShape::	ClearPll(void);
void TDrawShape::	ClearPlg(void);
void TDrawShape::	ClearAll(void);


};
//--------------------------------------------------------------------------

#endif
