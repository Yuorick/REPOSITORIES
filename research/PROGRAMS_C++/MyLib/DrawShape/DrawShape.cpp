#include "DrawShape.h"
#include <float.h>

#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define min(a, b)  (((a) < (b)) ? (a) : (b))

// ---------------------------------------------------------------------------
TDrawShape::TDrawShape(void) {
	Holst = NULL;
	// X0= Y0= 0;
	Box[0] = Box[1] = Box[2] = Box[3] = 0;
	// Xmin= Ymin= 0;
	memset(ColPts, 0, sizeof(ColPts));
	memset(ColPll, 0, sizeof(ColPll));
	memset(ColPlg, 0, sizeof(ColPlg));
	memset(RPts, 0, sizeof(RPts));
	memset(RPll, 0, sizeof(RPll));
	memset(RPlg, 0, sizeof(RPlg));

	dkScale = dkx = dky = 1.0; // коэфф.масштабирования

	int nPts = 0; // к-во точек
	int nPll = 0; // к-во полиний
	int nPlg = 0; // к-во полигонов
	IS_FRAME_NEED = true;
	IS_AXES_NEED = true;

}

// ---------------------------------------------------------------------------
TDrawShape::TDrawShape(TImage* Im) {
	Holst = Im;
	// X0= Y0= 0;
	Box[0] = Box[1] = Box[2] = Box[3] = 0;
	memset(ColPts, 0, sizeof(ColPts));
	memset(ColPll, 0, sizeof(ColPll));
	memset(ColPlg, 0, sizeof(ColPlg));
	memset(RPts, 0, sizeof(RPts));
	memset(RPll, 0, sizeof(RPll));
	memset(RPlg, 0, sizeof(RPlg));
	dkScale = dkx = dky = 1.0; // коэфф-ты масштабирования
	int nPts = 0; // к-во точек
	int nPll = 0; // к-во полиний
	int nPlg = 0; // к-во полигонов

	IS_FRAME_NEED = true;
	IS_AXES_NEED = true;

}

// ---------------------------------------------------------------------------
void TDrawShape::CalcFrameBox() {
	int i, k, n;
	double xmin = 10000000, xmax = -10000000, ymin = 10000000, ymax = -10000000;

	Box[0] = Box[1] = Box[2] = Box[3] = 0; // огр-щий фрейм-бокс  {xmin,ymin,xmax,ymax}

	// ТОЧКИ
	if (nPts > 0) {
		if (bPts[0]) {
			xmin = Pts[0].X;
			xmax = Pts[0].X;
			ymin = Pts[0].Y;
			ymax = Pts[0].Y;
		}
		for (i = 1; i < nPts; i++) {
			if (!bPts[i])
				continue;
			xmin = min(xmin, Pts[i].X);
			xmax = max(xmax, Pts[i].X);
			ymin = min(ymin, Pts[i].Y);
			ymax = max(ymax, Pts[i].Y);
		}
	}
	// ПОЛИЛИНИИ
	if (nPll > 0) {
		if (bPll[0]) {
			xmin = min(xmin, Pll[0].Box[0]);
			xmax = max(xmax, Pll[0].Box[2]);
			ymin = min(ymin, Pll[0].Box[1]);
			ymax = max(ymax, Pll[0].Box[3]);
		}
		for (i = 1; i < nPll; i++) {
			if (!bPll[i])
				continue;
			xmin = min(xmin, Pll[i].Box[0]);
			xmax = max(xmax, Pll[i].Box[2]);
			ymin = min(ymin, Pll[i].Box[1]);
			ymax = max(ymax, Pll[i].Box[3]);
		}
	}
	// ПОЛИГОНЫ
	if (nPlg > 0) {
		if (bPlg[0]) {
			xmin = min(xmin, Plg[0].Box[0]);
			xmax = max(xmax, Plg[0].Box[2]);
			ymin = min(ymin, Plg[0].Box[1]);
			ymax = max(ymax, Plg[0].Box[3]);
		}
		for (i = 1; i < nPlg; i++) {
			if (!bPlg[i])
				continue;
			xmin = min(xmin, Plg[i].Box[0]);
			xmax = max(xmax, Plg[i].Box[2]);
			ymin = min(ymin, Plg[i].Box[1]);
			ymax = max(ymax, Plg[i].Box[3]);
		}
	}
	Box[0] = xmin;
	Box[1] = ymin;
	Box[2] = xmax;
	Box[3] = ymax;

}

// ---------------------------------------------------------------------------
void TDrawShape::DrawPoints() {
	TRect rec;
	int i, n, k;
	double sdvX, sdvY;
	int W, H;
	int x1, y1, x2, y2;

	CalcFrameBox();

	W = (int)(dkScale * dkx * (double)(Box[2] - Box[0]));
	H = (int)(dkScale * dky * (double)(Box[3] - Box[1]));

	Holst->Width = W + 2 * Wframe;
	Holst->Height = H + 2 * Wframe;

	sdvX = -Box[0];
	sdvY = -Box[1];

	// --------------------------------------
	// TRect r(0,0,5000,2000);
	// Holst->Canvas->Brush->Color= clWhite;
	// Holst->Canvas->FillRect(r);
	// --------------------------------------

	for (n = 0; n < nPts; n++)
	{ if(!bPts[n]) continue;
		Holst->Canvas->Pen->Color = ColPts[n];
		Holst->Canvas->Brush->Color = ColPts[n];

		x1 = (Pts[n].X + sdvX) * dkScale * dkx;
		y1 = Holst->Height - (int)((Pts[n].Y + sdvY) * dkScale * dky);
		rec.left = x1 - RPts[i] + Wframe;
		rec.right = x1 + RPts[i] + Wframe;
		rec.top = y1 - RPts[i] - Wframe;
		rec.bottom = y1 + RPts[i] - Wframe;
		Holst->Canvas->Ellipse(rec.left, rec.top, rec.right, rec.bottom);
	}

}

// ------------------------------------------------------------------------------
void TDrawShape::DrawPolylines()
{ int i, n, k;
	double sdvX, sdvY;
	int W, H;
	int x1, y1, x2, y2;

	CalcFrameBox();

//	Box[0]= Box[1]= Box[2]= Box[3]= 0; // огр-щий фрейм-бокс  {xmin,ymin,xmax,ymax}


	W = (int)(dkScale * dkx * (double)(Box[2] - Box[0]));
	H = (int)(dkScale * dky * (double)(Box[3] - Box[1]));

	Holst->Width = W + 2 * Wframe;
	Holst->Height = H + 2 * Wframe;

	sdvX = -Box[0];
	sdvY = -Box[1];

	// --------------------------------------
	// TRect r(0,0,5000,2000);
	// Holst->Canvas->Brush->Color= clWhite;
	// Holst->Canvas->FillRect(r);
	// --------------------------------------

	for (n = 0; n < nPll; n++)
	{ if(!bPll[n]) continue;
		Holst->Canvas->Pen->Color = ColPll[n];
		Holst->Canvas->Brush->Color = ColPll[n];
		Holst->Canvas->Pen->Width = RPll[n];
		for(i = 0; i < Pll[n].NumParts - 1; i++)
		{ for(k = Pll[n].Parts[i]; k < Pll[n].Parts[i + 1]; k++)
			{ x1= (Pll[n].Points[k].X + sdvX) * dkScale * dkx;
				y1= Holst->Height - (int)((Pll[n].Points[k].Y + sdvY) * dkScale * dky);
				x2= (Pll[n].Points[k + 1].X + sdvX) * dkScale * dkx;
				y2= Holst->Height-(int)((Pll[n].Points[k+1].Y+sdvY)*dkScale*dky);

				Holst->Canvas->MoveTo(x1 + Wframe, y1 - Wframe);
				Holst->Canvas->LineTo(x2 + Wframe, y2 - Wframe);
			}
		}
		for(k= Pll[n].Parts[Pll[n].NumParts-1];k<Pll[n].NumPoints-1;k++)
		{ x1= (Pll[n].Points[k].X+sdvX)*dkScale*dkx;
			y1= Holst->Height - (int)((Pll[n].Points[k].Y + sdvY) * dkScale * dky);
			x2= (Pll[n].Points[k + 1].X + sdvX) * dkScale * dkx;
			y2= Holst->Height-(int)((Pll[n].Points[k+1].Y+sdvY)*dkScale*dky);
			Holst->Canvas->MoveTo(x1+ Wframe,y1- Wframe);
			Holst->Canvas->LineTo(x2+ Wframe,y2- Wframe);
		}
	}

	if (IS_FRAME_NEED) // нарисовать оси коор-т
	{ Holst->Canvas->Pen->Color= clBlack;
		Holst->Canvas->Pen->Width= 1;
		Holst->Canvas->Brush->Color= clBlack;

		x1= (Box[0]+sdvX)*dkScale*dkx; // xmin
		y1= Holst->Height-(int)((Box[1]+sdvY)*dkScale*dky); // ymin
		x2= (Box[2]+sdvX)*dkScale*dkx; // xmax
		y2= Holst->Height-(int)((Box[3]+sdvY)*dkScale*dky); // ymax

		// Рамка
		Holst->Canvas->MoveTo(x1 + Wframe, y1 - Wframe);
		Holst->Canvas->LineTo(x1 + Wframe, y2 - Wframe);
		Holst->Canvas->LineTo(x2 + Wframe, y2 - Wframe);
		Holst->Canvas->LineTo(x2 + Wframe, y1 - Wframe);
		Holst->Canvas->LineTo(x1 + Wframe, y1 - Wframe);
	}
	if (IS_AXES_NEED) // нарисовать оси коор-т
	{ // Сами оси
		x1 = sdvX * dkScale * dkx; // это начало координат
		y1 = Holst->Height - (int)(sdvY * dkScale * dky); // это начало координат

		Holst->Canvas->Pen->Color = clRed;
		Holst->Canvas->Brush->Color = clRed;

		Holst->Canvas->MoveTo(Wframe, y1 - Wframe);
		Holst->Canvas->LineTo(Holst->Width - Wframe, y1 - Wframe);

		Holst->Canvas->MoveTo(x1 + Wframe, Wframe);
		Holst->Canvas->LineTo(x1 + Wframe, Holst->Height - Wframe);
	}
}

// ------------------------------------------------------------------------------
void TDrawShape::DrawPolygons() {
	int i, n, k;
	double sdvX, sdvY;
	int W, H;
	int x1, y1, x2, y2;

	CalcFrameBox();

	W = (int)(dkScale * dkx * (double)(Box[2] - Box[0]));
	H = (int)(dkScale * dky * (double)(Box[3] - Box[1]));

	Holst->Width = W + 2 * Wframe;
	Holst->Height = H + 2 * Wframe;

	sdvX = -Box[0];
	sdvY = -Box[1];

	// --------------------------------------
	// TRect r(0,0,5000,2000);
	// Holst->Canvas->Brush->Color= clWhite;
	// Holst->Canvas->FillRect(r);
	// --------------------------------------

	for (n = 0; n < nPlg; n++) {
		if (!bPlg[n])
			continue;
		Holst->Canvas->Pen->Color = ColPlg[n];
		Holst->Canvas->Brush->Color = ColPlg[n];
		Holst->Canvas->Pen->Width = RPlg[n];
		for (i = 0; i < Plg[n].NumParts - 1; i++) {
			for (k = Plg[n].Parts[i]; k < Plg[n].Parts[i + 1]; k++) {
				x1 = (Plg[n].Points[k].X + sdvX) * dkScale * dkx;
				y1 = Holst->Height - (int)((Plg[n].Points[k].Y + sdvY) * dkScale * dky);
				x2 = (Plg[n].Points[k + 1].X + sdvX) * dkScale * dkx;
				y2 = Holst->Height - (int)((Plg[n].Points[k + 1].Y + sdvY)
					* dkScale * dky);

				Holst->Canvas->MoveTo(x1 + Wframe, y1 - Wframe);
				Holst->Canvas->LineTo(x2 + Wframe, y2 - Wframe);
			}
		}
		for (k = Plg[n].Parts[Plg[n].NumParts - 1]; k < Plg[n].NumPoints - 1; k++) {
			x1 = (Plg[n].Points[k].X + sdvX) * dkScale * dkx;
			y1 = Holst->Height - (int)((Plg[n].Points[k].Y + sdvY) * dkScale * dky);
			x2 = (Plg[n].Points[k + 1].X + sdvX) * dkScale * dkx;
			y2 = Holst->Height - (int)((Plg[n].Points[k + 1].Y + sdvY)
				* dkScale * dky);

			Holst->Canvas->MoveTo(x1 + Wframe, y1 - Wframe);
			Holst->Canvas->LineTo(x2 + Wframe, y2 - Wframe);
		}
	}

	/*
	if(IS_FRAME_NEED) // нарисовать оси коор-т
	{
	Holst->Canvas->Pen->Color  = clBlack;
	Holst->Canvas->Pen->Width  = 1;
	Holst->Canvas->Brush->Color= clBlack;

	x1= (Box[0]+sdvX)*dkScale*dkx;              // xmin
	y1= Holst->Height-(int)((Box[1]+sdvY)*dkScale*dky);              // ymin
	x2= (Box[2]+sdvX)*dkScale*dkx;              // xmax
	y2= Holst->Height-(int)((Box[3]+sdvY)*dkScale*dky);             // ymax

	// Рамка
	Holst->Canvas->MoveTo(x1+Wframe,y1-Wframe);
	Holst->Canvas->LineTo(x1+Wframe,y2-Wframe);
	Holst->Canvas->LineTo(x2+Wframe,y2-Wframe);
	Holst->Canvas->LineTo(x2+Wframe,y1-Wframe);
	Holst->Canvas->LineTo(x1+Wframe,y1-Wframe);
	}

	if(IS_AXES_NEED) // нарисовать оси коор-т
	{ // Сами оси
	x1= sdvX*dkScale*dkx;                          // это начало координат
	y1= Holst->Height- (int)(sdvY*dkScale*dky);    // это начало координат

	//Holst->Canvas->Pen->Color  = ColPlg;
	//Holst->Canvas->Brush->Color= ColPlg;
	}
	 */
}

// ------------------------------------------------------------------------------
void TDrawShape::DrawObjects() { // Установка размеров Image1
	/*
	int i;
	int xmin,xmax,ymin,ymax;
	int W,H;

	if(nPts>0)
	{	xmin= xmax=Pts[0].X; ymin= ymax=Pts[0].Y;
	for(i=1;i<nPts;i++)
	{ xmin= min(xmin,Pts[i].X); xmax= max(xmax,Pts[i].X);
	ymin= min(ymin,Pts[i].Y); ymax= max(ymax,Pts[i].Y);
	}
	}
	if(nPll>0)
	{ xmin= min(xmin,Pll[0].Box[0]); xmax= max(xmax,Pll[0].Box[2]);
	ymin= min(ymin,Pll[0].Box[1]); ymax= max(ymax,Pll[0].Box[3]);
	for(i=1;i<nPll;i++)
	{	xmin= min(xmin,Pll[i].Box[0]); xmax= max(xmax,Pll[i].Box[2]);
	ymin= min(ymin,Pll[i].Box[1]); ymax= max(ymax,Pll[i].Box[3]);
	}
	}
	if(nPlg>0)
	{ xmin= min(xmin,Plg[0].Box[0]); xmax= max(xmax,Plg[0].Box[2]);
	ymin= min(ymin,Plg[0].Box[1]); ymax= max(ymax,Plg[0].Box[3]);
	for(i=1;i<nPlg;i++)
	{	xmin= min(xmin,Plg[i].Box[0]); xmax= max(xmax,Plg[i].Box[2]);
	ymin= min(ymin,Plg[i].Box[1]); ymax= max(ymax,Plg[i].Box[3]);
	}
	}

	xmax= (int)((double)(xmax)*dkScale*dkx);
	xmin= (int)((double)(xmin)*dkScale*dkx);
	ymax= (int)((double)(ymax)*dkScale*dky);
	ymin= (int)((double)(ymin)*dkScale*dkx);

	W= xmax-xmin;
	H= ymax-ymin;
	Holst->Width= W+ 2*Wframe; Holst->Height= H+ 2*Wframe;

	//	Xmin= xmin;
	//	Ymin= ymin;
	 */
	// --------------------------------------
	// Очистка экрана

	// TRect r(0,0,5000,2000);

	Holst->Canvas->Brush->Color = clWhite;
	Holst->Canvas->FillRect(Holst->ClientRect);
	// --------------------------------------

	DrawPoints();
	DrawPolylines();
	DrawPolygons();

	// Рисование осей координат
	// Holst->Canvas->Pen->Color= clBlack;
	// Holst->Canvas->Pen->Width= 1;
	// Holst->Canvas->MoveTo(0, Holst->Height-(Y0));
	// Holst->Canvas->LineTo(Holst->Width, Holst->Height-(Y0));
	// Holst->Canvas->MoveTo(X0, 0);
	// Holst->Canvas->LineTo(X0, Holst->Height);

}

// ------------------------------------------------------------------------------
// Рисование функции времени на холсте
// ------------------------------------------------------------------------------
void __fastcall TDrawShape::DrawTimeFun(double *pF,
	// массив значений ф-ии времени - ось ординат
	double *pT, // массив замеров времени - ось абсцисс
	int N, // размер массивов
	double &dT, // масштабирование по абсциссе
	double &dF // масштабирование по ординате
		) {
	int i;
	TURPointXY *pPts = new TURPointXY[N];

	for (i = 0; i < N; i++) {
		pPts[i].Y = pF[i] * dF;
		pPts[i].X = pT[i] * dT;
	}

	// Инициализирую массив полилиний Pll
	// X0= Y0= 0;  				// начало координат  (pix)
	ColPll[0] = clRed; // цвет для рисования полилиний
	RPll[0] = 2; // толщина полилиний (pix)
	nPll = 1; // к-во рисуемых полилиний
	int iarrParts[1] = {
		0
	}; // массив сдвигов int

	// массив рисуемых полилиний
	TURPolyLine Pll0(1, N, iarrParts, pPts);
	Pll[0] = Pll0;
	nPll = 1;

	DrawPolylines();

	delete[]pPts;

}

// ------------------------------------------------------------------------------
// Возвращает первое значение массива большее (-DBL_MAX)
// irez содержит его индекс
double __fastcall TDrawShape::maxDoubleArr(double *parr, // массив чисел double
	int N, // размер массива
	int &irez // размер массивов
		) {
	int i;
	double rez = -DBL_MAX;
	for (i = 0; i < N; i++) {
		if (parr[i] > rez) {
			rez = parr[i];
			irez = i;
		}
	}
	return rez;
}

// ------------------------------------------------------------------------------
// Возвращает первое значение массива меньшее DBL_MAX
// irez содержит его индекс
double __fastcall TDrawShape::minDoubleArr(double *parr, int N, int &irez) {
	int i;
	double rez = DBL_MAX;
	for (i = 0; i < N; i++) {
		if (parr[i] < rez) {
			rez = parr[i];
			irez = i;
		}
	}
	return rez;
}

// ------------------------------------------------------------------------------
void __fastcall TDrawShape::WriteOneReport(double *mtxData,
	// матрица данных nRows x nCols
	const int nCols, // - к-во функций
	const int nRows, // - к-во значений абсцисс
	wchar_t *Names, // массив с именаими переменных - матрица nCols X lenName
	int lenName, // макс.длина имени переменной
	int nX, // № столбца переменной, явл-ся абсциссой
	int nY, // № столбца переменной, явл-ся ординатой
	double dX, // масштаб по абсциссе
	double dY // масштаб по ординате
		) {
	int i;
	// создание массива абсцисс и ординат
	TURPointXY *pPts = new TURPointXY[nRows];

	for (i = 0; i < nRows; i++) {
		pPts[i].X = mtxData[nCols * i + nX] * dX;
		pPts[i].Y = mtxData[nCols * i + nY] * dY;
	}
	// Инициализирую массив полилиний Pll
	// X0= Y0= 0;  				// начало координат  (pix)
	ColPll[0] = clRed; // цвет для рисования полилиний
	RPll[0] = 2; // толщина полилиний (pix)
	nPll = 1; // к-во рисуемых полилиний
	int iarrParts[1] = {
		0
	}; // массив сдвигов int

	// массив рисуемых полилиний
	TURPolyLine Pll0(1, nRows, iarrParts, pPts);
	Pll[0] = Pll0;

	nPll = 1;
	DrawPolylines();

	delete[]pPts;
}

// ------------------------------------------------------------------------------
void TDrawShape::DrawAxes(double xmin, double xmax, double ymin, double ymax) {
	int numParts = 2;
	int N = 4;
	int iarrParts[2] = {
		0, 2
	};
	TURPointXY pPts[4];

	pPts[0] = TURPointXY(xmin, 0);
	pPts[1] = TURPointXY(xmax, 0);
	pPts[2] = TURPointXY(0, ymin);
	pPts[3] = TURPointXY(0, ymax);
	// массив рисуемых полилиний
	TURPolyLine Pll0(1, N, iarrParts, pPts);
	Pll[0] = Pll0;
	nPll = 1;

	DrawPolylines();
}

// ------------------------------------------------------------------------------
// Создание  осей координат  со стрелками
void TDrawShape::DrawAxes(double xmin, double xmax, double ymin, double ymax,
	double Len) {
	TURPointXY pXL(xmin, 0);
	TURPointXY pXR(xmax, 0);
	TURPointXY pYB(ymin, 0);
	TURPointXY pYT(ymax, 0);

	TURPolyLine Axes = TURPolyLine::fncCreateAxes(pXL, pXR, pYB, pYT, Len);
	nPll = 1;
	Pll[0] = Axes;
	DrawPolylines();
}

// ------------------------------------------------------------------------------
// Создание осями коор-т со стрелками c нач.коор-т в точке XY0
void TDrawShape::DrawAxes(double xmin, double xmax, double ymin, double ymax,
	double Len, TURPointXY XY0 // начало координат
		) {
	TURPointXY pXL(xmin, 0);
	TURPointXY pXR(xmax, 0);
	TURPointXY pYB(ymin, 0);
	TURPointXY pYT(ymax, 0);

	TURPolyLine Axes1 = TURPolyLine::fncCreateAxes(pXL, pXR, pYB, pYT, Len);
	TURPolyLine Axes = Axes1.SdvigTransform(XY0);

	nPll = 1;
	Pll[0] = Axes;
	DrawPolylines();
}

// ------------------------------------------------------------------------------
// Создание  полилинии с секторами показывающими угол
void TDrawShape::CreateAngleMarks(TURPointXY P0, // вершина угла
	TURPointXY P1, // 1-я точка
	TURPointXY P2, // 2-я точка
	double D0, // расстояние до 1-го сектора
	double D1 // расстояние до 2-го сектора
		) {
	double fi0 = 0;
	double fi1 = 0;

	if (fabs(P1.X - P0.X) < 0.00001)
		fi0 = SGN(P1.Y - P0.Y) * M_PI / 2.0;
	else
		fi0 = atan2(P1.Y - P0.Y, P1.X - P0.X);

	if (fabs(P2.X - P0.X) < 0.00001)
		fi1 = SGN(P2.Y - P0.Y) * M_PI / 2.0;
	else
		fi1 = atan2(P2.Y - P0.Y, P2.X - P0.X);

	TURPolyLine Sect0 = TURPolyLine::fncCreateSector(P0, D0, fi0, fi1, 1500);
	TURPolyLine Sect1 = TURPolyLine::fncCreateSector(P0, D1, fi0, fi1, 1500);

	int iarrParts[2] = {
		0, 1500
	};
	const int N = 3000;
	const int iNumParts = 2;
	TURPointXY Pts[3000];

	memcpy(Pts, Sect0.Points, 1500 * sizeof(TURPointXY));
	memcpy(&Pts[1500], Sect1.Points, 1500 * sizeof(TURPointXY));
	TURPolyLine pllRez(2, N, iarrParts, Pts);

	nPll = 1;
	Pll[0] = pllRez;
	DrawPolylines();
}

// ------------------------------------------------------------------------------
double TDrawShape::SGN(double x) {
	if (x > 0)
		return 1;
	else {
		if (x < 0)
			return-1;
	}
	return 0;
}

// ------------------------------------------------------------------------------
void TDrawShape::ShowNormProbDistr(double A, double Sig2) {
	int i, N = 1000;
	double xmin = A - 5.0 * sqrt(Sig2);
	double xmax = A + 5.0 * sqrt(Sig2);
	double dx = (xmax - xmin) / ((double)N);
	double x = xmin;

	x = xmin;
	TURPolyLine pll(1, N); // создание простой пустой полининии
	for (i = 0; i < N; i++) {
		pll.Points[i].X = x;
		pll.Points[i].Y = exp(-((x - A) * (x - A) / 2.0) / Sig2) / sqrt
				(2.0 * M_PI * Sig2);
		x += dx;
	}

	nPll = 1;
	Pll[0] = pll;
	DrawPolylines();
}

// ------------------------------------------------------------------------------
// картинка с ФАР
void TDrawShape::PictFar() {
	int i;
	// первые 4 точки - излучатели
	// arrPoints[4] - начало координат - фазовый центр
	TURPointXY arrPts[5];
	Pts[0] = arrPts[0] = TURPointXY(0.0, -3.0 / 2.0);
	Pts[1] = arrPts[1] = TURPointXY(0.0, -1.0 / 2.0);
	Pts[2] = arrPts[2] = TURPointXY(0.0, 1.0 / 2.0);
	Pts[3] = arrPts[3] = TURPointXY(0.0, 3.0 / 2.0);
	Pts[4] = arrPts[4] = TURPointXY(100.0, 0.0);
	nPts = 5;

	// угол цели
	double valAl = 10.0 / 180.0 * M_PI;
	// линия фазового фронта проходящая через фазовый центр
	TURPointXY pntPhaseCntr(0.0, 0.0);
	double valxfront = 6.0;

	TURPointXY pntFront0(pntPhaseCntr.X - valxfront * sin(valAl),
		pntPhaseCntr.Y + valxfront * cos(valAl));
	TURPointXY pntFront1(-pntFront0.X, -pntFront0.Y);
	TURPolyLine pllFront(pntFront0, pntFront1);

	// точки пересечения лучей с линией фазового фронта
	TURPointXY arrPFront[5];
	for (i = 0; i < 4; i++) {
		arrPFront[i] = TURPointXY(pntPhaseCntr.X - arrPts[i].Y * cos(valAl) * sin
			(valAl), pntPhaseCntr.X + arrPts[i].Y * cos(valAl) * cos(valAl));
	}
	arrPFront[4] = pntPhaseCntr;

	// линии с сигналом от цели проходят через точки arrPoints
	TURPolyLine arrPll[5], arrSegmSdvig[5];
	double valx = 10.0;
	for (i = 0; i < 5; i++) {
		TURPointXY pnt2(valx, valx * tan(valAl) + arrPts[i].Y);
		arrPll[i] = TURPolyLine(arrPFront[i], pnt2);
		arrSegmSdvig[i] = TURPolyLine(arrPts[i], arrPFront[i]);
	}
	Pll[0] = arrPll[0]; // массив рисуемых полиний
	Pll[1] = arrPll[1];
	Pll[2] = arrPll[2];
	Pll[3] = arrPll[3];
	Pll[4] = arrPll[4];
	nPll = 5; // к-во рисуемых полиний

	// нормаль к ФАР в точке фазового центра
	TURPointXY pntNorm(10.0, 0.0);
	TURPolyLine pllNorm(pntNorm, pntPhaseCntr);

	// массив антенных модулей
	TURPolyLine arrPllAM[4];
	for (i = 0; i < 4; i++) {
		TURPointXY pntTopRight(0.0, arrPts[i].Y + 0.5);
		TURPointXY pntBottomLeft(-0.05, arrPts[i].Y - 0.5);
		arrPllAM[i] = TURPolyLine(pntTopRight, pntBottomLeft);
	}

	// нарисовать arrPts[5]
	DrawPoints();
	// нарисовать arrPll[5]
	DrawPolylines();
	// нарисовать pllNorm
	Pll[0] = pllNorm; // массив рисуемых полиний
	nPll = 1; // к-во рисуемых полиний
	DrawPolylines();

	TURPolyLine pllFar(arrPts[0], arrPts[4]);

	// нарисовать  arrPlgAM[4]
	Pll[0] = arrPllAM[0]; // массив рисуемых полиний
	Pll[1] = arrPllAM[1];
	Pll[2] = arrPllAM[2];
	Pll[3] = arrPllAM[3];
	nPll = 4; // к-во рисуемых полиний
	DrawPolylines();

	// обозначение угла между нормалью и лучом
	TURPointXY pnt0(0.0, 0.0);
	double Dist0 = 0.8, Dist1 = 0.85;
	CreateAngleMarks(pnt0, pntNorm, arrPll[4].Points[1], Dist0, Dist1);

	// обозначение угла между фронтом и ФАР
	CreateAngleMarks(pnt0, arrPts[0], pntFront1, Dist0, Dist1);
	CreateAngleMarks(pnt0, arrPts[3], pntFront0, Dist0, Dist1);

	// линия фронта волны
	Pll[0] = pllFront; // массив рисуемых полиний
	nPll = 1; // к-во рисуемых полиний
	DrawPolylines();

	// отрезки фазоовых сдвигов
	Pll[0] = arrSegmSdvig[0];
	Pll[1] = arrSegmSdvig[1];
	Pll[2] = arrSegmSdvig[2];
	Pll[3] = arrSegmSdvig[3];
	Pll[4] = arrSegmSdvig[4];
	nPll = 5; // к-во рисуемых полиний
	DrawPolylines();

	// фазовый центр
	Pts[0] = pntPhaseCntr;
	nPts = 1;
	DrawPoints();

}
// ------------------------------------------------------------------------------
TDrawShape TDrawShape:: operator = (TDrawShape dsh) { // TDrawShape tmp;
	int i;

	for (i = 0; i < dsh.nPts; i++)
		Pts[i] = dsh.Pts[i];
	for (i = 0; i < dsh.nPll; i++)
		Pll[i] = dsh.Pll[i];
	for (i = 0; i < dsh.nPlg; i++)
		Plg[i] = dsh.Plg[i];

	for (i = 0; i < dsh.nPts; i++)
		ColPts[i] = dsh.ColPts[i];
	for (i = 0; i < dsh.nPll; i++)
		ColPll[i] = dsh.ColPll[i];
	for (i = 0; i < dsh.nPlg; i++)
		ColPlg[i] = dsh.ColPlg[i];

	for (i = 0; i < dsh.nPts; i++)
		RPts[i] = dsh.RPts[i];
	for (i = 0; i < dsh.nPll; i++)
		RPll[i] = dsh.RPll[i];
	for (i = 0; i < dsh.nPlg; i++)
		RPlg[i] = dsh.RPlg[i];

	nPts = dsh.nPts; // к-во точек
	nPll = dsh.nPll; // к-во полиний
	nPlg = dsh.nPlg; // к-во полигонов
	// члены визуализации
	Holst = dsh.Holst; // орган управления TImage
	Wframe = dsh.Wframe; // ширина рамки
	dkScale = dsh.dkScale; // коэфф.масштабирования
	dkx = dsh.dkx; // коэфф.масштабирования
	dky = dsh.dky; // коэфф.масштабирования
	// X0= dsh.X0;  // начало координат
	// Y0= dsh.Y0;  // начало координат

	IS_AXES_NEED = true;

	Box[0] = dsh.Box[0];
	Box[1] = dsh.Box[1];
	Box[2] = dsh.Box[2];
	Box[3] = dsh.Box[3];

	return *this;
}
// ------------------------------------------------------------------------------
TDrawShape TDrawShape:: operator = (TDrawShape & dsh) { // TDrawShape tmp;
	int i;

	for (i = 0; i < dsh.nPts; i++)
		Pts[i] = dsh.Pts[i];
	for (i = 0; i < dsh.nPll; i++)
		Pll[i] = dsh.Pll[i];
	for (i = 0; i < dsh.nPlg; i++)
		Plg[i] = dsh.Plg[i];

	for (i = 0; i < dsh.nPts; i++)
		ColPts[i] = dsh.ColPts[i];
	for (i = 0; i < dsh.nPll; i++)
		ColPll[i] = dsh.ColPll[i];
	for (i = 0; i < dsh.nPlg; i++)
		ColPlg[i] = dsh.ColPlg[i];

	for (i = 0; i < dsh.nPts; i++)
		RPts[i] = dsh.RPts[i];
	for (i = 0; i < dsh.nPll; i++)
		RPll[i] = dsh.RPll[i];
	for (i = 0; i < dsh.nPlg; i++)
		RPlg[i] = dsh.RPlg[i];

	nPts = dsh.nPts; // к-во точек
	nPll = dsh.nPll; // к-во полиний
	nPlg = dsh.nPlg; // к-во полигонов
	// члены визуализации
	Holst = dsh.Holst; // орган управления TImage
	Wframe = dsh.Wframe; // ширина рамки
	dkScale = dsh.dkScale; // коэфф.масштабирования
	dkx = dsh.dkx; // коэфф.масштабирования
	dky = dsh.dky; // коэфф.масштабирования
	// X0= dsh.X0;  // начало координат
	// Y0= dsh.Y0;  // начало координат

	IS_AXES_NEED = true;

	Box[0] = dsh.Box[0];
	Box[1] = dsh.Box[1];
	Box[2] = dsh.Box[2];
	Box[3] = dsh.Box[3];

	return *this;
}
// ------------------------------------------------------------------------------
TDrawShape TDrawShape:: operator+(TURPolyLine pll) { // TDrawShape tmp;
	if (nPll == MAX_NPLL) { /* tmp= *this; */ return(*this);
	}
	else { // tmp= *this;
		Pll[nPll] = pll;
		nPll++;
		return(*this);
	}
}
// ------------------------------------------------------------------------------
TDrawShape TDrawShape:: operator += (TURPointXY pt) {
	if (nPts <= MAX_NPTS) {
		Pts[nPts] = pt;
		nPts++;
	}
	return(*this);
}
// ------------------------------------------------------------------------------
TDrawShape TDrawShape:: operator += (TURPolyLine pll) {
	if (nPll <= MAX_NPLL) {
		Pll[nPll] = pll;
		nPll++;
	}
	return(*this);
}
// ------------------------------------------------------------------------------
TDrawShape TDrawShape:: operator += (TURPolygon plg) {
	if (nPlg <= MAX_NPLG) {
		Plg[nPlg] = plg;
		nPlg++;
	}
	return(*this);
}

// ------------------------------------------------------------------------------
void TDrawShape::ClearPts(void) {
	nPts = 0; // к-во точек
}

// ------------------------------------------------------------------------------
void TDrawShape::ClearPll(void) {
	nPll = 0; // к-во полиномов
}

// ------------------------------------------------------------------------------
void TDrawShape::ClearPlg(void) {
	nPlg = 0; // к-во полигонов
}

// ------------------------------------------------------------------------------
void TDrawShape::ClearAll(void) {
	nPts = 0; // к-во точек
	nPll = 0; // к-во полиномов
	nPlg = 0; // к-во полигонов
}

// ------------------------------------------------------------------------------
void TDrawShape::Save(AnsiString fn) {
	int i, k, n, sz, m;
	char* cp;
	double x, y;
	FILE* pf;

	pf = fopen(fn.c_str(), "w+b");

	cp = (char*)(&nPts);
	sz = sizeof(int);
	for (i = 0; i < sz; i++)
		fputc((int)cp[i], pf);

	sz = sizeof(double);
	for (m = 0; m < nPts; m++) {
		cp = (char*)(&Pts[m].X);
		for (i = 0; i < sz; i++)
			fputc((int)cp[i], pf);
		cp = (char*)(&Pts[m].Y);
		for (i = 0; i < sz; i++)
			fputc((int)cp[i], pf);
	}

	// СОХРАНЕНИЕ ПОЛИНИЙ
	cp = (char*)(&nPll);
	sz = sizeof(int);
	for (i = 0; i < sz; i++)
		fputc((int)cp[i], pf);

	for (n = 0; n < nPll; n++) { // Сохранить Pll[n].NumParts
		cp = (char*)(&Pll[n].NumParts);
		sz = sizeof(int);
		for (i = 0; i < sz; i++)
			fputc((int)cp[i], pf);
		// Сохранить Pll[n].Parts[]
		for (i = 0; i < Pll[n].NumParts; i++) {
			cp = (char*)(&Pll[n].Parts[i]);
			sz = sizeof(int);
			for (i = 0; i < sz; i++)
				fputc((int)cp[i], pf);
		}
		// Сохранить Pll[n].NumPoints
		cp = (char*)(&Pll[n].NumPoints);
		sz = sizeof(int);
		for (i = 0; i < sz; i++)
			fputc((int)cp[i], pf);

		sz = sizeof(double);
		for (i = 0; i < Pll[n].NumParts - 1; i++) {
			for (k = Pll[n].Parts[i]; k < Pll[n].Parts[i + 1]; k++) {
				x = Pll[n].Points[k].X;
				y = Pll[n].Points[k].Y;
				cp = (char*)(&x);
				for (m = 0; m < sz; m++)
					fputc((int)cp[m], pf);
				cp = (char*)(&y);
				for (m = 0; m < sz; m++)
					fputc((int)cp[m], pf);
			}
		}
		for (k = Pll[n].Parts[Pll[n].NumParts - 1]; k < Pll[n].NumPoints; k++) {
			x = Pll[n].Points[k].X;
			y = Pll[n].Points[k].Y;
			cp = (char*)(&x);
			for (m = 0; m < sz; m++)
				fputc((int)cp[m], pf);
			cp = (char*)(&y);
			for (m = 0; m < sz; m++)
				fputc((int)cp[m], pf);
		}
	}

	cp = (char*)(&nPlg);
	sz = sizeof(int);
	for (i = 0; i < sz; i++)
		fputc((int)cp[i], pf);

	for (n = 0; n < nPlg; n++) { // Сохранить Plg[n].NumParts
		cp = (char*)(&Pll[n].NumParts);
		sz = sizeof(int);
		for (i = 0; i < sz; i++)
			fputc((int)cp[i], pf);
		// Сохранить Plg[n].Parts[]
		for (i = 0; i < Pll[n].NumParts; i++) {
			cp = (char*)(&Pll[n].Parts[i]);
			sz = sizeof(int);
			for (i = 0; i < sz; i++)
				fputc((int)cp[i], pf);
		}
		// Сохранить Plg[n].NumPoints
		cp = (char*)(&Pll[n].NumPoints);
		sz = sizeof(int);
		for (i = 0; i < sz; i++)
			fputc((int)cp[i], pf);

		sz = sizeof(double);
		for (i = 0; i < Plg[n].NumParts - 1; i++) {
			for (k = Plg[n].Parts[i]; k < Plg[n].Parts[i + 1]; k++) {
				x = Plg[n].Points[k].X;
				y = Plg[n].Points[k].Y;
				cp = (char*)(&x);
				for (m = 0; m < sz; m++)
					fputc((int)cp[m], pf);
				cp = (char*)(&y);
				for (m = 0; m < sz; m++)
					fputc((int)cp[m], pf);
			}
		}
		for (k = Plg[n].Parts[Plg[n].NumParts - 1]; k < Plg[n].NumPoints; k++) {
			x = Plg[n].Points[k].X;
			y = Plg[n].Points[k].Y;
			cp = (char*)(&x);
			for (m = 0; m < sz; m++)
				fputc((int)cp[m], pf);
			cp = (char*)(&y);
			for (m = 0; m < sz; m++)
				fputc((int)cp[m], pf);
		}
	}

	cp = (char*)(&bPts[0]);
	sz = sizeof(bPts);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);
	cp = (char*)(&bPll[0]);
	sz = sizeof(bPll);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);
	cp = (char*)(&bPlg[0]);
	sz = sizeof(bPlg);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);

	cp = (char*)(&ColPts[0]);
	sz = sizeof(ColPts);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);
	cp = (char*)(&ColPll[0]);
	sz = sizeof(ColPll);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);
	cp = (char*)(&ColPlg[0]);
	sz = sizeof(ColPlg);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);

	cp = (char*)(&RPts[0]);
	sz = sizeof(RPts);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);
	cp = (char*)(&RPll[0]);
	sz = sizeof(RPll);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);
	cp = (char*)(&RPlg[0]);
	sz = sizeof(RPlg);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);

	cp = (char*)(&Wframe);
	sz = sizeof(int);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);
	// cp= (char*)(&X0); sz= sizeof(int);
	// for(m=0;m<sz;m++) fputc((int)cp[m],pf);
	// cp= (char*)(&Y0); sz= sizeof(int);
	// for(m=0;m<sz;m++) fputc((int)cp[m],pf);

	cp = (char*)(&dkScale);
	sz = sizeof(double);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);
	cp = (char*)(&dkx);
	sz = sizeof(double);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);
	cp = (char*)(&dky);
	sz = sizeof(double);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);

	cp = (char*)(&Box[0]);
	sz = sizeof(Box);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);

	cp = (char*)(&IS_AXES_NEED);
	sz = sizeof(bool);
	for (m = 0; m < sz; m++)
		fputc((int)cp[m], pf);

	goto end;

end: ;
	fclose(pf); ;
}

// ------------------------------------------------------------------------------
void TDrawShape::Open(AnsiString fn) {
	int i, k, n, sz, m, size, pos;
	char* cp;
	double x, y;
	FILE* pf;
	TURPolyLine pll;

	int nparts, parts[100], npts;
	TURPointXY pts[5000];

	char* buf;

	pf = fopen(fn.c_str(), "rb");
	fseek(pf, 0, SEEK_END);
	size = ftell(pf);
	fseek(pf, 0, SEEK_SET);
	buf = new char[size];
	for (i = 0; i < size; i++)
		buf[i] = (char)fgetc(pf);
	fclose(pf);

	pos = 0;
	// СЧИТАТЬ ТОЧКИ Pts
	cp = (char*)(&nPts);
	sz = sizeof(int);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;

	sz = sizeof(double);
	for (m = 0; m < nPts; m++) {
		cp = (char*)(&Pts[m].X);
		for (i = 0; i < sz; i++)
			cp[i] = buf[pos + i];
		pos += sz;
		cp = (char*)(&Pts[m].Y);
		for (i = 0; i < sz; i++)
			cp[i] = buf[pos + i];
		pos += sz;
	}

	// СЧИТАТЬ ПОЛИЛИНИИ Pll
	cp = (char*)(&nPll);
	sz = sizeof(int);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;

	for (n = 0; n < nPll; n++) { // Считать Pll[n].NumParts
		cp = (char*)(&nparts);
		sz = sizeof(int);
		for (i = 0; i < sz; i++)
			cp[i] = buf[pos + i];
		pos += sz;
		// Считать Pll[n].Parts[]
		for (i = 0; i < nparts; i++) {
			cp = (char*)(&parts[i]);
			sz = sizeof(int);
			for (i = 0; i < sz; i++)
				cp[i] = buf[pos + i];
			pos += sz;
		}
		// Считать Pll[n].NumPoints
		cp = (char*)(&npts);
		sz = sizeof(int);
		for (i = 0; i < sz; i++)
			cp[i] = buf[pos + i];
		pos += sz;

		sz = sizeof(double);
		// Считывание точек полилиний
		for (i = 0; i < nparts - 1; i++) {
			for (k = parts[i]; k < parts[i + 1]; k++) {
				cp = (char*)(&x);
				for (i = 0; i < sz; i++)
					cp[i] = buf[pos + i];
				pos += sz;
				cp = (char*)(&y);
				for (i = 0; i < sz; i++)
					cp[i] = buf[pos + i];
				pos += sz;
				pts[k].X = x;
				pts[k].Y = y;
			}
		}
		for (k = parts[nparts - 1]; k < npts; k++) {
			cp = (char*)(&x);
			for (i = 0; i < sz; i++)
				cp[i] = buf[pos + i];
			pos += sz;
			cp = (char*)(&y);
			for (i = 0; i < sz; i++)
				cp[i] = buf[pos + i];
			pos += sz;
			pts[k].X = x;
			pts[k].Y = y;
		}
		// инициализируем Pll[n]
		TURPolyLine pll(nparts, npts, parts, pts);
		Pll[n] = pll;

	}

	// СЧИТАТЬ ПОЛИГОНЫ Plg
	cp = (char*)(&nPlg);
	sz = sizeof(int);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;

	for (n = 0; n < nPlg; n++) { // Считать Pll[n].NumParts
		cp = (char*)(&nparts);
		sz = sizeof(int);
		for (i = 0; i < sz; i++)
			cp[i] = buf[pos + i];
		pos += sz;
		// Считать Pll[n].Parts[]
		for (i = 0; i < nparts; i++) {
			cp = (char*)(&parts[i]);
			sz = sizeof(int);
			for (i = 0; i < sz; i++)
				cp[i] = buf[pos + i];
			pos += sz;
		}
		// Считать Pll[n].NumPoints
		cp = (char*)(&npts);
		sz = sizeof(int);
		for (i = 0; i < sz; i++)
			cp[i] = buf[pos + i];
		pos += sz;

		sz = sizeof(double);
		// Считывание точек полигонов
		for (i = 0; i < nparts - 1; i++) {
			for (k = parts[i]; k < parts[i + 1]; k++) {
				cp = (char*)(&x);
				for (i = 0; i < sz; i++)
					cp[i] = buf[pos + i];
				pos += sz;
				cp = (char*)(&y);
				for (i = 0; i < sz; i++)
					cp[i] = buf[pos + i];
				pos += sz;
				pts[k].X = x;
				pts[k].Y = y;
			}
		}
		for (k = parts[nparts - 1]; k < npts; k++) {
			cp = (char*)(&x);
			for (i = 0; i < sz; i++)
				cp[i] = buf[pos + i];
			pos += sz;
			cp = (char*)(&y);
			for (i = 0; i < sz; i++)
				cp[i] = buf[pos + i];
			pos += sz;
			pts[k].X = x;
			pts[k].Y = y;
		}
		// инициализируем Pll[n]
		TURPolygon plg(nparts, npts, parts, pts);
		Plg[n] = plg;
	}

	cp = (char*)(&bPts[0]);
	sz = sizeof(bPts);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;
	cp = (char*)(&bPll[0]);
	sz = sizeof(bPll);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;
	cp = (char*)(&bPlg[0]);
	sz = sizeof(bPlg);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;

	cp = (char*)(&ColPts[0]);
	sz = sizeof(ColPts);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;
	cp = (char*)(&ColPll[0]);
	sz = sizeof(ColPll);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;
	cp = (char*)(&ColPlg[0]);
	sz = sizeof(ColPlg);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;

	cp = (char*)(&RPts[0]);
	sz = sizeof(RPts);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;
	cp = (char*)(&RPll[0]);
	sz = sizeof(RPll);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;
	cp = (char*)(&RPlg[0]);
	sz = sizeof(RPlg);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;

	cp = (char*)(&Wframe);
	sz = sizeof(int);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;

	// cp= (char*)(&X0); sz= sizeof(int);
	// for(i=0;i<sz;i++) cp[i]= buf[pos+i]; pos+=sz;
	// cp= (char*)(&Y0); sz= sizeof(int);
	// for(i=0;i<sz;i++) cp[i]= buf[pos+i]; pos+=sz;

	cp = (char*)(&dkScale);
	sz = sizeof(double);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;
	cp = (char*)(&dkx);
	sz = sizeof(double);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;
	cp = (char*)(&dky);
	sz = sizeof(double);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;

	cp = (char*)(&Box[0]);
	sz = sizeof(Box);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;

	cp = (char*)(&IS_AXES_NEED);
	sz = sizeof(bool);
	for (i = 0; i < sz; i++)
		cp[i] = buf[pos + i];
	pos += sz;

	goto end;
end: ;
	delete[]buf; ;
}
// ------------------------------------------------------------------------------
