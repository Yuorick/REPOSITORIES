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

	TURPointXY  Pts[MAX_NPTS];    // ������ �������� �����
	TURPolyLine Pll[MAX_NPLL];    // ������ �������� �������
	TURPolygon  Plg[MAX_NPLG];    // ������ �������� ���������
	int nPts;         // �-�� �����
	int nPll;         // �-�� �������
	int nPlg;         // �-�� ���������

	double Box[4];      // �������������� �����-����  {xmin,ymin,xmax,ymax}


	bool bPts[MAX_NPTS];    // �������� ����������� �� ������
	bool bPll[MAX_NPLL];    // �������� ����������� �� ������
	bool bPlg[MAX_NPLG];    // �������� ����������� �� ������

	TColor ColPts[MAX_NPTS];   	// ����� �����
	TColor ColPll[MAX_NPLL];   	// ����� ���������
	TColor ColPlg[MAX_NPLG];   	// ����� ���������

	int RPts[MAX_NPTS];        	// ������� �����
	int RPll[MAX_NPLL];        	// ������� �������
	int RPlg[MAX_NPLG];        	// ������� ���������


	// ����� ������������
	TImage* Holst; 			// ����� ���������� TImage
	int Wframe;         // ������ �����

	double dkScale;     // �����.���������������
	double dkx,dky;     // ���.�����.���������-� �������� �� X � Y



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


	void Save(AnsiString fn);     // ��������� ������ � ���� ����� draw
	void Open(AnsiString fn);     // ������� ������ �� ����� draw


void __fastcall TDrawShape::DrawTimeFun
														(	double *pF,   //  ������ �������� �-�� ������� - ��� �������
															double *pT,   //  ������ ������� ������� - ��� �������
															int N,   			//  ������ ��������
															double &dT,   //  ��������������� �� ��������
															double &dF    //  ��������������� �� ��������
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
															(	double *mtxData, 	// ������� ������ nRows x nCols
																const int nCols, 	// - �-�� �������
																const int nRows, 	// - �-�� �������� �������
																wchar_t *Names, 		// ������ � �������� ���������� - ������� nCols X lenName
																int lenName, // ����.����� ����� ����������
																int nX,  	// � ������� ����������, ���-�� ���������
																int nY,  	// � ������� ����������, ���-�� ���������
																double dX,  	// ������� �� ��������
																double dY  	// ������� �� ��������
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
															TURPointXY XY0    // ������ ���������
														);
void TDrawShape::	 CreateAngleMarks
									 ( TURPointXY P0,  // ������� ����
										 TURPointXY P1,  // 1-� �����
										 TURPointXY P2,  // 2-� �����
										 double D0,      // ���������� �� 1-�� �������
										 double D1       // ���������� �� 2-�� �������
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
