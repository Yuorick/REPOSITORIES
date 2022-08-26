//---------------------------------------------------------------------------

#ifndef URPolygonZH
#define URPolygonZH

#include <stdio.h>
 class 	TURPointXY;
 class  TURPolygon;
//---------------------------------------------------------------------------

class TURPolygonZ
{
public:
   double Box[4] ;		// Bounding Box
	int NumParts ;
	int NumPoints ;
	int *Parts ;
	 TURPointXY *Points ;
	double ZRange[2] ;
	double *ZArray ;
	double MRange[2] ;
	double *MArray ;

	__fastcall ~TURPolygonZ() ;
	// �������� ������������
	 TURPolygonZ &operator=(const TURPolygonZ  &R)  ;
	 // ����������� �����������
	 TURPolygonZ () ;
	 TURPolygonZ (const TURPolygonZ &R) ;
	 TURPolygonZ(wchar_t*FileName);
	 TURPolygonZ( const int iNumParts, const int iNumPoints);
	 TURPolygonZ(const TURPolygon &R) ;


//
	int WriteToASCII(wchar_t*FileName);
	static void fillPartsOfPolygonFromASCII(wchar_t*FileName, const int quanParts
			,const int quanPoints,int *piarrParts) ;
	static void calcPointsOfPolygonFromASCII(wchar_t*FileName,int *quanParts
			,int *quanPoints);
	void calcZRange();
	void calcMRange();
	void calcBoundBox();
	double calcPartSq(const int n) ;
	 double calcVectSq()  ;
	 double calcLeng();
	 double calcPartLeng(const int n);
	 static double dist(TURPointXY*p0, TURPointXY*p1) ;
	bool ReadPolygonZFromASCII(wchar_t*FileName);
// ��������� ������� ���������� ������ � 4 ������� ����� (32 ���� ������)
//  input : chstr - ��������� �� ������ char[4]
// output: chstr - ��������� �� ������  char[4] c ���������� �������� ����������
// ������  char[0] = char[3] ; char[1] = char[2] ;  char[2] = char[1] ; char[3] = char[0] ;
static void ChangeByteOrder(int * pi0) ;


static void WriteDBASEFile(wchar_t *wchFileName,TURPolygonZ *purPlg, const int quantPlg) ;

static void WriteMainFile(wchar_t *wchFileName,TURPolygonZ *purPlg, const int quantPlg) ;
static void WriteIndexFile(wchar_t *wchFileName,TURPolygonZ *purPlg, const int quantPlg) ;
static void WriteSetSHPFiles(wchar_t *wchFileName,TURPolygon *purPlg, const int quantPlg) ;

static void ReadSHPFile(wchar_t *wchFileName,TURPolygonZ **ppurPlg,  int *pquantPlg)  ;
static void WriteSetSHPFiles(wchar_t *wchFileName,TURPolygonZ *purPlg, const int quantPlg) ;
};

#endif




