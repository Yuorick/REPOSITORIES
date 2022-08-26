//---------------------------------------------------------------------------

#ifndef URPolyLineZH
#define URPolyLineZH

#include <stdio.h>

class TURPointXY;
class TURPolyLine ;
//---------------------------------------------------------------------------

class TURPolyLineZ
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
	  TURPolyLineZ () ;

	__fastcall ~TURPolyLineZ() ;
	 // конструктор копирования
	  TURPolyLineZ  (const TURPolyLineZ  &R) ;
	  // оператор присваивания
      TURPolyLineZ  &operator=(const TURPolyLineZ   &R2) ;
//  параметр констр
   TURPolyLineZ(wchar_t*FileName);
//  параметр констр1
	TURPolyLineZ (const TURPolyLine &R) ;
//  параметр констр 2
   TURPolyLineZ( const int iNumParts, const int iNumPoints);
  // парам констр
TURPolyLineZ( const int iNumParts, const int iNumPoints, int *iarrParts
	  ,TURPointXY *arrPoints, double *arrZ);



	 double calcLeng();
	 double calcPartLeng(const int n);
	 static double dist(TURPointXY*p0, TURPointXY*p1) ;


	int WriteToASCII(wchar_t*FileName);
	static void fillPartsOfPolyLineFromASCII(wchar_t*FileName, const int quanParts
					,const int quanPoints,int *piarrParts) ;
   static void calcPointsOfPolyLineFromASCII(wchar_t*FileName,int *quanParts
					,int *quanPoints);
	void calcZRange();
	void calcMRange();
	void calcBoundBox();
	bool ReadPolyLyneZFromASCII(wchar_t*FileName);
// изменение порядка следования байтов в 4 байтном слове (32 бита массив)
//  input : chstr - указатель на массив char[4]
// output: chstr - указатель на массив  char[4] c измененным порядком следования
// байтов  char[0] = char[3] ; char[1] = char[2] ;  char[2] = char[1] ; char[3] = char[0] ;
static void ChangeByteOrder(int * pi0) ;


static void WriteDBASEFile(wchar_t *wchFileName,TURPolyLineZ *purPlg, const int quantPlg) ;

static void WriteMainFile(wchar_t *wchFileName,TURPolyLineZ *purPlg, const int quantPlg) ;
static void WriteIndexFile(wchar_t *wchFileName,TURPolyLineZ *purPlg, const int quantPlg) ;
static void WriteSetSHPFiles(wchar_t *wchFileName,TURPolyLineZ *purPlg, const int quantPlg) ;

static void ReadSHPFile(wchar_t *wchFileName,TURPolyLineZ **ppurPlg,  int *pquantPlg)  ;
};

#endif
