﻿//---------------------------------------------------------------------------

#ifndef UrPointXYH
#define UrPointXYH
#include "URFigure.h"


class TURPointXY : public TURFigure
{
public:

	double X ;
	double Y ;

	 TURPointXY() ;


	TURPointXY(const double x1,const double y1);
	 // конструктор копирования
		TURPointXY (const TURPointXY &R) ;

	  // оператор присваивания
        TURPointXY &operator=(const TURPointXY  &R2) ;

		TURPointXY operator+(TURPointXY  p2) ;

		TURPointXY operator-(TURPointXY  p2) ;

	  double operator*(TURPointXY  p2)  ;

		void ShowMe()	;
		static double calcVectS(TURPointXY P1,TURPointXY P2);

		static double dist(TURPointXY P1,TURPointXY P2);

		static TURPointXY  MultReal(const TURPointXY P, const double alf);

	  static TURPointXY ParamPoint(const TURPointXY P1,const TURPointXY P2,const double alf);
		double Norm();

		static void DoNorm(TURPointXY &P);

		static double  ScalMult(const TURPointXY P1,const TURPointXY P2);

	  static int PutPointsToCsvFile(wchar_t*FileName,TURPointXY *urpntP,
							double* pZ,const int lenArray,int * lenVars) ;

	  static int PutPointsToTxtFile(wchar_t*FileName,TURPointXY *urpntP,
							 double* pZ,const int lenArray,int * lenVars) ;

      static int PutPointsToCsvFile_(wchar_t*FileName,TURPointXY *urpntP,
                                                 const int lenArray);

// изменение порядка следования байтов в 4 байтном слове (32 бита массив)
//  input : chstr - указатель на массив char[4]
// output: chstr - указатель на массив  char[4] c измененным порядком следования
// байтов  char[0] = char[3] ; char[1] = char[2] ;  char[2] = char[1] ; char[3] = char[0] ;
static void ChangeByteOrder(int * pi0) ;


static void WriteDBASEFile(wchar_t *wchFileName,TURPointXY *purPnt, const int quantPoint) ;

static void WriteMainFile(wchar_t *wchFileName,TURPointXY *purPnt, const int quantPoint) ;

static void WriteIndexFile(wchar_t *wchFileName,TURPointXY *purPnt, const int quantPoint) ;

static void WriteSetSHPFiles(wchar_t *wchFileName,TURPointXY  *purPnt, const int quantPoint) ;

static void ReadSHPFile(wchar_t *wchFileName,TURPointXY  **ppurPnt,  int *pquantPoint) ;

static bool fncIsPointsEqual(const TURPointXY p0, const TURPointXY p1)  ;

TURPointXY    fncLinTrasform(double * arrMtxPer);

static double  calcAng( TURPointXY P1, TURPointXY P2);

static  double  calcCosAng( TURPointXY P1, TURPointXY P2) ;

static void fncCleanSetPoints(const double  valEps, TURPointXY **pparrPnt, int &lenArr);

TURPointXY   LinTransform(const double  valAng , const TURPointXY pntSdvig,const double valRastigenie );

TURPointXY SdvigTransform(const TURPointXY pntSdvig );

TURPointXY  LinTransform( const TURPointXY pntSdvig,double * arrMtxPer );

static void subtractEqualPoints(TURPointXY *parrPnts0 // массив точек, input
			, const int lenarr0 //  длина массива точек , input
			, const double VAlTolerance // точность
			, int *plenarr // к-во отличающихся точек
				);

virtual void createTargPointsArray( const int valTargCellSize, TURPointXY **ppTargPntArray, int *lenTargPntArray);

static double  calcDiam(TURPointXY *pPntArr, const int lenarr, int *pnum0,  int *pnum1);

virtual void   LinearTransformation(const double  valAng , const TURPointXY pntSdvig,const double valRastigenie );

virtual void WriteSetSHPFiles(wchar_t *wchFileName);

void sortPointsArray_X (TURPointXY *pPntarr, const int lenarr);

static void calcBoundBox__(TURPointXY *arrPoints, const int lenarr, double *arrBox);

static void swap(TURPointXY &pnt0, TURPointXY &pnt1);

};

int fncCmp( const void *a, const void *b);



//---------------------------------------------------------------------------
#endif
