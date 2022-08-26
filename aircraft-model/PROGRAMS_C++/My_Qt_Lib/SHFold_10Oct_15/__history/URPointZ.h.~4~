//---------------------------------------------------------------------------

#ifndef URPointZH
#define URPointZH
#include "URFigure.h"






//---------------------------------------------------------------------------
class  TYrRastr;
class TURPointXY ;
class TURPointZ : public TURFigure
{


public:
	// virtual ~ TURPointZ() ;
	double X  ; // X coordinate
	double Y  ;// Y coordinate
	double Z  ;// Z coordinate
	double M  ; // Measure
	 TURPointZ() ;
	 TURPointZ(const double x1,const double y1,const double z1) ;
	 TURPointZ(const double x1,const double y1,  TYrRastr *Rastr,const int i);
	 TURPointZ(const  TURPointXY P, TYrRastr *Rastr) ;
	 TURPointZ(const TURPointZ&p); // констр копир
	 TURPointZ operator = (TURPointZ p2) ; // оператор присваивания

	 static double Dist(const TURPointZ p1,const TURPointZ p2) ;
	 static double ScalMult(const TURPointZ p1,const TURPointZ p2);
	 static double Norm(const TURPointZ p1);
	 static TURPointZ Plus(const TURPointZ p1,const TURPointZ p2);
	 static TURPointZ Minus(const TURPointZ p1,const TURPointZ p2); // p1-p2
	 static TURPointZ MultReal(const TURPointZ p1,const double alf) ;
	 void showMe() ;
	 static TURPointZ ParamPoint(const TURPointZ p1,const TURPointZ p2,const double alf);
	 static double CosA(const TURPointZ p1,const TURPointZ p2) ;
	 static TURPointZ VectorProduct(const TURPointZ p1,const TURPointZ p2); // вект произв
	 static void DoNorm(TURPointZ &P);
// изменение порядка следования байтов в 4 байтном слове (32 бита массив)
//  input : chstr - указатель на массив char[4]
// output: chstr - указатель на массив  char[4] c измененным порядком следования
// байтов  char[0] = char[3] ; char[1] = char[2] ;  char[2] = char[1] ; char[3] = char[0] ;
	static void ChangeByteOrder(int * pi0) ;


	static void WriteDBASEFile(wchar_t *wchFileName,TURPointZ *purPnt, const int quantPoint) ;

	static void WriteMainFile(wchar_t *wchFileName,TURPointZ *purPnt, const int quantPoint) ;
	static void WriteIndexFile(wchar_t *wchFileName,TURPointZ *purPnt, const int quantPoint) ;
	static void WriteSetSHPFiles(wchar_t *wchFileName,TURPointZ  *purPnt, const int quantPoint) ;

	static void ReadSHPFile(wchar_t *wchFileName,TURPointZ  **ppurPnt,  int *pquantPoint) ;

	static int findIntertsectLine_And_Sphere(const TURPointZ PNtPos,const TURPointZ PNtVelo, const double VAlSphereRad
	, double *pT, TURPointZ *pPnt);

	virtual void   LinearTransformation(const double  valAng , const TURPointXY pntSdvig,const double valRastigenie );
};

#endif
