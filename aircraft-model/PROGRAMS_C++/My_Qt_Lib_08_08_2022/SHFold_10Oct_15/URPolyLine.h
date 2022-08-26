﻿//---------------------------------------------------------------------------

#ifndef URPolyLineH
#define URPolyLineH


#include "URFigure.h"
 class TURPointXY;
 class TURPolyLineZ ;
 class TArcEllipse ;
 class TArcParab ;
 class TSector ;
 class TCircle ;
 class THyperbola;
 class QLineF;
//---------------------------------------------------------------------------

class TURPolyLine : public TURFigure
{
public:
	double Box[4] ;		// Bounding Box
	int NumParts ;
	int NumPoints ;
	int *Parts ;
	TURPointXY *Points ;

        ~TURPolyLine() ;

  TURPolyLine();
 // конструктор копирования
 TURPolyLine (const TURPolyLine &R) ;
 // оператор присваивания
 TURPolyLine &operator=(const TURPolyLine  &R);
 // параметр констр создания из  TURPolyLineZ
 TURPolyLine (const TURPolyLineZ &R);
// парам констр
TURPolyLine(double *parrx,double *parry,const int quanPoints);

 // парам констр
 TURPolyLine(wchar_t*FileName) ;
// парам констр
 TURPolyLine(TURPointXY * pPoints,const int qPoints) ;
 // парам констр
 TURPolyLine ( const int iNumParts, const int iNumPoints) ;
 // парам констр
 TURPolyLine(const int iNumParts,const int iNumPoints,int *iarrParts
			,TURPointXY *arrPoints) ;
 // парам констр
 TURPolyLine( TURPolyLine *parrPln, const int lenarrPln);
  // парам констр
 TURPolyLine(  const int iNumParts,const int iNumPoints,int *iarrParts
                       ,double *arr);
 // параметр конструкто создания из математической линии эллипса
 TURPolyLine(  const TArcEllipse arcEll,const int iNumPoints)  ;
 // параметр конструкто создания из математической линии параболы
TURPolyLine (  const TArcParab arcPrbl,const int iNumPoints) ;
// параметр конструкто создания из математической линии сектора окружности
TURPolyLine(  const TSector Sect,const int iNumPoints);
// параметр конструкто создания из математической линии окружности
TURPolyLine(  const TCircle Circle,const int iNumPoints) ;

// параметр конструкто создания 2 точек
TURPolyLine( const TURPointXY  pnt1, const TURPointXY  pnt2) ;

TURPolyLine(const THyperbola hyperb,const int iNumPoints0, const double VAlDiap);

TURPolyLine( QLineF &qline);



//
int WriteToASCII(wchar_t*FileName) ;
bool ReadFromASCII(wchar_t*FileName);
virtual void calcBoundBox();
double calcLeng();
double calcPartLeng(const int n);
static double dist(TURPointXY*p0, TURPointXY*p1) ;

// изменение порядка следования байтов в 4 байтном слове (32 бита массив)
//  input : chstr - указатель на массив char[4]
// output: chstr - указатель на массив  char[4] c измененным порядком следования
// байтов  char[0] = char[3] ; char[1] = char[2] ;  char[2] = char[1] ; char[3] = char[0] ;
static void ChangeByteOrder(int * pi0) ;


static void WriteDBASEFile(wchar_t *wchFileName,TURPolyLine *purPlg, const int quantPlg) ;

static void WriteMainFile(wchar_t *wchFileName,TURPolyLine *purPlg, const int quantPlg) ;
static void WriteIndexFile(wchar_t *wchFileName,TURPolyLine *purPlg, const int quantPlg) ;
static void WriteSetSHPFiles(wchar_t *wchFileName,TURPolyLine *purPlg, const int quantPlg) ;

static void ReadSHPFile(wchar_t *wchFileName,TURPolyLine **ppurPlg,  int *pquantPlg)  ;

TURPolyLine   fncLinTransform(double * arrMtxPer );

static TURPolyLine fncCreateSector(const TURPointXY pointCentre, const double valR,
					const double valFi0,const double valFi1,const int NPoints) ;

static TURPolyLine fncCreateArrow(const TURPointXY pointBegin, const TURPointXY pointEnd
					,const double valLength,const double valAng);

static TURPolyLine fncCreateAxes(const TURPointXY pointBeginX, const TURPointXY pointEndX
									   ,const TURPointXY pointBeginY, const TURPointXY pointEndY
									  ,const double valLength);

 TURPolyLine   RastiagenieTransform( const double valRastigenie );

 TURPolyLine SdvigTransform(const TURPointXY pntSdvig ) ;

 TURPolyLine LinTransform(const double  valAng , const TURPointXY pntSdvig,const double valRastigenie ) ;

 static double Sign(const double x);

 static TURPolyLine   createFishNet(const double ValLeft, const double ValRight
   ,const double ValDown,const double ValUp, const double ValCellSizeWidth, const double ValCellSizeHeight );

static void   createTangGraph(wchar_t*FileName, const int NUmPoints );

static void createOddPolinomGraph(wchar_t*wchrPolinomPoliLine, int iNPolinom, double *arrCoeff, double valDiapazon )  ;

static void   createOddPolinomMinusZGraph(wchar_t* FileName, int iNPolinom, double *arrCoeff, double valDiapazon );

static void   createPolinomGraph(wchar_t* FileName, int iNPolinom, double *arrCoeff, double valDiapazon0, double valDiapazon1 );

static  TURPolyLine create_SinX_Div_X_Line(TURPointXY pntSdvig, double scalex, double scaley, int numRoot);

bool calcWayPoint( double VAlDist, TURPointXY &pntRez);

double LinearValueApprox(double x);

TURPolyLine nomalizationDiagr(double valMax, double *pXMax) ;

void stretchDiagrAlongXY(double valCoeffX, double valCoeffY);

TURPolyLine RepaireArg();

double calcDiagrWidth();

int calcGraphArgMax();

double LinearDerivApprox(double x);

static TURPolyLine createLineDiagram(
		 double *parrAimingPoints_X , int *piarrRepeatQuants_X
		 , const int iQuantAimingPoints_X) ;

static TURPolyLine createLineDiagram(
		 double *parrAimingPoints_X , double *parrRepeatQuants_X
		 , const int iQuantAimingPoints_X);

TURPolyLine   MultScalar(const double VAla );

 // извлечение из составного   полилинии  простой с номером n
TURPolyLine extractSimplePart(const int n);

int calcQuantApproximatingPoints(const double valTargCellSize) ;

void fillTargPointsArray( const double valTargCellSize,  TURPointXY *pPntArr) ;

int  WriteToASCII__(wchar_t*FileName);

/*******************************************************************************************************************************/
virtual int  calcQuantsNetPoints(const double  VAlCellSize) ;

virtual void createTargPointsArray(const int valTargCellSize, TURPointXY **ppTargPntArray, int *lenTargPntArray);

virtual void WriteSetSHPFiles(wchar_t *wchFileName);

virtual void createUnatedPointsArray(TURPointXY **ppunatedPoints, int *quantUnitedPoints) ;

virtual void   LinearTransformation(const double  valAng , const TURPointXY pntSdvig,const double valRastigenie );

virtual void   ConvexShell(TURPolygon *pPolgConv);

static TURPolyLine createUZP_Skaliga_AirPlane(const double VAlR,const double VAllProbMean,  const int NUmPoints);

static double fncScaligaPlaneLaw__(const double VAlR1, const double VAlProbMean, const double VAlMiss);

static double fncDiagrSinx_div_x__(double tet);

TURPolyLine   LinTransform(const double  valAng, const TURPointXY pntCentre
  , const TURPointXY pntSdvig,const double valRastigenie );

TURPolyLine   SimOtragenieY_and_flip();

void flip();

void flipPart(const int n) ;

TURPolyLine   setUpGraph( const TURPointXY pntDate00
,const double  valAng ,const double valCoeffX,const double valCoeffY );



};

 double fncOddPolinom1(double *arrCoeff, int lenarr, double valArg);
 double fncOddPolinom1MinusZ(double *arrCoeff, int lenarr, double valArg);
 double fncPolinom1(double *arrCoeff, int iNPolinom, double valArg);

#endif
