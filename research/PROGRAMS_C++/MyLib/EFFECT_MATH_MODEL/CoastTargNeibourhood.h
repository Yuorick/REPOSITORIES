//---------------------------------------------------------------------------

#ifndef CoastTargNeibourhoodH
#define CoastTargNeibourhoodH
#include "URMultiPoint.h"
#include "YrRastr.h"
#include "URPolygon.h"

class TURPointXY;
class TURFigure;
class TURPolygon;
class TURMultiPoint;
class TYrRastr;
class TCoastTargNeighbourhood
{
 public:
 double marrMtrxCorrFluct [4]; // флуктуац ошибки
 double marrMtrxCorrSyst [4]; // ситемат ошибки
 double mKillingRange;
 int mQuantShells ;
 TURMultiPoint mMultPntAim;
 TURMultiPoint mMultPntTarg;
 TYrRastr mRstrBiasAims;
 TYrRastr mRstrProbBias;
 TURPolygon mPlgAim;
 //TYrRastr mRstrAimPoints;

	__fastcall ~TCoastTargNeighbourhood();

	TCoastTargNeighbourhood();

	TCoastTargNeighbourhood(const TCoastTargNeighbourhood &R);

	TCoastTargNeighbourhood &operator=(const TCoastTargNeighbourhood  &R);

	TCoastTargNeighbourhood( double *arrMtrxCorr,
 const double VAlKillingRange, const int QUuantShells ,  TURPointXY *pPntArrAim,
 TURPointXY *pPntArrTarg, const int QUuantPointsTarg, const  int QUantPointsAims);

 TCoastTargNeighbourhood( double *arrMtrxCorr,double *arrMtrxCorrSyst,
 const double VAlKillingRange, const int QUuantShells ,  TURPointXY *pPntArrAim,
 TURPointXY *pPntArrTarg, const int QUuantPointsTarg, const  int QUantPointsAims);

 TCoastTargNeighbourhood(wchar_t *pwchSHP_PointsTargFile,
 const double VATargCellSize ,const double VAAimCellSize , double *arrMtrxCorr,
 const double VAlKillingRange, const int QUuantShells );

 TCoastTargNeighbourhood(TURFigure  *pTargFigure,
 const double VAlTargCellSize ,const double VAlAimCellSize , double *arrMtrxCorr,
 const double VAlKillingRange, const int QUuantShells );
/*
 TCoastTargNeighbourhood(TURFigure  *pTargFigure,
 const double VAlTargCellSize , TURPointXY *pPntArrAim,const  int QUantPointsAims, double *arrMtrxCorr,
 const double VAlKillingRange, const int QUuantShells );
*/
 TCoastTargNeighbourhood(TURFigure  *pTargFigure,
 const double VAlTargCellSize ,const double VAlAimCellSize , double *arrMtrxCorrFluct
 , double *arrMtrxCorrSyst, const double VAlKillingRange, const int QUuantShells , const double VAlCoeff );

TCoastTargNeighbourhood(TURFigure  *pTargFigure,
 const double VAlTargCellSize , TURPointXY *pPntArrAim,const  int QUantPointsAims
 , double *arrMtrxCorrFluct, double *arrMtrxCorrSyst,
 const double VAlKillingRange, const int QUuantShells );

TCoastTargNeighbourhood(TURFigure  *pTargFigure,
 const double VAlTargCellSize , const double VAlProbSystCellSize, TURPointXY *pPntArrAim,const  int QUantPointsAims
 , double *arrMtrxCorrFluct, double *arrMtrxCorrSyst,
 const double VAlKillingRange, const int QUuantShells );



 static void calcOptimalDestructionProb_For_GroupLinedPoints (const int  QAantIspit, double *parrMtrxCorr
	,const double VAlKillingRange, const int QUantShells,const int  In, double *parrTargX, const int LEnarrX
	, double *parrAimingPoints, int *piarrRepeatQuants, int *pQuantAimingPoints, TURPointXY *ppntArrAimingPoints_SKT
	, double *pvalProb, TURPointXY *pPntArrDensity);

 static bool calcOptimalArray_Of_AimPoints_ForGroupLinedPoints(double *arrElK
	,const double VAlKillingRange, const int QUantShells,const int  In, double *parrZ, const int LEnarrZ
	, double *parrAimingPoints, int *piarrRepeatQuants, int *pQuantAimingPoints, double &valObj, TURPointXY *pPntArrDensity);

bool calcOptimalArray_Of_AimPoints( TURPointXY *pPntArrAimingPoints
		, int *pQuantAimingPoints, int *piarrRepeatQuants, double *pvalObj) ;


 double calcMinObjectFunc_AlongWithLine(double *parr_L,double * parrX, double * parrXZv, double *parrRez) ;

 void calc_d2FGr_po_dx(double *parr_L,double * parrX, double *arr_d2FGr_po_dx);

 void  calcGradFGr(double *parr_L,double * parrX,double * arrGrad);

 double calcFi_i(double *parr_L,double * parrX, const int NUmi);

double calcFGr(double *parr_L,double * parrX);


static double _MAX_double( double x, double y)
{
	return (x>y)?x:y;
}

static int _MIN_int( int x, int y)
{
	return (x<y)?x:y;
}

static void calcAimingPoints_For_OpenManPower_MUS (double *arrElK,const double VAlKillingRange
		, const int QUantShells, TURPolygon plgTarg, double valDistAppPoint
		,TURPointXY *ppntArrAimingPoints, int *piarrRepeatAimingPoints
		, int* piQuantAimingPoints);

static void calcAimingPoints_For_OpenManPower_MUS_old(double *arrElK,const double VAlKillingRange
		, const int QUantShells, TURPolygon plgTarg, double valDistAppPoint
		,TURPointXY *ppntArrAimingPoints, int *piarrRepeatAimingPoints
		, int* piQuantAimingPoints);

void createMatrxL(const double VAlStepIntegr, double *parr_L);

double calcEfficiencyOfStrategy(double *parrStrategy) ;

double findNextAimingPoint(double *parr_L, double *parrStrategy, int *pinumNext);

double estimateStrategy_MonteCarlo(const int QUantIspit, int *piarrStrategy) ;

void applyKillingRange(const TURPointXY  pntFall, bool* barrHit);

bool calcOptimalArray_Of_AimPoints_With_SystMtrx( TURPointXY *pPntArrAimingPoints
		, int *pQuantAimingPoints, int *piarrRepeatQuants, double *pvalObj) ;

void createMatrxL_With_SystMtrx(const double VAlStepIntegr, double *parr_L);

void calcGradFGr_With_SystMtrx(double *parr_L,double * parrX,double * arrGrad);

double calcFi_jk(double *parr_L, double *parrX, const int j,const int  k) ;

double calcFGr_With_SystMtrx(double *parr_L  ,double * parrX) ;

double calcMinObjectFunc_AlongWithLine_With_SystMtrx(double *parr_L,double * parrX, double * parrXZv, double *parrRez);

void createMtrxExpSum(double *parr_L, double *parr_U, double *arrExpSum) ;

double findNextAimingPoint_With_SystMtrx(double *parr_L, double *parrStrategy, int *pinumNext);

double calcEfficiencyOfStrategy_With_SystMtrx( double *parrStrategy);

double calcEfficiencyOfStrategy_With_SystMtrx_Var1 (double *parrStrategy) ;



};
#endif
