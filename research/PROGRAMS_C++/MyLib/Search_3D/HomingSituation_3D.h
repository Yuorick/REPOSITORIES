//---------------------------------------------------------------------------

#ifndef HomingSituation_3DH
#define HomingSituation_3DH
#include "MainForm.h"
#include "Bomb.h"
#include "ShipTarg.h"
#include "BombTraj_3D.h"
#include "Filt.h"
#include "CntlFuncPar.h"
class  TBombTraj_3D;
class  TShipTarg;
class  TCntlFuncPar;
class    TURPointXY ;
class TURPolygon ;


class THomingSituation_3D {
public:

	// траектория АУБ
	TBombTraj_3D mBombTraj_3D;
	TCntlFuncPar mCntlFuncPar ;

	long double mTCur;

	// путь у папке с отчетом
	wchar_t *mpwcharrFoldReport;

	~THomingSituation_3D();
	// конструктор по умолчанию
	THomingSituation_3D();
	// конструктор копирования
	THomingSituation_3D(const THomingSituation_3D &R);

	// оператор присваивания
	THomingSituation_3D operator = (THomingSituation_3D R2);

	// парам конструктор
	THomingSituation_3D( TBombTraj_3D BombTraj_3D //
		, TCntlFuncPar CntlFuncPar
		,long double TCur
		, wchar_t *pwcharrFoldReport // путь у папке с отчетом
		);


	void fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs( wchar_t *wcharrPath, long double &valDHoriz );

	void fncCalcU( long double TCur,long double &valUGor, long double &valUVert);

	THomingSituation_3D fncFindOptimalControl_for_Dist_MethodPerebora_1Point(TCntlFuncPar &CntlFuncPar
	   ,const long double valShagT, long double& valDGorizOpt) ;

	static int fncSign (long double a);

	void fncMoveClas_TO_BlastPoint_AND_ShowGraphs( wchar_t *wcharrPath, long double &valDHoriz );
	void fncFindSetPointsOfApproachibility(const int nCircleCoef,const long double valShagT
		, TURPointXY **pparrPnt, int &lenArr);

	TURPolygon fncFindSetOfDopPoints(wchar_t *pwchOutFile, const int QuantVar, const long double valStepT );





};
#endif
