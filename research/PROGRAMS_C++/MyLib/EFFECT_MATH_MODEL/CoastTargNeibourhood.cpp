//---------------------------------------------------------------------------


#pragma hdrstop

#include "CoastTargNeibourhood.h"
#include <string.h>
#include <math.h>
#include <float.h>
#include "MatrixProccess.h"
#include "Detonator.h"
#include "Gauss.h"
#include "CalcCorMatrx.h"
#include "YrRastr.h"
#include "Gauss.h"
#include "URPointXY.h"
#include "ProbabilityTheory.h"
#include "GameTheory.h"
#include "Fight.h"
#include "LinOptimization.h"
#include "URFigure.h"
#include "URMultiPoint.h"
#include "URPolygon.h"
extern const double NODATA;
extern const double NODATA1;
//---------------------------------------------------------------------------

__fastcall TCoastTargNeighbourhood::~TCoastTargNeighbourhood()
{
 //	if(mMultPntAim.Points) delete [] mMultPntAim.Points ;
 //	mMultPntAim.Points = NULL ;
 //	if(mMultPntTarg.Points) delete []mMultPntTarg.Points ;
 //	mMultPntTarg.Points = NULL ;
}
//---------------------------------------------------------------------------


	TCoastTargNeighbourhood ::TCoastTargNeighbourhood()
{
 memset(marrMtrxCorrFluct , 0, 4 * sizeof(double));
 memset(marrMtrxCorrSyst, 0, 4 * sizeof(double));
 mKillingRange = 0.;
 mQuantShells = 0 ;

}

// ����������� �����������
 TCoastTargNeighbourhood ::TCoastTargNeighbourhood(const TCoastTargNeighbourhood &R)
 {
	mKillingRange  = R.mKillingRange ;
	mQuantShells = R.mQuantShells;
	mMultPntAim = R.mMultPntAim;
	mMultPntTarg = R.mMultPntTarg;
	mRstrBiasAims = R.mRstrBiasAims;
	mRstrProbBias = R.mRstrProbBias;
	mPlgAim = R.mPlgAim;
 }
 // �������� ������������
 TCoastTargNeighbourhood &TCoastTargNeighbourhood::operator=(const TCoastTargNeighbourhood  &R)
 {
	mKillingRange  = R.mKillingRange ;
	mQuantShells = R.mQuantShells;
	mMultPntAim = R.mMultPntAim;
	mMultPntTarg = R.mMultPntTarg;
	mRstrBiasAims = R.mRstrBiasAims;
	mRstrProbBias = R.mRstrProbBias;
	mPlgAim = R.mPlgAim;

	return *this ;
 }
//----------------------------------------------------------------------------------------------
TCoastTargNeighbourhood::TCoastTargNeighbourhood( double *arrMtrxCorr,
 const double VAlKillingRange, const int QUuantShells ,  TURPointXY *pPntArrAim,
 TURPointXY *pPntArrTarg, const int QUuantPointsTarg, const  int QUantPointsAims)
{
	mKillingRange  = VAlKillingRange ;
	mQuantShells = QUuantShells;
	memcpy(marrMtrxCorrFluct, arrMtrxCorr, 4 * sizeof(double));
	memset(marrMtrxCorrSyst, 0, 4 * sizeof(double));
	mMultPntAim = TURMultiPoint( pPntArrAim,  QUantPointsAims);
	mMultPntTarg = TURMultiPoint( pPntArrTarg, QUuantPointsTarg);

}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
TCoastTargNeighbourhood::TCoastTargNeighbourhood( double *arrMtrxCorr,double *arrMtrxCorrSyst,
 const double VAlKillingRange, const int QUuantShells ,  TURPointXY *pPntArrAim,
 TURPointXY *pPntArrTarg, const int QUuantPointsTarg, const  int QUantPointsAims)
{
	mKillingRange  = VAlKillingRange ;
	mQuantShells = QUuantShells;

	memcpy(marrMtrxCorrFluct, arrMtrxCorr, 4 * sizeof(double));
	memcpy(marrMtrxCorrSyst, arrMtrxCorrSyst, 4 * sizeof(double));

	mMultPntAim = TURMultiPoint( pPntArrAim,  QUantPointsAims);
	mMultPntTarg = TURMultiPoint( pPntArrTarg, QUuantPointsTarg);

}

//----------------------------------------------------------------------------------------------
TCoastTargNeighbourhood::TCoastTargNeighbourhood(wchar_t *pwchSHP_PointsTargFile,
 const double VAlTargCellSize ,const double VAlAimCellSize , double *arrMtrxCorr,
 const double VAlKillingRange, const int QUuantShells )
 {
		mKillingRange  = VAlKillingRange ;
		mQuantShells = QUuantShells;
		memcpy(marrMtrxCorrFluct, arrMtrxCorr, 4 * sizeof(double));
		memset(marrMtrxCorrSyst, 0, 4 * sizeof(double));

		// ����������� ���� ������� SHP � �-�� ��������
		TURFigure  Figure, *pFigure;
		int iShapeType = -1  // ���
		, iObjQuant = -1;       // �-��
		Figure.find_Objects_Type_And_Quant(pwchSHP_PointsTargFile, &iShapeType, &iObjQuant)  ;

		int lenObjArr =1;
		TURPointXY *pPnt0 =  (TURPointXY *)malloc(sizeof( TURPointXY)* lenObjArr);
		TURPointXY **ppPnt0 = &pPnt0;

		TURPolyLine *pPln0 =  (TURPolyLine *)malloc(sizeof( TURPolyLine)* lenObjArr);
		TURPolyLine **ppPln0 = &pPln0;

		TURPolygon *pPlg0 =  (TURPolygon *)malloc(sizeof( TURPolygon)* lenObjArr);
		TURPolygon **ppPlg0 = &pPlg0;

		TURMultiPoint MultiPoint00, *pMultiPoint00;
		///



		TURPolygon Polygon00, *pPolygon00;
		TURPolyLine PolyLine00, *pPolyLine00;
		TURPointXY  PointXY00(1.,1.), *pPointXY00;
		switch(iShapeType)
		{
			case 1:
			TURPointXY::ReadSHPFile(pwchSHP_PointsTargFile,ppPnt0,  &iObjQuant) ;
			MultiPoint00 =  TURMultiPoint(*ppPnt0, iObjQuant);
			pMultiPoint00 = &MultiPoint00;
			pFigure = (TURFigure*)(pMultiPoint00);
			break;

			case 3:
			TURPolyLine::ReadSHPFile(pwchSHP_PointsTargFile,ppPln0,  &iObjQuant) ;
			PolyLine00 =  TURPolyLine( *ppPln0, iObjQuant);
			pPolyLine00 = &PolyLine00;
			pFigure = (TURFigure*) (pPolyLine00);

			break;

			case 5:
			TURPolygon::ReadSHPFile(pwchSHP_PointsTargFile,ppPlg0,  &iObjQuant) ;
			Polygon00 = TURPolygon ( *ppPlg0, iObjQuant);
			pPolygon00 = &Polygon00 ;
			pFigure = (TURFigure*)(pPolygon00);

			break;

			default:
			ShowMessage(L"Error in Data SHP File");
			free(pPnt0);
			free(pPln0);
			free(pPlg0);
			Abort();

		}
		free(pPnt0);
		free(pPln0);
		free(pPlg0);


		int quantPointsTarg =1;
		TURPointXY *pTargPntArray =  (TURPointXY *)malloc(sizeof( TURPointXY)* quantPointsTarg);
		TURPointXY **ppTargPntArray = &pTargPntArray;

		pFigure->createTargPointsArray(VAlTargCellSize, ppTargPntArray, &quantPointsTarg);

		int quantPointsAims = 1.;
		TURPointXY *pAimPntArray =  (TURPointXY *)malloc(sizeof( TURPointXY)*  quantPointsAims);
		TURPointXY **ppAimPntArray = &pAimPntArray;
		pFigure->createAimPointsArray(VAlAimCellSize, ppAimPntArray, &quantPointsAims);

		mMultPntAim = TURMultiPoint( *ppAimPntArray,  quantPointsAims);
		mMultPntTarg = TURMultiPoint( *ppTargPntArray, quantPointsTarg);

		free(pAimPntArray);
		free(pTargPntArray);

 }

//----------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------
TCoastTargNeighbourhood::TCoastTargNeighbourhood(TURFigure  *pTargFigure,
 const double VAlTargCellSize ,const double VAlAimCellSize , double *arrMtrxCorr,
 const double VAlKillingRange, const int QUuantShells )
 {
		mKillingRange  = VAlKillingRange ;
		mQuantShells = QUuantShells;
		memcpy(marrMtrxCorrFluct, arrMtrxCorr, 4 * sizeof(double));
		memset(marrMtrxCorrSyst, 0, 4 * sizeof(double));
		int quantPointsTarg =1;
		int quantPointsAims = 1.;


		TURPointXY *pTargPntArray =  (TURPointXY *)malloc(sizeof( TURPointXY)* quantPointsTarg);
		TURPointXY **ppTargPntArray = &pTargPntArray;

		pTargFigure->createTargPointsArray(VAlTargCellSize, ppTargPntArray, &quantPointsTarg);


		TURPointXY *pAimPntArray =  (TURPointXY *)malloc(sizeof( TURPointXY)*  quantPointsAims);
		TURPointXY **ppAimPntArray = &pAimPntArray;
		pTargFigure->createAimPointsArray(VAlAimCellSize, ppAimPntArray, &quantPointsAims);

		mMultPntAim = TURMultiPoint( *ppAimPntArray,  quantPointsAims);
		mMultPntTarg = TURMultiPoint( *ppTargPntArray , quantPointsTarg);

		free(pAimPntArray);
		free(pTargPntArray);
 }

//----------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------
TCoastTargNeighbourhood::TCoastTargNeighbourhood(TURFigure  *pTargFigure,
 const double VAlTargCellSize ,const double VAlAimCellSize , double *arrMtrxCorrFluct
 , double *arrMtrxCorrSyst, const double VAlKillingRange, const int QUuantShells, const double VAlCoeff )
 {
		mKillingRange  = VAlKillingRange ;
		mQuantShells = QUuantShells;
		memcpy(marrMtrxCorrFluct, arrMtrxCorrFluct, 4 * sizeof(double));
		memcpy(marrMtrxCorrSyst, arrMtrxCorrSyst, 4 * sizeof(double));
		int quantPointsTarg =1;
		int quantPointsAims = 1.;

		TURPointXY *pTargPntArray =  (TURPointXY *)malloc(sizeof( TURPointXY)* quantPointsTarg);
		TURPointXY **ppTargPntArray = &pTargPntArray;

		pTargFigure->createTargPointsArray(VAlTargCellSize, ppTargPntArray, &quantPointsTarg);

		mPlgAim  = pTargFigure->createBuffPolygonInAccordanceWithCorrMtrx(arrMtrxCorrSyst, VAlCoeff);

	  //	plgAim.WriteSetSHPFiles(L"E:\\Ametist\\28-03-2018\\COAST\\Optimal_Strat_Rez\\plgAim.shp", &plgAim, 1);
		TURPolygon *pplgAim = &mPlgAim ;
		///
		TURFigure  Figure, *pFigure;
		pFigure = (TURFigure*)(pplgAim);

		TURPointXY *pAimPntArray =  (TURPointXY *)malloc(sizeof( TURPointXY)*  quantPointsAims);
		TURPointXY **ppAimPntArray = &pAimPntArray;
		pFigure->createAimPointsArray(VAlAimCellSize, ppAimPntArray, &quantPointsAims);

		mMultPntAim = TURMultiPoint( *ppAimPntArray,  quantPointsAims);
		mMultPntTarg = TURMultiPoint( *ppTargPntArray , quantPointsTarg);
		//mRstrAimPoints = TYrRastr( plgAim, VAlAimCellSize,  NODATA, mMultPntAim.Points[0]);

		free(pAimPntArray);
		free(pTargPntArray);
			///


  // �������� ������ ����������� ������� ������
				// ������� ������� ���� ������ �� ������ 3
				double arrMtrxInv[4] ={0.};
				InverseMtrx2(marrMtrxCorrSyst, arrMtrxInv);
				const double VAlDet = arrMtrxInv[3] * arrMtrxInv[0] - arrMtrxInv[1] * arrMtrxInv[1];
				double Box[4] = {0.};
				Box[0] = -sqrt(arrMtrxInv[3]/ VAlDet) * 3.;
				Box[1] = -sqrt(arrMtrxInv[0]/ VAlDet) * 3.;
				Box[2] = -Box[0];
				Box[3] = -Box[1];


		 TURPolygon Rect = TURPolygon::createRect(Box)  ;
		  TURPointXY pntZero(0.,0.);
		 mRstrProbBias = TYrRastr   ( Rect, VAlAimCellSize,  1., pntZero) ;



		 // ���������� ������� ������

		 const double valSqrtDet = sqrt(marrMtrxCorrSyst[3] * marrMtrxCorrSyst[0] - marrMtrxCorrSyst[1] * marrMtrxCorrSyst[1]);
		 double sum = 0.;
		 for (int i =0; i < mRstrProbBias.ncols * mRstrProbBias.nrows; i++)
		 {
			 TURPointXY pntTemp =  mRstrProbBias.getCellCentre(i) ;
			 double arrx[2] = {0.};
			 arrx[0] = pntTemp.X;
			 arrx[1] = pntTemp.Y;


			 double temp = calcYT_D_Y(arrx, arrMtrxInv, 2 );
			// if (temp >2.)
			// {
			//	mRstrProbBias.pflt_rastr[i] =  NODATA;
			//	continue;
			// }
			 mRstrProbBias.pflt_rastr[i] = 1./ (2. * M_PI)/ valSqrtDet* exp(-temp /2.) * VAlAimCellSize * VAlAimCellSize;
			 sum+=  mRstrProbBias.pflt_rastr[i];
		 }
		 for (int i =0; i < mRstrProbBias.ncols * mRstrProbBias.nrows; i++)
		 {
			 mRstrProbBias.pflt_rastr[i] = mRstrProbBias.pflt_rastr[i] / sum;
		 }
	 ///

		mPlgAim.calcBoundBox() ;
		TURPolygon plgBiasAim = mPlgAim.Buffer(Box[2],Box[3]);
	   //	TURPolygon plgBiasAim =  pTargFigure->createBuffPolygonInAccordanceWithCorrMtrx(marrMtrxCorrSyst, 3.);
		mRstrBiasAims = 	 TYrRastr ( plgBiasAim, VAlAimCellSize,  1., mMultPntAim.Points[0]);

 }

//----------------------------------------------------------------------------------------------

/*
//----------------------------------------------------------------------------------------------
TCoastTargNeighbourhood::TCoastTargNeighbourhood(TURFigure  *pTargFigure,
 const double VAlTargCellSize , TURPointXY *pPntArrAim,const  int QUantPointsAims, double *arrMtrxCorr,
 const double VAlKillingRange, const int QUuantShells )
{
	mKillingRange  = VAlKillingRange ;
	mQuantShells = QUuantShells;
	memcpy(marrMtrxCorrFluct, arrMtrxCorr, 4 * sizeof(double));
	int quantPointsTarg =1;

	TURPointXY *pTargPntArray =  (TURPointXY *)malloc(sizeof( TURPointXY)* quantPointsTarg);
	TURPointXY **ppTargPntArray = &pTargPntArray;

	pTargFigure->createTargPointsArray(VAlTargCellSize, ppTargPntArray, &quantPointsTarg);


	mMultPntAim = TURMultiPoint( pPntArrAim,  QUantPointsAims);
	mMultPntTarg = TURMultiPoint( *ppTargPntArray , quantPointsTarg);

	free(pTargPntArray);
}
   */
//----------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------
TCoastTargNeighbourhood::TCoastTargNeighbourhood(TURFigure  *pTargFigure,
 const double VAlTargCellSize , const double VAlProbSystCellSize, TURPointXY *pPntArrAim,const  int QUantPointsAims
 , double *arrMtrxCorrFluct, double *arrMtrxCorrSyst,
 const double VAlKillingRange, const int QUuantShells )
{
	mKillingRange  = VAlKillingRange ;
	mQuantShells = QUuantShells;
	memcpy(marrMtrxCorrFluct, arrMtrxCorrFluct, 4 * sizeof(double));
	memcpy(marrMtrxCorrSyst, arrMtrxCorrSyst, 4 * sizeof(double));
	int quantPointsTarg =1;

	TURPointXY *pTargPntArray =  (TURPointXY *)malloc(sizeof( TURPointXY)* quantPointsTarg);
	TURPointXY **ppTargPntArray = &pTargPntArray;

	pTargFigure->createTargPointsArray(VAlTargCellSize, ppTargPntArray, &quantPointsTarg);


	mMultPntAim = TURMultiPoint( pPntArrAim,  QUantPointsAims);
	mMultPntTarg = TURMultiPoint( *ppTargPntArray , quantPointsTarg);


  // �������� ������ ����������� ������� ������
				// ������� ������� ���� ������ �� ������ 3
				double arrMtrxInv[4] ={0.};
				InverseMtrx2(marrMtrxCorrSyst, arrMtrxInv);
				const double VAlDet = arrMtrxInv[3] * arrMtrxInv[0] - arrMtrxInv[1] * arrMtrxInv[1];
				double Box[4] = {0.};
				Box[0] = -sqrt(arrMtrxInv[3]/ VAlDet) * 3.;
				Box[1] = -sqrt(arrMtrxInv[0]/ VAlDet) * 3.;
				Box[2] = -Box[0];
				Box[3] = -Box[1];

		 TURPolygon Rect = TURPolygon::createRect(Box)  ;
		  TURPointXY pntZero(0.,0.);
		 mRstrProbBias = TYrRastr   ( Rect, VAlProbSystCellSize,  1., pntZero) ;

		 const double valSqrtDet = sqrt(marrMtrxCorrSyst[3] * marrMtrxCorrSyst[0] - marrMtrxCorrSyst[1] * marrMtrxCorrSyst[1]);
		 double sum = 0.;
		 for (int i =0; i < mRstrProbBias.ncols * mRstrProbBias.nrows; i++)
		 {
			 TURPointXY pntTemp =  mRstrProbBias.getCellCentre(i) ;
			 double arrx[2] = {0.};
			 arrx[0] = pntTemp.X;
			 arrx[1] = pntTemp.Y;


			 double temp = calcYT_D_Y(arrx, arrMtrxInv, 2 );

			 mRstrProbBias.pflt_rastr[i] = 1./ (2. * M_PI)/ valSqrtDet* exp(-temp /2.) * VAlTargCellSize * VAlTargCellSize;
			 sum+=  mRstrProbBias.pflt_rastr[i];
		 }
		 for (int i =0; i < mRstrProbBias.ncols * mRstrProbBias.nrows; i++)
		 {

			 mRstrProbBias.pflt_rastr[i] = mRstrProbBias.pflt_rastr[i] / sum;
		 }
	 ///

	free(pTargPntArray);
}

//----------------------------------------------------------------------------------------------



// ���������� ������� ����������� ����� ������������  ��� ��������� �������� ����

//OUTPUT:
//  parrAimingPoints -  ������ ��������� ����� ������������
//  piarrRepeatQuants - ������ ���������� ����� ������������
//  *pQuantAimingPoints - ����� �- �� ��������� ����� ������������
bool TCoastTargNeighbourhood::calcOptimalArray_Of_AimPoints( TURPointXY *pPntArrAimingPoints
		, int *pQuantAimingPoints, int *piarrRepeatQuants, double *pvalObj)

{
	int iNumArgMin = -1;

	// ������������ �������
	double *parr_L = new double  [ mMultPntTarg.NumPoints *  mMultPntAim.NumPoints];
	const double VAlStepIntegr = mKillingRange /10.;
	///
	createMatrxL(VAlStepIntegr, parr_L) ;

	///
	double *parrX = new double [mMultPntAim.NumPoints];
	memset(parrX, 0, sizeof(double) * mMultPntAim.NumPoints);

 for (int i=0; i < mMultPntAim.NumPoints; i++)
 {
	 parrX [i] = ((double)mQuantShells)/ ((double)mMultPntAim.NumPoints) ;
 }

	double *parrXZv = new double [mMultPntAim.NumPoints];
	memset(parrXZv, 0, sizeof(double) * mMultPntAim.NumPoints);

double valEps = 0.001 ;
	double valObj0 = calcFGr(parr_L, parrX);
	bool brez = false;
	for (int i = 0; i < 500; i++)
	{
	 // ����������   ��  ���
		 // ������� �������
			double *parr_f= new double [mMultPntAim.NumPoints];
			double *arrGrad= new double [mMultPntAim.NumPoints];
			double *arr_d2FGr_po_dx= new double [mMultPntAim.NumPoints * mMultPntAim.NumPoints];
			double *parrRez = new double [mMultPntAim.NumPoints];
			calcGradFGr(parr_L,  parrX, arrGrad);
			memcpy(parr_f, arrGrad,sizeof(double) * mMultPntAim.NumPoints);
			int numArgMin = -1;
			double valrez = MinDoubleArray(parr_f, mMultPntAim.NumPoints, &numArgMin) ;
			memset( parrXZv, 0, mMultPntAim.NumPoints * sizeof(double));
			parrXZv[numArgMin] =  (double)mQuantShells;
			///



			double *parrRez0 = new double [mMultPntAim.NumPoints];
			*pvalObj = calcMinObjectFunc_AlongWithLine(parr_L,  parrX, parrXZv, parrRez0);


		if (fabs((*pvalObj) - valObj0) < valEps)
			{
				brez = true;
				break;
			}
			memcpy(parrX,parrRez0, mMultPntAim.NumPoints * sizeof(double));
			valObj0 = *pvalObj;
			delete []parr_f;
			delete []arrGrad;
			delete []arr_d2FGr_po_dx;
			delete []parrRez;
			delete []parrRez0;

	}

	delete []parrXZv;

	int *piarr_q = new int [mMultPntAim.NumPoints];   // �-�� ���������� ����� �������� � ������� i
	memset( piarr_q , 0,mMultPntAim.NumPoints * sizeof(int));
 //	double *parr_z =  new double[mMultPntAim.NumPoints];
 //	memset(parr_z , 0,mMultPntAim.NumPoints * sizeof(double));
	double *parrStrategy =  new double [mMultPntAim.NumPoints];
	memset(parrStrategy , 0,mMultPntAim.NumPoints * sizeof(double));
	int iSum = 0;
	for (int i = 0; i < mMultPntAim.NumPoints; i++)
	{
		piarr_q[i] =  ((int)parrX[i]);
		iSum +=  piarr_q[i];
 //		parr_z [i] =   parrX[i] - ((double)piarr_q[i]);
		parrStrategy [i] = (double)piarr_q[i];
	}
	///

	int iJ = mQuantShells - iSum;
	int inumNext = -1;
	double valObjective = -1;

	for (int i = 0; i < iJ; i++)
	{

		*pvalObj  = findNextAimingPoint(parr_L,parrStrategy, &inumNext);
		parrStrategy[inumNext] += 1.;
	}
	///

	for (int i = 0; i < mMultPntAim.NumPoints; i++)
	{
		piarr_q[i] =  ((int)(parrStrategy[i] + 0.05));

	}

	*pQuantAimingPoints = 0;
	int icur = 0;
	for (int i = 0; i < mMultPntAim.NumPoints; i++)
	{
	  if (piarr_q[i] > 0)
	  {
		piarrRepeatQuants[icur] =  piarr_q[i];
		pPntArrAimingPoints [icur] = mMultPntAim.Points[i];
		icur++;

	  }
	}
	*pQuantAimingPoints = icur;
	delete []piarr_q;
	delete []parrX;
	delete []parrStrategy;
	delete []parr_L;
	return brez;
}


//-------------------------------------------------------------------------

void TCoastTargNeighbourhood::calcGradFGr(double *parr_L,double * parrX,double * arrGrad)
{
	memset(arrGrad, 0, mMultPntAim.NumPoints * sizeof(double));
	double *arrFi = new double [mMultPntTarg.NumPoints];
	for (int i =0; i < mMultPntTarg.NumPoints; i++)
	{
		 arrFi[i]  = calcFi_i(parr_L, parrX, i);

	}
	MtrxMultMatrx(arrFi,1, mMultPntTarg.NumPoints, parr_L,mMultPntAim.NumPoints, arrGrad) ;
	delete []arrFi;
}

void TCoastTargNeighbourhood::calc_d2FGr_po_dx(double *parr_L,double * parrX,double * arr_d2FGr_po_dx)
{
	memset(arr_d2FGr_po_dx, 0, mMultPntAim.NumPoints *mMultPntAim.NumPoints * sizeof(double));
	double *arrTemp = new double [mMultPntAim.NumPoints *mMultPntAim.NumPoints];
	double *arrTemp0 = new double [mMultPntAim.NumPoints *mMultPntAim.NumPoints];
	for (int i =0; i < mMultPntTarg.NumPoints; i++)
	{
		 double valTemp = calcFi_i(parr_L,  parrX, i);
		 MtrxMultMatrx(&parr_L [ i * mMultPntAim.NumPoints],mMultPntAim.NumPoints, 1, &parr_L [ i * mMultPntAim.NumPoints],mMultPntAim.NumPoints, arrTemp) ;
		 MatrxMultScalar(arrTemp, mMultPntAim.NumPoints, mMultPntAim.NumPoints, valTemp,arrTemp0);
		 MtrxSumMatrx(arr_d2FGr_po_dx, arrTemp0,mMultPntAim.NumPoints, mMultPntAim.NumPoints, arrTemp) ;
		 memcpy(arr_d2FGr_po_dx, arrTemp, mMultPntAim.NumPoints * mMultPntAim.NumPoints * sizeof(double));
	}

	delete []arrTemp;
	delete []arrTemp0;
}

double TCoastTargNeighbourhood::calcMinObjectFunc_AlongWithLine(double *parr_L,double * parrX, double * parrXZv, double *parrRez)
{
	double *arrTemp0 = new double [mMultPntAim.NumPoints];
	double *arrTemp1 = new double [mMultPntAim.NumPoints];
	double *arrTemp2 = new double [mMultPntAim.NumPoints];
	double valF =  calcFGr(parr_L,  parrX) ;
	memcpy(parrRez,  parrX, mMultPntAim.NumPoints * sizeof(double));
	int iC = 1000;
	for (int i = 1; i < iC; i++)
	{
	 double lamb = ((double)i) / ((double)iC -1.);
	 MatrxMultScalar(parrX, 1, mMultPntAim.NumPoints, (1. - lamb),arrTemp0);
	 MatrxMultScalar(parrXZv, 1, mMultPntAim.NumPoints,  lamb,arrTemp1);
	 MtrxSumMatrx(arrTemp0, arrTemp1 ,1, mMultPntAim.NumPoints, arrTemp2) ;
	 double valFTemp = calcFGr(parr_L ,arrTemp2);
		if(valFTemp >= valF)
		{
			break;
		}
		memcpy(parrRez, arrTemp2, mMultPntAim.NumPoints * sizeof(double));
		valF = valFTemp;
	}


	delete []arrTemp0;
	delete []arrTemp1;
	delete []arrTemp2;
	return valF;
}

double TCoastTargNeighbourhood::calcFi_i(double *parr_L,double * parrX, const int NUmi)
{
	double temp = ScalProduct(parrX , &parr_L[NUmi * mMultPntAim.NumPoints], mMultPntAim.NumPoints) ;

	return exp(temp);
}

double TCoastTargNeighbourhood::calcFGr(double *parr_L  ,double * parrX)
{
	double sum = 0.;
	for (int i =0; i < mMultPntTarg.NumPoints; i++)
	{
		 sum+= calcFi_i(parr_L,  parrX, i);
	}
	return sum;
}

//---------------------------------------------------------------------------

//-----------------------------------------------------------------------------------
// ���������� ������� ����� ������������  ��� ������� ��������
//INPUT:
// arrElK[4] - �������������� ������� ��������� ����� �������
// VAlKillingRange  - ������ ���������
// QUantShells  - �-�� ���������
// OUTPUT:
//  parrAimingPoints [QUantShells]-  ������ ��������� ����� ������������ �� ��� X
//  piarrRepeatQuants [QUantShells] - ������ ���������� ����� ������������ �� ��� X
//  *pQuantAimingPoints - �- �� ��������� ����� ������������ �� ��� X
//  parrAimingPoints [QUantShells]-  ������ ��������� ����� ������������ �� ��� Y
//  piarrRepeatQuants [QUantShells] - ������ ���������� ����� ������������ �� ��� Y
//  *pQuantAimingPoints - �- �� ��������� ����� ������������ �� ��� Y
//  mMultPntAim.NumPoints  - �-�� ������� ����� ������������
// Im  - �-�� ����������� ����� ���� �� ��� X (�-�� ����� �� �����)
void TCoastTargNeighbourhood::calcAimingPoints_For_OpenManPower_MUS_old(double *arrElK,const double VAlKillingRange
		, const int QUantShells, TURPolygon plgTarg, double valDistAppPoint
		,TURPointXY *ppntArrAimingPoints, int *piarrRepeatAimingPoints
		, int* piQuantAimingPoints)
{
	//if ( !( (mTarget.menumTargetType ==  OPEN_MANPOWER_LIE)
	//			||(mTarget.menumTargetType ==  OPEN_MANPOWER_STAND)
	//			||(mTarget.menumTargetType ==  BULLET_PROOF_LIE)
	 //			||(mTarget.menumTargetType ==  BULLET_PROOF_STAND)
		//		)
	 //	)
	 //	{
	 //		return;
	 //	}
	double *parrAimingPoints_X = new double  [QUantShells];
	double *parrAimingPoints_Y = new double  [QUantShells];
	int *piarrRepeatQuants_X = new int  [QUantShells];
	int *piarrRepeatQuants_Y = new int  [QUantShells];
	memset(parrAimingPoints_X, 0, QUantShells * sizeof(double));
	memset(piarrRepeatQuants_X, 0, QUantShells * sizeof(int));

	memset(parrAimingPoints_Y, 0, QUantShells * sizeof(double));
	memset(piarrRepeatQuants_Y, 0, QUantShells * sizeof(int));
	///
		plgTarg.calcBoundBox() ;

	// const double VAlDeep = plgTarg.Box[2] - plgTarg.Box[0] ;
	// const double VAlWidth = plgTarg.Box[3] - plgTarg.Box[1] ;
	 // �� ��� X
		double VAlDeep = plgTarg.Box[2] - plgTarg.Box[0];
		double valSig = sqrt(arrElK[0]);
	 double valPp = VAlDeep / valSig /2. + 2.;
	 int iTemp = valPp/ 2.;
	 int iPp =  (valPp/ 2. - ((double)iTemp)< 0.5)? 2 * iTemp + 1: 2 * iTemp + 3;
	 double valStepPp  =  VAlDeep / ((double)iPp - 1.);

	 for (int i =0; i < iPp; i++)
	 {
		parrAimingPoints_X [i] = plgTarg.Box[0] + ((double )i) *  valStepPp;

	 }

	 int iCur = 0;
	 for (int i =0; i < QUantShells; i++)
	 {
		piarrRepeatQuants_X[iCur]++;
		iCur++;
		if (iCur == iPp)
		{
			iCur = 0;
		}

	 }
	 *piQuantAimingPoints =  iPp;



	 ///

		// �� ��� Y
		const double VAlWidth =  plgTarg.Box[3] - plgTarg.Box[1];
		valSig = sqrt(arrElK[3]);
		double valC = (VAlWidth + 0.008 * valDistAppPoint)/ valSig/2. ;
	 iTemp = valC / 2.;
	 int iC =  (valC/ 2. - ((double)iTemp)< 0.5)? 2 * iTemp + 1: 2 * iTemp + 3;
	 double valStepC  =  VAlWidth / ((double)iC - 1.);

	 for (int i =0; i < iC; i++)
	 {
		parrAimingPoints_Y [i] = plgTarg.Box[1] + ((double )i) *  valStepC;
	 }

	 iCur = 0;
	 for (int i =0; i < QUantShells; i++)
	 {
		piarrRepeatQuants_Y[iCur]++;
		iCur++;
		if (iCur == iC)
		{
			iCur = 0;
		}

	 }
	// *piQuantAimingPoints_Y =  iC;
	 ///

	 // ������������ ������� ����� ������������
	 int nx = 0, ny = -1;
		for (int i =0; i < QUantShells; i++)
		{

		 nx = i % iPp;
		 if (nx == 0)
		 {
			 ny++;
		 }
		 if (ny == iC)
		 {
       ny = 0;
		 }
		 ppntArrAimingPoints [i] = TURPointXY (parrAimingPoints_X [nx], parrAimingPoints_Y [ny] );

		}
		*piQuantAimingPoints = QUantShells;
		for (int i =0; i < QUantShells; i++)
		{
		 piarrRepeatAimingPoints [i] = 1;
		 }

	delete []parrAimingPoints_X ;
	delete []parrAimingPoints_Y ;
	delete []piarrRepeatQuants_X;
	delete []piarrRepeatQuants_Y ;
}


//-----------------------------------------------------------------------------------
// ���������� ������� ����� ������������  ��� ������� ��������
//INPUT:
// arrElK[4] - �������������� ������� ��������� ����� �������
// VAlKillingRange  - ������ ���������
// QUantShells  - �-�� ���������
// OUTPUT:
//  parrAimingPoints [QUantShells]-  ������ ��������� ����� ������������ �� ��� X
//  piarrRepeatQuants [QUantShells] - ������ ���������� ����� ������������ �� ��� X
//  *pQuantAimingPoints - �- �� ��������� ����� ������������ �� ��� X
//  parrAimingPoints [QUantShells]-  ������ ��������� ����� ������������ �� ��� Y
//  piarrRepeatQuants [QUantShells] - ������ ���������� ����� ������������ �� ��� Y
//  *pQuantAimingPoints - �- �� ��������� ����� ������������ �� ��� Y
//  mMultPntAim.NumPoints  - �-�� ������� ����� ������������
// Im  - �-�� ����������� ����� ���� �� ��� X (�-�� ����� �� �����)
void TCoastTargNeighbourhood::calcAimingPoints_For_OpenManPower_MUS(double *arrElK,const double VAlKillingRange
		, const int QUantShells, TURPolygon plgTarg, double valDistAppPoint
		,TURPointXY *ppntArrAimingPoints, int *piarrRepeatAimingPoints
		, int* piQuantAimingPoints)
{

	double *parrAimingPoints_X = new double  [QUantShells];
	double *parrAimingPoints_Y = new double  [QUantShells];
	//int *piarrRepeatQuants_X = new int  [QUantShells];
 //	int *piarrRepeatQuants_Y = new int  [QUantShells];
	memset(parrAimingPoints_X, 0, QUantShells * sizeof(double));
//	memset(piarrRepeatQuants_X, 0, QUantShells * sizeof(int));

	memset(parrAimingPoints_Y, 0, QUantShells * sizeof(double));
 //	memset(piarrRepeatQuants_Y, 0, QUantShells * sizeof(int));
	///
		plgTarg.calcBoundBox() ;

	// const double VAlDeep = plgTarg.Box[2] - plgTarg.Box[0] ;
	// const double VAlWidth = plgTarg.Box[3] - plgTarg.Box[1] ;
	 // �� ��� X
		double VAlDeep = plgTarg.Box[2] - plgTarg.Box[0];
		double valSig = sqrt(arrElK[0]);
	 double valPp = VAlDeep / valSig /2. + 2.;
	 int iTemp = valPp/ 2.;
	 int iPp =  (valPp/ 2. - ((double)iTemp)< 0.5)? 2 * iTemp + 1: 2 * iTemp + 3;
	 double valStepPp  =  VAlDeep / ((double)iPp - 1.);

	 for (int i =0; i < iPp; i++)
	 {
		parrAimingPoints_X [i] = plgTarg.Box[0] + ((double )i) *  valStepPp;

	 }

	 ///

		// �� ��� Y
		const double VAlWidth =  plgTarg.Box[3] - plgTarg.Box[1];
		valSig = sqrt(arrElK[3]);
		double valC = (VAlWidth + 0.008 * valDistAppPoint)/ valSig/2. ;
	 iTemp = valC / 2.;

	 int iC =  (valC/ 2. - ((double)iTemp)< 0.5)? 2 * iTemp + 1: 2 * iTemp + 3;
	 // ��� ������� ��� ����, ����� �������� ������� �� 0
	 // � ��������  double valStepC  =  VAlWidth / ((double)iC - 1.);
	 // iC ������ ���� ������ >= 2
	 if (iC < 2 )
	 {
	   iC = 2;
	 }

	 double valStepC  =  VAlWidth / ((double)iC - 1.);

	 for (int i =0; i < iC; i++)
	 {
		parrAimingPoints_Y [i] = plgTarg.Box[1] + ((double )i) *  valStepC;
		}

	 *piQuantAimingPoints = _MIN_int( iC * iPp, QUantShells);

	 // ������������ ������� ����� ������������
	 int nx = 0, ny = -1;
		for (int i =0; i < *piQuantAimingPoints; i++)
		{

		 nx = i % iPp;
		 ny = i / iPp;

		 ppntArrAimingPoints [i] = TURPointXY (parrAimingPoints_X [nx], parrAimingPoints_Y [ny] );

		}

		for (int i =0; i < QUantShells; i++)
		{
			 int iCurr = i% (iC * iPp);
			 piarrRepeatAimingPoints[iCurr]++;
		}

	delete []parrAimingPoints_X ;
	delete []parrAimingPoints_Y ;

}
//-----------------------------------------------------------------------------
// �������� ��������� ������� ��� ������ �����������
void TCoastTargNeighbourhood::createMatrxL(const double VAlStepIntegr, double *parr_L)
{
for (int i = 0; i < mMultPntTarg.NumPoints; i++)
	{
		for (int j = 0; j < mMultPntAim.NumPoints; j++)
		{
		 double delX = mMultPntTarg.Points[i].X - mMultPntAim.Points[j].X;
		 double delY = mMultPntTarg.Points[i].Y - mMultPntAim.Points[j].Y;
		 double valp = TProbabilityTheory::calcIntegralNormalDensity_2D(delX, delY
			 ,marrMtrxCorrFluct,   mKillingRange, VAlStepIntegr);
			parr_L[i * mMultPntAim.NumPoints + j] = log(1. - valp)  ;
		}
	}
}
//---------------------------------------------------------
// ���������� ������������� �������� ��������� ��� ����������� ������
// INPUT:
// parrStrategy[mMultPntAim.NumPoints] - ��������� ������������� ���� � ���� �������
double TCoastTargNeighbourhood::calcEfficiencyOfStrategy( double *parrStrategy)
{
		// ������������ �������
	double *parr_L = new double  [mMultPntTarg.NumPoints  * mMultPntAim.NumPoints];
	const double VAlStepIntegr = mKillingRange/10.;
	createMatrxL( VAlStepIntegr,parr_L);

	double valreturn  = calcFGr(parr_L, parrStrategy) ;
	delete []parr_L;
	return valreturn;
}
//--------------------------------------------------------------------------------------------
// ���������� ��������� ����� �������������
// ������ ��������� �������� �������� parrStrategy [mMultPntAim.NumPoints]
// ������� ������� ����� ���������� ����� ������������
// ����������� ����� ������������ ������ �������������� ������� �������.
// �������� ������� ��������� ������ ������ ������������
// INPUT:
// parr_L - ������� ���������� ����������� �������
// parrStrategy - ���������, ��� ������ �������� ��������� �������� ������������ ���� � ������ ����� ������������
// OUTPUT:
// *pinumNext - ����� ����������� ����� ������������
// ���������� ����������� �������� ������� �������
double TCoastTargNeighbourhood::findNextAimingPoint(double *parr_L, double *parrStrategy, int *pinumNext)
{
	double valObj = 1000000000000.;
	for (int i = 0; i < mMultPntAim.NumPoints; i++)
	{
	 double *parrStrategyCur = new double [mMultPntAim.NumPoints];
	 memcpy(parrStrategyCur, parrStrategy, mMultPntAim.NumPoints * sizeof(double));
	 parrStrategyCur[i] += 1.;
	 double valreturn  = calcFGr(parr_L, parrStrategyCur) ;
	 if (valreturn < valObj )
	 {
		valObj =  valreturn;
		*pinumNext = i;
	 }
	 delete  []parrStrategyCur;
	}
	return valObj;
}

//-------------------------------------------------------------------------------------------------
double TCoastTargNeighbourhood::estimateStrategy_MonteCarlo(const int QUantIspit, int *piarrStrategy)
{
	 double Box[4] ={0.};
	 TURPointXY::calcBoundBox__(mMultPntTarg.Points, mMultPntTarg.NumPoints, Box);
		double valProbSum = 0;

		bool *barrHit = new bool[mMultPntTarg.NumPoints];
		memset(barrHit, 0, sizeof(bool) * mMultPntTarg.NumPoints);

		///
	 double arrF_Syst[4] = {0.} , arrMtrxLamb_Syst[4] = {0.}, arrPos00[2] = {0.};
	 CalcProperVectors2(marrMtrxCorrSyst, arrF_Syst, arrMtrxLamb_Syst) ;
	///

  	///
	 double arrF_Fluct[4] = {0.} , arrMtrxLamb_Fluct[4] = {0.};
	 CalcProperVectors2(marrMtrxCorrFluct, arrF_Fluct, arrMtrxLamb_Fluct) ;
	///


	// int iQuantValuePiksels =  rastrTargGSK.calcQuantValuablePiksels();
	//		rastrTargGSK.WriteMeAsFltFile(L"E:\\Ametist\\22-02-2018\\New\\rastrTargGSK.flt");
	for (int i =0; i < QUantIspit; i++)
	{
	 //	TYrRastr rastrTargCur = rastrTargGSK;

	 //	rastrTarg.WriteMeAsFltFile(L"E:\\PROJECTS_C++\\TARAN\\New\\rastrTarg.flt"); //!!!!!
	//	ppntArrCorrecting[0].WriteSetSHPFiles(L"E:\\PROJECTS_C++\\TARAN\\New\\ppntArrCorrecting.shp", ppntArrCorrecting,QUantShells);
	//	plgTarg0.WriteSetSHPFiles(L"E:\\PROJECTS_C++\\TARAN\\New\\plgTarg01.shp", &plgTarg0,1); //!!!!!
	//
	// ������������ ��������������� ������
	double arrDelPos_GSK_Syst[2] = {0.};
	getGaussVector(2, arrPos00,  arrF_Syst, arrMtrxLamb_Syst, arrDelPos_GSK_Syst);
	///
	memset(barrHit, 0, sizeof(bool) * mMultPntTarg.NumPoints);

		for (int j = 0; j < mMultPntAim.NumPoints; j++)
		{

			double arrAimPos_GSK [2] = {0.};
			arrAimPos_GSK[0]  = mMultPntAim.Points[j].X;
			arrAimPos_GSK [1] = mMultPntAim.Points[j].Y;

				for (int k = 0; k < piarrStrategy[j]; k++)
				{
						// ����������� ������� �������� � ������������� �������� � ����
			double  arrDelPos_GSK_Fluct[2] = {0.};



			getGaussVector(2,  arrAimPos_GSK, arrF_Fluct, arrMtrxLamb_Fluct, arrDelPos_GSK_Fluct);

			double arrPos_GSK[2] = {0.};
			 MtrxSumMatrx(arrDelPos_GSK_Fluct, arrDelPos_GSK_Syst,1, 2, arrPos_GSK) ;

			if (
			 (arrPos_GSK[0] > (Box[2]   + mKillingRange) )
			|| (arrPos_GSK[0] < (Box[0] - mKillingRange) )
			|| (arrPos_GSK[1] > (Box[3] + mKillingRange) )
			|| (arrPos_GSK[1] < (Box[1] - mKillingRange) )
			   )
			{
			   continue;
			}

			TURPointXY pntFall ( arrPos_GSK[0],  arrPos_GSK[1] );  // ����� �������
		 //	pntFall.WriteSetSHPFiles(L"E:\\Ametist\\22-02-2018\\NEW\\pntFallCur.shp", &pntFall,1);
		 //	pntFall.WriteSetSHPFiles(L"E:\\Ametist\\22-02-2018\\NEW\\pntAimCur.shp", &ppntArrCorrecting[j],1);

		   // ������� �������� �������� ������

			applyKillingRange( pntFall, barrHit);
	 //		rastrTargCur.WriteMeAsFltFile(L"E:\\Ametist\\22-02-2018\\New\\rastrTargCur.flt"); //!!!!!
			int iii = 0;
				}

	  }
	  double valSum = 0.;

		for (int n = 0; n <  mMultPntTarg.NumPoints; n++)
		{
		if (barrHit[n])
		{
			valSum += 1.;
		}
		}
		double valPTemp = valSum / ((double)mMultPntTarg.NumPoints );
		valProbSum +=  valPTemp;

		 //	rastrTargCur.WriteMeAsFltFile(L"E:\\Ametist\\22-02-2018\\COAST_1D\\plgTarg0.flt"); //!!!!!
	}
 //	*pvalProb =  valProbSum/((double)QUantIspit);


	delete []barrHit;
	return valProbSum/((double)QUantIspit);;
}

//--------------------------------------
void TCoastTargNeighbourhood::applyKillingRange(const TURPointXY  pntFall, bool* barrHit)
{
	for (int i = 0; i < mMultPntTarg.NumPoints ; i++)
	{
		 double valDelY = mMultPntTarg.Points[i].Y - pntFall.Y;
		 if (fabs(valDelY) >  mKillingRange)
		 {
			continue;
		 }
		 double valDelX  =  mMultPntTarg.Points[i].X - pntFall.X;
		 if (fabs(valDelX) >  mKillingRange)
		 {
			continue;
		 }
		 if (barrHit[i])
		 {
			 continue;
		 }
		 if (  (valDelX * valDelX +  valDelY * valDelY -mKillingRange * mKillingRange) <= 0.)
		 barrHit[i] = true;

	}
}

//----------------------------------------------------------------------------------
// ���������� ������� ����������� ����� ������������  ��� ��������� �������� ����

//OUTPUT:
//  parrAimingPoints -  ������ ��������� ����� ������������
//  piarrRepeatQuants - ������ ���������� ����� ������������
//  *pQuantAimingPoints - ����� �- �� ��������� ����� ������������
bool TCoastTargNeighbourhood::calcOptimalArray_Of_AimPoints_With_SystMtrx( TURPointXY *pPntArrAimingPoints
		, int *pQuantAimingPoints, int *piarrRepeatQuants, double *pvalObj)

{
	int iNumArgMin = -1;
	// �������� ������� parr_L  ����������� rstBiasAim.ncols * rstBiasAim.nrows  �����
			//  �  CoastTargNeighbourhood.mMultPntTarg.NumPoints ��������
			// �� ���� ��� ������ ��������������� (��������� �� ������������) ����� ������������  �� ������  rstBiasAim
			//  ����������� ����������� ��������� ������ ����� ����
			 double *parr_L  = new double [mRstrBiasAims.ncols * mRstrBiasAims.nrows * mMultPntTarg.NumPoints ];
			 const double VAlStepIntegr =  mKillingRange/ 10.;
			 createMatrxL_With_SystMtrx( VAlStepIntegr,parr_L);


	/*
	// ������������ �������
	double *parr_L = new double  [ mMultPntTarg.NumPoints *  mMultPntAim.NumPoints];
	const double VAlStepIntegr = mKillingRange /10.;
	///
	createMatrxL(VAlStepIntegr, parr_L) ;
	*/
	///
	double *parrX = new double [mMultPntAim.NumPoints];
	memset(parrX, 0, sizeof(double) * mMultPntAim.NumPoints);

 for (int i=0; i < mMultPntAim.NumPoints; i++)
 {
	 parrX [i] = ((double)mQuantShells)/ ((double)mMultPntAim.NumPoints) ;
 }

	double *parrXZv = new double [mMultPntAim.NumPoints];
	memset(parrXZv, 0, sizeof(double) * mMultPntAim.NumPoints);

double valEps = 0.001 ;

	double valObj0 = calcFGr_With_SystMtrx(parr_L, parrX);
	bool brez = false;
	for (int i = 0; i < 500; i++)
	{
	 // ����������   ��  ���
		 // ������� �������
			double *parr_f= new double [mMultPntAim.NumPoints];
			double *arrGrad= new double [mMultPntAim.NumPoints];
		   //	double *arr_d2FGr_po_dx= new double [mMultPntAim.NumPoints * mMultPntAim.NumPoints];
			double *parrRez = new double [mMultPntAim.NumPoints];
			calcGradFGr_With_SystMtrx(parr_L,  parrX, arrGrad);
			memcpy(parr_f, arrGrad,sizeof(double) * mMultPntAim.NumPoints);
			int numArgMin = -1;
			double valrez = MinDoubleArray(parr_f, mMultPntAim.NumPoints, &numArgMin) ;
			memset( parrXZv, 0, mMultPntAim.NumPoints * sizeof(double));
			parrXZv[numArgMin] =  (double)mQuantShells;
			///



			double *parrRez0 = new double [mMultPntAim.NumPoints];
			*pvalObj = calcMinObjectFunc_AlongWithLine_With_SystMtrx(parr_L,  parrX, parrXZv, parrRez0);


			if (fabs((*pvalObj) - valObj0) < valEps)
			{
				brez = true;
				break;
			}
			memcpy(parrX,parrRez0, mMultPntAim.NumPoints * sizeof(double));
			valObj0 = *pvalObj;
			delete []parr_f;
			delete []arrGrad;
			//delete arr_d2FGr_po_dx;
			delete []parrRez;
			delete []parrRez0;

	}

	delete []parrXZv;

	int *piarr_q = new int [mMultPntAim.NumPoints];   // �-�� ���������� ����� �������� � ������� i
	memset( piarr_q , 0,mMultPntAim.NumPoints * sizeof(int));
 //	double *parr_z =  new double[mMultPntAim.NumPoints];
 //	memset(parr_z , 0,mMultPntAim.NumPoints * sizeof(double));
	double *parrStrategy =  new double [mMultPntAim.NumPoints];
	memset(parrStrategy , 0,mMultPntAim.NumPoints * sizeof(double));
	int iSum = 0;
	for (int i = 0; i < mMultPntAim.NumPoints; i++)
	{
		piarr_q[i] =  ((int)parrX[i]);
		iSum +=  piarr_q[i];
 //		parr_z [i] =   parrX[i] - ((double)piarr_q[i]);
		parrStrategy [i] = (double)piarr_q[i];
	}
	///

	int iJ = mQuantShells - iSum;
	int inumNext = -1;
	double valObjective = -1;

	for (int i = 0; i < iJ; i++)
	{

		*pvalObj  = findNextAimingPoint_With_SystMtrx(parr_L,parrStrategy, &inumNext);
		parrStrategy[inumNext] += 1.;
	}
	///

	for (int i = 0; i < mMultPntAim.NumPoints; i++)
	{
		piarr_q[i] =  ((int)(parrStrategy[i] + 0.05));

	}

	*pQuantAimingPoints = 0;
	int icur = 0;
	for (int i = 0; i < mMultPntAim.NumPoints; i++)
	{
	  if (piarr_q[i] > 0)
		{
		piarrRepeatQuants[icur] =  piarr_q[i];
		pPntArrAimingPoints [icur] = mMultPntAim.Points[i];
		icur++;

	  }
	}
	*pQuantAimingPoints = icur;

	delete []piarr_q;
	delete []parrX;
	delete []parrStrategy;
	delete []parr_L;
	return brez;

}
//-------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ��������  ������� ���������� ��� ������ �����������
// mMultPntTarg.NumPoints  - ����� ���� �������  parr_L
// mRstrBiasAims.ncols * mRstrBiasAims.nrows   - ����� �������� �������  parr_L
void TCoastTargNeighbourhood::createMatrxL_With_SystMtrx(const double VAlStepIntegr, double *parr_L)
{
for (int i = 0; i <  mMultPntTarg.NumPoints; i++)
	{

		for (int j = 0; j <  mRstrBiasAims.ncols * mRstrBiasAims.nrows; j++)
		{
		 TURPointXY pntAimTemp = mRstrBiasAims.getCellCentre(j) ;
		double delX = mMultPntTarg.Points[i].X - pntAimTemp.X;
		double delY = mMultPntTarg.Points[i].Y - pntAimTemp.Y;
		double valp = TProbabilityTheory::calcIntegralNormalDensity_2D(delX, delY
		,marrMtrxCorrFluct,   mKillingRange, VAlStepIntegr);
		parr_L [ i * mRstrBiasAims.ncols * mRstrBiasAims.nrows+ j] = log(1. - valp)  ;
	 }
	}

}
//---------------------------------------------------------


void TCoastTargNeighbourhood::calcGradFGr_With_SystMtrx(double *parr_L,double * parr_U,double * arrGrad)
{
 //	int quantValuableBiasPoints = mRstrProbBias.calcQuantValuablePiksels() ;
	memset(arrGrad, 0, mMultPntAim.NumPoints * sizeof(double));
	double *arrExpSum = new double [mMultPntTarg.NumPoints * mRstrProbBias.ncols* mRstrProbBias.nrows];
	memset(arrExpSum, 0, mMultPntTarg.NumPoints * mRstrProbBias.ncols* mRstrProbBias.nrows * sizeof(double));
	createMtrxExpSum(parr_L, parr_U, arrExpSum);

	for (int i =0; i < mMultPntAim.NumPoints; i++)
	{
	   arrGrad[i] = 0.;
	   for (int j = 0; j < mMultPntTarg.NumPoints; j++)
	   {
			double *parr =  &parr_L [ j * mRstrBiasAims.ncols * mRstrBiasAims.nrows];
			for (int k = 0; k < mRstrProbBias.ncols * mRstrProbBias.nrows; k++)
			{
			 TURPointXY pntBiasOfPointFall = mRstrProbBias.getCellCentre(k) ;
			 TURPointXY pntCur(pntBiasOfPointFall.X + mMultPntAim.Points[i].X - mRstrProbBias.cellsize/2.
			 , pntBiasOfPointFall.Y + mMultPntAim.Points[i].Y - mRstrProbBias.cellsize/2.);
			///
			int iCur = mRstrBiasAims.getPixelNum(pntCur);   //����� ������� � ������ ����������
			if(iCur <0)
			{
			int iii = 0;
			}
			double valLn =  parr[iCur];
			arrGrad[i] += mRstrProbBias.pflt_rastr[k] * valLn * arrExpSum[ j * mRstrProbBias.ncols * mRstrProbBias.nrows +k];
			int iii = 0;
			}
	   }
	}
   //	MatrxMultScalar(arrGrad, 1, mMultPntAim.NumPoints, mRstrProbBias.cellsize * mRstrProbBias.cellsize,arrGrad);
	delete []arrExpSum;
}

 //-------------------------------------------------------------------------
double TCoastTargNeighbourhood::calcFi_jk(double *parr_L, double *parr_U, const int j,const int  k)
{
   TURPointXY pntBiasOfPointFall = mRstrProbBias.getCellCentre(k) ; // ���� ������ ����� ������������
   double sum = 0.;
   double *parr =  &parr_L [ j * mRstrBiasAims.ncols * mRstrBiasAims.nrows];
   for (int i =0; i < mMultPntAim.NumPoints; i++)
   {
	 // ����� ������������
	TURPointXY pntCur(pntBiasOfPointFall.X + mMultPntAim.Points[i].X - mRstrProbBias.cellsize/2.
	, pntBiasOfPointFall.Y + mMultPntAim.Points[i].Y - mRstrProbBias.cellsize/2.);
	///
	int iCur = mRstrBiasAims.getPixelNum(pntCur);
	if(iCur <0)
	{
		int iii = 0;
    }
	sum  +=  parr[iCur] * parr_U[i] ;
   }
	return exp(sum );
}

double TCoastTargNeighbourhood::calcFGr_With_SystMtrx(double *parr_L  ,double * parr_U)
{

	double *arrExpSum = new double [mMultPntTarg.NumPoints * mRstrProbBias.ncols* mRstrProbBias.nrows];
	memset(arrExpSum, 0, mMultPntTarg.NumPoints * mRstrProbBias.ncols* mRstrProbBias.nrows * sizeof(double));
	createMtrxExpSum(parr_L, parr_U, arrExpSum);

	double *parrTemp = new double [mRstrProbBias.ncols* mRstrProbBias.nrows];
	memset(parrTemp, 0,  mRstrProbBias.ncols* mRstrProbBias.nrows * sizeof(double));
	double valRez = 0.;
	for (int k = 0; k < mRstrProbBias.ncols* mRstrProbBias.nrows; k++)
	{
		double sum = 0.;

		for (int j = 0; j < mMultPntTarg.NumPoints; j++)
		{
		sum+= arrExpSum[ j *  mRstrProbBias.ncols* mRstrProbBias.nrows + k];
		}

		valRez +=  sum * mRstrProbBias.pflt_rastr[k];
	}

	delete []arrExpSum;
	delete []parrTemp;
	return valRez ;
}


double TCoastTargNeighbourhood::calcMinObjectFunc_AlongWithLine_With_SystMtrx(double *parr_L
	,double * parr_U, double * parr_UZv, double *parrRez)
{
	double *arrTemp0 = new double [mMultPntAim.NumPoints];
	double *arrTemp1 = new double [mMultPntAim.NumPoints];
	double *arrTemp2 = new double [mMultPntAim.NumPoints];
	double valF =  calcFGr_With_SystMtrx(parr_L,  parr_U) ;
	memcpy(parrRez,  parr_U, mMultPntAim.NumPoints * sizeof(double));
	int iC = 1000;
	for (int i = 1; i < iC; i++)
	{
	 double lamb = ((double)i) / ((double)iC -1.);
	 MatrxMultScalar(parr_U, 1, mMultPntAim.NumPoints, (1. - lamb),arrTemp0);
	 MatrxMultScalar(parr_UZv, 1, mMultPntAim.NumPoints,  lamb,arrTemp1);
	 MtrxSumMatrx(arrTemp0, arrTemp1 ,1, mMultPntAim.NumPoints, arrTemp2) ;
	 double valFTemp = calcFGr_With_SystMtrx(parr_L ,arrTemp2);
		if(valFTemp >= valF)
		{
			break;
		}
		memcpy(parrRez, arrTemp2, mMultPntAim.NumPoints * sizeof(double));
		valF = valFTemp;
	}
	delete []arrTemp0;
	delete []arrTemp1;
	delete []arrTemp2;
	return valF;
}

void TCoastTargNeighbourhood::createMtrxExpSum(double *parr_L, double *parr_U, double *arrExpSum)
{
	for (int j =0; j < mMultPntTarg.NumPoints;j++)
	{
		for (int k = 0; k < mRstrProbBias.ncols* mRstrProbBias.nrows; k++)
		{
		 arrExpSum[j * mRstrProbBias.ncols* mRstrProbBias.nrows + k]  = calcFi_jk(parr_L, parr_U, j, k);
		}
	}
}


//--------------------------------------------------------------------------------------------
// ���������� ��������� ����� �������������
// ������ ��������� �������� �������� parrStrategy [mMultPntAim.NumPoints]
// ������� ������� ����� ���������� ����� ������������
// ����������� ����� ������������ ������ �������������� ������� �������.
// �������� ������� ��������� ������ ������ ������������
// INPUT:
// parr_L - ������� ���������� ����������� �������
// parrStrategy - ���������, ��� ������ �������� ��������� �������� ������������ ���� � ������ ����� ������������
// OUTPUT:
// *pinumNext - ����� ����������� ����� ������������
// ���������� ����������� �������� ������� �������
double TCoastTargNeighbourhood::findNextAimingPoint_With_SystMtrx(double *parr_L, double *parrStrategy, int *pinumNext)
{
  /*	double valObj = 1000000000000.;
	for (int i = 0; i < mMultPntAim.NumPoints; i++)
	{
	 double *parrStrategyCur = new double [mMultPntAim.NumPoints];
	 memcpy(parrStrategyCur, parrStrategy, mMultPntAim.NumPoints * sizeof(double));
	 parrStrategyCur[i] += 1.;
	 double valreturn  = calcFGr_With_SystMtrx(parr_L, parrStrategyCur) ;
	 if (valreturn < valObj )
	 {
		valObj =  valreturn;
		*pinumNext = i;
	 }
	 delete  []parrStrategyCur;
	}
	return valObj; */
	double *arrGrad= new double [mMultPntAim.NumPoints];

	double *parrRez = new double [mMultPntAim.NumPoints];
	calcGradFGr_With_SystMtrx(parr_L,  parrStrategy, arrGrad);
	int iNumArgMin = -1;
	MinDoubleArray(arrGrad, mMultPntAim.NumPoints, &iNumArgMin) ;
	*pinumNext =  iNumArgMin;
	parrStrategy[ iNumArgMin] += 1.;
	double valreturn  = calcFGr_With_SystMtrx(parr_L, parrStrategy) ;
	parrStrategy[ iNumArgMin] -= 1.;
	return valreturn ;
}

// ���������� ������������� �������� ��������� ��� ����������� ������
// ��������� ������ ���� ��������� � �����
// �� ����, ����� ��������� ������ ���� ��������� �������� ������ mRstrBiasAims
// INPUT:
// parrStrategy[mMultPntAim.NumPoints] - ��������� ������������� ���� � ���� �������
double TCoastTargNeighbourhood::calcEfficiencyOfStrategy_With_SystMtrx( double *parrStrategy)
{
	// ������������ �������
	// �������� ������� parr_L  ����������� rstBiasAim.ncols * rstBiasAim.nrows  �����
	//  �  CoastTargNeighbourhood.mMultPntTarg.NumPoints ��������
	// �� ���� ��� ������ ��������������� (��������� �� �����������) ����� ������������  �� ������  rstBiasAim
	//  ����������� ����������� ��������� ������ ����� ����
	double *parr_L  = new double [mRstrBiasAims.ncols * mRstrBiasAims.nrows * mMultPntTarg.NumPoints ];
	const double VAlStepIntegr =  mKillingRange/ 10.;
	createMatrxL_With_SystMtrx( VAlStepIntegr,parr_L);
	double valreturn  = calcFGr_With_SystMtrx(parr_L, parrStrategy) ;
	delete []parr_L;
	return valreturn;
}

//----------------------------------------------------------------------------------

// ���������� ������������� �������� ��������� ��� ����������� ������
// ��������� �� �������  ���� ��������� � �����
// INPUT:
// parrStrategy[mMultPntAim.NumPoints] - ��������� ������������� ���� � ���� �������
// ������� ��������� ��������
double TCoastTargNeighbourhood::calcEfficiencyOfStrategy_With_SystMtrx_Var1 (double *parrStrategy)
{
	const int INumY = mRstrProbBias.ncols * mRstrProbBias.nrows ; // �-�� ����� ���������������� ������
	double valrez = 0.;
	double arrt[2] = {0.};
	const double VAlStepIntegr = mKillingRange / 10.;
	for (int k =0; k< INumY ; k++)
	{
	 TURPointXY pntYk =  mRstrProbBias.getCellCentre(k) ;
	 for (int j =0; j < mMultPntTarg.NumPoints; j++)
	 {
		 TURPointXY pntZjMinusYk(mMultPntTarg.Points[j].X - pntYk.X, mMultPntTarg.Points[j].Y - pntYk.Y);
		 double sum = 0.;
			for (int i = 0; i < mMultPntAim.NumPoints; i++)
			{
				arrt[0] =  pntZjMinusYk.X - mMultPntAim.Points[i].X;
				arrt[1] =  pntZjMinusYk.Y - mMultPntAim.Points[i].Y;
				double valp = TProbabilityTheory::calcIntegralNormalDensity_2D(arrt[0], arrt[1]
				,marrMtrxCorrFluct,   mKillingRange, VAlStepIntegr);
				sum +=  log(1. - valp)  *  parrStrategy[i] ;
			}
		 valrez += exp(sum) * mRstrProbBias.pflt_rastr[k];
	 }
	}
	return valrez;
}

//----------------------------------------------------------------------------------


#pragma package(smart_init)
