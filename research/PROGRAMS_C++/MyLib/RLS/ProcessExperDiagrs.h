//---------------------------------------------------------------------------

#ifndef ProcessExperDiagrsH
#define ProcessExperDiagrsH
class TMeasStand ;
class TURPointXY;
class TURPolyLine;
  class TProcessExperDiagrs
  {
    public:
	static void extractDiagrArray(double *parrData, int iNumRows, int iNumCol, double *parrOut) ;

	static  TURPolyLine extractLevelSubLine( TURPolyLine plnInp, double valLevel);

	static  void  normDiagr( TURPolyLine &plnInp, int &iarg);

	static double findOptCoeff ( TURPolyLine &plnInp, double valXTemp ) ;

	static 	double  approximateDiagr( TURPolyLine &plnInp, double eps,double *pvalCoeff, double *pvalX);

   static double findOptX ( TURPolyLine &plnInp, double valCoeff )  ;

   static double calcNeviazka(TURPolyLine &plnInp, double valCoeff, double valXTemp)  ;

   static void createGraphDiagrAndApproxDiagr(wchar_t *Fold,TURPolyLine &plnInp,  double valCoeff, double valX);

   static int extractFreqNum_FromFileName(wchar_t *FileName);

   static void extractMeasuresArrFromKaurovStandFile(double *parrBuff, const int NUmRows
   , int NUmDiagr, int iNumTriple, TMeasStand *parrMeas);

   static void drawAplitudeGraphs(wchar_t *wchFileName,double *parrBuff, int numDiagrBegin, int countDiagr
   , const int NUmBuffRows, const int NUmBuffCols, TURPointXY pntSdvig, double scalex, double  scaley);

   static void extractMeasuresArrFromKaurovStandFile(double *parrBuff, const int NUmRows
   , int NUmDiagr, int iNumAns, int lenAns,  TMeasStand *parrMeas);

   static void   drawApproximatedAplitudeGraphs(wchar_t *wchFileName1, double *arrDiagrSdvigX
	,double *arrDiagrMultCoeffY, double *arrDiagrMultCoeffX);


  };
#endif
