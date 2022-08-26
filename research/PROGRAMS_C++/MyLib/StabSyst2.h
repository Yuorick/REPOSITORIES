//---------------------------------------------------------------------------

#ifndef StabSyst2H
#define StabSyst2H
class TComp ;
class TStabSyst2
{
public:
	double mh ;        // ��� �� �������
	double msigm_w ;   //  ��� w
	double msigm_mmo ; //  ��� ���
	double msigm_bmo ; //  ��� ���
	double marrA [4] ;    // ������� �������� �� ���� ����;
	double marrB [2];    // � ������ ������ � ������ �������
	double marrC [2] ;   // �������  ����������
	double msigm_fl; // �������������� ������� ufl � ������ ����� � ��������� x1 = x1 + h * x2 + h * h /2 * w + ufl



	__fastcall ~TStabSyst2() ;

	TStabSyst2 () ;

	// ����������� �����������
	TStabSyst2  (const TStabSyst2  &R) ;
	// �������� ������������
	TStabSyst2  &operator=(const TStabSyst2   &R2) ;
	// ����� ������
	TStabSyst2(  const double h,const double sigm_w, const double sigm_mmo
	,const double sigm_bmo ) ;
	// ����� ������
	// ����� ������
	TStabSyst2(  const double h,const double sigm_w, const double sigm_mmo
	,const double sigm_bmo, const double sigm_fl) ;

	// ������� �����
	void StabSolutionKalm( double *arrK,double * arrP,TComp &rt1,TComp &rt2,
	double &valT1, double &valT2);

	double  SolvCharactEqKlm(const double c);

	static int  CalcProperVectors(double * arrKInp,double *arrV , double *arrLamb);
	void FltKlm( double *arrK,double * arrP);

	void OneStepFltrKlm( double *arrKInp,double *arrKOut,double * arrP);

	bool CalcAccuracyKlm_Real(const double sigm_w,const double sigm_mmo,const double sigm_bmo
	,double *arrK, double *arrPStab
	,const double psi1,const double psi2,double *arrKTol , double **pparrK,int * plenparrarrK) ;

	bool CalcFuncMaximumKlm_Real( const double sigm_w,const double sigm_mmo
	,const double sigm_bmo ,double *arrK, double *arrPStab
	, double &psi1, double &psi2,double *arrKTol , double **pparrK,int * plenparrarrK);

	bool OneStepRecalcTol( const double psi1,const double psi2,const double sigm_w
	,const double  sigm_mmo,const double sigm_bmo,  double *arrP
	,double *arrKInp,double *arrKOut,double *arrS,double *arrQ);

	static bool  CalcCoeff(double fi,double sigm,double *arrK,double *arrS);

	bool FindRoots(double *arrP
	, double *arrS,double *arrQ, TComp &rt1,TComp &rt2);

	void StabSolutionKalm_Real( double *arrK,double * arrP,TComp &rt1,TComp &rt2,
	double &valT1, double &valT2,double **pparrOut, int *plenparrOut, double **pparrWeight0
	, double **pparrWeight1,int * plenparrWeight);

	double  SolvCharactEqKlm_Real(const double mu) ;

	void FltKlm_Real(  double *arrK,double * arrP);

	void OneStepFltrKlm_Real(double *arrKInp,double *arrKOut,double * arrP);

	void CalcWeightArray(double * arrP, double **pparrWeight0
	,double **pparrWeight1, int *plenparrWeight) ;

	static double maxDoubleArr(double *parr, const int lenarr, int &irez);

	static void CreateShpFile(wchar_t*pstrNameFile, double *parrInf
	,const int lenarr, const int quantRows,const double h,double &scalex,double &scaley);

	static void CreateShpAxes(wchar_t *wchFileName,const double xmin,const double xmax
	,const double ymin,const double ymax) ;











	//******************************************************************************************************************//
	//************ ������� ����������� ����������        ******************************************************************************************************//
	//******************************************************************************************************************//
	//******************************************************************************************************************//
	//******************************************************************************************************************//




	,double **pparrWeight1,int * plenparrWeight);











	bool CalcFuncMaximum(double *arrK, double *arrPStab,const double fi1,const double fi2
	, double &psi1, double &psi2,double *arrKTol ) ;
	bool F_MaxForKlmFilt( const double sigm_u
	,const double sigm_bmo ,double *arrK, double *arrPStab
	, double *valPsi,double *valLamb,double * valMu,
	double *arrKTol , double **pparrK,int * plenparrarrK) ;

	bool CalcToleranceKalmFilt( const double sigm_u, const double sigm_bmo  ,double *arrK, double *arrPStab

	,double*  arrKTol , double **pparrK,int * plenparrarrK);


















      void 	ExtrapolStepGolubev_flAdd(const double hExtr,double *arrK,double *arrKExtr);












	bool CalcAccuracy(double *arrK, double *arrPStab,const double fi1,const double fi2
	,const double psi1,const double psi2,double *arrKTol ) ;

	int FSIGN(const double x);
	void RecalcExtendedMtrx(const double valS1,const double valS2,const double sigm_u,const double sigm_bmo , double *arrP
	, const double LambCur, const double MuCur,  const double FiCur ,double *arrKExtend ,double *arrKExtendOut);

	bool IsRootGood(const double sigm_w , const double LambCur, const double MuCur
	,const double FiCur ,double *arrKExtend,const double s0, const double s1 );
	void CalcMatrK(double x,int *iarrPar,double *arrKTemp);



};
#endif