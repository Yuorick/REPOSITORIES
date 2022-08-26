//---------------------------------------------------------------------------

#ifndef FloatRastrH
#define FloatRastrH
class TURPointXY;
class TFloatRastr
  {

	 public:
	int mCols;
	int mRows ;
	double mXllcorner ;
	double mYllcorner ;
	double mCellSize ;
	double mNodata ;
	char mcharrByteorder[30]  ;

	// массив с данными. хранятся сверху вниз , слева направо
	float *mparrFlt;


	~ TFloatRastr() ;
	TFloatRastr() ;
	TFloatRastr(wchar_t *wFileName) ;

	TFloatRastr(int nc,int nr,double xll,double yll,double cs,double nd,float *arrData) ;

	TFloatRastr(int nc,int nr,double xll,double yll,double cs,double nd,double *arrData) ;


	TFloatRastr(int nc,int nr,double xll,double yll,double cs,double nd,const bool irastr);

	TFloatRastr (wchar_t *wFileName,const TURPointXY pntP0,const double d0);
	// конструктор копирования
	TFloatRastr (const TFloatRastr &R) ;
	// оператор присваивания
	TFloatRastr &operator=(const TFloatRastr  &R2) ;

	static bool replace(char*str);

	static  void ReadHDRFile (wchar_t*wFileName,int *ncols, int *nrows, double *xllcorner
		 ,double *yllcorner, double *cellsize, double *nodata, char *pbyteorder)  ;

	double   get_X(const int i,const int j) ;
	double   get_Y(const int i,const int j) ;
	double   get_PiksVal(const double x, const double y) ;

	static int WriteMassiveInFltFile(wchar_t*FileName,float * parrZ, const int nrows,
					 const int ncols,const float xllcorner ,const float yllcorner,
					 const float cellsize,const float NODATA_value );

	void  WriteMeAsFltFile(wchar_t *wFileName);

	static int ChangeCellValue(wchar_t*wFileName,const TURPointXY urpntP,const double valZ) ;

	static  double  get_PiksVal(wchar_t*wFileName,const TURPointXY urpntP)  ;

	static TFloatRastr   SumOfRastrs(TFloatRastr *pRastr1, TFloatRastr *pRastr2);
	static TFloatRastr   MinusOfRastrs(TFloatRastr *pRastr1, TFloatRastr *pRastr2);
	static int YMIN(const int n, const int m);

	static int YMAX (const int n, const int m) ;
	bool IsPointInsideExtent(const double x,const double y ) ;
	void GetExtent(double &xmin,double &ymin,double &xmax,double &ymax) ;

	static int  ISign(const double a);
	static double YMINd(const double n, const double m) ;

	static double YMAXd (const double n, const double m) ;

	double  getValue(const int iNum);
	TURPointXY  getCellCentre(const int iNum) ;


/*
	 bool  TFloatRastr:: InitRastr();
	static bool FindParamsForCutRastr(const TFloatRastr *srcRastr, const double xmin
				 ,const double ymin,const double xmax, const double ymax
				  ,int *nrowfirst,int *ncolfirst, int *nrowscut,int *ncolscut) ;
	static void  ExportRastrExtent(TFloatRastr *dstRastr, TFloatRastr *srcRastr)  ;
	static bool  CutFltRastr(TFloatRastr *dstRastr, TFloatRastr *srcRastr,const int nrowfirst
					  ,const int ncolfirst,const int nrowscut, const int ncolscut);
	static bool CustomizeToGrid( TFloatRastr *srcRastr,  TFloatRastr *destRastr);
	static void  WriteRastrsMassiveInIkdFile(wchar_t *FileName,TFloatRastr *pRastr
				 , const int quanRastrs)  ;
	static double BilinearValue (double* arrSh,const double cellsize
			,const double xt,const double yt) ;
	static double NearestPointValue (double* arrSh,const double cellsize
			,const double  Nodata,const double xt,const double yt);
	static double LinApprox(double* arrSh,const double cellsize
			,const double Nodata,const double xt,const double yt) ;
	double InterpolateZ( const double x,const double y );

	public:
   void   ShowMe(int i)  ;

void   ModifyOrderOfDataArray();

bool DigCanyon(const TURPointXY P0,const TURPointXY P1);
TURPointXY OutPnt(const TURPointXY P0,const double kx,const double ky);




// void ModifyDEM(TURPointXY *arrP0,TURPointXY *arrP1, const int quantPoints) ;
static void AdjustDEM( wchar_t *wFileDEM  //исходный  FLT файлом DEM
						 ,wchar_t *wFileSetNull   // исходный  FLT файл с растром потоков типа setnull
						 ,TURPointXY *arrP0,TURPointXY *arrP1, const int quantPoints,TURPointXY *purpntOut);
static bool Link2Streams(const TURPointXY urpntUpper,const TURPointXY urpntLower
					,wchar_t *wFileDEM ,wchar_t *wFileSetNull
					,TURPointXY *urpntStart,TURPointXY *urpntEnd) ;
void CreateCluster(const TURPointXY P0, int **ppiarrNums , int *quantArrNums) ;
int FindSinkIndex (int **ppiarrNums, const int quantNums) ;

 bool IsLowerPointSuit(const TURPointXY P0,const TURPointXY P1 ) ;
 void CreateRasterFromCluster(int *piarrNums , const int quantNums,TFloatRastr &rstrRez)  ;

int InsertRastrInFltFile(wchar_t*FileName);




  */

  };
  #endif
