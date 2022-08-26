//---------------------------------------------------------------------------

#ifndef YrWriteH
#define YrWriteH

   class TYrWrite
  {
	 public:


	 static int PutPointsToCsvFile(wchar_t*FileName,float * pX, float * pY,
							 float* pZ,const int lenArray,int * lenVars);
	static int PutPointsToCsvFile(wchar_t*FileName,double* pX, double * pY,
							double* pZ,const int lenArray,int * lenVars) ;

	  static int PutPointsToTxtFile(wchar_t*FileName,float* pX, float * pY,
							 float* pZ,const int lenArray,int * lenVars) ;
	 static int PutPointsToTxtFile(wchar_t*FileName,double* pX, double * pY,
							 double* pZ,const int lenArray,int * lenVars) ;
	 static int WriteMassiveInFltFile(wchar_t*FileName,float * parrZ, const int nrows,
							 const int ncols,const float xllcorner ,const float yllcorner,
							 const float cellsize,const float NODATA_value ) ;
	 static int WriteMassiveInFltFile(wchar_t*FileName,double * parrZ, const int nrows,
							 const int ncols,const double xllcorner ,const double yllcorner,
							 const double cellsize,const double NODATA_value ) ;
	 static int PutSetOfMassivesToPointsTxtFiles(const wchar_t*FolderName,const wchar_t *FilesNames[],const int quanSet,float* pX
						   , float * pY,float* pTotalZ,const int lenArray) ;
	 static int PutSetOfMassivesToPointsTxtFiles(const wchar_t*FolderName,const wchar_t *TypesFileNames[]
								 ,const int lenSet,double* pX
						   , double * pY,double* pTotalZ,const int lenArray)  ;
	 static int WriteSetOfMassivesInIkfFile(const wchar_t*FileName,float * parrZ,const int quanMass
							 ,const int nrows,const int ncols,const float xllcorner
							 ,const float yllcorner, const float cellsize,const float nodat) ;
	 static int WriteMassiveInRastrTextFile(wchar_t*FileName,float * parrZ, const int nrows,
							 const int ncols,const float xllcorner ,const float yllcorner,
							 const float cellsize,const float NODATA_value );
	 static int WriteMassiveInRastrTextFile(wchar_t*FileName,double * parrZ, const int nrows,
							 const int ncols,const double xllcorner ,const double yllcorner,
							 const double cellsize,const double NODATA_value ) ;
	 static int WriteReportForFloatMassiveTXT(const wchar_t*FileName,float*parr,
										  const int ncols,const int nrows) ;
	 static int WriteReportForFloatMassiveTXT(const wchar_t*FileName,double*parr,
										  const int ncols,const int nrows)  ;
	 static int WriteSetOfMassivesInFltFiles(const wchar_t*FolderName
			 ,const wchar_t *TypesFileNames[] ,const int lenSet
						 ,float* pTotalZ, const int nrows, const int ncols,
						  const float xllcorner ,const float yllcorner
							 ,const float cellsize,const float nodata ) ;
	static int WriteSetOfMassivesInFltFiles(const wchar_t*FolderName
			 ,const wchar_t *TypesFileNames[] ,const int lenSet
						 ,double* pTotalZ, const int nrows, const int ncols
						  ,const double xllcorner ,const double yllcorner,
							 const double cellsize,const double nodata ) ;
	static int WriteHdrForIkdFile(const wchar_t*FileName,const int quanMass
							 ,const int nrows,const int ncols,const double xllcorner
							 ,const double yllcorner, const double cellsize
							 ,const double nodat) ;
	static int WriteReportForIntMassiveTXT(const wchar_t*FileName,int*parr,
										  const int ncols,const int nrows);
	 static int WriteReportForIntMassiveTXT(const wchar_t*FileName,long*plarr,
										  const int ncols,const int nrows) ;

static  int WriteReportForFloatMassiveTXT_(const wchar_t*FileName,double*parr,
										  const int ncols,const int nrows)   ;

static int WriteMassiveInFIleSCV(wchar_t*FileName,double *parrBuff, int iNumRows, int iNumCols
							 ,wchar_t *pwcharrRowNames,wchar_t *pwcharrColNames, int iLenName);

static int WriteReportForIntMassiveTXT_Comma(const wchar_t*FileName,double *parr,
										  const int ncols,const int nrows);

static int PutPointsToCsvFile(wchar_t*FileName,double* pX, double * pY,const int lenArray,int * lenVars);


  };

#endif
