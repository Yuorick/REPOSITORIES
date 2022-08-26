//---------------------------------------------------------------------------

#ifndef YrReadH
#define YrReadH
class QString;

   class TYrRead
  {
	 public:

	static  int YrReadCSV_(wchar_t *FileName, const int ncols,int *nrows,double *arrMass);

	static int YrReadCSV(wchar_t *FileName, const int ncols,int *nrows,int *iarrMass) ;

	static int YrCalcRows(wchar_t *FileName);

	static int YrReadColOfCSV(wchar_t *FileName,const int quanTitleRows, const int ncol,const int nrows,bool *barr) ;

	static int YrReadColOfCSV(wchar_t *FileName,const int quanTitleRows
	,const int ncol,const int nrows,double *darr) ;

	static	int YrReadColOfCSV(wchar_t *FileName,const int quanTitleRows
	, const int ncol,const int nrows,int *iarr)  ;

	static int YrReadTabCSV(wchar_t *FileName // файл с таблицей
	,const int quanTitleRows // к-во заголов строк
	,const int quanTitleCols // к-во заголов cтолбцов
	,const int nRowsTab // к-во строк содержательной части таблицы
	,const int nColsTab // к-во столбцов содержательной части таблицы
	,int *nrows  // к-во прочитанных строк
	,int *iarrMass  // массив в который принимается информация
	)  ;

	static int MergeTwoFiles(wchar_t *nameSrcFile1,const int quanTitleRows
	, wchar_t *nameSrcFile2,wchar_t *nameDstFile);

    static int YrReadCharRow(wchar_t *FileName,const int numberRow,char *pw);

    //-- Чтуние массива double из CSV таблицы в массив   iarrMass
    static int YrReadTabCSV_1(wchar_t *FileName // файл с таблицей
    ,const int quanTitleRows // к-во заголов строк
    ,const int quanTitleCols // к-во заголов cтолбцов
    ,const int nRowsTab // к-во строк содержательной части таблицы
    ,const int nColsTab // к-во столбцов содержательной части таблицы
    ,int *nrows  // к-во прочитанных строк
    ,double *parrMass  // массив в который принимается информация
    )  ;
    static int YrReadTabCSV_1(QString qstrFile // файл с таблицей
                      ,const int quanTitleRows // к-во заголов строк
                      ,const int quanTitleCols // к-во заголов cтолбцов
                      ,const int nRowsTab // к-во строк содержательной части таблицы
                      ,const int nColsTab // к-во столбцов содержательной части таблицы
                      ,int *nrows  // к-во прочитанных строк
                      ,double *parrMass  // массив в который принимается информация
                      );



	static int YrReadStrArrFromCSV(wchar_t *FileName // файл с таблицей
	,const int quanTitleRows // к-во заголов строк
	,const int quanTitleCols // к-во заголов cтолбцов
	,const int nRowsTab // к-во строк содержательной части таблицы
	,const int ncols// к-во столбцов массива strNames
	,int *nrows  // к-во прочитанных строк
	,char  *strNames // массив в который принимается информация
	) ;

	static int YrReadColWithCharFromCSV(wchar_t *FileName,const int quanTitleRows
	,const int ncol,const int nrows,const int lenword,char *sarr);

	static int ReadHdrFileFromFltFile(wchar_t*NameOfFltFile, int* nrows,
		 int* ncols, float* xllcorner , float* yllcorner,
		  float* cellsize, float* NODATA_value ) ;

	static int ReadHdrFileFromFltFile(wchar_t*NameOfFltFile, int* nrows,
		 int* ncols, double* xllcorner , double* yllcorner,
		 double* cellsize, double* NODATA_value )  ;

	static int ReadFltFile(wchar_t*NameOfFltFile
		,float* parrX, float* parrY,const float Nodata,float* parrZ);
	static int ReadFltFile(wchar_t*NameOfFltFile
		,double* parrX, double* parrY,const double Nodata,double* parrZ);

	static int ReadSetOfFltRastrs(const wchar_t*FolderName,const wchar_t *TypesFileNames[]
	,const int lenSet, const int nrows, const int ncols,const float xllcorner
	,const float yllcorner,const float cellsize
	,float* pTotalZ,float *pX,float*pY );

	static int ReadHdrOfSetFltRastrs(const wchar_t*FolderName,const wchar_t *TypesFileNames[]
	, int* nrows,int* ncols, float* xllcorner , float* yllcorner,
		 float* cellsize, float* NODATA_value ) ;

    static int ReadHdrFromGearFileForDriver(wchar_t*NameOfSCVFile, bool *bVelo
                                                     , double *pvalL, double *pvalR,  double *pvalPsi_f, double *pvalJ0
                                                     , double *pvalMresidual,double *pvalJLoad, double *pvalCx_om, double *pvalOm0
                                                      , double*pval_a , double*pval_c  , double*pval_l );



	static bool replace(wchar_t*str);

	static bool replace(char*str) ;

	static int CalcQuantMeasuresInKirnosLog(wchar_t*NameOfLogFile);

	static int ReadMeasuresInKirnosLog(wchar_t*NameOfLogFile,const int NumMeasures, double *parrMeasures) ;

	static void  ExtractMeasure (char *str, int *numRow, double *arrMeasCur) ;

	static void  CleaneMeasArrFromZer0(int *quantMeas,double *arrMeas);

	static int ReadDataFromKaurovTXT(wchar_t*NameOfDataFile, int *piNumRows, double *parrData);

	static int calcColCountFromKaurovTXT(wchar_t*NameOfDataFile);

	static int ReadDataFromKaurovTXT_(wchar_t*NameOfDataFile,const int NUmCols, int *piNumRows, double *parrData);

	static void  flip(double *parrData, int len);

    static int  YrCalcRows(QString qstrFile );

    static int calcColCountFromCSV(wchar_t*NameOfDataFile);

    static int ReadHdrFromGearFileForDriver(wchar_t*NameOfSCVFile, bool *bVelo
            , double *pvalL, double *pvalR,  double *pvalPsi_f, double *pvalJ0
            , double *pvalMresidual,double *pvalJLoad, double *pvalCx_om, double *pvalOm0)  ;

    static int ReadOrientationInputFile(wchar_t *wchInputFileName, int *pQuantBaseFreq
           ,double *arrBaseFreqData, int *QuantWorkFreq,double *arrWorkFreq
           , double* arrCavityAntenna, int *piPredeterminationType, double *arrPredeterminatedAngs);



							 
  };

#endif
