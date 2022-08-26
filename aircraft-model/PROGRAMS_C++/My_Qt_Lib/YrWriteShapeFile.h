//---------------------------------------------------------------------------

#ifndef YrWriteShapeFileH
#define YrWriteShapeFileH
class TURPointXY;
class TYrWriteShapeFile
{
public:
static void __fastcall CreateShpFile(wchar_t *wchFileName, double *parrInf, double *parrT
	 ,const int lenarr, double &scalex, double &scaley);
static double maxDoubleArr(double *parr, const int lenarr, int &irez) ;
static double minDoubleArr(double *parr, const int lenarr, int &irez) ;

static void WriteOneReport(wchar_t *wcharrPath  // путь к папке
								  ,double * parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,const int nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  ,const int nBuffRows //  - к-во точек
								  ,wchar_t *wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,const int lenName // максимальная длина имени переменной
								  ,const int numx  // номер переменной по оси Y
								  ,const int numy  // номер переменной по оси X
								  ,const double scalex  //  масштаб по оси Y
								  ,const double scaley  // масштаб по оси X
								   ) ;

static void createFileName( wchar_t *wcharrPath  // путь к папке
										,wchar_t *wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
										,const int lenName // максимальная длина имени переменной
										,const int numx  // номер переменной по оси X
										,const int numy  // номер переменной по оси Y
										,wchar_t *wcharrFileName  // w строка с именем
										);

static void  CreateShpAxes(wchar_t *wchFileName,const double xmin,const double xmax
	 ,const double ymin,const double ymax);

static void  ShowNormProbDistr(wchar_t *wchFileName,const double valA,const double valSigm2)  ;

static void CreateShpArrowedAxes(wchar_t *wchFileName,const double xmin,const double xmax
	 ,const double ymin,const double ymax,const double valLength) ;

static double  Sign_( const double x);

static void CreateAngleMarks(wchar_t *wchFileName, const TURPointXY Pnt0
	 , const TURPointXY Pnt1, const TURPointXY Pnt2,const double Dist0,const double Dist1 ) ;

static void  PictFar(wchar_t *wchFoldName);

static void CreateShpArrowedAxes(wchar_t *wchFileName,const double xmin,const double xmax
	 ,const double ymin,const double ymax,const double valLength, const TURPointXY pnt00);

static void WriteOneReport_Points(wchar_t *wcharrPath1  // путь к папке
                                                     ,double * parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                                     ,const int nBuffCols // - к-во переменных о корорых накоплена информация в буфере
                                                     ,const int nBuffRows //  - к-во точек
                                                     ,wchar_t *wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                                     ,const int lenName // максимальная длина имени переменной
                                                     ,const int numx  // номер переменной по оси X
                                                     ,const int numy  // номер переменной по оси Y
                                                     ,const double scalex  //  масштаб по оси X
                                                     ,const double scaley  // масштаб по оси Y
                                                      );


};


#endif
