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

static void WriteOneReport(wchar_t *wcharrPath  // ���� � �����
								  ,double * parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,const int nBuffCols // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,const int nBuffRows //  - �-�� �����
								  ,wchar_t *wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,const int lenName // ������������ ����� ����� ����������
								  ,const int numx  // ����� ���������� �� ��� Y
								  ,const int numy  // ����� ���������� �� ��� X
								  ,const double scalex  //  ������� �� ��� Y
								  ,const double scaley  // ������� �� ��� X
								   ) ;

static void createFileName( wchar_t *wcharrPath  // ���� � �����
										,wchar_t *wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
										,const int lenName // ������������ ����� ����� ����������
										,const int numx  // ����� ���������� �� ��� X
										,const int numy  // ����� ���������� �� ��� Y
										,wchar_t *wcharrFileName  // w ������ � ������
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


};


#endif
