//---------------------------------------------------------------------------

#ifndef ShipTargH
#define ShipTargH
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#include "PartData.h"
#include "ShipTraj.h"
#include "InitTargData.h"
#include "TargData.h"

class TPartData ;
class TTraject;
class TInitTargData ;
class TTargData;

//extern const int QUANT_TYPES_OF_TRAJECT = 3;
class TShipTarg
{
public:

// ������ �� ���������� �������� ����������  mquantParts<= 10
   int mquantParts ;
// ������ ����������, ����������� ������� ����������(�� ����� 10 ��)
   TPartData marrPartData[10];
   // ������� �����
   long double mT;
   // ����������
   TShipTraj mTraject ;

   // ������ �� ���������� �������
   TTargData mTargData;

	 // ��� ������
	 // �-��� ����� ����������
	int mQuantPntReport ;

	 //�������� ��������������� ���������������� ������ ��� ����� ������
	int mLenMemoryAlloc ;

	// ����� ������
	 double *mparrBuff    ;

	// ���� � ����� � �������
	wchar_t *mpwcharrFoldReport ;


     // ����������
	 ~TShipTarg() ;
	// ����������� �� ���������
	TShipTarg () ;
	// ����������� �����������
	TShipTarg  (const TShipTarg  &R) ;
	// �������� ������������
	TShipTarg  operator=(TShipTarg   R2) ;

	// ����� �����������
	TShipTarg ( const long double Bearing, const long double TargCourse
	, const long double TargZenitAng,  const long double V, const long double H ,
	const long double R
	,const long double valT
	,const int quantParts, TPartData *arrPartData, wchar_t* pwcharrFoldReport);;
		// ����� �����������
	TShipTarg (const TInitTargData InitData
	,const int quantParts, TPartData *arrPartData, wchar_t* pwcharrFoldReport);

	 // ����� �����������  3
 TShipTarg (const TInitTargData InitData,const int quantParts, TPartData *arrPartData,TTargData TargData, wchar_t* pwcharrFoldReport);

	bool recalcTrajPoint(const long double tNext);
	int  getNumCurrentPart(const long double valt) ;
	long double  getTimePartStarted(const long double valt) ;
	long double  getTimePartWillFinish(const long double valt) ;
	int getNumM�neuvreType(const long double valt);
	 void WriteReport() ;

	 void updateReportData() ;



}  ;
#endif
