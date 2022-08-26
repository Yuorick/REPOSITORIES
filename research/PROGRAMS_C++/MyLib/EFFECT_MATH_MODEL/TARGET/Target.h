//---------------------------------------------------------------------------

#ifndef TargetH
#define TargetH

#include "Traject.h"
#include "InitTargData.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "PlanePolygon.h"
#include "URPolygon.h"


#define MAX_POSSIBLE_QUANTS_PLANE_POLYGONS 20
enum enumTargetType{ GARPUN_V300, GARPUN_V700, JET_F16 ,JET_A_10A, DESTROYER, MOTORBOAT
, OPEN_MANPOWER_LIE, OPEN_MANPOWER_STAND, BULLET_PROOF_LIE, BULLET_PROOF_STAND, COVERED_MANPOWER_ENTRENCH, COVERED_MANPOWER_TRENCH
,MANPOWER_ARMOURED_CARRIER, MANPOWER_CAR, PLATOON_POINT, GROUP_POINT_COAST, GARPUN_V1000};

class TTraject;
class TInitTargData ;
class TURPolyLine;
class TPlanePolygon;

class TTarget
{
public:
	 // ����������
	 TTraject mTraject ;

	 // ������ ������������
	 TURPolyLine mplnGlagkovsky;
	 // ��� ����
	 enum enumTargetType  menumTargetType  ;
	 // ���:
		double mTargEPR ;

		// �������� ���������
	   //	TPlanePolygon *mpArrPlanePolygon ;// ������ ��������� �� ���������
	   TPlanePolygon mpArrPlanePolygon [MAX_POSSIBLE_QUANTS_PLANE_POLYGONS ];
	   int mLenArrPlanePolygon ; // ����� ������� ���������


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
	 ~TTarget() ;
	// ����������� �� ���������
	TTarget () ;
	// ����������� �����������
	TTarget  (const TTarget  &R) ;
	// �������� ������������
	TTarget  &operator=(const TTarget   &R2) ;
 /*
	// ����� �����������
TTarget (const double Bearing, const double TargCourse
	, const double TargZenitAng,  const double V, const double H
	,const double R, const double 	SigW, const double valT
	, TURPolyLine plnGlagkovsky
	,enumTargetType  enumTargetType,  wchar_t* pwcharrFoldReport);
		// ����� �����������
	TTarget (const TInitTargData InitData ,  const double 	SigW
	, TURPolyLine plnGlagkovsky
	,enumTargetType  EnumTargetType,  wchar_t* pwcharrFoldReport);

 TTarget (const TInitTargData InitData ,  const double 	SigW
	, TURPolyLine plnGlagkovsky,enumTargetType  EnumTargetType, const double TargEPR, wchar_t* pwcharrFoldReport);

TTarget (TTraject   Traject
	, TURPolyLine plnGlagkovsky
	,enumTargetType  EnumTargetType, const double TargEPR, wchar_t* pwcharrFoldReport);
 */
TTarget  (TTraject   Traject	,enumTargetType  EnumTargetType, const double TargEPR, wchar_t* pwcharrFoldReport);

void recalcTrajPoint(const double tNext);
	 void WriteReport() ;

void updateReportData() ;

void WriteReport(wchar_t *pwcharrPath)  ;

static void createVesselArrayOfPlanePolygons(const double VAlBoardX, const double VAlBoardY, const double VAlBoardH
  ,const double VAlStructureX, const double VAlStructureY, const double VAlStructureH
   ,const double VAlSdvgStructureX, TPlanePolygon *pPlanePolygonArr );




}  ;

#endif