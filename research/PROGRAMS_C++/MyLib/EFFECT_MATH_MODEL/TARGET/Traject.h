//---------------------------------------------------------------------------

#ifndef TrajectH
#define TrajectH
#include "InitTargData.h"

//#define LEN_TARG_VS  9
class TInitTargData;
class TTraject
{
public:
	double mTCur; // ������� �����
	double mTBegin ; // ����� ������ ����������

	double marrSigW[3]; // ��� ���� �������
	double marrVectSostGSK_Begin[9]; // ������ ��������� � ��������� ������ - ��� �������
	double marrVectSostGSK[9];     // ������ ��������� �� ������ mTCur
	double marrVectDeviations[9];   // ������ ���������� �� ������ mTCur
 //	THomingHead3D mHomingHead3D ;


	// ����������� �� ���������
	TTraject () ;
	// ����������� �����������
	TTraject  (const TTraject  &R) ;
	// �������� ������������
	TTraject  &operator=(const TTraject   &R2) ;

	// ����� �����������
	TTraject (const double TCur,  const double 	SigW
		,double *arrVectDeviations,double *arrVectSostGSK,double *arrVectSostGSK_Begin);

	TTraject (const double TCur,  double 	*arrSigW,	TInitTargData InitData );

	//
   TTraject (const double TCur,  const double 	SigW,TInitTargData InitData );

   TTraject (const double TBegin,  const double SigW ,double *arrVectSostGSK_Begin );

	int Sign(const double a);



	void recalcTrajPoint(const double tNext) ;

   void extrapolateTargVS(const double VAlTExtr, double *arrTargExtrapVS);

   void extrapolateTarg_BeginVS(const double VAlTExtr, double *arrTargExtrapVS);




}  ;
#endif
