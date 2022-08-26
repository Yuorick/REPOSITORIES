//---------------------------------------------------------------------------

#ifndef ShipTrajH
#define ShipTrajH
#include "PartData.h"

//---------------------------------------------------------------------------
class TPartData;
class THomingHead3D ;
class TShipTraj
{
public:
 // ��� ����������:
   int miType;


	long double   marrParams[25];

	long double   mTCur; // ������� �����
	long double   mTBegin ; // ����� ������ ����������
	long double   mT;        // ����������������� ���������� �� �������
	long double   mSigW; // ��� ���� �������
	long double   marrVectSostGSK_Begin[9]; // ������ ��������� � ��������� ������ - ��� �������
	long double   marrVectSostGSK[9];     // ������ ��������� �� ������ mTCur
	long double   marrVectDeviations[9];   // ������ ���������� �� ������ mTCur



	// ����������� �� ���������
	TShipTraj () ;
	// ����������� �����������
	TShipTraj  (const TShipTraj  &R) ;
	// �������� ������������
	TShipTraj  operator=(TShipTraj   R2) ;

	// ����� �����������
	TShipTraj (const long double   TCur, const int  iType , const long double   	SigW
	,long double   *arrParams 	,long double   *arrVectDeviations,long double   *arrVectSostGSK,long double   *arrVectSostGSK_Begin);
	// ����� �����������
	TShipTraj( TPartData PartData,long double   *arrS, const long double   TCur) ;
	//
	void  CreateMtrxTansf_SKPM_TO_GSK(long double   *arrS, long double   *arrCentre,long double   *arrN,long double   *arrOut);

	void FindFinalPoint_SKPM(const long double   valV,const long double   valR, long double   valT, long double   *arrSFin) ;

	bool FindNormVect(long double   *arrS,const long double   angGam, long double   *arrN) ;

	int Sign(const long double   a);

	void  FindCircleCentre(long double   *arrN,long double   *arrS, const long double   valR, long double   *arrOut) ;

	long double    norm(long double   *arr);

	void  VectMult(long double   *px, long double   *py, long double   *arrRez) ;

	int getPartNum(const long double   valT) ;

	void TrasformSSKPM_TO_GSK(const int numPart,long double   *arrSSKPM, long double   *arrSGSK) ;

	void calcAddDeviations_Gorka(const long double   h, long double   *arrOut) ;

	void recalcSys2(long double   *arrA, long double   *arrB, const long double   valNoise,long double   *arrVectDeviations) ;

	void reformVectDev(long double   *arrInp,long double   *arrOut) ;

	void createTransfMtrx(const long double   alf, long double   *arrMtrxOut);

	void recalcTrajPoint_2(const long double   tNext) ;

	void recalcTrajPoint_3(const long double   tNext)  ;



	static long double   getRand01( );

	bool recalcTrajPoint(const long double   tNext) ;

	void recalcTrajPoint_0(const long double   tNext) ;



	static long double   getGauss( const long double   a, const long double   sig) ;
    static long double    Rand_();

}  ;
#endif
