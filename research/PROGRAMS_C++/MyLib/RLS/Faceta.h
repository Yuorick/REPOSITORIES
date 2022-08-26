//---------------------------------------------------------------------------

#ifndef FacetaH
#define FacetaH
//#include "Far.h"

class TFaceta
{
public:

// ���������� ����������� � ������
 int m_n;
 // ���������� ����� ������������ � ������
 double  m_d;
 // ����� �����
 double	 mLambda ;
 //  ��� �����. ���� ����������
 double mSigEmitNoise ;
	//  ��� ���� �� ���������( 1 + delta)
	// ��� ����������� ��� � ����������  mSigEmitAmplFact*mSigEmitAmplFact
	// �� ������ ����������� � ����������� ���� �������
	// �������� ���������� ���� � ������ ����� �����
	// mSigEmitNoise*mSigEmitNoise  * m_n +  mSigEmitAmplFact  * mSigEmitAmplFact *  A * A * m_n
	// A - ������ ������� ������������� �����������
 double mSigEmitAmplFact ;


 __fastcall  TFaceta() ;
// ����������� �����������
__fastcall  TFaceta (const TFaceta &R2) ;
 // ����� ������
 __fastcall TFaceta(const int n,const double d,const double Lambda) ;
 // ����� ������
 __fastcall TFaceta(const int n);
 __fastcall TFaceta(const int n,const double d,const double Lambda
	 , const double SigEmitNoise , const double SigEmitAmplFact );

// �������� ������������
TFaceta   &operator=(const TFaceta  &R2) ;
double  fncFSource (const double valTetta);
double  fnc_dFSource_po_dTetta (const double valTetta);
double  fncFFaceta (const double valTetta);
double  fnc_dFFaceta_po_dTet (const double tet);



static   double fnc_dF3_po_tet(const double val_gam, const double tet) ;

static   double fnc_F3(const double val_gam, const double tet) ;

static   double fnc_F4(const double v,const  double u, const double v1,const  double u1) ;

double findDiagrWidth();

void createDiagrGraphs(wchar_t *wchFoldName1 );

 double findDiagrWidthApprox();

 double  fncFFacetaApprox (const double valTetta) ;

 double  fnc_dFFacetaApprox_po_dTet (const double tet);

 double fnc_d2FFacetaApprox_po_dTet2(const double tet);

 double calcNoiseSKZ();

 double calcAmplFactSKZ() ;

};
#endif