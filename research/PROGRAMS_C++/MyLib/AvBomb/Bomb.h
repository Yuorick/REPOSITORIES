//---------------------------------------------------------------------------

#ifndef BombH
#define BombH
//---------------------------------------------------------------------------

#include "HomingHead.h"
#include "Ebomb.h"

class THomingHead;
class TEbomb ;


// ����� ��������� �������������� ���
// ��� - ���� �������� � ������������� �������� ������ (������� �������������, ��� 152, � 9)
//  � ��������������� ��������
class TBomb
{
public:
 // 1.������� ����� ��������
	long double mMidSq;
// 2.������� ������� ����������� ��������
	long double mSurfFuzSq ;

// 3.�����
	long double mMass ;
// 4.���� ���� �����:
	long double mMaxAt ;
 // 5.���������� �������
	long double mLamb_k ;
 // 6.������� �������, � ��
	long double mSq_k ;
 // 7. ���
	THomingHead mHomingHead;
 // 8.���
	TEbomb mEbomb;
 // 9. ���� ������� ������� ������������ �� ���� �����
	long double mTStabSyst;


	// ����������� �� ���������
	TBomb () ;
	// ����������� �����������
	TBomb  (const TBomb  &R) ;

	// �������� ������������
	TBomb  operator=(TBomb   R2) ;
   // ���������������� �����������
	TBomb(
	 long double MidSq // ������� ����� ��������
	,long double SurfFuzSq // ������� ������� ����������� ��������
	,long double Mass// �����
	,long double MaxAt // ���� ���� �����
	,long double Lamb_k // ���������� �������
	,long double Sq_k // ������� �������, � ��
	,bool bHomeHead // ������� ������� ���
	,THomingHead HomingHead // ���
	,TEbomb Ebomb// ���
	);
	// ���������������� ����������� 2
   TBomb(THomingHead HomingHead );
	// ���������������� ����������� 3
   TBomb(THomingHead HomingHead,TEbomb Ebomb ) ;

   // ���������������� ����������� 4
   TBomb(THomingHead HomingHead,TEbomb Ebomb, const long double TStabSyst );



	long double fncCx( long double valMach);
	long double fncCy( long double valMach);
	long double fnc_dCx_po_dMach( long double valMach);
	long double fnc_dCy_po_dMach( long double valMach);

}  ;
#endif
