//---------------------------------------------------------------------------

#ifndef ArtComplexH
#define ArtComplexH
#include "ArtCannon.h"
class TArtCannon;
class TArtComplex
{
public:
// ��
	TArtCannon mArtCannon;
//  �������� ���������  �����
	double mSigU;
// �������� ���������  �������� �����
	double mSig_dU_po_dt;

	// ����������� �� ���������
	TArtComplex () ;
	// ����������� �����������
	TArtComplex  (const TArtComplex  &R) ;
	// �������� ������������
	TArtComplex  &operator=(const TArtComplex  &R2) ;

	// ����� �����������1
	TArtComplex (const TArtCannon ArtCannon, const double SigU, const double Sig_dU_po_dt)  ;

}  ;
#endif
