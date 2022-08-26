//---------------------------------------------------------------------------


#pragma hdrstop

#include "ArtComplex.h"


TArtComplex::TArtComplex()
{
	mArtCannon = TArtCannon ();
	mSigU  = 0.;
	mSig_dU_po_dt = 0.;

}



// ����������� �����������
 TArtComplex ::TArtComplex (const TArtComplex &R)
 {
	mArtCannon  =  R.mArtCannon;
	mSigU = R.mSigU;
	mSig_dU_po_dt= R.mSig_dU_po_dt ;
 }
 // �������� ������������
 TArtComplex &TArtComplex::operator=(const TArtComplex  &R)
 {
	mArtCannon  =  R.mArtCannon;
	mSigU = R.mSigU;
	mSig_dU_po_dt= R.mSig_dU_po_dt ;

	return *this ;
 }

  // ����� �����������1
 TArtComplex::TArtComplex (const TArtCannon ArtCannon, const double SigU, const double Sig_dU_po_dt)
 {
	mArtCannon  = ArtCannon ;
	mSigU = SigU ;
	mSig_dU_po_dt = Sig_dU_po_dt;
 }
#pragma package(smart_init)
