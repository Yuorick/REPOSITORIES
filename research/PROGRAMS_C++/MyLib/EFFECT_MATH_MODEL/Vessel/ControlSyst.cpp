//---------------------------------------------------------------------------


#pragma hdrstop

#include "ControlSyst.h"


TControlSyst::TControlSyst()
{
	//���� ����������
		 mFiltT  =0. ;
	// �������� ����
		 mSinsDelayT  =0. ;
		 mRzvT = 0.;



}

//---------------------------------------------------------------------------


// ����������� �����������
 TControlSyst ::TControlSyst (const TControlSyst &R)
 {
	mFiltT = R.mFiltT;
	mSinsDelayT = R.mSinsDelayT;
	mRzvT = R.mRzvT;
	mTblPereletSimulated = R.mTblPereletSimulated;

 }
 // �������� ������������
 TControlSyst &TControlSyst::operator=(const TControlSyst  &R)
 {
	mFiltT = R.mFiltT;
	mSinsDelayT = R.mSinsDelayT;
	mRzvT = R.mRzvT;
	mTblPereletSimulated = R.mTblPereletSimulated;
	return *this ;
 }

	// ����� �����������1
// TControlSyst::TControlSyst (const double FiltT, const double SinsDelayT)
// {
//	mFiltT  = FiltT;
//	mSinsDelayT = SinsDelayT ;
 //	mRzvT = 0.;
// }
	// ����� �����������2
 TControlSyst::TControlSyst (const double FiltT, const double SinsDelayT, const double RzvT)
 {
	mFiltT  = FiltT;
	mSinsDelayT = SinsDelayT ;
	mRzvT = RzvT;
 }

  TControlSyst::TControlSyst (const double FiltT, const double SinsDelayT
  , const double RzvT, const TTable_1D TblPereletSimulated)
 {
	mFiltT  = FiltT;
	mSinsDelayT = SinsDelayT ;
	mRzvT = RzvT;
	mTblPereletSimulated = TblPereletSimulated;
 }


#pragma package(smart_init)
