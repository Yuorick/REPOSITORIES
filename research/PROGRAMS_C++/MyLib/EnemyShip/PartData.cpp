//---------------------------------------------------------------------------


#pragma hdrstop
#include <string.h>
#include "PartData.h"




//---------------------------------------------------------------------------
TPartData::TPartData()
{
	miTypePart    = 0;
	mTimePart = 5. ;
	memset(marrData,0,10 * sizeof(double)) ;
	mSigW  = 33 ;





}

//---------------------------------------------------------------------------


// ����������� �����������
 TPartData ::TPartData (const TPartData &R)
 {

	miTypePart = R.miTypePart;
	mTimePart = R.mTimePart;
	mSigW = R.mSigW ;

	memcpy(marrData,R.marrData, 10 * sizeof(double)) ;


 }
 // �������� ������������
 TPartData TPartData::operator=(TPartData  R)
 {


	miTypePart = R.miTypePart;
	mTimePart = R.mTimePart;
	mSigW = R.mSigW ;
	memcpy(marrData,R.marrData, 10 * sizeof(double)) ;




	return *this ;
 }

  // ����� �����������
 TPartData::TPartData ( const int  iTypePart, const long double 	TimePart ,const long double SigW ,long double *arrData 	)
 {



	miTypePart = iTypePart ;
	mTimePart  =TimePart ;
	mSigW = SigW ;
	memcpy(marrData,arrData, 10 * sizeof(double)) ;




 }


//

#pragma package(smart_init)
