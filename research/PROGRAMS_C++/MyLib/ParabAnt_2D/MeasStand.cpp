//---------------------------------------------------------------------------


#pragma hdrstop
#include <string.h>


#include "MeasStand.h"
__fastcall TMeasStand::~TMeasStand()
{
	if(mpcmparrMeas) delete []mpcmparrMeas ;
	mpcmparrMeas = NULL ;

}
//---------------------------------------------------------------------------

TMeasStand::TMeasStand()
{
	// угол наклона оси антенны
	mRadAxeAnt =0.;
	mCountDgrs = 3;
	// массив измерений диаграмм в веере
	mpcmparrMeas = NULL;

	mpcmparrMeas = new TComp[mCountDgrs];
	if(mpcmparrMeas  == NULL)
	{
	// ShowMessage(L"Not memory for Parts") ;
	// return ;
	}
	memset( mpcmparrMeas, 0, 3 * sizeof(TComp)) ;
}


// Конструктор копирования
TMeasStand::TMeasStand (const TMeasStand &R)
 {
   mRadAxeAnt = R. mRadAxeAnt ;
   mCountDgrs = R.mCountDgrs;
   mpcmparrMeas = NULL;

	if(R.mpcmparrMeas != NULL)
	{
		mpcmparrMeas = new TComp[R.mCountDgrs];

		if(mpcmparrMeas == NULL)
		{
	  //	ShowMessage(L"Not memory for Parts") ;
	   //	Abort() ;
		}

		memcpy( mpcmparrMeas,R.mpcmparrMeas, R.mCountDgrs  * sizeof(TComp));
	}

}

  // оператор присваивания
  TMeasStand TMeasStand::operator=(TMeasStand  R)
{
   mRadAxeAnt = R. mRadAxeAnt ;
   mCountDgrs = R.mCountDgrs;
   mpcmparrMeas = NULL;

	if(R.mpcmparrMeas != NULL)
	{
		mpcmparrMeas = new TComp[R.mCountDgrs];

		if(mpcmparrMeas == NULL)
		{
	  //	ShowMessage(L"Not memory for Parts") ;
	   //	Abort() ;
		}

		memcpy( mpcmparrMeas,R.mpcmparrMeas, R.mCountDgrs  * sizeof(TComp));
	}

  return *this ;
}

__fastcall TMeasStand::TMeasStand(const double RadAxeAnt,  TComp *pcmparrMeas)
{
 mRadAxeAnt = RadAxeAnt;
 mCountDgrs =3;
 mpcmparrMeas = NULL;
 mpcmparrMeas = new TComp[3];
 memcpy(mpcmparrMeas, pcmparrMeas, sizeof(TComp) * 3);
}

__fastcall TMeasStand::TMeasStand(const double RadAxeAnt,  int CountDgrs)
{
 mRadAxeAnt = RadAxeAnt;
 mCountDgrs = CountDgrs;
 mpcmparrMeas = NULL;
 mpcmparrMeas = new TComp[CountDgrs];
 memset(mpcmparrMeas, 0, CountDgrs * sizeof(TComp));

}

#pragma package(smart_init)
