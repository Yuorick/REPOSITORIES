//---------------------------------------------------------------------------


#pragma hdrstop
#include "CntlFuncPar.h"
#include <string.h>





//---------------------------------------------------------------------------
TCntlFuncPar::TCntlFuncPar()
{
	miControlType    = 0;
	memset(miArrParams,0,NUM_INT_PARAMS * sizeof(int)) ;
	memset(mDblArrParams,0,NUM_LD_PARAMS * sizeof(long double)) ;
}

//---------------------------------------------------------------------------


// конструктор копирования
 TCntlFuncPar ::TCntlFuncPar (const TCntlFuncPar &R)
 {

	miControlType = R.miControlType;
	memcpy(miArrParams,R.miArrParams, NUM_INT_PARAMS* sizeof(int)) ;
	memcpy(mDblArrParams,R.mDblArrParams, NUM_LD_PARAMS * sizeof(long double)) ;

 }
 // оператор присваивания
 TCntlFuncPar TCntlFuncPar::operator=(TCntlFuncPar  R)
 {
	miControlType = R.miControlType;
	memcpy(miArrParams,R.miArrParams, NUM_INT_PARAMS* sizeof(int)) ;
	memcpy(mDblArrParams,R.mDblArrParams, NUM_LD_PARAMS * sizeof(long double)) ;
	return *this ;
 }

  // парам конструктор
 TCntlFuncPar::TCntlFuncPar ( const int  iControlType, int *iArrParams , long double *DblArrParams)
 {
	miControlType = iControlType;
	memcpy(miArrParams,iArrParams, NUM_INT_PARAMS* sizeof(int)) ;
	memcpy(mDblArrParams,DblArrParams, NUM_LD_PARAMS * sizeof(long double)) ;
 }


//



#pragma package(smart_init)
