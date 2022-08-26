//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include <stdlib.h>
#include <String.h>

#include "YrString.h"
#include <vcl.h>
#include <stdio.h>
double StrTo_Dbl_(String strng)
{
	wchar_t arrwch[300] = {0};
	wcscpy(arrwch,strng.w_str());
	int n = wcslen(arrwch);
	bool bNumber = false ;
	for (int i = 0; i < n; i++)
	{
		bNumber = false ;
	   for (int j = 0; j < 10; j++)
	   {

		 String s_22 = j ;
		 s_22 =  s_22 + L"j";
		 wchar_t arrwcht[3] = {0};
		 wcscpy(arrwcht,s_22.w_str());
		 if (arrwch[i] == arrwcht[0])
		 {
		   bNumber = true;
		   break;
		 }
		 s_22 = -1 ;
		 wcscpy(arrwcht,s_22.w_str());
		 if (arrwch[i] == arrwcht[0])
		 {
		   bNumber = true;
		   break;
		 }
	   }
	   //if (arrwch[i] == L'0')
	 //  {
	   //    bNumber = true;
      // }
	   if (!bNumber)
	   {
		 arrwch[i]= L'.';

		 break ;
	   }
	}
//double bb = _wtof(arrwch);
	return _wtof(arrwch);

}

#pragma package(smart_init)
