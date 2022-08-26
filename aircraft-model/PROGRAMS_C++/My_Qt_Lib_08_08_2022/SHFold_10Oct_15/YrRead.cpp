
#include "YrRead.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <QFile>
#include <QString>
#include <QStringList>
//---------------------------------------------------------------------------


 #include <math.h>
extern const double NODATA;
extern const double NODATA1 ;

extern const double DNIL  ;
extern const double DNIL1  ;
extern const double EPSSING   ;
extern const double EPS   ;


//---------------------------------------------------------------------------
//-- чтение CSV файла в DOUBLE массив  arrMass
//-- Вх инф:
//-- 1. FileName - файл
//-- 2.  ncols - к-во столбцов
//-- Вых инф:
//--  1. nrows - к-во прочит строк
//--  arrMass - массив в который произведена запись
int TYrRead::YrReadCSV_(wchar_t *FileName, const int ncols,int *nrows,double *arrMass)
   {
		FILE *fr;
	   if(( fr=_wfopen(FileName,L"r")) == NULL )
	   {
		   //ShowMessage(L"Not posiible to open .CSV file");

		   return 1;
       }
		char str[10000];

		int nrows1 =0;
		while (fgets(str,10001,fr))
			{

				if(str[0] == 0) break;
				int lenline = strlen(str);
				 str[lenline] = ';';
				 str[lenline + 1] = 0 ;
				 lenline++;
				for (int i =0; i < lenline; i++)
				{
				  if (str[i] == ',') str[i] = '.';


				}
				char str1[50] ={0};
				int i_temp = 0;
				int ncols_temp = 0;
				for (int i =0; i < lenline; i++)
				{
				   if (str[i]!=';')
				   {
					 str1[i_temp] = str[i];
					  i_temp++;
				   }
				   else
				   {
					   arrMass[nrows1*ncols + ncols_temp] = atof(str1);

					   ncols_temp++;
						i_temp = 0;
						for (int j = 0; j < 50; j++) str1[j] =0;

				   }

				}
			   nrows1++;
            }

		fclose(fr);
		*nrows = nrows1 ;
	   return 0 ;
   }

int TYrRead::YrReadCSV(wchar_t *FileName, const int ncols,int *nrows,int *iarrMass)
	// пергруженная под INT
   {
		FILE *fr;
		fr=_wfopen(FileName,L"r");
		char str[1000];
		fgets(str,1001,fr);

		int nrows1 =0;
		while (fgets(str,1001,fr))
			{

				if(str[0] == 0) break;
				int lenline = strlen(str);
				 str[lenline] = ';';
				 str[lenline + 1] = 0 ;
				 lenline++;
				for (int i =0; i < lenline; i++)
				{
				  if (str[i] == ',') str[i] = '.';


				}
				char str1[50] ={0};
				int i_temp = 0;
				int ncols_temp = 0;
				for (int i =0; i < lenline; i++)
				{
				   if (str[i]!=';')
				   {
					 str1[i_temp] = str[i];
					  i_temp++;
				   }
				   else
				   {
					   iarrMass[nrows1*ncols + ncols_temp] = atoi(str1);

					   ncols_temp++;
						i_temp = 0;
						for (int j = 0; j < 50; j++) str1[j] =0;

				   }

				}
			   nrows1++;
            }

		fclose(fr);
		*nrows = nrows1 ;
	   return 0 ;
   }

//-- Чтуние массива INT из CSV таблицы в массив   iarrMass
int TYrRead::YrReadTabCSV(wchar_t *FileName // файл с таблицей
					,const int quanTitleRows // к-во заголов строк
					,const int quanTitleCols // к-во заголов cтолбцов
					,const int nRowsTab // к-во строк содержательной части таблицы
					,const int nColsTab // к-во столбцов содержательной части таблицы
					,int *nrows  // к-во прочитанных строк
					,int *iarrMass  // массив в который принимается информация
					)

   {
		FILE *fr;
	   if((fr=_wfopen(FileName,L"r")) ==  NULL)
	   {
		   //ShowMessage(L"ERROR IN YrReadTabCSV\n Not possible to open FileName") ;
		   return 1;
       }
		char str[1000];
		for (int i = 0; i <quanTitleRows ; i++) fgets(str,1001,fr);



		int nrows1 =0;
	for ( int k = 0; k < nRowsTab ; k++)
	 {
		if (fgets(str,1001,fr))
		{
				if(str[0] == 0) break;
				int lenline = strlen(str);
				 str[lenline] = ';';
				 str[lenline + 1] = 0 ;
				 lenline++;
				for (int i =0; i < lenline; i++)
				{
				  if (str[i] == ',') str[i] = '.';


				}
				char str1[50] ={0};
				int i_temp = 0;
				int ncols_temp = 0;
				for (int i =0; i < lenline; i++)
				{
				   if (str[i]!=';')
				   {
					 str1[i_temp] = str[i];
					  i_temp++;
				   }
				   else
				   {
					   if( ( ncols_temp - quanTitleCols ) >= 0 )
					   iarrMass[nrows1*nColsTab  + ncols_temp- quanTitleCols ] = atoi(str1);

					   ncols_temp++;
						i_temp = 0;
						for (int j = 0; j < 50; j++) str1[j] =0;

				   }

				}
			   nrows1++;
		  }
		}

		fclose(fr);
		*nrows = nrows1 ;
	   return (nrows1 == nRowsTab)?0:1 ;
   }

int TYrRead::YrCalcRows(wchar_t *FileName)
	{
	  FILE *fr;
	 if( (fr=_wfopen(FileName,L"r")) == NULL )
	  {
		  //String s_22 =  L"Error in  YrCalcRows.\nNot possible to open file" ;
		  //ShowMessage(s_22 + FileName) ;
		  //Abort();
	  }
		//fr=_wfopen(FileName,L"r");
		char str[1000];

		int nrows1 =0;
		while (fgets(str,1001,fr))
			{
				if(str[0] == 0) break;

			   nrows1++;
            }

		fclose(fr);
	   return nrows1 ;
   }

int TYrRead::YrReadColOfCSV(wchar_t *FileName,const int quanTitleRows, const int ncol,const int nrows,bool *barr)
  // чтение столбца с номером ncol из CSV файла в bool массив barr  (ncol >=0 and <= ncols-1)
  // quanTitleRows - к-во строк заголовока
  {
	  FILE *fr;
	 if( (fr=_wfopen(FileName,L"r")) == NULL )
	  {
		  //ShowMessage(L"Error in YrReadColOfCSV bool\nNot possible to open file") ;
		  //Abort();
	  }
	   char str[1000];
	   char *st = str;
	   char *st1 = str ;
	   for (int i = 0; i < quanTitleRows; i++) 	fgets(str,1001,fr);

	   for (int i = 0; i < nrows; i++)
	   {

			fgets(str,1001,fr);
		 int k = strlen(str);
		 str[k] = ';';
		 str[k+1] = 0;

		 st = str;
		 for (int j = 0; j < ncol ; j++)
		 {
		   st1=	strchr(st,';');
		   st = st1+1;
		 }
		 char strrez[100] = {0};

		 for (int j = 0; j < 1000; j++)
		 {
		   if (st[j] == ';') break ;

		   strrez[j] =	st[j];
		 }
		 int irez = atoi(strrez) ;

		 barr[i] = ( irez == 0 )?false:true;
	   }
	   fclose(fr) ;
	  return 0 ;
  }

int TYrRead::YrReadColOfCSV(wchar_t *FileName,const int quanTitleRows, const int ncol,const int nrows,int *iarr)
  // чтение столбца с номером ncol из CSV файла в bool массив barr  (ncol >=0 and <= ncols-1)
  // quanTitleRows - к-во строк заголовока
  {
	  FILE *fr;
	if( (fr=_wfopen(FileName,L"r")) == NULL )
	  {
		  //ShowMessage(L"Error in  YrReadColOfCSV int.\nNot possible to open file") ;
		  //Abort();
	  }
	   char str[1000];
	   char *st = str;
	   char *st1 = str ;
	   for (int i = 0; i < quanTitleRows; i++) 	fgets(str,1001,fr);

	   for (int i = 0; i < nrows; i++)
	   {

			fgets(str,1001,fr);
		 int k = strlen(str);
		 str[k] = ';';
		 str[k+1] = 0;

		 st = str;
		 for (int j = 0; j < ncol ; j++)
		 {
		   st1=	strchr(st,';');
		   st = st1+1;
		 }
		 char strrez[100] = {0};

		 for (int j = 0; j < 1000; j++)
		 {
		   if (st[j] == ';') break ;

		   strrez[j] =	st[j];
		 }
		 iarr[i] = atoi(strrez) ;
	   }
	   fclose(fr) ;
	  return 0 ;
  }

int TYrRead::YrReadColOfCSV(wchar_t *FileName,const int quanTitleRows
		   ,const int ncol,const int nrows,double *darr)
  // чтение столбца с номером ncol из CSV файла в DOUBLE массив darr  (ncol >=0 and <= ncols-1)
  // quanTitleRows - к-во строк заголовока
  {
	  FILE *fr;
		if( (fr=_wfopen(FileName,L"r")) == NULL )
	  {
		//String s_22 =  L"Error in  YrReadColOfCSV double.\nNot possible to open file" ;
		  //ShowMessage(s_22 + FileName) ;
		  //Abort();
	  }
	   char str[1000];
	   char *st = str;
	   char *st1 = str ;
	   for (int i = 0; i < quanTitleRows; i++) 	fgets(str,1001,fr);

	   for (int i = 0; i < nrows; i++)
	   {

			fgets(str,1001,fr);
		 int k = strlen(str);
		 str[k] = ';';
		 str[k+1] = 0;

		 st = str;
		 for (int j = 0; j < ncol ; j++)
		 {
		   st1=	strchr(st,';');
		   st = st1+1;
		 }
		 char strrez[100] = {0};

		 for (int j = 0; j < 1000; j++)
		 {
		   if (st[j] == ';') break ;

		   strrez[j] =	st[j];
		 }
		 for (int k = 0; k < (int)strlen(strrez); k++)
		 {
			if( strrez[k]  == ',') strrez[k]  = '.' ;
		 }
		 darr[i] =(double) atof(strrez) ;


		// cout<<"i = "<< i<< "barr[i] = "<<barr[i] <<endl;
	   }
	   fclose(fr) ;
	  return 0 ;
  }

int TYrRead::YrReadColWithCharFromCSV(wchar_t *FileName,const int quanTitleRows
		   ,const int ncol,const int nrows,const int lenword,char *sarr)
  // чтение столбца WCHAR_T с номером ncol из CSV файла в WCHAR_T массив wharr
  // cдлиной  слова  lenword.  (ncol >=0 and <= ncols-1)
  // quanTitleRows - к-во строк заголовока
  {
	  FILE *fr;
	 if( (fr=_wfopen(FileName,L"r")) == NULL )
	  {
		  //ShowMessage(L"Error in  YrReadColWithWcharFromCSV.\nNot possible to open file") ;
		  //Abort();
	  }

	   char str[1000];
	   char *st = str;
	   char *st1 = str ;
	   for (int i = 0; i < quanTitleRows; i++) 	fgets(str,1001,fr);

	   for (int i = 0; i < nrows; i++)
	   {

			fgets(str,1001,fr);
		 int k = strlen(str);
		 str[k] = ';';
		 str[k+1] = 0;

		 st = str;
		 for (int j = 0; j < ncol ; j++)
		 {
		   st1=	strchr(st,';');
		   st = st1+1;
		 }
		 char strrez[100] = {0};

		 for (int j = 0; j < 1000; j++)
		 {
		   if (st[j] == ';')
		   {
			strrez[j] == '0' ;
			break ;
		   }

		   strrez[j] =	st[j];
		 }
		 if ( (int)strlen(strrez) > lenword )
		 {
			//ShowMessage(L"Error in YrReadColWithWcharFromCSV:\nlen word doesn't comply lenword") ;
			 //Abort() ;
		 }
		strcpy( &sarr[i * lenword] , strrez) ;


	   }
	   fclose(fr) ;
	  return 0 ;
  }

// слияние 2-х CSV файлов  nameSrcFile1    и  nameSrcFile2
//  nameSrcFile1 - имеет  quanTitleRows zagolovochnih strok
//  nameSrcFile2 - ne imeet zagolovochnih strok
// nameDstFile - rezultat sliania fajlov, imeet zagolovochnyu stroky iz nameSrcFile1
  int TYrRead::MergeTwoFiles(wchar_t *nameSrcFile1,const int quanTitleRows, wchar_t *nameSrcFile2,wchar_t *nameDstFile)
{
	  //String s_22 ;
	 int i1 = YrCalcRows(nameSrcFile1);

	 int i2 = YrCalcRows(nameSrcFile2);
	 if( (i1 - quanTitleRows) != i2)
	 {
	   //ShowMessage(L" Merged files have different row's quanttity. ERROR.");
	   return -1;
	 }
	FILE *fr1;
	FILE *fr2;
	FILE *fw;
	fr1=_wfopen(nameSrcFile1,L"r") ;
	 if(fr1 == NULL)
	 {
		//String s_22 = L" ERROR. MergeTwoFiles\n Not possible to open " ;
	   //ShowMessage(s_22+nameSrcFile1);
	   return -1;
	 }
   //	 //ShowMessage(L"Merged files 10 ") ;

	fr2=_wfopen(nameSrcFile2,L"r") ;

	 if(fr2 == NULL)
	 {
	   //String s_22 = L" ERROR. MergeTwoFiles\n Not possible to open " ;
	   //ShowMessage(s_22+nameSrcFile2);
	   return -1;
	 }
	////ShowMessage(L"Merged files 20 ") ;

	fw = _wfopen(nameDstFile,L"w") ;
	if(fw == NULL)
	 {
	   //String s_22 = L" ERROR. MergeTwoFiles\n Not possible to open " ;
	   //ShowMessage(s_22+nameDstFile);
	   return -1;
	 }
	const int num = 1000;
	 char  pch1[num] ={0}
		  ,pch2[num] ={0}
		  ,pch3[num] ={0} ;
	  for (int i = 0; i < quanTitleRows ; i++)
	  {
		 fgets(pch1,num,fr1);
		 fputs(pch1,fw) ;
		 memset( pch1,0,num);
	  }


	  for (int i = 0; i < i2; i++)
	  {


		  fgets(pch1,num,fr1);
		 //  s_22 =  pch1 ;
		//  //ShowMessage(s_22 +L"  после считывания ") ;

		// const int len0 = strlen(pch1) -1 ;
		  pch1[strlen(pch1) -1] = ';';

		  pch1[strlen(pch1)] = 0;
		 //   s_22 =  pch1 ;
		 // //ShowMessage(s_22 +L"  до цикла ") ;

		 /* for ( int j = len0; j >0; j--)
		  {
			if ((pch1[j -1] == ';' ) || (!isgraph(pch1[j -1] )) )
			  {
			   //	pch1[j] == 0 ;  изм 07 Oct 2011
					pch1[j] = 0 ;
				j--;
			   }
			else break ;

		  }  */
		 // s_22 =  pch1 ;
		 // //ShowMessage(s_22 +L"  после цикла ") ;
		  fgets(pch2,num,fr2);
		 // //ShowMessage(L"Merged files 50 ") ;
		  strcpy(pch3,pch1);
		 // //ShowMessage(L"Merged files 51 ") ;
		  strcat (pch3,pch2);
		 // //ShowMessage(L"Merged files 52 ") ;
		  fputs(pch3,fw) ;
		 // s_22 = L" i = ";
		 // //ShowMessage(s_22 + i +L" pch1 = " + pch1+ +L" pch2 = " + pch2 +L" pch3 = " + pch3);
		  memset( pch1,0,num);
		  memset( pch2,0,num);
		  memset( pch3,0,num);
	 }
	// //ShowMessage(L"Merged files 40 ") ;
	fclose(fr1);
	fclose(fr2);
	fclose(fw);
	return 0;
}

int TYrRead::YrReadCharRow(wchar_t *FileName,const int numberRow,char *pw)
  // чтение строки с номером numberRow из CSV файла в wchar_t * массив pw[30 *QUANT_CROPS]
  //
  {
	   FILE *fr;
	 if(  (fr=_wfopen(FileName,L"r")) == NULL)
	 {
	   //ShowMessage(L"Not possible to open file Predecessor.csv") ;
	   //Abort() ;
	 }
		char wstr[1000] = {0};
			char wstr_temp[30] ={0};
	 for (int i =0 ; i < numberRow; i++) fgets(wstr,1001,fr);
	 int ilen = strlen(wstr);
	 int it = 0;
	 int icur = 0 ;
	 for (int i = 0; i < ilen; i++)
	 {
		if( ( wstr[i] != ';' )&& ( wstr[i] != '\n' ) )
		 {
			 wstr_temp[it] =  wstr[i];
			 it++;
		 }
		 else
		 {
		   wstr_temp[it] = 0 ;
		   strcpy( &pw[icur * 30], wstr_temp );
		   it = 0;
		   icur++ ;
		   memset(wstr_temp,0,sizeof(char)*30) ;
         }
	 }


	   fclose(fr) ;




	  return icur ;
  }


//-----------------------------------------------------------------------------

//-- Чтуние массива double из CSV таблицы в массив   iarrMass
int TYrRead::YrReadTabCSV_1(wchar_t *FileName // файл с таблицей
                  ,const int quanTitleRows // к-во заголов строк
                  ,const int quanTitleCols // к-во заголов cтолбцов
                  ,const int nRowsTab // к-во строк содержательной части таблицы
                  ,const int nColsTab // к-во столбцов содержательной части таблицы
                  ,int *nrows  // к-во прочитанных строк
                  ,double *parrMass  // массив в который принимается информация
                  )

 {
      FILE *fr;
     if((fr=_wfopen(FileName,L"r")) ==  NULL)
     {

         int iii=0;
     }
      char str[1000];
      for (int i = 0; i <quanTitleRows ; i++)
      {
        if(NULL ==  fgets(str,1001,fr))
        {
            return -1;
        }
      }



      int nrows1 =0;
  for ( int k = 0; k < nRowsTab ; k++)
   {
      if (fgets(str,1001,fr))
      {
              if(str[0] == 0) break;
              int lenline = strlen(str);
               str[lenline] = ';';
               str[lenline + 1] = 0 ;
               lenline++;
              for (int i =0; i < lenline; i++)
              {
                if (str[i] == ',')
                {
                    str[i] = '.';
                }

              }
              char str1[50] ={0};
              int i_temp = 0;
              int ncols_temp = 0;
              for (int i =0; i < lenline; i++)
              {
                 if (str[i]!=';')
                 {
                   str1[i_temp] = str[i];
                    i_temp++;
                 }
                 else
                 {
                     if(( ( ncols_temp - quanTitleCols ) >= 0 )
                            && ( ( ncols_temp - quanTitleCols ) < nColsTab ))
                     parrMass[nrows1*nColsTab  + ncols_temp- quanTitleCols ] = atof(str1);

                     ncols_temp++;
                      i_temp = 0;
                      for (int j = 0; j < 50; j++) str1[j] =0;

                 }

              }
             nrows1++;
        }
      }

      fclose(fr);
      *nrows = nrows1 ;
     return (nrows1 == nRowsTab)?0:1 ;
}


//-- Чтуние массива double из CSV таблицы в массив   parrMass

int TYrRead::YrReadTabCSV_1(QString qstrFile // файл с таблицей
                  ,const int quanTitleRows // к-во заголов строк
                  ,const int quanTitleCols // к-во заголов cтолбцов
                  ,const int nRowsTab // к-во строк содержательной части таблицы
                  ,const int nColsTab // к-во столбцов содержательной части таблицы
                  ,int *nrows  // к-во прочитанных строк
                  ,double *parrMass  // массив в который принимается информация
                  )

 {

    QFile qFile;  
    qFile.setFileName(qstrFile);
    qFile.open(QIODevice::ReadOnly | QIODevice::Text);
    int nrows1 =0;
    QString line0;
    for (int i = 0; i <quanTitleRows ; i++)
    {
      line0 = qFile.readLine();

    }
    for (int i = 0; i < nRowsTab ; i++)
    {
      line0 = qFile.readLine();
      QStringList w =line0.split(";");
      for (int j = 0; j < nColsTab; j++)
      {
          QString line2 = w.at(quanTitleCols + j);
          QString y = ",";
          if(line2.indexOf(y) >= 0)
               {
               line2.replace(line2.indexOf(y), 1, '.');

               }



        parrMass[i * nColsTab + j] = line2.toDouble();
        double x = parrMass[i * nColsTab + j];
        int iii=0;
      }
      nrows1++;
    }
    qFile.close();
    *nrows = nrows1 ;
   return (nrows1 == nRowsTab)?0:1 ;

}

// INPUT:
// qstrFile - файл с таблицей
int  TYrRead::YrCalcRows(QString qstrFile )//
    {
    QFile qFile;
    qFile.setFileName(qstrFile);
    qFile.open(QIODevice::ReadOnly | QIODevice::Text);

    QString line0;
    qint64 i =0;
    char data[10000] = {0};
    qint64 maxSize = 2000000;
    for ( i = 0; i <2000000 ; i++)
    {
      //line0 = qFile.readLine();
      if (qFile.readLine(data,  maxSize) <= 0) break;


    }
    qFile.close();
    return i -1;


   }

//---------------------------------------------------------------------------
   //-- Чтуние массива строк  из CSV таблицы в массив   strNames с кодлич столбцов ncols
   // строка из CSV файла начинается с новой строки массива   strNames
   int TYrRead::YrReadStrArrFromCSV(wchar_t *FileName // файл с таблицей
					,const int quanTitleRows // к-во заголов строк
					,const int quanTitleCols // к-во заголов cтолбцов
					,const int nRowsTab // к-во строк содержательной части таблицы
					,const int ncols// к-во столбцов массива strNames
					,int *nrows  // к-во прочитанных строк
					,char  *strNames // массив в который принимается информация
					)

   {
		memset(strNames,0, ncols * nRowsTab) ;
		FILE *fr;
	   if((fr=_wfopen(FileName,L"r")) ==  NULL)
	   {
		   //ShowMessage(L"ERROR IN YrReadTabCSV\n Not possible to open FileName") ;
		   //Abort() ;
       }
		char str[1000];
		for (int i = 0; i <quanTitleRows ; i++) fgets(str,1001,fr);



		int nrows1 =0;
	for ( int k = 0; k < nRowsTab ; k++)
	 {
		if (fgets(str,1001,fr))
		{
		   if(str[0] == 0) break;
		   int lenline = strlen(str);
		   char *pstr = str ;
		   if (quanTitleCols > 0 )
		   {
			int ncols_temp = 0;

			 for (int i =0; i < lenline; i++)
			 {
				if (str[i] ==';')
				 {
				  ncols_temp++;
					if (ncols_temp == quanTitleCols)
					 {
					  pstr = &str[i +1 ] ;
					  break ;
					 }
				 }
			  }
		  }

		 strcpy( &strNames[ncols * nrows1], pstr) ;
		  nrows1++;
		}
		else break ;
	  }

		fclose(fr);
		*nrows = nrows1 ;
	   return (nrows1 == nRowsTab)?0:1 ;
   }
int TYrRead::ReadHdrFileFromFltFile(wchar_t*NameOfFltFile, int* nrows,							 int* ncols, float* xllcorner , float* yllcorner,
							 float* cellsize, float* NODATA_value )
 //
// NameOfFltFile - имя файла header которого читается. Расширение  должно быть.flt
// nrows   ncols  - к-во строк и столбцов соответственно
//  xllcorner  - координата x левого нижнего угла
//   yllcorner  - координата y левого нижнего угла
// cellsize    - шаг решетки
//  NODATA_value  - значение которым заполняется массив в случае отсутствия величины
{  int len = wcslen(NameOfFltFile) ;	if ( !( (NameOfFltFile[len - 1] == 't') && (NameOfFltFile[len - 2] == 'l') // проверка, что
	 && (NameOfFltFile[len - 3] == 'f') ) )  // указанный файл имеет расширение .flt
	{
	 //String St = NameOfFltFile ;
	  //ShowMessage(L"TYrRead::ReadHdrFileFromFltFile\nExtention of " + St + L"is wrong") ;
	  return 1 ;
	}
	wchar_t* HdrFileName = new wchar_t[len] ;
	wcscpy( HdrFileName, NameOfFltFile) ;
	 HdrFileName[len - 1] = 'r';
	 HdrFileName[len - 2] = 'd';
	 HdrFileName[len - 3] = 'h';
	FILE *fr ;	if ((fr = _wfopen(HdrFileName,L"r"))== NULL)	{
	 //String St =  HdrFileName ;
	 //ShowMessage(L"TYrRead::ReadHdrFileFromFltFile\n ERROR ! Not possible to open " +St) ;
	 return 1 ;
    }	int ncols1,nrows1;	double xllcorner1,yllcorner1,cellsize1,NODATA_value1 ;
    char  str_xll[50],str_yll[50],str_sell[50],str_nodat[50];
    /*fscanf(fr,"ncols         %d\n",&ncols1) ;
	  fscanf(fr,"nrows         %d\n",&nrows1) ;	  fscanf(fr,"xllcorner     %s\n",str_xll) ;	  fscanf(fr,"yllcorner     %s\n",str_yll) ;	  fscanf(fr,"cellsize      %s\n",str_sell ) ;	  fscanf(fr,"NODATA_value  %s\n",str_nodat  ) ;  */	  	  if	 (		( fscanf(fr,"ncols         %d\n",&ncols1)== 0)	 ||( fscanf(fr,"nrows         %d\n",&nrows1 ) == 0)	 || ( fscanf(fr,"xllcorner     %s\n",str_xll) == 0)	 || ( fscanf(fr,"yllcorner     %s\n",str_yll)== 0 )	 || ( fscanf(fr,"cellsize      %s\n",str_sell )== 0)	 || (fscanf(fr,"NODATA_value  %s\n",str_nodat  ) == 0 )	   )	   {		   //ShowMessage(L"Error in format of HDR File");
		   //Abort() ;
	   }	 TYrRead::replace(str_xll);	  TYrRead::replace(str_yll);	   TYrRead::replace(str_sell);		TYrRead::replace(str_nodat);	fclose(fr) ;	delete  [] HdrFileName ;	*ncols =ncols1;	*nrows =nrows1;	*xllcorner =atof(str_xll);	*yllcorner = atof(str_yll);	*cellsize = atof(str_sell);	*NODATA_value = atof(str_nodat);	return 0 ;
}

//--------------------------------------------------------------
int TYrRead::ReadHdrFileFromFltFile(wchar_t*NameOfFltFile, int* nrows,							 int* ncols, double* xllcorner , double* yllcorner,
							 double* cellsize, double* NODATA_value )
 // ПОЛИМОРФНАЯ ДЛЯ DOUBLE
// NameOfFltFile - имя файла header которого читается. Расширение  должно быть.flt
// nrows   ncols  - к-во строк и столбцов соответственно
//  xllcorner  - координата x левого нижнего угла
//   yllcorner  - координата y левого нижнего угла
// cellsize    - шаг решетки
//  NODATA_value  - значение которым заполняется массив в случае отсутствия величины
{  int len = wcslen(NameOfFltFile) ;	if ( !( (NameOfFltFile[len - 1] == 't') && (NameOfFltFile[len - 2] == 'l') // проверка, что
	 && (NameOfFltFile[len - 3] == 'f') ) )  // указанный файл имеет расширение .flt
	{
	 //String St = NameOfFltFile ;
	  //ShowMessage(L"TYrRead::ReadHdrFileFromFltFile\nExtention of " + St + L"is wrong") ;
	  return 1 ;
	}
	wchar_t* HdrFileName = new wchar_t[len] ;
	wcscpy( HdrFileName, NameOfFltFile) ;
	 HdrFileName[len - 1] = 'r';
	 HdrFileName[len - 2] = 'd';
	 HdrFileName[len - 3] = 'h';
	FILE *fr ;	if ((fr = _wfopen(HdrFileName,L"r"))== NULL)	{
	 //String St =  HdrFileName ;
	 //ShowMessage(L"TYrRead::ReadHdrFileFromFltFile\n ERROR ! Not possible to open " +St) ;
	 return 1 ;
    }	int ncols1,nrows1;
    double xllcorner1,yllcorner1,cellsize1,NODATA_value1 ;
    char  str_xll[50],str_yll[50],str_sell[50],str_nodat[50];
    /*fscanf(fr,"ncols         %d\n",&ncols1) ;
	  fscanf(fr,"nrows         %d\n",&nrows1) ;	  fscanf(fr,"xllcorner     %s\n",str_xll) ;	  fscanf(fr,"yllcorner     %s\n",str_yll) ;	  fscanf(fr,"cellsize      %s\n",str_sell ) ;	  fscanf(fr,"NODATA_value  %s\n",str_nodat  ) ;  */	  if	 (		( fscanf(fr,"ncols         %d\n",&ncols1)== 0)	 ||(  fscanf(fr,"nrows         %d\n",&nrows1 ) == 0)	 || ( fscanf(fr,"xllcorner     %s\n",str_xll) == 0)	 || ( fscanf(fr,"yllcorner     %s\n",str_yll)== 0 )	 || ( fscanf(fr,"cellsize      %s\n",str_sell )== 0)	 || ( fscanf(fr,"NODATA_value  %s\n",str_nodat  ) == 0 )	   )	   {		   //ShowMessage(L"Error in format of HDR File");
		   //Abort() ;
       }
    TYrRead::replace(str_xll);
    TYrRead::replace(str_yll);
    TYrRead::replace(str_sell);	 TYrRead::replace(str_nodat);		////String s_22= L" str_nodat = " ;		////ShowMessage(s_22 + str_nodat) ;	fclose(fr) ;	delete  [] HdrFileName ;	*ncols =ncols1;	*nrows =nrows1;	*xllcorner =atof(str_xll);	*yllcorner = atof(str_yll);	*cellsize = atof(str_sell);	*NODATA_value = atof(str_nodat);	return 0 ;
}


bool TYrRead::replace(wchar_t*str)  {	 int n = wcslen(str);
	 for (int i=0; i < n; i++)
	 {
		 if (str[i] == L',')
		 {
			str[i] = L'.';
			return true;
		 }
	 }
	 return false ;
  }

bool TYrRead::replace(char*str)  {	 int n = strlen(str);
	 for (int i=0; i < n; i++)
	 {
		 if (str[i] == ',')
		 {
			str[i] = '.';
			return true;
		 }
	 }
	 return false ;
  }


int TYrRead::ReadFltFile(wchar_t*NameOfFltFile							,float* parrX, float* parrY,const float Nodata,float* parrZ)
//  Функция предназначкена для чтения массива из файла типа .flt
//  результат работы функции хранится в массиве   parrZ , parrX и parrY
//   В файле . flt массив хранится в порядке слева направо, сверху вниз
//   У нас принято хранить слева направо снизу вверх (в порядке возрастания координат X и Y)
//  Т е  parrZ[0] - это левый нижний элемент
//   Массивы parrZ , parrX и parrY должны быть созданы в подпрограмме из которой происходит
//   обращение. Т е предварительно надо при помощи функции TYrRead::ReadHdrFileFromFltFile
// прочитать содержание .hdr и определить длину массива, а затем создать его. Если один из указателей parrX или parrY  ==NULL, то заполнение этих
// массивов не производится
{
   int nrows = 0,ncols = 0 ;
   float  xllcorner = 0,yllcorner =0, cellsize = 0,  NODATA_value = 0;
   TYrRead::ReadHdrFileFromFltFile(NameOfFltFile, &nrows,
							&ncols, &xllcorner , &yllcorner,
							  &cellsize, &NODATA_value )  ;

  float *Ztemp = new float [ nrows * ncols] ;
	if (Ztemp == NULL)
	{
	   //ShowMessage(L"TYrRead::ReadFltFile\nThere are not memory for Ztemp in TYrRead::ReadFltFile") ;
	   return 1 ;
	}

	FILE *fr ;
   if (	(fr=_wfopen(NameOfFltFile,L"rb")) == NULL)
   {
	   //String s_22 =  L"TYrRead::ReadFltFile\nNot possible to open file" ;
	   //ShowMessage(s_22+ NameOfFltFile) ;
	   return 1 ;
   }
	 ////String s_22 = L"ReadFltFile_nrows = " ;
	//  //ShowMessage(s_22 + nrows +L" ncols= " +  ncols );
	fread(Ztemp,sizeof(float),nrows * ncols,fr) ;
   //	TYrWrite::WriteReportForFloatMassiveTXT(L"C:\\zemproekt\\DebugMV\\New\\ReadFltFile_Ztemp.txt"
	//							 , Ztemp,4,3) ;

		fclose(fr);
	 int i1 = 0;
	 int j1 = 0;
   //  s_22 = L"ReadFltFile_Ztemp = " ;
	for ( i1 = 0; i1 < nrows; i1++)
	for ( j1 = 0; j1 < ncols; j1++)
	{
		parrZ[i1*ncols + j1]=   Ztemp[(nrows-1-i1)*ncols +j1] ;
	   //	if (parrZ[i1*ncols + j1] <= (NODATA_value +1)) parrZ[i1*ncols + j1] = Nodata;  изм 01.07.2011

  //	  изм 01.07.2011
     if ( ( (parrZ[i1*ncols + j1]/2 -0.5) < (NODATA_value/2) )  && ( (parrZ[i1*ncols + j1]/2 + 0.5) > (NODATA_value/2))  )
		parrZ[i1*ncols + j1] = Nodata;




	}

   	if( (parrX!=NULL) && (parrY!=NULL) )
	{
		for (int i = 0; i < nrows; i++)
		for (int j =0; j < ncols; j++)
		{
		  parrX[i*ncols + j] = xllcorner + cellsize/2 + j * cellsize ;
		  parrY[i*ncols + j] = yllcorner + cellsize/2 + i * cellsize ;

		}

	}

   delete [] Ztemp ;
   return 0 ;
}

 int TYrRead::ReadFltFile(wchar_t*NameOfFltFile
							,double* parrX, double* parrY,const double Nodata,double* parrZ)
 // ПОЛИМОРФНАЯ ДЛЯ DOUBLE
//  Функция предназначкена для чтения массива из файла типа .flt
//  результат работы функции хранится в массиве   parrZ , parrX и parrY
//   В файле . flt массив хранится в порядке слева направо, сверху вниз
//   У нас принято хранить слева направо снизу вверх (в порядке возрастания координат X и Y)
//  Т е  parrZ[0] - это левый нижний элемент
//   Массивы parrZ , parrX и parrY должны быть созданы в подпрограмме из которой происходит
//   обращение. Т е предварительно надо при помощи функции TYrRead::ReadHdrFileFromFltFile
// прочитать содержание .hdr и определить длину массива, а затем создать его. Если один из указателей parrX или parrY  ==NULL, то заполнение этих
// массивов не производится
{
   //	 //String s_22 = L"TYrRead::ReadFltFile \nNameOfFltFile  = " ;
	//  //ShowMessage(s_22 + NameOfFltFile);

   int nrows = 0,ncols = 0 ;
   double  xllcorner = 0,yllcorner =0, cellsize = 0,  NODATA_value = 0;
   TYrRead::ReadHdrFileFromFltFile(NameOfFltFile, &nrows,
							&ncols, &xllcorner , &yllcorner,
							  &cellsize, &NODATA_value )  ;

  float *Ztemp = new float [ nrows * ncols] ;
	if (Ztemp == NULL)
	{
	   //ShowMessage(L"TYrRead::ReadFltFile\nThere are not memory for Ztemp in TYrRead::ReadFltFile") ;
	   return 1 ;
	}

	FILE *fr ;
   if (	(fr=_wfopen(NameOfFltFile,L"rb")) == NULL)
   {
	   //String s_22 =  L"TYrRead::ReadFltFile\nNot possible to open file" ;
	   //ShowMessage(s_22+ NameOfFltFile) ;
	   return 1 ;
   }
   //	//String s_22 = L"ReadFltFile_nrows = " ;
	 // //ShowMessage(s_22 + nrows +L" ncols= " +  ncols + L"Nodata = " + Nodata + L"NODATA_value = " + NODATA_value);
	fread(Ztemp,sizeof(float),nrows * ncols,fr) ;
   //	TYrWrite::WriteReportForFloatMassiveTXT(L"C:\\zemproekt\\DebugMV\\New\\ReadFltFile_Ztemp.txt"
	//							 , Ztemp,4,3) ;

		fclose(fr);
	 int i1 = 0;
	 int j1 = 0;
   //  s_22 = L"ReadFltFile_Ztemp = " ;
	for ( i1 = 0; i1 < nrows; i1++)
	for ( j1 = 0; j1 < ncols; j1++)
	{
		parrZ[i1*ncols + j1]=   Ztemp[(nrows-1-i1)*ncols +j1] ;
	  //	if (parrZ[i1*ncols + j1] <= (NODATA_value +1)) parrZ[i1*ncols + j1] = Nodata;  изм 01.07.2011

  //	  изм 01.07.2011
   if (fabs( parrZ[i1*ncols + j1]/2 - NODATA_value/2) < 1.) parrZ[i1*ncols + j1] = Nodata;

	}

	if( (parrX!=NULL) && (parrY!=NULL) )
	{
		for (int i = 0; i < nrows; i++)
		for (int j =0; j < ncols; j++)
		{
		  parrX[i*ncols + j] = xllcorner + cellsize/2 + j * cellsize ;
		  parrY[i*ncols + j] = yllcorner + cellsize/2 + i * cellsize ;

		}

	}

   delete [] Ztemp ;
   return 0 ;
}
  
int TYrRead ::ReadSetOfFltRastrs(const wchar_t*FolderName,const wchar_t *TypesFileNames[]
			 ,const int lenSet, const int nrows, const int ncols,const float xllcorner
			 ,const float yllcorner,const float cellsize
			 ,float* pTotalZ,float *pX,float*pY )
 //Производится чтение набора растров в количестве lenset  в формате FLT в массив
 //  pTotalZ.
 // Входная информация:
 // 1. FolderName - путь к папке (с "хвостом") в которой храняться растры .FLT
 // например L"C:\\ZemProekt\\DebugMV\\Flt_Pr_MV\\NeprMV_";
 // 2. TypesFileNames[] - массв типов файлов. Например, {L"Z.",L"Zm.",L"p.",L"q.",L"pm.",L"qm."
 //	,L"rm.",L"sm.",L"tm."};
 // 3. lenSet - к-во растров с именами NeprMV_Z.flt, NeprMV_Zm.flt,NeprMV_p.flt и т.д.
 // 4.,5  nrows,  ncols  - к-во строк и столбцов в рвстре
 // 6.,7.  xllcorner, yllcorner - координата X,Y левого нижнего угла сетки
 //      Внимание! Значение привязано в центру квадрата(пиксела)
 // 8.  cellsize - шаг сетки(сторона пиксела)
 //  Выходная информация:
 // 1. pTotalZ - массив длиной  nrows*ncols*lenSet в который поочередно записываются  растры
 // 2. pX и  pY - массивы коородинат X,Y соответственно.
 //    Если pX или pY ==NULL, то эти массивы не вычисляются
 // Руководство по использованию:
 // 1. Надо прочитать припомощи функции  TYrRead ::ReadHdrOfSetFltRastrs
 // параметры  nrows,  ncols, xllcorner, yllcorner, cellsize
 //  TYrRead ::ReadHdrOfSetFltRastrs(FolderFltFilesInterpGeomVals
 //          ,GeomFilesName//String, &nrows, &ncols, &xllcorner , &yllcorner
 // 		 ,&cellsize, &NODATA_value ) ;
 // 2.Создать массивы типа float
 //  pTotalZ - размер 9 * nrows*ncols
 // pX,pY  - размер  nrows*ncols
 // 3. обратиться к настоящей функции.
 // после того, как эта функция отработала выставить указатели на массивы геометрических
 // величин в порядке Z ,Zm, p, q, pm, qm, rm, sm, tm с шагом  nrows*ncols

{
	 wchar_t FileName[300] ;

   wchar_t* ExtStr = L"flt" ;



   for (int i = 0; i < lenSet; i++)
   {
   wcscpy(FileName,FolderName) ;

   wcscat(FileName,TypesFileNames[i]) ;
   wcscat(FileName,ExtStr) ;
  // //String s_22  =  L"ReadSetOfFltRastrs:FileName = "  ;
   //	   //ShowMessage( s_22 + FileName) ;

	 TYrRead::ReadFltFile(FileName,NULL,NULL,NODATA,&pTotalZ[i * nrows * ncols]);

   }

	if ((pX != NULL) && (pY != NULL))
	{
	   for (int i = 0; i < nrows; i++)
	   for (int j =0 ;  j < ncols; j++)
	   {
		   pX[i * ncols + j] =  xllcorner +  cellsize/2 + j * cellsize;
		   pY[i * ncols + j] =  yllcorner +  cellsize/2 + i * cellsize;
	   }

	}

	return 0 ;
}

int TYrRead ::ReadHdrOfSetFltRastrs(const wchar_t*FolderName,const wchar_t *TypesFileNames[]
			 , int* nrows,int* ncols, float* xllcorner , float* yllcorner,
							 float* cellsize, float* NODATA_value )
 {
     wchar_t FileName[300] ;

	 wchar_t* ExtStr = L"flt" ;

	 wcscpy(FileName,FolderName) ;
	 wcscat(FileName,TypesFileNames[0]) ;
	 wcscat(FileName,ExtStr) ;
	 return	 TYrRead::ReadHdrFileFromFltFile(FileName, nrows,
							 ncols, xllcorner , yllcorner,
							  cellsize, NODATA_value ) ;

 }


// вычисление к-ва измерений записанных в файле .LOG кирноса
int TYrRead::CalcQuantMeasuresInKirnosLog(wchar_t*NameOfLogFile)
{
  int len = wcslen(NameOfLogFile) ;	if ( !( (NameOfLogFile[len - 1] == 'g') && (NameOfLogFile[len - 2] == 'o') // проверка, что
	 && (NameOfLogFile[len - 3] == 'l') ) )  // указанный файл имеет расширение .flt
	{
	 //String St = NameOfLogFile;
	  //ShowMessage(L"TYrRead::CalcQuantMeasuresInKirnosLog\nExtention of " + St + L"is wrong") ;
	  return -1 ;
	}

	FILE *fr ;	if ((fr = _wfopen(NameOfLogFile,L"r"))== NULL)	{
	 //String St =  NameOfLogFile ;
	 //ShowMessage(L"TYrRead::ReadHdrFileFromFltFile\n ERROR ! Not possible to open " +St) ;
	 return -1 ;
	}	char str[1000];	int numMeasures = 0;
   //	fseek
	//fread
	char *pch = NULL;
	char *pch0 = NULL;
	bool bMeasure = true;
	 for (int i = 0; i < 1000000; i++)
	 {

	   pch =	fgets(str,1001,fr);
	   if (!pch )break;

	   if(!strstr(pch, "airlow")) continue;

	  // char *pch0 =  strstr(pch, "vel");
	   pch0 = pch0 + 3;
	   int ivel0 = atoi(pch0);
		bMeasure = true;
	   for (int j=0; j < 3; j++)
	   {
		pch =	fgets(str,1001,fr);

			if(!strstr(pch, "airlow"))
			{
				bMeasure = false;
				break;
			}
			else
			{
			 char *pch1 =  strstr(pch, "vel");
			 pch1 = pch1 + 3;
			 int ivel1 = atoi(pch1);
			 if (ivel1 != ivel0)
			 {
			   bMeasure = false;
			   break;
			 }

			  // изменение 18 июня 2016
			 // double arrMeasCur[2] = {0.};
			//  int numRow = 0;
			//  ExtractMeasure (str, &numRow,arrMeasCur);
			 // if (sqrt(arrMeasCur[0] * arrMeasCur[0] + arrMeasCur[1] * arrMeasCur[1])< 0.000001)
			 // {
			 //	bMeasure = false;
			 //	break;
			 // }

			  ///
			}
	   }
	   if (bMeasure)
	   {
		numMeasures++;
	   }

		// //ShowMessage( str) ;
	 }


	 fclose(fr);	return numMeasures ;
}
// чтение измерений из файла .LOG кирноса
// предварительно надо обратиться к функции CalcQuantMeasuresInKirnosLog
// чтобы вычислить к-во измерений изаказать память для parrMeasures
// память  нужна 8* NumMeasures* sizeof(double)
//INPUT:
//NameOfLogFile - имя файла
//NumMeasures  -
// OUTPUT:
// parrMeasures -  массив измеренных сигналов
// порядок хранения информации
// S0.Re  S0.Im  S1.Re  S1.Im S2.Re  S2.Im  S3.Re  S3.Im
int TYrRead::ReadMeasuresInKirnosLog(wchar_t*NameOfLogFile,const int NumMeasures, double *parrMeasures)
{
	FILE *fr ;	if ((fr = _wfopen(NameOfLogFile,L"r"))== NULL)	{
	 //String St =  NameOfLogFile ;
	 //ShowMessage(L"TYrRead::ReadMeasuresInKirnosLog\n ERROR ! Not possible to open " +St) ;
	 return -1 ;
	}	char str[1000];
	char *pch = NULL;
	char *pch0 = NULL;
	bool bMeasure = true;
	int numCur = 0;
	long int iCurPosBegin = 0;
	long int iCurPosEnd = 0;
	 int i = 0 ;
	 for (i = 0; i < 10000; i++)
	 {
	   iCurPosBegin =  ftell(fr);
	   pch =	fgets(str,1001,fr);
	   if (!pch )break;
	   if(!strstr(pch, "airlow")) continue;

	   // проверка vel
	   char *pch0 =  strstr(pch, "vel");
	   pch0 = pch0 + 3;
	  int ivel0 = atoi(pch0);

	   for (int j=0; j < 3; j++)
	   {
		pch =	fgets(str,1001,fr);
		bMeasure = true;
		if(!strstr(pch, "airlow"))
		{
			bMeasure = false;
			break;
		}
		else
		{
		 char *pch1 =  strstr(pch, "vel");
		 pch1 = pch1 + 3;
		 int ivel1 = atoi(pch1);
		 if (ivel1 != ivel0)

		 {
           bMeasure = false;
		   break;
		 }
		}
	   }
	   if (!bMeasure)
	   {
		continue;
	   }
	  iCurPosEnd =  ftell(fr);
	  fseek(fr,(iCurPosBegin -  iCurPosEnd)* sizeof(char), SEEK_CUR);

	   for (int j=0; j< 4; j++)
	   {
		   pch =	fgets(str,1001,fr);
		   if (!pch )
		   {
			//ShowMessage(L"TYrRead::ReadMeasuresInKirnosLog\n ERROR 0!" ) ;
		   break;
		   }
		   double arrMeasCur[2] = {0.};
		   int numRow = 0;
		   ExtractMeasure (str, &numRow,arrMeasCur);
		   parrMeasures[8 *  numCur + 2 * numRow    ] = arrMeasCur[0];
		   parrMeasures[8 *  numCur + 2 * numRow + 1] = arrMeasCur[1];

	   }

	   numCur++;
	   if (numCur == NumMeasures)
	   {
		break;
	   }
	 }
    fclose(fr);
	return 0 ;
}



void  TYrRead::ExtractMeasure (char *str, int *numRow, double *arrMeasCur)
{
 char *pstr0 = strstr(str, "num") ;
 char *pstr1 = strstr(str, "Re") ;
 char *pstr2 = strstr(str, " Im") ;
 char *pstr3 = strchr(str, ']') ;

 pstr0 = pstr0 + 4;
 pstr0[1] =0;
 *numRow = atoi(pstr0);

  pstr1 = pstr1 + 3;
  pstr2[0] =0;
  arrMeasCur[0] = atof (pstr1);

  pstr2 = pstr2 + 3;
  pstr3[0] = 0;
  arrMeasCur[1] = atof (pstr2);



}

void  TYrRead::CleaneMeasArrFromZer0(int * quantMeas,double *arrMeas)
{
  bool bReady = true;
  int nC = * quantMeas;
  for (int i =0; i < nC; i++)
  {
	  bReady = true;
	 for (int j=0;  j < *quantMeas; j++)
	 {
	   double *pnter = &arrMeas [8*j] ;
	   for (int k = 0; k <4; k++)
	   {
		 if (sqrt(pnter [2 * k]* pnter [2 * k] + pnter [2 * k+1]* pnter [2 * k+1])< 0.0001)
		 {
			memcpy(&arrMeas [8*j], &arrMeas [8*j +8], (*quantMeas -1 -j) * 8 * sizeof(double));
			*quantMeas -= 1;
			bReady = false;
			break;
		 }
	   }
	   if( !bReady )
	   {
		   break;
	   }
	 }
	 if( bReady )
	   {
		   break;
	   }
  }
}

int TYrRead::ReadDataFromKaurovTXT(wchar_t*NameOfDataFile, int *piNumRows, double *parrData)
 {
		FILE *fr ;
	if ((fr = _wfopen(NameOfDataFile,L"r"))== NULL)	{
	 //String St =  NameOfDataFile;
	 //ShowMessage(L"TYrRead::ReadMeasuresInKirnosLog\n ERROR ! Not possible to open " +St) ;
	 return -1 ;
	}	char str[1000];	*piNumRows = 0;
	 char *pch;
	 int i = 0 ;
	 for (i = 0; i < 10000; i++)
	 {
	   pch =	fgets(str,1001,fr);
	   if (!pch )break;
	   replace(str) ;
	   char *pstr1 = str;
	   parrData[ 7 * i] = 	atof (pstr1);
	   for (int j =0; j < 6; j++)
	   {
		   pstr1 = strchr(pstr1, '\t') ;
		   for (int k =0; k < 20; k++)
		   {
			 if (pstr1 [0] == '\t')
			 {
			   pstr1++;
			   continue;
			 }
			 break;
		   }
		   parrData[ 7 * i +1 +j] = 	atof (pstr1);
		 }
		 (*piNumRows)++  ;

	   }


	fclose(fr);
	return 0 ;
 }

int TYrRead::calcColCountFromKaurovTXT(wchar_t*NameOfDataFile)
{
  	FILE *fr ;
	if ((fr = _wfopen(NameOfDataFile,L"r"))== NULL)	{
	 //String St =  NameOfDataFile;
	 //ShowMessage(L"TYrRead::ReadMeasuresInKirnosLog\n ERROR ! Not possible to open " +St) ;
	 return -1 ;
	}	char str[1000];
	 char *pch;
	 pch =	fgets(str,1001,fr);
	if (!pch)
	{
	return -1;
	}
	int iColCount = 0;
	if (str[0] != '\t')
	{
	  iColCount = 1;
	}
	char *pstr1 = str;
	int iReturn = -1;
	for (int i = 0; i < 10000; i++)
	{
	 pstr1 = strchr(pstr1, '\t');
	 if(! pstr1)
	 {
	   iReturn = iColCount;
	   break;
	  // return iColCount;
	 }
	 else
	 {
		iColCount++;
		for (int j =0; j < 20; j++)
		{
		if (pstr1[0] /*[j]*/ == '\t')
		{
		pstr1++;
		continue;
		}
		break;
		}
	 }
    }
    fclose(fr);
   return iReturn;
}



int TYrRead::ReadDataFromKaurovTXT_(wchar_t*NameOfDataFile,const int NUmCols, int *piNumRows, double *parrData)
 {
		FILE *fr ;
	if ((fr = _wfopen(NameOfDataFile,L"r"))== NULL)	{
	 //String St =  NameOfDataFile;
	 //ShowMessage(L"TYrRead::ReadMeasuresInKirnosLog\n ERROR ! Not possible to open " +St) ;
	 return -1 ;
	}	char str[1000];	*piNumRows = 0;
	 char *pch;
	 int i = 0 ;
	 for (i = 0; i < 10000; i++)
	 {
	   pch =	fgets(str,1001,fr);
	   if (!pch )break;
	   replace(str) ;
	   char *pstr1 = str;
	   parrData[ NUmCols * i] = 	atof (pstr1);
	   for (int j =0; j < (NUmCols-1); j++)
	   {
		   pstr1 = strchr(pstr1, '\t') ;
		   for (int k =0; k < 20; k++)
		   {
			 if (pstr1 [0] == '\t')
			 {
			   pstr1++;
			   continue;
			 }
			 break;
		   }
		   parrData[ NUmCols * i +1 +j] = 	atof (pstr1);
		 }


		 (*piNumRows)++  ;
       }

	fclose(fr);
	return 0 ;
 }

void  TYrRead::flip(double *parrData, int len)
{
	for (int i = 0; i < len /2; i++)
	{
	   double temp = parrData[i];
	   parrData[i] = parrData[len -1 -i] ;
	   parrData[len -1 -i] =  temp;
	}
}

int TYrRead::calcColCountFromCSV(wchar_t*NameOfDataFile)
{
    FILE *fr ;
    if ((fr = _wfopen(NameOfDataFile,L"r"))== NULL)	{
     //String St =  NameOfDataFile;
     //ShowMessage(L"TYrRead::ReadMeasuresInKirnosLog\n ERROR ! Not possible to open " +St) ;
     return -1 ;
    }	char str[1000];
     char *pch;
     pch =	fgets(str,1001,fr);
    if (!pch)
    {
    return -1;
    }
    int iColCount = 0;

    char *pstr1 = str;
    for (int i = 0; i < 10000; i++)
    {
     pstr1 = strchr(pstr1, ';');
     if(! pstr1)
     {
       break;
     }
     else
     {
        iColCount++;
        pstr1 = &(pstr1[1]);
     }
    }
    fclose(fr);
    return iColCount+1 ;
}

//--------------------------------------------------------------
int TYrRead::ReadHdrFromGearFileForDriver(wchar_t*NameOfSCVFile, bool *bVelo
            , double *pvalL, double *pvalR,  double *pvalPsi_f, double *pvalJ0
            , double *pvalMresidual,double *pvalJLoad, double *pvalCx_om, double *pvalOm0
             , double*pval_a , double*pval_c  , double*pval_l )
{
    int len = wcslen(NameOfSCVFile) ;
    if ( !( (NameOfSCVFile[len - 1] == 'v') && (NameOfSCVFile[len - 2] == 's') // проверка, что
     && (NameOfSCVFile[len - 3] == 'c') ) )  // указанный файл имеет расширение .flt
    {
     //String St = NameOfFltFile ;
      //ShowMessage(L"TYrRead::ReadHdrFileFromFltFile\nExtention of " + St + L"is wrong") ;
      return -1 ;
    }

    int quanTitleCols =4;

    double parrMass[12] = {0.};
    FILE *fr;
   if((fr=_wfopen(NameOfSCVFile,L"r")) ==  NULL)
   {
       //ShowMessage(L"ERROR IN YrReadTabCSV\n Not possible to open FileName") ;
       //Abort() ;
   }
    char str[1000];
    fgets(str,1001,fr);
    if (fgets(str,1001,fr))
    {
        if(str[0] == 0) return -1 ;
        int lenline = strlen(str);
         str[lenline] = ';';
         str[lenline + 1] = 0 ;
         lenline++;
        for (int i =0; i < lenline; i++)
        {
          if (str[i] == ',')
          {
              str[i] = '.';
          }

        }
        char str1[50] ={0};
        int i_temp = 0;
        int ncols_temp = 0;
        int jcur = 0;
        char *pch = str;
        for (int i =0; i < quanTitleCols; i++)
        {
         pch = strchr(pch,';');
         pch++;
        }

        char *pch1 =pch;
        for (int i =0; i < 12; i++)
        {
            pch1 = strchr(pch,';');
            pch1[0] = 0;
            parrMass[i] = atof(pch);
          pch = pch1 +1;

        }
    }
    fclose(fr);
    *bVelo = (parrMass[0]< 0.5)?false:true;
    *pvalL = parrMass[1];
    *pvalR = parrMass[2] ;
    *pvalPsi_f = parrMass[3] ;
    *pvalJ0 = parrMass[4] ;
    *pvalMresidual= parrMass[5] ;
    *pvalJLoad = parrMass[6] ;
    *pvalCx_om = parrMass[7] ;
    *pvalOm0 = parrMass[8] ;
    *pval_a = parrMass[9] ;
    *pval_c = parrMass[10] ;
    *pval_l = parrMass[11] ;
    return 0 ;
}

//------------------------------------------
int TYrRead::ReadHdrFromGearFileForDriver(wchar_t*NameOfSCVFile, bool *bVelo
            , double *pvalL, double *pvalR,  double *pvalPsi_f, double *pvalJ0
            , double *pvalMresidual,double *pvalJLoad, double *pvalCx_om, double *pvalOm0)
{
    int len = wcslen(NameOfSCVFile) ;
    if ( !( (NameOfSCVFile[len - 1] == 'v') && (NameOfSCVFile[len - 2] == 's') // проверка, что
     && (NameOfSCVFile[len - 3] == 'c') ) )  // указанный файл имеет расширение .flt
    {
     //String St = NameOfFltFile ;
      //ShowMessage(L"TYrRead::ReadHdrFileFromFltFile\nExtention of " + St + L"is wrong") ;
      return -1 ;
    }

    int quanTitleCols =4;

    double parrMass[9] = {0.};
    FILE *fr;
   if((fr=_wfopen(NameOfSCVFile,L"r")) ==  NULL)
   {
       //ShowMessage(L"ERROR IN YrReadTabCSV\n Not possible to open FileName") ;
       //Abort() ;
   }
    char str[1000];
    fgets(str,1001,fr);
    if (fgets(str,1001,fr))
    {
        if(str[0] == 0) return -1 ;
        int lenline = strlen(str);
         str[lenline] = ';';
         str[lenline + 1] = 0 ;
         lenline++;
        for (int i =0; i < lenline; i++)
        {
          if (str[i] == ',')
          {
              str[i] = '.';
          }

        }
        char str1[50] ={0};
        int i_temp = 0;
        int ncols_temp = 0;
        int jcur = 0;
        char *pch = str;
        for (int i =0; i < quanTitleCols; i++)
        {
         pch = strchr(pch,';');
         pch++;
        }

        char *pch1 =pch;
        for (int i =0; i < 9; i++)
        {
            pch1 = strchr(pch,';');
            pch1[0] = 0;
            parrMass[i] = atof(pch);
          pch = pch1 +1;

        }
    }
    fclose(fr);
    *bVelo = (parrMass[0]< 0.5)?false:true;
    *pvalL = parrMass[1];
    *pvalR = parrMass[2] ;
    *pvalPsi_f = parrMass[3] ;
    *pvalJ0 = parrMass[4] ;
    *pvalMresidual= parrMass[5] ;
    *pvalJLoad = parrMass[6] ;
    *pvalCx_om = parrMass[7] ;
    *pvalOm0 = parrMass[8] ;
    return 0 ;
}
//--------------------------------------------
int TYrRead::ReadOrientationInputFile(wchar_t *FileName, int *pQuantBaseFreq
      , double *arrBaseFreqData, int *pQuantWorkFreq,double *arrWorkFreq
     , double* arrCavityAntenna, int *piPredeterminationType, double *arrPredeterminatedAngs)
{
    FILE *fr;
    fr=_wfopen(FileName,L"r");
    char str[1000];

    //1
    fgets(str,1001,fr);
    if(str[0] == 0) return 1;
    char *pch = str;
    pch = strchr(pch,';')+1;
    char *pch1 = pch;
    pch1 = strchr(pch1,';');
    pch1[0] = 0;

    *pQuantBaseFreq = atoi(pch);

    //2
    fgets(str,1001,fr);
    if(str[0] == 0) return 2;
    pch = str;
    pch = strchr(pch,';')+1;
    pch1 = pch;
    pch1 = strchr(pch1,';');
    pch1[0] = 0;

    *pQuantWorkFreq = atoi(pch);

    //3
    int nrows1 = 0;
    int ncols = 3;
    fgets(str,1001,fr);
    for (int i = 0; i < *pQuantBaseFreq; ++i)
    {
        fgets(str,1001,fr);
        if(str[0] == 0) break;
        int lenline = strlen(str);
         str[lenline] = ';';
         str[lenline + 1] = 0 ;
         lenline++;
        for (int i =0; i < lenline; i++)
        {
          if (str[i] == ',') str[i] = '.';


        }
        char str1[50] ={0};
        int i_temp = 0;
        int ncols_temp = 0;
        for (int i =0; i < lenline; i++)
        {
           if (str[i]!=';')
           {
             str1[i_temp] = str[i];
              i_temp++;
           }
           else
           {
               arrBaseFreqData[nrows1*ncols + ncols_temp] = atof(str1);

               ncols_temp++;
                i_temp = 0;
                for (int j = 0; j < 50; j++) str1[j] =0;

           }

        }
       nrows1++;
    }


    //4
    int nrows2 = 0;
    int ncols2 = 1;
    fgets(str,1001,fr);
    for (int i = 0; i < *pQuantWorkFreq; ++i)
    {
        fgets(str,1001,fr);
        if(str[0] == 0) break;
        int lenline = strlen(str);
         str[lenline] = ';';
         str[lenline + 1] = 0 ;
         lenline++;
        for (int i =0; i < lenline; i++)
        {
          if (str[i] == ',') str[i] = '.';


        }
        char str1[50] ={0};
        int i_temp = 0;
        int ncols_temp = 0;
        for (int i =0; i < lenline; i++)
        {
           if (str[i]!=';')
           {
             str1[i_temp] = str[i];
              i_temp++;
           }
           else
           {
               arrWorkFreq[nrows2*ncols2 + ncols_temp] = atof(str1);

               ncols_temp++;
                i_temp = 0;
                for (int j = 0; j < 50; j++) str1[j] =0;

           }

        }
       nrows2++;
    }

    //5

    int ncols3 = 5;
    int nrows3 = 0;
    fgets(str,1001,fr);
    fgets(str,1001,fr);
    fgets(str,1001,fr);
    for (int i = 0; i < 1; ++i)
    {
        fgets(str,1001,fr);
        if(str[0] == 0) break;
        int lenline = strlen(str);
         str[lenline] = ';';
         str[lenline + 1] = 0 ;
         lenline++;
        for (int i =0; i < lenline; i++)
        {
          if (str[i] == ',') str[i] = '.';


        }
        char str1[50] ={0};
        int i_temp = 0;
        int ncols_temp = 0;
        for (int i =0; i < lenline; i++)
        {
           if (str[i]!=';')
           {
             str1[i_temp] = str[i];
              i_temp++;
           }
           else
           {
               arrCavityAntenna[nrows3*ncols3 + ncols_temp] = atof(str1);

               ncols_temp++;
                i_temp = 0;
                for (int j = 0; j < 50; j++) str1[j] =0;

           }

        }
       nrows3++;
    }

    //6
    fgets(str,1001,fr);
    fgets(str,1001,fr);
    if(str[0] == 0) return 1;
    pch = str;
    pch = strchr(pch,';');
    pch[0] = 0;

    int iNum = atoi(str);
    if (iNum ==0)
    {
      *piPredeterminationType = 0;
        return 0;
    }
    if (iNum ==2)
    {
      *piPredeterminationType = 3;
      fgets(str,1001,fr);
      if (strstr(pch,"mest"))
      {
          fgets(str,1001,fr);
          if(str[0] == 0) return 1;
          pch = str;
          pch = strchr(pch,';');
          pch[0] = 0;
          arrPredeterminatedAngs[0] = atof(str);
          fgets(str,1001,fr);
          fgets(str,1001,fr);
          if(str[0] == 0) return 1;
          pch = str;
          pch = strchr(pch,';');

          pch[0] = 0;
          arrPredeterminatedAngs[1] = atof(str);
          return 0;

      }
      else {
          fgets(str,1001,fr);
          if(str[0] == 0) return 1;
          pch = str;
          pch = strchr(pch,';');
          pch[0] = 0;
          arrPredeterminatedAngs[1] = atof(str);
          fgets(str,1001,fr);
          if(str[0] == 0) return 1;
          pch = str;
          pch = strchr(pch,';');

          pch[0] = 0;
          arrPredeterminatedAngs[0] = atof(str);
          return 0;
      }
      return 0;
    }
    if (iNum ==1)
    {
        fgets(str,1001,fr);
        if (strstr(pch,"mest"))
        {
            fgets(str,1001,fr);
            if(str[0] == 0) return 1;
            pch = str;
            pch = strchr(pch,';');

            pch[0] = 0;
            arrPredeterminatedAngs[0] = atof(str);

            arrPredeterminatedAngs[1] = 0.;
            *piPredeterminationType = 1;
            return 0;
        }
        else
        {
            fgets(str,1001,fr);
            if(str[0] == 0) return 1;
            pch = str;
            pch = strchr(pch,';');

            pch[0] = 0;
            arrPredeterminatedAngs[1] = atof(str);

            arrPredeterminatedAngs[0] = 0.;
            *piPredeterminationType = 2;
            return 0;
        }
    }


  return 0;
}

















































