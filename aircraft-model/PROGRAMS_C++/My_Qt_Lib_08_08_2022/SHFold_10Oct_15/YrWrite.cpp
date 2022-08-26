﻿
#include "YrWrite.h"
  #include <stdio.h>
#include <string.h>
#include <wchar.h>
//---------------------------------------------------------------------------

//#pragma package(smart_init)
extern const double NODATA;
extern const double NODATA1 ;
extern const double PI;
extern const double DNIL  ;
extern const double DNIL1  ;
extern const double EPSSING   ;
extern const double EPS   ;

 int TYrWrite::PutPointsToCsvFile(wchar_t*FileName,float * pX, float * pY,
							 float* pZ,const int lenArray,int * lenVars)
// Запись массива точек в файл .csv Десятичный разделитель - запятая, табуляция ;
{
   FILE  *fw ;
		int num = 0 ;
	 fw=_wfopen(FileName,L"w");
	 if(!fw) //ShowMessage (L"TYrWrite::PutPointsToCsvFile\nFile is not opened !") ;
	 fprintf(fw,"Object-ID;");
	 fprintf(fw,"X;");
	 fprintf(fw,"Y;");
	 fprintf(fw,"Z;\n");

	for (int i = 0; i < lenArray; i++)
	{
	   if (pZ[ i ] > NODATA)
	   {
		   char ch[100] ;

			sprintf(ch,"%d; %f; %f; %f;\n",num,pX[i], pY[i],pZ[i]) ;
			int lenNumberStr = strlen(ch) ;
			for (int n = 0 ; n < lenNumberStr ; n++)
			{
			   if (ch [n] == '.') ch [n] = ',';
			   
		   
			}
			fprintf(fw,"%s",ch);
			num++;

	   }
	} 

	 fclose(fw);
       *lenVars  = num ;
	return 0 ;
}

int TYrWrite::PutPointsToCsvFile(wchar_t*FileName,double* pX, double * pY,
							double* pZ,const int lenArray,int * lenVars)
// Запись массива точек в файл .csv Десятичный разделитель - запятая, табуляция ;
// полимарфная для double
{
   FILE  *fw ;
		int num = 0 ;
	 fw=_wfopen(FileName,L"w");
	 if(!fw) //ShowMessage (L"TYrWrite::PutPointsToCsvFile\nFile is not opened !") ;
	 fprintf(fw,"Object-ID;");
	 fprintf(fw,"X;");
	 fprintf(fw,"Y;");
	 fprintf(fw,"Z;\n");

	for (int i = 0; i < lenArray; i++)
	{
	   if (pZ[ i ] > NODATA)
	   {
		   char ch[100] ;

			sprintf(ch,"%d; %f; %f; %f;\n",num,pX[i], pY[i],pZ[i]) ;
			int lenNumberStr = strlen(ch) ;
			for (int n = 0 ; n < lenNumberStr ; n++)
			{
			   if (ch [n] == '.') ch [n] = ',';


			}
			fprintf(fw,"%s",ch);
			num++;

	   }
	}

	 fclose(fw);
       *lenVars  = num ;
	return 0 ;
}



int TYrWrite::PutPointsToTxtFile(wchar_t*FileName,float* pX, float * pY,
							 float* pZ,const int lenArray,int * lenVars)
// Запись массива в txt файл формата Point
// n x y z  для float
{



FILE  *fw ;

	 fw=_wfopen(FileName,L"w");
	 if(!fw)
	 {
	  //ShowMessage (L"TYrWrite::PutPointsToTxtFile\nFile is not opened !") ;
	  return  1 ;
	 }
	 fprintf(fw,"Point\n") ;
	  char ch[100] ;
	  int lenCh; //длина ch
	   int num = 0 ;
	for (int i = 0; i < lenArray; i++)
	{
	   if (pZ[ i ] > NODATA)
	   {


			 int num1 = num +1 ;
			 int ibegin,iend;
		sprintf(ch,"%d %f %f %f\n",num1,pX[i], pY[i],pZ[i]) ;
		lenCh = strlen(ch) ;
			for (int n = 0 ; n < lenCh ; n++)
			{
			   if (ch [n] == '.') ch [n] = ',';


			}
			fprintf(fw,"%s",ch);

			num++;

	   }
	}
	fprintf(fw,"END") ;

	 fclose(fw);
       *lenVars  = num ;
	return 0 ;

}

int TYrWrite::PutPointsToTxtFile(wchar_t*FileName,double* pX, double * pY,
							 double* pZ,const int lenArray,int * lenVars)
// Запись массива в txt файл формата Point
// n x y z  для double  ПОЛИМОРФНАЯ
{



FILE  *fw ;

	 fw=_wfopen(FileName,L"w");
	 if(!fw)
	 {
	  //ShowMessage (L"TYrWrite::PutPointsToTxtFile\nFile is not opened !") ;
	  return  1 ;
	 }
	 fprintf(fw,"Point\n") ;
	  char ch[100] ;
	  int lenCh; //длина ch
	   int num = 0 ;
	for (int i = 0; i < lenArray; i++)
	{
	   if (pZ[ i ] > NODATA)
	   {


			 int num1 = num +1 ;
			 int ibegin,iend;
		sprintf(ch,"%d %f %f %f\n",num1,pX[i], pY[i],pZ[i]) ;
		lenCh = strlen(ch) ;
			for (int n = 0 ; n < lenCh ; n++)
			{
			   if (ch [n] == '.') ch [n] = ',';


			}
			fprintf(fw,"%s",ch);

			num++;

	   }
	}
	fprintf(fw,"END") ;

	 fclose(fw);
       *lenVars  = num ;
	return 0 ;

}
int TYrWrite::PutSetOfMassivesToPointsTxtFiles(const wchar_t*FolderName,const wchar_t *TypesFileNames[]
								 ,const int lenSet,float* pX
						   , float * pY,float* pTotalZ,const int lenArray)
 // ДЛЯ FLOAT
// Запись набора массивов в текстовые файлы точек
// FolderName - путь к папке с файлами. Например:
// const wchar_t TxtFilesNotPreobrMVPoints[]  = L"C:\\ZemProekt\\DebugMV\\TxtFilesNotPreobrMVPoints\\NotPreobr_";
// NotPreobr_ - это общее начало в именах файлов
// TypesFileNames[] - это массив строк окончаний файлов. Например:
// const wchar_t * GeomFilesName//String[11] = {L"Z.",L"Zm.",L"p.",L"q.",L"pm.",L"qm."
//						,L"rm.",L"sm.",L"tm.",L"MCA.",L"MDA."};
//  Таким образом в папке TxtFilesNotPreobrMVPoints откроются
// 11 файлов  NotPreobr_Point_Z.txt
//, NotPreobr_Point_Zm.txt, NotPreobr_Point_p.txt,....,NotPreobr_Point_MDA.txt
//  lenSet - к-во записываемых файлов  (= 11)
// pX,pY - массивы координат X,Y ( длина lenArray)
// pTotalZ - массив с величинами  Z  ( длина lenArray * lenSet)
{
   wchar_t FileName[300] ;
   wchar_t * strType = L"Point_" ;
   wchar_t* ExtStr = L"txt" ;

   for (int i = 0; i < lenSet; i++)
   {
   wcscpy(FileName,FolderName) ;
   wcscat(FileName,strType) ;
   wcscat(FileName,TypesFileNames[i]) ;
   wcscat(FileName,ExtStr) ;
   int lenVars ;

	 TYrWrite::PutPointsToTxtFile(FileName,pX, pY,
							 &pTotalZ[i*lenArray], lenArray,&lenVars);

   }



	return 0 ;

}


int TYrWrite::PutSetOfMassivesToPointsTxtFiles(const wchar_t*FolderName,const wchar_t *TypesFileNames[]
								 ,const int lenSet,double* pX
						   , double * pY,double* pTotalZ,const int lenArray)
 // ДЛЯ DOUBLE   ПОЛИМОРФНАЯ
// Запись набора массивов в текстовые файлы точек
// FolderName - путь к папке с файлами. Например:
// const wchar_t TxtFilesNotPreobrMVPoints[]  = L"C:\\ZemProekt\\DebugMV\\TxtFilesNotPreobrMVPoints\\NotPreobr_";
// NotPreobr_ - это общее начало в именах файлов
// TypesFileNames[] - это массив строк окончаний файлов. Например:
// const wchar_t * GeomFilesName//String[11] = {L"Z.",L"Zm.",L"p.",L"q.",L"pm.",L"qm."
//						,L"rm.",L"sm.",L"tm.",L"MCA.",L"MDA."};
//  Таким образом в папке TxtFilesNotPreobrMVPoints откроются
// 11 файлов  NotPreobr_Point_Z.txt
//, NotPreobr_Point_Zm.txt, NotPreobr_Point_p.txt,....,NotPreobr_Point_MDA.txt
//  lenSet - к-во записываемых файлов  (= 11)
// pX,pY - массивы координат X,Y ( длина lenArray)
// pTotalZ - массив с величинами  Z  ( длина lenArray * lenSet)
{
   wchar_t FileName[300] ;
   wchar_t * strType = L"Point_" ;
   wchar_t* ExtStr = L"txt" ;
   for (int i = 0; i < lenSet; i++)
   {
   wcscpy(FileName,FolderName) ;
   wcscat(FileName,strType) ;
   wcscat(FileName,TypesFileNames[i]) ;
   wcscat(FileName,ExtStr) ;
   int lenVars ;

	 TYrWrite::PutPointsToTxtFile(FileName,pX, pY,
							 &pTotalZ[i*lenArray], lenArray,&lenVars);

   }



	return 0 ;

}

int TYrWrite::WriteMassiveInFltFile(wchar_t*FileName,float * parrZ, const int nrows,
							 const int ncols,const float xllcorner ,const float yllcorner,
							 const float cellsize,const float NODATA_value )
// Для FLOAT массивов
// Функция записывает массив высот parrZ  в файл .flt и формирует файл .hdr
// в формате FLT ESRI. FileName - имя файла с расширением .flt. nrows и  ncols - к-во строк и столбцов
// соответвенно.  xllcorner и yllcorner   - координата x и y левого нижнего угла,
// cellsize- шаг сетки,  NODATA_value - значение, которым заполняется parrZ,
//  если величина Z  не определена
// Массив  parrZ хранится в порядке слева направо и снизу вверх. То есть parrZ[0]
// это нижняя левая точка, а     parrZ[nrows * ncols ] это правая верхняя точка
// В ArcMap считается, что каждая вершина (точка) лежит в центре квадрата со стороной cellsize
// Таким обраэзом, в начале координат находится какбы не сама точка а нижний левый угол
// ее квадрата. Т е начало системы координат оказывается сдвинутым относительно левой нижней
// точки. Координаты левой нижней точки равны (cellsize/2,cellsize/2)
{
// 1. Формирование .hdr
	int len = wcslen(FileName) ;
	if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'l') // проверка, что
	 && (FileName[len - 3] == 'f') ) )  // указанный файл имеет расширение .flt
	{
	 //String St =  FileName ;
	  //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nExtention of " + St + L" is wrong") ;
	  return 1 ;
	}
	wchar_t* HdrFileName = new wchar_t[len] ;
	wcscpy( HdrFileName, FileName) ;
	 HdrFileName[len - 1] = L'r';
	 HdrFileName[len - 2] = L'd';
	 HdrFileName[len - 3] = L'h';
	 
	FILE *fw ;


	if ((fw = _wfopen(HdrFileName,L"w"))== NULL)
	{
	 //String St =  HdrFileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}
   fprintf(fw,"ncols         %i\n",ncols);
   fprintf(fw,"nrows         %i\n",nrows);
   fprintf(fw,"xllcorner     %f\n",xllcorner);
   fprintf(fw,"yllcorner     %f\n",yllcorner);
   fprintf(fw,"cellsize      %f\n",cellsize);
   fprintf(fw,"NODATA_value  %f\n",NODATA_value);
   fprintf(fw,"byteorder     LSBFIRST");
   delete [] HdrFileName ;
   fclose(fw);
 // 2. Запись .flt
   	FILE *fw1 ;
   if ((fw1 = _wfopen(FileName,L"wb"))== NULL)
	{
	 //String St =  FileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\n ERROR ! Not possible to open " +St) ;
	 return 1 ;
	}
 
	 float *Ztemp = new float [ nrows * ncols] ;
	if (Ztemp == NULL)
	{
	   //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nThere are not memory for Ztemp in TYrRead::ReadFltFile") ;
	   return 1 ;
	}

	for (int i = 0; i < nrows; i++)
	for (int j =0; j < ncols; j++)
	{

	 Ztemp[ ncols * i +j]= parrZ[(nrows - 1 -i)*ncols + j]   ;
	 

	}
  fwrite(Ztemp,sizeof(float),ncols * nrows,fw1) ;
		   //	TYrWrite::WriteReportForFloatMassiveTXT(L"C:\\zemproekt\\DebugMV\\New\\ZtempYrWrite.txt"
							//	 , Ztemp,4,3) ;
	  //	TYrWrite::WriteReportForFloatMassiveTXT(L"C:\\zemproekt\\DebugMV\\New\\parrZYrWrite.txt"
							  //	 , parrZ,4,3) ;
  delete [] Ztemp ;
  fclose(fw1);

  return 0 ;

}



int TYrWrite::WriteSetOfMassivesInIkfFile(const wchar_t*FileName,float * parrZ,const int quanMass
							 ,const int nrows,const int ncols,const float xllcorner
							 ,const float yllcorner, const float cellsize,const float nodat)

  {
		int len = wcslen(FileName) ;
	if ( !( (FileName[len - 1] == 'f') && (FileName[len - 2] == 'k') // проверка, что
	 && (FileName[len - 3] == 'i') ) )  // указанный файл имеет расширение .flt
	{
	 //String St =  FileName ;
	//  //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nExtention of " + St + L" is wrong") ;
	  return 1 ;
	}
	wchar_t* HdrFileName = new wchar_t[len] ;
	wcscpy( HdrFileName, FileName) ;
	 HdrFileName[len - 1] = 'h';

	FILE *fw ;


	if ((fw = _wfopen(HdrFileName,L"w"))== NULL)
	{
	 //String St =  HdrFileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}
   fprintf(fw,"ncols         %i\n",ncols);
   fprintf(fw,"nrows         %i\n",nrows);
   fprintf(fw,"xllcorner     %f\n",xllcorner);
   fprintf(fw,"yllcorner     %f\n",yllcorner);
   fprintf(fw,"cellsize      %f\n",cellsize);
   fprintf(fw,"NODATA_value  %f\n",nodat);
   fprintf(fw,"nmassives     %i\n",quanMass);

   delete [] HdrFileName ;
   fclose(fw);
 // 2. Запись .flt
	FILE *fw1 ;
   if ((fw1 = _wfopen(FileName,L"wb"))== NULL)
	{
	 //String St =  FileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\n ERROR ! Not possible to open " +St) ;
	 return 1 ;
	}

  fwrite(parrZ,sizeof(float),ncols * nrows *quanMass,fw1) ;


  fclose(fw1);

  return 0 ;
  }

  int TYrWrite::WriteHdrForIkdFile(const wchar_t*FileName,const int quanMass
							 ,const int nrows,const int ncols,const double xllcorner
							 ,const double yllcorner, const double cellsize,const double nodat)
// запись HDR файла для IKD файла
 {
     int len = wcslen(FileName) ;
	if ( !( (FileName[len - 1] == 'd') && (FileName[len - 2] == 'k') // проверка, что
	 && (FileName[len - 3] == 'i') ) )  // указанный файл имеет расширение .ikd
	{
	 //String St =  FileName ;
	  //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nExtention of !!!! " + St + L" is wrong") ;
	  return 1 ;
	}
	wchar_t* HdrFileName = new wchar_t[len] ;
	wcscpy( HdrFileName, FileName) ;
	 HdrFileName[len - 1] = 'h';

	FILE *fw ;


	if ((fw = _wfopen(HdrFileName,L"w"))== NULL)
	{
	 //String St =  HdrFileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}
   fprintf(fw,"ncols         %i\n",ncols);
   fprintf(fw,"nrows         %i\n",nrows);
   fprintf(fw,"xllcorner     %f\n",xllcorner);
   fprintf(fw,"yllcorner     %f\n",yllcorner);
   fprintf(fw,"cellsize      %f\n",cellsize);
   fprintf(fw,"NODATA_value  %f\n",nodat);
   fprintf(fw,"nmassives     %i\n",quanMass);

   delete [] HdrFileName ;
   fclose(fw);
   return 0 ;
 }

 int TYrWrite::WriteMassiveInRastrTextFile(wchar_t*FileName,float * parrZ, const int nrows,
							 const int ncols,const float xllcorner ,const float yllcorner,
							 const float cellsize,const float NODATA_value )
 // ДЛЯ МАССИВОВ FLOAT
// Запись массива parrZ в формат растра txt
//  nrows и  ncols - к-во строк и столбцов массива
// xllcorner и  yllcorner - координаты левого нижнего угла левого нижнего пиксела
// cellsize - размер пиксела(ячейки)
// NODATA_value - значение которым заполнен элемент массива  parrZ в случае
//               отсутствия информации
// Если имеется массив высот parrZ, которые заданы в точках pX и pY,
// то  xllcorner = pX[0] - cellsize/2;  yllcorner = pY[0] - cellsize/2;
// если нумерация массива начинается снизу, слева и идет слева направо, снизу вверх
// то есть в порядке возрастания географических координат в Северном полушарии
{
	// 1. Формирование .hdr
	int len = wcslen(FileName) ;
	if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'x') // проверка, что
	 && (FileName[len - 3] == 't') ) )  // указанный файл имеет расширение .flt
	{
	 //String St =  FileName ;
	  //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nExtention of " + St + L" is wrong") ;
	  return 1 ;
	}
	
	char ch[300] ;

	FILE *fw ;


	if ((fw = _wfopen(FileName,L"w"))== NULL)
	{
	 //String St =  FileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}
   fprintf(fw,"ncols         %i\n",ncols);
   fprintf(fw,"nrows         %i\n",nrows);
   fprintf(fw,"xllcorner     %f\n",xllcorner);
   fprintf(fw,"yllcorner     %f\n",yllcorner);
   fprintf(fw,"cellsize      %f\n",cellsize);
   fprintf(fw,"NODATA_value  %f\n",NODATA_value);

  // 2.Формирование массива высот
		   int  inumber = 0;
		   float temp = 0;
		   for (int i=0 ; i< nrows; i++)
		   {
			 for (int j = 0; j < ncols; j++)
			 {
				temp =  parrZ[(nrows -1 - i) * ncols + j] ;

				if(temp > (NODATA_value + 1)) inumber++;


				if (j != (ncols-1))
				{
				  sprintf(ch,"%f ",temp) ;
				}
				else
				{
				sprintf(ch,"%f\n",temp) ;
				}

			   int	lenCh = strlen(ch) ;
					for (int n = 0 ; n < lenCh ; n++)
					{
					   if (ch [n] == '.') ch [n] = ',';
					}
					fprintf(fw,"%s",ch);
			 }

		   }


			fclose(fw);



		 return inumber;


	}
 int TYrWrite::WriteMassiveInRastrTextFile(wchar_t*FileName,double * parrZ, const int nrows,
							 const int ncols,const double xllcorner ,const double yllcorner,
							 const double cellsize,const double NODATA_value )
 // ДЛЯ МАССИВОВ DOUBLE ПОЛИМОРФНАЯ
// Запись массива parrZ в формат растра txt
//  nrows и  ncols - к-во строк и столбцов массива
// xllcorner и  yllcorner - координаты левого нижнего угла левого нижнего пиксела
// cellsize - размер пиксела(ячейки)
// NODATA_value - значение которым заполнен элемент массива  parrZ в случае
//               отсутствия информации
// Если имеется массив высот parrZ, которые заданы в точках pX и pY,
// то  xllcorner = pX[0] - cellsize/2;  yllcorner = pY[0] - cellsize/2;
// если нумерация массива начинается снизу, слева и идет слева направо, снизу вверх
// то есть в порядке возрастания географических координат в Северном полушарии
{
	// 1. Формирование .hdr
	int len = wcslen(FileName) ;
	if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'x') // проверка, что
	 && (FileName[len - 3] == 't') ) )  // указанный файл имеет расширение .flt
	{
	 //String St =  FileName ;
	  //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nExtention of " + St + L" is wrong") ;
	  return 1 ;
	}

	char ch[300] ;

	FILE *fw ;


	if ((fw = _wfopen(FileName,L"w"))== NULL)
	{
	 //String St =  FileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}
   fprintf(fw,"ncols         %i\n",ncols);
   fprintf(fw,"nrows         %i\n",nrows);
   fprintf(fw,"xllcorner     %f\n",(float)xllcorner);
   fprintf(fw,"yllcorner     %f\n",(float)yllcorner);
   fprintf(fw,"cellsize      %f\n",(float)cellsize);
   fprintf(fw,"NODATA_value  %f\n",(float)NODATA_value);

  // 2.Формирование массива высот
		   int  inumber = 0;
		   float temp = 0;
		   for (int i=0 ; i< nrows; i++)
		   {
			 for (int j = 0; j < ncols; j++)
			 {
				temp = (float) parrZ[(nrows -1 - i) * ncols + j] ;

				if(temp > (NODATA_value + 1)) inumber++;


				if (j != (ncols-1))
				{
				  sprintf(ch,"%f ",temp) ;
				}
				else
				{
				sprintf(ch,"%f\n",temp) ;
				}

			   int	lenCh = strlen(ch) ;
					for (int n = 0 ; n < lenCh ; n++)
					{
					   if (ch [n] == '.') ch [n] = ',';
					}
					fprintf(fw,"%s",ch);
			 }

		   }


			fclose(fw);



		 return inumber;


	}

 int TYrWrite::WriteReportForFloatMassiveTXT(const wchar_t*FileName,float*parr,
										  const int ncols,const int nrows)
 // запись 2 мерного массива типа float в текстовый файл по строкам
 // разделитель между числами - ;
 // знак десятичного разделителя - . (десятичная точка)
 {

 char ch[300] ;

	FILE *fw ;


	if ((fw = _wfopen(FileName,L"w"))== NULL)
	{
	 //String St =  FileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}



			 int  inumber = 0;

		   for (int i=0 ; i< nrows; i++)
		   {
			 for (int j = 0; j < ncols; j++)
			 {



						inumber++;

				if (j != (ncols-1))
				{
				  sprintf(ch,"%f; ",parr[ i * ncols + j]) ;
				}
				else
				{
				sprintf(ch,"%f;\n",parr[ i * ncols + j]) ;
				}

			   /*int	lenCh = strlen(ch) ;
					for (int n = 0 ; n < lenCh ; n++)
					{
					   if (ch [n] == '.') ch [n] = ',';
					}  */
					 
						fprintf(fw,"%s",ch);
			 }

		   }


			fclose(fw);



		 return inumber;
 }

 int TYrWrite::WriteReportForFloatMassiveTXT(const wchar_t*FileName,double*parr,
										  const int ncols,const int nrows)
// ПОЛИМОРФНАЯ ДЛЯ DOUBLE
 // запись 2 мерного массива типа double в текстовый файл по строкам
 // разделитель между числами - ;
 // знак десятичного разделителя - . (десятичная точка)
 {

 char ch[300] ;

	FILE *fw ;


	if ((fw = _wfopen(FileName,L"w"))== NULL)
	{
	 //String St =  FileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}



			 int  inumber = 0;

		   for (int i=0 ; i< nrows; i++)
		   {
			 for (int j = 0; j < ncols; j++)
			 {



						inumber++;

				if (j != (ncols-1))
				{
				  sprintf(ch,"%f; ",parr[ i * ncols + j]) ;
				}
				else
				{
				sprintf(ch,"%f;\n",parr[ i * ncols + j]) ;
				}

			   /*int	lenCh = strlen(ch) ;
					for (int n = 0 ; n < lenCh ; n++)
					{
					   if (ch [n] == '.') ch [n] = ',';
					}  */

						fprintf(fw,"%s",ch);
			 }

		   }


			fclose(fw);



		 return inumber;
 }
  int WriteReportForFloatMassiveTXT(const wchar_t*FileName,bool *parr,
										  const int ncols,const int nrows)
// ПОЛИМОРФНАЯ ДЛЯ BOOL
 // запись 2 мерного массива типа double в текстовый файл по строкам
 // разделитель между числами - ;
 // знак десятичного разделителя - . (десятичная точка)
 {

 char ch[300] ;

	FILE *fw ;


	if ((fw = _wfopen(FileName,L"w"))== NULL)
	{
	 //String St =  FileName ;
   //	cout<<"WriteMassiveInFltFile\nERROR ! Not possible to open FileName"<<endl ;
	 return 1 ;
	}



			 int  inumber = 0;

		   for (int i=0 ; i< nrows; i++)
		   {
			 for (int j = 0; j < ncols; j++)
			 {



						inumber++;

				if (j != (ncols-1))
				{
				  sprintf(ch,"%i; ",parr[ i * ncols + j]) ;
				}
				else
				{
				sprintf(ch,"%i;\n",parr[ i * ncols + j]) ;
				}

						fprintf(fw,"%s",ch);
			 }

		   }


			fclose(fw);



		 return inumber;
 }

 int WriteReportForFloatMassiveTXT(const wchar_t*FileName,int *piarr,
										  const int ncols,const int nrows)
// ПОЛИМОРФНАЯ ДЛЯ INT
 // запись 2 мерного массива типа double в текстовый файл по строкам
 // разделитель между числами - ;
 // знак десятичного разделителя - . (десятичная точка)
 {

 char ch[300] ;

	FILE *fw ;


	if ((fw = _wfopen(FileName,L"w"))== NULL)
	{
	 //String St =  FileName ;
	//cout<<"WriteMassiveInFltFile\nERROR ! Not possible to open FileName"<<endl ;
	 return 1 ;
	}



			 int  inumber = 0;

		   for (int i=0 ; i< nrows; i++)
		   {
			 for (int j = 0; j < ncols; j++)
			 {



						inumber++;

				if (j != (ncols-1))
				{
				  sprintf(ch,"%i; ",piarr[ i * ncols + j]) ;
				}
				else
				{
				sprintf(ch,"%i;\n",piarr[ i * ncols + j]) ;
				}

						fprintf(fw,"%s",ch);
			 }

		   }


			fclose(fw);



		 return inumber;
 }


int WriteReportForFloatMassiveTXT(const wchar_t*FileName,long *piarr,
										  const int ncols,const int nrows)
// ПОЛИМОРФНАЯ ДЛЯ LONG
 // запись 2 мерного массива типа double в текстовый файл по строкам
 // разделитель между числами - ;
 // знак десятичного разделителя - . (десятичная точка)
 {

 char ch[300] ;

	FILE *fw ;


	if ((fw = _wfopen(FileName,L"w"))== NULL)
	{
	 //String St =  FileName ;
   //	cout<<"WriteMassiveInFltFile\nERROR ! Not possible to open FileName"<<endl ;
	 return 1 ;
	}



			 int  inumber = 0;

		   for (int i=0 ; i< nrows; i++)
		   {
			 for (int j = 0; j < ncols; j++)
			 {



						inumber++;

				if (j != (ncols-1))
				{
				  sprintf(ch,"%i; ",piarr[ i * ncols + j]) ;
				}
				else
				{
				sprintf(ch,"%i;\n",piarr[ i * ncols + j]) ;
				}

						fprintf(fw,"%s",ch);
			 }

		   }


			fclose(fw);



		 return inumber;
 }

  int TYrWrite::WriteSetOfMassivesInFltFiles(const wchar_t*FolderName
			 ,const wchar_t *TypesFileNames[] ,const int lenSet
						 ,float* pTotalZ, const int nrows, const int ncols
						  ,const float xllcorner ,const float yllcorner,
							 const float cellsize,const float nodata )
 // Записывет  lenSet массивов типа float во FLT файлы

  {

   wchar_t FileName[300] ;
  // wchar_t * strType = L"flt_" ;
   wchar_t* ExtStr = L"flt" ;



   for (int i = 0; i < lenSet; i++)
   {
   wcscpy(FileName,FolderName) ;
  // wcscat(FileName,strType) ;
   wcscat(FileName,TypesFileNames[i]) ;
   wcscat(FileName,ExtStr) ;

	 TYrWrite::WriteMassiveInFltFile(FileName,&pTotalZ[nrows*ncols*i],  nrows,
							  ncols,xllcorner, yllcorner,cellsize,  nodata );

   }



	return 0 ;

}

int TYrWrite::WriteSetOfMassivesInFltFiles(const wchar_t*FolderName
			 ,const wchar_t *TypesFileNames[] ,const int lenSet
						 ,double* pTotalZ, const int nrows, const int ncols
						  ,const double xllcorner ,const double yllcorner,
							 const double cellsize,const double nodata )
 // Записывет  lenSet массивов типа double во FLT файлы
 // ПОЛИМОРФНАЯ для double

  {
   wchar_t FileName[300] ;
  // wchar_t * strType = L"flt_" ;
   wchar_t* ExtStr = L"flt" ;

        //String	s_22 = L"TYrWrite::WriteSetOfMassivesInFltFiles,FolderName = " ;
  //	//ShowMessage(s_22 + FolderName) ;
//	s_22 = L"lenset =" ;
	  //	//ShowMessage(s_22 + lenSet + L"nrows = " + nrows + L" ncols = " + ncols) ;
   for (int i = 0; i < lenSet; i++)
   {
    // s_22 = L"i= ";
   //   //ShowMessage(s_22 + i) ;
   wcscpy(FileName,FolderName) ;
  // wcscat(FileName,strType) ;
   wcscat(FileName,TypesFileNames[i]) ;
   wcscat(FileName,ExtStr) ;
      //   s_22 = L"TYrWrite::WriteSetOfMassivesInFltFiles,FileName = " ;
	 // //ShowMessage(s_22 + FileName) ;
	 // //ShowMessage(L"Input TYrWrite::WriteSetOfMassivesInFltFiles") ;

	 TYrWrite::WriteMassiveInFltFile(FileName,&pTotalZ[nrows*ncols*i],  nrows,
							  ncols,xllcorner, yllcorner,cellsize,  nodata );

   }



	return 0 ;

}

 int TYrWrite::WriteMassiveInFltFile(wchar_t*FileName,double * parrZ, const int nrows,
							 const int ncols,const double xllcorner ,const double yllcorner,
							 const double cellsize,const double NODATA_value )
// ДЛЯ DOUBLE МАССИВОВ ПОЛИМОРФНАЯ
// Функция записывает массив высот parrZ  в файл .flt и формирует файл .hdr
// в формате FLT ESRI. FileName - имя файла с расширением .flt. nrows и  ncols - к-во строк и столбцов
// соответвенно.  xllcorner и yllcorner   - координата x и y левого нижнего угла,
// cellsize- шаг сетки,  NODATA_value - значение, которым заполняется parrZ,
//  если величина Z  не определена
// Массив  parrZ хранится в порядке слева направо и снизу вверх. То есть parrZ[0]
// это нижняя левая точка, а     parrZ[nrows * ncols ] это правая верхняя точка
// В ArcMap считается, что каждая вершина (точка) лежит в центре квадрата со стороной cellsize
// Таким обраэзом, в начале координат находится какбы не сама точка а нижний левый угол
// ее квадрата. Т е начало системы координат оказывается сдвинутым относительно левой нижней
// точки. Координаты центра левого нижнего  пиксела равны (cellsize/2,cellsize/2)
{

	 // //String s_22 = L"TYrWrite::WriteMassiveInFltFile Input,FileName = " ;
	 // //ShowMessage(s_22 + FileName) ;
	 //  s_22 = L" nrows,ncols,xllcorner, yllcorner,cellsize,NODATA_value= "  ;
	  //	 //ShowMessage(s_22 +  nrows+ L" ; " +ncols+ L" ; " +xllcorner+ L" ; " + yllcorner+ L" ; " +cellsize+ L" ; " +NODATA_value ) ;
// 1. Формирование .hdr
	int len = wcslen(FileName) ;
	if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'l') // проверка, что
	 && (FileName[len - 3] == 'f') ) )  // указанный файл имеет расширение .flt
	{
	 //String St =  FileName ;
	  //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nExtention of " + St + L" is wrong") ;
	  return 1 ;
	}
    wchar_t  HdrFileName [500] = {0};
	wcscpy( HdrFileName, FileName) ;
	 HdrFileName[len - 1] = L'r';
	 HdrFileName[len - 2] = L'd';
	 HdrFileName[len - 3] = L'h';
	//  s_22 = L"TYrWrite::WriteMassiveInFltFile Input,HdrFileName = " ;
	//  //ShowMessage(s_22 +HdrFileName) ;
	FILE *fw ;


	if ((fw = _wfopen(HdrFileName,L"w"))== NULL)
	{
	 //String St =  HdrFileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}
   fprintf(fw,"ncols         %i\n",ncols);
   fprintf(fw,"nrows         %i\n",nrows);
   fprintf(fw,"xllcorner     %f\n",xllcorner);
   fprintf(fw,"yllcorner     %f\n",yllcorner);
   fprintf(fw,"cellsize      %f\n",cellsize);
   fprintf(fw,"NODATA_value  %f\n",NODATA_value);
   fprintf(fw,"byteorder     LSBFIRST");

   fclose(fw);
 // 2. Запись .flt
   	FILE *fw1 ;
   if ((fw1 = _wfopen(FileName,L"wb"))== NULL)
	{
	 //String St =  FileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\n ERROR ! Not possible to open " +St) ;
	 return 1 ;
	}

	 float *Ztemp = new float [ nrows * ncols] ;
	if (Ztemp == NULL)
	{
	   //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nThere are not memory for Ztemp in TYrRead::ReadFltFile") ;
	   return 1 ;
	}

	for (int i = 0; i < nrows; i++)
	for (int j =0; j < ncols; j++)
	{

	 Ztemp[ ncols * i +j]=(float) parrZ[(nrows - 1 -i)*ncols + j]   ;


	}
  fwrite(Ztemp,sizeof(float),ncols * nrows,fw1) ;
		   //	TYrWrite::WriteReportForFloatMassiveTXT(L"C:\\zemproekt\\DebugMV\\New\\ZtempYrWrite.txt"
							//	 , Ztemp,4,3) ;
	  //	TYrWrite::WriteReportForFloatMassiveTXT(L"C:\\zemproekt\\DebugMV\\New\\parrZYrWrite.txt"
							  //	 , parrZ,4,3) ;
  delete [] Ztemp ;
  fclose(fw1);

  return 0 ;

}

 int TYrWrite::WriteReportForIntMassiveTXT(const wchar_t*FileName,int*parr,
										  const int ncols,const int nrows)
// ПОЛИМОРФНАЯ ДЛЯ DOUBLE
 // запись 2 мерного массива типа double в текстовый файл по строкам
 // разделитель между числами - ;
 // знак десятичного разделителя - . (десятичная точка)
 {

 char ch[300] ;

	FILE *fw ;


	if ((fw = _wfopen(FileName,L"w"))== NULL)
	{
	 //String St =  FileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}



			 int  inumber = 0;

		   for (int i=0 ; i< nrows; i++)
		   {
			 for (int j = 0; j < ncols; j++)
			 {



						inumber++;

				if (j != (ncols-1))
				{
				  sprintf(ch,"%d; ",parr[ i * ncols + j]) ;
				}
				else
				{
				sprintf(ch,"%d;\n",parr[ i * ncols + j]) ;
				}

			   /*int	lenCh = strlen(ch) ;
					for (int n = 0 ; n < lenCh ; n++)
					{
					   if (ch [n] == '.') ch [n] = ',';
					}  */

						fprintf(fw,"%s",ch);
			 }

		   }


			fclose(fw);



		 return inumber;
 }
 int TYrWrite::WriteReportForIntMassiveTXT(const wchar_t*FileName,long*plarr,
										  const int ncols,const int nrows)
// ПОЛИМОРФНАЯ ДЛЯ LONG
 // запись 2 мерного массива типа long в текстовый файл по строкам
 // разделитель между числами - ;
 // знак десятичного разделителя - . (десятичная точка)
 {

 char ch[300] ;

	FILE *fw ;


	if ((fw = _wfopen(FileName,L"w"))== NULL)
	{
	 //String St =  FileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}



			 int  inumber = 0;

		   for (int i=0 ; i< nrows; i++)
		   {
			 for (int j = 0; j < ncols; j++)
			 {



						inumber++;

				if (j != (ncols-1))
				{
				  sprintf(ch,"%d; ",plarr[ i * ncols + j]) ;
				}
				else
				{
				sprintf(ch,"%d;\n",plarr[ i * ncols + j]) ;
				}

			   /*int	lenCh = strlen(ch) ;
					for (int n = 0 ; n < lenCh ; n++)
					{
					   if (ch [n] == '.') ch [n] = ',';
					}  */

						fprintf(fw,"%s",ch);
			 }

		   }


			fclose(fw);



		 return inumber;
 }



 // запись 2 мерного массива типа double в текстовый файл подряд
 // разделитель между числами - ,
 // знак десятичного разделителя - . (десятичная точка)
  int TYrWrite::WriteReportForFloatMassiveTXT_(const wchar_t*FileName,double*parr,
										  const int ncols,const int nrows)

 {

 char ch[300] ;

	FILE *fw ;


	if ((fw = _wfopen(FileName,L"w"))== NULL)
	{
	 //String St =  FileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}



			 int  inumber = 0;
			int iChetchik = 0;
		   for (int i=0 ; i< nrows * ncols; i++)
		   {
			 //for (int j = 0; j < ncols; j++)
			 //{



				iChetchik++;

			 /*	if (j != (ncols-1))
				{
				  sprintf(ch,"%f, ",parr[ i * ncols + j]) ;
				}
				else
				{
				sprintf(ch,"%f;\n",parr[ i * ncols + j]) ;
				} */

						if (i  < nrows * ncols)
						{
						 if (iChetchik == 100)
						 {
						   sprintf(ch,"%f, \n",parr[ i ]) ;
						   iChetchik = 0;
						 }
						 else
						 {
						 sprintf(ch,"%f, ",parr[ i ]) ;
						 }
						}
						else
						{
						  sprintf(ch,"%f ",parr[ i ]) ;
                        }
						fprintf(fw,"%s",ch);
			 //}

			 }


			fclose(fw);
		 return inumber;
 }

// запись массива информации в  CSV файл
// INPUT:
// FileName
// parrBuff[ iNumRows * iNumCols] - массив с информацией
// iNumRows- к-во строк массива
// iNumCols - к-во столбцов массива
// pwcharrRowNames[ iLenName* iNumRows] - имена строк массива
// pwcharrColNames [ iLenName* iNumCols] - имена сьолбцов массива
int TYrWrite::WriteMassiveInFIleSCV(wchar_t*FileName,double *parrBuff, int iNumRows, int iNumCols
							 ,wchar_t *pwcharrRowNames,wchar_t *pwcharrColNames, int iLenName)
// Запись массива точек в файл .csv Десятичный разделитель - запятая, табуляция ;
{
	  FILE  *fw ;
	 fw=_wfopen(FileName,L"w");
    // if(!fw) //ShowMessage (L"TYrWrite::PutPointsToCsvFile\nFile is not opened !") ;

	 wchar_t wchTemp[500];
	 memset(wchTemp, 0, 500 * sizeof(wchar_t));
	 if(pwcharrRowNames)
	 {
	 wcscat(  wchTemp, L" ;");
	 fwprintf(fw,wchTemp);
	 }

	 if(pwcharrColNames)
	 {
	 for (int i = 0; i < iNumCols; i++)
	 {
	   wcscpy(wchTemp, &pwcharrColNames[ i *  iLenName]);
	   wcscat(  wchTemp, L";");
	   fwprintf(fw,wchTemp);
	 }
	 fwprintf(fw,L"\n");
     }


	 for (int i = 0; i <iNumRows; i++)
	 {
	   memset(wchTemp, 0, 500 * sizeof(wchar_t));
	   if(pwcharrRowNames)
	   {
	   wcscpy(wchTemp, &pwcharrRowNames[ i *  iLenName]);
	   wcscat(  wchTemp, L";");
       }
	   for (int j =0; j < iNumCols; j++)
	   {
		 wchar_t wch[100] ={0};
		 swprintf(wch, L"%f", parrBuff[ i *  iNumCols + j]);
         if (j != (iNumCols - 1))
         {
		 wcscat(  wch, L";");
         }
		 wcscat(  wchTemp,wch);
	   }
	   wcscat(  wchTemp,L"\n");
	   fwprintf(fw,wchTemp);
	 }
	 fclose(fw);


	return 0 ;
}



 int TYrWrite::WriteReportForIntMassiveTXT_Comma(const wchar_t*FileName,double *parr,
										  const int ncols,const int nrows)
// ПОЛИМОРФНАЯ ДЛЯ DOUBLE
 // запись 2 мерного массива типа double в текстовый файл по строкам
 // разделитель между числами - ,
 // знак десятичного разделителя - . (десятичная точка)
 {

 char ch[300] ;

	FILE *fw ;


	if ((fw = _wfopen(FileName,L"w"))== NULL)
	{
	 //String St =  FileName ;
	 //ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}




		   int  inumber = 0;

		   for (int i=0 ; i< nrows; i++)
		   {
			 for (int j = 0; j < ncols; j++)
			 {
				inumber++;

				if (j != (ncols-1))
				{
				  sprintf(ch,"%f, ",parr[ i * ncols + j]) ;
				}
				else
				{
				sprintf(ch,"%f,\n",parr[ i * ncols + j]) ;
				}


						fprintf(fw,"%s",ch);
			 }

		   }


			fclose(fw);



		 return inumber;
 }

int TYrWrite::WriteMassiveWithHeaderInFIleSCV(wchar_t*FileName,double *parrBuff, int iNumRows, int iNumCols
                              ,wchar_t *pwcharrHeader)
 {
    FILE  *fw ;
   fw=_wfopen(FileName,L"w");
   fwprintf(fw,pwcharrHeader);
   fwprintf(fw,L"\n");
   wchar_t wchTemp[500];
   memset(wchTemp, 0, 500 * sizeof(wchar_t));
   for (int i = 0; i <iNumRows; i++)
   {
     memset(wchTemp, 0, 500 * sizeof(wchar_t));

     for (int j =0; j < iNumCols; j++)
     {
       wchar_t wch[100] ={0};
       swprintf(wch, L"%f", parrBuff[ i *  iNumCols + j]);
       if (j != (iNumCols - 1))
       {
       wcscat(  wch, L";");
       }
       wcscat(  wchTemp,wch);
     }
     wcscat(  wchTemp,L"\n");
     fwprintf(fw,wchTemp);
   }
   fclose(fw);


  return 0 ;
 }






















