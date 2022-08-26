#include "Platform.h"

#include <string.h>
#include <wchar.h>
#include <stdio.h>
#include <time.h>

#include "Environment.h"

QPlatform::QPlatform()
{
 mpHidroRLS = NULL;
 memset(marrPosParams, 0, 6 * sizeof(double));
}

// конструктор копирования
  QPlatform :: QPlatform (const  QPlatform &R)
 {
   mpHidroRLS = R.mpHidroRLS;
   memcpy( marrPosParams, R.marrPosParams, 6 * sizeof(double));
 }

 // оператор присваивания
 QPlatform  &QPlatform::operator=( const QPlatform  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      mpHidroRLS = R.mpHidroRLS;
      memcpy( marrPosParams, R.marrPosParams, 6 * sizeof(double));

     return *this ;
 }

  // парам конструктор
 QPlatform:: QPlatform (QHidroRLS *pHidroRLS, const double *arrPosParams)

 {
    mpHidroRLS = pHidroRLS;
    memcpy( marrPosParams, arrPosParams, 6 * sizeof(double));
 }

 //
 //----------------------------


 int QPlatform::createInputDataReport(wchar_t*FileName, const bool bHeader)
 {
     int len = wcslen(FileName) ;

     if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'x') // проверка, что
      && (FileName[len - 3] == 't') ) )  // указанный файл имеет расширение .flt
     {
       return 1 ;
     }

     FILE *fw ;

     if ((fw = _wfopen(FileName,L"a"))== NULL)

     {
      return 1 ;
     }
 if (bHeader)
 {
    fprintf(fw,"***************************************\n");
    fprintf(fw,"  Дата и время формирования отчета\n");
    time_t t = time(NULL);
    struct tm* aTm = localtime(&t);
    fprintf(fw,"  Год = %04d\n",aTm->tm_year+1900);
    fprintf(fw,"  Mесяц = %02d\n",aTm->tm_mon+1);
    fprintf(fw,"  День = %02d\n",aTm->tm_mday);
    fprintf(fw,"  Время = %02d:%02d:%02d\n",aTm->tm_hour, aTm->tm_min, aTm->tm_sec);
    fprintf(fw,"***************************************\n");
    fprintf(fw,"***************************************\n");
    fprintf(fw,"***************************************\n");
 }
 fprintf(fw,"***************************************\n");
 fprintf(fw,"      Платформа   \n");
 fprintf(fw,"   1.параметры позиционирования в ПСК  \n");
 fprintf(fw,"  вектор параллакса (м):\n");
 fprintf(fw,"  X = %5.3f; Y  =  %5.3f; Z =  %5.3f;\n",marrPosParams[0], marrPosParams[1], marrPosParams[2]);
 fprintf(fw,"  углы (рад):\n");
 fprintf(fw,"  Betta = %5.4f; Eps  =  %5.4f; Alf =  %5.4f;\n",marrPosParams[3], marrPosParams[4], marrPosParams[5]);
 fprintf(fw,"        *Примечание*\n");
 fprintf(fw,"  *Betta - поворот относительно оси OZ ПСК по часовой стрелке*\n");
 fprintf(fw,"  *Eps - поворот относительно оси OX ПСК против часовой стрелки*\n");
 fprintf(fw,"  *Alf - поворот относительно оси OY ПСК по часовой стрелке*\n");

 fprintf(fw,"***************************************\n");
 fprintf(fw,"    2.устройство \n");


 fclose(fw);
 //mHidroRLS->createInputDataReport(FileName, false);
 createInheritedInputDataReport(FileName, false);

 }
 //----------------------------------
 int QPlatform::createInheritedInputDataReport(wchar_t*FileName, const bool bHeader)
 {
     int len = wcslen(FileName) ;

     if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'x') // проверка, что
      && (FileName[len - 3] == 't') ) )  // указанный файл имеет расширение .flt
     {
       return 1 ;
     }

     FILE *fw ;

     if ((fw = _wfopen(FileName,L"a"))== NULL)

     {
      return 1 ;
     }
 if (bHeader)
 {
    fprintf(fw,"***************************************\n");
    fprintf(fw,"  Дата и время формирования отчета\n");
    time_t t = time(NULL);
    struct tm* aTm = localtime(&t);
    fprintf(fw,"  Год = %04d\n",aTm->tm_year+1900);
    fprintf(fw,"  Mесяц = %02d\n",aTm->tm_mon+1);
    fprintf(fw,"  День = %02d\n",aTm->tm_mday);
    fprintf(fw,"  Время = %02d:%02d:%02d\n",aTm->tm_hour, aTm->tm_min, aTm->tm_sec);
    fprintf(fw,"***************************************\n");
    fprintf(fw,"***************************************\n");
    fprintf(fw,"***************************************\n");
 }
  fprintf(fw,"***************************************\n");
 fprintf(fw,"      3. тип платформы - стационарная   \n");

  fclose(fw);
 }




