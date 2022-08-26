#include "SubWaterBeam.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <QFile>
#include <QString>
#include "Table_1D.h"
#include "MatrixProccess.h"
#include <QVector>
#include <QPointF>
#include "Table_1D.h"
#include "CoordSystTrsf.h"
#include "BigMeasure.h"
#include "PeaceVess.h"

QSubWaterBeam::QSubWaterBeam()
{

}

int QSubWaterBeam::calcColCountFrom000(wchar_t*NameOf000File)
{
    FILE *fr ;
    if ((fr = _wfopen(NameOf000File,L"r"))== NULL)	{
     //String St =  NameOfDataFile;
     //ShowMessage(L"TYrRead::ReadMeasuresInKirnosLog\n ERROR ! Not possible to open " +St) ;
     return -1 ;
    }

    char str[1000];
    for(int i =0; i < 10000; ++i)
    {
        if(!fgets(str,1001,fr))
        {
            return 0;
        }
        if(strstr(str, "VELOCITY"))
        {
            break;
        }

    }
    int iColCount=0;
    for(int i =0; i < 10000; ++i)
    {
        if(!fgets(str,1001,fr))
        {
            break;
        }
        ++iColCount;

    }
    return iColCount;

}

//-------------------------------

int QSubWaterBeam::ReadDataFrom000(wchar_t*NameOfDataFile, int *piNumRows, double *parrData)
 {
    FILE *fr ;
    if ((fr = _wfopen(NameOfDataFile,L"r"))== NULL)	{
     //String St =  NameOfDataFile;
     //ShowMessage(L"TYrRead::ReadMeasuresInKirnosLog\n ERROR ! Not possible to open " +St) ;
     return -1 ;
    }

    char str[1000];
    *piNumRows = 0;
    for(int i =0; i < 10000; ++i)
    {
        if(!fgets(str,1001,fr))
        {
            return 0;
        }
        if(strstr(str, "VELOCITY"))
        {
            break;
        }

    }
    //!

    for(int i =0; i < 100000; ++i)
    {
        if(!fgets(str,1001,fr))
        {
            break;
        }
        else
        {
            (*piNumRows)++;
          char *p = str;
          p = strchr(str,'\t');
          p++;
          parrData[2 * i +1]= atof(p);
          p++;
          p = strchr(p,'\t');
          parrData[2 * i]= atof(p);

     }
    }
    return 0;

 }
//-------------------------------

int QSubWaterBeam::ReadDataFrom000_11_03_2022(wchar_t*NameOfDataFile, int *piNumRows, double *parrData)
 {
    FILE *fr ;
    if ((fr = _wfopen(NameOfDataFile,L"r"))== NULL)	{
     //String St =  NameOfDataFile;
     //ShowMessage(L"TYrRead::ReadMeasuresInKirnosLog\n ERROR ! Not possible to open " +St) ;
     return -1 ;
    }

    char str[1000];
    *piNumRows = 0;
    for(int i =0; i < 10000; ++i)
    {
        if(!fgets(str,1001,fr))
        {
            return 0;
        }
        if(strstr(str, "VELOCITY"))
        {
            break;
        }

    }
    //!

    for(int i =0; i < 100000; ++i)
    {
        if(!fgets(str,1001,fr))
        {
            break;
        }
        else
        {
            (*piNumRows)++;
          char *p = str;
          parrData[2 * i ]= atof(p);
          p = strchr(str,'\t');
          p++;
          p = strchr(p,'\t');
          parrData[2 * i +1]= atof(p);




     }
    }
    return 0;

 }

//-------------------------------

bool QSubWaterBeam::replace(wchar_t*str)
{	 int n = wcslen(str);
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

bool QSubWaterBeam::replace(char*str)
{	 int n = strlen(str);
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

//-------------------------------
void QSubWaterBeam::sortProfile(double *arrProfile, const int quantRows)
{
   qsort(arrProfile, quantRows, 2 * sizeof(double),fncCmp1);
}
//----------------------------------------------
//коррекция профиля скорости звука
// на вход подается отсортированный по глубине массив профиля звука
//что делает:
// 1. убрать узловые точки с расстоянием по глубине менее 1 мм
// 2. объединить соседние отрезки глубины с близким градиентом скорости
// 3. среди узловых точек с глубинами менее 6 м оставить единственную
//    глубина равна минимальной, а скорость равна средней скорости среди таких точек
void QSubWaterBeam::adjustProfile(double *parrData, const int iNumRows, int &iNumRowsNew)
{
    QVector <QPointF> vectPrfl(iNumRows);
    for (int i =0; i< iNumRows; ++i)
    {
      vectPrfl.replace(i,QPointF(parrData[i * 2], parrData[i * 2 + 1])) ;
    }
    // убрать отрицательные высоты
    for (int i =0; i< iNumRows; ++i)
    {
        bool bStop = true;
        for (int j = 0; j < (vectPrfl.length() ); ++j )
        {
            if (vectPrfl.at(j).x() <= 0.5)
            {
               vectPrfl.remove(j );
               bStop = false;
               continue;
            }

        }
      if (bStop)
      {
          break;
      }
    }
    //!

    // убрать нулевые скорости
    for (int i =0; i< iNumRows; ++i)
    {
        bool bStop = true;
        for (int j = 0; j < (vectPrfl.length() ); ++j )
        {
            if (vectPrfl.at(j).y() <= 100.)
            {
               vectPrfl.remove(j );
               bStop = false;
               continue;
            }

        }
      if (bStop)
      {
          break;
      }
    }
    //!


    //!// убрать узловые точки с расстоянием по глубине менее 1 мм
    for (int i =0; i< iNumRows; ++i)
    {
        bool bStop = true;
        for (int j = 0; j < (vectPrfl.length() -1); ++j )
        {
            if (fabs(vectPrfl.at(j).x() - vectPrfl.at(j +1).x()) <= 0.001)
            {
               vectPrfl.remove(j + 1);
               bStop = false;
               continue;
            }

        }
      if (bStop)
      {
          break;
      }
    }

    //!

    // объединить соседние отрезки глубины с близким градиентом скорости
    int iN0 = vectPrfl.length();
    for (int i =0; i< iN0; ++i)
    {
        bool bStop = true;
        for (int j = 0; j < (vectPrfl.length() -2); ++j )
        {
            double k0 = (vectPrfl.at(j).y() - vectPrfl.at(j +1).y())/(vectPrfl.at(j).x() - vectPrfl.at(j +1).x());
            double k1 = (vectPrfl.at(j+1).y() - vectPrfl.at(j +2).y())/(vectPrfl.at(j +1).x() - vectPrfl.at(j +2).x());
            if (fabs(k1 - k0) <= 0.0001)
            {
               vectPrfl.remove(j + 1);
               bStop = false;
               break;
            }

        }
      if (bStop)
      {
          break;
      }
    }
    //!


    int quant = 1;
    double sum = vectPrfl.at(0).y();
    // оставить первую точку, далее, убрать све точки до первой точки с глубиной >= valDepth м
    double valDepth =15.;

    for (int i =0; i< vectPrfl.length() -1; ++i)
    {

            if (vectPrfl.at(i).x() < valDepth)
            {
              sum += vectPrfl.at(i).y();
              quant++;

            }
            else
            {
             break;
            }

    }
    //!

    // исключаем quant точек начиная с превого
    vectPrfl.remove(1, quant);


    // изменяем 0-ю точку, делая ее средней скоростью
    double x0 = 0.;//vectPrfl.at(0).x();
    vectPrfl.replace(0, QPointF(x0, sum/ ((double)quant) ));

    iNumRowsNew = vectPrfl.length();
    memset(parrData, 0, 2 * iNumRows * sizeof(double));
    for (int i =0; i< iNumRowsNew; ++i)
    {
        parrData[i * 2] = vectPrfl.at(i).x();
        parrData[i * 2 + 1] = vectPrfl.at(i).y();
    }
}
//-------------------------------------------
//создание профиля скорости звука
bool QSubWaterBeam::createProfileTbl(wchar_t*NameOfDataFile, TTable_1D *ptbl, enumTypeOf_000 TypeOf_000)
{
    //создание профиля скорости звука

    int iNumRows = 0;
    double parrData[1000] = {0.}, parrData1[1000] = {0.};
    switch(TypeOf_000)
    {
    case VAR0:
        if(QSubWaterBeam::ReadDataFrom000(NameOfDataFile, &iNumRows, parrData)!=0)
        {
            return false;
        }
        break;

    case VAR1:
        if(QSubWaterBeam::ReadDataFrom000_11_03_2022(NameOfDataFile, &iNumRows, parrData)!=0)
        {
            return false;
        }
        break;

    default:

        break;
    }

     //1!

    // 2.сортировка по высоте
    QSubWaterBeam::sortProfile(parrData, iNumRows);
    if (iNumRows ==2) // это случай отладки
    {
       *ptbl =  TTable_1D (parrData,iNumRows);
       return true;
    }
    //2!

    //3. надо добавить самую глубокую точку
    parrData[iNumRows * 2]= parrData[(iNumRows-1) * 2] +25.;
    parrData[iNumRows * 2 + 1]= parrData[(iNumRows-1) * 2 + 1] - 0.06;
    ++iNumRows;

    int iNumRowsNew = -1;
    memcpy(parrData1, parrData,1000 * sizeof(double) );
    QSubWaterBeam::adjustProfile(parrData1,  iNumRows, iNumRowsNew);
   // 3!



    *ptbl =  TTable_1D (parrData1,iNumRowsNew);
    return true;
    // !


}
//----------------------------------------------
int fncCmp1( const void *a, const void *b)
{
   double *pa =  (double*)a;
   double *pb =  (double*)b;
   if (*pa < *pb)
   {
   return -1;
   }
   else
   {
       return 1;

   }

}

//-------------------------------------
// 1/sqrt(n*n- cos*cos)
double calcIntegral_I(const double VAlC0,const double VAlZn1
             ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta)
{
    double k = (VAlC2 -VAlC1)/(VAlZn2 -VAlZn1);
    if (fabs(k) < 0.00001)
    {
        double valn = VAlC0/VAlC1;
        double temp1= valn* valn -  VAlCosTetta * VAlCosTetta;
       double temp = 1./sqrt(temp1 );
       return (VAlZn2 -VAlZn1) *temp;
    }
    double t2 = fncSq1(VAlC0,VAlC2,VAlCosTetta);
    double t1 = fncSq1(VAlC0,VAlC1,VAlCosTetta);
    return -VAlC0/VAlCosTetta/VAlCosTetta /k
           *(t2 -t1);
}
//-------------------------------------
// n*n/pow(n*n- cos*cos,3/2)
double calcIntegral_I1(const double VAlC0,const double VAlZn1
             ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta)
{
    double k = (VAlC2 -VAlC1)/(VAlZn2 -VAlZn1);
    if (fabs(k) < 0.00001)
    {
        double valn = VAlC0/VAlC1;
       double temp = valn * valn/pow(valn* valn -  VAlCosTetta * VAlCosTetta,1.5 );
       return (VAlZn2 -VAlZn1) *temp;
    }
    double t2 = fncSq1(VAlC0,VAlC2,VAlCosTetta);
    double t1 = fncSq1(VAlC0,VAlC1,VAlCosTetta);
   return VAlC0/VAlCosTetta/VAlCosTetta * (VAlZn2 -VAlZn1)/(VAlC2 -VAlC1)
           *(1./t2 - 1./t1);
}
//-------------------------------------
// n*n/sqrt(n*n- cos*cos)
double calcIntegral_I2(const double VAlC0,const double VAlZn1
             ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta)
{
    double k = (VAlC2 -VAlC1)/(VAlZn2 -VAlZn1);
    if (fabs(k) < 0.00001)
    {
        double valn = VAlC0/VAlC1;
       double temp = valn * valn/sqrt(valn* valn -  VAlCosTetta * VAlCosTetta );
       return (VAlZn2 -VAlZn1) *temp;
    }
    double b = (VAlC2 -VAlC1)/(VAlZn2 -VAlZn1);
    double t2 = fncSq1(VAlC0,VAlC2,VAlCosTetta);
    double t1 = fncSq1(VAlC0,VAlC1,VAlCosTetta);

    double t =VAlC1/VAlC2 * (1. + t2)/(1. + t1);
   return -VAlC0/b * log(t);
}
//----------------------------------
double fncSq1(const double VAlC0,const double VAlCn,const double VAlCosTetta)
{
    double t = VAlCosTetta * VAlCn/VAlC0;
    return sqrt(1. - t * t);
}
//----------------------------------------
// вычисление горизоньтального расстояния луча от антенны
//INPUT:
//VAlza - глубина антенны
//VAlzm - глубина точки луча
//tblPrfl - профиль звука
//VAlCosTetta -косинус начального угла скольжения
//OUTPUT:
//возвращает горизонтадьное расстояние
double calcXHoriz(const double VAlza,const double VAlzm
         , TTable_1D &tblPrfl,const double VAlCosTetta)
{

   double val = calcSumIntegral_I( VAlza, VAlzm
                                  ,tblPrfl, VAlCosTetta,calcIntegral_I);
   return VAlCosTetta * val;
}
//-------------------------------------------
double calcSumIntegral_I( const double VAlza,const double VAlzm
    , TTable_1D &tblPrfl,const double VAlCosTetta
    ,double (*f)(const double VAlC0,const double VAlZn1
    ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta))
{
  const double VAlC0 = tblPrfl.calcValue(VAlza);
  const double VAlCend = tblPrfl.calcValue(VAlzm);
  int in0 = tblPrfl.getSegmentNum(VAlza);
  int in1 = tblPrfl.getSegmentNum(VAlzm);
  if (in0 == in1)
  {
      if(fabs(VAlza - VAlzm) <= 1.E-6)
      {
          return 0.;
      }
     return f( VAlC0,VAlza
               ,VAlzm,VAlC0,VAlCend, VAlCosTetta);
  }
  double sum = 0.;
  if (fabs(VAlza - tblPrfl.mparrArg[in0 +1]) > 1.E-4)
  {
    sum += f( VAlC0,VAlza
             ,tblPrfl.mparrArg[in0 +1],VAlC0,tblPrfl.mparrVal[in0 +1], VAlCosTetta);
  }
  for (int i = (in0 +1); i < in1; ++i)  {

      sum += f( VAlC0,tblPrfl.mparrArg[i]
               ,tblPrfl.mparrArg[i +1],tblPrfl.mparrVal[i],tblPrfl.mparrVal[i +1], VAlCosTetta);
  }
  if (fabs(VAlzm - tblPrfl.mparrArg[in1]) > 1.E-4)
  {
    sum += f( VAlC0,tblPrfl.mparrArg[in1]
             ,VAlzm,tblPrfl.mparrVal[in1],VAlCend, VAlCosTetta);
  }
  return sum;
}
//--------------------------------------------
// вычисление максимального значения скорости звука c(z)
// на отрезке [za;zm]
double calcMaxZ( const double VAlza,const double VAlzm, TTable_1D &tblPrfl)
{
    int in0 = tblPrfl.getSegmentNum(VAlza);
    int in1 = tblPrfl.getSegmentNum(VAlzm);
    double max = tblPrfl.calcValue(VAlza);
    for(int i = (in0+1); i<= in1;++i)
    {
        if(tblPrfl.mparrVal[i] > max)
        {
          max = tblPrfl.mparrVal[i] ;
        }
    }
    double temp = tblPrfl.calcValue(VAlzm);
    if(temp > max)
    {
      max = temp ;
    }
    return max;
}
//-----------------------------------------------
//вычисление критического начального угла скольжения Tetta
//при котором горизонтальная координата луча на заданной
// глубине будет максимальна
double calcTettaCrit( const double VAlza,const double VAlzm, TTable_1D &tblPrfl)
{
    double temp1 = tblPrfl.calcValue(VAlza);
    double temp2 = calcMaxZ( VAlza, VAlzm,tblPrfl);
    double temp = temp1/temp2;
    if (temp > 0.99999999)
    {
        temp =0.99999999;
    }
  return  acos(temp) ;
}

// вычисление критического радиуса распространения звуковой волны
//на заданной глубине
double calcXMCrit( const double VAlza,const double VAlzm, TTable_1D &tblPrfl)
{
  double tetta=  calcTettaCrit( VAlza, VAlzm,tblPrfl);
  return  calcXHoriz(VAlza, VAlzm,tblPrfl
                     ,cos(tetta +0.017)) ;
}
//-----------------------------------------------
//вычисление Tetta
//Tetta- угол скольжения луча, исходящего из антенны.
//INPUT:
//VAlza - глубина антенны
//VAlzm - глубина маяка
// VAlxm - гориз. расстояние от антенны до маяка
//tblPrfl - профиль скорости звука
//OUTPUT:
//valTetta - начальный угол скольжения
//возвращет true, если tetta существует и false
//в противном случае
bool calcTetta( const double VAlza,const double VAlzm,const double VAlxm
               , TTable_1D &tblPrfl, double &valTetta)
{
   // 1. проверка существования решения
    double valMaxX = calcXMCrit(  VAlza, VAlzm,tblPrfl);
    if(VAlxm >=valMaxX)
    {
        return false;
    }
    // 1!

    // 2. решение нелинейного уравнения
    double tetta0= atan((VAlzm- VAlza)/VAlxm);
    double delta = 0.00001;
    for(int i =0;i< 100; ++i)
    {
        const double VAlCosTetta = cos(tetta0);
        double fi = VAlxm - VAlCosTetta * calcSumIntegral_I(  VAlza,VAlzm
                ,tblPrfl, VAlCosTetta,calcIntegral_I);
        double dfi = calc_dfi_po_dtetta( VAlza, VAlzm, tetta0, tblPrfl);
        double step= fi/ dfi;
        tetta0 = tetta0 - 0.3 * step;
        if (fabs(step) < delta)
        {
          valTetta = tetta0;
          return true;
        }
    }
    return false;
}
//-------------------------------
//вычисление производной функции fi(tetta) при решении уравнения относительно tetta
//это производная горизонтальной координаты точки луча (Xs) по углу
//скольжения tetta
double calc_dfi_po_dtetta( const double VAlza,const double VAlzm,const double tetta0
               , TTable_1D &tblPrfl)
{
    return sin(tetta0)* calcSumIntegral_I(  VAlza,VAlzm
                                            ,tblPrfl, cos(tetta0),calcIntegral_I1);
}
//--------------------------------------------
//вычисление вектора градиента функции tetta
//по переменным xm, zm, za -
// горизонтальное расстояние от антенны до маяка, глубина маяка,
// глубина антенны
// INPUT:
//VAlza - глубина антенны
//VAlzm - глубина маяка
//tblPrfl - профиль скорости звука
//tetta0 - начальный угол скольжения
//OUTPUT:
//arrgrad[3] - вектор градиента функции tetta
// преполагается, что луч с нач углом вкольжения tetta0
// проходит через маяк !!!
void calc_gradTetta_po_dZ( const double VAlza,const double VAlzm
      ,const double tetta0, TTable_1D &tblPrfl, double *arrgrad)
{
   double temp =  1./calc_dfi_po_dtetta( VAlza, VAlzm, tetta0,tblPrfl);
   double val_nzs = tblPrfl.calcValue(VAlza)/tblPrfl.calcValue(VAlzm);
   double valCosTetta = cos(tetta0);
 arrgrad[0] = -temp;
 arrgrad[1] = temp * valCosTetta/sqrt(val_nzs * val_nzs - valCosTetta*valCosTetta);
 double val_dc0= tblPrfl.calc_d_po_dx(VAlza);
 const double VAlC0 = tblPrfl.calcValue(VAlza);
 arrgrad[2] = -temp *(1./ tan(tetta0) +valCosTetta *val_dc0/VAlC0
                      *calcSumIntegral_I(  VAlza, VAlzm
                          ,tblPrfl, valCosTetta
                          ,calcIntegral_I1));
}
//-----------------------------------
//вычисление времени достижения лучем маяка
// INPUT:
//VAlza - глубина антенны
//VAlzm - глубина маяка
//tblPrfl - профиль скорости звука
//tetta0 - начальный угол скольжения
double calc_t( const double VAlza,const double VAlzm
                                   ,const double tetta0, TTable_1D &tblPrfl)
{
   return calcSumIntegral_I(  VAlza,VAlzm
                              ,tblPrfl, cos(tetta0),calcIntegral_I2)/tblPrfl.calcValue(VAlza);
}
/*
//-----------------------------------
//вычисление времени достижения лучем маяка
// INPUT:
//VAlza - глубина антенны
//VAlzm - глубина маяка
//tblPrfl - профиль скорости звука
//tetta0 - начальный угол скольжения
bool calc_Tetta_and_t( const double VAlza,const double VAlzm,const double VAlxm
                       ,double &valTetta,double* arrGradTetta,double &t, TTable_1D &tblPrfl)
{
   if(! calcTetta( VAlza, VAlzm,VAlxm, tblPrfl, valTetta))
   {
           return false;
   }
   calc_gradTetta_po_dZ(  VAlza, VAlzm
         ,valTetta, tblPrfl, arrGradTetta);
    t = calcSumIntegral_I(  VAlza,VAlzm
                              ,tblPrfl, cos(valTetta),calcIntegral_I2)/tblPrfl.calcValue(VAlza);
    return true;
}
*/
//--------------------------------------------
//вычисление вектора градиента функции t - времени
//по переменным xm, zm, za -
// горизонтальное расстояние от антенны до маяка, глубина маяка,
// глубина антенны
// INPUT:
//VAlza - глубина антенны
//VAlzm - глубина маяка
//tblPrfl - профиль скорости звука
//tetta0 - начальный угол скольжения
// ARrGradTetta[3]- градиент функции tetta
//OUTPUT:
//arrgrad[3] - вектор градиента функции t
// преполагается, что луч с нач углом скольжения tetta0
// проходит через маяк !!!
/*
void calc_grad_t_po_dZ( const double VAlza,const double VAlzm
      ,const double tetta0, TTable_1D &tblPrfl,const double *ARrGradTetta, double *arrgrad)
{
    double arr_dU_po_dX[9] ={0.};
    arr_dU_po_dX[1]= arr_dU_po_dX[5] = 1.;
    memcpy(&arr_dU_po_dX[6], ARrGradTetta, 3 * sizeof(double));
    //!

    //
    double valC0 = tblPrfl.calcValue(VAlza);
    double val_nzs = valC0/tblPrfl.calcValue(VAlzm);
    double valCosTetta = cos(tetta0);
    double val_dc0 = tblPrfl.calc_d_po_dx(VAlza);
    double arr_dPsi_po_dU[3] = {0.};

    double valSumIntegral_I1 = calcSumIntegral_I(VAlza, VAlzm
                                   ,tblPrfl, valCosTetta
                                   ,calcIntegral_I1);

    // по Zs
    arr_dPsi_po_dU[0] = val_nzs * val_nzs/sqrt(val_nzs * val_nzs - valCosTetta*valCosTetta)/valC0;

    // По тетта
    arr_dPsi_po_dU[2] = -valCosTetta *sin(tetta0) /valC0 *valSumIntegral_I1;
    // по Za
    arr_dPsi_po_dU[1] = -1./(valC0 * sin(tetta0))
     -valCosTetta*valCosTetta* val_dc0/ (valC0*valC0)
      *valSumIntegral_I1;

    MtrxMultMatrx( arr_dPsi_po_dU,1, 3, arr_dU_po_dX,3, arrgrad) ;

}
*/

//--------------------------------------------
//вычисление вектора частных производных функции t - времени
//по переменным  zm, za, tetta-
//  глубина маяка,глубина антенны, угол скольжения
// INPUT:
//VAlza - глубина антенны
//VAlzm - глубина маяка
//tblPrfl - профиль скорости звука
//tetta0 - начальный угол скольжения
//OUTPUT:
//arrdt_po_db[3] - вектор градиента функции t
// преполагается, что луч с нач углом скольжения tetta0
void calc_dt_po_db( const double VAlza,const double VAlzm
      ,const double tetta0, TTable_1D &tblPrfl, double *arr_dt_po_db)
{

    //
    double valC0 = tblPrfl.calcValue(VAlza);
    double val_nzs = valC0/tblPrfl.calcValue(VAlzm);
    double valCosTetta = cos(tetta0);
    double val_dc0 = tblPrfl.calc_d_po_dx(VAlza);    

    double valSumIntegral_I1 = calcSumIntegral_I(VAlza, VAlzm
                                   ,tblPrfl, valCosTetta
                                   ,calcIntegral_I1);

    // по Zs
    arr_dt_po_db[0] = val_nzs * val_nzs/sqrt(val_nzs * val_nzs - valCosTetta*valCosTetta)/valC0;

    // По тетта
    arr_dt_po_db[2] = -valCosTetta *sin(tetta0) /valC0 *valSumIntegral_I1;
    // по Za
    arr_dt_po_db[1] = -1./(valC0 * sin(tetta0))
     -valCosTetta*valCosTetta* val_dc0/ (valC0*valC0)
      *valSumIntegral_I1;
}

//--------------------------------------------
//вычисление функции времени и  вектора градиента функции  времени
//по переменным xm, zm, za -
// горизонтальное расстояние от антенны до маяка, глубина маяка,
// глубина антенны
// INPUT:
//VAlza - глубина антенны Z[2]
//VAlzm - глубина маяка Z[1]
// VAlxm - расстояние по горизонтали от маяка до антенны Z[0]
//tblPrfl - профиль скорости звука

//OUTPUT:
// valTetta - угол скольжения
// arr_dTetta_po_dZ[3] - градиент угла скольжения
// val_t -время прохождения сигнала
//arr_dt_po_dZ[3] - вектор градиента функции t

bool calc_t_and_dt_po_dZ( const double VAlza,const double VAlzm,const double VAlxm
                          , TTable_1D &tblPrfl, double &valTetta, double *arr_dTetta_po_dZ
                          , double &val_t, double *arr_dt_po_dZ)

{
    if(!calcTetta( VAlza, VAlzm, VAlxm
                   ,tblPrfl, valTetta))
    {
        return false;
    }

    calc_gradTetta_po_dZ(  VAlza, VAlzm
          ,valTetta,tblPrfl, arr_dTetta_po_dZ);//!

    val_t = calc_t(VAlza, VAlzm,  valTetta, tblPrfl);

    double arr_db_po_dZ[9] ={0.};
    arr_db_po_dZ[1]= arr_db_po_dZ[5] = 1.;
    memcpy(&arr_db_po_dZ[6], arr_dTetta_po_dZ, 3 * sizeof(double));

   double arr_dt_po_db[3] ={0.};
    calc_dt_po_db( VAlza, VAlzm
          ,valTetta, tblPrfl, arr_dt_po_db);


   MtrxMultMatrx( arr_dt_po_db,1, 3, arr_db_po_dZ,3, arr_dt_po_dZ) ;

return true;

}
// вычисление длины дуги луча
double calc_CurveLength( const double VAlza,const double VAlzm
                                   ,const double tetta0,TTable_1D &tblPrfl)
{
    return calcSumIntegral_I(VAlza, VAlzm,tblPrfl, cos(tetta0),calcIntegral_I3);
}
//-------------------------------------------
// n/sqrt(n*n- cos*cos)
double calcIntegral_I3(const double VAlC0,const double VAlZn1
 ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta)
{
    double k = (VAlC2 -VAlC1)/(VAlZn2 -VAlZn1);
    if (fabs(k) < 0.00001)
    {
        double valn = VAlC1/VAlC0;
       double temp = valn/sqrt(valn* valn -  VAlCosTetta * VAlCosTetta );
       return (VAlZn2 -VAlZn1) *temp;
    }
   double b = (VAlC2 -VAlC1)/(VAlZn2 -VAlZn1);
   return VAlC0/(VAlCosTetta *b) * (asin(VAlCosTetta * VAlC2/VAlC0) - asin(VAlCosTetta * VAlC1/VAlC0) );
}
// вычисление текщего угла скольжения
double calc_CurTetta( const double VAlza,const double VAlzm
                                   ,const double tetta0, TTable_1D &tblPrfl)
{
    const double VAl_n = tblPrfl.calcValue(VAlza)/tblPrfl.calcValue(VAlzm);
    double valCos = cos(tetta0);
    return M_PI/2. - atan(valCos/ sqrt(VAl_n * VAl_n - valCos * valCos));
}
//------------------------------------
// Вычисление угла скольжения для варианта задания координат
// антенны и миаяка в ГСК
//INPUT:
//tblPrfl - профиль звука
//arrAntY[3] - координаты антенны в ГСК
//arrTrueBeaconPos[3]- координаты маяка в ГСК
// OUTPUT:
// valTetta - угол скольжения
// возвращает true, если задача решается
bool calcTetta_InpGSK(TTable_1D &tblPrfl, const double *arrAntY,const double *arrTrueBeaconPos
               ,  double &valTetta)
{
    //1. формирование вектора Z ф-ла (3.6)
    double arrZ[3] = {0.};
    arrZ[0] = sqrt((arrTrueBeaconPos[0] - arrAntY[0]) * (arrTrueBeaconPos[0] - arrAntY[0])
            + (arrTrueBeaconPos[1] - arrAntY[1]) * (arrTrueBeaconPos[1] - arrAntY[1]) );

    arrZ[1]= -arrTrueBeaconPos[2]; // глуб. маяка

    arrZ[2] = - arrAntY[2];  // глуб. антенны
    // 3!

    // 4. решение уравнения (3.8) (или (2.3))

    if(!calcTetta( arrZ[2],arrZ[1],arrZ[0],tblPrfl, valTetta))
    {
        return false;
    }
    return true;
}

//------------------------------------
// Вычисление времени прохождения сигнала для варианта задания координат
// антенны и миаяка в ГСК
//INPUT:
//tblPrfl - профиль звука
//arrAntY[3] - координаты антенны в ГСК
//arrTrueBeaconPos[3]- координаты маяка в ГСК
// OUTPUT:
// val_t - время
// возвращает true, если задача решается
bool calc_t_InpGSK(TTable_1D &tblPrfl, const double *arrAntY,const double *arrTrueBeaconPos
               ,  double &val_t)
{
    //1. формирование вектора Z ф-ла (3.6)
    double arrZ[3] = {0.};
    arrZ[0] = sqrt((arrTrueBeaconPos[0] - arrAntY[0]) * (arrTrueBeaconPos[0] - arrAntY[0])
            + (arrTrueBeaconPos[1] - arrAntY[1]) * (arrTrueBeaconPos[1] - arrAntY[1]) );

    arrZ[1]= -arrTrueBeaconPos[2]; // глуб. маяка

    arrZ[2] = - arrAntY[2];  // глуб. антенны
    // 3!

    // 4. решение уравнения (3.8) (или (2.3))
    double valTetta = 0.;
    if(!calcTetta_InpGSK(tblPrfl, arrAntY,arrTrueBeaconPos
                         ,  valTetta))
    {
        return false;
    }
    //5.вычисление времени подхода сигнала
    val_t =  calc_t( arrZ[2],arrZ[1],valTetta, tblPrfl);

    return true;

}
//------------------------------------
// Вычисление времени прохождения сигнала для варианта задания координат
// антенны и миаяка в ГСК
//INPUT:
//tblPrfl - профиль звука
//arrAntY[3] - координаты антенны в ГСК
//arrTrueBeaconPos[3]- координаты маяка в ГСК
// OUTPUT:
// val_t - время
// возвращает true, если задача решается
bool calc_t_and_dt_po_dYAnt_InpGSK(TTable_1D &tblPrfl, const double *arrAntY
                ,const double *arrTrueBeaconPos,  double &val_t, double* arrdt_po_dYAnt)
{
    //1. формирование вектора Z ф-ла (3.6)
    double arrZ[3] = {0.};
    arrZ[0] = sqrt((arrTrueBeaconPos[0] - arrAntY[0]) * (arrTrueBeaconPos[0] - arrAntY[0])
            + (arrTrueBeaconPos[1] - arrAntY[1]) * (arrTrueBeaconPos[1] - arrAntY[1]) );

    arrZ[1]= -arrTrueBeaconPos[2]; // глуб. маяка

    arrZ[2] = - arrAntY[2];  // глуб. антенны
    // 1!

    // 2. вычисление tetta и t
    double valTetta = 0., arr_dTetta_po_dZ[3] = {0.}, arr_dt_po_dZ[3] = {0.};
    if(!calc_t_and_dt_po_dZ( arrZ[2],arrZ[1],arrZ[0]
                             , tblPrfl, valTetta, arr_dTetta_po_dZ
                             , val_t, arr_dt_po_dZ))
    {
        return false;
    }
    //2!



    // 5. вчисление dZ_po_dYAnt ф-ла (3.11)
    double arr_dZ_po_dYAnt[9] = {0.};
    calc_dZ_po_dYAnt(arrAntY,arrTrueBeaconPos, arr_dZ_po_dYAnt);
    // 5!

    // вычисление dt_po_dYAnt
    MtrxMultMatrx(arr_dt_po_dZ,1, 3,  arr_dZ_po_dYAnt,3, arrdt_po_dYAnt);
    return true;

}
//------------------------------------
// Вычисление времени прохождения сигнала для варианта задания координат
// антенны и миаяка в ГСК
// маяк двигается равномерно и прямолинейно
//INPUT:
//tblPrfl - профиль звука
//arrSAnt[3] - координаты антенны в ГСК
//arrBeaconPos[3]- координаты маяка в ГСК на момент излучения
//arrBeaconVelo[3]- скорость маяка в ГСК
//если arrBeaconVelo==NULL то решается задача для неподвижного маяка !!
// OUTPUT:
// val_t - время
// возвращает true, если задача решается
bool calc_tZapr(TTable_1D &tblPrfl, const double *arrSAnt
                ,const double *arrBeaconPos,const double *arrBeaconVelo,  double &val_t)
{
    double arrdtZapr_po_dSBeacon[3] = {0.};
    // 0.
    if (arrBeaconVelo == NULL)
    {
        return calc_t_and_dt_po_dSBeacon_InpGSK(tblPrfl, arrSAnt
                        ,arrBeaconPos,  val_t, arrdtZapr_po_dSBeacon);
    }
    //1.вычисление начального приближения

    if(!calc_t_and_dt_po_dSBeacon_InpGSK(tblPrfl, arrSAnt
                    ,arrBeaconPos,  val_t,  arrdtZapr_po_dSBeacon))
    {
        return false;
    }
    // !1

    // 2. Решение нелинейного уравнения
    double arrBeaconPosCur[3] = {0.};

    double val_tZapr = 0.;
    bool breturn = false;
    for (int i = 0; i < 50; ++i)
    {
        arrBeaconPosCur[0] = arrBeaconPos[0] + val_t * arrBeaconVelo[0];
        arrBeaconPosCur[1] = arrBeaconPos[1] + val_t * arrBeaconVelo[1];
        arrBeaconPosCur[2] = arrBeaconPos[2] + val_t * arrBeaconVelo[2];
        calc_t_and_dt_po_dSBeacon_InpGSK(tblPrfl, arrSAnt
                        ,arrBeaconPosCur,  val_tZapr, arrdtZapr_po_dSBeacon);
        double fi = val_tZapr -val_t;
        double dfi_po_dt = ScalProduct(arrdtZapr_po_dSBeacon , arrBeaconVelo, 3) -1. ;
        double del_t = fi/ dfi_po_dt;
        val_t += -0.7 * del_t;

        if (fabs(del_t) <= 1.E-7)
        {
            breturn =true;
            break;
        }
    }
    return breturn;
}

//--------------------------------------
// ф-ла (3.11)
void calc_dZ_po_dYAnt(const double *arrYAnt,const double *arrTrueBeaconPos, double *arr_dZ_po_dYAnt)
{
    memset(arr_dZ_po_dYAnt, 0, 9 * sizeof(double));
    double r = sqrt((arrTrueBeaconPos[0] - arrYAnt[0]) * (arrTrueBeaconPos[0] - arrYAnt[0])
            +(arrTrueBeaconPos[1] - arrYAnt[1]) * (arrTrueBeaconPos[1] - arrYAnt[1]));
    arr_dZ_po_dYAnt [8] = -1.;
    arr_dZ_po_dYAnt[0] = -(arrTrueBeaconPos[0] - arrYAnt[0])/r;
    arr_dZ_po_dYAnt[1] =  -(arrTrueBeaconPos[1] - arrYAnt[1])/r;

}

//---------------------------------------------------------
// arr_dZ_po_dY[7*9] - матрица частных производных (Якоби)
// вектор-функции
// Z[0] = sqrt((Y3-Y0) * (Y3-Y0) +(Y4-Y1) * (Y4-Y1))
// Z[1] =  Y6  угол
// Z[2] =  Y7  угол
// Z[3] =  Y8  угол
// Z[4] =  Y3-Y0
// Z[5] =  Y4-Y1
// Z[6] =  Tetta(-Y[2],-Y[3], sqrt((Y3-Y0) * (Y3-Y0) +(Y4-Y1) * (Y4-Y1)))
void calc_dZ_po_dY(const double *arrY,double * arr_dTetta_po_dZ,double *arr_dZ_po_dY)
{
    memset(arr_dZ_po_dY, 0, 8 * 9 * sizeof(double));
    double r = sqrt((arrY[3] - arrY[0]) * (arrY[3] - arrY[0])
            +(arrY[4] - arrY[1]) * (arrY[4] - arrY[1]));

    arr_dZ_po_dY[0] = -(arrY[3] - arrY[0])/r;
    arr_dZ_po_dY[3] = -arr_dZ_po_dY[0];
    arr_dZ_po_dY[1] =  -(arrY[4] - arrY[1])/r;
    arr_dZ_po_dY[4] = -arr_dZ_po_dY[1];

    arr_dZ_po_dY [15] =arr_dZ_po_dY [25]=arr_dZ_po_dY [35] = 1.;

    arr_dZ_po_dY [4 * 9 ] = -1.;
    arr_dZ_po_dY [4 * 9 +3] = 1.;

    arr_dZ_po_dY [5 * 9 +1] = -1.;
    arr_dZ_po_dY [5 * 9 +4] = 1.;


    MtrxMultMatrx(arr_dTetta_po_dZ,1, 3, arr_dZ_po_dY,9, &arr_dZ_po_dY[6 *9]) ;
}
//---------------------------------------------------------
// arr_dz_po_dY[3*9] - матрица частных производных (Якоби)
// вектор-функции
// Z[0] = sqrt((Y3-Y0) * (Y3-Y0) +(Y4-Y1) * (Y4-Y1))
// Z[1] =  -Y5
// Z[2] =  -Y2
void calc_dz_po_dY(const double *arrY,double *arr_dz_po_dY)
{
    memset(arr_dz_po_dY, 0, 3 * 9 * sizeof(double));
    double r = sqrt((arrY[3] - arrY[0]) * (arrY[3] - arrY[0])
            +(arrY[4] - arrY[1]) * (arrY[4] - arrY[1]));

    arr_dz_po_dY[0] = -(arrY[3] - arrY[0])/r;
    arr_dz_po_dY[3] = -arr_dz_po_dY[0];
    arr_dz_po_dY[1] =  -(arrY[4] - arrY[1])/r;
    arr_dz_po_dY[4] = -arr_dz_po_dY[1];
    arr_dz_po_dY [14] =arr_dz_po_dY [20] = -1.;
}
//---------------------------------------------------------
// arr_dZ_po_dY[3*6] - матрица частных производных (Якоби)
// вектор-функции
// Z[0] = sqrt((Y3-Y0) * (Y3-Y0) +(Y4-Y1) * (Y4-Y1))
// Z[1] =  -Y5
// Z[2] =  -Y2
//
// arr_dZ_po_dY[3*6] - матрица частных производных (Якоби)
// вектор-функции
// Z[0] = sqrt((Y3-Y0) * (Y3-Y0) +(Y4-Y1) * (Y4-Y1))
// Z[1] = -Y5
// Z[2] = -Y2
//
void calc_dZ_3Vars_po_dY_6Vars(const double *arrYAnt,const double *arrTrueBeaconPos, double *arr_dZ_po_dY)
{
    memset(arr_dZ_po_dY, 0, 18 * sizeof(double));
    double r = sqrt((arrTrueBeaconPos[0] - arrYAnt[0]) * (arrTrueBeaconPos[0] - arrYAnt[0])
            +(arrTrueBeaconPos[1] - arrYAnt[1]) * (arrTrueBeaconPos[1] - arrYAnt[1]));
    arr_dZ_po_dY [11] =arr_dZ_po_dY [14] = -1.;
    arr_dZ_po_dY[0] = -(arrTrueBeaconPos[0] - arrYAnt[0])/r;
    arr_dZ_po_dY[3] = -arr_dZ_po_dY[0];
    arr_dZ_po_dY[1] =  -(arrTrueBeaconPos[1] - arrYAnt[1])/r;
    arr_dZ_po_dY[4] = -arr_dZ_po_dY[1];

}



//------------------------------------
// Вычисление времени прохождения сигнала для варианта задания координат
// антенны и миаяка в ГСК
//
//INPUT:
//tblPrfl - профиль звука
// arrEilers0[3] - палубные углы
//arrSAntPSK[3] - координаты антенны в ПСК (вектор параллакса)
// arrSVessGSK[3] - вектор положения центра корабля в ГСК
// arrSBeaconGSK[3] - вектор положения маяка в ГСК
// OUTPUT:
// val_t - время
// возвращает true, если задача решается
// arrdt_po_dX - градиент t по следующим 6- ти переменныч
// X[0], X[1], X[2] - вектор параллакса антенны в ПСК
// X[3], X[4], X[5] - вектор положения маяка в ГСК
bool calc_t_and_dt_po_dX(TTable_1D &tblPrfl,  double *arrEilers0,  double *arrSAntPSK
                , double *arrSVessGSK,const double *arrSBeaconGSK
                         ,double &val_t, double* arrdt_po_dX/*, QDataExchange &DataExchange*/)
{
    // персчет вектора положения антенны в ГСК
    double arrAntY[3] = {0.};
    // создание матрицы перехода из   КГСК в ПСК
     double arrEilers[3] = {0.},  arr_KGSK[3] = {0.};
    MatrxMultScalar(arrEilers0, 1, 3, -1.,arrEilers);
    double matrPereh_PSK_V_KGSK[9] = {0} ;
    QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilers, matrPereh_PSK_V_KGSK) ;
    // вектор положения в ПСК-центр тяжести
    // вычисление вектора положениея в ПСК

    MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrSAntPSK, 1, arr_KGSK) ;

    MtrxSumMatrx(arr_KGSK, arrSVessGSK,1, 3, arrAntY) ;
    // !


    // 2.
    //1. формирование вектора Z ф-ла (3.6)
    double arrZ[3] = {0.};
    arrZ[0] = sqrt((arrSBeaconGSK[0] - arrAntY[0]) * (arrSBeaconGSK[0] - arrAntY[0])
            + (arrSBeaconGSK[1] - arrAntY[1]) * (arrSBeaconGSK[1] - arrAntY[1]) );

    arrZ[1]= -arrSBeaconGSK[2]; // глуб. маяка

    arrZ[2] = - arrAntY[2];  // глуб. антенны
    // 1!

    // 2. вычисление tetta и t
    double valTetta = 0., arr_dTetta_po_dZ[3] = {0.}, arr_dt_po_dZ[3] = {0.};
    if(!calc_t_and_dt_po_dZ( arrZ[2],arrZ[1],arrZ[0]
                             , tblPrfl, valTetta, arr_dTetta_po_dZ
                             , val_t, arr_dt_po_dZ))
    {
        return false;
    }
    //2!



    double arr_dZ_po_dY[18] = {0.};
    calc_dZ_3Vars_po_dY_6Vars(arrAntY, arrSBeaconGSK, arr_dZ_po_dY);

    double arrdY_po_dX[36] = {0.};
    for (int i =0; i < 3; ++i)
        for (int j =0; j < 3; ++j)
        {
          arrdY_po_dX[i * 6 + j]  = matrPereh_PSK_V_KGSK[ i * 3 + j];
        }
    arrdY_po_dX[3 * 6 + 3] = arrdY_po_dX[4 * 6 + 4] = arrdY_po_dX[5 * 6 + 5] = 1.;

    double arrT0[6] = {0.};

    MtrxMultMatrx( arr_dt_po_dZ,1, 3, arr_dZ_po_dY,6, arrT0) ;
    MtrxMultMatrx( arrT0,1, 6, arrdY_po_dX,6, arrdt_po_dX) ;

   // DataExchange = QDataExchange(arrAntY,arrZ, arr_dTetta_po_dZ, valTetta);
    return true;

}

//--------------------------------------
// пересчет замера 3D в ГСК
// OUTPUT:
// arrSZv_GSK[3] - замер в ГСК
// valTZv - время излучения ответного сигнала
bool transfMeasure_to_GSK(TTable_1D &tblPrfl,QBigMeasure &BigMeasureInp
       ,double* arrAntPosParams,  double *arrSZv_GSK, double*valTZv)
{
    // 1. вычисление единичного вектора направления принятого
    // сигнала маяка в КГСК

    double arrV[3] = {0.}; // вектор измерений
    arrV[0] = BigMeasureInp.mqzv;
    arrV[1] =1.;
    arrV[2] = BigMeasureInp.mezv;

    double arr_e_kgsk[3] = {0.}, arrTemp0[6] = {0.};
    memcpy(&arrTemp0[3],&arrAntPosParams[3], 3 * sizeof(double) );

    QPeaceVess::recalcPositionFromGdgSphericalSK_to_KGSK(arrV,BigMeasureInp.marrMuZv
                                                            ,arrTemp0, arr_e_kgsk);
    // 1!

    // 2.Вычисление курсового угла направления на маяк относительно
    // судна в момент приема ответного сигнала и угла скольжения луча ответного сигнала

    double val_q_gsk = QPeaceVess::calcCourseAngle(arr_e_kgsk[0], arr_e_kgsk[1]); // курс угол
    double val_e_gsk = asin(arr_e_kgsk[2]); //угол места

    // 2!

    // 3.вычисление вектора положения антенны в момент приема в ГСК
    double arrV1[3] = {0.}, arrSAnt_KGSK[3] = {0.}, arrSAnt_GSK[3] = {0.};
    QPeaceVess::recalcPositionFromGdgSphericalSK_to_KGSK(arrV1,BigMeasureInp.marrMuZv
                                                            ,arrAntPosParams, arrSAnt_KGSK);
    MtrxSumMatrx(BigMeasureInp.marrSVessZv, arrSAnt_KGSK,1, 3, arrSAnt_GSK) ;
    // 3!

    // 4.вычисление вектора положения антенны в момент излучения в ГСК
    double  arrSAntWave_KGSK[3] = {0.}, arrSAntWave_GSK[3] = {0.};
    QPeaceVess::recalcPositionFromGdgSphericalSK_to_KGSK(arrV1,BigMeasureInp.marrMuWaveZv
                                                            ,arrAntPosParams, arrSAntWave_KGSK);
    MtrxSumMatrx(BigMeasureInp.marrSVessWaveZv, arrSAntWave_KGSK,1, 3, arrSAntWave_GSK) ;
    // 4!

    //5. начальное приближение времени ответного сигнала
    const double VAl_Delta_T = (BigMeasureInp.mTotvZv - BigMeasureInp.mTzaprZv - BigMeasureInp.mTobr);
    double val_tOtv = VAl_Delta_T/2.;
    // 5!

    //6. решение нелинейного уравнения относительно времени ответного сигнала
    double val_tZapr = 0., val_dtZapr_po_dtOtv = 0., arrGskBeacon[3] = {0.};
    bool breturn = false;
    for(int i =0; i < 100; ++i)
    {
       if(! calc_TZapr_and_dTZapr_po_dTOtv(tblPrfl,arrSAnt_GSK, arrSAntWave_GSK, val_q_gsk, -val_e_gsk
                    ,val_tOtv, arrGskBeacon,val_tZapr, val_dtZapr_po_dtOtv))
       {
           return false;
       }

       double fi = (val_tOtv + val_tZapr -VAl_Delta_T);
       double dfi = 1. + val_dtZapr_po_dtOtv;
       double dt = fi/dfi;
       val_tOtv -= 0.5 * dt;
       if (fabs(dt)<= 1.E-7)
       {
           breturn = true;
           break;
       }
    }
    if(!breturn)
    {
        return false;
    }
    memcpy(arrSZv_GSK, arrGskBeacon,3 * sizeof(double));

    *valTZv = BigMeasureInp.mTotvZv - val_tOtv;

 return true;
}

//--------------------------------------
// пересчет замера 3D в ГСК
// arrAntPosParams[6] - вектор позиционирования антенны в ПСК
bool transfMeasure_to_GSK(TTable_1D &tblPrfl,QBigMeasure &BigMeasureInp
       ,double* arrAntPosParams,const double *arrBeaconVelo
       ,  double *arrGskBeacon, double*valTZv)
{
    // 1. вычисление единичного вектора направления принятого
    // сигнала маяка в КГСК
    double arrV[3] = {0.}; // вектор измерений
    arrV[0] = BigMeasureInp.mqzv;
    arrV[1] =1.;
    arrV[2] = BigMeasureInp.mezv;

    double arr_e_kgsk[3] = {0.}, arrTemp0[6] = {0.};
    memcpy(&arrTemp0[3],&arrAntPosParams[3], 3 * sizeof(double) );

    QPeaceVess::recalcPositionFromGdgSphericalSK_to_KGSK(arrV,BigMeasureInp.marrMuZv
                                                            ,arrTemp0, arr_e_kgsk);
    // 1!

    // 2.Вычисление курсового угла направления на маяк относительно
    // судна в момент приема ответного сигнала и угла скольжения луча ответного сигнала

    double val_q_gsk = QPeaceVess::calcCourseAngle(arr_e_kgsk[0], arr_e_kgsk[1]); // курс угол
    double val_e_gsk = asin(arr_e_kgsk[2]); //угол места

    // 2!

    // 3.вычисление вектора положения антенны в момент приема в ГСК
    double arrV1[3] = {0.}, arrSAnt_KGSK[3] = {0.}, arrSAnt_GSK[3] = {0.};
    QPeaceVess::recalcPositionFromGdgSphericalSK_to_KGSK(arrV1,BigMeasureInp.marrMuZv
                                                            ,arrAntPosParams, arrSAnt_KGSK);
    MtrxSumMatrx(BigMeasureInp.marrSVessZv, arrSAnt_KGSK,1, 3, arrSAnt_GSK) ;
    // 3!

    // 4.вычисление вектора положения антенны в момент излучения в ГСК
    double  arrSAntWave_KGSK[3] = {0.}, arrSAntWave_GSK[3] = {0.};
    QPeaceVess::recalcPositionFromGdgSphericalSK_to_KGSK(arrV1,BigMeasureInp.marrMuWaveZv
                                                            ,arrAntPosParams, arrSAntWave_KGSK);
    MtrxSumMatrx(BigMeasureInp.marrSVessWaveZv, arrSAntWave_KGSK,1, 3, arrSAntWave_GSK) ;
    // 4!

    //5. начальное приближение времени ответного сигнала
    const double VAl_Delta_T = (BigMeasureInp.mTotvZv - BigMeasureInp.mTzaprZv - BigMeasureInp.mTobr);
    double val_tOtv = VAl_Delta_T/2.;
    // 5!

    //6. решение нелинейного уравнения относительно времени ответного сигнала
    double val_tZapr = 0., val_dtZapr_po_dtOtv = 0.;
    bool breturn = false;
    for(int i =0; i < 100; ++i)
    {

       if(! calc_TZapr_and_dTZapr_po_dTOtv_(tblPrfl,arrSAnt_GSK, arrSAntWave_GSK, arrBeaconVelo
                    ,BigMeasureInp.mTobr,val_q_gsk, -val_e_gsk
                    ,val_tOtv, arrGskBeacon,val_tZapr, val_dtZapr_po_dtOtv))
       {
           return false;
       }

       double fi = (val_tOtv + val_tZapr -VAl_Delta_T);
       double dfi = 1. + val_dtZapr_po_dtOtv;
       double dt = fi/dfi;
       val_tOtv -= 0.5 * dt;
       if (fabs(dt)<= 1.E-6)
       {
           breturn = true;
           break;
       }
    }
    if(!breturn)
    {
        return false;
    }   

    *valTZv = BigMeasureInp.mTotvZv - val_tOtv;

 return true;
}
//---------------------------------------------------
// вычисление момента излучения запросного сигнала и его производной по мементу ответного сигнала
// заданы 2 точки положения антенны в ГСК
// arrSAnt_GSK[3] и arrSAntWave_GSK[3]
// val_tOtv - время прохождения сигнала от маяка до точки arrSAnt_GSK
// val_q_gsk и val_tetta - углоые координаты луча в антенной сист. кординат
// требуется найти время прохождения луча от точки arrSAntWave_GSK до маяка
// и производную этого времени по val_tOtv
// INPUT:
// tblPrfl -профиль  скорости звука
// arrSAnt_GSK[3] - положение антенны в момент приема
// arrSAntWave_GSK[3] - положение антенны в момент излучения
// val_q_gsk, val_e_gsk - угловые координаты луча
// val_tOtv - время прохождения ответного сигнала
//OUTPUT:
//arrGskBeacon[3] -  положения маяка
// val_tZapr - время прохождения запросного сигнала
// val_dtZapr_po_dtOtv - производная
bool calc_TZapr_and_dTZapr_po_dTOtv(TTable_1D &tblPrfl, double* arrSAnt_GSK,double* arrSAntWave_GSK
               ,const double val_q_gsk, const double val_tetta
          ,const double val_tOtv, double *arrGskBeacon,double&val_tZapr,double& val_dtZapr_po_dtOtv)
{
    // 1. вычисление вектора положения луча в ГСК и его матрицы Якоби по времени
    //задача решается в плоскости распространения ответного луча (+)
   double arrt[3] = {0.},arr_dGskBeacon_po_dt[3] = {0.};
   calcBeamPosGSK_from_t_and_dBeamPosGSK_po_dt(tblPrfl,  arrSAnt_GSK
                         , val_q_gsk,  val_tetta
                              ,val_tOtv,arrGskBeacon ,arr_dGskBeacon_po_dt);


   // !1

    // 2. вычисление времени запросного сигнала и его матрицы Якоби по вектору координат маяка
    double arr_dt_po_darrSBeacon_GSK[3] ={0.};

    calc_t_and_dt_po_dSBeacon_InpGSK(tblPrfl, arrSAntWave_GSK
                    ,arrGskBeacon,  val_tZapr, arr_dt_po_darrSBeacon_GSK);

    // !2

    //3. вычисление производной времени запросного сигнала по времени ответного
    val_dtZapr_po_dtOtv = ScalProduct(arr_dt_po_darrSBeacon_GSK , arr_dGskBeacon_po_dt, 3) ;
    // !3

    return true;
}

//---------------------------------------------------
// вычисление момента излучения запросного сигнала и его производной по мементу ответного сигнала
// заданы 2 точки положения антенны в ГСК
// arrSAnt_GSK[3] и arrSAntWave_GSK[3]
// val_tOtv - время прохождения сигнала от маяка до точки arrSAnt_GSK
// val_q_gsk и val_tetta - углоые координаты луча в антенной сист. кординат
// требуется найти время прохождения луча от точки arrSAntWave_GSK до маяка
// и производную этого времени по val_tOtv
// INPUT:
// tblPrfl -профиль  скорости звука
// arrSAnt_GSK[3] - положение антенны в момент приема
// arrSAntWave_GSK[3] - положение антенны в момент излучения
// val_q_gsk, val_e_gsk - угловые координаты луча
// val_tOtv - время прохождения ответного сигнала
//OUTPUT:
//arrGskBeacon[3] -  положения маяка на момент ответа
// val_tZapr - время прохождения запросного сигнала
// val_dtZapr_po_dtOtv - производная
bool calc_TZapr_and_dTZapr_po_dTOtv_(TTable_1D &tblPrfl, double* arrSAnt_GSK
    ,double* arrSAntWave_GSK,const double *arrBeaconVelo,const double VAlTobr
    , const double val_q_gsk, const double val_tetta
    ,const double val_tOtv, double *arrGskBeacon,double&val_tZapr,double& val_dtZapr_po_dtOtv)
{
    // 1. вычисление вектора положения луча в ГСК и его матрицы Якоби по времени
    // задача решается в плоскости распространения ответного луча (+)
   double arrt[3] = {0.},arr_dGskBeacon_po_dt[3] = {0.};
   //  ошибка здесь!
   if(!calcBeamPosGSK_from_t_and_dBeamPosGSK_po_dt(tblPrfl,  arrSAnt_GSK
              ,val_q_gsk, val_tetta,val_tOtv,arrGskBeacon ,arr_dGskBeacon_po_dt))
   {
       return false;
   }

   // !1

    // 2. вычисление времени запросного сигнала и его матрицы Якоби по вектору координат маяка
    double arr_dt_po_darrSBeaconPriem_GSK[3] ={0.};
    double arrGskBeaconPriem[3] = {0.}; // положение маяка на момент приема
    memcpy(arrGskBeaconPriem, arrGskBeacon, 3 * sizeof(double));
    if (arrBeaconVelo != NULL)
    {
        arrGskBeaconPriem[0] -=  VAlTobr * arrBeaconVelo[0];
        arrGskBeaconPriem[1] -=  VAlTobr * arrBeaconVelo[1];
        arrGskBeaconPriem[2] -=  VAlTobr * arrBeaconVelo[2];
    }


    if(!calc_t_and_dt_po_dSBeacon_InpGSK(tblPrfl, arrSAntWave_GSK
                    ,arrGskBeaconPriem,  val_tZapr, arr_dt_po_darrSBeaconPriem_GSK))
    {
        return false;
    }

    // !2

    //3. вычисление производной времени запросного сигнала по времени ответного
    val_dtZapr_po_dtOtv = ScalProduct(arr_dt_po_darrSBeaconPriem_GSK , arr_dGskBeacon_po_dt, 3) ;
    // !3

    return true;
}


//-------------------------------------------------------------
bool calcDeepth_and_dDeepth_po_dt(TTable_1D &tblPrfl, const double VAlza, const double val_tetta
                                 ,const double val_t, double &val_zm, double &val_dz_po_dt)
{
    val_zm = tblPrfl.calcValue(10.)* sin(val_tetta) * val_t;
    double val_costetta =cos(val_tetta);
    double cza = tblPrfl.calcValue(VAlza);
    // вычисление начального приближения
    double zcur = VAlza;
    double valIntegral =0.;
    double step = 0.01;
    int Nc= tblPrfl.mparrArg[tblPrfl.mNumCols - 1]/step;
    bool breturn1 =false;
    for(int i=0; i < Nc; ++i)
    {
      val_zm =  VAlza + ((double)i)*step;
      double val_nz =cza/tblPrfl.calcValue(val_zm);
      double f = val_nz * val_nz /(cza * sqrt(val_nz * val_nz -val_costetta *val_costetta));
      valIntegral +=  step * f;
      if (valIntegral >= val_t)
      {
         breturn1 = true;
         break;
      }
    }
    if(!breturn1)
    {
        return false;
    }
 // Итерационный процесс
    bool breturn = false;
    double psi = 0., dpsi = 0.;

// !!!!
    for(int i = 0; i < 50; ++i)
    {
      psi =  calcSumIntegral_I( VAlza,val_zm, tblPrfl,val_costetta,calcIntegral_I2)/cza
              - val_t;
      double val_nz =cza/tblPrfl.calcValue(val_zm);
      dpsi = val_nz * val_nz /(cza * sqrt(val_nz * val_nz -val_costetta *val_costetta));
      double dz = psi/dpsi;
      val_zm -= 0.5 * dz;
      if(fabs(dz) < 0.001)
      {
          breturn =true;
          break;
      }

    }
    if(!breturn)
    {
        return false;
    }

    double val_nz =cza/tblPrfl.calcValue(val_zm);
    val_dz_po_dt =cza * sqrt(val_nz * val_nz -val_costetta *val_costetta)/( val_nz * val_nz);


    return true;
}

//-------------------------------------------
bool calcBeamPosGSK_from_t_and_dBeamPosGSK_po_dt(TTable_1D &tblPrfl, double* arrSAnt_GSK
                      ,const double val_q_gsk, const double val_tetta
                           ,const double val_t,double *arrBeamGskPos ,double *arr_dBeamGskPos_po_dt)
{
    // 1.вычисление глубины маяка и ее прозводной по времени
    double valZm = 0., val_dZm_po_dt = 0.;
    if(!calcDeepth_and_dDeepth_po_dt(tblPrfl, -arrSAnt_GSK[2],val_tetta
                                 ,val_t, valZm, val_dZm_po_dt))
    {
        return false;
    }
    // 1!

    //2.вычисление  YГСК_с_волной и ее производной по времени
    double val_yGskWave = cos(val_tetta)
            *calcSumIntegral_I( -arrSAnt_GSK[2],valZm
                , tblPrfl, cos(val_tetta),calcIntegral_I);
    double cz = tblPrfl.calcLinearValueApprox(valZm);
    double cza = tblPrfl.calcLinearValueApprox(-arrSAnt_GSK[2]);
    double val_dyGskWave_po_dt = cos(val_tetta) * cz * cz/cza;
    // 2!

    // 3. формирование SGskWave и Якобиана dSGskWave_po_dt
    double arrGskBeamWave[3] = {0.}, arr_dGskBeamWave_po_dt[3] = {0.};
    arrGskBeamWave[1] = val_yGskWave;
    arrGskBeamWave[2] = valZm;
    arr_dGskBeamWave_po_dt[1] = val_dyGskWave_po_dt;
    arr_dGskBeamWave_po_dt[2] = val_dZm_po_dt;
    // 3!

    // 4.матрица поворота
    double arrL[9] = {0.};
    arrL[0] = arrL[4] = cos(val_q_gsk);
    arrL[1] = -sin(val_q_gsk);
    arrL[3] = -arrL[1];
    arrL[8] = -1.;
    // 4!

     // 5. вычисление вектора положения луча в ГСК и его матрицы Якоби по времени
    double arrt[3] = {0.};
    MtrxTranspMultMatrx(arrL,3, 3, arrGskBeamWave,1, arrt) ;
    arrBeamGskPos[0] = arrt[0] + arrSAnt_GSK[0];
    arrBeamGskPos[1] = arrt[1] + arrSAnt_GSK[1];
    arrBeamGskPos[2] = arrt[2] ;

    MtrxTranspMultMatrx(arrL,3, 3, arr_dGskBeamWave_po_dt,1, arr_dBeamGskPos_po_dt) ;
    // 5!

   return true;
}

//---------------------------------
// пересчет замера 3D в АСПК
bool transfMeasure_to_ASPK(TTable_1D &tblPrfl,QBigMeasure &BigMeasureInp, double* arrAntPosParams,  double *arrSZv_ASPK, double*valTZv)
{
    // 1. пересчет замера в ГСК
    double arrSZv_GSK[3] = {0.}, arrKGSK[3] = {0.};
    if(!transfMeasure_to_GSK(tblPrfl,BigMeasureInp, arrAntPosParams,  arrSZv_GSK, valTZv))
    {
        return false;
    }
    // !1

    // 2. пересчет замера в КГСК
    MtrxMinusMatrx(arrSZv_GSK, BigMeasureInp.marrSVessZv,3, 1, arrKGSK);
    // !2

    // 3. пересчет замера в АСПК

    QPeaceVess::recalcVect_KGSK_INTO_GdgSobSK( arrKGSK
               , BigMeasureInp.marrMuZv, NULL
               ,  arrAntPosParams,arrSZv_ASPK,3 );
    // !3
return true;
}
//------------------------------------
// Вычисление времени прохождения сигнала для варианта задания координат
// антенны и миаяка в ГСК
//INPUT:
//tblPrfl - профиль звука
//arrAntY[3] - координаты антенны в ГСК
//arrTrueBeaconPos[3]- координаты маяка в ГСК
// OUTPUT:
// val_t - время
// возвращает true, если задача решается
bool calc_t_and_dt_po_dSBeacon_InpGSK(TTable_1D &tblPrfl, const double *arrSAnt
                ,const double *arrSBeacon,  double &val_t, double* arrdt_po_dSBeacon)
{
    //1. формирование вектора Z ф-ла (3.6)
    double arrZ[3] = {0.};
    arrZ[0] = sqrt((arrSBeacon[0] - arrSAnt[0]) * (arrSBeacon[0] - arrSAnt[0])
            + (arrSBeacon[1] - arrSAnt[1]) * (arrSBeacon[1] - arrSAnt[1]) );

    arrZ[1]= -arrSBeacon[2]; // глуб. маяка

    arrZ[2] = - arrSAnt[2];  // глуб. антенны
    // 1!

    // 2. вычисление tetta и t
    double valTetta = 0., arr_dTetta_po_dZ[3] = {0.}, arr_dt_po_dZ[3] = {0.};
    if(!calc_t_and_dt_po_dZ( arrZ[2],arrZ[1],arrZ[0]
                             , tblPrfl, valTetta, arr_dTetta_po_dZ
                             , val_t, arr_dt_po_dZ))
    {
        return false;
    }
    //2!


    // 5. вчисление dZ_po_dYAnt ф-ла (3.11)
    double arr_dZ_po_dSBeacon[9] = {0.};
    calc_dZ_po_dSBeacon(arrSAnt,arrSBeacon, arr_dZ_po_dSBeacon);
    // 5!

    // вычисление dt_po_dSBeacon
    MtrxMultMatrx(arr_dt_po_dZ,1, 3,  arr_dZ_po_dSBeacon,3, arrdt_po_dSBeacon);
    return true;

}

//--------------------------------------
//
void calc_dZ_po_dSBeacon(const double *arrSAnt
                         ,const double *arrSBeacon, double *arr_dZ_po_dSBeacon)
{
    memset(arr_dZ_po_dSBeacon, 0, 9 * sizeof(double));
    double r = sqrt((arrSBeacon[0] - arrSAnt[0]) * (arrSBeacon[0] - arrSAnt[0])
            +(arrSBeacon[1] - arrSAnt[1]) * (arrSBeacon[1] - arrSAnt[1]));
    arr_dZ_po_dSBeacon [5] = -1.;
    arr_dZ_po_dSBeacon[0] = (arrSBeacon[0] - arrSAnt[0])/r;
    arr_dZ_po_dSBeacon[1] =  (arrSBeacon[1] - arrSAnt[1])/r;

}
//------------------------------------------------
// Преобразование координат точки из ГСК в Антенную сферическую систему координат
//INPUT:
//tblPrfl - профиль скорости звука
//arrGSK_XYZ[3] - положение точки в ГСК
// arrSVessGSK[3] - положение корабля в ГСК
//arrMu[3] - вектор палубных углов
//arrAntPosParams[6] - вектор параметров позиционирования антенны
//OUTPUT:
//pval_q - курсовой угол сигнала в АСфСК
//pval_e - угол места сигнала в АСфСК
//pval_t - время прохождения сигнала
bool transf_GSK_XYZ_to_USBL3D(TTable_1D &tblPrfl, double *arrGSK_XYZ, double *arrSVessGSK, double *arrMu
  ,double* arrAntPosParams,  double *pval_q, double *pval_e,  double *pval_t)
{
    // 1. вычисление положения антенны в КГСК
    double arr_Ant_KGSK[3] = {0.};
    QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( arrAntPosParams, arrMu
                 , NULL,arr_Ant_KGSK,3);
    // !1

    // 2. вычисление положения антенны в ГСК
    double arr_Ant_GSK[3] = {0.};
    MtrxSumMatrx(arr_Ant_KGSK, arrSVessGSK,1, 3, arr_Ant_GSK) ;
    // !2

    // 3.Вычисление угла скольжения луча и угла места в ГСК
    double valTetta = 0., val_e_gsk = 0.;;
    if(!calcTetta_InpGSK(tblPrfl, arr_Ant_GSK,arrGSK_XYZ, valTetta))
    {
        return false;
    }
    val_e_gsk = -valTetta;
    // !3

    // 4. Вычисление курсового угла принятого луча в ГСК
    double val_q_gsk = QPeaceVess::calcCourseAngle(arrGSK_XYZ[0] - arr_Ant_GSK[0], arrGSK_XYZ[1] - arr_Ant_GSK[1]);
    // !4

    // 5.Формирование вектора единичного направления принятого сигнала в ГСК
    double arr_e_gsk[3] = {0.};
    arr_e_gsk[0] = cos(val_e_gsk) * sin(val_q_gsk);
    arr_e_gsk[1] = cos(val_e_gsk) * cos(val_q_gsk);
    arr_e_gsk[2] = sin(val_e_gsk) ;
    // !5

    // 6. Пересчет вектора единичного направления сигнала в АСПК
    double arr_e_aspk[3] = {0.};
    // пересчет вектора соcтояния из KGSK в собственную систему
    // если  lenarrKGSK == 6 , то пересчитывается положение и скорость
    // если   lenarrKGSK == 3   , то пересчитывается только положение
    // на вход подается вектор состоящий из 3 или 6  координат.
    // первые 3 координаты представляют из себя  положение точки в КГС
    // последние 3 координаты представляют из себя скорость точки в КГСК
    // INPUT:
    //arrKGSK[lenarrKGSK] - вектор положния (положения и скорости) цели
    //arrGdgPosParams[6]  - вектор положение гаджета в ПСК, вектор углов ориентации гаджета в ПСК
    // arrEilerCntrKP[3] -вектора углов ориентации корабля (Q, Tetta, Psi)
    // arrOmegaPSK0[3] - вектор угловых скоростей (Psi, Tet, Q)
    // OUTPUT:
    // arrSobSK- вектор положения (или положениея и скорости) цели в ПСК -гаджет
    double arrPosCur[6] = {0.};
    memcpy(&arrPosCur[3], &arrAntPosParams[3], 3 * sizeof(double) );
    QPeaceVess::recalcVect_KGSK_INTO_GdgSobSK(arr_e_gsk
               , arrMu, NULL,  arrPosCur,arr_e_aspk,3 );
    // !6

    // 7. Вычисление угловых координат единичного вектора направления принятого сигнала в АСПК
    *pval_q = QPeaceVess::calcCourseAngle(arr_e_aspk[0], arr_e_aspk[1]);
    *pval_e = asin(arr_e_aspk[2]);
    // !7

    // 8. Вычисление времени прохождения сигнала
    *pval_t = calc_t( -arr_Ant_GSK[2],-arrGSK_XYZ[2]
                      , valTetta, tblPrfl);
    // !8
    return true;

}
//------------------------------------
// Вычисление времени прохождения сигнала для варианта задания координат
// корабля и маяка в ГСК
//
//INPUT:
//tblPrfl - профиль звука
// arrEilers0[3] - палубные углы
//arrSAntPSK[3] - координаты антенны в ПСК (вектор параллакса)
// arrSVessGSK[3] - вектор положения центра корабля в ГСК
// arrSBeaconGSK[3] - вектор положения маяка в ГСК
// OUTPUT:
// val_t - время
// возвращает true, если задача решается
// arrdt_po_dX - градиент t по следующим 6- ти переменныч
// X[0], X[1], X[2] - вектор параллакса антенны в ПСК
// X[3], X[4], X[5] - вектор положения маяка в ГСК
bool calc_t_for_Vess_gsk(TTable_1D &tblPrfl,  double *arrEilers0,  double *arrSAntPSK
                , double *arrSVessGSK,const double *arrSBeaconGSK
                         ,double &val_t)
{
    // персчет вектора положения антенны в ГСК
    double arrAntY[3] = {0.};
    // создание матрицы перехода из   КГСК в ПСК
     double arrEilers[3] = {0.},  arr_KGSK[3] = {0.};
    MatrxMultScalar(arrEilers0, 1, 3, -1.,arrEilers);
    double matrPereh_PSK_V_KGSK[9] = {0} ;
    QCoordSystTrsf::calcMatr_PSK_v_KGSK_LeftRot( arrEilers, matrPereh_PSK_V_KGSK) ;
    // вектор положения в ПСК-центр тяжести
    // вычисление вектора положениея в ПСК

    MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrSAntPSK, 1, arr_KGSK) ;

    MtrxSumMatrx(arr_KGSK, arrSVessGSK,1, 3, arrAntY) ;
    // !

    if(!calc_t_InpGSK(tblPrfl, arrAntY,arrSBeaconGSK,  val_t))
    {
        return false;
    }

    return true;

}
//--------------------------------------
// пересчет угловых координат цели из антенной сферической системы координат
// в корабельную сферическую систему координат
// INPUT:
//tblPrfl - профиль скорости звука
//VAl_ASSK_q - КУ АССК
//VAl_ASSK_e - УМ в АССК
//arrEilersMu[3] - углы Эйлера ПСК
// arrAntPosParams[6] - вектор позиционирования антенны в ПСК
// OUTPUT:
//*pval_KGSK_q - КУ
//*pval_KGSK_e УМ
void calcSphericalAnglesKGSK(TTable_1D &tblPrfl, const double VAl_ASSK_q, const double VAl_ASSK_e, double *arrEilersMu
       ,double* arrAntPosParams, double *pval_KGSK_q, double *pval_KGSK_e)
{
    // 1. вычисление единичного вектора направления принятого
    // сигнала маяка в КГСК
    double arrV[3] = {0.}; // вектор измерений
    arrV[0] = VAl_ASSK_q;
    arrV[1] =1.;
    arrV[2] = VAl_ASSK_e;

    double arr_e_kgsk[3] = {0.}, arrTemp0[6] = {0.};
    memcpy(&arrTemp0[3],&arrAntPosParams[3], 3 * sizeof(double) );

    QPeaceVess::recalcPositionFromGdgSphericalSK_to_KGSK(arrV,arrEilersMu
                                                            ,arrTemp0, arr_e_kgsk);
    // 1!

    // 2.Вычисление курсового угла направления на маяк относительно
    // судна в момент приема ответного сигнала и угла скольжения луча ответного сигнала

    *pval_KGSK_q = QPeaceVess::calcCourseAngle(arr_e_kgsk[0], arr_e_kgsk[1]); // курс угол
    *pval_KGSK_e = asin(arr_e_kgsk[2]); //угол места

    // 2!


 return;
}

