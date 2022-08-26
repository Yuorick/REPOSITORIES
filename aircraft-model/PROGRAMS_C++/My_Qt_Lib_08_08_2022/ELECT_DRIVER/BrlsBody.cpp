#include "BrlsBody.h"
#include <string.h>
#include <dir.h>
# include <QString>


extern const double VAL_BRLS_Cx;
extern const double ARrBeginPlot_Furk2540[] =
{0., 0.37/ 2.,0.6/2., 1.1/2.
};

extern const int QUantColsReport_Wind = 8;

QBrlsBody::QBrlsBody()
{

    mQuantPlots = QUANT_PLOTS;
    mcs = 0.;

}
//---------------------------------------------------------------------------
/*
QBrlsBody::~QBrlsBody()
{
    if(marrPlots) delete marrPlots ;//delete []marrPlots ;
    marrPlots = NULL ;
    mQuantPlots = 0;
    if(marrBeginPlot) delete []marrBeginPlot;
    marrBeginPlot = NULL;
}
*/
//---------------------------------------------------------------------------


  // конструктор копирования
  QBrlsBody ::QBrlsBody (const QBrlsBody &R):QBrls(R)
  {
      mQuantPlots  = R.mQuantPlots ;
      memcpy( marrPlots,R.marrPlots, R.mQuantPlots * sizeof(QBeamPlot));
      memcpy( marrBeginPlot,R.marrBeginPlot, R.mQuantPlots * sizeof(double));
      mcs = R.mcs;

  }
//---------------------------------------------------------
  // оператор присваивания
QBrlsBody &QBrlsBody::operator=(const QBrlsBody  &R)
  {
      if(this == &R)
      {
          return *this;
      }
     QBrls:: operator= (R);
     mQuantPlots  = R.mQuantPlots ;
     memcpy( marrPlots,R.marrPlots, R.mQuantPlots * sizeof(QBeamPlot));
     memcpy( marrBeginPlot,R.marrBeginPlot, R.mQuantPlots * sizeof(double));
     mcs = R.mcs;

     return *this ;
  }
/*
//---------------------------------------------------------
QBrlsBody ::QBrlsBody( const double  Cv, const double   Cx, const double  VAlM
                       , const double  VAlLengthTotal,const double  VAlHeight
                       ,const double  Max_dOm_po_dt, const QSegment SgmBarrier
                       , const TPlane Plane, const QBeamPlot *arrPlots,const double* arrBeginPlot
                       ,const double cs)
    :QBrls( Cv, Cx, VAlM
                , VAlLengthTotal, VAlHeight
                , Max_dOm_po_dt, SgmBarrier
                , Plane)
{
   mQuantPlots = QUANT_PLOTS ;
   memcpy( marrBeginPlot, arrBeginPlot, QUANT_PLOTS * sizeof(double));
   memcpy( marrPlots,arrPlots, QUANT_PLOTS * sizeof(QBeamPlot));
   mcs = cs;

}
*/
//---------------------------------------------------------------------------

QBrlsBody ::QBrlsBody( const double  Cv, const double   Cx, const double  VAlM
                       , const double  VAlLengthTotal,const double  VAlHeight
                       ,const double  Max_dOm_po_dt, const QSegment SgmBarrier
                       , const TPlane Plane, const double  cs,const int ITypeOfBRLS)
    :QBrls( Cv, Cx, VAlM
                , VAlLengthTotal, VAlHeight
                , Max_dOm_po_dt, SgmBarrier
                , Plane)
{
    mcs = cs;
   switch (ITypeOfBRLS)
   {
   case 0: // furke 2540
       mQuantPlots = QUANT_PLOTS ;



       //marrPlots[0] = QBeamPlot(0);
       //marrPlots[1] = QBeamPlot(1);
       //marrPlots[2] = QBeamPlot(2);
      // marrPlots[3] = QBeamPlot(3);
       QBeamPlot::createBeamPlotFurke2540(mcs,0, (marrPlots[0]));
       QBeamPlot::createBeamPlotFurke2540(mcs,1, (marrPlots[1]));
       QBeamPlot::createBeamPlotFurke2540(mcs,2, (marrPlots[2]));
       QBeamPlot::createBeamPlotFurke2540(mcs,3, (marrPlots[3]));

       memcpy(marrBeginPlot, ARrBeginPlot_Furk2540, mQuantPlots * sizeof(double));

       marrPlots[0].marrPlgProfile[0].calcBoundBox();
       mHeight = marrPlots[0].marrPlgProfile[0].Box[3] - marrPlots[0].marrPlgProfile[0].Box[1];
       mJPayLoad = calcJY0();



       break;
   default:
       break;
   }

}

//---------------------------------------------------------------------------
// Вычисление полной массы
double QBrlsBody::calcMass()
{
    double valSum = 0.;
    for (int i =0; i < mQuantPlots-1; ++i)
    {
      valSum +=  marrPlots[i].cal_LinDensityMass() * (marrBeginPlot[i + 1] - marrBeginPlot[i]) * 2.;
    }
    valSum +=  marrPlots[mQuantPlots -1].cal_LinDensityMass() * (mLength/ 2. - marrBeginPlot[mQuantPlots -1]) * 2.;
    return valSum;
}
//----------------------------------------------

// вычисление момента инерции относительно точки вращения
//
double QBrlsBody::calcJY0()
{
    double valSum = 0.;
    for (int i =0; i < mQuantPlots-1; ++i)
    {
      valSum +=  marrPlots[i].calJY0(marrBeginPlot[i], marrBeginPlot[i + 1]) * 2.;

    }

    valSum +=  marrPlots[mQuantPlots -1].calJY0(marrBeginPlot[mQuantPlots -1], mLength/2.) * 2.;

}

//-----------------------------------------------------------------
// создание граифков участков БРЛС
void QBrlsBody::createBodyGraph(wchar_t *wchOutPutFold)
{


    for (int i = 0; i <mQuantPlots; ++i )
    {
                wchar_t wchFoldCur[500] = {0};
                wcscpy(wchFoldCur, wchOutPutFold);
                wcscat(wchFoldCur, L"\\Plot_");

                wchar_t wchtemp[20] ={0};
                swprintf(wchtemp, 20, L"%i", i +1);
                wcscat(wchFoldCur, wchtemp);
                _wmkdir(wchFoldCur);
                marrPlots[i].createGraph(wchFoldCur);


    }

}
//-----------------------------------------------------------------

void QBrlsBody::createSideView(const double scale_x,const double  scaley,wchar_t *wchOutPutFold)
{
    wchar_t wchFold[500] = {0};
    wcscpy(wchFold, wchOutPutFold);
    wcscat(wchFold, L"\\Side_View");
     _wmkdir(wchFold);
    for (int i = 0; i <mQuantPlots; ++i )
    {
        wchar_t wchFoldCur[500] = {0};
        wcscpy(wchFoldCur, wchFold);
        wcscat(wchFoldCur, L"\\SideViewPlot_");

        wchar_t wchtemp[20] ={0};
        swprintf(wchtemp, 20, L"%i", i +1);
        wcscat(wchFoldCur, wchtemp);
        _wmkdir(wchFoldCur);



                double a = marrBeginPlot[i];
                double b = 0.;
                if ( i == (mQuantPlots -1))
                {
                 b = mLength/2.;
                }
                else
                {
                  b =  marrBeginPlot[i +1];
                }
                marrPlots[i].createSideView(a, b, scale_x, scaley, wchFoldCur);


    }
}
//----------------------------------------------
void QBrlsBody::analysisWind(const double VAlOmega, const double VAlAccellerOmega, const double VAlWindV
         ,const double VAlStep,  double *pvalFSum
       , double *pvalMSum,double *arrBuff, int *pquantDoneSteps)
{

  *pquantDoneSteps = 0;
  const int NUmStep = int(mLength/2./VAlStep);
  double valJY0 = calcJY0();

  double valMWanted = VAlAccellerOmega * valJY0;
  *pvalFSum = calcIqu(VAlOmega,    VAlWindV
                       , VAlStep, mLength/2.);
  *pvalMSum = calcKqu(VAlOmega,  VAlAccellerOmega,  VAlWindV
                      , VAlStep,  mLength/2.);
  *pvalMSum = *pvalMSum-valMWanted;

  int iPlotNumCur = 0;


  // вычисление коэффицмента, входящего в SigmaZ для каждого отдельного подсечения
  // участка балки с номером 0
  double arrSigmaCoeff[QUANT_PLOTS] = {0.};
  marrPlots[0].calcCoeffMaxSigma_ZX(arrSigmaCoeff );

  for (int i = 0; i < NUmStep; ++i)
  {
      if( i == (NUmStep - 1))
      {
          int ddj = 0;
      }
      double valZCur = (double(i))*VAlStep;
      if ((iPlotNumCur < (mQuantPlots -1))&&(valZCur >= (marrBeginPlot[iPlotNumCur + 1]+ 0.00001)))
      {
          iPlotNumCur++;
          // вычисление коэффицмента, входящего в SigmaZ для каждого отдельного подсечения
          // участка балки с номером iPlotNumCur
          memset(arrSigmaCoeff, 0, 4 * sizeof(double));
          marrPlots[iPlotNumCur].calcCoeffMaxSigma_ZX( arrSigmaCoeff );
      }

      double valSectionR = calcRqu(VAlOmega,  VAlAccellerOmega,  VAlWindV
                                   , VAlStep,  valZCur);

      double valSectionR0 = calcRqu(VAlOmega,  VAlAccellerOmega,  VAlWindV
                                   , VAlStep,  mLength/2.);
      double valSectionM = valMWanted +valSectionR0 - valSectionR;



    arrBuff[i * QUantColsReport_Wind   ] = valZCur;
    arrBuff[i * QUantColsReport_Wind +1] = calc_Wind_qu( VAlOmega,  VAlWindV, valZCur);
    arrBuff[i * QUantColsReport_Wind +2]
            = -calcIqu( VAlOmega,   VAlWindV, VAlStep,   mLength/2.)
            +calcIqu( VAlOmega,   VAlWindV , VAlStep,   valZCur);
    arrBuff[i * QUantColsReport_Wind +3] = valSectionM;
    for (int j =0; j < marrPlots[iPlotNumCur].mquantPlg; ++j)
    {
      arrBuff[i * QUantColsReport_Wind +4 + j] =  fabs(valSectionM)
              *  arrSigmaCoeff[j]/ marrPlots[iPlotNumCur].marrSigma[j];
    }

    (*pquantDoneSteps)++;
  }
}

//--------------------------------------------------
//----------------------------------------------
void QBrlsBody::analysisIce(const double VAlPIce,const double VAlStep, double *pvalFSum
       , double *pvalMSum,double *arrBuff, int *pquantDoneSteps)
{

  *pquantDoneSteps = 0;
  const int NUmStep = int(mLength/2./VAlStep);

  marrPlots[0].marrPlgProfile[0].calcBoundBox();
  double val_w = marrPlots[0].marrPlgProfile[0].Box[2] - marrPlots[0].marrPlgProfile[0].Box[0];

  *pvalFSum = mLength/2. * VAlPIce * val_w;
  *pvalMSum = VAlPIce * val_w * mLength* mLength/8. ;


  int iPlotNumCur = 0;


  // вычисление коэффицмента, входящего в SigmaZ для каждого отдельного подсечения
  // участка балки с номером 0
  double arrSigmaCoeff[QUANT_PLOTS] = {0.},arrSigmaCoeff1[QUANT_PLOTS] = {0.};
  marrPlots[0].calcCoeffMaxSigma_ZY(arrSigmaCoeff );


  for (int i = 0; i < NUmStep; ++i)
  {
      double valZCur = (double(i))*VAlStep;
      if ((iPlotNumCur < (mQuantPlots -1))&&(valZCur >= (marrBeginPlot[iPlotNumCur + 1]+ 0.00001)))
      {
          iPlotNumCur++;
          // вычисление коэффицмента, входящего в SigmaZ для каждого отдельного подсечения
          // участка балки с номером iPlotNumCur
          memset(arrSigmaCoeff, 0, 4 * sizeof(double));
          marrPlots[iPlotNumCur].calcCoeffMaxSigma_ZY( arrSigmaCoeff );
      }


      double valSectionM =  VAlPIce * val_w/2.
              * (mLength/2. - valZCur)* (mLength/2. - valZCur);



    arrBuff[i * QUantColsReport_Wind   ] = valZCur;
    arrBuff[i * QUantColsReport_Wind +1] = -VAlPIce * val_w;
    arrBuff[i * QUantColsReport_Wind +2]
            = VAlPIce * val_w * (mLength/2. - valZCur);
    arrBuff[i * QUantColsReport_Wind +3] = valSectionM;
    for (int j =0; j < marrPlots[iPlotNumCur].mquantPlg; ++j)
    {
      arrBuff[i * QUantColsReport_Wind +4 + j] =  fabs(valSectionM)
              *  arrSigmaCoeff[j]/ marrPlots[iPlotNumCur].marrSigma[j];
    }

    (*pquantDoneSteps)++;
  }
}

//--------------------------------------------------
//
double QBrlsBody::calcIqu(const double VAlOmega,  const double VAlWindV
                       ,const double VAlStep,  const double x)
{
    double temp = VAlOmega * VAlOmega * x * x + 3. * VAlOmega * x * VAlWindV
            + 3. *  VAlWindV*  VAlWindV;
    return -temp * VAL_BRLS_Cx  * mHeight * ATM_RoN0 * x / 6.;
}

//--------------------------------------------------
//
double QBrlsBody::calcKqu(const double VAlOmega, const double VAlAccellerOmega, const double VAlWindV
                       ,const double VAlStep,  const double x)
{
    double temp = VAlOmega * VAlOmega * x * x/ 4. + 2./3. * VAlOmega * x * VAlWindV
            +   VAlWindV*  VAlWindV/ 2.;
    return -temp * VAL_BRLS_Cx  * mHeight * ATM_RoN0 * x * x / 2.;
}
//--------------------------------------------------
//
double QBrlsBody::calcRqu(const double VAlOmega, const double VAlAccellerOmega, const double VAlWindV
                       ,const double VAlStep,  const double x)
{
    double temp = VAlOmega * VAlOmega * x * x/ 4. +  VAlOmega  * VAlWindV * x
            +   3./2. * VAlWindV*  VAlWindV;
    return -temp * VAL_BRLS_Cx  * mHeight * ATM_RoN0 * x * x / 6.
            - x * calcIqu(VAlOmega,   VAlWindV, VAlStep, mLength/2.)   ;
}
//----------------------------------------------------------
//вычисление распределнной нагрузки при обдувке ветром
double QBrlsBody::calc_Wind_qu(const double VAlOmega, const double VAlWindV, const double VAlz)
{
   double temp = -(VAlOmega * VAlz + VAlWindV) * (VAlOmega * VAlz + VAlWindV)
           * VAL_BRLS_Cx  * mHeight * ATM_RoN0 / 2.;
   return temp;
}

//-------------------------------------------
// вычисление скорости града
double QBrlsBody::calc_Hail_Mass(const double VAlDiam)
{
    // удельная плотность льда
    double valIceRo = 925;
    double valHailM = 4./ 3. * M_PI *VAlDiam *VAlDiam * VAlDiam / 8. * valIceRo;

    return valHailM;
}
//-------------------------------------------
// вычисление скорости града
double QBrlsBody::calc_Hail_V(const double VAlDiam)
{
    // удельная плотность льда
    double valIceRo = 925;
    double valHailM = calc_Hail_Mass(VAlDiam);
    // 1.24 - плотность воздуха
    // 0.47 - Cx
    return sqrt(8. * valHailM * 9.81 / 0.47 / M_PI /VAlDiam /VAlDiam/1.24);
}
//-------------------------------------------
// вычисление эффективной высоты падения  града
double QBrlsBody::calc_Hail_EffectH(const double VAlDiam)
{
   double valHailV = calc_Hail_V( VAlDiam);
    return valHailV * valHailV / 9.81/ 2.;
}
//----------------------------------------------

//----------------------------------------------

void QBrlsBody::analysisHail_(const double VAlDiam,const double VAlHailZ0,const double VAlStep, double *pvalFSum
       , double *pvalMSum,double *arrBuff, int *pquantDoneSteps)
{


    double valHailM = calc_Hail_Mass(VAlDiam);
    const double VAlP = valHailM * 9.81;
    double valEndProgib =0., valEndDerivProgib=0.;
    calcProgib(VAlHailZ0,VAlHailZ0, VAlP,&valEndProgib, &valEndDerivProgib);





    double valHail_EffectH =  calc_Hail_EffectH( VAlDiam);
    double valBlowCoeff = 1. + sqrt(1. + 2. * valHail_EffectH/fabs(valEndProgib) );

    const double VAlP0 = VAlP * valBlowCoeff;


    *pvalFSum = VAlP0;
    *pvalMSum = VAlP0 * VAlHailZ0;
  *pquantDoneSteps = 0;
  const int NUmStep = int(mLength/2./VAlStep);

  int iPlotNumCur = 0;


  // вычисление коэффицмента, входящего в SigmaZ для каждого отдельного подсечения
  // участка балки с номером 0
  double arrSigmaCoeff[QUANT_PLOTS] = {0.};
  marrPlots[0].calcCoeffMaxSigma_ZY(arrSigmaCoeff );
  double valOmega1 = 0., valOmega2 =0.;
  for (int i = 0; i < NUmStep; ++i)
  {
      double valZCur = (double(i))*VAlStep;
      if ((iPlotNumCur < (mQuantPlots -1))&&(valZCur >= (marrBeginPlot[iPlotNumCur + 1]+ 0.00001)))
      {
          iPlotNumCur++;
          // вычисление коэффицмента, входящего в SigmaZ для каждого отдельного подсечения
          // участка балки с номером iPlotNumCur
          memset(arrSigmaCoeff, 0, 4 * sizeof(double));
          marrPlots[iPlotNumCur].calcCoeffMaxSigma_ZY( arrSigmaCoeff );
      }


      double valSectionM = 0.;
      if (valZCur <= VAlHailZ0)
      {
         valSectionM = VAlP0 * (VAlHailZ0 -valZCur);
      }

      calcProgib(VAlHailZ0,valZCur, VAlP0,&valOmega1, &valOmega2);


    arrBuff[i * QUantColsReport_Wind   ] = valZCur;
    arrBuff[i * QUantColsReport_Wind +1] = valSectionM;
    arrBuff[i * QUantColsReport_Wind +2] = valOmega2;
    arrBuff[i * QUantColsReport_Wind +3] = valOmega1;
    for (int j =0; j < marrPlots[iPlotNumCur].mquantPlg; ++j)
    {
      arrBuff[i * QUantColsReport_Wind +4 + j] =  fabs(valSectionM)
              *  arrSigmaCoeff[j]/ marrPlots[iPlotNumCur].marrSigma[j];
    }

    (*pquantDoneSteps)++;
  }
int iiu=0;
}

//-------------------------------------------
void QBrlsBody::calcProgib(const double VAlHailZ0,const double VAlZCur, const double VAlP
                           ,double *pvalProgib, double *pvalDerivProgib)
{

    if (VAlZCur <= VAlHailZ0)
    {
        calcProgibInside( VAlHailZ0, VAlZCur, VAlP
                         ,pvalProgib, pvalDerivProgib);
    }
    else
    {
        calcProgibInside( VAlHailZ0, VAlHailZ0, VAlP,pvalProgib, pvalDerivProgib);
        *pvalProgib -= (*pvalDerivProgib)* (VAlHailZ0 - VAlZCur);
    }

}
//---------------------------
void QBrlsBody::calcProgibInside(const double VAlHailZ0,const double VAlZCur, const double VAlP
                            ,double *pvalProgib, double *pvalDerivProgib)
{
    *pvalProgib = 0.;
    *pvalDerivProgib =0.;
    int iIntervalNum = calcIntervalNum(VAlZCur);

    double arr_k[QUANT_PLOTS] = {0.};
    for (int i = 0;i < iIntervalNum; ++i)
    {

      double valRo1 = marrPlots[i].calcCurveRad_PlaneZY(  1.);
      arr_k[i] = VAlP / valRo1;
    }

    double arrZ[QUANT_PLOTS +1] = {0.};
    memcpy(arrZ, marrBeginPlot,mQuantPlots *sizeof(double) );
    arrZ[iIntervalNum]= VAlZCur;
    for (int m = 0;m < iIntervalNum; ++m)
    {

        *pvalProgib += (
        arr_k[m]/2.* VAlHailZ0*(arrZ[m+1] * arrZ[m+1] - arrZ[m]* arrZ[m])
        -arr_k[m]/6.*(arrZ[m+1] * arrZ[m+1]* arrZ[m+1] -arrZ[m] * arrZ[m]* arrZ[m])
        +((*pvalDerivProgib) -arr_k[m]* VAlHailZ0*arrZ[m]
          +arr_k[m]/2.*arrZ[m]*arrZ[m])
                *(arrZ[m+1] -arrZ[m])
                );


        *pvalDerivProgib += (arr_k[m]* VAlHailZ0 *(arrZ[m +1] - arrZ[m])
    -arr_k[m] / 2. *(arrZ[m +1] * arrZ[m+1] -arrZ[m] * arrZ[m]));


    }


}

//-------------------------------------------
// вычисление номера интервала (плота)(номера с 1)
// которому принадлежит точка VAlHailZ0
int QBrlsBody::calcIntervalNum(const double VAlHailZ0)
{
    int ireturn =-1;
    ///
    double arrPlots[QUANT_PLOTS +1] = {0.};
    memcpy(arrPlots, marrBeginPlot, mQuantPlots * sizeof(double));
    arrPlots[mQuantPlots] = mLength/2.;


    for (int i =0; i < mQuantPlots;++i)
    {
       if ((arrPlots[i] - VAlHailZ0)* (arrPlots[i + 1] - VAlHailZ0) <= 0.)
       {
           ireturn = i +1;
           break;
       }
    }
    return ireturn;
}


