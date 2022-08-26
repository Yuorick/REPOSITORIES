#include "BeamPlot.h"
#include <string.h>
#include "UrPointXY.h"
#include "MatrixProccess.h"
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"


// исходные данные для полигона кожуха для фурке 2540
extern const int NumPartFurkCaske_2540 = 2;
extern const int NumPointsFurk_Caske = 78 ;

extern const double arrPointsCaske_Furk_2540[78 * 2] =
 {
     0.087680 , 0.162611
    , 0.098711 , 0.155737
    , 0.120799 , 0.153787
    , 0.138323 , 0.148001
    , 0.154691 , 0.137541
    , 0.164598 , 0.126905
    , 0.174280 , 0.111617
    , 0.180557 , 0.095936
    , 0.182285 , 0.075799
    , 0.178443 , 0.054908
    , 0.171053 , 0.039518
    , 0.156861 , 0.023033
    , 0.142791 , 0.014311
    , 0.130804 , 0.009365
    , 0.114132 , 0.005548
    , 0.096128 , 0.003318
    , 0.090052 , 0.000000
    , -0.090000 , 0.000000
    , -0.103528 , 0.007434
    , -0.112227 , 0.020573
    , -0.116945 , 0.035981
    , -0.120019 , 0.060431
    , -0.119975 , 0.075689
    , -0.114850 , 0.090668
    , -0.104919 , 0.103636
    , -0.082335 , 0.119703
    , -0.057216 , 0.133658
    , -0.020820 , 0.145244
    , -0.001426 , 0.148761
    , 0.013753 , 0.151513
    , 0.039104 , 0.155850
    , 0.057402 , 0.156603
    , 0.064115 , 0.161421
    , 0.087680 , 0.162611
    , 0.094685 , 0.005959
    , 0.095027 , 0.006119
    , 0.095387 , 0.006235
    , 0.095758 , 0.006305
    , 0.113610 , 0.008516
    , 0.129888 , 0.012243
    , 0.141415 , 0.017000
    , 0.154882 , 0.025347
    , 0.168513 , 0.041180
    , 0.175555 , 0.055848
    , 0.179251 , 0.075946
    , 0.177596 , 0.095233
    , 0.171589 , 0.110241
    , 0.162205 , 0.125058
    , 0.152748 , 0.135210
    , 0.137020 , 0.145261
    , 0.120186 , 0.150819
    , 0.098446 , 0.152739
    , 0.097982 , 0.152817
    , 0.097536 , 0.152966
    , 0.097119 , 0.153183
    , 0.086889 , 0.159557
    , 0.065151 , 0.158460
    , 0.059157 , 0.154157
    , 0.058784 , 0.153929
    , 0.058382 , 0.153757
    , 0.057960 , 0.153645
    , 0.057526 , 0.153595
    , 0.039421 , 0.152851
    , 0.014276 , 0.148549
    , -0.000889 , 0.145799
    , -0.020092 , 0.142317
    , -0.056017 , 0.130881
    , -0.080726 , 0.117153
    , -0.102805 , 0.101446
    , -0.112162 , 0.089228
    , -0.116967 , 0.075184
    , -0.117008 , 0.060615
    , -0.113991 , 0.036614
    , -0.109477 , 0.021871
    , -0.101428 , 0.009715
    , -0.089227 , 0.003010
    , 0.089284 , 0.003010
    , 0.094685 , 0.005959
   };


extern const int PartsCaske_Furk_2540[] = {0, 34 };
///

// толстый текстолит
extern const int NumPartTextFurk_2540 = 1;
extern const int NumPointsTextFurk_2540 = 5 ;

extern const double arrPoints_TextFurk_2540[] =
 {  0.089, 0.003
   ,0.089, 0.007
   ,-0.089, 0.007
   ,-0.089, 0.003
   ,0.089, 0.003
   };
extern const int Parts_TextFurk_2540[] = {0 };
//

// алюминиевый тавр на антенне
extern const int NumPartAlFurk_2540 = 1;
extern const int NumPointsAlFurk_2540 = 9 ;

extern const double arrPoints_AlFurk_2540[9 * 2] =
 {
    0.089, 0.
    ,-0.089, 0.
    ,-0.089, -0.0055
    ,-0.086, -0.0055
    ,-0.086, -0.003
    ,0.086, -0.003
    ,0.086, -0.0055
    ,0.089, -0.0055
    ,0.089, 0.

   };
extern const int Parts_AlFurk_2540[] = {0 };
//---------------------------------------------------------------------------
// алюминиевая пластина на приводе
extern const int NumPartDrvFurk_2540 = 1;
extern const int NumPointsDrvFurk_2540 = 5 ;

extern const double arrPoints_DrvFurk_2540[] =
 {
    0.055, -0.003
    ,-0.055, -0.003
    ,-0.055, -0.018
    ,0.055, -0.018
    ,0.055, -0.003

   };
extern const int Parts_DrvFurk_2540[] = {0 };
/*
// алюминиевый тавр на антенне
extern const int NumPartAlFurk_2540 = 1;
extern const int NumPointsAlFurk_2540 = 5 ;

extern const double arrPoints_AlFurk_2540[9 * 2] =
 {
    0.00001, -0.003
    ,-0.00001, -0.003
    ,-0.00001, -0.018
    ,0.00001, -0.018
    ,0.00001, -0.003

   };
extern const int Parts_AlFurk_2540[] = {0 };
//---------------------------------------------------------------------------
// алюминиевая пластина на приводе
extern const int NumPartDrvFurk_2540 = 1;
extern const int NumPointsDrvFurk_2540 = 5 ;

extern const double arrPoints_DrvFurk_2540[] =
 {
    0.00001, -0.003
    ,-0.00001, -0.003
    ,-0.00001, -0.018
    ,0.00001, -0.018
    ,0.00001, -0.003

   };
extern const int Parts_DrvFurk_2540[] = {0 };
*/
//---------------------------------------------------------------------------
// П.Г. Королев, Справочник Сопротивление материалов
// Издателдьское объединение "Вища школа" Головное издательство
// Киев, 1974

// В.З. Васильев.Краткий курс сопротивления материалов с основами теории упругости.
// Санкт-Петербург, Иван Федоров, 2001 (стр 250)

// уд. плотность дюралюминия кг/(м*м*м)
extern const double VAL_DURAL_DENS = 2800;
// уд. плотность текстолита кг/(м*м*м)
extern const double VAL_TEXTOLIT_DENS = 1300;

// макс норм напряжение дюралюминия 80 - 150 Па
extern const double VAL_DURAL_SIGMA = 120 * 1000000.;
// макс норм напряжение текстолита 30 - 45 Па
extern const double VAL_TEXTOLIT_SIGMA = 35 * 1000000.;

// модуль упругости дюралюминия ГПа (взято из  В.З. Васильев.)
extern const double VAL_DURAL_E = 71. * 1000000000.;
//  модуль упругости текстолита 6 - 10 МПа
extern const double VAL_TEXTOLIT_E = 8. * 1000000000.;


QBeamPlot:: QBeamPlot()
{        
        mquantPlg = MAX_QUANT_PROFILES;
        for (int i = 0; i < MAX_QUANT_PROFILES; ++i)
        {
          marrPlgProfile[i] =  TURPolygon() ;
        }
        memset(marrPseudoCentre, 0, 2 * sizeof(double));
        memset(marrMtrxPseudoInertia, 0, 4 * sizeof(double));
        mcs = 0.;
}

//-------------------------------------------------
  // конструктор копирования
  QBeamPlot ::QBeamPlot (const QBeamPlot &R)
  {
    mquantPlg = R.mquantPlg;

    memcpy( marrPlgProfile,R.marrPlgProfile, R.mquantPlg  * sizeof(TURPolygon));

    memcpy( marrRo,R.marrRo, R.mquantPlg  * sizeof(double));

    memcpy( marrSigma,R.marrSigma, R.mquantPlg  * sizeof(double));

    memcpy( marrE,R.marrE, R.mquantPlg  * sizeof(double));

    memcpy( marrPseudoCentre,R.marrPseudoCentre, 2  * sizeof(double));

    memcpy( marrMtrxPseudoInertia,R.marrMtrxPseudoInertia, 4  * sizeof(double));

    mcs = R.mcs;


  }
//--------------------------------------------------------
  // оператор присваивания
QBeamPlot &QBeamPlot::operator=(const QBeamPlot  &R)
  {
    mquantPlg = R.mquantPlg;

    memcpy( marrPlgProfile,R.marrPlgProfile, R.mquantPlg  * sizeof(TURPolygon));

    memcpy( marrRo,R.marrRo, R.mquantPlg  * sizeof(double));

    memcpy( marrSigma,R.marrSigma, R.mquantPlg  * sizeof(double));

    memcpy( marrE,R.marrE, R.mquantPlg  * sizeof(double));

    memcpy( marrPseudoCentre,R.marrPseudoCentre, 2  * sizeof(double));

    memcpy( marrMtrxPseudoInertia,R.marrMtrxPseudoInertia, 4  * sizeof(double));

    mcs = R.mcs;


     return *this ;
 }
//---------------------------------------------------------

QBeamPlot::QBeamPlot(  const int quantPlg, const TURPolygon *arrPlgProfile, const double *arrRo
                       , const double *arrSigma, const double *arrE, const double cs)
{

   mquantPlg =  quantPlg;
   memcpy( marrPlgProfile,arrPlgProfile, quantPlg  * sizeof(TURPolygon));
   memcpy( marrRo,arrRo, mquantPlg  * sizeof(double));
   memcpy( marrSigma,arrSigma, mquantPlg  * sizeof(double));
   memcpy( marrE,arrE, mquantPlg  * sizeof(double));
   calcPseudoInertiaMtrx( marrPseudoCentre, marrMtrxPseudoInertia);
   mcs = cs;


}
//---------------------------------------------------------
/*
QBeamPlot::QBeamPlot(const int INumPlotFurke2540)
{
    mquantPlg = 4 - INumPlotFurke2540;


 switch(INumPlotFurke2540)
 {
 case 0:
     // кожух
     marrPlgProfile[0] = TURPolygon ( NumPartFurkCaske_2540,NumPointsFurk_Caske,PartsCaske_Furk_2540
                                       ,arrPointsCaske_Furk_2540);
    marrRo [0]= VAL_TEXTOLIT_DENS;
     marrSigma [0] = VAL_TEXTOLIT_SIGMA;
     marrE [0] = VAL_TEXTOLIT_E;

     // толстый текстолит
     marrPlgProfile[1] = TURPolygon ( NumPartTextFurk_2540,NumPointsTextFurk_2540,Parts_TextFurk_2540
                                      ,arrPoints_TextFurk_2540);
     marrRo [1]= VAL_TEXTOLIT_DENS;
     marrSigma [1] = VAL_TEXTOLIT_SIGMA;
     marrE [1] = VAL_TEXTOLIT_E;
     // алюминиевый тавр на антенне
     marrPlgProfile[2] = TURPolygon ( NumPartAlFurk_2540,NumPointsAlFurk_2540,Parts_AlFurk_2540
                                      ,arrPoints_AlFurk_2540);
     marrRo [2]= VAL_DURAL_DENS;
     marrSigma [2] = VAL_DURAL_SIGMA;
     marrE [2] = VAL_DURAL_E;
     // алюминиевая пластина на приводе
     marrPlgProfile[3] = TURPolygon ( NumPartDrvFurk_2540,NumPointsDrvFurk_2540,Parts_DrvFurk_2540
                                      ,arrPoints_DrvFurk_2540);
     marrRo [3]= VAL_DURAL_DENS;
     marrSigma [3] = VAL_DURAL_SIGMA;
     marrE [3] = VAL_DURAL_E;
    break;
 case 1:
     // кожух
  marrPlgProfile[0] = TURPolygon ( NumPartFurkCaske_2540,NumPointsFurk_Caske,PartsCaske_Furk_2540
                                       ,arrPointsCaske_Furk_2540);
     marrRo [0]= VAL_TEXTOLIT_DENS;
     marrSigma [0] = VAL_TEXTOLIT_SIGMA;
     marrE [0] = VAL_TEXTOLIT_E;

     // толстый текстолит
     marrPlgProfile[1] = TURPolygon ( NumPartTextFurk_2540,NumPointsTextFurk_2540,Parts_TextFurk_2540
                                      ,arrPoints_TextFurk_2540);
     marrRo [1]= VAL_TEXTOLIT_DENS;
     marrSigma [1] = VAL_TEXTOLIT_SIGMA;
     marrE [1] = VAL_TEXTOLIT_E;
     // алюминиевый тавр на антенне
     marrPlgProfile[2] = TURPolygon ( NumPartAlFurk_2540,NumPointsAlFurk_2540,Parts_AlFurk_2540
                                      ,arrPoints_AlFurk_2540);
     marrRo [2]= VAL_DURAL_DENS;
     marrSigma [2] = VAL_DURAL_SIGMA;
     marrE [2] = VAL_DURAL_E;
    break;
 case 2:
     // кожух
   marrPlgProfile[0] = TURPolygon ( NumPartFurkCaske_2540,NumPointsFurk_Caske,PartsCaske_Furk_2540
                                       ,arrPointsCaske_Furk_2540);
     marrRo [0]= VAL_TEXTOLIT_DENS;
     marrSigma [0] = VAL_TEXTOLIT_SIGMA;
     marrE [0] = VAL_TEXTOLIT_E;

     // толстый текстолит
     marrPlgProfile[1] = TURPolygon ( NumPartTextFurk_2540,NumPointsTextFurk_2540,Parts_TextFurk_2540
                                      ,arrPoints_TextFurk_2540);
     marrRo [1]= VAL_TEXTOLIT_DENS;
     marrSigma [1] = VAL_TEXTOLIT_SIGMA;
     marrE [1] = VAL_TEXTOLIT_E;
    break;
 case 3:
     // кожух
   marrPlgProfile[0] = TURPolygon ( NumPartFurkCaske_2540,NumPointsFurk_Caske,PartsCaske_Furk_2540
                                       ,arrPointsCaske_Furk_2540);
     marrRo [0]= VAL_TEXTOLIT_DENS;
     marrSigma [0] = VAL_TEXTOLIT_SIGMA;
     marrE [0] = VAL_TEXTOLIT_E;

    break;

 default:
    break;
 }

}
//-----------------------------------------------
*/
//---------------------------------------------------------
void QBeamPlot::createBeamPlotFurke2540(const double cs , int INumPlotFurke2540, QBeamPlot &BeamPlot)
{
    BeamPlot.mquantPlg = 4 - INumPlotFurke2540;
// отладка тестирование
   /*
    const int NumPart = 1;
    const int NumPoints = 5 ;

     const double arrPoints[] =
     {
        0.2, 0.1
        ,-0.2, 0.1
        ,-0.2, -0.1
        ,0.2, -0.1
        ,0.2, 0.1

       };
    const int Parts[] = {0 };
  TURPolygon rect( NumPart,NumPoints,Parts,arrPoints);

  double  valAng = M_PI/ 6.;
  const TURPointXY pntSdvig(1., 3.);
  TURPolygon rect1  ;
   rect1 =   rect.LinTransform( valAng ,  pntSdvig,1. ) ;

   rect1.WriteSetSHPFiles(L"D:\\АМЕТИСТ_2019\\БРЛС_ПРОЧНОСТЬ\\TEST\\rect1.shp", &rect1, 1);
*/
 BeamPlot.mcs = cs;
 switch(INumPlotFurke2540)
 {
 case 0:
     // кожух
     BeamPlot.marrPlgProfile[0] = TURPolygon ( NumPartFurkCaske_2540,NumPointsFurk_Caske,PartsCaske_Furk_2540
                                     ,arrPointsCaske_Furk_2540);
     //BeamPlot.marrPlgProfile[0] = rect1;
     BeamPlot.marrRo [0]= VAL_TEXTOLIT_DENS;
     BeamPlot.marrSigma [0] = VAL_TEXTOLIT_SIGMA;
     BeamPlot.marrE [0] = VAL_TEXTOLIT_E;

     // толстый текстолит
     BeamPlot.marrPlgProfile[1] = TURPolygon ( NumPartTextFurk_2540,NumPointsTextFurk_2540,Parts_TextFurk_2540
                                      ,arrPoints_TextFurk_2540);
     BeamPlot.marrRo [1]= VAL_TEXTOLIT_DENS;
     BeamPlot.marrSigma [1] = VAL_TEXTOLIT_SIGMA;
     BeamPlot.marrE [1] = VAL_TEXTOLIT_E;
     // алюминиевый тавр на антенне
     BeamPlot.marrPlgProfile[2] = TURPolygon ( NumPartAlFurk_2540,NumPointsAlFurk_2540,Parts_AlFurk_2540
                                      ,arrPoints_AlFurk_2540);
     BeamPlot.marrRo [2]= VAL_DURAL_DENS;
     BeamPlot.marrSigma [2] = VAL_DURAL_SIGMA;
     BeamPlot.marrE [2] = VAL_DURAL_E;
     // алюминиевая пластина на приводе
     BeamPlot.marrPlgProfile[3] = TURPolygon ( NumPartDrvFurk_2540,NumPointsDrvFurk_2540,Parts_DrvFurk_2540
                                      ,arrPoints_DrvFurk_2540);
     BeamPlot.marrRo [3]= VAL_DURAL_DENS;
     BeamPlot.marrSigma [3] = VAL_DURAL_SIGMA;
     BeamPlot.marrE [3] = VAL_DURAL_E;
    break;
 case 1:
     // кожух
     BeamPlot.marrPlgProfile[0] = TURPolygon ( NumPartFurkCaske_2540,NumPointsFurk_Caske,PartsCaske_Furk_2540
                                     ,arrPointsCaske_Furk_2540);
    // BeamPlot.marrPlgProfile[0] = rect1;
     BeamPlot.marrRo [0]= VAL_TEXTOLIT_DENS;
     BeamPlot.marrSigma [0] = VAL_TEXTOLIT_SIGMA;
     BeamPlot.marrE [0] = VAL_TEXTOLIT_E;

     // толстый текстолит
     BeamPlot.marrPlgProfile[1] = TURPolygon ( NumPartTextFurk_2540,NumPointsTextFurk_2540,Parts_TextFurk_2540
                                      ,arrPoints_TextFurk_2540);
     BeamPlot.marrRo [1]= VAL_TEXTOLIT_DENS;
     BeamPlot.marrSigma [1] = VAL_TEXTOLIT_SIGMA;
     BeamPlot.marrE [1] = VAL_TEXTOLIT_E;
     // алюминиевый тавр на антенне
     BeamPlot.marrPlgProfile[2] = TURPolygon ( NumPartAlFurk_2540,NumPointsAlFurk_2540,Parts_AlFurk_2540
                                      ,arrPoints_AlFurk_2540);
     BeamPlot.marrRo [2]= VAL_DURAL_DENS;
     BeamPlot.marrSigma [2] = VAL_DURAL_SIGMA;
     BeamPlot.marrE [2] = VAL_DURAL_E;
    break;
 case 2:
     // кожух
    BeamPlot.marrPlgProfile[0] = TURPolygon ( NumPartFurkCaske_2540,NumPointsFurk_Caske,PartsCaske_Furk_2540
                                        ,arrPointsCaske_Furk_2540);
     // BeamPlot.marrPlgProfile[0] = rect1;
     BeamPlot.marrRo [0]= VAL_TEXTOLIT_DENS;
     BeamPlot.marrSigma [0] = VAL_TEXTOLIT_SIGMA;
     BeamPlot.marrE [0] = VAL_TEXTOLIT_E;

     // толстый текстолит
     BeamPlot.marrPlgProfile[1] = TURPolygon ( NumPartTextFurk_2540,NumPointsTextFurk_2540,Parts_TextFurk_2540
                                      ,arrPoints_TextFurk_2540);
     BeamPlot.marrRo [1]= VAL_TEXTOLIT_DENS;
     BeamPlot.marrSigma [1] = VAL_TEXTOLIT_SIGMA;
     BeamPlot.marrE [1] = VAL_TEXTOLIT_E;
    break;
 case 3:
     // кожух
     BeamPlot.marrPlgProfile[0] = TURPolygon ( NumPartFurkCaske_2540,NumPointsFurk_Caske,PartsCaske_Furk_2540
                                     ,arrPointsCaske_Furk_2540);
     // BeamPlot.marrPlgProfile[0] = rect1;
     BeamPlot.marrRo [0]= VAL_TEXTOLIT_DENS;
     BeamPlot.marrSigma [0] = VAL_TEXTOLIT_SIGMA;
     BeamPlot.marrE [0] = VAL_TEXTOLIT_E;

    break;

 default:
    break;
 }
 BeamPlot.calcPseudoInertiaMtrx(BeamPlot.marrPseudoCentre, BeamPlot.marrMtrxPseudoInertia);
}
//-----------------------------------------------

double QBeamPlot::calcNeutral_Y()
{
  double sumUp = 0., sumDown = 0.;
  for (int i =0; i < mquantPlg; ++i)
  {
    sumUp  +=  marrE[i] * marrPlgProfile[i].calcSy(mcs);
    sumDown += marrE[i] * fabs(marrPlgProfile[i].calcVectSq());
  }
  return sumUp / sumDown;
}
//-----------------------------------------------

double QBeamPlot::calcNeutral_X()
{
  double sumUp = 0., sumDown = 0.;
  for (int i =0; i < mquantPlg; ++i)
  {
    sumUp  +=  marrE[i] * marrPlgProfile[i].calcSx(mcs);
    sumDown += marrE[i] * fabs(marrPlgProfile[i].calcVectSq());
  }
  return sumUp / sumDown;
}

//-----------------------------------------------
// mtrxPseudoInertia = SUM(X * XT * E)
// arrPseudoCentre = SUM(E * X)/ SUM( E )
void QBeamPlot::calcPseudoInertiaMtrx( double *arrPseudoCentre, double *mtrxPseudoInertia)
{
   memset (mtrxPseudoInertia, 0, 4 * sizeof(double));
   arrPseudoCentre[1] =  calcNeutral_X();
   arrPseudoCentre[0] =  calcNeutral_Y();
   const TURPointXY PNtPseudoCentre(arrPseudoCentre[0], arrPseudoCentre[1]);

   double mtrxInertia[4] = {0.}, arrT0[4] = {0.},arrT1[4] = {0.};
   for (int i =0; i < mquantPlg; ++i)
   {
      marrPlgProfile[i].calcInertiaMtrx(PNtPseudoCentre,  mcs, mtrxInertia) ;
      MatrxMultScalar(mtrxInertia, 2, 2, marrE[i],arrT0);
      MtrxSumMatrx(arrT0, mtrxPseudoInertia,2, 2, arrT1) ;
      memcpy(mtrxPseudoInertia, arrT1, 4 * sizeof(double));
   }


}
//-----------------------------------------------
// нахождение радиуса кривизны в плоскости OXZ
// нагрузка действует в плоскости OXZ
double QBeamPlot::calcCurveRad_PlaneZX(const double VAlMomY)
{
  return marrMtrxPseudoInertia[0]/ VAlMomY ;
}

//-----------------------------------------------
// нахождение радиуса кривизны в плоскости OYZ
// нагрузка действует в плоскости OYZ
double QBeamPlot::calcCurveRad_PlaneZY( const double VAlMomX)
{

  return marrMtrxPseudoInertia[3]/ VAlMomX ;
}
//--------------------------------------------------
// вычисление плотности массы по длине
double QBeamPlot::cal_LinDensityMass()
{
    double sum = 0.;
    for (int i = 0; i < mquantPlg; ++i)
    {
       sum +=  fabs(marrPlgProfile[i].calcVectSq()) * marrRo[i];
    }
    return sum;
}
//--------------------------------------------------------------------
// вычиление момента инерции тела относительно орси вращения
//INPUT:
//cs - длина пикселя растра
//a - начало участка (расстояние от центра вращения)
//b - начало участка (расстояние от центра вращения)
double QBeamPlot::calJY0( const double a, const double b)
{

    double sum = 0.;
    for (int i = 0; i < mquantPlg; ++i)
    {
       double temp1 =  marrPlgProfile[i].calcJy(mcs, 0.)* marrRo[i] * (b - a);
       double temp2 =  marrRo[i] * fabs(marrPlgProfile[i].calcVectSq()) * (b * b * b - a * a * a)/3.;
       sum +=  (temp1 + temp2);
    }

    return sum;
}

//-----------------------------------
void QBeamPlot::createGraph(wchar_t *wchFold)
{
   for (int i = 0; i < mquantPlg; ++i)
   {
       wchar_t wchFileCur[500] = {0};
       wcscpy(wchFileCur, wchFold);
       wcscat(wchFileCur, L"\\Plg_");

       wchar_t wchtemp[20] ={0};
       swprintf(wchtemp, 20, L"%i", i +1);
       wcscat(wchFileCur, wchtemp);
       wcscat(wchFileCur, L".shp");
       marrPlgProfile[i].WriteSetSHPFiles(wchFileCur,&(marrPlgProfile[i]),1);

   }
   // оси  координат
   wchar_t wchAxesFileName0[300] ={0};
   wcscpy(  wchAxesFileName0,  wchFold);
   wcscat(wchAxesFileName0, L"\\AxesArr.shp");
   TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
   ,-10000., 10000.,30.) ;
   //double arrPseudoCentre[2] = {0.}, mtrxPseudoInertia[4] = {0.};
  // calcPseudoInertiaMtrx( arrPseudoCentre, mtrxPseudoInertia);

   double arrV[4] = {0.}, arrLamb[4] = {0.};
  //CalcProperVectors2(mtrxPseudoInertia,arrV , arrLamb)  ;
   CalcProperVectors2(marrMtrxPseudoInertia,arrV , arrLamb)  ;

   const double  valAng = atan2(arrV[2], arrV[0]);

  TURPointXY pointBeginX(-100.,0.), pointEndX(100.,0.),pointBeginY(0.,-100.), pointEndY(0., 100.);
  double valLength = 10.;
   TURPolyLine plnAxes = TURPolyLine::fncCreateAxes( pointBeginX,  pointEndX, pointBeginY,  pointEndY
                                         , valLength);
   TURPointXY pntSdvig(marrPseudoCentre[0], marrPseudoCentre[1]);
   TURPolyLine plnAxes1 = plnAxes.LinTransform( valAng ,  pntSdvig,1 ) ;

   wchar_t wchAxes1[500] ={0};
   wcscpy( wchAxes1, wchFold);
   wcscat( wchAxes1, L"\\Axes1.shp");

   plnAxes1.WriteSetSHPFiles(wchAxes1,&(plnAxes1),1);

}

//---------------------------------------------------------
void QBeamPlot::createSideView(const double a,const double  b
  ,const double scale_x,const double  scaley,wchar_t * wchFoldCur0)
{
    for (int i =0; i < mquantPlg; ++i)
    {
        wchar_t wchFile[500] = {0};
        wcscpy(wchFile, wchFoldCur0);
        wcscat(wchFile, L"\\SideViewPlg_");

        wchar_t wchtemp[20] ={0};
        swprintf(wchtemp, 20, L"%i", i +1);
        wcscat(wchFile, wchtemp);
        wcscat(wchFile,L".shp");
        TURPolygon plgSideView(5);
        marrPlgProfile[i].calcBoundBox();
        plgSideView.Points[0] = TURPointXY(a *scale_x, marrPlgProfile[i].Box[1] * scaley);
        plgSideView.Points[1] = TURPointXY(b *scale_x, marrPlgProfile[i].Box[1] * scaley);
        plgSideView.Points[2] = TURPointXY(b *scale_x, marrPlgProfile[i].Box[3] * scaley);
        plgSideView.Points[3] = TURPointXY(a *scale_x, marrPlgProfile[i].Box[3] * scaley);
        plgSideView.Points[4] = plgSideView.Points[0];
        plgSideView.WriteSetSHPFiles(wchFile, &plgSideView, 1);


    }

}

//-----------------------------------------------
//вычисление массива относительных максимальнных нормальных напряжений
void QBeamPlot::calcCoeffMaxSigma_ZX( double *arrSigmaCoeff )
{    

    double valRo1 = calcCurveRad_PlaneZX( 1.);

    for (int i =0; i < mquantPlg; ++i)
    {
      marrPlgProfile[i].calcBoundBox();
      double valMaxX = marrPlgProfile[i].Box[2];
      double valMinX = marrPlgProfile[i].Box[0];
      double valCritX = max__(fabs(valMaxX - marrPseudoCentre[0]), fabs(valMinX - marrPseudoCentre[0]));
      arrSigmaCoeff[i] = valCritX * marrE[i] / valRo1;
    }

}//-----------------------------------------------
//вычисление массива относительных максимальнных нормальных напряжений
void QBeamPlot::calcCoeffMaxSigma_ZY( double *arrSigmaCoeff )
{


    double valRo1 = calcCurveRad_PlaneZY(1.);


    for (int i =0; i < mquantPlg; ++i)
    {
      marrPlgProfile[i].calcBoundBox();
      double valMaxY = marrPlgProfile[i].Box[3];
      double valMinY = marrPlgProfile[i].Box[1];
      double valCritY = max__(fabs(valMaxY - marrPseudoCentre[1]), fabs(valMinY - marrPseudoCentre[1]));
      arrSigmaCoeff[i] = valCritY * marrE[i] / valRo1;
    }

}


//-----------------------------------
double QBeamPlot::max__(const double a, const double b)
{
    return (a > b)?a:b;
}
