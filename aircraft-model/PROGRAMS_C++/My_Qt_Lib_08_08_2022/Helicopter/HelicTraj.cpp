#include "HelicTraj.h"

#include <string.h>
#include <QFile>
#include <QDebug>
#include <math.h>
#include <float.h>

#include "MatrixProccess.h"

#include "Comp.h"
#include "HelicConstants.h"
#include "GaussKruger.h"




extern const long  double constArrGearSets[];
extern const long double constHelicMass;


extern const int QUantActualTaskPoints0 ;
extern long double ARRCourse0[];
extern long double ARRV0 [] ;
extern long double ARRCoord0[];
extern enumTypeOfManeuvre ARRTypeMove0[] ;
extern long double ARRYaw0[];
extern long double ARRRadius0[];
extern int IARRNumSetOfGears0[];
extern long double  ARRCourseEnd0[];
extern long double   ARRTHovering0[];

extern long double ARR_dPsi_po_dt0[];


extern const int QUantActualTaskPoints1 ;
extern long double ARRCourse1[];
extern long double ARRV1 [] ;
extern long double ARRCoord1[];
extern enumTypeOfManeuvre ARRTypeMove1[] ;
extern long double ARRYaw1[];
extern long double ARRRadius1[];
extern int IARRNumSetOfGears1[];
extern long double  ARRCourseEnd1[];
extern long double   ARRTHovering1[];

extern long double ARR_dPsi_po_dt1[];



THelicTraj::THelicTraj()
{
    //  вертолет
    mHelic =THelic ();
    // атмосфера
    mEnvironment = TEnvironment ();
    // полетное задание
    mFlyTask =TFlyTask ();
    //

    // вектор координат нулевой (первой по счету) точки полетного задания в ГСК
    memset(marr_GSKxyz_Pnt0, 0, 3 * sizeof(long double));
    ///

    // долгота и широта нулевой (первой по счету) точки полетного задания в ГСК
    mLongitude = 0.; // долгота
    mLatitude = 0.; // широта

    memset(marrPhaseVect, 0, QUantCurNZSKVarsVS * sizeof(long double));
    mTimeCur = 0.;

    fillMagDevTab();

    // номер текущей части траектории
     mNumPart = 0;;

     // признак иницуиализации текущей части траектории
     mbPartInit = false;


    // матрица перехода из текущей НЗСК в НЗСК
     memset(marrTransf_CurNZSK_to_NZSK, 0, 9 * sizeof(double));
    // матрица перехода из  НЗСК в текущую НЗСК
     memset(marrTransf_NZSK_to_CurNZSK, 0, 9 * sizeof(double));
    // вектор координат начала НЗСК
     memset(marrO_CurNZSK, 0, 3 * sizeof(double));

     //memset(miarr, 0, MAX_POSSIBLE_TRAJPART_QUANT* sizeof(int));


    // участок траектории
    mpPartTraj = NULL;
    //
    mStepIntegr = 0.01;
}

//---------------------------------------------------------------------------

// оператор присваивания
THelicTraj THelicTraj::operator=(THelicTraj  R)
{
    // 1.
    mHelic = R.mHelic;

    //2.
    mEnvironment = R.mEnvironment;

    // 3.
    mFlyTask = R.mFlyTask;

    // 4.
    memcpy(marr_GSKxyz_Pnt0, R.marr_GSKxyz_Pnt0, 3 * sizeof(long double));

    // 5, 6 долгота и широта нулевой (первой по счету) точки полетного задания в ГСК
    mLongitude = R.mLongitude; // долгота
    mLatitude = R.mLatitude; // ш

    // 7.
    memcpy(marrPhaseVect, R.marrPhaseVect, QUantCurNZSKVarsVS * sizeof(long double));

    // 8.
    mTimeCur = R.mTimeCur;

    // 9.
    mqvectTabMagDev = R.mqvectTabMagDev;
    // 10.
    mQuantRowsTabMagDev = R.mQuantRowsTabMagDev;
    //11.
    mQuantColsTabMagDev  = R.mQuantColsTabMagDev;


    // 12 номер текущей части траектории
    mNumPart = R.mNumPart;

     // 13 массив номеров матриц передаточных чисел для участков траектории
    //memcpy(miarr, R.miarr, MAX_POSSIBLE_TRAJPART_QUANT * sizeof(int));

     // 14 признак иницуиализации текущей части траектории
    mbPartInit = R.mbPartInit;


    // 15 матрица перехода из текущей НЗСК в НЗСК
    memcpy(marrTransf_CurNZSK_to_NZSK, R.marrTransf_CurNZSK_to_NZSK, 9 * sizeof(long double));

    // 16 матрица перехода из  НЗСК в текущую НЗСК
    memcpy(marrTransf_NZSK_to_CurNZSK, R.marrTransf_NZSK_to_CurNZSK, 9 * sizeof(long double));

    // 17 вектор координат начала НЗСК
    memcpy(marrO_CurNZSK, R.marrO_CurNZSK, 3 * sizeof(long double));

    // 18 участок траектории
   // mpPartTraj = R.mpPartTraj;;
    mTurnMove = R.mTurnMove;
    mLineMove = R.mLineMove;
    mHover = R.mHover;
    mRotating = R.mRotating;

    if (R.mpPartTraj == &(R.mLineMove))
    {
       mpPartTraj = &mLineMove;
    }
    if (R.mpPartTraj == &(R.mTurnMove))
    {
       mpPartTraj = &mTurnMove;
    }

    if (R.mpPartTraj == &(R.mHover))
    {
       mpPartTraj = &mHover;
    }
    if (R.mpPartTraj == &(R.mRotating))
    {
       mpPartTraj = &mRotating;
    }

    // 19
    mStepIntegr = R.mStepIntegr;


    return *this;
}

// конструктор копирования
THelicTraj::THelicTraj (const THelicTraj &R)
{
    // 1.
    mHelic = R.mHelic;

    //2.
    mEnvironment = R.mEnvironment;

    // 3.
    mFlyTask = R.mFlyTask;

    // 4.
    memcpy(marr_GSKxyz_Pnt0, R.marr_GSKxyz_Pnt0, 3 * sizeof(long double));

    // 5, 6 долгота и широта нулевой (первой по счету) точки полетного задания в ГСК
    mLongitude = R.mLongitude; // долгота
    mLatitude = R.mLatitude; // ш

    // 7.
    memcpy(marrPhaseVect, R.marrPhaseVect, QUantCurNZSKVarsVS * sizeof(long double));

    // 8.
    mTimeCur = R.mTimeCur;

    // 9.
    mqvectTabMagDev = R.mqvectTabMagDev;
    // 10.
    mQuantRowsTabMagDev = R.mQuantRowsTabMagDev;
    //11.
    mQuantColsTabMagDev  = R.mQuantColsTabMagDev;


    // 12 номер текущей части траектории
    mNumPart = R.mNumPart;

     // 13 массив номеров матриц передаточных чисел для участков траектории
    //memcpy(miarr, R.miarr, MAX_POSSIBLE_TRAJPART_QUANT * sizeof(int));

     // 14 признак иницуиализации текущей части траектории
    mbPartInit = R.mbPartInit;


    // 15 матрица перехода из текущей НЗСК в НЗСК
    memcpy(marrTransf_CurNZSK_to_NZSK, R.marrTransf_CurNZSK_to_NZSK, 9 * sizeof(long double));

    // 16 матрица перехода из  НЗСК в текущую НЗСК
    memcpy(marrTransf_NZSK_to_CurNZSK, R.marrTransf_NZSK_to_CurNZSK, 9 * sizeof(long double));

    // 17 вектор координат начала НЗСК
    memcpy(marrO_CurNZSK, R.marrO_CurNZSK, 3 * sizeof(long double));

    // 18 участок траектории
   // mpPartTraj = R.mpPartTraj;;
    mTurnMove = R.mTurnMove;
    mLineMove = R.mLineMove;
    mHover = R.mHover;
    mRotating = R.mRotating;

    if (R.mpPartTraj == &(R.mLineMove))
    {
       mpPartTraj = &mLineMove;
    }
    if (R.mpPartTraj == &(R.mTurnMove))
    {
       mpPartTraj = &mTurnMove;
    }

    if (R.mpPartTraj == &(R.mHover))
    {
       mpPartTraj = &mHover;
    }
    if (R.mpPartTraj == &(R.mRotating))
    {
       mpPartTraj = &mRotating;
    }

    // 19
    mStepIntegr = R.mStepIntegr;


}


//---------------------------------------------------------------------------

 // парам констр автономного создания - без использования интерфейса
// если ITYpeof_FLyTask = 0, то траектория на высоте 300 м в соответствии с РЛЭ
THelicTraj::THelicTraj(const TEnvironment Environment , const int ITYpeof_FLyTask
                       ,const bool bTypeofHelic, const bool bTypeofRotor)
{
    mHelic =  THelic(bTypeofHelic, bTypeofRotor, constHelicMass);

    mEnvironment =  Environment ;


 switch(ITYpeof_FLyTask)
 {
 case 0:
     mFlyTask = TFlyTask(QUantActualTaskPoints0, ARRCoord0, ARRV0,
                        ARRCourse0, ARRTypeMove0, ARRYaw0, ARRRadius0, IARRNumSetOfGears0
                         , ARRCourseEnd0,   ARRTHovering0, ARR_dPsi_po_dt0 );
 break;
 case 1:
     mFlyTask = TFlyTask(QUantActualTaskPoints1, ARRCoord1, ARRV1,
                        ARRCourse1, ARRTypeMove1, ARRYaw1, ARRRadius1, IARRNumSetOfGears1
                         , ARRCourseEnd1,    ARRTHovering1, ARR_dPsi_po_dt1 );
 break;

 default:
 break;

 }
    mStepIntegr = 0.01;

    mLatitude  = 55.666429;
    mLongitude = 37.947845;
   double valLatitude = static_cast<double>(mLatitude);
   double valLongitude = static_cast<double>(mLongitude);
   double t0 = valLatitude * M_PI/ 180.;
   double t1 = valLongitude* M_PI/ 180.;

   double valGaussX = 0., valGaussY = 0.;
  mhGeo2Gauss_( t0, t1, &valGaussX, &valGaussY );

   marr_GSKxyz_Pnt0[0] = valGaussX;
   marr_GSKxyz_Pnt0[1] = valGaussY;
   marr_GSKxyz_Pnt0[2] = mFlyTask.mparrPartTrajInputData[0].mTaskPointsZ.Y;

   mTimeCur = 0.;

   memset(marrPhaseVect, 0, QUantCurNZSKVarsVS * sizeof(long double));
   marrPhaseVect [0] = mFlyTask.mparrPartTrajInputData[0].mTaskPointsZ.X;
   marrPhaseVect [1] = mFlyTask.mparrPartTrajInputData[0].mTaskPointsZ.Y;// 15.10.18
   marrPhaseVect [2] = mFlyTask.mparrPartTrajInputData[0].mTaskPointsZ.Z;

   marrPhaseVect [3] = mFlyTask.mparrPartTrajInputData[0].mV * cos(mFlyTask.mparrPartTrajInputData[0].mCourse);
   marrPhaseVect [4] = 0.;
   marrPhaseVect [5] = mFlyTask.mparrPartTrajInputData[0].mV * sin(mFlyTask.mparrPartTrajInputData[0].mCourse);
   marrPhaseVect [11] = atan2l(-(mFlyTask.mparrPartTrajInputData[1].mTaskPointsZ.Z - mFlyTask.mparrPartTrajInputData[0].mTaskPointsZ.Z)
                 ,mFlyTask.mparrPartTrajInputData[1].mTaskPointsZ.X - mFlyTask.mparrPartTrajInputData[0].mTaskPointsZ.X)
           + mFlyTask.mparrPartTrajInputData[1].mYaw ;


   ///
   mNumPart = 1;

   initPartTraj(1);

   fillMagDevTab();

}



//--------------------------------------------------------------------------

void THelicTraj::initPartTraj(const int NUmPart)
{


   // mpPartTraj->mTimeCur = mTimeCur;
    mbPartInit = true;
    ///
    // пересчет направленгия ветра в ТНЗСК
    TEnvironment EnvironmentCur = mEnvironment;


    // это заготовка для точки центра вращения, котороая нужна для класса TTurnMove
    long double arrRotCentre_CurNZSK[3] = {0.};
    ///

   // номер набора передаточных чисел из файла "HelicConstants.cpp"
    int iNumOfGear= mFlyTask.mparrPartTrajInputData[NUmPart].mNumSetOfGears;
    ///



    switch(mFlyTask.mparrPartTrajInputData[NUmPart].menumType)
    {
    case (enLine):
        // mFlyTask.mparrPartTrajInputData[NUmPart].mCourse - это курсовой угол оси OX ТНЗСК

        // перевод вектора скорости ветра в ТНЗСК
        EnvironmentCur.mWind_Alf = mEnvironment.mWind_Alf - mFlyTask.mparrPartTrajInputData[NUmPart].mCourse;
        if (fabsl(EnvironmentCur.mWind_Alf) > M_PI)
         {
          EnvironmentCur.mWind_Alf -= SIGNUM(EnvironmentCur.mWind_Alf) * 2. * M_PI;
         }
        ///


        // матрица линейного преобразования координат из НЗСК в ТНЗСК
        create_transfMtrx_NZSK_to_CurNZSK(NUmPart);

        // матрица линейного преобразования координат из ТНЗСК в НЗСК
        create_transfMtrx_CurNZSK_to_NZSK(NUmPart);

        // вектор начала координат ТНЗСК в НЗСК
        create_vectO(NUmPart);

        // создание объекта класса TLineMove
        mLineMove = TLineMove( mHelic, EnvironmentCur, mFlyTask.mparrPartTrajInputData[NUmPart].mTaskPointsZ.Y
                    ,mFlyTask.mparrPartTrajInputData[NUmPart].mV, mFlyTask.mparrPartTrajInputData[NUmPart].mYaw, mTimeCur);


        // запоминание матрицы передаточыных чисел в классе
         memcpy(mLineMove.marrC0, &(constArrGearSets[48 *iNumOfGear]), 48 * sizeof(long double));
         mpPartTraj = &mLineMove;
        break;

    case(enTurn):
        // mFlyTask.mparrPartTrajInputData[NUmPart-1].mCourse - это курсовой угол оси OX ТНЗСК

        // перевод вектора скорости ветра в ТНЗСК
        EnvironmentCur.mWind_Alf = mEnvironment.mWind_Alf - mFlyTask.mparrPartTrajInputData[NUmPart -1].mCourse;
        if (fabsl(EnvironmentCur.mWind_Alf) > M_PI)
         {
          EnvironmentCur.mWind_Alf -= SIGNUM(EnvironmentCur.mWind_Alf) * 2. * M_PI;
         }
        ///

          // матрица линейного преобразования координат из НЗСК в ТНЗСК
        create_transfMtrx_NZSK_to_CurNZSK(NUmPart -1);

         // матрица линейного преобразования координат из ТНЗСК в НЗСК
        create_transfMtrx_CurNZSK_to_NZSK(NUmPart -1);

         // вектор начала координат ТНЗСК в НЗСК
        create_vectO(NUmPart -1);

        // координаты точки центра вращения в ТНЗСК
        arrRotCentre_CurNZSK[0] = 0.;
        arrRotCentre_CurNZSK[1] = mFlyTask.mparrPartTrajInputData[NUmPart].mTaskPointsZ.Y;
        arrRotCentre_CurNZSK[2] = mFlyTask.mparrPartTrajInputData[NUmPart].mRadius;
        ///

        //
        mTurnMove = TTurnMove( mHelic, EnvironmentCur, mFlyTask.mparrPartTrajInputData[NUmPart].mTaskPointsZ.Y
                    ,mFlyTask.mparrPartTrajInputData[NUmPart].mV, mFlyTask.mparrPartTrajInputData[NUmPart].mYaw
                   ,mFlyTask.mparrPartTrajInputData[NUmPart].mRadius
                               , arrRotCentre_CurNZSK, mTimeCur, mFlyTask.mparrPartTrajInputData[NUmPart].mCourseEnd);

        memcpy(mTurnMove.marrC0, &(constArrGearSets[48 *iNumOfGear]), 48 * sizeof(long double));// !!!!!!!!!!
        mpPartTraj = &mTurnMove;

        // курсовой угол движения ЛА на момент начала вращения
        mTurnMove.mCurrentCourseIntegrated = atan2l(marrPhaseVect[5], marrPhaseVect[3])-mFlyTask.mparrPartTrajInputData[NUmPart -1].mCourse;

     break;


    case (enSpotHovering):

        // mFlyTask.mparrPartTrajInputData[NUmPart-1].mCourse - это курсовой угол оси OX ТНЗСК

        // перевод вектора скорости ветра в ТНЗСК
        EnvironmentCur.mWind_Alf = mEnvironment.mWind_Alf - mFlyTask.mparrPartTrajInputData[NUmPart].mCourse;
        if (fabsl(EnvironmentCur.mWind_Alf) > M_PI)
         {
          EnvironmentCur.mWind_Alf -= SIGNUM(EnvironmentCur.mWind_Alf) * 2. * M_PI;
         }
        create_transfMtrx_NZSK_to_CurNZSK(NUmPart);
        create_transfMtrx_CurNZSK_to_NZSK(NUmPart);
        create_vectO(NUmPart);


        mHover = THover( mHelic, EnvironmentCur, mFlyTask.mparrPartTrajInputData[NUmPart].mTaskPointsZ.Y                    
                    , mFlyTask.mparrPartTrajInputData[NUmPart].mTHovering, mTimeCur);

         memcpy(mHover.marrC0, &(constArrGearSets[48 *iNumOfGear]), 48 * sizeof(long double));// !!!!!!!!!!
         mpPartTraj = &mHover;
        break;

    case (enRotating):

        mFlyTask.mparrPartTrajInputData[NUmPart].mCourse =0.;

        create_transfMtrx_NZSK_to_CurNZSK(NUmPart);
        create_transfMtrx_CurNZSK_to_NZSK(NUmPart);
        create_vectO(NUmPart);
        mRotating = TRotating( mHelic, mEnvironment, mFlyTask.mparrPartTrajInputData[NUmPart].mTaskPointsZ.Y
                    ,mFlyTask.mparrPartTrajInputData[NUmPart].m_dPsi_po_dt
                    , mFlyTask.mparrPartTrajInputData[NUmPart].mTHovering, mTimeCur, marrPhaseVect[11]);

         memcpy(mRotating.marrC0, &(constArrGearSets[48 *iNumOfGear]), 48 * sizeof(long double));// !!!!!!!!!!
         mpPartTraj = &mRotating;
        break;

    default:
        break;
    }

    memcpy(mpPartTraj->marrPhaseVect,marrPhaseVect, QUantCurNZSKVarsVS * sizeof(long double));
    ///

    // пересчет угла рыскания
    marrPhaseVect[11] = TPartHelicTraj::fnc_Minus_PI_Plus_PI( marrPhaseVect[11] );
    mpPartTraj->marrPhaseVect[11] = recalcPsi_NZSK_to_CurNZSK(marrPhaseVect[11]);

   ///  
    transf_xyzNZSK_to_xyzCurNZSK(marrPhaseVect, mpPartTraj->marrPhaseVect);
    MtrxMultMatrx(marrTransf_NZSK_to_CurNZSK,3, 3, &(marrPhaseVect[3]),1, &(mpPartTraj->marrPhaseVect[3])) ;

   mpPartTraj->mTimeCur = mTimeCur;

   mpPartTraj->doAirDensity();

}

//---------------------------------------------------------------------------
// пересчет угла рыскания из НЗСК в текущую НЗСК
long double THelicTraj::recalcPsi_NZSK_to_CurNZSK(long double valPsi_NZSK)
{
    // угол повороьа текущей НЗСК относительно НЗСК
  double valAlfRot = atan2(marrTransf_CurNZSK_to_NZSK[6], marrTransf_CurNZSK_to_NZSK[0]);

  double temp = valAlfRot + static_cast< double>(valPsi_NZSK);
  return static_cast<long double> (TPartHelicTraj::fnc_Minus_PI_Plus_PI( temp ));

}
//---------------------------------------------------------------------------


// пересчет угла рыскания из текущей НЗСК в  НЗСК
long double THelicTraj::recalcPsi_CurNZSK_to_NZSK(long double valPsi_CurNZSK)
{
    // угол повороьа текущей НЗСК относительно НЗСК
  double valAlfRot = atan2(marrTransf_CurNZSK_to_NZSK[6], marrTransf_CurNZSK_to_NZSK[0]); 
  double temp = -valAlfRot + static_cast< double>(valPsi_CurNZSK);
  return static_cast<long double> (TPartHelicTraj::fnc_Minus_PI_Plus_PI( temp ));

}


 //---------------------------------------------------------------------------
// формирование матрицы перехода из НЗСК в ТНЗСК для участка траектории с ном еровм INumPart
 // нумерация траектрий начинается с 1 !!!
void THelicTraj::create_transfMtrx_NZSK_to_CurNZSK(const int INumPart)
{

   memset(marrTransf_CurNZSK_to_NZSK, 0, 9 * sizeof(long double));
   marrTransf_CurNZSK_to_NZSK[0] = cos(mFlyTask.mparrPartTrajInputData[INumPart].mCourse);
   marrTransf_CurNZSK_to_NZSK[2] = -sin(mFlyTask.mparrPartTrajInputData[INumPart].mCourse);
   marrTransf_CurNZSK_to_NZSK[4] = 1.;
   marrTransf_CurNZSK_to_NZSK[6] = -marrTransf_CurNZSK_to_NZSK[2];
   marrTransf_CurNZSK_to_NZSK[8] = marrTransf_CurNZSK_to_NZSK[0];

}

//---------------------------------------------------------------------------
// формирование матрицы перехода из ТНЗСК в НЗСК для участка траектории с ном еровм INumPart
// нумерация траектрий начинается с 1 !!!
void THelicTraj::create_transfMtrx_CurNZSK_to_NZSK(const int INumPart)
{
  memset(marrTransf_NZSK_to_CurNZSK, 0, 9 * sizeof(long double));
  marrTransf_NZSK_to_CurNZSK[0] = cos(mFlyTask.mparrPartTrajInputData[INumPart].mCourse);
  marrTransf_NZSK_to_CurNZSK[2] = sin(mFlyTask.mparrPartTrajInputData[INumPart].mCourse);
  marrTransf_NZSK_to_CurNZSK[4] = 1.;
  marrTransf_NZSK_to_CurNZSK[6] = -marrTransf_NZSK_to_CurNZSK[2];
  marrTransf_NZSK_to_CurNZSK[8] = marrTransf_NZSK_to_CurNZSK[0];
}

// формирование ВЕКТОРА положения начала координат  ТНЗСК в НЗСК для участка траектории с номером INumPart
// нумерация траектрий начинается с 1 !!!
void THelicTraj::create_vectO(const int INumPart)
{
  marrO_CurNZSK[0] = mFlyTask.mparrPartTrajInputData[INumPart].mTaskPointsZ.X;//    mFlyTask.mparrTaskPointsZ[INumPart].X;
  marrO_CurNZSK[2] = mFlyTask.mparrPartTrajInputData[INumPart].mTaskPointsZ.Z;//mFlyTask.mparrTaskPointsZ[INumPart].Z;
  marrO_CurNZSK[1] = 0.;
}

// пересчет вектора из НЗСК в ТНЗСК
void THelicTraj::transf_xyzNZSK_to_xyzCurNZSK(long double *arrNZSKInp, long double *arrCurNZSKOut)
{
    long double arrT0[3] = {0.};
    MtrxMinusMatrx(arrNZSKInp, marrO_CurNZSK, 1, 3, arrT0);
    MtrxMultMatrx(marrTransf_NZSK_to_CurNZSK,3, 3, arrT0,1, arrCurNZSKOut) ;
}

//------------------------------------------------------------------

// пересчет вектора из TНЗСК в НЗСК
void THelicTraj::transf_xyzCurNZSK_to_xyzNZSK(long double *arrCurNZSKInp, long double *arrNZSKOut)
{
    long double arrT0[3] = {0.};
    MtrxMultMatrx(marrTransf_CurNZSK_to_NZSK,3, 3, arrCurNZSKInp,1, arrT0) ;
    MtrxSumMatrx(arrT0, marrO_CurNZSK, 1, 3, arrNZSKOut);
}


//-------------------------------------------------------------------



// загрузка таблицы магнитного склонения из файла русурса в qвектор mqvectTabMagDev
void THelicTraj::fillMagDevTab()
{
    mQuantRowsTabMagDev = 360;
    mQuantColsTabMagDev = 181;
    mqvectTabMagDev.resize(mQuantRowsTabMagDev *mQuantColsTabMagDev);

    QFile MagDev;

    std::vector<double> D(3);

    MagDev.setFileName(":/deviations.csv");
    MagDev.open(QIODevice::ReadOnly | QIODevice::Text);
    int icur =0;
    QString line0 = MagDev.readLine();
    while( !MagDev.atEnd())
    {
            QString line = MagDev.readLine();
            QStringList w =line.split(";");

            if (w.at(0).contains("S"))
            {
                D[0] = (-1.0) * w.at(0).right(2).toDouble();
            }
            else
            {
                D[0] = w.at(0).right(2).toDouble();
            }

            if (w.at(1).contains("W"))
            {
                D[1] = (-1.0) * w.at(1).right(3).toDouble();
            }
            else
            {
                D[1] = w.at(1).right(3).toDouble();
            }
            QString line2 = w.at(2);

             QString y = ",";
             if(line2.indexOf(y) >= 0)
                 {
                 line2.replace(line2.indexOf(y), 1, '.');

                 }


            D[2] = line2.toDouble();

            mqvectTabMagDev[icur] = D[2];// * M_PI/ 180.;

      icur++;
    }
    MagDev.close();
}
// вычислениме величины магнитного склонения
// INPUT:
// VAlLongtitude - долгота точки. Измеряется от 0 до 359 град в направлении на восток. всего 360 строк
//                с шагом в 1 град.
// VAlLatitude - широта точки. Измеряется от -90 до 90 град, от южного полюся к северному. всего 181 столбец
double THelicTraj::calcMagDevValue(const double VAlLongtitude, const double VAlLatitude0)
{
    const double VAlLatitude = VAlLatitude0 + 90.;
    int irow = VAlLongtitude;
    int irowNext = (irow + 1) % mQuantRowsTabMagDev;
    int icol = (int)(VAlLatitude );

    // отсечь случай icol= +-90 !!!!!!
    if (icol == mQuantColsTabMagDev)
    {
     return mqvectTabMagDev[irow * mQuantColsTabMagDev + icol];
    }
    // координаты точки внутри пиксела:
    double xt = VAlLatitude - ((double)icol);
    double yt = VAlLongtitude - ((double)irow);
    ///

    // формирование шаблона из 4-х окружающих точек
    double arrSh[4] = {0.} ;
    arrSh[0] = mqvectTabMagDev[irow * mQuantColsTabMagDev + icol];
    arrSh[1] = mqvectTabMagDev[irow * mQuantColsTabMagDev + icol +1];
    arrSh[2] = mqvectTabMagDev[irowNext * mQuantColsTabMagDev + icol];
    arrSh[3] = mqvectTabMagDev[irowNext * mQuantColsTabMagDev + icol +1];
    ///
    return BilinearValue (arrSh,1., xt, yt);
}


// билинейная интерполяция
double THelicTraj::BilinearValue (double* arrSh,const double cellsize
            ,const double xt,const double yt)
{
  double ax = xt/ cellsize ;
  double ay = yt/ cellsize ;
  double z4 = (1 - ay) *  arrSh[0] +  ay * arrSh[2] ;
  double z5 = (1 - ay) *  arrSh[1] +  ay * arrSh[3] ;
  return ( (1- ax) * z4 + ax * z5);
}
//----------------------------------------------------------------------------

long double THelicTraj::SIGNUM(long double x)
{
   return (x >= 0.)?1.:-1.;
}
 //-------------------------------------------------------------------------------------------
 //-----------------------------------------------------------------------------------
  // перемещение вертолета в ссответсвии с полетным заданием
  // на время полета VAlFlyTime
  // INPUT:
  // VAlFlyTime - время полета
  // VAlStIntegr - шаг интегрирования
 // LEnarrC - длина массива передаточных чисел = 44
  //OUTPUT:
  // arrBuff - буфер (массив) для накопления информации
  // Если arrBuff==NULL то информация не накапливается
  // *pquantSteps - к-во шагов интегрирования

  // вектор выходной информации arrOutputVect[25]
 // составляющие воздушной скорости в связанной системе координат; 3
 // составляющие путевой скорости в связанной системе координат; 3
 // приборная скорость;                                          1
 // бароинерциальная вертикальная скорость;                      1
 // геометрическая высота;                                       1
 // абсолютная барометрическая высота;                           1
 // относительная барометрическая высота;                         1
 // геодезическая высота;                                        1
 // бароинерциальная высота;                                     1
 // угол сноса;                                                  1
 // угол атаки;                                                  1
 // угол скольжения;                                              1
 // скорость и направление ветра;                                2
 // крен;                                                        1
 // тангаж;                                                      1
 // курс истинный;                                                1
  //магнитное склонение;                                          1
 // температура наружного воздуха;                               1
 // текущие координаты местоположения объекта в географической системе координат.
 // INPUT:
 //VAlFlyTime - ориентировочное значение полетного времени заведомо превосходящее реальное
 //VAl_OutPut_dt - шаг с которым выдается информация во внешнюю программу
 //VAlStIntegr - шаг интегрирования
 //OUTPUT:
 //arrBuff - массив с выходной информацией построения моих графиков
 //*pquantSteps - к-во реально сделанных шагов с интервалом VAl_OutPut_dt
 // arrOutSideInfo - массив с информацией для внешней мсодели с зарезервированной памятью с (1 + QUANT_COLS_OUT_VECT) столбцов
// QUantColsOutArr - к-во сооллбцов массива arrOutSideInfo (или число переменных)
 void THelicTraj::imitateMoving_and_collectData( const long double VAlFlyTime ,const long double  VAl_OutPut_dt
                                                 , double *arrBuff, const int QUantColsBuff
                                                 , int *pquantSteps, double *arrOutSideInfo)

  {
      // к-во столбцов вектора внешней модели, используемое для плстроения графиков
     // в первом столбце ( с нулевым номером) хранится время, далее хранится QUANT_COLS_OUT_VECT
     // элементов вектора
     const int  QUantColsOutArr = 1 + QUANT_COLS_OUT_VECT;
      // признак полета
      bool bFly = false;
      ///
   int quantSt = (int)(VAlFlyTime/ VAl_OutPut_dt) + 1;
      *pquantSteps = 0;
   //  заполнение буфера с информацией
      fillOutPutBuffRow(arrBuff);

      // заполнение буфера с для графического отображения выходного вектора
      arrOutSideInfo[ 0] = mTimeCur;
      fillOutPutVect( &(arrOutSideInfo[1]));
      (*pquantSteps)++;
   ///


    for (int i =1; i < quantSt; ++i)
    {
        //  перемещение вертолета на время VAl_OutPut_dt вперед
        long double valTEnd = mTimeCur + VAl_OutPut_dt;

        bool brez = goAhead_in_accordance_with_time(valTEnd);

        arrOutSideInfo[ i * QUantColsOutArr] = mTimeCur;


        fillOutPutBuffRow(&(arrBuff[ i * QUantColsBuff]));
        (*pquantSteps)++;
        if ( ((*pquantSteps) > 2) &&((marrPhaseVect[1] < 0.0005) || (!brez)))
        {
            break;
        }
    }

  }
 //----------------------------------------------------------
 // вектор выходной информации arrOutputVect[25]
// составляющие воздушной скорости в связанной системе координат; 3
// составляющие путевой скорости в связанной системе координат; 3
// приборная скорость;                                          1
// бароинерциальная вертикальная скорость;                      1
// геометрическая высота;                                       1
// абсолютная барометрическая высота;                           1
// относительная барометрическая высота;                         1
// геодезическая высота;                                        1
// бароинерциальная высота;                                     1
// угол сноса;                                                  1
// угол атаки;                                                  1
// угол скольжения;                                              1
// скорость и направление ветра;                                2
// крен;                                                        1
// тангаж;                                                      1
// курс истинный;                                                1
 //магнитное склонение;                                          1
// температура наружного воздуха;                               1
// текущие координаты местоположения объекта в географической системе координат. 2

 bool THelicTraj::goAhead_in_accordance_with_time(const long double VAlTend)
 {

     if (mpPartTraj->IsEndOfPart())
         {

         if (mNumPart == 11)
         {
          int jj =0;
         }
             if (mNumPart == (mFlyTask.mquantTaskPoint -1  ))
             {
                 return false;
             }
             mNumPart++;
             initPartTraj(mNumPart);
         }

     mpPartTraj->replaceYourSelf( VAlTend,  mStepIntegr);
     mTimeCur = VAlTend;
     transform_VS_from_CurNZSK_to_NZSK();

     return true;
 }

 //---------------------------------------------------------------
 // формирование вектора выходной информации для внешней программы
 // вектор выходной информации arrOutputVect[30]
 // составляющие воздушной скорости в связанной системе координат; 3
 // составляющие путевой скорости в ГСК (восточная, северная, модуль); 3
 // приборная скорость;                                          1
 // бароинерциальная вертикальная скорость;                      1
 // геометрическая высота;                                       1
 // абсолютная барометрическая высота;                           1
 // относительная барометрическая высота;                         1
 // геодезическая высота;                                        1
 // бароинерциальная высота;                                     1
 // угол сноса;                                                  1
 // угол атаки;                                                  1
 // угол скольжения;                                              1
 // скорость и направление ветра;                                2
 // крен;                                                        1
 // тангаж;                                                      1
 // курс истинный;                                                1
 //магнитное склонение;                                          1
 // температура наружного воздуха;                               1
 // текущие координаты местоположения объекта в географической системе координат. 2
 // модуль воздушной скорости                         1
 // вектор путевой скорости в СвСК                    3
 // боковая пергрузка                                   1
 void THelicTraj::fillOutPutVect( double *arrOutputVect0)
 {
     long double  arrOutputVect[QUANT_COLS_OUT_VECT ] ={0.};
     // 1.- вектор воздушной скорости в СвСК. arrOutputVect[0]- arrOutputVect[2]
     mpPartTraj->calcUaSvSK_Case_NZSK(arrOutputVect);
     ///

     // 2. составляющие путевой скорости в ГСК. arrOutputVect[3]- arrOutputVect[5]
     long double arrV_Ground[3];
     transf_xyzNZSK_to_xyzGSK(&(marrPhaseVect[3]), arrV_Ground);
     arrOutputVect[3] = arrV_Ground[0]; // восточная
     arrOutputVect[4] = arrV_Ground[1]; // северная
     arrOutputVect[5] = Norm3(arrV_Ground); // модуль

     // составляющие путевой скорости в связанной системе координат. arrOutputVect[3]- arrOutputVect[5]
     long double arrV_SvSK[3] = {0.};
     long double arrMtrxTransf_NZSK_SvSK[9] = {0.};
     mpPartTraj->calcMtrxTransf_CurNZSK_To_SvSK (arrMtrxTransf_NZSK_SvSK);
     MtrxMultMatrx(arrMtrxTransf_NZSK_SvSK,3, 3, &(mpPartTraj->marrPhaseVect[3]), 1,arrV_SvSK);



     ///

     // 3. Приборная скорость arrOutputVect[6]
     // https://en.wikipedia.org/wiki/Impact_pressure
    //https://ru.wikipedia.org/wiki/Приборная_скорость
      // вычисление числа Маха
     long double valTKel0 = 273.15 + mEnvironment.mAirT0;
     long double valTay = valTKel0  - 0.0065 * marrPhaseVect[1];
     long double valMach = mpPartTraj->calcMach(valTay,Norm3(arrOutputVect) );
     ///

     // вычисление относительной плотности воздуха
     long double valRelRo = mEnvironment.calcAirRelativeDensity(marrPhaseVect[1]);
     ///

     // вычисление CAS.
     long double temp0 = powl(1. + 0.2 * valMach * valMach, 7./2.);
     long double temp1 = powl(valRelRo *(temp0 - 1.) + 1., 2./7.);
     arrOutputVect[6] = ATM_AN0 * sqrtl( 5 *(temp1 -1.));
  ///

     // 4. бароинтерциальная вертикальная скорость. arrOutputVect[7]
     arrOutputVect[7] = marrPhaseVect[4];
     ///

     // 5. геометрическая высота;  arrOutputVect[8]
     arrOutputVect[8] = marrPhaseVect[1];

     ///
     // 6. абсолютная барометрическая высота;  arrOutputVect[9]
     arrOutputVect[9] = marrPhaseVect[1];
     ///


     // 7. относительная барометрическая высота;
     arrOutputVect[10] = marrPhaseVect[1];
     ///

     // 8. геодезическая высота;
     arrOutputVect[11] = marrPhaseVect[1];
     ///

     // 9.бароинерциальная высота;      arrOutputVect[12]
     arrOutputVect[12] = marrPhaseVect[1];
     ///

     // 10. Угол сноса (УС) — угол между векторами: воздушной скорости и путевой скорости.
     // Правый снос (+), когда самолет сносит относительно линии курса вправо, левый (-),
     // когда самолет сносит относительно линии курса влево.      arrOutputVect[13]
     // http://adminland.ru/crimea/books/sh_hb/part04.htm
     // модуль путевой скорости
     long double valModV = Norm3(arrV_SvSK);
     long double valModVa = Norm3(arrOutputVect);
      long double  arrWindV[3] = {0.};
      mEnvironment.createVectWindV(arrWindV);
      long double valModWindV = Norm3(arrWindV);


     if ((valModV <  0.0001) || (valModVa <  0.0001)|| (valModWindV < 0.01))
     {
         arrOutputVect[13] = 0.;
     }
     else
     {
         long double t = ScalProduct(arrOutputVect , arrV_SvSK, 3) /valModV/ valModVa;
         if (fabsl(t) > 1.)
         {
             t = SIGNUM(t);
         }
         arrOutputVect[13] = acosl(t) * 180./ M_PI;


     // определение знакаа угла сноса
     long double arrn[3] = {0.}, arrZenit[3] = {0.,1.,0.};
     OuterProduct(arrV_SvSK , arrZenit, arrn) ;

     if (ScalProduct(arrn , arrWindV, 3) > 0.)
     {
        arrOutputVect[13] = fabsl(arrOutputVect[13]) ;

     }
     else
     {
         arrOutputVect[13] = -fabsl(arrOutputVect[13]) ;
     }
     }

     ///

     // 11. угол атаки. arrOutputVect[14]
     // по поределению это  угол между продольной осью летательного аппарата (ЛА) OX и проекцией скорости
     // ЛА на плоскость OXY СвСК; Угол атаки следует считать положительным, если в СвСК Vy <0. То есть, это -arctg(Vy/Vx)
     double val_v = sqrtl(arrV_SvSK[0] * arrV_SvSK[0]  + arrV_SvSK[1] * arrV_SvSK[1]);
     if (val_v  < 0.000001)
     {
         arrOutputVect[14] = 0.;
     }
     else
     {
      arrOutputVect[14] = -asin(arrV_SvSK[1] / val_v) * 180. / M_PI;
     }

     ///


     // 12. угол скольжение. arrOutputVect[15]
     // по определению это  угол между вектором скорости  летательного аппарата (ЛА)  и плоскостью
     // OXY СвСК; Угол скольжение следует считать положительным, если в СвСК Vz >0.
     if ( valModV < 0.0001)
     {
       arrOutputVect[15] = 0.;
     }
     else
     {
         arrOutputVect[15] = SIGNUM(arrV_SvSK[2])
                 * acosl(sqrtl(arrV_SvSK[0] * arrV_SvSK[0] + arrV_SvSK[1] * arrV_SvSK[1])/ valModV)
                 * 180./ M_PI;

     }
     ///

     // 13. скорость и направление ветра; arrOutputVect[16] - arrOutputVect[17]
     arrOutputVect[16] = mEnvironment.mWind_V;
     arrOutputVect[17] = mEnvironment.mWind_Alf * 180./ M_PI;
     ///

     // 14. крен;  arrOutputVect[18]
     arrOutputVect[18] = mpPartTraj->marrPhaseVect[9] * 180. / M_PI;
     ///


     // 15. тангаж;    arrOutputVect[19]
     arrOutputVect[19] = mpPartTraj->marrPhaseVect[10] * 180. / M_PI;
     ///

     //16.  курс истинный; arrOutputVect[20]
     long double temp00 =  atan2l(arrV_Ground[0], arrV_Ground[1]) ;
     arrOutputVect[20] = temp00 * 180. / M_PI;
     ///

    // 17. магнитное склонение; arrOutputVect[21]
            // вычисление текущих значений координат в ГСК  !!!!!!
            long double arrGSKTemp[3] = {0.}, arrGSKCur[3] = {0.};
            transf_xyzNZSK_to_xyzGSK(marrPhaseVect, arrGSKTemp);
            MtrxSumMatrx(arrGSKTemp, marr_GSKxyz_Pnt0,1, 3, arrGSKCur) ;
            ///

            // вычисление текущих долготы и широты
             double valLongtitude = 0., valLatitude = 0.;
             mhGauss2Geo_( arrGSKCur[0],  arrGSKCur[1], &valLatitude, &valLongtitude );
             ///


       arrOutputVect[21] = calcMagDevValue( valLongtitude * 180. / M_PI,  valLatitude * 180. / M_PI);
     ///

     //18. температура наружного воздуха; arrOutputVect[22]     
     arrOutputVect[22] =mEnvironment.mAirT0  - 0.0065 * marrPhaseVect[1];
     ///

     // 19. текущие координаты местоположения объекта в географической системе координат.
     arrOutputVect[23]  = valLatitude * 180./ M_PI;
     arrOutputVect[24]  = valLongtitude * 180./ M_PI;

     // 20. модуль воздушной скорости
       arrOutputVect[25] = Norm3(arrOutputVect);

     // 21. вектор путевой скорости в СвСК arrOutputVect[26] - arrOutputVect[28]
       memcpy(&(arrOutputVect[26]), arrV_SvSK, 3 * sizeof(long double));

     // 22. боковая пергрузка
       arrOutputVect[29] = mpPartTraj->calcPeregrNorm();


     for (int i = 0; i < QUANT_COLS_OUT_VECT; ++i)
     {
         arrOutputVect0[i] = arrOutputVect[i];
     }

     for (int i =0; i < 8; ++i)
     {
     arrOutputVect0[i] = arrOutputVect0[i] * 3600. / 1000.;
     }
     arrOutputVect0[16] = arrOutputVect0[16] * 3600. / 1000.;

     for (int i = 25; i < 29; ++i)
     {
      arrOutputVect0[i] = arrOutputVect0[i] * 3600. / 1000.;
     }

 }



 //-----------------------------------------
 void THelicTraj::fillOutPutBuffRow(double *p)
 {
     p[0]= mTimeCur;


     // положение и скорость в ГСК
     long double arrVS_GSK[6] = {0.};
     transf_xyzNZSK_to_xyzGSK(marrPhaseVect, arrVS_GSK);
     transf_xyzNZSK_to_xyzGSK(&(marrPhaseVect[3]), &(arrVS_GSK[3]));

     for (int j =0; j < 6  ; ++j)
     {
      p[1 + j]  = arrVS_GSK[j];
     }
     //

     //  углы эйлера
     for (int j =0; j < 3 ; ++j)
     {
      p[7 + j]  =  marrPhaseVect[9 + j];
      p[7 + j] = TPartHelicTraj::fnc_Minus_PI_Plus_PI( p[7 + j] );
     }


 }
//------------------------------------------------------------------
//------------------------------------------------------------------
 // ВЫЧИСЛЕНИЕ вектора состояния вертолета в НЗСК
 //Задан вектор состояния вертолета в члене класса mpPartTraj
 // он состоит из 13 элементов
 // требуетмя вычислить вектор состояния вертолета в НЗСК
 // который состоит из 6 элементов - 3 положения и 3 скорости
 void THelicTraj::transform_VS_from_CurNZSK_to_NZSK()
 {
     transf_xyzCurNZSK_to_xyzNZSK(mpPartTraj->marrPhaseVect, marrPhaseVect);
     MtrxMultMatrx(marrTransf_CurNZSK_to_NZSK,3, 3, &(mpPartTraj->marrPhaseVect[3]) ,1, &(marrPhaseVect[3])) ;
     memcpy(&(marrPhaseVect[6]), &(mpPartTraj->marrPhaseVect[6]), 5 * sizeof(long double));

     marrPhaseVect[11] = recalcPsi_CurNZSK_to_NZSK(mpPartTraj->marrPhaseVect[11]);

     marrPhaseVect[12] = mpPartTraj->marrPhaseVect[12];


 }

 //------------------------------------------------------------------


 // перевод вектора из НЗСК в ГСК
 // сиситемы координат
 // ГСК - ось X на восток, Y на север, Z в зенит
 // НЗСК - ось X на север, Y  в зенит, Z на восток
 void THelicTraj::transf_xyzNZSK_to_xyzGSK(long double *arrNZSKInp, long double *arrGSKOut)
 {
     arrGSKOut[0] =  arrNZSKInp[2];
     arrGSKOut[1] =  arrNZSKInp[0];
     arrGSKOut[2] =  arrNZSKInp[1];
 }
 //------------------------------------------------------------------

 // перевод вектора из НЗСК в ГСК
 // сиситемы координат
 // ГСК - ось X на восток, Y на север, Z в зенит
 // НЗСК - ось X на север, Y  в зенит, Z на восток
 void THelicTraj::transf_xyzGSK_to_xyzNZSK(long double *arrGSKInp, long double *arrNZSKOut)
 {
     arrNZSKOut[0] =  arrGSKInp[1];
     arrNZSKOut[1] =  arrGSKInp[2];
     arrNZSKOut[2] =  arrGSKInp[0];
 }

 //------------------------------------------------------------------

 // перевод вектора из НЗСК в ГСК
 // сиситемы координат
 // ГСК - ось X на восток, Y на север, Z в зенит
 // НЗСК - ось X на север, Y  в зенит, Z на восток
 void THelicTraj::transf_xyzNZSK_to_xyzGSK( double *arrNZSKInp,  double *arrGSKOut)
 {
     arrGSKOut[0] =  arrNZSKInp[2];
     arrGSKOut[1] =  arrNZSKInp[0];
     arrGSKOut[2] =  arrNZSKInp[1];
 }
 //------------------------------------------------------------------

 // перевод вектора из НЗСК в ГСК
 // сиситемы координат
 // ГСК - ось X на восток, Y на север, Z в зенит
 // НЗСК - ось X на север, Y  в зенит, Z на восток
 void THelicTraj::transf_xyzGSK_to_xyzNZSK( double *arrGSKInp,  double *arrNZSKOut)
 {
     arrNZSKOut[0] =  arrGSKInp[1];
     arrNZSKOut[1] =  arrGSKInp[2];
     arrNZSKOut[2] =  arrGSKInp[0];
 }



