#include "FlyTask.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>



// параметры инициализации полетного задания вариант 0
// высота 300 м
// движение по прямоугольнику
// все задано в ГСК

extern const int QUantActualTaskPoints0 = 14;

extern  long double   ARRCourse0[] = {0.,-M_PI / 2., -M_PI / 2., 0., 0.
                                      ,M_PI /2. , M_PI/2.,  M_PI,  M_PI, -M_PI /2.
                                      , -M_PI /2., -M_PI /2., -M_PI /2., -M_PI /2.};//, -M_PI /2.,-M_PI /2.,-M_PI/ 2,-M_PI/2.,-M_PI/2 };

extern  long double   ARRV0 [] =     {0., 40.,   40.,   40.,   50.
                                    ,50., 50.,   50.,   50.,   50.,
                                     40., 40.,    10.,    0.};

 extern  long double   ARRCoord0[] = {
                                0.,    0.,   0.  //1. стоит на месте
                           ,-2225.,    0., 300.// 2.взлет с (0,0,0 ) на  скорость 40 на запад и движение до точки (-2225.,    0., 300.)
                           ,-2225.,    0., 300.//3. поворот на 90 град на север радиус 675 скорость 40
                           ,-2900., 2580., 300. // 4. движение на север до точки (-2900., 2580., 300.)
                          ,-2900., 2580., 300. //5. поворот на 90 град на восток радиус 820 скорость в конце 50

                          , 4180., 3400., 300. //6. движение на восток до точки (4180., 3400., 300.) скорость 50
                          , 4180., 3400., 300.//7.поворот на 90 скорость 50 на юг
                          , 5000., 820., 300. // 8.движение на юг до точки (5000., 820., 300.) скорость 50
                          , 5000., 820., 300. //9.поворот на 90 скорость 50 на запад
                          , 3000., 0., 300. // 10.движение на запад до точки (3000., 0., 300.) скорость 50

                          , 2000., 0., 300.// 11. торможение  на запад до точки (2000., 0., 300.)до скорости 40
                          ,  500.,    0., 100.  // 12.спуск (пикирование) на высоту 100 к точке (500.,    0., 100.)скорость 40
                          ,  50.,    0., 20.   // 13. спуск на высоту 20м  с торможением до 10 м на запад до точки (50., 0., 300.)
                          ,  0.,    0., 0.// 14. посадка по самолетному
                             };



extern enumTypeOfManeuvre ARRTypeMove0[] = {enLine, enLine, enTurn, enLine, enTurn
                                           ,enLine, enTurn, enLine, enTurn , enLine
                                           ,enLine, enLine, enLine, enLine};
extern long double   ARRYaw0[] = {0.,0.,0., 0.,0.
                                 ,0.,0.,0., 0.,0.
                                 ,0.,0., 0.,0.};
extern long double   ARRRadius0[] = {
                                             0.,100000000.,        675., 100000000.,       820.
                                    ,100000000.,      820.,  100000000.,       820., 100000000.
                                    ,100000000.,100000000.,  100000000.,  100000000.
                                   };
extern int IARRNumSetOfGears0[] = {-1, 2, 7,3, 8, 3, 8, 3, 8, 3, 2, 5, 1, 0};

extern long double   ARRCourseEnd0[] ={M_PI/ 2.,M_PI/ 2.,M_PI/ 2. + 2.*M_PI, M_PI/ 2., M_PI/ 2.
                                      ,M_PI/ 2.,M_PI/ 2.,           M_PI/ 2.,M_PI/ 2., M_PI/ 2.
                                      ,M_PI/ 2.,M_PI/ 2.,           M_PI/ 2.,M_PI/ 2.};

extern long double   ARRTHovering0[] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

// угловая скорость вращения
extern long double   ARR_dPsi_po_dt0[] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};



///


extern const int QUantActualTaskPoints1 = 26;

extern  long double   ARRCourse1[] = {0.,-M_PI / 2.,-M_PI / 2.,-M_PI / 2., -M_PI / 2.
                                      ,-M_PI / 4.,  -M_PI / 2., -M_PI / 2., -M_PI /2. +0.01, -M_PI / 2.
                                      , 0., 0., M_PI-0.001, 0. , M_PI - 0.001
                                      , 0. , M_PI / 2., M_PI / 2.,   M_PI/2. , M_PI/4.
                                      , M_PI/ 2., M_PI/ 2., M_PI, M_PI, M_PI
                                      , M_PI};

extern  long double   ARRV1 [] =     {0., 40., 40., 10., 0.
                                      , 0., 0., 0., 0.,40.
                                      , 40., 40., 40., 40., 40.
                                      , 40., 50.,    0., 0. ,0.
                                     , 0., 40., 40., 10., 0.
                                      ,0.};//,   40.,   40.,   50., 50., 50., 50., 50., 50.,   40.,   40., 10., 0.};//,     50. , 40.,   40.,  10.,   0.};

 extern  long double   ARRCoord1[] = {
                                0.,    0.,   0.  //0. стоит на месте
                           ,-2225.,    0., 100.// 1.взлет-кабрирование с (0,0,0 ) на  скорость 40 на запад и движение до точки (-2225.,    0., 100.)
                            ,-5000.,    50., 300.//2. кабрирование со смещением на север до скорости 40 на запад до точки (-5000.,    50., 300.)
                           ,-6000.,    50., 300. //3. торможение до скорости 10
                           ,-6047.6,    50., 300. //4. висение в течение 10 с

                           ,-6047.6,    50., 100. //5. вертикальное снижение до высоты 100 м c поворотом на 45 град и висение
                           , -6047.6,    50., 300.  //6. подъем на высоту 300 м с вращением 2 рад/с
                           , -6047.6,    50., 300.  //7. гашение угловой скорости до 0,1 рад/с на высоте 300 м
                           , -6047.6,    50., 300. // 8. остановка с углом рыскания 180 град в нзск (нос смотрит на запад) и висение
                          , -8000.,    0., 300. //9. разгон  до скорость 40 курс -90град

                           , -8675.,  675., 300. // 10.enTurn вираж с поворотом на север и подъемом до 300 и скорости до 40
                            ,-8675.,  3000., 300. // 11. движение на север до точки (-8000.,  3000., 300.)
                             ,-7325.,  3000., 300. //12. enTurn звено змейки 1
                            ,-5975.,  3000., 300. //13.enTurn звено змейки 2
                           ,-4625.,  3000., 300. //14.enTurn звено змейки 3

                           ,-3275.,  3000., 300. //15.enTurn звено змейки 4
                          , -2455., 3820., 300.   // 16. выход из змейки виражом с поворотом на восток на 90 град и скоростью 50    
                          ,-2200.2,3820., 300. // 17. торможение до 0
                          ,-2200.2,3820., 300. // 18. поворот на месте на -45 град
                          ,-2200.2,3820., 300. // 19. поворот на месте на 45 град

                           ,-2200.2,3820., 300. // 20. поворот на месте на -90 град
                          ,  -675,3820., 300.  // 22. разгон на восток  до 40
                           ,  0.,3145., 300.  // 23. поворот на юг
                           ,  0.,48., 100.  // 24. кабрирование с высоты 300 и скорости 40 на высоту 100  и скорость 10
                           ,  0.,0., 100.  // торможение на высоте 100 до 0

                           ,  0.,0., 0.
  };




 extern enumTypeOfManeuvre ARRTypeMove1[] = {enLine, enLine, enLine, enLine,  enLine
                                             , enSpotHovering, enRotating,  enRotating,enSpotHovering, enLine
                                             ,enTurn, enLine, enTurn, enTurn, enTurn
                                              , enTurn, enTurn, enLine,   enSpotHovering, enSpotHovering
                                            , enSpotHovering, enLine,  enTurn, enLine, enSpotHovering
                                             , enSpotHovering};

 extern long double   ARRYaw1[] = {0.,0., 0.,0., 0.
                                   , 0., 0.,0., 0., 0.
                                   , 0.,0.,0.,0., 0.
                                   ,0., 0., 0., 0., 0.
                                  ,0., 0., 0., 0.
                                  ,0.};

 extern long double   ARRRadius1[] = {0.,100000000.,100000000.,100000000.,100000000.
                                      ,100000000.,100000000.,100000000.,100000000.,100000000.
                                      , 675.,100000000., 675., -675., 675.
                                      , -675., 820.,100000000.,100000000.,100000000.
                                      ,100000000.,100000000., 675.,100000000.,100000000.
                                      ,100000000.};

extern int IARRNumSetOfGears1[] = {-1, 1, 3, 1, 9
                                   ,14, 11, 11, 12, 3/*8*/
                                   ,7, 3, 8, 8, 8
                                   , 8, 8, 9,   12, 12
                                   , 12, 2, 7, 1, 9
                                   , 9};//, 7};//, 7,3, 8, 3, 8, 3, 8, 3, 2, 5, 1, 0};

extern long double   ARRCourseEnd1[] ={0.,0.,0.,0.,0.
                                       ,0.,0., 0., 0., M_PI/ 2.
                                       , M_PI/ 2.,0., M_PI-0.001, -M_PI, M_PI - 0.001
                                       , -M_PI , M_PI/ 2., 0., 0., 0.
                                       , 0.,0.,  M_PI/ 2., 0., 0.
                                       , 0.};

extern long double   ARRTHovering1[] ={0.,0.,0.,0.,30.
                                       ,30.,30. ,19., 30., 0.
                                       , 0. ,0.,0., 0., 0.
                                       ,0., 0., 0., 40.,40.
                                      , 40.,0., 0., 0., 15.
                                       , 120.};

extern long double   ARR_dPsi_po_dt1[] ={0.,0.,0.,0.,0.
                                         ,0.,1. ,0.1, 0., 0.
                                         , 0. ,0.,0., 0., 0.
                                         ,0. , 0., 0., 30., 30.
                                        , 30.,0., 0., 0., 10.
                                         , 10.};


//-------------------------------------------------------------------------------
 //-------------------------------------------------------------------------------
 //-------- ОПИСАНИЕ КЛАССА TPointXYZ -----------------------------------------------------------------------
 //-------------------------------------------------------------------------------
 //-------------------------------------------------------------------------------


TPointXYZ::TPointXYZ()
{
 X =0.;
 Y = 0.;
 Z = 0.;
 M = 0.;

}
// Конструктор копирования
TPointXYZ::TPointXYZ (const TPointXYZ &R)
{
    X = R.X;
    Y = R.Y;
    Z = R.Z;
    M = R.M;
  }
// оператор присваивания
 TPointXYZ TPointXYZ::operator=(TPointXYZ  R)
{
     X = R.X;
     Y = R.Y;
     Z = R.Z;
     M = R.M;
    return *this ;
}

 // парам констр 1
 TPointXYZ::TPointXYZ(const long double   X_,const long double   Y_,const long double   Z_, const long double   M_)
 {
  X = X_;
  Y=Y_;
  Z=Z_;
  M=M_;
 }
 // парам констр 2
 TPointXYZ::TPointXYZ(const long double   X_,const long double   Y_,const long double   Z_)
 {
  X = X_;
  Y=Y_;
  Z=Z_;
  M=0.;
 }
 //-------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------
  //-------- ОПИСАНИЕ КЛАССА TPartTrajInputData-----------------------------------------------------------------------
  //-------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------


TPartTrajInputData::TPartTrajInputData()
{
    mTaskPointsZ  = TPointXYZ();
    mV = 0.;// одномерный массив горизонтальных скоростей
    mCourse = 0.;// одномерный массив курсовых углов
    menumType = enLine;
    mYaw = 0.;
    mRadius =0.;
    mNumSetOfGears = 0;
    mCourseEnd = 0.;
    mTHovering = 0.;
    m_dPsi_po_dt = 0.;
    //mPsi = 0.;

}
// Конструктор копирования
TPartTrajInputData::TPartTrajInputData (const TPartTrajInputData &R)
{
    mTaskPointsZ  = R.mTaskPointsZ;
    mV = R.mV; // одномерный массив горизонтальных скоростей
    mCourse = R.mCourse;// одномерный массив курсовых углов
    menumType = R.menumType;
    mYaw = R.mYaw;
    mRadius = R.mRadius;
    mNumSetOfGears = R.mNumSetOfGears;
    mCourseEnd = R.mCourseEnd;
    mTHovering = R.mTHovering;
    m_dPsi_po_dt = R.m_dPsi_po_dt;
    //mPsi = R.mPsi;
  }
// оператор присваивания
 TPartTrajInputData TPartTrajInputData::operator=(TPartTrajInputData  R)
{
     mTaskPointsZ  = R.mTaskPointsZ;
     mV = R.mV; // одномерный массив горизонтальных скоростей
     mCourse = R.mCourse;// одномерный массив курсовых углов
     menumType = R.menumType;
     mYaw = R.mYaw;
     mRadius = R.mRadius;
     mNumSetOfGears = R.mNumSetOfGears;
     mCourseEnd = R.mCourseEnd;
     mTHovering = R.mTHovering;
     m_dPsi_po_dt = R.m_dPsi_po_dt;
     //mPsi = R.mPsi;

    return *this ;
}

// парам констр 1
TPartTrajInputData::TPartTrajInputData(const  TPointXYZ   PointXYZ,const long double   V
  ,const long double   Course,const enumTypeOfManeuvre   TypeOfManeuvre,const long double   Yaw
               ,const long double   Radius, const int  NumSetOfGears,const long double   CourseEnd,const long double   THovering
                                       ,const long double   _dPsi_po_dt)
{
  mTaskPointsZ = PointXYZ;
  mV = V;
  mCourse = Course;
  menumType = TypeOfManeuvre;
  mYaw = Yaw;
  mRadius = Radius;
  mNumSetOfGears = NumSetOfGears;
  mCourseEnd = CourseEnd;
  mTHovering = THovering;
  //mPsi = Psi;
  m_dPsi_po_dt = _dPsi_po_dt;

}
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------



//-------------------------------------------------------------------------------
 //-------------------------------------------------------------------------------
 //-------- ОПИСАНИЕ КЛАССА TFlyTask -----------------------------------------------------------------------
 //-------------------------------------------------------------------------------
 //-------------------------------------------------------------------------------


TFlyTask::TFlyTask()
{
mquantTaskPoint = 0;
memset(mparrPartTrajInputData, 0, MAX_QUANT_TASK_POINTS * sizeof(TPartTrajInputData));
}


//---------------------------------------------------------------------------

 TFlyTask::~TFlyTask()
{

}
 //---------------------------------------------------------------------------

 //--------------------------

 // оператор присваивания
 TFlyTask TFlyTask::operator=(TFlyTask  R)
 {
     mquantTaskPoint  = R.mquantTaskPoint;
     memcpy(mparrPartTrajInputData , R.mparrPartTrajInputData , MAX_QUANT_TASK_POINTS *sizeof(TPartTrajInputData));
  return *this ;

 }

 // конструктор копирования
 TFlyTask::TFlyTask (const TFlyTask &R)
 {
     mquantTaskPoint  = R.mquantTaskPoint;
     memcpy(mparrPartTrajInputData , R.mparrPartTrajInputData , MAX_QUANT_TASK_POINTS *sizeof(TPartTrajInputData));

}

 // парам констр
 TFlyTask::TFlyTask( const int quantTaskPoint, long double   *parrCoord_GSK, long double   *parrV,
                  long double   *parrCourse, enumTypeOfManeuvre * parrenumTypes, long double   *parrYaw
                    , long double   *parrRadius, int *iarrNumSetOfGears, long double   *parrCourseEnd
                 ,long double   *parrTHovering,long double   *parr_dPsi_po_dt)
 {
    mquantTaskPoint = quantTaskPoint ;


    if (parrCoord_GSK && (quantTaskPoint >0 ))
    {

            for (int i =0; i < quantTaskPoint; i++)
            {
                mparrPartTrajInputData[i].mTaskPointsZ.X =  parrCoord_GSK[i * 3 + 1];
                mparrPartTrajInputData[i].mTaskPointsZ.Y = parrCoord_GSK[i * 3 + 2];
                mparrPartTrajInputData[i].mTaskPointsZ.Z =parrCoord_GSK[i * 3 ];
                mparrPartTrajInputData[i].mTaskPointsZ.M = 0.;

                mparrPartTrajInputData[i].mV = parrV[i];
                mparrPartTrajInputData[i].mYaw = parrYaw [i];
                mparrPartTrajInputData[i].mCourse = parrCourse[i];
                mparrPartTrajInputData[i].mRadius = parrRadius[i];
                mparrPartTrajInputData[i].menumType = parrenumTypes[i];
                mparrPartTrajInputData[i].mCourseEnd = parrCourseEnd[i];
                mparrPartTrajInputData[i].mTHovering = parrTHovering[i];
                mparrPartTrajInputData[i].mNumSetOfGears = iarrNumSetOfGears[i];

                mparrPartTrajInputData[i].m_dPsi_po_dt = parr_dPsi_po_dt[i];

            }

    }

 }





