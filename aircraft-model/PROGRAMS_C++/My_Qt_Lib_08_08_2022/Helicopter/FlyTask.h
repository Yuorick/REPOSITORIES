#ifndef FLYTASK_
#define FLYTASK_
#include "PartHelicTraj.h"

#define MAX_QUANT_TASK_POINTS 100
class TPointXYZ
{
public:

    long double   X  ; // X coordinate
    long double   Y  ;// Y coordinate
    long double   Z  ;// Z coordinate
    long double   M  ; // Measure
    TPointXYZ();
    TPointXYZ (const TPointXYZ &R);
    TPointXYZ operator=(TPointXYZ  R);
    TPointXYZ(const long double   X_,const long double   Y_,const long double   Z_, const long double   M_);
    TPointXYZ(const long double   X_,const long double   Y_,const long double   Z_);

};

class TPartTrajInputData
{
public:

    // для прямолинейного TLineMove: конечная точка участка траектории в НЗСК
     // для виража TTurnMove: конечная точка участка траектории в НЗСК
    // для висения THover: координаты ЦТ
    // для вращения TRotating: координаты ЦТ
    TPointXYZ mTaskPointsZ;

    //для прямолинейного TLineMove: скорость  в конечной точке
     // для виража TTurnMove:  скорость  в конечной точке
    // для висения THover: не используется (=0)
    // для вращения TRotating: не используется (=0)
    long double   mV;

    // для прямолинейного TLineMove: курсовой угол при движении по прямой
    // для виража TTurnMove:  курсовой угол конечной точки участка траектории
    // для висения THover: курсовой угол оси OX СвСК (= -угол рыскания в НЗСК)
                // например, надо, чтобы главная стрительная ось ЛА в конце маневра висения была направлена на восток
                // в ГСК курсовой угол равен 90 град.
                // следует иметь ввиду, что в этот момент угол рыскания в ГСК будет равен =-90 град
     // для вращения TRotating: не используется
    long double   mCourse;

    // тип маневра
    enumTypeOfManeuvre menumType;

    //для прямолинейного TLineMove: угол скольжения
     // для виража TTurnMove: угол скольжения
    // для висения THover: не используется
    // для вращения TRotating: не используется
    long double   mYaw;

    // для прямолинейного TLineMove: не используется
    // для виража TTurnMove:радиус виража
    // для висения THover: не используется
    // для вращения TRotating: не используется
    long double   mRadius;

    // номер набора передаточных чисел
    int mNumSetOfGears;

    // для прямолинейного TLineMove:не используется
     // для виража TTurnMove: величина изменения угла рыскания за время выполнения виража в ТНЗСК
    // для висения THover: не используется
    // для вращения TRotating: не используется
    long double   mCourseEnd;


     // для прямолинейного TLineMove:не используется
     // для виража TTurnMove: не исползуется
    // для висения THover: время висения
    // для вращения TRotating: время вращения
    long double   mTHovering;


     // для прямолинейного TLineMove:не используется
     // для виража: не используется
    // для висения THover: : не используется
    // для вращения TRotating: угловая скорость вращения
    long double m_dPsi_po_dt;



    TPartTrajInputData();
    TPartTrajInputData (const TPartTrajInputData &R);
    TPartTrajInputData operator=(TPartTrajInputData  R);
    TPartTrajInputData(const  TPointXYZ   PointXYZ,const long double   V
      ,const long double   Course,const enumTypeOfManeuvre   TypeOfManeuvre,const long double   Yaw
                   ,const long double   Radius, const int  NumSetOfGears,const long double   CourseEnd,const long double   THovering
                                           ,const long double   _dPsi_po_dt);
};
class TFlyTask
{
public:
    int mquantTaskPoint;
    TPartTrajInputData mparrPartTrajInputData[MAX_QUANT_TASK_POINTS];

    ~TFlyTask();
    TFlyTask();
    TFlyTask( const int quantTaskPoint, long double   *parrCoord_GSK, long double   *parrV,
                      long double   *parrCourse, enumTypeOfManeuvre * parrenumTypes, long double   *parrYaw
                        , long double   *parrRadius, int *iarrNumSetOfGears, long double   *parrCourseEnd
                     ,long double   *parrTHovering,long double   *parr_dPsi_po_dt);
    TFlyTask (const TFlyTask &R);
    TFlyTask operator=(TFlyTask  R);
};

#endif // FLYTASK_
