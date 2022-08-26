#ifndef MEASURMENTIMITATOR_H
#define MEASURMENTIMITATOR_H


class QDriverMeasure
{
public:
    double marrYzv[3];
    double marrKYzv[9];
    QDriverMeasure();
    QDriverMeasure (const  QDriverMeasure &R);
    QDriverMeasure  &operator=( const QDriverMeasure  &R);
    QDriverMeasure (const double  *arrYzv,const double  *arrKYzv);

};
//-----------------------------------------------------------
//-----------------------------------------------------------
//-----------------------------------------------------------
//-----------------------------------------------------------
class QMeasurmentImitator
{
public:
    // погрешность измерения тока Id, %
    double mErrPerCent_Id;
    // погрешность измерения тока Iqu, %
    double mErrPerCent_Iqu;
    // СКЗ ошибки измерения кглового положения
    double mSgmTetta;

    QMeasurmentImitator();

    QMeasurmentImitator (const  QMeasurmentImitator &R);

    QMeasurmentImitator  &operator=( const QMeasurmentImitator  &R);

    QMeasurmentImitator (const double  ErrPerCent_Id , const double  ErrPerCent_Iqu,
        const double  SgmTetta);

    void createMeasure(double *arrPfVect,  QDriverMeasure *measureCur);

    int createInputDataReport(wchar_t*FileName, const bool bHeader);


};

#endif // MEASURMENTIMITATOR_H
