#ifndef OPTCTRLSYST_H
#define OPTCTRLSYST_H



class QLinOptCtrlSyst
{
public:
    ~QLinOptCtrlSyst();
    QLinOptCtrlSyst();
    QLinOptCtrlSyst (const  QLinOptCtrlSyst &R);
    QLinOptCtrlSyst  &operator=( const QLinOptCtrlSyst  &R);
    QLinOptCtrlSyst (const   int QuantXVar,
         const int mQuantUVar,
         double *mpA,
         double *mpB,
         double *mpC,
         const double T0,
         double *marrPhVect0,
         const double mTCur,
         double *arrPhVect);

    int mQuantXVar;
    int mQuantUVar;
    double *mpA;
    double *mpB;
    double *mpC;
    double mT0;
    double *marrPhVect0;
    double mTCur;
    double *marrPhVect;
};

#endif // OPTCTRLSYST_H
