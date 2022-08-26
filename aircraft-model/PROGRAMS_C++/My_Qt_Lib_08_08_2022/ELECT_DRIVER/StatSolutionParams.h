#ifndef STATSOLUTIONPARAMS_H
#define STATSOLUTIONPARAMS_H
#include "ElectDriver.h"

class QStatSolutionParams
{
public:
    double marrStatPhVect[QVARS];
    double marrStatU[2];
    double marrGears[2 * QVARS];

    QStatSolutionParams();

    QStatSolutionParams  (const QStatSolutionParams &R);

    QStatSolutionParams &operator=(const QStatSolutionParams  &R);

    QStatSolutionParams(  const double *arrStatPhVect,const double* arrStatU,const double*arrGears);

};

#endif // STATSOLUTIONPARAMS_H
