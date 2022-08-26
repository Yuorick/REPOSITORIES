#ifndef SOLVER_LBL_H
#define SOLVER_LBL_H
#include "PosSolver.h"
#include <QVector>


class QSolver_LBL : public QPosSolver
{
public:
  // QSolver_LBL();

    QSolver_LBL (const  QSolver_LBL &R);

    QSolver_LBL  &operator=( const QSolver_LBL  &R);

    QSolver_LBL  (const QVector<QAbstractMeasure> vctMeasures
   ,const TTable_1D tblEstPrfl
,const double*arrGPSAprioriPosParams
, const int DimX, const int DimY);
};

#endif // SOLVER_LBL_H
