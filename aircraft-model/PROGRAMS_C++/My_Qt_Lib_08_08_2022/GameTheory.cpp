//---------------------------------------------------------------------------


#pragma hdrstop
#include <float.h>
#include <string.h>
#include "GameTheory.h"
#include "MatrixProccess.h"
#include "LinOptimization.h"

// решение матричной игры в смешанных стратегиях
// первый игрок выбирает строку матрицы i, второй столбец j
// первый максимизирует, второй минимизирует
//
// INPUT:
// arrMtrx - массив элементов матрицы
// NUmRows - число строк матрицы
// NUmCols - число столбцов  матрицы
// OUTPUT:
// arr_x - вектор смешанных стратегий 1-го игрока, размерность NUmRows
// arr_y - вектор смешанных стратегий второго игрока, размерность  NUmCols
// *fval  - цена игры
// если arr_x ==NULL, то вектор стратегий первого игрока не формируется
// если arr_y ==NULL, то вектор стратегий второго  игрока не формируется
void TGameTheory::solvMartrxGame(double *arrMtrx, const int NUmRows, const int NUmCols
									, double *arr_x, double *arr_y, double &fval)
{
	 if (arr_x != NULL)
	 {
		 // формирование исх данных для задачи лин программирования
			//1.1 матрица ограничений
			int NumArgMin = -1;
			double valMin = MinDoubleArray(arrMtrx, NUmRows * NUmCols, &NumArgMin) ;
			double val_add = (valMin <= 0.)? -valMin + 1.: 0.;
			double *parr_a = new double[NUmRows * NUmCols];
			MatrTransp(arrMtrx, NUmRows, NUmCols, parr_a);

			for (int i =0; i < NUmRows * NUmCols; i++)
			{
			 parr_a [i] = -(parr_a [i] + val_add);

			}
			///

			// 1.2 формирование вектора правой части ограничений
				double *parr_b = new double[NUmRows];
				for (int i = 0; i < NUmRows; i++)
				{
					parr_b[i] = -1.;
				}
				///

			// 1.3 формирование вектора целевой функции
				double *parr_f = new double[NUmRows];
				for (int i = 0; i < NUmRows; i++)
				{
					parr_f[i] = 1.;
				}
				///

			 // 1.4 формирование вектора ограничений снизу
					double *parr_lb = new double[NUmRows];
					memset(parr_lb, 0,NUmRows * sizeof(double));
					///

				 // 1.5 формирование вектора ограничений сверху
				 double *parr_ub = new double[NUmRows];
				 for (int i = 0; i < NUmRows; i++)
				 {
					 parr_ub[i] = DBL_MAX /2.;
				 }
				 ///

				 ///

				 // 1.7 вспомогат переменные
				 int  nrows_eq = 0; // к-во равенств
				 double a_eq[1] = {1.  // матрица равенств  не используется
													};
				 double b_eq[1] = {1.  // правая часть равенств  не используется
													 };

				int nvars = NUmRows;  // к-во переменных
				int nrows =  NUmCols; // к-во стограничений в виде неравенств рок в матрице
				int bvars =0;       // к-во булевых переменных
				int ix [1] ={0}; // вспомогательный вектор булевых переменных

				 LinNumericalSolver(  nvars, bvars, parr_f,  nrows,
						parr_a, parr_b, nrows_eq,	a_eq, b_eq, parr_lb,
							parr_ub,ix, arr_x, fval) ;

				for (int i = 0; i < nvars; i++)
				{
				 arr_x [i] =   arr_x [i] / fval;
				}
			fval = 1. / fval -val_add;
			delete parr_a;
			delete parr_b;
			delete parr_f;
			delete parr_lb ;
			delete parr_ub;
	 }

		if (arr_y != NULL)
	 {
			int nvars = NUmCols;  // к-во переменных
			int nrows =  NUmRows; // к-во ограничений в виде неравенств рок в матрице
		 // формирование исх данных для задачи лин программирования
			//1.1 матрица ограничений
			int NumArgMin = -1;
			double valMin = MinDoubleArray(arrMtrx, NUmRows * NUmCols, &NumArgMin) ;
			double val_add = (valMin <= 0.)? -valMin + 1.: 0.;
			double *parr_a = new double[NUmRows * NUmCols];
			memcpy(parr_a, arrMtrx, NUmRows * NUmCols * sizeof(double));
			for (int i =0; i < NUmRows * NUmCols; i++)
			{
			 parr_a [i] = (parr_a [i] + val_add);
			}
			///

			// 1.2 формирование вектора правой части ограничений
				double *parr_b = new double[nrows];
				for (int i = 0; i < nrows; i++)
				{
					parr_b[i] = 1.;
				}
				///

			// 1.3 формирование вектора целевой функции
				double *parr_f = new double[ nvars];
				for (int i = 0; i < nvars; i++)
				{
					parr_f[i] = -1.;
				}
				///

			 // 1.4 формирование вектора ограничений снизу
					double *parr_lb = new double[ nvars];
					memset(parr_lb, 0,NUmRows * sizeof(double));
					///

				 // 1.5 формирование вектора ограничений сверху
				 double *parr_ub = new double[ nvars];
				 for (int i = 0; i <  nvars; i++)
				 {
					 parr_ub[i] = DBL_MAX /2.;
				 }
				 ///
				 

				 // 1.7 вспомогат переменные
				 int  nrows_eq = 0; // к-во равенств
				 double a_eq[1] = {1.  // матрица равенств  не используется
													};
				 double b_eq[1] = {1.  // правая часть равенств  не используется

													 };

				int bvars =0;       // к-во булевых переменных
				int ix [1] ={0}; // вспомогательный вектор булевых переменных

				LinNumericalSolver(  nvars, bvars, parr_f,  nrows,
						parr_a, parr_b, nrows_eq,	a_eq, b_eq, parr_lb,
							parr_ub,ix, arr_y, fval) ;

				for (int i = 0 ; i < nvars; i++)
				{
				 arr_y [i] =   -arr_y [i] / fval;
				}
			 //	fval -=  val_add;
				fval = -fval;
				fval = 1. / fval -val_add;
			delete parr_a;
			delete parr_b;
			delete parr_f;
			delete parr_lb ;
			delete parr_ub;
	 }
}
//---------------------------------------------------------------------------

#pragma package(smart_init)
