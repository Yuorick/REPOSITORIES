//---------------------------------------------------------------------------

#ifndef GaussH
#define GaussH
	double   getRand01( );

        double   getGauss( const  double   a, const  double   sig) ;

	double    Rand_();

        void  CalcStatParams(const double valX, const int N
	, double &valAver,double  &valAverSquare);

	 void __fastcall getGaussVector_dim3(double *arrMean, double *arrF, double *arrMtrxLamb, double *arrOut);

	 void __fastcall MtrxMultMatrx_(double *parrA,int nRowsA, int nColsA, double * parrB,int nColsB, double *parrRez);

	 void __fastcall MtrxSumMatrx_(double *pA, double * pB,int nRows, int nCols, double *pRez);

	 void __fastcall MtrxTranspMultMatrx_(double *parrA,int nRowsA, int nColsA, double * parrB,int nColsB, double *parrRez) ;

	 void __fastcall getGaussVector(int iDim, double *arrMean, double *arrF, double *arrMtrxLamb, double *arrOut);

#endif
