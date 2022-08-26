//---------------------------------------------------------------------------


#pragma hdrstop

#include "Atmosphere.h"
	extern	const  double HT[5][3]=
	{   // высота темпер градиент
		 0.0, 289.00, -0.006328,
		 9324.0, 230.00, -1.0,
		12000.0, 221.50,  0.0,
		29600.0, 221.50,  0.0028,
		47000.0, 269.66,  0.0028
	};

//------------------------------------------------------------------------------
void fncCalcNormTemperature(const  double H, double &valHT,  double &valGradHT)
{ int i;

	 double gr0 = HT[0][2], a = -HT[0][2]/( HT[2][0] -HT[1][0]);

	if(H<=HT[0][0]) { valHT=289.0; valGradHT=HT[0][2]; return; }

	for(i=0;i<4;i++)
	{
	 if((H > HT[i][0])&&(H <= HT[i+1][0])) break;
	}

	switch(i)
	{ case 0:

			valHT= HT[0][1]+ (H-HT[0][0])* HT[0][2];
			valGradHT=HT[0][2];
		break;

		case 1:

		   valHT = HT[1][1] +  gr0*(H-HT[1][0]) + a * (H-HT[1][0]) * (H-HT[1][0]) / 2. ;
		   valGradHT = gr0 + a * (H-HT[1][0]);
		break;

		case 2:
			valHT= HT[2][1];
			valGradHT=0;
		break;

		case 3:

			valHT= HT[3][1]+ HT[3][2]*(H - HT[3][0]);
			valGradHT= HT[3][2];
		break;

		case 4:
			valHT=     HT[4][1];
			valGradHT= 0.0;
		break;

	}

}
#pragma package(smart_init)
