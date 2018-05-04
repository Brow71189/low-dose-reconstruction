#include "stdio.h"
#include <assert.h>
#include <stdlib.h>
#include "string.h"
#include "globals.h"
#include "fast_mpfr.h"
#include "fast_cache.h"
#include "communicator.hpp"

ssflt *histoB(nullptr);
unsigned long int *histoC(nullptr); 

#ifdef USE_MPFR
ssint *histoS(nullptr);
#endif

//Currently only counting encounters, no weighting with cctable
//also how to account for multiplicity of grayvalues


void create_Histo(
				ccflt **cctableB,
#ifdef USE_MPFR
                ccint **cctableS,
#endif
				unsigned short *qImg,
				unsigned short *Models
			)
{
	const int modelsize( globalModelSize );
	const int modelarea( globalModelArea );
	const int modelnum( globalModelNum );
	const int modelval( globalModelval ); //aka graylevels
	const int ptablewidth( globalPtableWidth );
	const int ptablesize( globalPtableSize );
	const int sympersf( globalsympersf );
	const int subframeCount( globalSubFrameCount );
	const int qImgCol( globalqImgCol );
	//int* const pixpos( pixposG );
	unsigned short * const inv_transformations( inv_transformationsG );
	delete[] histoB;
	histoB = new ssflt[ptablesize];
	memset(histoB,0,ptablesize*sizeof(ssflt));
	delete[] histoC;
	histoC = new unsigned long int[modelval];
	memset(histoB,0,modelval*sizeof(unsigned long int));
#ifdef USE_MPFR
	delete[] histoS;
	histoS = new ssint[ptablesize];
	memset(histoS,0,ptablesize*sizeof(ssint));
#endif	

	for(int m(0); m < modelnum ; ++m)
	{
		unsigned short* model( Models + m * modelarea );
		for(int hp(0); hp < globalhpMax; ++hp)
		{
			const int indpp( hp );
			const int indModel( globalhexPixels[indpp] );
			//assert(indModel >= 0);
			//assert(indModel < modelarea);
			//assert(indpp == pixpos[indModel]);
			const int lamda( model[indModel] ); //the current expectation value 
			//printf("m:%d hp:%d indModel:%d lamda:%d\n",m, hp,indModel,lamda);
			//fflush(stdout);
			++histoC[lamda];
			ccflt * cctableBMp( cctableB[m] );
#ifdef USE_MPFR		
			ccint * cctableSMp( cctableS[m] );
#endif	
			unsigned short* transfP = inv_transformations + (sympersf * indpp);
			
			for (int sym(0); sym < sympersf; ++sym) //aka lp && rot && mirror
			{
				const int indpp_data( *(transfP++) );
				//assert(indpp_data >= 0);
				//assert(indpp_data < modelarea);
				const int data_q( indpp_data % modelsize );
				const int data_r( indpp_data / modelsize );				
				unsigned short *datap( qImg + ( data_q * subframeCount + data_r * qImgCol ) );
				
				for( int sfC(0); sfC < subframeCount; ++sfC)
				{
					//assert(datap-qImg >= 0);
					//assert(datap-qImg < globalqImgSize);
					const int datav( *(datap++) );
					const int histpp( datav + lamda * ptablewidth);
					//assert(histpp >= 0);
					//assert(histpp < ptablesize);
					//assert(cctableBMp - cctableB[m] >= 0);
					//assert(cctableBMp - cctableB[m] < globalqImgCol);
					ssflt probvalB( *(cctableBMp++) );
													
#ifdef USE_MPFR				
					ssint probvalS( *(cctableSMp++) );
					qfr_add( histoB + histpp, histoS + histpp, probvalB, probvalS);
										
#else
					histoB[ histpp ] += probvalB;
					
#endif		
				}//end for sfC		
				
			}//end for sym
			
		}//end for hp

	}//end for m
	for(int ind(0); ind < ptablesize; ++ind)
	{
		const int lamda( ind/ptablewidth);
		const unsigned long hC( histoC[lamda]*qImgCol);		
#ifdef USE_MPFR
		
		const ssflt iC( hC ? (ssflt)(1.0 / (double)hC) : 1.0 );
		qfr_mul_d(histoB[ind],histoS[ind],iC);
#else		
		histoB[ind] /= ( hC ? (double)hC : 1.0 );	
#endif	
	
	}
}

void report_Histo(bool reporting)
{
	assert(histoB != nullptr);
#ifdef USE_MPFR	
	assert(histoS != nullptr);	
#endif	
	
	sendHisto(  histoB, 
#ifdef USE_MPFR
				histoS,
#endif
				reporting
				);	
}




void clear_Histo()
{
	delete[] histoB;
	histoB = nullptr;
	delete[] histoC;
	histoC = nullptr;
#ifdef USE_MPFR
	delete[] histoS;
	histoS = nullptr;
#endif			
}


