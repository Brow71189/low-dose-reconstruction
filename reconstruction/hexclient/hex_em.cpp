#include "stdio.h"
#include "string.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include <assert.h>

#include "globals.h"
#include "fast_cache.h"
#include "fast_mpfr.h"
#include "transformations.hpp"
#include "hex_em.hpp"

void em_from_cache	(	
						ccflt *cctableBM, ssflt **ssmB,								
#ifdef USE_MPFR								
						ccint *cctableSM, ssint **ssmS,
#endif
						unsigned short *qImg,
						//ssflt* modelweight,
						int *chx , int *chy, 
						ssflt *meanB, ssflt *wghtB,
#ifdef USE_MPFR						
						ssint *meanS, ssint *wghtS,
#endif						
						int chlistlen, int active_m
					)
{
	const int modelsize(globalModelSize);
	const int sympersf(globalsympersf);
	const int subframeCount(globalSubFrameCount);
	const int modelnum(globalModelNum);
	const int qImgCol(globalqImgCol);
	int* const pixpos(pixposG);
	unsigned short * const inv_transformations(inv_transformationsG);
	//The weights must be identical for all changes as they belong to the same active_m
	ssflt *meanBp( meanB );
	ssflt *wghtBp( wghtB );
	ssflt wghtBA( 0.0 );//identical weights
#ifdef USE_MPFR	
	ssint *meanSp( meanS );
	ssint *wghtSp( wghtS );
	ssint  wghtSA ( 0 );//identical weights	
#endif	
	const ssflt *const ssB ( ssmB[active_m] );
#ifdef USE_MPFR				
	const ssint *const ssS ( ssmS[active_m] );
#endif				
	double *contribution = new double[subframeCount];
	int *best_m = new int[subframeCount];
	memset(best_m,-1,subframeCount*sizeof(int));
#ifdef USE_MAX	
	int *best_sym = new int[subframeCount];
	memset(best_sym,-1,subframeCount*sizeof(int));
#endif	
	
	
	

	ccflt *cctableBMp( cctableBM );
#ifdef USE_MPFR
	ccint *cctableSMp( cctableSM );
#endif		
	
	
	//STEP 1 populate best_m[] and contribution[] and calculate wghtB/SA
	
	for( int sfC(0);sfC < subframeCount; ++sfC)
	{		
		ssflt bestB( 0.0 );
		ssflt allB (0.0);
#ifdef USE_MPFR				
		ssint bestS( 0 );
		ssint allS ( 0 );
#endif				
		bool first_m = true;
		
		for(int m = 0; m < modelnum; ++m)
		{
			ssflt thisB( *(ssmB[m]+sfC) );
#ifdef USE_MPFR
					
			ssint thisS( *(ssmS[m]+sfC) );
			qfr_add(&allB,&allS,thisB,thisS);
			if( first_m || qfr_greater(thisB,thisS,bestB,bestS) )
			{
				first_m = false;
				bestB = thisB;
				bestS = thisS;
				best_m[sfC] = m;
			}	
#else									
			allB +=  thisB;
			if( first_m || (thisB > bestB) )
			{
				first_m = false;
				bestB = thisB;
				best_m[sfC] = m;
			}	
#endif					
		}
#ifdef USE_MPFR
		ssflt contributionB(ssB[sfC]);
		ssint contributionS(ssS[sfC]);
		qfr_div(contributionB,contributionS,allB,allS);
		
		
		const double contr( contributionB * pow(2.0,contributionS) );//between 0 and 1
		//HalfLorentz rocks
		contribution[sfC] = use_weights ?  2.0/(1.0+100.0*pow(1.1-contr,2.0)) : 1.0; 		
		if(use_weights || (best_m[sfC]==active_m) )
		{	qfr_add(&wghtBA, &wghtSA, contribution[sfC], 0);}
						
#else				
		ssflt contributionB( ssB[sfC]/allB );
		///Gaussstyle
		//contribution[sfC] = use_weights ?  pow(2.0,-100.0*pow( (1.0-contributionB), 2.0)) : 1.0 ; 
		//HalfLorentz rocks	
		contribution[sfC] = use_weights ?  2.0/(1.0+100.0*pow(1.1-contributionB,2.0)) : 1.0; 
				
		if(use_weights || (best_m[sfC]==active_m) )
		{	wghtBA += contribution[sfC];} 
		
#endif				
		
	
		//Precache the indices of maximal cctable entries
		//TODO align the pointers with the outer sfC loop
#ifdef USE_MAX		
		int best_ccind(-1);
		//int tmp_best_sym(-1); 
		if(use_max)
		{
			bool first = true;
			for (int sym(0); sym < sympersf; ++sym)
			{
				const int ccind(sfC + subframeCount*sym);
#ifdef USE_MPFR
				if( first || qfr_greater(cctableBMp[ccind],cctableSMp[ccind],cctableBMp[best_ccind],cctableSMp[best_ccind]) )
#else
				if( first || cctableBMp[ccind] > cctableBMp[best_ccind] )
#endif					
				{
					first = false;
					best_ccind = ccind;
					//tmp_best_sym = sym;
				}
			}
			
			best_sym[sfC] = (best_ccind/*-sfC*/) / subframeCount;
			//assert(best_sym[sfC]== tmp_best_sym);
#if 0
#ifndef USE_MPFR			
			if( ssB[sfC]!=cctableBMp[best_ccind] )
			{
				printf("sfC: %d\tprobvalB: %E\tssB: %E\n",sfC, cctableBMp[best_ccind], ssB[sfC]);
				fflush(stdout);
			}
			assert( ssB[sfC]==cctableBMp[best_ccind] );
#endif
#endif					
		}
		
		
#endif							

	} //end for sfC
	
	for (int chnum(0); chnum < chlistlen; ++chnum)
	{
		//TODO make these arrays of chnum
		const int q(chx[chnum]);
		const int r(chy[chnum]);
		const int indpp( pixpos[q + r * modelsize]);
		*meanBp = 0.0;
		*wghtBp = wghtBA;
#ifdef USE_MPFR
		*meanSp = 0;
		*wghtSp = wghtSA; 
#endif
		cctableBMp = cctableBM;
#ifdef USE_MPFR
		cctableSMp = cctableSM;
#endif		
		unsigned short* transfP( inv_transformations + (sympersf * indpp) );			
		for (int sym(0); sym < sympersf; ++sym) //aka lp && rot && mirror
		{
			
			//TODO make these arrays of chnum
			const int indpp_data( *(transfP++) );
			const int data_q( indpp_data % modelsize );
			const int data_r( indpp_data / modelsize );
			const unsigned short * datap( qImg + ( data_q * subframeCount + data_r * qImgCol ) );	
			
			for( int sfC(0); sfC < subframeCount; ++sfC)
			{
				const bool skip(!use_weights && (best_m[sfC] != active_m) );			
#ifdef USE_MAX
				if( skip || ( use_max && (best_sym[sfC] != sym) )  )
#else				
				if( skip )
#endif				
				
				{   // if there are no weights in use  consider only the best subframe
					++datap;
					++cctableBMp;
#ifdef USE_MPFR
					++cctableSMp;
#endif					
					continue; 
				}
				//+globalDatamin offset is now applied once by master after receiving the averaged value
				const unsigned short datav( *(datap++) );
				ssflt probvalB( *(cctableBMp++) );
#ifdef USE_MPFR				
				ssint probvalS( *(cctableSMp++));
				
				if(datav != 0)
				{	//dont do that for datav == 0
					qfr_div(probvalB,probvalS,ssB[sfC],ssS[sfC]);
					//const double rel_prob( probvalB * pow(2.0,probvalS ) );
					qfr_add(meanBp, meanSp, contribution[sfC] * probvalB * datav , probvalS); //rel_prob,  0
				}	
#else
				
				if(datav != 0)
				{	
					const double rel_prob(probvalB/ssB[sfC]);
					*meanBp += contribution[sfC] * rel_prob * datav;
				}
#endif	
					
			}// end for sfC
		} //end for sym			
		
#ifdef USE_MPFR
		qfr_norm(meanBp,meanSp);
		qfr_norm(wghtBp,wghtSp);
		
		
		++meanSp;
		++wghtSp;
#endif	
		++meanBp;
		++wghtBp;
	
	} //end for chnum
#ifdef USE_MAX	
	delete[] best_sym;
#endif	
	delete[] best_m;
	delete[] contribution;

}					














void emodel_from_cache	(	
						ccflt **cctableB, ssflt **ssmB,	ssflt *EModelsB, 
						ssflt *wghtB,												
#ifdef USE_MPFR								
						ccint **cctableS, ssint **ssmS, ssint *EModelsS,
						ssint *wghtS,
#endif
						unsigned short *zipImg
						//const unsigned short *qImg,
						//const unsigned short *Models,
						//const ssflt* modelweight,
						//double pval, int fm, int lm
					)
{					
	const int subFrameArea(globalSubFrameArea);
	unsigned short *sf = new unsigned short[subFrameArea];
	
	const int modelsize(globalModelSize);
	const int modelarea(modelsize*modelsize);
	const int sympersf(globalsympersf);
	const int subframeCount(globalSubFrameCount);
	const int modelnum(globalModelNum);
	const int hpMax(globalhpMax);
	const int*const hexpix(globalhexPixels);
	const unsigned short * const inv_transformations(inv_transformationsG);
	ssflt *const allB = new ssflt[subframeCount];
	memset(allB,0,subframeCount * sizeof(ssflt));
#ifdef USE_MPFR
	ssint *const allS = new ssint[subframeCount];
	memset(allS,0,subframeCount * sizeof(ssint));
#endif				
	unsigned char *best_m = new unsigned char[subframeCount];
	memset(best_m,0,subframeCount * sizeof(unsigned char));
	for( int sfC(0); sfC < subframeCount; ++sfC)
	{
		bool first = true;
		
		double m1(0.0);
		
		
		for(int m = 0; m < modelnum; ++m)
		{
#ifdef USE_MPFR
			ssflt thisB( *(ssmB[m]+sfC) );
			ssint thisS( *(ssmS[m]+sfC) );
			qfr_add(allB+sfC,allS+sfC,thisB,thisS);
			const double m0( log2(thisB) + (double)thisS );	
#else									
			allB[sfC] += ( ( *(ssmB[m]+sfC) ) );
			const double m0( ( *(ssmB[m]+sfC) ) ); 
#endif					
			if(first)
			{
				m1 = m0;
				first = false;
			}
			else if (m0 > m1) //a new best fit
			{
				m1 = m0;
				best_m[sfC] = m;
			}
			else if (m0 == m1) // a tie
			{
				if( (sfC%2 == 1 ) )
				{	best_m[sfC] = m;}
				
			}
		
		}	
	}
					
	//for(int m(fm) ; m <= lm ; ++m)				
	{				
	
		//int max_errors = 5;
		unsigned short* zipPtr( zipImg ); ///OUCH This was inside the sfC loop
		for( int sfC(0); sfC < subframeCount; ++sfC)
		{
			//collect the distributed pixels of the current subframe
			//Doing that for every model again is inefficient, but probably
			//still better than reading qImg
			memset(sf, 0, sizeof(unsigned short)*subFrameArea);
			unsigned short val( 0 );
			unsigned short elements( 0 );
			do
			{
				val = *(zipPtr++);
				elements = *(zipPtr++);
				for(unsigned short el(0); el < elements; ++el)
				{	sf[ hexpix[ *(zipPtr++) ] ] = val; }
			}
			while( (val > 0) || (elements > 0) );
			
			for(int mm(0); mm < modelnum; ++mm)
			{
				int m (best_m[sfC]);
				if( (!use_weights) && (mm != m) ) //single run loop for best value only
				{	continue;}
				if(use_weights)//effectively regular looping
				{	m = mm;}
				ssflt *emodelB(EModelsB+m*modelarea);
			
#ifdef USE_MPFR
				ssint *emodelS(EModelsS+m*modelarea);
				
#endif		
			
				//calculate this models share in the frames total likelyhood
				ssflt ssB( *(ssmB[m]+sfC) );
				ssflt oursB( ssB );		
#ifdef USE_MPFR	
				ssint ssS( *(ssmS[m]+sfC) );	
				ssint oursS( ssS );
#endif				
								
#ifdef USE_MPFR
				ssflt contributionB(oursB);
				ssint contributionS(oursS);
				qfr_div(contributionB,contributionS,allB[sfC],allS[sfC]);
				const double contr( contributionB * pow(2.0,contributionS) );//between 0 and 1
				const double contribution = use_weights ?  2.0/(1.0+100.0*pow(1.1-contr,2.0)) : 1.0 ; 
#else				
				ssflt contributionB( oursB/allB[sfC] );
				//const double contribution( contributionB );//between 0 and 1
				const double contribution = use_weights ?  2.0/(1.0+100.0*pow(1.1-contributionB,2.0)) : 1.0 ; 
#endif			
			
			
				ssflt *emB = new ssflt[modelarea];
				memset(emB, 0, sizeof(ssflt)*modelarea);
				const ccflt * cctableBMp( cctableB[m] + sfC);
#ifdef USE_MPFR			
				ssint *emS = new ssint[modelarea];
				memset(emS, 0, sizeof(ssint)*modelarea);
				const ccint * cctableSMp( cctableS[m] + sfC);					
#endif			
				
#ifdef USE_MAX
				int best_ccind(-1);
				if(use_max)
				{
					bool first = true;
					for (int sym(0); sym < sympersf; ++sym)
					{
						int ccind = sfC + subframeCount*sym;
#ifdef USE_MPFR
						if( first || qfr_greater(cctableBMp[ccind],cctableSMp[ccind],cctableBMp[best_ccind],cctableSMp[best_ccind]) )
#else
						if( first || cctableBMp[ccind] > cctableBMp[best_ccind] )
#endif					
						{
							first = false;
							best_ccind = ccind;
						}
					}
				}

#endif			
			
				double rel_prob(1.0);
				for (int sym(0); sym < sympersf; ++sym) //aka lp && rot && mirror
				{
#ifdef USE_MAX				
					int ccind = sfC + subframeCount*sym;
					if(use_max)
					{
						
						if(ccind != best_ccind)
						{	continue;} //skip all non maximal values	
					}
					
					ssflt probvalB( cctableBMp[ccind] );			
#ifdef USE_MPFR				
					ssint probvalS( cctableSMp[ccind] );
					cctableSMp += subframeCount;
					qfr_div(probvalB,probvalS,ssB,ssS);
					rel_prob = (probvalB * pow(2.0,probvalS));//between 0 and 1
#else				
					rel_prob = (probvalB /= ssB);//between 0 and 1				
#endif //USE_MPFR
					
#else // !USE_MAX								
				
					ssflt probvalB( *(cctableBMp) );
					cctableBMp += subframeCount;
#ifdef USE_MPFR				
					ssint probvalS( *(cctableSMp) );
					cctableSMp += subframeCount;
					qfr_div(probvalB,probvalS,ssB,ssS);
					rel_prob = (probvalB * pow(2.0,probvalS));//between 0 and 1
#else				
					rel_prob = (probvalB /= ssB);//between 0 and 1				
#endif //USE_MPFR

#endif //USE MAX		
				

					//now add the transformed sf to emodel
					for(int hpsf(0); hpsf < hpMax; ++hpsf)
					{
						const unsigned short indsf( inv_transformations[ TRANSSYM(hpsf,sym)] );
						const int datav( sf[indsf] );
						if( datav != 0)
						{
							const int emind( hexpix[hpsf] );
#ifdef USE_MPFR						
							qfr_add(emB+emind, emS+emind, rel_prob * datav, 0);
#else						
							emB[emind] += (rel_prob * datav);
#endif						
						}				
					}//end for hpsf
				}//end for sym
			
			//emB/S are now the probability weighted average of the symmetrized subframe

				
#ifdef USE_MAX				
				double wB( use_max?1:contribution );
#else
				double wB(  contribution );
#endif

#ifdef USE_MPFR				
				qfr_add(wghtB+m, wghtS+m, wB, 0);
#else				
				wghtB[m] += wB;
#endif			
				for(int ind(0); ind < modelarea; ++ind)
				{
					if(emB[ind] > 0.0)
					{	
#ifdef USE_MPFR				
						qfr_add(emodelB+ind, emodelS+ind, wB*emB[ind], emS[ind]);
#else
						emodelB[ind] += (wB * emB[ind]); //ssB
#endif						
					}
				}
			
				delete[] emB;
#ifdef USE_MPFR			
				delete[] emS;
#endif				

				}//end for m
		
			}//end for sfC				
	}
	
#ifdef USE_MPFR
	//normalize all calculated qfr numbers		
	for(int m(0); m < modelnum; ++m)
	{	
		qfr_norm(wghtB+m, wghtS+m);
		for(int ind(0); ind < modelarea; ++ind)
		{
			qfr_norm(EModelsB + m * modelarea + ind, EModelsS + m * modelarea + ind );
		}	
	}
#endif					
	delete[] sf;
	delete[] best_m;
	//delete[] sfchk;
	//delete[] sum_contrB;
	delete[] allB;
#ifdef USE_MPFR
	//delete[] sum_contrS;
	delete[] allS;
	

#endif	
	
	
					
}					
					
