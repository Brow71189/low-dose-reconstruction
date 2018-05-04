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
#include "hexlikelyhood.hpp"
#include "transformations.hpp"

//is this one still a correct approximation with aberrations?
//probably not //TODO actually test that before expanding black from modelnum to 12*modelnum
void init_black(ssflt *blackB,
				//ssflt *blackL,
				//ssflt *blackK,
#ifdef USE_MPFR				
				ssint *blackS,
#endif				
				unsigned short *sModels,
				ssflt *lPtable,
				int fm, int lm )
{
	//const int modelnum(globalModelNum);
	const int modelarea(globalModelArea);
	const int shadownum( globalShadowNum );
	//memset(blackK,0,modelnum * sizeof(ssflt));
	for(int m(fm); m <= lm ; ++m)
	{
		unsigned short* model( sModels + (12*m) * modelarea * shadownum);
		ssflt sum(0.0);
		ssflt k(0.0);
		for(int hp(0); hp < globalhpMax; ++hp)
		{
			const int indModel ( globalhexPixels[hp] );	
			const double val( lPtable[ globalPtableWidth * model[ indModel ] ]); //datav = 0
			const double y( val - k );
			const double t( y + sum );
			k =	( t - sum ) - y;
			sum = t;	
		}
		//blackL[m] = sum;
#ifdef USE_MPFR
		ssint tmpS( 0 );
		qfr_exp2Z(sum, tmpS);
		
		blackB[m] = sum;
		blackS[m] = tmpS;
#else
		blackB[m] = (ssflt)exp2(sum);		
#endif		
	}
}

void update_black(ssflt *blackB,
				  //ssflt *blackL,
				  //ssflt *blackK,	
				  ssflt *nblackB,
				  //ssflt *nblackL,
				  //ssflt *nblackK,	
#ifdef USE_MPFR				
				ssint *blackS,
				ssint *nblackS,
#endif				
				ssflt *eitherPtable,
#ifndef MUL_PTBL				
				ssflt *InvPtable,	
#endif				
				int active_m,// int* chx, int* chy,
				int* newval, int* oldval, int chlistlen)
{
	const int ptablewidth(globalPtableWidth);
#ifdef MUL_PTBL	
	const int ptablesize(globalPtableSize);
#endif		
	ssflt probvalB( blackB[active_m] );
	//ssflt probvalL( blackL[active_m] );
	//ssflt probvalK( blackK[active_m] );
#ifdef USE_MPFR
	ssint probvalS( blackS[active_m] );
#endif
	
	for(int chnum(0);chnum < chlistlen; ++chnum)
	{
#ifdef MUL_PTBL
		const int mulPvalInd = oldval[chnum] * ptablewidth + newval[chnum] * ptablesize;
#else
		const int pwidthnewval(ptablewidth * newval[chnum]);
		const int pwidtholdval(ptablewidth * oldval[chnum]);
#endif
	
#ifdef MUL_PTBL
		probvalB *= eitherPtable[ mulPvalInd];// datav = 0
		//const ssflt x( eitherPtable[ mulPvalInd] );
#else
		probvalB *= ( eitherPtable[ pwidthnewval ] * InvPtable[ pwidtholdval ]  );// datav = 0
		//const ssflt x( eitherPtable[ pwidthnewval ] + InvPtable[ pwidtholdval ]  );
#endif
		/*const ssflt y( x - probvalK );
		const ssflt t( probvalL + y);
		probvalK = ( t - probvalL ) - y;
		probvalL = t;*/ 	
#ifdef USE_MPFR
		//probvalB = probvalL;
		//probvalS = 0.0;
		qfr_norm(&probvalB, &probvalS);									
#endif	
	}
	
	nblackB[active_m] = probvalB;
	//nblackL[active_m] = probvalL;
	//nblackK[active_m] = probvalK;
#ifdef USE_MPFR	
	nblackS[active_m] = probvalS;
#endif
	

}

//#if 0
void update_subsums_from_cache(	ccflt *cctableBM,
								//ccflt *cctableLM,
								//ccflt *cctableKM,
								ccflt *ncctableB,
								//ccflt *ncctableL,
								//ccflt *ncctableK,
								ssflt *nssmB,								
#ifdef USE_MPFR								
								ccint *cctableSM,
								ccint *ncctableS,
								ssint *nssmS,
#endif
								unsigned short *qImg,
#ifdef MUL_PTBL
								ssflt *mulPtbl,
#else
								ssflt *Ptable, ssflt *InvPtable,
#endif
								int *chx , int *chy, int chlistlen,
								int *newval, int *oldval)
{
#ifdef MUL_PTBL	
	const int ptablesize(globalPtableSize);
#endif	
#ifdef USE_MAX
	const bool usmx(use_max);
#endif	
	
	const int ptablewidth(globalPtableWidth);
	const int modelsize(globalModelSize);
	const int sympersf(globalsympersf);
	const int subframeCount(globalSubFrameCount);
	const int qImgCol(globalqImgCol);
	int* const pixpos(pixposG);
	unsigned short * const inv_transformations(inv_transformationsG);
	/*ssflt *nssmK = new ssflt[subframeCount];
	memset(nssmK,0,subframeCount*sizeof(ssflt));
	ssint *nssmKs = new ssint[subframeCount];
	memset(nssmKs,0,subframeCount*sizeof(ssint));*/
	for (int chnum(0); chnum < chlistlen; ++chnum)
	{
		//-funswitch-loops will optimize the innermost conditions
		const bool last_change(chnum == chlistlen - 1);
		const int q(chx[chnum]);
		const int r(chy[chnum]);
		const int indpp( pixpos[q + r * modelsize]);
#ifdef MUL_PTBL
		const int mulPvalInd = oldval[chnum] * ptablewidth + newval[chnum] * ptablesize;
#else
		const int pwidthnewval(ptablewidth * newval[chnum]);
		const int pwidtholdval(ptablewidth * oldval[chnum]);
#endif
		
		ccflt *ncctableBp( ncctableB );
		ccflt *cctableBMp( (chnum == 0) ? cctableBM : ncctableB );
		//ccflt *ncctableLp( ncctableL );
		//ccflt *ncctableKp( ncctableK );
		//first change initializes from old cache, all further from new cache
		//ccflt *cctableLMp( (chnum == 0) ? cctableLM : ncctableL );
		//ccflt *cctableKMp( (chnum == 0) ? cctableKM : ncctableK );
#ifdef USE_MPFR
		ccint *ncctableSp( ncctableS );
		ccint *cctableSMp( (chnum == 0) ? cctableSM : ncctableS);
#endif
		unsigned short* inv_transfP = inv_transformations + (sympersf * indpp);
		
		for (int sym(0); sym < sympersf; ++sym) //aka lp && rot && mirror
		{
			const int indpp_data( *(inv_transfP++) );
			const int data_q( indpp_data % modelsize );
			const int data_r( indpp_data / modelsize );				
			unsigned short *datap( qImg + ( data_q * subframeCount + data_r * qImgCol ) );

			static ssflt *const prob_changes = new ssflt[subframeCount];
			ssflt *prob_ch(prob_changes);
			for( int sfC(0); sfC < subframeCount; ++sfC)
			{
				const int datav( *(datap++) );
#ifdef MUL_PTBL
				*(prob_ch++) = mulPtbl[datav + mulPvalInd];
#else
				*(prob_ch++) = ( Ptable[ datav + pwidthnewval ] * InvPtable[ datav + pwidtholdval ]  );
#endif
			}
		
			prob_ch = prob_changes;	
			//In principle we only need to initialze for the last_change
			ssflt *nssmBp( nssmB );
			//ssflt *nssmKp( nssmK );
#ifdef USE_MPFR
			ssint *nssmSp( nssmS );
			//ssint *nssmKsp( nssmKs );
#endif
			
			for( int sfC(0); sfC < subframeCount; ++sfC)
			{
				//ssflt probvalL( *(cctableLMp++) );
				//ssflt probvalK( *(cctableKMp++) );
				ssflt probvalB( *(cctableBMp++) * *(prob_ch++) );
				//const int datav( *(datap++) );
#ifdef USE_MPFR				
				ssint probvalS( *(cctableSMp++) );
#endif				
/*
#ifdef MUL_PTBL
				probvalB *= mulPtbl[datav + mulPvalInd];
				//const ssflt x( mulPtbl[datav + mulPvalInd] );
#else
				probvalB *= ( Ptable[ datav + pwidthnewval ] * InvPtable[ datav + pwidtholdval ]  );
				//const ssflt x( Ptable[ datav + pwidthnewval ] + InvPtable[ datav + pwidtholdval ]  );
#endif
*/ 
				/*const ssflt y(x - probvalK);
				const ssflt t(y + probvalL);
				probvalK = (t - probvalL ) - y;
				probvalL = t;
				
				*(ncctableLp++) = probvalL;
				*(ncctableKp++) = probvalK;*/
				
#ifdef USE_MPFR
				qfr_norm(&probvalB, &probvalS);
				*(ncctableSp++) = (ccint)probvalS;
				*(ncctableBp++) = (ccflt)probvalB;
				if(last_change)
				{	
#ifdef USE_MAX		
				    if(usmx)
				    {					
						if( qfr_greater( probvalB, probvalS, *nssmBp, *nssmSp ) )
						{
							*nssmBp = probvalB;
							*nssmSp = probvalS;
						}
						++nssmBp;
						++nssmSp;					
					}
					else
#endif					
					{
						qfr_add(nssmBp, nssmSp, probvalB, probvalS);
						++nssmBp;
						++nssmSp;
						//qfr_kahan_add(nssmBp++, nssmSp++, probvalB, probvalS, nssmKp++, nssmKsp++ );
					}
					
				}
#else // NO_MPFR
				//ssflt probvalB = exp2(probvalL);
				if(last_change)
				{	
#ifdef USE_MAX
					*nssmBp = usmx?( (probvalB > *nssmBp) ? probvalB : *nssmBp ) : *nssmBp + probvalB;
					++nssmBp;
#else					
					/*const ssflt y( probvalB - *nssmKp );
					const ssflt t( y + *nssmBp );
					*(nssmKp++) = ( t - *nssmBp ) - y;
					*(nssmBp++) = t;*/
					*(nssmBp++) += probvalB;
#endif					
				}

				*(ncctableBp++) = (ssflt)probvalB;
#endif //MPFR			
			
			}//end for sfC								
			//delete[] prob_changes;
		} //end for sym
	} //end for chnum
	//delete[] nssmK;
	//delete[] nssmKs;
}

//#endif //0



int update_subsums_diff_cc(	ccflt *cctableBM,
								ccflt *ncctableB,
								ssflt *nssmB,								
#ifdef USE_MPFR								
								ccint *cctableSM,
								ccint *ncctableS,
								ssint *nssmS,
#endif
								unsigned short *qImg,
#ifdef MUL_PTBL
								ssflt *mulPtbl,
#else
								ssflt *Ptable, ssflt *InvPtable,
#endif
								short *diffmodel, unsigned short *smodel)
{
	
	//printf("entering update_subsums_diff_cc()\n");
	//fflush(stdout);
#ifdef MUL_PTBL	
	const int ptablesize(globalPtableSize);
#endif	
	const int ptablewidth(globalPtableWidth);
	const int modelsize(globalModelSize);
	const int modelarea(globalModelArea);
	const int shadownum(globalShadowNum); 
	const int sympersf(globalsympersf);
	const int subframeCount(globalSubFrameCount);
	const int qImgCol(globalqImgCol);
	const int difarea(12*shadownum*modelarea);
	const int lpMax(  globallpMax );
	//const int hpMax(  globalhpMax );
	int* const pixpos(pixposG);
#ifdef USE_MAX	
	const bool usmx = use_max;
#endif	
	unsigned short * const inv_translations(inv_translationsG);
	short *diffP(diffmodel);
	
	int num_changes(0);
	static ssflt *const prob_changes = new ssflt[subframeCount];
	
	static bool* fresh_blocks = new bool[12*shadownum];
		
	//for(int i = 0; i < 12*shadownum;++i)
	//{	fresh_blocks[i]=false; }
	std::fill_n(fresh_blocks,12*shadownum,true);
	
	
	for (int dif(0); dif < difarea; ++dif)
	{
		const short diffval = *(diffP++);
		if( diffval == 0 ){	continue;}
		++num_changes;
		const int ms( dif/modelarea	);
		const int shadow( ms / 12);
		//if(shadow == 1) continue; ///FIXME
		const int ind( dif % modelarea );
		const int mirror( (ms%12) / 6 );
		const int rot( ms % 6 );
		
		const int hp( pixpos[ind] );
		//assert(hp >= 0);
		//assert(hp <= hpMax);
		const unsigned short newval( smodel[dif] );
		if(newval < diffval)
		{
			printf("Shit: newval < diffval: %d < %d\n",newval, diffval);
		}
		const unsigned short oldval( newval - diffval);
		if( ptablewidth * newval > globalPtableSize || ptablewidth * oldval > globalPtableSize)
		{
			printf("bad newval: %d at dif: %d with ind: %d hp: %d ms: %d rot: %d mirror: %d\n", newval, dif, ind, hp, ms, rot, mirror);
			printf("or bad oldval: %d?\n", oldval);
			printf("and also bad diffval: %d\n",diffval);
			exit(0);
		}
		
		
		const bool fresh_block( fresh_blocks[ms] );
		//if(dif / (12*modelarea) == shadownum-1) ///HMM that should be the right way to do it
		{	fresh_blocks[ms] = false; } //after this time the updated values shall be used
		
		//printf("change in ms:%d at ind:%d  ov:%u -> nv:%u dv:%d fresh_block:%s\n",ms,ind,oldval,newval,diffval,(fresh_block?"true":"false"));
		//fflush(stdout);
		
		const int cc_start( CCIND(0,0,(rot+6*shadow),mirror) );
		//assert(cc_start + lpMax*subframeCount <= sfmrlG);
		ccflt *ncctableBp( ncctableB + cc_start);
		ccflt *cctableBMp( (fresh_block ? cctableBM : ncctableB) + cc_start );
#ifdef USE_MPFR
		ccint *ncctableSp( ncctableS + cc_start);
		ccint *cctableSMp( (fresh_block ? cctableSM : ncctableS) + cc_start );
#endif
		
		
		
#ifdef MUL_PTBL
		const int mulPvalInd = oldval * ptablewidth + newval* ptablesize;
#else
		const int pwidthnewval(ptablewidth * newval);
		const int pwidtholdval(ptablewidth * oldval);
		assert(pwidthnewval+ptablewidth-1 < globalPtableSize);
		assert(pwidtholdval+ptablewidth-1 < globalPtableSize);
#endif
		

		unsigned short* inv_translP( inv_translations + MOVE(hp,0) ); // MOVEHP0(hp)
		//assert(MOVE(hp,0)+lpMax <= hpMax * lpMax);
		
		for (int lp(0); lp < lpMax; ++lp)
		{
			const int ind_data( *(inv_translP++) );
			const int data_q( ind_data % modelsize );
			const int data_r( ind_data / modelsize );				
			unsigned short *datap( qImg + ( data_q * subframeCount + data_r * qImgCol ) );

			
			
				ssflt *prob_ch(prob_changes);
				for( int sfC(0); sfC < subframeCount; ++sfC)
				{
					const int datav( *(datap++) );
#ifdef MUL_PTBL
					*(prob_ch++) = mulPtbl[datav + mulPvalInd];
#else
					*(prob_ch++) = ( Ptable[ datav + pwidthnewval ] * InvPtable[ datav + pwidtholdval ]  );
#endif
				}

			
				prob_ch = prob_changes;	
				for( int sfC(0); sfC < subframeCount; ++sfC)
				{
					ssflt probvalB( *(cctableBMp++) * *(prob_ch++) );
#ifdef USE_MPFR				
					ssint probvalS( *(cctableSMp++) );
					qfr_norm(&probvalB, &probvalS);
					*(ncctableSp++) = (ccint)probvalS;
#endif			
					*(ncctableBp++) = (ccflt)probvalB;
				}//end for sfC
											
		} //end for lp
		//first_change = false;
	} //end for dif
	//delete[] prob_changes;
	//const bool no_change = first_change;
	
	//printf("updated cctable any_change: %s\n", (any_change?"true":"false") );
	//fflush(stdout);

	//update the subsums from the new cache values 
	
		
	int dead_updates(0);	
	for(int ms(0); ms < 12*shadownum; ++ms)
	{
		const int shadow = ms/12;
		//if(shadow == 1) continue; ///FIXME
		if(fresh_blocks[ms])
		{
			++dead_updates;
			const int mirror( (ms%12) / 6 );
			const int rot( ms % 6 );
			const int cc_start( CCIND(0,0,(rot+6*shadow),mirror) );
			const int block( lpMax * subframeCount);
			memcpy(ncctableB+cc_start,cctableBM+cc_start, block*sizeof(ccflt) );
#ifdef USE_MPFR		
			memcpy(ncctableS+cc_start,cctableSM+cc_start, block*sizeof(ccint) );
#endif		
		}

	}
	/*
	if( !( (dead_updates == 0) || (dead_updates == 12) ) )
	{
		printf("dead_updates: %d\n",dead_updates);
		fflush(stdout);
	}
	*/
	//assert( (dead_updates == 0) || (dead_updates == 12) );		

	ccflt const *cctableBMp( ncctableB );
#ifdef USE_MPFR
	ccint const *cctableSMp( ncctableS);
#endif
	const int lastsym ( shadownum * sympersf );
	for (int sym(0); sym < lastsym; ++sym) //aka lp && rot && mirror && shadow
	{
	
		ssflt *nssmBp( nssmB );
#ifdef USE_MPFR			
		ssint *nssmSp( nssmS );
#endif			
			
		for( int sfC(0); sfC < subframeCount; ++sfC)
		{	
			ssflt probvalB( *(cctableBMp++) );
#ifdef USE_MPFR				
			ssint probvalS( *(cctableSMp++) );
#endif							
						
#ifdef USE_MPFR
			
#ifdef USE_MAX		
				if(usmx)
				{					
					if( qfr_greater( probvalB, probvalS, *nssmBp, *nssmSp ) )
					{
						*nssmBp = probvalB;
						*nssmSp = probvalS;
					}
					++nssmBp;
					++nssmSp;					
				}
				else
#endif					
				{
					qfr_add(nssmBp, nssmSp, probvalB, probvalS);
					++nssmBp;
					++nssmSp;
				}
#else // NO_MPFR
			
#ifdef USE_MAX
				*nssmBp = usmx?( (probvalB > *nssmBp) ? probvalB : *nssmBp ) : *nssmBp + probvalB;
				++nssmBp;
#else					
				*(nssmBp++) += probvalB;
#endif					
							
#endif //MPFR			
		}//for sfC
	}//end for sym
	
	
	//printf("update_subsums_diff_cc() passed\n");
	//fflush(stdout);
	
	return (num_changes+6)/(12); //12 small changes are equivalent to one full traditional update
}




void init_subsums_from_black(	ssflt *blackB,
								//ssflt *blackL,
								//ssflt *blackK,
								ssflt **ssmB,
#ifdef USE_MPFR				
								ssint *blackS,
								ssint **ssmS,
#endif				
								ssflt *Ptbl0,
								unsigned short *zipImg,
								unsigned short *sModels,
								int fm, int lm ) //fm==-1 signals to just update lm
{
	const bool update( fm == -1 );
	fm = (update ? lm : fm);
	
	const int modelnum(globalModelNum);
	const int shadownum( globalShadowNum );
	const int subFrameCount(globalSubFrameCount);
	const int modelArea(globalModelArea);
	const int sympersf(globalsympersf);
	const int modsyms( (update ? 1 : modelnum) * sympersf);
	const int lpMax(globallpMax);
#ifdef USE_MAX
	const bool usmx = use_max;
#endif	
	unsigned short* zipPtr(zipImg);
	ssflt *const tmpcB( new ssflt[modsyms] );
	ssflt *const tmpcB2( new ssflt[modsyms] );
	{
		ssflt * tmpcBP2(  tmpcB2 );
		for(int m(fm) ; m<=lm ; ++m)
		{	tmpcBP2 = std::fill_n(tmpcBP2,sympersf,blackB[m]); }
	}
			
#ifdef USE_MPFR	
	ssint *const tmpcS( new ssint[modsyms] );
	ssint *const tmpcS2( new ssint[modsyms] );
	{
		ssint * tmpcSP2(  tmpcS2 );
		for(int m(fm) ; m<=lm ; ++m)
		{	tmpcSP2 = std::fill_n(tmpcSP2,sympersf,blackS[m]); }
	}
#endif	
	//work on subframes one by one, avoid the heavy cc caches
	for (int sfC(0); sfC < subFrameCount; ++sfC)
	{ 
		memcpy(tmpcB,tmpcB2,modsyms*sizeof(ssflt));
#ifdef USE_MPFR		
		memcpy(tmpcS,tmpcS2,modsyms*sizeof(ssint));
#endif		
/*		
		ssflt * tmpcBP(  tmpcB );
		ssint * tmpcSP( tmpcS );
		
		for(int m(fm) ; m<=lm ; ++m)
		{
			tmpcBP = std::fill_n(tmpcBP,sympersf,blackB[m]);
#ifdef USE_MPFR
			tmpcSP = std::fill_n(tmpcSP,sympersf,blackS[m]);
#endif			
		}
*/		
		//scan the different non-zero in the subframe
		unsigned short val(0);
		unsigned short elements(0);
		
		do
		{
			val = *(zipPtr++); //the currentvalue
			elements = *(zipPtr++); //its multiplicity
			for(unsigned short el(0); el < elements; ++el)
			{	
				
				
				ssflt *tmpcBP = (update ? tmpcB : (tmpcB + sympersf*fm) ); //rewind the temp cache Pointer //fm
#ifdef USE_MPFR				
				ssint *tmpcSP = (update ? tmpcS : (tmpcS + sympersf*fm) );		
#endif				
				const unsigned short indpp(*(zipPtr++) ); //hexpix index of val in Idata
				for(int m(fm); m <= lm; ++m)
				{
					//rewind transformations Pointer for every model
					/*** CK
					* NOTE sum over all inverse transforamtion is equivalent to sum
					* over all forward transformations. Numerically tested on 06.11.2014
					***/
					//unsigned short* transfP = transformationsG + (sympersf * indpp);
					const unsigned short* const model(sModels + (12*shadownum*m) * modelArea);
					
					for(int sym(0); sym < sympersf; ++sym)
					{
						
						const int lp( sym%lpMax ); //inv_translations
						const int sm( sym/lpMax ); //smodels with inverse rot and mirror
						const int sind(inv_translationsG[ MOVE(indpp,lp) ]); ///or inv_translationsG
						const unsigned short modelval( model[ sm*modelArea + sind] );
						/*
						if(*tmpcBP == 0.0L)
						{
							printf("Error: tmpcBP == 0\n");
							printf("m:%d   sym:%d   fm:%d   lm:%d   sf:%d\n", m, sym, fm, lm, sfC);
							printf("modsyms:%d    sympersf:%d   update: %s\n", modsyms, sympersf, update?"true":"false");
						}
						*/
						//assert(*tmpcBP != 0.0L);
						//assert(Ptbl0[modelval + val * globalModelval] != 0.0L);
						*tmpcBP *= Ptbl0[modelval + val * globalModelval];
						//assert(*tmpcBP != 0.0L);
						
#ifdef USE_MPFR					
						qfr_norm(tmpcBP, tmpcSP);
						//assert(*tmpcBP != 0.0L);
						++tmpcSP;
#endif					
						++tmpcBP;
					}
				}
			
			}
		}
		while( (val > 0) || (elements > 0) ); //finish reading a compressed frame
		
		const ssflt *tmpcBP( update ? tmpcB : (tmpcB + sympersf*fm) );
#ifdef USE_MPFR
		const ssint *tmpcSP( update ? tmpcS : (tmpcS + sympersf*fm) );
#endif
	
		for(int m(fm); m <= lm; ++m)
		{
			ssflt symsumB(0.0);	
#ifdef USE_MPFR	
			ssint symsumS(0);
#endif			
			for(int sym(0); sym < sympersf; ++sym)
			{	//sum over all symetric copies	
#ifdef USE_MPFR
#ifdef USE_MAX			
				if(usmx)
				{				
					if( qfr_greater(*tmpcBP, *tmpcSP, symsumB, symsumS) )
					{
						symsumB = *tmpcBP;
						symsumS = *tmpcSP;
					}
					++tmpcBP;
					++tmpcSP;
				}
				else
#endif			
				{	qfr_add(&symsumB, &symsumS, *(tmpcBP++), *(tmpcSP++) ); }				
#else // NO_MPFR
#ifdef USE_MAX
			
				symsumB = usmx?(((*tmpcBP)>symsumB) ? *tmpcBP : symsumB) : symsumB + *tmpcBP;
				++tmpcBP;	
#else			
				symsumB += *(tmpcBP++);
#endif	//USE_MAX
#endif	//MPFR		
			}
			
			//store to subsums and advance to next model
			if(update)
			{
				*((*ssmB) + sfC) = symsumB;
#ifdef USE_MPFR
				*((*ssmS) + sfC) = symsumS;
#endif		
			}
			else
			{
				*(ssmB[m] + sfC) = symsumB;
#ifdef USE_MPFR
				*(ssmS[m] + sfC) = symsumS;
#endif				
			}
		}
	}
	delete[] tmpcB;
	delete[] tmpcB2;
#ifdef USE_MPFR
	delete[] tmpcS;
	delete[] tmpcS2;
#endif		
}

void init_cache_from_black(
							ssflt *blackB,
							//ssflt *blackL,
							//ssflt *blackK,
							ccflt **cctableB,
							//ccflt **cctableL,
							//ccflt **cctableK,
							ssflt **ssmB,               
#ifdef USE_MPFR
							ssint *blackS,
							ccint **cctableS,
							ssint **ssmS,
					  
#endif
							ssflt *Ptbl0,
							unsigned short *zipImg,
							unsigned short *sModels,
							int fm, int lm )	
{
	//const int modelnum(globalModelNum);
	const int subFrameCount(globalSubFrameCount);
	const int modelArea(globalModelArea);
	const int sympersf(globalsympersf);
	const int shadownum( globalShadowNum );
	const int lpMax(globallpMax);
	unsigned short const* zipPtr(zipImg);
#ifdef USE_MAX	
	const bool usmx = use_max;
#endif	
	//unsigned short* inv_transformations = new unsigned short[globalhpMax * sympersf];
	//unsigned short* inv_transformations(inv_transformationsG);
	unsigned short const * const translations(translationsG);
	/*
	unsigned short *inv_transP(inv_transformations);
	//precache the inverse Transformations
	for(int hp(0); hp < globalhpMax; ++hp)//all hexpixels
	{
		for(int lp(0); lp < globallpMax; ++lp )//all latticepoints
		{
			for(int rot(0); rot < 6; ++rot) //six rotations
			{
				for(int mirror(0); mirror < 2; ++mirror)//two mirrors
				{
					*(inv_transP++) = (unsigned short)inv_transform_hp( hp, lp, rot, mirror );
				}	
			}
		}
	}
	*/
	for(int m(fm) ; m<=lm ; ++m)
	{
		//set ccpointers to model m
		ccflt *cctableBp( cctableB[m] );
		//ccflt *cctableLp( cctableL[m] );
		//ccflt *cctableKp( cctableK[m] );
		ssflt blackBv = blackB[m];
		//ssflt blackLv = blackL[m];
		//ssflt blackKv = blackK[m];
#ifdef USE_MPFR
		ccint *cctableSp( cctableS[m] );
		ssint blackSv = blackS[m];
#endif
		
		const int range(sympersf*subFrameCount);
		std::fill_n(cctableBp,range,blackBv); //Yeah 10% faster initialization
#ifdef USE_MPFR	
		std::fill_n(cctableSp,range,blackSv); //Yeah 10% faster initialization
#endif		
/*		
		//initialize all values with the (symetric) black values
		for(int sym(0); sym < sympersf; ++sym)
		{
			for (int sfC(0); sfC < subFrameCount; ++sfC)
			{
				*(cctableBp++) = (ccflt) blackBv;
				// *(cctableLp++) = (ccflt) blackLv;
				// *(cctableKp++) = (ccflt) blackKv;
#ifdef USE_MPFR			
				*(cctableSp++) = (ccint) blackSv;
#endif		
			}
		}
*/	
	}
	//scan the different non-zero in the subframe
	for (int sfC(0); sfC < subFrameCount; ++sfC)
	{
		unsigned short val(0);
		unsigned short elements(0);
		do
		{
			val = *(zipPtr++); //the currentvalue
			elements = *(zipPtr++); //its multiplicity
			for(unsigned short el(0); el < elements; ++el)
			{	
				const unsigned short indpp(*(zipPtr++) ); //hexpix index of val in Idata aka hp	
				for(int m(fm); m <= lm; ++m)
				{
					//rewind the transformations for every model
					//unsigned short* transfP = transformationsG + (sympersf * indpp);
					//set model m
					const unsigned short* model(sModels + (12*shadownum*m) * modelArea);
					//set ccpointers to model m
					ccflt* const cctableBM( cctableB[m] );
					//ssflt* const cctableLM( cctableL[m] );
					//ssflt* const cctableKM( cctableK[m] );
					
#ifdef USE_MPFR
					ccint* const cctableSM( cctableS[m] );
#endif
					for(int sym(0); sym < sympersf*shadownum; ++sym)
					{
						
						const int lp( sym % lpMax); // inv_translations
						const int sm( sym / lpMax); // smodels with inverse rot and mirror
						
						//const int sss( sfC + sym * subFrameCount );
						
						const int rot( sm %6 ); //
						const int mirror( (sm%12) / 6); //
						const int shadow( sm/12 ); 
						const int sss( CCIND(sfC,lp,(rot+6*shadow),mirror) );
						//assert( sss == sfC + sym * subFrameCount); //works out!
						const int sind( translations[ MOVE(indpp,lp)] );
						const unsigned short modelval( model[ sm*modelArea + sind ] );
						
						//get the value at the corresponding position inside the Model
						
						/*
						const unsigned short modelval2 = model[ transformationsG[ TRANSSYM(pixposG[sind],sym) ] ];
						if(modelval != modelval2)
						{
							printf("Shit happened\n");
							printf("modelval: %d != modelval2: %d\n",modelval,modelval2);
							printf("indpp: %d sym: %d   sm:%d   lp: %d\n", indpp, sym, sm, lp);
							
							
							int qs = sind % globalModelSize;
							int rs = sind / globalModelSize;
							printf("smodel sind: %d  q: %d r: %d\n", sind, qs, rs);
							
							fflush(stdout);
							abort();
						}
						
						///TODO test that assertion
						//assert(modelval == modelval2); 
						*/
						//update the Cache per model and sym for the non zero val from the subframe
						ssflt tmpB( cctableBM[sss] * Ptbl0[modelval + val * globalModelval]);
						//ssflt tmpL( cctableLM[sfC + sym * subFrameCount] );
						//ssflt tmpK( cctableKM[sfC + sym * subFrameCount] ); 
						
						/*const ssflt x( Ptbl0[modelval + val * globalModelval] );
						const ssflt y( x - tmpK );
						const ssflt t( tmpL + y);
						tmpK = (t - tmpL ) - y;
						tmpL = t;
						cctableKM[sfC + sym * subFrameCount] = (ccflt)tmpK;
						cctableLM[sfC + sym * subFrameCount] = (ccflt)tmpL;*/									
#ifdef USE_MPFR					
						ssint tmpS( cctableSM[sss] ); 
						qfr_norm(&tmpB,&tmpS);
						cctableSM[sss] = (ccint)tmpS;
#endif
						cctableBM[sss] = (ccflt)tmpB;																	
					}//end symmetries
				}//end models
			}//end non-zero esquence
		}
		while( (val!=0) || (elements!=0) ); //finish reading a compressed frame

		

		for(int m(fm) ; m <= lm ; ++m)
		{
			ssflt symsumB( 0.0 );
			//ssflt symsumK( 0.0 ); 
			ccflt* const cctableBM( cctableB[m] );
#ifdef USE_MAX			
			//ccflt* const cctableLM( cctableL[m] );
#endif			
#ifdef USE_MPFR
			ssint symsumS( 0 );
			//ssint symsumKs( 0 );
			ccint* const cctableSM( cctableS[m] );
#endif						
			for(int sym(0); sym < sympersf; ++sym)
			{
				const int sss( sfC + sym * subFrameCount );
#ifdef USE_MPFR
#ifdef USE_MAX				
				if(usmx)
				{
					if( qfr_greater(cctableBM[sss], 
						cctableSM[sss],
						symsumB, symsumS ) )
					//if(	cctableLM[sfC + sym * subFrameCount] > symsumS + log2(symsumB) )	
					{
						symsumB = cctableBM[sss];
						symsumS = cctableSM[sss];
					}					
				}	
				else
#endif				
				{	
					qfr_add(&symsumB, &symsumS, 
						cctableBM[sss],
						cctableSM[sss] );
				}		
				

#else //NO_MPFR				
#ifdef USE_MAX
				symsumB = usmx?(( cctableBM[sss] > symsumB )
						? cctableBM[sss] : symsumB) : 
						symsumB + cctableBM[sss];
#else //NO MAX				
				symsumB += cctableBM[sss];
				/*const ssflt val( cctableBM[sfC + sym * subFrameCount] ); 
				const ssflt y( val - symsumK);
				const ssflt t( symsumB + y);
				symsumK = ( t - symsumB ) - y;
				symsumB = t;*/
#endif //END MAX				
				
#endif	//END MPFR			
			}//end symmetries
		
			//store to subsums and advance to next model
			*(ssmB[m] + sfC) = symsumB;
#ifdef USE_MPFR
			*(ssmS[m] + sfC) = symsumS;
#endif				
		}//end models
	} //end subframes
	//delete[] inv_transformations; //ARGH
}



ssflt init_LPB_from_subsums(	ssflt **ssmB,
#ifdef USE_MPFR				
								ssint **ssmS,
#endif				
								ssflt* modelweight, int active_m)
{
	const int subframeCount(globalSubFrameCount);
	const int modelsp1(globalModelNum + 1);
	//const int modelnum(globalModelNum);
	ssflt totalLPBf( 0.0 );
	ssint totalLPBi( 0 );
	const bool usw(use_weights);
	//ssflt tot_K( 0.0 );
	//memset(cases,'\0',modelsp1*sizeof(int));
	for( int sfC(0); sfC < subframeCount; ++sfC)
	{
		//recalculate the whole sum
		bool model0(true);
		ssflt ELPB( 0.0 );
#ifdef USE_MPFR
		ssint ELPS( 0 );
#else		
		ssflt kahan( 0.0 );
#endif
		for( int m(0); m < modelsp1; ++m)
		{
			if( m == active_m ) //todo ++m here
			{
				continue;
			}

			if(usw)
			{
#ifdef USE_MPFR
				qfr_add(&ELPB, &ELPS,
							 (*( ssmB[m] + sfC )) * modelweight[m], *( ssmS[m] + sfC ) );
#else
				const ssflt val( (*( ssmB[m] + sfC )) * modelweight[m] );
				const ssflt y = (val - kahan);
				const ssflt t = (y + ELPB);
				kahan = (t - ELPB) - y;
				ELPB = t;
#endif
			}
			else 
			{
#ifdef USE_MPFR
				if(model0 || qfr_greater(*( ssmB[m] + sfC ), *( ssmS[m] + sfC ), ELPB, ELPS))
				{
					model0 = false;
					ELPB = *( ssmB[m] + sfC );
					ELPS = *( ssmS[m] + sfC );
				}
#else
				if( model0 || ( *( ssmB[m] + sfC ) > ELPB ) )
				{	
					model0 = false;
					ELPB = (*( ssmB[m] + sfC )); 
				}
#endif
			}	
			
		}
#ifdef USE_MPFR
		if(usw)
		{	qfr_norm(&ELPB,&ELPS); }
		totalLPBi += ELPS;		
#endif
		totalLPBf += log2(ELPB);
		/*const double dLPB( log2(ELPB) );
		const double y( dLPB - tot_K );
		const double t( y + totalLPBf );
		tot_K = (t - totalLPBf ) - y;
		totalLPBf = t;*/
	}
	return ( (totalLPBf + (ssflt)totalLPBi) * inv_pvalscaling);
}
/*
ssflt update_LPB_from_subsums(	ssflt *ssmBM,
								ssflt *nssmB,
								ssflt  *sflcB,
								ssflt *nsflcB,
#ifdef USE_MPFR				
								ssint *ssmSM,
								ssint *nssmS,
								ssint   *sflcS,
								ssint  *nsflcS,
#endif				
								ssflt weightM)
{
	const int subframeCount(globalSubFrameCount);
	ssflt totalLPB(0.0);
	//quick update + new - old
	for( int sfC(0); sfC < subframeCount; ++sfC)
	{
		ssflt ELPB( *sflcB );
#ifdef USE_MPFR
		ssint ELPS( *sflcS );
		ssflt tmpB (- *ssmBM);
		ssint tmpS (*ssmSM);
		qfr_add(&tmpB, &tmpS, *nssmB, *nssmS);
		tmpB *= weightM;
		qfr_add(&ELPB, &ELPS, tmpB, tmpS);
		qfr_norm(&ELPB, &ELPS);
		*nsflcB = ELPB;
		*nsflcS = ELPS;
		totalLPB += (log2(ELPB) + ELPS);
#else
		*nsflcB = ( ELPB += ( (*nssmB - *ssmBM) * weightM ) );
		totalLPB += log2(ELPB);
#endif
		++nsflcB;
		++nssmB;
		++ssmBM;
		++sflcB;
#ifdef USE_MPFR
		++nsflcS;
		++nssmS;
		++ssmSM;
		++sflcS;
#endif		
	
	}
	return (totalLPB * inv_pvalscaling);
}
*/

void update_cases(
                      ssflt **ssmB,
#ifdef USE_MPFR
                      ccint **ssmS,
#endif
                      int *cases
                )
{
	const int modelnum(globalModelNum);
	const int subFrameCount(globalSubFrameCount);
	memset(cases,0,modelnum*sizeof(int));
	for( int sfC(0); sfC < subFrameCount; ++sfC)
	{	
		ccflt highest_subsumB = *(ssmB[0] + sfC);	
#ifdef USE_MPFR
		ccint highest_subsumS = *(ssmS[0] + sfC);
#endif		
		int best_m(0);
		for (int m(1); m < modelnum; ++m)
		{
				ccflt subsumB = *(ssmB[m] + sfC);
#ifdef USE_MPFR
				ccint subsumS = *(ssmS[m] + sfC);
#endif
			if(
#ifdef USE_MPFR			
				qfr_greater(subsumB,subsumS,highest_subsumB,highest_subsumS)				
#else
				(subsumB > highest_subsumB)
#endif		
			)	
			{
				best_m = m;
				highest_subsumB = subsumB;	
#ifdef USE_MPFR
				highest_subsumS = subsumS;
#endif				
			}	
		}
		++(cases[best_m]);
	}
}


//Caution uses lPtable
ssflt hexlogprobability( const ssflt *lPtable,
			const unsigned short *sfData,
			const unsigned short *Model, //inv mirrored and rotated smodel
			const int lp) // const int mirror, const int rot
{
	const int hpMax( globalhpMax );
	const int lpMax( globallpMax );
	ssflt sum( 0.0 );
	ssflt kahan( 0.0 );
	int const* indModel( globalhexPixels );
	unsigned short const *itr( inv_translationsG + lp);
	for( int hpC( 0 ) ; hpC < hpMax ; ++hpC, ++indModel, itr+=lpMax) //++vali,
	{
		const ssflt val( lPtable[ sfData[ *itr ] + globalPtableWidth * Model[ *indModel ] ] );
		const ssflt y( val - kahan );
		const ssflt t( y + sum);
		kahan = ( t - sum ) - y;
		sum = t;	
	}
	return sum;
	/* ///Exact same timig as version above 
	ssflt sum( 0.0 );
	ssflt kahan( 0.0 );
	//int sf_sum( 0 );
	for(int hpC( 0 ); hpC < hpMax; ++hpC)
	{
		//const int indModel( globalhexPixels[hpC] );
		//const int indData( transformationsG[ TRANSFORM(hpC,lp,rot,mirror) ] );
		
		
		const int indModel = globalhexPixels[hpC];
		const int indData = inv_translationsG[ MOVE(hpC,lp) ];
		
		
		const int datav( sfData[ indData ] );
		//sf_sum += datav;
		const int modelv( Model[ indModel ] );
		const ssflt val( lPtable[ datav + globalPtableWidth * modelv ] );
		const ssflt y( val - kahan );
		const ssflt t( y + sum);
		kahan = ( t - sum ) - y;
		sum = t;
	}
	*/ 
	/*if(sum == 0)
	{
		printf("hexlogprobability: sf_sum == 0\n");
		abort();
	}*/
	
	//return sum;
}





ssflt init_hexlikelyhood(
                      ssflt **ssmB,
                      ccflt **cctableB,
                      //ccflt **cctableL,
                      //ccflt **cctableK,
#ifdef USE_MPFR
                      ssint **ssmS,
                      ccint **cctableS,
#endif //END MPFR
                      ssflt *lPtable,
                      unsigned short *zipImg,
                      unsigned short *sModels, ssflt *modelweight,
                      int fm, int lm)
{
	//const bool reinit_for_fun(::reinit_for_fun);
	//int max_errors = 20;
	//int bad_cc_entries = 0;
	const int subFrameArea(globalSubFrameArea);
	unsigned short *sf = new unsigned short[subFrameArea];
	//unsigned short *sf( (unsigned short*) alloca(subFrameArea*sizeof(unsigned short) ) ); //absolutely no difference
	//memset(sf0, 0, sizeof(unsigned short)*subFrameArea);
	ssflt totalLPBf(0.0);
	ssint totalLPBi(0);
	//ssflt tot_K(0.0);
	const int modelnum(globalModelNum);
	const int shadownum(globalShadowNum);
	const int subFrameCount(globalSubFrameCount);
	unsigned short* zipPtr = zipImg;
	const int modelArea(globalModelArea);
	//local handles for frequent global varibales in CCINDEX
	const int sfG = ::sfG;
	//const int sfmG = ::sfmG;
	//const int sfmrG = ::sfmrG;
	const int sfGlpM = ::sfGlpM;
	const int sfGlpMrot = ::sfGlpMrot;
	const bool usw(use_weights);
#ifdef USE_MAX
	const bool usmx(use_max);
#endif
	
	/*for(int m = 0; m< modelnum +1; ++m)
	{	memset(cctableK[m],0, sfmrlG*sizeof(ccflt) );}*/
	//memset(cases,0,(modelnum+1)*sizeof(int));
	//int counter = -1;
	
	for (int sfC(0); sfC < subFrameCount; ++sfC)
	{
		//collect the distributed pixels of the current subframe
		memset(sf, 0, sizeof(unsigned short)*subFrameArea);
		//memcpy(sf,sf0,sizeof(unsigned short)*subFrameArea);
		unsigned short val(0);
		unsigned short elements(0);
		//unsigned int sf_sum(0);
		do
		{
			val = *(zipPtr++);
			elements = *(zipPtr++);
			//sf_sum += val*elements;
			for(unsigned short el(0); el < elements; ++el)
			{	sf[globalhexPixels[*(zipPtr++)] ] = val;}
			
		}
		while( (val > 0) || (elements > 0) );
		/*if(sf_sum == 0)
		{
			printf("init__hexlikelyhood: sf_sum == 0 for sf: %d\n", sfC);	
		}*/
		
		ssflt ELPB(0.0);
		
#ifdef USE_MPFR
		ssint ELPS(0);
		//ssint ELPKs( 0 );
#else
		ssflt ELPK(0.0);
#endif //END MPFR
		bool model0(true);
		//int best_m(0);
		//int prev_fail(0);
		for (int m(0); m < modelnum; ++m)
		{
			const bool changed( m >= fm && m <=lm );
			ccflt *cctableBM = cctableB[m];
			//ccflt *cctableLM = cctableL[m];
			//ccflt *cctableKM = cctableK[m];
			ssflt subsumB( changed?0.0:*(ssmB[m] + sfC) );
			

#ifdef USE_MPFR
			ccint *cctableSM = cctableS[m];
			ssint subsumS( changed?0:*(ssmS[m] + sfC) );
			//ssint subsumKs( 0 );
//#else			
			//ssflt subsumK( 0.0 );		
#endif  //END MPFR
			if(changed)
			{
				//recalculate the cctableB/S and subsumB/S for that model
				for(int shadow = 0; shadow < shadownum; ++shadow) 
				{
					//if(shadow==1) continue; ///FIXME
					//int cc_order( CCIND(sfC,0,0,0) );
					for(int mirror(0); mirror < 2 ; ++mirror)
					{
						for(int rot(0); rot < 6; ++rot)
						{
							for(int lp(0); lp < globallpMax; ++lp)
							{
								//++counter;
								//FIXME extend cachtable to accomodate shadows
								const int ind( CCIND(sfC,lp,(rot+6*shadow),mirror) );
								/*
								if(ind != cc_order)
								{
									printf("ind: %d != cc_order: %d\n",ind , cc_order);
									printf("sfC: %d m: %d, mirror: %d, rot: %d, lp, %d\n", sfC, m, mirror, rot, lp);
									printf("sm: %d\n", (12*m + 6*mirror + rot) );
									printf("CCIND(sfC,0,0,0): %d\n" , CCIND(sfC,0,0,0));
									exit(0);
								}
								*/ 
								//assert(ind == cc_order);
								//cc_order += globalSubFrameCount;
								
								ssflt probval ( hexlogprobability( lPtable,
												sf,
												sModels + ( (12*m*shadownum + 12*shadow + (6*mirror + rot) ) * modelArea),
												lp ) ); //,0,0//mirror, rot										
#ifdef USE_MPFR
								ssint tmpS( 0 );
								qfr_exp2Z(probval, tmpS);
								/*
								if(reinit_for_fun)
								{
									if( abs((cctableBM[ind] / (ccflt)probval)-1.0)>1.0E-8 || (cctableSM[ind] != (ccint)tmpS) )
									{
										printf("CCTABLE entry error (old != new) : %lf != %lf\t || %d != %d\n", cctableBM[ind], (ccflt)probval, cctableSM[ind], (ccint)tmpS);
										printf("ind: %d sfC: %d m: %d, mirror: %d, rot: %d, lp, %d\n", ind, sfC, m, mirror, rot, lp);
										printf("sm: %d\n", (12*m + 6*mirror + rot) );
										printf("CCIND(sfC,0,0,0): %d\n" , CCIND(sfC,0,0,0));
										exit(0);
									}
									
								}
								*/
								/*
								if( ( fabs(1.0 - cctableBM[ind] / (ccflt)probval) < 0.000001 ) && (cctableSM[ind] == (ccint)tmpS) )
								{
									
									prev_fail = 0;
								}
								else
								{
									if(--max_errors > 0) 
									{
										printf("Match(%d) black_init != full_init  previous failures:%d\n", 20-max_errors,prev_fail);
										printf("Base (black!=full) : %lf != %lf\n",cctableBM[ind],(ccflt)probval);	
										printf("Shift (black!=full) : %d != %d\n",cctableSM[ind],tmpS);
										printf("at ind:%d \t sf: %d \t m:%d \t mirror: %d \t rot: %d \t lp: %d \n", ind, sfC, m, mirror, rot, lp); 	
									}
									
									
									++bad_cc_entries;
									++prev_fail;
								}
								*/
								/*
								if( (m==0) && ( cctableBM[ind]!=0.0 )
								{
									printf("dirty cctable entry: %lf\n", cctableBM[ind]);  
									printf("CCIND: %d\n", ind);
									printf("frame: %d\n",sfC);
									printf("model: %d\n",m);
									printf("shadow: %d\n",shadow);
									printf("mirror: %d\n",mirror);
									printf("rotation: %d\n",rot);
									printf("translation: %d\n",lp);
									printf("counter: %d\n", counter);
									if (cctableBM[ind]!=0.0) exit(0);
								}
								*/
								cctableBM[ind] = (ccflt)probval;
								cctableSM[ind] = (ccint)tmpS;
												
#ifdef USE_MAX						
								if(usmx)
								{							
									//probval is always >= 0.0
									if( qfr_greater(probval, tmpS, subsumB, subsumS) )
									{
										subsumB = probval;
										subsumS = tmpS;
									}						
								}
								else //DIRTY HACK to ensure identical initialization
						
#endif //END MAX						
								{	
									qfr_add(&subsumB, &subsumS, probval, tmpS);
									//qfr_norm(&subsumB, &subsumS);
								}
#else //NO MPFR
								probval = exp2(probval);
								/*
								if(reinit_for_fun)
								{
									if( cctableBM[ind] != (ccflt)probval )
									{
										printf("CCTABLE entry error (old != new) : %lf != %lf\n", cctableBM[ind], (ccflt)probval);
										printf("ind: %d  sfC: %d m: %d, mirror: %d, rot: %d, lp, %d\n", ind, sfC, m, mirror, rot, lp);
										printf("sm: %d\n", (12*m + 6*mirror + rot) );
										printf("CCIND(sfC,0,0,0): %d\n" , CCIND(sfC,0,0,0));
										exit(0);
									}	
								}
								*/ 
								cctableBM[ind] = (ccflt)probval;						
#ifdef USE_MAX
								subsumB = usmx?((probval > subsumB) ? probval : subsumB) : subsumB+probval;
#else //NO MAX						
								subsumB += probval;
								/*const ssflt y( probval - subsumK);
								const ssflt t( y + subsumB );
								subsumK = ( t - subsumB ) - y;
								subsumB = t;*/
#endif //END MAX					
#endif //END MPFR
							}
						}
					}
				}
				//store the recalculated subsumB/S
				*(ssmB[m] + sfC) = subsumB;
#ifdef USE_MPFR
				*(ssmS[m] + sfC) = subsumS;
#endif //END MPFR
			}//end changed 
			//subsumB/S are now either loaded or updated
			
			
#ifdef USE_MPFR
			if(usw)
			{
				ssflt subsB( subsumB * modelweight[m]);
				ssint subsS( subsumS);
				//qfr_norm(&subsB, &subsS);
				qfr_add(&ELPB, &ELPS, subsB, subsS);
			}
			else if(model0 || qfr_greater(subsumB, subsumS, ELPB, ELPS) ) 
			{
				model0 = false;
				ELPB = subsumB;
				ELPS = subsumS;
			}
			
#else //NO MPFR
			
			if(usw)
			{
				const ssflt x( subsumB * modelweight[m] );
				const ssflt y( x - ELPK );
				const ssflt t( y + ELPB );
				ELPK = ( t - ELPB ) - y;
				ELPB = t;
			}
			else if(model0 || (subsumB > ELPB))
			{	
				model0 = false;
				ELPB = subsumB;
			}
			
#endif //END MPFR

		
		} //end for m
		//calculate log-prob and add to overall sum
		//++cases[best_m];
#ifdef USE_MPFR
		if(usw)
		{	qfr_norm(&ELPB, &ELPS);}
		totalLPBi += ELPS; //integer contributions		
#endif
		totalLPBf += log2(ELPB);
		/*const ssflt dLPB( log2(ELPB) );
		const ssflt y( dLPB - tot_K );
		const ssflt t( y + totalLPBf );
		tot_K = ( t - totalLPBf ) - y;
		totalLPBf = t;*/
	}
	//printf("bad_cc_entries: %d\n", bad_cc_entries);
	delete[] sf;
	return ( (totalLPBf+(ssflt)totalLPBi) * inv_pvalscaling);
}

double cache_report( std::ofstream &reportfile,
                      ssflt *blackB,
                      //ssflt *blackL,
                      //ssflt *blackK,
                      ssflt **ssmB,
#ifdef USE_CCTABLE                      
                      ccflt **cctableB,
                      //ccflt **cctableL,
                      //ccflt **cctableK,
#endif                      
#ifdef USE_MPFR
                      ssint *blackS,
                      ssint **ssmS,
#ifdef USE_CCTABLE                       
                      ccint **cctableS,
#endif					  
#endif
					int* cases, bool reporting)
{
	
	update_cases(
				ssmB,
#ifdef USE_MPFR
				ssmS,
#endif
				cases);
	const int modelnum(globalModelNum);
	const int sfG(::sfG);

	//get the final state statistics
#ifdef USE_CCTABLE	
	double  minlog(0.0), maxlog(0.0), sumccB(0.0), sumccS(0.0);
	//double min_k(0.0), max_k(0.0), sum_k(0.0);
	bool first_cc( true );
#endif 	
	double minssmS(0.0),  maxssmS(0.0),  sumssB(0.0), sumssS(0.0);
	bool first_ss( true );

	for(int m (0); m < modelnum; ++m) 
	{
#ifdef USE_CCTABLE	
		ccflt *ccB( cctableB[m] );
		//ccflt *ccL( cctableL[m] );
		//ccflt *ccK( cctableK[m] );
#endif		
		ssflt *ssB( ssmB[m] );
#ifdef USE_MPFR
#ifdef USE_CCTABLE		
		ccint *ccS( cctableS[m] );
#endif		
		ssint *ssS( ssmS[m] );
#endif

#ifdef USE_CCTABLE			
		for( int i(0); i < sfmrlG; ++i) //cache table
		{
#ifdef USE_MPFR
			ssflt ccBv( log2(*(ccB)) + *(ccS++) );
#else
			ssflt ccBv( log2(*(ccB)) );
#endif
			//ssflt ccLv( *(ccL++) );
			//ssflt ccKv( *(ccK++) );
			int trash;
			sumccB += frexp( *(ccB++), &trash);
			//sumccL += ccLv;
			sumccS += ccBv;
			//sum_k += ccKv;
			
			minlog = ( first_cc || (minlog > ccBv) ) ? ccBv : minlog ;
			maxlog = ( first_cc || (maxlog < ccBv) ) ? ccBv : maxlog ;
			//minlog2 = ( first_cc || (minlog2 > ccLv) ) ? ccBv : minlog2 ;
			//maxlog2 = ( first_cc || (maxlog2 < ccLv) ) ? ccBv : maxlog2 ;
			//min_k = ( first_cc || (min_k > ccKv) ) ? ccKv : min_k ;
			//max_k = ( first_cc || (max_k < ccKv) ) ? ccKv : max_k ;
			
			first_cc = false;
		}
#endif
		for( int sfC(0); sfC < sfG; ++sfC) // subsums
		{
#ifdef USE_MPFR
			ssflt ssBv( log2( *(ssB) ) + *(ssS++) );
#else
			ssflt ssBv( log2(*(ssB)) );
#endif
			int trash;
			sumssB += frexp( *(ssB++), &trash);
			sumssS += ssBv;
			minssmS = ( first_ss || (minssmS > ssBv) ) ? ssBv : minssmS ;
			maxssmS = ( first_ss || (maxssmS < ssBv) ) ? ssBv : maxssmS ;
			first_ss = false;
		}
	}

	//some Idiot defined pvaloffset with negative sign, ok typicall values will hence be positive
#ifdef USE_CCTABLE	
	double best_ln_offset = pvaloffset+0.5*log(2.0)*(maxlog+minlog)/globalhpMax;
	double next_maxlog = 0.5*(maxlog-minlog);
	if(next_maxlog > CC_EXP_MAX)
	{
		best_ln_offset += (next_maxlog-CC_EXP_MAX)/globalhpMax;
	}
#endif
	std::stringstream report;
	//might be usefull for debugging
	
	
	report << "cases:";
	for(int m(0); m < modelnum; ++m) 
	{	report << "\t" << cases[m] << "\t";}
	report << "\n";
  
#ifdef USE_MPFR		
	report << "black propabilities Exp2:";
	for(int m(0); m < modelnum; ++m) 
	{	report << "\t" << blackS[m] << "\t";}
	report << "\n";
/*	
	report << "black propabilities Log2:";
	for(int m(0); m < modelnum; ++m) 
	{	report << "\t" << blackL[m] << "\t";}
	report << "\n";
*/	
#endif
	report << "black propabilities Base:";
	for(int m(0); m < modelnum; ++m) 
	{	report << "\t" << blackB[m];}
	report << "\n";
		
	
	

#ifdef USE_CCTABLE		
	report << "cctable: sumB = "<< sumccB << "\t Min ln2 = " << minlog << 
		"\t Max ln2 = " << maxlog  << "\t Sum ln2 = " << sumccS << "\n";
	/*report << "cctableL: Sum log2 = " << sumccL << "\t Min log2 = " <<  minlog2 << 
		"\t Max log2 = " << maxlog2  << "\n";*/
	/*report << "cctableK: Sum K = " << sum_k << "\t Min K = " <<  min_k << 
		"\t Max K = " << max_k  << "\n";*/		
#endif		
	report << "subsums: sumB = " << sumssB << "\t Min ln2 = " << minssmS <<
		"\t Max ln2 = " << maxssmS << "\t Sum ln2 = " << sumssS << "\n";

#ifdef USE_CCTABLE		
	if(pvaloffset_ok)
	{
		report << "current pvaloffset: " << pvaloffset << ", suggesting delta pvaloffset: " << best_ln_offset << "\n";
	}
	else
	{
		report << "Initial pvaloffset is UNKOWN, but should be changed by: " <<  best_ln_offset << "\n";
	}
#endif	
	if(reporting) {std::cout << report.str();}
	if(reportfile.is_open())
	{	reportfile << report.str();}
	
#ifdef USE_CCTABLE	
	return best_ln_offset;
#else
	return 0;
#endif	
}
//////////////////////////////////////////////TRASH SECTION

#if 0

void update_subsums_from_black(	ssflt *nblackB,
								//ssflt *nblackL,
								//ssflt *nblackK,
								ssflt *nssmB,
#ifdef USE_MPFR				
								ssint *nblackS,
								ssint *nssmS,
#endif				
								ssflt *Ptbl0,
								unsigned short *zipImg,
								unsigned short *sModels, int active_m )
{
	const int subFrameCount(globalSubFrameCount);
	const int modelArea(globalModelArea);
	const int sympersf(globalsympersf);
	const int lpMax(globallpMax);
	//const int ptablewidth(globalPtableWidth);
	const bool usmx(use_max);
	unsigned short* zipPtr(zipImg);
	
	ssflt *tmpcB = new ssflt[sympersf];
#ifdef USE_MPFR		
	ssint *tmpcS = new ssint[sympersf];
#endif	
	//work on subframes one by one, avoid the heavy cc caches
	for (int sfC(0); sfC < subFrameCount; ++sfC)
	{
		//initalize likelyhoods with black values
		
		ssflt *tmpcBP( tmpcB );
		/*ssflt *tmpcL = new ssflt[sympersf];
		ssflt *tmpcLP( tmpcL );
		ssflt *tmpcK = new ssflt[sympersf];
		ssflt *tmpcKP( tmpcK );*/
#ifdef USE_MPFR		
		
		ssint *tmpcSP( tmpcS );
#endif		 
		
		for(int sym(0); sym < sympersf; ++sym)
		{
			*(tmpcBP++) = nblackB[active_m];
			//*(tmpcLP++) = nblackL[active_m];
			//*(tmpcKP++) = nblackK[active_m];
#ifdef USE_MPFR			
			*(tmpcSP++) = nblackS[active_m];
#endif				
		}
		/*
		if(sfC%100==0)
		{	
			printf("tmpcBP and tmpCS for Frame %d were initialzed\n", sfC);
			fflush(stdout);
		}*/
		//scan the different non-zero in the subframe
		unsigned short val(0);
		unsigned short elements(0);
		//unsigned int sf_sum(0);
		do
		{
			//assert(zipPtr - zipImg < zip_img_cap);
			val = *(zipPtr++); //the currentvalue
			//assert(zipPtr - zipImg < zip_img_cap);
			elements = *(zipPtr++); //its multiplicity
			
			for(unsigned short el(0); el < elements; ++el)
			{	
				//assert(zipPtr - zipImg < zip_img_cap);
				const unsigned short indpp(*(zipPtr++) ); //pix index of val in Idata
				tmpcBP = tmpcB; //rewind the temp cache Pointer
				//tmpcLP = tmpcL;
				//tmpcKP = tmpcK;
#ifdef USE_MPFR				
				tmpcSP = tmpcS;		
#endif				
				/***
				//rewind tranformations Pointer for every model
				// use forward instead of inverse transformations sum is invariant
				// Numerically tested on 06.11.2014 by CK
				
				Now changed to sModels that require inv- transformations 19.04.2016
				
				***/
				//unsigned short* transfP = transformationsG + (sympersf * indpp);
				const unsigned short* model(sModels + active_m * modelArea);
				
				for(int sym(0); sym < sympersf; ++sym)
				{
					const int lp = sym%lpMax; //inv_translations
					const int sm = sym/lpMax; //smodels with inverse rot and mirror
					
					
						
					const int sind( inv_translationsG[ MOVE(indpp,lp)] ); /// or inv_translationsG
					const unsigned short modelval( model[ sm*modelArea + sind ] );
					
					
					//get the value at the corresponding position inside the Model
					//const unsigned short modelval = model[ *(transfP++) ];
					//update the temp Cache per model and sym for the non zero val from the subframe						
					*tmpcBP *= Ptbl0[modelval + val * globalModelval];
					/*const ssflt x( Ptbl0[modelval + val * globalModelval] );
					const ssflt y( x - *tmpcKP );
					const ssflt t( y + *tmpcLP );
					*tmpcKP = ( t - *tmpcLP ) - y;
					*tmpcLP = t;*/
#ifdef USE_MPFR					
					/*ssflt xb(*tmpcLP);
					ssint xs(0.0);
					qfr_exp2Z(xb,xs);
					*tmpcBP = xb;
					*tmpcSP = xs;*/
					qfr_norm(tmpcBP, tmpcSP);
					++tmpcSP;
/*#else					
					*tmpcBP = exp2(*tmpcLP);*/		
#endif					
					++tmpcBP;
					//++tmpcLP;
					//++tmpcKP;
				}
			}
		}
		while( (val > 0) || (elements > 0) ); //finish reading a compressed frame
		/*if(sf_sum == 0)
		{
			printf("update_subsums_from_black: sf_sum == 0 at sfC: %d\n", sfC);
			abort();
		}*/
	
	
		tmpcBP = tmpcB;
		//tmpcLP = tmpcL;
		//tmpcKP = tmpcK;
		ssflt symsumB( 0.0 );
		//ssflt symsumK( 0.0 );
#ifdef USE_MPFR
		tmpcSP = tmpcS;
		ssint symsumS( 0 );
		//ssint symsumKs( 0 );
#endif		
		
		for(int sym(0); sym < sympersf; ++sym)
		{	//sum over all symetric copies	
#ifdef USE_MPFR
#ifdef USE_MAX			
			if(usmx)
			{
				if( qfr_greater(*tmpcBP, *tmpcSP, symsumB, symsumS) )
				//if( *tmpcLP > symsumS + log2(symsumB) )
				{
					symsumB = *tmpcBP;
					symsumS = *tmpcSP;
				}
				++tmpcBP;
				//++tmpcLP;
				//++tmpcKP;
				++tmpcSP;			
			}			
			else
#endif			
			{	
				qfr_add(&symsumB, &symsumS, *(tmpcBP++), *(tmpcSP++) );
			}
#else //NO MPFR
#ifdef USE_MAX
			symsumB = usmx?((*tmpcBP > symsumB) ? *tmpcBP : symsumB) : symsumB + *tmpcBP;
			++tmpcBP;
#else			
			
			symsumB += *(tmpcBP++);
			/*const ssflt x( *(tmpcBP++) );
			const ssflt y( x - symsumK );
			const ssflt t( y + symsumB );
			symsumK = (t - symsumB) - y;
			symsumB = t;*/	
#endif
#endif	//MPFR		
		}
		
		//delete[] tmpcL;
		//delete[] tmpcK;
		//write symsum into subsums 
		*(nssmB++) = symsumB;
		
#ifdef USE_MPFR
		*(nssmS++) = symsumS;
#endif						
	}
	delete[] tmpcB;
#ifdef USE_MPFR	
	delete[] tmpcS;
#endif	
}

#endif
