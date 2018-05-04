#ifndef HEXLIKELYHOOD_HPP
#define HEXLIKELYHOOD_HPP

#include "globals.h"
#include "fast_cache.h"
#include <fstream>

void init_black(ssflt *blackB,
				//ssflt *blackL,
				//ssflt *blackK,
#ifdef USE_MPFR				
				ssint *blackS,
#endif				
				unsigned short *sModels,
				ssflt *lPtable,
				int fm = 0, int lm = (globalModelNum-1) );


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
				int* newval, int* oldval, int chlistlen);

void update_black_diff(ssflt *blackB,
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
				int active_m, 
				short *diffmodel, unsigned short *smodel);

                      
                      
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
								int fm = 0, int lm = (globalModelNum-1) );
							
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
							ssflt *Ptable0,
							unsigned short *zipImg,
							unsigned short *sModels,
							int fm = 0, int lm = (globalModelNum-1) );							
/*
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
								unsigned short *sModels, int active_m );
*/
//#if 0
void update_subsums_from_cache(	ccflt *cctableBM,
								//ccflt *cctbaleLM,
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
								//unsigned short *active_model,
#ifdef MUL_PTBL
								ssflt *mulPtbl,
#else
								ssflt *Ptable, ssflt *InvPtable,
#endif
								int *chx , int *chy, int chlistlen,
								int *newval, int *oldval);
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
								short *diffmodel, unsigned short *smodel);

								




ssflt init_LPB_from_subsums(	ssflt **ssmB,
#ifdef USE_MPFR				
								ssint **ssmS,
#endif				
								ssflt* modelweight, int active_m);
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
								ssflt weightM);
*/

void update_cases(
                      ssflt **ssmB,
#ifdef USE_MPFR
                      ccint **ssmS,
#endif
                      int *cases
                );


ssflt hexlogprobability(const ssflt *Ptable,
						const unsigned short *sfData,
						const unsigned short *Model,
						const int lp ); //const ,int mirror, const int rot

ssflt init_hexlikelyhood(
                      ssflt **ssmB,
                      ccflt **cctableB,
                      //ccflt **cctableL,
                      //ccflt **cctableK,
#ifdef USE_MPFR
                      ssint **ssmS,
                      ccint **cctableS,
#endif
                      ssflt *lPtable,
                      unsigned short *zipImg,
                      unsigned short *sModels, ssflt *modelweight,
                      int fm = 0, int lm = (globalModelNum-1) );


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
					int* cases, bool reporting);
#endif
