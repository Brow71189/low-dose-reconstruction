#ifndef HEX_EM_HPP
#define HEX_EM_HPP

#include "globals.h"
#include "fast_cache.h"
#include <fstream>



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
					);

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
						//double pval, int fm = 0, int lm = (globalModelNum-1)
					);




#endif
