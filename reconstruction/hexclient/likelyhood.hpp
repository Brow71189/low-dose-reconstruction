#ifndef LIKELYHOOD_HPP
#define LIKELYHOOD_HPP

#include "fast_cache.h"

/* functions for calculatin pixel probabilities and model likelyhoods  
 * 
 */

void hexagonal_trap(int &dx, int &dy, int &dz, const int &N);


void probability( double &base, int &probexp, const double *Ptable, const int &ptablewidth,
			 const unsigned char *sfData, const int &subframesize, const int &datamin,
			 const unsigned short *Model, const int &modelsize,
			 const int x0, const int y0);
			 
void hexprobability( double &base, int &probexp, const double *Ptable, const int &ptablewidth,
			 const unsigned char *sfData, const int &subframesize, const int &datamin,
			 const unsigned short *Model, const int &modelsize,
			 const int x0, const int y0);



void init_likelyhood( double &totalLPB, int modelnum, int subframeCount,
                      int oCount, int xoCount, int stepsize, double *ssmB, int *ssmS,
                      double *nssmB, int *nssmS, ccflt *cctablebase,
                      ccint  *cctableexp, ccflt *ncctablebase, ccint *ncctableexp,
                      double *lPtable, unsigned char *qImg,
                      int subframesize, int datamin, unsigned short *Models,
                      int modelsize, double *modelweight );
                      
void init_hexlikelyhood( double &totalLPB, int modelnum, int subframeCount,
                      int stepsize, double *ssmB, int *ssmS,
                      double *nssmB, int *nssmS, ccflt *cctablebase,
                      ccint  *cctableexp, ccflt *ncctablebase, ccint *ncctableexp,
                      double *lPtable, unsigned char *qImg,
                      int subframesize, int datamin, unsigned short *Models,
                      int modelsize, double *modelweight );

#endif
