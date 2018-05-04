#include "stdio.h"
#include "string.h"
#include <sstream>
#include <string>
#include <math.h>

#include "fast_cache.h"
#include "globals.h"
#include "fast_mpfr.h"
#include "likelyhood.hpp"


//Caution assumes lPtable is the ln of Ptable
void probability( double &base, int &probexp, const double *lPtable,
			 const unsigned char *sfData, const int &subframesize, const int &datamin,
			 const unsigned short *Model, const int &modelsize,
			 const int x0, const int y0)
{
	base = 0.0;
	int tmpexp;
    //probexp = 0;
	int ptablewidth = globalPtableWidth;
	for (int y = 0; y < modelsize; ++ y)
    {
        int yoffsmod = modelsize * y;
        int yoffsdata=  y + y0;
        if ((yoffsdata) >= subframesize)
        {
            yoffsdata -= subframesize;
        }
        yoffsdata *= subframesize;
        for (int x = 0; x < modelsize; ++x)
        {
            int xoffsdata = x + x0;
            if ( xoffsdata >= subframesize)
            {
                xoffsdata -= subframesize;
            }
            int datav = sfData[ xoffsdata + yoffsdata ];
            int modelv = Model[x + yoffsmod];
            base += lPtable[ (datav-datamin) + ptablewidth * modelv ];
        }

    }
	base = frexp( exp(base), &tmpexp);
	probexp = tmpexp;
}

void init_likelyhood(  double &totalLPB, int modelnum, int subframeCount,
					   int oCount, int xoCount, int stepsize, double *ssmB, int *ssmS,
					   double *nssmB, int *nssmS, ccflt *cctablebase,
                       ccint  *cctableexp, ccflt *ncctablebase, ccint *ncctableexp,
                       double *lPtable, unsigned char *qImg,
                       int subframesize, int datamin, unsigned short* Models,
                       int modelsize, double *modelweight )
{
	int modelarea    = modelsize * modelsize;
	int subframearea = subframesize * subframesize;
	int moCount      = oCount * modelnum;
	
	for (int sfC = 0; sfC < subframeCount; ++sfC)
	{
		double ELPB(0.0);
		int ELPS(0);
		double newsubsumB;
		int newsubsumS;
		for (int m = 0; m < modelnum; ++m)
		{			
			int ssmid = QSSMIND(sfC,m);
			double *subsumB = &ssmB[ssmid];
			int *subsumS = &ssmS[ssmid];
			*subsumB = 0.0;
			*subsumS = 0;
			for( int oC = 0; oC < oCount ; ++oC )
			{
				int yoC = oC / xoCount;
				int xoC = oC % xoCount;					
				int yo = yoC * stepsize;
				int xo = xoC * stepsize;
				
				//TODO: replace by something that is not specific for the present lattice
				int yoadj = ( xoC & 1 ) ? ( stepsize / 2 ) : 0;
				int probexp;
				double probval;
				probability( probval, probexp, lPtable, qImg + (sfC * subframearea), //0.0 is for log Pval
							 subframesize, datamin, Models + (m * modelarea), modelsize,
							 xo, yo + yoadj );				
				const unsigned int index = QQCCIND(sfC,m,oC);
				cctablebase[ index ] = probval;
				cctableexp [ index ] = probexp;
				ncctablebase[ index ] = probval;
				ncctableexp [ index ] = probexp;
				/*
				if (reporting && sfC==0 && probdebug==m)
				{
					printf("xoC,yoC %d,%d\tprob. %f << %d\n", xoC, yoC, probval, probexp);
				}
				*/ 
				qfr_add(*subsumB, *subsumS, probval, probexp);
			}
			nssmB[ssmid] = *subsumB;
			nssmS[ssmid] = *subsumS;
			qfr_mul_d(newsubsumB, newsubsumS, *subsumB, *subsumS, modelweight[m]);
			qfr_add(ELPB, ELPS, newsubsumB, newsubsumS);			
		} //end of summation over models (new)
		//calculate log-prob and add to overall sum		
		qfr_log(ELPB, ELPS);		
		totalLPB += ELPB;
	}
}


