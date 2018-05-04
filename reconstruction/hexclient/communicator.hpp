#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

#include "globals.h"
#include "fast_cache.h"

/* These functions read data from stdin in the format that is expected
 * from the Hex_Magic Java ImageJ Plugin
 */
//#define GETCHAR getc(jobfile?jobfile:stdin)

//short getshortBE(void);
//short get16bitint(void);
//unsigned short getushort(void);
double getdouble(FILE* input);
void senddouble(double);
ssflt getWeights(FILE* input,ssflt *weight);
void sendWeights(ssflt *weight);
int getModels(FILE* input, unsigned short* Models, int fm = 0, int lm = (globalModelNum-1) );
int getBeam(FILE* input, unsigned short* Beam);
int sendModels(unsigned short* Models, int fm = 0, int lm = (globalModelNum-1));
int sendEModels(ssflt* EModelsB, ssflt* wghtB,
#ifdef USE_MPFR
				ssint* EModelsS, ssint* wghtS, 
#endif
				bool reporting,
				int fm = 0, int lm = (globalModelNum-1));
int sendHisto(
				ssflt* histoB,
#ifdef USE_MPFR
				ssint* histoS,
#endif
				bool reporting
			 );
#ifdef MUL_PTBL
ssflt getPtable(FILE* input, ssflt *lPtable, ssflt *mulPtbl, ssflt *Ptbl0);
#else
ssflt getPtable(FILE* input, ssflt *lPtable, ssflt *Ptable, ssflt *InvPtable, ssflt *Ptbl0);
#endif

#ifdef USE_CCTABLE
int getIdata(FILE* input, unsigned short *qImg);
#else //NO CC_TABLE
int getIdata(FILE* input, unsigned short *zipImg, long &zipImgUsed);
#endif


int getZipIdata(FILE* input, unsigned short *zipImg, long &zipImgUsed);



//float getfloat()


#endif
