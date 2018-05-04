
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "string.h"
#include "stdlib.h"
#include "globals.h"
#include "fast_cache.h"
#include "communicator.hpp"
#include "zippedimg.hpp"
/*
//BigEndian Encoding
short getshortBE()
{
	short one = (short)getchar();
	short two = (short)getchar();
	return ( (one << 8) | (two) );
}

short get16bitint()
{
  short first = (short)getchar();
  first += ((short)getchar())*256;
  return first;
}

unsigned short getushort()
{
  unsigned short first = (unsigned short)getchar();
  first += ((unsigned short)getchar())*256;
  return first;
}
*/




double getdouble(FILE* input)  //we have to reverse the byte order from JAVA (at least for linux/gcc)
{
	
	
	union 
	{
		double dbl;
		char byts[sizeof(double)/sizeof(char)];  	
	};
	for(size_t i = sizeof(double)/sizeof(char); i > 0u; --i)
	{		byts[i-1] = getc(input);}
	return dbl;
}

void senddouble(double val)  //we have to reverse the byte order for JAVA (at least from linux/gcc)
{
  union 
	{
		double dbl;
		char byts[sizeof(double)/sizeof(char)];  	
	};
  dbl = val;
  
  for(size_t i = sizeof(double)/sizeof(char); i > 0u; --i)
	{	putchar(byts[i-1]);}
  
  /*
  char *ret = (char*)&val; //quick and dirty
  putchar(ret[7]);
  putchar(ret[6]);
  putchar(ret[5]);
  putchar(ret[4]);
  putchar(ret[3]);
  putchar(ret[2]);
  putchar(ret[1]);
  putchar(ret[0]);
  */
  

}

float getfloat(FILE* input)  //we have to reverse the byte order from JAVA (at least for linux/gcc)
{	
	union 
	{
		float flt;
		char byts[sizeof(float)/sizeof(char)];  	
	};
	for(size_t i = sizeof(float)/sizeof(char); i > 0u; --i)
	{		byts[i-1] = getc(input);}
	return flt;
}


/* getfloat is broken.  solve once we need it.
float getfloat()  //we have to reverse the byte order from JAVA (at least for linux/gcc)
{
  float ret;
  ((char*)&ret)[3]=getchar();
  ((char*)&ret)[2]=getchar();
  ((char*)&ret)[1]=getchar();
  ((char*)&ret)[0]=getchar();
  return ret;
}
*/

void sendfloat(float val)  //we have to reverse the byte order for JAVA (at least from linux/gcc)
{
	union 
	{
		float flt;
		char byts[sizeof(float)/sizeof(char)];  	
	};
	flt = val;
  
	for(size_t i = sizeof(float)/sizeof(char); i > 0u; --i)
	{	putchar(byts[i-1]);}
}



ssflt getWeights(FILE* input, ssflt *weight)
{
	double weightSum = 0.0;
	int modelnum(globalModelNum);
    globalCommand->assertCommand(input, "BeginBinary", modelnum * 8 );
    for(int q = 0; q < modelnum; ++q )
    {
		double mw = getdouble(jobfile?jobfile:stdin);
		//assert( (mw >= 0.0) && (mw <= 1.0) ); //scaled weights as of 17.09.2015
		weight[q] = (ssflt)mw;
		weightSum += (double)mw;
	}
    globalCommand->assertCommand(input, "EndBinary", weightSum );
    return weightSum;
}

void sendWeights(ssflt *weight)
{
    int modelnum(globalModelNum);
    for(int q = 0; q < modelnum; ++q )
    {
		printf("%lf\n", (double)weight[q]);
    }
    fflush(stdout);
}

int getModels(FILE* input, unsigned short* Models, int fm, int lm )
{
	const int modelarea(globalModelArea);
	globalCommand->assertCommand(input, "BeginBinary", 2 * (lm-fm+1) * modelarea );

	int checksum = 0;
    for (int m = fm; m <= lm; ++m)
    {
        for (int qr = 0; qr < modelarea; ++qr)
        {
            unsigned short val = 256 * getc(input);
            val += getc(input);
            checksum += val;
            Models[ m * modelarea + qr] = val;
            if(((int)val) >= globalModelval)
            {
				printf("Error: modelvalue %d excceeds modelval: %d\n", val, globalModelval);
				fflush(stdout);
				return -1;
			}
        }
    }
    globalCommand->assertCommand(input, "EndBinary", checksum );
    return checksum;
}

int getBeam(FILE* input, unsigned short* Beam)
{
	const int beamarea( globalBeamArea );
	globalCommand->assertCommand(input, "BeginBinary", 2 * beamarea );

	int checksum(0);
	unsigned short beam_min(0);
   
	for (int qr(0); qr < beamarea; ++qr)
	{
		unsigned short val (256 * getc(input));
		val += getc(input);
		checksum += val;
		Beam[qr] = val;
		if( (qr==0) || (val < beam_min) )
		{ beam_min = val;}
	}
    globalCommand->assertCommand(input, "EndBinary", checksum );
	
	for (int qr(0); qr < beamarea; ++qr)
	{
		Beam[qr] -= beam_min;
	}
	
	globalBeamNorm = checksum - beamarea*beam_min;	
    return checksum;
}




int sendModels(unsigned short* Models, int fm, int lm )
{
	int modelsize(globalModelSize);
	int modelarea = modelsize * modelsize;
	///HMMM maybe we should do some handshakes here. we will see	
    for (int m = fm; m <= lm; ++m)
    {
        printf("MODEL\t%d\tPIXELS%d\n", m, modelarea);
        fflush(stdout); //make sure master is aware of the awaing pixels
        for (int qr = 0; qr < modelarea; ++qr)
        {
            printf("%d\n", Models[ m * modelarea + qr]);
			fflush(stdout);
        }
    }
    fflush(stdout);//should be redundant, but safety first
    return 0;
}


int sendEModels(ssflt* EModelsB, ssflt* wghtB,
#ifdef USE_MPFR
				ssint* EModelsS, ssint* wghtS, 
#endif
				bool reporting,
				int fm, int lm)
{
	const int modelarea( globalModelArea );
	for (int m = fm; m <= lm; ++m)
    {
				
#ifdef USE_MPFR
		printf("EMODEL\t%d\tPIXELS\t%d\tWGHTB\t%lf\tWGHTS\t%d\n", m, modelarea, wghtB[m], wghtS[m]);
		
#else
		printf("EMODEL\t%d\tPIXELS\t%d\tWGHTB\t%E\n", m, modelarea, wghtB[m]);
#endif		
        fflush(stdout); //make sure master is aware of the awaing pixels
		if(!reporting)
		{		
			int ind( m * modelarea); 
			for (int qr = 0; qr < modelarea; ++qr)
			{
			
#ifdef USE_MPFR            
				printf("%lf\t%d\n", EModelsB[ind], EModelsS[ind]);	
#else
				printf("%E\n", EModelsB[ind]);	
#endif
				fflush(stdout);
				++ind;
			}
		}
		else
		{
			printf("hidding %d pixels\n",modelarea);
		}
	}
    fflush(stdout);//should be redundant, but safety first
	return 0;
}




#ifdef USE_CCTABLE

int getIdata(FILE* input, unsigned short *qImg)
{
	const int IdataSize(globalqImgSize);
	const int qImgCol(globalqImgCol);
	const int subframearea(globalSubFrameArea);
	const int subframesize(globalSubFrameSize);
	const int subframecount(globalSubFrameCount);
	const int datamin(globalDatamin);
	const int * const pixpos(pixposG);
	globalCommand->assertCommand(input, "BeginBinary", 3 * IdataSize/2 );
	int checksum = 0;
	int chsum = 0;
	pix_non_zero = 0;
	for (  int idx = 0; idx < IdataSize; ++idx)
	{
		int sfC = idx / subframearea;
		int ontheframe = idx % subframearea;
		
		

		
		if(pixpos[ontheframe] != -1)
		{
			unsigned short val = 256 * getc(input);
			val += getc(input);
			chsum += val;
			val -= datamin;
			const int x = ontheframe % subframesize;
			const int y = ontheframe / subframesize;
			const int ind = (sfC + x * subframecount + y * qImgCol);
			qImg[ind] = val;
			checksum += val;
			pix_non_zero += (val?1:0); 
			
			if(val > globalDatamax)
			{
				printf("Error datavalue %d exceeds limit: %d\n", val, globalDatamax);
				return -1;
			}
		}
		
	}
	globalCommand->assertCommand(input, "EndBinary", chsum );
	return checksum;
}

#else //NO_CCTABLE we only need the zipped image
int getIdata(FILE* input, unsigned short *zipImg, long &zipImgUsed)
{
	const int IdataSize(globalqImgSize);
	const int subframearea(globalSubFrameArea);
	//const int subframesize(globalSubFrameSize);
	const int subframecount(globalSubFrameCount);
	const int datamin(globalDatamin);
	const int datamax(globalDatamax);
	globalCommand->assertCommand(input, "BeginBinary", 3 * IdataSize / 2 );
	int checksum = 0;
	int chsum = 0;
	pix_non_zero = 0;
	unsigned short *zipPtr = zipImg;	
	unsigned short *sf = new unsigned short[subframearea];
	unsigned short *histo = new unsigned short[datamax];
	for(int sfC(0); sfC < subframecount; ++sfC)
	{
		memset(histo,'\0',sizeof(unsigned short)*datamax);
		memset(sf,'\0', sizeof(unsigned short)*subframearea);
		unsigned short sf_max_val(0);	
		for(int ind = 0; ind < subframearea; ++ind)
		{
			if(pixposG[ind]!=-1)
			{
				unsigned short val = 256 * getc(input);
				val += getc(input);
				chsum += val;
				val -= datamin;
				checksum += val;
				sf[ ind] = val;
				if(val)//skip the zeros
				{
					++pix_non_zero;
					++(histo[val-1]);
					sf_max_val = (val > sf_max_val) ? val : sf_max_val;
				}
			}
		}
		//assert((int)sf_max_val <= datamax - datamin);
		for(int vv(1); vv <=  sf_max_val; ++vv)
		{
			if(histo[vv-1])
			{
				*(zipPtr++) = vv; //value
				*(zipPtr++) = histo[vv-1]; //runlength
				//unsigned short chcksum(0);
				for(int ind = 0; ind < subframearea; ++ind)
				{
					if((pixposG[ind] != -1) && (sf[ind] == vv))
					{
						*(zipPtr++) = pixposG[ind]; //indices
						//++chcksum;
					}
				}
				//assert(chcksum == histo[vv-1]);
			}
		}
		*(zipPtr++) = 0; //value
		*(zipPtr++) = 0; //runlength
	}
	//pointer is finally one behind the actual data,-> correct length
	zipImgUsed = zipPtr - zipImg;
	
	globalCommand->assertCommand(input, "EndBinary", chsum );
	delete[] sf;
	delete[] histo;
	return checksum;

}
#endif

int getZipIdata(FILE* input, unsigned short *zipImg, long &zipImgUsed)
{
	int subframecount(globalSubFrameCount);
	int datamin(globalDatamin);
	globalCommand->assertCommand(input, "BeginBinary", 2 * zip_img_cap );
	int checksum = 0;
	int chsum = 0;
	long shorts_read = 0;
	pix_non_zero = 0;
	unsigned short* zipPtr = zipImg;
	for(int sfC(0); sfC < subframecount; ++sfC)
	{
		//printf("FRAME %d of %d\n",sfC, subframecount);
		//fflush(stdout);
		unsigned short count(0);
		unsigned short val0(0);
		unsigned short val(0);
		unsigned int sf_sum(0);
		do
		{
			val0 = 256 * getc(input);
			val0 += getc(input);
			chsum += val0;
			val = val0;
			if(val0 >= datamin)
			{	val = val0 - datamin;}
			//else
			//{	assert(val0 == 0);}
			//checksum += val;
			count = 256 * getc(input);
			count += getc(input);
			chsum += count;
			//printf("val: %d count: %d\n",val,count);
			//fflush(stdout);
			shorts_read += 2;
			
			if( (count > 0) && (val > 0) )
			{
				*(zipPtr++) = val;
				*(zipPtr++) = count;
				checksum += (val * count);
				sf_sum += (val * count);
			}
			
			for(unsigned short pos = 0; pos < count; ++pos)
			{
				unsigned short qr = 256 * getc(input);
				qr += getc(input);
				//printf("%d ",qr);
				//fflush(stdout);
				*(zipPtr++) = pixposG[qr];
				chsum += qr;
			}
			shorts_read += count;
			//printf("\n");
			//printf("remaining %ld\n", zip_img_cap - (zipPtr - zipImg));
			fflush(stdout);
		}
		while( (val > 0) || (count > 0) );
		*(zipPtr++) = 0; // 0 times 0 marks the end of a subframe
		*(zipPtr++) = 0;
		
		if(sf_sum == 0)
		{
			printf("getZipIdata: sf_sum == 0 at sfC: %d\n", sfC);
		}
	
	}
	zipImgUsed = zipPtr - zipImg;
	//printf("zipReading %ld of max. %ld\n", zipImgUsed, zip_img_cap);
	//fflush(stdout);
	//long padcount = 0;
	for(long padding = shorts_read; padding < zip_img_cap; ++padding)
	{
		getc(input);
		getc(input);
		//++padcount;
	}
	//printf("zipReading %ld extra shorts totalling %ld\n", padcount, padcount + zipImgUsed);
	//fflush(stdout);
	globalCommand->assertCommand(input, "EndBinary", chsum );
	return checksum;
}

#ifdef MUL_PTBL
ssflt getPtable(FILE* input, ssflt *lPtable, ssflt *mulPtbl, ssflt *Ptbl0)
#else
ssflt getPtable(FILE* input, ssflt *lPtable, ssflt *Ptable, ssflt* InvPtable, ssflt *Ptbl0)
#endif
{
	const ssflt InvLn2(1/log(2.0));
	double fchecksum(0.0);
	const int modelval(globalModelval);
	const int ptablewidth(globalPtableWidth);
    globalCommand->assertCommand(input, "BeginBinary", 8 * globalPtableSize );
    for (int datav = 0; datav < ptablewidth; ++datav)
    {
        for (int modelv = 0; modelv < modelval; ++modelv)
        {
            int ind = datav + ptablewidth * modelv;
            double pval = getdouble(jobfile?jobfile:stdin);
            //assert( isfinite(pval) );
            fchecksum += pval;
            pval *= pvalscaling;
            lPtable[ind] = (ssflt)pval*InvLn2;
#ifndef MUL_PTBL
            pval = exp(pval);
            Ptable[ind] = (ssflt)pval;
            InvPtable[ind] = 1.0/(ssflt)pval;
#endif
        }
    }
    globalCommand->assertCommand(input, "EndBinary", fchecksum );

#ifdef MUL_PTBL
	ssflt *mulPtblp( mulPtbl );
	//double min_log_jmp;
	//double max_log_jmp;
	//double sum(0.0); //should be zero in the end again!
	//bool first = true;
	for(int newv(0); newv < modelval; ++newv )
	{
		for(int oldv(0); oldv < modelval; ++ oldv)
		{
			for(int datav(0); datav < ptablewidth; ++datav)
			{
				//multiplying exps should be more accurate for floating point numbers than exping the difference
				*(mulPtblp++) = exp2( lPtable[datav + newv * ptablewidth]) * exp2( -lPtable[datav + oldv * ptablewidth] );
				//*mulPtblp = lPtable[datav + newv * ptablewidth] - lPtable[datav + oldv * ptablewidth];
				
				/*sum += *mulPtblp;
				if(first || *mulPtblp < min_log_jmp )
				{
					min_log_jmp = *mulPtblp;	
				}
				if(first || *mulPtblp > max_log_jmp )
				{
					max_log_jmp = *mulPtblp;
					first = false;
				}
				++mulPtblp;
				*/
				
				
			}
		}
	}
	//printf("log mulPtbl smallest: %lf biggest: %lf sum: %lf\n", min_log_jmp, max_log_jmp, sum);
#endif

	ssflt *Ptbl0P(Ptbl0);
	for(int datav(0); datav < ptablewidth; ++datav)
	{
		for(int newv(0); newv < modelval; ++newv )
		{
			*(Ptbl0P++) = exp2( lPtable[datav + newv * ptablewidth]) * exp2( -lPtable[ /*datav=0*/ newv * ptablewidth] );
			//*(Ptbl0P++) = lPtable[datav + newv * ptablewidth] - lPtable[ /*datav=0*/ newv * ptablewidth]; 
		}	
	}
    return fchecksum;
}


int sendHisto(ssflt* histoB, 
#ifdef USE_MPFR
				ssint* histoS,
#endif
				bool reporting
				)
{
	const int ptableSize ( globalPtableSize );			
	
	
	 //make sure master is aware of the awaing data
	if(!reporting)
	{					
#ifdef USE_MPFR
		printf("HISTO\tPOINTS\t%d\tBASE\tSHIFT\n", ptableSize);
	
#else
		printf("HISTO\tPOINTS\t%d\tBASE\n", ptableSize);
#endif	
		fflush(stdout);	
		for (int ind = 0; ind < ptableSize; ++ind)
		{
		
#ifdef USE_MPFR            
			printf("%lf\t%d\n", *(histoB++), *(histoS++));	
#else
			printf("%E\n", *(histoB++));	
#endif
			fflush(stdout);
		}
	}
	else
	{
		printf("\nhidding %d histogramm entries\n",ptableSize);
	}
    fflush(stdout);//should be redundant, but safety first
	return 0;
}

