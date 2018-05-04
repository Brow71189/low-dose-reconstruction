#include <stdio.h>
#include <assert.h>
#include "string.h"
#include "globals.h"
#include "zippedimg.hpp"


long write_zipImg(	unsigned short* zipImg,	unsigned short* qImg)				
{
	const int subFrameArea(globalSubFrameArea);
	unsigned short *sf = new unsigned short[subFrameArea];
	unsigned short *zipPtr = zipImg;	
	const int datamax(globalDatamax);
	unsigned short *histo = new unsigned short[datamax];
	const int subFrameCount(globalSubFrameCount);
	const int subFrameSize(globalSubFrameSize);
	const int qImgCol(globalqImgCol);
	
	for (int sfC(0); sfC < subFrameCount; ++sfC)
	{
		//clear the histogram
		memset(histo,'\0',sizeof(unsigned short) * datamax);
		memset(sf,'\0',sizeof(unsigned short) * subFrameArea);
		unsigned short* sfp(sf);
		unsigned short sf_max_val(0);
		//collect the distributed pixels of the current subframe and determine histogramm
		for(int ind(0); ind < subFrameArea; ++ind)
		{
			const int x(ind%subFrameSize);
			const int y(ind/subFrameSize);
			const unsigned short val = qImg[ sfC + x * subFrameCount + y * qImgCol ];
			*(sfp++) = val;
			if(val)//skip the zeros
			{
				++(histo[val-1]);
				sf_max_val = (val > sf_max_val) ? val : sf_max_val;
			}	
		}
		//assert(sf_max_val <= (unsigned short)globalDatamax );
		for(int vv(1); vv <=  sf_max_val; ++vv)
		{
			if(histo[vv-1])
			{
				//assert(zipPtr - zipImg < zip_img_cap - 2);
				*(zipPtr++) = vv; //value
				*(zipPtr++) = histo[vv-1]; //runlength
				unsigned short chcksum(0);
				for(int ind=0; ind<subFrameArea; ++ind )
				{
					if(pixposG[ind]!=-1)
					{
						if(sf[ind] == vv)
						{
							//assert(zipPtr - zipImg < zip_img_cap - 1);
							*(zipPtr++) = pixposG[ind]; //indices
							++chcksum;
						}
					}
				}
				//assert(chcksum == histo[vv-1]);
			}
		}
		//assert(zipPtr - zipImg < zip_img_cap - 2);
		*(zipPtr++) = 0; //value
		*(zipPtr++) = 0; //runlength
	}
	
	assert(zipPtr - zipImg < zip_img_cap);
	/*
	FILE * pFile;
	pFile = fopen ("zibImg.bin", "wb");
	fwrite (zipImg , sizeof(unsigned short), zipPtr - zipImg, pFile);
	fclose (pFile);
	*/
	delete[] sf;
	delete[] histo;
	return(zipPtr - zipImg);
}

void read_zipImg(	unsigned short* zipImg,	unsigned short* qImg)
{
	unsigned short* zipPtr = zipImg;
	const int subFrameCount(globalSubFrameCount);
	const int subFrameSize(globalSubFrameSize);
	const int qImgCol(globalqImgCol);

	for (int sfC(0); sfC < subFrameCount; ++sfC)
	{
		unsigned short val(0);
		unsigned short elements(0);
		unsigned int sf_sum(0);
		do
		{
			val = *(zipPtr++);
			elements = *(zipPtr++);
			sf_sum += val*elements;
			for(unsigned short el(0); el < elements; ++el)
			{	
				const int indqr = globalhexPixels[*(zipPtr++)];
				const int x = indqr % subFrameSize;
				const int y = indqr / subFrameSize;			
				qImg[ sfC + x * subFrameCount + y * qImgCol ] = val;		
			}
		}
		while( (val > 0) || (elements > 0) );
		if(sf_sum == 0)
		{
			printf("read_zipImg: sf_sum == 0 at sfC: %d\n", sfC);
			exit(-1);
		}
		
	}
}

unsigned short look_up_sf(int indpp)
{
	unsigned short res(0);
	unsigned short val(0);
	unsigned short elements(0);
	do
	{
		val = *(zip_img_sfP++);
		elements = *(zip_img_sfP++);
		for(unsigned short el(0); el < elements; ++el)
		{
			if( *(zip_img_sfP++) == indpp)
			{	res = val;}
		}
		
	} while ( (val>0) || (elements > 0) );
	return res;
}

