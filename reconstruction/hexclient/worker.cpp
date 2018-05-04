
#include "stdio.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "string.h"

#include "fast_mpfr.h"
#include "fast_cache.h"
#include "communicator.hpp"
#include "command.h"
#include "globals.h"
#include "likelyhood.hpp"

/***********************
 * PIVATE DECLARATIONS *
 ***********************/
/*
inline void
probability( double &base, int &probexp, const double *Ptable, const int &ptablewidth, const double pvaloffset,
			 const unsigned char *sfData, const int &subframesize, const int &datamin,
			 const unsigned short *Model, const int &modelsize,
			 const int x0, const int y0);
*/


/***************
 * DEFINITIONS *
 ***************/
/*
//Caution assumes log Ptable
inline void
probability( double &base, int &probexp, const double *Ptable, const int &ptablewidth, const double pvaloffset,
			 const unsigned char *sfData, const int &subframesize, const int &datamin,
			 const unsigned short *Model, const int &modelsize,
			 const int x0, const int y0)
{
	base = pvaloffset;
	int tmpexp;
    probexp = 0;
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
            base += Ptable[ (datav-datamin) + ptablewidth * modelv ];
        }

    }
	base = frexp( exp(base) , &tmpexp);
	probexp += tmpexp;
}
*/
int main()
{
    Command command("");
    FILE *log = 0; 
    globalCommand = &command;
    clock_t timer = clock();
    //Dirty Hack to ignore possible leading alignment bytes in a binary stream
    getchar();
    getchar();
    getchar();
    getchar();
    command.next(); //skip possibly truncated first Dummy() command 
    
    
	int threadnum      = 0; 		//default
    int impWidth       =-1;  		//width of every Frame
    int impHeight      =-1; 		//height of every Frame
    int modelsize      =-1; 		//width and height of every model
    int modelnum       =-1;  		//number of models
    int modelval       =-1;  		//scaling for gray values in model      
    int datamax        =-1;   		//brightes pixel in any Frame
    int datamin        =-1;   		//darkest pixel in any Frame
    int subframesize   =-1;
    int stepsize       =-1;		  //lattice of the graphene
    int subframeoffset = 0;	      //default
    int probdebug      =-1;      //default 
    bool reporting     = false;  //default
    {
		const char *head; 
		const char *args;
		do
		{
			command.next();
			head = command.getHead().c_str();
			args = command.getArgs().c_str();
			//TODO check readout
			if (strcmp(head, "ThreadNum") == 0)
			{
				sscanf( args, "%d", &threadnum);	
			}
			else if (strcmp(head, "ImpWidth") == 0 )
			{
				sscanf( args, "%d", &impWidth);	
			}
			else if (strcmp(head, "ImpHeight") == 0 )
			{
				sscanf( args, "%d", &impHeight);
			}
			else if (strcmp(head, "DataMin") == 0 )
			{
				sscanf( args, "%d", &datamin);
			}
			else if (strcmp(head, "DataMax") == 0 )
			{
				sscanf( args, "%d", &datamax);
			}
			else if (strcmp(head, "SubFrameSize") == 0 )
			{
				sscanf( args, "%d", &subframesize);
			}
			else if (strcmp(head, "ModelSize") == 0 )
			{
				sscanf( args, "%d", &modelsize);
				globalModelSize = modelsize;
			}
			else if (strcmp(head, "ModelNum") == 0 )
			{
				sscanf( args, "%d", &modelnum);
			}
			else if (strcmp(head, "ModelVal") == 0 )
			{
				sscanf( args, "%d", &modelval);
			}
			else if (strcmp(head, "StepSize") == 0 )
			{
				sscanf( args, "%d", &stepsize);
			}
			else if (strcmp(head, "SubFrameOffset") == 0 )
			{
				sscanf( args, "%d", &subframeoffset);
			}
			else if (strcmp(head, "ProbDebug") == 0 )
			{
				sscanf( args, "%d", &probdebug);
			}
			else if (strcmp(head, "Reporting") == 0 )
			{
				int val;
				sscanf( args, "%d", &val);
				reporting = (val == 1);
			}
			else if (strcmp( head, "Run") == 0) 
			{
				int val;
				sscanf( args, "%d", &val);
				if (val == 0)
				{
					return 0;
				}
			}
			else
			{	
				printf("unknown initialization command %s(%s)\n", head, args );
				return -1;
			}
		} while ( strcmp( head, "Run") ); 
	}
    if( (impWidth <= 0) ||
		(impHeight <= 0) ||
		(modelsize <= 0) ||
		(modelnum <= 0) ||
		(modelval <= 0) ||
		(datamin < 0) ||
		(datamax <= 0) ||
		(stepsize <= 0) ||
		(subframesize <= 0) )
	{
		printf("Missing or invalid input parameters.\n");
		return -1;
	}
    
    printf("Worker %i online\n", threadnum);
    if(reporting)
    {
        log = fopen("log.txt","w");
        printf("Frame Width = %d\tFrame Height = %d\n", impWidth, impHeight);
        printf("Data min = %d\t\tData max = %d\n", datamin, datamax);
        printf("Number of Models = %d\tLattice = %d\n", modelnum, stepsize);
        printf("Model Size = %d\t\tModel max = %d\n", modelsize, modelval);
        printf("subframesize = %d\tsubframeoffset = %d\n", subframesize, subframeoffset);
    }

    fflush(stdout);

    if ( (modelsize % stepsize) || 
		 (subframesize % stepsize) || 
		 (subframeoffset % stepsize) )
	{
		printf("Lattice mismatch for worker: %d\n", threadnum);
		if(log) fclose(log);
		return -1;
	}
    
    double modelweight[modelnum], old_modelweight[modelnum];
    getWeights(modelweight, modelnum);
    for(int q = 0; q < modelnum; ++q )
    {
		old_modelweight[q] = modelweight[q]; 
	}
    
    
    //allocate and read array of model data
       
    const int subframearea = subframesize * subframesize;
    const int modelarea = modelsize * modelsize;
    globalModelArea = modelarea;
    const int modelDataCount = modelnum * modelarea;
	
    unsigned short* Models = new unsigned short[modelDataCount];
    globalModels = Models;
    getModels(modelnum);
    
    int ptablewidth = 1+datamax-datamin;
    globalPtableWidth = ptablewidth;
    int ptablesize = ptablewidth*modelval;
    //Ptable is read as logvalues
    //good for initalization, will be exped before mainloop
    globalDatamax = datamax;
    globalDatamin = datamin;
    globalModelval = modelval;
    
    double* lPtable = new double[ptablesize];
    double* Ptable = new double[2*ptablesize];
    double *InvPtable =  Ptable + ptablesize;   
    getPtable(lPtable, ptablesize);
    
    //read our piece of data.
    
    int framestart = -1;
    int framestop = -1; 
    command.getCommand( "FrameStart", framestart);
    command.getCommand( "FrameStop", framestop);
    
    
    int framecount =1+framestop-framestart;

    
    const int xsfoCount 		= (impWidth-subframeoffset) / subframesize;
    const int ysfoCount 		= (impHeight-subframeoffset) / subframesize;
    const int sfoCount			= xsfoCount * ysfoCount;
    const int framearea        = subframearea * sfoCount;
    const int xoCount 			= subframesize / stepsize;
    const int yoCount 			= subframesize / stepsize;
    const int oCount           = xoCount * yoCount;
    const int moCount          = oCount * modelnum;
    const int subframeCount 	= framecount * sfoCount;
    const int qImgCount 	= subframeCount * subframearea;

    unsigned char *qImg = new unsigned char[qImgCount];

    const int impSize = impWidth * impHeight;
    const int IdataSize = framecount * impSize;
    if(reporting)
    {
        printf("IdataSize,subframeCount     %d,%d\n", IdataSize, subframeCount);
        printf("xsfoCount,ysfoCount         %d,%d\n", xsfoCount, ysfoCount);
        printf("xoCount,yoCount             %d,%d\n", xoCount, yoCount);
    }
    
    {
		command.assertCommand("BeginBinary", 1 * IdataSize );
		int checksum = 0;
		int chsum = 0;
		for (  int idx = 0; idx < IdataSize; ++idx)
		{
			int frame = idx / impSize;
			int ontheframe = idx % impSize;
			int framey = (ontheframe / impWidth) - subframeoffset;
			int framex = (ontheframe % impWidth) - subframeoffset;
			int ysfoC = framey / subframesize;
			int xsfoC = framex / subframesize;
			int y = framey % subframesize;
			int x = framex % subframesize;

			short val = 256 * getchar();
			val += getchar();
			
			chsum += val;
		   
			//Simply skip the pixels outside of our subframes
			if( (xsfoC >= 0) && (xsfoC < xsfoCount) && 
				(ysfoC >= 0) && (ysfoC < ysfoCount) &&
				(x >= 0) && (y >= 0))
			{
				qImg[ QIDIND(frame,xsfoC,ysfoC,x,y) ] = val;
				checksum += val;
			}
		}
		command.assertCommand("EndBinary", chsum );
		printf("%i Frames (%i to %i) checksum: %i\n",framecount,framestart+1,framestop+1,chsum);
		if(reporting)
		{
			printf("%d Pixels of ImageData imported\n", qImgCount);
		}
	}
    const int ccSize = subframeCount * moCount;
    ccflt  *cctablebase = new ccflt[ccSize];
    ccint  *cctableexp  = new ccint[ccSize];
    ccflt *ncctablebase = new ccflt[ccSize];
    ccint *ncctableexp  = new ccint[ccSize];
    ccflt *tmptablebase;
    ccint *tmptableexp;
	
    if(reporting)
    {
        printf("%d cctable elements\n", 4*ccSize);
    }
	//Base and Shift of subsums
    const int ssmCount = subframeCount * modelnum;
    double  *ssmB = new double[ssmCount];
    int     *ssmS = new int[ssmCount];
    double  *nssmB = new double[ssmCount];
    int    *nssmS = new int[ssmCount];
    double *tmpssmB;
    int *tmpssmS;
    
    if(reporting)
    {
        printf("%d ssm elements\n", 4*ssmCount);
        clock_t now = clock();
        double elapsed = double(now-timer)/CLOCKS_PER_SEC;
        printf("%fs for loading\n", elapsed);
        timer = now;
    }

    //calculate initial values B..base and S..shift aka exp2
    {
		//double ELPB = 0;
		//int ELPS = 0;
		double totalLPB = 0;
		init_likelyhood( totalLPB, modelnum, subframeCount, 
		                 oCount, xoCount, stepsize, ssmB, ssmS,
		                 nssmB, nssmS, cctablebase, cctableexp,
		                 ncctablebase, ncctableexp, lPtable, 
		                 qImg, subframesize,
		                 datamin, Models, modelsize, modelweight );
		/*
		for (int sfC = 0; sfC < subframeCount; ++sfC)
		{
			ELPB = 0;
			ELPS = 0;
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
					probability( probval, probexp, Ptable, ptablewidth, qImg + (sfC * subframearea), //0.0 is for log Pval
								 subframesize, datamin, Models + (m * modelarea), modelsize,
								 xo, yo + yoadj );				
					const unsigned int index = QQCCIND(sfC,m,oC);
					cctablebase[ index ] = probval;
					cctableexp [ index ] = probexp;
					ncctablebase[ index ] = probval;
					ncctableexp [ index ] = probexp;
					
					if (reporting && sfC==0 && probdebug==m)
					{
						printf("xoC,yoC %d,%d\tprob. %f << %d\n", xoC, yoC, probval, probexp);
					}
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
		*/ 
		//fflush(stdout);
		//Now Make regular Ptable probabilities for updating 
		{
			double * pval = Ptable;
			double * lpval = lPtable;
			for (int i = 0; i < ptablesize;  ++i)
			{
				*pval = exp(*lpval);
				*(pval + ptablesize) = 1 / *pval;
				++pval;
				++lpval;
			}
		}

		//return the first calculated value

		qfr_out_str(totalLPB);
		fflush(stdout);
		if(reporting)
		{
			fprintf(log, "Pval = %lf\n", totalLPB);
			fflush(log);
		}
		
	}

    int nccok[modelnum];
    int chthis[modelnum]; //redundant with nccok !=0 but +1 s¹ performance increas
    
    //double weightB[modelnum];
    //int weightS[modelnum];

    for (int z=0; z<modelnum; ++z)
    {
        nccok[z]=1;
    }

    //nccok: 0 - has changed;  1 - has not changed once, 2- has not changed twice

    //////////////////////////////////////////
    //here comes the main loop
    ///////////////////////////////
   
    if(reporting)
    {
        printf("%d ssm values calculated\n",ssmCount);
        clock_t now = clock();
        double elapsed = double(now-timer)/CLOCKS_PER_SEC;
        printf("%fs for initalization\n", elapsed);
        timer = now;
    }
    {
		int maxchlistlen = 16;
		int chlistlen = 0;
		int chx[maxchlistlen];
		int chy[maxchlistlen];
		int chm[maxchlistlen];
		int newval[maxchlistlen];
		int oldval[maxchlistlen];		
		int total_runs = 0;
		do
		{
			bool quit_now = false;
			bool update_weights = false;
			//bool model_guess = false;
			//bool guess_weights = false;
			bool sp_adjust = false;
			
			for (int m = 0; m < modelnum; ++m)
			{
				chthis[m] = 0;
				++nccok[m];
			}
						
			{
				const char * head;
				const char * args;
				chlistlen = 0;
				do
				{
					command.next();
					head = command.getHead().c_str();
					args = command.getArgs().c_str();
					if(reporting && strcmp(head,"Run") )
					{
						fprintf(log, "%s[%s]\n", head, args);
					}
					/*
					if (strcmp(head, "PixelGuess") == 0)
					{
						if( (chlistlen == 0) && !update_weights && !model_guess && !sp_adjust)
						{
							unsigned short m,x,y;
							sscanf( args, "%hu,%hu,%hu", &m, &x, &y);
							chx[chlistlen] = x;
							chy[chlistlen] = y;
							chm[chlistlen] = m;
							oldval[chlistlen] = Models[ MODIDX(m,x,y) ];
							++chthis[m];
							nccok[m] = 1; //actually not needed but faster updates
							++chlistlen;
							model_guess = true; //FIXME true
						}
						else
						{
							printf("Error: %s[%s] must be alone.\n", head, args);
						}
					}
					*/
					if (strcmp(head, "ModelPixel") == 0 /*&& !model_guess*/ && !update_weights)
					{
						if(chlistlen < maxchlistlen)
						{
							unsigned short m,x,y,nv;
							sscanf( args, "%hu,%hu,%hu,%hu", &m, &x, &y, &nv);
							chx[chlistlen] = x;
							chy[chlistlen] = y;
							chm[chlistlen] = m;
							newval[chlistlen] = nv;
							oldval[chlistlen] = Models[ MODIDX(m,x,y) ];
							++chthis[m];
							nccok[m] = 0;
							++chlistlen;
							sp_adjust = true;
						}
					}
					else if (strcmp(head,"NewWeights") == 0 /*&& !model_guess*/ && !sp_adjust && !update_weights)
					{
						getWeights(modelweight, modelnum);
						update_weights = true;
					}
					/*
					else if (strcmp(head, "GuessWeights") == 0 && !model_guess && !sp_adjust && !guess_weights)
					{
						for(int m = 0; m < modelnum; ++m)
						{
							weightB[m] = 0.0;
							weightS[m] = 0;
						}
						guess_weights = true;
						
					}
					*/ 
					else if (strcmp( head, "Run") == 0) 
					{
						int val;
						sscanf( args, "%d", &val);
						quit_now = (val == 0);
					}
					else if (strcmp( head, "UpdateModels") == 0)
					{
						getModels(modelnum);
						double totalLPB(0.0);
						init_likelyhood( totalLPB, modelnum, subframeCount, 
							oCount, xoCount, stepsize, ssmB, ssmS,
							nssmB, nssmS, cctablebase, cctableexp,
							ncctablebase, ncctableexp, lPtable, 
							qImg, subframesize,
							datamin, Models, modelsize, modelweight );
						for (int m=0; m<modelnum; ++m)
						{
							nccok[m]=2;
						}
						
						if(reporting)
						{
							fprintf(log, "Pval = %lf\n", totalLPB);
							fflush(log);
						}	
					}
					else
					{	
						printf("unknown/invalid command %s[%s]\n", head, args );	
					}
				} while ( strcmp( head, "Run") ); 
				fflush(log);
			}
			if (quit_now)
			{
				printf("\nWorker %d received Run(0) command.\n", threadnum);
				if (reporting)
				{
					fprintf(log,"Run[0]\n");
				}
				break;
			}
			
			double ELPB(0.0);
			int ELPS(0);
			double totalLPB(0.0);
			//The following are for averaging the Imagedata -> the actual expectation value of ModelPixel
			//double dataB(0), probB(0);
			//int dataS(0), probS(0);
			
			//loop for calculating modified log-prob value
			for( int sfC(0); sfC < subframeCount; ++sfC)
			{
				double newsubsumB;
				int newsubsumS;
				ELPB = 0.0;
				ELPS = 0;
				for (int m(0); m < modelnum; ++m) 
				{   
					int ssmid( QSSMIND(sfC,m) );
					double subsumB( ssmB[ ssmid ] );
					int subsumS( ssmS[ ssmid ] );
					
					/*if(model_guess && chthis[m])
					{
						qfr_add( probB, probS, subsumB, subsumS);
					}*/					
					
					if (nccok[m] < 2) //recent change
					{
						subsumB = 0.0;
						subsumS = 0;
						for (int oC(0); oC < oCount; ++oC)
						{									
							int yoC( oC / xoCount);
							int xoC( oC % xoCount);
							int index( QQCCIND(sfC,m,oC) );
							double probvalupd(cctablebase[ index ]);
							int probexpupd(cctableexp [ index ]);
							int yos(yoC * stepsize);
							int xos(xoC * stepsize);
							int yoadj( ( xoC & 1 ) ? ( stepsize / 2 ) : 0 );
							if(chthis[m]) //current change
							{
								for (int chnum(0); chnum < chlistlen; ++chnum)
								{		
									
									if( chm[chnum] != m)
									{
										continue;
									}
									
									int y( yos + chy[chnum] + yoadj);
									y -= ( y >= subframesize )?subframesize:0;
																	
									int x( xos + chx[chnum]);
									x -= ( x >= subframesize )?subframesize:0;
									
									int datav( qImg[ QQIND(sfC,x,y) ] - datamin );
									
									/*
									if(model_guess)
									{
										qfr_add( dataB, dataS, probvalupd * datav, probexpupd);	
									}
									else*/
									//{
										probvalupd *= ( ( Ptable[ datav + ptablewidth * newval[chnum] ] * 
												   InvPtable[ datav + ptablewidth * oldval[chnum] ] ) );				 		 			
									//}
								}
								int tmpexp(0);
								probvalupd  = frexp(probvalupd, &tmpexp );
								probexpupd  += (tmpexp);
								
							}
							ncctablebase[ index ] = probvalupd;
							ncctableexp [ index ] = probexpupd;
							qfr_add(subsumB, subsumS, probvalupd, probexpupd);	
						}						
						nssmB[ssmid] = subsumB;
						nssmS[ssmid] = subsumS;
					}
					/*if(guess_weights)
					{	qfr_add(weightB[m], weightS[m], subsumB, subsumS);}*/
					qfr_mul_d(newsubsumB, newsubsumS, subsumB, subsumS, modelweight[m]);
					qfr_add(ELPB, ELPS, newsubsumB, newsubsumS);
					
				} //endfor m
				//calculate log-prob and add to overall sum
				qfr_log(ELPB, ELPS);
				totalLPB += ELPB;
			}
			
			/*if(guess_weights)
			{
				double normB = 0.0;
				int normS = 0;
				for(int m = 0; m < modelnum; ++m)
				{
					qfr_add(normB, normS, weightB[m], weightS[m] );
				}				
				for(int m = 0; m < modelnum; ++m)
				{
					modelweight[m] = ldexp(weightB[m]/normB, weightS[m] - normS);
				}
			}*/
			
			/*if(model_guess)
			{
				dataB = ldexp( dataB / probB, dataS - probS);
				dataS = 0;
				totalLPB = (dataB * (double)modelval/( datamax - datamin) );
				if (reporting)
				{
					fprintf(log, "m,x,y,oldval,fit\t %d\t%d\t%d\t%d\t%lf\n", 
					chm[0], chx[0], chy[0], oldval[0], totalLPB);
				} 
			}*/
			
			if (reporting && ((total_runs % 20) != 19))
			{
				printf(".");
				fflush(stdout);
			}
			else
			{
				//senddouble(totalLPB);
				qfr_out_str(totalLPB);
				//if(guess_weights) {	sendWeights(modelweight, modelnum);}
				 
			}
			fflush(stdout);
			if(reporting)
			{
				fprintf(log, "%lf", totalLPB);
				//if(guess_weights) {	fprintf(log,"\tnewWeights");}
				fflush(log);	
			}
			
			
			int keep;
			if( !command.getCommand("Keep", keep) )
			{
				break;
			}
			if(reporting)
			{
				fprintf(log,"\t%s[%d]\n", command.getHead().c_str(), keep);
				fflush(log);
			}
			
			if ((chlistlen > 0) /*&& !model_guess*/)
			{
				if ((keep==1))
				{
					//exchange current and new cctables
					
					tmptablebase	= cctablebase;
					tmptableexp		= cctableexp;
					cctablebase		= ncctablebase;
					cctableexp		= ncctableexp;
					ncctablebase	= tmptablebase;
					ncctableexp		= tmptableexp;
					
					//exchange current and new subsums
					
					tmpssmB          = ssmB;
					tmpssmS          = ssmS;
					ssmB             = nssmB;
					ssmS             = nssmS;
					nssmB            = tmpssmB;
					nssmS            = tmpssmS;
					
					//apply changes to models
					for (int chnum=0; chnum<chlistlen; ++chnum)
					{
						Models[ MODIDX(chm[chnum],chx[chnum],chy[chnum]) ] = newval[chnum];
					}
				}
			}
			
			// In principle correlated weight and model changes are possible
			if (update_weights)
			{
				if (keep==1) //update the backup copy
				{
					for (int m = 0; m < modelnum; ++m)
					{
						old_modelweight[m] = modelweight[m];
					}
				}
				else //revert to backup copy
				{
					for (int m = 0; m < modelnum; ++m)
					{
						modelweight[m] = old_modelweight[m];
					}
				}
			}
			
			total_runs += (chlistlen + update_weights);
		} while(true);
		
		if(reporting)
		{
			printf("\n%d iterations finished in ", total_runs);
			clock_t now = clock();
			double elapsed = double(now-timer)/CLOCKS_PER_SEC;
			printf("%fs\t (Test/second = %f)\n", elapsed, total_runs/elapsed);
			timer = now;
		}
	}
     
    delete[] cctablebase;
    delete[] cctableexp;
    delete[] ncctablebase;
    delete[] ncctableexp;
    delete[] ssmB;
    delete[] ssmS;
    delete[] qImg;
    delete[] Ptable;
    delete[] Models;
    fclose(log);
    return 0;
}


