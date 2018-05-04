#ifdef USE_MPI
#include "mpi.h"
#endif
#include "stdio.h"
#include <stdlib.h>
#include <sched.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>

#include <time.h>
#include <assert.h>
#include "string.h"
#include <sstream>
#include <fstream>

#include "globals.h"
#include "fast_mpfr.h"
#include "fast_cache.h"
#include "communicator.hpp"
#include "command.h"
#include "hexlikelyhood.hpp"
#include "hex_em.hpp"
#include "transformations.hpp"
#include "zippedimg.hpp"
#include "histogram.hpp"
#include "abberation.hpp"


void sigint_handler(int )
{	
	//printf("SIGINT -> pushing_pack Run(0)\n");
	globalCommand->push_back("Run(0)");
}


int cpu_binding(int threadnum) {
	cpu_set_t mask;
	sched_getaffinity(0, sizeof(mask), &mask);
	unsigned int cpus = 0;
	while (CPU_ISSET(cpus, &mask)) {++cpus;} 
	//we know that there is always al least one cpu, but scan-build doesnt
	unsigned int rank = threadnum % ( (cpus>0)?cpus : 1 );
	unsigned int cpu = 0;
	cpu = (rank&1) ? ( cpus - 1 - (rank >> 1) ) : (rank >> 1);
	CPU_ZERO(&mask);
	CPU_SET(cpu, &mask);
	sched_setaffinity(0, sizeof(mask), &mask);
	return cpu;
}


void* save_realloc(void* ptr, size_t size, int lnr = -1)
{
	void* const nptr( realloc(ptr, size) );
	if(nptr != nullptr) 
	{	return nptr; }//everything is ok
	else
	{
		printf("ERROR: reallocation failed in %s at LINE: %d\n", __FILE__, lnr);
		fflush(stdout);
		free(ptr); //free memory of callie, its ptr is now dangling
		return nullptr;
	}
}

int main(int argc, char** argv)
{
	static Command command("init(1)\n");
	globalCommand = &command;
	signal(SIGINT, sigint_handler);
	int pid = (int)getpid(); 
#ifdef USE_MPI	
	int mpi_rank, mpi_size;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif
#ifdef FANCY_FEATURES
	bool fancyf(true);
#else	
	bool fancyf(false);
#endif
	if(argc > 1)
	{
		if( strcmp(argv[1],"-v")==0 )
		{
			printf("This is hexclient compiled at %s %s\n", __DATE__, __TIME__ );
			printf("by %s %s\n", 
#ifdef __clang__
		"CLANG"		
#else
#ifdef __ICC
		"ICC"
#else
		"GNUC"
#endif		
#endif						
			,__VERSION__ );
			printf("This software has not been released under any license\n");
			printf("You are not supposed to run it or trust the guy who gave it to you\n");
			printf("The source code is maintained by Christian.Kramberger-Kaplan@univie.ac.at\n");
			
#ifdef USE_CCTABLE			
			printf("ENABLED  ");
#else			
			printf("DISABLED ");
#endif
			printf("cachetable\n");

#ifdef USE_MPFR			
			printf("ENABLED  (ccflt:%lu  ccint:%lu) ", sizeof(ccflt), sizeof(ccint) );
#ifndef MPFR_TRUE
			printf(" (FAKE)");
#endif
#else			
			printf("DISABLED ");
#endif					
			printf("MPFR numbers\n");

#ifdef USE_MPFR
#ifdef MACRO_MPFR
			printf("MPFRSTYLE MACRO\n");
#else
			printf("MPFRSTYLE INLINE\n");
#endif
#endif					
			
#ifdef MUL_PTBL			
			printf("ENABLED  ");
#else			
			printf("DISABLED ");
#endif
			printf("3D Ptable\n");
			
			
#ifdef USE_MAX			
			printf("ENABLED  ");
#else			
			printf("DISABLED ");
#endif
			printf("picking best fit among symetric orbit of models\n");
			
#ifdef USE_MPI			
			printf("ENABLED  MPI with rank/size %d/%d\n", mpi_rank, mpi_size);
			MPI_Finalize();
#else
			printf("DISABLED MPI\n");		
#endif			
			printf("%s FANCY_FEATURES (polling Histogramms and EM)\n", fancyf ? "ENABLED" : "DISABLED");
			
			
	
			
			
			printf("MAXCHLISTLEN %d\n", MAXCHLISTLEN);
			
			
			return 0;
		}
		
		jobfile = fopen(argv[1], "r");
		bool success = ( jobfile != 0);
		if (!success)
		{
			printf("could not open %s\n", argv[1]);
			return -1;
		}
		//else {	printf("reading: %s\n", argv[1] );}
#ifdef USE_MPI	
		if(mpi_size > 1)
		{
			std::stringstream ss;
			ss << "out" << mpi_rank << "of" << mpi_size << ".txt";
			success = ( freopen(ss.str().c_str(), "w", stdout) != 0);
			if( !success)
			{
				printf("cannot write to %s\n", ss.str().c_str() );
				return -1;
			}
		}
#endif		
    }
#ifdef USE_MPI   
    else if(mpi_size > 1)
    {
		if(mpi_rank == 0)
		{
			printf("overcrowded stdin: %d processes tried to read from same instream\n", mpi_size);
		}
		return -1;
	}
#endif
    clock_t timer = clock();
    
    FILE *logfile = 0;
    std::ofstream reportfile;
    
	long zipImgUsed     = -1;
								    //TODO and also transmit zipped data
    long zipImgLen      =-1;
    int debug_level     =-1;         //default, only effective if reporting
	int threadnum       = 0; 		//default
	int modelsize       =-1; 		//width and height of every model
	int beamsize        =-1;
	int beamradius      =-1;
    int modelnum        =-1;  		//number of models
    int shadownum		= 1;		//every model is its own first shadow
    int modelval        =-1;  		//scaling for gray values in model
    int datamax         =-1;   		//brightes pixel in any Frame
    int datamin         =-1;   		//darkest pixel in any Frame
    int subframesize    =-1;
    int stepsize        =-1;		    //lattice of the graphene
    int wobble          = 0;			//include points up to that distance around lattice sites
    int framecount      =-1;
    bool use_beam       = false;
    bool reporting      = false;  //default
    bool use_init_black = false; //carefull default
    bool lattice_hopping = true; 
    {
		const char *head;
		const char *args;
		do
		{
			//printf("polling command ...\t");
			command.next(jobfile?jobfile:stdin);
			head = command.getHead().c_str();
			args = command.getArgs().c_str();
			//printf("%s[%s]\n",head,args);
			if (strcmp(head, "ThreadNum") == 0)
			{
				sscanf( args, "%d", &threadnum);
			}
			else if (strcmp(head, "FrameCount") == 0 )
			{
				sscanf( args, "%d", &framecount);
			}
			else if (strcmp(head, "ZipImgLen") == 0 )
			{
				sscanf( args, "%ld", &zipImgLen);
			}
			else if (strcmp(head, "DataMin") == 0 )
			{
				sscanf( args, "%d", &datamin);
			}
			else if (strcmp(head, "DataMax") == 0 )
			{
				sscanf( args, "%d", &datamax);
			}
			else if (strcmp(head, "ModelSize") == 0 )
			{
				//in principle there could be smaller models thean subframes
				//but that should be totally equivalent to same size overlapping subframes
				sscanf( args, "%d", &modelsize);
				globalModelSize = modelsize;
				subframesize = modelsize;
				globalSubFrameSize = subframesize;
			}
			else if (strcmp(head, "BeamSize") == 0 )
			{
				sscanf( args, "%d", &beamsize);
				globalBeamSize = beamsize;
				globalBeamRadius = beamsize/2-1;
			}
			else if (strcmp(head, "BeamRadius") == 0 )
			{
				sscanf( args, "%d", &beamradius);
				globalBeamRadius = beamradius;
			}
			else if (strcmp(head, "ModelNum") == 0 )
			{
				sscanf( args, "%d", &modelnum);
			}
			else if (strcmp(head, "ShadowNum") == 0 )
			{
				sscanf( args, "%d", &shadownum);
			}
			else if (strcmp(head, "ModelVal") == 0 )
			{
				sscanf( args, "%d", &modelval);
			}
			else if (strcmp(head, "PvalOffSet") == 0 )
			{
				pvaloffset_ok = (sscanf( args, "%lf", &pvaloffset) == 1);
			}
			else if (strcmp(head, "pvalscaling") == 0 )
			{
				sscanf( args, "%lf", &pvalscaling);
				inv_pvalscaling = 1.0/pvalscaling;
			}
			else if (strcmp(head, "BondLength") == 0 )
			{
				sscanf( args, "%d", &stepsize);
			}
			else if (strcmp(head, "Wobble") == 0 )
			{
				sscanf( args, "%d", &wobble);
			}
			else if (strcmp(head, "Translations") == 0 )
			{
				int val;
				sscanf( args, "%d", &val);
				lattice_hopping = (val == 1);
			}
			else if (strcmp(head, "Reporting") == 0 )
			{
				int val;
				sscanf( args, "%d", &val);
				reporting = (val == 1);
			}
			else if (strcmp(head, "use_weights") == 0 )
			{
				int val;
				sscanf( args, "%d", &val);
				use_weights = (val == 1);
			}
			else if (strcmp(head, "use_profile") == 0 )
			{
				int val;
				sscanf( args, "%d", &val);
				use_beam = (val == 1);
			}
			else if (strcmp(head, "use_max") == 0 )
			{
			
				int val;
				sscanf( args, "%d", &val);
#ifdef USE_MAX					
				use_max = (val == 1);
#else				
				if(val==1)
					{	printf("ERROR: compile with -DUSE_MAX to enable %s(%s)\n",head,args);}	
#endif				
			}
			else if (strcmp(head, "black_init") == 0 )
			{
				int val;
				sscanf( args, "%d", &val);
				use_init_black = (val == 1);
			}
			else if (strcmp(head, "real_log2") == 0 )
			{
				int val;
				sscanf( args, "%d", &val);
				real_log2 = (val == 1);
			}
			else if (strcmp(head, "debug_level") == 0 )
			{
				sscanf( args, "%d", &debug_level);// <0 or >2 will prompt the user
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
			else if (strcmp(head,"DummyForByteAlignment") == 0)
			{
				//nothing to do, simply ignore this one
			}
			else
			{
				printf("unknown initialization command %s(%s)\n", head, args );
				return -1;
			}
		} while ( strcmp( head, "Run") );
	}
    if( (beamradius == -1) || (beamradius >  beamsize/2 - 1) )
    {
		beamradius = beamsize/2 - 1; //will stay -1 if beamsize was not set
		globalBeamRadius = beamradius;
	}
    //THis might break many many things down the road
    if(stepsize == 0) lattice_hopping = false;
    
    if(
		(modelsize <= 0) ||
		(modelnum <= 0) ||
		(shadownum < 1) ||
		(modelval <= 0) ||
		(datamax <= datamin) ||
		(lattice_hopping && (stepsize <= 0) ) ||
		(framecount < 1) ||
		(use_beam && (beamsize <= 4) )
		)
	{
		printf("Missing or invalid input parameters.\n");
		return -1;
	}
	if(use_beam) fancyf = false; //fancy features cannot be combined with the beamprofile
    if(reporting)
    {
        printf("Data min = %d\t\tData max = %d\n", datamin, datamax);
        printf("frames = %d\t\tLattice = %d%s\n", framecount, stepsize,lattice_hopping?"":" OFF");
        printf("zipImg Size = %lu kB\n", sizeof(unsigned short) * zipImgLen/1024);
        printf("Model Size = %d\t\tNumber of Models = %d x %d \n", modelsize, modelnum, shadownum);
        printf("Model min = 0\t\tModel max = %d\n", modelval);
        printf("Beam Size = %d\t\tWobble = %d\n", beamsize, wobble);
    }

    fflush(stdout);

    if ( lattice_hopping && (modelsize % stepsize) )
	{
		printf("Lattice mismatch for worker: %d\n", threadnum);
		return -1;
	}

	int cpu = -1;
	if(!reporting)
	{
		cpu = cpu_binding(threadnum);
	}
    //cache for equivalent lattice sites
    int wc(0);
    for(int w(1); w <= wobble; ++w )
    {	wc += w;  }
	wc = 1 + 6*wc;
    int *hexPixels(nullptr);
    int *latticePoints(nullptr);
	int hpMax( (3 * subframesize * subframesize) / 4 );
	int lpMax = lattice_hopping?( wc * (subframesize * subframesize) / ( 4 * stepsize * stepsize ) ):wc;
    if(lpMax>hpMax)
    {	lpMax = hpMax;}//There could never be more translations than there are pixels
    
    globallpMax = lpMax;
	globalhpMax = hpMax;
    int * lpX = new int[lpMax];
	int * lpZ = new int[lpMax];
    globallpX = lpX;
	globallpZ = lpZ;

    globalModelNum = modelnum;
    globalShadowNum = shadownum;
    const int modelsp1 = modelnum + 1;
    ssflt *modelweight = new ssflt[modelsp1]; 
    ssflt *old_modelweight = new ssflt[modelsp1];
    int *cases = new int[modelnum];
    const int subframearea = subframesize * subframesize;
    globalSubFrameArea = subframearea;
    const int modelarea = modelsize * modelsize;
    globalModelArea = modelarea;
    const int modelDataCount = modelnum * modelarea;
    const int smodelDataCount = 12 * shadownum * modelnum * modelarea;
    const int beamarea = beamsize * beamsize;
    globalBeamArea = beamarea;
    const int sympersf = lpMax * 12;
    const int transElements = fancyf?(hpMax * sympersf):0;
    const int moveElements = hpMax * lpMax;
    const int mirrotElements = hpMax*12;
    globalsympersf = sympersf;
    transElementsG = transElements;
    moveElementsG = moveElements;
    mirrotElementsG = mirrotElements;

    int* pixpos = new int[modelarea];
    pixposG = pixpos;
    memset(pixpos, 0xFF, modelarea * sizeof(int));//should actually become (int)-1
	
    unsigned short* transformations = fancyf ?  (new unsigned short[transElements]) : nullptr ;
	transformationsG = transformations;
	unsigned short* inv_transformations = fancyf ? (new unsigned short[transElements]) : nullptr ;
	inv_transformationsG = inv_transformations;
	unsigned short* translations = new unsigned short[moveElements];
	translationsG = translations;
	unsigned short* inv_translations = new unsigned short[moveElements];
	inv_translationsG = inv_translations;
	unsigned short* mirrot = new unsigned short[mirrotElements];
	mirrotG = mirrot;
	unsigned short* inv_mirrot = new unsigned short[mirrotElements];
    inv_mirrotG = inv_mirrot;
    
    {
		//prefetch the equivalent lattice sites inside the hex subframe
		const int rH = subframesize/2;
		if (reporting)
		{
			printf("Pixels/LatticePoints per frame %d/%d\n", hpMax, lpMax );
		}
		hexPixels = new int[hpMax];
		globalhexPixels = hexPixels;
		latticePoints = new int[lpMax];
		globallatticePoints = latticePoints;

		int hpC(0);
		int lpC(0);
		//odd-r -> cube
		const int mx = rH - (rH) / 2;
		const int mz = rH;
		globalMx = mx;
		globalMz = mz;

		for(int dx = -rH; dx < rH ; ++dx)
		{
			int dz1 = -rH;
			int dz2 = rH - dx;
			if (dx < 0)
			{
				dz1 = -rH - dx;
				dz2 = rH;
			}
			for(int dz = dz1; dz < dz2; ++dz)
			{
				int dy = -dx - dz;
				//translate to center
				int x_data = dx + mx;
				int z_data = dz + mz;
				// cube -> odd-r
				int q_data = x_data + (z_data) / 2;
				int r_data = z_data;
				int ind = q_data + r_data * subframesize;
				pixpos[ind] = hpC;
				hexPixels[hpC++] = ind;

				int x = dx;
				int y = dy;
				int z = dz;
				if(lattice_hopping)
				{	hexagonal_trap(x,y,z,stepsize);}
				
				if(  abs(x)+abs(y)+abs(z) <= 2*wobble ) //point is equivalent to center
				{
					lpX[lpC] = dx;
					lpZ[lpC] = dz;
					if((dx == 0) && (dz == 0))
					{	global_lp0 = lpC;}
					latticePoints[lpC++] = ind;
					
				}
			}
		}
	}
	
	//OUCH prefecthing tranformations before receiving last data causes 
	//inefficiency if wobble > 0

    //allocate the expected amount of memory

    unsigned short* const Models = new unsigned short[modelDataCount];
    memset(Models,0,modelDataCount*sizeof(unsigned short)); 
    unsigned short* const sModels = new unsigned short[smodelDataCount];
    memset(sModels,0,smodelDataCount*sizeof(unsigned short)); 
    unsigned short* const Beam = new unsigned short[beamarea];
    memset(Beam,0,beamarea*sizeof(unsigned short)); 
	//add this to make the output the actual log2
#ifdef USE_MAX	
	double PvalR = framecount * (  (real_log2 ? ( ( modelarea * 0.75 * pvaloffset)  / log(2.0) ) : 0 ) - ( use_max?0:log2(sympersf*shadownum) ) );
#else
	double PvalR = framecount * (  (real_log2 ? ( ( modelarea * 0.75 * pvaloffset)  / log(2.0) ) : 0 ) - log2(sympersf*shadownum) );
#endif    
    int ptablewidth = 1+datamax-datamin;
    globalPtableWidth = ptablewidth;
    int ptablesize = ptablewidth*modelval;
    globalPtableSize = ptablesize;
    globalDatamax = datamax;
    globalDatamin = datamin;
    globalModelval = modelval;
    ssflt *lPtable = new ssflt[ptablesize];
    ssflt *Ptbl0 = new ssflt[ptablesize];

#ifdef MUL_PTBL
    ssflt *mulPtbl = new ssflt[modelval * modelval * ptablewidth];
#else
	ssflt *PtableL = new ssflt[ptablesize];
    ssflt *Ptable = PtableL;
    ssflt *InvPtableL = new ssflt[ptablesize];
    ssflt *InvPtable = InvPtableL;
#endif

    const int subframeCount 	= framecount;
    globalSubFrameCount = subframeCount;
    const int qImgCount 	   = subframeCount * subframearea;
    globalqImgSize = qImgCount;

    const int ccperm = sympersf * subframeCount;
    globalccperm = ccperm;
    const int qImgCol = subframesize * subframeCount;
    globalqImgCol = qImgCol;

    /*global constants for ccindex*/
    sfG = subframeCount;
	sfmG = 2 * sfG;
	sfmrG = 6 * sfmG * shadownum; ///WARNING
	sfmrlG = lpMax * sfmrG;
	/*global constants for inverse ordered ccindex*/
	//sfG
	sfGlpM = sfG * lpMax;
	sfGlpMrot = sfG * lpMax * 6 * shadownum; ///WARNING


#ifdef USE_CCTABLE
	int cc_elements = sfmrlG * modelsp1;
#else
	int cc_elements = 0;
	use_init_black = true; //black values are always required without cctable 
#endif

	int ss_elements = sfG * modelsp1;
	int ptbl_elements = ptablesize;
#ifdef MUL_PTBL
	ptbl_elements += modelval * 2 * ptablesize;
#else
	ptbl_elements += 3 * ptablesize;
#endif
	unsigned long cc_width = sizeof(ccflt);
	unsigned long ss_width = sizeof(ssflt);
#ifdef USE_MPFR
	cc_width += sizeof(ccint);
	ss_width += sizeof(ssint);
#endif
	zip_img_cap = (zipImgLen > 0) ? zipImgLen : subframeCount * ( 2 * (datamax + 1) + 3*subframearea/4 + 2);
	unsigned long mem_usage = (cc_elements * cc_width + ss_elements * ss_width +
			ptbl_elements * sizeof(ssflt) +
			qImgCount * sizeof(unsigned short) + (modelDataCount+smodelDataCount+beamarea) * sizeof(unsigned short)
			+ (2*transElements + 2*moveElements + 2*mirrotElements + zip_img_cap) * sizeof(unsigned short)
			);

#ifdef USE_CCTABLE
    unsigned short *qImg = new unsigned short[qImgCount];
#endif
	unsigned short *zipImg = (unsigned short *) malloc(zip_img_cap * sizeof(unsigned short));

#ifdef USE_CCTABLE
    //cachetable and bufferline for likelyhoods per subframe, model, translation,rotation and mirror
    ccflt ** const cctableB = new ccflt*[modelsp1];
    //ccflt ** const cctableL = new ccflt*[modelsp1];
    //ccflt ** const cctableK = new ccflt*[modelsp1];
#endif
    ssflt  ** const ssmB = new ssflt*[modelsp1];
#ifdef USE_MPFR
#ifdef USE_CCTABLE
	ccint ** const cctableS = new ccint*[modelsp1];
#endif
	ssint  ** const ssmS = new ssint*[modelsp1];
#endif
#ifdef USE_CCTABLE
    for(int m = 0; m < modelsp1; ++m )
    {
		cctableB[m] = new ccflt[sfmrlG];
		//std::fill_n(cctableB[m],(ccflt)(-42),sfmrlG);
		//cctableL[m] = new ccflt[sfmrlG];
		//cctableK[m] = new ccflt[sfmrlG];
#ifdef USE_MPFR
		cctableS[m] = new ccint[sfmrlG];
#endif
	}
	ccflt *cctableBM = cctableB[0];
	//ccflt *cctableLM = cctableL[0];
	//ccflt *cctableKM = cctableK[0];
    ccflt *ncctableB = cctableB[modelnum];
    //ccflt *ncctableL = cctableL[modelnum];
    //ccflt *ncctableK = cctableK[modelnum];
#endif

	for(int m = 0; m < modelsp1; ++m )
    {
		ssmB[m] = new ssflt[subframeCount];
#ifdef USE_MPFR
		ssmS[m] = new ssint[subframeCount];
#endif
	}

    ssflt *ssmBM = 0;
    ssflt *nssmB = ssmB[modelnum];

#ifdef USE_MPFR
#ifdef USE_CCTABLE
    ccint *cctableSM = cctableS[0];
    ccint *ncctableS = cctableS[modelnum];
#endif
    ssint *ssmSM = 0;
    ssint *nssmS = ssmS[modelnum];
#endif    
    ssflt totalLPB(0.0);

#ifndef USE_CCTABLE
	ssflt *nblackB = new ssflt[modelnum];
	//ssflt *nblackL = new ssflt[modelnum];
	//ssflt *nblackK = new ssflt[modelnum];
#ifdef USE_MPFR
	ssint *nblackS = new ssint[modelnum];
#endif
#endif

	ssflt *blackB = new ssflt[modelnum];
	//ssflt *blackL = new ssflt[modelnum];
	//ssflt *blackK = new ssflt[modelnum];
	for(int m(0); m < modelnum; ++m)
	{	
		blackB[m] = ssflt(0.0);
		//blackL[m] = ssflt(0.0);
		//blackK[m] = ssflt(0.0);	
	}
#ifdef USE_MPFR
	ssint *blackS = new ssint[modelnum];
#endif
    if(reporting)
    {
        printf("log2 Pval correction %lf\n", PvalR);
        printf("Bytes for cc_flt:%lu, ss/sfc_flt:%lu\n", sizeof(ccflt), sizeof(ssflt));
#ifdef USE_MPFR
		printf("Bytes for cc_int:%lu, ss/sfc_int:%lu\n", sizeof(ccint), sizeof(ssint));
#endif
#ifdef USE_CCTABLE
        printf("Idata            %d\n", qImgCount);
#endif
        printf("zipIdata       < %ld\n", zip_img_cap);

        printf("Transformations  %d\n", 2*(transElements + moveElements + mirrotElements));
        printf("Ptable           %d\n", ptbl_elements);
        printf("cc               %d\n", cc_elements);
        printf("ssm              %d\n", ss_elements);
        printf("sModels, Models  %d, %d\n", smodelDataCount, modelDataCount);
        printf("Beam             %d, %s\n", beamarea, (use_beam?"active":"unused") );
        printf("totaling       < %lu MB\n", 1+(mem_usage >> 20) );
	}
	printf("Worker%i: ", threadnum);
#ifdef USE_MPFR

#ifdef MPFR_TRUE
	printf("MPFR");
#else
	printf("MPFR (fake)");
#endif

#else
	printf("FLP");
#endif
#ifdef USE_MAX	
	printf(use_max?" MAX(%d)":" SUM(%d)", sympersf*shadownum);
#else
	printf(" SUM(%d)", sympersf*shadownum);
#endif
#ifdef USE_CCTABLE
	printf(" with cc_table");
#endif	
	if(use_init_black)
	{	printf( " (black init)");}

    printf("  ( < %lu MB) with pid: %d on cpu: %d weights: %s\n", (unsigned long) 1+(mem_usage >> 20), pid, cpu, (use_weights?"true":"false"));
    fflush(stdout);

	



    //Now receive the binary data for Model,Ptable Weights and Idata
    //Checks that every item is EXACTLY SEND ONCE
    {
		const char *head;
		const char *args;
		int IdataSum(0);
		int modelSum(0);
		int beamSum(0);
		double pTableSum(0.0);
		double weightSum(0.0);
		bool new_Models = false;
		bool new_Ptable = false;
		bool new_Weights = false;
		bool new_Idata = false;
		bool new_Beam = false;
#ifdef USE_CCTABLE
		bool new_zipIdata = false;
#endif
		do
		{
			command.next(jobfile?jobfile:stdin);
			head = command.getHead().c_str();
			args = command.getArgs().c_str();
			
			if( strcmp( head, "Run") == 0)
			{break;} //check for run BEFORE head might become a dangling pointer
			if ( strcmp(head, "updateModels") == 0 && !new_Models )
			{
				new_Models = true;
				modelSum = getModels(jobfile?jobfile:stdin, Models);
				if(reporting)
				{
					printf("Models checksum: %i\n", modelSum);
				}
			}
			else if ( strcmp(head, "updateBeamProfile") == 0 && !new_Beam )
			{
				new_Beam = true;
				beamSum = getBeam(jobfile?jobfile:stdin, Beam);
				if(reporting)
				{
					printf("Beam checksum: %i\n", beamSum);
				}
			}
			else if (strcmp(head, "updatePtable") == 0 && !new_Ptable)
			{
				new_Ptable = true;
#ifdef MUL_PTBL
				pTableSum = getPtable(jobfile?jobfile:stdin, lPtable, mulPtbl, Ptbl0);
#else
				pTableSum = getPtable(jobfile?jobfile:stdin, lPtable, Ptable, InvPtable, Ptbl0);
#endif

				if (reporting)
				{
					printf("Ptable checksum: %lf\n", pTableSum);
				}
			}
			else if (strcmp(head, "updateWeights") == 0 && !new_Weights)
			{
				new_Weights = true;
				weightSum = getWeights(jobfile?jobfile:stdin, modelweight);
				modelweight[modelnum] = 0.0;
				memcpy(old_modelweight, modelweight, sizeof(ssflt) * modelsp1 );
				if (reporting)
				{
					printf("Weights checksum %f\n", weightSum);
				}
			}

			else if (strcmp(head, "updateIdata") == 0 && !new_Idata)
			{
				new_Idata = true;

				//getIdata subtracts globalDatamin from all values
#ifdef USE_CCTABLE
				IdataSum = getIdata(jobfile?jobfile:stdin, qImg);
				zipImgUsed = write_zipImg(	zipImg, qImg);
#else //NO_CCTABLE
				IdataSum = getIdata(jobfile?jobfile:stdin, zipImg, zipImgUsed);
#endif
				zipImg = (unsigned short *)save_realloc(zipImg, zipImgUsed * sizeof(unsigned short), __LINE__);
				assert(zipImg);
				zip_img_cap = zipImgUsed;
				if(reporting)
				{
					printf("zipped Idata : %ld\t< %lu MB\n", zipImgUsed, 1+(zipImgUsed>>19) );
				}
				if(reporting)
				{
					printf("Idata checksum: %i\n", IdataSum);
				}
				datamax-= datamin;
				datamin = 0;
				globalDatamax = datamax;
				globalDatamin = datamin;
				//This is the place where we would switch from jobfile to stdin
				if(jobfile)
				{
					fclose(jobfile);
					jobfile = 0;
				}
			}

			else if (strcmp(head, "updateZipIdata") == 0 && !new_Idata)
			{
				new_Idata = true;
#ifdef USE_CCTABLE
				new_zipIdata = true;
#endif
				//geZiptIdata subtracts globalDatamin from all values

				IdataSum = getZipIdata(jobfile?jobfile:stdin, zipImg, zipImgUsed);
				zipImg = (unsigned short *)save_realloc(zipImg, zipImgUsed * sizeof(unsigned short), __LINE__);
				assert(zipImg);
				zip_img_cap = zipImgUsed;

				if(reporting)
				{
					printf("zipped Idata : %ld\t< %lu MB\n", zipImgUsed, 1+(zipImgUsed>>19));
					printf("Idata checksum: %i\n", IdataSum);
				}
				datamax-= datamin;
				datamin = 0;
				globalDatamax = datamax;
				globalDatamin = datamin;
				//This is the place where we would switch from jobfile to stdin
				if(jobfile)
				{
					fclose(jobfile);
					jobfile = 0;
				}
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
				printf("unknown/duplicate update request %s(%s)\n", head, args );
				return -1;
			}

		}while ( true );
		
		
		
#ifdef USE_CCTABLE
		if(new_zipIdata)
		{	read_zipImg( zipImg, qImg);}
#endif
		//report with one line to master, will be forwarded to console
		if( !( new_Models && new_Ptable && new_Weights && new_Idata && ( (!use_beam) || new_Beam))  )
		{
			printf("Missing input? Checksums: Idata %d, Models %d, pTable %f, Weights %f\n", IdataSum, modelSum, pTableSum, weightSum);
			if(reporting)
			{
				if(!new_Models) printf("no Models received\n");
				if(!new_Ptable) printf("no Ptable received\n");
				if(!new_Weights) printf("no Weights received\n");
				if(!new_Idata) printf("no ImageData received\n");
				if(use_beam && (!new_Beam) ) printf("no Beam Profile received\n");
			}
		}
		else //everythig should be all right
		{
			printf("Worker%i: Input validated\n", threadnum);
		}
		fflush(stdout);
	}
	
	//prefetch forward and backward tranformations only after receiving the input
	//master will feed the next worker meanwhile
	{
		unsigned short *transP = transformations;
		unsigned short *inv_transP = inv_transformations;
		unsigned short *moveP = translations;
		unsigned short *inv_moveP = inv_translations;
		unsigned short *mirrotP = mirrot;
		unsigned short *inv_mirrotP = inv_mirrot;
		
		///TODO check if the full transformations are still needed anywhere
		//split the loops for different memory chunks
		if(fancyf) //only fancy features require full transformations, they are incompatible with use of beam profile, though
		{
			for(int hp(0); hp < hpMax; ++hp)//all hexpixels
			{
				//Full transformations aligned to inverse
				for(int mirror(0); mirror < 2; ++mirror)//two mirrors
				{
					for(int rot(0); rot < 6; ++rot) //six rotations
					{
						for(int lp(0); lp < lpMax; ++lp)//all latticepoints
						{
							assert(transP - transformations == TRANSFORM(hp,lp,rot,mirror) );
							*(transP++) = (unsigned short)transform_hp( hp, lp, rot, mirror );
							*(inv_transP++) = (unsigned short)inv_transform_hp( hp, lp, rot, mirror );
						}
					}
				}
			}
		}
		for(int hp(0); hp < hpMax; ++hp)//all hexpixels
		{
			//Translations aligned to inv
			for(int lp(0); lp < lpMax; ++lp)//all latticepoints
			{
				assert(moveP - translations == MOVE(hp,lp));
				*(moveP++) = (unsigned short)transform_hp( hp, lp, 0, 0 );
				*(inv_moveP++) = (unsigned short)inv_transform_hp( hp, lp, 0, 0 );
			}
		}
		
		for(int hp(0); hp < hpMax; ++hp)//all hexpixels
		{
			//PointGroup aligned to inv
			for(int mirror(0); mirror < 2; ++mirror) //six inv rotations
			{
				for(int rot(0); rot < 6; ++rot)//two inv mirrors
				{
					assert(mirrotP - mirrot == MIRROT(hp,rot,mirror) );
					*(mirrotP++) = (unsigned short)transform_hp( hp, global_lp0, rot, mirror );
					*(inv_mirrotP++) = (unsigned short)inv_transform_hp( hp, global_lp0, rot, mirror );
				}
				
			}
		}
		
	
		
	}
	
	init_smodels(sModels,Models,Beam,use_beam);
	
	if(reporting && ( (debug_level < 0 ) || (debug_level > 2)))
	{
        FILE *tty = fopen("/dev/tty", "r");
        clock_t now = clock();
        double elapsed = double(now-timer)/CLOCKS_PER_SEC;

        printf("%fs for loading\n", elapsed);
        if (tty)
        {
			printf("Please choose a debug_level (0 ... 2): ");
			fflush(stdout);
			if(fscanf(tty,"%i",&debug_level)!=1 || (debug_level < 0) || (debug_level > 2))
			{
				debug_level = 0; //speed shall be the default
			}
			printf("\rdebug_level is %i: ", debug_level);
			switch(debug_level)
			{
				case 0:
					printf("speed mode\n");
					break;
				case 1:
					printf("full output\n");
					break;
				case 2:
					printf("full output & log.txt\n");
					break;
				default:
					printf("default in debug_level _> 1 (full_output)\n");
					debug_level = 1;
			}
			fclose(tty);
		}
		timer = clock();
    }

#ifdef USE_CCTABLE
	if(!use_init_black ) //|| true at this place fixes init_black issue => blavk cctable should be wrong PVAL and report seem fine though
	{
 
    //This variant is about half the performance, but may become
    //necessary at higher doses when black probabilities are no longer feasible
	totalLPB = init_hexlikelyhood(	ssmB,
									cctableB,
									//cctableL,
									//cctableK,
#ifdef USE_MPFR
									ssmS,
									cctableS,
#endif
									lPtable,
									zipImg,
									sModels, modelweight );


	
	}
	else //init with black values
	{
 
		init_black( blackB,
				    //blackL,
				    //blackK,
#ifdef USE_MPFR
					blackS,
#endif
					sModels,lPtable );

		init_cache_from_black(
								blackB,
								//blackL,
								//blackK,
								cctableB,
								//cctableL,
								//cctableK,
								ssmB,
#ifdef USE_MPFR
								blackS,
								cctableS,
								ssmS,

#endif
								Ptbl0,
								zipImg,
								sModels );
/*		
		printf("black cctables done\n");
		totalLPB = init_hexlikelyhood(	ssmB,
									cctableB,
									//cctableL,
									//cctableK,
#ifdef USE_MPFR
									ssmS,
									cctableS,
#endif
									lPtable,
									zipImg,
									sModels, modelweight );
		printf("asserted cctables done\n");
*/		

	
		
		
		
		
		totalLPB = init_LPB_from_subsums(	ssmB,
#ifdef USE_MPFR
											ssmS,
#endif
											modelweight, modelnum);
	}



#else //NO CCTABLE
	
	init_black( blackB,
				//blackL,
				//blackK,
#ifdef USE_MPFR
				blackS,
#endif
				sModels,lPtable );
								
	
	init_subsums_from_black(	blackB,
								//blackL,
								//blackK,
								ssmB,
#ifdef USE_MPFR
								blackS,
								ssmS,
#endif
								Ptbl0,
								zipImg, sModels);


	totalLPB = init_LPB_from_subsums(	ssmB,
#ifdef USE_MPFR
										ssmS,
#endif
										modelweight, modelnum);
#endif

	totalLPB += PvalR;
	//return the first calculated value and the timing
	{
		clock_t now = clock();
		double elapsed = double(now-timer)/CLOCKS_PER_SEC;
		printf("%lf\t%lf\n", totalLPB, elapsed);
		fflush(stdout);
	}
	if(reporting)
	{
		if(debug_level > 0)
		{
			std::stringstream ss;
			ss << "cachereport" << threadnum << ".txt";
			reportfile.open(ss.str().c_str());
		}
		cache_report(
				reportfile,
				blackB,
				//blackL,
				//blackK,
				ssmB,
#ifdef USE_CCTABLE
				cctableB,
				//cctableL,
				//cctableK,
#endif
#ifdef USE_MPFR
				blackS,
				ssmS,
#ifdef USE_CCTABLE
				cctableS,
#endif
#endif
				cases, reporting);


		if(debug_level > 0)
		{
			
			logfile = fopen("log.txt","wt");
			fprintf(logfile, "Pval = %lf\n", totalLPB);
			fflush(logfile);
		}
        clock_t now = clock();
        double elapsed = double(now-timer)/CLOCKS_PER_SEC;
        printf("%fs for initialization\n", elapsed);
        timer = clock();
    }

    ///////////////////////////////
    //here comes the main loop
    ///////////////////////////////
    const int maxchlistlen(MAXCHLISTLEN);
    {
		int num_changes(0);
		int chlistlen(0);
		int chx[maxchlistlen];
		int chy[maxchlistlen];
		int newval[maxchlistlen];
		int oldval[maxchlistlen];
		int total_runs = 0;
		int report_runs = 0;

		//local handles for frequent global varibales in CCINDEX
		const int sfG = ::sfG;
		bool quit_now = false;
		int fm(0), lm(modelnum-1);
		do
		{
			///Removed mutual exclusion of pixel and weight adjustment 17.09.2015
			///for testing symmetry weighting
			bool update_weights(false);
			bool sp_adjust(false);
			bool sp_em(false);
			bool model_em(false);
			bool any_em(false);
			bool collect_histo(false);
			
			bool reinitialize(false);
			int active_m (modelnum);
			{
				const char * head;
				const char * args;
				chlistlen = 0;
				num_changes = 0;
				do
				{
					
					//printf("fetching command ...\t");
					//fflush(stdout);
					command.next(jobfile?jobfile:stdin);
					head = command.getHead().c_str();
					args = command.getArgs().c_str();
					//printf("%s[%s]\n", head, args);
					if(strcmp( head, "Run") == 0)
					{
						int val;
						sscanf( args, "%d", &val);
						quit_now = (val != 1);
						break;
					}
					if(reporting && (debug_level > 1) && strcmp(head,"Run") )
					{
						fprintf(logfile, "%s[%s]\n", head, args);
					}

					if ( (strcmp(head, "ModelPixel") == 0) && !any_em && !reinitialize )
					{
						if(chlistlen < maxchlistlen)
						{
							unsigned short m,x,y,nv;
							sscanf( args, "%hu,%hu,%hu,%hu", &m, &x, &y, &nv);
							if( (active_m == modelnum) || (active_m == m) ) //we are not checking for a valid m here
							{
								active_m = m;
								chx[chlistlen] = x;
								chy[chlistlen] = y;
								newval[chlistlen] = nv;
								oldval[chlistlen] = Models[ x + y * modelsize + m * modelarea ];
								++chlistlen;
								sp_adjust = true;
								//printf("ModelPixel(%hu,%hu,%hu,%hu)\n",m,x,y,nv);
							}
							else
							{
								printf("Invalid request for ModelPixel: all Pixels must belong to the same Model!");
							}

						}
						else
						{
							printf("Invalid request for ModelPixel\n");
						}
					}
					else if (strcmp(head, "Pixel_EM") == 0)
					{
						if(fancyf && (chlistlen < maxchlistlen) && !sp_adjust && !update_weights && !reinitialize )
						{
							unsigned short m,x,y;
							sscanf( args, "%hu,%hu,%hu", &m, &x, &y);
							if( (active_m == modelnum) || (active_m == m) )
							{
								active_m = m;
								chx[chlistlen] = x;
								chy[chlistlen] = y;
								newval[chlistlen] = -1;
								oldval[chlistlen] = Models[ x + y * modelsize + m * modelarea ];
								++chlistlen;
								sp_em = true;
							}
							else
							{
								printf("Invalid request for Pixel_EM\n");
							}
						}
						else
						{
							printf("Pixel_EM is not possible with current settingd\n");
						}
					}
					else if (strcmp(head, "Model_EM") == 0)
					{
						if( fancyf && (chlistlen == 0) && !update_weights && !any_em && !reinitialize  )
						{
							if(args[0]!='\0') //either no arguments or fm,lm are valid syntax
							{
								sscanf( args, "%d,%d", &fm, &lm);
							}
							model_em = true;
						}
						else
						{	printf("Invalid or unsupported request for Model_EM[%d,%d]\n",fm,lm);}
					}
					else if (strcmp(head, "Total_EM") == 0)
					{
						if( fancyf && (chlistlen == 0) && !update_weights && !any_em && !reinitialize )
						{
							fm = 0;
							lm = modelnum-1;
							
							model_em = true;
						}
						else
						{	printf("Invalid or unsupported request for Total_EM\n");}
							
					}
					
					else if (strcmp( head, "PollCases" ) == 0)
					{
						if(!reporting)
						{
							update_cases(
											ssmB,
#ifdef USE_MPFR
											ssmS,
#endif
											cases);
							for(int m=0; m < modelnum; ++m)
							{
								printf("%d\n", cases[m]);
								fflush(stdout);
							}
						}
					}
					else if ( (strcmp(head,"updateWeights") == 0) && !any_em )
					{
						getWeights(jobfile?jobfile:stdin, modelweight);
						modelweight[modelnum] = 0.0;
						update_weights = true;
					}
					else if(strcmp(head,"use_weights") == 0)
					{
						int val;
						sscanf( args, "%d", &val);
						use_weights = (val == 1);
					}
					else if (strcmp(head, "use_max") == 0 )
					{
						
						int val;
						sscanf( args, "%d", &val);
#ifdef USE_MAX						
						use_max = (val == 1);
						reinitialize = true; //Rare case but maybe implement something like soft reinit ()
						PvalR = framecount * ( modelarea * 0.75 * pvaloffset / log(2.0) - (use_max?0:log2(sympersf*shadownum)));
#else
						if(val==1)
						{	printf("ERROR: compile with -DUSE_MAX to enable %s(%s)\n",head,args);}
#endif					
					}
					else if (strcmp(head, "black_init") == 0 )
					{
						int val;
						sscanf( args, "%d", &val);
						use_init_black = (val == 1);
					}
					
					else if (strcmp( head, "CacheReport") == 0 )
					{
						int wrapnum;
						sscanf( args, "%d", &wrapnum);
						time_t timer;
						time (&timer);
						reportfile << "wrapnum " << wrapnum << " time: " << ctime(&timer) << "\n";
						cache_report(reportfile,
									  blackB,
									  //blackL,
									  //blackK,
									  ssmB,
#ifdef USE_CCTABLE
									  cctableB,
									  //cctableL,
									  //cctableK,
#endif
#ifdef USE_MPFR
									  blackS,
									  ssmS,
#ifdef USE_CCTABLE
									  cctableS,
#endif
#endif
									  cases, reporting);
					}
					else if (strncmp( head, "update", 6) == 0)
					{
						if (strcmp( head, "updateModels") == 0)
						{
							
							if(args[0]!='\0') //either no arguments or fm,lm are valid syntax
							{
								sscanf( args, "%d,%d", &fm, &lm);
							}
							getModels(jobfile?jobfile:stdin, Models, fm, lm);
							reinitialize = true;
							init_smodels(sModels,Models,Beam,use_beam);
						}
						else if (strcmp( head, "updateIdata") == 0)
						{
#ifdef USE_CCTABLE
							getIdata(jobfile?jobfile:stdin, qImg);
#else
							zip_img_cap = (zipImgLen > 0) ? zipImgLen : subframeCount * ( 2 * (datamax + 1) + 3*subframearea/4 + 2) ;
							zipImg = (unsigned short*)save_realloc(zipImg, zip_img_cap * sizeof(unsigned short), __LINE__);
							assert(zipImg);
							getIdata(jobfile?jobfile:stdin, zipImg, zipImgUsed);
							zipImg = (unsigned short*)save_realloc(zipImg, zipImgUsed * sizeof(unsigned short), __LINE__);
							assert(zipImg);
							zip_img_cap = zipImgUsed;
#endif
							reinitialize = true;
						}
						else if (strcmp( head, "updateZipIdata") == 0)
						{
							zip_img_cap = (zipImgLen > 0) ? zipImgLen : subframeCount * ( 2 * (datamax + 1) + 3*subframearea/4 + 2) ;
							zipImg = (unsigned short*)save_realloc(zipImg, zip_img_cap * sizeof(unsigned short), __LINE__);
							assert(zipImg);
							getZipIdata(jobfile?jobfile:stdin, zipImg, zipImgUsed);
							zipImg = (unsigned short *)save_realloc(zipImg, zipImgUsed * sizeof(unsigned short), __LINE__);
							assert(zipImg);
							zip_img_cap = zipImgUsed;
#ifdef USE_CCTABLE
							read_zipImg( zipImg, qImg);
#endif
						}
						else if (strcmp( head, "updatePtable") == 0)
						{
#ifdef MUL_PTBL
							getPtable(jobfile?jobfile:stdin, lPtable, mulPtbl, Ptbl0);
#else
							getPtable(jobfile?jobfile:stdin, lPtable, Ptable, InvPtable, Ptbl0);
#endif
							reinitialize = true;
						}
						else if (strcmp( head, "updateHisto") == 0)
						{
							//collect_histo must be stand alone
#ifdef USE_CCTABLE						
							collect_histo = fancyf && (!(sp_em || model_em || sp_adjust || reinitialize || update_weights ));
							if(!collect_histo)
							{
								printf("ERROR: updateHisto must be a stand alone task, and cannot be used with a beam profile\n");
							}
#else
							{	printf("ERROR: compile with -DUSE_CCTABLE to enable %s(%s)\n",head,args);}
							collect_histo = false;
#endif						
						
						}
					}
					else
					{
						printf("unknown/invalid command %s[%s]\n", head, args );
						fflush(stdout);
						return -1;
					}
					any_em = (sp_em || model_em );
					fflush(stdout); //be sure errors are noticed
				} while ( true ); // OUCH strcmp( head, "Run") head may be dangling
				if(debug_level > 0)
				{	fflush(logfile);}
			}
			if (quit_now)
			{
				printf("Worker %d received Run(0) command.\n", threadnum);
				if (reporting && (debug_level > 1))
				{
					fprintf(logfile,"Run[0]\n");
				}
				break;
			}
			//double check if master tried to trick us by issuing another command AFTER "collectHisto()"
			if(collect_histo && (sp_em || model_em || sp_adjust || reinitialize || update_weights) )
			{
				collect_histo = false;
				printf("ERROR: updateHisto must be a stand alone task\n");
				fflush(stdout);
			}
#ifdef USE_CCTABLE		
			if (collect_histo /*|| (reporting && (debug_level>1) && (report_runs%50==49) )*/ )
			{
				
				create_Histo(
							cctableB,
#ifdef USE_MPFR
							cctableS,
#endif
							qImg,
							Models
							);
				report_Histo(reporting);
			
				clear_Histo();
				//if(collect_histo) continue;
				continue;
				//In principle we could combine collect_histo with other tasks, 
				//yet it would always reflect an outdated state
			}
#endif				
			/*
			bool force_keep = ( false && ((report_runs % 5000) == 24) );
			if(force_keep)
			{	
				printf(" ");
				fflush(stdout);
			}
			reinit_for_fun = ( false && ((report_runs % 5000) == 25) );
			*/
			if(reinitialize  ) // ||  reinit_for_fun
			{
				/*
				if(reinit_for_fun)
				{
					printf("PVAL: %lf\treinitalizing for fun...\n", totalLPB);
					fm = 0;
					lm = modelnum-1;
					init_smodels(sModels,Models,Beam,use_beam);
				}
				*/
				//printf("reinitialzing...\n");
#ifdef USE_CCTABLE
				
					if(!use_init_black)
					{
				 
					//This variant is about half the performance, but may become
					//necessary at higher doses when black probabilities are no longer feasible
					totalLPB = init_hexlikelyhood(	ssmB,
													cctableB,
													//cctableL,
													//cctableK,
#ifdef USE_MPFR
													ssmS,
													cctableS,
#endif
													lPtable,
													zipImg,
													sModels, modelweight,
													fm, lm );



					}
					else //init with black values
					{
				 
						init_black( blackB,
									//blackL,
									//blackK,
#ifdef USE_MPFR
									blackS,
#endif
									sModels,lPtable,
									fm, lm );

						init_cache_from_black(
									blackB,
									//blackL,
									//blackK,
									cctableB,
									//cctableL,
									//cctableK,
									ssmB,
#ifdef USE_MPFR
									blackS,
									cctableS,
									ssmS,

#endif
									Ptbl0,
									zipImg,
									sModels,
									fm, lm );

						totalLPB = init_LPB_from_subsums(	ssmB,
#ifdef USE_MPFR
															ssmS,
#endif
															modelweight, modelnum);
					}
				
#else //NO_CCTABLE

				init_black(	blackB,
							//blackL,
							//blackK,
#ifdef USE_MPFR
							blackS,
#endif
							sModels,lPtable,
							fm, lm );
				init_subsums_from_black(	blackB,
											//blackL,
											//blackK,
											ssmB,
#ifdef USE_MPFR
											blackS,
											ssmS,
#endif
											Ptbl0,
											zipImg,
											sModels,
											fm, lm );
				totalLPB = init_LPB_from_subsums(	ssmB,
#ifdef USE_MPFR
													ssmS,
#endif
													modelweight, modelnum);
#endif
				totalLPB += PvalR;
				if( ( reporting && (debug_level > 0) ) ) 
				{
					if(  (chlistlen != 0) ) //(!reinit_for_fun) &&
					{
						printf("WARNING: changelist (%d) gracefully ignored\n", chlistlen);
						fprintf(logfile, "WARNING: changelist gracefully ignored\n");
					}
					printf("new Pval = %lf\n", totalLPB);
					if(debug_level > 1)
					{
						fprintf(logfile, "Pval = %lf\n", totalLPB);
						fflush(logfile);
					}
				}
				//if(!reinit_for_fun)
				//{
				
					qfr_out_str(totalLPB);
					int keep(0);
					//always accept the last changes so stats reflect last Pval
					if( !command.getCommand(jobfile?jobfile:stdin,"Keep", keep) || (keep == 0))
					{
						printf("ERROR: did not receive \"Keep(1)\" after reinitialization");
						fflush(stdout);
						quit_now = true;
					}
					continue;
				//}
				/*
				active_m = modelnum;
				chlistlen = 0; //disregard all changes
				sp_em = false;
				model_em = false;
				any_em = false;
				*/
				
			}
			
				
			//set pointers for READING the cached values
			if(active_m != modelnum)
			{
				modelweight[modelnum] = modelweight[active_m];
#ifdef USE_CCTABLE
				cctableBM = cctableB[active_m];
				//cctableLM = cctableL[active_m];
				//cctableKM = cctableK[active_m];
#endif
				ssmBM = ssmB[active_m];
#ifdef USE_MPFR
#ifdef USE_CCTABLE
				cctableSM = cctableS[active_m];
#endif
				ssmSM = ssmS[active_m];
#endif
			}
			//clear the new subsums
			memset(nssmB, '\0', sizeof(ssflt)*sfG);
#ifdef USE_MPFR
			memset(nssmS, '\0', sizeof(ssint)*sfG);
#endif
			//apply tentative changes to models //needed for black mode
			bool any_changes(false);
			short *diffmodel(nullptr);
#ifdef USE_CCTABLE			
			unsigned short * sModel(nullptr);		
#endif			
			if(sp_adjust)
			{
				
				for (int chnum( 0 ); chnum<chlistlen; ++chnum)
				{
					Models[ chx[chnum] + chy[chnum] * modelsize + active_m * modelarea ] = newval[chnum];
				}
			
				diffmodel = update_smodels(sModels, active_m, chlistlen, chx , chy, newval, /*oldval,*/ any_changes);
				if(!any_changes)
				{
					printf(reporting?"*":"NCH\n");
					fflush(stdout);
					int keep;
					if(command.getCommand(jobfile?jobfile:stdin,"Keep", keep))
					{	
						if(keep != 1)
						{
							for (int chnum=0; chnum<chlistlen; ++chnum)
							{
								Models[ chx[chnum] + chy[chnum] * modelsize + active_m * modelarea ] = oldval[chnum];
							}	
							//simply exchange newval <-> oldval
							update_smodels(sModels, active_m, chlistlen, chx , chy, oldval, /*newval,*/ any_changes);
						}
						continue;
					}
					else
					{ 
						printf("effectively useless ModelPixel was not followed by any keep\n");
						break;
					}
				}
				
#ifdef USE_CCTABLE					
				sModel = sModels + (12 * shadownum * active_m * modelarea);
				//printf("update_smodels OK\n");
#endif			
			}
			
#ifndef USE_CCTABLE
			memcpy(nblackB,blackB,sizeof(ssflt)*modelsp1);
			//memcpy(nblackL,blackL,sizeof(ssflt)*modelsp1);
			//memcpy(nblackK,blackK,sizeof(ssflt)*modelsp1);
#ifdef USE_MPFR
			memcpy(nblackS,blackS,sizeof(ssint)*modelsp1);
#endif
#endif
			if(!any_em && active_m != modelnum)
			{
				//STEP 1
				//calculate the cc values in the cc buffer line
#ifdef USE_CCTABLE
				if(false && !use_beam)
				{		
					num_changes = chlistlen;
					update_subsums_from_cache(	cctableBM,
												//cctableLM,
												//cctableKM,
												ncctableB,
												//ncctableL,
												//ncctableK,
												nssmB,
#ifdef USE_MPFR
												cctableSM,
												ncctableS,
												nssmS,
#endif
												qImg,
#ifdef MUL_PTBL
												mulPtbl,
#else
												Ptable, InvPtable,
#endif
												chx , chy, chlistlen,
												newval, oldval);
					printf("update_subsums_from_cache   OK\n");
					exit(0);						
				}							
				else
				{			
					num_changes = 
					update_subsums_diff_cc(	cctableBM,
											ncctableB,
									nssmB,								
#ifdef USE_MPFR								
									cctableSM,
									ncctableS,
									nssmS,
#endif
									qImg,
#ifdef MUL_PTBL
									mulPtbl,
#else
									Ptable, InvPtable,
#endif
									diffmodel, sModel);
					//printf("update_subsums_diff_cc OK\n");
					//fflush(stdout);
					//exit(0);				
				}


													
#else  // NO_CCTABLE

				update_black(	 blackB,  //blackL,  blackK, 
								nblackB, //nblackL, nblackK,
#ifdef USE_MPFR
								blackS, nblackS,
#endif
#ifdef MUL_PTBL
								mulPtbl,
#else
								Ptable, InvPtable,
#endif
								active_m,// chx, chy, 
								newval, oldval, chlistlen);
				//printf("update_black OK\n");				

				init_subsums_from_black(	nblackB,
											&nssmB,
#ifdef USE_MPFR
											nblackS,
											&nssmS,
#endif
											Ptbl0,
											zipImg,
											sModels,
											-1, active_m ); //-1 triggers update on active_m


/*
				update_subsums_from_black(	nblackB,
											//nblackL,
											//nblackK,
											nssmB,
#ifdef USE_MPFR
											nblackS,
											nssmS,
#endif
											Ptbl0,
											zipImg,
											sModels, active_m );
				//printf("update_subsums_from_black  OK\n");							
*/

#endif
			}
			//STEP 2
			//calculate totalLPB from the updated newsubsums
			if(!any_em)
			{
				totalLPB = init_LPB_from_subsums(	ssmB,
#ifdef USE_MPFR
													ssmS,
#endif
													modelweight, active_m);
				
				totalLPB += PvalR;
				//printf("total_LPB_from_subsums OK %lf\n",totalLPB);
				if(debug_level > 0)
				{
					if (reporting && ((report_runs % 50) != 49))
					{
						printf(".");
						fflush(stdout);
					}
					else
					{
						qfr_out_str(totalLPB);
					}
					fflush(stdout);
				}
				if (!reporting)
				{
					qfr_out_str(totalLPB);
				}
			
				int keep(0);
				//always accept the last changes so stats reflect last Pval
				if( !command.getCommand(jobfile?jobfile:stdin,"Keep", keep) )
				{
					quit_now = true;
					keep = 1;
				}
				if(reporting && (debug_level > 1) )
				{
					fprintf(logfile, "%lf", totalLPB);
					fprintf(logfile,"\t%s[%d]\n", command.getHead().c_str(), keep);
					fflush(logfile);
					
				}
				/*
				if(force_keep)
				{	
					printf("f");
					fflush(stdout);
					keep = 1;
				}*/
				//printf("Keep[%d]\n",keep);
				//STEP 3 take or leave the freshly calculated values
				if(keep==1)
				{
					if ( sp_adjust )
					{
						
#ifdef USE_CCTABLE
						//accept the buffered ccache line
						cctableB[active_m] = ncctableB;
						//cctableL[active_m] = ncctableL;
						//cctableK[active_m] = ncctableK;
						cctableB[modelnum] = cctableBM;
						//cctableL[modelnum] = cctableLM;
						//cctableK[modelnum] = cctableKM;
						ncctableB = cctableBM;
						//ncctableL = cctableLM;
						//ncctableK = cctableKM;
						cctableBM = nullptr;
						//cctableLM = nullptr;
						//cctableKM = nullptr;
#else
						memcpy(blackB,nblackB,sizeof(ssflt)*modelsp1);
						//memcpy(blackL,nblackL,sizeof(ssflt)*modelsp1);
						//memcpy(blackK,nblackK,sizeof(ssflt)*modelsp1);
#endif
						//accept the buffered subsum line
						ssmB[active_m] = nssmB;
						ssmB[modelnum] = ssmBM;
						nssmB = ssmBM;
						ssmBM = nullptr;
#ifdef USE_MPFR
						//do the same for the exponenents aka shifts
#ifdef USE_CCTABLE
						cctableS[active_m] = ncctableS;
						cctableS[modelnum] = cctableSM;
						ncctableS = cctableSM;
						cctableSM = nullptr;
#else
#ifdef USE_MPFR
						memcpy(blackS,nblackS,sizeof(ssflt)*modelsp1);
#endif
#endif
						ssmS[active_m] = nssmS;
						ssmS[modelnum] = ssmSM;
						nssmS = ssmSM;
						ssmSM = nullptr;
#endif
					}
					if (update_weights)
					{	memcpy(old_modelweight, modelweight, sizeof(ssflt)*modelsp1);}
				}
				else //keep != 1
				{
					//revert changes to models
					if(sp_adjust)
					{
						for (int chnum=0; chnum<chlistlen; ++chnum)
						{
							Models[ chx[chnum] + chy[chnum] * modelsize + active_m * modelarea ] = oldval[chnum];
						}
						//printf("reverting models OK\n");
						//simply exchange newval <-> oldval
						update_smodels(sModels, active_m, chlistlen, chx , chy, oldval, /*newval,*/ any_changes);
						//revert_smodels(sModels, active_m);
						//printf("revert_smodels   OK\n");
					}
					if (update_weights)
					{
						memcpy(modelweight, old_modelweight, sizeof(ssflt)*modelsp1);
					}
				}	
			}
			else if(sp_em)
			{
#ifdef USE_CCTABLE	
				ssflt meanB[maxchlistlen];
				ssflt wghtB[maxchlistlen];
#ifdef USE_MPFR
				ssint meanS[maxchlistlen];
				ssint wghtS[maxchlistlen];		
#endif				
				
				em_from_cache	(	
									cctableBM, ssmB,								
#ifdef USE_MPFR								
									cctableSM, ssmS,
#endif
									qImg,
									//modelweight,
									chx , chy, 
									meanB, wghtB,
#ifdef USE_MPFR						
									meanS, wghtS,
#endif						
									chlistlen, active_m
								);				
				
				for(int q(0); q < chlistlen; ++q)
				{
#ifdef USE_MPFR					
					printf("MEANB\t%lf\tMEANS\t%d\tWGHTB\t%lf\tWGHTS\t%d\n", meanB[q], meanS[q], wghtB[q], wghtS[q]);
#else
					printf("MEANB\t%lf\tWGHTB\t%lf\n", meanB[q], wghtB[q]);
#endif					
				}
#else // !USE_CCTABLE				
					printf("ERROR single pixel EM step requires cctable, recompile with -DUSE_CCTABLE\n");
#endif				
				
				fflush(stdout);
			
			
			
			}
			else if (model_em)
			{
#ifdef USE_CCTABLE				
				ssflt *EModelsB = new ssflt[modelnum*modelarea];
				memset(EModelsB,0,modelnum*modelarea * sizeof(ssflt));
				ssflt *wghtB = new ssflt[modelnum];
				memset(wghtB,0,modelnum * sizeof(ssflt));
#ifdef USE_MPFR				
				ssint *EModelsS = new ssint[modelnum*modelarea];
				memset(EModelsS,0,modelnum*modelarea * sizeof(ssint));
				ssint *wghtS = new ssint[modelnum];
				memset(wghtB,0,modelnum * sizeof(ssint));		
#endif				
				
				emodel_from_cache(	cctableB, ssmB,	EModelsB, wghtB,												
#ifdef USE_MPFR								
									cctableS, ssmS, EModelsS, wghtS,
#endif
									zipImg//const unsigned short *qImg,
									//modelweight, totalLPB-PvalR, fm, lm
								);
					
				
				sendEModels(	EModelsB, wghtB,
#ifdef USE_MPFR
								EModelsS, wghtS, 
#endif
								reporting,
								fm, lm);
				delete[] EModelsB;
				delete[] wghtB;
#ifdef USE_MPFR
				delete[] EModelsS;
				delete[] wghtS;
#endif					
#else //!USE_CCTABLE
				printf("ERROR full model EM step requires cctable, recompile with -DUSE_CCTABLE\n");
				fflush(stdout);
#endif			
			}
			else
			{
				printf("ERROR unrecognized task requested!\n");
				fflush(stdout);
			}
			total_runs += ( (num_changes > 0 ? num_changes : chlistlen) + update_weights + (any_em?1:0));
			++report_runs;
		} while( !quit_now );

		if( reporting) 
		{
			clock_t now = clock();
			double elapsed = double(now-timer)/CLOCKS_PER_SEC;
			qfr_out_str(totalLPB);
			printf("%d iterations finished in ", total_runs);
			printf("%fs\t (Test/second = %f)\n", elapsed, total_runs/elapsed);
		}
		if(true && reporting) ///DEBUG reinitialize once again at the end
		{	
			//timer = now;
#ifdef USE_CCTABLE
			assert(cctableB[modelnum] == ncctableB);
			//assert(cctableL[modelnum] == ncctableL);
			//assert(cctableK[modelnum] == ncctableK);
#endif
			assert(ssmB[modelnum] == nssmB);
#ifdef USE_MPFR
#ifdef USE_CCTABLE
			assert(cctableS[modelnum] == ncctableS);
#endif
			assert(ssmS[modelnum] == nssmS);
#endif
			cache_report(	reportfile,
							blackB,
							//blackL,
							//blackK,
							ssmB,
#ifdef USE_CCTABLE
							cctableB,
							//cctableL,
							//cctableK,
#endif
#ifdef USE_MPFR
							blackS,
							ssmS,
#ifdef USE_CCTABLE
							cctableS,
#endif
#endif
							cases, reporting);
		
		
			if(reporting)
			{
				//reinit_for_fun = true;
#ifdef USE_CCTABLE
					init_smodels(sModels,Models,Beam,use_beam);
					//printf("init_smodels (re_fresh)  OK\n");
					if(!use_init_black)
					{
				 
						//This variant is about half the performance, but may become
						//necessary at higher doses when black probabilities are no longer feasible
						totalLPB = init_hexlikelyhood(	ssmB,
														cctableB,
														//cctableL,
														//cctableK,
#ifdef USE_MPFR
														ssmS,
														cctableS,
#endif
														lPtable,
														zipImg,
														sModels, modelweight );
					}
					else //init with black values
					{
				 
						init_black( blackB,
									//blackL,
									//blackK,
#ifdef USE_MPFR
									blackS,
#endif
									sModels,lPtable );

						init_cache_from_black(
									blackB,
									//blackL,
									//blackK,
									cctableB,
									//cctableL,
									//cctableK,
									ssmB,
#ifdef USE_MPFR
									blackS,
									cctableS,
									ssmS,

#endif
									Ptbl0,
									zipImg,
									sModels );
/*						
						printf("black cctables done\n");
						totalLPB = init_hexlikelyhood(	ssmB,
														cctableB,
														//cctableL,
														//cctableK,
#ifdef USE_MPFR
														ssmS,
														cctableS,
#endif
														lPtable,
														zipImg,
														sModels, modelweight );
						printf("asserted cctables done\n");			
*/
						totalLPB = init_LPB_from_subsums(	ssmB,
#ifdef USE_MPFR
															ssmS,
#endif
															modelweight, modelnum);
					}
				
#else //NO_CCTABLE

				init_black(	blackB,
							//blackL,
							//blackK,
#ifdef USE_MPFR
							blackS,
#endif
							sModels,lPtable );
				init_subsums_from_black(	blackB,
											//blackL,
											//blackK,
											ssmB,
#ifdef USE_MPFR
											blackS,
											ssmS,
#endif
											Ptbl0,
											zipImg,
											sModels );
				totalLPB = init_LPB_from_subsums(	ssmB,
#ifdef USE_MPFR
													ssmS,
#endif
													modelweight, modelnum);
#endif
				
				
				totalLPB += PvalR;
				printf("Pval = %lf\n", totalLPB);
				
				
				cache_report(reportfile,
					blackB,
					//blackL,
					//blackK,
					ssmB,
#ifdef USE_CCTABLE
					cctableB,
					//cctableL,
					//cctableK,
#endif
#ifdef USE_MPFR
					blackS,
					ssmS,
#ifdef USE_CCTABLE
					cctableS,
#endif
#endif
					cases, reporting);
				
				
			}// end if(true)
				
		}
	}

    for(int m = 0; m<modelsp1; ++m)
    {
#ifdef USE_CCTABLE
		delete[] cctableB[m];
		//delete[] cctableL[m];
		//delete[] cctableK[m];
#endif
		delete[] ssmB[m];
	}

#ifdef USE_MPFR
    for(int m = 0; m<modelsp1; ++m)
    {
#ifdef USE_CCTABLE
		delete[] cctableS[m];
#endif
		delete[] ssmS[m];
	}
#endif

	delete[] blackB;
	//delete[] blackL;
	//delete[] blackK;
	delete[] ssmB;
#ifdef USE_CCTABLE
    delete[] qImg;
    delete[] cctableB;
    //delete[] cctableL;
    //delete[] cctableK;
#else
	delete[] nblackB;
	//delete[] nblackL;
	//delete[] nblackK;
#endif


#ifdef USE_MPFR
#ifdef USE_CCTABLE
	delete[] cctableS;
#else	
	delete[] nblackS;
#endif

	delete[] ssmS;
	delete[] blackS;
#endif	

#ifdef MUL_PTBL
    delete[] mulPtbl;
#else
	delete[] PtableL;
    delete[] InvPtableL; 
#endif	
	delete[] Ptbl0; 
	
	
	free(zipImg);
	delete[] Models;
    delete[] sModels;
    delete[] hexPixels;
    delete[] pixpos;
    delete[] transformations;
    delete[] inv_transformations;
    delete[] translations;
    delete[] inv_translations;
    delete[] mirrot;
    delete[] inv_mirrot;
	delete[] latticePoints;
	delete[] lpX;
	delete[] lpZ;
	delete[] modelweight;
	delete[] old_modelweight;
	delete[] cases;
	close_smodels();
    if(debug_level > 0)
    {	fclose(logfile);}
	reportfile.close();
	
#ifdef USE_MPI	
	MPI_Finalize();
#endif	
    return 0;
}


