#include "stdio.h"
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <execinfo.h>
#include <sched.h> 

#include <time.h>
#include <assert.h>
#include <exception>
#include <math.h>
#include "string.h"
#include <sstream>
#include <fstream>

#include "globals.hpp"
#include "basis.hpp"
#include "command.hpp"
#include "communicator.hpp"
#include "optimizer.hpp"
#include "hex_optimizer.hpp"
#include "subframes.hpp"
#include "unitcell.hpp"
#include "xorshift1024star.hpp"


//private headers
void byebye(void);
void sigint_handler(int sig);
void segfault_handler(int sig);
void bad_signal_handler(int sig);
void abortfunc(mcheck_status mstatus);
int sane_input(void); // 0 everying alright, -1bad

void sigint_handler(int sig)
{	
	if(busy)
	{	
		if(reporting && !has_console)
		{	printf("Message ");} 
		if(reporting)
		{	printf("final report requested, exiting and deleting %s at next SIGINT\n", logname.c_str()); }
		fflush(stdout);
		config.fresh_config = true;
		busy = false;
	}
	else
	{	
		if(reporting && !has_console)
		{	printf("Message ");} 
		printf("quitting due to SIGINT, deleting %s\n", logname.c_str() );
		fflush(stdout);
		fclose(logfile);
		std::remove(logname.c_str());
		exit(sig);
	}
}

void segfault_handler(int sig)
{	
	fprintf(stderr,"ERROR worker%d frame: %d SEGMENTATION FAULT log file: %s\n",threadnum, frameNr, logname.c_str());
	
	fflush(stderr);
	fprintf(logfile,"ERROR worker%d frame: %d SEGMENTATION FAULT log file: %s\n",threadnum, frameNr, logname.c_str());
	fclose(logfile);
	exit(sig);
}

void bad_signal_handler(int sig)
{
	fprintf(stderr,"ERROR worker%d frame: %d terminated with signal: %d log file: %s\n",threadnum, frameNr, sig, logname.c_str());
	fflush(stderr);
	fprintf(logfile,"ERROR worker%d frame: %d terminated with signal: %d log file: %s\n",threadnum, frameNr, sig, logname.c_str());
	fclose(logfile);
	exit(sig);
}

void evil_signal_handler(int sig)
{
	fprintf(stderr,"EVIL ERROR worker%d frame: %d terminated with signal: %d log file: %s\n",threadnum, frameNr, sig, logname.c_str());
	fflush(stderr);
	fprintf(logfile,"EVIL ERROR worker%d frame: %d terminated with signal: %d log file: %s\n",threadnum, frameNr, sig, logname.c_str());
	fclose(logfile);
	exit(sig);
}

void abortfunc(mcheck_status mstatus)
{
	if(reporting && !has_console)
	{	printf("Message ");}
	fprintf(stderr,"MEMORY ERROR DETECTED type:%d %s\n", mstatus, mem_check_location.c_str());
	fflush(stderr);
	fprintf(logfile,"MEMORY ERROR DETECTED type:%d %s\n", mstatus, mem_check_location.c_str());
	if(reporting && !has_console)
	{	printf("Message ");}
	switch( mstatus )
	{
		case MCHECK_HEAD:
			fprintf(stderr,"MCHECK_HEAD memory access before chunk\n");
			fprintf(logfile,"MCHECK_HEAD memory access before chunk\n");
		break;
		
		case MCHECK_TAIL:
			fprintf(stderr,"MCHECK_TAIL memory access after chunk\n");
			fprintf(logfile,"MCHECK_TAIL memory access after chunk\n");
		break;
		
		case MCHECK_FREE:
			fprintf(stderr,"MCHECK_FREE double free\n");
			fprintf(logfile,"MCHECK_FREE double free\n");
		break;
		
		case MCHECK_OK:
			fprintf(stderr,"MCHECK_OK this should not have been reported\n");
			fprintf(logfile,"MCHECK_OK this should not have been reported\n");
		break;
		case MCHECK_DISABLED:
			fprintf(stderr,"MCHECK_DISABLED initialization of mcheck() has failed\n");
			fprintf(logfile,"MCHECK_DISABLED initialization of mcheck() has failed\n");
		break;
		default:
			fprintf(stderr,"UNSPECIFIED STATUS\n");
			fprintf(logfile,"UNSPECIFIED STATUS\n");
	}
	fclose(logfile);
	fprintf(stderr,"ERROR worker%d frame: %d log file: %s\n",threadnum, frameNr, logname.c_str());
	fflush(stderr);
	exit(mstatus);
}

int sane_input(void) // 0 everying alright, -1bad
{
	int result( 0 );
	
	if (frameNr < 1) 
	{ 
		printf("Message ERROR: frameNr: %d < 1\n", frameNr);
		result = -1;
	} 
	if (modelsize < 1) 
	{	
		printf("Message ERROR: modelsize: %d < 1\n", modelsize);
		result = -1;
	}
	if ( (bondlength_CC < 1) || (bondlength_CC%2 != 0) ) 
	{
		printf("Message ERROR: bondlength: %d < 1 || %d%%2 != 0\n", bondlength_CC, bondlength_CC);
		result = -1;
	}
	if (modelsize%bondlength_CC != 0) 
	{	
		printf("Message ERROR: modelsize%%bondlength: %d%%%d != 0\n", modelsize, bondlength_CC);
		result = -1;
	}
	if (solidsize%bondlength_CC != 0) 
	{	
		printf("Message ERROR: solidsize%%bondlength: %d%%%d != 0\n", solidsize, bondlength_CC);
		result = -1;
	}
	if (solidsize<1 ) 
	{	
		printf("Message ERROR: solidsize: %d < 1\n", solidsize);
		result = -1;
	}
	if (solidsize > modelsize ) 
	{	
		printf("Message ERROR: solidsize: %d > modelsize: %d\n", solidsize,modelsize);
		result = -1;
	}
	if ( (bondlength_Mini < 1) || (bondlength_Mini%2 != 0) ) 
	{
		printf("Message ERROR: minilength: %d < 1 || %d%%2 != 0\n", bondlength_Mini, bondlength_Mini);
		result = -1;
	}
	/*
	if (modelsize%bondlength_Mini != 0) 
	{	
		printf("Message ERROR: modelsize%%minilength: %d%%%d != 0\n", modelsize, bondlength_Mini);
		result = -1;
	}
	*/
	if ((impWidth * impHeight != framelen) || (impWidth <= 0) || (impHeight <= 0) ) //should be impossible
	{	
		printf("Message ERROR: frame dimensions: %d x %d != %d\n", impWidth, impHeight, framelen);
		result = -1;
	}
	if (fov < 0.0) 
	{	
		printf("Message ERROR: field of view: %lf <= 0.0\n",fov);
		result = -1;
	}
	//These should be impossible
	if (ea <= 0.0)
	{
		printf("Message ERROR: ellipse_A: %lf <= 0.0\n",ea);
		result = -1;
	}
	if (eb <= 0.0)
	{
		printf("Message ERROR: ellipse_B: %lf <= 0.0\n",eb);
		result = -1;
	}
	return result;
}


int main(int argc, char* argv[])
{
	mem_check_location.assign("not available"); //default before the first call to memory_check()
	int mcheck_result = 1;
#ifndef NOCATCH
	mcheck_result = mcheck( &abortfunc );	
#endif	
	signal(SIGHUP, &bad_signal_handler); //parent has died
	signal(SIGINT, &sigint_handler); //ctrl c from keyboard
	signal(SIGQUIT, &bad_signal_handler); //ctrl z from keyboard
	signal(SIGILL, &bad_signal_handler); //illeagl instruction
	signal(SIGABRT, &bad_signal_handler); //abort() e.g. exception aka SIGIOT
	signal(SIGFPE, &bad_signal_handler);  //floating point exception
	signal(SIGKILL, &bad_signal_handler); //could not be catched
	if(mcheck_result != 0) //check if we did install mcheck
	{	
		signal(SIGSEGV, &segfault_handler);
	} 
	else //we can still catch the signal
	{
		signal(SIGSEGV, &bad_signal_handler); 
	}
	signal(SIGPIPE, &bad_signal_handler); //broken pipe
	signal(SIGALRM, &evil_signal_handler); //timer expired
	signal(SIGTERM, &bad_signal_handler); // terminate()
	signal(SIGUSR1, &evil_signal_handler); //user signal1
	signal(SIGUSR2, &evil_signal_handler); //user signal2
	//SIGCHLD //we have no children
	/*
	Dont mess with stop and continue signals
	-> SIGSTOP would be a perfect silent killer 
	*/
	signal(SIGBUS, &bad_signal_handler); //physical memory error
	signal(SIGPOLL, &bad_signal_handler); //pollable avent aka SIGIO
	signal(SIGPROF, &bad_signal_handler); //profiler timing expired
	signal(SIGSYS, &bad_signal_handler); // bad argument to routine
	signal(SIGTRAP, &bad_signal_handler); //trace/breakpoint trap
	signal(SIGVTALRM, &evil_signal_handler); //virtual alram clock expired
	signal(SIGXCPU, &evil_signal_handler);   // cpu time limit exceeded
	signal(SIGXFSZ, &bad_signal_handler); //max file size exceed
	
	signal(SIGSTKFLT, &bad_signal_handler); //stack fault
	
	signal(SIGIO, &bad_signal_handler); // I/O now possible
	signal(SIGPWR, &bad_signal_handler); //Power failure
	
	
	Command command("init(1)");
	globalCommand = &command;
	
	
	if(argc > 1)
	{
		if( strcmp(argv[1],"-v")==0 )
		{
			printf("This is hexsampler compiled at %s %s\n", __DATE__, __TIME__ );
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
			puts("This software has not been released under any license");
			puts("You are not supposed to use it or even trust the guy who gave it to you");
			puts("The source code is maintained by Christian.Kramberger-Kaplan@univie.ac.at");	
			puts("compiled settings defined in globals.hpp");
			printf("SAMPLE_RATE: %d\n", SAMPLE_RATE);
			printf("MIN_ITER: %d\tMAX_ITER: %d\tMAX_ITER_LIMIT: %d\n", MIN_ITER, MAX_ITER, MAX_ITER_LIMIT);
#ifdef NOCATCH
			puts("CATCH exceptions in main is DISABLED");			
#else
			puts("CATCH exceptions in main is ENABLED");
			printf("mcheck is %sAVAILABLE\n", mcheck_result?"NOT ":"");
#ifdef ACTIVE_CHECKS
			puts("active memory checking is ENABLED");
#else
			puts("active memory checking is DISABLED");	
#endif
			
#endif			
			return 0;
		}
		if( strcmp(argv[1],"-h")==0 )
		{
			puts("hexsampler identifies the lattice parameters of graphene,");
			puts("resamples a 16bit 2D array into hexagonal pixels and slices it into supercells (aka models)");
			printf("excluded areas can be marked with %d\n", dirt_val);
			puts("hexsampler reads commands and binary data (BigEndian) from stdin and returns text formated data to stdout");
			puts("the comand is the first token of a line and followed by one argument and a newline.");
			puts("spaces, tabs or () or [] maybe used to sepparate command and argument");
			puts("lines starting with # will be ignorred");
			puts("The following commands are used for I/O and are symmetric:\n");
			
			puts("Run(0|1)\t\t exits | validates the necessary input and launches the job");
			puts("Task(%d)\t\t 0 .. optimize, 1 .. reslice models, 2 .. poll models, 3 .. debug");
			puts("ThreadNum(%d)\t\t sets id of this worker");
			puts("Frame(%d)\t\t sets id of the image to process");
			puts("ImpWidth(%d)\t\t width of the image in pixels");
			puts("ImpHeight(%d)\t\t height of the image in pixels");
			puts("updateFrameData()\t signal begin image data");
			puts("BeginPixels %d\t\t header for text pixels argument is number of pixels");
			puts("EndPixels %d\t\t footer for text pixels argument is checksum");
			puts("HexImgScale(%lf)\t canvas scaling for rotated hex pixel image");
			puts("Resolution(%d)\t\t resolution for sampling of rhombic unitcell");
			puts("BondLength(%d)\t\t hex pixels per C-C bond");
			puts("ModelSize(%d)\t\t diameter in hex pixels");
			puts("MiniLength(%d)\t\t super resolution in smaller hex pixels per bondlength");
			puts("FieldofView(%lf)\t ImpWidth in nm");
			puts("OffsetX(%lf)\t\t origin pixels of lattice (hollow site)");
			puts("OffsetY(%lf)\t\t origin pixels of lattice (hollow site)");
			puts("HexEdgeLen(%lf)\t\t length of an edge of a hex pixel in pixels");
			puts("Tilt(%lf)\t\t tilt angle of graphene basis vectors");
			puts("RotateTilt(%d)\t\t rotate basis vectors in 60° steps");
			puts("EllipseA(%d)\t\t long axis of first Bragg reflex ellipsoid in pixels");
			puts("EllipseB(%d)\t\t short axis of first Bragg reflex ellipsoid in pixels");
			puts("EllipsePhi(%d)\t\t orientation of long axis");
			puts("OffsetHexQ(%d)\t\t offset in unitcells for slicing models");
			puts("OffsetHexR(%d)\t\t offset in unitcells for slicing models");
			puts("NumSF(%d)\t\t number of complete models");
			puts("Raw_Avg(%lf)\t\t average counts of pixels in rhombic unitcell");
			puts("Hex_Avg(%lf)\t\t average counts of hex pixels");
			puts("Hex_Low(%lf)\t\t average counts of hex pixels at hollow sites");
			puts("Hex_High(%lf)\t\t average counts of hex pixels at atoms");
			puts("Contrast(%d)\t\t contrast of the lattice");
			puts("MinHexContrast(%lf)\t contrast required to employ hex merit");
			puts("Hex_Merit(%lf)\t\t advanced black magic on super resolution hexagonal unitcell");
			puts("Hex_Shape(%d)\t\t the core of the advanced black magic");
			puts("Smoothing(%d)\t\t smoothing passes before evaluating merits");
			puts("Merit(%lf)\t\t black magic on rhombic unitcell (backup)");
			puts("Mirror_Q(%d)\t\t measure of mirror symmetry in rhombic unitcell");
			puts("UseMirror(0|1)\t\t switches if mirror symmetry should enter merit and hex_merit");
			puts("Position_Q(%lf)\t\t measure of localization of atoms in rhombic unitcell");
			puts("UsePosition(0|1)\t switches if position_Q enters merit");
			puts("Stability(%d)\t\t 0,1,2,3 widens the search range");
			puts("PeakPoints(%d)\t\t repetitions to level fluctuations in sampling");
			puts("ReportUC(0|1)\t\t switches reporting of rhombic unitcell, position map and super resolution hexagon");
			puts("ReportHexImg(0|1)\t switches reporting of hex pixeled overview image");
			puts("ReportSF(0|1)\t\t switches reporting the sliced models");
			puts("Reporting(0|1)\t\t switches rudimentary extra output");
			puts("HasConsole(0|1)\t\t switches verboose extra output");
			
			puts("\nThe following commands and keywords can be used for networking and are asymmetric\n");
			
			puts("LAUNCH\t\t\t at the beginning of a line confirms valid input");
			puts("REPORT\t\t\t signals that a report is ready");
			puts("Report(1)\t\t signals that the report can now be received");
			puts("Elapsed %lf\t\t tells elapsed time in seconds and also marks the end of a report");
			puts("sumUC(%d)\t\t signals an averaged rhombic unitcell argument is pixel area");
			puts("symUC(%d)\t\t signals during optimization a rhombic localization map");
			puts("         \t\t or the final rhombic unitcell rotated by 60° argument is pixel area");
			puts("miniUC(%d)\t\t signals a super resolution hexagonal unitcell argument is pixel area");
			puts("sfSum(%d)\t\t signals an averaged model argument is pixel area");
			puts("hexImg(%d)\t\t signals a resampled overview image argument is pixel area");
			puts("sfStack(%d)\t\t signals all models are send argument is pixel area per model");
			//puts("requestingTarget\t worker will receive another target unitcell after the report");
			
			puts("updateUCTarget(0|1)\t master will or will not update the target unitcell");
			puts("SendPixels(1)\t\t pixels of one image may be printed as single column text");
			puts("             \t\t every sliced model in sfStack has to be requested sepparately");
			puts("Message\t\t\t may preceed any message that should be shown to the user");
			puts("ERROR\t\t\t preceeds a serious error");
			puts("WARNING\t\t\t preceeds any warnings");
			puts("BeginBinary(%d)\t\t marks start of binary stream argument is number of bytes");
			puts("EndBinary(%d)\t\t argument is int32 checksum, may overflow");
			return 0;
		}
		
		printf("unknown commandline option %s\n", argv[1]);
		puts("-v shows version info and compiled FLAGS");
		puts("-h lists the supported commands and keywords");
		
		return 0;
    }
	int pid = (int)getpid();
	std::stringstream lognamestream;
	lognamestream << "log_" << pid << ".txt";
	logname = lognamestream.str();
	logfile = fopen(logname.c_str(),"wt");
	fprintf(logfile,"logfile created\n");fflush(logfile);
#ifndef NOCATCH    
    try
    { 
#endif		
		int frame_counter(0);
		while( read_input() != -1 ) // -1 from Run(0)
		{
			if(	sane_input() != 0) 
			{	
				printf("Message WARNING: skipping frame %d\n", frameNr);
				fflush(stdout);
				continue; //ignore task and await fresh input
			}
			
			if(frame_counter == 0)
			{
				srand(time(0));
				for(int i(0); i<threadnum; ++i)
				{
					srand((unsigned int)rand());
				}
			}
			init_xorshift1024star();//uses lower bits of rand()
			
			
			timer = clock(); //we have read and validated the input
			elapsed = 0.0;
			busy = false;
			config.fresh_config = false;
			bondlength = &bondlength_CC; //the default behaviour
			++frame_counter;  
			printf("LAUNCH %d on frame %d task: %d log: %s\n",frame_counter,frameNr,op_mode,logname.c_str());fflush(stdout);
			fprintf(logfile,"LAUNCH %d on frame %d task: %d\n",frame_counter,frameNr,op_mode);fflush(logfile);
			if(reporting && !has_console)
			{	printf("LAUNCH %d on frame %d log: %s\n",frame_counter,frameNr,logname.c_str());fflush(stdout);}
			switch(op_mode)
			{
				case 0:
					busy = true;
					optimize_merit(true); // set tilt,hel,phi,excent,offset_X & offset_Y
					//refesh and resample
					optimize_subframes(true,true); //set Num_sf,offsetQ,offsetR
					break;
				case 1: //rather quick no busy mechanism
					//fallthrough
				case 2: //rather quick no busy mechanism
					if(retune_origin)
					{	
						center_offset();
						floating_offset = true;
						ruling_merit = 1;
						optimize_hexagons(0,0); //retune offset_X and offset_Y
						ruling_merit = 0;
						offset_XS = offset_X;
						offset_YS = offset_Y;
						double m_std( -1.0 );
						merit = mean_merit(m_std,false);
					}
					else
					{
						floating_offset = false;
						ruling_merit = 0;
						double m_std( -1.0 );
						merit = mean_merit(m_std,false);
						ruling_merit = 1;
						double hm_std( -1.0 );
						hexagonal_merit = mean_merit(hm_std,false); //just reflect the current state
					}
					//refresh and resample
					
					optimize_subframes(true,true); //set Num_sf,offsetQ,offsetR	
					break;
				case 3: //DEBUG does a short run using all features
					busy = true;
					optimize_merit(false);
					optimize_subframes(true,true);
					break;
/*				
				case 4: //simply reflect the current state
					floating_offset = false;
					ruling_merit = 0;
					{
						double m_std( -1.0 );
						merit = mean_merit(m_std,false);//with statistics
						ruling_merit = 1;
						double hm_std( -1.0 );
						hexagonal_merit = mean_merit(hm_std,false); //with statistics
					}
					optimize_subframes(false,true); //use existing Num_sf, offsetQ, offsetR
					break;
*/				
				default:
					if(reporting && !has_console)
					{	printf("Message "); } 
					printf("ERROR unrecognzed op_mode %d in %s line: %d\n", op_mode, __FILE__, __LINE__);
					fflush(stdout);
					exit(0);	
			}
			busy = false; //the next report is the final one
			config.fresh_config = true;//and it will always be send 
			report();
			free_globals();
		}
#ifndef NOCATCH
	}
	catch (std::exception &e)
	{	//We never want to see that output //We are not catching any segfaults with gcc
		fprintf(logfile, "top level exception in main: %s\n",e.what());
		fclose(logfile);
		printf("Message ERROR top level exception in main\n");
		printf("Message keeping logfile %s\n", logname.c_str());
		puts( e.what() );
		fprintf(stderr, "top level exception in main: %s\n",e.what());
		fflush(stdout);
		fflush(stderr);
		return -1;
	}
#endif
	fflush(stdout);
	fflush(stderr);
	fclose(logfile);
	std::remove(logname.c_str()); //a clean regular exit
	return 0;
}


