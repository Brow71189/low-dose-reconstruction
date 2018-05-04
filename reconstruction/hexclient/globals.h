#ifndef GLOBALS_H
#define GLOBALS_H
#include "command.h"
#include "stdio.h"

#ifdef USE_MPFR
//Comment the next line to have fake mpfr behavior with exponent always 0
#define MPFR_TRUE
#endif
#define MAXCHLISTLEN 1024
#define MAXSHADOWS 32

/*
 * Global Variables are global handles/copies to the local variables in speed optimized main()
 *
 */
extern FILE* jobfile;
 
extern int pix_non_zero;
extern long zip_img_cap;
extern unsigned short* zip_img_sfP;
extern int globalModelSize, globalModelArea;
extern int globalBeamSize, globalBeamArea, globalBeamRadius, globalBeamNorm;
extern int globalSubFrameSize, globalSubFrameArea;
extern int globalModelNum, globalSubFrameCount;
extern int globalShadowNum;

extern int globalhpMax,globallpMax, global_lp0;
extern int globalMx, globalMz;
extern int globalccpersf;
extern int globalsymperm;
extern int globalccperm;
extern int globalsympersf;
extern int* globalhexPixels;
extern int* globallatticePoints;
extern int* globallpX;
extern int* globallpZ;


/*constants for ccindex*/
extern int sfG;
extern int sfmG;
extern int sfmrG;
extern int sfmrlG;


extern int sfGlpM;
extern int sfGlpMrot;

extern int* pixposG;
extern unsigned short* transformationsG;
extern unsigned short* inv_transformationsG;
extern unsigned short* translationsG;
extern unsigned short* inv_translationsG;
extern unsigned short* mirrotG;
extern unsigned short* inv_mirrotG;

extern int transElementsG;
extern int moveElementsG;
extern int mirrotElementsG;


extern int globalDatamin, globalDatamax, globalModelval;
extern int globalPtableWidth, globalPtableSize;
extern int globalqImgSize, globalqImgCol;

extern double pvaloffset;
extern double pvalscaling;
extern double inv_pvalscaling;
extern bool pvaloffset_ok;
extern bool use_weights;
extern bool real_log2;

//extern bool reinit_for_fun;  //used for DEBUGING 

#ifdef USE_MAX
extern bool use_max;
#endif
extern Command *globalCommand;


#endif
