
#include "globals.h"

FILE* jobfile = 0;

Command *globalCommand = 0;
int globalModelSize(-1), globalModelArea(-1);
int globalBeamSize(-1), globalBeamArea(-1), globalBeamRadius(-1), globalBeamNorm(-1);
int globalSubFrameSize(-1), globalSubFrameArea(-1);
int globalDatamax(-1), globalDatamin(-1), globalModelval(-1);
int globalPtableWidth(-1);
int globalPtableSize(-1);
int globalModelNum(-1);
int globalSubFrameCount(-1);
int globalqImgSize(-1);
int globalqImgCol(-1);

int globalhpMax(-1);
int globallpMax(-1);
int global_lp0(-1);
int globalccpersf(-1);
int globalccperm(-1);
int globalsymperm(-1);
int globalsympersf(-1);
int globalShadowNum(1);//the first shadow is the model itself

int globalMx(-1);
int globalMz(-1);

int pix_non_zero(0);
long zip_img_cap(0);
unsigned short* zip_img_sfP = 0;

//constants for ccindex
int sfG(-1);
int sfmG(-1);
int sfmrG(-1);
int sfmrlG(-1);
int sfGlpM(-1);
int sfGlpMrot(-1);


int* pixposG = 0;
unsigned short* transformationsG = nullptr; 
unsigned short* inv_transformationsG = nullptr;
unsigned short* translationsG = nullptr;
unsigned short* inv_translationsG = nullptr;
unsigned short* mirrotG = nullptr;
unsigned short* inv_mirrotG = nullptr;

int transElementsG = -1;
int moveElementsG = -1;
int mirrotElementsG = -1;

int* globalhexPixels = 0;
int* globallatticePoints = 0;
int* globallpX = 0;
int* globallpZ = 0;

double pvaloffset(0.0);
double pvalscaling(1.0);
double inv_pvalscaling(1.0);
bool pvaloffset_ok(false);
bool use_weights(true);
bool real_log2(false);

//bool reinit_for_fun(false);

#ifdef USE_MAX
bool use_max(false);
#endif
