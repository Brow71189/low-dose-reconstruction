
#ifndef FAST_CACHE_H
#define FAST_CACHE_H

#include "globals.h"

typedef double ssflt;
typedef int ssint;



#ifdef USE_MPFR

// in mpfr_mode we want to use smaller types
//maybe long double can handle outliers better? No it cannot
//Only removing them with Outlier_Remover beforehand does the trick
typedef double ccflt; 
typedef int ccint;
//tolerance level for exponent adjustment
#define CC_EXP_MAX 1000000
#else
//there is no reason to treat cc and ss differently
//tolerance level for exponent adjustment
#define CC_EXP_MAX 800
typedef ssflt ccflt;
typedef ssint ccint;
#endif

//Macros beat inline functions
//WARNING beware of memcpy and pointerarithmetics in hexworker.cpp when changing these!




//#define TRANSLATE(hp,lp) (TRANSFORM(hp,lp,0,0))
//forward order
//#define TRANSFORM(hp,lp,rot,mirror) (mirror + 2 * rot + 12 * lp + globalsympersf * hp)
//#define CCINDEX(m,sfC,lp,rot,mirror) (sfC + sfG * mirror + sfmG * rot + sfmrG * lp + sfmrlG * m)
//#define CCIND(sfC,lp,rot,mirror) (sfC + sfG * mirror + sfmG * rot + sfmrG * lp)
//#define MIRROT(hp,rot,mirror) (mirror + 2 * rot + 12 * hp)



//inverse order
#define MIRROT(hp,rot,mirror) (rot + 6 * mirror + 12 * hp)
#define TRANSFORM(hp,lp,rot,mirror) (lp + rot*globallpMax + mirror * 6 * globallpMax + globalsympersf * hp)
//                                sfC + sfG * lp + rot * sfG*globallpMax + mirror * 6*sfG*globallpMax 
#define CCIND(sfC,lp,rot,mirror) (sfC + sfG * lp + rot * sfGlpM + mirror * sfGlpMrot)
#define CCIND00(rot,mirror) (rot * sfGlpM + mirror * sfGlpMrot) 

//same ordering for forward and inverse
#define MOVE(hp,lp) (lp + hp * globallpMax)
#define MOVEHP0(hp) (hp*globallpMax) 
#define TRANSSYM(hp,sym) (sym + globalsympersf * hp)
#define SSMID(m,sfC) (sfC + m * sfG)
#define CCINDEX(m,sfC,lp,rot,mirror) ( CCIND(sfC,lp,rot,mirror) + sfmrlG * m)

#endif 
