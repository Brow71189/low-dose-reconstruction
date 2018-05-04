#ifndef FAST_MPFR_H
#define FAST_MPFR_H

#include <math.h>
#include "fast_cache.h"
#include <assert.h>

#ifdef USE_MPFR
/*
inline __attribute__((always_inline))
void qfr_norm(ssflt &xb, ssint &xs)
{
#ifdef MPFR_TRUE	
	int tmp;
	xb = frexp(xb, &tmp);
	xs += (ssint)tmp;
#else
	assert(xs == 0);
#endif
}
*/
#if defined(MACRO_MPFR) && defined(MPFR_TRUE)
#define qfr_norm(xb,xs) \
		{ \
			int tmp; \
			*(xb) = frexp(*(xb), &tmp); \
			*(xs) += (ssint)tmp; \
		}		 
#else //MACRO_MPFR
inline __attribute__((always_inline))
void qfr_norm(ssflt *xb, ssint *xs) //pointer version
{
#ifdef MPFR_TRUE	
	//if(*xb != 0.0)
	//{
		int tmp;
		*xb = frexp(*xb, &tmp);
		*xs += (ssint)tmp;
	/*}
	else
	{
		*xs = 0;
	}*/
#else
	assert(*xs == 0);
#endif
}

#endif //MACRO_MPFR
#ifndef MPFR_TRUE
inline __attribute__((always_inline))
void qfr_assert_norm(const ssint xs)  //Ja Ja ineffizient, aber wurscht weil inline
{
	assert(xs == 0);
}
#define qfr_assert_norm(b,s) qfr_assert_norm((s))
#else //MPFR_TRUE
inline __attribute__((always_inline))
void qfr_assert_norm(const ssflt xb)  //Ja Ja ineffizient, aber wurscht weil inline
{
	assert( (xb >= 0.5) && (xb < 1.0) );
}
#define qfr_assert_norm(b,s) qfr_assert_norm((b))
#endif //MPFR_TRUE
/*
inline __attribute__((always_inline))
void qfr_assert_norm(const ssflt xb, const ssint xs)  //Ja Ja ineffizient, aber wurscht weil inline
{
#ifdef MPFR_TRUE
	assert( (xb >= 0.5) && (xb < 1.0) );
#else
	assert(xs == 0);
#endif
}
*/
inline __attribute__((always_inline))
bool qfr_greater(const ssflt xb, const ssint xs, const ssflt yb, const ssint ys)
{
//qfr numbers cannot be negative
#ifdef MPFR_TRUE
	if(yb <= 0.0)
		return true;
	if(xb <= 0.0) //always reject underflow
		return false;
	int s;
	frexp(xb/yb,&s);
	return (s + xs - ys > 0);
#else
	assert(xs == 0 && ys == 0);
	return xb > yb;
#endif	
}

/*
inline __attribute__((always_inline))
void qfr_add(ssflt &xb, ssint &xs, const ssflt yb, const ssint ys)
{
#ifdef MPFR_TRUE	
	
	//bool big_xs( (xs >= ys) && (xb != 0) );
	//xb = big_xs ? ( xb + ldexp(yb, ys - xs) ) : yb + ( (xb==0) ? 0 : ldexp(xb, xs - ys) ) ;
	//xs = big_xs ? xs : ys ;
	if( (xs >= ys) && (xb != 0) )
	{
		xb += ldexp(yb, ys - xs);
	}
	else
	{
		xb = yb + ( (xb==0) ? 0 : ldexp(xb, xs - ys) );
		xs = ys;
	}
#else
	xb += yb;
	assert(xs == 0);
#endif
	
}
*/
#if defined(MACRO_MPFR) && defined(MPFR_TRUE)
#define qfr_add(xb,xs,yb,ys) \
	if( (*(xs) >= (ys)) && (*(xb) != 0) ) \
	{ \
		*(xb) += ldexp((yb), (ys) - *(xs)); \
	} \
	else \
	{ \
		*(xb) = (yb) + ( (*(xb)==0) ? 0 : ldexp(*(xb), *(xs) - (ys)) ); \
		*(xs) = (ys); \
	}

#else //MACRO_MPFR
inline __attribute__((always_inline))
void qfr_add(ssflt *xb, ssint *xs, const ssflt yb, const ssint ys) //pointer version
{
#ifdef MPFR_TRUE	
	
	//bool big_xs( (*xs >= ys) && (*xb != 0) );
	//xb = big_xs ? ( *xb + ldexp(yb, ys - *xs) ) : yb + ( (*xb==0) ? 0 : ldexp(*xb, *xs - ys) ) ;
	//xs = big_xs ? *xs : ys ;
	if( (*xs >= ys) && (*xb != 0) )
	{
		*xb += ldexp(yb, ys - *xs);
	}
	else
	{
		*xb = yb + ( (*xb==0) ? 0 : ldexp(*xb, *xs - ys) );
		*xs = ys;
	}
#else
	*xb += yb;
	assert(*xs == 0);
#endif //MPFR_TRUE
	
}
#endif //MACRO_MPFR

inline __attribute__((always_inline))
void qfr_sub(ssflt *xb, ssint *xs, const ssflt yb, const ssint ys) //pointer version
{
#ifdef MPFR_TRUE	
	
	//bool big_xs( (*xs >= ys) && (*xb != 0) );
	//xb = big_xs ? ( *xb + ldexp(yb, ys - *xs) ) : yb + ( (*xb==0) ? 0 : ldexp(*xb, *xs - ys) ) ;
	//xs = big_xs ? *xs : ys ;
	if( (*xs >= ys) && (*xb != 0) )
	{
		*xb -= ldexp(yb, ys - *xs);
	}
	else
	{
		*xb = -yb + ( (*xb==0) ? 0 : ldexp(*xb, *xs - ys) );
		*xs = ys;
	}
#else
	*xb -= yb;
	assert(*xs == 0);
#endif
	
}
#if 0 //no use of qfr_kahan
inline __attribute__((always_inline))
void qfr_kahan_add(ssflt *xb, ssint *xs, ssflt yb, ssint ys, ssflt *kb, ssint *ks) //pointer version
{
#ifdef MPFR_TRUE	
#if 1
	//simple adding without kahan algorithm
	qfr_add(xb,xs,yb,ys);
#else //1	
	/**	we assume that x and y are non negative, k can have either sign
		if y is bigger than k we do a regular kahan summation step, 
		otherwise we simply reduce k by y and leave the sum x unchanged
		In Principle we should have effectively long double precision  
	**/
	qfr_norm(&yb,&ys);
	if( (*kb < 0.0) || (*ks > ys) || 
		( (*ks == ys) && (yb > *kb) ) )
	{
		//do the kahan summation
		
		//ssflt yb2(yb);
		//ssint ysi(ys);
		if(*kb != 0.0)
		{	
			qfr_sub(&yb,&ys,*kb,*ks);
			//qfr_norm(&yb,&ys);	
		} 
		const ssflt zb(yb);       
		/*if(zb < 0)
		{
			printf("\nERROR z < 0 from z = y-k (with y > k!) \n");
			printf("zb: %lf zs: %d\n", yb, ys);
			printf("xb: %lf xs: %d yb: %lf ys: %d kb: %lf ks: %d\n", *xb, *xs, yb2, ysi, *kb, *ks);
			abort();
		}*/
		
		const ssint zs(ys);       // z = y1 = y0-k
		if(*xb != 0.0)
		{	
			qfr_add(&yb,&ys,*xb,*xs);
			qfr_norm(&yb,&ys);
		}
		
		const ssflt tb(yb);
		const ssflt ts(ys);       // t = y2 = y1+x = y0-k+x 
		if(*xb != 0.0)
		{	
			qfr_sub(&yb,&ys,*xb,*xs);
			qfr_norm(&yb,&ys);
		}
		
		*xb = tb;	*xs = ts;     // x = y3 = y2-x = y1+x-x = y0-k+x-x
		
		if(zb != 0.0)
		{	
			qfr_sub(&yb,&ys,zb,zs);
			qfr_norm(&yb,&ys);
		}
		*kb = yb;	*ks = ys;     // k = y4 = y3 - z = y0-k+x-x-y0+k
		
	}
	else //simply reduce kahan correction
	{
		//if even that does not matter we should be fine anyways
		qfr_sub(kb,ks,yb,ys); 
		qfr_norm(kb,ks);	 
	}	
#endif //1	
#else //USE_MPFR
	const ssflt z( yb - *kb );
	const ssflt t( z + *xb );
	*kb = ( t - *xb ) - z;
	*xb = t;
	
	assert(*ks == 0);
	assert(ys == 0);
	assert(*xs == 0);
#endif
	
}

#endif //no use of qfr_kahan

#if defined(MACRO_MPFR) && defined(MPFR_TRUE)
#define qfr_mul_d(xb,xs,z) \
	int xi; \
	(xb) = frexp((xb)*(z), &xi); \
	(xs) += xi; \

#else //MACRO_MPFR
inline __attribute__((always_inline))
void qfr_mul_d(ssflt &xb, ssint &xs, const ssflt z)
{
#ifdef MPFR_TRUE	
	int xi;
	xb = frexp(xb*z, &xi);
	xs += xi;
#else
	xb *= z;
	assert(xs == 0);
#endif
}
#endif //MACRO_MPFR


inline __attribute__((always_inline))
void qfr_div(ssflt &xb, ssint &xs, const ssflt yb, const ssint ys)
{
#ifdef MPFR_TRUE	
	int tmps;
	xb = frexp(xb/yb, &tmps);
	xs += (tmps - ys);
#else
	xb /= yb;
	assert(xs == 0);
#endif
}

inline __attribute__((always_inline))
void qfr_log2( ssflt &xb, ssint &xs)
{
#ifdef MPFR_TRUE	
	//xb = (xb > 0) ? log2(xb) + xs : log2(-xb) - xs;
	xb = log2(xb) + xs;
	xs = 0;
#else
	xb = log2(xb);
	assert(xs == 0);
#endif
}

inline __attribute__((always_inline))
void qfr_exp2( ssflt &xb, ssint &xs)
{
	//static const double ln2E = log2(exp(1.0));
	long double l2( ldexp((long double)xb,xs) );// * ln2E;
#ifdef MPFR_TRUE		
	xs = ceil(l2);
	xb = exp2(l2-xs);
#else
	xb = l2;
	assert(xs == 0);
#endif
}


#if defined(MACRO_MPFR) && defined(MPFR_TRUE)
#define qfr_exp2Z(xb,xs) \
	ssflt l2((xb)); \
	(xs) = (ssint)ceil((l2)); \
	(xb) = exp2(l2-(xs));

#else //MACRO_MPFR

inline __attribute__((always_inline))
void qfr_exp2Z( ssflt &xb, ssint &xs)
{
	//here xs is assumed to be zero;
#ifdef MPFR_TRUE
	ssflt l2( xb ); //const
	xs = (ssint)ceil(l2);
	xb = exp2(l2-xs);
#else
	xb = exp2(xb);
	assert(xs == 0);
#endif
}
#endif //MACRO_MPFR

#endif //USE_MPFR

inline __attribute__((always_inline))
void qfr_out_str(const ssflt &xb)
{
	printf("PVAL\t%lf\n", xb);
	fflush(stdout);
}

inline __attribute__((always_inline))
void qfr_out_str(const ssflt &xb, const ssint xs)
{
#ifdef MPFR_TRUE
	printf("PVAL\t%lf\n", ldexp(xb,xs));
#else
	printf("PVAL\t%lf\n", xb);
	assert(xs == 0);
#endif
	fflush(stdout);
}













#endif
