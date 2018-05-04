#ifndef BASIS_HPP
#define BASIS_HPP

#include "globals.hpp"
//needed for inlining
extern double ra1x,a1x;
extern double ra1y,a1y;
extern double ra2x,a2x;
extern double ra2y,a2y;
extern double p1x,p1y,p2x,p2y;
extern double ridet,idet,pidet;


void set_basis(double new_hel, double new_tilt, bool distort = true);


void set_uc_spacings(void);

//applies elliptic distortion to basis, this is the default state
//void ellipsify_basis(void);
//undoes elliptic distortions, has to be called manually and should be undone by calling ellipsify_basis after use
//this is useful for caculating undistorted euclidian distances
//void correct_basis(void);

//sets the basis vectors to match neighboring pixels

void distort_basis(bool distort = true);
//DEBUG helper
bool get_basis_distortion(void);

void set_limits(bool hard);

//void add_hel(double d_hel);

//void reset_hel(double new_hel);

//void reset_tilt(double new_tilt);

//void rotate_basis(double phi);

void center_offset(void);

void set_offset(double new_dx, double new_dy);

//void shift_offset(double dx, double dy);

//double a1(double x, double y);
//No difference to inline functions as long as basis vectors are still global variables
//#define MACRO_BASIS 

#ifdef MACRO_BASIS

#define a1(x,y) ( ( a2y*(x-offset_X) - a2x*(y-offset_Y) ) * idet )
#define p1(x,y) ( ( p2y*(x-offset_X) - p2x*(y-offset_Y) ) * pidet )
#define ra1(x,y) ( ( ra2y*(x-offset_X) - ra2x*(y-offset_Y) ) * ridet )
#define a2(x,y) ( ( a1x*(y-offset_Y) - a1y*(x-offset_X) ) * idet )
#define p2(x,y) ( ( p1x*(y-offset_Y) - p1y*(x-offset_X) ) * pidet )
#define ra2(x,y) ( ( ra1x*(y-offset_Y) - ra1y*(x-offset_X) ) * ridet )
#define a1d(x,y) ( ( a2y*x - a2x*y ) * idet )
#define p1d(x,y) ( ( p2y*x - p2x*y ) * pidet )
#define ra1d(x,y) ( ( ra2y*x - ra2x*y ) * ridet )
#define a2d(x,y) ( ( a1x*y - a1y*x ) * idet )
#define p2d(x,y) ( ( p1x*y - p1y*x ) * pidet )
#define ra2d(x,y) ( ( ra1x*y - ra1y*x ) * ridet )
#define x_a(a1,a2) ( offset_X + a1x*a1 + a2x*a2 )
#define x_p(p1,p2) ( offset_X + p1x*p1 + p2x*p2 )
#define x_ra(ra1,ra2) ( offset_X + ra1x*ra1 + ra2x*ra2 )
#define y_a(a1,a2) ( offset_Y + a1y*a1 + a2y*a2 )
#define y_p(p1,p2) ( offset_Y + p1y*p1 + p2y*p2 )
#define y_ra(ra1,ra2) ( offset_Y + ra1y*ra1 + ra2y*ra2 )
#define xd_a(a1,a2) ( a1x*a1 + a2x*a2 )
#define xd_p(p1,p2) ( p1x*p1 + p2x*p2 )
#define xd_ra(ra1,ra2) ( ra1x*ra1 + ra2x*ra2 )
#define yd_a(a1,a2) ( a1y*a1 + a2y*a2 )
#define yd_p(p1,p2) ( p1y*p1 + p2y*p2 )
#define yd_ra(ra1,ra2) ( ra1y*ra1 + ra2y*ra2 )

#else

inline __attribute__((always_inline))
double a1(double x, double y)
{
	return ( a2y*(x-offset_X) - a2x*(y-offset_Y) ) * idet;
}


inline __attribute__((always_inline))
double p1(double x, double y)
{
	return ( p2y*(x-offset_X) - p2x*(y-offset_Y) ) * pidet;
}


inline __attribute__((always_inline))
double ra1(double x, double y)
{
	return ( ra2y*(x-offset_X) - ra2x*(y-offset_Y) ) * ridet;
}


inline __attribute__((always_inline))
double a2(double x, double y)
{
	return ( a1x*(y-offset_Y) - a1y*(x-offset_X) ) * idet;
}


inline __attribute__((always_inline))
double p2(double x, double y)
{
	return ( p1x*(y-offset_Y) - p1y*(x-offset_X) ) * pidet;
}


inline __attribute__((always_inline))
double ra2(double x, double y)
{
	return ( ra1x*(y-offset_Y) - ra1y*(x-offset_X) ) * ridet;
}


inline __attribute__((always_inline))
double a1d(double x, double y)
{
	return ( a2y*x - a2x*y ) * idet;
}


inline __attribute__((always_inline))
double p1d(double x, double y)
{
	return ( p2y*x - p2x*y ) * pidet;
}


inline __attribute__((always_inline))
double ra1d(double x, double y)
{
	return ( ra2y*x - ra2x*y ) * ridet;
}


inline __attribute__((always_inline))
double a2d(double x, double y)
{
	return ( a1x*y - a1y*x ) * idet;
}


inline __attribute__((always_inline))
double p2d(double x, double y)
{
	return ( p1x*y - p1y*x ) * pidet;
}


inline __attribute__((always_inline))
double ra2d(double x, double y)
{
	return ( ra1x*y - ra1y*x ) * ridet;
}


inline __attribute__((always_inline))
double x_a(double a1, double a2)
{
	return offset_X + a1x*a1 + a2x*a2;
}


inline __attribute__((always_inline))
double x_p(double p1, double p2)
{
	return offset_X + p1x*p1 + p2x*p2;
}


inline __attribute__((always_inline))
double x_ra(double ra1, double ra2)
{
	return offset_X + ra1x*ra1 + ra2x*ra2;
}


inline __attribute__((always_inline))
double y_a(double a1, double a2)
{
	return offset_Y + a1y*a1 + a2y*a2;
}


inline __attribute__((always_inline))
double y_p(double p1, double p2)
{
	return offset_Y + p1y*p1 + p2y*p2;
}


inline __attribute__((always_inline))
double y_ra(double ra1, double ra2)
{
	return offset_Y + ra1y*ra1 + ra2y*ra2;
}


inline __attribute__((always_inline))
double xd_a(double a1, double a2)
{
	return a1x*a1 + a2x*a2;
}


inline __attribute__((always_inline))
double xd_p(double p1, double p2)
{
	return p1x*p1 + p2x*p2;
}


inline __attribute__((always_inline))
double xd_ra(double ra1, double ra2)
{
	return ra1x*ra1 + ra2x*ra2;
}


inline __attribute__((always_inline))
double yd_a(double a1, double a2)
{
	return a1y*a1 + a2y*a2;
}


inline __attribute__((always_inline))
double yd_p(double p1, double p2)
{
	return p1y*p1 + p2y*p2;
}


inline __attribute__((always_inline))
double yd_ra(double ra1, double ra2)
{
	return ra1y*ra1 + ra2y*ra2;
}

#endif //NO MACRO_BASIS

#endif




