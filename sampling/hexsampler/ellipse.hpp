#ifndef ELLIPSE_HPP
#define ELLIPSE_HPP

#include "globals.hpp"

extern double eax;
extern double eay;
extern double ebx;
extern double eby;
extern double idet_ellipse;

void set_ellipse(double new_excent, double new_phi);

extern double stretch_b;

//sets ea, eb and phi in proper format for communication with master
//double ead(double x, double y);
inline __attribute__((always_inline))
double ead(const double x, const double y)
{
	return ( eby*x - ebx*y ) * idet_ellipse;
}
//double ebd(double x, double y);
inline __attribute__((always_inline))
double ebd(const double x, const double y)
{
	return ( eax*y - eay*x ) * idet_ellipse;
}
//double x_ed(double eac, double ebc);
inline __attribute__((always_inline))
double x_ed(const double eac, const double ebc)
{
	return eax*eac + ebx*ebc;
}
//double y_ed(double eac, double ebc);
inline __attribute__((always_inline))
double y_ed(const double eac, const double ebc)
{
	return eay*eac + eby*ebc;
}



// x,y image coordinates cx,cy deellipsoized coordinates
//void correct_ellipse(double &x, double &y);
inline __attribute__((always_inline))
void correct_ellipse(double &x, double &y)
{
	//eac,ebc vector componenents along a and b axis;
	const double eac( ead(x,y));
	const double ebc( ebd(x,y) * stretch_b);
	//transform back
	x = x_ed(eac,ebc);
	y = y_ed(eac,ebc);
}
//void ellipsify(double &x, double &y)
inline __attribute__((always_inline))
void ellipsify(double &x, double &y)
{
	//eac,ebc vector componenents along a and b axis;
	const double eac( ead(x,y));
	const double ebc( ebd(x,y) * excent);
	//transform back
	x = x_ed(eac,ebc);
	y = y_ed(eac,ebc);
}


#endif
