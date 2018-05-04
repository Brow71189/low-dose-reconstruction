#include "globals.hpp"
#include "ellipse.hpp"
#include "basis.hpp"
#include <math.h>
#include <assert.h>

//local use only
void update_idet_ellipse(void);
void update_ea_eb_phi(void);

double eax = 1.0;
double eay = 0.0;
double ebx = 0.0;
double eby = 1.0;
double idet_ellipse = 1.0;

double stretch_b(1.0);

void update_ea_eb_phi(void)
{
	ea = 2 * impWidth/(3 * hel * (*bondlength) * sqrt(3.0) ); //scaling according to FFT
	eb = ea/excent;

}

void set_ellipse(double new_excent, double new_phi)
{
	phi = new_phi;
	excent = new_excent;
	update_ea_eb_phi();
	stretch_b = 1.0/excent;
	//Note the inverted y-axis
	eax = cos(phi);
	eay = -sin(phi);
	ebx = cos(phi+M_PI_2l); 
	eby = -sin(phi+M_PI_2l);
	
	update_idet_ellipse();
	set_basis(hel, tilt);
	//if(new_hel <= 0.0)
	//{	new_hel = 2 * impWidth/(3 * ea * bondlength_CC * sqrt(3.0) );}
	
	/*
	if(reporting)
	{
		//printf("a1: (%lf,%lf)\ta2: (%lf,%lf)\tlen: %lf\n", a1x, a1y, a2x, a2y, len);
		
		double eac = ead(impWidth/2,impHeight/2);
		double ebc = ebd(impWidth/2,impHeight/2);
		printf("imgcenter: ea=%lf\teb=%lf\n", eac, ebc);
		double x = x_ed(eac, ebc);
		double y = y_ed(eac, ebc);
		printf("imgcenter: x=%lf\ty=%lf\n", x, y);
		correct_ellipse(x, y);
		printf("de-elipsified: x=%lf\ty=%lf\n",x,y);	
	}
	*/
	
}

void update_idet_ellipse(void)
{
	idet_ellipse = 1.0 / ( eax*eby - eay*ebx );// ( a1x*a2y - a1y*a2x );
}
/* These are now inline functions
double ead(double x, double y)
{
	return ( eby*x - ebx*y ) * idet_ellipse;
}

double ebd(double x, double y)
{
	return ( eax*y - eay*x ) * idet_ellipse;
}

double x_ed(double eac, double ebc)
{
	return eax*eac + ebx*ebc;
}

double y_ed(double eac, double ebc)
{
	return eay*eac + eby*ebc;
}

void correct_ellipse(double &x, double &y)
{
	//eac,ebc vector componenents along a and b axis;
	double eac( ead(x,y));
	double ebc( ebd(x,y) * stretch_b);
	//transform back
	x = x_ed(eac,ebc);
	y = y_ed(eac,ebc);
}

void ellipsify(double &x, double &y)
{
	//eac,ebc vector componenents along a and b axis;
	double eac( ead(x,y));
	double ebc( ebd(x,y) * excent);
	//transform back
	x = x_ed(eac,ebc);
	y = y_ed(eac,ebc);
}
 
*/
