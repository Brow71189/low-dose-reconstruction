#include "globals.hpp"
#include "basis.hpp"
#include "ellipse.hpp"
#include <math.h>
#include <assert.h>

//local use only
void update_idet(void);

void ellipsify_basis(void);
//void correct_basis(void);
double a1x = 0.0;
double a1y = 0.0;
double a2x = 0.0;
double a2y = 0.0;
double ra1x = 0.0;
double ra1y = 0.0;
double ra2x = 0.0;
double ra2y = 0.0;
double p1x = 0.0;
double p1y = 0.0;
double p2x = 0.0;
double p2y = 0.0;

double idet = 1.0;
double ridet = 1.0;
double pidet = 1.0;
bool is_ellipsified = false;


void set_basis(double new_hel, double new_tilt, bool distort)
{
	hel = new_hel;
	tilt = new_tilt;
	
	const double real_tilt( tilt + rotate_tilt*M_PI/3.0 );
	const double len( 3 * (*bondlength) * hel ); 
	//Note the tilt is flipped since screen vs. image coordinates
	//have an inverted y-axis
	a1x = len * cos(-real_tilt-M_PI/6.0);
	a1y = len * sin(-real_tilt-M_PI/6.0);
	a2x = len * cos(-real_tilt+M_PI/6.0);
	a2y = len * sin(-real_tilt+M_PI/6.0);
	ra1x = len * cos(-real_tilt-M_PI/6.0 + M_PI/3.0); //==a2x
	ra1y = len * sin(-real_tilt-M_PI/6.0 + M_PI/3.0); //==a2y
	ra2x = len * cos(-real_tilt+M_PI/6.0 + M_PI/3.0);
	ra2y = len * sin(-real_tilt+M_PI/6.0 + M_PI/3.0);
	
	const double plen( sqrt(3.0) * hel );
	p1x = plen * cos(-real_tilt);
	p1y = plen * sin(-real_tilt);
	p2x = plen * cos(-real_tilt+M_PI/3.0);
	p2y = plen * sin(-real_tilt+M_PI/3.0);
	is_ellipsified = false;
	if(distort)
	{	ellipsify_basis();} //does also update idet
	else
	{	update_idet();}
	set_offset(offset_X, offset_Y); //better to alwasy do that, we might need recentering
	/*//Rectangular cell dor debuging
	a1x = sqrt(3.0)*len;
	a1y = 0;
	a2x = 0;
	a2y = len;
	*/
	
	//if(distort)
	//{	set_offset(offset_X, offset_Y);}//the new basis might require a different periodic centered offset
	/*if(reporting)
	{
		//printf("a1: (%lf,%lf)\ta2: (%lf,%lf)\tlen: %lf\n", a1x, a1y, a2x, a2y, len);
		
		double a1len = a1(impWidth/2,impHeight/2);
		double a2len = a2(impWidth/2,impHeight/2);
		printf("imgcenter: a1=%lf\ta2=%lf\n", a1len, a2len);
		double x = x_a(a1len, a2len);
		double y = y_a(a1len, a2len);
		printf("imgcenter: x=%lf\ty=%lf\n", x, y);
		
		
		x = x_a(floor(a1len), floor(a2len));
		y = y_a(floor(a1len), floor(a2len));
		
		a1len -= floor(a1len);
		a2len -= floor(a2len);
		
		printf("center home: x=%lf\ty=%lf\n",x,y);
		printf("loc vec: a1=%lf\ta2=%lf\n",a1len,a2len);
		
	}
	*/
	 
}

void ellipsify_basis(void)
{
	//assert(is_ellipsified == false);
	ellipsify(a1x, a1y);
	ellipsify(a2x, a2y);
	ellipsify(ra1x, ra1y);
	ellipsify(ra2x, ra2y);
	ellipsify(p1x, p1y);
	ellipsify(p2x, p2y);
	//here idet and ridet may become different
	update_idet();
	is_ellipsified = true;	
}

void distort_basis(bool distort)
{
	if(distort == is_ellipsified)
	{	return;} //easy case nothing to do
	else if(distort)
	{	ellipsify_basis();	} //apply elliptic distorions
	else
	{	set_basis(hel,tilt,false);} //reset basis and skip distortions
}

bool get_basis_distortion(void)
{
	return is_ellipsified;
}

/*void correct_basis(void)
{
	assert(is_ellipsified == true);
	correct_ellipse(a1x, a1y);
	correct_ellipse(a2x, a2y);
	correct_ellipse(ra1x, ra1y);
	correct_ellipse(ra2x, ra2y);
	correct_ellipse(p1x, p1y);
	correct_ellipse(p2x, p2y);
	//here idet and ridet may become different
	update_idet();
	is_ellipsified = false;
	
}*/

void set_uc_spacings(void)
{
	uc_gridX =  (int)(0.5*( fabs(a1x) + fabs(a2x) ));
	uc_gridY =  (int)(0.5*( fabs(a1y) + fabs(a2y) ));
	if(uc_gridX == 0) {	uc_gridX = 1;}
	if(uc_gridY == 0) { uc_gridY = 1;}
	uc_rangeX = impWidth / uc_gridX + 1; //We have to account for incomplete  uc at the edges
	uc_rangeY = impHeight / uc_gridY + 1; //in both directions
	uc_len = uc_rangeX * uc_rangeY;
	
	ruc_gridX =  (int)(0.5*( fabs(ra1x) + fabs(ra2x) ));
	ruc_gridY =  (int)(0.5*( fabs(ra1y) + fabs(ra2y) ));
	if(ruc_gridX == 0) { ruc_gridX = 1;}
	if(ruc_gridY == 0) { ruc_gridY = 1;}
	ruc_rangeX = impWidth / ruc_gridX + 1; //We have to account for incomplete  uc at the edges
	ruc_rangeY = impHeight / ruc_gridY + 1; //in both directions
	ruc_len = ruc_rangeX * ruc_rangeY;
	
	//printf("a1: (%2.2lf,%2.2lf) a2: (%2.2lf,%2.2lf) ra1: (%2.2lf,%2.2lf) ra2: (%2.2lf,%2.2lf) \n", a1x, a1y, a2x, a2y, ra1x, ra1y, ra2x, ra2y );
	//printf("uc_len: %d uc_gridX,Y: %d,%d   ", uc_len, uc_gridX, uc_gridY);
	//printf("ruc_len: %d ruc_gridX,Y: %d,%d\n", ruc_len, ruc_gridX, ruc_gridY);
		
}

//Only rarely needed, dont call during single stepped hex_optimize() as it affects both x and y
void center_offset(void)
{
	double a1len = a1(clean_X,clean_Y);
	double a2len = a2(clean_X,clean_Y);
	a1len = floor(a1len + 0.5);
	a2len = floor(a2len + 0.5);
	double off_x = x_a(a1len, a2len);
	double off_y = y_a(a1len, a2len);
	offset_X = off_x;
	offset_Y = off_y;
	
	//assert( ( (new_dx == 0.0) && (new_dy == 0.0) ) || ( (new_dx != 0.0) && (new_dy != 0.0) ) );
	
	
	if( ( ( a1len != 0.0) || (a2len != 0.0) ) && reporting && has_console)
	{
		printf("recentered offsets %lf,%lf\n", offset_X, offset_Y);
		fflush(stdout);
	}
	 
}

/* //outdated functions
void add_hel(double d_hel)
{
	scale_hel( (hel+d_hel)/hel );
}

void reset_hel(double new_hel)
{
	scale_hel( new_hel/hel );
}

void scale_hel(double r)
{
	a1x *= r;
	a1y *= r;
	a2x *= r;
	a2y *= r;
	hel *= r;
	idet /= r; 
}

void reset_tilt(double new_tilt)
{
	tilt = new_tilt;
	double len = 3 * (*bondlength_CC) * hel;
	//Note the tilt is flipped since screen/image coordinates
	//have an inverted y-axis
	a1x = len * cos(-tilt-M_PI/6);
	a1y = len * sin(-tilt-M_PI/6);
	a2x = len * cos(-tilt+M_PI/6);
	a2y = len * sin(-tilt+M_PI/6);
	update_idet();
	
}

void rotate_basis(double phi)
{
	double new_tilt = tilt + phi;
	reset_tilt(new_tilt);
}
*/

void update_idet(void)
{
	idet =  1.0 / (  a1x *  a2y -  a1y *  a2x );
	ridet = 1.0 / ( ra1x * ra2y - ra1y * ra2x );
	pidet =  1.0 / (  p1x *  p2y -  p1y *  p2x );	
}

void set_offset(double offx, double offy)
{
	offset_X = offx;
	offset_Y = offy;
	//center_offset(); //a1 and a2 originate as close as possiple from clean X/Y
}


/*void shift_offset(double dx, double dy)
{
	dx += offset_X;
	dy += offset_Y;
	//offset_X -= floor(offset_X / abs(a1x + a2x)) * abs(a1x + a2x);
	//offset_Y -= floor(offset_Y / abs(a2y - a1y)) * abs(a2y - a1y);
	double a1len = a1d(dx,dy); 
	double a2len = a2d(dx,dy);
	a1len -= floor(a1len);
	a2len -= floor(a2len);
	offset_X = xd_a(a1len,a2len);
	offset_Y = yd_a(a1len,a2len);
	//	printf("new shifted offsets: %lf,%lf\n", offset_X, offset_Y);	
}
*/
/* These are moved to the header for inlining
double a1(double x, double y)
{
	return ( a2y*(x-offset_X) - a2x*(y-offset_Y) ) * idet;
}

double a2(double x, double y)
{
	return ( a1x*(y-offset_Y) - a1y*(x-offset_X) ) * idet;
}

double a1d(double x, double y)
{
	return ( a2y*x - a2x*y ) * idet;
}

double a2d(double x, double y)
{
	return ( a1x*y - a1y*x ) * idet;
}

double x_a(double a1, double a2)
{
	return offset_X + a1x*a1 + a2x*a2;
}

double y_a(double a1, double a2)
{
	return offset_Y + a1y*a1 + a2y*a2;
}

double xd_a(double a1, double a2)
{
	return a1x*a1 + a2x*a2;
}

double yd_a(double a1, double a2)
{
	return a1y*a1 + a2y*a2;
}
*/

void set_limits(bool hard)
{
	if(hard) //also reset global serach box
	{
		switch( search_type )
		{
			case 0: //original intent scan hel and tilt for better merit
				box_hel_min = hel;
				box_hel_max = hel;
				box_tilt_min = tilt;
				box_tilt_max = tilt;
			break;
			case 1: // abuse of variablenames scan instead excent and phi for better merit
				box_hel_min = excent;
				box_hel_max = excent;
				box_tilt_min = phi;
				box_tilt_max = phi;
			break;
			default:
				printf("Message ERROR unrecognized search type %d in %s line: %d\n",search_type, __FILE__, __LINE__);
				fflush(stdout);
		}	
		
		switch( stability )
		{
			case 0: //thight and dangerous
				box_hel_min *= 0.995; //-0.5%
				box_hel_max *= 1.005; //+0.5%
				box_tilt_min -= M_PI/360; //-0.5°
				box_tilt_max += M_PI/360; //+0.5°
				peak_rank = 2; //start looking for peak right at the beginning with 5x5 grid
			break;
			case 2: //more carefull
				box_hel_min *= 0.96; //-4%
				box_hel_max *= 1.04; //+4%
				box_tilt_min -= M_PI/45; //-4°
				box_tilt_max += M_PI/45; //+4°
				peak_rank = 4; //start looking for peak with 17 x 17 grid
			break;
			case 3: //extensive just to be sure
				box_hel_min *= 0.92; //-8%
				box_hel_max *= 1.08; //+8%
				box_tilt_min -= M_PI/22.5; //-8°
				box_tilt_max += M_PI/22.5; //+8°
				peak_rank = 5; //start looking for peak with 33 x 33 grid
			break;
			default:
				if( reporting)
				{
					printf("Unrecognized stability level %d, defaulting to 1 (reasonable)\n", stability);
					fflush(stdout);
					stability = 1;
				} 
			/*fall through*/
			case 1: //reasonable
				box_hel_min *= 0.98; //-2%
				box_hel_max *= 1.02; //+2%
				box_tilt_min -= M_PI/90; //-2°
				box_tilt_max += M_PI/90; //+2°
				peak_rank = 3; //start looking for peak with 9 x 9 grid
		}
	}
	//reset the peak search area
	
	peak_hel_min = box_hel_min;
	peak_hel_max = box_hel_max;
	box_hel_range = box_hel_max - box_hel_min; 
	peak_hel_range = box_hel_range;
	peak_tilt_min = box_tilt_min;
	peak_tilt_max = box_tilt_max;
	box_tilt_range = box_tilt_max - box_tilt_min; 
	peak_tilt_range = box_tilt_range;
	phint_min = 0;
	phint_max = mg_len - 1;
	ptint_min = 0;
	ptint_max = mg_len - 1;
	
/*		if(reporting)
	{
		printf("hel : min,avg.,max\t%lf,%lf,%lf\n", box_hel_min, hel, box_hel_max);
		printf("tilt: min,avg.,max\t%lf,%lf,%lf\n", box_tilt_min, tilt, box_tilt_max);
	}*/
}

