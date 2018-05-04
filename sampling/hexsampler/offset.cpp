#include "offset.hpp"
#include "globals.hpp"
#include "basis.hpp"
#include "ellipse.hpp"
#include <math.h>
#include <assert.h>

#define THIRD (1.0/3.0)
#define TWOTHIRD (2.0/3.0)

void update_rsteps()
{
	//assert(uc_stepX >= 0 && uc_stepY >= 0);
	//assert(uc_stepX < uc_res && uc_stepY < uc_res);
	const double da1( (double)uc_stepX/(double)uc_res );
	const double da2( (double)uc_stepY/(double)uc_res );
	const double dx( xd_a(da1,da2) );
	const double dy( yd_a(da1,da2) );
	const double dra1(	ra1d(dx,dy) );
	const double dra2(	ra2d(dx,dy) );
	ruc_stepX = ((int)(dra1*uc_res) + uc_res) % uc_res; 
	ruc_stepY = ((int)(dra2*uc_res) + uc_res) % uc_res;
	//assert(ruc_stepX >= 0 && ruc_stepY >= 0);
	//assert(ruc_stepX < uc_res && ruc_stepY < uc_res);
}

int uc_sum_val(int x, int y)
{
	
	x = (x + uc_stepX) % uc_res;
	y = (y + uc_stepY) % uc_res;
	const int ind( x + uc_res * y );
	return uc_sum[ind];	
}

int ruc_sum_val(int x, int y)
{
	x = (x + ruc_stepX) % uc_res;
	y = (y + ruc_stepY) % uc_res;
	const int ind( x + uc_res * y);
	return ruc_sum[ind];	
}

//map corners on the center
int uc_pos_val(int x, int y)
{
	//assert(uc_stepX != -1 && uc_stepY != -1);
	x = (x + uc_stepX + 3*uc_res/2) % uc_res;
	y = (y + uc_stepY + 3*uc_res/2) % uc_res;
	const int ind( x + uc_res * y );
	return uc_pos[ind];	
}






double mirror(void)
{
	double ss(0.0); //max mirror quality is 1
	//threshold should capture atoms and bonds
	const int treshold( (int)( 0.5 *data_max + 0.5 *data_min ) );
	const int min_val ( (int)( 0.25*data_max + 0.75*data_min) ); 
	int norm(0);
	for(int ind=0; ind < uc_area; ++ind)
	{
		//coordinates in symetrized cell
		int xs = (ind % uc_res); 
		int ys = (ind / uc_res);
		//add up shifted values from sum
		const int val( uc_sum_val(xs,ys) ); //original
		if( val > treshold)
		{
			const int val0( val - min_val );
			const int val1( uc_sum_val(ys,xs) - min_val ); //first mirror
			const int val2( uc_sum_val(uc_res-1-xs,uc_res-1-ys) - min_val ); //second mirror
			const int val3( uc_sum_val(uc_res-1-ys,uc_res-1-xs) - min_val ); //both mirrors;
			const int sum( val0 + val1 + val2 + val3 ); //sum of mirror images
			//every contribution is 0 for totally symetric or 1 for totally antisymetric
			const double ds( ( (double)( abs(sum - 4*val0) + abs(sum - 4*val1) + abs(sum - 4*val2) + abs(sum - 4*val3) )) / ( 6 * sum) );	  
			ss += ( ( (ds > 1.0) || (ds < 0.0) ) ? 1.0 : ds );
			++norm;
		}
	}
	//assert(norm > 0);
	return ( 1.0-ss/(double)norm );
}

double rmirror(void)
{
	double ss(0.0); //max mirror quality is 1
	//threshold should capture atoms and bonds
	const int treshold( (int)( 0.5 *rdata_max + 0.5 *rdata_min ) );
	const int min_val ( (int)( 0.25*rdata_max + 0.75*rdata_min) ); 
	int norm(0);
	for(int ind=0; ind < uc_area; ++ind)
	{
		//coordinates in symetrized cell
		int xs = (ind % uc_res); 
		int ys = (ind / uc_res);
		//add up shifted values from sum
		const int val( ruc_sum_val(xs,ys) ); //original
		if( val > treshold)
		{
			const int val0( val - min_val );
			const int val1( ruc_sum_val(ys,xs) - min_val ); //first mirror
			const int val2( ruc_sum_val(uc_res-1-xs,uc_res-1-ys) - min_val ); //second mirror
			const int val3( ruc_sum_val(uc_res-1-ys,uc_res-1-xs) - min_val ); //both mirrors;
			const int sum( val0 + val1 + val2 + val3 ); //sum of mirror images
			//every contribution is 0 for totally symetric or 1 for totally antisymetric
			const double ds( ( (double)( abs(sum - 4*val0) + abs(sum - 4*val1) + abs(sum - 4*val2) + abs(sum - 4*val3) )) / ( 6 * sum) );	  
			ss += ( ( (ds > 1.0) || (ds < 0.0) ) ? 1.0 : ds );
			++norm;
		}
	}
	//assert(norm > 0);
	return ( 1.0-ss/(double)norm );
}



double gauss_fit( void )
{
	double res(0.0);
	
	for(int ind(0); ind < uc_area; ++ind)
	{
		const int x( ind % uc_res );
		const int y( ind / uc_res );
		const double val( uc_sum_val(x,y) );
		const double gauss( uc_gauss[ind] );
		const double dres( val-gauss );
		res += dres*dres;	
	}
	return sqrt(res)/data_sum;
}

double rgauss_fit( void )
{
	double rres(0.0);
	for(int ind(0); ind < uc_area; ++ind)
	{
		const int x( ind % uc_res );
		const int y( ind / uc_res );
		const double rval( ruc_sum_val(x,y) );
		const double rgauss( ruc_gauss[ind] );
		const double drres( rval-rgauss );
		rres += drres*drres;	
	}
	return sqrt(rres)/rdata_sum;
}


double find_steps()
{
	bool first( true );
	double best( 0.0 );
	double worst( 0.0 );
	int sxM( 0 );
	int syM( 0 );
	double *const mcache( new double[uc_area] );
	memset(mcache, 0, uc_area*sizeof(double) );
	double mnorm(0.0);
	for(int sx(0); sx<uc_res; ++sx )
	{
		uc_stepX = sx;
		for(int sy(0); sy<uc_res; ++sy )
		{
			uc_stepY = sy;
			update_rsteps();
			double m( 1000*gauss_fit() ); //typical range 10^-3
			const double rm( 1000*rgauss_fit() ); //typical range 10^-3
			m+=rm;
			mnorm += m;
			mcache[ sx + sy * uc_res ] = m;
			if ( (m>worst) || first)
			{	worst = m;}
			if( (m < best) || first)
			{
				best = m;
				sxM = sx;
				syM = sy;
			}
			first = false;
		}
	}
	mnorm /= (double)uc_area;
	mnorm -= worst; //is negative
	for(int ind( 0 ); ind<uc_area; ++ind )
	{   //invert and set average to 1000
		uc_pos[ind] = (int)(1000 * (mcache[ind] - worst) / mnorm );
	}
	const int pos_bg_level( 1000 * (0.5*(best+worst) - worst) / mnorm);
	//assert(floating_offset == false);
	if(floating_offset) //auto correct display for offset
	{
		uc_stepX = sxM;
		uc_stepY = syM;
	}
	else //always evaluate and show the unitcells as sampled
	{
		uc_stepX = 0;
		uc_stepY = 0;
	}
	update_rsteps();//actually not yet needed
	next_offset_a1 =  (double)sxM / (double)(uc_res); 
	next_offset_a2 =  (double)syM / (double)(uc_res);
	//the offsets are in ellipsified frame coordinates
	pending_offsets = true; 
	
	double fxM( 0.0 );
	double fyM( 0.0 );
	double xM2( 0.0 );
	double yM2( 0.0 );
	
	double w( 0.0 );
	//double mm( 0.0 );
	const double idelta( 1.0/(worst-best) );
	//int peak_area( 0 );
	for(int sx=-uc_res/2; sx<(uc_res+1)/2; ++sx )
	{
		const int locX( (uc_stepX + sx + uc_res) % uc_res ); //uc_stepX
		
		for(int sy=-uc_res/2; sy<(uc_res+1)/2; ++sy )
		{
			const int locY( (uc_stepY + sy + uc_res) % uc_res ); //uc_stepY
			const double m0 ( mcache[ locX + locY * uc_res] );
			const double mr0( (worst-m0)*idelta );
			if( (mr0 > 0.5) ) //critical for scaling of varX and varY
			{
				const double mr( mr0 - 0.5); //affects std but not the center
				w += mr;
				const double dx( ((double)sx)/((double)uc_res) );
				const double dy( ((double)sy)/((double)uc_res) );
				fxM += dx*mr;
				fyM += dy*mr;
				xM2 += dx*dx*mr;
				yM2 += dy*dy*mr;	 
			}
			else
			{	uc_pos[locX + locY * uc_res] = pos_bg_level;} //nullify all disregarded offsets
		}
	}
	fxM = fxM/w;
	fyM = fyM/w;
	const double varX( fabs( xM2/w -fxM*fxM ) );
	const double varY( fabs( yM2/w - fyM*fyM ) );
	const double area(  30*(varX + varY) ); //scaled area varX and varY can become as small as 0.01 and hence area >~0.55
	distort_basis(true);
	
	delete[] mcache;
	//if(reporting)
	//{	printf("stdX: %lf  stdY: %lf   msym / area: %lf / %lf = %lf pts: %d\n", sqrt(varX), sqrt(varY), msym, area, msym/area, peak_area );}
	mirror_quality = mirror(); //msym
	mirror_quality *= rmirror();
	position_quality = 1.0/area;
	return  ( use_mirror ? (mirror_quality*position_quality) : position_quality );
}

bool apply_a1_a2(void)
{
	if(pending_offsets)
	{
		
		double fx( next_offset_a1 );
		double fy( next_offset_a2 );
		
		//fx -= floor(fx); 
		//fy -= floor(fy);
		if(floating_offset) //dont float the display without floating offset
		{
			uc_stepX = ( (int)(fx*(double)(uc_res) ) + uc_res ) % uc_res;
			uc_stepY = ( (int)(fy*(double)(uc_res) ) + uc_res ) % uc_res;
			update_rsteps();
		}
		
		distort_basis(true);
		const double offx = x_a(fx,fy);
		const double offy = y_a(fx,fy);
		offset_X = offx;
		offset_Y = offy;
		
		//set_offset(offx, offy);
		
		if(reporting)
		{
			printf("a1,a2: %lf,%lf -> X0,Y0:  %lf,%lf\n",next_offset_a1,next_offset_a2,offset_X,offset_Y);
			fflush(stdout);
		}
		
		pending_offsets = false;
		next_offset_a1 = 0.0;
		next_offset_a2 = 0.0;
		return true; 
	}
	return false;
}

