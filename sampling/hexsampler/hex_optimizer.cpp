#include "globals.hpp"
#include "hex_optimizer.hpp"
#include "optimizer.hpp"
#include "unitcell.hpp"
#include "hexagons.hpp"
#include "offset.hpp"
#include "basis.hpp"
#include "ellipse.hpp"
#include "communicator.hpp"
#include "xorshift1024star.hpp"
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <signal.h>

//utility for optimize_hexagons
void do_a_step_X(double *h_merit,double *x_now, double *x_step, bool *x_bounded, int *x_counter);
void do_a_step_Y(double *h_merit,double *y_now, double *y_step, bool *y_bounded, int *y_counter);
void do_a_step_T(double *h_merit,double *t_now, double *t_step, bool *t_bounded, int *t_counter);
void do_a_step_H(double *h_merit,double *h_now, double *h_step, bool *h_bounded, int *h_counter);
void do_a_step_P(double *h_merit,double *p_now, double *p_step, bool *p_bounded, int *p_counter);
void do_a_step_E(double *h_merit,double *e_now, double *e_step, bool *e_bounded, int *e_counter);
//void step_error(double *h_merit,double *x_now, double *x_step, bool *x_bounded, int *x_counter);
void (*do_a_step[])(double *h_merit,double *now, double *step, bool *bounded, int *counter) =
{&do_a_step_X,&do_a_step_Y,&do_a_step_T,&do_a_step_H,&do_a_step_P,&do_a_step_E}; //,&step_error

//debug helper
bool check_state(double *const*states);

static const std::string lbls[] = {"X0,Y0","tilt,hel","phi,excent","all","tilt,phi","hel,excent"};

double optimize_hexagons(int lvl, int rank)
{
	if( (lvl < -1) || (lvl > 5) )
	{	
		if( !(reporting && has_console) )
		{	printf("Message ");}
		printf("ERROR unknown lvl in %s at line: %d\n", __FILE__, __LINE__);
		fflush(stdout);
		exit(0);
	}
	//assert(ruling_merit == 1); //can happen if called inside eval_point
	fprintf(logfile,"optimize_hexagons lvl: %d, rank: %d ... ", lvl, rank); fflush(logfile);
	double h_std(0.0);
	//TODO assert( config.exists() );
	//double h_start( config.exists() ? hexagonal_merit : get_hex_merit(offset_X,offset_Y,tilt,hel,phi,excent,h_std) );
	floating_offset = true;
	double h_start( get_hex_merit(offset_X,offset_Y,tilt,hel,phi,excent,h_std,true) );
	
	int trials(0);
	while(apply_p1_p2() && (++trials < 3) )
	{	h_start = get_hex_merit(offset_X,offset_Y,tilt,hel,phi,excent,h_std,true);}

	if(lvl == -1) //no optimization is attempted only p1 and p2 are fixed
	{	
		fprintf(logfile,"done\n"); fflush(logfile);
		return h_start;
	}

	double h_merit ( h_start );
	
	int total_steps( 0 );
	int previous_steps; //( total_steps )
	center_offset();
	center_offset();
	center_offset();
	double offX_now( offset_X );
	double offY_now( offset_Y );
	double tlt_now( tilt );
	double hel_now( hel );
	double phi_now( phi );
	double exc_now( excent ); 
	double *const  now[6] = {&offX_now, &offY_now, &tlt_now, &hel_now, &phi_now, &exc_now};
	//const double start[6] = { offX_now,  offY_now,  tlt_now,  hel_now,  phi_now,  exc_now};
	const double fine( 0.2*pow(0.7,rank) );
	double offX_step( fine * (0.8+0.4*rand_d() ) * hel * ((double)bondlength_CC/(double)bondlength_Mini) );
	double offY_step( fine * (0.8+0.4*rand_d() ) * hel * ((double)bondlength_CC/(double)bondlength_Mini) );
	double tlt_step( fine * (0.08+0.04*rand_d() )  ); //made this 2 times bigger, much better for tiny frames
	double hel_step( fine * (0.08+0.04*rand_d() ) * hel );
	double phi_step( fine * (0.24+0.12*rand_d() ) ); //made this 6 times bigger, much better for tiny frames 
	double exc_step( fine * (0.04+0.02*rand_d() ) ); 
	if(rand_d() > 0.5)
	{	offX_step = -offX_step;}
	if(rand_d() > 0.5)
	{	tlt_step  = -tlt_step;}
	if(rand_d() > 0.5)
	{	phi_step  = -phi_step;}
	if(rand_d() > 0.5)
	{	offY_step = -offY_step;}
	if(rand_d() > 0.5)
	{	hel_step  = -hel_step;}
	if(rand_d() > 0.5)
	{	exc_step  = -exc_step;}
	
	double *const step[6] = {&offX_step, &offY_step, &tlt_step, &hel_step, &phi_step, &exc_step};
	bool offX_bounded( false );
	bool offY_bounded( false );
	bool tlt_bounded( false );
	bool hel_bounded( false );
	bool phi_bounded( false );
	bool exc_bounded( false ); 
    bool *const bounded[6] = {&offX_bounded, &offY_bounded, &tlt_bounded, &hel_bounded, &phi_bounded, &exc_bounded};
    int offX_counter( 0 );
	int offY_counter( 0 );
	int tlt_counter( 0 );
	int hel_counter( 0 );
	int phi_counter( 0 );
	int exc_counter( 0 ); 
    int *const counter[6] = {&offX_counter, &offY_counter, &tlt_counter, &hel_counter, &phi_counter, &exc_counter};
	
	if(reporting)
	{	
		if( (lvl >= 0) && (lvl <= 2) )
		{	printf("%s: %lf, %lf  hex_merit: %lf hex_shape: %lf mom2: %lf\n", lbls[lvl].c_str(), *now[2*lvl], *now[2*lvl+1], h_merit, hex_shape, moment2);}
		else if(lvl == 3)
		{	
			printf("optimizing hexagons hex_merit: %lf hex_shape: %lf mom2: %lf\n", h_merit, hex_shape, moment2);
			printf("X0 : %lf Y0 : %lf tlt: %lf hel: %lf phi: %lf exc: %lf\n",offset_X, offset_Y, tilt, hel, phi, excent);	
		}
		//printf("X0,Y0: %lf,%lf\n", offset_X, offset_Y);
		fflush(stdout);		
	}
	//this clumsy switch silences scan-build
	const int i_max( 6 );
	bool do_is[i_max] = {false, false, false, false, false, false}; //do nothing
	switch(lvl)
	{
		case 0:
			do_is[0] = true;
			do_is[1] = true;
			break;
		case 1:
			do_is[2] = true;
			do_is[3] = true;
			break;
		case 2:
			do_is[4] = true;
			do_is[5] = true;
			break;
		case 3: //do all is
			for(int i(0); i < i_max; ++i)
			{	do_is[i] = true;}
			
			break;
		case 4:
			do_is[2] = true;
			do_is[4] = true;
			if(reporting)
			{	
				printf("%s: %lf, %lf  hex_merit: %lf hex_shape: %lf mom2: %lf\n", lbls[lvl].c_str(), *now[2], *now[4], h_merit, hex_shape, moment2);
				fflush(stdout);
			}
			break;
		case 5:
			do_is[3] = true;
			do_is[5] = true;
			if(reporting)
			{	
				printf("%s: %lf, %lf  hex_merit: %lf hex_shape: %lf mom2: %lf\n", lbls[lvl].c_str(), *now[3], *now[5], h_merit, hex_shape, moment2);
				fflush(stdout);	
			}
			break;	
		default:
		if( !(reporting && has_console) )
		{	printf("Message ");}
		printf("ERROR unknown lvl in %s at line: %d\n", __FILE__, __LINE__);
		fflush(stdout);
		exit(0);		
	}
	
	
	do
	{	
		previous_steps = total_steps;
		int j(0);
		double activity(0.0); //measured in frame pixels
		double active(0.0);
		for(int i(0); i < i_max; ++i)
		{
			if(!do_is[i])
			{	continue;} //skip them
			switch(i) //cases for measuring the step
			{
				case 0: //offX
				case 1: //offY
					active = *step[i];
					break;
				case 2: //tilt
				case 5: //excent
					active = *step[i] * 10 * ea; 
					break;	
				case 3: //hel
					active = *step[i] * 10 * ea / hel; 
					break;	
				case 4: //phi
					active = *step[i] * 10 * ea * (ea - eb);
					break;
				default:
				if( !(reporting && has_console) )
				{	printf("Message ");}
				printf("ERROR unrecognized optimization parameter %d in %s line: %d\n",search_type, __FILE__, __LINE__);
				fflush(stdout);
				exit(0);
			}
			active = fabs(active);
			if( 	(active > activity) || 
					( (active == activity) && (rand_d() > 0.5) )
			  )
			{
				activity = active;
				j = i;
			}
		}
		//printf("activity: %lf j: %d\n",activity, j);
		if( (activity > ( hel/(double)impWidth)) /*&& (j >= i_min) && (j < i_max)*/ )  //YES such fine steps make a difference
		{
			do_a_step[j]( &h_merit, now[j], step[j], bounded[j], counter[j] );
			//these may have been changed by center_offset()
			//offX_now = offset_X;
			//offY_now = offset_Y;
			
			++total_steps;
			 
			
			//check_state(now);
			/*
			if(!check_state(now))
			{
				printf("Incoherrent bounded state detected with j: %d\n", j);
			}*/
		}
		//}
	}
	while( (total_steps > previous_steps) && (total_steps < 120) );
	
	
	//this section should totally ensure a coherent but not necessarily correct state,
	//irrespective of the order or acceptance of the last trials
	/*
	offset_X = offX_now;
	offset_Y = offY_now;
	offset_XS = offX_now;
	offset_YS = offY_now;
	tilt = tlt_now;
	hel = hel_now;
	phi = phi_now;
	excent = exc_now;		
	set_ellipse(excent,phi);			
	*/		
	if(reporting)
	{	
		if(lvl >= 0 && lvl <= 2)
		{	printf("%s: %lf, %lf  hex_merit: %lf hex_shape: %lf mom2: %lf trials: %d x %lu\n", lbls[lvl].c_str(), *now[2*lvl],*now[2*lvl+1],h_merit,hex_shape,moment2,total_steps,top_points);}
		else if (lvl == 3)
		{	
			printf("X0 : %lf Y0 : %lf tlt: %lf hel: %lf phi: %lf exc: %lf\n",offset_X, offset_Y, tilt, hel, phi, excent);	
			printf("optimized hexagons hex_merit: %lf hex_shape: %lf mom2: %lf trials: %d x %lu\n", h_merit, hex_shape, moment2, total_steps,top_points);
		}
		else if (lvl == 4)
		{	printf("%s: %lf, %lf  hex_merit: %lf hex_shape: %lf mom2: %lf trials: %d x %lu\n", lbls[lvl].c_str(), *now[2],*now[4],h_merit,hex_shape,moment2,total_steps,top_points);}
		else if (lvl == 5)
		{	printf("%s: %lf, %lf  hex_merit: %lf hex_shape: %lf mom2: %lf trials: %d x %lu\n", lbls[lvl].c_str(), *now[3],*now[5],h_merit,hex_shape,moment2,total_steps,top_points);}
		fflush(stdout);		
	}
	fprintf(logfile,"done\n"); fflush(logfile);
	//check_memory();
	return (hexagonal_merit = h_merit);	
}


double get_hex_merit(double offx, double offy, double t, double h, double p, double e, double &std, bool quick)
{
	//assert(ruling_merit == 1);
	/*******SWITCHING******/
	/*
	bondlength = &bondlength_Mini;
	const double bl_scaling(((double)bondlength_CC)/((double)bondlength_Mini));
	const int original_offset_Q( offset_Q );
	const int original_offset_R( offset_R );
	const int original_modelsize( modelsize );
	const int original_hexImpWidth( hexImpWidth );
	const int original_hexImpHeight( hexImpHeight );
	const int original_coh_len( coh_len );
	
	hexImpWidth = (int) (((double)hexImpWidth) / bl_scaling);
	hexImpHeight = (int)(((double)hexImpHeight) / bl_scaling);
	coh_len = (int) ( coh_len / bl_scaling);
	if(hexImpWidth % 2 == 1) { ++hexImpWidth; }
	if(hexImpHeight % 2 == 1) { ++hexImpHeight; }
	modelsize = 2 * bondlength_Mini;
	offset_Q = 0;
	offset_R = 0;
	*/ 
	set_large_hexes();
	/***********************/
	const double bl_scaling(((double)bondlength_CC)/((double)bondlength_Mini));
	h *= bl_scaling;
	/**********************/
	
	offset_X = offx;
	offset_Y = offy;
	tilt = t;	
	hel = h; //these changes are not yet effective
	set_ellipse(e,p); //updates all vector stuff
	
	double avg(0.0);
	double avg_hsp(0.0);
	double avg_mom2(0.0);
	double avg2(0.0);
	const int last_pt ((int)top_points);
	int pts(0);
	
	const bool old_floating_offset( floating_offset );
	 //determine p1 and p2 only once at the end, otherwise the first floating merit might outperform all others
	do
	{
		floating_offset = old_floating_offset && (quick || (pts == 0));
		double hsp(0.0);
		const double val( sample_hexagons(hsp,false) );
		avg_hsp += hsp;
		avg_mom2 += moment2;
		avg += val;
		avg2 += val*val;
		++pts;
		 
		if( (val < 0.0) || (val < 0.9 * hexagonal_merit) )
		{	
			cancel_p1_p2();
			break;
		}
	}
	while( (!quick) && (pts < last_pt) );
	floating_offset = old_floating_offset;
	avg /= (double)pts;
	std = sqrt( fabs(avg2/(double)pts - avg*avg) );
	avg_hsp /= (double)pts;
	hex_shape = avg_hsp;
	avg_mom2 /= (double)pts;
	moment2 = avg_mom2;
	
	/*******SWITCHING_BACK*/
	/*
	bondlength = &bondlength_CC;
	
	hexImpWidth = original_hexImpWidth;
	hexImpHeight = original_hexImpHeight;
	coh_len = original_coh_len;
	offset_Q = original_offset_Q;
	offset_R = original_offset_R;
	modelsize = original_modelsize;
	*/
	set_small_hexes();
	/***********************/
	hel /= bl_scaling;
	/***********************/
	return avg;
}

void do_a_step_X(double *h_merit,double *x_now, double *x_step, bool *x_bounded, int *counter)
{	
	double hm_std(0.0);
	const double x_hm( get_hex_merit(*x_now + *x_step,offset_Y,tilt,hel,phi,excent,hm_std));
	if(x_hm <= *h_merit)
	{
		//stay switch direction and decrease stepsize unless at startup
		cancel_p1_p2();
		*x_bounded |= (*counter > 0);
		*x_step = - *x_step;
		set_offset(*x_now, offset_Y);
		if(*x_bounded) {*x_step *= 0.5;}
	}
	else
	{
		apply_p1_p2();
		*h_merit = x_hm;
		*x_now = offset_X;
		if(reporting && has_console)
		{
			printf("new X0 : %lf(%lf) hex_merit: %lf(%lf) hex_shape: %lf mom2: %lf\n",*x_now, *x_step, *h_merit, hm_std, hex_shape, moment2);
			fflush(stdout);
		}
		if(!*x_bounded && (*x_step < hel))
		{	*x_step *= 1.5;}
	}
	++(*counter);
	//assert( p1_and_p2_applied );
}

void do_a_step_Y(double *h_merit, double *y_now, double *y_step, bool *y_bounded, int *counter)
{	
	double hm_std(0.0);
	const double y_hm( get_hex_merit(offset_X, *y_now + *y_step,tilt,hel,phi,excent,hm_std));
	if(y_hm <= *h_merit) 
	{
		cancel_p1_p2();
		*y_bounded |= (*counter > 0);
		*y_step = -*y_step; 
		if(*y_bounded) {	*y_step *= 0.5;}
		set_offset(offset_X, *y_now);
	}
	else
	{
		apply_p1_p2();
		*h_merit = y_hm;
		*y_now = offset_Y;
		if(reporting && has_console)
		{
			printf("new Y0 : %lf(%lf) hex_merit: %lf(%lf) hex_shape: %lf mom2: %lf\n", *y_now, *y_step, *h_merit, hm_std, hex_shape, moment2);
			fflush(stdout);
		}
		if(!*y_bounded && (*y_step < hel))
		{	*y_step *= 1.5;}
	}
	++(*counter);
	//assert( p1_and_p2_applied );
}

void do_a_step_T(double *h_merit,double *t_now, double *t_step, bool *t_bounded, int *counter)
{	
	double hm_std(0.0);
	const double t_hm( get_hex_merit(offset_X,offset_Y,*t_now + *t_step,hel,phi,excent,hm_std));
	if(t_hm <= *h_merit)
	{
		//stay switch direction and decrease stepsize unless at startup
		cancel_p1_p2();
		*t_bounded |= (*counter > 0);
		*t_step = -*t_step;
		set_basis(hel, *t_now);
		if(*t_bounded) {*t_step *= 0.5;}
	}
	else
	{
		apply_p1_p2();
		*h_merit = t_hm;
		*t_now = tilt; 
		if(reporting && has_console)
		{
			printf("new tlt: %lf(%lf) hex_merit: %lf(%lf) hex_shape: %lf mom2: %lf\n",*t_now, *t_step, *h_merit, hm_std, hex_shape, moment2);
			fflush(stdout);
		}
		if(!*t_bounded && (*t_step < 0.005)) // tilt is limited to < 1°
		{	*t_step *= 1.5;}
	}
	++(*counter);
	//assert( p1_and_p2_applied );
}

void do_a_step_H(double *h_merit, double *h_now, double *h_step, bool *h_bounded, int *counter)
{	
	double hm_std(0.0);
	const double h_hm( get_hex_merit(offset_X,offset_Y,tilt,*h_now+*h_step,phi,excent,hm_std));
	if(h_hm <= *h_merit)
	{
		//stay switch direction and decrease stepsize unless at startup
		cancel_p1_p2();
		*h_bounded |= (*counter > 0);
		*h_step = -*h_step;
		set_basis(*h_now, tilt);
		if(*h_bounded) {*h_step *= 0.5;}
	}
	else
	{
		apply_p1_p2();
		*h_merit = h_hm;
		*h_now = hel;
		if(reporting && has_console)
		{
			printf("new hel: %lf(%lf) hex_merit: %lf(%lf) hex_shape: %lf mom2: %lf\n",*h_now, *h_step, *h_merit, hm_std, hex_shape, moment2);
			fflush(stdout);
		}
		if(!*h_bounded && (*h_step < 0.01*hel)) //step in hel is limited to 1.5% 
		{	*h_step *= 1.5;}
	}
	++(*counter);
	//assert( p1_and_p2_applied );
}

void do_a_step_P(double *h_merit,double *p_now, double *p_step, bool *p_bounded, int *counter)
{	
	double hm_std(0.0);
	const double p_hm( get_hex_merit(offset_X,offset_Y,tilt,hel,*p_now+*p_step,excent,hm_std));
	if(p_hm <= *h_merit)
	{
		//stay switch direction and decrease stepsize unless at startup
		cancel_p1_p2();
		*p_bounded |= (*counter > 0);
		*p_step = -*p_step;
		set_ellipse(excent, *p_now);
		if(*p_bounded) {*p_step *= 0.5;}
	}
	else
	{
		apply_p1_p2();
		*h_merit = p_hm;
		*p_now = phi;
		if(reporting && has_console)
		{
			printf("new phi: %lf(%lf) hex_merit: %lf(%lf) hex_shape: %lf mom2: %lf\n",*p_now, *p_step, *h_merit, hm_std, hex_shape, moment2);
			fflush(stdout);
		}
		if(!*p_bounded && (*p_step < 0.005)) // tilt is limited to < 1°
		{	*p_step *= 1.5;}
	}
	++(*counter);
	//assert( p1_and_p2_applied );
}

void do_a_step_E(double *h_merit,double *e_now, double *e_step, bool *e_bounded, int *counter)
{	
	double hm_std(0.0);
	const double e_hm( get_hex_merit(offset_X,offset_Y,tilt,hel,phi,*e_now + *e_step,hm_std));
	if(e_hm <= *h_merit)
	{
		//stay switch direction and decrease stepsize unless at startup
		cancel_p1_p2();
		*e_bounded |= (*counter > 0);
		*e_step = -*e_step;
		set_ellipse(*e_now, phi);
		if(*e_bounded) {*e_step *= 0.5;}
	}
	else
	{
		apply_p1_p2();
		*h_merit = e_hm;
		*e_now = excent;
		if(reporting && has_console)
		{
			printf("new exc: %lf(%lf) hex_merit: %lf(%lf) hex_shape: %lf mom2: %lf\n",*e_now, *e_step, *h_merit, hm_std, hex_shape, moment2);
			fflush(stdout);
		}
		if(!*e_bounded && (*e_step < 0.01)) // excentricity is sensitive
		{	*e_step *= 1.5;}
	}
	++(*counter);
	//assert( p1_and_p2_applied );
}
/*
void step_error(double *h_merit,double *e_now, double *e_step, bool *e_bounded, int *counter)
{
	printf("OUCH somebody called step_error\n");
	fflush(stdout);
	raise(SIGABRT);
}
*/

bool check_state(double *const*states)
{
	bool ok( true );
	 //these are not strict because of periodic boundary conitions
	if (*states[0] != offset_X)
	{
		printf("error X0 : %lf != offset_X: %lf\n", *states[0], offset_X);
		ok = false;
	}
	if (*states[1] != offset_Y)
	{
		printf("error Y0 : %lf != offset_Y: %lf\n", *states[1], offset_Y);
		ok = false;
	}
	if (*states[2] != tilt)
	{
		printf("error tlt: %lf != tilt: %lf\n", *states[2], tilt);
		ok = false;
	}
	if (*states[3] != hel)
	{
		printf("error hel: %lf != hexedgelen: %lf\n", *states[3], hel);
		ok = false;
	}
	if (*states[4] != phi)
	{
		printf("error phi: %lf != ellipse phi: %lf\n", *states[4], phi);
		ok = false;
	}
	if (*states[5] != excent)
	{
		printf("error exc: %lf != excent: %lf\n", *states[5], excent);
		ok = false;
	}
	if(!ok)
	{	fflush(stdout);}
	
	return ok;
	
}

//TRASH SECTION
		/*
		//if( (h_merit < 0.0) && (total_steps == 0) )
		if( false && ( (next_offset_p1 != 0.0) || (next_offset_p2 != 0.0) || (h_merit < 0.0) ) )
		{
			int search_count( 0 );
			// rmax is 1.5 C-C bondlength
			const double xda( 0.5 * xd_a(1.0,1.0) );
			const double yda( 0.5 * yd_a(1.0,1.0) );
			
			const double rmax = sqrt( xda*xda + yda*yda );
			double x0( (0.5 * (double)impWidth -  rmax) + rmax*rand_d() );
			double y0( (0.5 * (double)impHeight - rmax) + rmax*rand_d() );
			
			if(config.exists()) //That should always be true but we can afford saftey here
			{
				if(!p1_and_p2_applied) //p1 and p2 are expensive and only created once at start up
				{
					x0 = x_p( next_offset_p1, next_offset_p2);
					y0 = y_p( next_offset_p1, next_offset_p2);
					if(reporting)
					{
						printf("p1,p2: %.0lf,%.0lf -> X0,Y0:  %lf,%lf\n",next_offset_p1,next_offset_p2,x0,y0);
						fflush(stdout);
					}
					
					next_offset_p1 = 0.0;
					next_offset_p2 = 0.0;
					p1_and_p2_applied = true;
				}
				else if(pending_offsets) // use cheap a1 and a2 for later refinements
				{
					x0 = x_a( next_offset_a1, next_offset_a2);
					y0 = y_a( next_offset_a1, next_offset_a2);
					if(reporting)
					{
						printf("a1,a2: %lf,%lf -> X0,Y0:  %lf,%lf\n",next_offset_a1,next_offset_a2,x0,y0);
						fflush(stdout);
					}
					next_offset_a1 = 0.0;
					next_offset_a2 = 0.0;
					pending_offsets = false;
				}
			}
			
			h_merit = get_hex_merit(x0, y0, tilt, hel, phi, excent, h_std,true);
			if(h_merit < 0)
			{	
				const double r0 = 0.15 * rmax;
				double radius = r0;
				const double a0 = 2*M_PI * rand_d();
				double alpha = a0;
				
				double best_h_merit(h_merit);
				double best_offset_X(offset_X);
				double best_offset_Y(offset_Y);
				do
				{
					int N_pts = (int) (2*M_PI * radius / r0);
					double d_alpha = (2*M_PI / (double) N_pts);
					alpha += 0.5 * d_alpha; //avoid hitting 2Pi exactlly
					do
					{
						++search_count;
						
						h_merit = get_hex_merit(x0 + radius * cos(alpha),y0 + radius * sin(alpha),tilt,hel,phi,excent,h_std,true);
						if(h_merit > best_h_merit)
						{	
							best_h_merit = h_merit;
							best_offset_X = offset_X;
							best_offset_Y = offset_Y;
							offX_now = offset_X;
							offY_now = offset_Y;
							if(reporting && has_console)
							{
								printf("new X0,Y0: %lf %lf hex_merit: %lf hex_shape: %lf\n", best_offset_X, best_offset_Y, h_merit, hex_shape);
							}
						}
						else //that should not be needed
						{
							h_merit  = best_h_merit;
							offX_now = best_offset_X;
							offY_now = best_offset_Y;
						}
						//if(reporting && has_console)
						//{
						//	printf("circular search r: %lf a: %lf x: %lf y: %lf h_merit: %lf\n", radius, alpha, offset_X, offset_Y, h_merit);
						//	fflush(stdout);
						//}
						alpha += d_alpha * (0.9 + 0.2 * rand_d() );
					}
					while( (alpha < (a0 + 2*M_PI) ) && (h_merit < 0.0) );
					alpha -= (2*M_PI);
					radius += r0 * (0.9 + 0.2 * rand_d()); 
				}
				//better to waste a little time on trying realy hard here
				//non hexagonal merit will take for about one hour
				while( (radius < 3*rmax) && (h_merit < 0.0) );		
				if(h_merit > 0)
				{
					if(reporting)
					{	printf("circular search on frame %d succeded in %d trials hex_merit: %lf\n",frameNr ,search_count, h_merit);}
				}
				else
				{
					printf("Message WARNING frame: %d circular search gave up after %d trials\n", frameNr, search_count);
				}
			}
			else
			{
				offX_now = offset_X;
				offY_now = offset_Y;
				//if(reporting && has_console)
				//{
				//	printf("lucky X0,Y0: %lf %lf hex_merit: %lf hex_shape: %lf\n", offset_X,offset_Y, h_merit, hex_shape);
				//}
			}
		}
		p1_and_p2_applied = true; //they will be outdated by optimization anyways
		*/
