#include "globals.hpp"
#include "optimizer.hpp"
#include "hex_optimizer.hpp"
#include "unitcell.hpp"
#include "hexagons.hpp"
#include "offset.hpp"
#include "basis.hpp"
#include "ellipse.hpp"
#include "communicator.hpp"
#include <math.h>
#include <algorithm>
#include <assert.h>

/* private headers  */
//double sampled_mean(double &variance);
//double sampled_mean5(double &variance);

//updates the config creates an initial config if it does not exist
//force will always acceppt the first new_merit
int step_search(int &rank, const bool force_first, const bool do_first = true);

double grid_search(int &rank);



//returns true if the averaged coords are accepted and applied 
//bool peak_stats(int &rank, int &peak_sample, double &var_t, double &var_h);
//utility for grid_search
void eval_point(int h_ind, int t_ind, double &bh, double &bt, const int rank);

//utility for eval_point
bool comp_merit_ind(int i, int j);

/*double sampled_mean(double &variance)
{
	//not supposed to be used without any sampled merits
	if(sampled_merits.empty())
	{
		variance = 0.0;
		return -1.0; 	
	}
	double sum (0.0);
	double sum2 (0.0);
	for(std::map<int,double>::iterator it = sampled_merits.begin();
		it != sampled_merits.end(); ++it)
	{
		const double val(it->second);
		sum += val;
		sum2 += val*val;
	}
	const int size( sampled_merits.size() );
	const double mean(sum/size);
	variance = fabs( sum2/size - mean*mean );
	return mean;	
}*/

/*double sampled_mean5(double &variance)
{
	//not supposed to be used without any top5
	if(top5.empty())
	{
		variance = 0.0;
		return -1.0; 	
	}
	double sum (0.0);
	double sum2 (0.0);
	const int size( top5.size() );
	for(int i(0); i < size; ++i)
	{
		const double val( sampled_merits[ top5[i] ] );
		sum += val;
		sum2 += pow(val, 2);
	}
	const double mean( sum/size );
	variance = fabs( sum2/size - pow(mean,2) );
	return mean;	
}*/

bool comp_merit_ind(int i, int j)
{
	return sampled_merits[i] > sampled_merits[j];
}

void eval_point(int h_ind, int t_ind, const int rank)
{
	//const bool last_fl_off(floating_offset);
	floating_offset = (stability > 0);
	//cancel_p1_p2();
	const int ind( t_ind + h_ind * mg_len );
	if(sampled_merits.count(ind) != 0)
	{	return; } //silently skip known points
	const double h( box_hel_min  + h_ind * box_hel_range  / (mg_len-1) );
	const double t( box_tilt_min + t_ind * box_tilt_range / (mg_len-1) );
	switch(search_type)
	{
		case 0: //hel and tilt
			set_basis(h, t);
		break;
		case 1: //excent and phi
			set_ellipse(h, t);
		break;
		default:
			if(reporting && !has_console)
			{	printf("Message ");} 
			printf("ERROR unrecognized search type %d in %s line: %d\n",search_type, __FILE__, __LINE__);
			fflush(stdout);
			exit(0);
	}
	
	
	
	double any_std(-1.0);
	double new_merit( mean_merit(any_std,true) ); //quick test
	//cancel_p1_p2();//simply never use roughly tested p1 and p2, 
	//the point should already be preoptimzed by step_search
	
	///THIS CAN GO WRONG WITH AN ACTIVE UC_TARGET
	/*if(ruling_merit == 1) //we are not interested in p1 and p2, better to cancel them
	{	cancel_p1_p2();}
	if( (ruling_merit == 0) && sampled_merits.empty() && apply_p1_p2() )
	{	//In this special case we should apply the offset immediately
		new_merit = mean_merit(any_std,true);
	}*/
	
	//a good candidate for top5
	//always test the first
	//always test a new tentative winner
	//also test a candidate if we are in peakmode, he is better than the mediam, and at least 90% of the best 
	double top_merit( (ruling_merit==0) ? merit : hexagonal_merit );
	double bottom_merit(  top_merit );
	if(!(sampled_merits.empty()||top5.empty()))
	{
		top_merit = sampled_merits[top5.front()];
		bottom_merit = sampled_merits[top5.back()];
	}
	/*
	if(		sampled_merits.empty() || 
			(new_merit > top_merit) || 
			(	(rank >= peak_rank) && 
				(new_merit >= bottom_merit ) && 
				(new_merit > 0.9*top_merit ) 
			) ) 
	*/
	if(sampled_merits.empty() || (new_merit >= bottom_merit ) || (new_merit > 0.95*top_merit ))
	{
		/*if(reporting && has_console)
		{	
			switch(search_type)
			{
				case 0:
				printf("testing tilt: %lf hel: %lf merit: %lf\n", tilt, hel, new_merit);
				break;
				case 1:
				printf("testing phi: %lf exc: %lf merit: %lf\n", phi, excent, new_merit);
				break;
				default:
					if(reporting && !has_console)
					{	printf("Message ");}
					printf("ERROR unrecognized search type %d in %s line: %d\n",search_type, __FILE__, __LINE__);
					fflush(stdout);
					exit(0);
			}
			fflush(stdout);
		}*/
		//get some statistics for promising points
		
		new_merit = (new_merit + (double)top_points * mean_merit(any_std,false) ) / ( 1.0 + (double)top_points );	
	}
	
	sampled_merits[ind] = new_merit;
	if( (new_merit >= bottom_merit) || (top5.size() < top_points) ) //dont rely on floating point equivalency
	{
		top5.push_back(ind);
		std::sort (top5.begin(),top5.end(), comp_merit_ind );
		while (top5.size() > top_points) {	top5.pop_back();}
		if(top5.front() == ind)//a new best value
		{
			double m_std( any_std );
			double hm_std( any_std );
			const double ruling_contrast( contrast );
			//cancel_p1_p2();
			if(floating_offset)
			{
				if(ruling_merit == 0)
				{	apply_a1_a2();}
				else
				{	apply_p1_p2();}
			}
			//move sampling offset
			offset_XS = offset_X;
			offset_YS = offset_Y;
			floating_offset = false;
			
			if(ruling_merit == 0)
			{
				merit = new_merit;
				ruling_merit = 1;
				hexagonal_merit = mean_merit(hm_std,false);
				ruling_merit = 0;
			}
			else
			{
				hexagonal_merit = new_merit;
				ruling_merit = 0;
				merit = mean_merit(m_std,false);
				ruling_merit = 1;
			}
			contrast = ruling_contrast;
			
			config.save();
			
			if(reporting)
			{	
				switch(search_type)
				{
					case 0:
					printf("new tilt: %lf hel: %lf merit: %lf(%lf), hex_merit: %lf(%lf)\n", tilt, hel, merit, m_std , hexagonal_merit, hm_std);
					break;
					case 1:
					printf("new  phi: %lf exc: %lf merit: %lf(%lf), hex_merit: %lf(%lf)\n", phi, excent, merit, m_std, hexagonal_merit, hm_std);
					break;
					default:
						if(reporting && !has_console)
						{	printf("Message ");} 
						printf("ERROR unrecognized search type %d in %s line: %d\n",search_type, __FILE__, __LINE__);
						fflush(stdout);
						exit(0);
				}
				fflush(stdout);
			}
		}
	}
	cancel_p1_p2();//that should not be needed!?, cheap anyways
	floating_offset = (stability > 0);
	//check_memory();
}

double grid_search(int &rank)
{
	if(rank==max_rank)
	{	
		printf("Message WARNING grid_search is not supposed to be called if rank is already max_rank!\n");
		return -1.0;
	}
	double old_merit( (ruling_merit==0)?merit:hexagonal_merit );
	fprintf(logfile,"grid_search search_type: %d, rank: %d ... ", search_type, rank); fflush(logfile);
	
	
	int dh( (rank > 0) ? (mg_len >> rank) : (mg_len-1) ); 
	int dt( (rank > 0) ? (mg_len >> rank) : (mg_len-1) ); 
	
	int h0( phint_min );
	int t0( ptint_min );
	
	if(rank == 0) //sample the initial 9 points and increase rank 
	{
		//center
		eval_point(h0 + dh/2   , t0 + dt/2 , rank); // initialization and report to master
		//4 direct neighbors
		
		eval_point(h0 + dh    , t0 + dt/2  , rank);
		eval_point(h0         , t0 + dt/2  , rank);
		eval_point(h0 + dh/2  , t0         , rank);
		eval_point(h0 + dh/2  , t0 + dt    , rank);
		//d diagonal neighbors
		eval_point(h0         , t0         , rank); 
		eval_point(h0 + dh    , t0         , rank);
		eval_point(h0         , t0 + dt    , rank);
		eval_point(h0 + dh    , t0 + dt    , rank);
	
	
		++rank;
		dh /= 2;
		dt /= 2;
	}
	
	dh/=2;
	dt/=2;
	
	for(int h( phint_min ) ; busy && (h <= phint_max); h += dh )
	{
		for(int t( ptint_min ) ; busy && (t <= ptint_max); t += dt )
		{	eval_point( h, t , rank);	}	
	}

	++rank;
	
	config.load(); //fetch/restore the best config so far
	if(rank >= peak_rank) //We have at least a 19x19 grid
	{
		
		//tighten peak area to have the five biggest values in the interior
		int nptmax( ptint_min );
		int nptmin( ptint_max );	
		int nphmax( phint_min );
		int nphmin( phint_max );
		int size( top5.size() );
		
		for(int i(0); i < size ; ++i)
		{
			const int ddt ( ( (i==0)&&(stability > 0) ) ? 2*dt : dt ); //just in case highest would lie on edge of peakbox
			const int ddh ( ( (i==0)&&(stability > 0) ) ? 2*dh : dh ); //just in case highest would lie on edge of peakbox
			const int ind( top5[i] );
			const int t( ind % mg_len );
			const int h( ind / mg_len );
			
/*			if(reporting)
			{	
				const double m( sampled_merits[ind] );
				printf("%d\t(%d,%d) with %lf\n", i, t, h, m);
			}*/
			if(nptmax < t + dt ) {nptmax = t + ddt;}
			if(nptmin > t - dt ) {nptmin = t - ddt;}
			if(nphmax < h + dh ) {nphmax = h + ddh;}
			if(nphmin > h - dh ) {nphmin = h - ddh;}		
		}
		
		if( (nptmin != ptint_min) && (nptmin >= 0) )
		{ 
			ptint_min = nptmin;
			peak_tilt_min = box_tilt_min + ptint_min * box_tilt_range / (mg_len-1);
		}
		
		if( (nptmax != ptint_max) && (nptmax < mg_len) )
		{
			ptint_max = nptmax;
			peak_tilt_max = box_tilt_min + ptint_max * box_tilt_range / (mg_len-1);
		}
		peak_tilt_range = peak_tilt_max - peak_tilt_min;
		
		if( (nphmin != phint_min) && (nphmin >= 0) )
		{
			phint_min = nphmin;
			peak_hel_min = box_hel_min + phint_min * box_hel_range / (mg_len-1);
		}
		
		if( (nphmax != phint_max) && (nphmax < mg_len) )
		{
			phint_max = nphmax;
			peak_hel_max = box_hel_min + phint_max * box_hel_range / (mg_len-1);
		}
		peak_hel_range = peak_hel_max - peak_hel_min;
		
		if( rank > peak_rank ) //We dont yet care about close to edge at iniializing the peak box 
		{
			if( (ptint_min == 0) || (phint_min == 0) || (ptint_max == mg_len - 1) || (phint_max == mg_len - 1))
			{
				edge_incidence = true;
				if(reporting)
				{
					switch(search_type)
					{
						case 0:
						printf("%speak box at edge of search box tilt: %lf hel: %lf\n", reporting?"":"Message WARNING ", tilt, hel);
						break;
						case 1:
						printf("%speak box at edge of search box phi: %lf excent: %lf\n", reporting?"":"Message WARNING ", phi, excent);
						break;
						default:
							if(reporting && !has_console)
							{	printf("Message ");} 
							printf("Message ERROR unrecognized search type %d in %s line: %d\n",search_type, __FILE__, __LINE__);
							fflush(stdout);
							exit(0);
					}
					fflush(stdout);
				}
			}
		}
	}
	fprintf(logfile,"done\n"); fflush(logfile);
	return ( ( (ruling_merit==0) ? merit : hexagonal_merit) - old_merit);
}

/*bool peak_stats(int &rank, int &peak_sample, double &var_t, double &var_h)
{
	double sum_h ( 0.0 ),  sum_t ( 0.0 );
	double sum_h2( 0.0 ),  sum_t2( 0.0 );
	double sum_m ( 0.0 );
	peak_sample = 0;
	
	const int size( top5.size() );
	//int best_ind = top5.front(); //dont ever do that without any top5
	for(int i(0); i < size; ++i)
	{
		const int ind(top5[i]);
		if(sampled_merits.count(ind) == 0)
		{	continue; } 
		const int t( ind % mg_len );
		const int h( ind / mg_len );	
		double m( sampled_merits[ind] );
		const double t_val( box_tilt_min + t * box_tilt_range / (mg_len-1) );
		const double h_val( box_hel_min  + h * box_hel_range  / (mg_len-1) );
		
		sum_m  += m;
		sum_t  += m * t_val;
		sum_t2 += m * pow(t_val,2);
		sum_h  += m * h_val;
		sum_h2 += m * pow(h_val,2);
			
		//if(reporting)
		//{
		//	printf("(%d)peak\ttilt: %lf\thel: %lf\tmerit: %lf\n", peak_sample, t_val, h_val, m);
		//}
		++peak_sample;	
	}
	const double mean_t( sum_t/sum_m );
	var_t = fabs( (sum_t2/sum_m - pow(mean_t,2)) / peak_sample );
	const double mean_h( sum_h/sum_m );
	var_h = fabs( (sum_h2/sum_m - pow(mean_h,2)) / peak_sample );
	
	//apply the averaged tilt/phi and hel/ellipse
	
	switch(search_type)
	{
		case 0: //hel and tilt
			set_basis(mean_h, mean_t);
		break;
		case 1: //excent and phi
			set_ellipse(mean_h, mean_t);
		break;
		default:
		    if(reporting && !has_console)
			{	printf("Message ");} 
			printf("ERROR unrecognized search type %d in %s line: %d\n",search_type, __FILE__, __LINE__);
			fflush(stdout);
			exit(0);
	}
	//merit = sampleImg(uc_res);	//twice the effort
	//std::sort (top5.begin(),top5.end(), comp_merit_ind );
	double std( -1.0 );
	merit = mean_merit(std,false);
	
	if( top5.empty() || (merit > (0.5*sampled_merits[top5.front()] + 0.5 * sum_m/peak_sample ) ) )
	{
		if( reporting)
		{
			printf("new tilt: %lf hel: %lf phi: %lf exc: %lf merit: %lf\n", tilt, hel, phi, excent, merit);
			fflush(stdout);
		}
		config.save();
		return true;
	}
	else
	{
		config.load(); //dismiss averaged h,t,merit
		return false;
	}
}*/

void optimize_merit(const bool lasting)
{		
	int iter(0);
	int hex_iter(0);
	int min_iter(MIN_ITER); //at least MIN_ITER of each type final is always search_type 1
	int max_iter(MAX_ITER); //at max MAX_ITER to MAX_ITER_LIMIT of each type final is always search_type 1
	bool another_iter = true;
	double hyper_merit( merit );
	double hyper_hex_merit( hexagonal_merit );
	stability = init_stability;
	int hex_rank(0);
	do //iter loop cyles through search_types before incrementing
	{
		edge_incidence = false;
		if(reporting && (search_type == 0) )
		{
			printf("MERIT optimization %d of (min: %d, max: %d)\n", (iter+1), min_iter, max_iter);
			fflush(stdout);
		}
		sampled_merits.clear();
		top5.clear();
		//cancel_p1_p2();
		
		if( search_type == 0 ) //measure improvement in full cycles 
		{	
			int runs( step_search(hex_rank, true, hex_iter == 0) ); 
			if(runs != -1) //-1 means no good at all
			{
				if (stability >= runs)  
				{	stability -= runs;}
				else
				{	stability = 0;}
				
				
			}
			if(skip_grid) 
			{	return;}
			/*if( (ruling_merit == 1) )
			{
				if( (runs > 0) && (++iter < max_iter) )
				{	continue;}
				else
				{	break;}
			}*/
			hyper_merit = merit;
			hyper_hex_merit = hexagonal_merit;
		}
		
		int rank(0);
		bool got_better( false );
		double old_peak_tilt_range; //( peak_tilt_range )
		double old_peak_hel_range; //( peak_hel_range )
		int t_steps = 4;
		int h_steps = 4;
		bool better_merit( false );
		
		
		set_limits(true);	
		do //going up the ranks and narrowing the search grid
		{	
			report();
			if(reporting)
			{
				printf("%s (%dx%d) lvl: %d  ", (rank >= peak_rank) ? "peak" : "grid", t_steps, h_steps,rank);
				switch(search_type)
				{
					case 0: //hel and tilt
					printf(" %lf < tilt < %lf  %lf < hel < %lf  merit: %lf  hex_merit: %lf\n", 
						peak_tilt_min, peak_tilt_max, peak_hel_min, peak_hel_max, merit, hexagonal_merit);
					break;
					case 1:
					printf(" %lf < phi < %lf  %lf < excent < %lf  merit: %lf hex_merit: %lf\n", 
						peak_tilt_min, peak_tilt_max, peak_hel_min, peak_hel_max, merit, hexagonal_merit);
					break;
					default:
						if(!has_console)
						{	printf("\nMessage ");}
						else
						{	printf("\n");} 
						printf("ERROR unrecognized search type %d in %s line: %d\n",search_type, __FILE__, __LINE__);
						fflush(stdout);
						exit(0);
				}
				fflush(stdout);
			}
			old_peak_tilt_range = peak_tilt_range;
			old_peak_hel_range = peak_hel_range;
			double delta_merit = grid_search(rank);
			if(ruling_merit == 0)
			{
				better_merit = ( (delta_merit/merit > 0.01) || (rank == 0) );
			}
			else
			{
				better_merit = ( (delta_merit/hexagonal_merit > 0.01) || (rank == 0) );	
			}
			
			got_better = better_merit || (old_peak_tilt_range > peak_tilt_range) || (old_peak_hel_range > peak_hel_range);
			if(rank > 0)
			{
				const int res ( mg_len >> rank );
				t_steps = 2 * ( (ptint_max - ptint_min) / res );
				h_steps = 2 * ( (phint_max - phint_min) / res );
				//check if the peak search box is shrinking as the resolution grows 10*10 is performance limit
				got_better &= (t_steps * h_steps <= 10*10);
			}
			else
			{	t_steps = 4;	h_steps = 4;}
		} 
		while ( lasting && busy && (rank < max_rank) && ( got_better || (rank < peak_rank) ) ); //(rank <= peak_rank)
		if( ( (iter+1) == max_iter) && ( max_iter < MAX_ITER_LIMIT ) && edge_incidence && (hexagonal_merit > 0) )
		{
			++max_iter;
			printf("%sraising max_iter for frame: %d to %d \n", reporting ? "" : "Message ", frameNr, max_iter);
			fflush(stdout);
			
		}				
		
		if(search_type == 1) //the last defined search_type
		{	
			++iter;
			edge_incidence = false; //only care about recent cases
			(stability > 0) && (--stability);
		}
		search_type = (search_type + 1) % 2; //2 supported search_types
													  //(!stable_position) || 
		//we are happy as long as either one goes up //MAYBE define hyper_merit as product?
		const bool better_hyper_merit ( ( (merit-hyper_merit)/merit > 0.01) || ( (hexagonal_merit-hyper_hex_merit)/hexagonal_merit > 0.01)  );
		another_iter = ( lasting || (search_type != 0) ) && busy && (iter < max_iter) && (  (search_type == 1) || better_hyper_merit || (iter < min_iter) ); 
		
		if(another_iter) 
		{	report();}
	}
	while( another_iter );
	//Final polishing
	//retune offset as it may be stuck somewhere if the the optimization was cancelled
	{
		config.load(); //just in case there was an SIGINT during optimization
		ruling_merit = 1;
		floating_offset = true;
		optimize_hexagons(0, hex_rank); 
		optimize_hexagons(3, hex_rank); // 0 -> 3, lets see if step search is in overall better than grid_search
		center_offset();
		offset_XS = offset_X;
		offset_YS = offset_Y;
		floating_offset = false; //the offsets ought to be perfect AFTER optimization, but who knows ...
		double m_std( -1.0 );
		double hm_std( -1.0 );
		hexagonal_merit = mean_merit(hm_std,false);
		double const hex_contrast( contrast );
		ruling_merit = 0;
		merit = mean_merit(m_std,false);//get some statistics
		ruling_merit = 1;
		contrast = hex_contrast;
		config.save();
		if(reporting)
		{
			printf("optimized merit: %lf(%lf) hex_merit: %lf(%lf)\n", merit, m_std, hexagonal_merit, hm_std);
			printf("optimized X0,Y0: %lf,%lf tilt,hel: %lf,%lf phi,exc: %lf,%lf\n", offset_X, offset_Y, tilt, hel, phi, excent);
			//offsets will be retuned one last time before cutting the subframes
			fflush(stdout);
		}
	}
}

double mean_merit(double &std, const bool quick)
{
	double mean( 0.0 );
	
	if(ruling_merit == 1)
	{
		mean = get_hex_merit(offset_X,offset_Y,tilt,hel,phi,excent,std,quick);
	}
	else
	{
		int pts(0);
		const int last_ind (top_points);
		double sum(0.0);
		double sum2(0.0);
		do
		{
			++pts;
			const double m( sampleImg(uc_res) ); 
			sum += m;
			sum2 += pow(m,2);
			if (m < 0.9 * merit)
			{	break;}
		}
		while( (!quick) && (pts < last_ind) );
		mean = sum/(double)pts;
		std = sqrt( fabs( sum2/(double)pts - pow(mean,2) ) );
	}
	return mean;
}

int step_search(int &hex_rank, const bool force_first, const bool do_first)
{
	bool start = !config.exists();
	suppress_no_mini_cells = false;
	double m_std( -1.0 );
	double hm_std( -1.0 );
	fprintf(logfile,"step_search start: %s, force_first: %s, do_first: %s\n", start?"true":"false", force_first?"true":"false", do_first?"true":"false"); fflush(logfile);
	double new_merit;//(merit);
	const double old_hex_merit ( hexagonal_merit );
	double new_hex_merit; //( hexagonal_merit ); //scan-build
	floating_offset = true;
	if( start ) //create the starting point
	{
		if(reporting && has_console)
		{
			printf("starting X0,Y0: %lf,%lf tilt,hel: %lf,%lf phi,exc: %lf,%lf\n", offset_X, offset_Y, tilt, hel, phi, excent);
			fflush(stdout);
		}
		ruling_merit = 0;
		new_merit = mean_merit(m_std,true);//quick mode
		merit = new_merit;
		if( (offset_X <= 1.0) || (offset_X >= impWidth-1) || (offset_Y <= 1.0) || (offset_Y >= impHeight-1) )
		{	//quite inaccurate but our best first guess, hex_merit would fail outside the frame 
			apply_a1_a2(); //a1 and a2 are inaccurate but more robust
			center_offset();
		}
		ruling_merit = 1;
		new_hex_merit = mean_merit(hm_std,true); //quick_mode
		hexagonal_merit = new_hex_merit;
		
		apply_p1_p2(); //Also done during optimization
		center_offset();
		offset_XS = offset_X;
		offset_YS = offset_Y;
		
		config.save();
		//dont add the quick point with poor statistics to the top5 list
		report();
		
	}
	if( reporting )
	{	
		printf("%s X0,Y0: %lf,%lf tlt,hel: %lf,%lf phi,exc: %lf,%lf merit: %lf hex_merit: %lf\n",
			start?"initial":"current", offset_X, offset_Y, tilt, hel, phi, excent, merit, hexagonal_merit); 
		fflush(stdout);
	}
	int old_ruling_merit( ruling_merit );
	if(do_first)
	{
		ruling_merit = 1;
		new_hex_merit = optimize_hexagons(0,hex_rank);
		hexagonal_merit = new_hex_merit;
		
		if(force_first) //accept anyways
		{
			offset_XS = offset_X;
			offset_YS = offset_Y;
			ruling_merit = 0;
			new_merit = mean_merit(m_std,false);
			merit = new_merit;
			ruling_merit = 1; //just to be safe
			config.save();
		}
	}
	int runs(0);
	//bool stable_coh_len( false );
	const int hex_rank_reset( hex_rank );
	//cancel_p1_p2();
	if( busy /*&& !stable_coh_len*/) //while
	{
		//the contrast is not even strictly monotone from hollow to ring to bond to atom
		//stable_coh_len = false;
		hex_rank = hex_rank_reset;
		if( (hexagonal_merit < 0.0) || (contrast < min_hex_contrast) )
		{
			//start = false; //scan-build
			if(old_ruling_merit == 1)
			{	
				printf("Message frame: %d bad situation hex_merit: %lf contrast: %lf hex_shape: %lf\n", frameNr, hexagonal_merit, contrast, hex_shape);
				(stability < 2) && (++stability); //dont go to 3 automatically
			}
			ruling_merit = 0; 
			search_pattern = 0;
			return 0;
			///It seems this does not realy help either
			/*
			small_coh_len/=3;
			if(reporting)
			{	printf("smaller coh_len: %d\n", small_coh_len);}
			if(small_coh_len < (int) (12 * (*bondlength) ) )
			{
				small_coh_len = (int) (12 * (*bondlength) );
				stable_coh_len = true;
				if(reporting)
				{
					printf("switching to classic merit minimal small_coh_len: %d\n", small_coh_len);
				}
				ruling_merit = 0; 
				search_pattern = 0;
				return 0; //we did not do any actually useful
			}
			coh_len = small_coh_len 
			if(busy){	optimize_hexagons(0, hex_rank);} //retune offsets with smaller coh_len
			*/ 
		}
		else
		{
			if(old_ruling_merit == 0)
			{	printf("Message frame: %d reactivating hex mode hex_merit: %lf contrast: %lf hex_shape: %lf\n", frameNr, hexagonal_merit, contrast, hex_shape);}
			search_pattern = 1;
			ruling_merit = 1;
			//stable_coh_len = true; //dont make endless loop
			/* 
			small_coh_len = small_coh_len * 3 / 2;
			if(small_coh_len > max_coh_len)
			{	
				small_coh_len = max_coh_len;
				stable_coh_len = true;
				if(reporting)
				{	printf("stable coh_len: %d\n", small_coh_len);}
			}
			else if(reporting)
			{	printf("bigger coh_len: %d\n", small_coh_len);}
			coh_len = small_coh_len;
			*/ 
		}
		++runs;
		if(start)
		{
			if(reporting)
			{	printf("hexagonal single step search rank: %d\n", hex_rank); }
			if(busy){	optimize_hexagons(1, hex_rank);} //1 .. sets new tilt and hel 
			if(busy){	optimize_hexagons(2, hex_rank);} //2 .. sets new phi and excent
			++hex_rank;
			if(busy){	optimize_hexagons(4, hex_rank);} //4 .. sets new tilt and phi
			if(busy){	optimize_hexagons(5, hex_rank);} //5 .. sets new hel and excent
			++hex_rank;
			if(busy){	optimize_hexagons(0, hex_rank);} //0 .. sets new offset_X and offset_Y 
		}
		if(busy && (config.any_changes() || !start) )
		{	
			++hex_rank;
			optimize_hexagons(3, hex_rank);
		} //3 .. all at the same time, tilt, hel,excent and phi are locked as long as offsets are not bounded
		//check for changes in primary data
		if(	config.any_changes() ) 	
		{
			//move sampling offset
			offset_XS = offset_X;
			offset_YS = offset_Y;
			floating_offset = false; //false is more honest
			ruling_merit = 0;
			const double ruling_contrast(contrast);
			merit = mean_merit(m_std,false);//get some statistics
			ruling_merit = 1;
			contrast = ruling_contrast;
			config.save();
			report();
		}
		else //no improvement we are done
		{	
			runs = -1;
			hexagonal_merit = old_hex_merit;
		}
		++hex_rank;
		config.load(); //retrieve the redundant data
	}
	floating_offset = (stability > 0); //spare the effort at fine grid/peak search
	fprintf(logfile,"step_search done\n"); fflush(logfile);
	return (hex_rank>1)?runs:0; //dont decrease stability after first run
}


