#include "subframes.hpp"
#include "globals.hpp"
#include "hexagons.hpp"
#include "ellipse.hpp"
#include "assert.h"
#include "math.h"
//#include <set>
//#include <deque>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>

//private headers
void init_subframes(int most_sf);
//eiter returns position on the frame or -1 ouside
int check_bounds(int dx, int dz, bool apply);
//returns number of complete subframes with offset q0,r0
int mark_subframes(int q0, int r0, bool apply, bool fill_sf);
//checks if the next subframe at ind would be clean and increases next_sf if so
//apply decides if the hexpixels would actually be collected
bool check_and_add_sf(const int ind, int &next_sf, const bool apply, const bool fill_sf);

//void hexagonal_trap(int &lx, int &ly, int &lz, const int N);

//key are positions of possible subframes
//value counts colliding other possible frames
std::map< int, int > sf_collisions;
int monkey_picks(0);
int full_area(0);
//int the_sublattice;
//std::map<int,int> sf_sublattice;
//std::map<int,int> lattice_chart;
int last_x,last_z;
//initiates/refreshes the count of colliding possible subframe locations
void update_sf_collisions(void);
//picks a frame and eliminates all other colliding possible locations
//for_real realizes the location, other wise only discounts neighboring collisions 
void pick_a_sf(int ind, bool for_real);
//picks an least overlapping but not yet totally free location
int most_lonely_sf(void);
//removes possible locations until none of them is any longer colliding
int find_best_tiling(void);

void init_subframes(const int most_sf, bool resample)
{
	const int modelarea( modelsize * modelsize);
	if(resample)
	{
		if(subframes != nullptr)
		{
			for(int sf(0); sf < Num_sf; ++sf)
			{	delete[] subframes[sf];}
			delete[] subframes;
		}
		Num_sf = most_sf;
		subframes = new int*[Num_sf];
		for(int sf(0); sf < Num_sf; ++sf)
		{	
			subframes[sf] = new int[modelarea];
			memset(subframes[sf], 0, modelarea*sizeof(int));
		}
	}
	memcpy(hexarea, hexagons, hexframelen*sizeof(int));
	delete[] sf_sum;
	sf_sum = new int[modelarea];
	//assert(sf_sum != nullptr);
	memset(sf_sum, 0, modelarea*sizeof(int));	
}

int check_bounds(const int dx, const int dz)		
{		
	const int dq( dx+(dz-(dz&1))/2 );
	const int dr( dz );
	const int q( hexImpWidth/2 + dq );
	const int r( hexImpHeight/2 + dr );
	if( (q < 0) || (q >= hexImpWidth) ||
	(r < 0) || (r >= hexImpHeight) )
	{	return -1;}
	
	const int pos( q + r * hexImpWidth );
	if(hexagons[pos]==-1)
	{	return -1;}
	return pos;
}		

bool check_and_add_sf(const int ind, int &next_sf, const bool apply, const bool fill_sf)
{	
	const int q( ind % hexImpWidth - hexImpWidth/2);
	const int r( ind / hexImpWidth - hexImpHeight/2);
	const int cx( q - (r - (r&1) )/2);
	const int cz( r );
	
	bool accepted( true );
	const bool ald( allow_light_dirt);
	
	const int rHH( modelsize/2 );
	bool trivial( *bondlength == rHH );
	const int rH(rHH+2); //FIXME why do we need that safety border at all?
						//I guess it is because sampling is asymmetric around a pixel
						//hence the top and left edge have less intensity	
	if( !apply || sf_collisions.empty() ) // || true //for double checking
	{
		int ldc(0); //light_dirt_counter
		
		//long sum( 0 );
		for(int dx(-rH); dx < rH ; ++dx)
		{
			int dz1( -rH );
			int dz2( rH - dx );
			if (dx < 0)
			{
				dz1 = -rH - dx;
				dz2 = rH;
			}
			for(int dz( dz1 ); dz < dz2; ++dz)
			{
				//const int dy( -dx - dz);
				const int px( cx + dx );
				const int pz( cz + dz );
				const int qh(px + (pz - (pz&1) )/2 + hexImpWidth/2);
				const int rh(pz + hexImpHeight/2);
				if( (qh < 0) || (qh >= hexImpWidth) ||
					(rh < 0) || (rh >= hexImpHeight) )
				{	//might occur if the hexedImg actually cuts the frame
					accepted = false;
					//dz = dz2; //quit for dz //scan-build
					dx = rH; //quit for dx
					//printf("o");fflush(stdout);
					break;
				}	
				const int pos(qh + rh * hexImpWidth);
				/*
				if(pos < 0 || pos >= hexframelen)
				{
					printf("check_and_add_sf\tpos:%d qh:%d rh:%d hiW:%d hiH:%d hfl:%d\n",pos,qh,rh,hexImpWidth,hexImpHeight,hexframelen);
					fflush(stdout);
					exit(0);
				}
				*/
				//const int hval( hexagons[ pos ]);
				//sum += (long)hval;
				
				const int hmsk( hex_mask[pos] );
				
				if( (ald && (hmsk != 1) && (hmsk != light_dirt_val) ) ||
					( (!ald) && (hmsk != 1) ) )
				{	
					accepted = false;
					//dz = dz2; //quit for dz//scan-build
					dx = rH; //quit for dx
					//printf("d");fflush(stdout);
					break;
				}
				
				
				if(ald && ( hmsk == light_dirt_val ) )
				{
					++ldc;
				}
				
			}
		}
		if(ald && accepted)
		{
			const double light_dirt_coverage( ((double)ldc)/(3.0*rH*rH) );
			accepted = ( (light_dirt_coverage >= min_light_dirt) && 
					   (light_dirt_coverage <= max_light_dirt) );
		}
		/*
		if(accepted)
		{	
			printf("acceptable pos: %d acceppted: %s\n",ind,accepted?"true":"false");
			fflush(stdout);
			exit(0);	
		}
		*/ 
		/*
		if(accepted && (sum == 0))
		{
			if(!reporting)
			{	printf("Message ");}
			printf("WARNING skipped empty subframe centered at %d,%d\n", ind % hexImpWidth, ind / hexImpWidth);
			fflush(stdout);
		}
		accepted &= (sum > 0);
		*/
		if( (!trivial) && (!apply) && accepted) //initial scan
		{
			sf_collisions[ind] = 0;
			//sf_sublattice[ind] = the_sublattice;
		}

	}
	
	if(apply && accepted ) //accepted may or may not have been rechecked
	{	
		const int rH( modelsize/2 );
		const int bl( *bondlength ); 
		const int mcx( rH - (rH - (rH&1) )/2 );
		const int mcz( rH);
		int *const subf( fill_sf ? subframes[next_sf] : nullptr );
		
		for(int dx( -rH); dx < rH ; ++dx)
		{
			int dz1( -rH );
			int dz2( rH - dx );
			if (dx < 0)
			{
				dz1 = -rH - dx;
				dz2 = rH;
			}
			for(int dz( dz1); dz < dz2; ++dz)
			{
				const int px( cx + dx );
				const int pz( cz + dz );
				const int qh(px + (pz - (pz&1) )/2 + hexImpWidth/2);
				const int rh(pz + hexImpHeight/2);
				const int pos( qh + rh * hexImpWidth );
				const int mdx( mcx + dx);
				const int mdz( mcz + dz);
				const int mq( mdx + (mdz - (mdz&1) )/2 );
				const int mr( mdz );
				const int mpos( mq + mr * modelsize);
				const int hex_pix( hexagons[ pos ] );
				sf_sum[mpos] += hex_pix;
				if(fill_sf)
				{
					if(highlight_sf)
					{
						hexarea[ pos ] += hex_pix;
						hexarea[ mq + mr * hexImpWidth] += (2*hex_pix+Num_sf/2)/Num_sf; //DEBUG display sum in top left corner
					}
					subf[mpos] = hex_pix;
					hex_avg += (double)hex_pix;
					int xt(dx), yt(-dx-dz), zt(dz);
					hexagonal_trap(xt,yt,zt,bl);
					
					int mm = abs(xt);
					int nn = abs(xt);
					if(abs(yt) > mm) {	mm = abs(yt);}
					if(abs(zt) > mm) {	mm = abs(zt);} 
					if(abs(yt) < nn) {	nn = abs(yt);}
					if(abs(zt) < nn) {	nn = abs(zt);}
					if( mm == 0 )
					{	hex_low += (double)hex_pix;}
					else if( (mm == bl) && (nn == 0) )
					{	hex_high += (double)hex_pix;}
				}
			}
		}			
	}
	if(accepted)
	{	++next_sf;}
	return accepted;
}

int mark_subframes(const int ddx, const int ddz, const bool apply, const bool fill_sf)
{
	//check_memory();
	const int dx( (*bondlength) * (2*ddx - ddz));
	const int dz( (*bondlength) * (2*ddz - ddx));
	const int rH( (soft_sf?solidsize:modelsize)/2 );
	
	int new_sf( 0 );
	int good_sf( 0 );
	int radius( 0 );
	if(!apply || sf_collisions.empty() ) //classic search pattern
	{
		do
		{
			new_sf = 0;
			for(int x(-radius); x<=radius ;++x)
			{
				for(int z(-radius); z<= radius;
						z += ( (x==-radius)||(x==+radius) ) ? 1 : (2*radius) )
				{	
					const int ddx = rH*(2*x-z);
					const int ddz = rH*(2*z-x);
					const int pos( check_bounds(dx + ddx, dz + ddz) );
					if(pos != -1)
					{	
						check_and_add_sf(pos, good_sf, apply, fill_sf); //false == apply
						++new_sf;
						//printf("pos: %d good_sf: %d\n",pos,good_sf);
						//fflush(stdout);
					}
				}	
			}
			++radius;	
		}
		while(new_sf > 0);
		return good_sf;
	}
	else // apply
	{
		const std::map< int,int >::iterator itb = sf_collisions.begin();
		const std::map< int,int >::iterator ite = sf_collisions.end();
		std::map< int, int >::iterator it;
		for(it = itb; it != ite; ++it)
		{
			//assert(it->second == 0); // assert we dont apply any colliding subframes
			//const int pos(it->first);
			check_and_add_sf(it->first, good_sf, apply, fill_sf); // true == apply 
		}
		return good_sf;
	}
}

int optimize_subframes(bool refresh, bool resample)
{
	//check_memory();
	double hsp(0.0);
	const int bl( modelsize/(2* (*bondlength) ) );
	if(resample)
	{	//skip calculating uc_mini and hex_merit
		set_ellipse(excent,phi);
		allow_light_dirt = (max_light_dirt > 0.0);
		sample_hexagons(hsp,true); //true = skip_uc_mini and hence dont call back optimize_subframes
		soft_sf = true;
	}
	/*else
	{
		assert(!allow_light_dirt);
		asser(!soft_sf);
	}*/ 
	
	int most_sf( Num_sf );
	if( refresh )
	{
		sf_collisions.clear();
		bool first( true );
		//the_sublattice = 0;
		for(int dx(0); dx < bl; ++dx)
		{
			for(int dz(0); dz < bl; ++dz)
			{
				
				const int sf( mark_subframes( dx, dz, false, false) ); //false just count dont apply
				//++the_sublattice;
				//printf("bl: %d offset_Q: %d  offset_R: %d  sf: %d  most_sf: %d\n",bl,dx,dz,sf,most_sf);
				//fflush(stdout);
				if(first || (sf > most_sf) )
				{
					first = false;
					most_sf = sf;
					offset_Q = dx;
					offset_R = dz;
				}
			}
		}
	}
	if(resample)
	{
		hex_avg = 0.0;
		hex_low = 0.0;
		hex_high = 0.0;
	}
	//assert(most_sf > 0);
	if(most_sf > 0) //otherwise spare the effort
	{	
		int simple_counting = most_sf;
		const int real_size( soft_sf?solidsize:modelsize );
		if(2*(*bondlength) != real_size)
		{	most_sf = find_best_tiling();}
		if( most_sf < simple_counting )
		{
			most_sf = simple_counting;
			sf_collisions.clear();
		}
		init_subframes(most_sf,resample);
		
		mark_subframes(offset_Q,offset_R, true, resample);
		if(resample)
		{
			hex_low  /= ( (double)(most_sf * bl * bl)); 
			hex_high /= ( (double)(2 * most_sf * bl * bl));
			hex_avg  /= ( (double)(most_sf * modelsize * modelsize * 3 / 4) );
		}
	}
	else
	{	Num_sf = 0;} //be sure we dont ever try to read the subframes
	if(resample) 
	{	config.save();}
	//assert( sf_sum != nullptr );
	allow_light_dirt = false;
	soft_sf = false;
	return most_sf;
}

//update the count of colliding subframes

void update_sf_collissions(void)
{
	//const int bl (*bondlength);
	//const int rH( 2*(modelsize/(2*bl)-1) ); //need to search within double radius but not reach next touching position
	//assert(rH > 0); //not supposed to be called at all in trivial case
	
	const std::map< int,int >::iterator itb = sf_collisions.begin();
	const std::map< int,int >::iterator ite = sf_collisions.end();
	std::map< int,int >::iterator it,it2;
	const int real_size (soft_sf?solidsize:modelsize);
	for(it = itb; it != ite; ++it)
	{
		//wipe all earlier collision counts
		it->second = 0; 
	}
	//stupid but working
	for(it = itb; it != ite; ++it)
	{
		const int ind = it->first;
		const int q( ind % hexImpWidth - hexImpWidth/2);
		const int r( ind / hexImpWidth - hexImpHeight/2);
		const int cx( q - (r - (r&1) )/2);
		const int cz( r );
		
		for(it2 = std::next(it,1); it2 != ite; ++it2)
		{
			const int ind2 = it2->first;
			const int q2( ind2 % hexImpWidth - hexImpWidth/2);
			const int r2( ind2 / hexImpWidth - hexImpHeight/2);
			const int cx2( q2 - (r2 - (r2&1) )/2);
			const int cz2( r2 );
			
			const int dx( cx - cx2 );
			const int dz( cz - cz2 );
			const int dy( -(dx+dz) );
			
			if( (abs(dx)+abs(dy)+abs(dz)) < 2*real_size )
			{
				it->second++;
				it2->second++;
			}
			
		}
	} 	
}

void pick_a_sf(const int ind, const bool for_real)
{
	//assert(sf_collisions.count(ind)); //check for valid location 
	int collisions( sf_collisions[ind] );
	//assert(collisions>0);       //dont use this function for trivial cases
	//const int bl (*bondlength); 
	const int real_size (soft_sf?solidsize:modelsize);
	//const int rH( 2*(real_size/(2*bl)-1) ); //need to search within double radius
	
	//assert(rH > 0); //not supposed to be called at all in trivial case
	
	const int q( ind % hexImpWidth - hexImpWidth/2);
	const int r( ind / hexImpWidth - hexImpHeight/2);
	const int cx( q - (r - (r&1) )/2);
	const int cz( r );
	last_x = cx;
	last_z = cz;
	std::map<int,int>::iterator itb = sf_collisions.begin();
	std::map<int,int>::iterator ite = sf_collisions.end();
	std::map<int,int>::iterator it2 = itb;
	
	while(it2 != ite)
	{
		
		const int ind2 = it2->first; 
		
		if(ind2 == ind)//skip home
		{	
			++it2;
			continue;
		} 
		
		const int q2( ind2 % hexImpWidth - hexImpWidth/2);
		const int r2( ind2 / hexImpWidth - hexImpHeight/2);
		const int cx2( q2 - (r2 - (r2&1) )/2);
		const int cz2( r2 );
		
		const int dx( cx - cx2 );
		const int dz( cz - cz2 );
		const int dy( -(dx+dz) );
		
		
		if(for_real)
		{
			if( (abs(dx)+abs(dy)+abs(dz)) < 2*real_size )
			{
				pick_a_sf(ind2, false);
				it2 = sf_collisions.erase(it2);
				//assert(collisions-1 == sf_collisions[ind]);
				if( --collisions == 0 )
				{	break;}
			}
			else
			{
				++it2;
			}
		}
		else
		{
			if( (abs(dx)+abs(dy)+abs(dz)) < 2*real_size )
			{
				sf_collisions[ind2]--;
				sf_collisions[ind]--;
				//assert(collisions-1 == sf_collisions[ind]);
				if( --collisions == 0 )
				{	break;}
			}
			++it2;
			
		}
		
	}
	//assert(sf_collisions[ind].size()==0);
	sf_collisions[ind]=0;
	
	
	
	
	
	//update_sf_collissions(); //bulk update		
	//printf("remaining neigbohors: %d\n", sf_collisions[ind]);
	//fflush(stdout);
}

int most_lonely_sf(void)
{
	const std::map< int, int >::iterator itb = sf_collisions.begin();
	const std::map< int, int >::iterator ite = sf_collisions.end();
	std::map< int, int >::iterator it;
	int smallest = 0; // we are looking for the smallest value that is not 0
	//int sublattice_count = -1;
	int ind = -1; // key belonging to the smallest value
	bool first( true );
	bool entire_frame(false);
	//const bool get_corner = (last_x == -1) && (last_y == -1);
	std::vector<int> candidates;
	candidates.reserve(128); //should totally suffice
	std::vector<int> next_candidates;
	next_candidates.reserve(128); 
	const int real_size (soft_sf?solidsize:modelsize);
	for(it = itb; it != ite; ++it)
	{
		const int collisions (it->second);
		if(collisions < 0)
		{
			printf("neg. collisions at ind: %d\tcollisions: %d\n", it->first, collisions);
			fflush(stdout);
			exit(0);
		}
		if(collisions > 0) // >0
		{
			const int pos( it->first );
			/*
			const int q( pos % hexImpWidth - hexImpWidth/2);
			const int r( pos / hexImpWidth - hexImpHeight/2);
			const int cx( q - (r - (r&1) )/2);
			const int cz( r );
			const int cy( -cx - cz);
			const int last_y(-last_x -last_z);
				
				
			const int dist( (abs(cx-last_x)+abs(cz-last_z)+abs(cy-last_y))/2 );
			*/
			if(first || (collisions <= smallest) )
			{
				if(collisions < smallest)
				{
					candidates.clear();
				}
				first = false;
				smallest = collisions;
				ind = pos;
				candidates.push_back(ind);	
			}
		}
	}
	// no upper limit, but a break as soon as the entrie_frame is connected
	for(int range(1); (candidates.size() > 1); ++range ) 
	{
		bool first(true);
		unsigned int lowest(0);
		unsigned int max_radius(0);
		for(unsigned int i(0); i < candidates.size(); ++i)
		{
			const int loc = candidates[i];
			const int q( loc % hexImpWidth - hexImpWidth/2);
			const int r( loc / hexImpWidth - hexImpHeight/2);
			const int cx( q - (r - (r&1) )/2);
			const int cz( r );
			const unsigned int radius( (abs(cx) + abs(cz) + abs(cx+cz))/2 );
			unsigned int close_by(0);
			
			for(it = itb; it != ite; ++it) //sf_collisions
			{
				const int loc2 = it->first;
				const int q2( loc2 % hexImpWidth - hexImpWidth/2);
				const int r2( loc2 / hexImpWidth - hexImpHeight/2);
				const int cx2( q2 - (r2 - (r2&1) )/2);
				const int cz2( r2 );
			
				const int dx( cx - cx2 );
				const int dz( cz - cz2 );
				const int dy( -(dx+dz) );
				const int dist = (abs(dx)+abs(dy)+abs(dx))/2; 
				if(	(dist <= (real_size + range * real_size)) )
				{
					++close_by;
				}	
			} //sf_collisions	
			if( first || (close_by <= lowest) )
			{
				if(first || (close_by < lowest) )
				{
					next_candidates.clear();
				}
				else if (radius > max_radius) //close_by == lowest, but this candidate is further from center
				{
					next_candidates.clear();
				}
				
				first = false;
				lowest = close_by;
				max_radius = radius;
				next_candidates.push_back(loc);
				//assert(sf_collisions.count(loc));
			}
		
		}// end candidates	
		if(!next_candidates.empty())
		{	
			candidates = next_candidates;
			next_candidates.clear();
			//assert(sf_collisions.count(candidates.front()));	
		}
		
		if( (lowest == sf_collisions.size() ) && (candidates.size() > 1) )
		{
			entire_frame = true;
			++full_area;
			break;
		}
	}//end range
	int any = -1;//-1 no further colliding subframes
	if(!candidates.empty())
	{	
		any = candidates.back(); 
		//assert(sf_collisions.count(any));
		//What to do in that case? just stay go for the monky picks?
		if(entire_frame)
		{
			any = *(std::max_element(candidates.begin(),candidates.end()));
			//case is handled
			full_area = 0;
			//entire_frame = false; //scan-build
		}
		else if(candidates.size() > 1) //may be as good as a monkey
		{
			++monkey_picks;
			any = candidates[rand()%candidates.size()];	
			//assert(sf_collisions.count(any));
		}
	}
	
	return any;	
}

int find_best_tiling(void)
{
	monkey_picks = 0;
	full_area = 0;
	/*
	const int bl (*bondlength);
	const int real_size( soft_sf?solidsize:modelsize );
	if( real_size == (2*bl) )
	{
		//confirmed with solidsize and modelsize 29.02.2016
		assert(false); //not supposed to be called at all in trivial case
		return sf_collisions.size();
	}
	*/ 
	last_x = -1;
	last_z = -1;
	//lattice_chart.clear();
	update_sf_collissions();
	int ind( most_lonely_sf() ); //probably in some corner
	while( ind != -1   ) //-1 no more collissions
	{
		//printf("picking sf at: %d\t overlapp: %d\t total candidates: %lu\n", ind, sf_collisions[ind], sf_collisions.size());
		//fflush(stdout);
		pick_a_sf(ind,true); //solidfy this one and dismiss all other colliding candidates
		ind = most_lonely_sf();
	}
	 
	if(full_area > 0)
	{
		if(!reporting)
		{	printf("Message WARNING ");}
		printf("frame: %d indecisive tiling heuristics %d\n", frameNr, full_area);
		fflush(stdout);
		full_area = 0;	
	}
	
	
	if(monkey_picks > 0)
	{
		if(!reporting)
		{	printf("Message WARNING ");}
		printf("frame: %d needed %d monkey picks for tiling\n", frameNr, monkey_picks);
		fflush(stdout);
		monkey_picks = 0;	
	}
	
	return sf_collisions.size();
}
