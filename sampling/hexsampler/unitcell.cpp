#include "globals.hpp"
#include "basis.hpp"
#include "ellipse.hpp"
#include "unitcell.hpp"
#include "offset.hpp"
#include "xorshift1024star.hpp"
#include <assert.h>
// private headers
inline __attribute__((always_inline))
int ind00_a1_a2(double a1, double a2); //checks one corner
inline __attribute__((always_inline))
int ind0_a1_a2(double a1, double a2); //checks all four corners
inline __attribute__((always_inline))
int rind00_a1_a2(double a1, double a2); //checks one corner
inline __attribute__((always_inline))
int rind0_a1_a2(double a1, double a2); //checks all four corners
inline __attribute__((always_inline))
int ind1_a1_a2(double a1, double a2);
inline __attribute__((always_inline))
void init_sampling(void);
inline __attribute__((always_inline))
int periodic_ind(int x, int y);
inline __attribute__((always_inline))
int mirrored_ind(int ind);
inline __attribute__((always_inline))
void smooth_uc(int* uc, int* sm);

inline __attribute__((always_inline))
int periodic_ind(int x, int y)
{
	x = (x+uc_res) % uc_res;
	y = (y+uc_res) % uc_res;
	return (x + y*uc_res);
}

inline __attribute__((always_inline))
void smooth_uc(int* uc, int* sm)
{
	for(int y(0); y < uc_res; ++y)
	{
		for(int x(0); x < uc_res; ++x)
		{
			int sum( 3+1*uc[x+y*uc_res] ); //add 3 to center for correct rounding
			int ind( periodic_ind(x+1,y));
			sum += uc[ind];
			ind = periodic_ind(x-1,y);
			sum += uc[ind];
			ind = periodic_ind(x,y+1);
			sum += uc[ind];
			ind = periodic_ind(x,y-1);
			sum += uc[ind];
			ind = periodic_ind(x+1,y-1);
			sum += uc[ind];
			ind = periodic_ind(x-1,y+1);
			sum += uc[ind];
			sm[x+y*uc_res] = sum/7;
		}		
	}
}


void init_sampling(void)
{
	uc_stepX = -1;
	uc_stepY = -1;
	ruc_stepX = -1;
	ruc_stepY = -1;
	
	distort_basis(true);//ensure elliptic distortions x,y are actuall pixels
	//pending_offsets = false;
	if( unitcells != nullptr )
	{
		int **uc(unitcells);
		for(int raw_ind = 0; raw_ind < uc_len ; ++raw_ind)
		{	delete[] *(uc++);}
		delete[] unitcells;
		unitcells = nullptr;
		delete[] dirty;
		dirty = nullptr;	
	}
	
	if( runitcells != nullptr )
	{
		int **ruc(runitcells);
		for(int raw_ind = 0; raw_ind < ruc_len ; ++raw_ind)
		{	delete[] *(ruc++);}
		delete[] runitcells;
		runitcells = nullptr;
		delete[] rdirty;
		rdirty = nullptr;	
	}
	set_uc_spacings();
	
	
	dirty = new bool[uc_len];
	unitcells = new int*[uc_len];
	//is there a way to use memset for these?
	for(int ind( 0 ) ; ind < uc_len; ++ind )
	{
		dirty[ind] = false;
		unitcells[ind] = nullptr;
	}
	rdirty = new bool[ruc_len];
	runitcells = new int*[ruc_len];
	for(int rind( 0 ) ; rind < ruc_len; ++rind )
	{
		rdirty[rind] = false;
		runitcells[rind] = nullptr;
	}
	return;
}

inline __attribute__((always_inline))
int ind00_a1_a2(double a1, double a2)
{
	const int x( (int) ( floor(xd_a(a1, a2) + offset_XS) ) ); //FIX now cast to int AFTER addition
	const int y( (int) ( floor(yd_a(a1, a2) + offset_YS) ) ); //FIX now cast to int AFTER addition
	//if( (x<0) || (y<0) || (x>=impWidth) || (y>=impHeight) )
	//{	return -1;}
	return ( ( (x<0) || (y<0) || (x>=impWidth) || (y>=impHeight) ) ? -1 : ( (x/uc_gridX) + (y/uc_gridY) * uc_rangeX ) );
}

inline __attribute__((always_inline))
int rind00_a1_a2(double a1, double a2)
{
	const int x( (int) ( floor(xd_ra(a1, a2) + offset_XS ) ) ); //FIX now cast to int AFTER addition
	const int y( (int) ( floor(yd_ra(a1, a2) + offset_YS ) ) ); //FIX now cast to int AFTER addition
	//if( (x<0) || (y<0) || (x>=impWidth) || (y>=impHeight) )
	//{	return -1;}
	return ( ( (x<0) || (y<0) || (x>=impWidth) || (y>=impHeight) ) ? -1 : ( (x/ruc_gridX) + (y/ruc_gridY) * ruc_rangeX ) );
}

inline __attribute__((always_inline))
int ind0_a1_a2(double a1, double a2)
{
	//define the four corners
	const double a1l( floor(a1) );
	const double a1h( a1l + 1.0 );//ceil(a1);
	const double a2l( floor(a2) );
	const double a2h( a2l + 1.0 ); //ceil(a2)
	//check if any of them lies outside the frame
	return ( ( ind00_a1_a2(a1l, a2h) < 0 ) ||
			 ( ind00_a1_a2(a1h, a2l) < 0 ) ||
			 ( ind00_a1_a2(a1h, a2h) < 0 ) ) 
			? -1 : ind00_a1_a2(a1l, a2l);
}

inline __attribute__((always_inline))
int rind0_a1_a2(double a1, double a2)
{
	//define the four corners
	const double a1l( floor(a1) );
	const double a1h( a1l + 1.0 ); //ceil(a1);
	const double a2l( floor(a2) );
	const double a2h( a2l + 1.0 ); //ceil(a2);
	//check if any of them lies outside the frame
	return ( ( rind00_a1_a2(a1l, a2h) < 0 ) ||
			 ( rind00_a1_a2(a1h, a2l) < 0 ) ||
			 ( rind00_a1_a2(a1h, a2h) < 0 ) ) 
			? -1 : rind00_a1_a2(a1l, a2l);
}

inline __attribute__((always_inline))
int ind1_a1_a2(double a1len, double a2len)
{
	const double fa1( a1len - floor(a1len) );
	const double fa2( a2len - floor(a2len) );
	//strict periodic condition 0<= ia1,ia2 < uc_res
	const int ia1( ( (int)floor(fa1 * uc_res) + uc_res ) % uc_res );
	const int ia2( ( (int)floor(fa2 * uc_res) + uc_res ) % uc_res );
	//ia1 = (ia1+uc_res) % uc_res; //now inside initialization
	//ia2 = (ia2+uc_res) % uc_res; //now inside initialization
	//if(ia1 == uc_res) {ia1 = 0;}
	//if(ia2 == uc_res) {ia2 = 0;}
	return (ia1 + uc_res * ia2);
}

double sampleImg(int resolution)
{
	//assert(resolution > 0);
	uc_res = resolution;
	uc_area = resolution*resolution;
	init_sampling();
	const int *fr(frame);
	
	const int mskWidth(impWidth/mask_scaling);
	//int cell_check(0);
	//int cell_rcheck(0);
	
	for(int raw_ind(0); raw_ind < framelen ; ++raw_ind)
	{
		const int val( *(fr++) );
		
		const int dv0( 1 + (val-sample_rate/2)/sample_rate);
		
		const int x( raw_ind % impWidth);
		const int y( raw_ind / impWidth);
		const int flg( mask[ x/mask_scaling + y/mask_scaling * mskWidth ] );
		const double cx( (double)x - offset_XS );
		const double cy( (double)y - offset_YS );
		//if(true || (cx*cx+cy*cy < coh_radius2)) //only sample the central circle
		//{}
			for(int v(0); v < val; v+=dv0) //very clean sampling where a pixel may even be shared among adjacent unitcells
			{
				int dv( ( (v+dv0) > val ) ? (val-v) : dv0  );
				if(dv==0) {	break;}
				const double subx( rand_d() );
				const double suby( rand_d() ); 
				
				{
					const double  a1len(  a1d(cx+subx,cy+suby) );
					const double  a2len(  a2d(cx+subx,cy+suby) );
					const int uc_key(  ind0_a1_a2( a1len,  a2len) );
					
					if( (uc_key != -1) && (!dirty[uc_key]) ) 
					//the corrsponding unitcell is enitrely inside the frame and not yet known to be dirty
					{
						//assert(uc_key < uc_len);
						if( (unitcells[uc_key] == nullptr) && /*(val < dirt_val)*/ (flg&1) ) //on graphene
						{
							unitcells[uc_key] = new int[uc_area];
							memset(unitcells[uc_key], 0, uc_area*sizeof(int));
							//++cell_check;
						}
						 
						if( !(flg==1) ) //( (val >= dirt_val) ) //anything else than graphene
						{
							dirty[uc_key] = true;
							if(unitcells[uc_key] != nullptr) //could spare this check
							{
								delete[] unitcells[uc_key];
								unitcells[uc_key] = nullptr;
								//--cell_check;
							}
							 
						} 
						
						if( !dirty[uc_key] )
						{	
							const int loc_ind( ind1_a1_a2( a1len, a2len) );
							unitcells[uc_key][loc_ind] += dv;
						}
						//else assert( unitcells[uc_key] == nullptr );
					}
				}
				//TODO make that optional
				{
					const double ra1len( ra1d( cx+subx, cy+suby) ); //
					const double ra2len( ra2d( cx+subx, cy+suby) );
					const int ruc_key( rind0_a1_a2(ra1len, ra2len) );
					
					if( (ruc_key != -1) && (!rdirty[ruc_key]) ) //the corrsponding unitcell is enitrely inside the frame
					{
						//assert(ruc_key < ruc_len);
						if(  (runitcells[ruc_key] == nullptr) && (flg&1)/*(val < dirt_val)*/ ) //on graphene
						{
							runitcells[ruc_key] = new int[uc_area];
							memset(runitcells[ruc_key], 0, uc_area*sizeof(int));
							//++cell_rcheck;
						}
						 
						if( !(flg==1) ) //( (val >= dirt_val) ) //anything else than graphene
						{
							rdirty[ruc_key] = true; 
							if(runitcells[ruc_key] != nullptr) //could spare that check
							{
								delete[] runitcells[ruc_key];
								runitcells[ruc_key] = nullptr;
								//--cell_rcheck;
							}
							
						} 
						
						if( !rdirty[ruc_key] )
						{	
							const int loc_ind( ind1_a1_a2( ra1len, ra2len) );
							runitcells[ruc_key][loc_ind] += dv;
						}
						//else assert( runitcells[ruc_key] == nullptr );
					}
				}
				if(!(flg==1))//(val >= dirt_val) //anything else than graphene
				{	break;}
			}
		
	}
	
	delete[] uc_sum; //save to call with nullptr
	uc_sum = new int[uc_area];
	delete[] ruc_sum; //save to call with nullptr
	ruc_sum = new int[uc_area];
	delete[] uc_rough; //save to call with nullptr
	uc_rough = new int[uc_area];
	//delete[] uc_stat; //save to call with nullptr
	//uc_stat = new double[uc_area];
	delete[] uc_pos; //save to call with nullptr
	uc_pos = new int[uc_area];
	delete[] uc_gauss; //save to call with nullptr
	uc_gauss = new int[uc_area];
	delete[] ruc_gauss; //save to call with nullptr
	ruc_gauss = new int[uc_area];
		
	memset(uc_pos,0,uc_area*sizeof(int));
	//memset(uc_rough,0,uc_area*sizeof(int)); //should not be needed
	memset(uc_sum,0,uc_area*sizeof(int));
	memset(ruc_sum,0,uc_area*sizeof(int));
	//memset(uc_stat,0,uc_area*sizeof(double));
	
	int **uc( unitcells );
	int good_cells( 0 );
	
	for(int ind(0); ind < uc_len ; ++ind)
	{
		if( *uc != nullptr ) 
		{
			++good_cells;
			int *s( uc_sum );
			//double *stat( uc_stat );
			int *a( *uc );
			
			for(int ind0(0); ind0 < uc_area; ++ind0)
			{	
				const int val( *(a++) );
				
				*(s++) += val;
				//*(stat++) += (double)(val*val);
			}	
		}
		++uc;
	}
	//assert(cell_check == good_cells);
	
	int **ruc( runitcells );
	int good_rcells( 0 );
	for(int ind(0); ind < ruc_len ; ++ind)
	{
		if( *ruc != nullptr ) 
		{
			++good_rcells;
			int *s( ruc_sum );
			int *a( *ruc );
			for(int ind0(0); ind0 < uc_area; ++ind0)
			{	
				const int val( *(a++) );
				*(s++) += val;
			}	
		}
		++ruc;
	}
	//assert(cell_rcheck == good_rcells);
	
	/*if(reporting)
	{	printf("good uc: %d\n", good_cells);fflush(stdout);}*/
	//if(smooth_passes==0)
	//{	memcpy(uc_rough,ruc_sum,uc_area*sizeof(int));}
	for(int s(0); s < smooth_passes; ++s)
	{
		memcpy(uc_rough,ruc_sum,uc_area*sizeof(int));
		smooth_uc(uc_rough,ruc_sum); //ruc_sum becomes smoothed version of uc_rough
	}
	
	
	//if(smooth_passes==0)
	//{	memcpy(uc_rough,uc_sum,uc_area*sizeof(int));}
	for(int s(0); s < smooth_passes; ++s)
	{
		memcpy(uc_rough,uc_sum,uc_area*sizeof(int));
		smooth_uc(uc_rough,uc_sum); //uc_sum becomes smoothed version of uc_rough
	}
	 data_sum = 0;
	rdata_sum = 0;
	//double variance( 0.0 );
	//double stddev( 0.0 ); 
	double sum(0.0);
	double sum2(0.0);
	int *s(uc_sum);
	int *rs(ruc_sum);
	//double *stat( uc_stat);
	bool first(true);
	const int total_cells( good_cells + good_rcells); 
	for(int ind0(0); ind0 < uc_area; ++ind0)
	{	
		const int val( *(s++) );
		
		data_sum += val;
		const int rval( *(rs++) );
		rdata_sum += rval;
		const double mean( ((double)val)/(total_cells) );
		//const int sval(val+rval);
		sum  += (double)(mean);
		sum2 += (double)(mean*mean);
		
		//*stat =  fabs(*stat/(total_cells) - mean*mean) ;//total variance
		//variance += *stat;
		//stddev += sqrt(*stat); //total stddev
		//++stat;
		if( first || (rval > rdata_max) )
		{	rdata_max = (double)rval;}
		if( first || (rval < rdata_min) )
		{	rdata_min = (double)rval;}
		if( first || (val > data_max) )
		{	data_max = (double)val;}
		if( first || (val < data_min) )
		{
			data_min = (double)val;
			first = false;
		}
	} 
	/*
	printf("good_cells: %d  good_rcells: %d\n", good_cells, good_rcells);
	fflush(stdout);
	assert(uc_area > 0);
	assert(good_cells > 0);
	assert(good_rcells > 0);
	assert(data_sum > 0);
	assert(rdata_sum > 0);
	assert(false);*/
	data_avg = data_sum/uc_area;
	raw_avg = data_avg/(double)good_cells;
	rdata_avg = rdata_sum/uc_area;
	rraw_avg = rdata_avg/(double)good_rcells;
	uc_avg = 0.5*(raw_avg + rraw_avg);
	//const double avg (data_sum/uc_area);
/*	//const double sum2 ( (double) lsum2 );
	//printf("data_sum %lf\n",sum);
	//printf("sum2 %lf\n",sum2);
	//printf("<sum2> %lf\n",sum2/uc_area);
	//printf("avg %lf\n",avg);
	//printf("avg2 %lf\n", avg*avg);*/
	const double avg( sum/uc_area ); //data_avg + rdata_avg
	contrast = sqrt( fabs(sum2/(double)uc_area - avg*avg ) ) / avg;
	//assert(contrast < 0.9);
	double new_merit ( 500.0*contrast );
	//DONT CALCULATE hex_avg here, the bondlength of hex_low and hex_high cannot be known
	//hex_avg = data_sum / ( (double)good_cells * 3 * (*bondlength)* (*bondlength) ); //geometric expectation of average counts on hexagons
	//double new_merit ( 10000.0*( (sum2/(double)uc_area - avg*avg ) / (avg*avg) ) ); //(std/mean)^2 can become as good as 100 for exp data
/*	//printf("rvar %lf\n", rvar);
	//const double new_merit (rvar); //measure of contrast in summed unitcell, should be independend of type of noise
	// These work well for generated Poisson Noise, but totally fail for experimental data                     
	//const double new_merit( data_sum/(variance + 0.001*data_sum) ); 
	// data_sum*idet/( stddev + 0.001 * good_cells ) //works very well, but initially small r is favored and really flat peak
	// data_sum/(variance + 0.001*data_sum) //works even better, and favours smaller r in the beginning peak is much steeper
	// data_sum/(stddev*stddev)+ 0.001*data_sum) //works great, but very shallow maximum */
	
	distort_basis(false);// cancel elliptic distortions x,y are now ideal coordinates
	static const double ca1[6] = {THIRD, TWOTHIRD, -THIRD   , TWOTHIRD, 4*THIRD,   THIRD};
	static const double ca2[6] = {THIRD, TWOTHIRD,  TWOTHIRD, -THIRD  ,   THIRD, 4*THIRD};
	const double ilen( 1.0/(3 * (*bondlength) * hel));
	double cx[6];
	double cy[6];
	for(int i(0); i < 6; ++i )
	{
		cx[i] = ilen*xd_a(ca1[i],ca2[i]);
		cy[i] = ilen*yd_a(ca1[i],ca2[i]);
	}
	const double sigmaf( 1.0/(2.35) ); //FWHM = a0;
	const double sigmax( sigmaf*ilen*xd_a(THIRD,THIRD) );
	const double sigmay( sigmaf*ilen*yd_a(THIRD,THIRD) ); //should be identical to sigmax
	const double nis( -1.0/( sigmax*sigmax + sigmay*sigmay ) );
	int *gauss( uc_gauss );
	int *rgauss( ruc_gauss );
	for(int ind(0); ind < uc_area; ++ind)
	{
		const int ix( ind % uc_res );
		const int iy( ind / uc_res );
		double aa1( ((double)ix)/uc_res );
		double aa2( ((double)iy)/uc_res );
		
		//aa1 -= floor(aa1);
		//aa2 -= floor(aa2);
		
		const double ax( ilen*xd_a(aa1,aa2) );
		const double ay( ilen*yd_a(aa1,aa2) );
		double val(0.0);
		
		for(int i(0); i < 6; ++i)
		{
			const double acx( ax - cx[i] );
			const double acy( ay - cy[i] );
			
			const double dd( acx*acx + acy*acy );
			val += ( exp( dd*nis ) );
		}
		*(gauss++) = (int) (data_min + (data_max-data_min)*val);
		*(rgauss++) = (int) (rdata_min + (rdata_max-rdata_min)*val);
	}
	const double shapeQ( find_steps() ); //find_steps() restores elliptic distortions by calling apply_offsets()
	return use_position?new_merit*shapeQ:new_merit;
}
