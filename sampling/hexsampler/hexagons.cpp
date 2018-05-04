#include "hexagons.hpp"
#include "basis.hpp"
#include "ellipse.hpp"
#include "unitcell.hpp"
#include "subframes.hpp"
#include "graphene.hpp"
#include "correlate_target.hpp"
#include "globals.hpp"
#include "communicator.hpp"
#include "xorshift1024star.hpp"
#include <math.h>
#include <assert.h>

//private headers
void init_hexagons(void);
//returns 0 if the pixel belongs to the inital frame and -1 otherwise;
int check_interior(int q, int r);
//returns 1 if the pixel hexagons is inside e.g. != -1 
int check_hexpixel(int x, int z);
//returns home index of hit hexagon

double get_graphene_correlation(void);

//type .. 0 no defect
//type .. 1 transposes a square in the center 
void create_defect(int type);

inline __attribute__((always_inline))
int hex_ind_a1_a2(double a1, double a2);

//calculates the average intensity on valid hexagonal pixels
void set_hex_avg(void);
double eval_uc_mini(double &phs);
double get_hex_shape(void);
void translate_uc_mini(int dx, int dz);
void translate_uc_norm(double* const uc_mini_norm, double* uc_mini_norm2, const int tx, const int tz);

//center of hexagons are at corners and (1/3,1/3) and (2/3,2/3) of a subunitcell
//we could also use a simple loop, but maybe we want to test more special geometries
const int sites_a1[] = { 0, 1, 0, 1 };
const int sites_a2[] = { 0, 0, 1, 1 };
static_assert( sizeof(sites_a1)==sizeof(sites_a2), "Dimesion mismatch of sites_a1[] and sites_a2[]" );
//cartesians coordinates inside the subunitcell
double sites_x[sizeof(sites_a1)/sizeof(sites_a1[0])], sites_y[sizeof(sites_a2)/sizeof(sites_a2[0])];

static int num_MC(-1);

void init_hexagons(void)
{
	distort_basis(true);
	//memset(hexagons, -1, hexframelen*sizeof(int) ); //-1 is OxFFFF
	//memset(hex_mask, 0, hexframelen*sizeof(int) );
	for(size_t i(0); i<( sizeof(sites_a1)/sizeof(sites_a1[0])) ; ++i )
	{
		sites_x[i] = xd_p( (double)sites_a1[i], (double)sites_a2[i]); 
		sites_y[i] = yd_p( (double)sites_a1[i], (double)sites_a2[i]); 
	}
	int *hgs( hexagons );
	int *hmk( hex_mask );
	for(int qr( 0 ); qr < hexframelen; ++qr )
	{
		const int q( qr % hexImpWidth );
		const int r( qr / hexImpWidth );
		const int flg( check_interior(q,r) ); 
		*(hgs++) = flg; //0 inside and -1 outside
		*(hmk++) = (flg==0)?1:32; //1 inside and 32 outside
	}
	return;	
}

inline __attribute__((always_inline))
int check_interior( int qO, int rO)
{
	const int q( qO - hexImpWidth/2);
	const int r( rO - hexImpHeight/2);
	const double a1len ( (double)(q - (r - (r&1) )/2) );
	const double a2len ( (double)r );
	const double x ( x_p( a1len, a2len ) );
	const double y ( y_p( a1len, a2len ) );
	/*//sufficiently inclusive test for debugging
	if( (floor(x) >= -1-1.5*hel) && (ceil(x) < (double)impWidth + 1 + 1.5*hel) &&
		(floor(y) >= -1-1.5*hel) && (ceil(y) < (double)impHeight + 1 + 1.5*hel) )
		{	return 0;}*/
	
	//check that the center of the hex would be deep enough inside the image
	/*if( (floor(x) >= (1.0 + 1.5*hel) ) && ( ceil(x) < (double)impWidth - 1.0 - 1.5*hel) &&
		(floor(y) >= (1.0 + 1.5*hel) ) && ( ceil(y) < (double)impHeight - 1.0 - 1.5*hel) )
	{	return 0;}*/
	return ( (floor(x) >= (1.0 + 1.5*hel) ) && ( ceil(x) < (double)impWidth - 1.0 - 1.5*hel) &&
		(floor(y) >= (1.0 + 1.5*hel) ) && ( ceil(y) < (double)impHeight - 1.0 - 1.5*hel) ) ? 0 : -1;
}


inline __attribute__((always_inline))
int hex_ind_a1_a2(double a1, double a2)
{
	//containing unitcell
	const int a1l( (int)floor(a1) );
	const int a2l( (int)floor(a2) );
	//offset in unitcell
	const double a1f( a1-(double)a1l);
	const double a2f( a2-(double)a2l);
	//relative cartesian around unitcells
	const double x( xd_p( a1f, a2f) );
	const double y( xd_p( a1f, a2f) ); 
	//find closest site
	double d2_min( 0.0 );
	int i_min( 0 ); //index of closest site aka center of hexagon
	for(int i(0); i < 4; ++i) //sizeof is unsigned and would spoil unrolling
	{
		const double dx( x-sites_x[i] );
		const double dy( y-sites_y[i] );
		const double d2( dx*dx + dy*dy );
		if( (d2 < d2_min) || (i==0) )
		{
			d2_min = d2;
			i_min = i;
		}
	}
	const int ix( a1l + sites_a1[i_min] );
	const int iz( a2l + sites_a2[i_min] );
	const int q( ix + (iz - (iz&1) )/2 + hexImpWidth/2);
	const int r( iz + hexImpHeight/2); 
	return ( ( (q < 0) || (q >= hexImpWidth) || (r < 0) || (r >= hexImpHeight) ) ? -1 : (q + r * hexImpWidth) );
}


inline __attribute__((always_inline))
int check_hexpixel(int x, int z)
{
	const int q( x + (z - (z&1) )/2 + hexImpWidth/2) ;
	const int r( z + hexImpHeight/2);
	if( (q < 0) || (q >= hexImpWidth) ||
		(r < 0) || (r >= hexImpHeight) )
	{ return 0;}
	const int qr( q + r * hexImpWidth);
	return ( hexagons[qr] != -1  ) ? 1 : 0; 
}

void set_hex_avg(void)
{
	const int* hexfr( hexagons );
	const int* hexmsk( hex_mask );
	long good_hexes ( 0 );
	long hex_sum ( 0 );
	for(int ind(0); ind < hexframelen; ++ind)
	{
		int hval = *(hexfr++);
		int hmsk = *(hexmsk++);
		if( hmsk == 1 )
		{
			++good_hexes;
			hex_sum += hval;
		} 
	}
	hex_avg = (double)hex_sum/(double)good_hexes;	
}

double sample_hexagons(double &phs, bool skip_mini)
{
	//int trials (0);
	//do
	//{
		init_hexagons();
		const int * fr(frame);
		const bool ald(allow_light_dirt);
		const int mskWidth(impWidth/mask_scaling);
		for(int raw_ind = 0; raw_ind < framelen ; ++raw_ind)
		{
			const int val( *(fr++) );
			
			
			int dv0( 1 + (val-sample_rate/2)/sample_rate);
			
			const int x( raw_ind % impWidth); //no hex frame here
			const int y( raw_ind / impWidth); //no hex frame here
			const int flg( mask[x/mask_scaling + (y/mask_scaling) * mskWidth ] );
			
			//Dont use offset_XS and offset_YS here as the are not updated inside optimize_subframes()
			double cx( (double)x - offset_X);
			double cy( (double)y - offset_Y);
			
				int v(0);
				do
				{
					const int dv( ( (v+dv0) > val ) ? (val-v) : dv0  );
					//if(dv==0) //Hmm is that check worth the rare cases of data==0?
					//{	continue;}
					const double subx( rand_d() );
					const double suby( rand_d() );
					const double xs( cx + subx );
					const double ys( cy + suby );
					const double p1len( p1d(xs,ys) );
					const double p2len( p2d(xs,ys) );
					const int hexind( hex_ind_a1_a2(p1len,p2len) );
					if(hexind == -1) 
					{	break;	}	//break;
					hex_mask[hexind] = flg;
					if( (!ald) && (flg!=1) ) 
					{
						break;		//break;
					}
					else if( ald && (flg!=1) && (flg!=light_dirt_val) )
					{	
						break;		//break;
					}
					hexagons[hexind] += dv;
				}
				while( (v+=dv0) < val);
		}
		
		if(defect_type != 0)
		{
			create_defect(defect_type);
		}
		
		
		set_hex_avg();
		if(skip_mini)
		{	
			/*
			//nasty hacky check that the offset is correctly atop hollow site // seems to be ok
			const bool old_floating_offset (floating_offset);
			const int old_modelsize (modelsize);
			const int old_uc_mini_len( uc_mini_len );
			int *const old_uc_mini( uc_mini );
			uc_mini = nullptr;
			floating_offset = true;
			modelsize = 2 * (*bondlength);
			eval_uc_mini(phs);
			floating_offset = old_floating_offset;
			modelsize = old_modelsize;
			delete[] uc_mini;
			uc_mini_len = old_uc_mini_len;
			uc_mini = old_uc_mini;
			if( (next_offset_p1 != 0) || (next_offset_p2 != 0) )
			{
				if(!reporting)
				{	printf("Message ");}
				printf("ERROR detected wrong center while sampling hex pixels p1,p2: %.0lf,%.0lf != 0,0\n", next_offset_p1, next_offset_p2);
				fflush(stdout);
				apply_p1_p2();
				continue;
			}
			*/
			return -1000.0;
		} //signal no hex_merit was calculated and dont change uc_mini
		else
		{  return eval_uc_mini(phs);}
		
	//}
	//while(++trials < 5);
	//return -1000.0;
}	
	
double eval_uc_mini(double &phs)
{	
	//assert(*bondlength == bondlength_Mini);
	const int rH( *bondlength );
	
	uc_mini_len = 4*rH*rH;
	//assert(uc_mini_len > 0);
	
	
	//assert( uc_mini != nullptr );
	//assert(uc_mini_len == modelsize * modelsize);
	num_MC = optimize_subframes(true,false) ; // dummy optimize offsets (on 1x1 grid) but DONT call back sample_hexagons  
	//assert( sf_sum != nullptr );
	if(num_MC > 0)
	{
		delete[] uc_mini;
		uc_mini = new int[uc_mini_len];
		memcpy( uc_mini, sf_sum, uc_mini_len * sizeof(int) );
		if(uc_target == nullptr)
		{	phs = get_hex_shape();} //will suggest next p1 and p2 if floating_offset is true
		else
		{	
			
			
			moment2 = 0.0; // mark it as not calculated
			phs = 100*get_graphene_correlation();
		}
	}
	else 
	{
		if(!suppress_no_mini_cells )
		{
			if(!reporting){	printf("Message ");}
			printf("WARNING there were not any complete hexagonal cells found\n");
			fflush(stdout);
			suppress_no_mini_cells = true;
		}
		phs = -1000.0;
		//contrast = 0.0; simply keep the one from sample_unitcells
		return -1000.0;
	}
	double hl(0.0);//, hl2(0.0); //hex_low
	double hh(0.0);//, hh2(0.0); //hex_high
	//double hr(0.0);//, hr2(0.0); //ring around center
	//double hb(0.0);//, hb2(0.0);//remaining pixels on the bonds
	int hlc(0);
	int hhc(0);
	//int hrc(0);
	//int hbc(0);
	
	int rl(0); //radius around hole
	int rh(0); //radius around atom
	int rr(0); //width of ring around hole
	
	const int mcx( rH - (rH - (rH&1) )/2 );
	const int mcz( rH );
	while(  ( (1 + rl )     + (1 + rr )     + ( rh + 1) )     < ( rH + 1 ) )
	{
		if( ( (1 + (++rl) ) + (1 + rr )     + ( rh + 1) )     >= ( rH + 1 ) ) break;
		if( ( (1 + rl )     + (1 + rr )     + ( (++rh) + 1) ) >= ( rH + 1 ) ) break;
		if( ( (1 + rl )     + (1 + (++rr) ) + ( rh + 1) )     >= ( rH + 1 ) ) break;
	}
	
	double avg(0.0);
	double std(0.0);
	
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
			const int px( mcx + dx );
			const int pz( mcz + dz );
			const int qh( px + (pz - (pz&1) )/2 );
			const int rh( pz );
			const int pos( qh + rh * 2 * rH );
			const double hex_pix( (double)uc_mini[pos]/( num_MC ) );
			//simply use the last average here
			
			{
				const double dev( hex_pix-hex_avg);
				avg += dev;
				std += dev*dev;
			}
			const int xt(dx), yt(-dx-dz), zt(dz);
			int mm( abs(xt) );
			int nn( abs(xt) );
			if(abs(yt) > mm) {	mm = abs(yt);}
			if(abs(zt) > mm) {	mm = abs(zt);} 
			if(abs(yt) < nn) {	nn = abs(yt);}
			if(abs(zt) < nn) {	nn = abs(zt);}
			
			if( mm <= rl ) //central position
			{	
				hl  += hex_pix;
			    ++hlc;	
			}
			else if( mm <= rl + (rr + 1) ) //a ring around the center
			{
				//hr += hex_pix;
				//++hrc; 	
			}
			else if( mm >= 4*nn) //closer to the atom 
			{
				hh += hex_pix;
				++hhc;
			}
			/*
			else //must be inside an concave shape around the center of the bond
			{
				hb  += hex_pix;
				++hbc;	
			}
			*/
		}
	}
	
	if(hlc != 0) {	hl /= (double)(hlc);}// hl2 /= (double)(hlc);	}
	//if(hrc != 0) {  hr /= (double)(hrc);}// hr2 /= (double)(hrc);	}
	//if(hbc != 0) {  hb /= (double)(hbc);}// hb2 /= (double)(hbc); }
	if(hhc != 0) {  hh /= (double)(hhc);}// hh2 /= (double)(hhc);}
	
	const double ucm_pts(0.75*(double)uc_mini_len);
	hex_low = hl;    //approximation during runtime
	hex_high = hh;   //approximation during runtime
	avg /= ucm_pts;  //approximation during runtime
	std = sqrt( fabs(std/ucm_pts - avg*avg) );
	hex_avg += avg; 
	//can become as much ~ 0.1
	
	contrast = std/hex_avg;  //approximation during runtime
	return 500 * contrast * phs;
	///lets try to always return a positive hex_merit, if that works we can also eliminate much of the code above
	//assert(contrast < 0.9); passes this cannot cause the issue visible in ImageJ
	/*
	const double std_l( sqrt( fabs( hl2 - pow(hl,2) ) / (double)(hlc) ) );
	const double std_r( sqrt( fabs( hr2 - pow(hr,2) ) / (double)(hrc) ) );
	const double std_b( sqrt( fabs( hb2 - pow(hb,2) ) / (double)(hbc) ) );
	const double std_h( sqrt( fabs( hh2 - pow(hh,2) ) / (double)(hhc) ) );
	*/
	/*
	double next_hex_merit(-1000.0);
	
	if(uc_target == nullptr)
	{
		if( (hlc != 0) && (hhc != 0) && (hrc!=0) && (hbc!=0) )
		{
			if((hh > hb ) && (hb > hr ) && (hr > hl))
			{	
				next_hex_merit = phs * 500 * contrast;
			}
			else
			{
				next_hex_merit = 100.0 * (  (hh-hl) / (hh+hl) - 1.0 );
			}
		}
		// should become 1000 for perfect contrast realistically ~100
		// 0 for no contrast
		// <0 for wrong phase;
		else if( (hlc != 0) && (hhc != 0) )
		{	
			//assert(rH<=2);
			const double phi( (phs>0.0) ? (0.5+phs) : (0.5 + atan( phs )/M_PI) );
			next_hex_merit = 1000.0 * phi * (hh-hl) / (hh+hl);
		}
	}
	else //uc_target != nullptr
	{
		next_hex_merit = 500 * contrast * phs;
	}
	return next_hex_merit;
	*/ 
	
}


double get_hex_shape(void)
{
	
	bool first = true;
	const int numSym( use_mirror?12:6 );
	const double numSymd ( (double)numSym );
	//double best_moment (-1000.0 );
	double best_hsp( -1000.0 );
	//double best_cross( -1000.0 );
	//double best_cross_P( -1000.0 );
	//double best_shape( -1000.0 );
	double best_loc( -1000.0 );
	//double best_d_std(1000.0);
/*	
	double moment0 (-1000.0 );
	double hsp0( -1000.0 );
	double cross0( -1000.0 );
	double cross_P0( -1000.0 );
	double shape0( -1000.0 );
	double loc0( -1000.0 );
	double d_std0(1000.0);
*/	
	int best_dxm( 0 );
	int best_dzm( 0 );
	//double bg_lvl( 0.0 );
	double pts( 0.75 * (double)uc_mini_len );
	//int const * ucm(uc_mini);
	
	
	//for(int ind(0); ind < uc_mini_len; ++ind)
	//{	bg_lvl += (double)*(ucm++);}
	//const double inorm( pts/bg_lvl );
	double * uc_mini_norm( new double[uc_mini_len] );
	double * uc_mini_norm2 ( new double[uc_mini_len] );
	double * uc_mini_sym( new double[uc_mini_len] );
	double * sym_St( new double[ numSym * uc_mini_len] );
	double * cross_St ( new double[numSym] );
	
	//double * ucm_norm(uc_mini_norm);
	
	const double nMC( num_MC );
	for(int ind(0); ind < uc_mini_len; ++ind)
	{	uc_mini_norm[ind] = ( ((double)uc_mini[ind])/nMC);} //zeros at border dont matter here
	for(int i(0); i < smooth_passes; ++i)
	{	
		double * const uc_mini_rough( uc_mini_norm );//capture
		uc_mini_norm = mini_smooth(uc_mini_rough); //allocates and returns a new array
		delete[] uc_mini_rough;//release
	} 
	double n_std;
	double n_avg( get_hex_stats( uc_mini_norm, n_std) );
	//const double vol(  );
	//const double scaling(1.0/vol);
	normalize_hex( uc_mini_norm, n_avg, n_std);
	/*{	///DEBUG section
		const double std0( get_hex_std(uc_mini_norm) ); 
		assert(fabs(std0-1.0) < 0.0001);
		double n0_std;
		double n0_avg( get_hex_stats( uc_mini_norm, n0_std) );
		assert(fabs(std0-1.0) < 0.0001);
		assert(fabs(n0_avg) < 0.0001);
	}*/
	
	
	const double ivar = 1.0/( pts ); //std == 1.0 in pow(std,2)
	//distort_basis(false);
	const double hbl( hel * (*bondlength) );
	const double ir_max2( 1.0/( 3 * hbl * hbl ) );
	
	const int rH(*bondlength);
	const int mcx( rH - (rH - (rH&1) )/2);
	const int mcz( rH);
	
	const int dxm0( (int)next_offset_p1 );
	const int dzm0( (int)next_offset_p2 );
	
	for(int dxm( -rH); dxm < rH ; ++dxm)
	{
		int dzm1( -rH );
		int dzm2( rH - dxm );
		if (dxm < 0)
		{
			dzm1 = -rH - dxm;
			dzm2 = rH;
		}
		for(int dzm( dzm1); dzm < dzm2; ++dzm)
		{
			
			memset(uc_mini_sym,0,uc_mini_len * sizeof(double));
			memset(sym_St,0, numSym * uc_mini_len * sizeof(double));
			
			if(!floating_offset) //one time loop
			{
				dxm = dxm0;
				dzm = dzm0;
			}
			
			if( (dxm != 0) || (dzm != 0) )
			{	translate_uc_norm(uc_mini_norm,uc_mini_norm2, -dxm, -dzm);	} 
			else
			{	memcpy(uc_mini_norm2, uc_mini_norm, uc_mini_len * sizeof(double));}
			
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
					int px(  dx  );
					int pz(  dz  );
					int py( -px - pz );
					hexagonal_trap(px,py,pz,rH);
					px += mcx;
					pz += mcz;
					
					const int qh( px + (pz - (pz&1) )/2 );
					const int rh( pz );
					const int pos( qh + rh * 2 * rH );
					const double norm_val(	uc_mini_norm2[pos] );
					
					int rx( dx);
					int rz( dz);
					int ry(-rx -rz);
					int symNr( 0 );
					for(int r(0); r < 6; ++r)
					{
						//Rotations
						if(r > 0)
						{
							const double swap( rx );
							rx = -ry;
							ry = -rz;
							rz = -swap;
						}
						int mrx( rx ); // + dxm
						int mrz( rz ); // + dzm
						int mry( -mrx - mrz);
						hexagonal_trap(mrx,mry,mrz,rH);
						mrx += mcx;
						mrz += mcz;
						
						const int qr( mrx + ( mrz - (mrz&1) ) / 2 );
						const int rr( mrz );
						
						const int rpos( qr + rr * 2 * rH );
						uc_mini_sym[rpos] += norm_val;
						
						const int symStpos (rpos + (symNr++) * uc_mini_len); //FIXME transpose sym_St 
						sym_St[symStpos] = norm_val;
						if(use_mirror)
						{
							//axis mirroring //point mirroring has det +1 in 2D
							const int ix = rx;
							const int iz = ry;
							int mix( ix  );
							int miz( iz  );
							int miy( -mix - miz );
							
							hexagonal_trap(mix,miy,miz,rH);
							mix += mcx;
							miz += mcz;
							
							const int qi(mix + ( miz - (miz&1) ) / 2);
							const int ri(miz);
							const int ipos( qi + ri * 2 * rH );
							uc_mini_sym[ipos] += norm_val;
							
							const int symStipos (ipos + (symNr++) * uc_mini_len); //FIXME transpose sym_St
							sym_St[symStipos] = norm_val;
						}
					}
				}
			}
			for(int ind(0); ind < uc_mini_len; ++ind)
			{	uc_mini_sym[ind] /= numSymd;}
		
			double s_var( 0.0 );
			//double d_var( 0.0 );
			double cross( 0.0 );
			double moment( 0.0 );
			double loc (0.0);
			const double * const ucms( uc_mini_sym );
			memset(cross_St,0, numSym * sizeof(double) );
			
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
					const int px( dx + mcx ); 
					const int pz( dz + mcz); 
					const int qh( px + (pz - (pz&1) )/2 );
					const int rh( pz );
					const int pos( qh + rh * 2 * rH );
					const double sym_val(	ucms[pos] );
					const double norm_val(	uc_mini_norm2[pos] );
					//const double d_val( sym_val - norm_val );
					s_var += sym_val*sym_val;
					//d_var += pow(d_val, 2);
					cross += sym_val * norm_val;
					
					for( int s( 0 ); s < numSym; ++s )
					{
						const double sym_val_s ( sym_St[ pos + s * uc_mini_len] );
						cross_St[s] += sym_val_s * norm_val;
					}
					//moment should not depend directly on dxm and dzm
					///Shape and contrast are competing if data is centered around 0
					{
						const double xdp( xd_p(dx,dz) ); 
						const double ydp( yd_p(dx,dz) );
						moment += norm_val * ( xdp*xdp + ydp*ydp ) * ir_max2; // (norm_val-avg)
					}
					if( uc_target != nullptr )
					{
						const double tar_val( uc_target[pos] );
						loc += tar_val * norm_val;
					}
				}
			}
			//moment *= (10 * scaling); //as if the unsigned average would be 1.0
			//Magic Number 10.0 makes hex_shape roughly Position_Q
			moment *= (10/pts);  //normalize moment to area 
			const double s_std( sqrt( s_var/pts ) ); //- pow(s_avg,2)
			cross /= (pts * s_std); //*std == 1
			//double cross_P( 1.0 );
			
			double min_cP (1.0);
			for( int s( 0 ); s < numSym; ++s )
			{
				cross_St[s] *= ivar;
				//cross_P *= (  cross_St[s] ); //maybe we want to inspect cross_St some later time
				if( s==0 || cross_St[s] < min_cP)
				{	min_cP = cross_St[s];}
			}
			/*
			if(min_cP > 0)
			{
				cross_P = pow(cross_P,1.0/numSymd); //geometric mean of cross [0..1]
			}
			*/
			/* 
			if(min_cP < 0 ) //keep constant at low level,use purely mom2
			{
				min_cP = 0.1; 
			}
			*/ 
			const double cross_P( (1.0+min_cP)/2);
			
			//const double d_std( sqrt( d_var/pts ) ); // - pow(d_avg,2)
			//const double new_shape( (d_std>0.0) ? 1.0/d_std : 1000.0); // std == 1.0
			
			//const double hsp( cross * new_shape * moment * cross_P );
			// 10 * moment
			const double hsp( moment * cross * cross_P ); 
			if( first || ( (uc_target==nullptr) ? (hsp > best_hsp) : ( loc/pts > best_loc) ) )
			//if( first || (hsp > best_hsp) )
			{
				first = false;
				best_loc = loc/pts;
				//best_cross_P = cross_P;
				//best_shape = new_shape;
				//best_moment = moment;
				//best_cross = cross;
				//best_d_std = d_std;
				best_hsp = hsp;
				moment2 = moment;
				best_dxm = dxm;
				best_dzm = dzm;	
			}
/*			
			if( (dxm == dxm0) && (dzm == dzm0) )
			{
				loc0 = loc/pts;
				cross_P0 = cross_P;
				shape0 = new_shape;
				moment0 = moment;
				cross0 = cross;
				d_std0 = d_std;
				hsp0 = hsp;
			}
*/
			
			if(!floating_offset)//one time loop
			{	break;}
		}//dzm
		if(!floating_offset)//one time loop
		{	break;}
	}//dxm
	//They will always be zero without floating_offset
	if(floating_offset)
	{
		
		translate_uc_mini(-best_dxm,-best_dzm); //only good for display, no longer used for calculation
		//assert( p1_and_p2_applied );
		next_offset_p1 = (double)best_dxm;
		next_offset_p2 = (double)best_dzm;
		p1_and_p2_applied = false;
	}
	//distort_basis(true); //just in case we want to apply next_offset_p1/2 somewhere else
	
#if 0	
	if(reporting && has_console)
	{
		if((best_dxm != dxm0) || (best_dzm != dzm0) || (best_hsp < 0.005 * hex_shape) )
		{
			printf("loc: %lf dxm,dzm: %d,%d corr: %lf shape: %lf moment: %lf cross_P: %lf d_std: %lf hex_shape: %lf\n", 
			best_loc, best_dxm, best_dzm, best_cross, best_shape, best_moment, best_cross_P, best_d_std, best_hsp);
			///They look reasonable
			/*if((best_dxm != dxm0) || (best_dzm != dzm0))
			{
				printf("loc0: %lf dx0,dz0: %d,%d corr0: %lf shape0: %lf moment0: %lf cross_P0: %lf d_std0: %lf hex_shape0: %lf\n", 
				loc0, dxm0, dzm0, cross0, shape0, moment0, cross_P0, d_std0, hsp0);
			}*/
			///These seem to be correct
			printf("tmp X0,Y0: %lf,%lf tlt,hel: %lf,%lf phi,exc: %lf,%lf merit: %lf hex_merit: %lf\n",
				offset_X, offset_Y, tilt, hel*bondlength_Mini/bondlength_CC , phi, excent, merit, hexagonal_merit);
			///Also fine
			//printf("num_MC: %d floating_offset: %s ruling_merit: %d\n", num_MC, floating_offset?"true":"false", ruling_merit);
			/*double norm_std = get_hex_std( uc_mini_norm ); 
			double norm2_std = get_hex_std( uc_mini_norm2 );
			double sym_std = get_hex_std(uc_mini_sym);
			printf("norm_std: %lf norm2_std: %lf sym_std: %lf\n",norm_std, norm2_std, sym_std);
			for(int s(0); s < numSym; ++s)
			{
				const double std_ss( get_hex_std(sym_St+s*uc_mini_len) );
				printf("std_%d: %lf ", s, std_ss);
			}
			puts("\n");*/
			fflush(stdout);
			if(best_hsp < 0.005 * hex_shape)
			{	exit(0);}
		}
	}
#endif	
	delete[] uc_mini_norm;
	delete[] uc_mini_norm2;
	delete[] uc_mini_sym;
	delete[] cross_St;
	delete[] sym_St;
	//MAGIC number 10.0 makes good target_correlation ~ 1.0
	///we now use match to entire graphene sheet for bad data that requires a target at all
	//const double target_cc( (uc_target ? (10.0*correlate_target(-best_dxm,-best_dzm)) : 1.0 )); //need to test the signs here but rarely !=0
	return best_hsp;// * target_cc;
	//return best_hsp;
}


double* mini_smooth( double * const rough)
{
	double *smooth = new double[uc_mini_len];
	memset(smooth, 0 , uc_mini_len * sizeof(double) );
	
	const bool median_filter(false); //tested but actually better to median filter the input once beforehand
	const int rH(*bondlength);
	const int mcx( rH - (rH - (rH&1) )/2 );
	const int mcz( rH);
	
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
			const int px0( mcx + dx );
			const int pz0( mcz + dz );
			const int qh0( px0 + (pz0 - (pz0&1) )/2 );
			const int rh0( pz0 );
			const int pos0( qh0 + rh0 * 2 * rH );
			double med( rough[pos0] );
			double sum( med ); 
			double smaller[7],bigger[7];
			int sp( -1 ), bp( -1 ), ep( 0 );
			//visit all 6 neighbors
			for(int n(0); n < 6; ++n)
			{
				int pxn( dx + cube_neighbors_X[n] );
				int pzn( dz + cube_neighbors_Z[n] );
				int pyn( -pxn - pzn);
				hexagonal_trap(pxn,pyn,pzn,rH);
				pxn += mcx;
				pzn += mcz;
				const int qhn( pxn + (pzn - (pzn&1) )/2 );
				//const int rhn( pzn );
				const int posn( qhn + pzn * 2 * rH );
				const double val = rough[posn];
				//simple sum for averaging
				sum += val;
				///median algorithm for small sets N = 7
				if(median_filter)
				{
					if(val > med)
					{
						bigger[++bp] = val;
					}
					else if (val < med)
					{
						smaller[++sp] = val;
					}
					else
					{
						++ep;//no need to store identical values
					}
					//check if we are unbalanced
					if( (abs(sp - bp) - ep) > 0 )
					{
						//easy case just push one of the equal above or below	
						if(ep > 1) 
						{
							//less smaller values
							if(sp > bp)
							{
								--ep;
								smaller[++sp] = med;
							}
							else //less bigger values
							{
								--ep;
								bigger[++bp] = med;
							}
						}
						else //single outdated median
						{
							//less smaller values
							if(sp < bp)
							{
								smaller[++sp] = med;
								int	minbi(0);	
								double minb = bigger[0];
								
								for(int i(1); i <= bp; ++i)
								{
									if(bigger[i] < minb)
									{
										minbi = i; //rember index of min bigger value
										minb = bigger[i]; //min bigger value 
									}
								}
								if(minbi < bp)
								{
									//move all behind the min one to the left
									memmove(bigger+minbi,bigger+minbi+1,(bp-minbi)*sizeof(int));
								}
								--bp; //and "forget" the rightmost
								med = minb;
							}
							else //less bigger values
							{
							
								bigger[++bp] = med;
								int	maxsi(0);	
								double maxs = smaller[0];
								
								for(int i(1); i <= sp; ++i)
								{
									if(smaller[i] > maxs)
									{
										maxsi = i; //rember index of max smaller value
										maxs = smaller[i]; //max smaller value 
									}
								}
								
								if(maxsi < sp)
								{
									//move all behind the max one to the left
									memmove(smaller+maxsi,smaller+maxsi+1,(sp-maxsi)*sizeof(int));
								}
								--sp; //and "forget" the rightmost
								med = maxs;
							}
						}
					}
				
				}//endif median filter
			}	
			//write averaged rough[] to smooth[]
			if(median_filter)
			{
				smooth[pos0] = med;
			}
			else
			{
				smooth[pos0] = (sum/7.0);
			}
		}
	}
	
	return smooth;
}

void translate_uc_mini(int tx, int tz)
{
	const int rH(*bondlength);
	const int mcx( rH - (rH - (rH&1) )/2 );
	const int mcz( rH);
	
	int *mini_clone = new int[uc_mini_len];
	memcpy(mini_clone, uc_mini, uc_mini_len * sizeof(int) );
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
			const int px0( mcx + dx );
			const int pz0( mcz + dz );
			const int qh0( px0 + (pz0 - (pz0&1) )/2 );
			const int rh0( pz0 );
			const int pos0( qh0 + rh0 * 2 * rH );
			int pxt( dx + tx );
			int pzt( dz + tz );
			int pyt( -pxt - pzt );
			hexagonal_trap(pxt,pyt,pzt, rH);
			pxt += mcx;
			pzt += mcz;
			const int qht( pxt + (pzt - (pzt&1) )/2 );
			const int rht( pzt );
			const int post( qht + rht * 2 * rH );
			uc_mini[post] = mini_clone[pos0];
		}
	}
	delete[] mini_clone;
	
	return;
}

void translate_uc_norm(double* const uc_mini_norm, double* uc_mini_norm2, const int tx, const int tz)
{
	const int rH(*bondlength);
	const int mcx( rH - (rH - (rH&1) )/2 );
	const int mcz( rH);
	memset(uc_mini_norm2, 0, uc_mini_len * sizeof(double)); //make sure areas outside the hexagon are 0
	//double norm(0.0);
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
			const int px0( mcx + dx );
			const int pz0( mcz + dz );
			const int qh0( px0 + (pz0 - (pz0&1) )/2 );
			const int rh0( pz0 );
			const int pos0( qh0 + rh0 * 2 * rH );
			int pxt( dx + tx );
			int pzt( dz + tz );
			int pyt( -pxt - pzt );
			hexagonal_trap(pxt,pyt,pzt, rH);
			pxt += mcx;
			pzt += mcz;
			const int qht( pxt + (pzt - (pzt&1) )/2 );
			const int rht( pzt );
			const int post( qht + rht * 2 * rH );
			//assert(uc_mini_norm[pos0] != 0);
			//assert(uc_mini_norm2[post] == 0);
			uc_mini_norm2[post] = uc_mini_norm[pos0];
			//norm += ( uc_mini_norm2[post] = uc_mini_norm[pos0] );
		}	
	}
	//assert( fabs(norm) < 0.00001 );
	
}

//Mean is assumed to be zero
double get_hex_std(double * const uc_mini_norm)
{
	const int rH(*bondlength);
	const int mcx( rH - (rH - (rH&1) )/2 );
	const int mcz( rH);
	double var(0.0);
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
			const int px0( mcx + dx );
			const int pz0( mcz + dz );
			const int qh0( px0 + (pz0 - (pz0&1) )/2 );
			const int rh0( pz0 );
			const int pos0( qh0 + rh0 * 2 * rH );
			const double ucmn0( uc_mini_norm[pos0] );
			var += ucmn0*ucmn0;
		}
	}
	
	return sqrt(var/3.0)/(double)rH;
}

double get_hex_stats(const double * const uc_mini_norm, double& std)
{
	const int rH(*bondlength);
	const int mcx( rH - (rH - (rH&1) )/2 );
	const int mcz( rH);
	double var(0.0);
	double avg(0.0);
	double pts( 3 * rH * rH );
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
			const int px0( mcx + dx );
			const int pz0( mcz + dz );
			const int qh0( px0 + (pz0 - (pz0&1) )/2 );
			const int rh0( pz0 );
			const int pos0( qh0 + rh0 * 2 * rH );
			const double val( uc_mini_norm[pos0] );
			avg += val;
			var += val*val;	
		}
	}
	avg /= pts;
	std = sqrt(fabs( var/pts - avg*avg ));
	return avg;
}

double normalize_hex(double * const uc_mini_norm, const double avg, const double std)
{
	const int rH(*bondlength);
	const int mcx( rH - (rH - (rH&1) )/2 );
	const int mcz( rH);
	double vol = 0.0;
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
			const int px0( mcx + dx );
			const int pz0( mcz + dz );
			const int qh0( px0 + (pz0 - (pz0&1) )/2 );
			const int rh0( pz0 );
			const int pos0( qh0 + rh0 * 2 * rH );
			const double val( uc_mini_norm[pos0] );
			const double nval( (val-avg)/std );
			uc_mini_norm[pos0] = nval;
			vol += fabs(nval);
		}
	}
	return vol;
}


bool apply_p1_p2(void)
{
	if( (!p1_and_p2_applied) && ( (next_offset_p1 != 0) || (next_offset_p2 != 0) ) )
	{
		//distort_basis(false); ///at least with 6 pix per a0 It does not matter at all up to 20%
		const double offx = x_p( next_offset_p1, next_offset_p2);
		const double offy = y_p( next_offset_p1, next_offset_p2);
		
		//set_offset(offx,offy);
		offset_X = offx;
		offset_Y = offy;
		
		if(reporting)
		{
			printf("p1,p2: %.0lf,%.0lf -> X0,Y0:  %lf,%lf\n",next_offset_p1,next_offset_p2,offset_X,offset_Y);
			fflush(stdout);
		}
		//translate_uc_mini(-next_offset_p1,-next_offset_p2);
		cancel_p1_p2();
		//distort_basis(true); ///at least with 6 pix per a0 It does not matter at all up to 20%
		return true;
	
	}
	p1_and_p2_applied = true;
	return false;			
}

void cancel_p1_p2(void)
{
	next_offset_p1 = 0.0;
	next_offset_p2 = 0.0;
	p1_and_p2_applied = true;				
}

void create_defect(int type)
{
	const int hiw2(hexImpWidth/2);
	const int hih2(hexImpHeight/2);
	const int rr(*bondlength * 3 / 2 );
	const int mark(1); //scaling for modified values
	switch(type)
	{
		case 1://transpose central square
		for(int x = hiw2-rr; x < hiw2+rr; ++x )
		{
			for(int y = hih2-rr; y < hih2+rr; ++y )
			{
				const int sr( x + y * hexImpWidth );
				const int dt( y + x * hexImpWidth );
				if(sr <= dt)
				{
					const int hv( hexagons[sr] );
					const int hm( hex_mask[sr] );
					hexagons[sr] = mark*hexagons[dt];
					hex_mask[sr] = hex_mask[dt];
					hexagons[dt] = mark*hv;
					hex_mask[dt] = hm;
				}		
			}
		}
		break;
		case 2://turn central square inside out
		for(int x = hiw2-rr; x < hiw2+rr; ++x )
		{
			int x2 ( hexImpWidth-x - ( (x < hiw2) ? rr : -rr) );
			/*
			if(x < hiw2)
			{
				x2 = hexImpWidth - rr - x;
			}
			else
			{
				x2 = hexImpWidth + rr - x;
			}
			*/
			for(int y = hih2-rr; y < hih2+rr; ++y )
			{
				const int sr( x + y * hexImpWidth );
				int y2( hexImpHeight-y - ( (y < hih2) ? rr : -rr ) );
				/*
				if(y < hih2)
				{
					y2 = hexImpHeight - rr - y;
				}
				else
				{
					y2 = hexImpHeight + rr - y;
				}
				*/
				
				const int dt( x2 + y2 * hexImpWidth );
				if(sr <= dt)
				{
					const int hv( hexagons[sr] );
					const int hm( hex_mask[sr] );
					hexagons[sr] = mark*hexagons[dt];
					hex_mask[sr] = hex_mask[dt];
					hexagons[dt] = mark*hv;
					hex_mask[dt] = hm;
				}		
			}
		}
		break;
		case 0: //do nothinng
		break;
		//falltrough
		default:
		break;	
	}
	return;
}


double get_graphene_correlation()
{
	int dx( 0 );
	int dz( 0 );
	const double graphene_match = match_graphene(dx,dz);	
	//may become useful one day
	if( floating_offset && ( (dx!=0) || (dz!=0) ) )
	{
		translate_uc_mini(-dx,-dz); 
		next_offset_p1 = (double)dx;
		next_offset_p2 = (double)dz;
		p1_and_p2_applied = false;	
	}
	return graphene_match;
}
