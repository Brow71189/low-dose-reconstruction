#include "globals.hpp"
#include "subframes.hpp"
#include <math.h>
//#include <assert.h>

double correlate_target(const int dx, const int dz)
{
	//const int* hexfr( hexagons );
	const int* hmsk( hex_mask );
	long good_hexes( 0 );
	long hex_sum( 0 );
	long hex_sum2( 0 );
	double tar_sum( 0.0 );
	double tar_sum2( 0.0 );
	double cross( 0.0 );
	//const int rH( modelsize/2 );
	//const int mcx( rH - (rH - (rH&1) )/2 );
	//const int mcz( rH);
	/*
	double hl(0.0);// hex_low
	double hh(0.0);// hex_high
	//double hr(0.0);// ring around center
	//double hb(0.0);// remaining pixels on the bonds
	int hlc(0);
	int hhc(0);
	//int hrc(0);
	//int hbc(0);
	int rl(0); //radius around hole
	int rh(0); //radius around atom
	int rr(0); //width of ring around hole
	while(  ( (1 + rl )     + (1 + rr )     + ( rh + 1) )     < ( rH + 1 ) )
	{
		if( ( (1 + (++rl) ) + (1 + rr )     + ( rh + 1) )     >= ( rH + 1 ) ) break;
		if( ( (1 + rl )     + (1 + rr )     + ( (++rh) + 1) ) >= ( rH + 1 ) ) break;
		if( ( (1 + rl )     + (1 + (++rr) ) + ( rh + 1) )     >= ( rH + 1 ) ) break;
	}
	*/
	
	for(int ind(0); ind < hexframelen; ++ind)
	{
		if( *(hmsk++) == 1 ) //only consider graphene
		{
			
			const int q( ind % hexImpWidth - hexImpWidth/2 ); //
			const int r( ind / hexImpWidth - hexImpHeight/2 ); //
			const int x( q - (r - (r&1) )/2 + dx);
			const int z( r + dz);
			const int qh1( x + (z - (z&1) )/2 );
			const int rh1( z );
			if( (qh1 >= 0) && (qh1 < hexImpWidth) &&
			    (rh1 >= 0) && (rh1 < hexImpHeight) )
			{    
				const int hval( hexagons[ind] ); //typically only 10% to 25% of enire bounding box are graphene
				const int pos1( qh1 + rh1 * hexImpWidth ); //by dx and dz shifted position
				const int pos0( (hex_trap[pos1]) ); 			//equivalent to pos1 but inside uc_target
				const double tval( uc_target[pos0] );
				++good_hexes;
				hex_sum += (double)hval;
				hex_sum2 += (double)(hval*hval);
				tar_sum += tval;
				tar_sum2 += tval*tval;
				cross += tval * (double)hval;
#if 0		
				int mm( abs(x) );
				int nn( abs(x) );
				if(abs(y) > mm) {	mm = abs(y);}
				if(abs(z) > mm) {	mm = abs(z);} 
				if(abs(y) < nn) {	nn = abs(y);}
				if(abs(z) < nn) {	nn = abs(z);}
				
				if( mm <= rl ) //central position
				{	
					hl  += hval;
					++hlc;	
				}
				/*else if( mm <= rl + (rr + 1) ) //a ring around the center
				{
					hr += hval;
					++hrc; 	
				}*/
				else if( (mm >= 4*nn) && (mm > rl + (rr + 1) ) ) //closer to the atom 
				{
					hh += hval;
					++hhc;
				}
				/*else //must be inside an concave shape around the center of the bond
				{
					hb  += hval;
					++hbc;	
				}*/
#endif		
			}
		} 
	}
	/*
	if(hlc != 0) {	hl /= (double)(hlc);}// hl2 /= (double)(hlc);	}
	//if(hrc != 0) {  hr /= (double)(hrc);}// hr2 /= (double)(hrc);	}
	//if(hbc != 0) {  hb /= (double)(hbc);}// hb2 /= (double)(hbc); }
	if(hhc != 0) {  hh /= (double)(hhc);}// hh2 /= (double)(hhc);}
	
	hex_low = hl;    //approximation during runtime
	hex_high = hh;   //approximation during runtime
	*/
	
	const double h_avg( (double)hex_sum/(double)good_hexes );
	const double h_std( sqrt(fabs(hex_sum2/(double)good_hexes - h_avg*h_avg )) );
	const double t_avg( (double)tar_sum/(double)good_hexes );
	const double t_std( sqrt(fabs(tar_sum2/(double)good_hexes - t_avg*t_avg )) );
	return  (cross/(double)good_hexes - h_avg*t_avg) / (h_std * t_std);
	
	//contrast = (hex_high-hex_low) / (hex_low + hex_high + hex_avg); //approximation 
	
	//contrast = h_std/h_avg;
	//return  phs * 100 * contrast;					  	
}
