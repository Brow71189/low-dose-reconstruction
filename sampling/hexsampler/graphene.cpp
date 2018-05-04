#include "graphene.hpp"
#include "hexagons.hpp"
#include "basis.hpp"
#include "ellipse.hpp"
#include "unitcell.hpp"
#include "subframes.hpp"
#include "correlate_target.hpp"
#include "globals.hpp"
#include "communicator.hpp"
#include "xorshift1024star.hpp"
#include <math.h>
#include <assert.h>

inline __attribute__((always_inline))
double get_match(int dx, int dz);


bool init_graphene()
{
	if(uc_target == nullptr)
	{	return false;} //maybe issue some error or warning?
	/*double min = 1;
	double max = -1;
	double sum = 0;
	double sum2 = 0;
	int npts = 0;*/
	int * old_bondlength = bondlength;
	const int old_uc_mini_len = uc_mini_len;
	
	bondlength = &bondlength_Mini;
	uc_mini_len = 4 * bondlength_Mini * bondlength_Mini;
	double std;
	double avg = get_hex_stats(uc_target,std);
	//printf("uc_target: avg: %lf   std: %lf\n", avg,std);
	
	int *gr (large_graphene);
	//bool always_zero = true;
	for(int ind = 0; ind < large_hexframelen; ++ind)
	{
		/*
		const int q = ind%large_hexImpWidth - large_hexImpWidth/2;
		const int r = ind/large_hexImpWidth - large_hexImpHeight/2;
		int x ( q - (r - (r&1) )/2 );
		int z = r;
		int y(-x-z);
		hexagonal_trap(x,y,z,bondlength_Mini);
		const int tq( x + (z-(z&1))/2 + bondlength_Mini);
		const int tr =  z + bondlength_Mini;
		const int tind( tq + 2 * tr * bondlength_Mini);
		*/
		const int tind(large_trap[ind]);
		//assert( (tind >= 0) && (tind < 4*bondlength_Mini*bondlength_Mini) );
		//always_zero &= (uc_target[tind] == 0);
		*gr++ = (int)(16*(uc_target[tind]-avg)/std + 0.5);
		/*
		if( (ind==0) || (val < min))
		{	min = val;}
		if( (ind==0) || (val > max))
		{	max = val;}
		*/
		//uc_target is normalized to avg=0 and sigma=1
		//graphene has usually a spread of 4 sigma 
		//so graphene values should be between -16 and +16
		/*++npts;
		sum+=val;
		sum2+= val*val;*/
	
	
	}
	//assert(!always_zero);
	//assert(max != min);
	//printf("graphene min: %lf max:%lf\n",min,max);
	//sum/=npts;
	//sum2/=npts;
	//const double stdg(sqrt(sum2-sum*sum));
	//printf("graphene avg: %lf std: %lf\n", sum, stdg);
	
	uc_mini_len = old_uc_mini_len;
	bondlength = old_bondlength;
	
	return true;
}

double match_graphene(int &deltaX, int &deltaZ)
{
	double best_match = -1.1;
	double worst_match = 1.1;
	const int rH(bondlength_Mini);
	int pts = 0;
	if(!floating_offset)
	{
		best_match = get_match(deltaX = 0,deltaZ = 0);
	}
	else
	{
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
				const double match( get_match(dx,dz) );
				++pts;
				//assert( (match >= -1.0001) && (match <= 1.0001) );
				if( match > best_match )
				{
					best_match = match;
					deltaX = dx;
					deltaZ = dz;
				}
				if(match < worst_match)
				{
					worst_match = match;
				}
			}
		}
	}
	//worst match would never be -1 so lets try best - 2*worst
	return (best_match >= 0) ? (best_match*best_match) :-pow(-best_match,0.5); //experimental
	
	//differences are a bit noisy
	//return best_match * (best_match - worst_match); //should typically be between 0 and 2
}

double get_match(int dx, int dz)
{
	int npts( 0 );
	double sumG( 0.0 );
	double sumGG( 0.0 );
	double sumT( 0.0 );
	double sumTT( 0.0 );
	double sumGT( 0.0 );
	for(int indg(0); indg < large_hexframelen; ++indg)
	{
		if( hex_mask[indg] == 1 )
		{
			const int q = indg%large_hexImpWidth - large_hexImpWidth/2;
			const int r = indg/large_hexImpWidth - large_hexImpHeight/2;
			const int x ( q - (r - (r&1) )/2 - dx);
			const int z (r - dz);
			const int qt( x + (z-(z&1))/2 + large_hexImpWidth/2);
			const int rt( z + large_hexImpHeight/2);
			if( (qt>=0) && (qt < large_hexImpWidth) &&
				(rt>=0) && (rt < large_hexImpHeight))
			{	
				//Hmm, we might use smaller graphene and large_trap instead
				const int indt( qt + rt * large_hexImpWidth);
				const double gval( large_graphene[indt] ); 
				const double tval( large_hexagons[indg] );
				++npts;
				sumG  += gval;
				sumGG += gval*gval;
				sumT  += tval;
				sumTT += tval * tval;
				sumGT += gval * tval;
			}
		}
	}
	//assert(npts>0);
	const double  crG( (sumGG - sumG*sumG/npts)/npts );
	const double  crT( (sumTT - sumT*sumT/npts)/npts );
	const double crGT( (sumGT - sumG*sumT/npts)/npts );
	const double corr( crGT*pow(crG*crT,-0.5) ); 
	//printf("MESSAGE: sumG: %lf  sumGG: %lf sumT: %lf sumTT: %lf sumGT: %lf\n", sumG, sumGG, sumT, sumTT, sumGT); 
	//printf("MESSAGE: varG: %lf  varT: %lf varGT: %lf corr: %lf npts: %d\n", crG, crT, crGT, corr, npts);
	
	return corr;
	
	
}


