#include "assert.h"
#include "globals.h"
#include "transformations.hpp" //inline hexagonal trap
#include "hexlikelyhood.hpp" 

int transform_qr(const int q, const int r, const int lp, const int rot, const int mirror)
{
	//odd-r -> cube 
	const int x(q - (r - (r&1))/2 );
	const int z(r);
	
	//offset around center
	int xs(x - globalMx);
	int zs(z - globalMz);
	
	//translation by lattice point
	const int xt( globallpX[lp] );
	const int zt( globallpZ[lp] );
	xs += xt;
	zs += zt;
	int ys = (-xs - zs);
	
	//rotation
	if(rot > 0)
	{
		int rotxyz[3] = {xs, ys, zs};
		xs = rotxyz[rot % 3];
		ys = rotxyz[(1+rot) % 3];
		if(rot&1) //positive and odd number
		{
			xs = -xs;
			ys = -ys;
		}
		zs = (-xs - ys);
	}
	//mirroring
	if(mirror&1)
	{
		//axis mirror //point inversion is no mirror in 2D
		zs = ys;
		ys = (-xs - zs);
	}

	//periodic hexagon
	
	
	
	hexagonal_trap(xs, ys, zs, globalMz);
	
	//undo offset by center
	xs += globalMx;
	zs += globalMz;
	
	const int qt = xs + (zs - (zs&1))/2;
	//int rt = zs;
	//assert(q+r*globalModelSize == inv_transform_qr(qt, zs, lp, rot, mirror) );//passes
	return (qt + zs * globalModelSize); //zs = rt
}

int inv_transform_qr(const int q, const int r, const int lp, const int rot, const int mirror)
{
	//odd-r -> cube 
	const int x(q - (r >> 1));
	const int z(r);
	
	//offset around center
	int xs(x - globalMx);
	int zs(z - globalMz);
	int ys = (-xs - zs);
	
	//mirroring //its own inverse
	if(mirror&1)
	{
		//axis mirror point inversion is no mirror in 2D
		zs = ys;
		ys = (-xs - zs);
	}
	
	//rotation
	if(rot > 0)
	{
		const int invrot (6- (rot % 6) ); //inv rotation
		const int rotxyz[3] = {xs, ys, zs}; 
		xs = rotxyz[invrot % 3];
		ys = rotxyz[(1+invrot) % 3];
		if(invrot&1) //odd number
		{
			xs = -xs;
			ys = -ys;
		}
		zs = (-xs - ys);
	}
	
	
	//translation by lattice point
	const int xt( globallpX[lp] );
	const int zt( globallpZ[lp] );
	xs -= xt; //negate translation
	zs -= zt; //negate translation
	ys = (-xs - zs);
	
	//periodic hexagon
	hexagonal_trap(xs, ys, zs, globalMz);
	
	//undo offset by center
	xs += globalMx;
	zs += globalMz;
	
	const int qt = xs + (zs >> 1);
	//int rt = zs;
	
	return (qt + zs * globalModelSize); //zs = rt
}

int transform_hp(const int hp, const int lp, const int rot, const int mirror)
{
	const int ind = globalhexPixels[hp];
	const int q = ind % globalModelSize;
	const int r = ind / globalModelSize;
	const int result = transform_qr( q, r, lp, rot, mirror);
	/*
	const int qt = result % globalModelSize;
	const int rt = result / globalModelSize;
	
	int check = inv_transform_qr( qt, rt, lp, rot, mirror);
	assert(check == ind);
	assert(hp == pixposG[ind]);
	*/ 
	return result;
}

int inv_transform_hp(const int hp, const int lp, const int rot, const int mirror)
{
	const int ind = globalhexPixels[hp];
	const int q = ind % globalModelSize;
	const int r = ind / globalModelSize;
	return inv_transform_qr( q, r, lp, rot, mirror);
}


int transform_qrsym(const int q, const int r, const int sym)
{
	const int lp = sym / 12;
	const int rs = sym % 12;
	const int rot = (rs >> 1);
	const int mirror = (rs&1);
	return transform_qr( q, r, lp, rot, mirror);
}

int inv_transform_qrsym(const int q, const int r, const int sym)
{
	const int lp = sym / 12;
	const int rs = sym % 12;
	const int rot = (rs >> 1);
	const int mirror = (rs&1);
	return inv_transform_qr( q, r, lp, rot, mirror);
}
/*
unsigned short *hextrap = 0;

void init_hextrap()
{
	static bool first_time = true;
	assert(first_time);
	first_time = false;
	
	hextrap = new unsigned short[9*globalModelArea];
	
	int rH = 3*globalModelSize/2;
	int dH = 3*globalModelSize;
	int H = globalModelSize/2;
	
	for(int(r); r<dH; ++r)
	{
		for(int(q); q < dH; ++q )
		{
			//odd-r -> cube 
			int x(q - (r >> 1));
			int z(r);
			
			//offset around center
			int xs(x - rH);
			int zs(z - rH);
			int ys(-xs-zs);
			
			hexagonal_trap(xs, ys, zs, H);
			// undo offset by center
			xs += rH;
			ys += rH;
			
			int qt (xs + (zs >> 1) );
			int rt (zs);
			
			hextrap[q+r*dH] = (qt+rt*H);
		}
		
	}
}
*/

