#ifndef SUBFRAMES_HPP
#define SUBFRAMES_HPP

#include "assert.h"
//#include "stdio.h"

int optimize_subframes(bool refresh, bool resample = true);


//tested version 10.12.2015 18:00
inline __attribute__((always_inline))
void hexagonal_trap(int &lx, int &ly, int &lz, const int N)
{
	//assert(lx+ly+lz == 0);
	//assert(N > 0);
	bool at_home;
	const int N2 ( 2*N );
	//apply periodic boundary conditions until we are home
	do
	{
		at_home = true;
		if(lx >= N)
		{
			at_home = false;
			const int s( 1 + (lx-1)/N2 );
			lx -= s*N2;
			ly += s*N;
			lz += s*N;
		}
		else if (lx < -N)
		{
			at_home = false;
			const int s( 1 - (lx+1)/N2 );
			lx += s*N2;
			ly -= s*N;
			lz -= s*N;
		}

		if(ly > N)
		{
			at_home = false;
			const int s( 1 + (ly-1)/N2 );
			lx += s*N;
			ly -= s*N2;
			lz += s*N;
		}
		else if (ly <= -N)
		{
			at_home = false;
			const int s( 1 - (ly+1)/N2 );
			lx -= s*N;
			ly += s*N2;
			lz -= s*N;
		}
		//if x and y are ok z will be fine after next shift if s==1
		if(lz >= N)
		{
			
			const int s( 1 + (lz-1)/N2 );
			at_home &= (s==1);
			lx += s*N;
			ly += s*N;
			lz -= s*N2;
		}
		else if (lz < -N)
		{
			const int s( 1 - (lz+1)/N2 );
			at_home &= (s==1);
			lx -= s*N;
			ly -= s*N;
			lz += s*N2;
		}
	} while (!at_home);
	//assert(lx < N);
	//assert(lx >= -N);
	//assert(ly <= N);
	//assert(ly > -N);
	//assert(lz < N);
	//assert(lz >= -N);
}


#endif
