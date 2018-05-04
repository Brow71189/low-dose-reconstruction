#include <cstdlib>
#include <climits>
#include "xorshift1024star.hpp"


unsigned long states[16];
const double range_H( 1.0/((double)ULONG_MAX) ); //prceision loss of ~ 10bit approx 1000 identical values in a period
int p_xor;

void init_xorshift1024star(void)
{		
		for(int i = 0; i < 16; ++i)
		{	//we only use the lower 16 bit of standard rand();
			unsigned long p0 = ( (unsigned long)(rand() & 0x0000FFFF) );
			unsigned long p1 = ( (unsigned long)(rand() & 0x0000FFFF) ) << 16;
			unsigned long p2 = ( (unsigned long)(rand() & 0x0000FFFF) ) << 32;
			unsigned long p3 = ( (unsigned long)(rand() & 0x0000FFFF) ) << 48;
			
			states[i] = ( p0 | p1 | p2 | p3 );
		}
		while( 0 != (p_xor = rand()%16) ){};
}

double rand_d(void)
{
	long s0( states[ p_xor ] );
	long s1( states[ p_xor = ( (p_xor+1) & 15) ] );
	s1 ^= (s1 << 31);
	s1 ^= (s1 >> 11);
	s0 ^= (s0 >> 30);
	return ( (double)  ( ( states[p_xor] = (s0 ^ s1) ) * (unsigned long)(1181783497276652981)  ) ) * range_H;
}
