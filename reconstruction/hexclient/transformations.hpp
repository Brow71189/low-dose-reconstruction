#ifndef TRANSFORMATIONS_HPP
#define TRANSFORMATIONS_HPP


//takes q,r coordinates and returns qt+rt*modelsize
int transform_qr(const int q, const int r, const int lp, const int rot, const int mirror);

//takes qt,rt coordinates and returns q+r*modelsize
int inv_transform_qr(const int q, const int r, const int lp, const int rot, const int mirror);

//takes an indexd hex pixel position and returns qt+rt*modelsize
int transform_hp(const int hp, const int lp, const int rot, const int mirror);

//takes hp and returns 
int inv_transform_hp(const int hp, const int lp, const int rot, const int mirror);

//takes q,r index and numbered symmetry and returns qt+rt*modelsize
int transform_qrsym(const int q, const int r, const int sym);

//takes  qt, rt and numbered symmetry and returns q+r*modelsize
int inv_transform_qrsym(const int q, const int r, const int sym);
/*
extern unsigned short *hextrap;
void init_hextrap();
*/
//ugly but we have to do that in lining in a header 
/* Wraps the cube hex position (&dx,&dy,&dz) back onto
 * a truncated hex around (0,0,0) with radius N
 */
inline __attribute__((always_inline))
void hexagonal_trap(int &lx, int &ly, int &lz, const int N)
{
	bool at_home(true);
	//apply periodic boundary conditions until we are home
	do
	{
		at_home = true;
		if(lx >= N)
		{
			at_home = false;
			lx -= 2*N;
			ly += N;
			lz += N;
		}
		else if (lx < -N)
		{
			at_home = false;
			lx += 2*N;
			ly -= N;
			lz -= N;
		}

		if(ly > N)
		{
			at_home = false;
			lx += N;
			ly -= 2*N;
			lz += N;
		}
		else if (ly <= -N)
		{
			at_home = false;
			lx -= N;
			ly += 2*N;
			lz -= N;
		}
		//if x and y are ok z will be fine after next shift
		if(lz >= N)
		{
			//at_home = false;
			lx += N;
			ly += N;
			lz -= 2*N;
		}
		else if (lz < -N)
		{
			//at_home = false;
			lx -= N;
			ly -= N;
			lz += 2*N;
		}

	} while (!at_home);
}

#endif
