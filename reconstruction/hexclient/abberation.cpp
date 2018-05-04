#include "abberation.hpp"
#include "string.h"
#include "globals.h"
#include "fast_cache.h"
#include "transformations.hpp"
#include "assert.h"
#include "stdio.h"
#include "math.h"

int *corrmodels = nullptr;
//int *ncrrmodels = nullptr;
short *diffmodel = nullptr;
unsigned short *lbeam = nullptr;
//unsigned short *nsmodels = nullptr;
unsigned short *rmodels = nullptr;
//unsigned short *nrmodels = nullptr;

double *rotor = nullptr;
int *corners = nullptr;

bool luse_beam = false;
bool lany_changes = false;

int	brH;
int	bcq;
int	bcr;
int	bcx;
int	bcz;
int	bcy;

void make_corrmodels(const int fm, const int lm );
void build_corrmodels(const int fm, const int lm );
void apply_corrmodels( unsigned short *smodels, const int fm, const int lm );
int setBeamCenter(const int beamsize);
void init_rotor(void);
void cast_shadow(unsigned short *smodel, unsigned short *model);
void addto_rotor(int qr, int x, int z, int shadow, int offset, double share);
void roll_out_diffs(unsigned short *rmodel0);

//integrates scattered intensity from beamprobe at dx0,dz0 in rmodel
int get_corrval(	const int ind0, const int modelsize,
						const int beamsize, const int beamarea, const int beamradius,
						unsigned short * const rmodel);

inline __attribute__((always_inline))
void put_beam(	const int ind1, const int delta, const int modelsize,
const int beamsize, const int beamarea, const int beamradius,
				int *const crmodel);


int setBeamCenter(const int beamsize)
{
	brH = beamsize/2;
	bcq = brH;
	bcr = brH;
	bcx = bcq - (bcr - (bcr&1) ) / 2;
	bcz = bcr;
	bcy = -bcx - bcz;
	return (bcq + bcr * beamsize);
}


void close_smodels()
{
	delete[] diffmodel;
	diffmodel = nullptr;
	delete[] rmodels;
	rmodels = nullptr;
	delete[] rotor;
	rotor = nullptr;
	delete[] corners;
	corners = nullptr;
	
	if(luse_beam)
	{
		delete[] corrmodels;
		corrmodels = nullptr;
		//delete[] ncrrmodels;
		delete[] lbeam;
		lbeam = nullptr;
		luse_beam = false;
	}
}



//smodels are inv_transformed models that can be compared to non transformed data
void init_smodels( unsigned short *smodels, unsigned short *models, unsigned short *beam, bool use_beam)
{
	//printf("init_models entered\n");
	//fflush(stdout);
	const int modelarea( globalModelArea );
	const int modelnum( globalModelNum );
	const int shadownum( globalShadowNum );
	const int hpMax( globalhpMax );
	const int smdlsize(12 * shadownum * modelnum * modelarea);
	
	luse_beam = use_beam;
	
	//unsigned short *model( models );
	
	
	//assert(smodel != nullptr);
	//assert( model != nullptr);
	//assert(globalhexPixels != nullptr);
	//assert(inv_mirrotG != nullptr);
	//printf("local initializers passed\n");
	//fflush(stdout);
	
	if(rotor == nullptr)
	{
		rotor = new double[4*shadownum*modelarea];
		corners = new int[4*shadownum*modelarea];
	}
	memset(rotor,0,4*shadownum*modelarea*sizeof(double));
	memset(corners,-1,4*shadownum*modelarea*sizeof(int)); 
	init_rotor();
	
	
	for(int m(0); m < modelnum; ++m)
	{
		memcpy((smodels + m*12*modelarea*shadownum), (models+m*modelarea), modelarea * sizeof(unsigned short));
		cast_shadow( (smodels + m*12*modelarea*shadownum), (models+m*modelarea) );
		
		
		for(int s = 0; s < shadownum; ++s)
		{
			//unsigned short *smodelo( smodels + 12*modelarea* (m*shadownum + s-1) );
			unsigned short *smodel0( smodels + 12*modelarea* (m*shadownum + s) );
			
			for(int sm(1); sm < 12; ++sm) //the first shadow was already cast
			{
				unsigned short *smodel( smodel0 + sm * modelarea );
				const int rot( sm % 6 );
				const int mirror( sm / 6);
					
				//assert(MIRROT(m,rot,mirror) == mirrot_check++); //m for testing different hp
				//printf("m: %d, sm: %d \n",m,sm);
				fflush(stdout);
				for(int hp( 0 ); hp < hpMax; ++hp)
				{
					const int indhome( globalhexPixels[hp] );
					//assert(indhome >= 0);
					//assert(indhome < modelarea);
					const int indsym( inv_mirrotG[ MIRROT(hp,rot,mirror) ] );
					//assert(indsym >= 0);
					//assert(indsym < modelarea);
					
					smodel[indsym] = smodel0[indhome];
					/*
					if(s==0)
					{
						assert(smodel0[indhome]==smodel1[indhome]);
					}
					else if(s==1)
					{
						assert(smodel0[indhome]==smodelo[indhome]);
					}
					*/
					
				}
				//assert(smodel-smodels <= smdlsize);
				
				//printf("put model:%d shadow:%d k:%d\n",m,s,sm);
				
			}
		}
		
		//model += modelarea;
		//assert(model-models <= modelnum * modelarea);
	}
	
	if(diffmodel == nullptr )
	{
		//delete[] diffmodel;
		diffmodel = new short[12*shadownum*modelarea];
		
	}
	memset(diffmodel,0,12*shadownum*modelarea*sizeof(short));
	
	
	if(luse_beam)
	{
		
		if(shadownum > 1)
		{
			printf("abberated beams and continous rotations are not yet implemented\n");
			exit(0);
		}	
		
		const int beamarea( globalBeamArea );
		const int beamsize( globalBeamSize );
		
		setBeamCenter(beamsize);
		if(lbeam == nullptr)
		{
			//delete[] lbeam;
			lbeam = new unsigned short[beamarea];
		}
		memcpy(lbeam,beam,beamarea*sizeof(unsigned short));
		
		if(rmodels == nullptr)
		{
			//delete[] rmodels;
			rmodels = new unsigned short[smdlsize];
		}
		memcpy( rmodels, smodels, smdlsize*sizeof(unsigned short));
		
		if(corrmodels == nullptr)
		{
			//delete[] corrmodels;
			corrmodels = new int[smdlsize];
		
			memset(corrmodels,0,smdlsize*sizeof(int));
			build_corrmodels(0,modelnum-1);
		}
		
		apply_corrmodels(smodels,0,modelnum-1);	
	}
/*	
	//debug section for dphi==0.0
	for(int m = 0; m < modelnum; ++m)
	{
		for(int s = 1; s < shadownum; ++s)
		{
			const int home = 12*modelarea*s;
			assert( memcmp(smodels+m*12*modelarea*shadownum,smodels+m*12*modelarea*shadownum+home,
			12*modelarea*sizeof(unsigned short))==0 );	
		}
	}
*/	
	
}

void init_rotor()
{
	//const int modelarea( globalModelArea );
	const int shadownum( globalShadowNum );
	const int hpMax( globalhpMax );
	const int subframesize( globalModelSize );
	const int mx( globalMx );
	const int mz( globalMz );
	
	double dphi = -M_PIl/(3*shadownum); //*shadownum
	//dphi = 0.0;
	
	
	for(int s = 1; s < shadownum; ++s)
	{
		double phi = s*dphi;
		
		const double ca = cos(phi);
		const double nca = (1-ca);
		const double sa = sin(phi)/sqrt(3);
		const double diag = (1+2*ca)/3;
		const double ca3 = nca/3;
		const double p = ca3+sa;
		const double m = ca3-sa;
		
		for(int hp = 0; hp < hpMax; ++hp)
		{
			//const int hp2 ((3+hp+hpMax/2)%hpMax);
			const int ind( globalhexPixels[hp] );
			const int q = ind%subframesize;
			const int r = ind/subframesize;
			const int x = q - (r-(r&1))/2 - mx;
			const int z = r - mz;
			const int y = -x-z;
			
			
			double xphi = ((double)x)*diag + ((double)y)*m    + ((double)z)*p;
			//const double yphi = x*p    + y*diag + z*m;
			double zphi = ((double)x)*m    + ((double)y)*p    + ((double)z)*diag;
			//printf( "%lf = %d*%lf + %d*%lf + %d*%lf\n", zphi,x,m,y,p,z,diag ); 
			//const double yphi = -xphi-zphi;
			//10^-6 of a pixel should never matter
			const int xl = (int)floor(xphi+1E-6);
			const int xh = xl+1;
			const int zl = (int)floor(zphi+1E-6);
			const int zh = zl+1;
#if 0			
			{												
				const int chk_qr = ind;//inv_mirrotG[ MIRROT(hp,s,0) ];
				const int cq = chk_qr%subframesize;
				const int cr = chk_qr/subframesize;
				const int xm = cq - (cr-(cr&1))/2 - mx;
				const int zm = cr - mz;
				const int ym = -xm-zm;
				
				
				
				int cx(xl), cy(-xl-zl), cz(zl);
				
				hexagonal_trap(cx, cy, cz, subframesize/2);
	
	
				int x_data = cx + mx;
				int z_data = cz + mz;
				// cube -> odd-r
				int q_data = x_data + (z_data) / 2;
				int r_data = z_data;
				int indd = q_data + r_data * subframesize;
				
				static int counter(1);
				if(indd != chk_qr)
				{
					printf("OUCH, faulty rotation (q,r): (%d,%d) (x,y,z): (%d,%d,%d)\n",q,r,x,y,z);
					printf("(%d,%d) (%d,%d,%d) != (%d,%d) (%d,%d,%d)\n", 
							cq,cr,xm,ym,zm,			q_data,r_data, x_data-mx, z_data-mz, mx+mz-x_data-z_data) ;
					printf("raw (%lf,%lf,%lf)\n",xphi,-xphi-zphi,zphi);
					printf("hp: %d, shadow/rot: %d\n ", hp,s);
					printf("diag: %lf\n", diag);
					printf("m: %lf\tp: %lf\n", m, p);
					printf("nca: %lf\n",nca);
					printf("ca3: %lf\n",ca3);
					printf("sa: %lf\n", sa);
					if(--counter < 1)	exit(0);
				}
			
			}
#endif		
			const double fx = fmax(0.0,xphi-xl);
			const double fz = fmax(0.0,zphi-zl);
			
			const double ucx = (fx+fz)-1.0; //relative horizontal -1 at left and +1 at right
			const double ucy = (fx-fz); //relative vertical -1 at botom and +1 at top 
/*			
			assert(ucx <= 1.0);
			assert(ucx >= -1.0);
			
			assert(ucy <= 1.0);
			assert(ucy >= -1.0);
*/			
			
			double cshare = ( 1.0 + (fabs(ucy)/0.5) ) / 3.0;
			double left_share    = (1.0-fabs(ucx+1.0)*(11.0/12.0)) /( 3 * cshare );
			double right_share   = (1.0-fabs(ucx-1.0)*(11.0/12.0)) /( 3 * cshare );
			if(left_share < 0.0) left_share = 0.0;
			if(right_share < 0.0) right_share = 0.0;
			double central_share = 1.0 - left_share - right_share;
			double top_share = (1.0-(1-ucy)*0.5)*central_share;
			double bottom_share = central_share - top_share;
#if 0			
			if( fabs(1.0-left_share) > 1E-6 )
			{
				printf("OUCH, too small left share: %lf\t right share:%lf\n",left_share,right_share);
				printf("OUCH, too big top share: %lf\t bottom share:%lf\n",top_share,bottom_share);
				printf("fx: %lf \t fz: %lf\n", fx, fz);
				printf("ucx: %lf \t ucy: %lf\n", ucx, ucy);
				printf("start (q,r): (%d,%d) (x,y,z): (%d,%d,%d)\n",q,r,x,y,z);
				printf("raw (%lf,%lf,%lf)\n",xphi,-xphi-zphi,zphi);
				printf("hp: %d, shadow/rot: %d\n ", hp,s);
				printf("diag: %lf\n", diag);
				printf("m: %lf\tp: %lf\n", m, p);
				printf("nca: %lf\n",nca);
				printf("ca3: %lf\n",ca3);
				printf("sa: %lf\n", sa);
				
				exit(0);
			}
#endif			
/*			
			assert(central_share <= 1.0);
			assert(central_share >= 0.0);
			
			
			assert(left_share <= 1.0);
			assert(left_share >= 0.0);
			
			assert(right_share <= 1.0);
			assert(right_share >= 0.0);
			
			assert(top_share <= 1.0);
			assert(top_share >= 0.0);
			
			assert(bottom_share <= 1.0);
			assert(bottom_share >= 0.0);
			
			assert(fabs(left_share + right_share + bottom_share + top_share - 1) < 0.001);
*/			
			if(left_share > 0.0)	{addto_rotor(ind, xl,zl,s, 0, left_share);}
			if(right_share > 0.0)	{addto_rotor(ind, xh,zh,s, 1, right_share);}
			if(top_share > 0.0)		{addto_rotor(ind, xh,zl,s, 2, top_share);}
			if(bottom_share > 0.0)	{addto_rotor(ind, xl,zh,s, 3, bottom_share);}
		}
	
	
	}
	
	
	
	
}

void addto_rotor(int qr, int x, int z, int shadow, int offset, double share)
{
	
	int y = -x-z;
	const int modelarea( globalModelArea );
	const int subframesize( globalModelSize );
	const int mx( globalMx );
	const int mz( globalMz );
	
	hexagonal_trap(x, y, z, subframesize/2);
	
	int x_data = x + mx;
	int z_data = z + mz;
	// cube -> odd-r
	int q_data = x_data + (z_data) / 2;
	int r_data = z_data;
	int indd = q_data + r_data * subframesize;
	rotor[shadow*4*modelarea+4*qr+offset] += share;
	corners[shadow*4*modelarea+4*qr+offset] = indd;
}



void cast_shadow( unsigned short *smodel, unsigned short *model)
{
	const int modelarea( globalModelArea );
	const int shadownum( globalShadowNum );
	const int hpMax( globalhpMax );
	
	
	static double *tmp_model = new double[modelarea]; //memory leak
	
	
	for(int s = 1; s < shadownum; ++s)
	{
		const int home = 12*modelarea*s;
		std::fill_n(tmp_model,modelarea,0.0);
		for(int hp = 0; hp < hpMax; ++hp)
		{
			const int ind( globalhexPixels[hp] );
			double signal = model[ind];
			//double weight_chk = 0.0;
			//bool any_dest = false;
			//int dest_count = 0;
			for(int offset = 0; offset < 4; ++offset)
			{
				const int plc = s*4*modelarea+4*ind+offset;
				const int dest( corners[plc] );
				const double weight( rotor[plc] ); 
				if(dest != -1)
				{
					//++dest_count;
					//any_dest = true;
					tmp_model[dest]+=signal*weight;
					//weight_chk += weight;
					//assert(dest == ind);
				}
				
			}
/*			
			if(!any_dest)
			{
				printf("no destinations (%d) from (%d,%d)\n", dest_count, ind%64, ind/64);
			}
			
			
			if(fabs(weight_chk-1.0)>1E-6)
			{
				printf("weight_chk: %lf!=1.0 count: %d from (%d,%d)\n", weight_chk , dest_count, ind%64, ind/64);
				
				exit(0);
			}		 
*/	
		}
		for(int hp = 0; hp < hpMax; ++hp)
		{
			const int ind( globalhexPixels[hp] );
			const unsigned short nval = (unsigned short)(tmp_model[ind]+0.5);
			if(diffmodel != nullptr)
			{	diffmodel[home+ind]=(nval-smodel[home+ind]);}
			smodel[home+ind] = nval;
		}
	}
/*	//debug section for dphi==0.0 and only the first orientation
	for(int s = 1; s < shadownum; ++s)
	{
		const int home = 12*modelarea*s;
		assert( memcmp(smodel,smodel+home,modelarea*sizeof(unsigned short))==0 );
		if(diffmodel != nullptr)
		{	assert(	memcmp(diffmodel,diffmodel+home,modelarea*sizeof(short))==0 );}
	}
*/
	
}

void make_corrmodels(const int fm, const int lm)
{
	//printf("making corrmodels\n");
	//fflush(stdout);
	//const bool reinit_for_fun(::reinit_for_fun);
	
	const int modelarea( globalModelArea );
	const int shadownum( globalShadowNum );
	//const int chunksize( (lm-fm+1)*12*modelarea );
	const int chunkhome( fm*12*shadownum*modelarea);
	
	const int modelsize( globalModelSize );
	const int beamarea( globalBeamArea );
	const int beamsize( globalBeamSize );
	const int beamradius( globalBeamRadius );
	//const int beamnorm( globalBeamNorm );
	const int hpMax( globalhpMax );
	
	unsigned short *rmodel( rmodels + chunkhome);
	int  *corrmodel(  corrmodels + chunkhome);
	//int *ncrrmodel( ncrrmodels + chunkhome);

	//memcpy(ncrrmodel,corrmodel, chunksize*sizeof(int));


	
	for(int sm(12*shadownum*fm); sm < 12*shadownum*(lm+1); ++sm)
	{
		for(int hp(0); hp < hpMax; ++hp)
		{
			const int ind0( globalhexPixels[hp] );
			

			
			
			const int corrval( get_corrval(	ind0, modelsize, 
											beamsize, beamarea, beamradius,
											rmodel ) );
			//if(reinit_for_fun)
			//{	assert(corrmodel[ind0] == corrval);}
			corrmodel[ind0] = corrval;
			//delta function check
			/*
			if(corrval/beamnorm != rmodel[ind0])
			{
				
				//printf("delta beam failure in make_corrmodels()\n");
				const int crr(corrval/beamnorm);
				const int m( sm/12 );
				const int rmir( sm%12 );
				printf("rmodel:%u != crr:%d(raw: %d) in m:%d and sm:%d at hp:%d ind0:%d q,r: %d,%d \n",rmodel[ind0],crr,corrval,m,rmir,hp,ind0,q0,rz0);
				exit(0);
			}
			*/
			
		}
		rmodel += modelarea;
		corrmodel += modelarea;
	}
	//printf("make_corrmodels range: %d\n",(int)((corrmodel-corrmodels)/modelarea));
}

void build_corrmodels(const int fm, const int lm)
{
	//printf("building corrmodels\n");
	//fflush(stdout);
	
	
	const int modelarea( globalModelArea );
	//const int chunksize( (lm-fm+1)*12*modelarea );
	const int shadownum( globalShadowNum ); 
	const int chunkhome( fm*12*shadownum*modelarea);
	
	const int modelsize( globalModelSize );
	const int beamarea( globalBeamArea );
	const int beamsize( globalBeamSize );
	const int beamradius( globalBeamRadius );
	//const int beamnorm( globalBeamNorm );
	const int hpMax( globalhpMax );
	
	unsigned short *rmodel( rmodels + chunkhome);
	int  *corrmodel(  corrmodels + chunkhome);
	//int *ncrrmodel( ncrrmodels + chunkhome);

	//memcpy(ncrrmodel,corrmodel, chunksize*sizeof(int));


	
	for(int sm(12*shadownum*fm); sm < 12*shadownum*(lm+1); ++sm)
	{
		for(int hp(0); hp < hpMax; ++hp)
		{
			const int ind0( globalhexPixels[hp] );
			const int rmval( rmodel[ind0] );
			put_beam(	ind0, rmval, modelsize,
							beamsize, beamarea, beamradius,
							corrmodel);
			
		}
		rmodel += modelarea;
		corrmodel += modelarea;
	}
	//printf("make_corrmodels range: %d\n",(int)((corrmodel-corrmodels)/modelarea));
}




void apply_corrmodels(unsigned short *smodels, const int fm, const int lm) //fills diffmodel
{
	//printf("applying_corrmodels\n");
	//fflush(stdout);
	//const int smdlsize(12 * globalModelNum * globalModelArea);
	const int modelarea( globalModelArea );
	const int shadownum( globalShadowNum );
	const int chunksize( (lm-fm+1)*12*shadownum * modelarea );
	const int chunkhome( fm*12*shadownum*modelarea);
	
	const int beamnorm( globalBeamNorm );

	unsigned short *smodP( smodels+chunkhome );
	//unsigned short *nsmodP( nsmodels+chunkhome );
	//memcpy(nsmodP, smodP, chunksize*sizeof(unsigned short));
	
	short *diffP( diffmodel );
	int *corrP( corrmodels+chunkhome );
	
	lany_changes = false;
	for(int ind(0); ind < chunksize; ++ind)
	{
		const unsigned short val ((*(corrP++)+beamnorm/2)/beamnorm); //correctly rounded int division
		if(fm==lm) //-funswitch loops, only write diff when updating a single model
		{	
			const short diffval( val - *smodP );
			*(diffP++) = diffval;
			lany_changes = (lany_changes || (diffval != 0));
		}
		*(smodP++) = val;
		
		
		/*
		assert(diffval <= globalModelval);
		assert(diffval >= -globalModelval);
		*/
		//delta function check
		/*
		if(rmodels[ chunkhome + ind] != val)
		{
			const int m( (chunkhome+ind)/(12*globalModelArea));
			const int sm( (ind-m*12*globalModelArea)/globalModelArea) ;
			const int ind1(ind%globalModelArea);
			printf("rmodel:%u != corr:%u(%d/%d) in m:%d and sm:%d at ind1:%d ind:%d\n", 
			rmodels[ chunkhome + ind], val, *(corrP-1), beamnorm, m, sm, ind1, ind);
			//assert(diffval==0);
			exit(0);
		}
		*/
		
		
	}
	
}

short *update_smodels(	unsigned short *smodels, int active_m, int chlistlen, 
						int *chx , int *chy, int *newval, /*int *oldval,*/ bool &any_changes)
{
	
	const int shadownum( globalShadowNum );
	const int modelsize( globalModelSize );
	const int modelarea( globalModelArea );
	const int beamsize( globalBeamSize );
	const int beamarea( globalBeamArea );
	const int beamradius( globalBeamRadius );
	const int chunkhome( active_m*12*shadownum*modelarea );
	const bool use_beam( luse_beam );
	memset(diffmodel, 0, 12*shadownum*modelarea*sizeof(short));
	lany_changes = false;
	
	for(int chnum(0); chnum < chlistlen; ++chnum)
	{
		const int ind0(chx[chnum]+chy[chnum]*modelsize);
		const int hp0( pixposG[ind0] );
		const int nval( newval[chnum] );
		
		unsigned short *rmodel( (use_beam?rmodels:smodels) + chunkhome );
		//assert(memcmp(rmodel,rmodel+12*modelarea, 12*modelarea*sizeof(unsigned short))==0);
		short *diffm( diffmodel );
		int *crmodel( corrmodels + chunkhome);
		unsigned short *rmodel0( rmodel );
		const short diffval( nval - rmodel[ ind0 ] );
		if(diffval > nval)
		{	printf("OUCH overflow with old: %d new: %d diff: %d\n", rmodel[ind0], nval, diffval);}
		lany_changes = lany_changes || (diffval != 0);

		for(int sm(0); sm < 12; ++sm)
		{
			
			const int rot( sm % 6 );
			const int mirror( sm / 6);
			const int ind1( inv_mirrotG[ MIRROT(hp0,rot,mirror) ] );
			if(use_beam)
			{
				put_beam(	ind1, diffval, modelsize,
							beamsize, beamarea, beamradius,
							crmodel);
			}
			else
			{
				diffm[ind1] = diffval;
				diffm  += modelarea;
			}
			rmodel[ ind1 ] = nval;
			
			
			rmodel += modelarea;
			crmodel += modelarea;		
		}
			
		if( shadownum>1 )
		{	
			if(!use_beam)
			{
				cast_shadow(rmodel0,rmodel0); //casts all shadows at once
				roll_out_diffs(rmodel0);
			}
			else
			{
				printf("abberated beams and continous rotations are not yet implemented\n");
				exit(0);
			}	
		} 		
	}
	
	
	
	
	if(use_beam)
	{
		apply_corrmodels(smodels, active_m, active_m );	//fills diffmodels
	}
	any_changes = lany_changes;
	//assert(memcmp(diffmodel,diffmodel+12*modelarea, 12*modelarea*sizeof(short))==0);
/*	
	//debug section for dphi==0.0
	for(int s = 1; s < shadownum; ++s)
	{
		const int home = 12*modelarea*s;
		assert( memcmp(smodels+chunkhome,smodels+chunkhome+home,12*modelarea*sizeof(unsigned short))==0 );
		if(diffmodel != nullptr)
		{	assert(	memcmp(diffmodel,diffmodel+home,12*modelarea*sizeof(short))==0 );}
	}
*/
	
	return diffmodel;
}						


void roll_out_diffs(unsigned short *rmodel0)
{
	const int shadownum( globalShadowNum );
	const int modelarea( globalModelArea );
	const int hpMax( globalhpMax );
	for(int s = 1; s < shadownum; ++s)
	{
		const int chunkhome = s*12*modelarea;
		for(int hp( 0 ); hp < hpMax; ++hp)
		{
			const int indhome( chunkhome + globalhexPixels[hp]);
			const short diffval = diffmodel[indhome];
			const unsigned short val = rmodel0[indhome]; 
			if(diffval != 0)
			{
				for(int sm(1); sm < 12; ++sm) 
				{
					unsigned short *smodel( rmodel0 + chunkhome + sm * modelarea );
					short *diffshadow( diffmodel + chunkhome + sm * modelarea );
					const int rot( sm % 6 );
					const int mirror( sm / 6);
						
					const int indsym( inv_mirrotG[ MIRROT(hp,rot,mirror) ] );
			
					smodel[indsym] = val;
					diffshadow[indsym] = diffval;	
				}
			}
			
		}
	}
}




inline __attribute__((always_inline))
int get_corrval(	const int ind0, const int modelsize,
const int beamsize, const int beamarea, const int beamradius,
				unsigned short *const rmodel)
{
	int corrval(0);
	
	const int gmrH( globalMz );
	const int q0( ind0 % modelsize );
	//absolute cube coordinates of the smodel
	const int rz0( ind0 / modelsize );//r0 = z0 = rz0
	const int x0( q0 - ( rz0 - (rz0&1) ) /2);

	//relative cube coordinates to center of smodel
	const int dx0(x0 - globalMx);
	const int dz0(rz0 - globalMz);
	
	unsigned short const *beam( lbeam );
	for(int ind1(0); ind1 < beamarea; ++ind1)//scanning the beam
	{
		const int beamval( *(beam++) );
		if(beamval == 0)
		{	continue;}
		const int q1( ind1 % beamsize );
		const int rz1( ind1 / beamsize );//r0 = z0 = rz0
		const int x1( q1 - ( rz1 - (rz1&1) ) /2);
		//cube coordinates relative to center of the beam
		const int dx1(x1 - bcx);
		const int dz1(rz1 - bcz);
		const int dy1(-dx1 - dz1);
		if( (abs(dx1) + abs(dy1) + abs(dz1)) > 2*beamradius)
		{	continue;} //there cannot be contributions from outside the beam radius
		
		//relative cube to pixel covered by the beam
		int bx (dx0 + dx1);
		int bz (dz0 + dz1);
		int by (-bx-bz);

		hexagonal_trap(bx, by, bz, gmrH);

		//absolute cube coords of pixel covered by the beam
		bx += globalMx;
		bz += globalMz;

		const int qt = bx + (bz - (bz&1))/2;
		const int ind2 (qt + bz * modelsize);
		const int pixval( rmodel[ind2] );
		corrval += beamval*pixval;
	}
	return corrval;
}				


inline __attribute__((always_inline))
void put_beam(	const int ind0, const int delta, const int modelsize,
const int beamsize, const int beamarea, const int beamradius,
				int *const crmodel)
{
	const int gmrH( globalMz );
	const int q0( ind0 % modelsize );
	//absolute cube coordinates of the smodel
	const int rz0( ind0 / modelsize );//r0 = z0 = rz0
	const int x0( q0 - ( rz0 - (rz0&1) ) /2);

	//relative cube coordinates to center of smodel
	const int dx0(x0 - globalMx);
	const int dz0(rz0 - globalMz);
	
	unsigned short const *beam( lbeam ); 
	for(int ind1(0); ind1 < beamarea; ++ind1)//scanning the beam
	{
		const int beamval( *(beam++) );
		if(beamval == 0)
		{	continue;}
		const int q1( ind1 % beamsize );
		const int rz1( ind1 / beamsize );//r0 = z0 = rz0
		const int x1( q1 - ( rz1 - (rz1&1) ) /2);
		//inverted cube coordinates relative to center of the beam
		const int dx1(-x1 + bcx);//this line is different
		const int dz1(-rz1 + bcz);//this line is different
		const int dy1(-dx1 - dz1);
		if( (abs(dx1) + abs(dy1) + abs(dz1)) > 2*beamradius)
		{	continue;} //there cannot be contributions from outside the beam radius
		
		//relative cube to pixel covered by the beam
		int bx (dx0 + dx1);
		int bz (dz0 + dz1);
		int by (-bx-bz);

		hexagonal_trap(bx, by, bz, gmrH);

		//absolute cube coords of pixel covered by the beam
		bx += globalMx;
		bz += globalMz;
	
		//scale beam with diff and add/subtract it from the crmodel
		const int qt = bx + (bz - (bz&1))/2;
		const int ind2 (qt + bz * modelsize);
		crmodel[ind2] += (delta * beamval);
	}	
}				
				
				
				
