#include "globals.hpp"
#include "basis.hpp"
#include "ellipse.hpp"
#include "offset.hpp"
#include "subframes.hpp"
#include "graphene.hpp"
#include "math.h"
#include <assert.h>

//#include <assert.h>

Command *globalCommand = 0;

std::map<int,double> sampled_merits;
std::vector<int> top5;
FILE *file_lock = nullptr;
int fd_lock = -1;
bool locking = false;
FILE *logfile = nullptr;
std::string logname, mem_check_location;
size_t top_points = 5;
const int max_rank( 15 ); //30bit bit for 32768x32768 search grid
const int mg_len( 1 + (2 << (max_rank-1)) );
int phint_min = 0;
int phint_max = mg_len - 1;
int ptint_min = 0;
int ptint_max = mg_len - 1;
int stability = 1; //0 .. quick and narrowsearch, 1 .. tested default, 2 .. safe, 3 .. extensive 
int init_stability = 1;
int peak_rank = 2;
int search_type = 0; // 0 .. hel and tilt, 1 .. excent and phi, 2 .. TODO 4d search


double min_light_dirt = 0.0;
double max_light_dirt = 0.0;
bool allow_light_dirt = false;
int light_dirt_val = 1;

int op_mode = 0;
clock_t timer;    
const double a0 = 0.142; //C-C bondlength in nm
             
int threadnum      = 0; //default
int frameNr        =-1; //our frame number, not relevant
int modelsize      =-1; //desired subframe/modelsize
int small_modelsize = -1;
int large_modelsize = -1;
int mask_scaling = 1;

int solidsize	   =-1;	 		
int impWidth       =-1; //frame width
int impHeight      =-1; //frame height
int hexImpWidth    =-1; //hexed frame width //TODO turn these to pointers
int hexImpHeight   =-1; //hexed frame height
int small_hexImpWidth    =-1; //hexed frame width
int small_hexImpHeight   =-1; //hexed frame height 
int large_hexImpWidth    =-1; //hexed frame width
int large_hexImpHeight   =-1; //hexed frame height   		
int bondlength_CC     = -1; //hexes per cc
int bondlength_Mini   = -1;
int *bondlength    = &bondlength_CC;
int offset_Q       = 1; //origin column for scanning subframes //TODO turn these to int*
int offset_R       = 1;	//origin row for scanning subframes 
int small_offset_Q = 1;
int small_offset_R = 1;
int large_offset_Q = 1;
int large_offset_R = 1;


int Num_sf         = 0; //Number of identified subframes

int *frame		    = nullptr; //array holding our raw array
unsigned char *mask = nullptr;	
int  **unitcells    = nullptr;
int **runitcells    = nullptr;
int *uc_sum		    = nullptr;
int *ruc_sum        = nullptr;
int *uc_rough	    = nullptr;
int *uc_mini        = nullptr;
double *uc_target   = nullptr;
int smooth_passes   = 2;	
int *uc_pos		    = nullptr;
int *uc_gauss	    = nullptr;
int *ruc_gauss      = nullptr;
double sigmaf       = (1/2.35); //FWHM ~ a0;
//double *uc_stat     = nullptr;
bool  *dirty        = nullptr;
bool *rdirty        = nullptr;
int defect_type     = 0;

int *hexagons       = nullptr;//TODO turn these to int**
int *hex_mask       = nullptr;
int *hexarea        = nullptr;
int *hex_trap       = nullptr;
int *graphene       = nullptr;

int *small_hexagons       = nullptr;
int *small_hex_mask       = nullptr;
int *small_hexarea        = nullptr;
int *small_trap           = nullptr;
int *small_graphene       = nullptr;

int *large_hexagons       = nullptr;
int *large_hex_mask       = nullptr;
int *large_hexarea        = nullptr;
int *large_trap           = nullptr;
int *large_graphene       = nullptr;



int **subframes     = nullptr;
int *sf_sum         = nullptr;

const int cube_neighbors_X[6] = {  1,  1,  0,  0, -1, -1};
const int cube_neighbors_Z[6] = {  0, -1,  1, -1,  1,  0};
const int cube_neighbors_Y[6] = { -1,  0, -1,  1,  0,  1};

int framelen	   = -1; //-1..nullpointer, otherwise length of frame;
int hexframelen    = -1; //-1..nullpointer, otherwise length of hexframe;
int small_hexframelen = -1;
int large_hexframelen = -1;

int uc_mini_len    = -1; //-1..nullpointer, otherwise length of uc_mini;

double data_avg    = 0.0; //integral
double data_sum    = 0.0; //integral
double data_min    = 0.0; //integral
double data_max    = 0.0; //integral
double raw_avg     = 0.0; //per unit cell
double rdata_avg   = 0.0; //integral
double rdata_sum   = 0.0; //integral
double rdata_min   = 0.0; //integral
double rdata_max   = 0.0; //integral
double rraw_avg    = 0.0; //per unit cell
double hex_avg     = 0.0; //per hex pixel
double hex_low     = 0.0;
double hex_high    = 0.0;
double hex_shape   = 0.0;
double uc_avg      = 0.0;
double contrast    = 0.0; //std/mean
double moment2     = 0.0;
double min_hex_contrast = 0.0125;

double fov		   =-1.0; //field_of_view aka impWidth in nm

double offset_X    = 0.0; //Origin of unitcell
double offset_Y    = 0.0; //Origin of unitcell
double offset_XS   = 0.0; //fixed Origin of unitcell for stable sampling
double offset_YS   = 0.0; //fixed Origin of unitcell for stable sampling
double next_offset_a1 = 0.0; //determined by find_steps() for the currently sampled unitcells
double next_offset_a2 = 0.0; //determined by find_steps() for the currently sampled unitcells
double next_offset_p1 = 0.0; //determined by get_hex_shape() for currently sampled uc_mini
double next_offset_p2 = 0.0; //determined by get_hex_shape() for currently sampled uc_mini
double offset_delta   = 0.0; //set inside apply_offset according to next_offset_a1 and next_offset_a2
double clean_X = 0.0;
double clean_Y = 0.0;
//int coh_len = -1;
//int small_coh_len = -1;
//int large_coh_len = -1;
//int max_coh_len = -1;

bool pending_offsets = false;
bool p1_and_p2_applied = true;

double box_hel_min = 0;
double box_hel_max = 100;
double box_hel_range = 100;
double box_tilt_min = -M_PI/30;
double box_tilt_max = M_PI/30;
double box_tilt_range = M_PI/15;

double peak_hel_min = box_hel_min;
double peak_hel_max = box_hel_max;
double peak_hel_range = box_hel_range;
double peak_tilt_min = box_tilt_min;
double peak_tilt_max = box_tilt_max;
double peak_tilt_range = box_tilt_range;


double elapsed = 0.0; //timing for sampling task
double merit = -1.0; //merit of optimized sampling -1 .. we did not even try
double hexagonal_merit = -1000.0;
double position_quality = -1.0;
double mirror_quality = -1.0;
//double coh_radius2 = -1.0;

//Vector basis of graphene unitcell
double hel		   =-1.0; //hex_edge_len = 2*impWidth/(3*r_mean*bondlength)
double tilt		   = 0.0; //radians for rotating the unitcell

//Ellipse
double ea = 1.0;     //communication
double eb = 1.0;     //communication
double excent = 1.0; //actual parameter
double phi = 0.0;    //actual parameter

int uc_res = 60;
//int uc_report_res = 60;
int sample_rate = SAMPLE_RATE;
int uc_area   = 3600;
int  uc_stepX = 0;
int ruc_stepX = 0;
int  uc_stepY = 0;
int ruc_stepY = 0;
int  uc_gridX = 1;
int ruc_gridX = 1;
int  uc_gridY = 1;
int ruc_gridY = 1;
int  uc_rangeX = 1;
int ruc_rangeX = 1;
int  uc_rangeY = 1;
int ruc_rangeY = 1;
int  uc_len    = -1;
int ruc_len    = -1;
int dirt_val   = 0xFFFF;//65535;

bool suppress_no_mini_cells = false;

bool skip_grid       = false;
bool retune_origin   = false;
bool has_console     = true;  //has to be set false by master
bool reporting       = true;  //extra output
bool report_uc       = false;
bool report_hexes    = false;
bool report_sf       = false;
bool volatile busy 	 = false;  //we are currently working on a frame
bool edge_incidence  = false;
bool use_mirror      = true; //definitely good in conjuction with hex_merit
bool use_position    = true;
bool floating_offset = true;
bool highlight_sf    = false;
bool soft_sf         = false;
int rotate_tilt      = 0;
int ruling_merit     = 1; // 0 .. merit, 1 hexagonal_merit
int search_pattern   = 1; // 0 .. tilt and hel, 1 .. phi and excent
int report_num       = 0;

bool small_hexes = true;

void set_large_hexes()	//actually smaller pixels
{	
	/*******SWITCHING******/
	assert(small_hexes);
	small_hexes = false;
	bondlength = &bondlength_Mini;
	hexImpWidth = large_hexImpWidth;
	hexImpHeight = large_hexImpHeight;
	modelsize = large_modelsize; 
	//coh_len = small_coh_len * (int)(bondlength_Mini/bondlength_CC);
	small_offset_Q = offset_Q;
	small_offset_R = offset_R;
	offset_Q = large_offset_Q;
	offset_R = large_offset_R;
	modelsize = large_modelsize;
	hex_mask = large_hex_mask;
	hexagons = large_hexagons;
	hexarea = large_hexarea;
	hex_trap = large_trap;
	graphene = large_graphene;
	hexframelen = large_hexframelen;
	/***********************/
}

void set_small_hexes() //actually bigger pixels
{
	/*******SWITCHING_BACK*/
	assert(!small_hexes);
	small_hexes = true;
	bondlength = &bondlength_CC;
	hexImpWidth = small_hexImpWidth;
	hexImpHeight = small_hexImpHeight;
	//coh_len = small_coh_len;
	large_offset_Q = offset_Q;
	large_offset_R = offset_R;
	offset_Q = small_offset_Q;
	offset_R = small_offset_R;
	modelsize = small_modelsize;
	hex_mask = small_hex_mask;
	hexagons = small_hexagons;
	hexarea = small_hexarea;
	hex_trap = small_trap;
	graphene = small_graphene;
	hexframelen = small_hexframelen;
	/***********************/
}

void init_small_hexes()
{
	const double bl_scaling(((double)bondlength_CC)/((double)bondlength_Mini));
	
	small_hexImpWidth = hexImpWidth;
	small_hexImpHeight = hexImpHeight;
	small_hexframelen = small_hexImpWidth * small_hexImpHeight;
	small_modelsize = modelsize;
	large_hexImpWidth  = (int) (((double)hexImpWidth) / bl_scaling);
	large_hexImpHeight = (int)(((double)hexImpHeight) / bl_scaling);
	if(large_hexImpWidth % 2 == 1) { ++large_hexImpWidth; }
	if(large_hexImpHeight % 2 == 1) { ++large_hexImpHeight; }
	large_hexframelen = large_hexImpWidth * large_hexImpHeight;
	large_modelsize = 2 * bondlength_Mini;
	delete[] small_hexagons;
	small_hexagons = new int[small_hexframelen];
	delete[] small_hexarea;
	small_hexarea = new int[small_hexframelen];
	delete[] small_hex_mask;
	small_hex_mask = new int[small_hexframelen];
	//delete[] small_graphene;
	//small_graphene = new int[small_hexframelen];
	
	delete[] large_hexagons;
	large_hexagons = new int[large_hexframelen];
	delete[] large_hexarea;
	large_hexarea = new int[large_hexframelen];
	delete[] large_hex_mask;
	large_hex_mask = new int[large_hexframelen];
	delete[] large_trap;
	large_trap = new int[large_hexframelen];
	delete[] large_graphene;
	large_graphene = new int[large_hexframelen];
	
	const int rH( bondlength_Mini );
	const int mcx( rH - (rH - (rH&1) )/2 );
	const int mcz( rH);
	int *ltrp(large_trap);
	for(int ind(0); ind < large_hexframelen; ++ind)
	{
		const int q( ind % large_hexImpWidth - large_hexImpWidth/2);
		const int r( ind / large_hexImpWidth - large_hexImpHeight/2);
		int x( q - (r - (r&1) )/2); // + dx
		int z( r ); //+ dz
		int y( -x -z);
		hexagonal_trap(x,y,z,rH);
		const int px0( mcx + x );
		const int pz0( mcz + z );
		const int qh0( px0 + (pz0 - (pz0&1) )/2 );
		const int rh0( pz0 );
		const int pos0( qh0 + rh0 * 2 * rH );
		*(ltrp++) = pos0;
	}
	init_graphene();//actually large graphene
	
	small_hexes = false;
	set_small_hexes();
}



Config::Config(void)
{	
	_uc_sum 	= nullptr;
	_ruc_sum	= nullptr;
	_uc_pos 	= nullptr;
	_uc_gauss 	= nullptr;
	_uc_mini 	= nullptr;
	empty 		= true;
	//recent_save = false;
}


Config::~Config(void)
{	
	clear();
}

bool Config::exists(void)
{
	return(!empty);
}

bool Config::any_changes(void) 
{ 	
	//merits are noisy they do not count here
	return( !( (_hel == hel) && (_tilt == tilt) && (_phi == phi) &&
			   (_excent == excent) && (_offset_X == offset_X) && (_offset_Y == offset_Y) &&
			   (_offset_Q == offset_Q) && (_offset_R == offset_R)	
		     )	  );
}

void Config::save(void)
{
	//if(reporting) printf("saving config\n");
	//apply_offsets(); //anyways always done in find_steps
	fresh_config = true; //true until send to master
	empty       = false;
	_hel 		= hel; //hex_edge_len = 2*impWidth/(3*r_mean*bondlength)
	_tilt 		= tilt; //radians for rotating the unitcell
	_phi        = phi;
	_excent     = excent;
	_hex_shape  = hex_shape;
	_hex_avg    = hex_avg;
	_uc_avg     = uc_avg;
	_hex_low    = hex_low;
	_hex_high   = hex_high;
	_contrast   = contrast;
	_uc_res     = uc_res;
	_uc_stepX   = uc_stepX;
	_ruc_stepX  = ruc_stepX;
	_uc_stepY   = uc_stepY;
	_ruc_stepY  = ruc_stepY;
	_offset_X 	= offset_X; //Origin of unitcell
	_offset_Y 	= offset_Y; //Origin of unitcell
	_offset_XS 	= offset_XS; //conservative Origin of unitcell
	_offset_YS 	= offset_YS; //conservative Origin of unitcell
	_merit 		= merit;
	_next_offset_a1   = next_offset_a1;
	_next_offset_a2   = next_offset_a2;
	_next_offset_p1   = next_offset_p1;
	_next_offset_p2   = next_offset_p2;
	_hexagonal_merit  = hexagonal_merit;
	_position_quality = position_quality;
	_mirror_quality   = mirror_quality;
	_moment2		  = moment2;	
	
	_offset_Q 	= offset_Q; //origin column for scanning subframes
	_offset_R 	= offset_R;	//origin row for scanning subframes 
	_Num_sf 	= Num_sf; //Number of identified subframes	
	if(uc_sum != nullptr) //check for instantiaion of global config
	{
		delete[] _uc_sum;
		_uc_sum = new int[uc_area];
		memcpy(_uc_sum,uc_sum,uc_area*sizeof(int));
	}
	if(ruc_sum != nullptr) //check for instantiaion of global config
	{
		delete[] _ruc_sum;
		_ruc_sum = new int[uc_area];
		memcpy(_ruc_sum,ruc_sum,uc_area*sizeof(int));
	}
	if(uc_pos != nullptr)
	{
		delete[] _uc_pos;
		_uc_pos = new int[uc_area];
		memcpy(_uc_pos,uc_pos,uc_area*sizeof(int));
	}
	if(uc_gauss != nullptr)
	{
		delete[] _uc_gauss;
		_uc_gauss = new int[uc_area];
		memcpy(_uc_gauss,uc_gauss,uc_area*sizeof(int));
	}
	if(uc_mini != nullptr)
	{
		delete[] _uc_mini;
		_uc_mini = new int[uc_mini_len];
		memcpy(_uc_mini,uc_mini,uc_mini_len*sizeof(int));
	}
}



void Config::load(void)
{
	//set_basis(_hel, _tilt);
	hel         = _hel;
	tilt        = _tilt;
	set_ellipse(_excent,_phi); //relies on current hel
	uc_res      = _uc_res;
	uc_stepX    = _uc_stepX;
	ruc_stepX   = _ruc_stepX;
	uc_stepY    = _uc_stepY;
	ruc_stepY   = _ruc_stepY;
	hex_shape   = _hex_shape;
	moment2     = _moment2;
	hex_avg     = _hex_avg;
	hex_low     = _hex_low;
	hex_high    = _hex_high;
	uc_avg      = _uc_avg;
	contrast    = _contrast;
	uc_area     = uc_res*uc_res;
	offset_X 	= _offset_X; //Origin of unitcell
	offset_Y 	= _offset_Y; //Origin of unitcell
	offset_XS 	= _offset_XS; //conservative Origin of unitcell
	offset_YS 	= _offset_YS; //conservative Origin of unitcell
	merit 		= _merit;
	next_offset_a1   = _next_offset_a1;
	next_offset_a2   = _next_offset_a2;
	next_offset_p1   = _next_offset_p1;
	next_offset_p2   = _next_offset_p2;
	hexagonal_merit  = _hexagonal_merit;
	position_quality = _position_quality;
	mirror_quality   = _mirror_quality;
	offset_Q 	= _offset_Q; //origin column for scanning subframes
	offset_R 	= _offset_R;	//origin row for scanning subframes 
	Num_sf 		= _Num_sf; //Number of identified subframes	
	if( (_uc_sum != nullptr) && (uc_sum != nullptr) )
	{
		delete[] uc_sum;
		uc_sum = new int[uc_area];
		memcpy(uc_sum, _uc_sum, uc_area*sizeof(int));
	}
	if( (_ruc_sum != nullptr) && (ruc_sum != nullptr) )
	{
		delete[] ruc_sum;
		ruc_sum = new int[uc_area];
		memcpy(ruc_sum, _ruc_sum, uc_area*sizeof(int));
	}
	if( (_uc_pos != nullptr) && (uc_pos != nullptr) )
	{
		delete[] uc_pos;
		uc_pos = new int[uc_area];
		memcpy(uc_pos, _uc_pos, uc_area*sizeof(int));
	}
	if( (_uc_gauss != nullptr) && (uc_gauss != nullptr) )
	{
		delete[] uc_gauss;
		uc_gauss = new int[uc_area];
		memcpy(uc_gauss, _uc_gauss, uc_area*sizeof(int));
	}
	if( (_uc_mini != nullptr) && (uc_mini != nullptr) )
	{
		delete[] uc_mini;
		uc_mini = new int[uc_mini_len];
		memcpy(uc_mini, _uc_mini, uc_mini_len*sizeof(int));
	}
}

void Config::clear(void)
{
	delete[] _uc_sum;
	delete[] _ruc_sum;
	delete[] _uc_pos;
	delete[] _uc_gauss;
	delete[] _uc_mini;
	_uc_sum = nullptr;
	_ruc_sum = nullptr;
	_uc_pos = nullptr;
	_uc_gauss = nullptr;
	_uc_mini = nullptr;
	empty = true;
	//check_memory();	
}

Config config = Config();

void free_globals(void)
{
#ifndef NOCATCH	
	try
	{
#endif		
		fprintf(logfile,"clearing globals ... "); fflush(logfile);
		delete[] frame;
		frame = nullptr;
		delete[] mask;
		mask = nullptr;
		//TODO no longer delete them once they are aliases
		/*
		delete[] hexagons;
		hexagons = nullptr;
		delete[] hex_mask;
		hex_mask = nullptr;
		delete[] hexarea;
		hexarea = nullptr;
		*/
		delete[] small_hexagons;
		small_hexagons = nullptr;
		delete[] small_hex_mask;
		small_hex_mask = nullptr;
		delete[] small_hexarea;
		small_hexarea = nullptr;
		
		delete[] large_hexagons;
		large_hexagons = nullptr;
		delete[] large_hex_mask;
		large_hex_mask = nullptr;
		delete[] large_hexarea;
		large_hexarea = nullptr;
		delete[] large_trap;
		large_trap = nullptr;
		delete[] large_graphene;
		large_graphene = nullptr;
		
		delete[] uc_sum;
		uc_sum = nullptr;
		delete[] ruc_sum;
		ruc_sum = nullptr;
		delete[] uc_pos;
		uc_pos = nullptr;
		//delete[] uc_stat;
		//uc_stat = nullptr;
		delete[] uc_rough;
		uc_rough = nullptr;
		delete[] uc_gauss;
		uc_gauss = nullptr;
		delete[] sf_sum;
		sf_sum = nullptr;
		delete[] uc_mini;
		uc_mini = nullptr;
		delete[] uc_target;
		uc_target = nullptr;
		
		sampled_merits.clear();
		top5.clear();
		if(subframes != nullptr)
		{
			for(int sf(0); sf < Num_sf; ++sf)
			{	delete[] subframes[sf];}
			delete[] subframes;
			subframes = nullptr;
		}
		
		if(unitcells != nullptr)
		{
			int **uc(unitcells);
			for(int raw_ind = 0; raw_ind < uc_len ; ++raw_ind)
			{	delete[] *(uc++);}
			delete[] unitcells;
			delete[] dirty;
			unitcells = nullptr;
			dirty = nullptr;
		}
		if(runitcells != nullptr)
		{
			int **ruc(runitcells);
			for(int raw_ind = 0; raw_ind < ruc_len ; ++raw_ind)
			{	delete[] *(ruc++);}
			delete[] runitcells;
			delete[] rdirty;
			runitcells = nullptr;
			rdirty = nullptr;
		}
		min_light_dirt = 0.0;
		max_light_dirt = 0.0;
		light_dirt_val = 1;
		allow_light_dirt = false;
		bondlength = &bondlength_CC;
		floating_offset = true;
		highlight_sf = false;
		soft_sf = false;
		retune_origin = false;
		skip_grid = false;
		defect_type = 0;
		ruling_merit = 1;
		search_pattern = 1;
		search_type = 0;
		framelen = -1;
		hexframelen = -1;
		small_hexframelen = -1;
		large_hexframelen = -1;
		mask_scaling = -1;
		Num_sf = 0;
		report_num = 0;
		 uc_len = -1;
		ruc_len = -1;
		uc_mini_len = -1;
		//coh_len = -1;
		//max_coh_len = -1;
		//coh_radius2 = -1.0;
		config.clear();
		//check_memory();
		fprintf(logfile,"done\n"); fflush(logfile);
#ifndef NOCATCH	
	}catch (std::exception &e)
	{	
		if(reporting && !has_console)
		{	printf("Message ");}
		printf("ERROR exception in free_globals()\n");
		fflush(stdout);
		throw(e);
	}
#endif	
}
