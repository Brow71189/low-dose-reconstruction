#ifndef GLOBALS_HPP
#define GLOBALS_HPP
#include "mcheck.h"
#include "stdio.h"
#include "string.h"
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <time.h>
#include "command.hpp"

#define SAMPLE_RATE 256 // possibly any difference to 256!?
#define MIN_ITER 1
#define MAX_ITER 1
#define MAX_ITER_LIMIT 2
//does not cost much and will always converge at a convex shape of valid hexagonal pixels
//#define LAST_PASS 4
//#define MISS_HEX_W 0

#define THIRD (1.0/3.0)
#define TWOTHIRD (2.0/3.0)

#if 0
#define check_memory() \
	{\
		std::stringstream mcl_stream;\
		mcl_stream << " file: " << __FILE__ << " line: " << __LINE__;\
		mem_check_location = mcl_stream.str();\
		mcheck_check_all();\
	}
	
#define ACTIVE_CHECKS 
	
#else

#define check_memory()

#endif

extern Command *globalCommand;
void free_globals(void);
extern FILE *file_lock;
extern int fd_lock;
extern bool locking;
extern FILE *logfile;
extern std::string logname;
extern std::string mem_check_location;
extern std::map<int,double> sampled_merits;
extern std::vector<int> top5;
extern size_t top_points;
//extern double *const merit_grid;
extern const int max_rank;
extern const int mg_len;
extern int stability;
extern int init_stability;
extern int peak_rank; 
extern const double a0; //C-C bondlength in nm

extern bool suppress_no_mini_cells;
 
extern clock_t timer;    
         
extern int threadnum; //default
extern int frameNr; //our frame number, not relevant
extern int modelsize; //desired subframe/modelsize
extern int solidsize; 		
extern int impWidth; //frame width
extern int impHeight; //frame height
extern int hexImpWidth; // hexed frame width
extern int hexImpHeight; //hexed frame height  	 		
extern int small_hexImpWidth; // hexed frame width
extern int small_hexImpHeight; //hexed frame height
extern int large_hexImpWidth; // hexed frame width
extern int large_hexImpHeight; //hexed frame height
extern int mask_scaling;

extern int *bondlength; //hexes per cc
extern int bondlength_CC;
extern int bondlength_Mini;
extern int offset_Q; //origin column for scanning subframes
extern int offset_R; //origin row for scanning subframes
extern int Num_sf;//Number of identified subframes	 

extern int *frame; //array holding our raw array
extern unsigned char *mask; //array holding our mask
extern int  **unitcells;
extern int **runitcells;
extern int *hexagons;
extern int *hex_mask;
extern int *hexarea;
extern int *hex_trap;
extern int *graphene;

extern int *small_hexagons;
extern int *small_hex_mask;
extern int *small_hexarea;
extern int *small_trap;
extern int *small_graphene;

extern int *large_hexagons;
extern int *large_hex_mask;
extern int *large_hexarea;
extern int *large_trap;
extern int *large_graphene;

extern int *uc_sum;
extern int *ruc_sum;
extern int *uc_rough;
extern int *uc_mini;
extern double *uc_target; 
extern int smooth_passes;
extern int *uc_pos;
extern int *uc_gauss;
extern int *ruc_gauss;
extern double sigmaf; //inverse sigma of atomic resolution in units of a0;
//extern double *uc_stat;
extern bool  *dirty;
extern bool *rdirty;
extern int defect_type;

extern int **subframes;
extern int *sf_sum;

extern const int cube_neighbors_X[6];
extern const int cube_neighbors_Z[6];
extern const int cube_neighbors_Y[6];


//extern bool uc_valid;
extern int framelen; //-1..nullpointer, otherwise length of array;
extern int hexframelen;
extern int small_hexframelen;
extern int large_hexframelen;
extern int uc_mini_len;
extern int sample_rate;
extern int op_mode;
extern int report_num;


extern double min_light_dirt;
extern double max_light_dirt;
extern bool allow_light_dirt;
extern int light_dirt_val;
extern double fov; //field_of_view aka impWidth in nm

extern double offset_X; //Origin of unitcell
extern double offset_Y; //Origin of unitcell
extern double offset_XS; //fixed offset for stable sampling
extern double offset_YS; //fixed offset for stable sampling
extern double clean_X; //center of clean area
extern double clean_Y; //center of clean area
extern double next_offset_a1; //tentative guess
extern double next_offset_a2; //tentative guess
extern double next_offset_p1; //tentative guess
extern double next_offset_p2; //tentative guess
extern double offset_delta; //cartesian length of tetative guesss
extern bool pending_offsets;  //if there are tentative guesses
extern bool p1_and_p2_applied; 

extern int search_type; // 0 .. tilt and hel, 1 .. phi and excent
extern double box_hel_min; //actually hel or excent
extern double box_hel_max; //actually hel or excent
extern double box_hel_range; //actually hel or excent
extern double box_tilt_min;  //actually tilt or phi
extern double box_tilt_max;   //actually tilt or phi
extern double box_tilt_range;  //actually tilt or phi

extern double peak_hel_min; //actually hel or excent
extern double peak_hel_max; //actually hel or excent
extern double peak_hel_range; //actually hel or excent
extern int 	  phint_min;
extern int    phint_max;	
extern double peak_tilt_min; //actually tilt or phi
extern double peak_tilt_max; //actually tilt or phi
extern int    ptint_min;
extern int    ptint_max;
extern double peak_tilt_range; //actually tilt or phi


extern double elapsed; //timing for sampling task
extern double merit; //merit of optimized sampling .. -1 we did not even try
extern double hexagonal_merit;
extern double contrast;
extern double min_hex_contrast;
extern double position_quality;
extern double mirror_quality;
extern double rdata_avg, data_avg;
extern double rdata_sum, data_sum;
extern double rdata_max, data_max;
extern double rdata_min, data_min;
extern double rraw_avg, raw_avg;
extern double uc_avg;
extern double hex_avg;
extern double hex_low;
extern double hex_high;
extern double hex_shape;
extern double moment2;
//extern int coh_len; // coherent length in pixels
//extern int small_coh_len; //TODO make only this public
//extern int large_coh_len;

//extern int max_coh_len;
//extern double coh_radius2; //measure of coh_len^2 in raw_data pixels

//Vector basis of graphene unitcell
extern double hel; //hex_edge_len = 2*impWidth/(3*r_mean*bondlength)
extern double tilt; //radians for rotating the unitcell 

//Ellipse
extern double ea;      //communication with master
extern double eb;      //communication with master
extern double excent;  //actuall parameter
extern double phi;     //actuall parameter

extern int uc_res; //resolution of unitcell
//extern int uc_report_res;
extern int uc_area;
extern int ruc_stepX, uc_stepX;
extern int ruc_stepY, uc_stepY;
extern int ruc_gridX, uc_gridX;
extern int ruc_gridY, uc_gridY;
extern int ruc_rangeX, uc_rangeX;
extern int ruc_rangeY, uc_rangeY;
extern int ruc_len, uc_len;

extern bool skip_grid;
extern int dirt_val;
extern bool retune_origin;
extern bool has_console;
extern bool reporting;  //extra output
extern bool report_uc; //send frames
extern bool report_hexes;
extern bool report_sf;
extern volatile bool busy;
extern bool edge_incidence;
extern bool use_mirror;
extern bool use_position;
extern bool floating_offset;
extern bool highlight_sf;
extern bool soft_sf;
extern int rotate_tilt;

extern int ruling_merit; //0 .. regular 1.. hex_merit
extern int search_pattern; // 0 .. grid/peak 1 .. linear trials

void set_large_hexes();  //for hex_optimizer get_hex_merit 
void set_small_hexes();  //for hex_optimizer get_hex_merit
void init_small_hexes(); //for communicator  read_input

class Config
{
public:
	Config();
	~Config();
	void save(void);
	void load(void);
	void clear(void);
	bool exists(void);
	bool any_changes(void);
	bool volatile fresh_config;
private:
	bool empty;
	double _hel; //hex_edge_len = 2*impWidth/(3*r_mean*bondlength)
	double _tilt; //radians for rotating the unitcell
	double _phi;  //radians for rotating long axis of ellipse in FFT
	double _excent; //excentrity >= 1 of the ellipse
	double _offset_X; //Origin of unitcell
	double _offset_Y; //Origin of unitcell
	double _offset_XS; //Origin of unitcell
	double _offset_YS; //Origin of unitcell
	double _next_offset_a1;
	double _next_offset_a2;
	double _next_offset_p1;
	double _next_offset_p2;
	
	
	double _hex_shape;
	double _hex_avg;
	double _hex_low;
	double _hex_high;
	double _uc_avg;
	double _contrast;
	double _merit;
	double _hexagonal_merit;
	double _position_quality;
	double _mirror_quality;
	double _moment2;
	int _uc_res;
	int _uc_stepX;
	int _ruc_stepX;
	int _uc_stepY;
	int _ruc_stepY;
	int _offset_Q; //origin column for scanning subframes
	int _offset_R;	//origin row for scanning subframes 
	int _Num_sf; //Number of identified subframes
	int *_uc_sum;
	int *_ruc_sum;
	int *_uc_pos;
	int *_uc_gauss;
	int *_uc_mini;
};

extern Config config,reset;

#endif
