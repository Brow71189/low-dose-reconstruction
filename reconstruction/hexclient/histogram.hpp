#ifndef HISTOGRAM
#define HISTOGRAM

void create_Histo(
				ccflt **cctableB,
#ifdef USE_MPFR
                ccint **cctableS,
#endif
				unsigned short *qImg,
				unsigned short *Models
			);
			
void report_Histo(bool reporting);


void clear_Histo();







#endif

