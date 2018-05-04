#ifndef ABBERATION_HPP
#define ABBERATION_HPP


void init_smodels( unsigned short *smodels, unsigned short *models, unsigned short *beam, bool use_beam = true);

short *update_smodels(	unsigned short *smodels, int active_m, int chlistlen,
						int *chx , int *chy, int *newval, /*int *oldval,*/ bool &any_changes);
/*
void revert_smodels(unsigned short *smodels, int active_m);
*/
void close_smodels();





#endif
