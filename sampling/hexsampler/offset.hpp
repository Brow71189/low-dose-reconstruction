#ifndef OFFSET_HPP
#define OFFSET_HPP

#include "globals.hpp"



//looks up uc_sum with current uc_stepX and uc_stepY
int uc_sum_val(int x, int y);

//the same for symmetrized unitcell
int uc_pos_val(int x, int y);

int ruc_sum_val(int x, int y);

//double lattice_match(double sx, double sy);

double gauss_fit(void);
double rgauss_fit(void);

//writes symetrized unitcell with current uc_stepX/Y
//and returns a measure of symetry in uc_sum
double mirror(void);
double rmirror(void);

//returns the asymmetry under 90°rotation in uc_sym
//double asym(void);

//search for best fit between uc_avg and uc_gauss returns symetry parameter of uc_pos
double find_steps(void);

//applies next_offset_a1 and next_offset_a2
bool apply_a1_a2(void);

#endif
