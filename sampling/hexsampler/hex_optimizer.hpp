#ifndef HEX_OPTIMIZER_HPP
#define HEX_OPTIMIZER_HPP

#include "globals.hpp"


//modfies offset_X and offset_Y to optimize hex_merit from sample_hexagons
//sampling offset are NOT UPDATED automatically, caller ought to verify them first
//no config is saved caller has to verify the changes before commiting
double optimize_hexagons(int mode, int rank); //-1 .. none, 0 .. offsets, 1 .. tilt & hel 2 .. phi & excent 
double get_hex_merit(double offx, double offy, double t, double h, double p, double e, double &std, bool quick = false);
#endif
