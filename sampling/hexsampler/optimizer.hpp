#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP

#include "globals.hpp"


//samples the frame top_points times and returns the mean
double mean_merit(double &std, const bool quick = false);

void optimize_merit(const bool lasting);




#endif
