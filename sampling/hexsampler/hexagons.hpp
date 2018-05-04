#ifndef HEXAGONS_HPP
#define HEXAGONS_HPP




double sample_hexagons(double &tmp_hex_shape, bool skip_mini = false);

bool apply_p1_p2(void);
void cancel_p1_p2(void);

//applies 2D smoothing to squared shaped double arrays with uc_mini_len pixels
double* mini_smooth( double * const rough);

//returns std assumes that the mean is zero
double get_hex_std(double * const uc_mini_norm);

double get_hex_stats(const double * const uc_mini_norm, double& std);

//normalizes the hex image so that avg=0 and std=1.0
//returns the new unsigned volume
double normalize_hex(double * const uc_mini_norm, const double avg, const double std);

#endif
