#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

#include "globals.hpp"


/* These functions read data from FILE* in the format that is expected
 * from the Hex_Sampler Java ImageJ Plugin or "< job.bin"
 */


int read_input(void); //0 .. normal, 1 .. reinit, -1 .. quit now

void report(void); //update the results and report_results

void report_results(void); //print all findings
/*
void sendshort(unsigned short val);
*/


#endif
