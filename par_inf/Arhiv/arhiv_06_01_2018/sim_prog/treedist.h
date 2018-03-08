#ifndef TREE_DIST_H_INCLUDED
#define TREE_DIST_H_INCLUDED

#include "cons.h"
#include <string>

float phylip_dist(std::string & Probe_mut_1,
								  std::string & Probe_mut_2, int dist_type, int seed);

#endif // TREE_DIST_H
