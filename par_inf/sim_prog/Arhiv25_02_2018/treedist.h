////////////////////////////////////////////////////////////////////////////////
//We are using treedist.cpp from PHYLogeny Inference Package (PHYLIP) to adopt
//two standard tree metrics: 1) Robinson–Foulds Distance (SYMMETRIC)
//                           2) Branch Score Distance (BSD)
//See http://evolution.genetics.washington.edu/phylip/doc/treedist.html
//for details
////////////////////////////////////////////////////////////////////////////////

//==============================================================================
// include guard
#ifndef TREE_DIST_H_INCLUDED
#define TREE_DIST_H_INCLUDED

//==============================================================================
// include dependencies
#include <string> // need for parameters of the function phylip_dist()

//==============================================================================
//function definition

//==================================
// phylip_dist() returns distances between two trees.
// @probe_mut_1 and @probe_mut_2  are  Newick notations of two compared trees.
// @dist_type is a setting for the type of tree distance:
// dist_type = 0 - Robinson – Foulds Distance (SYMMETRIC);
// dist_type = 1 - Branch Score Distance (BSD).
// @seed is for temporal file name:
// we need to simulate unique file to apply PHYLIP algorithm.
float phylip_dist(std::string & probe_mut_1,
								  std::string & probe_mut_2,
									int dist_type,
									int seed);


#endif // TREE_DIST_H
