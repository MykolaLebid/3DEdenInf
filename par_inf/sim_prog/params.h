// Copyright 2017 Lebid Mykola all rights reserved
#ifndef PARAMS_H_INCLUDED
#define PARAMS_H_INCLUDED

#include <cmath>
//---- exactly one method should be defined------
//#define GILLESPIE
//#define FASTER_KMC
#define NORMAL  // standard simulation method described in the paper (non-KMC)

//-----------------------------------------------

//#define PUSHING // if defined, cells can push away other cells as they grow
//#define CONST_BIRTH_RATE
// if defined, birth rate does not depend on the
//no. of empty neighbors as long as there is at least one

// --exactly one of these should be defined as the replication neighborhood--
#define VON_NEUMANN_NEIGHBOURHOOD // 6 neighbors
//#define VON_NEUMANN_NEIGHBOURHOOD_QUADRATIC
// 6 neighbors but non linear growth rate (12n-n^2)/36 where
// n = no. of empty neighbors
//#define MOORE_NEIGHBOURHOOD // 26 neighbors
//--------------------------------------------------------
// if defined, simulate treatment after reaching given size
//#define MAKE_TREATMENT_N
// if defined, simulate treatment after reaching given time
//#define MAKE_TREATMENT_T

const float gama_res = 0;
// these are rates per daughter cell. Rates per diploid exome
// will be 2x higher (these values are given in the paper)
extern float gama; //5e-8 ;
//#define MIGRATION_MATRIX

// used when MIGRATION_MATRIX is not defined
//const float migr=0.00000 ;
// used only when MIGRATION_MATRIX is defined
//const float migr[2][2]={{0,0} // before treatment: WT/resistant
//                        ,{0,1e-5}}; // after treatment: WT/resistant
// if defined then cells die on surface only upon treatment,
// if not then cells die also in the volume
//#define DEATH_ON_SURFACE ;
//#define CORE_IS_DEAD  // when set, core cells are set to dead
//#define SHOW_ONLY_DRIVERS
// if defined, no mechanics is simulated
// (this speeds up everything but looses info on spatial distribution)
//#define NO_MECHANICS

const float timescale=1./log(2.) ; // calculates the timescale factor from the division time [days] (first number, here 1.)
//extern float death0, growth0;   // before treatment

// if death on surface:
//float death1=/*0.1*/0.99, growth1=0.0 ;   // after treatment
// if death in volume
//const float death1=1., growth1=0.5 ; // after treatment

extern float driver_adv;


extern float driver_prob; // driver probability per haploid genome (should be 2e-5)
//const float driver_balance=1 ; // 1==drivers decrease death, 0==drivers increase growth, intermediate values are also possible
//const float driver_migr_adv=0,  max_migr=1e-5 ; // maximum migration prob. is taken into account only when driver_migr_adv>0
//const int driver_mode = 0 ; // 0== drivers affect bd only, 1==drivers affect simultaneously bd and migr, 2==drivers affect bd xor migr with equal prob.

const float cutoff=0.1 ;
//#ifdef __MAIN
//int max_size=int(1e6) ; // this is ignored when MAKE_TREATMENT_T is defined
//float time_to_treat=10*(365./12) ; // this is ignored when MAKE_TREATMENT_N is defined
//#endif

const int _resol=1 ; // spatial resolution of sampling [cells]
const int _bins=10000 ; // max number of bins
//#define PAUSE_WHEN_MEMORY_LOW 10000 ; // if defined, the program pauses if there is less than PAUSE_WHEN_MEMORY_LOW MB of free memory




#endif  // PARAMS_H_INCLUDED
