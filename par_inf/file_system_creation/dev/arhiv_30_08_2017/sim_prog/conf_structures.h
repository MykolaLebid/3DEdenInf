#ifndef ETALON_CONF_STRUCTURES_H
#define ETALON_CONF_STRUCTURES_H

#include <string>
#include <iostream>

//---------------conf file---------------------------begin
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

//---------------conf file-----------------------------end


//Mykola begin
//Parameter descriptions
///////////////////////////////////////////////////////////////////////////////
struct Evolution {
    float driver_adv;
    float mutation_rate;
    float driver_mutation_rate;
    // boundary for number of cells in population
    int stop_time_cell_num;
};


struct DeltaArea {
    int del_x;
    int del_y;
    int del_z;
};
// structure ProbeArea
// Define Rectangular cuboid for probe cells
struct Area {
    int min_x;  int max_x;
    int min_y;  int max_y;
    int min_z;  int max_z;
};

// structure dynamic
struct Borders {
    int min_x;  int max_x;
    int min_y;  int max_y;
    int min_z;  int max_z;
};

struct Probe{
    int is_random;
    Area area;
    DeltaArea delta_area;
    Borders population_borders;
    int cell_num;
};

struct ComparisonResults {
    float dist = 10000;
    int num_mut = - 1;
    int e_num_mut = -1;
    int delta_mut = - 1;
};

struct ProbeComparisonParameters {
    std::string comparison_file_name;
    std::string e_mut_file_name;
    std::string cataloge_with_e_mut_file;
    int e_seed;
    Probe probe_pars;  // for etalon
    Evolution evo_pars;  // for etalon
    ComparisonResults comparison_results;
    unsigned int dist_type;
    float max_dist_probe_trees;
    float threshold_frac_cell_diff;  // threshold frac cell number discrepancy
    float threshold_frac_mut_diff;  // see above (for mutation)
    int piece_num; // number of probe pieces (attempts)
};


struct Parameters {
    std::string par_file_name;
    std::string file_name;
    std::string catalog_file_name;
    std::string log_file;

    int seed;  // (RAND)
    int threshold_simulation_num;  // number of attempts to make simulation
                                   // that is important for small birth rate
    bool is_simulation_for_etalon_probe; // there is 0 - for trials
    Probe probe;
    Evolution evolution;
    ProbeComparisonParameters etalon;

		void GetFromFile(std::string & file_name_);
		void SetToFile(std::string & file_name_);
		void PrintToConsole();
		boost::property_tree::ptree GetPTree();
		void GetFromPTree(boost::property_tree::ptree pt);

		private:
		boost::property_tree::ptree p_tree;

};
///////////////////////////////////////////////////////////////////////////////


#endif // ETALON_CONF_STRUCTURES_H
