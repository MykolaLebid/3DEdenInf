#ifndef ETALON_CONF_STRUCTURES_H
#define ETALON_CONF_STRUCTURES_H

//#include <string>
#include <iostream>
//
//---------------conf file---------------------------begin
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/lexical_cast.hpp>
//---------------conf file-----------------------------end
//Mykola begin
//Parameter descriptions
////////////////////////////////////////////////////////////////////////////////
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

//struct ProbeComparisonParameters {
//    std::string comparison_file_name;
//    std::string e_mut_file_name;
//    std::string catalogue_with_e_mut_file;
//    int e_seed;
//    Probe probe_pars;  // for etalon
//    Evolution evo_pars;  // for etalon
//    ComparisonResults comparison_results;
//    unsigned int dist_type;
//    float max_dist_probe_trees;
//    float threshold_frac_cell_diff;  // threshold frac cell number discrepancy
//    float threshold_frac_mut_diff;  // see above (for mutation)
//    int piece_num; // number of probe pieces (attempts)
//};
//

struct BasicSimParameters{
    std::string related_path_to_par_file;
    std::string par_file_name;

    std::string etalon_result_catalogue;
    std::string etalon_result_file_name;

    int seed;  // (RAND)
    int threshold_simulation_num;  // number of attempts to make simulation
                                   // that is important for small birth rate
    Probe probe;
    Evolution evolution;
 //   ProbeComparisonParameters etalon;
};
//------------------------------------------------------------------------------
struct LogScaleSearchSpace {
		float driver_advantage; // that is a log scale of a search range
														// Example: LogScaleSearchSpace A;
														// then we begin to search a parameter of driver_advantage in a range
														// [true_value/10^(A.driver_advantage);true_value*10^(A.driver_advantage)]
		float mutation_rate;
		float driver_mutation_rate;
};

struct AlgorithmParameters {
		int final_iteration;
		float final_epsilon;
		float p_acc_min;
		int parameter_population_number;
		int parameter_population_num_alpha;
};

struct ComparisonParameters {
		int dist_type;
		float threshold_frac_cell_diff;
		float threshold_frac_mut_diff;
		int piece_num;
};

struct BasicSettings {
		std::string current_result_file_name;
		std::string log_file_name;
};

struct InferenceParameters{
		LogScaleSearchSpace log_scale_search_space;
    AlgorithmParameters algorithm_parameters;
		ComparisonParameters comparison_parameters;
		BasicSettings basic_settings;
};


struct AllSimParameters{
		BasicSimParameters basic_sim_parameters;
		InferenceParameters inference_parameters;
};

//------------------------------------------------------------------------------
class ParameterStructClass{
public:
		void GetFromFile(std::string & file_name);
		void SetToFile(std::string & file_name);
		void PrintToConsole();
		void GetFromPTree(boost::property_tree::ptree pt);
		boost::property_tree::ptree getPTree();

private:
		BasicSimParameters basic_sim_parameters;
		InferenceParameters inference_parameters;
		boost::property_tree::ptree p_tree;
};
////////////////////////////////////////////////////////////////////////////////

//initializes parameters for an etalon simulation
void GetEtalonPartParameters(BasicSimParameters & basic_sim_parameters);

//initializes parameters for an inference simulation
void GetInferencePartParameters(BasicSimParameters & basic_sim_parameters,
																InferenceParameters & inference_parameters);


// initializes parameters for inference algorithm

void GetAllParameters(BasicSimParameters & basic_sim_parameters,
											InferenceParameters & inference_parameters);

std::string GetSystemRelevantSlash();
std::string GetSystemRelevantInit();
std::string GetSystemRelevantTwoUpperRoot();
std::string GetSystemRelevantOneUpperRoot();

std::string GetProgName();
std::string GetABCProgName();
std::string GetCatalogueProgName();

void err(char *reason);
void err(const char *reason); // exits program with correspondent reason
bool FileExists(const char *name);

#endif // ETALON_CONF_STRUCTURES_H
