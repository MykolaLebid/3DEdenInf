#ifndef CONF_STRUCTURES_H
#define CONF_STRUCTURES_H

#include <string>
#include <iostream>

//---------------conf file---------------------------begin
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/lexical_cast.hpp>
//---------------conf file-----------------------------end


//Parameter descriptions
////////////////////////////////////////////////////////////////////////////////
struct Evolution {
private:
public:
	float driver_adv = 0;
	float mutation_rate = 0;
	float driver_mutation_rate = 0;
	int stop_time_cell_num = 0; // boundary for number of cells in population
	int simulation_regime = 0; // 0 - well_mixed; 1 - 3D; 2 - accurate 3D

	Evolution() : driver_adv(0), mutation_rate(0), driver_mutation_rate(0),
								stop_time_cell_num(0), simulation_regime(0) {};

	Evolution(float driver_adv_ex,
						float mutation_rate_ex,
						float driver_mutation_rate_ex,
						int 	stop_time_cell_num_ex,
						int 	simulation_regime_ex):
		driver_adv(driver_adv_ex), mutation_rate(mutation_rate_ex),
		driver_mutation_rate(driver_mutation_rate_ex),
		stop_time_cell_num(stop_time_cell_num_ex),
		simulation_regime(simulation_regime_ex){};

	friend bool operator==(const Evolution& lhs, const Evolution& rhs)
	{
			if (lhs.mutation_rate == rhs.mutation_rate &&
					lhs.driver_mutation_rate == rhs.driver_mutation_rate &&
					lhs.driver_adv == rhs.driver_adv &&
					lhs.stop_time_cell_num == rhs.stop_time_cell_num &&
					lhs.simulation_regime == rhs.simulation_regime) return true;
			else return false;
	}
	Evolution& operator=(const Evolution& right) {
			if (this == &right) {
					return *this;
			};
			driver_adv = right.driver_adv;
			mutation_rate = right.mutation_rate;
			driver_mutation_rate = right.driver_mutation_rate;
			stop_time_cell_num = right.stop_time_cell_num;
      simulation_regime =right.simulation_regime;
			return *this;
	};
//////////////////////////////////////////////////////////////////////////////
//TODO Delete
//    Evolution operator+(const Evolution& arg){
//			return Evolution(driver_adv + arg.driver_adv,
//								       mutation_rate + arg.mutation_rate,
//											 driver_mutation_rate + arg.driver_mutation_rate,
//								       stop_time_cell_num + arg.stop_time_cell_num);
//
//    };
//    Evolution operator/(const float div) {
//			return Evolution(driver_adv/div,mutation_rate/div,
//											 driver_mutation_rate/div, stop_time_cell_num/div);
//    };
//////////////////////////////////////////////////////////////////////////////
};



//struct DeltaArea {
//    int del_x;
//    int del_y;
//    int del_z;
//};
//
//// structure ProbeArea
//// Define Rectangular cuboid for probe cells
//struct Area {
//    int min_x;  int max_x;
//    int min_y;  int max_y;
//    int min_z;  int max_z;
//};

// structure dynamic
struct Borders {
    int min_x;  int max_x;
    int min_y;  int max_y;
    int min_z;  int max_z;
};


struct TreeCells{
	int use;
	int in_balk;
	int cell_num;
	std::string etalon_result_t_c_file_name;
} ;

struct BalkCells{
	int use;
	int needle_diam;
	std::string etalon_result_b_c_file_name;
};

struct Probes{
    //int is_random;
    TreeCells tree_cells;
    BalkCells balk_cells;
    //DeltaArea delta_area;
    Borders population_borders;
    //int cell_num;

};


struct ComparisonResults {
public:
		float balk_dist;
    float tree_dist;
		float dist;
    Evolution evolution;

    ComparisonResults() : balk_dist(10000), tree_dist(10000), dist(10000){};

    ComparisonResults (float balk_dist_ex, float tree_dist_ex,
											 float dist_ex, Evolution evo_ex):
    	balk_dist(balk_dist_ex), tree_dist(tree_dist_ex), dist(dist_ex),
    	evolution(evo_ex.driver_adv, evo_ex.mutation_rate,
								evo_ex.driver_mutation_rate, evo_ex.stop_time_cell_num,
								evo_ex.simulation_regime){};

//		friend bool operator==(const ComparisonResults& lhs, const ComparisonResults& rhs)
//    {
//				if ((lhs.evolution == rhs.evolution) &&
//				    (lhs.balk_dist == rhs.balk_dist) &&
//						(lhs.tree_dist == rhs.tree_dist) &&
//						(lhs.dist == rhs.dist)) return true;
//				else return false;
//    };
		ComparisonResults& operator=(const ComparisonResults & arg) {
        if (this == &arg) {
            return *this;
        };
				evolution = arg.evolution;
		    dist = arg.dist;
        tree_dist = arg.tree_dist;
        balk_dist = arg.balk_dist;
				return *this;
    };

    ComparisonResults operator+(const ComparisonResults& arg){
			return ComparisonResults(balk_dist + arg.balk_dist,
								               tree_dist + arg.tree_dist,
								               dist + arg.dist,
								               evolution);
    };
    ComparisonResults operator/(const float div) {
//    	balk_dist = balk_dist/div;
//    	tree_dist = tree_dist/div;
//    	dist= dist/div;
			return ComparisonResults(balk_dist/div, tree_dist/div, dist/div, evolution);
    };


    friend std::istream & operator>>(std::istream & s_in,
													    ComparisonResults & comparison_results) {
				if (!s_in.eof()){
						s_in >> comparison_results.evolution.driver_adv >>
										comparison_results.evolution.mutation_rate >>
								    comparison_results.evolution.driver_mutation_rate >>
										comparison_results.balk_dist >>
								    comparison_results.tree_dist >>
								    comparison_results.dist;
				};
				return s_in;
    };



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
    //File System info
    std::string related_path_to_par_file;
    std::string par_file_name;
    std::string etalon_result_catalogue;

    int seed;  // (RAND)
    int threshold_simulation_num;  // number of attempts to make simulation
                                   // that is important for small birth rate
    Probes probes; // approach to make probes
    Evolution evolution;
		//float treatment = 0;

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
		int tree_dist_type;
		int balk_dist_type;
		float proportion_of_tree_dist;
		float threshold_frac_cell_diff;
		float threshold_frac_mut_diff;
		int piece_num;
		bool par_scale_to_one;
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
		void getFromFile(std::string & file_name);
		void setToFile(std::string & file_name);
		void printToConsole();
		void getFromPTree(boost::property_tree::ptree pt);
		boost::property_tree::ptree getPTree();
		std::string getSearchSpaceParsString();
private:
		BasicSimParameters basic_sim_parameters;
		InferenceParameters inference_parameters;
		boost::property_tree::ptree p_tree;
};
////////////////////////////////////////////////////////////////////////////////


//initializes parameters for an etalon simulation
void getEtalonPartParameters(BasicSimParameters & basic_sim_parameters);

//initializes parameters for an inference simulation
void getInferencePartParameters(BasicSimParameters & basic_sim_parameters,
																InferenceParameters & inference_parameters);

// initializes parameters for inference algorithm

void getAllParameters(BasicSimParameters & basic_sim_parameters,
											InferenceParameters & inference_parameters);

std::string getSystemRelevantSlash();
std::string getSystemRelevantInit();
std::string getSystemRelevantTwoUpperRoot();
std::string getSystemRelevantOneUpperRoot();

std::string getProgName();
std::string getABCProgName();
std::string getCatalogueProgName();
std::string getCatalogueABCProgName();
std::string getSearchSpaceParsString();

void err(char *reason);
void err(const char *reason); // exits program with correspondent reason
bool fileExists(const char *name);
bool fileRemoveIfExists(const char *name);
#endif // CONF_STRUCTURES_H
