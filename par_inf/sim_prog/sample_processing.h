// Copyright 2017 Lebid Mykola all rights reserved

#ifndef SAMPLE_PROCESSING_H_INCLUDED
#define SAMPLE_PROCESSING_H_INCLUDED

//#include <stdio.h>
//#include <math.h>

#include <fstream>
#include <iostream>
//#include <stdexcept>  // out_of_range

#include <vector>

#include "conf_structures.h"
#include "classes.h"

//class Genotype;
//struct Cell;


//extern std::vector<Genotype*> genotypes;

// contains: max, min borders in x, y, z directions;
//           sweeps in x (del_1), y (del_2), z (del_3)
//           directions of the evolutionary population
//struct PopulationRange {
//    int min_x;  int max_x;
//    int min_y;  int max_y;
//    int min_z;  int max_z;
//    int del_1;  int del_2;  int del_3;
//};


// takes all cells in Rectangular cuboid @param Probe Pr and
// print them in a file with a name @param char *name;
// creates prob;
// write description of parameters or better to write structure for parameters
//float TakeProbe(char *name, DeltaProbeArea & pop_range,
//                 ofstream & file_tree_comparison,
//                 float max_dist_for_prob);
//float TakeProbe(Parameters & Pars);
void doesComparativeAnalysis(AllSimParameters &, SimData &);

void GetComparisonResultsFromFile(const AllSimParameters &,
															    std::vector<ComparisonResults> &);
void MakeAverageFile(const AllSimParameters &);

// takes an probe (sample) to build a tree
//
//void setEtalonProbe(char *file_name_part, ProbeArea &,
//                    unsigned int threshold_num_cells);
//void InitProbeArea(BasicSimParameters & pars);
void setEtalonProbe(BasicSimParameters &, SimData &);
///////////////////////////////////////////////////////////////////////////////
// structure CompareCells
// supports exclusion algorithm
// adds number of cells in unique genotype after exclusion
struct CompareGenotypesEq {
    inline bool operator() ( Genotype * a, Genotype * b){
        unsigned int size_a= a->sequence.size();
        unsigned int size_b= b->sequence.size();
        unsigned int last_a = -1;
        unsigned int last_b = -1;
        if (size_a !=0 ) last_a = a->sequence.at(size_a-1);
        if (size_b !=0 ) last_b = b->sequence.at(size_b-1);
        if ((size_a == size_b) && (last_a == last_b)) {
            a->number++; //b->number=2;
            return true;
        } else return false;
    };
};/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// structure included_genotype
// support unique algo for exclusion probe_cells by deletion the same cells
struct included_genotype {
inline bool operator() ( Genotype * a, Genotype * b){
    unsigned int size_a= a->sequence.size();
    unsigned int size_b= b->sequence.size();
    int last_a = -1;
    int last_b = -1;
    if (size_a !=0 ) last_a = a->sequence.at(size_a-1);
    if (size_b !=0 ) last_b = b->sequence.at(size_b-1);
    if ((size_a == size_b) && (last_a == last_b)) {
        a->number++; //b->number=2;
        return true;
    } else return false;
};
};/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//Structure MutGenotypeConnection
// we use this  node structure after:
//1) sorting genotype by length @function UniqueGenotypes();
//2) after exclusions of inclusions @function GenotypesWithoutInclusions()//
//example:
//Genotype B= probe_genotypes[6].sequence=" p 51 p 789 p 3465 p 12456"
//Genotype A= probe_genotypes[7].sequence=" p 51 p 789 p 3465 d 4570 p 10020"
// -> mut_genotypes_connections[6].num_genotype=7
// -> mut_genotypes_connections[6].num_mutation=2
//Stub {-1, -1}; // Stub ofstd::vector mut_genotypes_connections
struct Connect{
  // genotype B which we attach to A
  int num_genotype;
  // the place in line genotype A  where we attach the genotype B
  int num_mutation;
};
struct MutGenotypeConnection {
  // genotype A to which we attach correspondent (smaller) genotype B
  int num_genotype;
  // the place in line genotype A  where we attach the genotype B
  int num_mutation;
	std::vector <Connect> vec_index_attached_genotypes;//---------------------
																						// we save all indexes of attached
                                            // genotypes init in function
																						// GetNumChildrenCurNode(,)
};/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Structure MutationNode
// we use this node structure for:
// 1) comparison of trees;
// 2) all mutations is located in one vector and they are arranged
//    by the depth of mute
struct MutationNode {
    int abs_mut; // absolute number of this mutation (-1 for the root mut)
    unsigned int depth; // distance to the root mutation (begins with 0!!!)
    int abs_father_mut; // (-1 for the root mut)
    int index_father_mut;  //index in mutation_vector(zero_mut=-1)
    unsigned int num_children_mut;
    unsigned int number_cells;  // number of cells
                                // with this mutation in the end of genotypes
    std::vector <unsigned int> vec_child_places; // in vector mutation_vector
		//float frec;
};
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Structure MutationNode
// we use this structure for:
// 1) comparison of trees;
// 2) hub is mutation node or set of subsequent nodes with one of the following properties
//    a) node contains a cell
//    b) node has multiple ancestors
// 3) all hubs is located in one vector and they are arranged in a mut_hub vector
// by the depth of correspondent mutations
struct MutationHub {
	int abs_mut; // abs_number of the last mutation in a seq of subsequent nodes
  std::vector <unsigned int> label_vec; // we can have several labels adjusted to the hub
																	 // (each label is a cell)
																	 // size of label_vec should be number_cells
  unsigned int number_cells;  // number of cells in the hub
  unsigned int depth; // distance to the root mutation (begins with 0!!!)
  unsigned int node_number; // number of nodes in a Hub
  int abs_father_hub; // (-1 for the root mut)
	int index_father_hub;  // index in mut_hub_vector(zero_mut=-1)
	unsigned int num_children_hubs; //

	//vector <unsigned int> vec_child_places; // in vector mutation_hub_vector
};
///////////////////////////////////////////////////////////////////////////////

struct ProbePiece {
	std::vector <Cell> tree_cell_vector;
  std::vector <Cell> vaf_cell_vector;
 // int indicator = 0;
  ComparisonResults results;
};
///////////////////////////////////////////////////////////////////////////////
// ProbeMutProcessor
// Class take vector of sells, all sim genotypes and work with mutations
class TreeMutAnalyzer{
  public:
  //Constructor:
  //@param const vector <Cell> & - Cells in Probe
  TreeMutAnalyzer();
	//TreeMutAnalyzer(const std::vector <Cell> &);
  //TreeMutAnalyzer(const std::vector <Cell> &, unsigned int threshold_numb_cells);

  //genotypes came from simulations
  TreeMutAnalyzer(const std::vector <Cell> & _probe_cells,
									const std::vector <Genotype *> & _genotypes);
  //Destructor:
  ~TreeMutAnalyzer();

  //we use only mutation_vector() for these methods
  unsigned int get_cell_number() const;
  unsigned int get_mut_number() const;
  unsigned int get_tree_mut_depth() const; // max
  unsigned int get_tree_mut_degree() const;

  //vector <MutationNode> GetMutationNode(unsigned int index) const;
  std::vector <unsigned int> get_depth_mut_vec() const;
  std::vector <unsigned int> get_degree_mut_vec() const;
  //finds mean num of muts in probe cells
  float get_mean_mut_num_cell() const;
  float get_mean_pass_b_w_non_elem_nodes() const;
  //TODO change to private
  std::string get_str_newick_mut_tree_format() const;
  /////////////////////////////////////////////////////////////////////////////
  private:
  // Member variables

 //  main cells from outside:
  std::vector <Cell> probe_sim_cells;

// main cells () after reformation of gen: from (Genotype *) to (unsigned int)
  std::vector <SelectedCell> probe_cells;
//  std::vector <ProbeCell> raw_probe_cells;
  std::vector <MutationNode> mutation_vector;
  std::vector <MutationHub> mut_hub_vector;
  //genotype of main cells
  std::vector <Genotype *> probe_genotypes;
  //all genotypes from simulations
  std::vector <Genotype *> genotypes;
  //@vector: see description of struc MutGenotypeConnection:
  std::vector <MutGenotypeConnection> mut_genotypes_connections; ///////////////
  // @vector: needs for sub-genotypes
  // example: A and B in probe_genotypes
  // Genotype B= " p 51 p 789 p 3465 "
  // Genotype A= " p 51 p 789 p 3465 d 4570 p 10020"
  // then B is a sub-genotype of A
  std::vector <Genotype *> excluded_probe_genotypes;
  /////////////////////////////////////////////////////////////////////////////


  // Support functions for constructor initialization:
  void Constructor_init();  // main body of constructor. There are
                           // sequence of functions in Constructor_init():
  /////////////////////////////////////////////////////////////////////////////
  void setCellsGenotypes();    //  saves sorted vector of genotypes to
                               //  probe_genotypes (sorting was made
                               //  indirectly by sorting of cells by
                               //  lengths of genotype in constructor)
  void uniqueProbeCells(); // left cells with diff.genotypes
  void uniqueGenotypes(); // leaves unique genotypes and
                          // write new number of cells w.r.t. a genotype

  void genotypesWithoutInclusions();//after 1) void UniqueGenotypes();
  /////////////////////////////////////////////////
	/// \brief @function creatMutGenotypesConnections()
  /// needs in constructor after 1) void UniqueGenotypes() and
  /// 2) void GenotypesWithoutInclusions();
  /// creats vector mut_genotypes_connections;
  /// this vector describes tree like structure with preprocessed
  /// probe_genotypes vector
  /// for the longest genotype we find the longest intersection
  /// with small and remember the place in mut_genotypes_connections
  /// \return void
  ///
  /////////////////////////////////////////////////
  void creatMutGenotypesConnections();
  /////////////////////////////////////////////////////////////////////////////
  void creatMutVector();
  int GetNumCellsForCurNode(int index_cur_genotype_mut, int abs_mut);
  int GetFatherMutIndexForCurNode( int index_cur_probe_genotype,
                                   int index_cur_genotype_mut,
                                   int initial_place);
  MutationNode GetCurrentMutNode( int index_cur_probe_genotype,
                                  int index_cur_genotype_mut,
                                  int initial_place);
  int GetNumChildrenCurNode( int index_cur_probe_genotype,
                             int index_cur_genotype_mut);//and attach
  void InitZeroMutationNode(MutationNode& zero_mutation);
  void FillVecChildPlaces();
	///////////////////////////////////////////////////////////////////////////

	//form mut_hub_vector
	///////////////////////////////////////////////////////////////////////////
  void creatMutHubVector();
  void ChangeDepth(std::vector <MutationNode> & ,
													 unsigned int);
	void SetIndexFatherHubVec();
	//TODO rewrite following two functions (They are similar)
	void AdjustLabels2HubVecUsingMut();
  void AdjustLabels2HubVecUsingHub();
  ///////////////////////////////////////////////////////////////////////////

	std::vector <std::vector<int>> GetSPLMatrix() const;// Splitted path lengths
  unsigned int GetDistBetweenCells(unsigned int, unsigned int) const;
  ///////////////////////////////////////////////////////////////////////////
  //Add functions

  //string get_string_format_mut(unsigned int number);
  //Accessor functions
  int get_number_of_intrsections(Genotype * ,Genotype *);

  // File functions
  friend std::ostream & operator<<(std::ostream & output, const
                                     TreeMutAnalyzer & ProbeMutProcessor_);
  friend std::istream & operator>>(std::istream & s_in,
                                     TreeMutAnalyzer & ProbeMutProcessor_);
  // different tree metrics
  friend float t_dist_deg_vs_dep(const TreeMutAnalyzer &,
																 const TreeMutAnalyzer &);  // dist_type = 0

  // dist_type = 1
  friend float t_dist_deg_vs_mean_pass_non_elem (const TreeMutAnalyzer &,
																		             const TreeMutAnalyzer &);
  // dist_type = 2
  friend float t_dist_deg_vs_mean_pass_to_cell (const TreeMutAnalyzer &,
																		            const TreeMutAnalyzer &);

	friend float t_dist_mean_pass_to_cell(const TreeMutAnalyzer &,
																				const TreeMutAnalyzer &);
  // additional functions for tree metrics
  friend float chi_square_depth(const TreeMutAnalyzer &,
																const TreeMutAnalyzer &);

  friend float chi_square_degree(const TreeMutAnalyzer &,
																 const TreeMutAnalyzer &);

  friend float dist_mean_pass_b_w_non_elem_nodes(const TreeMutAnalyzer &,
                                                 const TreeMutAnalyzer &);
  friend float dist_mean_pass_to_cell(const TreeMutAnalyzer &,
																			const TreeMutAnalyzer &);

	friend float dist_branch_score(const TreeMutAnalyzer &,  // dist of Kuhner and
																 const TreeMutAnalyzer &, int seed); // Felsenstein (1994)
	// The Branch Score Distance uses branch lengths,
	// and can only be calculated when the trees have lengths on all branches
	friend float dist_robinson_and_foulds(const TreeMutAnalyzer &,
																				const TreeMutAnalyzer &, int seed);	// friend float

	friend float nodal_distance (const TreeMutAnalyzer &,
												       const TreeMutAnalyzer &);

  friend float splitted_nodal_distance (const TreeMutAnalyzer &,
																				const TreeMutAnalyzer &);
	// dist_type = 8
	friend float dist_splitted_nodal (const TreeMutAnalyzer &,
																		const TreeMutAnalyzer &);
	// dist_type = 9
	friend float dist_nodal (const TreeMutAnalyzer &,
													 const TreeMutAnalyzer &);
	// dist_type = 10
	friend float t_dist_deg_vs_dist_nodal(const TreeMutAnalyzer &,
													              const TreeMutAnalyzer &);

	// dist_type = 11
	friend float t_dist_deg (const TreeMutAnalyzer & Probe_mut_1,
									         const TreeMutAnalyzer & Probe_mut_2);
	// dist_type = 12
	friend float t_dist_control (const TreeMutAnalyzer & Probe_mut_1,
															 const TreeMutAnalyzer & Probe_mut_2);
	// dist_type = 13
	friend float t_dist_depth(const TreeMutAnalyzer & Probe_mut_1,
									          const TreeMutAnalyzer & Probe_mut_2);
	// dist_type = 14
	friend float t_dist_depth_cum_kolmogorov_smirnov(const TreeMutAnalyzer & Probe_mut_1,
															                     const TreeMutAnalyzer & Probe_mut_2);
	// dist_type = 15
	friend float t_dist_depth_cum_kuiper(const TreeMutAnalyzer & Probe_mut_1,
																			 const TreeMutAnalyzer & Probe_mut_2);
	// dist_type = 16
	friend float t_dist_depth_cum_cramer_mises(const TreeMutAnalyzer & Probe_mut_1,
																				     const TreeMutAnalyzer & Probe_mut_2);
	// dist_type = 17
	friend float t_dist_deg_cum_kolmogorov_smirnov(const TreeMutAnalyzer & Probe_mut_1,
															                   const TreeMutAnalyzer & Probe_mut_2);
	// dist_type = 18
	friend float t_dist_deg_cum_kuiper(const TreeMutAnalyzer & Probe_mut_1,
																		 const TreeMutAnalyzer & Probe_mut_2);
	// dist_type = 19
	friend float t_dist_deg_cum_cramer_mises(const TreeMutAnalyzer & Probe_mut_1,
																				   const TreeMutAnalyzer & Probe_mut_2);
	// dist_type = 20
	friend float t_dist_deg_depth_sum_cum_cramer_mises(const TreeMutAnalyzer & Probe_mut_1,
																				             const TreeMutAnalyzer & Probe_mut_2);
	// dist_type = 21
	friend float t_dist_deg_mean(const TreeMutAnalyzer & Probe_mut_1,
															 const TreeMutAnalyzer & Probe_mut_2);
  // dist_type = 22
  friend float t_dist_depth_mean(const TreeMutAnalyzer & Probe_mut_1,
															   const TreeMutAnalyzer & Probe_mut_2);
  // dist_type = 23
  friend float dist_zhao_michor(const TreeMutAnalyzer & Probe_mut_1,
																const TreeMutAnalyzer & Probe_mut_2);

protected:
  /////////////////////////////////////////////////////////////////////////////
  //Support connector functions
  void setProbeCells(const std::vector <Cell> & _outside_probe_cells);

  /////////////////////////////////////////////////////////////////////////////
  //Support functions for file writing
  void SaveProbeCells(const char *);
  void SaveGenotypes(char *, std::vector<Genotype*> &);
  void SaveMutGenotypeConnections(char *);
  void SaveMutVector(char *);
  void SaveMutHubVector(char *);

  void SaveGraph(const char * catalog_name, //create a tex file
								 const char * file_name);   //with simple dot graph construction

  void SaveHubGraph(const char * catalog_name, //create a tex file
								    const char * file_name);   //with simple dot graph
								                               //hub construction

  void SaveNewickTreeFormat(char *);
	// recursive function for tree Newick format
	std::string GetCurrentMutationString(int mut_vec_index) const;
	std::string FormCellString(int num_cells) const; // for tree Newick format
  /////////////////////////////////////////////////////////////////////////////
  int current_cell_index = 0;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
class BalkVAFAnalyzer{
public:
	//Constructor
	BalkVAFAnalyzer();
	BalkVAFAnalyzer(unsigned int _num_of_cells, unsigned int _num_of_muts,
									const std::vector<unsigned int> & _num_mut_vec);
	BalkVAFAnalyzer(const std::vector <Genotype*> & selected_cell_genotypes);
	BalkVAFAnalyzer(const std::vector <Cell> & cell_vec,
							    const std::vector <Genotype*> & cell_genotypes);
//void PrintToFileCumulativeDist(std::ostream & output);

	//get methods
	unsigned int getCellNum()	const {return num_of_cells;};
	unsigned int getMutNum()	const {return num_of_muts;};
	std::vector <float> getCumulativeVafVec()	const {return cumulative_vaf_vec;};
	std::vector <unsigned int> getNumMutVec()	const {return num_mut_vec;};

	//void free();
	//Destructor
	~BalkVAFAnalyzer();
private:
	// Member variables
	unsigned int num_of_cells;
	unsigned int num_of_muts;
	std::vector <unsigned int> num_mut_vec;
	std::vector <float> cumulative_vaf_vec;

	// Support functions
	// void AddOneGenotype(Cell one_cell);
	void initCumulativeVAFVec();
};

float KolmogorovSmirnovDist(const BalkVAFAnalyzer & balk_vaf_1,
														const BalkVAFAnalyzer & balk_vaf_2);
float KuiperDist(const BalkVAFAnalyzer & balk_vaf_1,
								 const BalkVAFAnalyzer & balk_vaf_2);
float CramerMisesDist(const BalkVAFAnalyzer & balk_vaf_1,
											const BalkVAFAnalyzer & balk_vaf_2);
std::ostream & operator<<(std::ostream & output,
																 const BalkVAFAnalyzer & balk_vaf_analyzer);
std::istream & operator>>(std::istream & s_in,
																	 BalkVAFAnalyzer & balk_vaf_analyzer);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class CellMutAnalyzer{
public:
	//Constructor
	CellMutAnalyzer();
	CellMutAnalyzer(const std::vector <Cell> & cell_vec,/*vector<Genotype*> & class_genotypes,*/
								                             /*unsigned int total_number_of_SNV,*/
									const std::vector<Genotype*> & cell_genotypes);
  void initSelectedCellVec(const std::vector <Cell> & _cell_vec);

	friend std::ostream & operator<<(std::ostream & output,
																     const CellMutAnalyzer & cell_mut_analyzer);
	friend std::istream & operator>>(std::istream& s_in,
                                     CellMutAnalyzer & cell_mut_analyzer);
	std::vector <Genotype> getGenotypeVec() const {return genotype_vec;};
	std::vector <SelectedCell> getCellVec() const { return cell_vec;};
	void CleanAll(){
		genotype_vec.clear();
		cell_vec.clear();
	};
	void setGenotypeVec (std::vector <Genotype> & _genotype_vec){
		genotype_vec = std::move(_genotype_vec);
	};
	  //Support connector functions
	void setCellVec(const std::vector <SelectedCell> & _outside_probe_cells);

	void MutationFilter(const unsigned int left_border,
												const unsigned int right_border);
	std::vector <unsigned int> GetMutCellFr();
	std::vector <unsigned int> GetCellMutFr();
	friend float ZhaoMichorDist(CellMutAnalyzer cell_mut_probe_1,
							                  CellMutAnalyzer cell_mut_probe_2,
							                  unsigned int left_border,
							                  unsigned int right_border);
	unsigned int GetCellNum(){return cell_vec.size();};
	 //Destructor
	~CellMutAnalyzer();
private:
		// Member variables
	std::vector <Genotype>   genotype_vec;
  std::vector <SelectedCell> cell_vec;

		// Support functions
};


////////////////////////////////////////////////////////////////////////////////
#endif // SAMPLE_PROCESSING_H_INCLUDED
