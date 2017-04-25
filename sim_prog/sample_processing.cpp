// Copyright 2017 Lebid Mykola all rights reserved


#include "sample_processing.h"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>


#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <memory>
#include <numeric>

#include "classes.h"
//#include "params.h"
//#include "params.h"
#include <assert.h>

using namespace std;

void ProbeMutProcessor::UniqueProbeCells() {
  vector<Cell>::iterator it = unique (probe_cells.begin(), probe_cells.end(),
                                      CompareCellsEq());
  probe_cells.resize(std::distance(probe_cells.begin(),it));
};


unsigned int ProbeMutProcessor::get_cell_number() const {
  int cell_num = 0;
  for (unsigned int i=0; i<mutation_vector.size(); i++) {
    cell_num += mutation_vector[i]. number_cells;
  };
  return cell_num;
};

unsigned int ProbeMutProcessor::get_mut_number() const {
 return mutation_vector.size();
};

void ProbeMutProcessor::UniqueGenotypes() {
  vector<Genotype *>::iterator it = unique (probe_genotypes.begin(),
                                            probe_genotypes.end(),
                                            CompareGenotypesEq());
  probe_genotypes.resize( std::distance(probe_genotypes.begin(),it) );
};

void ProbeMutProcessor::SetCellsGenotypes() {
    for (unsigned int i=0;i<probe_cells.size(); i++) {
            probe_genotypes.push_back(genotypes[probe_cells[i].gen]);
            probe_genotypes[i]->number = 1;
    };
};
void ProbeMutProcessor::GenotypesWithoutInclusions() {
  for(unsigned int i = 0; i < (probe_genotypes.size()-1); i++) {
    bool delete_indicator = false;
    for(unsigned int j = i + 1; (j < probe_genotypes.size()) &&
                       (delete_indicator == false); j++) {
      bool is_included = includes( probe_genotypes[j]->sequence.begin(),
                                   probe_genotypes[j]->sequence.end(),
                                   probe_genotypes[i]->sequence.begin(),
                                   probe_genotypes[i]->sequence.end());
      if (is_included) {
        excluded_probe_genotypes.push_back(probe_genotypes[i]);
        probe_genotypes.erase(probe_genotypes.begin() + i);
        delete_indicator = true;
        i--;
      };
    };
  };
};


int ProbeMutProcessor::get_number_of_intrsections(Genotype * a, Genotype * b) {
    int min_ab = (a->sequence.size() > b->sequence.size()) ?
                  b->sequence.size() : a->sequence.size();
    int counter = -1;
    for(int i=0;i<min_ab;i++) {
        if (a->sequence[i] == b->sequence[i]) {
          counter++;
        } else {
          break;
        };
    };
    return counter;
}

void ProbeMutProcessor::CreatMutGenotypesConnections() {
MutGenotypeConnection last={-1, -1}; // Stub of vec mut_genotypes_connections
mut_genotypes_connections.push_back(last);
for(int i = (probe_genotypes.size() - 2); i >= 0; i --) {
    // the maximal number of intersections
    // of mutations between genotype with index i and other
    // genes with indexes [i+1,probe_genotypes.size()-1]
    int i_max_intersection = -1;
    ///////////////////////////////////////////////////////////////////////////
    int num_genotype = - 1;
    for(int j = probe_genotypes.size() - 1; j >= i + 1; j--){
        //can return -1
        int intersect_num = get_number_of_intrsections(probe_genotypes[i],
                                                       probe_genotypes[j]);
        if (i_max_intersection < intersect_num) {
                num_genotype = j;
                i_max_intersection = intersect_num;
        };
    };
    MutGenotypeConnection con = {num_genotype, i_max_intersection};
    mut_genotypes_connections.push_back(con);

};
//order correction
reverse(mut_genotypes_connections.begin(),mut_genotypes_connections.end());
};

void ProbeMutProcessor::SaveMutGenotypeConnections(char* name) {
  char name_file[256];
  sprintf(name_file,"%s/Probes_mut_genotypes_connections.dat",name);
  FILE *f=fopen(name_file, "w") ;
  if (f==NULL) err("err");
  int size_mgc = mut_genotypes_connections.size();
  for(int i = 0 ;i < size_mgc; i++) {
    fprintf

    (f,"num_of_attachted_genotype=%d num_genotype=%d num_mutation=%d\n",
            i, mut_genotypes_connections[i].num_genotype,
           mut_genotypes_connections[i].num_mutation );
  };
};

void ProbeMutProcessor::SaveMutVector(char* name) {
  char name_file[256];
  sprintf(name_file,"mutation_vector_%s.dat",name);
  FILE *f=fopen(name_file, "w") ;
  if (f==NULL) err("err");

  int size_mv = mutation_vector.size();
  for(int i = 0 ;i < size_mv; i++) {
    fprintf(f,"%d  abs_mut=%d depth=%d abs_father_mut=%d index_father_mut=%d ",
            i,
            mutation_vector[i].abs_mut, mutation_vector[i].depth,
            mutation_vector[i].abs_father_mut,
            mutation_vector[i].index_father_mut);
    fprintf(f,"num_children_mut=%d number_cells=%d\n",
            mutation_vector[i].num_children_mut,
            mutation_vector[i].number_cells);
  };
};

void ProbeMutProcessor::SaveGraph(char * name){
    char name_file[256];
    sprintf(name_file,"%s/mutation_graph.txt",name);
    FILE *f=fopen(name_file, "w") ;
    if (f==NULL) err("err");
    int size_mv = mutation_vector.size();
    fprintf(f,"digraph graphname {\n");
    for(int i = 0 ;i < size_mv; i++) {
      unsigned int sun = mutation_vector[i].abs_mut;
      unsigned int father = mutation_vector[i].abs_father_mut;
      if (sun!=0){
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if ((father & RESISTANT_PM)==0) {
          if ((father & DRIVER_PM)==0) {
            fprintf(f,"p%u",father);
          } else {
            fprintf(f,"d%u",(father & (~DRIVER_PM)));
          };
        } else {
          fprintf(f,"r%u",(father & (~RESISTANT_PM)));
        };
       //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if ((sun & RESISTANT_PM)==0) {
          if ((sun & DRIVER_PM)==0) {
            fprintf(f," -> p%u",sun);
          } else {
            fprintf(f," -> d%u",(sun & (~DRIVER_PM)));
          };
        } else {
          fprintf(f," -> r%u",(sun & (~RESISTANT_PM)));
        };
        fprintf(f,"\n");
      };
    };
    fprintf(f,"}\n");
}


// calculates number of cells with current mutations
// can be optimize with "find" in sorted vector
int ProbeMutProcessor::GetNumCellsForCurNode(int current_genotype_mut,
                                             int abs_mut) {
  int sum = 0;
  // searches in probe_genotypes (without inclusion elements)
  // can be optimize with "find" in sorted vector
  unsigned int bound = abs(current_genotype_mut + 1);
  for(unsigned int k=0; k < probe_genotypes.size(); k++){
    if (probe_genotypes.at(k)->sequence.size() > bound) break;
    if (probe_genotypes.at(k)->sequence.size() == bound) {
      unsigned int x = abs(abs_mut);
			if ((probe_genotypes.at(k)->sequence.at(current_genotype_mut) == x) &&
					(abs_mut>-1)) sum += probe_genotypes.at(k)->number;
    };
  };
  // searches in excluded_probe_genotypes (in inclusion elements)
  // can be optimize with "find" in sorted vector
  for(unsigned int k = 0; (k < excluded_probe_genotypes.size()); k++){
    if (excluded_probe_genotypes.at(k)->sequence.size() > bound) break;
    if (excluded_probe_genotypes.at(k)->sequence.size() == bound) {
      unsigned int x = abs(abs_mut);
      if ((excluded_probe_genotypes.at(k)->sequence.at(current_genotype_mut) == x) &&
          (abs_mut > -1)) sum += excluded_probe_genotypes.at(k)->number;
		};
  };
  return sum;
};
// calculates number of children of current mutation in current genotype
int ProbeMutProcessor::GetNumChildrenCurNode(int index_cur_probe_genotype,
                                             int index_cur_genotype_mut) {
  // equals at list 1 if in current genotype mut is not the last one
  unsigned int bound = index_cur_genotype_mut + 1;
  int sum = (probe_genotypes[index_cur_probe_genotype]->
             sequence.size() > bound) ? 1 : 0;
//Structure MutGenotypeConnection
// we use this  node structure after:
//1) sorting genotype by length @function UniqueGenotypes();
//2) after exclusions of inclusions @function GenotypesWithoutInclusions()//
//example:
//Genotype B= probe_genotypes[6].sequence=" p 51 p 789 p 3465 p 12456"
//Genotype A= probe_genotypes[7].sequence=" p 51 p 789 p 3465 d 4570 p 10020"
// -> mut_genotypes_connections[6].num_genotype=7
// -> mut_genotypes_connections[6].num_mutation=2
//Stub {-1, -1}; // Stub of vector mut_genotypes_connections
  for(int i = (mut_genotypes_connections.size() - 1) ; i >= 0 ; i--) {
    bool is_genotype_consistent = (mut_genotypes_connections[i].num_genotype ==
                                   index_cur_probe_genotype);
    bool is_mut_in_genotype_consistent =
           ( mut_genotypes_connections[i].num_mutation ==
             ( index_cur_genotype_mut) );
    if ( is_genotype_consistent && is_mut_in_genotype_consistent) sum++;
  };

  return sum;
};

// Attention!!! depth begin with 0 !!!
void ProbeMutProcessor::InitZeroMutationNode(MutationNode& zero_mutation) {
 zero_mutation.abs_mut= - 1;
 zero_mutation.abs_father_mut = - 1;
 zero_mutation.depth = 0;//don't forget about shift in depth
 zero_mutation.index_father_mut= - 1;
 //Calculate number of cells seq without mutations (only basic mutation)
 zero_mutation.number_cells = 0;
 // explanations: that is important to have this row of conditions
 // (separation "of" conditions: the second condition can give out of range)!!!
 if (probe_genotypes.size() > 0) {
   if (probe_genotypes[0] -> sequence.size() == 0) {
     zero_mutation.number_cells = probe_genotypes[0]->number;
   } else {
    // explanations: read previous comment
     if (excluded_probe_genotypes.size()>0) {
       if (excluded_probe_genotypes[0]->sequence.size() == 0) {
        zero_mutation.number_cells = excluded_probe_genotypes[0]->number;
       };
     };
   };
 };
 //Calculate number of children of zero mutation
 zero_mutation.num_children_mut = 0;
 for (unsigned int i = 0; i < mut_genotypes_connections.size(); i++) {
    if (mut_genotypes_connections[i].num_genotype == -1) {
        zero_mutation.num_children_mut++;
    };
 };

};


int ProbeMutProcessor::GetFatherMutIndexForCurNode(
                        int index_cur_probe_genotype,
                        int index_cur_genotype_mut,
                        int initial_place) {
  if (index_cur_genotype_mut>initial_place) {
    return (mutation_vector.size() - 1);
  } else {
    int num_genotype = mut_genotypes_connections.
                         at(index_cur_probe_genotype).num_genotype;
    int num_mut_in_genotype = mut_genotypes_connections.
                         at(index_cur_probe_genotype).num_mutation;
    if (num_genotype!=-1 && num_genotype!=-1) {
      int abs_num_father_mut = probe_genotypes.at(num_genotype)->
                                 sequence.at(num_mut_in_genotype);
      auto it = std::find_if(mutation_vector.begin(), mutation_vector.end(),
                          [=](const MutationNode& mut_node) {
                       return (mut_node.abs_father_mut == abs_num_father_mut);
      });
      return (std::distance(mutation_vector.begin(), it)-1);
    } else {
        return 0;
    };
  };
};
//@function
//@param i_probe_genotype: index of genotype
//@param i_genotype_mut: index of mutation
//@param attachment_location: index of attachment of i_probe_genotype
MutationNode ProbeMutProcessor::GetCurrentMutNode(
                                     int index_cur_probe_genotype,
                                     int index_cur_genotype_mut,
                                     int initial_place) {
  MutationNode current_mutation;
  // determines abs_mut and depth:
  current_mutation.abs_mut = probe_genotypes.at(index_cur_probe_genotype)->
                           sequence.at(index_cur_genotype_mut);
  //calibrates w.r.t. zero mutation
  current_mutation.depth = index_cur_genotype_mut + 1;//
  /////////////////////////////////////////////////////////////////////////////
  // determines abs_father_mut and index_father_mut
  //&&  mut_genotypes_connections.at(i_probe_genotype).num_mutation == -1
  if ( index_cur_genotype_mut == 0  ) {
    current_mutation.abs_father_mut = -1;
    current_mutation.index_father_mut = 0;
  } else {
    current_mutation.abs_father_mut =
    probe_genotypes.at(index_cur_probe_genotype)->
      sequence.at(index_cur_genotype_mut - 1);

    current_mutation.index_father_mut = GetFatherMutIndexForCurNode(
                                         index_cur_probe_genotype,
                                         index_cur_genotype_mut,
                                         initial_place);
  };
  current_mutation.number_cells=GetNumCellsForCurNode(index_cur_genotype_mut,
                                                    current_mutation.abs_mut);
  current_mutation.num_children_mut=GetNumChildrenCurNode(
                                      index_cur_probe_genotype,
                                      index_cur_genotype_mut);
  return current_mutation;
};


// @function, small description:
//  of attachment of genotype
void ProbeMutProcessor::CreatMutVector() {
  MutationNode zero_mutation;
  InitZeroMutationNode(zero_mutation);
  mutation_vector.push_back(zero_mutation);
  for (int i = probe_genotypes.size() - 1; i >= 0; i--) {
    // CAN BE -1!!!!!!!
    int attachment_location = mut_genotypes_connections.at(i).num_mutation;
    //int initial_place = (attachment_location != -1) ? attachment_location : 0;
    unsigned int initial_place = abs(attachment_location + 1);
    for(unsigned int j = initial_place;
        j < probe_genotypes.at(i)->sequence.size();
        j++) {
      // adjoins mutation in a sequence to the main mutation
      MutationNode current_mutation =
      GetCurrentMutNode(i, j, initial_place);
      mutation_vector.push_back(current_mutation);
    };
  };
};

void ProbeMutProcessor::ConstructorInit() {
  SetCellsGenotypes();// SaveGenotypes("BeforeUniqueGenotypes();",probe_genotypes);
  UniqueGenotypes(); //SaveProbeCells("after_unique_genotypes");
  UniqueProbeCells();
  GenotypesWithoutInclusions();
  CreatMutGenotypesConnections();
  CreatMutVector();
};

ProbeMutProcessor::ProbeMutProcessor(const vector <Cell> & _probe_cells,
                                     unsigned int threshold_numb_cells) {

  raw_probe_cells=_probe_cells;
  probe_cells=_probe_cells;
  probe_cells.erase(probe_cells.begin() + threshold_numb_cells);
  // sorts cells by the length of genotype
  stable_sort (probe_cells.begin(), probe_cells.end(), CompareCells());
  ConstructorInit();
};

ProbeMutProcessor::ProbeMutProcessor(const vector <Cell> & _probe_cells) {
    raw_probe_cells=_probe_cells;
    probe_cells=_probe_cells;
    // sorts cells by the length of genotype
    //string file_name("before_sort");
    //SaveProbeCells(file_name.c_str());
    stable_sort (probe_cells.begin(), probe_cells.end(), CompareCells());
    //file_name = "after_sort";
    //SaveProbeCells(file_name.c_str());
    ConstructorInit();
};

ProbeMutProcessor::ProbeMutProcessor (const vector <Cell> &,
                   const vector <Genotype *> & genotype_vector) {

};


ProbeMutProcessor::ProbeMutProcessor() {
  //  cout <<"forgot to put parameters"<<endl;
};

ProbeMutProcessor::~ProbeMutProcessor() {

};

void ProbeMutProcessor::SaveProbeCells(const char *name) {
  char name_file[256];
  sprintf(name_file,"Probes_sort_%s.dat",name);
  FILE *f=fopen(name_file, "w") ;
  if (f==NULL) err("err");
  fprintf(f,"probe_cell_size = %u \n",probe_cells.size());
  for (unsigned int i=0;i<probe_cells.size();i++) {
    Genotype *g=genotypes[probe_cells[i].gen] ;
    fprintf(f,"i=%u  ", i);
    for (unsigned int j=0;j<g->sequence.size();++j) {
//mykola: we cut of additional bits for relevant number
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if ((g->sequence[j] & RESISTANT_PM)==0) {
            if ((g->sequence[j] & DRIVER_PM)==0) {
                fprintf(f,"p %u  ",g->sequence[j]);
            } else {
                fprintf(f,"d %u  ",(g->sequence[j] & (~DRIVER_PM)));
            };
            } else {
                fprintf(f,"r %u  ",(g->sequence[j] & (~RESISTANT_PM)));
            };
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      };
    fprintf(f,"%d",g->number); fprintf(f,"&&&\n");
  };
  fprintf(f,"the end");
}

void ProbeMutProcessor::SaveGenotypes(char *name,
                                       vector<Genotype *> & probe_genotypes_)
{
    char name_file[256];

    sprintf(name_file,"%s",name);
    FILE *f=fopen(name_file, "w") ;
    if (f==NULL) err("err");

    for (unsigned int i=0; i<probe_genotypes_.size(); i++) {

    Genotype *g= probe_genotypes_[i];

    fprintf(f,"Genotype_num = %d  ",i);

    for (unsigned int j=0; j<g->sequence.size(); j++) {
//mykola: we cut of additional bits for relevant number
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if ((g->sequence[j] & RESISTANT_PM)==0) {
            if ((g->sequence[j] & DRIVER_PM)==0) {
                fprintf(f,"p %u  ",g->sequence[j]);
            } else {
                fprintf(f,"d %u  ",(g->sequence[j] & (~DRIVER_PM)));
            };
            } else {
                fprintf(f,"r %u  ",(g->sequence[j] & (~RESISTANT_PM)));
            };
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      }
      fprintf(f,"%d",g->number); fprintf(f,"\n");
  }

}

std::ostream & operator<<(std::ostream & output_out,
                          const ProbeMutProcessor & probe_mut_processor) {
    for(unsigned int i=0; i<probe_mut_processor.mutation_vector.size(); i++)  {
        output_out<< probe_mut_processor.mutation_vector[i].abs_mut
        <<"\t"<<probe_mut_processor.mutation_vector[i].depth
        <<"\t"<<probe_mut_processor.mutation_vector[i].abs_father_mut
        <<"\t"<<probe_mut_processor.mutation_vector[i].index_father_mut
        <<"\t"<<probe_mut_processor.mutation_vector[i].num_children_mut
        <<"\t"<<probe_mut_processor.mutation_vector[i].number_cells<<"\t";
    };
    return output_out;
}

std::istream& operator >> (std::istream& s_in,
                         ProbeMutProcessor & probe_mut_processor) {
probe_mut_processor.mutation_vector.clear();
do { MutationNode A;
     s_in >> A.abs_mut>>A.depth >> A.abs_father_mut>>A.index_father_mut>>
            A.num_children_mut>>A.number_cells;
    if (!s_in.eof()) probe_mut_processor.mutation_vector.push_back(A);
} while (!s_in.eof());
return s_in;
};

unsigned int ProbeMutProcessor::get_tree_mut_depth() const {
  unsigned int max_depth = 0;
  for(auto it = mutation_vector.begin(); it < mutation_vector.end(); ++it)
    max_depth = (max_depth < it->depth) ? it->depth : max_depth;
  return max_depth;
};

vector <unsigned int> ProbeMutProcessor::get_depth_mut_vec() const {
  unsigned int tree_mut_depth = get_tree_mut_depth();
  vector <unsigned int> depth_mut_vector(tree_mut_depth + 1, 0);
  unsigned int mut_number = get_mut_number();
  for (unsigned int i=0; i < mut_number; i++)
    depth_mut_vector[mutation_vector[i].depth]++;
  return depth_mut_vector;
};

float DeltaSquare(unsigned int num_1, unsigned int sum_1,
                  unsigned int num_2, unsigned int sum_2) {
  float frac_1 = (1.0 * num_1)/(1.0 * sum_1);
  float frac_2 = (1.0 * num_2)/(1.0 * sum_2);
  return pow( frac_1 - frac_2, 2.0);
};

float HalfFracSum(unsigned int num_1, unsigned int sum_1,
                  unsigned int num_2, unsigned int sum_2) {
  float half_frac_1 = (1.0 * num_1)/(2.0 * sum_1);
  float half_frac_2 = (1.0 * num_2)/(2.0 * sum_2);
  return (half_frac_1 + half_frac_2);
};
float chi_square_sum(const vector <unsigned int> vec_1,
                     const vector <unsigned int> vec_2) {
  vector <unsigned int> vec1(vec_1);
  vector <unsigned int> vec2(vec_2);
  unsigned int sum1 = std::accumulate(vec1.begin(),vec1.end(), 0);
  unsigned int sum2 = std::accumulate(vec2.begin(),vec2.end(), 0);
  unsigned int max_depth = 0;
 //makes size of vec1 and vec2 equal to max of both
  if (vec1.size() < vec2.size()) {
    max_depth =vec2.size();
    unsigned int diff_langth = vec2.size() - vec1.size();
    vec1.insert(vec1.end(), diff_langth, 0);
  } else {
    max_depth = vec1.size();
    unsigned int diff_langth = vec1.size() - vec2.size();
    vec2.insert(vec2.end(), diff_langth, 0);
  };

  float chi_square = 0;
  for (unsigned int i=0; i < max_depth; i++) {
    float del_sq=DeltaSquare(vec1[i], sum1,
                             vec2[i], sum2);
    if ((vec1[i]!=0) || (vec2[i]!=0)) {
      float den=HalfFracSum(vec1[i], sum1,
                            vec2[i], sum2);
      chi_square += del_sq/den;
    };
  };
  return chi_square;
};

float chi_square_depth(const ProbeMutProcessor & probe_1,
                          const ProbeMutProcessor & probe_2) {
  vector <unsigned int> depth_mut_vec_pr_1 = probe_1.get_depth_mut_vec();
  vector <unsigned int> depth_mut_vec_pr_2 = probe_2.get_depth_mut_vec();
  return chi_square_sum (depth_mut_vec_pr_1,
                         depth_mut_vec_pr_2);
};

unsigned int ProbeMutProcessor::get_tree_mut_degree() const {
  unsigned int max_degree = 0;
  for(auto it = mutation_vector.begin(); it < mutation_vector.end(); ++it)
    max_degree = (max_degree < it->num_children_mut) ?
                                it->num_children_mut : max_degree;
  return max_degree;
};

vector <unsigned int> ProbeMutProcessor::get_degree_mut_vec() const {
  unsigned int tree_mut_degree = get_tree_mut_degree();
  vector <unsigned int> degree_mut_vector(tree_mut_degree + 1, 0);
  unsigned int mut_number = get_mut_number();
  for (unsigned int i=0; i < mut_number; i++)
    degree_mut_vector[mutation_vector[i].num_children_mut]++;
  return degree_mut_vector;
};

float chi_square_degree (const ProbeMutProcessor & probe_1,
												 const ProbeMutProcessor & probe_2) {
  vector <unsigned int> degree_mut_vec_pr_1 = probe_1.get_degree_mut_vec();
  vector <unsigned int> degree_mut_vec_pr_2 = probe_2.get_degree_mut_vec();
  return chi_square_sum (degree_mut_vec_pr_1,
                         degree_mut_vec_pr_2);
}

float ProbeMutProcessor::get_mean_mut_num_cell() const {
  unsigned int sum = 0;
  for(const MutationNode & m_n: mutation_vector){
		if (m_n.number_cells>0) sum += m_n.number_cells*m_n.depth;
  };
  return (1.*sum/get_cell_number());
}


float ProbeMutProcessor::get_mean_pass_b_w_non_elem_nodes() const {
  unsigned int sum_len = 0;
  for(const MutationNode & m_n: mutation_vector){
		if (m_n.number_cells>0) {
		  // lengths b/w the first_index and the first_index of non elementary
		  // mut nodes on the pass to a cell
		  unsigned int first_index = m_n.depth;
		  //bool is_first_index_fund = false;
		  unsigned int last_index = m_n.depth;
		  bool is_last_index_fund = false;
		  MutationNode current_mut_node = m_n;
		  for(unsigned int i=m_n.depth; i > 0; i--) {
			  if (current_mut_node.num_children_mut>1) {
						if (!is_last_index_fund) {
							last_index = current_mut_node.depth;
							is_last_index_fund = true;
			      };
						first_index = current_mut_node.depth;
			  };
	      if (current_mut_node.index_father_mut>-1) current_mut_node=
	        mutation_vector.at(current_mut_node.index_father_mut);
		  }
		  sum_len += (last_index - first_index) * m_n.number_cells;
		}
  }
  return (1.*sum_len/get_cell_number());
}

float dist_mean_pass_b_w_non_elem_nodes(const ProbeMutProcessor & Probe_mut_1,
                                        const ProbeMutProcessor & Probe_mut_2) {
	float mean_sum = (Probe_mut_1.get_mean_pass_b_w_non_elem_nodes() +
	                  Probe_mut_2.get_mean_pass_b_w_non_elem_nodes())/2;
	float abs_diff = abs( Probe_mut_1.get_mean_pass_b_w_non_elem_nodes() -
											  Probe_mut_2.get_mean_pass_b_w_non_elem_nodes() );
	return (abs_diff/mean_sum);
};

float dist_mean_pass_to_cell(const ProbeMutProcessor & Probe_mut_1,
														 const ProbeMutProcessor & Probe_mut_2) {
	float mean_sum = (Probe_mut_1.get_mean_mut_num_cell() +
	                  Probe_mut_2.get_mean_mut_num_cell())/2;
	float abs_diff = abs( Probe_mut_1.get_mean_mut_num_cell() -
											  Probe_mut_2.get_mean_mut_num_cell() );
	return (abs_diff/mean_sum);
};


float t_dist_deg_vs_dep (const ProbeMutProcessor & Probe_mut_1,
                         const ProbeMutProcessor & Probe_mut_2) {
    float chi_square_dep = chi_square_depth(Probe_mut_1, Probe_mut_2);
    float chi_square_deg = chi_square_degree(Probe_mut_1, Probe_mut_2);
    return (chi_square_dep + chi_square_deg);
};
float t_dist_deg_vs_mean_pass_non_elem(const ProbeMutProcessor & Probe_mut_1,
                               const ProbeMutProcessor & Probe_mut_2) {
	float chi_square_deg = chi_square_degree(Probe_mut_1, Probe_mut_2);
	float dist_mean = dist_mean_pass_b_w_non_elem_nodes(Probe_mut_1, Probe_mut_2);
  return (chi_square_deg + dist_mean);
};

float t_dist_deg_vs_mean_pass_to_cell (const ProbeMutProcessor & Probe_mut_1,
                                       const ProbeMutProcessor & Probe_mut_2) {
	float chi_square_deg = chi_square_degree(Probe_mut_1, Probe_mut_2);
	float dist_mean = dist_mean_pass_to_cell(Probe_mut_1, Probe_mut_2);
  return (chi_square_deg + dist_mean);
};

float t_dist_mean_pass_to_cell (const ProbeMutProcessor & Probe_mut_1,
																const ProbeMutProcessor & Probe_mut_2) {
	float dist_mean = dist_mean_pass_to_cell(Probe_mut_1, Probe_mut_2);
  return dist_mean;
};


bool IsCellFittedToPiece(const int cell_index,
                         const int piece_index,
                         const int piece_min_x,
                         const Area & probe_area){
    // x min is shifted for the different probes
    bool is_cell_x_in_range = (cells[cell_index].x >= piece_min_x) &&
                              (cells[cell_index].x <= probe_area.max_x);
    bool is_cell_y_in_range = (cells[cell_index].y >= probe_area.min_y) &&
                              (cells[cell_index].y <= probe_area.max_y);
    bool is_cell_z_in_range = (cells[cell_index].z >= probe_area.min_z) &&
                              (cells[cell_index].z <= probe_area.max_z);
    return (is_cell_x_in_range && is_cell_y_in_range && is_cell_z_in_range);
};

int GetMinxProbePiece(const Parameters & pars, const int index) {
    //makes piece_min_x 0 or bigger
    int del_x = pars.probe_pars.del_area.del_x;
    /*probe_area.max_x - probe_area.min_x;*/
    int del_y = pars.probe_pars.del_area.del_y;
    //  projection of vector with coords (?,i*del_y,0)
    //  shift vector (max_x,0,0) to (?, del_y, 0) and find ? = piece_min_x
    float square_rute = pow(pow(del_x, 2.0) - pow(index*del_y, 2.0), 0.5);  /*10*/
    int piece_min_x = del_x - int (square_rute);
    return piece_min_x;
};


void InitProbePieceVector (vector <ProbePiece> & probe_piece_vector,
                           const Parameters & pars){
  int piece_num = pars.etalon.piece_num;
  for (int i = 0; i < piece_num; i++) {
      ProbePiece local_probe_piece;
      local_probe_piece.indicator = 0;
      local_probe_piece.cell_vector.clear();
      local_probe_piece.piece_min_x_vec = GetMinxProbePiece(pars, i);
      probe_piece_vector.push_back(local_probe_piece);
  };
};

void GetVecOfProbePiece(vector <ProbePiece> & probe_piece_vector,
                        const Parameters & pars) {
  for (unsigned int i=0; i < cells.size(); i++) {  // conditional loop
    for (int j=0; j < pars.etalon.piece_num; j++) {
      if (IsCellFittedToPiece(i, j,
                              probe_piece_vector.at(j).piece_min_x_vec,
                              pars.probe_pars.area)
					&& (probe_piece_vector.at(j).indicator <
          pars.probe_pars.cell_num) ) {
        probe_piece_vector.at(j).cell_vector.push_back(cells.at(i));
        probe_piece_vector.at(j).indicator++;
        break;
      };
    };

    bool is_search_finished = true;
		for (int j=0; j < pars.etalon.piece_num; j++) {
      if (probe_piece_vector.at(j).indicator <
            pars.probe_pars.cell_num) is_search_finished = false;
    };
    if (is_search_finished) break;
  };
};

void GetVecOfRandProbePiece(vector <ProbePiece> & probe_piece_vector,
                            const Parameters & pars) {
  for (int i = 0; i < pars.etalon.piece_num; i++) {
		for (int j = 0; j < pars.probe_pars.cell_num; j++) {
			unsigned int index = trunc((cells.size()-1)*_drand48());
			probe_piece_vector.at(i).cell_vector.push_back(cells.at(index));
		};
  };
};


void SetComparisonResults2File(const Parameters & pars,
                               const vector <ProbePiece> & probe_piece_vector ) {
 //opens file_tree_comparison.dat
  ofstream comparison_file;
  string full_comparison_file_name = pars.catalog_file_name + "/" +
                                     pars.etalon.comparison_file_name;
  comparison_file.open(full_comparison_file_name.c_str(), std::ios::app);

  for (int i=0; i<pars.etalon.piece_num; i++) {
    if (probe_piece_vector.at(i).results.dist < pars.etalon.max_dist_probe_trees) {
      comparison_file.precision(5);
      comparison_file.setf(std::ios::fixed, std:: ios::floatfield);
      comparison_file << pars.evo_pars.driver_adv     <<"\t"<<
        pars.evo_pars.mutation_rate                   <<"\t"<<
        pars.evo_pars.driver_mutation_rate            <<"\t"<<
        probe_piece_vector.at(i).results.dist         <<"\t"<<
        probe_piece_vector.at(i).results.delta_mut    <<"\t"<<
        probe_piece_vector.at(i).results.num_mut            <<endl;
//      cout<<"all info succ_prob:  "<<pars.evo_pars.driver_adv <<
//        " "<<pars.evo_pars.mutation_rate                      <<
//        " "<<pars.evo_pars.driver_mutation_rate               <<
//        " "<<probe_piece_vector.at(i).results.dist            <<
//        " "<<probe_piece_vector.at(i).results.delta_mut       <<
//        " "<<probe_piece_vector.at(i).results.num_mut         <<endl;
    } else {
 //     cout<<"dist("<<i<<")=" << probe_piece_vector.at(i).results.dist <<endl;
    };
  };
  comparison_file.close();
}

ComparisonResults GetComparisonResults( const ProbeMutProcessor & probe,
                                        const ProbeMutProcessor & e_probe,
                                        const Parameters & pars) {
  ComparisonResults c_results;
  cout<<"successful probe"<<endl;
  switch (pars.etalon.dist_type) {
		case 0: {
			c_results.dist = t_dist_deg_vs_dep(probe, e_probe);
			break;
		 }
		case 1: {
			c_results.dist = t_dist_deg_vs_mean_pass_non_elem(probe, e_probe);
      break;
		}
		case 2: {
			c_results.dist = t_dist_deg_vs_mean_pass_to_cell (probe, e_probe);
      break;
		}
		case 3: {
			c_results.dist = t_dist_mean_pass_to_cell (probe, e_probe);
      break;
		}

     default:
     err("Mistake in distance type (see setting file)");
	};
  c_results.e_num_mut = e_probe.get_mut_number();
  c_results.num_mut   = probe.get_mut_number();
  c_results.delta_mut = abs(c_results.e_num_mut - c_results.num_mut);

  cout<<"tree dist is "                << c_results.dist      << endl;
  cout<<"Number of etalon mutations "  << c_results.e_num_mut << endl;
  cout<<"Number of probe mutations "   << c_results.num_mut   << endl;
  cout<<"Difference in mut number is " << c_results.delta_mut << endl;

  return c_results;
};


bool IsProbePieceCompatible2EPorobe(  const ProbeMutProcessor & probe,
                                      const ProbeMutProcessor & e_probe,
                                      const Parameters & pars) {
  int probe_num_cells = probe.get_cell_number();
  int e_probe_num_cells = e_probe.get_cell_number();
  int delta_size_cells = abs(probe_num_cells - e_probe_num_cells);
  float size_level = e_probe_num_cells * pars.etalon.threshold_frac_cell_diff;

  if (delta_size_cells <= size_level) {
    cout<<"successful probe (there are enough cells ) e_probe_cell= "
        <<e_probe.get_cell_number()<<"; probe_cell ="
        << probe.get_cell_number() <<endl;

    int probe_mut_num     =  probe.get_mut_number();
    int e_probe_mut_num   =  e_probe.get_mut_number();
    int delta_probe_e_mut =  abs(e_probe_mut_num - probe_mut_num);

    if (delta_probe_e_mut < e_probe_mut_num *
                            pars.etalon.threshold_frac_mut_diff) {
      return true;
    } else {
      cout<<"There are " << probe_mut_num <<
          " muts and the demand is " << e_probe_mut_num <<" muts"<<endl;
      return false;
    };

  } else {
    cout<<"there are "<< probe_num_cells <<
          " and the demand is " << e_probe_num_cells <<" cells"<<endl;
    return false;
  };
};


void ProcessProbePieceVector(vector <ProbePiece> & probe_piece_vector,
                       Parameters & pars) {
  string full_name_file_mut_vec (pars.etalon.cataloge_with_e_mut_file + "/" +
                                 pars.etalon.e_mut_file_name);
  ifstream file_mut_vec;
  file_mut_vec.open(full_name_file_mut_vec.c_str());
  if (file_mut_vec.is_open()) {
    ProbeMutProcessor e_probe;  // etalon probe
    file_mut_vec >> e_probe;
    for (int i=0; i < pars.etalon.piece_num; i++) {
      ProbeMutProcessor probe(probe_piece_vector.at(i).cell_vector);
			std::cout<<"piece_num = " << i << std::endl;
			std::cout<<"mean mut number of cells in e_probe = " <<
				e_probe.get_mean_mut_num_cell() << std::endl;
			std::cout<<"get_mean_pass_b_w_non_elem_nodes for e_probe = " <<
				e_probe.get_mean_pass_b_w_non_elem_nodes()<<std::endl;
			std::cout<<"mean mut number of cells in piece_probe = " <<
				probe.get_mean_mut_num_cell() << std::endl;
			std::cout<<"get_mean_pass_b_w_non_elem_nodes for piece_probe = " <<
				probe.get_mean_pass_b_w_non_elem_nodes()<<std::endl;
      bool is_compatible =
        IsProbePieceCompatible2EPorobe(probe, e_probe, pars);
      if (is_compatible){
        probe_piece_vector.at(i).results =
          GetComparisonResults(probe, e_probe, pars);
      };
    };
    SetComparisonResults2File(pars, probe_piece_vector);
    file_mut_vec.close();
  } else {
    file_mut_vec.close();
  };
};

void DoesInferenceAnalyze(Parameters & pars) {
  vector <ProbePiece> probe_piece_vector;
  InitProbeArea(pars);
  InitProbePieceVector(probe_piece_vector, pars);
  if (pars.probe_pars.is_random) {
		GetVecOfRandProbePiece(probe_piece_vector, pars);
  } else {
    GetVecOfProbePiece(probe_piece_vector, pars);
  };
  ProcessProbePieceVector(probe_piece_vector, pars);
};
// cells form area @par probe_area and write them in cell_vector
void GetEtalonCellVector(Parameters pars,
                         vector <Cell> & cell_vector) {
  int j = 0;//@j counter for fitted cells
  //for_loop: Find fitted cells (number less then threshold_num_cells)
  for (unsigned int i=0; i < cells.size() &&
                         j < pars.probe_pars.cell_num; ++i) {
    bool is_x_in_range =
      (cells[i].x >= pars.probe_pars.area.min_x) &&
      (cells[i].x <= pars.probe_pars.area.max_x);
    bool is_y_in_range =
      (cells[i].y >= pars.probe_pars.area.min_y) &&
      (cells[i].y <= pars.probe_pars.area.max_y);
    bool is_z_in_range =
      (cells[i].z >= pars.probe_pars.area.min_z) &&
      (cells[i].z <= pars.probe_pars.area.max_z);

    if ( is_x_in_range && is_y_in_range && is_z_in_range) {
            cell_vector.push_back(cells[i]);
            j++;
    };
  };
  //begin delete this
  //cout<<"number of sells in simulation "<<cells.size()<<endl;
  //end delete this
};

void GetRandomEtalonCellVector(Parameters pars ,
                               vector <Cell> & cell_vector) {
	for (int i=0; i < pars.probe_pars.cell_num; i++) {
		int index = trunc((cells.size()-1)*_drand48());
		cell_vector.push_back(cells.at(index));
	};
};

void small_test(ProbeMutProcessor& probe_mut_processor,
                int size_vector_cells) {
    int num_of_cells_e_prob = probe_mut_processor.get_cell_number();
    cout<<"num_of_cells_e_prob="<<num_of_cells_e_prob<<endl;
    cout<<"size_vector_cells="<<size_vector_cells<<endl;
    assert(size_vector_cells == num_of_cells_e_prob);
};

//Inits probe area. Takes cuboid with del_x, del_y, del_z sides in
// x,y,z direction. IN X DIRECTION IT GOES FROM BORDER SIDE MAX_X
void InitProbeArea(Parameters & pars) {
    pars.probe_pars.area.min_x = pars.probe_pars.population_borders.max_x -
                                 pars.probe_pars.del_area.del_x;
    pars.probe_pars.area.max_x = pars.probe_pars.population_borders.max_x;
    pars.probe_pars.area.min_y = - int (pars.probe_pars.del_area.del_y*1.0/2);
    pars.probe_pars.area.max_y = int (pars.probe_pars.del_area.del_y*1.0/2);
    pars.probe_pars.area.min_z = - int (pars.probe_pars.del_area.del_z*1.0/2);
    pars.probe_pars.area.max_z = int (pars.probe_pars.del_area.del_z*1.0/2);
};

void SetProbeResultsToFile(ProbeMutProcessor & probe_mut_processor,
                           Parameters & pars) {
  string full_pass_to_file =
         pars.catalog_file_name + "/" + pars.file_name;
  ofstream etalon_file_mut_vec;
  etalon_file_mut_vec.open (full_pass_to_file.c_str());
  //cout<<file_name_part<<"/etalon_mut_vec.dat is opened"<<endl;
  etalon_file_mut_vec << probe_mut_processor;
  etalon_file_mut_vec.close();
};

void SetEtalonProbe(Parameters pars) {
  InitProbeArea(pars);
  vector <Cell> cell_vector;
  cell_vector.clear();  // vector preparation
  if (pars.probe_pars.is_random) {  // init cell_vector
		GetRandomEtalonCellVector(pars, cell_vector);
  } else {
		GetEtalonCellVector(pars, cell_vector);
  };
  if (cell_vector.size() > 0) {
    ProbeMutProcessor probe_mut_processor(cell_vector);
    //small_test(probe_mut_processor, cell_vector.size());
    SetProbeResultsToFile(probe_mut_processor, pars);
  } else {
    cout<<"there is no cells in the etalon probe"<< endl;
  };
};

void AnalyzeFrec::add_one_genotype(Cell one_cell) {
  Genotype *g=genotypes[one_cell.gen];
  for (unsigned int j=0;j<g->sequence.size();j++) {
//mykola: we cut of additional bits for relevent number
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if ((g->sequence[j] & RESISTANT_PM)==0) {
            if ((g->sequence[j] & DRIVER_PM)==0) {
                vec_num_mutations.at(g->sequence[j])++;
               // cout<<"g->sequence["<<j<<"]=" <<g->sequence[j]<<endl;
            } else {
                    vec_num_mutations.at((g->sequence[j] & (~DRIVER_PM)))++;
                    //cout<<"g->sequence["<<j<<"] & (~DRIVER_PM) ="
                    //<<(g->sequence[j] & (~DRIVER_PM))<<endl;
            };
         } else {
                vec_num_mutations.at((g->sequence[j] & (~RESISTANT_PM)))++;
            };
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 };
}

AnalyzeFrec::AnalyzeFrec(vector <Cell> &A , unsigned int total_number_of_SNV,
                         unsigned int num_of_div) {
  vec_num_mutations.resize(total_number_of_SNV,0);
  //cout<<"vec_num_mutations size="<<vec_num_mutations.size()<<endl;
  for (unsigned int i=0;i<A.size();i++) add_one_genotype(A.at(i));
  frec_in_div.resize(num_of_div,0);
  for (unsigned int i=0; i < total_number_of_SNV; i++) {
    int detector = 0;
    float proportion = (vec_num_mutations.at(i)*1.0)/(A.size()*1.0)*
                        num_of_div;
    for (unsigned int j=0; j < num_of_div; j++ ){
      if ((proportion>=j) && (proportion<j+1)) detector=j;
    };
    frec_in_div.at(detector)++;
  };
};

AnalyzeFrec::~AnalyzeFrec()
{

};








