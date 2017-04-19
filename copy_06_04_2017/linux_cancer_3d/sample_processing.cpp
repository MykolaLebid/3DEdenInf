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
  for (int i=0; i<mutation_vector.size(); i++) {
    cell_num += mutation_vector[i]. number_cells;
  };
  return cell_num;
};

unsigned int ProbeMutProcessor::get_mut_number() const {
  mutation_vector.size();
};

void ProbeMutProcessor::UniqueGenotypes() {
  vector<Genotype *>::iterator it = unique (probe_genotypes.begin(),
                                            probe_genotypes.end(),
                                            CompareGenotypesEq());
  probe_genotypes.resize( std::distance(probe_genotypes.begin(),it) );
};

void ProbeMutProcessor::SetCellsGenotypes() {
    for (int i=0;i<probe_cells.size(); i++) {
            probe_genotypes.push_back(genotypes[probe_cells[i].gen]);
            probe_genotypes[i]->number = 1;
    };
};
void ProbeMutProcessor::GenotypesWithoutInclusions() {
  for(int i = 0; i < (probe_genotypes.size()-1); i++) {
    bool delete_indicator = false;
    for(int j = i + 1; (j < probe_genotypes.size()) &&
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
    fprintf(f,"num_of_attachted_genotype=%d num_genotype=%d num_mutation=%d\n",
            i, mut_genotypes_connections[i].num_genotype,
           mut_genotypes_connections[i].num_mutation );
  };
};

void ProbeMutProcessor::SaveMutVector(char* name) {
  char name_file[256];
  sprintf(name_file,"%mutation_vector_%s.dat",name);
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
  for(int k=0; k < probe_genotypes.size(); k++){
    if (probe_genotypes.at(k)->sequence.size() >
        current_genotype_mut + 1) break;
    if (probe_genotypes.at(k)->sequence.size() ==
        current_genotype_mut + 1) {
        if (probe_genotypes.at(k)->sequence.at(current_genotype_mut) ==
            abs_mut)
        sum+=probe_genotypes.at(k)->number;
    };
  };
  // searches in excluded_probe_genotypes (in inclusion elements)
  // can be optimize with "find" in sorted vector
  for(int k = 0; (k < excluded_probe_genotypes.size()); k++){
    if (excluded_probe_genotypes.at(k)->sequence.size() >
        current_genotype_mut + 1) break;
    if (excluded_probe_genotypes.at(k)->sequence.size() ==
          current_genotype_mut + 1) {
      if (excluded_probe_genotypes.at(k)->sequence.at(current_genotype_mut) ==
          abs_mut)
          sum+=excluded_probe_genotypes.at(k)->number;
      };
  };
  return sum;
};
// calculates number of children of current mutation in current genotype
int ProbeMutProcessor::GetNumChildrenCurNode(int index_cur_probe_genotype,
                                             int index_cur_genotype_mut) {
  // equals at list 1 if in current genotype mut is not the last one
  int sum = (probe_genotypes[index_cur_probe_genotype]->
             sequence.size() > index_cur_genotype_mut + 1) ? 1 : 0;
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
 for (int i = 0; i < mut_genotypes_connections.size(); i++) {
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

//    unsigned int index_sum = 0;
//    int current_genotype =
//        mut_genotypes_connections.at(index_cur_probe_genotype).num_genotype;
//    for(int i = (probe_genotypes.size() - 1);
//        i > current_genotype;
//        i--)
//    index_sum+= (probe_genotypes.at(i)->sequence.size() -
//                 mut_genotypes_connections.at(i).num_mutation);
//
//    index_sum+= mut_genotypes_connections.at(current_genotype).num_mutation -
//                mut_genotypes_connections.at(current_genotype).num_mutation ;
//    return index_sum;
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
  ///////////////////////////////////////////////////////////////////////////////
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
    int initial_place = attachment_location + 1;
    for(int j = initial_place;
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
  SetCellsGenotypes(); // SaveProbeCells("after_write_cells_genotypes")
  UniqueGenotypes(); //SaveProbeCells("after_unique_genotypes");
  UniqueProbeCells(); //SaveProbeCells("after_UniqueProbeCells");
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
    SaveProbeCells("before_sort");
    stable_sort (probe_cells.begin(), probe_cells.end(), CompareCells());
    SaveProbeCells("after_sort");
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

void ProbeMutProcessor::SaveProbeCells(char *name) {
  char name_file[256];
  sprintf(name_file,"Probes_sort_%s.dat",name);
  FILE *f=fopen(name_file, "w") ;
  if (f==NULL) err("err");
  fprintf(f,"probe_cell_size = %u \n",probe_cells.size());
  for (int i=0;i<probe_cells.size();i++) {
    Genotype *g=genotypes[probe_cells[i].gen] ;
    fprintf(f,"i=%u  ", i);
    for (int j=0;j<g->sequence.size();++j) {
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

    for (int i=0;i<probe_genotypes_.size();i++) {

    Genotype *g= probe_genotypes_[i];

    fprintf(f,"Genotype_num = %d  ",i);

    for (int j=0;j<g->sequence.size();j++) {
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
    for(int i=0;i<probe_mut_processor.mutation_vector.size();i++)  {
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
  int max_depth = 0;
  for(auto it = mutation_vector.begin(); it < mutation_vector.end(); ++it)
    max_depth = (max_depth < it->depth) ? it->depth : max_depth;
  return max_depth;
};

vector <unsigned int> ProbeMutProcessor::get_depth_mut_vec() const {
  unsigned int tree_mut_depth = get_tree_mut_depth();
  vector <unsigned int> depth_mut_vector(tree_mut_depth + 1, 0);
  unsigned int mut_number = get_mut_number();
  for (int i=0; i < mut_number; i++)
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

float my_chi_square_depth(const ProbeMutProcessor & probe_1,
                          const ProbeMutProcessor & probe_2) {
  vector <unsigned int> depth_mut_vec_pr_1 = probe_1.get_depth_mut_vec();
  vector <unsigned int> depth_mut_vec_pr_2 = probe_2.get_depth_mut_vec();
  unsigned int max_depth = 0;
//makes size of depth_mut_vec_Pr_1 depth_mut_vec_Pr_2 equal to max of both
  if (depth_mut_vec_pr_1.size() < depth_mut_vec_pr_2.size()) {
    max_depth = depth_mut_vec_pr_2.size();
    unsigned int diff_langth = depth_mut_vec_pr_2.size() -
                                 depth_mut_vec_pr_1.size();
    depth_mut_vec_pr_1.insert(depth_mut_vec_pr_1.end(), diff_langth, 0);
  } else {
    max_depth = depth_mut_vec_pr_1.size();
    unsigned int diff_langth = depth_mut_vec_pr_1.size() -
                                 depth_mut_vec_pr_2.size();
    depth_mut_vec_pr_2.insert(depth_mut_vec_pr_2.end(), diff_langth, 0);
  };

  float chi_square = 0;
  for (int i=0; i < max_depth; i++){
    float del_sq=DeltaSquare(depth_mut_vec_pr_1[i], probe_1.get_mut_number(),
                             depth_mut_vec_pr_2[i], probe_2.get_mut_number());
    if ((depth_mut_vec_pr_1[i]!=0) || (depth_mut_vec_pr_2[i]!=0)) {
      float den=HalfFracSum(depth_mut_vec_pr_1[i], probe_1.get_mut_number(),
                            depth_mut_vec_pr_2[i], probe_2.get_mut_number());
      chi_square += del_sq/den;
    };
  };
  return chi_square;
};

unsigned int ProbeMutProcessor::get_tree_mut_degree() const {
  int max_degree = 0;
  for(auto it = mutation_vector.begin(); it < mutation_vector.end(); ++it)
    max_degree = (max_degree < it->num_children_mut) ?
                                it->num_children_mut : max_degree;
  return max_degree;
};

vector <unsigned int> ProbeMutProcessor::get_degree_mut_vec() const {
  unsigned int tree_mut_degree = get_tree_mut_degree();
  vector <unsigned int> degree_mut_vector(tree_mut_degree + 1, 0);
  unsigned int mut_number = get_mut_number();
  for (int i=0; i < mut_number; i++)
    degree_mut_vector[mutation_vector[i].num_children_mut]++;
  return degree_mut_vector;
};

float my_chi_square_degree (const ProbeMutProcessor & probe_1,
                            const ProbeMutProcessor & probe_2) {
  vector <unsigned int> degree_mut_vec_pr_1 = probe_1.get_degree_mut_vec();
  vector <unsigned int> degree_mut_vec_pr_2 = probe_2.get_degree_mut_vec();

  unsigned int max_degree = 0;
//makes size of depth_mut_vec_Pr_1 depth_mut_vec_Pr_2 equal to max of both
  if (degree_mut_vec_pr_1.size() < degree_mut_vec_pr_2.size()) {
    max_degree = degree_mut_vec_pr_2.size();
    unsigned int diff_langth = degree_mut_vec_pr_2.size() -
                                 degree_mut_vec_pr_1.size();
    degree_mut_vec_pr_1.insert(degree_mut_vec_pr_1.end(), diff_langth, 0);
  } else {
    max_degree = degree_mut_vec_pr_1.size();
    unsigned int diff_langth = degree_mut_vec_pr_1.size() -
                                 degree_mut_vec_pr_2.size();
    degree_mut_vec_pr_2.insert(degree_mut_vec_pr_2.end(), diff_langth, 0);
  };

  float chi_square = 0;
  for (int i=0; i < max_degree; i++){
    float del_sq=DeltaSquare(degree_mut_vec_pr_1[i], probe_1.get_mut_number(),
                             degree_mut_vec_pr_2[i], probe_2.get_mut_number());
    if ((degree_mut_vec_pr_1[i]!=0) || (degree_mut_vec_pr_2[i]!=0)) {
      float den=HalfFracSum(degree_mut_vec_pr_1[i], probe_1.get_mut_number(),
                            degree_mut_vec_pr_2[i], probe_2.get_mut_number());
      chi_square += del_sq/den;
    };
  };
  return chi_square;
}



float mean_pass_to_cell_dist(const ProbeMutProcessor & Probe_mut_1,
                             const ProbeMutProcessor & Probe_mut_2) {

};




float tree_dist (const ProbeMutProcessor & Probe_mut_1,
                 const ProbeMutProcessor & Probe_mut_2) {
    float chi_square_depth = my_chi_square_depth(Probe_mut_1, Probe_mut_2);
    float chi_square_degree = my_chi_square_degree(Probe_mut_1, Probe_mut_2);
    return (chi_square_depth + chi_square_degree);
};

float tree_dist_degree_depth_opt (ProbeMutProcessor & Probe_mut_1,
                                  ProbeMutProcessor & Probe_mut_2) {
    float chi_square_degree = my_chi_square_degree(Probe_mut_1, Probe_mut_2);
}

//void SetPieceMinXVector(const Parameters pars/*const int piece_num,
//                        const Area & probe_area*/,
//                        vector<int> & piece_min_x_vec) {
//  //initiates piece shifts for x direction under y changes
//  for(int i = 0; i < pars.etalon.piece_num; i++) {
//    int piece_min_x = 0;
//    //makes piece_min_x 0 or bigger
//    int del_x = pars.etalon.probe_pars.del_area.del_x;
//    /*probe_area.max_x - probe_area.min_x;*/
//    int del_y = pars.etalon.probe_pars.del_area.del_y;
//    //  projection of vector with coords (?,i*del_y,0)
//    //  shift vector (max_x,0,0) to (?, del_y, 0) and find ? = piece_min_x
//    float square_rute = (pow(del_x, 2.0) - pow(i*del_y/*10*/, 2.0), 0.5);
//    piece_min_x = del_x - int (square_rute);
//
//    piece_min_x_vec.push_back(piece_min_x);
//    //vector <Cell> vector_cells_i;
//    //v_vector_cells.push_back(vector_cells_i);
//  };
////int indicator_vec=0;
////for(int j=0;j<piece_num;j++)
////    cout <<"min_x["<<j<<"]="<<piece_min_x_vec[j]
////         <<"max_x["<<j<<"]="<<pop_range.max_x
////         <<"min_y["<<j<<"]="<<-int(pop_range.del_2/2)
////         <<"max_y["<<j<<"]="<<int(pop_range.del_2/2)
////         <<"min_z["<<j<<"]="<<(-int(pop_range.del_3/2)+int(10*j))
////         <<"max_z["<<j<<"]="<<(int(pop_range.del_3/2)+int(10*j))
////         <<endl;
//};

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

//void InitVectorsOfCellVectors(vector <vector <Cell>>& v_vector_cells,
//                              const int piece_num) {
//  for (int i = 0; i < piece_num; i++) {
//      vector <Cell> local_vector_cells;
//      v_vector_cells.push_back(local_vector_cells);
//  };
//};

int GetMinxProbePiece(const Parameters & pars, const int index) {
    int piece_min_x = 0;
    //makes piece_min_x 0 or bigger
    int del_x = pars.probe_pars.del_area.del_x;
    /*probe_area.max_x - probe_area.min_x;*/
    int del_y = pars.probe_pars.del_area.del_y;
    //  projection of vector with coords (?,i*del_y,0)
    //  shift vector (max_x,0,0) to (?, del_y, 0) and find ? = piece_min_x
    float square_rute = (pow(del_x, 2.0) - pow(index*del_y/*10*/, 2.0), 0.5);
    piece_min_x = del_x - int (square_rute);
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

//void GetVecOfProbeVectors(const int piece_num, const int probe_cell_num,
//                      Area & probe_area, vector<int> &indicator_vec,
//                      vector <vector <Cell>> &v_vector_cells,
//                      vector<int> &piece_min_x_vec) {
//  for (int i=0; i < cells.size(); i++) {  //conditional loop
//    for (int j=0; j < piece_num; j++) {
//      if (indicator_vec[j] >= probe_cell_num) break;
//      //if (indicator_vec>= probe_cell_num) break;
//      if (IsCellFittedToPiece(i,j, piece_min_x_vec.at(j), probe_area)) {
//        v_vector_cells.at(j).push_back(cells.at(i));
//        indicator_vec[j]++;
//          // indicator_vec++;
//      };
//
//    };
//  };
//};

void GetVecOfProbePiece(vector <ProbePiece> & probe_piece_vector,
                        const Parameters & pars) {
  for (int i=0; i < cells.size(); i++) {  // conditional loop
    for (int j=0; j < pars.etalon.piece_num; j++) {

      if (probe_piece_vector.at(j).indicator >=
          pars.probe_pars.cell_num) break;

      if (IsCellFittedToPiece(i, j,
                              probe_piece_vector.at(j).piece_min_x_vec,
                              pars.probe_pars.area)) {
        probe_piece_vector.at(j).cell_vector.push_back(cells.at(i));
        probe_piece_vector.at(j).indicator++;
      };
    };

    bool is_search_finished = true;
    for (unsigned int j=0;j < pars.etalon.piece_num; j++) {
      if (probe_piece_vector.at(j).indicator <
            pars.probe_pars.cell_num) is_search_finished = false;
    };

    if (is_search_finished) break;
  };
};

//bool IsProbeAppropriate(Parameters pars,
//                        ProbeMutProcessor & e_probe,
//                        vector <Cell> & vector_probe_cells) {
//int probe_num_cells = vector_probe_cells.size();
//int e_probe_num_cells = e_probe.get_cell_number();
//int delta_size_cells = abs(probe_num_cells - e_probe_num_cells);
//
//
//if (delta_size_cells<=e_probe_num_cells*pars.etalon.threshold_frac_cell_diff) {
//    cout<<"successful probe (there are enough cells )"<<endl;
//    ProbeMutProcessor probe(vector_probe_cells);
//    int probe_mut_num = probe.get_mut_number();
//    int e_probe_mut_num = e_probe.get_mut_number();
//    int delta_probe_e_mut = abs(e_probe_mut_num - probe_mut_num);
//    if (delta_probe_e_mut < e_probe_mut_num *
//                            pars.etalon.threshold_frac_mut_diff) {
//    };
//
//} else {
//    cout<<"there is only "<<vector_probe_cells.size()<<
//          " and the demand is " << e_probe_num_cells <<" cells"<<endl;
//};
//
//
//};
//void CheckVecOfProbeVectors(vector <vector <Cell>> & v_vector_cells,
//                            const Parameters & pars,
//                            ProbeMutProcessor & e_probe) {
//    for (int i=0; i< pars.etalon.piece_num; i++){
//     if (IsProbeAppropriate(pars, e_probe , v_vector_cells.at(i))) {
//        int size_v_vec_cells = v_vector_cells.at(i).size();
//        int delta_size_cells = abs(size_v_vec_cells - e_prob_cell_num );
//        if (delta_size_cells <= e_prob_cell_num * level_num_cells){
//            cout<<"successful probe"<<endl;
//            ProbeMutProcessor prob_mut(v_vector_cells.at(i));
//            int num_prob_mut = prob_mut.mutation_vector.size();
//            int delta_prob_etalon_mut = abs( e_prob_mut_num - num_prob_mut);
//            if (delta_prob_etalon_mut<level_num_muts* e_prob_mut_num){
//
//                v_delta_mut.push_back(delta_prob_etalon_mut);
//                float dist = tree_dist(prob_mut,e_prob);
//                v_dist.push_back(dist);
//
//                v_n_muts.push_back(num_prob_mut);
//
//                if ( delta_size_cells < min_mut_diff) {
//                    min_mut_diff = delta_size_cells;
//                    index_min_mut_diff = i;
//                };
//            };
//
//        } else {cout<<"the is only "<<v_vector_cells.at(i).size()<<
//                " cells in a probe with number "<<i<<
//                " and the demand is " << e_prob_cell_num <<" cells"<<endl;
//        };
//
//      };
//    for (int i=0; i<piece_num; i++){
//        int size_v_vec_cells = v_vector_cells.at(i).size();
//        int delta_size_cells = abs(size_v_vec_cells - e_prob_cell_num );
//        float level_num_cells=0.1;
//        float level_num_muts=0.2;
//        if (delta_size_cells <= e_prob_cell_num * level_num_cells){
//            cout<<"successful probe"<<endl;
//            ProbeMutProcessor prob_mut(v_vector_cells.at(i));
//            int num_prob_mut = prob_mut.mutation_vector.size();
//            int delta_prob_etalon_mut = abs( e_prob_mut_num - num_prob_mut);
//            if (delta_prob_etalon_mut<level_num_muts* e_prob_mut_num){
//
//                v_delta_mut.push_back(delta_prob_etalon_mut);
//                float dist = tree_dist(prob_mut,e_prob);
//                v_dist.push_back(dist);
//
//                v_n_muts.push_back(num_prob_mut);
//
//                if ( delta_size_cells < min_mut_diff) {
//                    min_mut_diff = delta_size_cells;
//                    index_min_mut_diff = i;
//                };
//            };
//
//        } else {cout<<"the is only "<<v_vector_cells.at(i).size()<<
//                " cells in a probe with number "<<i<<
//                " and the demand is " << e_prob_cell_num <<" cells"<<endl;
//        };
//
//    };

//};

void SetComparisonResults2File(const Parameters & pars,
                               const vector <ProbePiece> & probe_piece_vector ) {
 //opens file_tree_comparison.dat
  ofstream comparison_file;
  string full_comparison_file_name = pars.catalog_file_name + "/" +
                                     pars.etalon.comparison_file_name;
  comparison_file.open(full_comparison_file_name.c_str());

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
      cout<<"all info succ_prob:  "<<pars.evo_pars.driver_adv <<
        " "<<pars.evo_pars.mutation_rate                      <<
        " "<<pars.evo_pars.driver_mutation_rate               <<
        " "<<probe_piece_vector.at(i).results.dist            <<
        " "<<probe_piece_vector.at(i).results.delta_mut       <<
        " "<<probe_piece_vector.at(i).results.num_mut         <<endl;
    } else {
      cout<<"dist("<<i<<")=" << probe_piece_vector.at(i).results.dist <<endl;
    };
  };
}

ComparisonResults GetComparisonResults( const ProbeMutProcessor & probe,
                                        const ProbeMutProcessor & e_probe,
                                        const Parameters & pars) {
  ComparisonResults c_results;

  cout<<"successful probe"<<endl;

  c_results.dist      = tree_dist(probe, e_probe);
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
    cout<<"successful probe (there are enough cells )"<<endl;

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
    cout<<"there are only "<< probe_num_cells <<
          " and the demand is " << e_probe_num_cells <<" cells"<<endl;
    return false;
  };
};


void DistComRes(vector <ProbePiece> & probe_piece_vector,
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
      bool is_compatible =
        IsProbePieceCompatible2EPorobe(probe, e_probe, pars);
      if (is_compatible){
        probe_piece_vector.at(i).results =
          GetComparisonResults(probe, e_probe, pars);
      };
    };
    SetComparisonResults2File(pars , probe_piece_vector);
    file_mut_vec.close();
  } else {
    file_mut_vec.close();
  };
};

//void CompareWithEtalonProbe( Area & probe_area,
//                             ProbeMutProcessor & e_prob,
//                             ofstream & file_tree_comparison,
//                             Parameters parameters) {
//void CompareWithEtalonProbe( Parameters & pars) {
  //int e_prob_cell_num = e_prob.get_cell_number();
  //int e_prob_mut_num = e_prob.get_mut_number();
  //int piece_num = pars.etalon.piece_num;  // divide circle in @piece_num peaces


  //
  //vector <vector <Cell>> v_vector_cells;
  //InitVectorsOfCellVectors(v_vector_cells, piece_num);

  //vector <int> indicator_vec(piece_num, 0);
  //vector <int> piece_min_x_vec;

//  SetPieceMinXVector(pars /*piece_num, probe_area*/, piece_min_x_vec);

//  GetVecOfProbeVectors( /*piece_num, e_prob_cell_num,
//                        probe_area,*/ indicator_vec,
//                        v_vector_cells, piece_min_x_vec);

//
//    int min_mut_diff = e_prob_mut_num;
//    int index_min_mut_diff= - 1;
//
//    vector<ComparisonResults> v_comparison_results;
//
//    for (int i=0; i<piece_num; i++){
//        int size_v_vec_cells = v_vector_cells.at(i).size();
//        int delta_size_cells = abs(size_v_vec_cells - e_prob_cell_num );
//        float level_num_cells=0.1;
//        float level_num_muts=0.2;
//        if (delta_size_cells <= e_prob_cell_num * level_num_cells){
//            ComparisonResults comparison_results;
//            cout<<"successful probe"<<endl;
//            ProbeMutProcessor prob_mut(v_vector_cells.at(i));
//            int num_prob_mut = prob_mut.mutation_vector.size();
//            int delta_prob_etalon_mut = abs( e_prob_mut_num - num_prob_mut);
//            if (delta_prob_etalon_mut<level_num_muts* e_prob_mut_num){
//                comparison_results.delta_mut = delta_prob_etalon_mut;
//                comparison_results.dist = tree_dist(prob_mut,e_prob);
//                comparison_results.num_mut = num_prob_mut;
//                v_comparison_results.push_back(comparison_results);
//
//                if ( delta_size_cells < min_mut_diff) {
//                    min_mut_diff = delta_size_cells;
//                    index_min_mut_diff = i;
//                };
//            };
//
//        } else {cout<<"the is only "<<v_vector_cells.at(i).size()<<
//                " cells in a probe with number "<<i<<
//                " and the demand is " << e_prob_cell_num <<" cells"<<endl;
//        };
//
//    };
//
//    int size_of_succ_prob = v_delta_mut.size();
//
//    //cout << "size_of_succ_prob="<< size_of_succ_prob<<endl;
//    for(int i=0; i<size_of_succ_prob; i++)
//    {
//        if (v_dist.at(i)<max_dist_for_prob){
//            file_tree_comparison.precision(5);
//            file_tree_comparison.setf(std::ios::fixed, std:: ios::floatfield);
//            file_tree_comparison
//            <<driver_adv<<"\t"<<gama<<"\t"<<driver_prob<<"\t"<<v_dist.at(i)
//            <<"\t"<<v_delta_mut.at(i)
//            <<"\t"<<v_n_muts.at(i)//e_prob.mutation_vector.size()
//            <<endl;
//           // cout<<"size_of_succ_prob  "<<driver_adv<<" "<<gama<<" "<<driver_prob<<" "<<v_dist.at(i)
//           // <<" "<<v_delta_mut.at(i)
//           // <<" "<<v_n_muts.at(i)//e_prob.mutation_vector.size()
//           // <<endl;
//        } else {
//            cout<<"v_dist.at("<<i<<")=" <<v_dist.at(i)<<endl;
//        };
//
//    };
//
//
//
//        //int min_diff_in_mut =  e_prob_mut_num;
//        //int min_diff_in_cells = e_prob_cell_num;
//        //for (int j=0; j<piece_num;j++){
//        //    if v_vector_cells.at(j).size()
//        //}

//}


//};


//float TakeProbe(char *name,
//                 Area & pop_range,
//                 ofstream & file_tree_comparison,
//                 Parameters Pars) {
void DoesInferenceAnalyze(Parameters & pars) {
  vector <ProbePiece> probe_piece_vector;
  InitProbeArea(pars);
  InitProbePieceVector (probe_piece_vector, pars);
  GetVecOfProbePiece(probe_piece_vector, pars);
  DistComRes(probe_piece_vector, pars);

  int i;
  cin>>i;
//          //begin * Takes info from file with etalon results
//          ofstream file_with_etalon_probe;
//          char name_file[256];
//          sprintf(name_file,"%s/file_tree_comparison.dat",
//                  pars.catalog_file_name);
//          file_tree_comparison.open(name_file, ios::ate | ios::in);
//          if (file_tree_comparison.is_open()==false)
//            file_tree_comparison.open(name_file);
//          //end *
//          file_tree_comparison.close();

  //begin * Open file with etalon data
  //char full_name_file_mut_vec[256];


    //int num_of_cells_e_prob=0;
    //for (int i=0; i<e_prob.mutation_vector.size();i++)
    //num_of_cells_e_prob+=e_prob.mutation_vector[i].number_cells;
    //cout<<"download num_of_cells_e_prob="<<num_of_cells_e_prob<<endl;
//    CompareWithEtalonProbe(pop_range, e_prob,
//                           file_tree_comparison,
//                           Pars);


        //Graph !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //e_prob.SaveGraph(name);
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


/*
        if (vector_cells.size()>0){

            ProbeMutProcessor A(vector_cells);
            float dist = tree_dist(A,e_prob);

            int num_of_cells_A=0; for (int i=0; i<A.mutation_vector.size();i++) num_of_cells_A+=A.mutation_vector[i].number_cells;



            file_tree_comparison
            <<"driver_adv="<<driver_adv<<"; "<<"dist="<<dist<<";N_c_A=" << num_of_cells_A
            <<";N_c_e_prob="<<num_of_cells_e_prob
            <<";A_m_v_s="<<A.mutation_vector.size()
            <<"e_prob_m_v_s="<<e_prob.mutation_vector.size()
            <<endl;



        } else{cout<<"vector_cells is empty"<<endl;
               file_tree_comparison<<"vector_cells is empty"<<endl;

               //cout<<"!!!!!!!!!!!!!!!!!!!!!"<<endl;
               //file_tree_comparison<<"Pr.x_max="<<Pr.x_max<<"; Pr.y_max="<<Pr.y_max<<"; Pr.z_max="<<Pr.z_max
               //<<"Pr.x_min="<<Pr.x_min<<"; Pr.y_min="<<Pr.y_min<<"; Pr.z_min="<<Pr.z_min<<endl;
               //file_tree_comparison<<"j="<<j<<endl;
            };
*/
};
// cells form area @par probe_area and write them in cell_vector
void GetEtalonCellVector(Parameters pars,
                         vector <Cell> & cell_vector) {
  int j = 0;//@j counter for fitted cells
  //for_loop: Find fitted cells (number less then threshold_num_cells)
  for (int i=0; i < cells.size() && j < pars.probe_pars.cell_num; ++i) {
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
  GetEtalonCellVector(pars, cell_vector);  // init cell_vector
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
  for (int j=0;j<g->sequence.size();j++) {
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
  for (int i=0;i<A.size();i++) add_one_genotype(A.at(i));
  frec_in_div.resize(num_of_div,0);
  for (int i=0; i < total_number_of_SNV; i++) {
    int detector = 0;
    float proportion = (vec_num_mutations.at(i)*1.0)/(A.size()*1.0)*
                        num_of_div;
    for (int j=0; j < num_of_div; j++ ){
      if ((proportion>=j) && (proportion<j+1)) detector=j;
    };
    frec_in_div.at(detector)++;
  };

  // for (int i=0; (i < num_of_div) && (i<200); i++)
  // cout<<"num_of_div["<<i<<"]="<<frec_in_div.at(i)<<endl;
};

AnalyzeFrec::~AnalyzeFrec()
{

};








