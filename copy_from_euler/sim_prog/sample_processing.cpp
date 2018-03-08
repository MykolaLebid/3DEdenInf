// Copyright 2017 Lebid Mykola all rights reserved


#include "sample_processing.h"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>


#include <vector>
#include <algorithm>
#include <iterator>
#include <memory>
#include <numeric>

#include "classes.h"

#include <assert.h>
#include "treedist.h"

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
  if (f==NULL) err("err SaveMutGenotypeConnections ");
  int size_mgc = mut_genotypes_connections.size();
  for(int i = 0 ;i < size_mgc; i++) {
    fprintf(f,"num_of_attachted_genotype=%d num_genotype=%d num_mutation=%d\n",
            i, mut_genotypes_connections[i].num_genotype,
           mut_genotypes_connections[i].num_mutation );
    int number_of_children = mut_genotypes_connections[i].
																vec_index_attached_genotypes.size();
		fprintf(f,"children with numbers:");

    for(int j=0;j < number_of_children;j ++ )
			fprintf(f,"attached genotype %d - place of attachment -  %d , \n",
							mut_genotypes_connections[i].vec_index_attached_genotypes[j].num_genotype,
							mut_genotypes_connections[i].vec_index_attached_genotypes[j].num_mutation);
  };
};

void ProbeMutProcessor::SaveMutVector(char* name) {
  char name_file[256];
  sprintf(name_file,"mutation_vector_%s.dat",name);
  FILE *f=fopen(name_file, "w") ;
  if (f==NULL) err("err SaveMutVector");

  int size_mv = mutation_vector.size();
  for(int i = 0 ;i < size_mv; i++) {
    fprintf(f,"%d  abs_mut=%d depth=%d abs_father_mut=%d index_father_mut=%d ",
            i,
            mutation_vector[i].abs_mut, mutation_vector[i].depth,
            mutation_vector[i].abs_father_mut,
            mutation_vector[i].index_father_mut);
    fprintf(f,"num_children_mut=%d number_cells=%d",
            mutation_vector[i].num_children_mut,
            mutation_vector[i].number_cells);
    for(unsigned int j=0; j < mutation_vector[i].vec_child_places.size(); j++) {
			fprintf(f,"children: %d", mutation_vector[i].vec_child_places[j]);
    };
		fprintf(f,"\n");

  };
};

void ProbeMutProcessor::SaveMutHubVector (char *name) {

  char name_file[256];
  sprintf(name_file,"mutation_hub_vector_%s.dat",name);
  FILE *f = fopen(name_file, "w");
  if (f==NULL) err("err SaveMutHubVector");

  int size_mv = mut_hub_vector.size();
  for(int i = 0 ;i < size_mv; i++) {
    fprintf(f,"%d  abs_mut_hub=%d depth=%d node_number=%d abs_father_hub=%d ", i,
            mut_hub_vector[i].abs_mut, mut_hub_vector[i].depth,
            mut_hub_vector[i].node_number, mut_hub_vector[i].abs_father_hub);
		fprintf(f,"index_father_hub=%d ",mut_hub_vector[i].index_father_hub);
    fprintf(f,"num_children_hubs=%d ", mut_hub_vector[i].num_children_hubs);
    fprintf(f," number_cells = %d ", mut_hub_vector[i].number_cells);

    //for(unsigned int j=0; j < mut_hub_vector[i].vec_child_places.size(); j++) {
		//	fprintf(f,"children: %d", mut_hub_vector[i].vec_child_places[j]);
    //};
		fprintf(f,"\n");

  };


};

string GetStringFromAbsMutNum(const int abs_mut_num){
	char buffer [50];
	int n = 0 ;

  if ((abs_mut_num & RESISTANT_PM)==0) {
			if ((abs_mut_num & DRIVER_PM)==0) {
				  n=sprintf (buffer, "p%d", abs_mut_num);
          string str(buffer);
          str = str;
          return str;
			} else {
					n=sprintf (buffer, "d%d", (abs_mut_num & (~DRIVER_PM)));
          string str(buffer);
          str = str;
          return str;};
	 } else {
					n=sprintf (buffer, "r%d * %d", (abs_mut_num & (~RESISTANT_PM)),abs_mut_num);
          string str(buffer);
          str = str;
          return str;};
}

string ProbeMutProcessor::FormCellString(int mut_vec_index)  const {
	int num_cells = mutation_vector[mut_vec_index].number_cells;

	string str/*=")"*/;
	if( num_cells > 0 ) {
		int current_cell_index = 0;
		for (int i = 0; i < mut_vec_index; i++)
			current_cell_index += mutation_vector[i].number_cells;

		for (int i = 0; i < num_cells; i ++ ) {
			char buffer [30];
			int n = sprintf (buffer, "%d", current_cell_index + i);
			string str1(buffer);
			str = str + str1 + ":1.0";
			if (i != (num_cells-1))
				str = str + ",";
			//std::cout<<"local str="<<str<<std::endl;
		};
	};
	//current_cell_index = current_cell_index + num_cells;
	return str;
};
string ProbeMutProcessor::GetCurrentMutationString(int mut_vec_index) const {
	/*unsigned int abs_mut = mutation_vector[mut_vec_index].abs_mut;*/
	string str ="):1.0" /*= GetStringFromAbsMutNum(abs_mut)*/;

	//if( mutation_vector[mut_vec_index].num_children_mut > 0 ) {
	str = FormCellString(mut_vec_index) + str;
	int num_cells = mutation_vector[mut_vec_index].number_cells;
  int mut_num = mutation_vector[mut_vec_index].num_children_mut;
  if (num_cells > 0 && mut_num > 0) str = "," + str;
	for (unsigned int i = 0; i <mut_num; i++) {
		unsigned int child_index = mutation_vector[mut_vec_index].vec_child_places[i];
		str = GetCurrentMutationString(child_index) + str;
		if (i != (mut_num-1)) {
			str = "," + str;
		};
	};
		str="("+str;
	//} else {
	//	str = FormCellString(mut_vec_index);
	//	str="("+str;
	//};
	return str;
};

string ProbeMutProcessor::get_str_newick_mut_tree_format() const {
	string str=";"/*="0d;"*/;
	str = "):1.0" + str;

  str = FormCellString(0) + str;
  int num_cells = mutation_vector[0].number_cells;
  int mut_num = mutation_vector[0].num_children_mut;
  if (num_cells > 0 && mut_num > 0) str = "," + str;
	for(unsigned int i = 0 ; i < mut_num; i++) {
		unsigned int child_index = mutation_vector[0].vec_child_places[i];
		str = GetCurrentMutationString(child_index) + str;
		if (i != (mut_num-1)) {
				str = "," + str;
		};
	};
	str = "(" + str;
	//cout<<str<<endl;
	//cout<<str.length()<<endl;

  return str;
};

void ProbeMutProcessor::SaveNewickTreeFormat(char *name) {
  char name_file[256];
  sprintf(name_file,"newick_%s.dat",name);
  FILE *f=fopen(name_file, "w") ;
  if (f==NULL) err("err err with newick");
	//vector <bool> is_written(size_mv,false);
	string str=";"/*"0d;"*/;
	str = ")" + str;

  str = FormCellString(0) + str;
	for(unsigned int i = 0 ; i < mutation_vector[0].num_children_mut; i++) {
		unsigned int child_index = mutation_vector[0].vec_child_places[i];
		str = GetCurrentMutationString(child_index) + str;
		if ((i != (mutation_vector[0].num_children_mut-1)) ||
					(mutation_vector[0].number_cells != 0) ) {
				str = "," + str;
		};
	};
	str = "(" + str;
	const char * c = str.c_str();
	fprintf(f,c);
	fclose(f);
};

void ProbeMutProcessor::SaveGraph(const char * catalog_name,
																	const char * file_name) {
    char name_file[256];
    sprintf(name_file,"%s/%sGraph.info",catalog_name, file_name);
    FILE *f=fopen(name_file, "w") ;
    if (f==NULL) err("err");
    int size_mv = mutation_vector.size();
    fprintf(f,"digraph graphname {\n");
    for(int i = 0 ;i < size_mv; i++) {
      int sun = mutation_vector[i].abs_mut;
      int father = mutation_vector[i].abs_father_mut;
      if (sun!=-1){
				if (father != -1){
						fprintf(f,"%s",GetStringFromAbsMutNum(father).c_str());
				} else fprintf(f,"init_mut");
			  if (sun != -1) {
						fprintf(f," -> %s",GetStringFromAbsMutNum(sun).c_str());
			  } else fprintf(f," -> init_mut");

        fprintf(f,"\n");
      };
			if (mutation_vector[i].number_cells > 0){
						fprintf(f,"%s -> in_%s_%d_cells",
													GetStringFromAbsMutNum(sun).c_str(),
													GetStringFromAbsMutNum(sun).c_str(),
											    mutation_vector[i].number_cells);
			fprintf(f,"\n");
			};

    };
    fprintf(f,"}\n");
}


void ProbeMutProcessor::SaveHubGraph(const char * catalog_name,
																		 const char * file_name) {
    char name_file[256];
    sprintf(name_file,"%s/%sHubGraph.info",catalog_name, file_name);
    FILE *f=fopen(name_file, "w") ;
    if (f==NULL) err("err");
    int size_mv = mut_hub_vector.size();
    fprintf(f,"digraph graphname {\n");
    for(int i = 0 ;i < size_mv; i++) {
			int sun = mut_hub_vector[i].abs_mut;
			int father = mut_hub_vector[i].abs_father_hub;
      if (sun!=-1) {
				if (father != -1){
						fprintf(f,"%s_%d",GetStringFromAbsMutNum(father).c_str(),
										mut_hub_vector[mut_hub_vector[i].index_father_hub].node_number);
				} else fprintf(f,"init_mut");
			  if (sun != -1) {
						fprintf(f," -> %s_%d",GetStringFromAbsMutNum(sun).c_str(),
										mut_hub_vector[i].node_number);
			  } else fprintf(f," -> init_mut");

        fprintf(f,"\n");
      };
			if (mut_hub_vector[i].number_cells > 0){
						fprintf(f,"%s_%d -> in_%s_%d_cells",
													GetStringFromAbsMutNum(sun).c_str(),
													mut_hub_vector[i].node_number,
													GetStringFromAbsMutNum(sun).c_str(),
											    mut_hub_vector[i].number_cells);
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
    if ( is_genotype_consistent && is_mut_in_genotype_consistent) {
			 Connect A;
			 A.num_genotype = i;
		 	 A.num_mutation = index_cur_genotype_mut;
			 mut_genotypes_connections[index_cur_probe_genotype].
			 vec_index_attached_genotypes.push_back(A);
			 sum++;
		};
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
 // Calculate number of children of zero mutation
 // and save in vec_child_places
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
  current_mutation.number_cells = GetNumCellsForCurNode(index_cur_genotype_mut,
                                                    current_mutation.abs_mut);
  current_mutation.num_children_mut = GetNumChildrenCurNode(
                                      index_cur_probe_genotype,
                                      index_cur_genotype_mut);
  return current_mutation;
};

void ProbeMutProcessor::FillVecChildPlaces() {
	for (unsigned int i = 0 ; i < mutation_vector.size(); i++) {
		for (unsigned int j = i + 1 ; j < mutation_vector.size(); j++) {
			if (mutation_vector[i].abs_mut == mutation_vector[j].abs_father_mut)
				mutation_vector[i].vec_child_places.push_back(j);
		};
	};
}

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
  FillVecChildPlaces();
};

void ProbeMutProcessor::ChangeDepth(
												vector <MutationNode> & local_mut_vec,
												unsigned int index) {

	local_mut_vec.at(index).depth--;

	for (unsigned int i = 0; i < local_mut_vec.at(index).num_children_mut; i++) {
		unsigned int next_index = local_mut_vec[index].vec_child_places[i];
		ChangeDepth(local_mut_vec, next_index);
	};
};

void ProbeMutProcessor::SetIndexFatherHubVec() {
	//current_mut_hub.abs_father_hub = mut_vec.at(i).abs_father_mut;
	//mut_hub_vector[0].index_father_hub = - 1;
	for (unsigned int i = 0; i < mut_hub_vector.size(); i++) {
		for (unsigned int j = 0; j < mut_hub_vector.size(); j++)
			if (mut_hub_vector[i].abs_father_hub == mut_hub_vector[j].abs_mut)
				mut_hub_vector[i].index_father_hub = j;
	};
};

void ProbeMutProcessor:: CreatMutHubVector() {
	vector <MutationNode> mut_vec(mutation_vector);
  MutationHub current_mut_hub;
  current_mut_hub.node_number = 0;
  current_mut_hub.abs_father_hub = -1;
  current_mut_hub.index_father_hub = 0;
 	bool save_lock = false;

  for (unsigned int i = 0; i < mut_vec.size(); i ++) {

			if (save_lock == false) {
				  current_mut_hub.abs_father_hub = mut_vec.at(i).abs_father_mut;
					save_lock = true;
			};
			if (mut_vec.at(i).number_cells == 0 &&
					mut_vec.at(i).num_children_mut == 1) {
				ChangeDepth(mut_vec, i);
				current_mut_hub.node_number++;
			} else {
				current_mut_hub.abs_mut = mut_vec.at(i).abs_mut;
				current_mut_hub.depth = mut_vec.at(i).depth;
				current_mut_hub.number_cells = mut_vec.at(i).number_cells;
				current_mut_hub.num_children_hubs = mut_vec.at(i).num_children_mut;
				current_mut_hub.node_number ++;

				mut_hub_vector.push_back(current_mut_hub);
				current_mut_hub.node_number = 0;

				save_lock = false;
			};
  };
  SetIndexFatherHubVec();
  AdjustLabels2HubVecUsingMut();
};

void ProbeMutProcessor::AdjustLabels2HubVecUsingMut() {
	vector <tuple<unsigned int, unsigned int, unsigned int>> hub_vec_number_mut;

	for(unsigned int i=0; i < mut_hub_vector.size(); i++) {
	  unsigned int gereral_mut_num = 0;
	  unsigned int j = i;
		while (mut_hub_vector[j].abs_father_hub!=-1) {
			gereral_mut_num += mut_hub_vector[j].node_number;
			j = mut_hub_vector[j].index_father_hub;
		};

    if ((mut_hub_vector[j].abs_father_hub == -1) &&
				(mut_hub_vector[j].index_father_hub==0)){
						gereral_mut_num += mut_hub_vector[j].node_number;
				};

    tuple<unsigned int ,unsigned int, unsigned int> A;
    A = make_tuple(i , mut_hub_vector[i].abs_mut, gereral_mut_num);

		hub_vec_number_mut.push_back(A);
		//hub_vec_number_mut[i] = gereral_mut_num;
	};

  std::sort(hub_vec_number_mut.begin(), hub_vec_number_mut.end(),
						[](const std::tuple<unsigned int ,unsigned int, unsigned int> &left,
							 const std::tuple<unsigned int ,unsigned int, unsigned int> &right) {
    return get<2>(left) < get<2>(right);
	});


  int label_index = 0;
	for(unsigned int i=0; i < hub_vec_number_mut.size(); i++ ) {
		//cout<<i<<" index in mut_hub_vector = "<< get<0>(hub_vec_number_mut[i])<<"; absolute mut = "<<
		//get<1>(hub_vec_number_mut[i]) <<"; gereral_hub_num = "<< get<2>(hub_vec_number_mut[i])<<endl;
		for(unsigned int j=0; j < mut_hub_vector[get<0>(hub_vec_number_mut[i])].number_cells; j++) {
			mut_hub_vector[get<0>(hub_vec_number_mut[i])].label_vec.push_back(label_index);
		  label_index++;
		};
	};


//	for(unsigned int i=0; i<mut_hub_vector.size(); i++) {
//		cout<<"for element with index "<< i << " in mut_hub_vector ";
//		for(unsigned int j=0; j<mut_hub_vector[i].label_vec.size(); j++)
//			cout<<" we have labels " << mut_hub_vector[i].label_vec[j];
//	  cout<<endl;
//	};
};


void ProbeMutProcessor::AdjustLabels2HubVecUsingHub() {
	vector <tuple<unsigned int ,unsigned int, unsigned int>> hub_vec_number;

	for(unsigned int i=0; i < mut_hub_vector.size(); i++) {
	  unsigned int gereral_hub_num = 0;
	  unsigned int j = i;
		while (mut_hub_vector[j].abs_father_hub!=-1) {
			gereral_hub_num += 1;
			j = mut_hub_vector[j].index_father_hub;
		};

    if ((mut_hub_vector[j].abs_father_hub == -1) &&
				(mut_hub_vector[j].index_father_hub==0)){
						gereral_hub_num += 1;
				};

    std::tuple<unsigned int ,unsigned int, unsigned int> A;
    A = make_tuple(i , mut_hub_vector[i].abs_mut, gereral_hub_num);
		hub_vec_number.push_back(A);
		//hub_vec_number_mut[i] = gereral_mut_num;
	};

  std::sort(hub_vec_number.begin(), hub_vec_number.end(),
						[](const std::tuple<unsigned int ,unsigned int, unsigned int> &left,
							 const std::tuple<unsigned int ,unsigned int, unsigned int> &right) {
    return get<2>(left) < get<2>(right);
	});


  int label_index = 0;
	for(unsigned int i=0;i <hub_vec_number.size();i++) {
//		cout<<i<<" index in mut_hub_vector = " <<get<0>(hub_vec_number[i])<<"; absolute mut = "<<
//		get<1>(hub_vec_number[i]) <<"; gereral_hub_num = "<< get<2>(hub_vec_number[i])<<endl;
		for(unsigned int j=0; j <mut_hub_vector[get<0>(hub_vec_number[i])].number_cells;j++){
			mut_hub_vector[get<0>(hub_vec_number[i])].label_vec.push_back(label_index);
		  label_index++;
		};
	};


//	for(unsigned int i=0; i<mut_hub_vector.size(); i++) {
//		cout<<"for element with index "<< i << " in mut_hub_vector ";
//		for(unsigned int j=0; j<mut_hub_vector[i].label_vec.size(); j++)
//			cout<<" we have labels " << mut_hub_vector[i].label_vec[j];
//	  cout<<endl;
//	};
};

vector <vector<int>> ProbeMutProcessor::GetSPLMatrix() const {
	unsigned int n = get_cell_number();
	vector<vector<int>> dist_matrix(n, vector<int>(n, -1));
	for(unsigned int i_1 = 0; i_1 < mut_hub_vector.size(); i_1++) {
		for (unsigned int i_2 = 0; i_2 < mut_hub_vector[i_1].number_cells; i_2++) {
			for (unsigned int j_1 = 0; j_1 < mut_hub_vector.size(); j_1++) {
				for (unsigned int j_2 = 0; j_2 < mut_hub_vector[j_1].number_cells; j_2++) {
					dist_matrix[mut_hub_vector[i_1].label_vec[i_2]]
										 [mut_hub_vector[j_1].label_vec[j_2]] =
						GetDistBetweenCells(i_1, j_1);
				};
			};
		};
	};
//	for (unsigned int i = 0; i < n; i ++) {
//		for (unsigned int j = 0; j < n; j ++) {
//			cout<<"A["<<i<<"]["<<j<<"]=" << dist_matrix[i][j]<<";";
//		};
//		cout<<endl;
//  };
	return dist_matrix;
};

unsigned int ProbeMutProcessor::GetDistBetweenCells(unsigned int j_1,
																										unsigned int j_2) const {
	vector <MutationHub> first_vec;
	vector <MutationHub> second_vec;

	while (mut_hub_vector[j_1].abs_father_hub != -1) {
		first_vec.push_back(mut_hub_vector[j_1]);
		j_1 = mut_hub_vector[j_1].index_father_hub;
	};

	if ((mut_hub_vector[j_1].abs_father_hub == -1) &&
			(mut_hub_vector[j_1].index_father_hub == 0))
			first_vec.push_back(mut_hub_vector[j_1]);

	while (mut_hub_vector[j_2].abs_father_hub != -1) {
		second_vec.push_back(mut_hub_vector[j_2]);
		j_2 = mut_hub_vector[j_2].index_father_hub;
	};

	if ((mut_hub_vector[j_2].abs_father_hub == -1) &&
			(mut_hub_vector[j_2].index_father_hub == 0))
			second_vec.push_back(mut_hub_vector[j_2]);

  //for(unsigned int i = 0; i < first_vec.size(); i++)
	//	cout<<"first_vec["<<i<<"]"<< first_vec[i].

	unsigned int min_length = 0;
	if (first_vec.size() < second_vec.size()) { min_length = first_vec.size();
	} else min_length = second_vec.size();

	std::reverse(first_vec.begin(),first_vec.end());
	std::reverse(second_vec.begin(),second_vec.end());

	unsigned int length = 0;
	for (unsigned int i = 0; i < min_length; i++) {
		if (first_vec[i].abs_mut == second_vec[i].abs_mut) {
			length++;
		};
	};

	unsigned int dist_ = 0;
	for (unsigned int i = length; i < first_vec.size(); i++) {
		dist_+=first_vec[i].node_number;
	};

	return dist_;
};

void ProbeMutProcessor::ConstructorInit() {
  SetCellsGenotypes();// SaveGenotypes("BeforeUniqueGenotypes();",probe_genotypes);
  UniqueGenotypes(); //SaveProbeCells("after_unique_genotypes");
  UniqueProbeCells();
  GenotypesWithoutInclusions();
  CreatMutGenotypesConnections();
  CreatMutVector();
	CreatMutHubVector();
//	SaveMutGenotypeConnections((char *)"catalog");
//  SaveMutVector((char *)"catalog");
//  SaveGraph((char *)"catalog","my");
//  SaveHubGraph((char *)"catalog","my_hub");
//	SaveProbeCells("catalog");
//  SaveGenotypes("catalog", probe_genotypes);
//  SaveNewickTreeFormat((char *)"catalog");

 // SaveMutHubVector ((char *)"catalog");
};

ProbeMutProcessor::ProbeMutProcessor(const vector <Cell> & _probe_cells,
                                     unsigned int threshold_numb_cells) {
  raw_probe_cells = _probe_cells;
  probe_cells = _probe_cells;
  probe_cells.erase(probe_cells.begin() + threshold_numb_cells);
  // sorts cells by the length of genotype
  stable_sort (probe_cells.begin(), probe_cells.end(), CompareCells());
  ConstructorInit();
};

ProbeMutProcessor::ProbeMutProcessor(const vector <Cell> & _probe_cells) {
    raw_probe_cells = _probe_cells;
    probe_cells = _probe_cells;
    // sorts cells by the length of genotype
    //string file_name("before_sort");
    //SaveProbeCells(file_name.c_str());
    stable_sort (probe_cells.begin(), probe_cells.end(), CompareCells());
    //file_name = "after_sort";
    //SaveProbeCells(file_name.c_str());
    ConstructorInit();
};

ProbeMutProcessor::ProbeMutProcessor (const vector <Cell> & _probe_cells,
                   const vector <Genotype *> & genotype_vector) {
	probe_genotypes = genotype_vector;
	raw_probe_cells = _probe_cells;
	probe_cells = _probe_cells;
    // sorts cells by the length of genotype
    //string file_name("before_sort");
    //SaveProbeCells(file_name.c_str());
	stable_sort (probe_cells.begin(), probe_cells.end(), CompareCells());
    //file_name = "after_sort";
    //SaveProbeCells(file_name.c_str());
	UniqueGenotypes(); //SaveProbeCells("after_unique_genotypes");
  UniqueProbeCells();
  GenotypesWithoutInclusions();
  CreatMutGenotypesConnections();
  CreatMutVector();
	CreatMutHubVector();
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

	probe_mut_processor.FillVecChildPlaces();
	probe_mut_processor.CreatMutHubVector();

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
  return chi_square_sum ( degree_mut_vec_pr_1,
                          degree_mut_vec_pr_2 );
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

float t_dist_deg (const ProbeMutProcessor & Probe_mut_1,
									const ProbeMutProcessor & Probe_mut_2) {
    float chi_square_deg = chi_square_degree(Probe_mut_1, Probe_mut_2);
    return (chi_square_deg);
};

float t_dist_control (const ProbeMutProcessor & Probe_mut_1,
											const ProbeMutProcessor & Probe_mut_2){
  return 1;
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

float dist_robinson_and_foulds(const ProbeMutProcessor & Probe_mut_1,
															 const ProbeMutProcessor & Probe_mut_2, int seed) {
  string tree_str_1 = Probe_mut_1.get_str_newick_mut_tree_format();
	string tree_str_2 = Probe_mut_2.get_str_newick_mut_tree_format();

  return phylip_dist(tree_str_1, tree_str_2, 0, seed);
};


float dist_branch_score (const ProbeMutProcessor & Probe_mut_1,
												 const ProbeMutProcessor & Probe_mut_2, int seed) {
  string tree_str_1 = Probe_mut_1.get_str_newick_mut_tree_format();
	string tree_str_2 = Probe_mut_2.get_str_newick_mut_tree_format();
  return phylip_dist(tree_str_1, tree_str_2, 1, seed);
};

float nodal_distance (const ProbeMutProcessor & Probe_mut_1,
											const ProbeMutProcessor & Probe_mut_2) {

  vector <vector <int>> PL_matrix_1 = Probe_mut_1.GetSPLMatrix(); // Path_Lengths_matrix

  //AdjustLabels2HubVecUsingHub();
  vector <vector <int>> PL_matrix_2 = Probe_mut_2.GetSPLMatrix(); // Path_Lengths_matrix

	unsigned int nodal_dist = 0;
	unsigned int Matrix_Size = PL_matrix_1.size();
	for (unsigned int i=0; i < Matrix_Size; i++)
		for (unsigned int j=0; j < Matrix_Size; j++) {
			PL_matrix_1 [i] [j] = PL_matrix_1 [i] [j] + PL_matrix_1 [j] [i];
			PL_matrix_2 [i] [j] = PL_matrix_2 [i] [j] + PL_matrix_2 [j] [i];
			if (i<j) nodal_dist += abs(PL_matrix_1 [i] [j] - PL_matrix_2 [i] [j]);
		};
	//cout <<"nodal_dist = " << nodal_dist<<endl;
  return pow(nodal_dist, 0.5)/(Matrix_Size * Matrix_Size);
};

float splitted_nodal_distance (const ProbeMutProcessor & Probe_mut_1,
												       const ProbeMutProcessor & Probe_mut_2) {
  //AdjustLabels2HubVecUsingHub();
  vector <vector <int>> PL_matrix_1 = Probe_mut_1.GetSPLMatrix(); // Splitted_Path_Lengths

  vector <vector <int>> PL_matrix_2 = Probe_mut_2.GetSPLMatrix(); // Splitted_Path_Lengths

	unsigned int splitted_nodal_dist = 0;
	unsigned int Matrix_Size = PL_matrix_1.size();
	//dcout<<"PL_matrix_1.size()="<<PL_matrix_1.size()<< endl;
	for (unsigned int i=0; i < Matrix_Size; i++) {
		for (unsigned int j=0; j < Matrix_Size; j++) {
			splitted_nodal_dist += abs(PL_matrix_1 [i] [j] - PL_matrix_2 [i] [j]) ;
    };
	};
//  cout <<"matrix_1"<<endl;
//	for (unsigned int i=0; i < Matrix_Size; i++) {
//		for (unsigned int j=0; j < Matrix_Size; j++) {
//			cout<<" ["<<i<<"] ["<<j<<"]=" << PL_matrix_1 [i] [j]<<"; ";
//		}
//		cout<<endl;
//	}
//
//	cout <<"matrix_2"<<endl;
//	for (unsigned int i=0; i < Matrix_Size; i++) {
//		for (unsigned int j=0; j < Matrix_Size; j++) {
//			cout<<"["<<i<<"] ["<<j<<"]=" << PL_matrix_2 [i] [j]<<"; ";
//		}
//		cout<<endl;
//	}
//
//
//
//    cout <<"splitted_nodal_dist = " << pow(splitted_nodal_dist,0.5)
//                                       /(Matrix_Size*Matrix_Size)<<endl;
  return pow(splitted_nodal_dist, 0.5) / (Matrix_Size * Matrix_Size);
};

float dist_splitted_nodal (const ProbeMutProcessor & Probe_mut_1,
													 const ProbeMutProcessor & Probe_mut_2 ) {
	vector <vector <int>> PL_matrix_1 = Probe_mut_1.GetSPLMatrix(); // Splitted_Path_Lengths
  vector <vector <int>> PL_matrix_2 = Probe_mut_2.GetSPLMatrix(); // Splitted_Path_Lengths

	unsigned int max_length_1 = 0;
	unsigned int max_length_2 = 0;

	unsigned int matrix_size = PL_matrix_1.size();

	for(unsigned int i = 0; i < matrix_size; i++) {
		for(unsigned int j = 0; j < matrix_size; j++) {
				if (max_length_1 < PL_matrix_1[i][j]) max_length_1 = PL_matrix_1[i][j];
				if (max_length_2 < PL_matrix_2[i][j]) max_length_2 = PL_matrix_2[i][j];
		};
	};

  vector <unsigned int> splitted_lengths_vector_1 (max_length_1 + 1, 0);
  vector <unsigned int> splitted_lengths_vector_2 (max_length_2 + 1, 0);

	for(unsigned int i = 0; i < matrix_size; i++) {
		for(unsigned int j = 0; j < matrix_size; j++) {
			splitted_lengths_vector_1[PL_matrix_1[i][j]]++;
			splitted_lengths_vector_2[PL_matrix_2[i][j]]++;
		};
	};

  return chi_square_sum(splitted_lengths_vector_1,
												splitted_lengths_vector_2);
};

float dist_nodal (const ProbeMutProcessor & Probe_mut_1,
									const ProbeMutProcessor & Probe_mut_2) {
	vector <vector <int>> PL_matrix_1 = Probe_mut_1.GetSPLMatrix(); // Splitted_Path_Lengths
  vector <vector <int>> PL_matrix_2 = Probe_mut_2.GetSPLMatrix(); // Splitted_Path_Lengths

	unsigned int max_length_1 = 0;
	unsigned int max_length_2 = 0;

	unsigned int matrix_size = PL_matrix_1.size();

	for(unsigned int i = 0; i < matrix_size; i++) {
		for(unsigned int j = i; j < matrix_size; j++) {
			PL_matrix_1 [i] [j] = PL_matrix_1 [i] [j] + PL_matrix_1 [j] [i];
			PL_matrix_2 [i] [j] = PL_matrix_2 [i] [j] + PL_matrix_2 [j] [i];
			if (max_length_1 < PL_matrix_1[i][j]) max_length_1 = PL_matrix_1[i][j];
			if (max_length_2 < PL_matrix_2[i][j]) max_length_2 = PL_matrix_2[i][j];
		};
	};

  vector <unsigned int> splitted_lengths_vector_1 (max_length_1 + 1, 0);
  vector <unsigned int> splitted_lengths_vector_2 (max_length_2 + 1, 0);

	for(unsigned int i = 0; i < matrix_size; i++) {
		for(unsigned int j = i; j < matrix_size; j++) {
			splitted_lengths_vector_1[PL_matrix_1[i][j]]++;
			splitted_lengths_vector_2[PL_matrix_2[i][j]]++;
		};
	};

  return chi_square_sum(splitted_lengths_vector_1,
												splitted_lengths_vector_2);
};

float t_dist_deg_vs_dist_nodal(const ProbeMutProcessor & Probe_mut_1,
															 const ProbeMutProcessor & Probe_mut_2) {

    float chi_square_splitted_nodal = dist_splitted_nodal(Probe_mut_1, Probe_mut_2);
    float chi_square_deg = chi_square_degree(Probe_mut_1, Probe_mut_2);
    return (chi_square_splitted_nodal + chi_square_deg);

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

int GetMinxProbePiece(const BasicSimParameters & pars, const int index) {
    //makes piece_min_x 0 or bigger
    int del_x = pars.probe.delta_area.del_x;
    /*probe_area.max_x - probe_area.min_x;*/
    int del_y = pars.probe.delta_area.del_y;
    //  projection of vector with coords (?,i*del_y,0)
    //  shift vector (max_x,0,0) to (?, del_y, 0) and find ? = piece_min_x
    float square_rute = pow(pow(del_x, 2.0) - pow(index*del_y, 2.0), 0.5);  /*10*/
    int piece_min_x = del_x - int (square_rute);
    return piece_min_x;
};


void InitProbePieceVector (const AllSimParameters & pars,
													 vector <ProbePiece> & probe_piece_vector){
  int piece_num = pars.inference_parameters.comparison_parameters.piece_num;
  for (int i = 0; i < piece_num; i++) {
      ProbePiece local_probe_piece;
      local_probe_piece.indicator = 0;
      local_probe_piece.cell_vector.clear();
      local_probe_piece.piece_min_x_vec =
				GetMinxProbePiece(pars.basic_sim_parameters, i);
      probe_piece_vector.push_back(local_probe_piece);
  };
};

void GetVecOfProbePiece(const AllSimParameters & pars,
												vector <ProbePiece> & probe_piece_vector) {
	int piece_number = pars.inference_parameters.comparison_parameters.piece_num;
	int cell_number = pars.basic_sim_parameters.probe.cell_num;
  for (unsigned int i=0; i < cells.size(); i++) {  // conditional loop
    for (int j=0; j < piece_number; j++) {
      if (IsCellFittedToPiece(i, j,
                              probe_piece_vector.at(j).piece_min_x_vec,
                              pars.basic_sim_parameters.probe.area)
					&& (probe_piece_vector.at(j).indicator <
          cell_number) ) {
        probe_piece_vector.at(j).cell_vector.push_back(cells.at(i));
        probe_piece_vector.at(j).indicator++;
        break;
      };
    };

    bool is_search_finished = true;
		for (int j=0; j < piece_number; j++) {
      if (probe_piece_vector.at(j).indicator <
            cell_number) is_search_finished = false;
    };
    if (is_search_finished) break;
  };
};

void GetVecOfRandProbePiece(const AllSimParameters & pars,
														vector <ProbePiece> & probe_piece_vector) {
	int piece_number = pars.inference_parameters.comparison_parameters.piece_num;
  int cell_number = pars.basic_sim_parameters.probe.cell_num;
  for (int i = 0; i < piece_number; i++) {
		for (int j = 0; j < cell_number; j++) {
			unsigned int index = trunc((cells.size()-1)*_drand48());
			probe_piece_vector.at(i).cell_vector.push_back(cells.at(index));
		};
  };
};


void SetComparisonResults2File(const AllSimParameters & pars,
                               const vector <ProbePiece> & probe_piece_vector) {
 //opens file_tree_comparison.dat
  ofstream comparison_file;

  string full_comparison_file_name =
				  pars.basic_sim_parameters.related_path_to_par_file
					+ pars.basic_sim_parameters.etalon_result_catalogue
					+ GetSystemRelevantSlash()
					+ pars.inference_parameters.basic_settings.current_result_file_name;

  comparison_file.open(full_comparison_file_name.c_str(), std::ios::app);
	int piece_number = pars.inference_parameters.comparison_parameters.piece_num;
  for (int i=0; i < piece_number; i++) {
   // if (probe_piece_vector.at(i).results.dist < pars.inference_parameters.comparison_parameters.max_dist_probe_trees) {
      comparison_file.precision(5);
      comparison_file.setf(std::ios::fixed, std:: ios::floatfield);
      comparison_file << pars.basic_sim_parameters.evolution.driver_adv <<"\t"<<
        pars.basic_sim_parameters.evolution.mutation_rate               <<"\t"<<
        pars.basic_sim_parameters.evolution.driver_mutation_rate        <<"\t"<<
        probe_piece_vector.at(i).results.dist                           <<"\t"<<
        probe_piece_vector.at(i).results.delta_mut                      <<"\t"<<
        probe_piece_vector.at(i).results.num_mut                        <<endl;
//      cout<<"all info succ_prob:  "<<pars.evolution.driver_adv <<
//        " "<<pars.evolution.mutation_rate                      <<
//        " "<<pars.evolution.driver_mutation_rate               <<
//        " "<<probe_piece_vector.at(i).results.dist            <<
//        " "<<probe_piece_vector.at(i).results.delta_mut       <<
//        " "<<probe_piece_vector.at(i).results.num_mut         <<endl;
//    } else {
 //     cout<<"dist("<<i<<")=" << probe_piece_vector.at(i).results.dist <<endl;
//    };
  };
  comparison_file.close();
}

ComparisonResults GetComparisonResults( const AllSimParameters & pars,
																			  const ProbeMutProcessor & probe,
                                        const ProbeMutProcessor & e_probe) {

  ComparisonResults c_results;

	int seed = pars.basic_sim_parameters.seed;
  switch (pars.inference_parameters.comparison_parameters.dist_type) {
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
		case 4: {
			c_results.dist = dist_robinson_and_foulds (probe, e_probe, seed);
      break;
		}
		case 5: {
			c_results.dist = dist_branch_score (probe, e_probe, seed);
      break;
		}
		case 6: {
			c_results.dist = nodal_distance (probe, e_probe);
			break;
		}
		case 7: {
			c_results.dist = splitted_nodal_distance (probe, e_probe);
			break;
		}
		case 8: {
			c_results.dist = dist_nodal (probe, e_probe);
			break;
		}
		case 9: {
			c_results.dist = dist_splitted_nodal (probe, e_probe);
			break;
		}
		case 10: {
			c_results.dist = t_dist_deg_vs_dist_nodal(probe, e_probe);
			break;
		}
		case 11: {
			c_results.dist = t_dist_deg(probe, e_probe);
		  break;
		}
		case 12: {
			c_results.dist = t_dist_control(probe, e_probe);
			break;
		}

     default:
     err("Mistake in distance type (see setting file)");
   };
	//probe.SaveGraph("pic","mut_graph");

//	cout<<"0_case = " << t_dist_deg_vs_dep(probe, e_probe)<<endl;;
//	cout<<"1_case = " << t_dist_deg_vs_mean_pass_non_elem(probe, e_probe)<<endl;
//	cout<<"2_case = " << t_dist_deg_vs_mean_pass_to_cell(probe, e_probe)<<endl;;
//	cout<<"3_case = " << t_dist_mean_pass_to_cell (probe, e_probe)<<endl;
//  cout<<"4_case = " << dist_robinson_and_foulds (probe, e_probe)<<endl;
//	cout<<"5_case = " << dist_branch_score (probe, e_probe)<<endl;
//  cout<<"6_case = " << nodal_distance (probe, e_probe)<<endl;
//  cout<<"7_case = " << splitted_nodal_distance (probe, e_probe)<<endl;
//  cout<<"8_case = " << dist_nodal (probe, e_probe)<<endl;
//  cout<<"9_case = " << dist_splitted_nodal (probe, e_probe)<<endl;
//  cout<<"10_case = " << t_dist_deg_vs_dist_nodal(probe, e_probe)<<endl;


  c_results.e_num_mut = e_probe.get_mut_number();
  c_results.num_mut   = probe.get_mut_number();
  c_results.delta_mut = abs(c_results.e_num_mut - c_results.num_mut);

  //cout<<"tree dist is "                << c_results.dist      << endl;
  //cout<<"Number of etalon mutations "  << c_results.e_num_mut << endl;
  //cout<<"Number of probe mutations "   << c_results.num_mut   << endl;
  //cout<<"Difference in mut number is " << c_results.delta_mut << endl;

  return c_results;
};


bool IsProbePieceCompatible2EProbe(const AllSimParameters  & pars,
																	 const ProbeMutProcessor & probe,
																	 const ProbeMutProcessor & e_probe) {
  int probe_num_cells = probe.get_cell_number();
  int e_probe_num_cells = e_probe.get_cell_number();
  int delta_size_cells = abs(probe_num_cells - e_probe_num_cells);
  float size_level = e_probe_num_cells *
				pars.inference_parameters.comparison_parameters.threshold_frac_cell_diff;

  if (delta_size_cells <= size_level) {
		// cout<<"successful probe (there are enough cells ) e_probe_cell= "
    //     <<e_probe.get_cell_number()<<"; probe_cell ="
    //     << probe.get_cell_number() <<endl;

    int probe_mut_num     =  probe.get_mut_number();
    int e_probe_mut_num   =  e_probe.get_mut_number();
    int delta_probe_e_mut =  abs(e_probe_mut_num - probe_mut_num);

    if (delta_probe_e_mut < e_probe_mut_num *
                            pars.inference_parameters.
															comparison_parameters.threshold_frac_mut_diff) {
      return true;
    } else {
    	std::cout
    	<<pars.inference_parameters.comparison_parameters.threshold_frac_cell_diff<<std::endl;
      std::cout<<probe_mut_num<< "-" <<  e_probe_mut_num <<" muts"<<std::endl;

      //cout<<"There are " << probe_mut_num <<
      //   " muts and the demand is " << e_probe_mut_num <<" muts"<<endl;
      return false;
    };

  } else {
      std::cout<<probe_num_cells<<"-"<<  e_probe_num_cells <<" cells"<<std::endl;
    //cout<<"there are "<< probe_num_cells <<
          //" and the demand is " << e_probe_num_cells <<" cells"<<endl;
    return false;
  };
};

void ProcessProbePieceVector(const AllSimParameters & pars,
														 vector <ProbePiece> & probe_piece_vector){
  std::string path(pars.basic_sim_parameters.related_path_to_par_file +
									 pars.basic_sim_parameters.etalon_result_catalogue +
									 GetSystemRelevantSlash() +
									 pars.basic_sim_parameters.etalon_result_file_name);
  ifstream file_mut_vec;
  file_mut_vec.open(path.c_str());
  if (file_mut_vec.is_open()) {
    ProbeMutProcessor e_probe;  // etalon probe
    file_mut_vec >> e_probe;
    int piece_number = pars.inference_parameters.comparison_parameters.piece_num;
    for (int i=0; i < piece_number; i++) {
      ProbeMutProcessor probe(probe_piece_vector.at(i).cell_vector);
			//std::cout<<"piece_num = " << i << std::endl;
			//std::cout<<"mean mut number of cells in e_probe = " <<
			//	e_probe.get_mean_mut_num_cell() << std::endl;
			//std::cout<<"get_mean_pass_b_w_non_elem_nodes for e_probe = " <<
			//	e_probe.get_mean_pass_b_w_non_elem_nodes()<<std::endl;
			//std::cout<<"mean mut number of cells in piece_probe = " <<
			//	probe.get_mean_mut_num_cell() << std::endl;
			//std::cout<<"get_mean_pass_b_w_non_elem_nodes for piece_probe = " <<
			//	probe.get_mean_pass_b_w_non_elem_nodes()<<std::endl;
      bool is_compatible =
        IsProbePieceCompatible2EProbe(pars, probe, e_probe);
      if (is_compatible) {
        probe_piece_vector.at(i).results = GetComparisonResults(pars, probe, e_probe);
      };
    };
    SetComparisonResults2File(pars, probe_piece_vector);
    file_mut_vec.close();
  } else {
    file_mut_vec.close();
  };

};

void DoesInferenceAnalyze(AllSimParameters & pars) {
  vector <ProbePiece> probe_piece_vector;
  InitProbeArea(pars.basic_sim_parameters);
  InitProbePieceVector(pars, probe_piece_vector);
  if (pars.basic_sim_parameters.probe.is_random) {
		GetVecOfRandProbePiece(pars, probe_piece_vector);
  } else {
    GetVecOfProbePiece(pars, probe_piece_vector);
  };
  ProcessProbePieceVector(pars, probe_piece_vector);
};
// cells form area @par probe_area and write them in cell_vector
void GetEtalonCellVector(const BasicSimParameters & pars,
                         vector <Cell> & cell_vector) {
  int j = 0;//@j counter for fitted cells
  //for_loop: Find fitted cells (number less then threshold_num_cells)
  for (unsigned int i=0; i < cells.size() &&
                         j < pars.probe.cell_num; ++i) {
    bool is_x_in_range =
      (cells[i].x >= pars.probe.area.min_x) &&
      (cells[i].x <= pars.probe.area.max_x);
    bool is_y_in_range =
      (cells[i].y >= pars.probe.area.min_y) &&
      (cells[i].y <= pars.probe.area.max_y);
    bool is_z_in_range =
      (cells[i].z >= pars.probe.area.min_z) &&
      (cells[i].z <= pars.probe.area.max_z);

    if ( is_x_in_range && is_y_in_range && is_z_in_range) {
            cell_vector.push_back(cells[i]);
            j++;
    };
  };
  //begin delete this
  //cout<<"number of sells in simulation "<<cells.size()<<endl;
  //end delete this
};

void GetRandomEtalonCellVector(const BasicSimParameters & pars,
                               vector <Cell> & cell_vector) {
	for (int i=0; i < pars.probe.cell_num; i++) {
		int index = trunc((cells.size()-1)*_drand48());
		cell_vector.push_back(cells.at(index));
	};
};

void small_test(ProbeMutProcessor & probe_mut_processor,
                int size_vector_cells) {
    int num_of_cells_e_prob = probe_mut_processor.get_cell_number();
    //cout<<"num_of_cells_e_prob="<<num_of_cells_e_prob<<endl;
    //cout<<"size_vector_cells="<<size_vector_cells<<endl;
    assert(size_vector_cells == num_of_cells_e_prob);
};

//Inits probe area. Takes cuboid with del_x, del_y, del_z sides in
// x,y,z direction. IN X DIRECTION IT GOES FROM BORDER SIDE MAX_X
void InitProbeArea(BasicSimParameters & pars) {
    pars.probe.area.min_x = pars.probe.population_borders.max_x -
                                 pars.probe.delta_area.del_x;
    pars.probe.area.max_x = pars.probe.population_borders.max_x;
    pars.probe.area.min_y = - int (pars.probe.delta_area.del_y*1.0/2);
    pars.probe.area.max_y = int (pars.probe.delta_area.del_y*1.0/2);
    pars.probe.area.min_z = - int (pars.probe.delta_area.del_z*1.0/2);
    pars.probe.area.max_z = int (pars.probe.delta_area.del_z*1.0/2);
};

void SetProbeResultsToFile(const BasicSimParameters & pars,
													 ProbeMutProcessor & probe_mut_processor) {
  string path = pars.related_path_to_par_file +
								pars.etalon_result_catalogue +
								GetSystemRelevantSlash() +
								pars.etalon_result_file_name;
  ofstream etalon_file_mut_vec;
  etalon_file_mut_vec.open (path.c_str());
  //cout<<file_name_part<<"/etalon_mut_vec.dat is opened"<<endl;
  etalon_file_mut_vec << probe_mut_processor;
  etalon_file_mut_vec.close();
};

void SetEtalonProbe(BasicSimParameters & pars) {
  InitProbeArea(pars);
  vector <Cell> cell_vector;
  cell_vector.clear();  // vector preparation
  if (pars.probe.is_random) {  // init cell_vector
		GetRandomEtalonCellVector(pars, cell_vector);
  } else {
		GetEtalonCellVector(pars, cell_vector);
  };
  if (cell_vector.size() > 0) {
    ProbeMutProcessor probe_mut_processor(cell_vector);
    //small_test(probe_mut_processor, cell_vector.size());
    SetProbeResultsToFile(pars, probe_mut_processor);
  } else {
    //cout<<"there is no cells in the etalon probe"<< endl;
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
  for (unsigned int i=0; i<A.size(); i++) add_one_genotype(A.at(i));
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











