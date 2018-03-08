// Copyright 2017 Lebid Mykola all rights reserved
#include "sample_processing.h"

#include <set>
#include <numeric>

#include "classes.h"
#include "treedist.h"


void TreeMutAnalyzer::UniqueProbeCells() {
  std::vector<Cell>::iterator it = unique (probe_cells.begin(), probe_cells.end(),
                                           CompareCellsEq());
  probe_cells.resize(std::distance(probe_cells.begin(),it));
};



unsigned int TreeMutAnalyzer::get_cell_number() const {
  int cell_num = 0;
  for (unsigned int i=0; i<mutation_vector.size(); i++) {
    cell_num += mutation_vector[i]. number_cells;
  };
  return cell_num;
};

unsigned int TreeMutAnalyzer::get_mut_number() const {
 return mutation_vector.size();
};

void TreeMutAnalyzer::UniqueGenotypes() {
  std::vector<Genotype *>::iterator it = unique (probe_genotypes.begin(),
                                                 probe_genotypes.end(),
                                                 CompareGenotypesEq());
  probe_genotypes.resize(std::distance(probe_genotypes.begin(),it));
};

void TreeMutAnalyzer::SetCellsGenotypes() {
    for (unsigned int i=0;i<probe_cells.size(); i++) {
            probe_genotypes.push_back(genotypes[probe_cells[i].gen]);
            probe_genotypes[i]->number = 1;
    };
};
void TreeMutAnalyzer::GenotypesWithoutInclusions() {
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


int TreeMutAnalyzer::get_number_of_intrsections(Genotype * a, Genotype * b) {
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

void TreeMutAnalyzer::CreatMutGenotypesConnections() {
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

void TreeMutAnalyzer::SaveMutGenotypeConnections(char* name) {
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

void TreeMutAnalyzer::SaveMutVector(char* name) {
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

void TreeMutAnalyzer::SaveMutHubVector (char *name) {

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

std::string GetStringFromAbsMutNum(const int abs_mut_num){
	char buffer [50];
	//int n = 0 ;

  if ((abs_mut_num & RESISTANT_PM)==0) {
			if ((abs_mut_num & DRIVER_PM)==0) {
				  int n = sprintf (buffer, "p%d", abs_mut_num);
          std::string str(buffer);
          str = str;
          return str;
			} else {
					int n = sprintf (buffer, "d%d", (abs_mut_num & (~DRIVER_PM)));
          std::string str(buffer);
          str = str;
          return str;};
	 } else {
					int n = sprintf (buffer, "r%d * %d", (abs_mut_num & (~RESISTANT_PM)),abs_mut_num);
          std::string str(buffer);
          str = str;
          return str;};
}


std::string TreeMutAnalyzer::FormCellString(int mut_vec_index)  const {
	int num_cells = mutation_vector[mut_vec_index].number_cells;

	std::string str/*=")"*/;
	if( num_cells > 0 ) {
		int current_cell_index = 0;
		for (int i = 0; i < mut_vec_index; i++)
			current_cell_index += mutation_vector[i].number_cells;

		for (int i = 0; i < num_cells; i ++ ) {
			char buffer [30];
			int n = sprintf (buffer, "%d", current_cell_index + i);
			std::string str1(buffer);
			str = str + str1 + ":1.0";
			if (i != (num_cells-1))
				str = str + ",";
			//std::cout<<"local str="<<str<<std::endl;
		};
	};
	//current_cell_index = current_cell_index + num_cells;
	return str;
};
std::string TreeMutAnalyzer::GetCurrentMutationString(int mut_vec_index) const {
	/*unsigned int abs_mut = mutation_vector[mut_vec_index].abs_mut;*/
	std::string str ="):1.0" /*= GetStringFromAbsMutNum(abs_mut)*/;

	//if( mutation_vector[mut_vec_index].num_children_mut > 0 ) {
	str = FormCellString(mut_vec_index) + str;
	int num_cells = mutation_vector[mut_vec_index].number_cells;
  int mut_num = mutation_vector[mut_vec_index].num_children_mut;
  if (num_cells > 0 && mut_num > 0) str = "," + str;
	for (int i = 0; i <mut_num; i++) {
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

std::string TreeMutAnalyzer::get_str_newick_mut_tree_format() const {
	std::string str=";"/*="0d;"*/;
	str = "):1.0" + str;

  str = FormCellString(0) + str;
  int num_cells = mutation_vector[0].number_cells;
  int mut_num = mutation_vector[0].num_children_mut;
  if (num_cells > 0 && mut_num > 0) str = "," + str;
	for(int i = 0 ; i < mut_num; i++) {
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

void TreeMutAnalyzer::SaveNewickTreeFormat(char *name) {
  char name_file[256];
  sprintf(name_file,"newick_%s.dat",name);
  FILE *f=fopen(name_file, "w") ;
  if (f==NULL) err("err err with newick");
	//vector <bool> is_written(size_mv,false);
	std::string str=";"/*"0d;"*/;
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

void TreeMutAnalyzer::SaveGraph(const char * catalog_name,
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


void TreeMutAnalyzer::SaveHubGraph(const char * catalog_name,
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
int TreeMutAnalyzer::GetNumCellsForCurNode(int current_genotype_mut,
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
int TreeMutAnalyzer::GetNumChildrenCurNode(int index_cur_probe_genotype,
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
void TreeMutAnalyzer::InitZeroMutationNode(MutationNode& zero_mutation) {
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


int TreeMutAnalyzer::GetFatherMutIndexForCurNode(
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
MutationNode TreeMutAnalyzer::GetCurrentMutNode(
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

void TreeMutAnalyzer::FillVecChildPlaces() {
	for (unsigned int i = 0 ; i < mutation_vector.size(); i++) {
		for (unsigned int j = i + 1 ; j < mutation_vector.size(); j++) {
			if (mutation_vector[i].abs_mut == mutation_vector[j].abs_father_mut)
				mutation_vector[i].vec_child_places.push_back(j);
		};
	};
}

// @function, small description:
//  of attachment of genotype
void TreeMutAnalyzer::CreatMutVector() {
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

void TreeMutAnalyzer::ChangeDepth(
												std::vector <MutationNode> & local_mut_vec,
												unsigned int index) {

	local_mut_vec.at(index).depth--;

	for (unsigned int i = 0; i < local_mut_vec.at(index).num_children_mut; i++) {
		unsigned int next_index = local_mut_vec[index].vec_child_places[i];
		ChangeDepth(local_mut_vec, next_index);
	};
};

void TreeMutAnalyzer::SetIndexFatherHubVec() {
	//current_mut_hub.abs_father_hub = mut_vec.at(i).abs_father_mut;
	//mut_hub_vector[0].index_father_hub = - 1;
	for (unsigned int i = 0; i < mut_hub_vector.size(); i++) {
		for (unsigned int j = 0; j < mut_hub_vector.size(); j++)
			if (mut_hub_vector[i].abs_father_hub == mut_hub_vector[j].abs_mut)
				mut_hub_vector[i].index_father_hub = j;
	};
};

void TreeMutAnalyzer:: CreatMutHubVector() {
	std::vector <MutationNode> mut_vec(mutation_vector);
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

void TreeMutAnalyzer::AdjustLabels2HubVecUsingMut() {
	std::vector <std::tuple<unsigned int, unsigned int, unsigned int>> hub_vec_number_mut;

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

    std::tuple<unsigned int ,unsigned int, unsigned int> A;
    A = std::make_tuple(i , mut_hub_vector[i].abs_mut, gereral_mut_num);

		hub_vec_number_mut.push_back(A);
		//hub_vec_number_mut[i] = gereral_mut_num;
	};

  std::sort(hub_vec_number_mut.begin(), hub_vec_number_mut.end(),
						[](const std::tuple<unsigned int ,unsigned int, unsigned int> &left,
							 const std::tuple<unsigned int ,unsigned int, unsigned int> &right) {
    return std::get<2>(left) < std::get<2>(right);
	});


  int label_index = 0;
	for(unsigned int i=0; i < hub_vec_number_mut.size(); i++ ) {
		//cout<<i<<" index in mut_hub_vector = "<< get<0>(hub_vec_number_mut[i])<<"; absolute mut = "<<
		//get<1>(hub_vec_number_mut[i]) <<"; gereral_hub_num = "<< get<2>(hub_vec_number_mut[i])<<endl;
		unsigned int current_index = std::get<0>(hub_vec_number_mut[i]);
		unsigned int current_cell_num = mut_hub_vector[current_index].number_cells;
		for(unsigned int j=0; j < current_cell_num; j++) {
			mut_hub_vector[std::get<0>(hub_vec_number_mut[i])].label_vec.push_back(label_index);
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


void TreeMutAnalyzer::AdjustLabels2HubVecUsingHub() {
	std::vector <std::tuple<unsigned int ,unsigned int, unsigned int>> hub_vec_number;

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
    A = std::make_tuple(i, mut_hub_vector[i].abs_mut, gereral_hub_num);
		hub_vec_number.push_back(A);
		//hub_vec_number_mut[i] = gereral_mut_num;
	};

  std::sort(hub_vec_number.begin(), hub_vec_number.end(),
						[](const std::tuple<unsigned int ,unsigned int, unsigned int> &left,
							 const std::tuple<unsigned int ,unsigned int, unsigned int> &right) {
    return std::get<2>(left) < std::get<2>(right);
	});


  int label_index = 0;
	for(unsigned int i=0;i <hub_vec_number.size();i++) {
//		cout<<i<<" index in mut_hub_vector = " <<get<0>(hub_vec_number[i])<<"; absolute mut = "<<
//		get<1>(hub_vec_number[i]) <<"; gereral_hub_num = "<< get<2>(hub_vec_number[i])<<endl;
		unsigned int current_index = std::get<0>(hub_vec_number[i]);
		unsigned int num_cells = mut_hub_vector[current_index].number_cells;
		for(unsigned int j=0; j < num_cells;j++){
			mut_hub_vector[current_index].label_vec.push_back(label_index);
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

std::vector <std::vector<int>> TreeMutAnalyzer::GetSPLMatrix() const {
	unsigned int n = get_cell_number();
	std::vector<std::vector<int>> dist_matrix(n, std::vector<int>(n, -1));
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

unsigned int TreeMutAnalyzer::GetDistBetweenCells(unsigned int j_1,
																										unsigned int j_2) const {
	std::vector <MutationHub> first_vec;
	std::vector <MutationHub> second_vec;

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

void TreeMutAnalyzer::Constructor_init() {
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

TreeMutAnalyzer::TreeMutAnalyzer(const std::vector <Cell> & _probe_cells,
                                     unsigned int threshold_numb_cells) {
  raw_probe_cells = _probe_cells;
  probe_cells = _probe_cells;
  probe_cells.erase(probe_cells.begin() + threshold_numb_cells);
  // sorts cells by the length of genotype
  std::stable_sort (probe_cells.begin(), probe_cells.end(), CompareCells());

  Constructor_init();
};

TreeMutAnalyzer::TreeMutAnalyzer(const std::vector <Cell> & _probe_cells) {
    raw_probe_cells = _probe_cells;
    probe_cells = _probe_cells;
    // sorts cells by the length of genotype
    //string file_name("before_sort");
    //SaveProbeCells(file_name.c_str());
    std::stable_sort (probe_cells.begin(), probe_cells.end(), CompareCells());
    //file_name = "after_sort";
    //SaveProbeCells(file_name.c_str());
    Constructor_init();
};

TreeMutAnalyzer::TreeMutAnalyzer (const std::vector <Cell> & _probe_cells,
                   const std::vector <Genotype *> & genotype_vector) {
	probe_genotypes = genotype_vector;
	raw_probe_cells = _probe_cells;
	probe_cells = _probe_cells;
    // sorts cells by the length of genotype
    //string file_name("before_sort");
    //SaveProbeCells(file_name.c_str());
	std::stable_sort (probe_cells.begin(), probe_cells.end(), CompareCells());
    //file_name = "after_sort";
    //SaveProbeCells(file_name.c_str());
	UniqueGenotypes(); //SaveProbeCells("after_unique_genotypes");
  UniqueProbeCells();
  GenotypesWithoutInclusions();
  CreatMutGenotypesConnections();
  CreatMutVector();
	CreatMutHubVector();
};


TreeMutAnalyzer::TreeMutAnalyzer() {
  //  cout <<"forgot to put parameters"<<endl;
};

TreeMutAnalyzer::~TreeMutAnalyzer() {

};

void TreeMutAnalyzer::SaveProbeCells(const char *name) {
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

void TreeMutAnalyzer::SaveGenotypes(char *name,
																		std::vector<Genotype *> & probe_genotypes_)
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
                          const TreeMutAnalyzer & probe_mut_processor) {
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
                         TreeMutAnalyzer & probe_mut_processor) {
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

unsigned int TreeMutAnalyzer::get_tree_mut_depth() const {
  unsigned int max_depth = 0;
  for(auto it = mutation_vector.begin(); it < mutation_vector.end(); ++it)
    max_depth = (max_depth < it->depth) ? it->depth : max_depth;
  return max_depth;
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

float chi_square_sum(const std::vector <unsigned int> vec_1,
                     const std::vector <unsigned int> vec_2) {
  std::vector <unsigned int> vec1(vec_1);
  std::vector <unsigned int> vec2(vec_2);
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


void GetCumulativeDists (const std::vector <unsigned int> vec_1,
												 const std::vector <unsigned int> vec_2,
												 std::vector <float> & comulative_vec_1,
												 std::vector <float> & comulative_vec_2) {
  std::vector <unsigned int> vec1(vec_1);
  std::vector <unsigned int> vec2(vec_2);

  unsigned int mut_num_1 = std::accumulate(vec1.begin(),vec1.end(), 0);
  unsigned int mut_num_2 = std::accumulate(vec2.begin(),vec2.end(), 0);

  // the length of step
  if (vec1.size() < vec2.size()) {
    unsigned int diff_langth = vec2.size() - vec1.size();
    vec1.insert(vec1.end(), diff_langth, 0);
  } else {
    unsigned int diff_langth = vec1.size() - vec2.size();
    vec2.insert(vec2.end(), diff_langth, 0);
  };

	float cumulative_num_mutations_1 = 0;
	for (auto current: vec1) {
			cumulative_num_mutations_1 += (1.0 * current)/(1.0 * mut_num_1);
			comulative_vec_1.push_back(cumulative_num_mutations_1);
	};

	float cumulative_num_mutations_2 = 0;
	for (auto current: vec2) {
			cumulative_num_mutations_2 += (1.0 * current)/(1.0 * mut_num_2);
			comulative_vec_2.push_back(cumulative_num_mutations_2);
	};

};


float KuiperDist(const std::vector <unsigned int> vec_1,
								 const std::vector <unsigned int> vec_2) {
	std::vector <float> comulative_vec_1, comulative_vec_2;
	GetCumulativeDists (vec_1, vec_2, comulative_vec_1, comulative_vec_2);
	unsigned int size_of_vecs = comulative_vec_1.size();
	float max_up = 0;
	float max_down = 0;
	for (unsigned int i = 0; i < size_of_vecs; i ++){
			float current = comulative_vec_1[i] - comulative_vec_2[i];
			if  (max_up < current) max_up = current;
			if  (max_down > current) max_down = current;
 	};

	float max_diff = max_up - max_down;

	float multiplier = sqrt((1.0*size_of_vecs)/2.0);
	return multiplier*max_diff;
};
float KolmogorovSmirnovDist(const std::vector <unsigned int> vec_1,
								            const std::vector <unsigned int> vec_2){
	std::vector <float> comulative_vec_1, comulative_vec_2;
	GetCumulativeDists (vec_1, vec_2, comulative_vec_1, comulative_vec_2);
	unsigned int size_of_vecs = comulative_vec_1.size();

	float max_diff = 0;
	for (unsigned int i = 0; i < size_of_vecs; i ++){
			float current = std::abs(comulative_vec_1[i] - comulative_vec_2[i]);
			if  (max_diff < current) max_diff = current;
 	};

	float multiplier = sqrt((1.0*size_of_vecs)/2.0);
	return multiplier*max_diff;
};

float CramerMisesDist(const std::vector <unsigned int> vec_1,
											const std::vector <unsigned int> vec_2) {
	std::vector <float> comulative_vec_1, comulative_vec_2;
	GetCumulativeDists (vec_1, vec_2, comulative_vec_1, comulative_vec_2);
	unsigned int size_of_vecs = comulative_vec_1.size();

	float sum_sq_diff = 0;
	for (unsigned int i = 0; i < size_of_vecs; i ++){
			float diff = std::abs(comulative_vec_1[i] - comulative_vec_2[i]);
			float sq_diff = diff * diff;
			sum_sq_diff += sq_diff;
	};
	return sum_sq_diff;
};



float chi_square_depth(const TreeMutAnalyzer & probe_1,
                          const TreeMutAnalyzer & probe_2) {
  std::vector <unsigned int> depth_mut_vec_pr_1 = probe_1.get_depth_mut_vec();
  std::vector <unsigned int> depth_mut_vec_pr_2 = probe_2.get_depth_mut_vec();
  return chi_square_sum (depth_mut_vec_pr_1,
                         depth_mut_vec_pr_2);
};



unsigned int TreeMutAnalyzer::get_tree_mut_degree() const {
  unsigned int max_degree = 0;
  for(auto it = mutation_vector.begin(); it < mutation_vector.end(); ++it)
	max_degree = (max_degree < it->num_children_mut) ?
                                it->num_children_mut : max_degree;
  return max_degree;
};



std::vector <unsigned int> TreeMutAnalyzer::get_degree_mut_vec() const {
  unsigned int tree_mut_degree = get_tree_mut_degree();
  std::vector <unsigned int> degree_mut_vector(tree_mut_degree + 1, 0);
  unsigned int mut_number = get_mut_number();
  for (unsigned int i=0; i < mut_number; i++)
    degree_mut_vector[mutation_vector[i].num_children_mut]++;
  return degree_mut_vector;
};


std::vector <unsigned int> TreeMutAnalyzer::get_depth_mut_vec() const {
  unsigned int tree_mut_depth = get_tree_mut_depth();
  std::vector <unsigned int> depth_mut_vector(tree_mut_depth + 1, 0);
  unsigned int mut_number = get_mut_number();
  for (unsigned int i=0; i < mut_number; i++)
    depth_mut_vector[mutation_vector[i].depth]++;
  return depth_mut_vector;
};

float chi_square_degree (const TreeMutAnalyzer & probe_1,
												 const TreeMutAnalyzer & probe_2) {
  std::vector <unsigned int> degree_mut_vec_pr_1 = probe_1.get_degree_mut_vec();
  std::vector <unsigned int> degree_mut_vec_pr_2 = probe_2.get_degree_mut_vec();
  return chi_square_sum ( degree_mut_vec_pr_1,
                          degree_mut_vec_pr_2 );
}

float TreeMutAnalyzer::get_mean_mut_num_cell() const {
  unsigned int sum = 0;
  for(const MutationNode & m_n: mutation_vector){
		if (m_n.number_cells>0) sum += m_n.number_cells*m_n.depth;
  };
  return (1.*sum/get_cell_number());
}


float TreeMutAnalyzer::get_mean_pass_b_w_non_elem_nodes() const {
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

float dist_mean_pass_b_w_non_elem_nodes(const TreeMutAnalyzer & Probe_mut_1,
                                        const TreeMutAnalyzer & Probe_mut_2) {
	float mean_sum = (Probe_mut_1.get_mean_pass_b_w_non_elem_nodes() +
	                  Probe_mut_2.get_mean_pass_b_w_non_elem_nodes())/2;
	float abs_diff = abs( Probe_mut_1.get_mean_pass_b_w_non_elem_nodes() -
											  Probe_mut_2.get_mean_pass_b_w_non_elem_nodes() );
	return (abs_diff/mean_sum);
};

float dist_mean_pass_to_cell(const TreeMutAnalyzer & Probe_mut_1,
														 const TreeMutAnalyzer & Probe_mut_2) {
	float mean_sum = (Probe_mut_1.get_mean_mut_num_cell() +
	                  Probe_mut_2.get_mean_mut_num_cell())/2;
	float abs_diff = abs( Probe_mut_1.get_mean_mut_num_cell() -
											  Probe_mut_2.get_mean_mut_num_cell() );
	return (abs_diff/mean_sum);
};


float t_dist_deg_vs_dep (const TreeMutAnalyzer & Probe_mut_1,
                         const TreeMutAnalyzer & Probe_mut_2) {
    float chi_square_dep = chi_square_depth(Probe_mut_1, Probe_mut_2);
    float chi_square_deg = chi_square_degree(Probe_mut_1, Probe_mut_2);
    return (chi_square_dep + chi_square_deg);
};

float t_dist_deg (const TreeMutAnalyzer & Probe_mut_1,
									const TreeMutAnalyzer & Probe_mut_2) {
    float chi_square_deg = chi_square_degree(Probe_mut_1, Probe_mut_2);
    return (chi_square_deg);
};

float t_dist_depth(const TreeMutAnalyzer & Probe_mut_1,
									 const TreeMutAnalyzer & Probe_mut_2) {
    float chi_square_dep = chi_square_depth(Probe_mut_1, Probe_mut_2);
    return (chi_square_dep);
};


// dist_type = 14
float t_dist_depth_cum_kolmogorov_smirnov(const TreeMutAnalyzer & Probe_mut_1,
															                  const TreeMutAnalyzer & Probe_mut_2){
  std::vector <unsigned int> depth_mut_vec_pr_1 = Probe_mut_1.get_depth_mut_vec();
  std::vector <unsigned int> depth_mut_vec_pr_2 = Probe_mut_2.get_depth_mut_vec();
  return KolmogorovSmirnovDist(depth_mut_vec_pr_1, depth_mut_vec_pr_2);

};

// dist_type = 15
float t_dist_depth_cum_kuiper(const TreeMutAnalyzer & Probe_mut_1,
																				 const TreeMutAnalyzer & Probe_mut_2){
  std::vector <unsigned int> depth_mut_vec_pr_1 = Probe_mut_1.get_depth_mut_vec();
  std::vector <unsigned int> depth_mut_vec_pr_2 = Probe_mut_2.get_depth_mut_vec();
  return KuiperDist(depth_mut_vec_pr_1, depth_mut_vec_pr_2);

};

// dist_type = 16
float t_dist_depth_cum_cramer_mises(const TreeMutAnalyzer & Probe_mut_1,
																				 const TreeMutAnalyzer & Probe_mut_2){
  std::vector <unsigned int> depth_mut_vec_pr_1 = Probe_mut_1.get_depth_mut_vec();
  std::vector <unsigned int> depth_mut_vec_pr_2 = Probe_mut_2.get_depth_mut_vec();

  return CramerMisesDist(depth_mut_vec_pr_1, depth_mut_vec_pr_2);

};

	// dist_type = 17
float t_dist_deg_cum_kolmogorov_smirnov(const TreeMutAnalyzer & Probe_mut_1,
															                  const TreeMutAnalyzer & Probe_mut_2){
  std::vector <unsigned int> deg_mut_vec_pr_1 = Probe_mut_1.get_degree_mut_vec();
  std::vector <unsigned int> deg_mut_vec_pr_2 = Probe_mut_2.get_degree_mut_vec();

  return KolmogorovSmirnovDist(deg_mut_vec_pr_1, deg_mut_vec_pr_2);

};
	// dist_type = 18
float t_dist_deg_cum_kuiper(const TreeMutAnalyzer & Probe_mut_1,
																				 const TreeMutAnalyzer & Probe_mut_2) {

  std::vector <unsigned int> deg_mut_vec_pr_1 = Probe_mut_1.get_degree_mut_vec();
  std::vector <unsigned int> deg_mut_vec_pr_2 = Probe_mut_2.get_degree_mut_vec();

  return  KuiperDist(deg_mut_vec_pr_1, deg_mut_vec_pr_2);
}

	// dist_type = 19
float t_dist_deg_cum_cramer_mises(const TreeMutAnalyzer & Probe_mut_1,
																				 const TreeMutAnalyzer & Probe_mut_2){
  std::vector <unsigned int> deg_mut_vec_pr_1 = Probe_mut_1.get_degree_mut_vec();
  std::vector <unsigned int> deg_mut_vec_pr_2 = Probe_mut_2.get_degree_mut_vec();

  return  CramerMisesDist(deg_mut_vec_pr_1, deg_mut_vec_pr_2);
};

	// dist_type = 20
float t_dist_deg_depth_sum_cum_cramer_mises(const TreeMutAnalyzer & Probe_mut_1,
																				 const TreeMutAnalyzer & Probe_mut_2){
  std::vector <unsigned int> deg_mut_vec_pr_1 = Probe_mut_1.get_degree_mut_vec();
  std::vector <unsigned int> deg_mut_vec_pr_2 = Probe_mut_2.get_degree_mut_vec();

  std::vector <unsigned int> depth_mut_vec_pr_1 = Probe_mut_1.get_depth_mut_vec();
  std::vector <unsigned int> depth_mut_vec_pr_2 = Probe_mut_2.get_depth_mut_vec();


  return  CramerMisesDist(depth_mut_vec_pr_1, depth_mut_vec_pr_2) +
					CramerMisesDist(deg_mut_vec_pr_1, deg_mut_vec_pr_2);
};



float diff_mean_vectors(std::vector <unsigned int> vec_1,
												std::vector <unsigned int> vec_2){

	float sum_1 = 0;
	for (unsigned int i = 0; i < vec_1.size(); i++) {
			sum_1 += vec_1[i] * i;
	};
	float mean_1 = sum_1/(std::accumulate(vec_1.begin(),
																			  vec_1.end(), 0) * 1.0);

	float sum_2 = 0;
	for (unsigned int i = 0; i < vec_2.size(); i++) {
			sum_2 += vec_2[i] * i;
	};
	float mean_2 = sum_2/(std::accumulate(vec_2.begin(),
																			  vec_2.end(), 0) * 1.0);

	return std::abs( mean_1 - mean_2 );

}

float t_dist_deg_mean(const TreeMutAnalyzer & Probe_mut_1,
											const TreeMutAnalyzer & Probe_mut_2){
  std::vector <unsigned int> deg_mut_vec_pr_1 = Probe_mut_1.get_degree_mut_vec();
  std::vector <unsigned int> deg_mut_vec_pr_2 = Probe_mut_2.get_degree_mut_vec();
  return diff_mean_vectors(deg_mut_vec_pr_1, deg_mut_vec_pr_2);
};


float t_dist_depth_mean(const TreeMutAnalyzer & Probe_mut_1,
											  const TreeMutAnalyzer & Probe_mut_2){
  std::vector <unsigned int> depth_mut_vec_pr_1 = Probe_mut_1.get_depth_mut_vec();
  std::vector <unsigned int> depth_mut_vec_pr_2 = Probe_mut_2.get_depth_mut_vec();

	return diff_mean_vectors(depth_mut_vec_pr_1, depth_mut_vec_pr_2);
};



float t_dist_deg_vs_mean_pass_non_elem(const TreeMutAnalyzer & Probe_mut_1,
                               const TreeMutAnalyzer & Probe_mut_2) {
	float chi_square_deg = chi_square_degree(Probe_mut_1, Probe_mut_2);
	float dist_mean = dist_mean_pass_b_w_non_elem_nodes(Probe_mut_1, Probe_mut_2);
  return (chi_square_deg + dist_mean);
};

float t_dist_deg_vs_mean_pass_to_cell (const TreeMutAnalyzer & Probe_mut_1,
                                       const TreeMutAnalyzer & Probe_mut_2) {
	float chi_square_deg = chi_square_degree(Probe_mut_1, Probe_mut_2);
	float dist_mean = dist_mean_pass_to_cell(Probe_mut_1, Probe_mut_2);
  return (chi_square_deg + dist_mean);
};

float t_dist_mean_pass_to_cell (const TreeMutAnalyzer & Probe_mut_1,
																const TreeMutAnalyzer & Probe_mut_2) {
	float dist_mean = dist_mean_pass_to_cell(Probe_mut_1, Probe_mut_2);
  return dist_mean;
};

float dist_robinson_and_foulds(const TreeMutAnalyzer & Probe_mut_1,
															 const TreeMutAnalyzer & Probe_mut_2, int seed) {
  std::string tree_str_1 = Probe_mut_1.get_str_newick_mut_tree_format();
	std::string tree_str_2 = Probe_mut_2.get_str_newick_mut_tree_format();

  return phylip_dist(tree_str_1, tree_str_2, 0, seed);
};


float dist_branch_score (const TreeMutAnalyzer & Probe_mut_1,
												 const TreeMutAnalyzer & Probe_mut_2, int seed) {
  std::string tree_str_1 = Probe_mut_1.get_str_newick_mut_tree_format();
	std::string tree_str_2 = Probe_mut_2.get_str_newick_mut_tree_format();
  return phylip_dist(tree_str_1, tree_str_2, 1, seed);
};

float nodal_distance (const TreeMutAnalyzer & Probe_mut_1,
											const TreeMutAnalyzer & Probe_mut_2) {

  std::vector <std::vector <int>> PL_matrix_1 = Probe_mut_1.GetSPLMatrix(); // Path_Lengths_matrix

  //AdjustLabels2HubVecUsingHub();
  std::vector <std::vector <int>> PL_matrix_2 = Probe_mut_2.GetSPLMatrix(); // Path_Lengths_matrix

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


float splitted_nodal_distance (const TreeMutAnalyzer & Probe_mut_1,
												       const TreeMutAnalyzer & Probe_mut_2) {
  //AdjustLabels2HubVecUsingHub();
  std::vector <std::vector <int>> PL_matrix_1 = Probe_mut_1.GetSPLMatrix(); // Splitted_Path_Lengths

  std::vector <std::vector <int>> PL_matrix_2 = Probe_mut_2.GetSPLMatrix(); // Splitted_Path_Lengths

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

float dist_splitted_nodal (const TreeMutAnalyzer & Probe_mut_1,
													 const TreeMutAnalyzer & Probe_mut_2 ) {
	std::vector <std::vector <int>> PL_matrix_1 = Probe_mut_1.GetSPLMatrix(); // Splitted_Path_Lengths
  std::vector <std::vector <int>> PL_matrix_2 = Probe_mut_2.GetSPLMatrix(); // Splitted_Path_Lengths

	unsigned int max_length_1 = 0;
	unsigned int max_length_2 = 0;

	unsigned int matrix_size = PL_matrix_1.size();

	for(unsigned int i = 0; i < matrix_size; i++) {
		for(unsigned int j = 0; j < matrix_size; j++) {
				if (max_length_1 < PL_matrix_1[i][j]) max_length_1 = PL_matrix_1[i][j];
				if (max_length_2 < PL_matrix_2[i][j]) max_length_2 = PL_matrix_2[i][j];
		};
	};

  std::vector <unsigned int> splitted_lengths_vector_1 (max_length_1 + 1, 0);
  std::vector <unsigned int> splitted_lengths_vector_2 (max_length_2 + 1, 0);

	for(unsigned int i = 0; i < matrix_size; i++) {
		for(unsigned int j = 0; j < matrix_size; j++) {
			splitted_lengths_vector_1[PL_matrix_1[i][j]]++;
			splitted_lengths_vector_2[PL_matrix_2[i][j]]++;
		};
	};

  return chi_square_sum(splitted_lengths_vector_1,
												splitted_lengths_vector_2);
};

float dist_nodal (const TreeMutAnalyzer & Probe_mut_1,
									const TreeMutAnalyzer & Probe_mut_2) {
	std::vector <std::vector <int>> PL_matrix_1 = Probe_mut_1.GetSPLMatrix(); // Splitted_Path_Lengths
  std::vector <std::vector <int>> PL_matrix_2 = Probe_mut_2.GetSPLMatrix(); // Splitted_Path_Lengths

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

  std::vector <unsigned int> splitted_lengths_vector_1 (max_length_1 + 1, 0);
  std::vector <unsigned int> splitted_lengths_vector_2 (max_length_2 + 1, 0);

	for(unsigned int i = 0; i < matrix_size; i++) {
		for(unsigned int j = i; j < matrix_size; j++) {
			splitted_lengths_vector_1[PL_matrix_1[i][j]]++;
			splitted_lengths_vector_2[PL_matrix_2[i][j]]++;
		};
	};

  return chi_square_sum(splitted_lengths_vector_1,
												splitted_lengths_vector_2);
};

float t_dist_deg_vs_dist_nodal(const TreeMutAnalyzer & Probe_mut_1,
															 const TreeMutAnalyzer & Probe_mut_2) {

    float chi_square_splitted_nodal = dist_splitted_nodal(Probe_mut_1, Probe_mut_2);
    float chi_square_deg = chi_square_degree(Probe_mut_1, Probe_mut_2);
    return (chi_square_splitted_nodal + chi_square_deg);

};

//float dist_zhao_michor(const TreeMutAnalyzer & Probe_mut_1,
//											 const TreeMutAnalyzer & Probe_mut_2) {
//
//};
//bool IsCellFittedToPiece(const int cell_index,
//                         const int piece_index,
//                         const int piece_min_x,
//                         const Area & probe_area){
//    // x min is shifted for the different probes
//    bool is_cell_x_in_range = (cells[cell_index].x >= piece_min_x) &&
//                              (cells[cell_index].x <= probe_area.max_x);
//    bool is_cell_y_in_range = (cells[cell_index].y >= probe_area.min_y) &&
//                              (cells[cell_index].y <= probe_area.max_y);
//    bool is_cell_z_in_range = (cells[cell_index].z >= probe_area.min_z) &&
//                              (cells[cell_index].z <= probe_area.max_z);
//    return (is_cell_x_in_range && is_cell_y_in_range && is_cell_z_in_range);
//};

//int GetMinxProbePiece(const BasicSimParameters & pars, const int index) {
//    //makes piece_min_x 0 or bigger
//    int del_x = pars.probe.delta_area.del_x;
//    /*probe_area.max_x - probe_area.min_x;*/
//    int del_y = pars.probe.delta_area.del_y;
//    //  projection of vector with coords (?,i*del_y,0)
//    //  shift vector (max_x,0,0) to (?, del_y, 0) and find ? = piece_min_x
//    float square_rute = pow(pow(del_x, 2.0) - pow(index*del_y, 2.0), 0.5);  /*10*/
//    int piece_min_x = del_x - int (square_rute);
//    return piece_min_x;
//};


//void InitProbePiece (const AllSimParameters & pars,
//										 ProbePiece & probe_piece_vector){
//  int piece_num = pars.inference_parameters.comparison_parameters.piece_num;
//  for (int i = 0; i < piece_num; i++) {
//      ProbePiece local_probe_piece;
//      local_probe_piece.indicator = 0;
//      local_probe_piece.cell_vector.clear();
//      local_probe_piece.piece_min_x_vec =
//				GetMinxProbePiece(pars.basic_sim_parameters, i);
//      probe_piece_vector.push_back(local_probe_piece);
//  };
//};

//void GetVecOfProbePiece(const AllSimParameters & pars,
//												vector <ProbePiece> & probe_piece_vector) {
//	int piece_number = pars.inference_parameters.comparison_parameters.piece_num;
//	int cell_number = pars.basic_sim_parameters.probe.cell_num;
//  for (unsigned int i=0; i < cells.size(); i++) {  // conditional loop
//    for (int j=0; j < piece_number; j++) {
//      if (IsCellFittedToPiece(i, j,
//                              probe_piece_vector.at(j).piece_min_x_vec,
//                              pars.basic_sim_parameters.probe.area)
//					&& (probe_piece_vector.at(j).indicator <
//          cell_number) ) {
//        probe_piece_vector.at(j).cell_vector.push_back(cells.at(i));
//        probe_piece_vector.at(j).indicator++;
//        break;
//      };
//    };
//
//    bool is_search_finished = true;
//		for (int j=0; j < piece_number; j++) {
//      if (probe_piece_vector.at(j).indicator <
//            cell_number) is_search_finished = false;
//    };
//    if (is_search_finished) break;
//  };
//};

//void GetVecOfRandProbePiece(const AllSimParameters & pars,
//														vector <ProbePiece> & probe_piece_vector) {
//	int piece_number = pars.inference_parameters.comparison_parameters.piece_num;
//  int cell_number = pars.basic_sim_parameters.probe.cell_num;
//  for (int i = 0; i < piece_number; i++) {
//		for (int j = 0; j < cell_number; j++) {
//			unsigned int index = trunc((cells.size()-1)*_drand48());
//			probe_piece_vector.at(i).cell_vector.push_back(cells.at(index));
//		};
//  };
//};


void GetComparisonResultsFromFile(const AllSimParameters & pars,
															    std::vector<ComparisonResults> & comparison_results){
  std::string full_comparison_file_name =
				  pars.basic_sim_parameters.related_path_to_par_file
					+ pars.basic_sim_parameters.etalon_result_catalogue
					+ getSystemRelevantSlash()
					+ pars.inference_parameters.basic_settings.current_result_file_name;
	std::ifstream comparison_file;

  comparison_file.open(full_comparison_file_name.c_str());
	if (comparison_file.is_open()) {
		do {
			ComparisonResults current;
			comparison_file >> current;
			if (!comparison_file.eof()) comparison_results.push_back(current);
		} while (!comparison_file.eof());
	};
};

std::vector <ComparisonResults> GetAvarageVec(const AllSimParameters & pars,
															           const std::vector<ComparisonResults> & comparison_results){
  std::vector <ComparisonResults> averaging_vec;
//
  ComparisonResults averaging = comparison_results[0];
  int number_of_elements = 1;
  //bool is_new = true;
  for (int i=0; i < (comparison_results.size() - 1); i++) {
	  if (comparison_results[i].evolution == comparison_results[i+1].evolution) {
				averaging = averaging + comparison_results[i+1];
//				std::cout << averaging.evolution.driver_adv <<";"<< averaging.evolution.mutation_rate<<";"<<
//				        averaging.evolution.driver_mutation_rate<<std::endl;
				number_of_elements++;

	  } else {
//				int ze;
//	  	  std::cin>>ze;
//	  	  std::cout << "HI"<<std::endl;
//				std::cout << averaging.evolution.driver_adv <<";"<< averaging.evolution.mutation_rate<<";"<<
//				        averaging.evolution.driver_mutation_rate<<std::endl;


	  	  averaging = averaging/number_of_elements;
//				std::cout << averaging.evolution.driver_adv <<";"<< averaging.evolution.mutation_rate<<";"<<
//				        averaging.evolution.driver_mutation_rate<<";";
//				std::cout << averaging.balk_dist <<";"<< averaging.tree_dist<<";"<<
//				        averaging.dist<<std::endl;

	  	  averaging_vec.push_back(averaging);
				averaging.balk_dist = comparison_results[i+1].balk_dist;
				averaging.tree_dist = comparison_results[i+1].tree_dist;
				averaging.dist = comparison_results[i+1].dist;
				averaging.evolution.mutation_rate = comparison_results[i+1].evolution.mutation_rate;
				averaging.evolution.driver_mutation_rate = comparison_results[i+1].evolution.driver_mutation_rate;
				averaging.evolution.driver_adv = comparison_results[i+1].evolution.driver_adv;
				number_of_elements = 1;
		};

		if ((i+1) == comparison_results.size() - 1) {
				averaging = averaging/number_of_elements;
				averaging_vec.push_back(averaging);
		};
  };
  return averaging_vec;
};

void SetAvarageComparisonVectorToFile(const AllSimParameters & pars,
															        const std::vector<ComparisonResults> & comparison_results){
  std::string full_comparison_file_name =
				  pars.basic_sim_parameters.related_path_to_par_file
					+ pars.basic_sim_parameters.etalon_result_catalogue
					+ getSystemRelevantSlash()
					+ "avarage" + pars.inference_parameters.basic_settings.current_result_file_name ;

	std::ofstream comparison_file;
  comparison_file.open(full_comparison_file_name.c_str());
	comparison_file.precision(5);
	comparison_file.setf(std::ios::fixed, std:: ios::floatfield);

  for(ComparisonResults current:comparison_results){
      comparison_file << current.evolution.driver_adv <<"\t"<<
												 current.evolution.mutation_rate<<"\t"<<
												 current.evolution.driver_mutation_rate        <<"\t"<<
        current.balk_dist                      <<"\t"<<
        current.tree_dist                      <<"\t"<<
        current.dist                           <<"\t"<< std::endl;
  };
  comparison_file.close();
};

void MakeAverageFile(const AllSimParameters & pars){
				std::vector <ComparisonResults> comparison_results;
				GetComparisonResultsFromFile(pars,comparison_results);
	      std::vector <ComparisonResults> average_vec = GetAvarageVec(pars,comparison_results);
				SetAvarageComparisonVectorToFile(pars,average_vec);
};

void SetComparisonResults2File(const AllSimParameters & pars,
															 const ProbePiece & probe_piece/*
                               const vector <ProbePiece> & probe_piece_vector*/) {
//opens file_tree_comparison.dat
  std::ofstream comparison_file;

  std::string full_comparison_file_name =
				  pars.basic_sim_parameters.related_path_to_par_file
					+ pars.basic_sim_parameters.etalon_result_catalogue
					+ getSystemRelevantSlash()
					+ pars.inference_parameters.basic_settings.current_result_file_name;

  comparison_file.open(full_comparison_file_name.c_str(), std::ios::app);
	//int piece_number = pars.inference_parameters.comparison_parameters.piece_num;
  //for (int i=0; i < piece_number; i++) {

    if (probe_piece.results.dist < 10000) {
      comparison_file.precision(5);
      comparison_file.setf(std::ios::fixed, std:: ios::floatfield);
      comparison_file << pars.basic_sim_parameters.evolution.driver_adv <<"\t"<<
        pars.basic_sim_parameters.evolution.mutation_rate               <<"\t"<<
        pars.basic_sim_parameters.evolution.driver_mutation_rate        <<"\t"<<
        probe_piece.results.balk_dist                      <<"\t"<<
        probe_piece.results.tree_dist                      <<"\t"<<
        probe_piece.results.dist                           <<"\t"<< std::endl;
        //probe_piece.results.tree_delta_mut                 <<"\t"<<
        //probe_piece.results.tree_num_mut                   <<

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
//  };
  comparison_file.close();
};




float CoefficientProbePieceCompatible2EProbe(const AllSimParameters  & pars,
																	 const TreeMutAnalyzer & probe,
																	 const TreeMutAnalyzer & e_probe) {
  int probe_num_cells = probe.get_cell_number();
  int e_probe_num_cells = e_probe.get_cell_number();
  int delta_size_cells = abs(probe_num_cells - e_probe_num_cells);
  float relative_delta_cell = (delta_size_cells * 1.0)/(e_probe_num_cells*1.0);
  float cell_number_factor;
  //if (relative_delta_cell  < pars.inference_parameters.comparison_parameters.threshold_frac_cell_diff) cell_number_factor = 1;
	/*else*/ cell_number_factor = pow(relative_delta_cell + 1, 2.0);

	int probe_mut_num     =  probe.get_mut_number();
	int e_probe_mut_num   =  e_probe.get_mut_number();
	int delta_probe_e_mut =  abs(e_probe_mut_num - probe_mut_num);
	float relative_delta_mut = (delta_probe_e_mut * 1.0)/(e_probe_mut_num*1.0);

  float mut_number_factor;
  //if (relative_delta_mut < pars.inference_parameters.comparison_parameters.threshold_frac_mut_diff) mut_number_factor = 1;
	/*else*/ mut_number_factor = relative_delta_mut + 1;

	return cell_number_factor * mut_number_factor;
};
void GetTreeComparisonResults( const AllSimParameters & pars,
															 const TreeMutAnalyzer & probe,
                               const TreeMutAnalyzer & e_probe,
                               ComparisonResults & c_results) {
	int seed = pars.basic_sim_parameters.seed;
  switch (pars.inference_parameters.comparison_parameters.tree_dist_type) {
		case 0: {
			c_results.tree_dist = t_dist_deg_vs_dep(probe, e_probe);
		 } break;
		case 1: {
			c_results.tree_dist = t_dist_deg_vs_mean_pass_non_elem(probe, e_probe);

		} break;
		case 2: {
			c_results.tree_dist = t_dist_deg_vs_mean_pass_to_cell (probe, e_probe);

		} break;
		case 3: {
			c_results.tree_dist = t_dist_mean_pass_to_cell (probe, e_probe);

		} break;
		case 4: {
			c_results.tree_dist = dist_robinson_and_foulds (probe, e_probe, seed);

		} break;
		case 5: {
			c_results.tree_dist = dist_branch_score (probe, e_probe, seed);

		} break;
		case 6: {
			c_results.tree_dist = nodal_distance (probe, e_probe);

		} break;
		case 7: {
			c_results.tree_dist = splitted_nodal_distance (probe, e_probe);

		} break;
		case 8: {
			c_results.tree_dist = dist_nodal (probe, e_probe);

		} break;
		case 9: {
			c_results.tree_dist = dist_splitted_nodal (probe, e_probe);

		} break;
		case 10: {
			c_results.tree_dist = t_dist_deg_vs_dist_nodal(probe, e_probe);

		} break;
		case 11: {
			c_results.tree_dist = t_dist_deg(probe, e_probe);

		} break;
		case 12: {
			c_results.tree_dist = CoefficientProbePieceCompatible2EProbe(pars,probe, e_probe);

		} break;

		case 13: {
			c_results.tree_dist = t_dist_depth(probe, e_probe);
		} break;

		case 14: {
			c_results.tree_dist = t_dist_depth_cum_kolmogorov_smirnov(probe, e_probe);
		} break;

		case 15: {
			c_results.tree_dist = t_dist_depth_cum_kuiper(probe, e_probe);
		} break;

		case 16: {
			c_results.tree_dist = t_dist_depth_cum_cramer_mises(probe, e_probe);
		} break;

		case 17: {
			c_results.tree_dist = t_dist_deg_cum_kolmogorov_smirnov(probe, e_probe);
		} break;

		case 18: {
			c_results.tree_dist = t_dist_deg_cum_kuiper(probe, e_probe);
		} break;

		case 19: {
			c_results.tree_dist = t_dist_deg_cum_cramer_mises(probe, e_probe);
		} break;

		case 20: {
			c_results.tree_dist = t_dist_deg_depth_sum_cum_cramer_mises(probe, e_probe);
		} break;

		case 21: {
			c_results.tree_dist = t_dist_deg_mean(probe, e_probe);
		} break;


		case 22: {
			c_results.tree_dist = t_dist_depth_mean(probe, e_probe);
		} break;


     default:
     err("Mistake in distance type (see setting file)");
   };
//  c_results.tree_dist *= CoefficientProbePieceCompatible2EProbe(pars,probe, e_probe);


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

//  c_results.tree_e_num_mut = e_probe.get_mut_number();
//  c_results.tree_num_mut   = probe.get_mut_number();
//  c_results.tree_delta_mut = abs(c_results.tree_e_num_mut - c_results.tree_num_mut);

  //cout<<"tree dist is "                << c_results.dist      << endl;
  //cout<<"Number of etalon mutations "  << c_results.e_num_mut << endl;
  //cout<<"Number of probe mutations "   << c_results.num_mut   << endl;
  //cout<<"Difference in mut number is " << c_results.delta_mut << endl;

//  return c_results;
};



void GetBulkComparisonResults(const AllSimParameters & pars,
															const BalkVAFAnalyzer & probe,
															const BalkVAFAnalyzer & e_probe,
															ComparisonResults & c_results){
	switch (pars.inference_parameters.comparison_parameters.balk_dist_type){
		case 0: {
			c_results.balk_dist = KolmogorovSmirnovDist(probe, e_probe);
			break;
		}
		case 1: {
			c_results.balk_dist = KuiperDist(probe, e_probe);
		  break;
		}

		case 2: {
			c_results.balk_dist = KuiperDist(probe, e_probe);
		  break;
		}

		default:
		err("Mistake in distance type (see setting file)");
	}

};

void GetCellMutComparisonResults( const AllSimParameters & pars,
															    CellMutAnalyzer & probe,
                                  CellMutAnalyzer & e_probe,
                                  ComparisonResults & c_results) {
	switch (pars.inference_parameters.comparison_parameters.balk_dist_type){
		case 100:{
			//probe.MutationFilter(2,6);
			//e_probe.MutationFilter(2,6);
			c_results.balk_dist = ZhaoMichorDist(probe, e_probe, 2, round(e_probe.GetCellNum()*0.15));
			break;
		}
		default:
		err("Mistake in distance type (see setting file)");
	}


};

void ProcessCellMutProbePiece(const AllSimParameters & pars,
										ProbePiece & probe_piece) {
		std::string path_tree(pars.basic_sim_parameters.related_path_to_par_file +
									 pars.basic_sim_parameters.etalon_result_catalogue +
									 getSystemRelevantSlash() + "cell_mut" +
									 pars.basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name);
		std::ifstream file_tree_mut_vec_stream;
		file_tree_mut_vec_stream.open(path_tree.c_str());
		if (file_tree_mut_vec_stream.is_open()) {
			CellMutAnalyzer e_probe;  // etalon probe
			file_tree_mut_vec_stream >> e_probe;
			CellMutAnalyzer probe(probe_piece.tree_cell_vector, genotypes);

//      vector <unsigned int> A = probe.GetMutCellFr();
//      vector <unsigned int> B = probe.GetCellMutFr();
//      for(unsigned int a:A) std::cout<<a<<";";
//      std::cout <<std::endl;
//      for(unsigned int b:B) std::cout<<b<<";";
//      std::cout <<std::endl;
//
//
//      vector <unsigned int> A1 = e_probe.GetMutCellFr();
//      vector <unsigned int> B1 = e_probe.GetCellMutFr();
//      for(unsigned int a:A1) std::cout<<a<<";";
//      std::cout <<std::endl;
//      for(unsigned int b:B1) std::cout<<b<<";";
//      std::cout <<std::endl;


			GetCellMutComparisonResults(pars, probe, e_probe, probe_piece.results);
//			std::cout<<"probe_piece.results.balk_dist="<<probe_piece.results.balk_dist<<"\n";
		} else {
			probe_piece.results.balk_dist = -10;
		}
		file_tree_mut_vec_stream.close();
};

void ProcessProbePiece(const AllSimParameters & pars,
											 ProbePiece & probe_piece) {
  if (pars.basic_sim_parameters.probes.tree_cells.use) {
		std::string path_tree(pars.basic_sim_parameters.related_path_to_par_file +
									 pars.basic_sim_parameters.etalon_result_catalogue +
									 getSystemRelevantSlash() +
									 pars.basic_sim_parameters.probes.tree_cells.etalon_result_t_c_file_name);
		std::ifstream file_tree_mut_vec_stream;
		file_tree_mut_vec_stream.open(path_tree.c_str());
		if (file_tree_mut_vec_stream.is_open()) {
			TreeMutAnalyzer e_probe;  // etalon probe
			file_tree_mut_vec_stream >> e_probe;
			//int piece_number = pars.inference_parameters.comparison_parameters.piece_num;
			//for (int i=0; i < piece_number; i++) {
			TreeMutAnalyzer probe(probe_piece.tree_cell_vector);
			//std::cout<<"piece_num = " << i << std::endl;
			//std::cout<<"mean mut number of cells in e_probe = " <<
			//	e_probe.get_mean_mut_num_cell() << std::endl;
			//std::cout<<"get_mean_pass_b_w_non_elem_nodes for e_probe = " <<
			//	e_probe.get_mean_pass_b_w_non_elem_nodes()<<std::endl;
			//std::cout<<"mean mut number of cells in piece_probe = " <<
			//	probe.get_mean_mut_num_cell() << std::endl;
			//std::cout<<"get_mean_pass_b_w_non_elem_nodes for piece_probe = " <<
			//	probe.get_mean_pass_b_w_non_elem_nodes()<<std::endl;
			//bool is_compatible =
      //  IsProbePieceCompatible2EProbe(pars, probe, e_probe);
			//if (is_compatible){
        //probe_piece_vector.at(i).results = GetComparisonResults(pars, probe, e_probe);
			GetTreeComparisonResults(pars, probe, e_probe, probe_piece.results);
			//};
		};
		file_tree_mut_vec_stream.close();
  } else {
		probe_piece.results.tree_dist = 0;
  };

  if (pars.basic_sim_parameters.probes.balk_cells.use) {
		std::string path_balk(pars.basic_sim_parameters.related_path_to_par_file +
									 pars.basic_sim_parameters.etalon_result_catalogue +
									 getSystemRelevantSlash() +
									 pars.basic_sim_parameters.probes.balk_cells.etalon_result_b_c_file_name);
		std::ifstream file_balk_vaf_stream;
		file_balk_vaf_stream.open(path_balk.c_str());
		if (file_balk_vaf_stream.is_open()) {
			BalkVAFAnalyzer e_probe;
			file_balk_vaf_stream >> e_probe;
			BalkVAFAnalyzer probe(probe_piece.balk_cell_vector, genotypes);
			if (pars.inference_parameters.comparison_parameters.balk_dist_type < 100) {
				GetBulkComparisonResults(pars, probe, e_probe, probe_piece.results);
			} else {
				ProcessCellMutProbePiece(pars, probe_piece);
			};
		};
  } else {
		probe_piece.results.balk_dist = 0;
  };
  probe_piece.results.dist = probe_piece.results.tree_dist + probe_piece.results.balk_dist;
  SetComparisonResults2File(pars, probe_piece/*probe_piece_vector*/);
};




// can be generalized with number of div and transfer points in area by Rotation of all points
void GetBalkCellVector(const BasicSimParameters & pars,
											 std::vector <Cell> & balk_cell_vector) {
  unsigned int squared_needle_diam = pars.probes.balk_cells.needle_diam *
																		 pars.probes.balk_cells.needle_diam;
  //for_loop: Find fitted cells
  for (unsigned int i=0; i < cells.size(); i++) {
		// cell in cylinder of diameter pars.probes.balk_cells.needle_diam from max.x to 0
		unsigned int squared_cell_y = cells[i].y*cells[i].y;
		unsigned int squared_cell_z = cells[i].z*cells[i].z;
		if( ((cells[i].x >= 0) && (cells[i].x <= pars.probes.population_borders.max_x))/*x direction*/ &&
				((squared_cell_y + squared_cell_z) * 4 <= squared_needle_diam) /*y and z direction*/
			)
		{
			 balk_cell_vector.push_back(cells[i]);
    };
  };
};


void GetTreeCellVectorFromBalkVector(const BasicSimParameters & pars,
														         const std::vector <Cell> & balk_cell_vector,
																		 std::vector <Cell> & tree_cell_vector) {
  //for_loop: Find fitted cells (number less then threshold_num_cells)
	for (int i=0; i < pars.probes.tree_cells.cell_num; i++) {
		int index = trunc((balk_cell_vector.size()-1)*_drand48());
		tree_cell_vector.push_back(balk_cell_vector.at(index));
	};
  //begin delete this
  //cout<<"number of sells in simulation "<<cells.size()<<endl;
  //end delete this
};

void GetRandomTreeCellVector(const BasicSimParameters & pars,
														 std::vector <Cell> & tree_cell_vector) {
	for (int i=0; i < pars.probes.tree_cells.cell_num; i++) {
		int index = trunc((cells.size()-1)*_drand48());
		tree_cell_vector.push_back(cells.at(index));
	};
};

void SmallTest(TreeMutAnalyzer & probe_mut_processor,
                int size_vector_cells) {
    int num_of_cells_e_prob = probe_mut_processor.get_cell_number();
    //cout<<"num_of_cells_e_prob="<<num_of_cells_e_prob<<endl;
    //cout<<"size_vector_cells="<<size_vector_cells<<endl;
    assert(size_vector_cells == num_of_cells_e_prob);
};

//Inits probe area. Takes cuboid with del_x, del_y, del_z sides in
// x,y,z direction. IN X DIRECTION IT GOES FROM BORDER SIDE MAX_X
//void InitProbeArea(BasicSimParameters & pars) {
//    pars.probe.area.min_x = pars.probe.population_borders.max_x -
//                                     pars.probe.delta_area.del_x;
//    pars.probe.area.max_x = pars.probe.population_borders.max_x;
//    pars.probe.area.min_y = - int (pars.probe.delta_area.del_y*1.0/2);
//    pars.probe.area.max_y = int (pars.probe.delta_area.del_y*1.0/2);
//    pars.probe.area.min_z = - int (pars.probe.delta_area.del_z*1.0/2);
//    pars.probe.area.max_z = int (pars.probe.delta_area.del_z*1.0/2);
//};

void SetBalkResultsToFile(const BasicSimParameters & pars,
													const BalkVAFAnalyzer & balk_vaf_analyzer) {
  std::string path = pars.related_path_to_par_file +
								pars.etalon_result_catalogue +
								getSystemRelevantSlash() +
								pars.probes.balk_cells.etalon_result_b_c_file_name;
	std::ofstream etalon_balk_results;
  etalon_balk_results.open (path.c_str());
  etalon_balk_results << balk_vaf_analyzer;
};

void SetTreeCellResultsToFile(const BasicSimParameters & pars,
													    TreeMutAnalyzer & probe_mut_processor) {
  std::string path = pars.related_path_to_par_file +
								pars.etalon_result_catalogue +
								getSystemRelevantSlash() +
								pars.probes.tree_cells.etalon_result_t_c_file_name;
  std::ofstream etalon_file_mut_vec;
  etalon_file_mut_vec.open (path.c_str());
  //cout<<file_name_part<<"/etalon_mut_vec.dat is opened"<<endl;
  etalon_file_mut_vec << probe_mut_processor;
  etalon_file_mut_vec.close();

};

void SetCellMutResultsToFile(const BasicSimParameters & pars,
														 CellMutAnalyzer & cell_mut_analyzer) {
  std::string path = pars.related_path_to_par_file +
				 pars.etalon_result_catalogue +
				 getSystemRelevantSlash() +
				 "cell_mut" +
				 pars.probes.tree_cells.etalon_result_t_c_file_name;
  std::ofstream etalon_file_mut_vec;
  etalon_file_mut_vec.open (path.c_str());
  etalon_file_mut_vec << cell_mut_analyzer;
  etalon_file_mut_vec.close();
};

void DoesInferenceAnalyze(AllSimParameters & pars) {
  //vector <ProbePiece> probe_piece_vector;
  //InitProbeArea(pars.basic_sim_parameters);
  //InitProbePieceVector(pars, probe_piece_vector);
  ProbePiece probe_piece;
  probe_piece.balk_cell_vector.clear();
  GetBalkCellVector(pars.basic_sim_parameters, probe_piece.balk_cell_vector); // can be generalized

  std::vector <Cell> tree_cell_vector;
	tree_cell_vector.clear();
  if (pars.basic_sim_parameters.probes.tree_cells.use){
		// init tree_cell_vector
		if (pars.basic_sim_parameters.probes.tree_cells.in_balk == false) {
			GetRandomTreeCellVector(pars.basic_sim_parameters, probe_piece.tree_cell_vector);
		} else {
			GetTreeCellVectorFromBalkVector(pars.basic_sim_parameters,
																			probe_piece.balk_cell_vector,
																			probe_piece.tree_cell_vector);
		};
  };
  ProcessProbePiece(pars, probe_piece);
//  if (pars.basic_sim_parameters.probe.is_random) {
//		GetVecOfRandProbePiece(pars, probe_piece_vector);
//  } else {
//    GetVecOfProbePiece(pars, probe_piece_vector);
//  };
};

void SetEtalonProbe(BasicSimParameters & pars) {
  std::vector <Cell> balk_cell_vector;
  balk_cell_vector.clear();
	//InitProbeArea(pars);
	GetBalkCellVector(pars, balk_cell_vector);

  std::vector <Cell> tree_cell_vector;
	tree_cell_vector.clear();
  if (pars.probes.tree_cells.use){
		// init tree_cell_vector
		if (pars.probes.tree_cells.in_balk == false) {
			GetRandomTreeCellVector(pars, tree_cell_vector);
		} else {
			GetTreeCellVectorFromBalkVector(pars, balk_cell_vector, tree_cell_vector);
		};

  };

  if (pars.probes.tree_cells.use && (balk_cell_vector.size() > 0)) {
		BalkVAFAnalyzer balk_vaf_analyzer(balk_cell_vector, genotypes);
		SetBalkResultsToFile(pars, balk_vaf_analyzer);
  };

  if (tree_cell_vector.size() > 0) {
    TreeMutAnalyzer tree_mut_analyzer(tree_cell_vector);
    CellMutAnalyzer cell_mut_analyzer(tree_cell_vector, genotypes);

    //SmallTest(tree_mut_analyzer, tree_cell_vector.size());//!!!
    SetTreeCellResultsToFile(pars, tree_mut_analyzer);
    SetCellMutResultsToFile(pars, cell_mut_analyzer);
  } else {
    //cout<<"there is no cells in the etalon probe"<< endl;
  };
};
BalkVAFAnalyzer::BalkVAFAnalyzer(){
		num_of_cells = 0;
		num_of_muts = 0;
    num_mutations_vec.clear();
    cumulative_vaf_vec.clear();
};

BalkVAFAnalyzer::BalkVAFAnalyzer(std::vector <Cell> & cell_vec,
												         std::vector <Genotype*> & cell_genotypes){
	num_of_cells = cell_vec.size();
  std::set<unsigned int> cell_gene_index_set;
  //inserts cell index of gen to set of gen indexes
  for(Cell cell: cell_vec) cell_gene_index_set.insert(cell.gen);
  //goes through set of gen indexes and gather mutations in set of muts
	std::set<unsigned int> cell_gen_mut_set;
	for(unsigned int cell_gene_index: cell_gene_index_set){
	  std::vector <unsigned int> sequence = cell_genotypes.at(cell_gene_index)->sequence;
		cell_gen_mut_set.insert(sequence.begin(),sequence.end());
	};
	num_of_muts = cell_gen_mut_set.size();
  //gathers number of cells with correspondent mut
  std::vector <unsigned int> cell_gen_mut_num_vec;
	for(unsigned int cell_mut: cell_gen_mut_set){
		unsigned int number = 0;
		for(Cell cell: cell_vec){
				std::vector <unsigned int> sequence = cell_genotypes.at(cell.gen)->sequence;
				for (unsigned int mut:sequence) if (mut == cell_mut) number++;
		};
		cell_gen_mut_num_vec.push_back(number);
	};
  //sorting vector of number of cell with some mut
  std::sort(cell_gen_mut_num_vec.begin(), cell_gen_mut_num_vec.end(),
						[](const unsigned int &left,const unsigned int &right){return left < right;});
	num_mutations_vec.resize(num_of_cells, 0);
	for(auto cell_gen_mut: cell_gen_mut_num_vec) num_mutations_vec.at(cell_gen_mut-1)++;
	InitCumulativeVAFVec();
};

void BalkVAFAnalyzer::InitCumulativeVAFVec(){
	float cumulative_num_mutations = 0;
	//creates cumulative_vaf_vec
	for (auto num_mutations: num_mutations_vec) {
			cumulative_num_mutations += (1.0 * num_mutations)/(1.0 * num_of_muts);
			cumulative_vaf_vec.push_back(cumulative_num_mutations);
	};
};

void BalkVAFAnalyzer::PrintToFileCumulativeDist(std::ostream & output){
		float frec = 0;
    for(auto x:cumulative_vaf_vec) {
				frec += 1.0 / (1.0 * num_of_cells);
				output << frec <<"\t"<< x << std::endl;
    };
};


float KuiperDist(const BalkVAFAnalyzer & balk_vaf_1,
								 const BalkVAFAnalyzer & balk_vaf_2){
	unsigned int step_num_1 = balk_vaf_1.get_cell_num();
	unsigned int step_num_2 = balk_vaf_2.get_cell_num();
	// the length of step
	float delta_1 = 1.0 / (1.0 * step_num_1);
	float delta_2 = 1.0 / (1.0 * step_num_2);
	// proportions of lengths
	float fraction_1 = delta_1/delta_2;
	float fraction_2 = delta_2/delta_1;
  // vectors of cumulative mut frequencies (distribution of muts in cells)
  std::vector<float> cum_vaf_vec_1 = balk_vaf_1.get_cumulative_vaf_vec();
	std::vector<float> cum_vaf_vec_2 = balk_vaf_2.get_cumulative_vaf_vec();

	float max_up = 0;
	float max_down = 0;
	for (unsigned int i = 0; i < step_num_1; i ++){
			float current = cum_vaf_vec_1[i] - cum_vaf_vec_2 [trunc(i*fraction_1)];
			if  (max_up < current) max_up = current;
			if  (max_down > current) max_down = current;
 	};

	for (unsigned int i = 0; i < step_num_2; i ++){
			float current = cum_vaf_vec_1 [trunc(i*fraction_2)] - cum_vaf_vec_2[i];
			if  (max_up < current) max_up = current;
			if  (max_down > current) max_down = current;
	};

	float max_diff = max_up - max_down;
//	std::cout << "max_up = " << max_down << "max_down = " << max_down <<std::endl;
//	std::cout <<"n="<<step_num_1 <<";m="<<step_num_2<<"; " <<sqrt((1.0*step_num_1*step_num_2)/(step_num_1*1.0 + step_num_2*1.0)) <<"; " <<max_diff << std::endl;
//	std::cout << "mut_num_1=" << balk_vaf_1.get_mut_num()
//						<< "; mut_num_2=" << balk_vaf_2.get_mut_num()<<std::endl;
	float multiplier = sqrt((1.0*step_num_1*step_num_2)/(step_num_1*1.0 + step_num_2*1.0));
	return multiplier*max_diff;
};


float KolmogorovSmirnovDist(const BalkVAFAnalyzer & balk_vaf_1,
												    const BalkVAFAnalyzer & balk_vaf_2){
	unsigned int step_num_1 = balk_vaf_1.get_cell_num();
	unsigned int step_num_2 = balk_vaf_2.get_cell_num();
	// the length of step
	float delta_1 = 1.0 / (1.0 * step_num_1);
	float delta_2 = 1.0 / (1.0 * step_num_2);
	// proportions of lengths
	float fraction_1 = delta_1/delta_2;
	float fraction_2 = delta_2/delta_1;
  // vectors of cumulative mut frequencies (distribution of muts in cells)
  std::vector<float> cum_vaf_vec_1 = balk_vaf_1.get_cumulative_vaf_vec();
	std::vector<float> cum_vaf_vec_2 = balk_vaf_2.get_cumulative_vaf_vec();

	float max_diff = 0;
	for (unsigned int i = 0; i < step_num_1; i ++){
			float current = std::abs(cum_vaf_vec_1[i] - cum_vaf_vec_2 [trunc(i*fraction_1)]);
			if  (max_diff < current) max_diff = current;
	};

	for (unsigned int i = 0; i < step_num_2; i ++){
			float current = std::abs(cum_vaf_vec_2[i] - cum_vaf_vec_1 [trunc(i*fraction_2)]);
			if  (max_diff < current) max_diff = current;
	};
//	std::cout <<"n="<<step_num_1 <<";m="<<step_num_2<<"; " <<sqrt((1.0*step_num_1*step_num_2)/(step_num_1*1.0 + step_num_2*1.0)) <<"; " <<max_diff << std::endl;
//	std::cout << "mut_num_1=" << balk_vaf_1.get_mut_num()
//						<< "; mut_num_2=" << balk_vaf_2.get_mut_num()<<std::endl;
  float multiplier = sqrt((1.0*step_num_1*step_num_2)/(step_num_1*1.0 + step_num_2*1.0));
	return multiplier*max_diff;
};




float CramerMisesDist(const BalkVAFAnalyzer & balk_vaf_1,
									    const BalkVAFAnalyzer & balk_vaf_2) {
	unsigned int step_num_1 = balk_vaf_1.get_cell_num();
	unsigned int step_num_2 = balk_vaf_2.get_cell_num();
	// the length of step
	float delta_1 = 1.0 / (1.0 * step_num_1);
	float delta_2 = 1.0 / (1.0 * step_num_2);
	// proportions of lengths
	float fraction_1 = delta_1/delta_2;
	float fraction_2 = delta_2/delta_1;
  // vectors of cumulative mut frequencies (distribution of muts in cells)
  std::vector<float> cum_vaf_vec_1 = balk_vaf_1.get_cumulative_vaf_vec();
	std::vector<float> cum_vaf_vec_2 = balk_vaf_2.get_cumulative_vaf_vec();

	float sum_sq_diff_1 = 0;
	for (unsigned int i = 0; i < step_num_1; i ++){
			float diff = std::abs(cum_vaf_vec_1[i] - cum_vaf_vec_2[trunc(i*fraction_1)]);
			float sq_diff = diff * diff;
			sum_sq_diff_1 += sq_diff;
	};

	float sum_sq_diff_2 = 0;
	for (unsigned int i = 0; i < step_num_2; i ++){
			float diff = std::abs(cum_vaf_vec_2[i] - cum_vaf_vec_1[trunc(i*fraction_2)]);
			float sq_diff = diff * diff;
			sum_sq_diff_2 += sq_diff;
	};

	float sq_sum = (step_num_1*1.0 + step_num_2*1.0) * (step_num_1*1.0 + step_num_2*1.0);
	float multiplier = (1.0*step_num_1*step_num_2)/sq_sum;
	return multiplier*(sum_sq_diff_1 + sum_sq_diff_2);
};

std::vector <unsigned int> CellMutAnalyzer::GetCellMutFr(){
	//vector <unsigned int> CellFr
	std::vector <unsigned int> cell_gen_mut_num_vec = GetMutCellFr();
  //std::cout<<"size"<< cell_gen_mut_num_vec.size()<<"\n";
 	std::vector <unsigned int> num_mutations_vec;
 	num_mutations_vec.resize(cell_vec.size() + 1, 0);
 	//std:cout<<"num_mutations_vec.size()="<<num_mutations_vec.size()<<std::endl;
	for(auto cell_gen_mut: cell_gen_mut_num_vec) num_mutations_vec.at(cell_gen_mut)++;
	return num_mutations_vec;
};

std::vector <unsigned int> CellMutAnalyzer::GetMutCellFr(){
  //goes through set of gen indexes and gather mutations in set of muts
	std::set<unsigned int> gen_mut_set;

	for(Genotype current_genotype: genotype_vec){
	  std::vector <unsigned int> sequence = current_genotype.sequence;
		gen_mut_set.insert(sequence.begin(),sequence.end());
	};

  //gathers number of cells with correspondent mut
  std::vector <unsigned int> cell_gen_mut_num_vec;

	for (unsigned int current_gen_mut: gen_mut_set){
		unsigned int cell_num_of_current_gen_mut = 0;

		//calculate number of cells with current_gen_mut
		for (Cell current_cell: cell_vec) {
			//find element of genotype_vec with current_cell_index
			std::vector <unsigned int> genotype_mut_seq;
			for (Genotype current_genotype : genotype_vec)
				if (current_genotype.index == current_cell.gen) {
							genotype_mut_seq =current_genotype.sequence;
							break;
				};
			unsigned int is_in_cell_genotype = 0;
				for (unsigned int current_mut:genotype_mut_seq)
					if (current_mut == current_gen_mut) {is_in_cell_genotype = 1; break;}
			cell_num_of_current_gen_mut += is_in_cell_genotype;
		};
		cell_gen_mut_num_vec.push_back(cell_num_of_current_gen_mut);
	};

  //sorting vector of number of cell with some mut
  std::sort(cell_gen_mut_num_vec.begin(), cell_gen_mut_num_vec.end(),
						[](const unsigned int &left,const unsigned int &right){return left < right;});
	return cell_gen_mut_num_vec;
};

void CellMutAnalyzer::MutationFilter(const unsigned int left_border,
										                 const unsigned int right_border){
  //goes through set of gen indexes and gather mutations in set of muts
	std::set<unsigned int> gen_mut_set;

	for(Genotype current_genotype: genotype_vec){
	  std::vector <unsigned int> sequence = current_genotype.sequence;
		gen_mut_set.insert(sequence.begin(),sequence.end());
	};

//  unsigned int inc = 0;
	for (unsigned int current_gen_mut: gen_mut_set){
			unsigned int cell_num_of_current_gen_mut = 0;
			//calculate number of cells with current_gen_mut
			for (Cell current_cell: cell_vec) {
			   //find element of genotype_vec with current_cell_index
			   std::vector <unsigned int> genotype_mut_seq;
			   for (Genotype current_genotype : genotype_vec)
					 if (current_genotype.index == current_cell.gen) {
							genotype_mut_seq = current_genotype.sequence;
							break;
					 };
			   unsigned int is_in_cell_genotype = 0;
				 for (unsigned int current_mut:genotype_mut_seq)
					 if (current_mut == current_gen_mut) {is_in_cell_genotype = 1; break;};
         cell_num_of_current_gen_mut += is_in_cell_genotype;
			};

			if ( (cell_num_of_current_gen_mut > right_border) ||
					 (cell_num_of_current_gen_mut < left_border ) ) {
					for (Genotype & current_genotype : genotype_vec)
						for (unsigned i = 0; i < current_genotype.sequence.size(); i++)
						  if (current_genotype.sequence.at(i) == current_gen_mut)
						    current_genotype.sequence.erase(current_genotype.sequence.begin() + i);
						//std::cout<<inc<<" - " << cell_num_of_current_gen_mut <<"\n";
			};
	};
};

float ZhaoMichorDist(CellMutAnalyzer cell_mut_probe_1,
							       CellMutAnalyzer cell_mut_probe_2,
							       unsigned int left_border,
							       unsigned int right_border){
	std::vector <unsigned int> cell_mut_fr_1 = cell_mut_probe_1.GetCellMutFr();
	std::vector <unsigned int> cell_mut_fr_2 = cell_mut_probe_2.GetCellMutFr();

	unsigned int max_size = 0;
	//std::cout<< cell_mut_fr_1.size()<<" - " <<cell_mut_fr_2.size()<<std::endl;
		// the length of step
	if (cell_mut_fr_1.size() < cell_mut_fr_2.size()) {
    max_size = cell_mut_fr_2.size();
    unsigned int diff = cell_mut_fr_2.size() - cell_mut_fr_1.size();
    cell_mut_fr_1.insert(cell_mut_fr_1.end(), diff, 0);
  } else {
    max_size = cell_mut_fr_1.size();
    int diff = cell_mut_fr_1.size() - cell_mut_fr_2.size();
    cell_mut_fr_2.insert(cell_mut_fr_2.end(), diff, 0);
  };

  //std::cout<< cell_mut_fr_1.size()<<" - " <<cell_mut_fr_2.size()<<std::endl;


  int sum_diff	=	0;
  if (left_border<=right_border)
		for (unsigned int i = left_border; ((i < max_size) && (i<= right_border)); i++) {
			sum_diff += (cell_mut_fr_2.at(i) - cell_mut_fr_1.at(i)) *
								  (cell_mut_fr_2.at(i) - cell_mut_fr_1.at(i));
//			std::cout<<"i="<<i<<";c_m_1="<<cell_mut_fr_1.at(i)<<
//			";c_m_2="<<cell_mut_fr_2.at(i)<<
//			";sq_diff="<<(cell_mut_fr_2.at(i) - cell_mut_fr_1.at(i)) *
//								  (cell_mut_fr_2.at(i) - cell_mut_fr_1.at(i))<<"\n";
		};
	return sum_diff;
};

std::ostream & operator<<(std::ostream & output,
                          const BalkVAFAnalyzer & balk_vaf_analyzer) {
    output << balk_vaf_analyzer.get_cell_num() << "\t" <<
							balk_vaf_analyzer.get_mut_num() << "\t";
		for (float x : balk_vaf_analyzer.get_num_mutations_vec()) output << x	<< "\t";
    return output;
}

std::istream & operator>>(std::istream & s_in,
													BalkVAFAnalyzer & balk_vaf_analyzer) {
	balk_vaf_analyzer.CleanAll();

	s_in >> balk_vaf_analyzer.num_of_cells;
	s_in >> balk_vaf_analyzer.num_of_muts;

	float current;
	std::vector <unsigned int> current_num_mut_vec;
	do {
		s_in >> current;
    if (!s_in.eof()) current_num_mut_vec.push_back(current);
	} while (!s_in.eof());
	balk_vaf_analyzer.set_num_mutations_vec(current_num_mut_vec);
	balk_vaf_analyzer.InitCumulativeVAFVec();
	return s_in;
};

void BalkVAFAnalyzer::CleanAll(){
		num_of_cells = 0;
    num_of_muts = 0;
    num_mutations_vec.clear();
    cumulative_vaf_vec.clear();
};

BalkVAFAnalyzer::~BalkVAFAnalyzer()
{

};



CellMutAnalyzer::CellMutAnalyzer(){
    genotype_vec.clear();
    cell_vec.clear();
};

CellMutAnalyzer::CellMutAnalyzer(std::vector <Cell> & cell_vector,
																 std::vector<Genotype*> & cell_genotypes) {
  std::set<unsigned int> cell_gene_index_set;
  //inserts cell index of gen to set of gen indexes
  for(Cell cell: cell_vector) cell_gene_index_set.insert(cell.gen);
  //goes through set of gen indexes and gather mutations in set of muts

  std::vector <Genotype> genotype_vector;
	for(unsigned int cell_gene_index: cell_gene_index_set){
	  std::vector <unsigned int> sequence = cell_genotypes.at(cell_gene_index)->sequence;
		Genotype A(0,0);
		A.index = cell_gene_index;
		A.sequence = std::move(sequence);
		genotype_vector.push_back(A);
	};

  genotype_vec = std::move(genotype_vector);
	cell_vec = cell_vector;
};

std::ostream & operator<<(std::ostream & output,
                          const CellMutAnalyzer & cell_mut_analyzer) {
		std::vector <Genotype> genotype_vec = cell_mut_analyzer.GetGenotypeVec();
		std::vector <Cell> cell_vec = cell_mut_analyzer.GetCellVec();

    output << cell_vec.size() << "\t";
		for (Cell x : cell_vec) output << x.gen	<< "\t";
		output << genotype_vec.size() << "\t";
		for (unsigned int i = 0; i < genotype_vec.size(); i++) {
				output << genotype_vec.at(i).index << "\t";
				output << genotype_vec.at(i).sequence.size()  << "\t";
				for (unsigned int x:genotype_vec.at(i).sequence) output << x << "\t";
		};
    return output;
}

std::istream & operator>>(std::istream & s_in,
													CellMutAnalyzer & cell_mut_analyzer) {
	cell_mut_analyzer.CleanAll();
  unsigned int num_of_cells;

	s_in >> num_of_cells;
  std::vector <Cell> cell_vec;

  unsigned int current_cell_gen;
	for (unsigned int i = 0 ; i < num_of_cells; i++) {
			s_in >> current_cell_gen;
			Cell cell;
			cell.gen = current_cell_gen;
			cell.lesion = -1;
			cell.x = -1;
			cell.y = -1;
			cell.z = -1;
			cell_vec.push_back(cell);
	};

	std::vector <Genotype> genotype_vec;

	unsigned int num_of_genotypes;
	s_in >> num_of_genotypes;

	unsigned int current_mut_num;
	for (unsigned int i = 0 ; i < num_of_genotypes; i++) {
			int index;
			s_in >> index;
			s_in >> current_mut_num;
			std::vector <unsigned int> mut_index_genotype_vec;
			unsigned int current_mut;
			for (unsigned int j = 0 ; j < current_mut_num; j++) {
				s_in >> current_mut;
			 	mut_index_genotype_vec.push_back(current_mut);
			};
			Genotype A(0,0);
			A.sequence = mut_index_genotype_vec;
			A.index = index;
			genotype_vec.push_back(A);
	};

	cell_mut_analyzer.SetCellVec(cell_vec);
	cell_mut_analyzer.SetGenotypeVec(genotype_vec);
	return s_in;
};

CellMutAnalyzer::~CellMutAnalyzer()
{

};

