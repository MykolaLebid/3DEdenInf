#ifndef SAMPLE_PROCESSING_H_INCLUDED
#define SAMPLE_PROCESSING_H_INCLUDED

#include <stdio.h>
#include <math.h>
#include <vector>
#include "classes.h"
#include <iostream>
#include <fstream>


//#include "params.h"

#include <stdexcept>// out_of_range


extern vector<Genotype*> genotypes;


// structure for taking all sells in Rectangular cuboid
struct Probe{
    int x_min; int x_max;
    int y_min; int y_max;
    int z_min; int z_max;
};


struct Min_max_diff{
int min_x;int max_x;
int min_y;int max_y;
int min_z;int max_z;
int del_1; int del_2;
int del_3;
};
// function: 1) takes all cells in Rectangular cuboid @param Probe Pr and
//              print them in a file with a name @param char *name
//           2) create prob
float take_probe(char *name, Min_max_diff & info_mmd,  ofstream & file_tree_comparison, float driver_adv, float gama, float driver_prob, float max_dist_for_prob);
void take_etalon_probe(char *name, Probe & Pr, unsigned int number_of_cells);



// structure Compare_cells
// support sort algorithm for probe_cells by length of genotype sequence (of mutations)
struct Compare_cells{

    inline bool operator() (const Cell & a,const Cell & b)
    {

        unsigned int size_a = genotypes.at(a.gen)->sequence.size();
        unsigned int size_b = genotypes.at(b.gen)->sequence.size();

        if (size_a != size_b || size_a == 0){
            return (size_a < size_b);

        } else return (genotypes[a.gen]->sequence.at(size_a-1) <
                     genotypes[b.gen]->sequence.at(size_b-1));
    };



};
///////////////////////////////////////////////////////////////////////////////////////

// structure Compare_cells
// support unique algorithm for sorted probe_cells by deletion the same cells
struct Compare_cells_eq{

    inline bool operator() (const Cell & a, const Cell & b)
    {
        unsigned int size_a= genotypes.at(a.gen)->sequence.size();
        unsigned int size_b= genotypes.at(b.gen)->sequence.size();

        unsigned int last_a = -1;
        unsigned int last_b = -1;
        if (size_a !=0 ) last_a = genotypes.at(a.gen)->sequence.at(size_a-1);
        if (size_b !=0 ) last_b = genotypes.at(b.gen)->sequence.at(size_b-1);

        if ((size_a == size_b) && (last_a == last_b)) {
            return true;
        } else return false;
    };
};
////////////////////////////////////////////////////////////////////////////////////////


// structure Compare_cells
// support exclusion algorithm
struct Compare_genotypes_eq{

    inline bool operator() ( Genotype * a, Genotype * b)
    {

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
};
////////////////////////////////////////////////////////////////////////////////////////

// structure included_genotype
// support unique algorithm for exclusion probe_cells by deletion the same cells///////////
struct included_genotype{

    inline bool operator() ( Genotype * a, Genotype * b)
    {

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
};
////////////////////////////////////////////////////////////////////////////////////////






//Structure Prob_branch_connaction
// we use this  node structure after:
//1) sorting genotype by length @function unique_sorted_genotypes();
//2) after exclusions of inclusions @function genotypes_without_inclusions()///////////
struct Mut_genotype_connection{
        unsigned int num_genotype;// genotype A to which we attach correspondent (smaller) genotype B
        unsigned int num_mutation;// the place in line genotype A  where we attach the a genotype B
} ;
//example:
//Genotype B= probe_genotypes[6].sequence="  p 51 p 789 p 3465 p 12456"
//Genotype A= probe_genotypes[7].sequence=" p 51 p 789 p 3465 d 4570 p 10020"
// -> mut_genotypes_connections[6].num_genotype=7
// -> mut_genotypes_connections[6].num_mutation=2
///////////////////////////////////////////////////////////////////////////////////////


//Structure Abs_mutation
// we use this  node structure for:
//1) comparison of trees;
//2) all mutations is located in one line and they are arranged by the depth of mute///////////
struct Abs_mutation{
        unsigned int abs_mut;// absolute number of this mutation
        unsigned int depth;// distance to the main mutation
        unsigned int abs_father_mut;
        unsigned int index_father_mut;//index in a mutation_vector
        unsigned int num_of_child_mut;
        unsigned int number_cells; //number of cells with this mutation
} ;

///////////////////////////////////////////////////////////////////////////////////////




// Probe_mutations
// Class take vector of sells and work with mutation
///////////////////////////////////////////////////////////////////////////////////////
class Probe_mutations{

public:

    //Constructor


    Probe_mutations(const vector <Cell> &, unsigned int );

    Probe_mutations(const vector <Cell> &);
    Probe_mutations();
    //@param const vector <Cell> & - Cells in Probe

    //Destructor
    ~Probe_mutations();

    // Accessor functions

    //Should be deleted later++++++++++++++++
    void Save_probe_cells(char *);
    void save_genotypes(char *, vector<Genotype*> &);
    void save_mut_genotypes_connections(char*);
    void save_mutation_vector(char*);
    void save_graph(char *);//create a tex file with simple dot graph construction

    //+++++++++++++++++++++++++++++++++++++++

    // Member variables

    vector <Cell> probe_cells;
    vector <Genotype *> probe_genotypes;
    vector <Mut_genotype_connection> mut_genotypes_connections;
    vector <Genotype *> excluded_probe_genotypes;
    vector <Abs_mutation> mutation_vector;

private:


 // Support functions



    //Mutator functions
    void unique_sorted_cells();
    void unique_sorted_genotypes();
    void genotypes_without_inclusions();//after 1) void unique_sorted_genotypes();
    void creat_mut_genotypes_connections();// after 1) void unique_sorted_genotypes() and 2) void genotypes_without_inclusions();

    void init_zero_abs_mutation(Abs_mutation& );// support function for void creat_mutation_vector();
    void creat_mutation_vector();

    //Add functions
    void write_cells_genotypes();
    string get_string_format_mut(unsigned int number);
    //Accessor functions
    int get_number_of_intrsections(Genotype * ,Genotype *);

    //File functions

    //
    friend std::ostream & operator<<(std::ostream & output_out, const Probe_mutations & Probe_mutations_);


    friend std::istream& operator>>(std::istream& s_in, Probe_mutations & Probe_mutations_);


    friend float tree_dist (Probe_mutations & Probe_mut_1, Probe_mutations & Probe_mut_2);

// dist tree  function

    friend float my_chi_square_depth (Probe_mutations & Probe_mut_1, Probe_mutations & Probe_mut_2);
    friend float my_chi_square_degree (Probe_mutations & Probe_mut_1, Probe_mutations & Probe_mut_2);

};



class Anal_frec{

public:

    //Constructor

    Anal_frec(vector <Cell> & A,unsigned int total_number_of_SNV, unsigned int num_of_div);


    //Destructor
    ~Anal_frec();

    // Member variables
    unsigned int num_of_div;
    unsigned int total_number_of_SNV;

    vector <int> vec_num_mutations;
    vector <int> frec_in_div;

    //Support functions
    void add_one_genotype(Cell one_cell);
};


#endif // SAMPLE_PROCESSING_H_INCLUDED
