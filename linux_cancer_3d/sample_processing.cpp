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




void Probe_mutations::unique_sorted_cells()
{

  vector<Cell>::iterator it;
  // using predicate comparison:

  it = unique (probe_cells.begin(), probe_cells.end(), Compare_cells_eq());

  probe_cells.resize( std::distance(probe_cells.begin(),it) );


}

void Probe_mutations::unique_sorted_genotypes()
{
  vector<Genotype *>::iterator it;
  // using predicate comparison:

  it = unique (probe_genotypes.begin(), probe_genotypes.end(), Compare_genotypes_eq());

  probe_genotypes.resize( std::distance(probe_genotypes.begin(),it) );


}


void Probe_mutations::write_cells_genotypes()
{
    for (int i=0;i<probe_cells.size(); i++)
    {
            probe_genotypes.push_back(genotypes[probe_cells[i].gen]);
            probe_genotypes[i]->number=1;

    }
}



void Probe_mutations::genotypes_without_inclusions()
{

 //cout << "we are in Probe_mutations::genotypes_without_inclusions()" << endl;
 //cout << "probe_genotypes.size()" << probe_genotypes.size() << endl;


 for(int i=0; i<probe_genotypes.size()-1;i++)
 {

   bool delete_indicator=false;
   for(int j=i+1;
   (
        (j<probe_genotypes.size())   &&
        (delete_indicator == false)

    );j++ )
   if (includes(probe_genotypes[j]->sequence.begin(),probe_genotypes[j]->sequence.end(),
                probe_genotypes[i]->sequence.begin(),probe_genotypes[i]->sequence.end()
               )
      ){
                    excluded_probe_genotypes.push_back(probe_genotypes[i]);

                    probe_genotypes.erase (probe_genotypes.begin()+i);
                    delete_indicator = true;
                    i=i-1;

       };

 };

}


int Probe_mutations::get_number_of_intrsections(Genotype * a, Genotype * b)
{
    int min_ab = (a->sequence.size() > b->sequence.size()) ? b->sequence.size() : a->sequence.size();

    int counter = 0;

    for(int i=0;i<min_ab;i++) {
            if (a->sequence[i] == b->sequence[i]) counter++;
            else break;
    };

    return counter;
}

void Probe_mutations::creat_mut_genotypes_connections()
{
    for(int i = probe_genotypes.size()-2; i>=0;i--)
    {
        int i_max_intersection = 0;// the maximal number of intersection of mutations between genotype with index i and other gens with indexes [i+1,probe_genotypes.size()-1]
        int num_genotype = -1;

        for(int j = probe_genotypes.size()-1; j>=i+1; j--){

            int intersect_num = get_number_of_intrsections(probe_genotypes[i],probe_genotypes[j]);
            if (i_max_intersection < intersect_num)
            {
                    num_genotype = j;
                    i_max_intersection = intersect_num;
            };

        };

        Mut_genotype_connection A;
        A.num_genotype = num_genotype;
        if (i_max_intersection!=-1) A.num_mutation = i_max_intersection;
        mut_genotypes_connections.push_back(A);

    };

    reverse(mut_genotypes_connections.begin(),mut_genotypes_connections.end());

    Mut_genotype_connection last;
    last.num_genotype = -1;
    last.num_mutation = 0;

    mut_genotypes_connections.push_back(last);
}


void Probe_mutations::save_mut_genotypes_connections(char* name)
{

    char name_file[256];

    sprintf(name_file,"%s/Probes_mut_genotypes_connections.dat",name);
    FILE *f=fopen(name_file, "w") ;
    if (f==NULL) err("err");

    int size_mgc = mut_genotypes_connections.size();
    for(int i = 0 ;i < size_mgc; i++)
    {
            fprintf(f,"num_of_attachted_genotype=%d num_genotype=%d num_mutation=%d \n", i, mut_genotypes_connections[i].num_genotype, mut_genotypes_connections[i].num_mutation );
    };
}


void Probe_mutations::save_mutation_vector(char* name)
{

    char name_file[256];

    sprintf(name_file,"%s/mutation_vector.dat",name);
    FILE *f=fopen(name_file, "w") ;
    if (f==NULL) err("err");

    int size_mv = mutation_vector.size();
    for(int i = 0 ;i < size_mv; i++)
    {
        fprintf(f,"%d  abs_mut=%d depth=%d abs_father_mut=%d index_father_mut=%d num_of_child_mut=%d number_cells=%d\n", i,
                mutation_vector[i].abs_mut, mutation_vector[i].depth, mutation_vector[i].abs_father_mut,
                mutation_vector[i].index_father_mut,mutation_vector[i].num_of_child_mut,mutation_vector[i].number_cells);

    };

}

void Probe_mutations::save_graph(char * name)
{
    char name_file[256];

    sprintf(name_file,"%s/mutation_graph.txt",name);
    FILE *f=fopen(name_file, "w") ;
    if (f==NULL) err("err");

    int size_mv = mutation_vector.size();
    fprintf(f,"digraph graphname {\n");
    for(int i = 0 ;i < size_mv; i++)
    {
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

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        fprintf(f,"\n");
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//        fprintf(f,"%d -> %d;\n",
//                mutation_vector[i].abs_father_mut,
//                mutation_vector[i].abs_mut);
        };
    };
    fprintf(f,"}\n");
}


void Probe_mutations::init_zero_abs_mutation(Abs_mutation& zero_mutation)
{
    zero_mutation.abs_mut= 0;
    zero_mutation.abs_father_mut = 0;
    zero_mutation.depth = 0;
    zero_mutation.index_father_mut= 0;

//Calculate number of cells seq without mutations
    zero_mutation.number_cells = 0;
//for(int i=0 ; probe_genotypes[i]->sequence.size()<=0 && (i+1<probe_genotypes.size());i++) zero_mutation.number_cells+=probe_genotypes[i]->number;

    if (probe_genotypes.size()>0)
        if (probe_genotypes[0]-> sequence.size()==0) zero_mutation.number_cells+= probe_genotypes[0]->number;


//    for(int i=0 ; excluded_probe_genotypes[i]->sequence.size()<=0 && (i+1<excluded_probe_genotypes.size());i++) zero_mutation.number_cells+=excluded_probe_genotypes[i]->number;
    if (excluded_probe_genotypes.size()>0)
        if(excluded_probe_genotypes[0]->sequence.size() == 0) zero_mutation.number_cells+=excluded_probe_genotypes[0]->number;

//////////////////////////////////////////////////

//Calculate number of children of zero mutation
    zero_mutation.num_of_child_mut = 0;//?????????????????????
    for(int i=0; i<mut_genotypes_connections.size();i++) if (mut_genotypes_connections[i].num_genotype==-1) zero_mutation.num_of_child_mut+=1;
//////////////////////////////////////////////////

}


void Probe_mutations::creat_mutation_vector()
{
    Abs_mutation zero_mutation;
    init_zero_abs_mutation(zero_mutation);

    mutation_vector.push_back(zero_mutation);


    int size_probe_genotypes=probe_genotypes.size();
    for (int i=size_probe_genotypes-1; i>=0 ;i--){


        int i_genotype_size_sequence = probe_genotypes.at(i)->sequence.size();
        int begin_gen_mutation = mut_genotypes_connections.at(i).num_mutation;
        for(int j = begin_gen_mutation; j < i_genotype_size_sequence; j++)
        {


        // we will adjoin mutation in a sequence to the main mutation
            Abs_mutation current_mutation;
            current_mutation.abs_mut = probe_genotypes.at(i)->sequence.at(j);//
            current_mutation.depth = j + 1;

            if ( mut_genotypes_connections.at(i).num_mutation==0 && j==0){
                current_mutation.abs_father_mut = 0;
                current_mutation.index_father_mut = 0;
            } else {
                current_mutation.abs_father_mut = probe_genotypes.at(i)->sequence.at(j-1);

                if (j>begin_gen_mutation) current_mutation.index_father_mut =  mutation_vector.size()-1;
                else{
                    unsigned int index_father_mut_ = 0;

                    int i_i;
                    for(i_i=size_probe_genotypes-1;i_i > mut_genotypes_connections.at(i).num_genotype; i_i--)
                        index_father_mut_+=(probe_genotypes.at(i_i)->sequence.size() - mut_genotypes_connections.at(i_i).num_mutation);

                    index_father_mut_+= mut_genotypes_connections.at(i).num_mutation - mut_genotypes_connections.at(i_i).num_mutation ;
                    current_mutation.index_father_mut = index_father_mut_;
                };


            };
        // Calculate number of cells with current mutations

            current_mutation.number_cells = 0;
            for(int k=0 ; k < probe_genotypes.size(); k++){
                    if (probe_genotypes.at(k)->sequence.size() > j + 1) break;
                    if (probe_genotypes.at(k)->sequence.size() == j + 1)
                    {
                       if (probe_genotypes.at(k)->sequence.at(j) == current_mutation.abs_mut) current_mutation.number_cells+=probe_genotypes.at(k)->number;
                    }

                    //if (probe_genotypes.at(k)->sequence.at(j) == current_mutation.abs_mut &&
                    //    probe_genotypes.at(k)->sequence.size() == current_mutation.depth)
                    //  cerr<<"i="<<i<<endl;
                    //  cerr<<"k="<<k<<" probe_genotypes.size()="<<probe_genotypes.size() <<endl;
                    //  cerr<<"j="<<j<<" probe_genotypes.at(k)->sequence.size()="<<probe_genotypes.at(k)->sequence.size()<<endl;
                    //  cerr<< "current_mutation.abs_mut" << current_mutation.abs_mut<<endl;

            };
            for(int k=0 ; (k < excluded_probe_genotypes.size()); k++){
                if (excluded_probe_genotypes.at(k)->sequence.size()>j + 1) break;

                    if (excluded_probe_genotypes.at(k)->sequence.size() == j+ 1)
                    {
                        if (excluded_probe_genotypes.at(k)->sequence.at(j) == current_mutation.abs_mut) current_mutation.number_cells+=excluded_probe_genotypes.at(k)->number;
                    }


            };
        ////////////////////////////////////////////////////////////////////////////////////////////////

        //Calculate number of children of zero mutations  probe_genotypes[i]->sequence[0]

            current_mutation.num_of_child_mut = (probe_genotypes[i]->sequence.size() > j + 1) ? 1 : 0;
            for(int p=mut_genotypes_connections.size(); p>=0; p--)
            if (mut_genotypes_connections[p].num_genotype == i &&
                mut_genotypes_connections[p].num_mutation == j + 1) current_mutation.num_of_child_mut+=1;

        mutation_vector.push_back(current_mutation);
        };



    };

    //int sum_cells=0;
    //int pr_sum_cells=0;
    //int ex_sum_cells=0;
    //for(int i=0; i<mutation_vector.size();i++) sum_cells+=mutation_vector[i].number_cells;
    //for(int i=0; i<probe_genotypes.size();i++) pr_sum_cells+=probe_genotypes[i]->number;
    //for(int i=0; i<excluded_probe_genotypes.size();i++) ex_sum_cells+=excluded_probe_genotypes[i]->number;

   // printf("La La La... %d and %d and %d and %d \n", sum_cells,probe_cells.size(), pr_sum_cells, ex_sum_cells) ;
}


Probe_mutations::Probe_mutations(const vector <Cell> &_probe_cells, unsigned int threshold_numb_cells)
 {
    probe_cells=_probe_cells;

    probe_cells.erase(probe_cells.begin() + threshold_numb_cells);

    stable_sort (probe_cells.begin(), probe_cells.end(), Compare_cells());
    //unique_sorted_cells();
    write_cells_genotypes();
    unique_sorted_genotypes();
    genotypes_without_inclusions();
    creat_mut_genotypes_connections();
    creat_mutation_vector();



 }


Probe_mutations::Probe_mutations(const vector <Cell> & _probe_cells)
{
    //cout<< "we are in Probe_mutations::Probe_mutations(const vector <Cell> & _probe_cells) "<< endl;
    probe_cells = _probe_cells;
    stable_sort (probe_cells.begin(), probe_cells.end(), Compare_cells());
    //cout<<"after stable_sort"<<endl;
   //unique_sorted_cells();
    write_cells_genotypes();
   // cout<<"after write_cells_genotypes();"<<endl;

    unique_sorted_genotypes();
   // cout<<"after unique_sorted_genotypes();"<<endl;

    genotypes_without_inclusions();
   // cout<<"after genotypes_without_inclusions();"<<endl;

    creat_mut_genotypes_connections();
   // cout<<"after creat_mut_genotypes_connections();"<<endl;

    creat_mutation_vector();
   // cout<<"aftercreat_mutation_vector();"<<endl;



}

Probe_mutations::Probe_mutations()
{

}


Probe_mutations::~Probe_mutations()
{

}



void Probe_mutations::Save_probe_cells(char *name)
{
    char name_file[256];

    sprintf(name_file,"%s/Probes_sort.dat",name);
    FILE *f=fopen(name_file, "w") ;
    if (f==NULL) err("err");

    for (int i=0;i<probe_cells.size();i++) {

    Genotype *g=genotypes[probe_cells[i].gen] ;

    for (int j=0;j<g->sequence.size();j++) {
//mykola: we cut of additional bits for relevent number
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

void Probe_mutations::save_genotypes(char *name, vector<Genotype *> & probe_genotypes_)
{
    char name_file[256];

    sprintf(name_file,"%s",name);
    FILE *f=fopen(name_file, "w") ;
    if (f==NULL) err("err");

    for (int i=0;i<probe_genotypes_.size();i++) {

    Genotype *g= probe_genotypes_[i];

    fprintf(f,"Genotype_num = %d  ",i);

    for (int j=0;j<g->sequence.size();j++) {
//mykola: we cut of additional bits for relevent number
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

std::ostream & operator<<(std::ostream & output_out, const Probe_mutations & Probe_mutations_)
{

        for(int i=0;i<Probe_mutations_.mutation_vector.size();i++)
        output_out<< Probe_mutations_.mutation_vector[i].abs_mut
        <<"\t"<<Probe_mutations_.mutation_vector[i].depth
        <<"\t"<<Probe_mutations_.mutation_vector[i].abs_father_mut
        <<"\t"<<Probe_mutations_.mutation_vector[i].index_father_mut
        <<"\t"<<Probe_mutations_.mutation_vector[i].num_of_child_mut
        <<"\t"<<Probe_mutations_.mutation_vector[i].number_cells<<"\t";
        return output_out;
}

std::istream& operator>>(std::istream& s_in, Probe_mutations & Probe_mutations_)
{
        Probe_mutations_.mutation_vector.clear();

        do {

            Abs_mutation A;

            s_in >> A.abs_mut>>A.depth >> A.abs_father_mut
            >>A.index_father_mut>>A.num_of_child_mut>>A.number_cells;


            if (!s_in.eof()) Probe_mutations_.mutation_vector.push_back(A);
        } while (!s_in.eof());

        return s_in;
}



float my_chi_square_depth(Probe_mutations & Probe_mut_1, Probe_mutations & Probe_mut_2)
{
    //cout<<"my_chi_square_depth"<<endl;
    vector <int> depth_vector_Pr_1, depth_vector_Pr_2;
/////////////////////////////////////////////////////////////////////////

    int max_depth_1=0;
    int size_mut_vec_1 = Probe_mut_1.mutation_vector.size();
    for (int i=0; i < size_mut_vec_1; i++)
        max_depth_1 = (max_depth_1 < Probe_mut_1.mutation_vector[i].depth) ? Probe_mut_1.mutation_vector[i].depth : max_depth_1;

    int max_depth_2=0;
    int size_mut_vec_2 = Probe_mut_2.mutation_vector.size();
    for (int i=0; i < size_mut_vec_2; i++)
        max_depth_2 = (max_depth_2 < Probe_mut_2.mutation_vector[i].depth) ? Probe_mut_2.mutation_vector[i].depth : max_depth_2;

    int max_depth = (max_depth_1<= max_depth_2) ? max_depth_2 : max_depth_1;

    depth_vector_Pr_1.reserve(max_depth+1);
    depth_vector_Pr_2.reserve(max_depth+1);

    for (int i=0; i <= max_depth; i++)  depth_vector_Pr_1[i]=depth_vector_Pr_2[i]=0;

    for (int i=0; i < size_mut_vec_1; i ++) depth_vector_Pr_1[Probe_mut_1.mutation_vector[i].depth]+=1;
    for (int i=0; i < size_mut_vec_2; i ++) depth_vector_Pr_2[Probe_mut_2.mutation_vector[i].depth]+=1;


    float chi_square=0;
    //cout<< "max_depth_1=" <<max_depth_1<<"\n";
    //cout<< "max_depth=" <<max_depth<<"\n";
    //cout<< "max_depth_2=" <<max_depth_2<<"\n";

    // cout<<"size_mut_vec_1="<<size_mut_vec_1<<
    // "size_mut_vec_2="<<size_mut_vec_2<<endl;

    for (int i=0; i<= max_depth; i++){

     //cout<< "depth_vector_Pr_1[i]" <<depth_vector_Pr_1[i]<<"\n";
        float delta_square = pow((depth_vector_Pr_1[i]*1.0)/(1.0*size_mut_vec_1) - (depth_vector_Pr_2[i]*1.0)/(1.0*size_mut_vec_2),2.0);
    // cout<< "(depth_vector_Pr_1[i]*1.0)/(1.0*size_mut_vec_1)=" <<(depth_vector_Pr_1[i]*1.0)/(1.0*size_mut_vec_1)<<"\n";
        if ((depth_vector_Pr_1[i]!=0) || (depth_vector_Pr_2[i]!=0)) chi_square+= delta_square/((depth_vector_Pr_1[i]*1.0)/(2*1.0*size_mut_vec_1) + (depth_vector_Pr_2[i]*1.0)/(2*1.0*size_mut_vec_2)); //if (depth_vector_Pr_1[i]!=depth_vector_Pr_2[i])
    };
    //cout<<"depth_chi_square="<<chi_square<<"\n";
    return chi_square;


}



float my_chi_square_degree (Probe_mutations & Probe_mut_1, Probe_mutations & Probe_mut_2)
{
    //cout<<"my_chi_square_degree"<<endl;

    vector <int> degree_vector_Pr_1, degree_vector_Pr_2;
/////////////////////////////////////////////////////////////////////////

    int max_degree_1 = 0;
    int size_mut_vec_1 = Probe_mut_1.mutation_vector.size();
    for (int i=0; i < size_mut_vec_1; i++)
        max_degree_1 = (max_degree_1 < Probe_mut_1.mutation_vector[i].num_of_child_mut) ? Probe_mut_1.mutation_vector[i].num_of_child_mut : max_degree_1;

    int max_degree_2 = 0;
    int size_mut_vec_2 = Probe_mut_2.mutation_vector.size();
    for (int i=0; i < size_mut_vec_2; i++)
        max_degree_2 = (max_degree_2 < Probe_mut_2.mutation_vector[i].num_of_child_mut) ? Probe_mut_2.mutation_vector[i].num_of_child_mut : max_degree_2;

    int max_degree = (max_degree_1<= max_degree_2) ? max_degree_2 : max_degree_1;

    degree_vector_Pr_1.reserve(max_degree+1);
    degree_vector_Pr_2.reserve(max_degree+1);

    for (int i=0; i <= max_degree; i++)  degree_vector_Pr_1[i]=degree_vector_Pr_2[i]=0;

    for (int i=0; i < size_mut_vec_1; i ++) degree_vector_Pr_1[Probe_mut_1.mutation_vector[i].num_of_child_mut]+=1;
    for (int i=0; i < size_mut_vec_2; i ++) degree_vector_Pr_2[Probe_mut_2.mutation_vector[i].num_of_child_mut]+=1;


    float chi_square=0;

    for (int i=0; i<= max_degree; i++){


        float delta_square= pow((degree_vector_Pr_1[i]*1.0)/(1.0*size_mut_vec_1) - (degree_vector_Pr_2[i]*1.0)/(1.0*size_mut_vec_2),2.0);
      // cout<< "(degree_vector_Pr_1[i]*1.0)/(1.0*size_mut_vec_1 )=" <<(degree_vector_Pr_1[i]*1.0)/(1.0*size_mut_vec_1)<<"\n";
      // cout<< "(degree_vector_Pr_2[i]*1.0)/(1.0*size_mut_vec_2)=" <<(degree_vector_Pr_2[i]*1.0)/(1.0*size_mut_vec_2)<<"\n";
      //  cout<< "delta_square="<<delta_square<<"\n";

        if ((degree_vector_Pr_1[i]!=0) || (degree_vector_Pr_2[i]!=0))
            chi_square+= delta_square/((degree_vector_Pr_1[i]*1.0)/(2.0*1.0*size_mut_vec_1) + (degree_vector_Pr_2[i]*1.0)/(2.0*1.0*size_mut_vec_2)); // cout<< "chi_square="<<chi_square<<"\n";
    };

    //cout<<"degree_chi_square="<<chi_square<<"\n";
     return chi_square;


}




float tree_dist (Probe_mutations & Probe_mut_1, Probe_mutations & Probe_mut_2)
{
    float chi_square_depth = my_chi_square_depth(Probe_mut_1, Probe_mut_2);
    float chi_square_degree = my_chi_square_degree(Probe_mut_1, Probe_mut_2);

    return (chi_square_depth + chi_square_degree);

}



void init_vector_cells(vector <Cell> & vector_cells, Min_max_diff & info_mmd,Probe_mutations & e_prob,
                       ofstream & file_tree_comparison,
                       float driver_adv,
                       float gama,
                       float driver_prob,
                       float max_dist_for_prob)
{
    int num_etalon_mut = e_prob.mutation_vector.size();

    int num_of_cells_e_prob=0;
    for (int i=0; i<e_prob.mutation_vector.size();i++)
        num_of_cells_e_prob+=e_prob.mutation_vector[i].number_cells;



    int dimension=1;//divide circle in @dimension peaces
    vector <vector <Cell> > v_vector_cells;

    vector<int> indicator;// vector of indicators for break from conditional loop
    vector<short int> min_x;

    for(int i=0;i<dimension;i++) {//init all
        indicator.push_back(0);
        short int m_x=0;
        if (
            (info_mmd.max_x-info_mmd.del_1)>=0
            && ((pow(info_mmd.max_x-info_mmd.del_1,2.0) - pow(i*10,2.0),0.5)>=0)
            )
        {
            m_x = int (pow (pow(info_mmd.max_x-info_mmd.del_1,2.0) - pow(i*10,2.0),0.5));
        } else {
            m_x=0;
        };

        min_x.push_back(m_x);
        vector <Cell> vector_cells_i;

        v_vector_cells.push_back(vector_cells_i);
    };



        //int indicator=0;


    //for(int j=0;j<dimension;j++)
    //    cout <<"min_x["<<j<<"]="<<min_x[j]
    //         <<"max_x["<<j<<"]="<<info_mmd.max_x
    //         <<"min_y["<<j<<"]="<<-int(info_mmd.del_2/2)
    //         <<"max_y["<<j<<"]="<<int(info_mmd.del_2/2)
    //         <<"min_z["<<j<<"]="<<(-int(info_mmd.del_3/2)+int(10*j))
    //         <<"max_z["<<j<<"]="<<(int(info_mmd.del_3/2)+int(10*j))
    //         <<endl;


    for (int i=0; i<cells.size(); i++) //conditional loop
    {

        for(int j=0; j<dimension;j++){

            if (indicator[j] >= num_of_cells_e_prob) break;
                //if (indicator >= num_of_cells_e_prob) break;
            if (
                //(_drand48()<0.5) &&
                (cells[i].x >= min_x.at(j)) && (cells[i].x <= info_mmd.max_x) &&
                (cells[i].y >= -int(info_mmd.del_2/2)) && (cells[i].y <= int(info_mmd.del_2/2)) &&
                (cells[i].z >= (-int(info_mmd.del_3/2)+int(10*j))) && (cells[i].z <= (int(info_mmd.del_3/2)+int(10*j)))
            ){
                    v_vector_cells.at(j).push_back(cells.at(i));
                    indicator[j]++;
                   // indicator++;
            };

        };
    };

    int min_mut_diff = num_etalon_mut;
    int index_min_mut_diff= - 1;

    vector <float> v_dist;
    vector <int> v_delta_mut;
    vector <int> v_n_muts;

    for (int i=0; i<dimension; i++){
        int size_v_vec_cells = v_vector_cells.at(i).size();
        int delta_size_cells = abs( size_v_vec_cells - num_of_cells_e_prob );

        float level_num_cells=0.1;
        float level_num_muts=0.2;
        if (delta_size_cells <= num_of_cells_e_prob*level_num_cells){
        cout<<"successful probe"<<endl;
            Probe_mutations prob_mut(v_vector_cells.at(i));
            int num_prob_mut = prob_mut.mutation_vector.size();
            int delta_prob_etalon_mut = abs(num_etalon_mut - num_prob_mut);
            if (delta_prob_etalon_mut<level_num_muts*num_etalon_mut){

                v_delta_mut.push_back(delta_prob_etalon_mut);
                float dist = tree_dist(prob_mut,e_prob);
                v_dist.push_back(dist);

                v_n_muts.push_back(num_prob_mut);

                if ( delta_size_cells < min_mut_diff) {
                    min_mut_diff = delta_size_cells;
                    index_min_mut_diff = i;
                };
            };

        } else {cout<<"the is only "<<v_vector_cells.at(i).size()<<
                " cells in a probe with number "<<i<<
                " and the demand is " << num_of_cells_e_prob <<" cells"<<endl;
        };

    };

    int size_of_succ_prob = v_delta_mut.size();

    //cout << "size_of_succ_prob="<< size_of_succ_prob<<endl;
    for(int i=0; i<size_of_succ_prob; i++)
    {
        if (v_dist.at(i)<max_dist_for_prob){
            file_tree_comparison.precision(5);
            file_tree_comparison.setf(std::ios::fixed, std:: ios::floatfield);
            file_tree_comparison
            <<driver_adv<<"\t"<<gama<<"\t"<<driver_prob<<"\t"<<v_dist.at(i)
            <<"\t"<<v_delta_mut.at(i)
            <<"\t"<<v_n_muts.at(i)//e_prob.mutation_vector.size()
            <<endl;
           // cout<<"size_of_succ_prob  "<<driver_adv<<" "<<gama<<" "<<driver_prob<<" "<<v_dist.at(i)
           // <<" "<<v_delta_mut.at(i)
           // <<" "<<v_n_muts.at(i)//e_prob.mutation_vector.size()
           // <<endl;
        } else {
            cout<<"v_dist.at("<<i<<")=" <<v_dist.at(i)<<endl;
        };

    };



        //int min_diff_in_mut = num_etalon_mut;
        //int min_diff_in_cells = num_of_cells_e_prob;
        //for (int j=0; j<dimension;j++){
        //    if v_vector_cells.at(j).size()
        //}

}

float take_probe(char *name,Min_max_diff & info_mmd,  ofstream & file_tree_comparison, float driver_adv , float gama, float driver_prob, float  max_dist_for_prob)
{
////Open file with etalon data//////
    Probe_mutations e_prob;
    ifstream file_mut_vec;

    float dist;

    char name_file[256];
    sprintf(name_file, "%s/etalon_mut_vec.dat", name);
    file_mut_vec.open(name_file);
    if(file_mut_vec) {
        file_mut_vec>>e_prob;
        file_mut_vec.close();

        int num_of_cells_e_prob=0;
        for (int i=0; i<e_prob.mutation_vector.size();i++)
        num_of_cells_e_prob+=e_prob.mutation_vector[i].number_cells;
        //cout<<"download num_of_cells_e_prob="<<num_of_cells_e_prob<<endl;

////////////////////////////////////


        vector <Cell> vector_cells;



        init_vector_cells(vector_cells, info_mmd, e_prob, file_tree_comparison, driver_adv,gama,driver_prob, max_dist_for_prob);


        //Graph !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //e_prob.save_graph(name);
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


/*
        if (vector_cells.size()>0){

            Probe_mutations A(vector_cells);
            dist = tree_dist(A,e_prob);

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
    } else {

        cout<<"file "<<name_file<< " doesn't exists"<<endl;
        file_mut_vec.close();
        int i;
        cin >> i;
        //file_tree_comparison<<"file 'name_of_dir/etalon_mut_vec.dat' doesn't exists"<<endl;
        return 0;
    };


  return dist;
}

void init_etalon_vector_of_cells(vector <Cell> &vector_cells, Probe& Pr,unsigned int threshold_num_cells)
{
    //@j counter for fitted cells
    int j=0;
    //for_loop: Find fitted cells (number less then threshold_num_cells)
    for (int i=0; i<cells.size() && j<threshold_num_cells;i++)
    {
        if ((cells[i].x >= Pr.x_min) && (cells[i].x <= Pr.x_max) &&
            (cells[i].y >= Pr.y_min) && (cells[i].y <= Pr.y_max) &&
            (cells[i].z >= Pr.z_min) && (cells[i].z <= Pr.z_max)
            ){
                vector_cells.push_back(cells[i]);
                j++;
            }

    };

}

void small_test(Probe_mutations& A,int size_vector_cells)
{
        int num_of_cells_e_prob=0;
        for (int i=0; i<A.mutation_vector.size();i++)
            num_of_cells_e_prob+=A.mutation_vector[i].number_cells;
        cout<<"num_of_cells_e_prob="<<num_of_cells_e_prob<<endl;
        cout<<"size_vector_cells="<<size_vector_cells<<endl;
        assert(size_vector_cells == num_of_cells_e_prob);
}

void take_etalon_probe(char *name, Probe& Pr, unsigned int threshold_num_cells){

    // take an etalon probe (sample) to build a tree
    vector <Cell> vector_cells;
    vector_cells.clear();

    init_etalon_vector_of_cells(vector_cells, Pr, threshold_num_cells);

    if (vector_cells.size()>0){
        Probe_mutations A(vector_cells);

        //cout<<"size of probe="<<vector_cells.size()<<endl;
        //cout <<"threshold_num_cells="<<threshold_num_cells<<"; vector_cells.size()="<<vector_cells.size()<<endl;
        small_test(A, vector_cells.size());

        //save data in etalon_mut_vec.dat
        char name_file[256]; sprintf(name_file,"%s/etalon_mut_vec.dat",name);
        ofstream etalon_file_mut_vec; etalon_file_mut_vec.open (name_file);
        //cout<<name<<"/etalon_mut_vec.dat is opened"<<endl;
        etalon_file_mut_vec << A; etalon_file_mut_vec.close();

    } else { cout<<"there is no cells in the etalon probe"<< endl;};
}

void Anal_frec::add_one_genotype(Cell one_cell)
{
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
                    //cout<<"g->sequence["<<j<<"] & (~DRIVER_PM) =" <<(g->sequence[j] & (~DRIVER_PM))<<endl;

            };
         } else {
                vec_num_mutations.at((g->sequence[j] & (~RESISTANT_PM)))++;

            };

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 };
}

Anal_frec::Anal_frec(vector <Cell> &A , unsigned int total_number_of_SNV, unsigned int num_of_div)
{
    vec_num_mutations.resize(total_number_of_SNV,0);
    //cout<<"vec_num_mutations size="<<vec_num_mutations.size()<<endl;

    for (int i=0;i<A.size();i++) add_one_genotype(A.at(i));

    frec_in_div.resize(num_of_div,0);

    for (int i=0; i < total_number_of_SNV; i++)
    {
       int detector = 0;

       float proportion = (vec_num_mutations.at(i)*1.0)/(A.size()*1.0)*num_of_div;
       for (int j=0; j < num_of_div; j++ )
       {
            if ((proportion>=j) && (proportion<j+1)) detector=j;
       };

       frec_in_div.at(detector)++;
    };

    //for (int i=0; (i < num_of_div) && (i<200); i++) cout<<"num_of_div["<<i<<"]="<<frec_in_div.at(i)<<endl;


};

Anal_frec::~Anal_frec()
{

};








