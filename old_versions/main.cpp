//i/o includes
#include <iostream>
//#include <string>
#include <stdio.h>
#include <stdlib.h>

//---------------conf file--------------------------------
//begin
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
//end
//---------------conf file--------------------------------


//initializes parameters
boost::property_tree::ptree GetPars(std::string par_file_name) {

  boost::property_tree::ptree pt;
  boost::property_tree::info_parser::read_info(par_file_name, pt);
  // begin* Name a directory with etalon_mut_vec.dat (etalon mutation vector)
	return pt;
};

void CreateInferenceConfFile(boost::property_tree::ptree & pt, std::string full_way_to_config_file) {
  boost::property_tree::info_parser::write_info(full_way_to_config_file, pt);
};



int main()
{
    std::string par_file_name("model_inference_info.info");
    std::string full_way_to_config_file("go.dat");

    //boost::property_tree::ptree pt = GetPars(par_file_name);

    //CreateInferenceConfFile(pt, full_way_to_config_file);

  char command_line[256];
  for(unsigned int i=0; i <= 1;i ++) {
    sprintf(command_line,"cd dist_%d",i);
    system (command_line);
    sprintf(command_line,"cd fit_adv");
    system (command_line);
    sprintf(command_line,"ls");
    system (command_line);
    boost::property_tree::ptree pt = GetPars(par_file_name);
    std::string full_way_to_config_file("go.dat");
    CreateInferenceConfFile(pt, full_way_to_config_file);
  };

//
//
//    sprintf(command_line,"cd sex");
//#endif
//    //cout <<"command_line[256]="<<command_line<<endl;
//    system (command_line);
//#if defined __linux
//    sprintf(command_line,"mkdir sex_1");
//#else
//    sprintf(command_line,"mkdir sex/sex_1 sex/sex_2 sex/sex_3");
//#endif
//		system (command_line);
//    //CreateInferenceConfFile(pars);
    std::cout << "Hello world!" << std::endl;
    return 0;
}



