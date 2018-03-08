//i/o includes
#include <iostream>
#include <fstream>
#include <string>
//#include <stdio.h>
//#include <stdlib.h>

//---------------conf file--------------------------------
//begin
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
//end
//---------------conf file--------------------------------
using namespace std;

//#include <boost/filesystem.hpp>

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

void CreateDirectoriesLinux() {

	std::string dist_str("dist_"), fit_adv_str("/fit_adv/"),
		mut_rate_str("/mut_rate/"), dr_mut_rate_str("/dr_mut_rate/"),
		test_str("test"), catalogue_str("statistics/"), space_str(" ");
		system ("mkdir statistics");
		//std::string dir_path(".\\dist_0\\fit_adv\\");

	for(unsigned int i = 0; i <= 12; i++) {
			//std::string dir_name_1("dist_");
			char number_i[256];
			sprintf(number_i,"%d",i);

			std::string s_i(number_i);


			std::string fit_adv_dir_path =catalogue_str + dist_str + s_i ;
			std::string mut_rate_dir_path =catalogue_str + dist_str + s_i ;
			std::string dr_mut_rate_dir_path =catalogue_str + dist_str + s_i ;

			std::string command_line_fit_adv ("mkdir ");
			std::string command_line_mut_rate ("mkdir ");
			std::string command_line_dr_mut_rate ("mkdir ");

			command_line_fit_adv += fit_adv_dir_path;
			command_line_mut_rate += mut_rate_dir_path;
			command_line_dr_mut_rate += dr_mut_rate_dir_path;

			std::cout << command_line_fit_adv << std::endl;
			std::cout << command_line_mut_rate << std::endl;
			std::cout << command_line_dr_mut_rate << std::endl;

			system(command_line_fit_adv.c_str());
			system(command_line_mut_rate.c_str());
			system(command_line_dr_mut_rate.c_str());
	};

	for(unsigned int i = 0; i <= 12;i ++) {

			//std::string dir_name_1("dist_");
			char number_i[256];
			sprintf(number_i,"%d",i);

			std::string s_i(number_i);


			std::string fit_adv_dir_path =catalogue_str + dist_str + s_i + fit_adv_str;
			std::string mut_rate_dir_path =catalogue_str + dist_str + s_i + mut_rate_str ;
			std::string dr_mut_rate_dir_path =catalogue_str + dist_str + s_i + dr_mut_rate_str ;

			std::string command_line_fit_adv ("mkdir ");
			std::string command_line_mut_rate ("mkdir ");
			std::string command_line_dr_mut_rate ("mkdir ");

			command_line_fit_adv += fit_adv_dir_path;
			command_line_mut_rate += mut_rate_dir_path;
			command_line_dr_mut_rate += dr_mut_rate_dir_path;

			std::cout << command_line_fit_adv << std::endl;
			std::cout << command_line_mut_rate << std::endl;
			std::cout << command_line_dr_mut_rate << std::endl;

			system(command_line_fit_adv.c_str());
			system(command_line_mut_rate.c_str());
			system(command_line_dr_mut_rate.c_str());

	};


	for(unsigned int i = 0; i <= 12;i ++) {
		for (unsigned int j = 0; j <10;j ++){
			//std::string dir_name_1("dist_");
			char number_i[256],number_j[256];
			sprintf(number_i,"%d",i);
			sprintf(number_j,"%d",j);
			std::string s_i(number_i);
			std::string s_j(number_j);

			std::string fit_adv_dir_path =catalogue_str + dist_str + s_i + fit_adv_str + test_str + s_j;
			std::string mut_rate_dir_path =catalogue_str + dist_str + s_i + mut_rate_str + test_str + s_j;
			std::string dr_mut_rate_dir_path =catalogue_str + dist_str + s_i + dr_mut_rate_str + test_str + s_j;

			std::string command_line_fit_adv ("mkdir ");
			std::string command_line_mut_rate ("mkdir ");
			std::string command_line_dr_mut_rate ("mkdir ");

			command_line_fit_adv += fit_adv_dir_path;
			command_line_mut_rate += mut_rate_dir_path;
			command_line_dr_mut_rate += dr_mut_rate_dir_path;

			std::cout << command_line_fit_adv << std::endl;
			std::cout << command_line_mut_rate << std::endl;
			std::cout << command_line_dr_mut_rate << std::endl;

			system(command_line_fit_adv.c_str());
			system(command_line_mut_rate.c_str());
			system(command_line_dr_mut_rate.c_str());
		};
	};






};

void CreateParFilesFromExistedOne() {

	std::string dist_str("dist_"), fit_adv_str("/fit_adv/"),
		mut_rate_str("/mut_rate/"), dr_mut_rate_str("/dr_mut_rate/"),
		test_str("test"),slash_str("/");
	std::string par_file_name("model_inference_info"),underscore("_"), filename_extension(".info");
	std::string old_par_file_name = par_file_name + filename_extension;

	//	std::string dir_path(".\\dist_0\\fit_adv\\");

	for(unsigned int i = 0; i <= 12;i ++) {
		char number_i[256];
		sprintf(number_i,"%d",i);
		std::string s_i(number_i);

		std::string fit_adv_dir_path = dist_str + s_i + fit_adv_str  + slash_str;
		std::string mut_rate_dir_path = dist_str + s_i + mut_rate_str  + slash_str;
		std::string dr_mut_rate_dir_path = dist_str + s_i + dr_mut_rate_str  + slash_str;

		std::string full_old_par_file_name_in_fit_adv = fit_adv_dir_path + old_par_file_name;
		std::string full_old_par_file_name_in_mut_rate = mut_rate_dir_path  + old_par_file_name;
		std::string full_old_par_file_name_in_dr_mut_rate = dr_mut_rate_dir_path + old_par_file_name;

		boost::property_tree::ptree pt_fit_adv = GetPars(full_old_par_file_name_in_fit_adv);
		boost::property_tree::ptree pt_mut_rate = GetPars(full_old_par_file_name_in_mut_rate);
		boost::property_tree::ptree pt_dr_mut_rate = GetPars(full_old_par_file_name_in_dr_mut_rate);

		for (unsigned int j = 1; j <= 10; j++) {
			char number_j[256];
			sprintf(number_j, "%d", j);
      std::string s_j(number_j);

      std::string full_new_par_file_name_in_fit_adv= fit_adv_dir_path +
        par_file_name + underscore + s_j + filename_extension;
			std::string full_new_par_file_name_in_mut_rate = mut_rate_dir_path +
			  par_file_name + underscore + s_j + filename_extension;
			std::string full_new_par_file_name_in_dr_mut_rate = dr_mut_rate_dir_path +
			  par_file_name + underscore + s_j + filename_extension;

			std::string new_directory = test_str + s_j;
			pt_fit_adv.put ("basic_settings.file_name", new_directory.c_str());
			pt_mut_rate.put ("basic_settings.file_name", new_directory.c_str());
			pt_dr_mut_rate.put ("basic_settings.file_name", new_directory.c_str());
			pt_fit_adv.put ("Parameters.ProbeComparisonParameters.cataloge_with_e_mut_file", new_directory.c_str());
			pt_mut_rate.put ("Parameters.ProbeComparisonParameters.cataloge_with_e_mut_file", new_directory.c_str());
			pt_dr_mut_rate.put ("Parameters.ProbeComparisonParameters.cataloge_with_e_mut_file", new_directory.c_str());

			CreateInferenceConfFile(pt_fit_adv, full_new_par_file_name_in_fit_adv);
			CreateInferenceConfFile(pt_mut_rate, full_new_par_file_name_in_mut_rate);
			CreateInferenceConfFile(pt_dr_mut_rate, full_new_par_file_name_in_dr_mut_rate);

			std::cout << full_new_par_file_name_in_fit_adv << std::endl;
			std::cout << full_new_par_file_name_in_mut_rate<< std::endl;
			std::cout << full_new_par_file_name_in_dr_mut_rate << std::endl;

		};
	};

};

void DuplicateEtalonMutVectorFileToAllTests() {

   //cp -r dir1/* dir2


	std::string dist_str("dist_"), fit_adv_str("/fit_adv/"),
		mut_rate_str("/mut_rate/"), dr_mut_rate_str("/dr_mut_rate/"),slash_str("/"),
		test_str("test"), etalon_mut_vec_str("etalon_mut_vec.dat ");
	//	std::string dir_path(".\\dist_0\\fit_adv\\");
	for(unsigned int i = 0; i <= 12;i ++) {
		for (unsigned int j = 1; j <= 10;j ++){
			//std::string dir_name_1("dist_");
			char number_i[256],number_j[256];
			sprintf(number_i,"%d",i);
			sprintf(number_j,"%d",j);
			std::string s_i(number_i);
			std::string s_j(number_j);

			std::string fit_adv_dir_path = dist_str + s_i + fit_adv_str + test_str + s_j;
			std::string mut_rate_dir_path = dist_str + s_i + mut_rate_str + test_str + s_j;
			std::string dr_mut_rate_dir_path = dist_str + s_i + dr_mut_rate_str + test_str + s_j;

			std::string command_line_fit_adv ("cp ");
			std::string command_line_mut_rate ("cp ");
			std::string command_line_dr_mut_rate ("cp ");

			command_line_fit_adv += (etalon_mut_vec_str + fit_adv_dir_path);
			command_line_mut_rate += (etalon_mut_vec_str + mut_rate_dir_path);
			command_line_dr_mut_rate += (etalon_mut_vec_str + dr_mut_rate_dir_path);

			std::cout << command_line_fit_adv << std::endl;
			std::cout << command_line_mut_rate << std::endl;
			std::cout << command_line_dr_mut_rate << std::endl;

			system(command_line_fit_adv.c_str());
			system(command_line_mut_rate.c_str());
			system(command_line_dr_mut_rate.c_str());
		};
	};


}

void DuplicateLogsAndSomeDataFiles() {

	//cp -r dir1/* dir2
 	std::string dist_str("dist_"), fit_adv_str("/fit_adv/"),
		mut_rate_str("/mut_rate/"), dr_mut_rate_str("/dr_mut_rate/"),
		test_str("test"), catalogue_str("statistics/"), space_str(" "),slash_str("/"),
		file_log_name_str("_file_tree_comparison.dat "), file_logs("logs.dat ");


  //etalon_mut_vec_str("etalon_mut_vec.dat ");
	//	std::string dir_path(".\\dist_0\\fit_adv\\");
	for(unsigned int i = 0; i <= 12; i++) {
	 if (i!=4 && i!=5){
		for (unsigned int j = 0; j <= 10; j++){
			//std::string dir_name_1("dist_");
			char number_i[256],number_j[256];
			sprintf(number_i,"%d",i);
			sprintf(number_j,"%d",j);
			std::string s_i(number_i);
			std::string s_j(number_j);

			std::string fit_adv_dir_path = dist_str + s_i + fit_adv_str + test_str + s_j;
			std::string mut_rate_dir_path = dist_str + s_i + mut_rate_str + test_str + s_j;
			std::string dr_mut_rate_dir_path = dist_str + s_i + dr_mut_rate_str + test_str + s_j;



		  //copy 9_file_tree_comparison.dat, 19_file_tree_comparison.dat, ..., 99_file_tree_comparison.dat
      for(unsigned int k = 1; k <= 10; k++) {
					char number_k[256];

					sprintf(number_k,"%d", (k * 10 - 1));
					std::string s_k(number_k);
          std::string command_line_fit_adv ("cp ");
			    std::string command_line_mut_rate ("cp ");
			    std::string command_line_dr_mut_rate ("cp ");
					command_line_fit_adv += (fit_adv_dir_path + slash_str + s_k + file_log_name_str + catalogue_str + fit_adv_dir_path);
					command_line_mut_rate += (mut_rate_dir_path + slash_str + s_k + file_log_name_str + catalogue_str + mut_rate_dir_path);
					command_line_dr_mut_rate += (dr_mut_rate_dir_path + slash_str + s_k + file_log_name_str + catalogue_str + dr_mut_rate_dir_path);

					system(command_line_fit_adv.c_str());
					system(command_line_mut_rate.c_str());
					system(command_line_dr_mut_rate.c_str());
      };
        //copy logs.dat
			std::string command_line_fit_adv ("cp ");
			std::string command_line_mut_rate ("cp ");
			std::string command_line_dr_mut_rate ("cp ");
			command_line_fit_adv += (fit_adv_dir_path + slash_str + file_logs + catalogue_str + fit_adv_dir_path);
			command_line_mut_rate += (mut_rate_dir_path + slash_str + file_logs + catalogue_str + mut_rate_dir_path);
			command_line_dr_mut_rate += (dr_mut_rate_dir_path + slash_str + file_logs + catalogue_str + dr_mut_rate_dir_path);
			system(command_line_fit_adv.c_str());
			system(command_line_mut_rate.c_str());
			system(command_line_dr_mut_rate.c_str());

			std::cout << command_line_fit_adv << std::endl;
			std::cout << command_line_mut_rate << std::endl;
			std::cout << command_line_dr_mut_rate << std::endl;
		};
	 };
	};


}


void COne() {

	std::string dist_str("dist_"), fit_adv_str("/fit_adv/"),
		mut_rate_str("/mut_rate/"), dr_mut_rate_str("/dr_mut_rate/"),
		test_str("test"),slash_str("/");
	std::string par_file_name("model_inference_info"),underscore("_"), filename_extension(".info");
	std::string old_par_file_name = par_file_name + filename_extension;

	//	std::string dir_path(".\\dist_0\\fit_adv\\");

	for(unsigned int i = 0; i <= 12;i ++) {
		char number_i[256];
		sprintf(number_i,"%d",i);
		std::string s_i(number_i);

		std::string fit_adv_dir_path = dist_str + s_i + fit_adv_str  + slash_str;
		std::string mut_rate_dir_path = dist_str + s_i + mut_rate_str  + slash_str;
		std::string dr_mut_rate_dir_path = dist_str + s_i + dr_mut_rate_str  + slash_str;

		std::string full_old_par_file_name_in_fit_adv = fit_adv_dir_path + old_par_file_name;
		std::string full_old_par_file_name_in_mut_rate = mut_rate_dir_path  + old_par_file_name;
		std::string full_old_par_file_name_in_dr_mut_rate = dr_mut_rate_dir_path + old_par_file_name;

		boost::property_tree::ptree pt_fit_adv = GetPars(full_old_par_file_name_in_fit_adv);
		boost::property_tree::ptree pt_mut_rate = GetPars(full_old_par_file_name_in_mut_rate);
		boost::property_tree::ptree pt_dr_mut_rate = GetPars(full_old_par_file_name_in_dr_mut_rate);


			char number_j[256];
			sprintf(number_j, "%d", 0);
      std::string s_j(number_j);

      std::string full_new_par_file_name_in_fit_adv= fit_adv_dir_path +
        par_file_name + underscore + s_j + filename_extension;
			std::string full_new_par_file_name_in_mut_rate = mut_rate_dir_path +
			  par_file_name + underscore + s_j + filename_extension;
			std::string full_new_par_file_name_in_dr_mut_rate = dr_mut_rate_dir_path +
			  par_file_name + underscore + s_j + filename_extension;

			std::string new_directory = test_str + s_j;
			pt_fit_adv.put ("basic_settings.file_name", new_directory.c_str());
			pt_mut_rate.put ("basic_settings.file_name", new_directory.c_str());
			pt_dr_mut_rate.put ("basic_settings.file_name", new_directory.c_str());
			pt_fit_adv.put ("Parameters.ProbeComparisonParameters.cataloge_with_e_mut_file", new_directory.c_str());
			pt_mut_rate.put ("Parameters.ProbeComparisonParameters.cataloge_with_e_mut_file", new_directory.c_str());
			pt_dr_mut_rate.put ("Parameters.ProbeComparisonParameters.cataloge_with_e_mut_file", new_directory.c_str());

			CreateInferenceConfFile(pt_fit_adv, full_new_par_file_name_in_fit_adv);
			CreateInferenceConfFile(pt_mut_rate, full_new_par_file_name_in_mut_rate);
			CreateInferenceConfFile(pt_dr_mut_rate, full_new_par_file_name_in_dr_mut_rate);

			std::cout << full_new_par_file_name_in_fit_adv << std::endl;
			std::cout << full_new_par_file_name_in_mut_rate<< std::endl;
			std::cout << full_new_par_file_name_in_dr_mut_rate << std::endl;


	};

};

struct Data{
	int iteration_num;
	int trials;
	float	driver_adv_true;
  float r_err_driver_adv;
  float samp_var_driver_adv;
  float mut_rate_true;
	float r_err_mut_rate;
	float samp_var_mut_rate;
	float driver_mut_rate_true;
	float r_err_driver_mut_rate;
	float samp_var_driver_mut_rate;
};

struct Results{
	float r_err_driver_adv;
	float r_err_mut_rate;
	float r_err_driver_mut_rate;
	float abs_min;
	unsigned int time_mean_hit;
	unsigned int time_min_hit;
};


Results ReadFromFile(std::string cataloge_path) {// 0 - fit_ad, 1 - mut_rate, 2 - dr_mut_rate
	std::ifstream log_tree_dist_comparison;
	char name_file[256];
	sprintf(name_file,"%s/logs.dat",cataloge_path.c_str());
	log_tree_dist_comparison.open(name_file);
	log_tree_dist_comparison.precision(5);
	log_tree_dist_comparison.setf( std::ios::fixed, std:: ios::floatfield);
	std::vector <Data> vec_Data;

	Data A;
  while(!log_tree_dist_comparison.eof()) {
		log_tree_dist_comparison >> A.iteration_num>>
		A.trials>>
		A.driver_adv_true>>
		A.r_err_driver_adv>>
		A.samp_var_driver_adv>>
		A.mut_rate_true>>
		A.r_err_mut_rate>>
		A.samp_var_mut_rate>>
		A.driver_mut_rate_true>>
		A.r_err_driver_mut_rate>>
		A.samp_var_driver_mut_rate;
		vec_Data.push_back(A);
  };

  Results results={0,0,0,0,0,0};
  std::cout<<"vec_size"<< vec_Data.size()<<std::endl;
	for (unsigned int i = 0; i < vec_Data.size(); i++ ) {
		results.r_err_driver_adv 			+= 	std::abs(vec_Data[i].r_err_driver_adv);
		results.r_err_mut_rate   			+=	std::abs(vec_Data[i].r_err_mut_rate);
		results.r_err_driver_mut_rate += 	std::abs(vec_Data[i].r_err_driver_mut_rate);

	};
  unsigned int v_size = vec_Data.size();
	results.r_err_driver_adv 			= results.r_err_driver_adv/(1.0* v_size);
	results.r_err_mut_rate   			= results.r_err_mut_rate/(1.0* v_size);
	results.r_err_driver_mut_rate = results.r_err_driver_mut_rate/(1.0 * v_size);

	std::cout<<"mean_r_err_driver_adv = "<< results.r_err_driver_adv<<std::endl;
  std::cout<<"mean_r_err_mut_rate = "<< results.r_err_mut_rate<<std::endl;
  std::cout<<"mean_err_driver_mut_rate = "<< results.r_err_driver_mut_rate<<std::endl;

//	float abs_min = 10;
//
//	for (unsigned int i = 0; i < vec_Data.size(); i++ ) {
//		if ((results.r_err_driver_adv ! = 0) && (abs_min> std::abs(vec_Data[i].r_err_driver_adv))
//						abs_min = std::abs(vec_Data[i].r_err_driver_adv);
//	  if  (results.r_err_mut_rate != 0)
//		if  (results.r_err_mut_rate != 0)
//	};
//


	log_tree_dist_comparison.close();
  return results;
};

void CreateSummary() {

	//cp -r dir1/* dir2
 	std::string s_fit_av_str("fit_ad_summary.dat"), s_mut_rate_str("mut_rat_summary.dat"),
							s_dr_mut_rate_str("dr_mut_rat_summary.dat");

	char name_file[256];
	/*file_tree_comparison.dat*/
	sprintf(name_file,"%s", s_fit_av_str.c_str());
	remove(name_file);
	sprintf(name_file,"%s", s_mut_rate_str.c_str());
	remove(name_file);
	sprintf(name_file,"%s", s_dr_mut_rate_str.c_str());
	remove(name_file);

 	std::string dist_str("dist_"), fit_adv_str("/fit_adv/"),
		mut_rate_str("/mut_rate/"), dr_mut_rate_str("/dr_mut_rate/"),
		test_str("test"), catalogue_str("statistics/"), space_str(" "),slash_str("/"),
		file_log_name_str("_file_tree_comparison.dat "), file_logs("logs.dat ");


  //etalon_mut_vec_str("etalon_mut_vec.dat ");
	//	std::string dir_path(".\\dist_0\\fit_adv\\");
	for(unsigned int i = 0; i <= 12; i++) {
	 if (i!=4 && i!=5){
		for (unsigned int j = 0; j <= 10; j++){
			//std::string dir_name_1("dist_");
			char number_i[256],number_j[256];
			sprintf(number_i,"%d",i);
			sprintf(number_j,"%d",j);
			std::string s_i(number_i);
			std::string s_j(number_j);

			std::string fit_adv_dir_path = dist_str + s_i + fit_adv_str + test_str + s_j;
			std::string mut_rate_dir_path = dist_str + s_i + mut_rate_str + test_str + s_j;
			std::string dr_mut_rate_dir_path = dist_str + s_i + dr_mut_rate_str + test_str + s_j;

			Results A_fit_mut_rate = ReadFromFile(fit_adv_dir_path);
			Results	A_mut_rate_dir_path = ReadFromFile(mut_rate_dir_path);
			Results A_dr_mut_rate_dir_path = ReadFromFile (dr_mut_rate_dir_path);

			std::ofstream s_fit_av_file;
			s_fit_av_file.open(s_fit_av_str.c_str(), std::ios_base::app);
			s_fit_av_file.precision(5);
			s_fit_av_file.setf( std::ios::fixed, std:: ios::floatfield);
			s_fit_av_file << i <<"\t"<< A_fit_mut_rate.r_err_driver_adv << std::endl;
			s_fit_av_file.close();

			std::ofstream mut_rate_file;
			mut_rate_file.open(s_mut_rate_str.c_str(),std::ios_base::app);
			mut_rate_file.precision(5);
			mut_rate_file.setf( std::ios::fixed, std:: ios::floatfield);
			mut_rate_file << i <<"\t"<< A_mut_rate_dir_path.r_err_mut_rate << std::endl;
			mut_rate_file.close();

			std::ofstream dr_mut_rate_file;
			dr_mut_rate_file.open(s_dr_mut_rate_str.c_str(),std::ios_base::app);
			dr_mut_rate_file.precision(5);
			dr_mut_rate_file.setf( std::ios::fixed, std:: ios::floatfield);
			dr_mut_rate_file << i <<"\t"<< A_dr_mut_rate_dir_path.r_err_driver_mut_rate << std::endl;
			dr_mut_rate_file.close();
		};
	 };
	};



};

int main()
{
	//std::string par_file_name("model_inference_info.info");
  //CreateDirectoriesLinux();
	//CreateParFilesFromExistedOne();
	//DuplicateLogsAndSomeDataFiles();

	//std::string A("001");
	//ReadFromFile(A);
	CreateSummary();

	//DuplicateEtalonMutVectorFileToAllTests();
	//CreateInferenceConfFile(pt, full_way_to_config_file);

//  char command_line_0[256], command_line_1[256], command_line_2[256],
//       command_line_3[256], command_line_4[256];
//
//	sprintf(command_line_2,"cd fit_adv &");
//	sprintf(command_line_3,"ls");



//  for(unsigned int i=0; i <= 0;i ++) {
//    sprintf(command_line_0,"dist_%d\\fit_ad\\", i);
//		system ("dir");
//    std::string catalog_path(command_line_0);
//    std::string full_way_to_par_file_name= catalog_path + par_file_name;
//    std::cout<< full_way_to_par_file_name<<std::endl;
//		boost::property_tree::ptree pt = GetPars(par_file_name);
//		std::string full_way_to_config_file("go.info");
//C:\\Users\\lebidm\\Desktop\\Work\\My_Articles\\article_spatial_inference_Martin_Novak\\3DEdenInf\\FSInference
//  };




    //std::cout << "Hello world!" << std::endl;
    return 0;
}



