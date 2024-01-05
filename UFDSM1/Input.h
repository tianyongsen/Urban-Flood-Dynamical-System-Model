#pragma once
#pragma once
#include<string>
#include<fstream>
#include <memory>
#include<vector>
#include<map>
//just as  class function
class Input
{
private:
	std::string ele_file= "";
	std::string roughness_file = "";
	std::string rain_file = "";
	std::ifstream  infile;

	//×´Ì¬
	bool   file_opened = false;
	//²Ù×÷
	void read_ele_or_roughness_file(std::string& file, std::vector<std::vector<double>>& matrix);
	void read_raindata_file(std::string &rain_file, std::vector<std::vector<std::pair<int, int>>>& aimcell,
		std::vector<std::map<double, double>>& rainmap);
	//
	bool open(std::string &file);
	inline void close();
public:
	Input();
	void init(std::string& ele_file,std::string & roughness_file,std::string& rain_file);
	~Input();
	//
	void get_Celldata( std::vector<std::vector<double>> & ele_matrix, 
		std::vector<std::vector<double>>&roughness_matrix);                //transfer celldata from file to cellfield
	void get_Raindata(std::vector<std::vector<std::pair<int, int>>> &aimcell, std::vector<std::map<double, double>>  &rainmap);
};





