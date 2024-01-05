#include<sstream>
#include<string>
#include<utility>
#include"Input.h"

Input::Input()
{	}
Input::~Input()
{
	close();
}
void  Input::init(std::string& ele_file_, std::string& roughness_file_, std::string& rain_file_)
{
	ele_file = ele_file_;
	roughness_file = roughness_file_;
	rain_file = rain_file_;
}

void Input::get_Celldata(std::vector<std::vector<double>>& ele_matrix,
	std::vector<std::vector<double>>&roughness_matrix)
	//get celldata from file to deliver to  cellfield  
{	
	read_ele_or_roughness_file(ele_file, ele_matrix);                           
	read_ele_or_roughness_file(roughness_file, roughness_matrix);
}
void Input::read_ele_or_roughness_file(std::string& file, std::vector<std::vector<double>>& matrix)
// function: read the file storing the elevation information 
// input:  elevation file name ;   //output:
{
	close();
	open(file);
	infile.seekg(std::ios_base::beg);           //定位到文件开头
	std::string line;
	double rows = 0.;
	double cols = 0.;
	while (!infile.eof())
	{
		std::getline(infile, line);     //默认终止字符为换行字符
		if (line[0] == 'n' && line[1] == 'c' && line[2] == 'o')		//找到标头  "nclos"
		{
			std::istringstream iss_0(line);
			std::string temp;
			iss_0 >> temp;          //去掉标头
			iss_0 >> cols;             //读取列数			

		}
		if (line[0] == 'n' && line[1] == 'r' && line[2] == 'o')		//找到标头  "nrows"
		{
			std::istringstream iss_1(line);
			std::string temp;
			iss_1 >> temp;          //去掉标头
			iss_1 >> rows;             //读取行数			
		}
		if (line[0] == 'N')
		{
			break;                      //定位到数据上一行
		}
	}
	int k = 0;
	
	while (!infile.eof())
	{
		std::getline(infile, line);     //读取数据
		std::istringstream iss(line);
		//可添加数据检验功能 ,之后再说吧*** 假定输入是规范的。产生的celldata也是规范的
		double data = 0.;
		matrix.push_back({});    //先添加一行
		while (iss >> data)
		{
			matrix[k].push_back(data);
		}
		k++;       //下一行
	}
	//check the data numbers
	if (rows != matrix.size() or cols != matrix[0].size())
		throw "the  rows or cols of cellddata is not corectly corresping to the celldata matrix.";
	close();
}
void Input::get_Raindata(std::vector<std::vector<std::pair<int, int>>>& aimcell, std::vector<std::map<double, double>>& rainmap)
{	
	read_raindata_file(rain_file, aimcell, rainmap);
}
void Input::read_raindata_file(std::string &rain_file, std::vector<std::vector<std::pair<int, int>>>&  aimcell, 
	std::vector<std::map<double, double>>& rainmap)
{
	aimcell.clear();   //清空
	rainmap.clear();  
	close(); 
	open(rain_file);
	infile.seekg(std::ios_base::beg);           //定位到文件开头
	std::string line;
	int index = 0;
	while (!infile.eof())
	{
		std::getline(infile, line);     //默认终止字符为换行字符
		if (line[0] == 'a' && line[1] == 'i' && line[2] == 'm')		//找到标头  [rain] or [inflow]
		{
			//读取到标头之后，下面的格式比如按照规定，要不然报错。以下代码便是依据这个规则写的。
			//规则：
			//aim_cell_index - 1 - 1
			//	time
			//	0 1 2 3 4 5 6 7 8 300
			//	intensity
			//	0 0 6.66667 6.66667 6.66667 0 0 0 0 0
			std::istringstream iss_0(line);
			std::string temp;
			iss_0 >> temp;    //删去无用的标头
			int row = -2, col = -2;
			aimcell.push_back({});
			while (iss_0 >> row && iss_0 >> col)    //一读读两个 ，要不然失败
			{
				aimcell[index].push_back(std::make_pair(row, col));
			}		
			std::getline(infile, line);
			if (line[0] != 't' || line[1] != 'i' || line[2] !='m' || line[3] != 'e')
			{
				throw " the format in raindata file is wrong, check the time";
			}
			else {
				std::string line_time, line_intensity;
				std::getline(infile, line_time);      //time 				
				std::istringstream iss_time(line_time);

				std::getline(infile, line_intensity); //标识行
				if (line_intensity[0] != 'i' || line_intensity[1] != 'n' || line_intensity[2] != 't' || line_intensity[3] != 'e')
				{
					throw "the format in raindata file is wrong ,check the inrensity";
				}
				else {
					std::getline(infile, line_intensity);    //数据行
					std::istringstream iss_iten(line_intensity);
					double  t = 0., inten = 0.;
					rainmap.push_back({});
					while (iss_time >> t && iss_iten >> inten)
					{
						rainmap[index].insert(std::make_pair(t, inten));   //添加数据
					}
				}
			}
			index++;
		}		
	}
	close();
}

bool Input::open(std::string &file)
{
	close();
	infile.open(file);
	if (!infile)
		throw "cannot open the input file";
	file_opened = true;
	return true;
}
void Input::close()
{
	if (file_opened)
		infile.close();
}