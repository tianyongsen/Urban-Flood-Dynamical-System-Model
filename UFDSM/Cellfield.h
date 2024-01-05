#pragma once
#include<vector>
#include<string>
#include<unordered_map>
#include<map>
#include <time.h>
#include"Cell.h"
#include"Input.h"
#include"Rain.h"
//#include"Four_tuple.h"
//#include"NSFD_solver.h"
typedef  std::vector<std::vector<double>>  Matrix;
struct Four_tuple;
size_t Four_tuple_hash(const Four_tuple& ft);	
class NSFD_solver;     //for  友元
//class Cell;
class Cellfield
{
private:
	std::vector<std::vector<Cell*>>          _cells;   //只考虑最简单的矩阵元胞域，邻域关于服从矩阵行列规则	
	int							            		         _rows=0;  //行数
	int														_cols = 0;  //列数
	std::vector<Rain*>                               _rains;    //降雨  
	Input												     _input;	

	////元胞边界编号 和边上的通量 。 有特殊约定：边界编号由低值指向高值...
	//... 矩阵元胞区域，具有自然整齐结构，这种结构与向量（矩阵）结构对应，可以充分利用这种对应关系...
	//...简化算法逻辑和计算复杂度。这种整齐结构和对应关系也会被NSFD求解器继承。
	std::vector<std::vector<double> >      _links ;     //这里面蕴涵了对矩阵边界的特殊映射 
	//std::unordered_map<const Four_tuple, double , hash_name> _emap;   //用这个觉得程序太慢，试试直接用矩阵

	void construct_cells();
	void construct_rains();
	void construct_links();
	//void construct_emap(); 


public:
	friend NSFD_solver;    //暂时不考虑严格封装，先用友元取得便利
	void construct_Cellfiled(std::string& ele_file, std::string& roughness_file, std::string& rain_file);	
	~Cellfield();	
};







