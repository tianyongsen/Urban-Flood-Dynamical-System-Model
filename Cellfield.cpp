#include "Cellfield.h"
#include"Default_Parameters.h"
#include"Input.h"
#include<iostream>

void Cellfield::construct_Cellfiled(std::string &ele_file, std::string &roughness_file,std::string& rain_file)
{
	//I construct they, I manage they.  
	_input.init(ele_file, roughness_file, rain_file);
	construct_cells();       
	construct_rains();	
	construct_links();
}

void Cellfield::construct_cells()
{
	Matrix  ele_matrix = {};
	Matrix  roughness_matrix = {};
	_input.get_Celldata(ele_matrix, roughness_matrix);
	const int row = ele_matrix.size();     
	const int col = ele_matrix[0].size(); 
	_cells.reserve(row);                   //pre apply the memory 
	for (int i = 0; i < row; i++)     
	{
		_cells.push_back({});  
		_cells[i].reserve(col); 
		for (int j = 0; j < col; j++)
		{
			double n = 0.;      
			//give the given manning's coef based on the roughness map 
			if (roughness_matrix[i][j] == 0)
				n = Default_Parameters::Default_general_n;
			else
				n = Default_Parameters::Default_road_n;				
			Cell* cell = nullptr;
			//new ,,,,,remember to release
			cell = new Cell(ele_matrix[i][j],Default_Parameters::Default_h0, Default_Parameters::Default_area,n);
			_cells[i].push_back(cell);
		}		
   }
	_rows = _cells.size();
	_cols = _cells[0].size();
}

void Cellfield::construct_rains()
{	
	std::vector<std::map<double, double>>  rainmap = {};
	std::vector<std::vector<std::pair<int,int>>>      rain_aimCell_index= {};
	_input.get_Raindata(rain_aimCell_index,rainmap);                          //read 
	for (int i=0;i< rain_aimCell_index.size();i++)
	{
		Rain* rain(new Rain(rain_aimCell_index[i], rainmap[i]));  //construct ,remember to release
		_rains.push_back(rain);
	}                
	
	
	
}

void Cellfield::construct_links()
{
	int  maxcols = _cols+_cols-1;     //每行存储的边数  

	for (int i = 0; i < _rows-1; i++)
	{
		_links.push_back({});                  //申请行
		_links[i].reserve(maxcols);         //申请空间
		for (int j = 0; j < maxcols; j++)
		{
			_links[i].push_back(0.);           //初始化为0
		}		
	}
	//最后一行特殊处理
	_links.push_back({});
	for (int j = 0; j < _cols - 1; j++)    //最后一行仅有这些元素
	{
		_links[_rows - 1].push_back(0.);  
	}
}

//void Cellfield::construct_emap()
//{	
//	if (_rows == 1)     //如果只有一行  
//	{
//		for (int j = 0; j < _cols-1; j++)
//		{
//			_emap.insert(std::make_pair(Four_tuple(0, j,0, j + 1), 0.));
//		}
//		return;
//	}
//	if (_cols == 1)         //如果只有一列
//	{
//		for (int i = 0; i < _rows - 1; i++)
//		{
//			_emap.insert(std::make_pair(Four_tuple(i, 0, i+1, 0), 0.));   
//		}
//		return;
//	}	
//	//多行多列 
//	for (int i = 0; i < _rows; i++)
//	{
//		for (int j = 0; j < _cols; j++)
//		{
//			if (i == 0)
//			{
//				if (j != _cols - 1)
//				{
//					_emap.insert(std::make_pair(Four_tuple(i, j, i, j + 1), 0.));
//					_emap.insert(std::make_pair(Four_tuple(i, j, i+1, j), 0.));
//				}
//				else 
//					_emap.insert(std::make_pair(Four_tuple(i, j, i+1, j), 0.));
//			}
//			else if (i == _rows - 1)
//			{
//				if (j != _cols - 1)	
//					_emap.insert(std::make_pair(Four_tuple(i, j, i, j + 1), 0.));   
//			}
//			else   //平凡情况
//			{
//				if (j != _cols - 1)
//				{
//					_emap.insert(std::make_pair(Four_tuple(i, j, i, j + 1), 0.));
//					_emap.insert(std::make_pair(Four_tuple(i, j, i+1, j ), 0.));
//				}
//				else
//					_emap.insert(std::make_pair(Four_tuple(i, j, i+1, j), 0.));					
//			}				
//		}
//	}
//	
//}

Cellfield::~Cellfield()
{
	//delete all elements already applied
	for (auto &row : _cells)
	{
		for (auto &cell : row)
		{
			delete cell;      //correspoing to new cell 
			cell = nullptr;
		}
	}

	for (auto &rain : _rains)
	{
		delete rain;     
		rain = nullptr;
	}
	_cells.clear();
	_rains.clear();
}

//能不能定义一种循环结构，然后每个循环结构做一些操作。
