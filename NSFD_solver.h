#pragma once
#include<vector>
#include<unordered_map>
#include"Cellfield.h"
#include"Rules_set.h"  
#include"Cell.h"
#include"ReportOutput.h"
//class Cell;
//class Cellfield;
//class Rules_set;
//struct Four_tuple;
	
class NSFD_solver
{
private:
	Cellfield&  _cellfield;          //元胞空间 
	Rules_set   _rules;              //规则空间
	//时间
	const double     _tau = 0.;   //标准步长
	const double     _q = 0.;
	const double     _end_time = 0.;
	double              _pt=0.;        //非标准步长 
	double              _ct=0.;        //当前时间 标准时间
	void NSFD_reserve_result(std::vector<std::pair<int, int>>& report_index,std::vector<std::vector<double>>& report_h);
	void NSFD_reserve_result(ReportOutput& report);
	void NSFD_advance_step();
	void cal_gij_generalCells();
	void cal_specialCells();    
	void update();		
	void cal_link_aij();	
	void  cal_sum_g(); 
public:
	NSFD_solver(Cellfield& cellfield_,const double& time_step,const double & end_time,const double q) :
		_cellfield(cellfield_), _tau(time_step),_end_time(end_time),_q(q) {};
	void run_NSFD_solver(std::vector<int>& report_time, std::vector<std::pair<int,int>>& report_index,
		std::vector<std::vector<double>> & report_h); 
	void run_NSFD_solver(ReportOutput& report);    //main 函数里面主调用程序 


};






