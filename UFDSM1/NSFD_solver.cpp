#include<cmath>
#include "NSFD_solver.h"

//debug
#include<iostream>

void NSFD_solver::run_NSFD_solver(std::vector<int>& report_time, std::vector<std::pair<int, int>>& report_index,
	std::vector<std::vector<double>>& report_h)
{
	//根据给定步长计算时间步长
	_pt = (1 - std::exp(-1 * _q * _tau)) / _q;  //时间步长(1 - e ^ (-_qh)) / _q;
	_ct = 0.;        //当前时间
	//做个循环时间序列 python循环思维	
	int k = 0;
	const double end_time = _end_time + _tau;
	while(_ct<end_time)
	{
		if (k<report_time.size()&&std::abs(_ct - report_time[k])<0.1 * _tau)
		{
			NSFD_reserve_result(report_index, report_h);
			k++;
			std::cout <<_ct<<std::endl; 
		}
		NSFD_advance_step();      //推进一个步长 		
		_ct = _ct + _tau;                                    //更新时间		
	}
}

void NSFD_solver::run_NSFD_solver(ReportOutput& report)
{
	//根据给定步长计算时间步长
	_pt = (1 - std::exp(-1 * _q * _tau)) / _q;  //时间步长(1 - e ^ (-_qh)) / _q;
	_ct = 0.;        //当前时间
	//做个循环时间序列 python循环思维	
	const double end_time = _end_time + _pt;    //这里用\phi(dt) 而不是dt  
	  
	int num = 0;             //记录是第几次迭代
	while (_ct < end_time)
	{
		//报告
		if(report.if_report(num))       //判断是否应该报告 
		{
			NSFD_reserve_result(report);   //输出结果  		
			std::cout << _ct << std::endl;   //输出报告时刻 
		}
		//计算
		NSFD_advance_step();      //推进一个步长 		
		num++;                         //num=0 有_ct=0, num=n时，有_ct=n*_pt。当 num=report_intval时，报告的实际时刻为report_intval*_pt.
		_ct = _ct + _pt;                                    //更新时间		
		
	}

}

void  NSFD_solver::NSFD_advance_step()
{
	//step1:calculate  the outflow of the general cells  (the gij of general cells)		
	cal_gij_generalCells();
	//step2: calculate the outflow of the special cells 
	cal_specialCells();  
	//step3: upadate the next state  
	update(); 
	//some other work
}
void NSFD_solver::update()
{
	//更新状态


	//处理降到所有元胞的降雨 
	double rain = 0.;
	for (auto& rainclass : _cellfield._rains)
	{
		if (rainclass->is_forAllCells())
			rain += rainclass->get_presentRain();
	}		
	//计算总出流通量绝对值，后面要加负号
	cal_sum_g();
//======================================
//更新元胞状态，包括对流和源项。其中源项便是与雨元胞的对流
//======================================
	//int i = 0 ,j=0;
	for (auto& rowcell : _cellfield._cells)
	{
		//i++;
		//j = 0;
		for (auto& cell : rowcell)
		{
			cell->_ht = cell->_ht - _pt * cell->_sum_g + rain;		
		/*	j++;
			if (cell->_ht < 0.)
			{
				std::cout << "the ct is :" << _ct << i << "," << j <<"  the rain is "<<rain
					<<" the sumrain"<< _cellfield._rains[0]->get_sumRain() << std::endl;
				throw  "wrong";
			}*/
		}
		
	}
	//====================================
	//处理个别元胞的源项
	//====================================
	for (auto& rainclass : _cellfield._rains)
	{
		if (!rainclass->is_forAllCells())  //不是对所有元胞
		{
			for (auto& aim : rainclass->_aimcell)
			{
				_cellfield._cells[aim.first][aim.second]->_ht+= rainclass->get_presentRain();   //加入源项入流 
			}
		}				
	}
}
void NSFD_solver::cal_gij_generalCells()
{		
	// loop general cells
	const int rows = _cellfield._rows;
	const int cols = _cellfield._cols;
	auto&  cells = _cellfield._cells;           //换个名字
	//cal aij for every link (just cal once )
	cal_link_aij();		
	//cal gij for every cell ( based on the aij on every link)
	auto& links = _cellfield._links;   //换个名字

	double outflow = 0.;
	double a1 = 0.;
	double a2 = 0.;
	double a3 = 0.;
	double a4 = 0.;
	//已经计算出来了aij, 并按照矩阵结构顺序定了正负。下面基于aij计算gij，对于等于0的边，不做处理...
	//...每个边有且仅有修正一次...
	//但是被访问和判断两次。这个两次就不作处理了。
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			//已经规定了正方向，与正方向相同为正，与正方向相反为负值
			outflow = 0.;   //重置为0
			a1 = 0.;          
			a2 = 0.;
			a3 = 0.;
			a4 = 0.;
			if (i == 0)
			{
				if (j ==0 )
				{					
					a1 = links[i][j];                
					a2 = links[i][j+cols-1];   
					if (a1 >0)    //outflow 
						outflow -= a1;    //对元胞i 而言出流为负值
					if (a2 > 0)   //outflow
						outflow -= a2;
					if (a1 > 0)
						links[i][j] =cells[i][j]->_ht * a1 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if(a2>0)
						links[i][j + cols - 1] = cells[i][j]->_ht * a2 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
				}
				else if(j==cols-1)
				{						
					a1 = links[i][j- 1];
					a2 = links[i][j+cols-1];
					if (a1 <0)    //对(i,j)而言是出流   对（i,i-1）而言是入流
						outflow -= -a1;      
					if (a2 >0)   
						outflow -= a2;
					if (a1 <0)
						links[i][j - 1] = cells[i][j]->_ht * a1 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);   //修正值
					if (a2 > 0)
						links[i][j + cols - 1] = cells[i][j]->_ht * a2 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);  //修正值
				}
				else
				{
					a1 = links[i][j - 1];
					a2 = links[i][j];
					a3 = links[i][j +cols-1];
					if (a1 <0)    //对(i,j)而言是出流   对（i,i-1）而言是入流
						outflow -= -a1; 
					if (a2 >0)
						outflow -= a2;
					if (a3 > 0)
						outflow -= a3;
					if (a1 < 0)
						links[i][j - 1] =cells[i][j]->_ht * a1 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if (a2 > 0)
						links[i][j + 1] = cells[i][j]->_ht * a2 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if(a3>0)
						links[i][j + cols - 1] = cells[i][j]->_ht * a3 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
				}
			}
			else if (i == rows - 1)
			{
				if (j == 0)
				{
					a1 = links[i][j];     
					a2 = links[i-1][j+cols-1]; 
					if (a1 > 0)    //outflow 
						outflow -= a1;
					if (a2 < 0)   //
						outflow -= -a2;
					if (a1 > 0)
						links[i][j] = cells[i][j]->_ht * a1 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if (a2 < 0)
						links[i - 1][j] = cells[i][j]->_ht * a2 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
				}
				else if (j == cols - 1)
				{
					a1 = links[i][j-1];
					a2 = links[i-1][j+cols-1];
					if (a1 < 0)    //对(i,j)而言是出流   对（i,i-1）而言是入流
						outflow -= -a1;
					if (a2 < 0)
						outflow -= -a2;
					if (a1 < 0)
						links[i][j - 1] =cells[i][j]->_ht * a1 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if (a2 < 0)
						links[i - 1][j + cols - 1] = cells[i][j]->_ht * a2 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
				}
				else
				{
					a1 = links[i][j - 1];
					a2 = links[i][j];
					a3 = links[i-1][j+cols-1];
					if (a1 < 0)    //对(i,j)而言是出流   对（i,i-1）而言是入流
						outflow -= -a1;
					if (a2 > 0)
						outflow -= a2;
					if (a3 < 0)
						outflow -= -a3;
					if (a1 < 0)
						links[i][j - 1] = cells[i][j]->_ht * a1 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if (a2 > 0)
						links[i][j] = cells[i][j]->_ht * a2 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if (a3 < 0)
						links[i - 1][j + cols - 1] = cells[i][j]->_ht * a3 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
				}
			}
			else
			{
				if (j == 0)
				{
					a1 = links[i][j];
					a2 = links[i-1][j +cols- 1];
					a3 = links[i][j+cols-1];
					if (a1 > 0)    //outflow 
						outflow -= a1;
					if (a2 < 0)   //
						outflow -= -a2;
					if (a3 > 0)
						outflow -= a3;
					if (a1 > 0)
						links[i][j] = cells[i][j]->_ht * a1 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if (a2 < 0)
						links[i - 1][j + cols - 1] = cells[i][j]->_ht * a2 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if (a3 > 0)
						links[i][j + cols - 1] = cells[i][j]->_ht * a3 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
				}
				else if (j == cols - 1)
				{
					a1 = links[i][j-1];
					a2 = links[i-1][j+cols-1];
					a3 = links[i][j+cols-1];
					if (a1 < 0)    //对(i,j)而言是出流   对（i,i-1）而言是入流
						outflow -= -a1;
					if (a2 < 0)
						outflow -= -a2;
					if (a3 > 0)
						outflow -= a3;
					if (a1 < 0)
						links[i][j - 1] =cells[i][j]->_ht * a1 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if (a2 < 0)
						links[i - 1][j + cols - 1] = cells[i][j]->_ht * a2 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if (a3 > 0)
						links[i][j + cols - 1] = cells[i][j]->_ht * a3 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
				}
				else
				{
					a1 = links[i][j ];
					a2 = links[i][j - 1];
					a3 = links[i][j +cols- 1];
					a4 = links[i-1][j+cols - 1];
					if (a1 > 0)    //对(i,j)而言是出流   对（i,i-1）而言是入流
						outflow -= a1;
					if (a2 < 0)
						outflow -= -a2;
					if (a3 >0)
						outflow -= a3;
					if (a4 < 0)
						outflow -= -a4;
					if (a1 > 0)
						links[i][j] = cells[i][j]->_ht * a1 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if (a2 < 0)
						links[i][j - 1] = cells[i][j]->_ht * a2 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if (a3 > 0)
						links[i][j + cols - 1] = cells[i][j]->_ht * a3 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
					if (a4 < 0)
						links[i - 1][j + cols - 1] = cells[i][j]->_ht * a4 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);
				}

			}
		}
	}
}
void NSFD_solver::cal_specialCells()
{
	//对 rains遍历 计算入流量 
	for (auto& rain:_cellfield._rains)
	{
		rain->cal_presentRain(_ct, _tau);   //得到入流 ,存在了rain 里面
	}
}
void NSFD_solver::cal_link_aij()
{	
	auto& links = _cellfield._links;   //换个名字
	//int  row = _cellfield._links.size();
	int  col = _cellfield._cols;    //列数
	
	//每个元胞只计算向右的邻域边界和向下的邻域边界 最右边界和最下边界特殊处理 
	for (int i = 0; i < _cellfield._rows - 1; i++)
	{
		for (int j = 0; j < _cellfield._cols-1; j++)   //对非最右和最下元胞做一般处理
		{	
			auto& cell1 = _cellfield._cells[i][j];    //(i,j)
			auto& cell2 = _cellfield._cells[i][j + 1];   //向右的邻域
			auto& cell3 = _cellfield._cells[i+1][j];    //向下的邻域
			//根据元胞矩阵编号到边矩阵的映射关系计算对应边的通量
			links[i][j]= Rules_set::rule_manning_formula(*cell1, *cell2); //计算  i,j--> i,j+1
			links[i][j + col - 1] = Rules_set::rule_manning_formula(*cell1, *cell3);   //计算i,j---> i+1,j
		}
	}
	//对最右边界特殊处理
	for (int i = 0; i < _cellfield._rows - 1; i++)
	{
		auto& cell1 = _cellfield._cells[i][col-1];     //最右边界
		auto& cell2 = _cellfield._cells[i+1][col-1];  // 其邻域 
		links[i][col-1 + col - 1] = Rules_set::rule_manning_formula(*cell1, *cell2);  //计算i,j---> i+1,j
	}
	//对最下边界特殊处理
	for (int j = 0; j < col-1;j++)
	{
		auto& cell1 = _cellfield._cells[_cellfield._rows-1][j];     //最右边界
		auto& cell2 = _cellfield._cells[_cellfield._rows - 1][j+1];  // 其邻域 
		links[_cellfield._rows-1][j] = Rules_set::rule_manning_formula(*cell1, *cell2);  //计算i,j---> i+1,j
	}	
}	
void NSFD_solver::NSFD_reserve_result(std::vector<std::pair<int, int>>& report_index, std::vector<std::vector<double>>& report_h)
{
	report_h.push_back({});                                //申请下一行
	report_h.back().reserve(report_index.size());    //预留空间
	for (auto& id : report_index)
	{
		report_h.back().push_back(_cellfield._cells[id.first][id.second]->_ht);   //存储该时间点结果。 
	}
	//报告降雨和源项信息
	for (auto& rain : _cellfield._rains)
	{
		report_h.back().push_back(rain->get_presentRain());     //单位m
		report_h.back().push_back(rain->get_sumRain());     //单位mm
	}
	
	//报告源项信息
}
void NSFD_solver::NSFD_reserve_result(ReportOutput& report)
{
	auto& cells = _cellfield._cells;
	if (report.if_Allcells_report())        //判断是否报告所有元胞
	{
		const int rows = _cellfield._rows;
		const int cols = _cellfield._cols;		
		report.report_h.resize(rows);                //申请行空间
		for(auto& re: report.report_h)
			re.resize(cols);									//申请列空间
		for (int i=0;i<_cellfield._rows;i++)
		{			
			for (int j = 0; j < _cellfield._cols; j++)
			{				
				report.report_h[i][j] = cells[i][j]->_ht;
			}			
		}

		//将结果输出到txt文件中
		report.writeResult2file_every();   //对全部元胞场的话，每个时刻一个文件 
	}
	else
	{
		report.report_h.resize(1);                //申请行空间，仅有一个
		report.report_h.resize(report.get_report_index().size());									//申请列空间
		int k = 0;
		for (auto& id : report.get_report_index())
		{
			report.report_h[0][k] = _cellfield._cells[id.first][id.second]->_ht;   //存储该时间点结果。 
			k++;
		}
		report.writeResult2file_single();
	}
}


void NSFD_solver:: cal_sum_g()
{
	double sum_g = 0.;
	const int rows = _cellfield._rows;
	const int cols = _cellfield._cols;
	auto& cells = _cellfield._cells;           //换个名字
	auto& links = _cellfield._links;         //换个名字
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			//下面计算通量更新
			if (i == 0)
			{
				if (j == 0)
					sum_g = links[i][j] + links[i][j + cols - 1];
				else if (j == cols - 1)
					sum_g = -links[i][j - 1] + links[i][j + cols - 1];
				else
					sum_g = -links[i][j - 1] + links[i][j] + links[i][j + cols - 1];
			}
			else if (i == rows - 1)
			{
				if (j == 0)
					sum_g = links[i][j] - links[i - 1][j + cols - 1];
				else if (j == cols - 1)
					sum_g = -links[i][j - 1] - links[i - 1][j + cols - 1];
				else
					sum_g = links[i][j] - links[i][j - 1] - links[i - 1][j + cols - 1];
			}
			else
			{
				if (j == 0)
					sum_g = links[i][j] - links[i - 1][j + cols - 1] + links[i][j + cols - 1];
				else if (j == cols - 1)
					sum_g = -links[i][j - 1] - links[i - 1][j + cols - 1] + links[i][j + cols - 1];
				else
					sum_g = links[i][j] - links[i][j - 1] - links[i - 1][j + cols - 1] + links[i][j + cols - 1];
			}
			
			cells[i][j]->_sum_g = sum_g;     //记录单位时间内的出流深度，即h的单位时间变化率
		}
	}


}
