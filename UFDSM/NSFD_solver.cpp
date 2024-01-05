#include<cmath>
#include "NSFD_solver.h"

//debug
#include<iostream>

void NSFD_solver::run_NSFD_solver(std::vector<int>& report_time, std::vector<std::pair<int, int>>& report_index,
	std::vector<std::vector<double>>& report_h)
{
	//���ݸ�����������ʱ�䲽��
	_pt = (1 - std::exp(-1 * _q * _tau)) / _q;  //ʱ�䲽��(1 - e ^ (-_qh)) / _q;
	_ct = 0.;        //��ǰʱ��
	//����ѭ��ʱ������ pythonѭ��˼ά	
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
		NSFD_advance_step();      //�ƽ�һ������ 		
		_ct = _ct + _tau;                                    //����ʱ��		
	}
}

void NSFD_solver::run_NSFD_solver(ReportOutput& report)
{
	//���ݸ�����������ʱ�䲽��
	_pt = (1 - std::exp(-1 * _q * _tau)) / _q;  //ʱ�䲽��(1 - e ^ (-_qh)) / _q;
	_ct = 0.;        //��ǰʱ��
	//����ѭ��ʱ������ pythonѭ��˼ά	
	const double end_time = _end_time + _pt;    //������\phi(dt) ������dt  
	  
	int num = 0;             //��¼�ǵڼ��ε���
	while (_ct < end_time)
	{
		//����
		if(report.if_report(num))       //�ж��Ƿ�Ӧ�ñ��� 
		{
			NSFD_reserve_result(report);   //������  		
			std::cout << _ct << std::endl;   //�������ʱ�� 
		}
		//����
		NSFD_advance_step();      //�ƽ�һ������ 		
		num++;                         //num=0 ��_ct=0, num=nʱ����_ct=n*_pt���� num=report_intvalʱ�������ʵ��ʱ��Ϊreport_intval*_pt.
		_ct = _ct + _pt;                                    //����ʱ��		
		
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
	//����״̬


	//����������Ԫ���Ľ��� 
	double rain = 0.;
	for (auto& rainclass : _cellfield._rains)
	{
		if (rainclass->is_forAllCells())
			rain += rainclass->get_presentRain();
	}		
	//�����ܳ���ͨ������ֵ������Ҫ�Ӹ���
	cal_sum_g();
//======================================
//����Ԫ��״̬������������Դ�����Դ���������Ԫ���Ķ���
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
	//�������Ԫ����Դ��
	//====================================
	for (auto& rainclass : _cellfield._rains)
	{
		if (!rainclass->is_forAllCells())  //���Ƕ�����Ԫ��
		{
			for (auto& aim : rainclass->_aimcell)
			{
				_cellfield._cells[aim.first][aim.second]->_ht+= rainclass->get_presentRain();   //����Դ������ 
			}
		}				
	}
}
void NSFD_solver::cal_gij_generalCells()
{		
	// loop general cells
	const int rows = _cellfield._rows;
	const int cols = _cellfield._cols;
	auto&  cells = _cellfield._cells;           //��������
	//cal aij for every link (just cal once )
	cal_link_aij();		
	//cal gij for every cell ( based on the aij on every link)
	auto& links = _cellfield._links;   //��������

	double outflow = 0.;
	double a1 = 0.;
	double a2 = 0.;
	double a3 = 0.;
	double a4 = 0.;
	//�Ѿ����������aij, �����վ���ṹ˳�����������������aij����gij�����ڵ���0�ıߣ���������...
	//...ÿ�������ҽ�������һ��...
	//���Ǳ����ʺ��ж����Ρ�������ξͲ��������ˡ�
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			//�Ѿ��涨������������������ͬΪ�������������෴Ϊ��ֵ
			outflow = 0.;   //����Ϊ0
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
						outflow -= a1;    //��Ԫ��i ���Գ���Ϊ��ֵ
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
					if (a1 <0)    //��(i,j)�����ǳ���   �ԣ�i,i-1������������
						outflow -= -a1;      
					if (a2 >0)   
						outflow -= a2;
					if (a1 <0)
						links[i][j - 1] = cells[i][j]->_ht * a1 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);   //����ֵ
					if (a2 > 0)
						links[i][j + cols - 1] = cells[i][j]->_ht * a2 / (cells[i][j]->_ht * cells[i][j]->_area - _pt * outflow);  //����ֵ
				}
				else
				{
					a1 = links[i][j - 1];
					a2 = links[i][j];
					a3 = links[i][j +cols-1];
					if (a1 <0)    //��(i,j)�����ǳ���   �ԣ�i,i-1������������
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
					if (a1 < 0)    //��(i,j)�����ǳ���   �ԣ�i,i-1������������
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
					if (a1 < 0)    //��(i,j)�����ǳ���   �ԣ�i,i-1������������
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
					if (a1 < 0)    //��(i,j)�����ǳ���   �ԣ�i,i-1������������
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
					if (a1 > 0)    //��(i,j)�����ǳ���   �ԣ�i,i-1������������
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
	//�� rains���� ���������� 
	for (auto& rain:_cellfield._rains)
	{
		rain->cal_presentRain(_ct, _tau);   //�õ����� ,������rain ����
	}
}
void NSFD_solver::cal_link_aij()
{	
	auto& links = _cellfield._links;   //��������
	//int  row = _cellfield._links.size();
	int  col = _cellfield._cols;    //����
	
	//ÿ��Ԫ��ֻ�������ҵ�����߽�����µ�����߽� ���ұ߽�����±߽����⴦�� 
	for (int i = 0; i < _cellfield._rows - 1; i++)
	{
		for (int j = 0; j < _cellfield._cols-1; j++)   //�Է����Һ�����Ԫ����һ�㴦��
		{	
			auto& cell1 = _cellfield._cells[i][j];    //(i,j)
			auto& cell2 = _cellfield._cells[i][j + 1];   //���ҵ�����
			auto& cell3 = _cellfield._cells[i+1][j];    //���µ�����
			//����Ԫ�������ŵ��߾����ӳ���ϵ�����Ӧ�ߵ�ͨ��
			links[i][j]= Rules_set::rule_manning_formula(*cell1, *cell2); //����  i,j--> i,j+1
			links[i][j + col - 1] = Rules_set::rule_manning_formula(*cell1, *cell3);   //����i,j---> i+1,j
		}
	}
	//�����ұ߽����⴦��
	for (int i = 0; i < _cellfield._rows - 1; i++)
	{
		auto& cell1 = _cellfield._cells[i][col-1];     //���ұ߽�
		auto& cell2 = _cellfield._cells[i+1][col-1];  // ������ 
		links[i][col-1 + col - 1] = Rules_set::rule_manning_formula(*cell1, *cell2);  //����i,j---> i+1,j
	}
	//�����±߽����⴦��
	for (int j = 0; j < col-1;j++)
	{
		auto& cell1 = _cellfield._cells[_cellfield._rows-1][j];     //���ұ߽�
		auto& cell2 = _cellfield._cells[_cellfield._rows - 1][j+1];  // ������ 
		links[_cellfield._rows-1][j] = Rules_set::rule_manning_formula(*cell1, *cell2);  //����i,j---> i+1,j
	}	
}	
void NSFD_solver::NSFD_reserve_result(std::vector<std::pair<int, int>>& report_index, std::vector<std::vector<double>>& report_h)
{
	report_h.push_back({});                                //������һ��
	report_h.back().reserve(report_index.size());    //Ԥ���ռ�
	for (auto& id : report_index)
	{
		report_h.back().push_back(_cellfield._cells[id.first][id.second]->_ht);   //�洢��ʱ������� 
	}
	//���潵���Դ����Ϣ
	for (auto& rain : _cellfield._rains)
	{
		report_h.back().push_back(rain->get_presentRain());     //��λm
		report_h.back().push_back(rain->get_sumRain());     //��λmm
	}
	
	//����Դ����Ϣ
}
void NSFD_solver::NSFD_reserve_result(ReportOutput& report)
{
	auto& cells = _cellfield._cells;
	if (report.if_Allcells_report())        //�ж��Ƿ񱨸�����Ԫ��
	{
		const int rows = _cellfield._rows;
		const int cols = _cellfield._cols;		
		report.report_h.resize(rows);                //�����пռ�
		for(auto& re: report.report_h)
			re.resize(cols);									//�����пռ�
		for (int i=0;i<_cellfield._rows;i++)
		{			
			for (int j = 0; j < _cellfield._cols; j++)
			{				
				report.report_h[i][j] = cells[i][j]->_ht;
			}			
		}

		//����������txt�ļ���
		report.writeResult2file_every();   //��ȫ��Ԫ�����Ļ���ÿ��ʱ��һ���ļ� 
	}
	else
	{
		report.report_h.resize(1);                //�����пռ䣬����һ��
		report.report_h.resize(report.get_report_index().size());									//�����пռ�
		int k = 0;
		for (auto& id : report.get_report_index())
		{
			report.report_h[0][k] = _cellfield._cells[id.first][id.second]->_ht;   //�洢��ʱ������� 
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
	auto& cells = _cellfield._cells;           //��������
	auto& links = _cellfield._links;         //��������
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			//�������ͨ������
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
			
			cells[i][j]->_sum_g = sum_g;     //��¼��λʱ���ڵĳ�����ȣ���h�ĵ�λʱ��仯��
		}
	}


}
