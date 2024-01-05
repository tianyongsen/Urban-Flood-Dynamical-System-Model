#include<cmath>
#include "Rules_set.h"

double Rules_set::rule_manning_formula(const Cell& cell1, const Cell& cell2)
{
	//已经规定了正方向，与正方向相同为正，与正方向相反为负值
	double H1 = cell1._ht + cell1._ele;
	double H2 = cell2._ht + cell2._ele;
	double n = (cell1._n + cell2._n)/2.;    //取平均   
	if (H1 == H2)
		return 0.;
	else if (H1 > H2)
	{
		if (cell1._ht <= 0)
			return 0;           //无出流  作截止
		return Default_Parameters::Default_length* std::pow(cell1._ht, 5. / 3.)* std::sqrt((H1 - H2) / Default_Parameters::Default_length)/n;
	}
	else
	{
		if (cell2._ht <= 0)
			return 0;           //无出流  作截止
		return -Default_Parameters::Default_length * std::pow(cell2._ht, 5. / 3.) * std::sqrt((H2 - H1) / Default_Parameters::Default_length) / n;
	}	
}