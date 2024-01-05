#include "Cell.h"
Cell::Cell(const double& elevation, const double& h0, const double& area, const double& roughness):
	_ele(elevation),_ht(h0), _area(area),_n(roughness)
{   }