#pragma once
#include "Net_Fragmentation.h"
#include <functional>
#include <cmath>

class Integral
{

private:
	vector <double> Sigma;
	vector <double> Alfa;
public:
    enum Integration_Scheme_Type
    {
        // Rectangle method
        Gauss1 = 1,
        // Parabola method
        Simpson = 2,
        Gauss4 = 4,
        Gauss5 = 5,
    };
	Integral(Integration_Scheme_Type Type);
	double Integration(const function<double(const double& P)>& F, const double& Start, const double& End, unsigned int Number_Nodes = 2, double multiplier = 1);
	
};
