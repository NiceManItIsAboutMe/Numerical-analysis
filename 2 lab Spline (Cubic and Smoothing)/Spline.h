#pragma once
#include "Net_Fragmentation.h"
#include <fstream>
class Spline
{
public :
	virtual void fill(const vector <double>& x, const vector <double>& y)=0;
	virtual void fill(double f(double), const vector <double>& x)=0;
	virtual void fill(double f(double), double start, double end, unsigned int amount, double multiplier = 1)=0;
	virtual void create_coef() = 0;
	virtual void clear_coef() = 0;
	virtual vector <double> interpolate(double point, vector<double> x = vector<double>(), vector <double> y = vector <double>())=0;
	virtual vector<vector<double>> interpolate(vector<double> point)=0;
};