#pragma once
#include "Spline.h"
class Cubic_Spline: public Spline
{
private:
	vector <double> x;
	vector <double> y;
	vector <double> h; // h=x[k]-x[k-1]
	vector <double> l; // l=(y[k]-y[k-1])/h[k]
	vector <double> delta; //delta[1] = - h[2]/(2*(h[1]+h[2]))
	//delta[k-1] = - h[k]/(2*h[k-1] + 2*h[k] + h[k-1]*delta[k-2]);
	vector <double> lambda; //lambda[1] = 1.5*(l[2] - l[1])/(h[1]+h[2]);
	//  lambda[k-1] = (3*l[k] - 3*l[k-1] - h[k-1]*lambda[k-2]) / (2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]);
	vector <double> c;      //  c[k-1] = delta[k-1]*c[k] + lambda[k-1];
	vector <double> d; //d[k] = (c[k] - c[k - 1]) / (3 * h[k]);
	vector <double> b;  //b[k] = l[k] + (2 * c[k] * h[k] + h[k] * c[k - 1]) / 3

public:
	void readfile(string File,string Path="");
	void fill(const vector <double>& x,const vector <double>& y) override;
	void fill(double f(double), const vector <double> &x) override;
	void fill(double f(double), double start, double end, unsigned int amount, double multiplier=1) override;
	void create_coef() override;
	void clear_coef() override;
	vector <double> interpolate(double point,vector<double> x=vector<double> (),vector <double> y=vector <double> ()) override;
	vector<vector<double>> interpolate(vector<double> point) override;
};