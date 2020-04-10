#pragma once
#include "Spline.h"
class Smoothing_Spline :public Spline
{
private:
	//параметр сглаживания
	double SMOOTH=0;
	//точки для сплайна
	vector<double> x;
	vector<double> y;
	//коэффициенты разложения по базису
	vector<double> alpha;
	//переход на мастер элемент [-1, 1]: 
	//Seg_Num - номер сегмента, Х - абсцисса, Ksi - координата на мастер-элементе
	void Transition_To_Master_Element(int Seg_Num, const double& X, double& Ksi);
	//базисные функции на [-1, 1]:
	//Number - номер функции, Ksi - координата на мастер-элементе
	double Basis_Function(int Number, const double& Ksi);
	//производные базисных функций на [-1, 1]:
	//Number - номер функции, Ksi - координата на мастер-элементе
	double Der_Basis_Function(int Number, const double& Ksi);
public:
	void fill(const vector <double>& x, const vector <double>& y) override;
	void fill(double f(double), const vector <double>& x) override;
	void fill(double f(double), double start, double end, unsigned int amount, double multiplier = 1) override;
	void Smooth(double SMOOTH);
	void create_coef() override;
	void clear_coef() override;
	vector <double> interpolate(double point, vector<double> x = vector<double>(), vector <double> y = vector <double>()) override;
	vector<vector<double>> interpolate(vector<double> point) override;
};