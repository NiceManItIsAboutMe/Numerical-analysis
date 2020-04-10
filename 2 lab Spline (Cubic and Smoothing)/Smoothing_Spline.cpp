#include "Smoothing_Spline.h"
#include <functional>

void Smoothing_Spline::Transition_To_Master_Element(int Seg_Num, const double& X, double& Ksi) 
{
	Ksi = 2.0 * (X - x[Seg_Num]) / (x[Seg_Num + 1] - x[Seg_Num]) - 1.0;
}

double Smoothing_Spline::Basis_Function(int Number, const double& Ksi) 
{
	switch (Number)
	{
	case 1: {return 0.5 * (1 - Ksi); break; }
	case 2: {return 0.5 * (1 + Ksi); break; }
	default: {throw std::exception("Error in the basis function number..."); break; }
	}
}

double Smoothing_Spline::Der_Basis_Function(int Number, const double& Ksi) 
{
	switch (Number)
	{
	case 1: {return -0.5; break; }
	case 2: {return  0.5; break; }
	default: {throw std::exception("Error in the basis function derivative number..."); break; }
	}
}

void  Smoothing_Spline::fill(const vector<double>& x, const vector<double>& y)
{
	this->x = x;
	this->y = y;
}

void  Smoothing_Spline::fill(double f(double), const vector<double>& x)
{
	this->x = x;
	for (const auto i : x)
		y.push_back(f(i));
}

void  Smoothing_Spline::fill(double f(double), double start, double end, unsigned int amount, double multiplier)
{
	x = Fragmentation::adaptive_fragmentation(start, end, amount, multiplier);
	for (const auto i : x)
		y.push_back(f(i));
}

void Smoothing_Spline::Smooth(double SMOOTH)
{
	this->SMOOTH = SMOOTH;
}

void Smoothing_Spline::create_coef()
{
	//число отрезков разбиения
	int Num_Segments = x.size() - 1;

	//коэффициенты разложения по базису
	alpha.resize(Num_Segments + 1);

	//диагонали матрицы
	std::vector <double> a, b, c;
	a.resize(Num_Segments + 1);
	b.resize(Num_Segments + 1);
	c.resize(Num_Segments + 1);

	//процедура для ассемблирования СЛАУ: 
	//Num_Segment - номер отрезка, P - точка, F_Val - значение табличной функции в точке, w - вес  
	std::function<void(int Num_Segment, const double & X, const double& F_Val, const double& w)>
		Assembling = [&](int i, const double & X, const double& F_Val, const double& w)
	{
		double  Ksi;
		//переход на мастер-элемент
		Transition_To_Master_Element(i, X, Ksi);
		//вычисление значений базисных функций на мастер-элементе
		double f1 = Basis_Function(1, Ksi);
		double f2 = Basis_Function(2, Ksi);

		//внесение вкладов в СЛАУ
		b[i] += (1.0 - SMOOTH) * w * f1 * f1;
		b[i + 1] += (1.0 - SMOOTH) * w * f2 * f2;
		a[i + 1] += (1.0 - SMOOTH) * w * f1 * f2;
		c[i] += (1.0 - SMOOTH) * w * f2 * f1;
		alpha[i] += (1.0 - SMOOTH) * w * f1 * F_Val;
		alpha[i + 1] += (1.0 - SMOOTH) * w * f2 * F_Val;
	};

	//сборка СЛАУ по сетке: сумма вкладов от каждого сегмента разбиения
	for (int i = 0; i < Num_Segments; i++)
	{
		//добавление узла сетки в СЛАУ
		double W = 1.0;
		Assembling(i, this->x[i], y[i], W);
		Assembling(i, this->x[i + 1], y[i + 1], W);

		//вклад от сглаживания по первой производной
		double h = x[i + 1] - x[i];
		b[i] += 1.0 / h * SMOOTH;
		b[i + 1] += 1.0 / h * SMOOTH;
		a[i + 1] -= 1.0 / h * SMOOTH;
		c[i] -= 1.0 / h * SMOOTH;
	}

	//метод прогонки: прямой ход
	for (int j = 1; j < Num_Segments + 1; j++)
	{
		//диагональ
		b[j] -= a[j] / b[j - 1] * c[j - 1];
		//правая часть         
		alpha[j] -= a[j] / b[j - 1] * alpha[j - 1];
	}

	//метод прогонки: обратный ход
	alpha[Num_Segments] /= b[Num_Segments];
	for (int j = Num_Segments - 1; j >= 0; j--)
		alpha[j] = (alpha[j] - alpha[j + 1] * c[j]) / b[j];
}

void Smoothing_Spline::clear_coef()
{
}

vector<double> Smoothing_Spline::interpolate(double P, vector<double> x, vector<double> y)
{
	if (!x.empty() || !y.empty())
	{
		this->x = x; this->y = y; create_coef();
	}
	if (this->x.empty() || this->y.empty())
		throw exception("x or y empty");
	if (alpha.empty())
		create_coef();
	vector <double> Res(3);
	//машинный ноль
	double eps = 1e-7;
	//число отрезков
	int Num_Segments = this->x.size() - 1;
	//координата точки
	double X = P;

	//поиск отрезка, которому принадлежит точка
	for (int i = 0; i < Num_Segments; i++)
	{
		if (X > this->x[i] && X < this->x[i + 1] ||
			fabs(X - this->x[i]) < eps ||
			fabs(X - this->x[i + 1]) < eps)
		{
			//длина отрезка
			double h = this->x[i + 1]- this->x[i];
			//переход на местер-элемент, Ksi - координата на мастер-элементе
			double Ksi;
			Transition_To_Master_Element(i, X, Ksi);
			//вычисляем значение сплайна и производной по базисным функциям
			Res[0] = alpha[i] * Basis_Function(1, Ksi) +
				alpha[i + 1] * Basis_Function(2, Ksi);
			Res[1] = (alpha[i] * Der_Basis_Function(1, Ksi) +
				alpha[i + 1] * Der_Basis_Function(2, Ksi)) * 2.0 / h;
			Res[2] = 0.0;
			return Res;
		}
	}
	throw std::exception("The point wasn't found in the segments...");
}

vector<vector<double>> Smoothing_Spline::interpolate(vector<double> point)
{
	if (this->x.empty() || this->y.empty())
		throw exception("x or y empty");
	if (alpha.empty())
		create_coef();
	int k;
	vector <vector<double>> res(point.size(), vector<double>(3));
	for (int i = 0; i < point.size(); i++)
		res[i] = Smoothing_Spline::interpolate(point[i]);
	return res;
}
