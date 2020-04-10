#include "Cubic_Spline.h"

void Cubic_Spline::readfile(string File, string Path)
{
	ifstream in(Path + File);
	if (!in.is_open())
		throw exception("file is not open");
	int x, y;
	while (!in.eof())
	{
		in >> x >> y;
		this->x.push_back(x);
		this->y.push_back(y);
	}
	in.close();
}

void Cubic_Spline::fill(const vector<double>& x,const vector<double>& y)
{
	this->x = x;
	this->y = y;
}

void Cubic_Spline::fill(double f(double),const vector<double> &x)
{
	this->x = x;
	for (const auto i : x)
		y.push_back(f(i));
}

void Cubic_Spline::fill(double f(double), double start, double end, unsigned int amount, double multiplier)
{
	x=Fragmentation::adaptive_fragmentation(start, end, amount, multiplier);
	for (const auto i : x)
		y.push_back(f(i));
}

void Cubic_Spline::create_coef()
{
	int k = 0;
	int N = x.size();
	h.resize(N, 0);
	l.resize(N, 0);
	delta.resize(N, 0);
	lambda.resize(N, 0);
	c.resize(N, 0);
	d.resize(N, 0);
	b.resize(N, 0);
	for (k = 1; k < N; k++) {
		h[k] = x[k] - x[k - 1];
		if (h[k] == 0) {
			throw exception("x[k]==x[k-1]");
		}
		l[k] = (y[k] - y[k - 1]) / h[k];
	}
	delta[1] = -h[2] / (2 * (h[1] + h[2]));
	lambda[1] = 1.5 * (l[2] - l[1]) / (h[1] + h[2]);
	for (k = 3; k < N; k++) {
		delta[k - 1] = -h[k] / (2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]);
		lambda[k - 1] = (3 * l[k] - 3 * l[k - 1] - h[k - 1] * lambda[k - 2]) /
			(2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]);
	}
	for (k = N-1; k >= 2; k--) {
		c[k - 1] = delta[k - 1] * c[k] + lambda[k - 1];
	}
	// ñ[0]==0 c[N-1]==0
	for (k = 1; k < N; k++) {
		d[k] = (c[k] - c[k - 1]) / (3 * h[k]);
		b[k] = l[k] + (2 * c[k] * h[k] + h[k] * c[k - 1]) / 3;
	}
}

void Cubic_Spline::clear_coef()
{
	h.clear();
	l.clear();
	delta.clear();
	lambda.clear();
	c.clear();
	d.clear();
	b.clear();
}

vector <double> Cubic_Spline::interpolate(double point, vector<double> x, vector<double> y)
{
	if (!x.empty() || !y.empty())
	{this->x = x; this->y = y; create_coef();}
	if (this->x.empty() || this->y.empty())
		throw exception("x or y empty");
	if (h.empty())
		create_coef();
	int k;
	for (k = 1; k < this->x.size(); k++)
		if (point >= this->x[k - 1] && point <= this->x[k])
		{
			vector <double> res(3);
			res[0] = this->y[k] + b[k] * (point - this->x[k]) + c[k] * pow(point - this->x[k], 2) + d[k] * pow(point - this->x[k], 3);
			res[1] = b[k] + 2 * c[k] * (point - this->x[k]) + 3 * d[k] * pow(point - this->x[k], 2);
			res[2] = 2 * c[k] + 6 * d[k] * (point - this->x[k]);
			return res;

		}
	throw exception("point does not belong to the specified interval");
	
}

vector<vector<double>> Cubic_Spline::interpolate(vector<double> point)
{
	if (this->x.empty() || this->y.empty())
		throw exception("x or y empty");
	if (h.empty())
		create_coef();
	int k;
	vector <vector<double>> res(point.size(),vector<double>(3));
	for (int i = 0; i < point.size(); i++)
		res[i] = Cubic_Spline::interpolate(point[i]);
	return res;
}
