#include "Cubic_Spline.h"
#include "Smoothing_Spline.h"
#include <iostream>
double X(double x)
{
	return x;
}
double square_x(double x)
{
	return x * x;
}
double cube_x(double x)
{
	return x * x * x;
}
double quadruple_x(double x)
{
	return x * x * x * x;
}
double foo(double x)
{
	return x * sin(10000 * x);
}
int main()
{
	try
	{
		Fragmentation a;
		// создание таблицы значений 
		vector<double> v1 = a.adaptive_fragmentation(0, 1, 10, 1);
		/*for (const auto i : v1)
			cout << i << endl;*/

		Cubic_Spline S;
		S.fill(X, v1);
		vector<double>point = { 0.04,0.12,0.26,0.35,0.41,0.52,0.67,0.79,0.84,0.94,0.97 };
		vector<vector<double>> result = S.interpolate(point);

		for (int i = 0; i < result.size(); i++)
		{
			cout << "x=\t" << point[i] << "\ty=\t" << result[i][0] << "\tdy/dx=\t" << result[i][1] << "\t\td^2y/dx^2=\t" << result[i][2] << endl;
		}
		/*
		Smoothing_Spline SP;
		SP.Smooth(0.0);
		SP.fill(foo, v1);
		vector<double>point2 = { -99,-90,-40,-30,-20,-10,-1,0,1,10,20,30,40,90,99 };
		vector<vector<double>> result2 = SP.interpolate(point2);
		cout << "/////////////////////////////////////////////////////////////////////////////////////////\n";
		for (int i = 0; i < result2.size(); i++)
		{
			cout << "x=\t" << point2[i] << "\ty=\t" << result2[i][0] << "\tdy/dx=\t" << result2[i][1] << "\t\td^2y/dx^2=\t" << result2[i][2] << endl;
		}
		/*
		ofstream fout("Result.txt");
		for (int i=0;i<result.size();i++)
		{
			fout << "x=\t" << point[i]<<"\ty=\t"<< result[i][0]<<"\tdy/dx=\t" << result[i][1]<<"\t\td^2y/dx^2=\t"<<result[i][2]<< endl;
		}
		fout.close();
		*/
	}
	catch (exception & err)
	{
		cout << err.what() << endl;
	}
}