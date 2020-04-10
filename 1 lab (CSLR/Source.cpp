#include "Sparse_Matrix_CSlR.h"
#include <chrono>
int main()
{
	setlocale(LC_ALL, "Russian");
	try {
		CSLR a("\\фпми\\численные методы\\лаба 1 CSIR\\лаба 1 CSIR\\Matrix1\\");
		/*for (const auto& it : a.iptr)
			cout << it << "\t";
		cout << endl;*/
		a.in_tight_format();
		auto t1 = chrono::high_resolution_clock::now();
		a.mult_on_vector(a.di);
		auto t2 = chrono::high_resolution_clock::now();
		double time = chrono::duration <double>(t2 - t1).count();
		auto t3 = chrono::high_resolution_clock::now();
		a.mult_on_vector_tight_form(a.di);
		auto t4 = chrono::high_resolution_clock::now();
		double time2 = chrono::duration <double>(t4 - t3).count();
		cout << "CSLR:\t" << time << endl << "tight form:\t" << time2 << endl;
		/*if (a.check_vector(a.mult_on_vector(a.di), a.mult_on_vector_tight_form(a.di)))
			cout << "ok2\n";
		if (fabs(a.tight_format[1][5] - a.at(1, 5)) < EPS)
			cout << "ok3\n";*/
		int temp = 1;
		double y2;
		vector <double> Y;
		vector <double> X(a.di.size(), 0);
		X[0] = 1;
		// 5.1

		/*while (true)
		{
			y2 = 0;
			Y = a.mult_on_vector(X);
			for (auto it : Y)
				y2 += it * it;
			y2 = sqrt(y2);
			if (y2 > DBL_MAX)
				throw exception("y2 Out of range");
			if (X[0]> DBL_MAX)
				throw exception("X[0] Out of range");
			cout << "итерация N " << temp << " y2=" << y2 << " x1=" << X[0]<<endl;
			temp++;
			X[0] *= 10;
		}*/

		//5.2
		while (true)
		{
			y2 = 0;
			Y = a.mult_on_vector(X);
			for (auto it : Y)
				y2 += it * it;
			y2 = sqrt(y2);
			if (y2 ==0)
				throw exception("y2 Underflow");
			if (X[0] ==0)
				throw exception("X[0] Underflow");
			cout << "итерация N " << temp << " y2=" << y2 << " x1=" << X[0] << endl;
			temp++;
			X[0] /= 10;
		}
	}
	catch (const std::exception & err)
	{
		cout << err.what() << endl;
	};
}