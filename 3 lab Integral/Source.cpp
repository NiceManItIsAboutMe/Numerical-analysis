#include "Integral.h"
#include <iostream>
int main()
{
	try {
		
		function<double(const double& x)> f =
			[](const double& x) { return sin(x) ; };
		
		function<double(const double& x)> F =
			[](const double& x) { return  -cos(x); };

		//квадратурная формула Гаусс-1	
		Integral Gauss1(Integral::Gauss1);
		Integral Gauss4(Integral::Gauss4);
		Integral Gauss5(Integral::Gauss5);
		Integral Simpson(Integral::Simpson);

		//начало и конец отрезка интегрирования
		double Begin = 0;
		double End = 1;

		//число узлов !!!
		const int Num_Nodes = 5; // на 1 больше чем сегментов
		const int Multiplier = 1;

		//точное значение интеграла (ф. Ньютона-Лейбница)
		double I_True = F(End) - F(Begin);

		//численное значение интеграла
		double resG1 = Gauss1.Integration(f, Begin, End, Num_Nodes,Multiplier );
		double resG1H2 = Gauss1.Integration(f, Begin, End, (Num_Nodes - 1) * 2 + 1, Multiplier);
		double resG4 = Gauss4.Integration(f, Begin, End, Num_Nodes, Multiplier);
		double resG4H2 = Gauss4.Integration(f, Begin, End, (Num_Nodes - 1) * 2 + 1, Multiplier);
		double resG5 = Gauss5.Integration(f, Begin, End, Num_Nodes, Multiplier);
		double resG5H2 = Gauss5.Integration(f, Begin, End, (Num_Nodes - 1) * 2 + 1, Multiplier);
		double resSimpson = Simpson.Integration(f, Begin, End, Num_Nodes, Multiplier);
		double resSimpsonH2 = Simpson.Integration(f, Begin, End, (Num_Nodes - 1) * 2 + 1, Multiplier);
		std::cout << std::scientific;
		std::cout << "h = " << (End - Begin) / (Num_Nodes -1) << std::endl;
		std::cout << "Gauss1 = " << resG1 << std::endl;
		std::cout << "Simpson = " << resSimpson << std::endl;
		std::cout << "Gauss4 = " << resG4 << std::endl;
		std::cout << "Gauss5 = " << resG5 << std::endl;
		cout << "I_TRUE = " << I_True << endl;
		std::cout << "|Gauss1 - I_True| = " << fabs(resG1 - I_True) << std::endl;	
		std::cout << "|GaussSimpson - I_True| = " << fabs(resSimpson - I_True) << std::endl;
		std::cout << "|Gauss4 - I_True| = " << fabs(resG4 - I_True) << std::endl;
		std::cout << "|Gauss5 - I_True| = " << fabs(resG5 - I_True) << std::endl;
		cout << "///////////////////////////////////////////" << endl;
		cout << "Gauss 1 k = " << log2(fabs(1 + ((resG1H2 - resG1) / (I_True - resG1H2)))) << endl;
		cout << "Simpson k = " << log2(fabs(1 + ((resSimpsonH2 - resSimpson) / (I_True - resSimpsonH2)))) << endl;
		cout << "Gauss 4 k = " << log2(fabs(1 + ((resG4H2 - resG4) / (I_True - resG4H2)))) << endl;
		cout << "Gauss 5 k = " << log2(fabs(1 + ((resG5H2 - resG5) / (I_True - resG5H2)))) << endl;
	
	}
	catch (exception & err)	
	{
		cout << err.what() << endl;
	}
}