#include "Integral.h"

Integral::Integral(Integration_Scheme_Type Type)
{
	switch (Type)
	{
		//схема метода Гаусс-1
	case Gauss1:
	{
		Sigma = { 0 };
		Alfa = { 2 };	
		break;
	}
	case Simpson:
	{
		Sigma = { 0,
				  1,
				 -1 };
		Alfa = { 4. / 3,
				 1. / 3,
				 1. / 3 };
		break;
	}
	case Gauss4:
	{
		Sigma = { sqrt((3 - 2 * sqrt(6. / 5)) / 7) ,
				 -sqrt((3 - 2 * sqrt(6. / 5)) / 7),
				  sqrt((3 + 2 * sqrt(6. / 5)) / 7),
				 -sqrt((3 + 2 * sqrt(6. / 5)) / 7) };
		Alfa = { (18 + sqrt(30)) / 36 ,
				 (18 + sqrt(30)) / 36,
				 (18 - sqrt(30)) / 36,
				 (18 - sqrt(30)) / 36 };
		break;
	}
	// gauss5
	default:
	{
		Sigma = { 0,
				  sqrt(5 - 2 * sqrt(10. / 7)) / 3 ,
				 -sqrt(5 - 2 * sqrt(10. / 7)) / 3,
				  sqrt(5 + 2 * sqrt(10. / 7)) / 3,
				 -sqrt(5 + 2 * sqrt(10. / 7)) / 3 };
		Alfa = { 128. / 225,
				 (322 + 13 * sqrt(70)) / 900 ,
				 (322 + 13 * sqrt(70)) / 900,
				 (322 - 13 * sqrt(70)) / 900,
				 (322 - 13 * sqrt(70)) / 900 };
		break;
	}
	}
}

double Integral::Integration(const function<double(const double& P)>& F, const double& Start, const double& End, unsigned int Number_Nodes, double multiplier)
{
	if (Number_Nodes < 2 || multiplier < 1 || Start <= End)
		exception("Incorrect Integration parametrs");
	//результат (квадратурная сумма)
	double Result = 0.0,temp;
	vector <double> h;
	vector <double> x = Fragmentation::adaptive_fragmentation(Start, End,h, Number_Nodes, multiplier);
	//сумма по всем сегментам разбиения
	for (int i = 0; i < Number_Nodes-1; i++)
	{
		//сумма по узлам интегрирования
		for (int Integ_Point = 0; Integ_Point < Sigma.size(); Integ_Point++)
		{
			//переход с мастер-элемента [-1, 1]
			temp = x[i] + (1 + Sigma[Integ_Point]) * h[i] / 2.0;
			Result += Alfa[Integ_Point] * F(temp)*h[i];
		}
		
	}
	//формируем результат с учётом якобиана на отрезке [-1, 1]
	return Result/2.0;
}	
