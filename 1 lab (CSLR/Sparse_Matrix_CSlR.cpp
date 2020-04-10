#include "Sparse_Matrix_CSlR.h"

CSLR::CSLR(const std::string& PATH) {
	//размер матрицы
	int N = 0;
	std::ifstream Reader(PATH + "size.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File size.bin was not found...");
	Reader.read((char*)&N, sizeof(N));
	Reader.close();
	//di
	Reader.open(PATH + "di.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File di.bin was not found...");
	di.resize(N);
	for (int i = 0; i < N; i++)
		Reader.read((char*)&di[i], sizeof(double));
	Reader.close();
	//iptr
	Reader.open(PATH + "iptr.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File iptr.bin was not found...");
	iptr.resize(N + 1);
	for (int i = 0; i < N + 1; i++)
		Reader.read((char*)&iptr[i], sizeof(int));
	Reader.close();
	//jptr
	Reader.open(PATH + "jptr.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File jptr.bin was not found...");
	int JPTR_SIZE = iptr[N] - 1;
	jptr.resize(JPTR_SIZE);
	for (int i = 0; i < JPTR_SIZE; i++)
		Reader.read((char*)&jptr[i], sizeof(int));
	Reader.close();
	//autr
	Reader.open(PATH + "autr.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File autr.bin was not found...");
	autr.resize(JPTR_SIZE);
	for (int i = 0; i < JPTR_SIZE; i++)
		Reader.read((char*)&autr[i], sizeof(int));
	Reader.close();
	// altr
	Reader.open(PATH + "altr.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File altr.bin was not found...");
	altr.resize(JPTR_SIZE);
	for (int i = 0; i < JPTR_SIZE; i++)
		Reader.read((char*)&altr[i], sizeof(int));
	Reader.close();


}

vector<double> CSLR::mult_on_vector(vector<double> &V)
{
	//размер матрицы
	int n = di.size();
	vector <double> Res(n);
	//инициализация результата через умножения вектора на диагональ
	for (int i = 0; i < n; i++) Res[i] = V[i] * di[i];
	//проход по всем строкам и столбцам с учётом формата
	for (int i = 0; i < n; i++)
		for (int j = iptr[i] - 1; j < iptr[i + 1] - 1; j++)
		{
			Res[i] += V[jptr[j] - 1] * altr[j];
			Res[jptr[j] - 1] += V[i] * autr[j];
		}
	return Res;

}

vector<vector<double>> CSLR::in_tight_format()
{
	int size = di.size();
	tight_format.resize(size);
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			tight_format[i].push_back(CSLR::at(i,j));

	return tight_format;
}

vector<double> CSLR::mult_on_vector_tight_form(vector<double> &V,vector<vector<double>> matrix)
{
	if (matrix == vector<vector<double>>())
		matrix = tight_format;
	int size = matrix.size();
	vector <double> result(size, 0);
	if (V.size() != size)
		throw exception("Vector size != length of the column");
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			result[i] += tight_format[i][j] * V[i];
	return result;
}



double CSLR::at(int Row, int Column)
{
	//диагональный элемент
	if (Row == Column) return di[Row];
	//нижний треугольник
	if (Row > Column)
	{
		//число элементов в строке
		int Num = iptr[Row + 1] - iptr[Row];
		if (Num == 0) return 0;
		//ищем столбец Column в массиве jptr
		for (int ind = iptr[Row] - 1, k = 0; k < Num; ind++, k++)
			if (Column == (jptr[ind] - 1)) return altr[iptr[Row] - 1 + k];
	}
	//верхний треугольник
	if (Row < Column)
	{
		//число элементов в столбце
		int Num = iptr[Column + 1] - iptr[Column];
		if (Num == 0) return 0;
		//ищем строку Row в массиве jptr
		for (int ind = iptr[Column] - 1, k = 0; k < Num; ind++, k++)
			if (Row == (jptr[ind] - 1)) return autr[iptr[Column] - 1 + k];
	}
}

bool CSLR::check_vector(const vector<double>& V1, const vector<double>& V2)
{
	if (V1.size() != V2.size())
		return false;
	for (int i = 0; i < V1.size(); i++)
		if (fabs(V1[i] - V2[i]) > EPS)
			return false;
	return true;
}

bool CSLR::check_matrix(const vector<vector<double>>& V1, const vector<vector<double>>& V2)
{
	if (V1.size() != V2.size())
		return false;
	for (int i = 0; i < V1.size(); i++)
	{
		if (V1[i].size() != V2[i].size())
			return false;
		for (int j = 0; j < V1.size(); j++)
			if (fabs(V1[i][j] - V2[i][j]) > EPS)
				return false;
	}
	return true;
}
