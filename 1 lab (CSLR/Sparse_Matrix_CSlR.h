#pragma once
#define EPS 1e-16
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <exception>
using namespace std;

class CSLR
{
	//private:
public:
	vector <vector<double>> tight_format;
	vector <int>  iptr, jptr;
	vector <double> altr, autr, di;

public:
	CSLR(const std::string& PATH);
	vector <double> mult_on_vector(vector<double> &V);
	vector <vector <double>> in_tight_format();
	vector <double> mult_on_vector_tight_form(vector<double>& V,vector<vector<double>> matrix = vector<vector<double>>());
	double at(int Row, int Column);
	bool check_vector(const vector <double>& V1, const vector <double>& V2);
	bool check_matrix(const vector <vector<double>>& V1, const vector <vector<double>>& V2);
};