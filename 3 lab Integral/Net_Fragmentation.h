#pragma once
#include <exception>
#include <vector>
#include <algorithm>
using namespace std;
class Fragmentation
{
public:
static vector <double> regular_fragmentation(double start,double end,unsigned int amount);
static vector <double> regular_fragmentation(double start, double end, unsigned int amount,vector<double> &step);
static vector <double> regular_fragmentation(double start, double end, double step);
static	vector <double> adaptive_fragmentation(double start, double end, unsigned int amount, double multiplier);
static	vector <double> adaptive_fragmentation(double start, double end,vector<double> &step, unsigned int amount, double multiplier=1);

};