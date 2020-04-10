#pragma once
#include <exception>
#include <vector>
#include <algorithm>
using namespace std;
class Fragmentation
{
public:
static vector <double> regular_fragmentation(double start,double end,unsigned int amount);
static	vector <double> adaptive_fragmentation(double start, double end, unsigned int amount, double multiplier);

};