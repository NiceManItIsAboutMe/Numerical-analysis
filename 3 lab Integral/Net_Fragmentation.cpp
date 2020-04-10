	#include "Net_Fragmentation.h"

	vector <double> Fragmentation::regular_fragmentation(double start, double end,unsigned int amount)
	{
		if ( end <= start || amount < 3||fabs(end-start)/(amount)==0)
			throw exception("incorrect input regular_fragmentation parameters");
		double step = fabs(end - start) / (amount-1);
		vector <double> res(1,start);
		for (int i = 0; i < amount-1; i++)
			res.push_back(res[i] + step);
		return res;
	}

	vector<double> Fragmentation::regular_fragmentation(double start, double end, unsigned int amount, vector<double> &step)
	{
		if (end <= start  || fabs(end - start) / (amount) == 0)
			throw exception("incorrect input regular_fragmentation parameters");
		step.push_back(fabs(end - start) / (amount - 1));
		vector <double> res(1, start);
		for (int i = 0; i < amount - 1; i++)
			res.push_back(res[i] + step[0]);
		step.resize(res.size()-1, fabs(end - start) / (amount - 1));
		return res;
	}

	vector<double> Fragmentation::regular_fragmentation(double start, double end, double step)
	{
		if (end <= start || step<=0)
			throw exception("incorrect input regular_fragmentation parameters");
		vector <double> res(1, start);
		for (int i = 0; res[i]<end; i++)
			res.push_back(res[i] + step);
		res.pop_back();
		res.push_back(end);
		return res;
	}


	vector <double> Fragmentation::adaptive_fragmentation(double start, double end, unsigned int amount, double multiplier)
	{
		if (end <= start || multiplier < 1 )
			throw exception("incorrect input adaptive_fragmentation parameters");
		if (multiplier == 1)
			return Fragmentation::regular_fragmentation(start, end, amount);
			vector <double> h(1, 0),res(1,start);
			h.push_back(fabs(end - start) / pow(multiplier, amount - 2));
			for (int i = 1; i < amount; i++)
			{
				h.push_back(h[i] * multiplier);
				res.push_back(h[i] + start);
			}
	 return res;

	}

	vector<double> Fragmentation::adaptive_fragmentation(double start, double end, vector<double> &step, unsigned int amount, double multiplier)
	{
		if (end <= start || multiplier < 1)
			throw exception("incorrect input adaptive_fragmentation parameters");
		if (multiplier == 1)
		{
			return Fragmentation::regular_fragmentation(start, end, amount,step);
		}
		vector <double> h(1, 0), res(1, start);
		h.push_back(fabs(end - start) / pow(multiplier, amount - 2));
		for (int i = 1; i < amount; i++)
		{
			h.push_back(h[i] * multiplier);
			res.push_back(h[i] + start);
			step.push_back(res[i] - res[i - 1]);
		}
		return res;
	}
