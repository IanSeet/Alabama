struct fitParam
{
	double weight = 1;
	unsigned char type;
	virtual double fitFunc(vector<vector<map<float, int> > > & pol, ofstream &ofs){}
	virtual double fitFunc(vector<energyConfig> &vcg, double threshold, ofstream &ofs){}
	virtual double fitFunc(vector<double> &expenditureList, double threshold, ofstream &ofs){}
	virtual void parse(istringstream &iss){}
};

struct fitTraj : fitParam
{
	double start, end, depth;
	bool isGreater = 0;
	int target;
	fitTraj(){type = 0;}
	double fitFunc(vector<vector<map<float, int> > > & pol, ofstream &ofs)
	{
		ofs << "fitTraj\n";
		ofs << start << ' ' << end << ' ' << depth << '\n';
		int a = start * pol.size(), b = end * pol.size(); if (a < 0) a = 0; if (b > pol.size()) b = pol.size();
		double loss = 0.0;
		double f, g, h = 0;
		for (int i = a; i < b; i++)
		{
			ofs << "Step no: " << i << '\n';
			map<float, int>::iterator it;
			for (it = pol[i][target].begin(); it != pol[i][target].end(); ++it)
			{
				f = it->first - depth;
				if (isGreater)
				{
					if (f > 0) f = 0;
				}
				else
				{
					if (f < 0) f = 0;
				}
				g = f*f*it->second;
				h += it->second;
				loss += g;
				ofs << it->first << ' ' << it->second << ' ' << g << '\n';
			}
		}
		if (std::isnan(loss)) return intmax;
		return loss/h;
	}
	void parse(istringstream &iss)
	{
		iss >> weight >> start >> end >> depth >> isGreater >> target;
	}
};

struct fitRestBound : fitParam
{
	double start, end, equil, bound;
	int target;
	fitRestBound(){type = 1;}
	double fitFunc(vector<vector<map<float, int> > > & pol, ofstream &ofs)
	{
		ofs << "fitRestBound\n";
		ofs << start << ' ' << end << ' ' << equil << ' ' << bound << '\n';
		int a = start * pol.size(), b = end * pol.size(); if (a < 0) a = 0; if (b > pol.size()) b = pol.size();
		double loss = 0.0;
		double f, g, h = 0;
		for (int i = a; i < b; i++)
		{
			ofs << "Step no: " << i << '\n';
			map<float, int>::iterator it;
			for (it = pol[i][target].begin(); it != pol[i][target].end(); ++it)
			{
				f = it->first - equil; if (fabs(f) > bound) f = fabs(f) - bound; else f = 0;
				g = f*f*it->second;
				loss += g;
				h += it->second;
				ofs << it->first << ' ' << it->second << ' ' << g << '\n';
			}
		}
		if (std::isnan(loss)) return intmax;
		return loss/h;
	}
	void parse(istringstream &iss)
	{
		iss >> weight >> start >> end >> equil >> bound >> target;
	}
};

struct fitRestMean : fitParam
{
	double start, end, equil;
	int target;
	fitRestMean(){type = 2;}
	double fitFunc(vector<vector<map<float, int> > > & pol, ofstream &ofs)
	{
		ofs << "fitRestMean\n";
		ofs << start << ' ' << end << ' ' << equil << '\n';
		int a = start * pol.size(), b = end * pol.size(); if (a < 0) a = 0; if (b > pol.size()) b = pol.size();
		double loss = 0.0;
		double f = 0, g = 0;
		for (int i = a; i < b; i++)
		{
			ofs << "Step no: " << i;
			map<float, int>::iterator it;
			for (it = pol[i][target].begin(); it != pol[i][target].end(); ++it)
			{
				f += it->first*it->second;
				g += it->second;
				ofs << it->first << ' ' << it->second << ' ' << g << '\n';
			}
		}
		f /= g;
		loss = (f - equil)*(f - equil);
		if (std::isnan(loss)) return intmax;
		return loss;
	}
	void parse(istringstream &iss)
	{
		iss >> weight >> start >> end >> equil >> target;
	}
};

struct fitPot : fitParam
{
	double timeStart, timeEnd;
	fitPot(){type = 10;}
	double fitFunc(vector<energyConfig> &vcg, double threshold, ofstream &ofs)
	{
		ofs << "fitPot\n";
		int a = timeStart * vcg.size(), b = timeEnd * vcg.size(); if (a < 0) a = 0; if (b >= vcg.size()) b = vcg.size() - 1;
		double f = vcg[a].V - vcg[b].V; f /= vcg.size();
		ofs << vcg[a].V << ' ' << vcg[b].V << '\n';
		if (std::isnan(f)) return intmax;
		return f*f;
	}
	void parse(istringstream &iss)
	{
		iss >> weight >> timeStart >> timeEnd;
	}
};

struct fitPotVar : fitParam
{
	fitPotVar(){type = 12;}
	double fitFunc(vector<energyConfig> &vcg, double threshold, ofstream &ofs)
	{
		ofs << "fitPot\n";
		double meanV = 0, varV = 0;
		if (vcg.size() == 0) return 0;
		for (int i = 0; i < vcg.size(); i++)
		{
			meanV += vcg[i].V;
		}
		meanV /= vcg.size();
		for (int i = 0; i < vcg.size(); i++)
		{
			double temp = (vcg[i].V - meanV)*(vcg[i].V - meanV);
			if (threshold >= temp) temp = 0; else temp -= threshold;
			varV += temp;
		}
		varV /= vcg.size();
		return varV;
	}
	void parse(istringstream &iss)
	{
		iss >> weight;
	}
};

struct fitKin : fitParam
{
	fitKin(){type = 11;}
	double fitFunc(vector<energyConfig> &vcg, double threshold, ofstream &ofs)
	{
		ofs << "fitKin\n";
		double meanKE = 0;
		for (int i = 0; i < vcg.size(); i++)
		{
			meanKE += vcg[i].TKE + vcg[i].RKE;
			ofs << vcg[i].TKE + vcg[i].RKE << '\n';
		}
		meanKE /= vcg.size();
		double f = meanKE - threshold;
		if (std::isnan(f)) return intmax;
		return f*f;
	}
	void parse(istringstream &iss)
	{
		iss >> weight;
	}
};

struct fitExp : fitParam
{
	double timeFrac;
	fitExp(){type = 20;}
	double fitFunc(vector<double> &expenditureList, double threshold, ofstream &ofs)
	{
		ofs << "fitExp\n";
		double f = 0;
		for (int i = 0; i < expenditureList.size(); i++)
		{
			if (threshold >= 0)
			{
				double g = expenditureList[i] - threshold;
				if (g < 0) g = 0;
				f += g;
			}
			else
			{
				f += expenditureList[i];
			}
			ofs << expenditureList[i] << '\n';
		}
		if (std::isnan(f)) return intmax;
		f /= expenditureList.size();
		return f;
	}
	void parse(istringstream &iss)
	{
		iss >> weight >> timeFrac;
	}
};

vector <fitParam*> fitParamList;
double fitness = 0;

double fitFunction(vector<vector<map<float, int> > > & pol, vector<energyConfig> &vcg, vector<double> &expList, double threshold, ofstream &ofs)
{
	ofs << fitParamList.size() << '\n';
	double ffo, fft = 0;
	vector<double> vd(fitParamList.size());
	for (int i = 0; i < fitParamList.size(); i++)
	{
		ofs << "Param No: " << i << '\n';
		fitParam * fp = fitParamList[i];
		if (fp->type < 10)
		{
			ffo = fp->fitFunc(pol, ofs);
		}
		else if (fp->type < 20)
		{
			ffo = fp->fitFunc(vcg, threshold, ofs);
		}
		else
		{
			ffo = fp->fitFunc(expList, threshold, ofs);
		}
		vd[i] = ffo * fp->weight;
		fft += vd[i];
	}
	for (int i = 0; i < fitParamList.size(); i++)
	{
		ofs << "fitParam " << i << " output: " << vd[i] << '\n';
	}
	ofs << "total output: " << fft << '\n';
	fitness = fft;
	return fft;
}

void parseFit(ifstream &ifs)
{
	string s; int i;
	while (getline(ifs, s))
	{
		istringstream iss0(s.c_str());
		iss0 >> i;
		cout << s << '\n' << i << '\n';
		switch(i)
		{
			fitParam *f;
			case 0:
				f = new fitTraj();
				f->parse(iss0);
				fitParamList.push_back(f);
				break;
			case 1:
				f = new fitRestBound();
				f->parse(iss0);
				fitParamList.push_back(f);
				break;
			case 2:
				f = new fitRestMean();
				f->parse(iss0);
				fitParamList.push_back(f);
				break;
			case 10:
				f = new fitPot();
				f->parse(iss0);
				fitParamList.push_back(f);
				break;
			case 11:
				f = new fitKin();
				f->parse(iss0);
				fitParamList.push_back(f);
				break;
			case 12:
				f = new fitPotVar();
				f->parse(iss0);
				fitParamList.push_back(f);
				break;
			case 20:
				f = new fitExp();
				f->parse(iss0);
				fitParamList.push_back(f);
				break;
		}
		iss0.clear();
	}
}
