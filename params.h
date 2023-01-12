struct intParam
{
	float s, e, step;
	float LJ4e, LJ4es6, LJ4es12, LJ24es6, LJ48es12, wcaConst;
	inline intParam(){}
	inline void preCalc()
	{
		LJ4e = 4*e; LJ4es6 = LJ4e*pow(s, 6); LJ4es12 = LJ4e*pow(s, 12); LJ24es6 = LJ4e*6*pow(s, 6); LJ48es12 = LJ4e*12*pow(s, 12);
	}
	inline intParam(double _s, double _e, double _step)
	{
		if (wca == 1)
		{
			s = _s + _step; e = _e; step = 0; wcaConst = s*s*pow(2, 0.3333333);
		}
		else
		{
			s = _s; e = _e; step = _step;
		}
		preCalc();
	}
};

struct LJparam
{
	float mass, s, e, step;
	string displayAs;
	inline LJparam(){}
	inline LJparam(double _mass, double _s, double _e, double _step){mass = _mass * massMult, s = _s; e = _e; step = _step; displayAs = "H";}
	inline LJparam(double _mass, double _s, double _e, double _step, string &_displayAs){mass = _mass * massMult, s = _s; e = _e; step = _step; displayAs = _displayAs;}
};

map<char, LJparam> paramTable;

inline void assignLJparam()
{
	string C = "C", O = "O", U = "U";
	LJparam LJA(1, 0.4, 3.328, 0.3, C);
	LJparam LJB(0.12, 0.39, 0.2, 0);
	LJparam LJS(0.34, 0.28, 0.8, 0.21, O);
	LJparam LJY(0, 0, 0, 0, U);
	LJparam LJZ(1000, 0, 0, 0);
	LJparam LJH(10000, 0, 0, 0);
	paramTable['A'] = LJA;
	paramTable['B'] = LJB;
	paramTable['C'] = LJB;
	paramTable['D'] = LJA;
	paramTable['F'] = LJS;
	paramTable['S'] = LJS;
	paramTable['T'] = LJS;
	paramTable['Z'] = LJZ;
	paramTable['Y'] = LJY;
	paramTable['H'] = LJH;
}

map<int, intParam> intTable;

inline void populateIntTable()
{
	set<char> us;
	vector<char> v;
	for (int i = 0; i < RBlistSize; i++)
	{
		for (int j = 0; j < RBlist[i].massListSize; j++)
		{
			char c = RBlist[i].massList[j].type;
			if (us.count(c) == 0)
			{
				us.insert(c);
				v.push_back(c);
			}
		}
	}
//	cout << v.size() << '\n';
	for (int i = 0; i < v.size(); i++) for (int j = i; j < v.size(); j++)
	{
		//cout << v[i] << ' ' << v[j] << '\n';
		int t = (min(v[i], v[j]) << 8) + max(v[i], v[j]);
		LJparam &LJ1 = paramTable[v[i]], &LJ2 = paramTable[v[j]];
		if (v[i] == 'Z' || v[j] == 'Z')
		{
		}
		if (v[i] == 'F' && v[j] == 'F')
		{
			intParam ip(0, 0, 0);
			intTable[t] = ip;
		}
		else if (v[i] == 'D' && v[j] == 'C')
		{
			intParam ip(0, 0, 0);
			intTable[t] = ip;
		}
		else if (v[i] == 'E' && v[j] == 'C')
		{
			intParam ip(0, 0, 0);
			intTable[t] = ip;
		}
		else if (v[i] == 'T' && v[j] == 'C')
		{
			intParam ip(0, 0, 0);
			intTable[t] = ip;
		}
		else
		{
			intParam ip(LJ1.s + LJ2.s, sqrt(LJ1.e * LJ2.e), LJ1.step + LJ2.step);
			intTable[t] = ip;
		}
	}
}



