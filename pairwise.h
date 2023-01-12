#include "params.h"
//const float LJepsilon = 3.328, LJsigma = 0.81, LJoffset = 0.6, LJcutoff = 4.0, LJm = -3.02e9, LJc = 1.23e10;//Lennard-Jones parameters
//const float LJepsilon = 0.4, LJsigma = 0.4, LJoffset = 0.3;
//const float LJ4e = 4*LJepsilon, LJ4es6 = LJ4e*pow(LJsigma, 6), LJ4es12 = LJ4e*pow(LJsigma, 12), LJ24es6 = LJ4e*6*pow(LJsigma, 6), LJ48es12 = LJ4e*12*pow(LJsigma, 12); //Pre-calculated LJ parameters for faster multiplication
int maxi, maxj, maxk;

struct gridCell
{
	vector<vector<int> > list;
	vector<short> map; vector<short> revmap;
	inline gridCell(){}
	inline void purge()
	{
		list.clear(); map.clear(); revmap.clear();
	}
	inline void print()
	{
		for (int i = 0; i < list.size(); i++)
		{
			cout << revmap[i] << ": ";
			for (int j = 0; j < list[i].size(); j++) cout << list[i][j] << ' ';
			cout << '\n';
		}
	}
};

struct LJinteraction
{
	rigidBody *rb1, *rb2;
	int b1, b2;
	intParam *ip;
	bool umbrella;
	inline LJinteraction(){}
	inline LJinteraction(rigidBody *_rb1, rigidBody *_rb2, int _b1, int _b2, intParam *_ip, bool _umbrella) : rb1(_rb1), rb2(_rb2), b1(_b1), b2(_b2), ip(_ip), umbrella(_umbrella)
	{
		//cout << "umbrellaInit " << umbrella << '\n';
		//cout << ' ' << rb1 - RBlistMain << ' ' << b1 << ' ' << rb2 - RBlistMain << ' ' << b2 << ' ' << ip->step << '\n';
	}
	inline double LJpairwise(bool calcForce, double &LJumbrellaEnergy) //Pairwise modified Lennard-Jones potential
	{
		//cout << "umbrella " << umbrella << '\n';
		vector3d &v1 = rb1->massList[b1].decomp, &v2 = rb2->massList[b2].decomp, vec = v2 - v1;
		if (wca)
		{
			double d2 = vec * vec;
			//cout << d2 << ' ' << ip->wcaConst << '\n';
			if (d2 >= ip->wcaConst) return 0;
			if (d2 < 0.001) d2 = 0.001;
			d2 = 1/d2; double d6 = d2*d2*d2, d12 = d6*d6;
			double v = ip->LJ4es12*d12 - ip->LJ4es6*d6 + ip->e;
			if (calcForce && !umbrella)
			{
				double f = (ip->LJ24es6*d6 - ip->LJ48es12*d12)*d2;
				vec *= f;
				rb1->force += vec; rb2->force -= vec;
				rb1->massList[b1].moment += vec;
				rb2->massList[b2].moment -= vec;
				return v;
			}
			else if (!umbrella) return v;
			else LJumbrellaEnergy += v;
		}
		const float LJm = -3.02e9, LJc = 1.23e10;//Lennard-Jones parameters
		double mag = vec.mag(), d = mag - ip->step, f, v; vec /= mag;
		if (d > LJcutoff) return 0;
		if (calcForce && !umbrella)
		{
			if (d <= 0)
			{
				f = LJm*(d + ip->s);
				v = LJm*(d + ip->s) + LJc;
			}
			else
			{
				d = 1/d;
				double d2 = d*d, d6 = d2*d2*d2, d12 = d6*d6;
				v = ip->LJ4es12*d12 - ip->LJ4es6*d6;
				f = (ip->LJ24es6*d6 - ip->LJ48es12*d12)*d;
			}
			vec *= f;
			rb1->force += vec; rb2->force -= vec;
			rb1->massList[b1].moment += vec;
			rb2->massList[b2].moment -= vec; 
			return v;
		}
		else if (!umbrella)
		{
			if (d <= 0) return LJm*(d + ip->s) + LJc;
			else
			{
				d = pow(ip->s/d, 6);
				return ip->LJ4e*d*(d - 1);
			}
		}
		else
		{
			if (d <= 0) LJumbrellaEnergy += LJm*(d + ip->s) + LJc;
			else
			{
				d = pow(ip->s/d, 6);
				LJumbrellaEnergy += ip->LJ4e*d*(d - 1);
			}
			return 0;
		}
	}
	/*inline void print()
	{
		cout << "LJinteraction:\n" ;
		cout << rb1 - &RBlist[0] << ' ' << b1 << '\n';
		cout << rb2 - &RBlist[0] << ' ' << b2 << '\n';
	}*/
};

inline double LJpairwise(rigidBody &rb1, rigidBody &rb2, int b1, int b2, bool calcForce) //Pairwise modified Lennard-Jones potential
{
	const float LJm = -3.02e9, LJc = 1.23e10;//Lennard-Jones parameters
	vector3d &v1 = rb1.massList[b1].decomp, &v2 = rb2.massList[b2].decomp, vec = v2 - v1;
	int t = (min(rb1.massList[b1].type, rb2.massList[b2].type) << 8) + max(rb1.massList[b1].type, rb2.massList[b2].type);
	intParam &ip = intTable[t];
	double mag = vec.mag(), d = mag - ip.step, f, v; vec /= mag;
	//if (useGrid && d > LJcutoff) return 0;
	if (calcForce)
	{
		if (d <= 0)
		{
			f = LJm*(d + ip.s);
			v = LJm*(d + ip.s) + LJc;
		}
		else
		{
			d = 1/d;
			double d2 = d*d, d6 = d2*d2*d2, d12 = d6*d6;
			v = ip.LJ4es12*d12 - ip.LJ4es6*d6;
			f = (ip.LJ24es6*d6 - ip.LJ48es12*d12)*d;
		}
		vec *= f;
		rb1.force += vec; rb2.force -= vec;
		rb1.massList[b1].moment += vec;
		rb2.massList[b2].moment -= vec; 
		return v;
	}
	else
	{
		if (d <= 0) return LJm*(d + ip.s) + LJc;
		else
		{
			d = pow(ip.s/d, 6);
			return ip.LJ4e*d*(d - 1);
		}
	}
}

vector<vector<vector<gridCell> > > cellList;
vector<LJinteraction> LJinteractionList;
vector<int> filledCells;
set<int> suppressLJ, umbrellaLJ;
vector<pair<int, int> > suppressed, umbrellaed;
vector<float> LJstrength;

inline int suppressf(int i, int j)
{
	return i + (j<<16);
}

inline void assignGridCell(int i, int j)
{
	int chash = RBlist[i].findCell(j);
	int a = (chash>>24), b = (chash<<12)>>24, c = (chash<<24)>>24;
	if (a > maxi) a = maxi;
	if (b > maxj) b = maxj;
	if (c > maxk) c = maxk;
	if (a < 0) a = 0;
	if (b < 0) b = 0;
	if (c < 0) c = 0;
	chash = (a<<24) + (b<<12) + c;
	gridCell &g = cellList[a][b][c];
	if (g.list.size() == 0)
	{
		//cout << "Regrid: " << a << ' ' << b << ' ' << c << '\n';
		g.map.assign(RBlistSize, -1);
		filledCells.push_back(chash);
	}
	if (g.map[i] == -1)
	{
		g.map[i] = g.list.size(); g.revmap.push_back(i);
		vector<int> v; v.assign(1, j);
		g.list.push_back(v);
	}
	else g.list[g.map[i]].push_back(j);
}

inline void populateList(gridCell &g, gridCell &neighbour)
{
	if (&g != &neighbour)
	{
		for (int i = 0; i < g.list.size(); i++) for (int j = 0; j < neighbour.list.size(); j++)
		if (g.revmap[i] != neighbour.revmap[j] && suppressLJ.count(suppressf(g.revmap[i], neighbour.revmap[j])) == 0) 
		{
			bool umbrella = 0;
			//cout << g.revmap[i] << ' ' << neighbour.revmap[j] << ' ' << suppressf(g.revmap[i], neighbour.revmap[j]) << ' ' << umbrellaLJ.count(suppressf(g.revmap[i], neighbour.revmap[j])) << '\n';
			if (isMD && umbrellaLJ.count(suppressf(g.revmap[i], neighbour.revmap[j]))) umbrella = 1;
			for (int i2 = 0; i2 < g.list[i].size(); i2++) for (int j2 = 0; j2 < neighbour.list[j].size(); j2++)
			{
				char c1 = RBlist[g.revmap[i]].massList[g.list[i][i2]].type, c2 = RBlist[neighbour.revmap[j]].massList[neighbour.list[j][j2]].type; int t;
				vector3d v = RBlist[g.revmap[i]].massList[g.list[i][i2]].decomp - RBlist[neighbour.revmap[j]].massList[neighbour.list[j][j2]].decomp;
				if (c1 != 'Z' && c2 != 'Z' && v.mag() < verletRadius)
				{
					if (c1 > c2) t = (c2<<8) + c1; else t = (c1<<8) + c2;
					//cout << (int)c1 << ' ' << (int)c2 << ' ' << t << '\n';
					LJinteraction LJi(&RBlist[g.revmap[i]], &RBlist[neighbour.revmap[j]], g.list[i][i2], neighbour.list[j][j2], &intTable[t], umbrella);
					//cout << LJinteractionList.size() << ' ' << &RBlist[g.revmap[i]] << '\n';
					LJinteractionList.push_back(LJi);
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < g.list.size(); i++) for (int j = i + 1; j < g.list.size(); j++)
		if (suppressLJ.count(suppressf(g.revmap[i], g.revmap[j])) == 0)
		{
			bool umbrella = 0;
			//cout << g.revmap[i] << ' ' << g.revmap[j] << ' ' << suppressf(g.revmap[i], g.revmap[j]) << ' ' << umbrellaLJ.count(suppressf(g.revmap[i], g.revmap[j])) << '\n';
			if (isMD && umbrellaLJ.count(suppressf(g.revmap[i], g.revmap[j]))) umbrella = 1;
			for (int i2 = 0; i2 < g.list[i].size(); i2++) for (int j2 = 0; j2 < g.list[j].size(); j2++)
			{
				char c1 = RBlist[g.revmap[i]].massList[g.list[i][i2]].type, c2 = RBlist[g.revmap[j]].massList[g.list[j][j2]].type; int t;	
				vector3d v = RBlist[g.revmap[i]].massList[g.list[i][i2]].decomp - RBlist[g.revmap[j]].massList[g.list[j][j2]].decomp;				
				if (c1 != 'Z' && c2 != 'Z' && v.mag() < verletRadius)
				{
					if (c1 > c2) t = (c2<<8) + c1; else t = (c1<<8) + c2;
					//cout << (int)c1 << ' ' << (int)c2 << ' ' << t << '\n';
					//cout << LJinteractionList.size() << ' ' << &RBlist[g.revmap[i]] << '\n';
					LJinteraction LJi(&RBlist[g.revmap[i]], &RBlist[g.revmap[j]], g.list[i][i2], g.list[j][j2], &intTable[t], umbrella);
					LJinteractionList.push_back(LJi);
				}
			}
		}
	}
}

inline void findNeighbours(int chash, gridCell &g)
{
	int i, j;
	int a = (chash>>24), b = (chash<<12)>>24, k = (chash<<24)>>24;
	//cout << chash << ' ' << a << ' ' << b << ' ' << k << '\n';
	//cout << maxi << ' ' << maxj << ' ' << maxk << '\n';
	if (a == 0 && b == 0)
	{
		if (k < maxk)
		{
			for (i = 0; i < 2 && i <= maxi; i++) for (j = 0; j < 2 && j <= maxj; j++)
			{
				gridCell &n = cellList[i][j][k + 1];
				if (n.revmap.size() != 0) populateList(g, n);
			}
		}
		for (i = 0; i < 2 && i <= maxi; i++) for (j = 0; j < 2 && j <= maxj; j++)
		{
			gridCell &n = cellList[i][j][k];
			if (n.revmap.size() != 0) populateList(g, n);
		}
	}
	else if (a == 0)
	{
		if (k < maxk)
		{
			for (i = 0; i < 2 && i <= maxi; i++) for (j = b - 1; j <= b + 1 && j <= maxj; j++)
			{
				gridCell &n = cellList[i][j][k + 1];
				if (n.revmap.size() != 0) populateList(g, n);
			}
		}
		for (i = 0; i < 2 && i <= maxi; i++) for (j = b; j <= b + 1 && j <= maxj; j++)
		{
			gridCell &n = cellList[i][j][k];
			if (n.revmap.size() != 0) populateList(g, n);
		}
	}
	else if (b == 0)
	{
		if (k < maxk)
		{
			for (i = a - 1; i <= a + 1 && i <= maxi; i++) for (j = 0; j < 2 && j <= maxj; j++)
			{
				gridCell &n = cellList[i][j][k + 1];
				if (n.revmap.size() != 0) populateList(g, n);
			}
		}
		for (i = a; i <= a + 1 && i <= maxi; i++) for (j = 0; j < 2 && j <= maxj; j++)
		{
			gridCell &n = cellList[i][j][k];
			if (n.revmap.size() != 0) populateList(g, n);
		}
		gridCell &n = cellList[a - 1][1][k];
		if (n.revmap.size() != 0) populateList(g, n);
	}
	else
	{
		if (k < maxk)
		{
			for (i = a - 1; i <= a + 1 && i <= maxi; i++) for (j = b - 1; j <= b + 1 && j <= maxj; j++)
			{
				//cout << i << ' ' << j << ' ' << k << ' ' << cellList[i][j].size() << '\n';
				gridCell &n = cellList[i][j][k + 1];
				if (n.revmap.size() != 0) populateList(g, n);
			}
		}
		//cout << "step1\n";
		for (i = a; i <= a + 1 && i <= maxi; i++) for (j = b; j <= b + 1 && j <= maxj; j++)
		{
			//cout << i << ' ' << j << ' ' << k << '\n';
			gridCell &n = cellList[i][j][k];
			if (n.revmap.size() != 0) populateList(g, n);
		}
		//cout << "step2\n";
		if (b < maxj)
		{
			gridCell &n = cellList[a - 1][b + 1][k];
			if (n.revmap.size() != 0) populateList(g, n);
		}
		//cout << "step3\n";
	}
}

inline void initDecompAll()
{
	vector3d berth(4, 4, 4);
	LJinteractionList.clear();
	minDim.assign(Big, Big, Big); maxDim.assign(Smol, Smol, Smol); int i, j;
	for (i = 0; i < filledCells.size(); i++)
	{
		cellList[filledCells[i]>>24][(filledCells[i]<<12)>>24][(filledCells[i]<<24)>>24].purge();
	}
	filledCells.clear();
	//cout << RBlistSize << '\n';
	for (i = 0; i < RBlistSize; i++) RBlist[i].minDecompose();
	maxDim += berth; minDim -= berth;
	boxDim = maxDim - minDim; boxDim /= verletRadius;
	maxi = (int)boxDim.i; maxj = (int)boxDim.j; maxk = (int)boxDim.k;
	//cout << maxi << ' ' << maxj << ' ' << maxk << '\n'; return;
	cellList.resize(maxi + 1);
	for (i = 0; i < cellList.size(); i++)
	{
		cellList[i].resize(maxj + 1);
		for (j = 0; j < cellList[i].size(); j++) cellList[i][j].resize(maxk + 1);
	}
	for (i = 0; i < RBlistSize; i++) for (j = 0; j < RBlist[i].massListSize; j++)
	{
		assignGridCell(i, j);
	}
	for (i = 0; i < RBlistSize; i++) for (j = 0; j < RBlist[i].massListSize; j++)
	{
		RBlist[i].massList[j].oldPosition = RBlist[i].massList[j].decomp;
	}
	updateCellList = 1;
}

inline void regrid()
{
	//cout << "Regridding... " << currStep << '\n';
	LJinteractionList.clear();
	int i, j;
	for (i = 0; i < filledCells.size(); i++)
	{
		cellList[filledCells[i]>>24][(filledCells[i]<<12)>>24][(filledCells[i]<<24)>>24].purge();
	}
	filledCells.clear();
	for (i = 0; i < RBlistSize; i++) for (j = 0; j < RBlist[i].massListSize; j++)
	{
		assignGridCell(i, j);
	}
	for (i = 0; i < RBlistSize; i++) for (j = 0; j < RBlist[i].massListSize; j++)
	{
		RBlist[i].massList[j].oldPosition = RBlist[i].massList[j].decomp;
	}
	//cout << "Finished regridding\n";
}

inline void printGrid()
{
	for (int i = 0; i < filledCells.size(); i++) cellList[filledCells[i]>>24][(filledCells[i]<<12)>>24][(filledCells[i]<<24)>>24].print();
}

inline double LJpotential(bool calcForce) //Modified Lennard-Jones potential with cell list
{
	int i, j, k, l; double total = 0; LJumbrellaEnergy = 0;
	bool useGrid = 1;
	if (!useGrid)
	{
		for (i = 0; i < RBlistSize; i++) for (k = i + 1; k < RBlistSize; k++)
		{
			if (suppressLJ.count(suppressf(i, k)) == 0)
			{
				for (j = 0; j < RBlist[i].massListSize; j++) for (l = 0; l < RBlist[k].massListSize; l++) total += LJpairwise(RBlist[i], RBlist[k], j, l, calcForce);
			}
		}
	}
	else
	{
		if (updateCellList)
		{
			int t = clock();
			regrid(); updateCellList = 0;
			for (i = 0; i < filledCells.size(); i++) 
			{
				findNeighbours(filledCells[i], cellList[filledCells[i]>>24][(filledCells[i]<<12)>>24][(filledCells[i]<<24)>>24]);
			}
			updateCellList = 0;
			cellListTime += clock() - t;
		}
		for (i = 0; i < LJinteractionList.size(); i++)
		{
			total += LJinteractionList[i].LJpairwise(calcForce, LJumbrellaEnergy);
		}
	}
	LJcontrib = total;
	return total;
}

inline double electrostatic(bool calcForce) //since most particles are uncharged, it makes more sense to keep a list of the few charged ones rather than do pairwise interactions
{
	double total = 0.0, d, d2, v;
	vector3d vec;
	for (int i = 0; i < RBlistSize; i++) for (int j = 0; j < RBlist[i].chargeListSize; j++) for (int k = i + 1; k < RBlistSize; k++) for (int l = 0; l < RBlist[k].chargeListSize; l++)
	{
		vec = RBlist[k].chargeList[l].decomp - RBlist[i].chargeList[j].decomp;
		d2 = vec.square();
		if (d2 != 0)
		{
			d = sqrt(d2); v = CAdj/d*RBlist[i].chargeList[j].charge*RBlist[k].chargeList[l].charge; total += v;
			if (calcForce)
			{
				double u = -v/d2; vec *= u;
				RBlist[i].force += vec; RBlist[k].force -= vec;
				RBlist[i].chargeList[j].moment += vec; vec *= -1;
				RBlist[k].chargeList[l].moment += vec;
			}
		}
	}
	return total;
}
