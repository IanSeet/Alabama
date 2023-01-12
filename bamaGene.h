void mutate(vector<gene>& geneList, float rate)
{
	int a = rate*1000, k = 1;
	double intDisp = 0.05, rigidDisp = 0.02;
	for (int i = 0; i < geneList.size(); i++)
	{
		int c = rand()%1000;
		if (a >= c)
		{
			//cout << "mutate: " << i << ' ' << a << ' ' << c << '\n';
			if (geneList[i].type == 0)
			{
				int pos = rand()%2;
				if (!pos) pos = -1;
				geneList[i].shift *= 1.0 + pos*intDisp;
			}
			else
			{
				int pos1 = rand()%3 - 1, pos2 = rand()%3 - 1, pos3 = rand()%3 - 1;
				vector3d v(pos1*rigidDisp, pos2*rigidDisp, pos3*rigidDisp);
				geneList[i].disp += v;
				//geneList[i].disp.printtocout();
			}
		}
	}
}

void printGenes(vector<gene>& geneList, ofstream &ofs)
{
	ofs << geneList.size() << '\n';
	for (int i = 0; i < geneList.size(); i++)
	{
		geneList[i].print(ofs);
	}
}

void parseGenes(vector<gene>& geneList, ifstream &ifs)
{
	int a;
	ifs >> a;
	cout << "a: " << a << '\n';
	if (a > 1000000) return;
	geneList.clear(); geneList.resize(a);
	for (int i = 0; i < geneList.size(); i++)
	{
		geneList[i].parse(ifs);
	}
}

struct parent
{
	vector<gene> geneList;
	double fitness;
	int idx;

	bool operator < (const parent &otr) const 
	{
	        return fitness < otr.fitness;
	}
	
	void print(ofstream &ofs)
	{
		ofs << "index: " << idx << '\n';
		ofs << "fitness: " << fitness << '\n';
		ofs << "geneList: " << geneList.size() << '\n';
		for (int i = 0; i < geneList.size(); i++) geneList[i].printVerbose(ofs);
	}

	void printSparse(ofstream &ofs)
	{
		ofs << "index: " << idx << '\n';
		ofs << "fitness: " << fitness << '\n';
		//ofs << "geneList: " << geneList.size() << '\n';
		//for (int i = 0; i < geneList.size(); i++) geneList[i].printVerbose(ofs);
	}
};

void crossover(parent &s1, parent &s2, vector<vector<gene> > &offspringList, int target)
{
	offspringList[target].resize(s1.geneList.size());
	for (int j = 0; j < s1.geneList.size(); j++)
	{
		char c = rand()%2;
		if (c) offspringList[target][j] = s1.geneList[j];
		else offspringList[target][j] = s2.geneList[j];
	}	
}

void printParentList(vector<vector<parent> > & parentList, int gen, ofstream &ofs)
{
	//cout << "size: " << parentList[gen].size() << ' ' << gen << '\n';
	ofs << "size: " << parentList[gen].size() << '\n';
	for (int i = 0; i < parentList[gen].size(); i++)
	{
		ofs << "parentNo: " << i << '\n';
		parentList[gen][i].printSparse(ofs);
	}
	ofs.close();
}

void printFitnessMatrix(vector<parent> &temp, ofstream &ofs)
{
	for (int i = 0; i < temp.size(); i++)
	{
		ofs << temp[i].fitness << ' ';
		for (int j = 0; j < temp[i].geneList.size(); j++)
		{
			if (temp[i].geneList[j].type == 0)
			{
				ofs << temp[i].geneList[j].shift << ' ';
			}
			else
			{
				temp[i].geneList[j].disp.printNoNewline(ofs);
			}
		}
		ofs << '\n';
	}
}

void parseGenAux(string &genAux, genAlgoParam &gp)
{
	ifstream genaux(genAux.c_str());
	genaux >> gp.generations >> gp.popsize >> gp.mutrate >> gp.maxthreads;
	genaux >> gp.steps >> gp.runs >> gp.eqsteps >> gp.roteqsteps;
	genaux >> gp.isConden;
	//cout << "gp.isConden: " << gp.isConden << '\n';
}

int probeGenActList(string s)
{
	ifstream ifs(s.c_str()); int i = 0;
	while (getline(ifs, s)) i++;
	return i;
}

void genAlgo(int generations, int popSize, int procid, int numprocs)
{
	genAlgoParam gp(1, 1, 0.5, 1, 100000, 1, 50000, 0);
	string genauxName("genAux");
	parseGenAux(genauxName, gp);
	string defFolder = "genr", genName;
	vector<vector<parent> > parentList;
	vector<vector<vector<gene> > > offspringList;
	vector<gene> seedList;
	int limit, elite;
	bool isElite;
	vector<int> fittest;
	int genActListSize = probeGenActList("genAct");
	if (procid == 0)
	{
		string systemCall = "./RigorMortis inputGen seed";
		int status = system(systemCall.c_str());
		ifstream seed("seed");
		parentList.resize(generations);
		offspringList.resize(generations + 1);
		fittest.resize(generations);
		for (int i = 0; i < generations + 1; i++)
		{
			offspringList[i].resize(popSize);
		}
		parseGenes(seedList, seed);
		for (int i = 0; i < popSize; i++)
		{
			offspringList[0][i] = seedList;
		}
		
		const double threshold = 0.25, elitism = 0.05;
		limit = threshold * popSize; if (limit == 0) limit = 1;
		elite = elitism * popSize; if (elite > limit) elite = limit;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (int j = 0; j < generations; j++)
	{
		genName = defFolder + to_string(j);
		if (procid == 0)
		{
			mkdir(genName.c_str(), 0700);

			chdir(genName.c_str());
			for (int i = 0; i < popSize; i++)
			{
				string gfilename = genName + '_' + to_string(i);
				ofstream ofs(gfilename.c_str());
				if (j == 0) mutate(offspringList[0][i], gp.mutrate);
				printGenes(offspringList[j][i], ofs);
			}
			chdir("..");
		}
	
		MPI_Barrier(MPI_COMM_WORLD);
		
		for (int i = 0; i < popSize; i++)
		{
			if (i % numprocs != procid) continue;
			string gfilename = genName + '_' + to_string(i);
			if (i < elite) isElite = 1; else isElite = 0;
			string systemCall = "./RigorMortis inputGen gen " + gfilename + ' ' + genName + ' ' + to_string(i) + ' ' + to_string(isElite);
			int status = system(systemCall.c_str());
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (procid == 0)
		{
			vector<parent> temp(popSize);
			
			chdir(genName.c_str());
			for (int i = 0; i < popSize; i++)
			{
				string gfilename = genName + '_' + to_string(i);
				string ofilename =  gfilename;
				cout << "ofilename: " << ofilename << '\n';
				ifstream ifs(ofilename.c_str());
				parseGenes(temp[i].geneList, ifs);
				temp[i].idx = i;
				ifstream ifsf((gfilename + "fit").c_str());
				ifsf >> temp[i].fitness;
			}
			
			sort(temp.begin(), temp.end());
			cout << "parents sorted\n";
			string fitMatrixName = genName + "fitMatrix";
			ofstream ofsfm(fitMatrixName.c_str());
			printFitnessMatrix(temp, ofsfm);
			parentList[j].resize(limit);
			
			for (int i = 0; i < limit; i++) parentList[j][i] = temp[i];
			fittest[j] = parentList[j][0].idx;
			
			string parentListName = genName + "parents";
			ofstream ofsp(parentListName.c_str());
			//cout << "parentListSize: " << parentList[j].size() << '\n';
			printParentList(parentList, j, ofsp);
			string fitListName = genName + "fitList";
			ofstream ofsfl(fitListName.c_str());
			double meanFitness = 0;
			for (int i = 0; i < limit; i++)
			{
				meanFitness += parentList[j][i].fitness;
			}
			meanFitness /= limit;
			ofsfl << meanFitness << '\n';
			for (int i = 0; i < limit; i++)
			{
				ofsfl << parentList[j][i].fitness << '\n';
			}
			chdir("..");
			
			vector<vector<gene> > &g = offspringList[j + 1];
			for (int i = 0; i < elite; i++)
			{
				g[i] = parentList[j][i].geneList;
			}
			for (int i = elite; i < popSize; i++)
			{
				int parent1 = rand()%limit, parent2 = rand()%limit;
				crossover(parentList[j][parent1], parentList[j][parent2], g, i);
			}
			for (int i = elite; i < popSize; i++)
			{
				mutate(g[i], gp.mutrate);
			}
			cout << "crossover done\n";
			string systemCall = "python purge.py " + to_string(j) + ' ' + to_string(popSize) + ' ' + to_string(genActListSize);
			int status = system(systemCall.c_str());
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
	}

	
	if (procid == 0)
	{
		ofstream ofs("fittest");
		for (int i = 0; i < fittest.size(); i++)
		{
			ofs << fittest[i] << '\n';
		}
	}
}
