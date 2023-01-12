struct simSet
{
	vector<simulation> simList;
	double fitness = 0;
	string name;
};

struct parent
{
	vector<gene> geneList;
	double fitness;
	int idx;
	simSet * sim;

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

struct lossParam
{
	double start, end, depth, equil;
	int target;
};

struct genParam
{
	vector<simSet*> vecsim;
	vector<lossParam*> vecloss;
	int steps, runs, roteqsteps, eqsteps, randNo;
	string folderName;
	bool writeGene = 1;
	genParam(int _steps, int _runs) : steps(_steps), runs(_runs){eqsteps = 50000, roteqsteps = 0;}
	genParam(int _steps, int _runs, int _eqsteps, int _roteqsteps):
	steps(_steps), runs(_runs), eqsteps(_eqsteps), roteqsteps(_roteqsteps){}
};

static void *genF(void *p)
{
	genParam *sp = (genParam*)p;
	int steps = sp->steps, runs = sp->runs;
	const int minsteps = 20000, eqsteps = sp->eqsteps, areq = sp->eqsteps, roteqsteps = sp->roteqsteps, write = steps/100;
	bool roteq = 0; if (sp->roteqsteps > 0) roteq = 1;
	string s, fout, ftype = "top";
	int seed = time(NULL) + sp->randNo;
	srand(seed);
	for (int j = 0; j < sp->vecsim.size(); j++)
	{
		for (int k = 0; k < sp->vecsim[j]->simList.size(); k++)
		{
			simulation &sim = sp->vecsim[j]->simList[k];
			cout << "starting simulation " << sim.name << ' ' << &sim << '\n';
			fout = sim.name;
			sim.initialise();
			sim.findObservables();
			//cout << "simulation initialised " << sim.name << '\n';
		//	s = fout + "Gene";
		//	ofstream ofsg(s.c_str());
		//	sim.printGenes(ofsg);
		//	cout << "genes printed " << sim.name << '\n';
			cout << sim.SP() << '\n';
			sim.printContributions(0);
			sim.SDmin(minsteps);
			cout << "min complete\n";
			sim.equilibriate(eqsteps, eqsteps/10);
			if (roteq)
			{
				sim.roteq(roteqsteps, roteqsteps/10);
			}			
			s = fout + ".mol2";
			ofstream ofsm(s.c_str());
			sim.printTopolMol2(ofsm, 0, 1);
			ofsm.close();
			sim.optConfig.init(sim.RBlist, sim.RBlistSize);
			sim.optConfig.assign();
			sim.combinedEnergies.resize(steps/write + 1);
			//parseInputCommon(ftype, fout, sim);
			sim.maxRuns = runs;
			
			s = fout + ".top";
			string s2 = fout + ".conf";
			ofstream ofst; ofst.open(s.c_str());
			sim.printTop(ofst);
			ofst.close();
			ofst.open(s2.c_str());
			sim.printConfig(ofst);
			ofst.close();
			
			s = fout + "Lossf";
			ofstream ofsl(s.c_str());
			for (int i = 0; i < runs; i++)
			{
				sim.runNo = i + 1;
				sim.optConfig.overwrite(); sim.catpInit();
				sim.areq(eqsteps, eqsteps/10);
				sim.langevinIntegrator(steps, write, 0);
				for (int j = 0; j < sim.MDtrajectory.size(); j++) sim.combinedEnergies[j].add(sim.MDtrajectory[j]);
				sim.expenditureList.push_back(sim.expenditure);
				//ofsl << "runNo: " << i << '\n'; 
				sim.findObservables();
			}
			//cout << "lossF: " << sim.lossFunction(0.25, 0.75, -30, 90.0, 0.1, 0, sim.expenditureList.size() - 1, ofsl) << '\n';
			sim.fitFunction(sim.printObservableList, sim.combinedEnergies, sim.expenditureList, 0, ofsl);
			s = fout + ".xyz";
			ofstream ofs(s.c_str());
			for (int i = 0; i < sim.MDtrajectory.size(); i++)
			{
				sim.printConfigxyz(i, sim.MDtrajectory, ofs); 
			}
			s = fout + "Energy";
			ofstream ofs2(s.c_str());
			for (int i = 0; i < sim.combinedEnergies.size(); i++)
			{
				sim.combinedEnergies[i].divide(runs);
				ofs2 << i << ' ' << sim.combinedEnergies[i].TKE << ' ' << sim.combinedEnergies[i].RKE << ' ' << sim.combinedEnergies[i].V << ' ' <<
				sim.combinedEnergies[i].TKE + sim.combinedEnergies[i].RKE + sim.combinedEnergies[i].V << ' ' << sim.combinedEnergies[i].expenditure << '\n';
			}
			s = fout + "Observables";
			ofstream ofs3(s.c_str());
			sim.printObservables(ofs3, writePOL);
			s = fout + "Expenditures";
			ofstream ofs4(s.c_str());
			sim.printExpenditures(ofs4);
			sp->vecsim[j]->fitness += sim.fitness;
		}
		if (sp->writeGene)
		{
			fout  = sp->vecsim[j]->name;
			/*s = fout + "Out";
			//ofstream ofsg(s.c_str());
			//sp->vecsim[j]->simList[0].printGenes(ofsg);
			cout << "genes printed " << '\n';*/
			s = fout + "fit";
			ofstream ofsf(s.c_str()); ofsf << sp->vecsim[j]->fitness;
		}
	}
	return NULL;
}

void printParentList(vector<vector<parent> > & parentList, int gen, ofstream &ofs)
{
	ofs << "size: " << parentList[gen].size() << '\n';
	for (int i = 0; i < parentList[gen].size(); i++)
	{
		ofs << "parentNo: " << i << '\n';
		parentList[gen][i].print(ofs);
	}
}

void parseGenActParams(ifstream &ifs, vector<vector<int> > &vgp, vector<string> &vgal, vector<string> &vs)
{
	string s;
	while (getline(ifs, s))
	{
		istringstream iss(s); int a; cout << "iss: " << s << '\n';
		iss >> s; vs.push_back(s);
		iss >> s; vgal.push_back(s);
		vector<int> vi;
		while (iss >> a)
		{
			vi.push_back(a);
		}
		vgp.push_back(vi);
	}
}

void genAlgo(simulation &seed, genAlgoParam &gap)
{
	bool isConden = gap.isConden;
	const int generations = gap.generations, popSize = gap.popsize;
	vector<vector<parent> > parentList(generations);
	vector<vector<vector<gene> > > offspringList(generations);
	vector<simSet*> fittest(generations);
	for (int i = 0; i < generations; i++)
	{
		offspringList[i].resize(popSize);
	}
	vector<vector<simSet> > vs(generations);

	vector<vector<int> > genActParams;
	vector<string> genFitList, genAuxList;
	ifstream ifsgap("genAct");
	parseGenActParams(ifsgap, genActParams, genAuxList, genFitList);

	for (int i = 0; i < generations; i++) vs[i].resize(popSize);
	
	const int maxthreads = gap.maxthreads;
	vector<pthread_t> genThreads(maxthreads);
	vector<genParam> vp;
	const double threshold = 0.25;
	int limit = threshold * popSize; if (limit == 0) limit = 1;
	string defFolder = "gen";
	
	seed.active.clear();
	cout << "parseFit complete\n";
	
	for (int j = 0; j < generations; j++)
	{
		string genName = defFolder + to_string(j);
		mkdir(genName.c_str(), 0700);
		vp.clear();
		for (int k = 0; k < maxthreads; k++)
		{
			genParam gp(gap.steps, gap.runs, gap.eqsteps, gap.roteqsteps);
			gp.folderName = genName;
			gp.randNo = rand()%10000;
			for (int i = k; i < popSize; i += maxthreads)
			{
				cout << j << ' ' << i << " new offspring\n";
				vs[j][i].simList.resize(genActParams.size());
				for (int l = 0; l < vs[j][i].simList.size(); l++)
				{
					seed.geneCopy(vs[j][i].simList[l]);
					ifstream ifsGen(genFitList[l].c_str());
					vs[j][i].simList[l].parseFit(ifsGen);
					string activeAux = genAuxList[l];
					parseRMaux(activeAux);
					vs[j][i].simList[l].overrideParameters(1);
					for (int m = 0; m < genActParams[l].size(); m++)
					{
						//cout << "gap: " << genActParams[l].size() << ' ' << m << ' ' << genActParams[l][m] << '\n';
						vs[j][i].simList[l].active.insert(genActParams[l][m]);
					}
				}
				if (j > 0)
				{
					for (int l = 0; l < vs[j][i].simList.size(); l++)
					{
						vs[j][i].simList[l].geneList = offspringList[j - 1][i];
					}
				}
				else
				{
					for (int l = 0; l < vs[j][i].simList.size(); l++)
					{
						vs[j][i].simList[l].makeGenes(isConden);
					}
				}
				for (int l = 0; l < vs[j][i].simList.size(); l++)
				{
					simulation &sim = vs[j][i].simList[l];
					sim.name = seed.name + "_c" + to_string(i) + "_a" + to_string(l);
					//cout << vs[i].name << '\n';
					sim.mutate(gap.mutrate);
					sim.phenotype(isConden);
					sim.initDecompAll();
					cout << "mutation complete\n";
				}
				gp.vecsim.push_back(&vs[j][i]);
			}
			vp.push_back(gp);
			//cout << k << ' ' << vp[k].vecsim[0] << '\n';
		}
		cout << "initalising threads\n";
		chdir(genName.c_str());
		for (int i = 0; i < maxthreads; i++)
		{
			int rc = pthread_create(&genThreads[i], NULL, genF, &vp[i]);
			//cout << &vs[i] << '\n';
	  		if (rc)
			{
				cout << "Error: unable to create thread " << rc << '\n';
				exit(-1);
	  		}
		}
		for (int i = 0; i < maxthreads; i++) 
		{
			pthread_join(genThreads[i], NULL);
		}
		cout << "threads joined\n";
		
		vector<parent> temp(popSize);
		for (int i = 0; i < popSize; i++)
		{
			temp[i].geneList = vs[j][i].simList[0].geneList;
			temp[i].fitness = vs[j][i].fitness;
			temp[i].idx = i;
			temp[i].sim = &(vs[j][i]);
		}
		sort(temp.begin(), temp.end());
		cout << "parents sorted\n";
		parentList[j].resize(limit);
		
		for (int i = 0; i < limit; i++) parentList[j][i] = temp[i];
		fittest[j] = parentList[j][0].sim;
		
		string parentListName = genName + "parents";
		ofstream ofsp(parentListName.c_str());
		printParentList(parentList, j, ofsp);
		chdir("..");
		
		vector<vector<gene> > &g = offspringList[j];
		for (int i = 0; i < popSize; i++)
		{
			int parent1 = rand()%limit, parent2 = rand()%limit;
			crossover(parentList[j][parent1], parentList[j][parent2], g, i);
		}
		cout << "crossover done\n";
	}
	ofstream ofs("fittest");
	for (int i = 0; i < fittest.size(); i++)
	{
		ofs << fittest[i]->fitness << '\n';
	}
}

void pyGenAlgo(simulation &seed, genAlgoParam &gap, string &geneListName, string &folderName, bool isElite, bool isTop)
{
	bool isConden = gap.isConden;
	cout << geneListName << ' ' << folderName << ' ' << isElite << '\n';
	const int popSize = gap.popsize;

	vector<vector<int> > genActParams;
	vector<string> genFitList, genAuxList;
	ifstream ifsgap("genAct");
	parseGenActParams(ifsgap, genActParams, genAuxList, genFitList);

	seed.active.clear();
	simSet ss;
	ss.name = geneListName;
	ss.simList.resize(genActParams.size());
	genParam gp(gap.steps, gap.runs, gap.eqsteps, gap.roteqsteps);
	gp.folderName = folderName;
	gp.randNo = rand()%10000;
	cout << "rand set\n";
	//cout << "simListSize " << ss.simList.size() << '\n';
	for (int i = 0; i < ss.simList.size(); i++)
	{
		simulation &sim = ss.simList[i];
		cout << "isTop: " << isTop << '\n';
		if (isTop) seed.deepCopy(sim); else seed.geneCopy(sim);
		ifstream ifsGen(genFitList[i].c_str());
		sim.parseFit(ifsGen);
		string activeAux = genAuxList[i];
		parseRMaux(activeAux);
		ss.simList[i].overrideParameters(1);
		for (int m = 0; m < genActParams[i].size(); m++)
		{
			cout << "gap: " << genActParams[i].size() << ' ' << m << ' ' << genActParams[i][m] << '\n';
			ss.simList[i].active.insert(genActParams[i][m]);
		}
		cout << "active inserted\n";
	}
	chdir(folderName.c_str());
	for (int i = 0; i < ss.simList.size(); i++)
	{
		simulation &sim = ss.simList[i];
		ifstream ifs(geneListName.c_str());
		sim.makeGenes(isConden);
		sim.parseGenes(ifs);
		cout << "genes read\n";
		sim.name = geneListName + "_a" + to_string(i);
		sim.phenotype(isConden);
		cout << "phenotyping done\n";
		sim.initDecompAll();
	}
	gp.vecsim.push_back(&ss);
	genF((void*)&gp);
	chdir("..");
	cout << geneListName << " done\n";
}

void readGenAlgo(simulation &seed, genAlgoParam &gap, string &geneListName)
{
	bool isConden = gap.isConden;
	cout << geneListName << '\n';
	const int popSize = gap.popsize;

	vector<vector<int> > genActParams;
	vector<string> genFitList, genAuxList;
	ifstream ifsgap("genAct");
	parseGenActParams(ifsgap, genActParams, genAuxList, genFitList);

	seed.active.clear();
	simSet ss;
	ss.name = geneListName;
	ss.simList.resize(genActParams.size());
	genParam gp(gap.steps, gap.runs, gap.eqsteps, gap.roteqsteps);
	gp.randNo = rand()%10000;
	cout << "rand set\n";
	//cout << "simListSize " << ss.simList.size() << '\n';
	for (int i = 0; i < ss.simList.size(); i++)
	{
		simulation &sim = ss.simList[i];
		seed.geneCopy(sim);
		ifstream ifsGen(genFitList[i].c_str());
		sim.parseFit(ifsGen);
		string activeAux = genAuxList[i];
		parseRMaux(activeAux);
		ss.simList[i].overrideParameters(1);
		for (int m = 0; m < genActParams[i].size(); m++)
		{
			//cout << "gap: " << genActParams[l].size() << ' ' << m << ' ' << genActParams[l][m] << '\n';
			ss.simList[i].active.insert(genActParams[i][m]);
		}
		cout << "active inserted\n";
	}
	for (int i = 0; i < ss.simList.size(); i++)
	{
		simulation &sim = ss.simList[i];
		ifstream ifs(geneListName.c_str());
		sim.makeGenes(isConden);
		sim.parseGenes(ifs);
		cout << "genes read\n";
		sim.name = seed.name + "_a" + to_string(i);
		sim.phenotype(isConden);
		cout << "phenotyping done\n";
		sim.initDecompAll();
	}
	gp.vecsim.push_back(&ss);
	gp.writeGene = 0;
	genF((void*)&gp);
	cout << geneListName << " done\n";
}
