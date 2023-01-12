static void *multiThreadF(void *p)
{
	string *sptr = (string*)p;
	istringstream iss(*sptr);
	cout << *sptr << '\n';
	string s, fin, ftype, calc; int no, steps, write, runs = 1;
	string RMaux = "RMaux", fout = "Default", active = "Default";
	iss >> no >> s >> fin >> ftype >> calc >> steps >> write >> runs >> RMaux >> active >> fout;
	int seed = time(NULL) + no;
	srand(seed);
	parseRMaux(RMaux);
	ifstream ift(fin.c_str());
	vector3d v(0, 0, 0); quaternion q(1, 0, 0, 0);
	simulation &sim = vecsim[no];
	sim.overrideParameters(0);
	if (active != "Default") parseActive(active, fout, sim);
	triple counter(0, 0, 0);
	cout << "Parameters overridden...\n";
	if (s == "new") parseFrag(ift, v, 0, q, sim, 1, counter, steps);
	else if (s == "super") parseSuper(ift, sim, steps);
	else if (s == "top")
	{
		string stop = fin + ".top";
		string sconf = fin + ".conf";
		ifstream if1(stop.c_str()), if2(sconf.c_str());
		sim.parseTop(if1, if2); //return;
	}
	else if (s == "mtop")
	{
		sim.parseMultiTop(ift);
	}
	sim.initialise();
	sim.findObservables();
	//cout << "Observables found...\n";
	RBlistMain = sim.RBlist;
	cout << "Parsing completed...\n";
	if (fout == "Default") fout = fin;
	sim.name = fin;
	if (calc == "SP")
	{
		cout << sim.SP() << '\n';
		sim.printContributions(0);
		parseInputCommon(ftype, fout, sim);
	}
	else if (calc == "SD")
	{
		sim.SDmin(steps);
		parseInputCommon(ftype, fout, sim);
	}
	else if (calc == "MD")
	{
		sim.SDmin(minsteps);
		sim.equilibriate(eqsteps, eqsteps/10);
		if (roteq)
		{
			sim.roteq(roteqsteps, roteqsteps/10);
		}
		sim.areq(eqsteps, eqsteps/10);
		parseInputCommon(ftype, fout, sim);
		sim.langevinIntegrator(steps, write, 0);
		s = fout + ".xyz";
		ofstream ofs(s.c_str());
		for (int i = 0; i < sim.MDtrajectory.size(); i++)
		{
			sim.printConfigxyz(i, sim.MDtrajectory, ofs); 
		}
		s = fout + "Energy";
		ofstream ofs2(s.c_str());
		for (int i = 0; i < sim.MDtrajectory.size(); i++)
		{
			ofs2 << i << ' ' << sim.MDtrajectory[i].TKE << ' ' << sim.MDtrajectory[i].RKE << ' ' << sim.MDtrajectory[i].V << ' ' <<
			sim.MDtrajectory[i].TKE + sim.MDtrajectory[i].RKE + sim.MDtrajectory[i].V << ' ' << sim.MDtrajectory[i].expenditure << ' ' <<
			sim.MDtrajectory[i].TKE + sim.MDtrajectory[i].RKE + sim.MDtrajectory[i].V - sim.MDtrajectory[i].expenditure << '\n';
		}
		s = fout + "Observables";
		ofstream ofs3(s.c_str());
		sim.printObservables(ofs3, writePOL);
	}
	else if (calc == "MDR")
	{
		sim.SDmin(minsteps);
		sim.equilibriate(eqsteps, eqsteps/10);
		if (roteq)
		{
			sim.roteq(roteqsteps, roteqsteps/10);
		}
		sim.optConfig.init(sim.RBlist, sim.RBlistSize);
		sim.optConfig.assign();
		sim.combinedEnergies.resize(steps/write + 1);
		parseInputCommon(ftype, fout, sim);
		sim.maxRuns = runs;
		for (int i = 0; i < runs; i++)
		{
			sim.runNo = i + 1;
			sim.optConfig.overwrite(); sim.catpInit();
			sim.areq(eqsteps, eqsteps/10);
			sim.langevinIntegrator(steps, write, 0);
			for (int j = 0; j < sim.MDtrajectory.size(); j++) sim.combinedEnergies[j].add(sim.MDtrajectory[j]);
			sim.expenditureList.push_back(sim.expenditure);
			//cout << "lossF: " << sim.lossFunction(0.25, 0.75, -30, 90.0, 0.1, 0, sim.expenditureList.size() - 1) << '\n';
			sim.findObservables();
		}
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
	}
	else if (calc == "MDC")
	{
		sim.SDmin(minsteps);
		sim.equilibriate(minsteps, minsteps/10);
		if (roteq)
		{
			sim.roteq(roteqsteps, roteqsteps/10);
		}
		sim.areq(eqsteps, eqsteps/10);
		sim.chain = runs;
		sim.langevinIntegrator(steps, write, 1);
		parseInputCommon(ftype, fout, sim);
		s = fout + ".xyz";
		ofstream ofs(s.c_str());
		for (int i = 0; i < sim.MDtrajectory.size(); i++)
		{
			sim.printConfigxyz(i, sim.MDtrajectory, ofs); 
		}
		s = fout + "Energy";
		ofstream ofs2(s.c_str());
		for (int i = 0; i < sim.MDtrajectory.size(); i++)
		{
			ofs2 << i << ' ' << sim.MDtrajectory[i].TKE << ' ' << sim.MDtrajectory[i].RKE << ' ' << sim.MDtrajectory[i].V << ' ' <<
			sim.MDtrajectory[i].TKE + sim.MDtrajectory[i].RKE + sim.MDtrajectory[i].V << ' ' << sim.MDtrajectory[i].expenditure << '\n';
		}
		s = fout + "Observables";
		ofstream ofs3(s.c_str());
		sim.printObservables(ofs3, writePOL);
		s = fout + "Expenditures";
		ofstream ofs4(s.c_str());
		sim.printExpenditures(ofs4);
	}
	else if (calc == "MDCR")
	{
		sim.maxRuns = runs;
		sim.optConfig.init(sim.RBlist, sim.RBlistSize);
		sim.optConfig.assign();
		for (int i = 0; i < runs; i++)
		{
			sim.optConfig.overwrite();
			sim.SDmin(minsteps);
			sim.equilibriate(eqsteps, eqsteps/10);
			if (roteq)
			{
				sim.roteq(roteqsteps, roteqsteps/10);
			}
			sim.areq(eqsteps, eqsteps/10);
			sim.combinedEnergies.resize(steps/write + 1);
			parseInputCommon(ftype, fout, sim);
			sim.runNo = i + 1;
			sim.catpInit();
			sim.chain = 10;
			sim.langevinIntegrator(steps, write, 1);
			for (int j = 0; j < sim.MDtrajectory.size(); j++) sim.combinedEnergies[j].add(sim.MDtrajectory[j]);
		}
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
	}
}

vector<pthread_t> MDthreads;

void printObservablesMulti(ofstream &ofs, map<string, int> &observableMatrix)
{
	map<string, int>::iterator it;
	for (it = observableMatrix.begin(); it != observableMatrix.end(); ++it)
	{
		for (int j = 0; j < it->first.length(); j++) ofs << (int)((unsigned char)it->first[j]) << ' ';
		ofs << it->second << '\n';
	}
}

void combinePOLMulti(vector<simulation> &vecsim, vector<vector<map<float, int> > > &pol)
{
	map<float, int>::iterator it;
	for (int i = 0; i < vecsim.size(); i++)
	{
		vector<vector<map<float, int> > > &polv = vecsim[i].printObservableList;
		if (!i)
		{
			pol.resize(polv.size());
			for (int j = 0; j < pol.size(); j++) pol[j].resize(polv[j].size());
		}
		for (int j = 0; j < polv.size(); j++) for (int k = 0; k < polv[j].size(); k++)
		{
			for (it = polv[j][k].begin(); it != polv[j][k].end(); ++it)
			{
				if (pol[j][k].count(it->first) == 0) pol[j][k][it->first] = it->second;
				else pol[j][k][it->first] += it->second;
			}
		}
	}
}

void printPOLMulti(ofstream &ofs, vector<vector<map<float, int> > > &printObservableList)
{
	
	for (int i = 0; i < printObservableList.size(); i++)
	{
		ofs << "Timestep No: " << i << '\n';
		for (int j = 0; j < printObservableList[i].size(); j++)
		{
			//cout << "pol: " << i << ' ' << j << '\n';
			ofs << "Observable No: " << j << '\n';
			map<float, int>::iterator it;
			for (it = printObservableList[i][j].begin(); it != printObservableList[i][j].end(); ++it)
			{
				ofs << it->first << ": " << it->second << '\n';
			}
			ofs << '\n';
		}
		ofs << '\n';
	}
}
	
void printInputMulti()
{
	string s = "gjf", name;
	unordered_map <string, vector<int> >::iterator it;
	for (it = ums.begin(); it != ums.end(); ++it)
	{
		vector<int> &v = it->second; name = it->first;
		simulation &sim = vecsim[v[0]];
		parseInputCommon(s, name, sim);
		if (sim.MDtrajectory.size() != 0)
		{
			vector<energyConfig> vcg(sim.MDtrajectory.size());
			for (int i = 0; i < v.size(); i++) for (int j = 0; j < vcg.size(); j++)
			{
				vcg[j].add(vecsim[v[i]].MDtrajectory[j]);
			}
			
			s = name + ".xyz";
			ofstream ofs(s.c_str());
			for (int i = 0; i < sim.MDtrajectory.size(); i++)
			{
				sim.printConfigxyz(i, sim.MDtrajectory, ofs); 
			}
			s = name + "Energy";
			ofstream ofs2(s.c_str());
			for (int i = 0; i < vcg.size(); i++)
			{
				vcg[i].divide(v.size());
				ofs2 << i << ' ' << vcg[i].TKE << ' ' << vcg[i].RKE << ' ' << vcg[i].V << ' ' << vcg[i].TKE + vcg[i].RKE + vcg[i].V << 
				' ' << vcg[i].expenditure << '\n';
			}
			if (!writePOL)
			{
				map<string, int> combinedObservables;
				map<string, pair<double, double> >::iterator it;
				for (int i = 0; i < v.size(); i++)
				{
					map<string, pair<double, double> > &om = vecsim[v[i]].observableMatrix;
					for (it = om.begin(); it != om.end(); ++it)
					{
						if (combinedObservables[it->first] == 0) combinedObservables[it->first] = it->second.first;
						else combinedObservables[it->first] += it->second.first;
					}
				}
				s = name + "Observables";
				ofstream ofs3(s.c_str());
				printObservablesMulti(ofs3, combinedObservables);
			}
			else
			{
				vector<vector<map<float, int> > > pol;
				combinePOLMulti(vecsim, pol);
				s = name + "Observables";
				ofstream ofs3(s.c_str());
				printPOLMulti(ofs3, pol);
			}
			for (int i = 1; i < v.size(); i++) for (int j = 0; j < vecsim[v[i]].expenditureList.size(); j++) sim.expenditureList.push_back(vecsim[v[i]].expenditureList[j]);
			s = name + "Expenditures";
			ofstream ofs4(s.c_str());
			//sim.printExpenditures(ofs4);
		}
	}
}

void parseInputMulti(ifstream &ifs)
{
	string s; int i = 0, sum = 0;
	vector<string> vs;
	while (getline(ifs, s))
	{
		ostringstream oss; oss << sum;
		string x = oss.str(), y; x += " ";
		vs.push_back(x + s);
		istringstream iss(s);
		iss >> y >> y;
		if (ums.count(y) == 0)
		{
			vector<int> v; v.assign(1, sum); ums[y] = v;
		}
		else ums[y].push_back(sum);
		sum++;
	}
	MDthreads.resize(sum);
	vecsim.resize(sum); //sum--;
	for (int i = 0; i < sum; i++)
	{
		int rc = pthread_create(&MDthreads[i], NULL, multiThreadF, &vs[i]);
		//cout << &vs[i] << '\n';
  		if (rc)
		{
			cout << "Error: unable to create thread " << rc << '\n';
			exit(-1);
  		}
	}
	for (int i = 0; i < sum; i++) 
	{
		pthread_join(MDthreads[i], NULL);
	}
	printInputMulti();
}
