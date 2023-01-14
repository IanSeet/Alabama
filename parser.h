vector3d empty(0, 0, 0);
#include "oldParser.h"
#include "conParser.h"
unordered_map <string, vector<int> > ums;

void parseInputCommon(string &ftype, string &fin, simulation &sim)
{
	sim.bondGen();
	cout << "bondgen complete\n";
	string s, s2;
	ofstream ofs;
	bool changeType = 1;
	//bool b = sim.colouring(3, sim.adjList2);
	//cout << "colouring done\n";
	switch (ftype[0])
	{
		case 'g':
			s = fin + ".gjf";
			ofs.open(s.c_str());
			sim.togjf(ofs);
			break;
		case 'p':
			s = fin + ".pdb";
			ofs.open(s.c_str());
			sim.topdb(ofs);
			break;
		case 't':
			s = fin + ".top";
			s2 = fin + ".conf";
			ofs.open(s.c_str());
			sim.printTop(ofs);
			ofs.close();
			ofs.open(s2.c_str());
			sim.printConfig(ofs);
			break;
		case 'r':
			s = fin + "_rigid.mol2";
			ofs.open(s.c_str());
			sim.printTopolMol2(ofs, 1, changeType);
			ofs.close();
			s = fin + "_bridge.mol2";
			ofs.open(s.c_str());
			sim.printTopolMol2(ofs, 2, changeType);
			ofs.close();
			s = fin + "_tether.mol2";
			ofs.open(s.c_str());
			sim.printTopolMol2(ofs, 3, changeType);
			cout << "printed rigid and bridge\n";
	}
	s = fin + ".mol2";
	ofstream ofs2(s.c_str());
	sim.printTopolMol2(ofs2, 0, changeType);
	ofs2.close();
}

void MDCexpenditure(ofstream &ofs, simulation &sim)
{
	for (int i = 0; i < sim.MDtrajectory.size(); i++)
	{
		ofs << i << ' ' << sim.MDtrajectory[i].TKE << ' ' << sim.MDtrajectory[i].RKE << ' ' << sim.MDtrajectory[i].V << ' ' <<
		sim.MDtrajectory[i].TKE + sim.MDtrajectory[i].RKE + sim.MDtrajectory[i].V << ' ' << sim.MDtrajectory[i].expenditure << '\n';
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

void parseActive(string &Active, string &Active2, simulation &sim)
{
	ifstream active(Active.c_str());
	string s; int x;
	while (getline(active, s))
	{
		istringstream iss(s);
		iss >> s;
		if (s == Active2)
		{
			while (iss >> x) sim.active.insert(x);
			break;
		}
	}
}

void parseInput(ifstream &ifs)
{
	string s, fin, ftype, calc, fout; int idx = 0; bool isTop = 0;
	int cnt = cntOffset, multiple = 0;
	while (getline(ifs, s))
	{
		multiple++;
	}
	ifs.clear(); ifs.seekg(0);
	vecsim.resize(multiple);
	while (getline(ifs, s))
	{
		cout << s << '\n';
		istringstream iss(s);
		int steps, write, runs = 1;
		string RMaux = "RMaux", fout = "Default", active = "Default";
		iss >> s >> fin >> ftype >> calc >> steps >> write >> runs >> RMaux >> active >> fout;
		parseRMaux(RMaux);
		ifstream ift(fin.c_str());
		vector3d v(0, 0, 0); quaternion q(1, 0, 0, 0);
		simulation &sim = vecsim[idx]; idx++;
		sim.overrideParameters(0);
		if (active != "Default") parseActive(active, fout, sim);
		triple counter(0, 0, 0);
		cout << "Parameters overridden...\n";
		bool fileExist = 1;
		
		if (s == "new") 
		{
			if (ift.fail()) fileExist = 0;
			else parseFrag(ift, v, 0, q, sim, 1, counter, steps);
		}
		else if (s == "super")
		{
			if (ift.fail()) fileExist = 0;
			else parseSuper(ift, sim, steps);
		}
		else if (s == "conden") 
		{
			if (ift.fail()) fileExist = 0;
			parseFragCondensed(ift, v, 0, q, sim, 1, counter, steps);
			if (calc != "G2") sim.decondenseIntList();
			cout << "intListSize " << sim.intListSize << '\n';
		}
		else if (s == "top")
		{
			string stop = fin + ".top";
			string sconf = fin + ".conf";
			ifstream if1(stop.c_str()), if2(sconf.c_str());
			if (if1.fail() || if2.fail()) fileExist = 0;
			else sim.parseTop(if1, if2); //return;
			isTop = 1;
		}
		else if (s == "mtop")
		{
			if (ift.fail()) fileExist = 0;
			else sim.parseMultiTop(ift);
		}
		
		if (fileExist)
		{
			//cout << "fexist\n";
			sim.initialise();
			//cout << "initialised\n";
			sim.findObservables();
			//cout << "Observables found...\n";
			RBlistMain = sim.RBlist;
			cout << "Parsing completed...\n";
			if (fout == "Default") fout = fin;
			if (multiple > 1  && cntOffset >= 0)
			{
				ostringstream oss; oss << cnt;
				fout += "_" + oss.str();
			}
			sim.name = fin;
			sim.resetExpend();
			if (calc == "SP")
			{
				cout << sim.SP() << '\n';
				sim.printContributions(0);
				parseInputCommon(ftype, fout, sim);
			}
			else if (calc == "test")
			{
				/*ifstream ifs("seed");
				sim.parseGenes(ifs);
				sim.phenotype(0);*/
				
				for (int i = 0; i < sim.RBlistSize; i++)
				{
					cout << "curr RB: " << i << '\n';
					sim.RBlist[i].recenter();
				}
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
					if (writePOL) sim.findObservables();
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
				s = fout + "Dipoles";
				ofstream ofs5(s.c_str());
				sim.printExtDipVectorsFrame(ofs5, 0, sim.MDtrajectory);
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
			else if (calc == "G2")
			{
				genAlgoParam gap(1, 1, 0.5, 1, 100000, 1, 50000, 0);
				string genauxName("genAux");
				parseGenAux(genauxName, gap);
				genAlgo(sim, gap);
			}
			else if (calc == "PG")
			{
				genAlgoParam gap(1, 1, 0.5, 1, 100000, 1, 50000, 0);
				string genauxName("genAux");
				parseGenAux(genauxName, gap);
				cout << "isConden: " << gap.isConden << '\n';
				cout << "genAux parsed\n";
				cout << "intListSize: " << sim.intListSize << '\n';
				cout << "cilsize: " << sim.condensedIntList.size() << '\n';
				if (args.size() > 2) 
				{
					cout << args.size() << ' ' << args[1] << ' ' << args[2] << '\n';
					if (args[2] == "seed")
					{
						sim.makeGenes(gap.isConden);
						ofstream ofseed("seed");
						sim.printGenes(ofseed);
					}
					else
					{
						cout << args[3] << ' ' << args[4] << ' ' << args[6] << '\n';
						bool isElite = 0; if (args[6] == "1") isElite = 1;
						//cout << "isElite: " << args[5] << ' ' << isElite << '\n';
						pyGenAlgo(sim, gap, args[3], args[4], isElite, isTop);
					}
				}
			}
			else if (calc == "RG")
			{
				genAlgoParam gap(1, 1, 0.5, 1, 100000, 1, 50000, 0);
				string genauxName("genAux");
				parseGenAux(genauxName, gap);
				readGenAlgo(sim, gap, args[2]);
			}
			cnt++;
		}
	}
}
