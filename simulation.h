class simulation
{
	//Simulation units: Energy - kT at 298K, Temperature - kelvin, Distance - 0.5nm, Charge - elementary charge
	public:
	int currStep, RBlistSize, dummyListSize, intListSize, conIntListSize, maxTimeStep, chain;
	double contribs[6]; double LJcontrib;
	bool printStep, updateCellList; //If a particle breaches the boundaries of a cell, trigger update
	vector3d minDim, maxDim, boxDim; //Smallest and largest values of x, y and z reached. Needed for grid cell calculations.
	rigidBody * RBlist;
	inline void haltAll()
	{
		for (int i = 0; i < RBlistSize; i++) RBlist[i].halt();
	}
	config optConfig;
	vector<double> expenditureList;
	vector<config> SDtrajectory, MDtrajectory, eqtrajectory; //Stores trajectory
	vector<energyConfig> combinedEnergies; //Averaged energies from multiple runs
	double expenditure; //Energy expenditure
	vector<constrainAxisTime> constrainAxisTimeList;
	vector<bond> bondList;
	
	struct trobs
	{
		vector<int> vi;
		int last;
		trobs()
		{
			last = -1;
		}
	};
	
	vector<interaction*> observables; unordered_map<float, trobs> triggerObs;
	vector< vector<unsigned char> > finalObservables;
	vector<vector<vector<unsigned char> > > finalObservableList;
	vector<vector<map<float, int> > > printObservableList;
	vector<vector<float> > polAverage;
	vector<map<float, int> > poltemp;
	map<string, pair<double, double> > observableMatrix;

	double interactionTime, integratorTime, liTime, cellListTime;
	interaction ** interactionList;
	vector<interaction*> condensedIntList;

	inline int min(int a, int b) {if (a < b) return a; else return b;}
	inline int max(int a, int b) {if (a > b) return a; else return b;}

	double offScale, offPhaseFactor, pauseRotFactor;
	double LJumbrellaEnergy;
	vector3d externalDipole, rotationAxis, initialDipole, prevDipole, offsetDipole, initOffDipole, prevOffDipole; bool clockwise;
	bool isEquil, isMin, isMD;
	string name; int runNo, maxRuns;
	vector<pair<int, int> > topSpecList;
	vector<vector<pair<int, int> > > mtopSpecList;
	unordered_set<int> active; bool isActive;
	double sinDipoleAngle; bool circular;
	double loss;
	int intervalCount;
	double sumTKE = 0, sumRKE = 0;
	int totalKEsample = 0;
	double expectedTKE, expectedRKE;

	map<int, int> massMap;
	vector<pair<int, int> > revMassMap;

	unordered_map<vector3d, vector<int> > symMap;
	
	#include "extDip.h"

	void findObservables()
	{
		observables.clear();
		intervalCount = 0;
		//printObservableList.clear();
		for (int i = 0; i < intListSize; i++)
		{
			if (interactionList[i]->obs != 0)
			{
				double temp;
				if (interactionList[i]->obs < 0) temp = interactionList[i]->obs + interactionList[i]->phase;
				else temp = interactionList[i]->obs - interactionList[i]->phase;
				temp *= 2*PI;
				observables.push_back(interactionList[i]);
				temp = (roundf(temp*1000) + 0.0)/1000;
				if (triggerObs.count(temp) == 0)
				{
					trobs v; v.vi.push_back(observables.size() - 1);
					triggerObs[temp] = v;
				}
				else
				{
					triggerObs[temp].vi.push_back(observables.size() - 1);
				}
			}
		}
		/*unordered_map<float, trobs>::iterator it;
		cout << "iterator: ";
		for (it = triggerObs.begin(); it != triggerObs.end(); ++it) cout << it->first << ' ' << it->second.vi.size() << '\n';
		cout << "obssize " << observables.size() << '\n';*/
		finalObservables.resize(observables.size());
		poltemp.clear(); poltemp.resize(observables.size());
	}
	
	inline void populateMatrix(int w, int step)
	{
		vector<float> vc; vc.resize(observables.size());
		double umbrella = 0;
		for (int i = 0; i < vc.size(); i++)
		{
			vc[i] = observables[i]->observable();
			umbrella += observables[i]->umbrella();
		}
		if (writePOL)
		{
			if (runNo > 1)
			{
				for (int i = 0; i < vc.size(); i++)
				{
					float k = roundf(vc[i]);
					map<float, int> &m = printObservableList[intervalCount][i];
					if (m.count(k) == 0) m[k] = 1;
					else m[k]++;
				}
				if (w == 0)
				{
					intervalCount++;
				}
			}
			else
			{	
				for (int i = 0; i < vc.size(); i++)
				{
					float k = roundf(vc[i]);
					if (poltemp[i].count(k) == 0) poltemp[i][k] = 1;
					else poltemp[i][k]++;
				}
				if (w == 0)
				{
					printObservableList.push_back(poltemp);
					poltemp.clear(); poltemp.resize(vc.size());
				}
			}
		}
		else
		{
			string s;
			for (int i = 0; i < vc.size(); i++)
			{
				char c = (unsigned char)((int)(roundf(vc[i]))<<24);
				s += c;
			}
			double LJumbExp = 1;
			if (LJumbrellaEnergy > 100) LJumbExp = 0;
			else LJumbExp =  exp(-LJumbrellaEnergy);
			if (observableMatrix.count(s) == 0)
			{
				observableMatrix[s].first = LJumbExp;
				observableMatrix[s].second = exp(-umbrella);
			}
			else observableMatrix[s].first += LJumbExp;
		}	
	}
	
	inline void recordFinalObservables()
	{
		int timeFrac;
		if (chain == 1) timeFrac = currStep;
		else
		{
			int chainLength = maxTimeStep/chain, remainder = currStep%chainLength;
			timeFrac = remainder;
		}
		if (triggerObs.count(sinDipoleAngle) == 1)
		{
			trobs &v = triggerObs[sinDipoleAngle];
			if (v.last == -1 || timeFrac - v.last > 0.05*maxTimeStep)
			{
				v.last = timeFrac;
				//cout << "FOUND OBSERVABLE\n";
				for (int i = 0; i < v.vi.size(); i++)
				{
					finalObservables[v.vi[i]].push_back(observables[v.vi[i]]->observable());
					//cout << v.vi[i] << '\n';
				}
			}
		}
	}

	inline void findFinalObservables()
	{
		unordered_map<float, trobs>::iterator it;
		for (it = triggerObs.begin(); it != triggerObs.end(); ++it)
		{
			it->second.last = -1;
		}
		finalObservableList.push_back(finalObservables);
		for (int i = 0; i < finalObservables.size(); i++) finalObservables[i].clear();
	}

	inline void catpInit()
	{
		for (int i = 0; i < constrainAxisTimeList.size(); i++)
		{
			constrainAxisTimeList[i].init();
		}
		externalDipole = initialDipole; prevDipole = initialDipole; expenditure = 0;
	}

	inline double constrainAxisTimePotential(bool calcForce, double adv, double &expend)
	{
		double total = 0.0; expend = 0.0;
		for (int i = 0; i < constrainAxisTimeList.size(); i++)
		{
			constrainAxisTimeList[i].fix(adv);
			expend += constrainAxisTimeList[i].expend;
			total += constrainAxisTimeList[i].potential();
			if (calcForce) constrainAxisTimeList[i].force();
		}
		return total;
	}

	inline double interactionPotential(bool calcForce, bool printStep)
	{
		double total = 0.0;
		for (int i = 0; i < 6; i++) contribs[i] = 0;
		for (int i = 0; i < intListSize; i++)
		{
			//cout << i << ' ' << intListSize << ' ' << interactionList[i]->type << '\n';
			if (isActive && interactionList[i]->active != 0 && active.count(interactionList[i]->active) == 0) {}
			else
			{	
				if (interactionList[i]->type >= 5) interactionList[i]->fix(printStep, isMin, isEquil);
				else interactionList[i]->fix(printStep, isMin, isEquil, sinDipoleAngle);
				contribs[interactionList[i]->type] += interactionList[i]->potential();
				if (calcForce) interactionList[i]->force();
			}
		}
		for (int i = 0; i < 6; i++) total += contribs[i];
		return total;
	}
	#include "pairwise.h" //Lennard-Jones and electrostatic potentials + cell list
	#include "nativeGene.h"
	
	vector<vector<pair<int, int> > > adjList;
	vector<vector<int> > adjList2;
	vector<vector<pair<int, int> > > combinedAdjList;

	inline void decompAll()
	{
		for (int i = 0; i < RBlistSize; i++)
		{
			RBlist[i].decompose();
			RBlist[i].force.assign(0, 0, 0); RBlist[i].torque.assign(0, 0, 0); RBlist[i].resetMoments();
		}
	}

	inline double customInteractions(){return 0.0;} //placeholder for custom interactions
	
	inline void resetExpend()
	{
		for (int i = 0; i < intListSize; i++)
		{
			interactionList[i]->resetExpend();
			interactionList[i]->repointer(RBlist, &expenditure);
		}
	}

	inline double totalPotential(double t)
	{
		double expend;
		double out = interactionPotential(0, 0) + constrainAxisTimePotential(0, t, expend) 
				+ LJpotential(0) + electrostatic(0) + customInteractions();
		return out;
	}

	inline void printContributions(double t)
	{
		double expend;
		cout << "constraint: " << contribs[0] << '\n';
		cout << "constrainAxis: " << contribs[1] << '\n';
		cout << "constrainAxisTime: " << constrainAxisTimePotential(0, t, expend) << '\n';
		cout << "bond: " << contribs[2] << '\n';
		cout << "angle: " << contribs[3] << '\n';
		cout << "dihedral: " << contribs[4] << '\n';
		cout << "external dipole: " << contribs[5] << '\n';
		cout << "dispersion: " << LJcontrib << '\n';
		cout << "electrostatic: " << electrostatic(0) << '\n';
	}

	inline double totalForce(double t, bool printStep)
	{
		double expend = 0;
		double iP = interactionPotential(1, printStep);
		double cATP = constrainAxisTimePotential(1, t, expend), eP = electrostatic(1);
		double ljP = LJpotential(1);
		for (int i = 0; i < RBlistSize; i++) RBlist[i].calcTorque();
		expenditure += expend;
		return iP + cATP + ljP + eP + customInteractions();
	}

	inline double SP() //Single point energy
	{
		initDecompAll();
		return totalPotential(0);
	}

	void printBonds()
	{
		cout << "bonds: " << bondList.size() << '\n';
		for (int i = 0; i < bondList.size(); i++) bondList[i].printtocout();
	}
	void initialise()
	{
		for (int i = 0; i < intListSize; i++)
		{
			if (interactionList[i]->type == 5 && interactionList[i]->subtype == 3)
			{
				//cout << i << '\n';
				staccatoDipole * sd = (staccatoDipole*)interactionList[i];
				sd->findAxes(sinDipoleFactor, chain, rotationAxis, initialDipole, clockwise);
			}
		}
	}
	#include "nativeParse.h"
	
	void decondenseIntList()
	{
		cout << "cilsize " << condensedIntList.size() << '\n';
		if (condensedIntList.size() == 0) return;
		if (intListSize > 0)
		{
			return;
			for (int i = condensedIntList.size(); i < intListSize; i++)
			{
				delete interactionList[i];
			}
			delete [] interactionList;
		}
		bondGen();
		combineAdjLists();
		vector<interaction*> nangle;
		for (int i = 0; i < condensedIntList.size(); i++)
		{
			//cout << condensedIntList.size() << ' ' << i << '\n';
			if (condensedIntList[i]->type == 2)
			{
				bond *b = (bond*)condensedIntList[i];
				if (!b->isChain && b->ground > 0)
				{
					pair<int, int> p1, p2;
					p1 = b->px[0];
					p2 = b->px[1];
					int hashf1 = ((p1.first + 1) << 16) + p1.second, hashf2 = ((p2.first + 1) << 16) + p2.second;
					
					vector3d t1 = RBlist[b->px[0].first].massList[b->px[0].second].decomp, 
					t2 = RBlist[b->px[1].first].massList[b->px[1].second].decomp;
					
					for (int j = 0; j < combinedAdjList[massMap[hashf1]].size(); j++)
					{
						if (combinedAdjList[massMap[hashf1]][j].first != massMap[hashf2])
						{
							pair <int, int> p = revMassMap[combinedAdjList[massMap[hashf1]][j].first];
							vector3d v = RBlist[p.first].massList[p.second].decomp;
							double angleGround = calcAngle(v, t1, t2)*180/PI;
							//cout << "bangletypeFirst: " << b->angleType.first << '\n';
							if (b->angleType.first >= 0) angleGround = b->angleType.first;
							angle * a = new angle(p, b->px[0], b->px[1], angleGround, 210, RBlist);
							a->isGene = b->isGene;
							nangle.push_back(a);
						}
					}
					
					for (int j = 0; j < combinedAdjList[massMap[hashf2]].size(); j++)
					{
						if (combinedAdjList[massMap[hashf2]][j].first != massMap[hashf1])
						{
							pair <int, int> p = revMassMap[combinedAdjList[massMap[hashf2]][j].first];
							vector3d v = RBlist[p.first].massList[p.second].decomp;
							double angleGround = calcAngle(t1, t2, v)*180/PI;
							//cout << "bangletypeSecond: " << b->angleType.second << '\n';
							if (b->angleType.second >= 0) angleGround = b->angleType.second;
							angle * a = new angle(b->px[0], b->px[1], p, angleGround, 210, RBlist);
							a->isGene = b->isGene;
							nangle.push_back(a);
						}
					}
					
					rigidBody A = RBlist[b->px[0].first], B = RBlist[b->px[1].first];
					if (A.scaffold || B.scaffold)
					{
						b->ground = dist(t1, t2);
					}
				}
			}
		}
		for (int i = 0; i < RBlistSize; i++)
		{
			if (RBlist[i].scaffold)
			{
				for (int j = 0; j < RBlist[i].massList.size(); j++)
				{
					pair<int, int> p(i, j);
					constraint * c = new constraint(p, RBlist[i].massList[j].decomp, 0.0, 13000, RBlist); nangle.push_back(c);
				}
			}
		}
		intListSize = nangle.size() + condensedIntList.size();
		interactionList = new interaction * [intListSize];
		for (int i = 0; i < condensedIntList.size(); i++)
		{
			interactionList[i] = condensedIntList[i];
		}
		for (int i = 0; i < nangle.size(); i++)
		{
			interactionList[i + condensedIntList.size()] = nangle[i];
		}
		//cout << "intListPtr after: " << interactionList << '\n';
		//cout << "intListSize after: " << intListSize << '\n';
	}
	//#include "thread.h"

	void mirrorSym(vector3d &reflect)
	{
		const double threshold = 1e3;
		map<int, int>::iterator it;
		vector3d systemCenter(0, 0, 0);
		vector<vector3d> adjustedCoords(revMassMap.size());
		for (it = massMap.begin(); it != massMap.end(); ++it)
		{
			pair<int, int> &p = revMassMap[it->second];
			Mass &m = RBlist[p.first].massList[p.second];
			vector3d coords = m.initial + RBlist[p.first].center;
			systemCenter += coords;
		}
		systemCenter /= revMassMap.size();
		for (it = massMap.begin(); it != massMap.end(); ++it)
		{
			pair<int, int> &p = revMassMap[it->second];
			Mass &m = RBlist[p.first].massList[p.second];
			vector3d coords = m.initial + RBlist[p.first].center;
			coords -= systemCenter; coords.round(threshold);
			adjustedCoords[it->second] = coords;
		}
		for (int i = 0; i < adjustedCoords.size(); i++)
		{
			vector3d &coords = adjustedCoords[i];
			double project = coords * reflect;
			if (project < 0) coords -= 2*project*reflect;
			if (symMap.count(coords) == 0)
			{
				vector<int> v(1); v[0] = i; symMap[coords] = v;
			}
			else symMap[coords].push_back(i);
		}
	}

	#include "langevin.h" //langevin thermostat
	#include "top.h"
	#include "costf.h"
	
	void expectedKEfinder()
	{
		int nonDummies = RBlistSize;
		for (int i = 0; i < RBlistSize; i++) if (RBlist[i].isDummy) nonDummies--;
		expectedTKE = 1.5*nonDummies; expectedRKE = expectedTKE;
	}
	
	void KEfinder()
	{
		expectedKEfinder();
		sumTKE /= totalKEsample; sumRKE /= totalKEsample;
		cout << "Total TKE excess: " << sumTKE << ' ' << sumTKE - expectedTKE << '\n';
		cout << "Total RKE excess: " << sumRKE << ' ' << sumRKE - expectedRKE << '\n';
	}
	
	simulation()
	{
		printStep = 0; expenditure = 0; intListSize = 0; runNo = 1; maxRuns = 1;
		updateCellList = 1; chain = 1;
		assignLJparam();	
		interactionTime = 0, integratorTime = 0, liTime = 0, cellListTime = 0;
		interactionList = NULL; RBlist = NULL;
		initialDipole = defDipole; rotationAxis = defAxis; clockwise = defDirection; 
		initialDipole.norm(); rotationAxis.norm();
		externalDipole = initialDipole; prevDipole = initialDipole;
		offScale = offScaleGen; offPhaseFactor = offPhaseFactorGen; pauseRotFactor = pauseRotFactorGen;
		initOffDipole = defOffDipole;
		offsetDipole = initOffDipole; prevOffDipole = offsetDipole;
		isEquil = 0; isMin = 0; isMD = 0; isActive = 0; sinDipoleAngle = 0; circular = 0;
	}
	
	void overrideParameters(bool preserveRB)
	{
		if (!preserveRB)
		{
			printStep = 0; expenditure = 0; runNo = 1; maxRuns = 1;
			updateCellList = 1; chain = 1;
			assignLJparam();
			interactionTime = 0, integratorTime = 0, liTime = 0, cellListTime = 0;	
			interactionList = NULL; RBlist = NULL;
			intListSize = 0;
			isEquil = 0; isMin = 0; isMD = 0; isActive = 0;
		}
		initialDipole = defDipole; rotationAxis = defAxis; clockwise = defDirection;
		if (rotateDip != 0)
		{
			quaternion q(rotationAxis, rotateDip);
			initialDipole.rotate(q);
			
		}
		initialDipole.norm(); rotationAxis.norm();
		externalDipole = initialDipole; prevDipole = initialDipole;
		offScale = offScaleGen; offPhaseFactor = offPhaseFactorGen; pauseRotFactor = pauseRotFactorGen;
		initOffDipole = defOffDipole;
		offsetDipole = initOffDipole; prevOffDipole = offsetDipole;
		sinDipoleAngle = 0; circular = Fcircular;
	}
	
	void deepCopy(simulation &target)
	{
		target = *this;
		target.RBlist = new rigidBody[RBlistSize];
		for (int i = 0; i < RBlistSize; i++)
		{
			target.RBlist[i] = RBlist[i]; target.RBlist[i].repointer(&target.minDim, &target.maxDim, &target.updateCellList);
		}
		target.interactionList = new interaction*[intListSize];
		string fname = "temp"; ofstream ofs(fname.c_str());
		for (int i = 0; i < intListSize; i++)
		{
			int type = stoi(to_string(interactionList[i]->type) + to_string(interactionList[i]->subtype));
			ofs << i << ' ' << type << ' ' << intListSize << '\n';
			switch (type)
			{
				case 20:
					target.interactionList[i] = new bond;
					break;
				case 30:
					target.interactionList[i] = new angle;
					break;
				case 31:
					target.interactionList[i] = new angleAxis;
					*((angleAxis*)target.interactionList[i]) = *((angleAxis*)interactionList[i]);
					break;
				case 40:
					target.interactionList[i] = new dihedral;
					break;
				case 0:
					target.interactionList[i] = new constraint;
					*((constraint*)target.interactionList[i]) = *((constraint*)interactionList[i]);
					break;
				case 10:
					target.interactionList[i] = new constrainAxis;
					*((constrainAxis*)target.interactionList[i]) = *((constrainAxis*)interactionList[i]);
					target.interactionList[i]->redist();
					break;
				case 50:
					target.interactionList[i] = new constrainAxisTime;
					*((constrainAxisTime*)target.interactionList[i]) = *((constrainAxisTime*)interactionList[i]);
					target.interactionList[i]->redist();
					break;
				case 51:
					target.interactionList[i] = new extDipole;
					*((extDipole*)target.interactionList[i]) = *((extDipole*)interactionList[i]);
					target.interactionList[i]->redist();
					target.interactionList[i]->repointer2(&target.externalDipole, &target.prevDipole, &target.expenditure, &target.sinDipoleAngle);
					break;
				case 52:
					target.interactionList[i] = new offDipole;
					*((offDipole*)target.interactionList[i]) = *((offDipole*)interactionList[i]);
					target.interactionList[i]->redist();
					target.interactionList[i]->repointer2(&target.offsetDipole, &target.prevOffDipole, &target.expenditure, &target.sinDipoleAngle);
					break;
				case 53:
					target.interactionList[i] = new staccatoDipole;
					*((staccatoDipole*)target.interactionList[i]) = *((staccatoDipole*)interactionList[i]);
					target.interactionList[i]->redist();
					target.interactionList[i]->repointer2(&target.externalDipole, &target.prevDipole, &target.expenditure, &target.sinDipoleAngle);
					break;
			}
			if (type > 10 && type < 50 && type != 31) *(target.interactionList[i]) = *(interactionList[i]);
			target.interactionList[i]->repointer(target.RBlist, &target.expenditure);
			interactionList[i]->print(ofs);
			target.interactionList[i]->print(ofs);
			//cout << "original: " << (this->interactionList[i])->ground << ' '  << (this->interactionList[i])->k << ' ' << (this->interactionList[i])->d << '\n';
			//cout << "copy: " << (target.interactionList[i])->ground << ' '  << (target.interactionList[i])->k << ' ' << (target.interactionList[i])->d << '\n';
		}
	}
	void geneCopy(simulation &target)
	{
		target = *this;
		target.RBlist = new rigidBody[RBlistSize];
		for (int i = 0; i < RBlistSize; i++)
		{
			target.RBlist[i] = RBlist[i]; target.RBlist[i].repointer(&target.minDim, &target.maxDim, &target.updateCellList);
		}
		target.condensedIntList.resize(this->condensedIntList.size());
		string fname = "temp"; ofstream ofs(fname.c_str());
		for (int i = 0; i < condensedIntList.size(); i++)
		{
			int type = stoi(to_string(condensedIntList[i]->type) + to_string(condensedIntList[i]->subtype));
			ofs << i << ' ' << type << ' ' << condensedIntList.size() << '\n';
			switch (type)
			{
				case 20:
					target.condensedIntList[i] = new bond;
					break;
				case 30:
					target.condensedIntList[i] = new angle;
					break;
				case 31:
					target.condensedIntList[i] = new angleAxis;
					*((angleAxis*)target.condensedIntList[i]) = *((angleAxis*)condensedIntList[i]);
					break;
				case 40:
					target.condensedIntList[i] = new dihedral;
					break;
				case 0:
					target.condensedIntList[i] = new constraint;
					*((constraint*)target.condensedIntList[i]) = *((constraint*)condensedIntList[i]);
					break;
				case 10:
					target.condensedIntList[i] = new constrainAxis;
					*((constrainAxis*)target.condensedIntList[i]) = *((constrainAxis*)condensedIntList[i]);
					target.condensedIntList[i]->redist();
					break;
				case 50:
					target.condensedIntList[i] = new constrainAxisTime;
					*((constrainAxisTime*)target.condensedIntList[i]) = *((constrainAxisTime*)condensedIntList[i]);
					target.condensedIntList[i]->redist();
					break;
				case 51:
					target.condensedIntList[i] = new extDipole;
					*((extDipole*)target.condensedIntList[i]) = *((extDipole*)condensedIntList[i]);
					target.condensedIntList[i]->redist();
					target.condensedIntList[i]->repointer2(&target.externalDipole, &target.prevDipole, &target.expenditure, &target.sinDipoleAngle);
					break;
				case 52:
					target.condensedIntList[i] = new offDipole;
					*((offDipole*)target.condensedIntList[i]) = *((offDipole*)condensedIntList[i]);
					target.condensedIntList[i]->redist();
					target.condensedIntList[i]->repointer2(&target.offsetDipole, &target.prevOffDipole, &target.expenditure, &target.sinDipoleAngle);
					break;
				case 53:
					target.condensedIntList[i] = new staccatoDipole;
					*((staccatoDipole*)target.condensedIntList[i]) = *((staccatoDipole*)condensedIntList[i]);
					target.condensedIntList[i]->redist();
					target.condensedIntList[i]->repointer2(&target.externalDipole, &target.prevDipole, &target.expenditure, &target.sinDipoleAngle);
					break;
			}
			if (type > 10 && type < 50 && type != 31) *(target.condensedIntList[i]) = *(condensedIntList[i]);
			target.condensedIntList[i]->repointer(target.RBlist, &target.expenditure);
			condensedIntList[i]->print(ofs);
			target.condensedIntList[i]->print(ofs);
		}
		target.intListSize = 0; target.interactionList = NULL;
		//cout << "intListSize " << target.intListSize << '\n';
		//cout << "intListPtr " << target.interactionList << '\n';
	}
	
	~simulation()
	{
		//cout << "deleting\n";
		/*cout << intListSize << ' ' << interactionList << '\n';
		for (int i = 0; i < intListSize; i++) 
		{
			cout << i << ' ' << interactionList[i] << '\n';
			delete interactionList[i];
		}*/
		//delete [] interactionList;
		//cout << "deleted intList\n";
		//delete [] RBlist;
		//cout << "deleted\n";
	}
};
