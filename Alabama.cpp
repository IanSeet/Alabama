#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <map>
#include <set>
#include <pthread.h>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>
#include <queue>
using namespace std;
bool writetocout = 1;
const double PI = 3.14159265359, C = 112.1483235, D = 1.0, CAdj = C/D, Big = 2e128, Smol = -2e128, sqrt3 = sqrt(3), intmax = 2147483647; //Pi, coulomb and dielectric constants
double h = 1e-3, h2, mu = 1.5, minh = 1e-4; //step size for integrator
const double scale = 0.1;
double LJcutoff = 4.0, verletRadius = 4.3, verletSkin = (verletRadius - LJcutoff)/2; //LJ potential cutoff distance
bool wca = 0, noise = 1, multi = 0, roteq = 0, dipoleOn = 1, sinDipole = 0, defDirection = 1, expSum = 0; //default direction is anticlockwise. expSum = 1 when expenditures are summed across each writeout
int minsteps = 20000, eqsteps = 20000, roteqsteps = 1e5, cntOffset = 0;
double massMult = 1, kT = 1.0, stepMult = 1; //massMult multiplies mass and stepMult multiplies the timestep by the specified constant
double wcaCutoff = 1.6, wcaRadius = 1.9, wcaSkin = (wcaRadius - wcaCutoff)/2;  //WCA potential constants
bool expFront = 0; double dipoleMass = 0; //these were for debugging the external dipole, do not modify
double sinDipoleFactor = 1, rotateDip = 0, takeObs = 0;
#include "quaternion.h" //Quaternion and 3d vector functions
vector3d defDipole, defOffDipole, defAxis;
double offScaleGen = 1.0, offPhaseFactorGen = 0.0; bool writePOL = 0;
double phaseObs = 0.0; bool inactivateDum = 1, Fcircular = 0;
#include "Eigen/Eigenvalues" //Finds eigenvalues and eigenvectors of a hermitian 3x3 matrix
using namespace Eigen;
#include "rigidBody.h" //Rigid body functions
rigidBody * RBlistMain;
#include "interactionHarmonic.h" //Constraint, bond, angle and dihedral interactions
vector<string> args;
struct triple
{
	int a, b, c;
	inline triple(){}
	inline triple(int _a, int _b, int _c){a = _a; b = _b; c = _c;}
};

struct genAlgoParam
{
	double mutrate;
	int generations, popsize, maxthreads;
	int steps, runs, eqsteps, roteqsteps;
	bool isConden = 0;
	genAlgoParam();
	genAlgoParam(int _gens, int _popsize, double _mutrate, int _maxthreads, int _steps, int _runs, int _eqsteps, int _roteqsteps):
	generations(_gens), popsize(_popsize), mutrate(_mutrate), maxthreads(_maxthreads), steps(_steps), runs(_runs), eqsteps(_eqsteps), roteqsteps(_roteqsteps){}
};

struct gene
{
	vector<pair<int, int> > vecMass;
	vector<int> vecInt;
	int type;
	double shift;
	vector3d disp;
	gene(){shift = 1; disp.assign(0, 0, 0);}
	
	void print(ofstream &ofs)
	{
		ofs << type << ' ' <<  shift << '\n';
		disp.print(ofs);
		ofs << vecMass.size() << '\n';
		for (int j = 0; j < vecMass.size(); j++)
		{
			ofs << vecMass[j].first << ' ' << vecMass[j].second << '\n';
		}
		ofs << vecInt.size() << '\n';
		for (int j = 0; j < vecInt.size(); j++)
		{
			ofs << vecInt[j] << '\n';
		}
	}
	
	void printVerbose(ofstream &ofs)
	{
		ofs << "type-shift-disp\n";
		ofs << type << ' ' <<  shift << '\n';
		disp.print(ofs);
		ofs << "vecMass: " << vecMass.size() << '\n';
		for (int j = 0; j < vecMass.size(); j++)
		{
			ofs << vecMass[j].first << ' ' << vecMass[j].second << '\n';
		}
		ofs << "vecInt: " << vecInt.size() << '\n';
		for (int j = 0; j < vecInt.size(); j++)
		{
			ofs << vecInt[j] << '\n';
		}
	}
	
	void parse(ifstream &ifs)
	{
		ifs >> type >>  shift;
		disp.read(ifs);
		int a;
		ifs >> a; vecMass.resize(a);
		for (int j = 0; j < vecMass.size(); j++)
		{
			ifs >> vecMass[j].first >> vecMass[j].second;
		}
		ifs >> a; vecInt.resize(a);
		for (int j = 0; j < vecInt.size(); j++)
		{
			ifs >> vecInt[j];
		}
	}
};

inline int massf(rigidBody * RBlist, vector<pair<int, int> > &vp, int i)
{
	return RBlist[vp[i].first].specList[vp[i].second]; 
}

#include "bamaGene.h"

int main(int argc, char *argv[])
{
	int provided;
	int ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	int procid, numprocs;
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	if (argc > 2)
	genAlgo(atoi(argv[1]), atoi(argv[2]), procid, numprocs);
	ierr = MPI_Finalize();
}
