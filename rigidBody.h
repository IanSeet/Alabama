struct state 
{
	vector3d center, momentum, force, torque;
	quaternion orient, angmom;
	inline state(){}
	inline state(vector3d& _center, vector3d& _force, vector3d& _torque, quaternion& _orient)
	{
		center = _center; force = _force; torque = _torque; orient = _orient; momentum.assign(0, 0, 0); angmom.assign(1, 0, 0, 0);
	}
	inline void assign(vector3d& _center, vector3d& _force, vector3d& _torque, quaternion& _orient)
	{
		center = _center; force = _force; torque = _torque; orient = _orient; momentum.assign(0, 0, 0); angmom.assign(1, 0, 0, 0);
	}
	inline state(vector3d& _center, vector3d& _force, vector3d& _torque, quaternion& _orient, vector3d& _momentum, quaternion& _angmom)
	{
		center = _center; force = _force; torque = _torque; orient = _orient; momentum = _momentum, angmom = _angmom;
	}
	inline void assign(vector3d& _center, vector3d& _force, vector3d& _torque, quaternion& _orient, vector3d& _momentum, quaternion& _angmom)
	{
		center = _center; force = _force; torque = _torque; orient = _orient; momentum = _momentum, angmom = _angmom;
	}
	inline void reposition(vector3d &_center, quaternion &_orient)
	{
			center += _center;
			orient = _orient * orient;
	}
	void printtocout()
	{
		cout << "printing state...\n";
		center.printtocout();
		orient.printtocout();
	}
};

//#include "vertexInteraction.h"

struct Mass
{
	float mass = 0; char type = 'Y'; bool isGene = 0;
	vector3d initial, decomp, moment, oldPosition;
	inline Mass(){}
	inline Mass(vector3d &_init){initial = _init;}
	inline void print(ofstream &ofs)
	{
		ofs << mass << ' ' << type << '\n';
		initial.print(ofs);
	}
	inline void parse(ifstream &ifs)
	{
		ifs >> mass >> type;
		initial.read(ifs);
	}
};

struct cuboid
{
	vector3d dim, center, offset; //dimensions, xyz coordinates of center of rotation, offset of (0, 0, 0) corner from center of rotation when system is unrotated
	quaternion orient; //quaternion defining orientation
	double spacing, mass; //spacing between point masses during decomposition
	vector<vector3d> decomp; //decompose into point masses for dispersion calculations
	inline cuboid(){}
	inline cuboid(vector3d& _dim, vector3d& _center, quaternion& _orient) //assumes center of rotation is center of mass of uniform density
	{
		dim = _dim; center = _center; orient = _orient; spacing = 1;
		offset.assign(-dim.i/2, -dim.j/2, -dim.k/2);
		decomp.resize((dim.i + 1)*(dim.j + 1)*(dim.k + 1));
	}
	inline cuboid(vector3d& _dim, vector3d& _center, vector3d& _offset, quaternion& _orient, double &_spacing, double _mass) //for custom center of rotation
	{
		dim = _dim; center = _center; orient = _orient; offset = _offset;
		decomp.resize((dim.i + 1)*(dim.j + 1)*(dim.k + 1));
		spacing = _spacing; mass = _mass;
	}
	inline void rotate(quaternion& r)
	{
		orient = r*orient;
	}
	inline void decompose()
	{
		double i, j, k; int l = 0;
		for (i = 0; i <= dim.i; i+= spacing) for (j = 0; j <= dim.j; j += spacing) for (k = 0; k <= dim.k; k += spacing)
		{
			vector3d v3d(i, j, k); v3d += offset;
			v3d.rotate(orient); v3d += center;
			decomp[l] = v3d; l++;
		}
	}
	inline int num()
	{
		double i, j, k, tot = 0;
		for (i = 0; i <= dim.i; i += spacing) for (j = 0; j <= dim.j; j += spacing) for (k = 0; k <= dim.k; k += spacing)
		{
			tot++;
		}
		return tot;
	}
	inline void decompose(vector<Mass>& massList, int &idx)
	{
		double i, j, k;
		for (i = 0; i <= dim.i; i += spacing) for (j = 0; j <= dim.j; j += spacing) for (k = 0; k <= dim.k; k += spacing)
		{
			vector3d v3d(i, j, k); v3d += offset;
			v3d.rotate(orient); v3d += center;
			massList[idx].initial = v3d; idx++;
		}
	}
};

struct Charge
{
	float charge;
	vector3d initial, decomp, moment;
	inline Charge(){}
	inline void print(ofstream &ofs)
	{
		ofs << charge << '\n';
		initial.print(ofs);
	}
	inline void parse(ifstream &ifs)
	{
		ifs >> charge;
		initial.read(ifs);
	}
};

struct rigidBody
{
	vector3d center; //xyz coordinates of center

	//static variables
	vector<Mass> massList; int massListSize;
	vector<Charge> chargeList; int chargeListSize;
	bool isDummy, isGene = 0; //if isDummy do not apply integrator
	vector3d principalAxes[3];
	double MItensor[3][3], totalMass, I[3]; //MItensor = Moment of inertia tensor, I = moments of inertia along principal axes
	double Tstokes, Tbrownian, Rchi[3], Rstokes[3], Rbrownian[3], h2m; //drag
	double bondThreshold = -1;
	//positional variables
	//dynamic variables
	vector3d force, torque, momentum; //angmom = angular momentum
	quaternion orient, qtorque, angmom, chi; vector3d overAngmom; //quaternion defining orientation and torque
	double rotMatrix[3][3]; //Rotation matrix associated with orientation quaternion
	bool overrideAngmom = 0, scaffold = 0;
	//bond adjacency list
	vector<vector<int> > adjList;
	//thread
	vector3d *minDim, *maxDim; bool * updateCellList; //int chash;
	unordered_map<vector3d, int> um;
	vector<int> specList;
	void repointer(vector3d *_minDim, vector3d *_maxDim, bool * _ucl)
	{
		minDim = _minDim; maxDim = _maxDim; updateCellList = _ucl;
		//cout << "repointer: ";
		//minDim->printtocout(); maxDim->printtocout();
	}
	inline int findCell(int i)
	{
		vector3d v = massList[i].decomp - *minDim; v /= verletRadius;
		if (v.i < 0) v.i = 0;
		if (v.j < 0) v.j = 0;
		if (v.k < 0) v.k = 0;
		return (int)v.i + ((int)v.j<<12) + ((int)v.k<<24);
	}
	void fillUm()
	{
		for (int i = 0; i < massList.size(); i++)
		{
			//massList[i].initial.printtocout();
			um[massList[i].initial] = i;
		}
	}

	struct PQnode
	{
		double k; int a, b;
		inline PQnode(){}
		inline PQnode(double _k, int _a, int _b): k(_k), a(_a), b(_b){}
	};
	
	struct comparePQnode
	{
		bool operator()(const PQnode &a, const PQnode &b)
		{return a.k > b.k;}
	};

	void MST() //populate adjacency list
	{
		if (scaffold) return;
		adjList.resize(massListSize);
		//cout << bondThreshold << '\n';
		if (bondThreshold <= 0) //bondThreshold = 1.55;
		{
			priority_queue<PQnode, vector<PQnode>, comparePQnode> PQ;
			vector<bool> visited(massListSize);
			for (int i = 0; i < massListSize; i++)
			{
				if (massList[i].type != 'Y') visited[i] = 0; else visited[i] = 1;
			}
			for (int i = 0; i < massListSize; i++) if (!visited[i])
			{
				visited[i] = 1;
				for (int j = 0; j < massListSize; j++) if (!visited[j])
				{
					double currdist = dist(massList[i].initial, massList[j].initial);
					PQnode pqn(currdist, i, j);
					PQ.push(pqn);
				}
				break;
			}
			while (!PQ.empty())
			{
				PQnode shortest = PQ.top(); PQ.pop();
				if (!visited[shortest.b])
				{
					visited[shortest.b] = 1;
					adjList[shortest.a].push_back(shortest.b);
					for (int j = 0; j < massListSize; j++) if (!visited[j])
					{
						double currdist = dist(massList[shortest.b].initial, massList[j].initial);
						PQnode pqn(currdist, shortest.b, j);
						PQ.push(pqn);
					}
				}
			}
		}
		/*{
			int i, j, curr = 0, next;
			double currmin = intmax;
			vector<double> mindist(massListSize);
			vector<int> neighbour(massListSize); vector<bool> visited(massListSize);
			mindist[0] = 0; neighbour[0] = 0; visited[0] = 1;
			for (i = 1; i < massListSize; i++)
			{
				mindist[i] = intmax; neighbour[i] = 0;
				if (massList[i].type != 'Y') visited[i] = 0; else visited[i] = 1;
			}
			for (i = 1; i < massListSize; i++)
			{
				for (j = 1; j < neighbour.size(); j++)
				{
					if (!visited[j])
					{
						double currdist = dist(massList[curr].initial, massList[j].initial);
						//initial[curr].printtocout(); initial[j].printtocout();
						//cout << curr << ' ' << j << ' ' << currdist << '\n';
						if (currdist < mindist[j])
						{
							neighbour[j] = curr;
							mindist[j] = currdist;
							if (currdist < currmin)
							{
								currmin = currdist;
								next = j;
							}
						}
					}
				}
				currmin = intmax;
				for (j = 1; j < neighbour.size(); j++)
				{
					if (!visited[j])
					{
						if (mindist[j] < currmin)
						{
							currmin = mindist[j];
							next = j;
						}
					}
				}
				adjList[neighbour[next]].push_back(next);
				curr = next;
				visited[next] = 1;
			}
		}*/
		else
		{
			for (int i = 0; i < massListSize - 1; i++)
			{
				if (massList[i].type != 'Y')
				{
					for (int j = i + 1; j < massListSize; j++)
					{
						if (massList[j].type != 'Y')
						{
							double currdist = dist(massList[i].initial, massList[j].initial);
							//cout << "currdist: " << currdist << '\n';
							if (currdist < bondThreshold)
							{
								adjList[i].push_back(j);
							}
						}
					}
				}
			}
		}
	}
	void printAdjList()
	{
		for (int i = 0; i < adjList.size(); i++)
		{
			cout << i << ": ";
			for (int j = 0; j < adjList[i].size(); j++) cout << adjList[i][j] << ' ';
			cout << '\n';
		}
	}
	vector3d centerMass(vector3d &trueCM) //gives xyz coordinates of center of mass relative to center of rotation assuming uniform density
	{
		vector3d CM;
		for (int i = 0; i < massListSize; i++) CM += massList[i].mass * massList[i].initial;
		CM /= totalMass;
		trueCM = CM; trueCM.rotate(orient);
		return CM;
	}
	void findMItensor()
	{
		double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0, MIsphere = 0;
		for (int i = 0; i < massListSize; i++)
		{
			//cout << massListSize << ' ' << i << ' ' << massList[i].mass << '\n';
			if (massList[i].mass > 0)
			{
				MIsphere = scale*pow(massList[i].mass, 5.0/3.0);
				//cout << massList[i].initial.i << ' ' << massList[i].initial.j << ' ' << massList[i].initial.k << '\n';
				xx += MIsphere + massList[i].mass*(massList[i].initial.j*massList[i].initial.j + massList[i].initial.k*massList[i].initial.k);
				yy += MIsphere + massList[i].mass*(massList[i].initial.i*massList[i].initial.i + massList[i].initial.k*massList[i].initial.k);
				zz += MIsphere + massList[i].mass*(massList[i].initial.i*massList[i].initial.i + massList[i].initial.j*massList[i].initial.j);
				xy -= massList[i].mass*massList[i].initial.i*massList[i].initial.j;
				xz -= massList[i].mass*massList[i].initial.i*massList[i].initial.k;
				yz -= massList[i].mass*massList[i].initial.j*massList[i].initial.k;
				/*cout << xx << ' ' << xy << ' ' << xz << '\n';
				cout << xy << ' ' << yy << ' ' << yz << '\n';
				cout << xz << ' ' << yz << ' ' << zz << '\n';*/
			}
		}
		MItensor[0][0] = xx; MItensor[0][1] = xy; MItensor[0][2] = xz;
		MItensor[1][0] = xy; MItensor[1][1] = yy; MItensor[1][2] = yz;
		MItensor[2][0] = xz; MItensor[2][1] = yz; MItensor[2][2] = zz;
	}
	void findI()
	{
		//for (int i = 0; i < 3; i++) {for (int j = 0; j < 3; j++) cout << MItensor[i][j] << ' '; cout << '\n';}
		//cout << "center: "; center.printtocout();
		Matrix3d principalMatrix(3, 3);
		for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
		{
			principalMatrix(i, j) = MItensor[i][j];
		}
		EigenSolver<MatrixXd> es(principalMatrix);
		vector<Vector3cd> vv(3);
		for (int i = 0; i < 3; i++)
		{
			vv[i] = es.eigenvectors().col(i);
		}
		for (int i = 0; i < 3; i++)
		{
			principalAxes[i].i = vv[i][0].real();
			principalAxes[i].j = vv[i][1].real();
			principalAxes[i].k = vv[i][2].real();
		}
		//cout << "before norm\n";
		for (int i = 0; i < 3; i++)
		{
			//principalAxes[i].printtocout();
			principalAxes[i].norm();
		}
		//cout << "after norm\n";
		/*for (int i = 0; i < 3; i++)
		{
			principalAxes[i].printtocout();
		}*/
		for (int i = 0; i < 3; i++)
		{
			I[i] = es.eigenvalues()[i].real();
		}
		//cout << "MI: " << I[0] << ' ' << I[1] << ' ' << I[2] << '\n';
		//cout << principalAxes[0]*principalAxes[0] << ' ' << principalAxes[1]*principalAxes[1] << ' ' << principalAxes[2]*principalAxes[2] << '\n';
		//cout << principalAxes[0]*principalAxes[1] << ' ' << principalAxes[0]*principalAxes[2] << ' ' << principalAxes[2]*principalAxes[1] << '\n';
		//cout << "Total Mass: " << totalMass << '\n';
	}
	void findDragConstants(double mu, double h)
	{
		h2m = h/(2*totalMass); int i;
		for (i = 0; i < 3; i++) Rchi[i] = 1/(4*I[i]);
		double radius = 0.5*cbrt(totalMass);
		double Tgamma = 6*PI*mu*radius/totalMass;
		Tstokes = exp(-Tgamma*h); Tbrownian = sqrt(totalMass*kT*(1 - exp(-2*Tgamma*h)))*sqrt3;
		double Rgamma = mu*8*PI*radius*radius*radius*(1/I[0] + 1/I[1] + 1/I[2]), M = 0;
		if (overrideAngmom)
		{
			Tgamma = 0; Rgamma = 0;
		}
		for (i = 0; i < 3; i++) M += 1/I[i];
		M = 4/M;
		for (i = 0; i < 3; i++) Rstokes[i] = exp(-Rgamma*M*h/(4*I[i]));
		for (i = 0; i < 3; i++)
		{
			Rbrownian[i] = sqrt(4*kT*I[i]*(1 - exp(-M*Rgamma*h/(2*I[i]))))*sqrt3; 
		}
		//cout << Tgamma << ' ' << Rgamma << '\n';
		//for (int i = 0; i < 3; i++) cout << Rbrownian[i] << ' ' << Rstokes[i] << '\n';
	}
	double momentInertia(vector3d axis) //Finds moment of inertia about an arbitrary axis.
	{
		axis.norm(); return axis * (MItensor * axis);
	}
	void convertTorque()
	{
		quaternion qtemp(0, torque.i, torque.j, torque.k);
		vector<vector<double> > Sorient;
		Sfier(orient, Sorient);
		qtorque = Sorient * qtemp; qtorque *= 2;
	}
	inline void rotate(quaternion& r)
	{
		orient = r*orient;
	}
	inline void decompose()
	{
		QtoRM(orient, rotMatrix);
		for (int i = 0; i < massListSize; i++)
		{
			massList[i].decomp = matrixMult3d(rotMatrix, massList[i].initial); massList[i].decomp += center;
			double disp = dist(massList[i].decomp, massList[i].oldPosition);
			if (disp > verletSkin) *updateCellList = 1;
		}
		for (int i = 0; i < chargeListSize; i++)
		{
			chargeList[i].decomp = matrixMult3d(rotMatrix, chargeList[i].initial); chargeList[i].decomp += center;
		}
	}
	inline void minDecompose()
	{
		QtoRM(orient, rotMatrix);
		for (int i = 0; i < massListSize; i++)
		{
			massList[i].decomp = matrixMult3d(rotMatrix, massList[i].initial); massList[i].decomp += center;
			if (massList[i].decomp.i < minDim->i) minDim->i = massList[i].decomp.i;
			if (massList[i].decomp.j < minDim->j) minDim->j = massList[i].decomp.j;
			if (massList[i].decomp.k < minDim->k) minDim->k = massList[i].decomp.k;
			if (massList[i].decomp.i > maxDim->i) maxDim->i = massList[i].decomp.i; 
			if (massList[i].decomp.j > maxDim->j) maxDim->j = massList[i].decomp.j;
			if (massList[i].decomp.k > maxDim->k) maxDim->k = massList[i].decomp.k;
			//massList[i].cell = findCell(i);
		}
		for (int i = 0; i < chargeListSize; i++)
		{
			chargeList[i].decomp = matrixMult3d(rotMatrix, chargeList[i].initial); chargeList[i].decomp += center;
		}
	}
	void printtocout()
	{
		cout << "decomp:\n";
		for (int i = 0; i < massListSize; i++) massList[i].decomp.printtocout();
		cout << "charges:\n";
		for (int i = 0; i < chargeListSize; i++) chargeList[i].decomp.printtocout();
	}
	void recenter()
	{
		vector3d trueCM, CM = centerMass(trueCM), xaxis(1, 0, 0), yaxis(0, 1, 0); center += trueCM; //CM.printtocout();
		for (int i = 0; i < massListSize; i++)
		{
			massList[i].initial -= CM;
		}
		for (int i = 0; i < chargeListSize; i++)
		{
			chargeList[i].initial -= CM;
		}
		findMItensor(); findI();
		vector3d xcross = cross(xaxis, principalAxes[0]); xcross.norm(); //cout << "xcross: "; xcross.printtocout();
		double xdot = principalAxes[0].i; if (xdot > 1) xdot = 1; if (xdot < -1) xdot = -1; double theta = acos(xdot); //cout << "theta: " << theta << '\n';
		quaternion qx = convert(xcross, theta); yaxis.rotate(qx); //cout << "yaxis: "; yaxis.printtocout();
		vector3d ycross = cross(yaxis, principalAxes[1]); ycross.norm(); //cout << "ycross: "; ycross.printtocout();
		double ydot = yaxis * principalAxes[1]; if (ydot > 1) ydot = 1; if (ydot < -1) ydot = -1; theta = acos(ydot); //cout << "theta: " << ydot << ' ' << theta << '\n';
		quaternion qy = convert(ycross, theta);
		quaternion qf = qy * qx; qf.norm();
		quaternion antiqf = qf; antiqf.conjugate();
		for (int i = 0; i < massListSize; i++)
		{
			massList[i].initial.rotate(antiqf);
		}
		for (int i = 0; i < chargeListSize; i++)
		{
			chargeList[i].initial.rotate(antiqf);
		}
		orient = orient * qf;
	}
	rigidBody()
	{
		force.assign(0, 0, 0);
		torque.assign(0, 0, 0);
		momentum.assign(0, 0, 0);
		qtorque.assign(0, 0, 0, 0);
		angmom.assign(0, 0, 0, 0);
		isDummy = 0;
	}
	rigidBody(vector<cuboid> &components, quaternion& _orient, vector3d& _center, vector<char> &ptype, vector3d *_minDim, vector3d *_maxDim, bool *_updateCellList) 
	: minDim(_minDim), maxDim(_maxDim), updateCellList(_updateCellList)
	{
		isDummy = 0;
		orient = _orient; center = _center;
		massListSize = 0;
		for (int i = 0; i < components.size(); i++)
		{
			int t = components[i].num();
			massListSize += t;
		}
		massList.resize(massListSize);
		int k = 0;
		for (int i = 0; i < components.size(); i++)
		{
			int idx = k;
			components[i].decompose(massList, idx);
			for (int j = k; j < idx; j++)
			{
				massList[j].mass = components[i].mass;
				massList[j].type = ptype[i];
			}
			k = idx;
		}
		totalMass = 0; for (int i = 0; i < massListSize; i++) totalMass += massList[i].mass;
		//for (int i = 0; i < 3; i++) {freezeRot[i] = 1; freezeCart[i] = 1;}
		force.assign(0, 0, 0);
		torque.assign(0, 0, 0);
		momentum.assign(0, 0, 0);
		qtorque.assign(0, 0, 0, 0);
		angmom.assign(0, 0, 0, 0);
		overrideAngmom = 0;
	}
	void resetMoments()
	{
		for (int i = 0; i < massListSize; i++) massList[i].moment.assign(0, 0, 0);
		for (int i = 0; i < chargeListSize; i++) chargeList[i].moment.assign(0, 0, 0);
	}
	inline void calcTorque()
	{
		for (int i = 0; i < massListSize; i++)
		{
			vector3d v = massList[i].decomp - center;
			torque += cross(v, massList[i].moment);
		}
		for (int i = 0; i < chargeListSize; i++)
		{
			vector3d v = chargeList[i].decomp - center;
			torque += cross(v, chargeList[i].moment);
		}
		quaternion conj = conjugate(orient);
		torque.rotate(conj);
	}
	inline rigidBody& operator = (const state& st)
	{
		this->center = st.center; this->momentum = st.momentum; this->angmom = st.angmom; this->force = st.force; this->torque = st.torque; this->orient = st.orient;
		return *this;
	}
//	inline void purgeFrozen() {force *= freezeCart; torque *= freezeRot;};
	void print(ofstream &ofs)
	{
		ofs << "decomp:\n";
		for (int i = 0; i < massListSize; i++) massList[i].decomp.print(ofs);
		cout << "charges:\n";
		for (int i = 0; i < chargeListSize; i++) chargeList[i].decomp.print(ofs);
	}
	void printInitial(ofstream &ofs)
	{
		center.print(ofs);
		orient.print(ofs);
		ofs << massListSize << '\n';
		for (int i = 0; i < massListSize; i++) massList[i].initial.print(ofs);
		ofs << chargeListSize << '\n';
		for (int i = 0; i < chargeListSize; i++) {ofs << chargeList[i].charge << ' '; chargeList[i].initial.print(ofs);}
	}
	inline void halt()
	{
		momentum.assign(0, 0, 0);
		angmom.assign(0, 0, 0, 0);
	}
};

struct baseConfig {double TKE, RKE, V, expenditure;};

struct config : baseConfig
{
	vector<state> vs;
	int RBlistSize; rigidBody * RBlist;
	config(){}
	config(rigidBody *_RBlist, int _RBlistSize){RBlist = _RBlist; RBlistSize = _RBlistSize;}
	void init(rigidBody *_RBlist, int _RBlistSize){RBlist = _RBlist; RBlistSize = _RBlistSize;}
	void overwrite()
	{
		if (vs.size() == RBlistSize)
		{
			for (int i = 0; i < RBlistSize; i++) RBlist[i] = vs[i];
		}
	}
	void reposition(vector3d &center, quaternion &orient)
	{
		for (int i = 0; i < vs.size(); i++) vs[i].reposition(center, orient);
	}
	void partialOverwrite(int offset, int run)
	{
		for (int i = 0; i < run; i++)
		{
			//vs[i].printtocout();
			RBlist[offset + i] = vs[i];
		}
	}
	void assignState(rigidBody& rb, state &s)
	{
		s.assign(rb.center, rb.force, rb.torque, rb.orient, rb.momentum, rb.angmom);
	}
	void assign()
	{
		vs.resize(RBlistSize);
		for (int i = 0; i < RBlistSize; i++) assignState(RBlist[i], vs[i]);
	}
	inline void overwriteEnergy(baseConfig &ec)
	{
		TKE = ec.TKE; RKE = ec.RKE, V = ec.V; expenditure = ec.expenditure;
	}
	void print(ofstream &ofs)
	{
		for (int i = 0; i < RBlistSize; i++)
		{
			vs[i].center.print(ofs);
			vs[i].momentum.print(ofs);
			vs[i].orient.print(ofs);
			vs[i].angmom.print(ofs);
		}
	}
	void read(ifstream &ifs)
	{
		vs.resize(RBlistSize);
		//cout << "config: " << RBlistSize << '\n';
		for (int i = 0; i < RBlistSize; i++)
		{
		//	cout << i << '\n';
			vs[i].center.read(ifs);
			vs[i].momentum.read(ifs);
			vs[i].orient.read(ifs);
			vs[i].angmom.read(ifs);
			//vs[i].printtocout();
		}
	}
};

struct energyConfig : baseConfig
{
	energyConfig(){TKE = 0; RKE = 0; V = 0; expenditure = 0;}
	inline void add(baseConfig & c)
	{
		TKE += c.TKE; RKE += c.RKE; V += c.V; expenditure += c.expenditure;
	}
	void divide(double d)
	{
		TKE /= d; RKE /= d; V /= d; expenditure /= d;
	}
};
