double Trot, Ttrans, V, totV, totTrot, totTtrans;

inline double rotationalKE(rigidBody &rb) //rotational KE calculated using formula from paper
{
	int i, j;
	double M[4][4], aa[4], ai[4], aj[4], ak[4], ii[4], ij[4], ik[4], jj[4], jk[4], kk[4], u;
	aa[0] = rb.orient.a*rb.orient.a; ai[0] = rb.orient.a*rb.orient.i; aj[0] = rb.orient.a*rb.orient.j; ak[0] = rb.orient.a*rb.orient.k;
	ii[0] = rb.orient.i*rb.orient.i; ij[0] = rb.orient.i*rb.orient.j; ik[0] = rb.orient.i*rb.orient.k;
	jj[0] = rb.orient.j*rb.orient.j; jk[0] = rb.orient.j*rb.orient.k;
	kk[0] = rb.orient.k*rb.orient.k;
	for (i = 1; i < 4; i++)
	{
		u = rb.I[i - 1];
		aa[i] = aa[0]/u; ai[i] = ai[0]/u; aj[i] = aj[0]/u; ak[i] = ak[0]/u;
		ii[i] = ii[0]/u; ij[i] = ij[0]/u; ik[i] = ik[0]/u;
		jj[i] = jj[0]/u; jk[i] = jk[0]/u;
		kk[i] = kk[0]/u;
	}
	M[0][0] = ii[1] + jj[2] + kk[3]; M[0][1] = -ai[1] + jk[2] - jk[3]; M[0][2] = -ik[1] - aj[2] + ik[3]; M[0][3] =  ij[1] - ij[2] - ak[3];
	M[1][1] = aa[1] + kk[2] + jj[3]; M[1][2] = ak[1] - ak[2] - ij[3]; M[1][3] =  -aj[1] - ik[2] + aj[3];
	M[2][2] = kk[1] + aa[2] + ii[3]; M[2][3] =  -jk[1] + ai[2] - ai[3];
	M[3][3] = jj[1] + ii[2] + aa[3];
	for (i = 0; i < 4; i++) for (j = 0; j < i; j++) M[i][j] = M[j][i];
	quaternion q = M * rb.angmom;
	return dot(rb.angmom, q)/8;
}

inline double translationalKE(rigidBody &rb)
{
	return 0.5*rb.momentum.square()/rb.totalMass;
}

inline int randomDistribution() {return rand()%6;}

inline void psi(rigidBody &rb, double h, int l, quaternion qi, quaternion pi, quaternion &qf, quaternion &pf)
{
	//Function psi(q, pi) -> (Q, PI) from the paper
	double chih; //chi * timestep
	quaternion sq = Smult(l, qi); quaternion sp = Smult(l, pi); //Smult returns the result of Sl * q, where q is a quaternion. It is defined in line 355 of quaternion.h.
	chih = dot(pi, sq) * rb.Rchi[l - 1] * h; //Rchi is defined in line 106 of the rigidBody.h file (under the function findDragConstants)
	double sint = sin(chih), cost = cos(chih);
	//double sint = chih, cost = 1 - chih*chih/2;
	sq *= sint; sp *= sint; qi *= cost; pi *= cost;
	qf = qi + sq;
	pf = pi + sp;
}

inline void psit(rigidBody &rb, double h, const quaternion &qi, const quaternion &pi, quaternion &qf, quaternion &pf)
{
	//Composite function psi
	psi(rb, h/2, 3, qi, pi, qf, pf);
	psi(rb, h/2, 2, qf, pf, qf, pf);
	psi(rb, h, 1, qf, pf, qf, pf);
	psi(rb, h/2, 2, qf, pf, qf, pf);
	psi(rb, h/2, 3, qf, pf, qf, pf);
}

inline vector3d Tnoise(double Tbrownian)
{
	//Translational noise
	vector3d out;
	int k = randomDistribution();
	if (k > 3)
	{
		if (k == 4) out.i = Tbrownian;
		else out.i = -Tbrownian;
	}
	k = randomDistribution();
	if (k > 3)
	{
		if (k == 4) out.j = Tbrownian;
		else out.j = -Tbrownian;
	}
	k = randomDistribution();
	if (k > 3)
	{
		if (k == 4) out.k = Tbrownian;
		else out.k = -Tbrownian;
	}
	return out;
}

inline quaternion piDampingNoise(rigidBody &rb, quaternion &q, quaternion &pi)
{
	//Rotational damping, Rstokes is defined in line 113 of rigidBody.h
	quaternion damping, noise(0, 0, 0, 0); int i, j, k;
	double M[4][4]; for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) M[i][j] = 0;
	for (i = 0; i < 3; i++)
	{
		quaternion t = Smult(i + 1, q);
		double arr[4];
		//Finding the outer product of Sl*q
		arr[0] = t.a * rb.Rstokes[i]; arr[1] = t.i * rb.Rstokes[i]; arr[2] = t.j * rb.Rstokes[i]; arr[3] = t.k * rb.Rstokes[i];
		M[0][0] += arr[0] * t.a; M[0][1] += arr[0] * t.i; M[0][2] += arr[0] * t.j; M[0][3] += arr[0] * t.k;
		M[1][1] += arr[1] * t.i; M[1][2] += arr[1] * t.j; M[1][3] += arr[1] * t.k;
		M[2][2] += arr[2] * t.j; M[2][3] += arr[2] * t.k;
		M[3][3] += arr[3] * t.k;
		k = randomDistribution();
		if (k > 3)
		{
			if (k == 4) t *= rb.Rbrownian[i];
			else t *= -rb.Rbrownian[i];
			noise += t;
		}
	}
	for (j = 0; j < 4; j++) for (k = 0; k < j; k++) M[j][k] = M[k][j];
	damping = M * pi; damping += noise;
	return damping;
}

inline quaternion piDamping(rigidBody &rb, quaternion &q, quaternion &pi)
{
	//Rotational damping, Rstokes is defined in line 113 of rigidBody.h
	quaternion damping, noise(0, 0, 0, 0); int i, j, k;
	double M[4][4]; for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) M[i][j] = 0;
	for (i = 0; i < 3; i++)
	{
		quaternion t = Smult(i + 1, q);
		double arr[4];
		//Finding the outer product of Sl*q
		arr[0] = t.a * rb.Rstokes[i]; arr[1] = t.i * rb.Rstokes[i]; arr[2] = t.j * rb.Rstokes[i]; arr[3] = t.k * rb.Rstokes[i];
		M[0][0] += arr[0] * t.a; M[0][1] += arr[0] * t.i; M[0][2] += arr[0] * t.j; M[0][3] += arr[0] * t.k;
		M[1][1] += arr[1] * t.i; M[1][2] += arr[1] * t.j; M[1][3] += arr[1] * t.k;
		M[2][2] += arr[2] * t.j; M[2][3] += arr[2] * t.k;
		M[3][3] += arr[3] * t.k;
	}
	for (j = 0; j < 4; j++) for (k = 0; k < j; k++) M[j][k] = M[k][j];
	damping = M * pi;
	return damping;
}

inline quaternion piDampingPrint(rigidBody &rb, quaternion &q, quaternion &pi)
{
	//Rotational damping, Rstokes is defined in line 113 of rigidBody.h
	quaternion damping, noise(0, 0, 0, 0); int i, j, k;
	double M[4][4]; for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) M[i][j] = 0;
	for (i = 0; i < 3; i++)
	{
		quaternion t = Smult(i + 1, q);
		double arr[4];
		//Finding the outer product of Sl*q
		cout << "t: "; t.printtocout();
		cout << "Stokes: " << rb.Rstokes[i] << '\n';
		arr[0] = t.a * rb.Rstokes[i]; arr[1] = t.i * rb.Rstokes[i]; arr[2] = t.j * rb.Rstokes[i]; arr[3] = t.k * rb.Rstokes[i];
		M[0][0] += arr[0] * t.a; M[0][1] += arr[0] * t.i; M[0][2] += arr[0] * t.j; M[0][3] += arr[0] * t.k;
		M[1][1] += arr[1] * t.i; M[1][2] += arr[1] * t.j; M[1][3] += arr[1] * t.k;
		M[2][2] += arr[2] * t.j; M[2][3] += arr[2] * t.k;
		M[3][3] += arr[3] * t.k;
	}
	for (j = 0; j < 4; j++) for (k = 0; k < j; k++) M[j][k] = M[k][j];
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++) cout << M[i][j] << ' ' ;
		cout << '\n';
	}
	damping = M * pi;
	return damping;
}

inline void SDmin(int timesteps)
{
	resetExpend();
	isMin = 1; isEquil = 0, isMD = 0;
	double t2;
	double h = 5e-4, h2 = h/2;
	int f = timesteps/10; if (f == 0) f++;
	SDtrajectory.resize(timesteps/f + 1, config(RBlist, RBlistSize));
	int i, j; vector3d x1;
	for (i = 0; i < RBlistSize; i++) RBlist[i].findDragConstants(1.5, h);
	quaternion q1, pi2, pi3;
	decompAll();
	V = totalForce(0, 0);
	cout << "time: " << 0 << " potential: " << V << '\n';
	printContributions(0);
	for (i = 0; i < RBlistSize; i++)
	{
		if (!RBlist[i].isDummy || !inactivateDum)
		{
			rigidBody &rb = RBlist[i];
			rb.convertTorque();
			rb.momentum += rb.force*h2;
			rb.angmom = rb.angmom + rb.qtorque*h2;
		}
	}
	for (j = 1; j <= timesteps; j++)
	{
		int w = j%f;
		currStep = j;
		t2 = clock();
		for (i = 0; i < RBlistSize; i++)
		{
			if (!RBlist[i].isDummy || !inactivateDum)
			{
				rigidBody &rb = RBlist[i];
				x1 = rb.center + (rb.h2m*rb.momentum);
				psit(rb, h2, rb.orient, rb.angmom, q1, pi2);
				rb.momentum = rb.momentum*rb.Tstokes;
				pi3 = piDamping(rb, q1, pi2);
				rb.center = x1 + (rb.h2m*rb.momentum);
				psit(rb, h2, q1, pi3, rb.orient, rb.angmom);
			}
		}
		decompAll();
		integratorTime += clock() - t2;
		t2 = clock();
		V = totalForce(0, 0);
		interactionTime += clock() - t2;
		t2 = clock();
		Trot = 0; Ttrans = 0;
		for (i = 0; i < RBlistSize; i++)
		{
			if (!RBlist[i].isDummy || !inactivateDum)
			{
				rigidBody &rb = RBlist[i]; rb.convertTorque();
				rb.momentum += rb.force*h;
				rb.angmom = rb.angmom + rb.qtorque*h;
				Trot += rotationalKE(rb);
				Ttrans += translationalKE(rb);
			}
		}
		liTime += clock() - t2;
		if (w == 0) 
		{
			SDtrajectory[j/f].assign(); SDtrajectory[j/f].V = V; SDtrajectory[j/f].TKE = Ttrans; SDtrajectory[j/f].RKE = Trot;
			if (writetocout)
			{
				cout << "time: " << j*h << " potential: " << V << " Rotational KE: " << Trot << " Translational KE: " << Ttrans << " TE: " << Trot + Ttrans + V << '\n';
				printContributions(0);
			}
		}
	}
}

inline void equilibriate(int timesteps, int f)
{
	resetExpend();
	isMin = 0; isEquil = 1, isMD = 0;
	maxTimeStep = timesteps;
	totV = 0; totTrot = 0; totTtrans = 0;
	int i, j; vector3d x1;
	double ratio = 0, t2;
	for (i = 0; i < RBlistSize; i++) RBlist[i].findDragConstants(mu, h);
	quaternion q1, pi2, pi3;
	decompAll(); V = totalForce(0, 0);
	for (i = 0; i < RBlistSize; i++)
	{
		if (!RBlist[i].isDummy || !inactivateDum)
		{
			rigidBody &rb = RBlist[i];
			rb.convertTorque();
			rb.momentum += rb.force*h2;
			rb.angmom = rb.angmom + rb.qtorque*h2;
		}
	}
	for (j = 0; j <= timesteps; j++)
	{
		int w = j%f; if (w == 0) printStep = 1; else printStep = 0;
		currStep = j;
		t2 = clock();
		for (i = 0; i < RBlistSize; i++)
		{
			if (!RBlist[i].isDummy || !inactivateDum)
			{
				rigidBody &rb = RBlist[i];
				x1 = rb.center + (rb.h2m*rb.momentum);
				psit(rb, h2, rb.orient, rb.angmom, q1, pi2);
				if (noise)
				{
					rb.momentum = rb.momentum*rb.Tstokes + Tnoise(rb.Tbrownian);
					pi3 = piDampingNoise(rb, q1, pi2);
				}
				else
				{
					rb.momentum = rb.momentum*rb.Tstokes;
					pi3 = piDamping(rb, q1, pi2);
				}
				rb.center = x1 + (rb.h2m*rb.momentum);
				psit(rb, h2, q1, pi3, rb.orient, rb.angmom);
			}
		}
		decompAll();
		integratorTime += clock() - t2;
		t2 = clock();
		V = totalForce(j*ratio, printStep);
		interactionTime += clock() - t2;
		t2 = clock();
		Trot = 0; Ttrans = 0;
		for (i = 0; i < RBlistSize; i++)
		{
			if (!RBlist[i].isDummy || !inactivateDum)
			{
				rigidBody &rb = RBlist[i]; rb.convertTorque();
				rb.momentum += rb.force*h;
				rb.angmom = rb.angmom + rb.qtorque*h;
				Trot += rotationalKE(rb);
				Ttrans += translationalKE(rb);
			}
		}
		liTime += clock() - t2;
		if (w == 0) 
		{
			if (writetocout)
			{
				cout << "time: " << j*h << " potential: " << V << " Rotational KE: " << Trot << " Translational KE: " << Ttrans << " TE: " << Trot + Ttrans + V << '\n';
				printContributions(0);
			}
		}
		totV += V; totTrot += Trot; totTtrans += Ttrans;
	}
	cout << totV/timesteps << ' ' << totTrot/timesteps << ' ' << totTtrans/timesteps << '\n';
}

inline void roteq(int timesteps, int f)
{
	resetExpend();
	isMin = 0; isEquil = 1, isMD = 0;
	maxTimeStep = timesteps;
	totV = 0; totTrot = 0; totTtrans = 0;
	double stepV = 0, stepRot = 0, stepTrans = 0;
	eqtrajectory.resize(timesteps/f + 1, config(RBlist, RBlistSize));
	int i, j; vector3d x1;
	double ratio = 0, t2;
	for (i = 0; i < RBlistSize; i++) RBlist[i].findDragConstants(mu, h);
	quaternion q1, pi2, pi3;
	decompAll(); V = totalForce(0, 0);
	for (i = 0; i < RBlistSize; i++)
	{
		if (!RBlist[i].isDummy || !inactivateDum)
		{
			rigidBody &rb = RBlist[i];
			rb.convertTorque();
			rb.momentum += rb.force*h2;
			rb.angmom = rb.angmom + rb.qtorque*h2;
		}
	}
	for (j = 0; j <= timesteps; j++)
	{
		int w = j%f; if (w == 0) printStep = 1; else printStep = 0;
		currStep = j;
		t2 = clock();
		for (i = 0; i < RBlistSize; i++)
		{
			if (!RBlist[i].isDummy || !inactivateDum)
			{
				rigidBody &rb = RBlist[i];
				x1 = rb.center + (rb.h2m*rb.momentum);
				psit(rb, h2, rb.orient, rb.angmom, q1, pi2);
				if (noise)
				{
					rb.momentum = rb.momentum*rb.Tstokes + Tnoise(rb.Tbrownian);
					pi3 = piDampingNoise(rb, q1, pi2);
				}
				else
				{
					rb.momentum = rb.momentum*rb.Tstokes;
					pi3 = piDamping(rb, q1, pi2);
				}
				rb.center = x1 + (rb.h2m*rb.momentum);
				psit(rb, h2, q1, pi3, rb.orient, rb.angmom);
			}
		}
		decompAll();
		integratorTime += clock() - t2;
		t2 = clock();
		rotateED();
		V = totalForce(j*ratio, printStep);
		interactionTime += clock() - t2;
		t2 = clock();
		Trot = 0; Ttrans = 0;
		for (i = 0; i < RBlistSize; i++)
		{
			if (!RBlist[i].isDummy || !inactivateDum)
			{
				rigidBody &rb = RBlist[i]; rb.convertTorque();
				rb.momentum += rb.force*h;
				rb.angmom = rb.angmom + rb.qtorque*h;
				Trot += rotationalKE(rb);
				Ttrans += translationalKE(rb);
			}
		}
		liTime += clock() - t2;
		stepV += V; stepRot += Trot; stepTrans += Ttrans;
		if (w == 0) 
		{
			if (j == 0)
			{
				eqtrajectory[j/f].assign(); eqtrajectory[j/f].V = stepV; eqtrajectory[j/f].TKE = stepTrans; eqtrajectory[j/f].RKE = stepRot;
			}
			else
			{
				eqtrajectory[j/f].assign(); eqtrajectory[j/f].V = stepV/f; eqtrajectory[j/f].TKE = stepTrans/f; eqtrajectory[j/f].RKE = stepRot/f;
			}
			eqtrajectory[j/f].expenditure = 0;
			stepV = 0; stepRot = 0; stepTrans = 0;
			if (writetocout)
			{
				cout << "time: " << j*h << " potential: " << V << " Rotational KE: " << Trot << " Translational KE: " << Ttrans << " TE: " << Trot + Ttrans + V << '\n';
				printContributions(0);
			}
		}
		totV += V; totTrot += Trot; totTtrans += Ttrans;
	}
	cout << totV/timesteps << ' ' << totTrot/timesteps << ' ' << totTtrans/timesteps << '\n';
}

inline void areq(int timesteps, int f)
{
	resetExpend();
	externalDipole = initialDipole; prevDipole = initialDipole;
	offsetDipole = initOffDipole; prevOffDipole = initOffDipole;
	isMin = 0; isEquil = 0, isMD = 0;
	maxTimeStep = timesteps;
	totV = 0; totTrot = 0; totTtrans = 0;
	int i, j; vector3d x1;
	double ratio = 0, t2;
	for (i = 0; i < RBlistSize; i++) RBlist[i].findDragConstants(mu, h);
	quaternion q1, pi2, pi3;
	decompAll(); V = totalForce(0, 0);
	for (i = 0; i < RBlistSize; i++)
	{
		if (!RBlist[i].isDummy || !inactivateDum)
		{
			rigidBody &rb = RBlist[i];
			rb.convertTorque();
			rb.momentum += rb.force*h2;
			rb.angmom = rb.angmom + rb.qtorque*h2;
		}
	}
	for (j = 0; j <= timesteps; j++)
	{
		int w = j%f; if (w == 0) printStep = 1; else printStep = 0;
		currStep = j;
		t2 = clock();
		for (i = 0; i < RBlistSize; i++)
		{
			if (!RBlist[i].isDummy || !inactivateDum)
			{
				rigidBody &rb = RBlist[i];
				x1 = rb.center + (rb.h2m*rb.momentum);
				psit(rb, h2, rb.orient, rb.angmom, q1, pi2);
				if (noise)
				{
					rb.momentum = rb.momentum*rb.Tstokes + Tnoise(rb.Tbrownian);
					pi3 = piDampingNoise(rb, q1, pi2);
				}
				else
				{
					rb.momentum = rb.momentum*rb.Tstokes;
					pi3 = piDamping(rb, q1, pi2);
				}
				rb.center = x1 + (rb.h2m*rb.momentum);
				psit(rb, h2, q1, pi3, rb.orient, rb.angmom);
			}
		}
		decompAll();
		integratorTime += clock() - t2;
		t2 = clock();
		V = totalForce(j*ratio, printStep);
		interactionTime += clock() - t2;
		t2 = clock();
		Trot = 0; Ttrans = 0;
		for (i = 0; i < RBlistSize; i++)
		{
			if (!RBlist[i].isDummy || !inactivateDum)
			{
				rigidBody &rb = RBlist[i]; rb.convertTorque();
				rb.momentum += rb.force*h;
				rb.angmom = rb.angmom + rb.qtorque*h;
				Trot += rotationalKE(rb);
				Ttrans += translationalKE(rb);
			}
		}
		liTime += clock() - t2;
		if (w == 0) 
		{
			if (writetocout)
			{
				cout << "time: " << j*h << " potential: " << V << " Rotational KE: " << Trot << " Translational KE: " << Ttrans << " TE: " << Trot + Ttrans + V << '\n';
				printContributions(0);
			}
		}
		totV += V; totTrot += Trot; totTtrans += Ttrans;
	}
	cout << totV/timesteps << ' ' << totTrot/timesteps << ' ' << totTtrans/timesteps << '\n';
}

inline void langevinIntegrator(int timesteps, int f, bool isChain)
{
	resetExpend();
	externalDipole = initialDipole; prevDipole = initialDipole;
	offsetDipole = initOffDipole; prevOffDipole = initOffDipole;
	isMin = 0; isEquil = 0, isMD = 1;
	maxTimeStep = timesteps; expenditure = 0;
	totV = 0; totTrot = 0; totTtrans = 0;
	double stepV = 0, stepRot = 0, stepTrans = 0;
	MDtrajectory.resize(timesteps/f + 1, config(RBlist, RBlistSize));
	int i, j; vector3d x1;
	double ratio = 1, t2;
	double stepTKE, stepRKE;
	for (i = 0; i < RBlistSize; i++)
	{
		rigidBody &rb = RBlist[i];
		rb.findDragConstants(mu, h);
		/*if (rb.overrideAngmom)
		{
			rb.overAngmom *= 2*PI;
			rb.overAngmom = vmult(rb.overAngmom, rb.I);
			rb.angmom = omegaConvert(rb.overAngmom, rb.orient);
		}*/
	}
	//double initRKE = rotationalKE(RBlist[RBlistSize - 1]);
	quaternion q1, pi2, pi3;
	decompAll(); V = totalForce(0, 0);
	double prevExpenditure = 0;
	for (i = 0; i < RBlistSize; i++)
	{
		if (!RBlist[i].isDummy || !inactivateDum)
		{
			rigidBody &rb = RBlist[i];
			rb.convertTorque();
			rb.momentum += rb.force*h2;
			rb.angmom = rb.angmom + rb.qtorque*h2;
		}
	}
	for (j = 0; j <= timesteps; j++)
	{
		int w = j%f; if (w == 0) printStep = 1; else printStep = 0;
		currStep = j;
		t2 = clock();
		for (i = 0; i < RBlistSize; i++)
		{
			if (!RBlist[i].isDummy || !inactivateDum)
			{
				rigidBody &rb = RBlist[i];
				x1 = rb.center + (rb.h2m*rb.momentum);
				psit(rb, h2, rb.orient, rb.angmom, q1, pi2);
				if (noise)
				{
					rb.momentum = rb.momentum*rb.Tstokes + Tnoise(rb.Tbrownian);
					pi3 = piDampingNoise(rb, q1, pi2);
				}
				else
				{
					rb.momentum = rb.momentum*rb.Tstokes;
					if (rb.overrideAngmom) pi3 = pi2;
					else pi3 = piDamping(rb, q1, pi2);
				}
				rb.center = x1 + (rb.h2m*rb.momentum);
				psit(rb, h2, q1, pi3, rb.orient, rb.angmom);
			}
		}
		decompAll();
		integratorTime += clock() - t2;
		t2 = clock();
		if (dipoleOn) rotateED();
		V = totalForce(j*ratio, printStep);
		interactionTime += clock() - t2;
		t2 = clock();
		Trot = 0; Ttrans = 0;
		for (i = 0; i < RBlistSize; i++)
		{
			if (!RBlist[i].isDummy || !inactivateDum)
			{
				rigidBody &rb = RBlist[i]; rb.convertTorque();
				rb.momentum += rb.force*h;
				rb.angmom = rb.angmom + rb.qtorque*h;
				stepRKE = rotationalKE(rb);
				stepTKE = translationalKE(rb);
				Trot += stepRKE;
				Ttrans += stepTKE;
				sumRKE += stepRKE;
				sumTKE += stepTKE;
			}
		}
		populateMatrix(w, f);
		liTime += clock() - t2;
		if (!expSum)
		{
			stepV += V; stepRot += Trot; stepTrans += Ttrans;
		}
		else
		{
			stepV = V; stepRot = Trot; stepTrans = Ttrans;
		}
		if (w == 0) 
		{
			if (!expSum)
			{
				if (j == 0)
				{
					MDtrajectory[j/f].assign(); MDtrajectory[j/f].V = stepV; MDtrajectory[j/f].TKE = stepTrans; MDtrajectory[j/f].RKE = stepRot;
				}
				else
				{
					MDtrajectory[j/f].assign(); MDtrajectory[j/f].V = stepV/f; MDtrajectory[j/f].TKE = stepTrans/f; MDtrajectory[j/f].RKE = stepRot/f;
				}
			}
			else
			{
				MDtrajectory[j/f].assign(); MDtrajectory[j/f].V = stepV; MDtrajectory[j/f].TKE = stepTrans; MDtrajectory[j/f].RKE = stepRot;
			}
			MDtrajectory[j/f].expenditure = expenditure;
			stepV = 0; stepRot = 0; stepTrans = 0;
			if (writetocout)
			{
				cout << "-----\n";
				cout << name << " run " << runNo << " of " << maxRuns <<": " << (j + 0.0)/timesteps*100.0 << "%\n";
				//cout << "ExtDipole: "; externalDipole.printtocout();
				cout << "time: " << j*h << " potential: " << V << " Rotational KE: " << Trot << " Translational KE: " << Ttrans << " TE: " << Trot + Ttrans + V << '\n';
				cout << "expenditure: " << expenditure << '\n';
				printContributions(j*ratio);
			}
		}
		recordFinalObservables();
		if (isChain && j != 0 && j%(timesteps/chain) == 0)
		{
			if (expenditureList.size() == 0) expenditureList.push_back(expenditure);
			else expenditureList.push_back(expenditure - prevExpenditure);
			prevExpenditure = expenditure;
		}
		if (isChain && j != 0 && j%(timesteps/chain) == 0)
		{
			findFinalObservables();
		}
		totV += V; totTrot += Trot; totTtrans += Ttrans;
	}
	if (!isChain) findFinalObservables();
	cout << totV/timesteps << ' ' << totTrot/timesteps << ' ' << totTtrans/timesteps << '\n';
	//cout << initRKE << ' ' << rotationalKE(RBlist[RBlistSize - 1]) << '\n';
}
