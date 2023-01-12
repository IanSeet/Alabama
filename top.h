void printRigid(ofstream &ofs)
{
	ofs << RBlistSize << '\n';
	for (int i = 0; i < RBlistSize; i++)
	{
		rigidBody &rb = RBlist[i];
		ofs << rb.massListSize << '\n';
		for (int j = 0; j < rb.massListSize; j++)
		{
			rb.massList[j].print(ofs);
		}
		ofs << rb.specList.size() << '\n';
		for (int j = 0; j < rb.specList.size(); j++)
		{
			ofs << rb.specList[j] << '\n';
		}
		ofs << rb.chargeListSize << '\n';
		for (int j = 0; j < rb.chargeListSize; j++)
		{
			rb.chargeList[j].print(ofs);
		}
		for (int j = 0; j < 3; j++) rb.principalAxes[j].print(ofs);
		for (int j = 0; j < 3; j++) ofs << rb.I[j] << ' ';
		ofs << '\n';
		ofs << rb.isDummy << ' ' << rb.isGene << ' ' << rb.scaffold << '\n';
	}
}

void printConfig(ofstream &ofs)
{
	config finalConfig(RBlist, RBlistSize); finalConfig.assign(); finalConfig.print(ofs);
}

void printInteractions(ofstream &ofs)
{
	ofs << intListSize << '\n';
	for (int i = 0; i < intListSize; i++)
	{
		interactionList[i]->print(ofs);
	}
}

void printSuppressed(ofstream &ofs)
{
	ofs << suppressed.size() << '\n';
	for (int i = 0; i < suppressed.size(); i++) ofs << suppressed[i].first << ' ' << suppressed[i].second << '\n';
	ofs << '\n' << umbrellaed.size() << '\n';
	for (int i = 0; i < umbrellaed.size(); i++) ofs << umbrellaed[i].first << ' ' << umbrellaed[i].second << '\n';
}

void printSpecial(ofstream &ofs)
{
	ofs << topSpecList.size() << '\n';
	for (int i = 0; i < topSpecList.size(); i++) ofs << topSpecList[i].first << ' ' << topSpecList[i].second << '\n';
}

int parseRigid(ifstream &ifs, ifstream &ifsConfig, bool ind, int offset, vector3d &center, quaternion &orient)
{
	int x;
	ifs >> x;
	cout << "RBlistSize: " << x << '\n';
	if (ind) 
	{
		RBlistSize = x;
		RBlist = new rigidBody[RBlistSize];
	}
	for (int i = offset; i < x + offset; i++)
	{
		//cout << "i: " << i << '\n';
		rigidBody &rb = RBlist[i];
		rb.minDim = &minDim;
		rb.maxDim = &maxDim;
		rb.updateCellList = &updateCellList;
		ifs >> rb.massListSize;
		//cout << "massListSize: " << rb.massListSize << '\n';
		rb.massList.resize(rb.massListSize);
		for (int j = 0; j < rb.massListSize; j++)
		{
			rb.massList[j].parse(ifs);
		}
		int specListSize;
		ifs >> specListSize;
		rb.specList.resize(specListSize);
		//cout << "specListSize : " << rb.specList.size() << '\n';
		for (int j = 0; j < rb.specList.size(); j++)
		{
			ifs >> rb.specList[j];
			//cout << "specListNo: " << rb.specList[j] << '\n';
		}
		ifs >> rb.chargeListSize;
		//cout << "chargeListSize : " << rb.chargeListSize << '\n';
		rb.chargeList.resize(rb.chargeListSize);
		for (int j = 0; j < rb.chargeListSize; j++)
		{
			rb.chargeList[j].parse(ifs);
		}
		for (int j = 0; j < 3; j++) rb.principalAxes[j].read(ifs);
		for (int j = 0; j < 3; j++) ifs >> rb.I[j];
		ifs >> rb.isDummy >> rb.isGene >> rb.scaffold;
		//cout << "dummy: " << rb.isDummy << ' ' << rb.isGene << ' ' << rb.scaffold << '\n';
		rb.totalMass = 0; for (int i = 0; i < rb.massListSize; i++) rb.totalMass += rb.massList[i].mass;
		quaternion q(1, 0, 0, 0); rb.orient = q;
		//cout << "hi\n";
		rb.decompose();
		//rb.recenter();
		rb.MST();
		//cout << "rb done\n";
	}
	cout << "rb top parsed\n";
	config readConfig(RBlist, x);
	readConfig.read(ifsConfig);
	if (ind)
	{
		colour.resize(RBlistSize);
		readConfig.overwrite();
	}
	else
	{
		readConfig.reposition(center, orient);
		readConfig.partialOverwrite(offset, x);
	}
	return x;
}

void parseInteractions(ifstream &ifs, bool ind, int offsetRB, int offsetInt, vector3d &center, quaternion &rotate, double phase, double offPhase)
{
	int x;
	ifs >> x;
	if (ind)
	{
		intListSize = x;
		interactionList = new interaction*[intListSize];
	}
	cout << "intListSize: " << intListSize << '\n';
	for (int i = offsetInt; i < x + offsetInt; i++)
	{
		int type; ifs >> type;
		//cout << i << ' ' << type << '\n';
		if (type == 51)
		{
			interactionList[i] = new extDipole;
			extDipole * ptr = (extDipole*)interactionList[i];
			ptr->parse(ifs, RBlist, &externalDipole, &prevDipole, &expenditure, offsetRB, rotationAxis);
			if (phase != 0)
			{
				ptr->deflection += phase * PI/180;
				if (ptr->deflection < 0) ptr->deflection += 2*PI;
				ptr->q = convert(rotationAxis, ptr->deflection);
				ptr->qrev = convert(rotationAxis, -ptr->deflection);
				//ptr->q.printtocout();
			}
		}
		else if (type == 52)
		{
			interactionList[i] = new offDipole;
			offDipole * ptr = (offDipole*)interactionList[i];
			ptr->parse(ifs, RBlist, &offsetDipole, &prevOffDipole, &expenditure, offsetRB, rotationAxis);
			if (offPhase != 0)
			{
				ptr->deflection += offPhase * PI/180;
				if (ptr->deflection < 0) ptr->deflection += 2*PI;
				ptr->q = convert(rotationAxis, ptr->deflection);
				ptr->qrev = convert(rotationAxis, -ptr->deflection);
				//ptr->q.printtocout();
			}
		}
		else if (type == 53)
		{
			interactionList[i] = new staccatoDipole;
			staccatoDipole * ptr = (staccatoDipole*)interactionList[i];
			ptr->parse(ifs, RBlist, &externalDipole, &prevDipole, &expenditure, offsetRB, rotationAxis, &sinDipoleAngle);
			if (phase != 0)
			{
				//ptr->start += (phase/360 - phaseObs)*2*PI; ptr->halt += (phase/360 - phaseObs)*2*PI;
				if (ptr->halt < 0)
				{
					ptr->start += (phase/360 - phaseObs)*2*PI; ptr->halt += (phase/360 - phaseObs)*2*PI;
					ptr->deflection += phase * PI/180;
					if (ptr->deflection < 0) ptr->deflection += 2*PI;
					ptr->q = convert(rotationAxis, ptr->deflection);
					ptr->qrev = convert(rotationAxis, -ptr->deflection);
				}
				else
				{
					ptr->start -= (phase/360 - phaseObs)*2*PI; ptr->halt -= (phase/360 - phaseObs)*2*PI;
					ptr->deflection -= phase * PI/180;
					if (ptr->deflection < 0) ptr->deflection -= 2*PI;
					ptr->q = convert(rotationAxis, ptr->deflection);
					ptr->qrev = convert(rotationAxis, -ptr->deflection);
				}
				//ptr->q.printtocout();
			}
		}
		else
		{
			bool rot = 0;
			switch (type)
			{
				case 2:
					interactionList[i] = new bond;
					break;
				case 30:
					interactionList[i] = new angle;
					break;
				case 31:
					interactionList[i] = new angleAxis;
					break;
				case 4:
					interactionList[i] = new dihedral;
					break;
				case 0:
					interactionList[i] = new constraint; rot = 1;
					break;
				case 1:
					interactionList[i] = new constrainAxis; rot = 1;
					break;
				case 50:
					interactionList[i] = new constrainAxisTime; rot = 1;
					break;
			}
			if (rot && !ind) interactionList[i]->parse(ifs, RBlist, offsetRB, center, rotate);
			else interactionList[i]->parse(ifs, RBlist, offsetRB);
			//cout << "interaction parsed\n";
			if (type == 2)
			{
				//cout << i << '\n';
				bond *temp = (bond*)interactionList[i];
				bondList.push_back(*temp);
				adjList2[temp->rbx[0] - RBlist].push_back(temp->rbx[1] - RBlist);
				adjList2[temp->rbx[1] - RBlist].push_back(temp->rbx[0] - RBlist);
				//bondList[bondList.size() - 1].printtocout();
			}
			interactionList[i]->phase = phase/360 - phaseObs;
		}
		if (interactionList[i]->active > 0) isActive = 1;									   
	}
}

void parseSuppressed(ifstream &ifs, int offset)
{
	int n, i, j;
	ifs >> n;
	for (int k = 0; k < n; k++)
	{
		ifs >> i >> j; i += offset; j += offset;
		suppressLJ.insert(suppressf(i, j));
		suppressLJ.insert(suppressf(j, i));
		pair<int, int> p(i, j);
		suppressed.push_back(p);
	}
	ifs >> n;
	for (int k = 0; k < n; k++)
	{
		ifs >> i >> j; i += offset; j += offset;
		umbrellaLJ.insert(suppressf(i, j));
		umbrellaLJ.insert(suppressf(j, i));
		pair<int, int> p(i, j);
		umbrellaed.push_back(p);
	}
}

void parseSpecial(ifstream &ifs, int offset, int fileno)
{
	int n, i, j;
	ifs >> n;
	for (int k = 0; k < n; k++)
	{
		ifs >> i >> j;
		pair<int, int> p(i + offset, j);
		topSpecList.push_back(p);
		mtopSpecList[fileno].push_back(p);
	}
}

void printColour(ofstream &ofs)
{
	for (int i = 0; i < colour.size(); i++) ofs << colour[i] << '\n';
}

void parseColour(ifstream &ifs)
{
	for (int i = 0; i < RBlistSize; i++) ifs >> colour[i];
}

void parseColourFrag(ifstream &ifs, vector<int> &colourList, int k)
{
	int x;
	for (int i = 0; i < k; i++)
	{
		ifs >> x; colourList.push_back(x);
	}
}

void printTop(ofstream &ofs)
{
	ofs << RBlistSize << ' ' << intListSize << '\n';
	printRigid(ofs);
	ofs << '\n';
	printInteractions(ofs);
	ofs << '\n';
	printSuppressed(ofs);
	ofs << '\n';
	printSpecial(ofs);
	ofs << '\n';
	printColour(ofs);
	ofs << '\n';
}

void parseTop(ifstream &ifs, ifstream &ifsConfig)
{
	string s; getline(ifs, s);
	vector3d null(0, 0, 0); quaternion q(1, 0, 0, 0);
	parseRigid(ifs, ifsConfig, 1, 0, null, q); cout << "Rigid parsed\n";
	adjList2.resize(RBlistSize);
	parseInteractions(ifs, 1, 0, 0, null, q, 0, 0); cout << "Interactions parsed\n";
	parseSuppressed(ifs, 0); cout << "Suppressed parsed\n";
	mtopSpecList.resize(1);
	parseSpecial(ifs, 0, 0); cout << "Special parsed\n";
	parseColour(ifs); cout << "Colour parsed\n";
	populateIntTable();
	initDecompAll(); cout << "InitDecompAll\n";
}

void parseTopFrag(ifstream &ifs, ifstream &ifsConfig, int offsetRB, int offsetInt, vector3d &center, quaternion &orient, int fileno, double phase, double offPhase,
vector<int> &colourList)
{
	string s; getline(ifs, s);
	int k = parseRigid(ifs, ifsConfig, 0, offsetRB, center, orient); cout << "Rigid parsed\n";
	parseInteractions(ifs, 0, offsetRB, offsetInt, center, orient, phase, offPhase); cout << "Interactions parsed\n";
	parseSuppressed(ifs, offsetRB);
	parseSpecial(ifs, offsetRB, fileno);
	parseColourFrag(ifs, colourList, k); cout << "Colour parsed\n";
}

void summa(ifstream &ifs, int &_a, int &_b)
{
	int a, b;
	ifs >> a >> b; _a += a; _b += b;
}

void parseTopInteraction(ifstream &ifs, string &s, int &intIdx)
{
	istringstream iss(s);
	//cout << s << '\n';
	string aux, aux2; bool observe = 0, min = 0, equil = 0, umb = 0, active = 0, isGene = 0;
	double ground = 0, k = 0, divisor;
	vector<pair<int, int> > vp(4), vt(4);
	iss >> s;
	if (s == "gen")
	{
		isGene = 1;
		iss >> s;
		//cout << s << '\n';
	}
	iss >> aux;
	if (aux == "obs") {observe = 1;}
	else if (aux == "obsv") {iss >> observe;}
	else
	{
		if (aux == "act")
		{
			iss >> active >> aux;
		}
		ground = atof(aux.c_str());
		iss >> k;
		if (!iss.eof())
		{
			iss >> aux2;
			//cout << aux2 << '\n';
			if (aux2 == "min") min = 1;
			else if (aux2 == "equil") equil = 1;
			else if (aux2 == "umb") {observe = 1; umb = 1;}
		}
	}
	if (s == "constraint")
	{
		string s2;
		getline(ifs, s); istringstream iss2(s);
		vector3d locus;
		iss2 >> s2;
		if (s2 == "curr")
		{
			iss2 >> vp[0].first >> vp[0].second;
			vt[0] = mtopSpecList[vp[0].first][vp[0].second];
			vt[0].second++;
			vt[0].second = massf(RBlist, vt, 0);
			locus = RBlist[vt[0].first].massList[vt[0].second].decomp;
		}
		else
		{
			vp[0].first = atof(s2.c_str());
			iss2 >> vp[0].second >> locus.i >> locus.j >> locus.k;
			vp[0].second++;
			vt[0].second = massf(RBlist, vt, 0);
		}
		constraint * c = new constraint(vp[0], locus, ground, k, RBlist); interactionList[intIdx] = c; intIdx++;
	}
	if (s == "constrainAxis")
	{
		getline(ifs, s); istringstream iss2(s);
		vector3d locus;
		iss2 >> vp[0].first >> vp[0].second >> locus.i >> locus.j >> locus.k;
		vt[0] = mtopSpecList[vp[0].first][vp[0].second];
		vt[0].second++;
		vt[0].second = massf(RBlist, vt, 0);
		constrainAxis * c = new constrainAxis(vt[0], locus, ground, k, RBlist); interactionList[intIdx] = c; intIdx++;
	}
	/*if (s == "constrainAxisTime")
	{
		getline(ifs, s); istringstream iss2(s);
		vector3d startaxis, endaxis;
		double t = 0;
		double duration, offset = 0, on = -1, off = -1;
		iss2 >> vp[0].first >> vp[0].second >> startaxis.i >> startaxis.j >> startaxis.k >> endaxis.i >> endaxis.j >> endaxis.k >> duration >> offset >> on >> off;
		if (duration <= 1)
		{
			duration *= steps;
		}
		if (offset <= 1)
		{
			offset *= steps;
		}		
		if (on <= 1)
		{
			on *= steps;
		}
		if (off <= 1)
		{
			off *= steps;
		}
		vp[0].second++;
		//cout << vp[0].first << ' ' << vp[0].second << '\n';
		vt[0] = mtopSpecList[vp[0].first][vp[0].second];
		constrainAxisTime c(vt[0], startaxis, endaxis, ground, k, (int)duration, t, (int)offset, (int)on, (int)off, RBlist); constrainAxisTimeList.push_back(c);
	}*/
	else if (s == "bond" || s == "chain")
	{
		bool bisChain = 0;
		if (s == "chain") bisChain = 1;
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second >> vp[1].first >> vp[1].second;
		for (int i = 0; i < 2; i++)
		{
			vt[i] = mtopSpecList[vp[i].first][vp[i].second];
			vt[i].second++;
			vt[i].second = massf(RBlist, vt, i);
			//cout << vt[i].first << ' ' << vt[i].second << '\n';
		}
		bond * b = new bond(vt[0], vt[1], ground, k, RBlist); interactionList[intIdx] = b; intIdx++;
		b->isChain = bisChain;
		bond temp = *b;
		bondList.push_back(temp);
		adjList2[temp.rbx[0] - RBlist].push_back(temp.rbx[1] - RBlist);
		adjList2[temp.rbx[1] - RBlist].push_back(temp.rbx[0] - RBlist);
	}
	else if (s == "angle")
	{
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second >> vp[1].first >> vp[1].second >> vp[2].first >> vp[2].second;
		for (int i = 0; i < 3; i++)
		{
			vt[i] = mtopSpecList[vp[i].first][vp[i].second];
			vt[i].second++;
			vt[i].second = massf(RBlist, vt, i);
		}
		angle * a = new angle(vt[0], vt[1], vt[2], ground, k, RBlist); a->obs = observe; interactionList[intIdx] = a; intIdx++;
	}
	else if (s == "gaussAngle")
	{
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second >> vp[1].first >> vp[1].second >> vp[2].first >> vp[2].second;
		for (int i = 0; i < 3; i++)
		{
			vt[i] = mtopSpecList[vp[i].first][vp[i].second];
			vt[i].second++;
			vt[i].second = massf(RBlist, vt, i);
		}
		gaussAngle * a = new gaussAngle(vt[0], vt[1], vt[2], ground, k, RBlist); a->obs = observe; interactionList[intIdx] = a; intIdx++;
	}
	else if (s == "angleAxis")
	{
		getline(ifs, s); istringstream iss2(s);
		vector3d axis;
		iss2 >> vp[0].first >> vp[0].second >> vp[1].first >> vp[1].second >> axis.i >> axis.j >> axis.k;
		for (int i = 0; i < 2; i++)
		{
			vt[i] = mtopSpecList[vp[i].first][vp[i].second];
			vt[i].second++;
			vt[i].second = massf(RBlist, vt, i);
		}
		angleAxis * a = new angleAxis(vt[0], vt[1], axis, ground, k, RBlist); a->obs = observe; interactionList[intIdx] = a; intIdx++;
	}
	else if (s == "dihedral")
	{
		iss >> divisor;
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second >> vp[1].first >> vp[1].second >> vp[2].first >> vp[2].second >> vp[3].first >> vp[3].second;
		for (int i = 0; i < 4; i++)
		{
			vt[i] = mtopSpecList[vp[i].first][vp[i].second];
			vt[i].second++;
			vt[i].second = massf(RBlist, vt, i);
		}
		dihedral * d = new dihedral(vt[0], vt[1], vt[2], vt[3], ground, k, divisor, RBlist); interactionList[intIdx] = d; intIdx++;
	}
	interactionList[intIdx - 1]->min = min;
	interactionList[intIdx - 1]->equil = equil;
	interactionList[intIdx - 1]->umb = umb;
	interactionList[intIdx - 1]->active = active;
	interactionList[intIdx - 1]->isGene = isGene;
}

void parseMultiTop(ifstream &ifs)
{
	double a, b, c, phase = 0, offPhase = 0; int i = 0;
	string s;
	RBlistSize = 0, intListSize = 0; int baseils = 0;
	vector<int> offsetCountRB, offsetCountInt, colourList;
	getline(ifs, s);
	cout << "probing...\n";
	while (true)
	{
		if (s[0] != ';')
		{
			if (s[0] == '!') break;
			i++;
			istringstream iss(s); iss >> s; s += ".top";
			ifstream ifsf(s.c_str());
			offsetCountRB.push_back(RBlistSize);
			offsetCountInt.push_back(intListSize);
			summa(ifsf, RBlistSize, intListSize);
	//cout << RBlistSize << ' ' << intListSize << '\n';
		}
		getline(ifs, s);
	}
	mtopSpecList.resize(i);
	baseils = intListSize;
	getline(ifs, s);
	if (s[0] != '!')
	{
		while (true) //probe super interactions
		{
			if (s[0] != ';')
			{
				if (s[0] == '!') break;
				if (s[0] > 96 && s[0] < 123) intListSize++;
			}
			getline(ifs, s);
		}
	}
	ifs.clear();
	ifs.seekg(0, ios::beg);
	//cout << intListSize << ' ' << baseils << '\n';
	cout << "probed...\n";
	RBlist = new rigidBody[RBlistSize]; interactionList = new interaction*[intListSize];
	adjList2.resize(RBlistSize);
	getline(ifs, s);
	i = 0;
	while (true)
	{
		if (s[0] != ';')
		{
			if (s[0] == '!') break;
			istringstream iss(s); iss >> s >> a >> b >> c;
			vector3d center(a, b, c);
			iss >> a >> b >> c;
			vector3d euler(a, b, c);
			quaternion orient = eulerConv(euler);
			iss >> phase >> offPhase;
			cout << "Parsing " << s << '\n';
			string stop = s + ".top", sconf = s + ".conf";
			ifstream ifsf(stop.c_str()), ifsc(sconf.c_str());
			//cout << offsetCountRB[i] << ' ' << offsetCountInt[i] << '\n';
			parseTopFrag(ifsf, ifsc, offsetCountRB[i], offsetCountInt[i], center, orient, i, phase, offPhase, colourList);
			cout << "Parsed " << s << '\n';
			i++;
		}
		getline(ifs, s);
	}
	/*cout << mtopSpecList.size() << '\n';
	for (int i = 0; i < mtopSpecList.size(); i++)
	{
		cout << i << ": \n";
		for (int j = 0; j < mtopSpecList[i].size(); j++) cout << mtopSpecList[i][j].first << ' ' << mtopSpecList[i][j].second << '\n';
	}*/
	colour.resize(colourList.size());
	for (int i = 0; i < colourList.size(); i++) colour[i] = colourList[i];
	populateIntTable();
	getline(ifs, s);
	if (s[0] != '!')
	{
		while (true) //stitch together disparate structures with new interactions
		{
			if (s[0] != ';')
			{
				if (s[0] == '!') break;
				//cout << s << '\n';
				parseTopInteraction(ifs, s, baseils);
			}
			getline(ifs, s);
		}
	}
	initDecompAll();
	cout << "done parsing\n";
}
