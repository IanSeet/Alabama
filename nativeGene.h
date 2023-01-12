vector<gene> geneList;

void mutate(float rate)
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
				geneList[i].disp.printtocout();
			}
		}
	}
}

void mutate(vector<gene> &gList, float rate)
{
	int a = rate*1000, k = 3;
	double unitDisp = 0.02;
	gList.clear(); gList.resize(geneList.size());
	for (int i = 0; i < gList.size(); i++)
	{
		int c = rand()%1000, pos = rand()%2;
		if (a < c)
		{
			if (gList[i].type == 0)
			{
				gList[i].shift *= 1.0 + pos*unitDisp*(c%k + 1);
			}
			else
			{
				vector3d v(rand()%(k + 1)*unitDisp, rand()%(k + 1)*unitDisp, rand()%(k + 1)*unitDisp);
				geneList[i].disp += pos*v;
				gList[i].disp.printtocout();
			}
		}
	}
}

void makeGenes(bool isConden)
{
	if (isConden)
	{
		for (int i = 0; i < condensedIntList.size(); i++)
		{
			if (condensedIntList[i]->type && !condensedIntList[i]->obs && condensedIntList[i]->isGene)
			{
				gene g;
				g.type = 0;
				g.vecInt.push_back(i);
				geneList.push_back(g);
			}
		}
	}
	else
	{
		decondenseIntList();
		cout << "intListSize: " << intListSize << '\n';
		for (int i = 0; i < intListSize; i++)
		{
			if (interactionList[i]->type && !interactionList[i]->obs && interactionList[i]->isGene)
			{
				gene g;
				g.type = 0;
				g.vecInt.push_back(i);
				geneList.push_back(g);
			}
		}
		cout << "interactions genefied\n";
	}
	for (int i = 0; i < RBlistSize; i++)
	{
		cout << i << ' ' << RBlistSize << '\n';
		if (RBlist[i].isGene)
		{
			cout << "hi\n";
			for (int j = 0; j < RBlist[i].massList.size(); j++)
			{
				cout << j << ' ' << RBlist[i].massList.size() << '\n';
				gene g;
				g.type = 1;
				pair<int, int> p(i, j);
				g.vecMass.push_back(p);
				geneList.push_back(g);
			}
		}
	}
	cout << "rb genefied\n";
}

void printGenes(ofstream &ofs)
{
	ofs << geneList.size() << '\n';
	for (int i = 0; i < geneList.size(); i++)
	{
		geneList[i].print(ofs);
	}
}

void phenotype(bool isConden)
{
	cout << "phenotyping...\n";
	//cout << "cilsize: " << condensedIntList.size() << ' ' << intListSize << '\n';
	if (!isConden) decondenseIntList();
	for (int i = 0; i < geneList.size(); i++)
	{
		//cout << "i: " << i << ' ' << geneList.size() <<  ' ' << isConden << ' ' << geneList[i].type << '\n';
		if (geneList[i].type == 0)
		{
			gene &g = geneList[i];
			if (isConden)
			{
				for (int i = 0; i < g.vecInt.size(); i++)
				{
					condensedIntList[g.vecInt[i]]->ground *= g.shift;
				}
			}
			else
			{
				//cout << "intListSize " << intListSize << ' ' << g.vecInt[0] << '\n'; 
				for (int i = 0; i < g.vecInt.size(); i++)
				{
					interactionList[g.vecInt[i]]->ground *= g.shift;
				}
			}
		}
		else
		{
			gene &g = geneList[i];
			for (int i = 0; i < g.vecMass.size(); i++)
			{
				int a = g.vecMass[i].first, b = g.vecMass[i].second;
				rigidBody &rb = RBlist[a];
				rb.massList[b].initial += g.disp;
			}
		}
	}
	cout << "applyPhysChanges done\n";
	if (isConden) decondenseIntList();
	for (int i = 0; i < RBlistSize; i++)
	{
		RBlist[i].recenter();
	}
}

void parseGenes(ifstream &ifs)
{
	int a;
	ifs >> a; geneList.clear(); geneList.resize(a);
	for (int i = 0; i < geneList.size(); i++)
	{
		geneList[i].parse(ifs);
	}
}
