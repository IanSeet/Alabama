double r3o2 = sqrt(3)/2;

struct complex 
{
	double re, im;
	double m, theta;
	complex(){re = 0; im = 0;}
	complex(double _re, double _im){re = _re; im = _im;}
	void mag() {m = sqrt(re*re + im*im);}
	complex& operator += (complex &C)
	{
		this->re += C.re; this->im += C.im; return *this;
	}
	complex& operator -= (complex &C)
	{
		this->re -= C.re; this->im -= C.im; return *this;
	}
	complex& operator *= (double d)
	{
		this->re *= d; this->im *= d; return *this;
	}
	complex& operator /= (double d)
	{
		this->re /= d; this->im /= d; return *this;
	}
	void convertExp()
	{
		mag(); double d = re/m;
		if (im >= 0) theta = acos(d);
		else theta = -acos(d);
	}
	void convertVec()
	{
		re = m*cos(theta); im = m*sin(theta);
	}
	void pwr(double d)
	{
		convertExp();
		theta *= d; m = pow(m, d);
		convertVec();
	}
	void invert()
	{
		convertExp();
		theta = -theta; m = 1/m;
		convertVec();
	}
	void cubeRoot()
	{
		convertExp();
		theta /= 3; m = cbrt(m);
		convertVec();
	}
	void print()
	{
		cout << "real: " << re << " imaginary: " << im << '\n';
	}
};
complex r120(-0.5, r3o2), r240(-0.5, -r3o2);

complex operator + (complex &c1, complex &c2)
{
	complex c(c1.re + c2.re, c1.im + c2.im);
	return c;
}

complex operator - (complex &c1, complex &c2)
{
	complex c(c1.re - c2.re, c1.im - c2.im);
	return c;
}

complex operator * (complex &c1, complex &c2)
{
	complex c(c1.re*c2.re - c1.im*c2.im, c1.im*c2.re + c1.re*c2.im);
	return c;
}

char cubicSolver(double a, double b, double c, double d, complex *out)
{
	//double D = 18*a*b*c*d - 4*b*b*b*d + b*b*c*c - 4*a*c*c*c - 27*a*a*d*d;
	double D0 = b*b - 3*a*c, D1 = 2*b*b*b - 9*a*b*c + 27*a*a*d;
	complex C(0, 0);
	double root = D1*D1 - 4*D0*D0*D0;
	//cout << D0 << ' ' << D1 << ' ' << root << '\n';
	if (D0 == 0)
	{
		if (D1 < 0) C.re = cbrt(D1);
		else 
		{
			C.re = cbrt(D1);
			return 0;
		}
		return 1;
		//cout << "case 1: " << C.re << ' ' << C.im << '\n';
	}
	else
	{
		if (root < 0)
		{
			C.re = D1/2;
			C.im = sqrt(-root)/2;
			C.cubeRoot();
			//cout << "case 2: " << C.re << ' ' << C.im << '\n';
		}
		else
		{
			C.re = cbrt((D1 + sqrt(root))/2);
			//cout << "case 3: " << C.re << ' ' << C.im << '\n';
			if (C.re == 0) C.re = cbrt((D1 - sqrt(root))/2);
		}
		complex I = C; I.invert(); I *= D0;
		double mult = -1/(3*a); //cout << mult << '\n';
		out[0] = C + I; out[0].re += b; out[0] *= mult; //out[0].print(); //C.print(); I.print(); //C.mag(); cout << C.m << '\n';
		complex C120 = C * r120;
		I = C120; I.invert(); I *= D0; //C120.print(); I.print(); //C.mag(); cout << C.m << '\n';
		out[1] = C120 + I; out[1].re += b; out[1] *= mult;
		complex C240 = C * r240;
		I = C240; I.invert(); I *= D0; //C240.print(); I.print(); //C.mag(); cout << C.m << '\n';
		out[2] = C240 + I; out[2].re += b; out[2] *= mult;
		return 2;
	}
}

double determinant(double M[][2])
{
	//cout << "det: " << M[0][0] << ' ' << M[0][1] << ' ' << M[1][0] << ' ' << M[1][1] << ' ' << M[0][0] * M[1][1] - M[0][1] * M[1][0] << '\n';
	return M[0][0] * M[1][1] - M[0][1] * M[1][0];
}

bool matrix2Dsolver(double M[][2], double v[], double out[])
{
	if (v[0] == 0 && v[1] == 0) return false;
	double det = determinant(M);
	if (det == 0) return false;
	det = 1/det;
	out[0] = det*(M[1][1]*v[0] - M[0][1]*v[1]);
	out[1] = det*(M[0][0]*v[1] - M[1][0]*v[0]);
	return true;
}

bool symEigenSolver(double M[][3], double EVal[], double EVec[][3])
{
	double a = -1, b = M[0][0] + M[1][1] + M[2][2];
	double xy2 = M[0][1] * M[0][1], xz2 = M[0][2] * M[0][2], yz2 = M[1][2] * M[1][2];
	double c = xy2 + xz2 + yz2 - M[0][0]*M[1][1] - M[1][1]*M[2][2] - M[0][0]*M[2][2];
	double xyxzyz = M[0][1] * M[0][2] * M[1][2], d = M[0][0]*M[1][1]*M[2][2] + 2*xyxzyz - M[0][0]*yz2 - M[1][1]*xz2 - M[2][2]*xy2;
	complex out[3];
	//cout << a << ' ' << b << ' ' << c << ' ' << d << '\n';
	char o = cubicSolver(a, b, c, d, out);
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) EVec[i][j] = 0;
//	for (int i = 0; i < 3; i++) cout << out[i].re << ' ' << out[i].im << '\n';
//	cout << (int)o << '\n';
	if (o == 3)
	{
		cout << "M is not hermitian\n"; return 0;
	}
	else if (o == 0)
	{
		for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
		if (i == j) EVec[i][j] = 1; else EVec[i][j] = 0;
		for (int i = 0; i < 3; i++) EVal[i] = M[0][0];
	}
	else if (o == 1)
	{
		double l = out[1].re, S[2][2], v[2], s[2];
		vector3d z, x(1, 0, 0), y;
		for (int i = 0; i < 3; i++)
		{
			int r = (i + 1)%3, p = (i + 2)%3;
			S[0][0] = M[i][r]; S[0][1] = M[i][p]; v[0] = l - M[i][i];
			S[1][0] = M[r][r] - l; S[1][1] = M[r][p]; v[1] = -M[i][r];
			if (matrix2Dsolver(S, v, s))
			{
				EVec[0][i] = 1; EVec[0][r] = s[0]; EVec[0][p] = s[1];
				double mag = 0;
				for (int j = 0; j < 3; j++) mag += (EVec[0][j])*(EVec[0][j]);
				mag = sqrt(mag);
				if (mag != 0)
				{
					z.i = EVec[0][0] / mag;
					z.j = EVec[0][1] / mag;
					z.k = EVec[0][2] / mag;
				}
				break;
			}
		}
		x = x - z * (x * z); y = cross(x, z);
		copy(z, &(EVec[0][0]));
		copy(x, &(EVec[1][0]));
		copy(y, &(EVec[2][0]));
		EVal[0] = out[1].re;
		EVal[1] = out[0].re; EVal[2] = out[0].re;
	}
	else
	{		
		vector3d x, y, z;
		double l = out[0].re, l2 = out[1].re, S[2][2], v[2], s[2]; bool done = 0;
		for (int i = 0; i < 3; i++)
		{
			int r, p;
			if (i == 0) {r = 1; p = 2;}
			if (i == 1) {r = 0; p = 2;}
			if (i == 2) {r = 0; p = 1;}
			S[0][0] = M[0][r]; S[0][1] = M[0][p]; v[0] = -M[0][i];
			S[1][0] = M[1][r]; S[1][1] = M[1][p]; v[1] = -M[1][i];
			if (i < 2) v[i] += l;
			if (i != 1) S[1][0] -= l;
			if (i != 0) S[0][0] -= l;
			//cout << "S00: " << M[0][r] << ' ' << l << ' ' << S[0][0] << '\n';
			if (determinant(S) == 0 || v[0] == 0 && v[1] == 0)
			{
				EVec[0][i] = 1; EVec[0][r] = 0; EVec[0][p] = 0;
			}
			else if (matrix2Dsolver(S, v, s))
			{
				EVec[0][i] = 1; EVec[0][r] = s[0]; EVec[0][p] = s[1];
				double mag = 0;
				for (int j = 0; j < 3; j++) mag += (EVec[0][j])*(EVec[0][j]);
				mag = sqrt(mag);
				if (mag != 0) for (int j = 0; j < 3; j++) EVec[0][j] /= mag;
			}	
			for (int i2 = 0; i2 < 3; i2++)
			{
				//cout << "cycle: " << i << ' ' << i2 << '\n';
				int r2, p2;
				if (i2 == 0) {r2 = 1; p2 = 2;}
				if (i2 == 1) {r2 = 0; p2 = 2;}
				if (i2 == 2) {r2 = 0; p2 = 1;}
				S[0][0] = M[0][r2]; S[0][1] = M[0][p2]; v[0] = -M[0][i2];
				S[1][0] = M[1][r2]; S[1][1] = M[1][p2]; v[1] = -M[1][i2];
				if (i2 < 2) v[i2] += l2;
				if (i2 != 1) S[1][0] -= l2;
				if (i2 != 0) S[0][0] -= l2;
				if (determinant(S) == 0 || v[0] == 0 && v[1] == 0)
				{
					EVec[1][i2] = 1; EVec[1][r2] = 0; EVec[1][p2] = 0;
					x.assign(&EVec[0][0]); y.assign(&EVec[1][0]);
					z = cross(x, y);
					vector3d zout = matrixMult3d(M, z); double eigenz;
					if (zout.i != 0) eigenz = zout.i/z.i;
					else if (zout.j != 0) eigenz = zout.j/z.j;
					else if (zout.k != 0) eigenz = zout.k/z.k;
					//cout << "zout: "; zout.printtocout();
					//cout << eigenz << ' ' << out[2].re << '\n';
					//cout << "NullAssign: "; y.printtocout(); cout << x* y << '\n';
					if (x * y < 0.001 && x * y > -0.001 && fabs(eigenz - out[2].re)/eigenz < 0.01) 
					{
						done = 1;
						break;
					}
				}
				else if (matrix2Dsolver(S, v, s))
				{
					EVec[1][i2] = 1; EVec[1][r2] = s[0]; EVec[1][p2] = s[1];
					double mag = 0;
					for (int j = 0; j < 3; j++) mag += (EVec[1][j])*(EVec[1][j]);
					mag = sqrt(mag);
					if (mag != 0) for (int j = 0; j < 3; j++) EVec[1][j] /= mag;
					x.assign(&EVec[0][0]); y.assign(&EVec[1][0]);
					z = cross(x, y);
					vector3d zout = matrixMult3d(M, z); double eigenz;
					if (zout.i != 0) eigenz = zout.i/z.i;
					else if (zout.j != 0) eigenz = zout.j/z.j;
					else if (zout.k != 0) eigenz = zout.k/z.k;
					//cout << "zout: "; zout.printtocout();
					//cout << eigenz << ' ' << out[2].re << '\n';
					//cout << "Assign: "; x.printtocout(); y.printtocout(); cout << x * y << '\n';
					if (x * y < 0.001 && x * y > -0.001 && fabs(eigenz - out[2].re)/eigenz < 0.01) 
					{
						done = 1;
						break;
					}
				}
			}
			if (done) break;
		}
		if (!done)
		{
			cout << "Warning, could not find principal axes\n";
			for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) EVec[i][j] = 0;
			for (int i = 0; i < 3; i++) EVec[i][i] = 1;
			return 0;
		}
		else
		{
			EVec[2][0] = z.i; EVec[2][1] = z.j; EVec[2][2] = z.k;
			if (x * y != 0)
			{
				y = y - (x * y) * x;
				EVec[1][0] = y.i; EVec[1][1] = y.j; EVec[1][2] = y.k;
			}
		}
		//for (int i = 0; i < 2; i++)
		//{
		//	for (int j = 0; j < 3; j++) cout << EVec[i][j] << ' ';
		//	cout << '\n';
		//}
		for (int i = 0; i < 3; i++) EVal[i] = out[i].re;
	}
	return 1;
}
