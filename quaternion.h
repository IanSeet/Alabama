struct vector3 {double i, j, k;};

struct quaternion
{
	double i, j, k, a; //a is real part
	inline quaternion()
	{
		a = 1; i = 0; j = 0; k = 0;
	}
	inline quaternion(double _a, double _i, double _j, double _k)
	{
		a = _a; i = _i; j = _j; k = _k;
	}
	inline quaternion(double _i, double _j, double _k)
	{
		a = 0.0; i = _i; j = _j; k = _k;
	}
	inline quaternion(vector3 &v)
	{
		a = 0.0; i = v.i; j = v.j; k = v.k;
	}
	inline quaternion(double *arr)
	{
		a = arr[0]; i = arr[1]; j = arr[2]; k = arr[3];
	}
	inline quaternion(vector3 &v, double theta) //converts axis-angle rotation to quaternion
	{
		a = cos(theta/2);
		i = v.i*sin(theta/2); j = v.j*sin(theta/2); k = v.k*sin(theta/2);
	}
	inline void assign(double _a, double _i, double _j, double _k)
	{
		a = _a; i = _i; j = _j; k = _k;
	}
	inline void assign(vector3 &v)
	{
		a = 0.0; i = v.i; j = v.j; k = v.k;
	}
	inline void convert(double x, double y, double z, double theta) //converts axis-angle rotation to quaternion
	{
		a = cos(theta/2);
		i = x*sin(theta/2); j = y*sin(theta/2); k = z*sin(theta/2);
	}
	inline double axisx() {return i/sqrt(i*i + j*j + k*k);} //Find xyz coordinates of axis of rotation
	inline double axisy() {return j/sqrt(i*i + j*j + k*k);}
	inline double axisz() {return k/sqrt(i*i + j*j + k*k);}
	inline double anglerot() {return 2*atan2(sqrt(i*i + j*j + k*k), a);} //Find angle of rotation in radians about aformentioned axis
	inline double toEuler(char c) //Euler angle conversion to radians (Tait-Bryan); 0 = yaw, 1 = pitch, 2 = roll
	{
		if (c == 0) return atan2(2*(a*i + j*k), 1 - 2*(i*i + j*j))*180/PI;
		if (c == 1) return asin(2*(a*j - i*k))*180/PI;
		return atan2(2*(a*k + i*j), 1 - 2*(j*j + k*k))*180/PI;
	}
	inline void conjugate()
	{
		i = -i; j = -j; k = -k;
	}
	inline quaternion& operator += (quaternion &q)
	{
		this->i += q.i; this->j += q.j; this->k += q.k; this->a += q.a; return *this;
	}
	inline quaternion& operator -= (quaternion &q)
	{
		this->i -= q.i; this->j -= q.j; this->k -= q.k; this->a -= q.a; return *this;
	}
	inline quaternion& operator *= (double d) //scalar multiplication
	{
		this->i *= d; this->j *= d; this->k *= d; this->a *= d; return *this;
	}
	inline quaternion& operator /= (double d) //scalar division
	{
		this->i /= d; this->j /= d; this->k /= d; this->a /= d; return *this;
	}
	inline double mag() {return sqrt(a*a + i*i + j*j + k*k);} //magnitude
	inline void norm() //normalise
	{
		double x = mag();
		if (x != 0) *this /= x;
	}
	inline void print(ofstream &ofs)
	{
		ofs << a << ' ' << i << ' ' << j << ' ' << k << '\n';
	}
	inline void read(ifstream &ifs)
	{
		ifs >> a >> i >> j >> k;
	}
	inline void printtocout()
	{
		cout << a << ' ' << i << ' ' << j << ' ' << k << '\n';
	}
	inline void clear()
	{
		a = 1; i = 0; k = 0; j = 0;
	}
};

inline void assign(double o[4], quaternion &q)
{
	o[0] = q.a; o[1] = q.i; o[2] = q.j; o[3] = q.k;
}

inline quaternion conjugate(quaternion &q)
{
	quaternion out(q.a, -q.i, -q.j, -q.k); return out;
}

inline quaternion operator * (const quaternion &q1, const quaternion &q2) //hamilton product
{
	quaternion q;
	q.a = q1.a*q2.a - q1.i*q2.i - q1.j*q2.j - q1.k*q2.k;
	q.i = q1.a*q2.i + q1.i*q2.a + q1.j*q2.k - q1.k*q2.j;
	q.j = q1.a*q2.j - q1.i*q2.k + q1.j*q2.a + q1.k*q2.i;
	q.k = q1.a*q2.k + q1.i*q2.j - q1.j*q2.i + q1.k*q2.a;
//	q.norm();
	return q;
}

inline quaternion operator + (const quaternion &q1, const quaternion &q2)
{
	quaternion q(q1.a + q2.a, q1.i + q2.i, q1.j + q2.j, q1.k + q2.k);
	return q;
}

inline quaternion operator - (const quaternion &q1, const quaternion &q2)
{
	quaternion q(q1.a - q2.a, q1.i - q2.i, q1.j - q2.j, q1.k - q2.k);
	return q;
}

inline quaternion operator * (const vector<vector<double> > &M, const quaternion &q)
{
	double arr[4];
	if (M.size() == 4 && M[0].size() == 4)
	{
		for (int i = 0; i < 4; i++)
		{
			arr[i] = M[i][0]*q.a + M[i][1]*q.i + M[i][2]*q.j + M[i][3]*q.k; //cout << arr[i] << '\n';
		}
	}
	quaternion out(arr);
	return out;
}

inline quaternion operator * (const double M[][4], const quaternion &q)
{
	quaternion out(M[0][0]*q.a + M[0][1]*q.i + M[0][2]*q.j + M[0][3]*q.k,
			M[1][0]*q.a + M[1][1]*q.i + M[1][2]*q.j + M[1][3]*q.k,
			M[2][0]*q.a + M[2][1]*q.i + M[2][2]*q.j + M[2][3]*q.k,
			M[3][0]*q.a + M[3][1]*q.i + M[3][2]*q.j + M[3][3]*q.k);
	return out;
}

inline quaternion operator * (const quaternion &q1, double k) //scalar multiplication - possibly unstable
{
	quaternion q(q1.a*k, q1.i*k, q1.j*k, q1.k*k); return q;
}

inline quaternion operator * (double k, const quaternion &q1) //scalar multiplication - possibly unstable
{
	quaternion q(q1.a*k, q1.i*k, q1.j*k, q1.k*k); return q;
}

inline double dot(const quaternion &q1, const quaternion &q2)
{
	return q1.a*q2.a + q1.i*q2.i + q1.j*q2.j + q1.k*q2.k;
}

struct vector3d : vector3 //3D vector
{
	inline vector3d()
	{
		i = 0; j = 0; k = 0;
	}
	inline vector3d(double _i, double _j, double _k)
	{
		i = _i; j = _j; k = _k;
	}
	inline vector3d(double *arr)
	{
		i = arr[0]; j = arr[1]; k = arr[2];
	}
	inline void assign(double _i, double _j, double _k)
	{
		i = _i; j = _j; k = _k;
	}
	inline void assign(vector3d &v)
	{
		i = v.i; j = v.j; k = v.k;
	}
	inline void assign(double *arr)
	{
		i = arr[0]; j = arr[1]; k = arr[2];
	}
	inline vector3d& operator += (const vector3d & v3d)
	{
		this->i += v3d.i; this->j += v3d.j; this->k += v3d.k; return *this;
	}
	inline vector3d& operator -= (const vector3d & v3d)
	{
		this->i -= v3d.i; this->j -= v3d.j; this->k -= v3d.k; return *this;
	}
	inline vector3d& operator *= (double d) //scalar multiplication
	{
		this->i *= d; this->j *= d; this->k *= d; return *this;
	}
	inline vector3d& operator *= (bool *arr) //purge variables
	{
		if (!arr[0]) this->i = 0;
		if (!arr[1]) this->j = 0;
		if (!arr[2]) this->k = 0;
		return *this;
	}
	inline vector3d& operator /= (double d) //scalar division
	{
		this->i /= d; this->j /= d; this->k /= d; return *this;
	}
	inline vector3d& operator = (quaternion& q) //assign quaternion to vector
	{
		this->i = q.i; this->j = q.j; this->k = q.k; return *this;
	}
	inline double square() {return i*i + j*j + k*k;}
	inline double mag() {return sqrt(square());} //magnitude
	inline void norm() //normalise
	{
		double x = mag();
		if (x != 0) *this /= x;
	}
	inline vector3d cross(vector3d &v1, vector3d &v2) //cross product
	{
		vector3d v3d(v1.j*v2.k - v1.k*v2.j, v1.k*v2.i - v1.i*v2.k, v1.i*v2.j - v1.j*v2.i); return v3d;
	}
	inline vector3d add(const vector3d & v1, const vector3d & v2)
	{
		vector3d v3d(v1.i + v2.i, v1.j + v2.j, v1.k + v2.k);
		return v3d;
	}
	inline void rotate(quaternion& q)
	{
		//double x = mag();
		vector3d r(q.i, q.j, q.k), w(i * q.a, j * q.a, k * q.a), v = r.add(r.cross(r, *this), w);
		vector3d t = r.cross(r, v); t *= 2;
		*this += t;
		//double y = mag(); if (y != 0) *this *= x/y;
	}
	inline void printtocout()
	{
		cout << i << ' ' << j << ' ' << k << '\n';
	}
	inline void print(ofstream &ofs)
	{
		ofs << i << ' ' << j << ' ' << k << '\n';
	}
	inline void printNoNewline(ofstream &ofs)
	{
		ofs << i << ' ' << j << ' ' << k << ' ';
	}
	inline void read(ifstream &ifs)
	{
		ifs >> i >> j >> k;
	}
	inline void round(double threshold)
	{
		int ir = roundf(i*threshold)/threshold;
		int jr = roundf(j*threshold)/threshold;
		int kr = roundf(k*threshold)/threshold;
		i = ir + 0.0; j = jr + 0.0; k = kr + 0.0;
	}
};

inline vector3d operator + (const vector3d & v1, const vector3d & v2)
{
	vector3d v3d(v1.i + v2.i, v1.j + v2.j, v1.k + v2.k);
	return v3d;
}

inline vector3d operator - (const vector3d & v1, const vector3d & v2)
{
	vector3d v3d(v1.i - v2.i, v1.j - v2.j, v1.k - v2.k);
	return v3d;
}

inline double operator * (const vector3d & v1, const vector3d & v2) //dot product
{
	return v1.i*v2.i + v1.j*v2.j + v1.k*v2.k;
}

inline bool operator == (const vector3d& v, const vector3d &v2) //compare vectors
{
	if (v2.i == v.i && v2.j == v.j && v2.k == v.k) return 1;
	else return 0;
}

inline vector3d vmult (const vector3d & v1, const vector3d & v2) //direct multiplication
{
	vector3d v(v1.i*v2.i, v1.j*v2.j, v1.k*v2.k);
	return v;
}

inline vector3d operator * (const vector3d & v1, double k) //scalar multiplication
{
	vector3d v(v1.i*k, v1.j*k, v1.k*k); return v;
}

inline vector3d operator * (double k, const vector3d & v1) //scalar multiplication
{
	vector3d v(v1.i*k, v1.j*k, v1.k*k); return v;
}

inline vector3d operator * (const vector3d &v, double matrix[][3]) //matrix multiplication
{
	double arr[3];
	for (int i = 0; i < 3; i++) arr[i] = v.i * matrix[0][i] + v.j * matrix[1][i] + v.k * matrix[2][i];
	vector3d t(arr); return t;
}

inline vector3d operator * (double matrix[][3], const vector3d &v) //matrix multiplication
{
	double arr[3];
	for (int i = 0; i < 3; i++) arr[i] = v.i * matrix[i][0] + v.j * matrix[i][1] + v.k * matrix[i][2];
	vector3d t(arr); return t;
}

inline vector3d operator / (const vector3d & v1, double k) //scalar division - possibly unstable
{
	vector3d v(v1.i/k, v1.j/k, v1.k/k); return v;
}

inline vector3d operator / (const vector3d & v1, double *v2) //scalar division - possibly unstable
{
	//cout << v2[0] << '\n';
	vector3d v(v1.i/v2[0], v1.j/v2[1], v1.k/v2[2]); return v;
}

inline vector3d cross(vector3d &v1, vector3d &v2) //cross product
{
	vector3d v3d(v1.j*v2.k - v1.k*v2.j, v1.k*v2.i - v1.i*v2.k, v1.i*v2.j - v1.j*v2.i); return v3d;
}

inline void copy(vector3d &v, double *arr)
{
	arr[0] = v.i; arr[1] = v.j; arr[2] = v.k;
}

inline void copy(quaternion &q, double *arr)
{
	arr[0] = q.a; arr[1] = q.i; arr[2] = q.j; arr[3] = q.k;
}

inline double dist(vector3d &v1, vector3d &v2) //euclidean distance
{
	vector3d v = v1 - v2; return v.mag();
}

inline double calcAngle(vector3d &v1, vector3d &v2, vector3d &v3)
{
	vector3d a1 = v1 - v2, a2 = v3 - v2;
	double a1mag = a1.mag(), a2mag = a2.mag();
	if (a1mag == 0 || a2mag == 0) return 0;
	a1 /= a1mag; a2 /= a2mag;
	double d = a1 * a2; if (d > 1) d = 1; if (d < -1) d = -1;
	return acos(d);
}

inline quaternion eulerConv(vector3d euler) //convert Tait-Bryan Euler angles in degrees to quaternion form; i = yaw, j = pitch, k = roll
{
	euler /= (360/PI);
	double cy = cos(euler.i);
	double sy = sin(euler.i);
	double cp = cos(euler.j);
	double sp = sin(euler.j);
	double cr = cos(euler.k);
	double sr = sin(euler.k);
	quaternion q;
	q.a = cy * cp * cr + sy * sp * sr;
	q.i = cy * cp * sr - sy * sp * cr;
	q.j = sy * cp * sr + cy * sp * cr;
	q.k = sy * cp * cr - cy * sp * sr;
	return q;
}

inline quaternion convert(vector3d &v, double theta) //converts axis-angle rotation to quaternion
{
	double t2 = theta/2, cost2 = cos(t2), sint2 = sin(t2);
	quaternion q(cost2, v.i*sint2, v.j*sint2, v.k*sint2);
	return q;
}

inline quaternion convert(vector3d &v, double cost2, double sint2) //converts axis-angle rotation to quaternion
{
	quaternion q(cost2, v.i*sint2, v.j*sint2, v.k*sint2);
	return q;
}

inline double minDistSquare(vector3d &l, vector3d &p) //returns square of minimum distance assuming equation of line l passes through origin and that line vector is normalised
{
	double dot = p*l;
	return (p*p) - dot*dot;
}

inline void matrixMult(vector<vector<double> > &M1, vector<vector<double> > &M2, vector<vector<double> > &output)
{
	if (M1.size() == 0 || M1[0].size() != M2.size()) return;
	int i, j, k; output.resize(M1.size()); for (i = 0; i < output.size(); i++) output[i].assign(M2[0].size(), 0);
	for (i = 0; i < M1.size(); i++) for (j = 0; j < M2[i].size(); j++) for (k = 0; k < M1[i].size(); k++)
	{
		output[i][j] += M1[i][k] * M2[k][j];
	}
}

inline void transpose(vector<vector<double> > &M, vector<vector<double> > &Mt)
{
	if (M.size() == 0 || M[0].size() == 0) return;
	int i, j; Mt.resize(M[0].size()); for (i = 0; i < Mt.size(); i++) Mt[i].resize(M.size());
	for (i = 0; i < M.size(); i++) for (j = 0; j < M[i].size(); j++) Mt[j][i] = M[i][j];
}

inline void Sfier(quaternion &q, vector<vector<double> > &M)
{
	M.resize(4);
	M[0].resize(4); M[0][0] = q.a; M[0][1] = -q.i; M[0][2] = -q.j; M[0][3] = -q.k;
	M[1].resize(4); M[1][0] = q.i; M[1][1] = q.a; M[1][2] = -q.k; M[1][3] = q.j;
	M[2].resize(4); M[2][0] = q.j; M[2][1] = q.k; M[2][2] = q.a; M[2][3] = -q.i;
	M[3].resize(4); M[3][0] = q.k; M[3][1] = -q.j; M[3][2] = q.i; M[3][3] = q.a;
}

inline quaternion Smult(int l, const quaternion &q)
{
	if (l == 1)
	{
		quaternion out(-q.i, q.a, q.k, -q.j); return out;
	}
	else if (l == 2)
	{
		quaternion out(-q.j, -q.k, q.a, q.i); return out;
	}
	else
	{
		quaternion out(-q.k, q.j, -q.i, q.a); return out;
	}
}

inline void QtoRM(quaternion &q, double RM[][3])
{
	double aa = q.a*q.a, ab = q.a*q.i, ac = q.a*q.j, ad = q.a*q.k;
	double bb = q.i*q.i, bc = q.i*q.j, bd = q.i*q.k;
	double cc = q.j*q.j, cd = q.j*q.k;
	double dd = q.k*q.k;
	RM[0][0] = aa + bb - cc - dd; RM[0][1] = 2*(bc - ad); RM[0][2] = 2*(bd + ac);
	RM[1][0] = 2*(bc + ad); RM[1][1] = aa - bb + cc - dd; RM[1][2] = 2*(cd - ab);
	RM[2][0] = 2*(bd - ac); RM[2][1] = 2*(cd + ab); RM[2][2] = aa - bb - cc + dd;
}

inline vector3d matrixMult3d(double M[][3], vector3d &v)
{
	vector3d out(M[0][0]*v.i + M[0][1]*v.j + M[0][2]*v.k,
			M[1][0]*v.i + M[1][1]*v.j + M[1][2]*v.k,
			M[2][0]*v.i + M[2][1]*v.j + M[2][2]*v.k);
	return out;
}

inline void matrixMult3dSafe(double M[][3], vector3d &v, vector3d &o)
{
	o.i = M[0][0]*v.i + M[0][1]*v.j + M[0][2]*v.k;
	o.j = M[1][0]*v.i + M[1][1]*v.j + M[1][2]*v.k;
	o.k = M[2][0]*v.i + M[2][1]*v.j + M[2][2]*v.k;
}

bool gluInvertMatrix(double m[16], double invOut[16])
{
    double inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    return true;
}

inline quaternion matrix4dInvertMult(double *I, quaternion &q)
{
	double O[16], out[4] = {0, 0, 0, 0}, in[4];
	assign(in, q);
	if (gluInvertMatrix(I, O))
	{
		for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) out[i] += O[i*4 + j] * in[j];
	}
	quaternion qout(out); return qout;
}

inline quaternion omegaConvert(vector3d &omega, quaternion &orient)
{
	double M[16], Inv[16], q0 = orient.a, q1 = orient.i, q2 = orient.j, q3 = orient.k;
	M[0] = q0; M[1] = q1; M[2] = q2; M[3] = q3;
	M[4] = -q1; M[5] = q0; M[6] = q3; M[7] = -q2;
	M[8] = -q2; M[9] = -q3; M[10] = q0; M[11] = q1;
	M[12] = -q3; M[13] = q2; M[14] = -q1; M[15] = q0;
	omega *= 2;
	quaternion qin(omega), qout = matrix4dInvertMult(M, qin);
	return qout;
}

namespace std
{
	template <> struct hash<vector3d>
	{
		inline size_t operator () (const vector3d& k) const
		{
			return ((hash<double>()(k.i) ^ (hash<double>()(k.j) << 1)) >> 1) ^ (hash<double>()(k.k) << 1);
		}
	};
}
