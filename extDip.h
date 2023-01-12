inline void rotateED() {rotateEDPrime(); rotateOD();}

inline void constRot(vector3d &targetDipole, vector3d &initDipole, vector3d &_prevDipole, 
const bool clockwise, const double timeFrac, double _offScale)
{
	float mult = 1; if (!clockwise) mult = -1;
	quaternion q(rotationAxis, mult*(timeFrac - offPhaseFactor)*2*chain*PI*stepMult*_offScale);
	sinDipoleAngle = mult*(timeFrac - offPhaseFactor)*2*chain*PI*stepMult; sinDipoleAngle = (roundf(sinDipoleAngle*1000) + 0.0)/1000;
	if (printStep) cout << "ED: " << sinDipoleAngle << '\n';
	_prevDipole = targetDipole;
	targetDipole = initDipole; targetDipole.rotate(q); targetDipole.norm();
}

inline void sinRot(vector3d &targetDipole, vector3d &initDipole, vector3d &_prevDipole, 
const bool clockwise, const bool circular, const double timeFrac, double _offScale)
{
	double angleConst = PI*sinDipoleFactor, angle = timeFrac*2*chain*PI;
	float mult = 1; if (!clockwise) mult = -1;
	double x = timeFrac; if (x > 0.5 && x <= 1) x = 1 - x; if (x > 1) x -= 1;
	quaternion q;
	if (circular)
	{
		quaternion q1(rotationAxis, mult*(x - offPhaseFactor)*chain*PI*stepMult*2*sinDipoleFactor*_offScale); q = q1;
		sinDipoleAngle = mult*(x - offPhaseFactor)*chain*PI*stepMult*2*sinDipoleFactor; sinDipoleAngle = (roundf(sinDipoleAngle*1000) + 0.0)/1000;
	}
	else
	{
		double offset = angleConst*_offScale*mult*(1 - cos(offPhaseFactor*2*chain*PI))/2;
		quaternion q1(rotationAxis, angleConst*_offScale*mult*(1 - cos(angle))/2 - offset); q = q1;
		sinDipoleAngle = angleConst*mult*(1 - cos(angle))/2 - offset; sinDipoleAngle = (roundf(sinDipoleAngle*1000) + 0.0)/1000;
	}
	if (printStep) cout << "ED: " << sinDipoleAngle << '\n';
	_prevDipole = targetDipole;
	targetDipole = initDipole; targetDipole.rotate(q); targetDipole.norm();
}

inline void pauseRot(vector3d &targetDipole, vector3d &initDipole, vector3d &_prevDipole, 
const bool clockwise, const double timeFrac, double freq, double _offScale)
{
	float mult = 1; if (!clockwise) mult = -1;
	
	double period = 1.0/freq, tf = timeFrac - offPhaseFactor, seg = tf/period;
	int disc = tf/period; double remainder = seg - disc;
	
	double rem = 0.5;
	if (remainder < 0.5) rem = remainder; rem *= 2;
	double adjustedTimeFrac = (disc + rem)*period;
	
	quaternion q(rotationAxis, mult*(adjustedTimeFrac)*2*chain*PI*stepMult*_offScale);
	sinDipoleAngle = mult*(adjustedTimeFrac)*2*chain*PI*stepMult; sinDipoleAngle = (roundf(sinDipoleAngle*1000) + 0.0)/1000;
	if (printStep) cout << "PED: " << sinDipoleAngle << '\n';
	_prevDipole = targetDipole;
	targetDipole = initDipole; targetDipole.rotate(q); targetDipole.norm();
}

inline void pauseSinRot(vector3d &targetDipole, vector3d &initDipole, vector3d &_prevDipole, 
const bool clockwise, const bool circular, const double timeFrac, double freq, double _offScale)
{
	double angleConst = PI*sinDipoleFactor, angle = timeFrac*2*chain*PI;
	float mult = 1; if (!clockwise) mult = -1;
	
	double period = 1/freq, tf = timeFrac - offPhaseFactor, seg = tf/period;
	int disc = tf/period; double remainder = seg - disc;
	
	double rem = 0.5;
	if (remainder < 0.5) rem = remainder; rem *= 2;
	double adjustedTimeFrac = (disc + rem)*period;
	
	double x = adjustedTimeFrac; if (x > 0.5 && x <= 1) x = 1 - x; if (x > 1) x -= 1;
	quaternion q;
	if (circular)
	{
		quaternion q1(rotationAxis, mult*(x - offPhaseFactor)*chain*PI*stepMult*2*sinDipoleFactor*_offScale); q = q1;
		sinDipoleAngle = mult*(x - offPhaseFactor)*chain*PI*stepMult*2*sinDipoleFactor; sinDipoleAngle = (roundf(sinDipoleAngle*1000) + 0.0)/1000;
	}
	else
	{
		double offset = angleConst*_offScale*mult*(1 - cos(offPhaseFactor*2*chain*PI))/2;
		quaternion q1(rotationAxis, angleConst*_offScale*mult*(1 - cos(angle))/2 - offset); q = q1;
		sinDipoleAngle = angleConst*mult*(1 - cos(angle))/2 - offset; sinDipoleAngle = (roundf(sinDipoleAngle*1000) + 0.0)/1000;
	}
	if (printStep) cout << "SPED: " << sinDipoleAngle << '\n';
	_prevDipole = targetDipole;
	targetDipole = initDipole; targetDipole.rotate(q); targetDipole.norm();
}

inline void rotateEDPrime()
{
	double timeFrac = ((double)currStep)/maxTimeStep + offPhaseFactor;
	if (pauseRotFactor > 0)
	{
		if (!sinDipole)
		{
			pauseRot(externalDipole, initialDipole, prevDipole, clockwise, timeFrac, pauseRotFactor, 1);
		}
		else
		{
			pauseSinRot(externalDipole, initialDipole, prevDipole, clockwise, circular, timeFrac, pauseRotFactor, 1);
		}
	}
	else
	{
		if (!sinDipole)
		{
			constRot(externalDipole, initialDipole, prevDipole, clockwise, timeFrac, 1);
		}
		else
		{
			sinRot(externalDipole, initialDipole, prevDipole, clockwise, circular, timeFrac, 1);
		}
	}
}

inline void rotateOD()
{
	double timeFrac = ((double)currStep)/maxTimeStep + offPhaseFactor;
	if (pauseRotFactor > 0)
	{
		if (!sinDipole)
		{
			pauseRot(offsetDipole, initOffDipole, prevOffDipole, clockwise, timeFrac, pauseRotFactor, offScale);
		}
		else
		{
			pauseSinRot(offsetDipole, initOffDipole, prevOffDipole, clockwise, circular, timeFrac, pauseRotFactor, offScale);
		}
	}
	else
	{
		if (!sinDipole)
		{
			constRot(offsetDipole, initOffDipole, prevOffDipole, clockwise, timeFrac, offScale);
		}
		else
		{
			sinRot(offsetDipole, initOffDipole, prevOffDipole, clockwise, circular, timeFrac, offScale);
		}
	}
}
