#ifndef MISCELLANEOUS_C
#define MISCELLANEOUS_C

#include "Navi.h"

//#define assert(b)
BOOL assert(BOOL b)
{
	int res;

	if(b)
	{
		res = 1;
	}
	else
	{
		res = 0;
	}
	return res;
}

// determine the sign of 'val' with the sensitivity of 'eps'
int signE(double val, double eps)
{
	int s;

	if(val<-eps)
	{
		s = -1;
	}
	else if(val>eps)
	{
		s = 1;
	}
	else
	{
		s = 0; 
	}
	return s;
}

// set double value 'val' between range 'minVal' and 'maxVal'
double range(double val, double minVal, double maxVal)
{
	double res;

	if(val<minVal)
	{ 
		res = minVal; 
	}
	else if(val>maxVal)	
	{ 
		res = maxVal; 
	}
	else				
	{ 
		res = val;
	}
	return res; 
}

int nextIndex(int curIdx, int maxIdx)
{
	int nIdx = curIdx + 1;

	if(nIdx>=maxIdx)
	{
		nIdx = 0;
	}
	return nIdx;
}

double asinEx(double x)
{
	return asin(range(x, -1.0, 1.0));
}

double atan2Ex(double y, double x)
{
	double res;

	if((sign(y)==0) && (sign(x)==0))
	{
		res = 0.0;
	}
	else
	{
		res = atan2(y, x);
	}
	return res;
}

#endif /*MISCELLANEOUS_C*/
