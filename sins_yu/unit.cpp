#include "navi.h"

double r2d(double rad)
{
	return rad/NAV_deg;
}

double d2r(double deg)
{
	return deg*NAV_deg;
}

double r2dm(double rad)
{
	double d, m;

	d = rad/NAV_deg;
	m = (d-(int)d)*60;
	return (int)d*100.0+m;
}

double dm2r(double deg)
{
	int d;
	double m;

	d = (int)(deg/100.0);
	m = deg-d*100.0; 
	return d*NAV_deg+m*NAV_min;
}

double r2dms(double rad)
{
	double d, m, s;

	d = rad/NAV_deg;
	m = (d-(int)d)*60; 
	s = (m-(int)m)*60;
	return (int)d*10000.0+(int)m*100.0+s;
}

double dms2r(double deg)
{
	int d, m;
	double s;

	d = (int)(deg/10000.0);
	m = (int)((deg-d*10000.0)/100.0); 
	s = deg-d*10000.0-m*100.0;
	return d*NAV_deg+m*NAV_min+s*NAV_sec;
}

double r2m(double rad)
{
	return rad/NAV_mil;
}

double m2r(double m)
{
	return m*NAV_mil;
}

