#include "navi.h"

CVect3::CVect3(void)
{
	return;
}

CVect3::CVect3(double val)
{
	i = j = k = val;
	return;
}

CVect3::CVect3(double xx, double yy, double zz)
{
	i = xx, j = yy, k = zz;
	return;
}

CVect3::CVect3(CMat3 &Cnb)
{
	*this = CVect3( asinEx(Cnb.e21),		/* pitch: -pi/2 -- pi/2 */
		atan2Ex(-Cnb.e20, Cnb.e22),			/* roll: -pi -- pi */
		atan2Ex(-Cnb.e01, Cnb.e11) );		/* azimuth: -pi -- pi */
	return;
}

CVect3::CVect3(CQuat &qnb)
{
	*this = CVect3(CMat3(qnb));
	return;
}

CVect3 CVect3::operator+(CVect3 &v)
{ 
	return CVect3(i+v.i, j+v.j, k+v.k); 
}

CVect3 CVect3::operator-(CVect3 &v)
{ 
	return CVect3(i-v.i, j-v.j, k-v.k);
}

CVect3& CVect3::operator+=(CVect3 &v)
{ 
	i += v.i, j += v.j, k += v.k;
	return *this;
}

CVect3& CVect3::operator-=(CVect3 &v)
{ 
	i -= v.i, j -= v.j, k -= v.k;
	return *this;
}

CVect3 CVect3::operator*(CVect3 &v)
{ 
	return CVect3(j*v.k-k*v.j, k*v.i-i*v.k,	i*v.j-j*v.i);
}

CVect3 CVect3::operator*(double f)
{ 
	return CVect3(i*f, j*f, k*f);
}

CVect3 CVect3::operator*(CMat3 &m)
{
	return CVect3(
		i*m.e00 + j*m.e10 + k*m.e20,
		i*m.e01 + j*m.e11 + k*m.e21,
		i*m.e02 + j*m.e12 + k*m.e22
	);
}

CVect3 operator*(double f, CVect3 &v)
{ 
	return CVect3(v.i*f, v.j*f, v.k*f);
}

CVect3& CVect3::operator*=(double f)
{ 
	i *= f, j *= f, k *= f;
	return *this;
}

CVect3 CVect3::operator/(double f)
{ 
	return CVect3(i/f, j/f, k/f);
}

CVect3& CVect3::operator/=(double f)
{
	assert( sign(f)!=0 );
	i /= f, j /= f, k /= f;
	return *this;
}

CVect3 operator-(CVect3 &v)
{ 
	return CVect3(-v.i, -v.j, -v.k);
}

double& CVect3::operator()(int i)
{
	double *pd = &this->i;
	return pd[i];
}

double operator!(CVect3 &v)
{
	return sqrt(v.i*v.i + v.j*v.j + v.k*v.k);
}

double VDot(CVect3 v1, CVect3 &v2)
{
	return (v1.i*v2.i + v1.j*v2.j + v1.k*v2.k);
}

/*
CQuat RV2Q(CVect3 &v)
{
	CQuat qtmp;
	double v_2 = (!v)/2.0;
	if(fabs(v_2)<1.0e-20)
	{
		qtmp.q0 = 1, qtmp.q1 = qtmp.q2 = qtmp.q3 = 0;
	}
	else
	{
		double tmp = sin(v_2)/v_2*0.5;
		qtmp.q0 = cos(v_2), 
			qtmp.q1 = tmp*v.i, qtmp.q2 = tmp*v.j, qtmp.q3 = tmp*v.k;
	}
	return qtmp;
}
*/
/*
CQuat RV2Q(CVect3 &v)
{
	double v_2 = (!v)/2.0, f;
	if(v_2<1.0e-4)
	{
		double v2=v_2*v_2;
		f = (1.0 - v2/6.0 + v2*v2/120.0)*0.5;
	}
	else
	{
		f = sin(v_2)/v_2*0.5;
	}
	return CQuat(cos(v_2), f*v.i, f*v.j, f*v.k);
}
*/

CQuat RV2Q(CVect3 &v)
{
#define F1	(   2 * 1)		// define: Fk=2^k*k! 
#define F2	(F1*2 * 2)
#define F3	(F2*2 * 3)
#define F4	(F3*2 * 4)
#define F5	(F3*2 * 5)
	double n2 = v.i*v.i+v.j*v.j+v.k*v.k, c, f;
	if(n2<NAV_deg*NAV_deg)	// 0.017^2 
	{
		double n4=n2*n2;
		c = 1.0 - n2*(1.0/F2) + n4*(1.0/F4);
		f = 0.5 - n2*(1.0/F3) + n4*(1.0/F5);
	}
	else
	{
		double n_2 = sqrt(n2)/2.0;
		c = cos(n_2);
		f = sin(n_2)/n_2*0.5;
	}
	return CQuat(c, f*v.i, f*v.j, f*v.k);
}

CMat3 P2Cen(CVect3 &v)
{
	double si = sin(v.i), ci = cos(v.i), 
		sj = sin(v.j), cj = cos(v.j);
	return CMat3(
		-sj, -si*cj,  ci*cj,  
		 cj, -si*sj,  ci*sj,  
		 0,   ci,     si      );	//Cen
}
