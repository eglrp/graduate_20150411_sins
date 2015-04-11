#include "navi.h"

CQuat::CQuat(void)
{
	q0 = 1.0, q1 = q2 = q3 = 0.0;
	return;
}

CQuat::CQuat(double val)
{
	q0 = val, q1 = q2 = q3 = 0.0;
	return;
}

CQuat::CQuat(double pch, double rll, double yaw)
{
	*this = CQuat(CVect3(pch,rll,yaw));
	return;
}

CQuat::CQuat(double qq0, double qq1, double qq2, double qq3)
{
	q0 = qq0, q1 = qq1, q2 = qq2, q3 = qq3;
	return;
}

CQuat::CQuat(CVect3 &att)
{
	*this = CQuat(CMat3(att));
	return;
}

CQuat::CQuat(CMat3 &Cnb)
{
	double tmp;

	tmp = 1.0 + Cnb.e00 - Cnb.e11 - Cnb.e22;
	q1 = sqrt(fabs(tmp))/2.0;
	tmp = 1.0 - Cnb.e00 + Cnb.e11 - Cnb.e22;
	q2 = sqrt(fabs(tmp))/2.0;
	tmp = 1.0 - Cnb.e00 - Cnb.e11 + Cnb.e22;
	q3 = sqrt(fabs(tmp))/2.0;
	tmp = 1.0 - q1*q1 - q2*q2 - q3*q3;
	q0 = sqrt(fabs(tmp));

	if(Cnb.e21 - Cnb.e12 < 0)	/* sign decision */
	{
		q1 = -q1;
	}
	if(Cnb.e02 - Cnb.e20 < 0)
	{
		q2 = -q2;
	}
	if(Cnb.e10 - Cnb.e01 < 0)
	{
		q3 = -q3;
	}

	*this = Normlize(*this);
}

CQuat CQuat::operator+(CVect3 &v)
{
	return RV2Q(CVect3(0.0)-v)*(*this);
//	return RV2Q(-v)*(*this);  // why NO in bc31!
}

CQuat& CQuat::operator+=(CVect3 &v)
{
	return *this=RV2Q(CVect3(0.0)-v)*(*this);
//	return RV2Q(-v)*(*this);  // why NO in bc31!
}

CQuat CQuat::operator-(CVect3 &v)
{
	return RV2Q(v)*(*this);
}

CQuat& CQuat::operator-=(CVect3 &v)
{
	return *this=RV2Q(v)*(*this);
}

CVect3 CQuat::operator-(CQuat &quat)
{
	CQuat dq;
	
	dq = quat*~(*this);
	if(dq.q0<0)
	{
		dq.q0=-dq.q0, dq.q1=-dq.q1, dq.q2=-dq.q2, dq.q3=-dq.q3;
	}
	double n2 = acos(dq.q0), f;
	if( sign(n2)!=0 )
	{
		f = 2.0/(sin(n2)/n2);
	}
	else
	{
		f = 2.0;
	}
	return CVect3(dq.q1,dq.q2,dq.q3)*f;
}

CQuat CQuat::operator*(CQuat &quat)
{
	CQuat qtmp;
	qtmp.q0 = q0*quat.q0 - q1*quat.q1 - q2*quat.q2 - q3*quat.q3;
	qtmp.q1 = q0*quat.q1 + q1*quat.q0 + q2*quat.q3 - q3*quat.q2;
	qtmp.q2 = q0*quat.q2 + q2*quat.q0 + q3*quat.q1 - q1*quat.q3;
	qtmp.q3 = q0*quat.q3 + q3*quat.q0 + q1*quat.q2 - q2*quat.q1;
	return qtmp;
}

CQuat& CQuat::operator*=(CQuat &quat)
{
	return (*this=*this*quat);
}

CVect3 CQuat::operator*(CVect3 &v)
{
	CQuat qtmp;
	CVect3 vtmp;
	qtmp.q0 =         - q1*v.i - q2*v.j - q3*v.k;
	qtmp.q1 = q0*v.i           + q2*v.k - q3*v.j;
	qtmp.q2 = q0*v.j           + q3*v.i - q1*v.k;
	qtmp.q3 = q0*v.k           + q1*v.j - q2*v.i;
	vtmp.i = -qtmp.q0*q1 + qtmp.q1*q0 - qtmp.q2*q3 + qtmp.q3*q2;
	vtmp.j = -qtmp.q0*q2 + qtmp.q2*q0 - qtmp.q3*q1 + qtmp.q1*q3;
	vtmp.k = -qtmp.q0*q3 + qtmp.q3*q0 - qtmp.q1*q2 + qtmp.q2*q1;
	return vtmp;
}

CQuat CQuat::operator*(double f)
{ 
	return CQuat(q0*f, q1*f, q2*f, q3*f);
}

double& CQuat::operator()(int i)
{
	double *pd = &this->q0;
	return pd[i];
}

CVect3 Q2RV(CQuat &q)
{
	CQuat dq;
	dq = q;
	if(dq.q0<0)
	{
		dq.q0=-dq.q0, dq.q1=-dq.q1, dq.q2=-dq.q2, dq.q3=-dq.q3;
	}
	double n2 = acos(dq.q0), f;
	if( sign(n2)==0 )
	{
		f = 2.0/(sin(n2)/n2);
	}
	else
	{
		f = 2.0;
	}
	return CVect3(dq.q1,dq.q2,dq.q3)*f;
}

CQuat operator~(CQuat &q)
{
	return CQuat(q.q0, -q.q1, -q.q2, -q.q3);
}

double operator!(CQuat &q)
{
	return sqrt(q.q0*q.q0+q.q1*q.q1+q.q2*q.q2+q.q3*q.q3);
}

CQuat Normlize(CQuat &q)
{
	double tmp=!q;
	if( sign(tmp)==0 )
	{
		return CQuat(1.0);
	}
	else
	{
		return CQuat(q.q0/tmp, q.q1/tmp, q.q2/tmp, q.q3/tmp);
	}
}
