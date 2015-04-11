#include "navi.h"

CMat3::CMat3(void)
{
//	memset(&e00, 0, 9*sizeof(double));
	e00 = e01 = e02 = e10 = e11 = e12 = e20 = e21 = e22 = 0.0;
	return;
}

CMat3::CMat3(double val)
{
	e01 = e02 = e10 = e12 = e20 = e21 = 0.0;
	e00 = e11 = e22 = val;
	return;
}

CMat3::CMat3(double xx, double yy, double zz)
{
	e01 = e02 = e10 = e12 = e20 = e21 = 0.0;
	e00 = xx, e11 = yy, e22 = zz;
	return;
}

CMat3::CMat3(double xx, double xy, double xz, 
			 double yx, double yy, double yz,
			 double zx, double zy, double zz )
{
	e00 = xx;	e01 = xy;	e02 = xz;
	e10 = yx;	e11 = yy;	e12 = yz;
	e20 = zx;	e21 = zy;	e22 = zz;
	return;
}

CMat3::CMat3(CVect3 &att)
{
	double
		si = sin(att.i), ci = cos(att.i),
		sj = sin(att.j), cj = cos(att.j),
		sk = sin(att.k), ck = cos(att.k);
	e00 =  cj*ck - si*sj*sk;	e01 =  -ci*sk;	e02 = sj*ck + si*cj*sk;
	e10 =  cj*sk + si*sj*ck;	e11 =  ci*ck;	e12 = sj*sk - si*cj*ck;
	e20 = -ci*sj;				e21 =  si;		e22 = ci*cj;
	return;
}

CMat3::CMat3(CQuat &qnb)
{
	double 
		q11 = qnb.q0*qnb.q0, q12 = qnb.q0*qnb.q1, q13 = qnb.q0*qnb.q2, q14 = qnb.q0*qnb.q3, 
		q22 = qnb.q1*qnb.q1, q23 = qnb.q1*qnb.q2, q24 = qnb.q1*qnb.q3,     
		q33 = qnb.q2*qnb.q2, q34 = qnb.q2*qnb.q3,  
		q44 = qnb.q3*qnb.q3;
    e00 = q11+q22-q33-q44,  e01 = 2*(q23-q14),     e02 = 2*(q24+q13),
	e10 = 2*(q23+q14),      e11 = q11-q22+q33-q44, e12 = 2*(q34-q12),
	e20 = 2*(q24-q13),      e21 = 2*(q34+q12),     e22 = q11-q22-q33+q44 ;
	return;
}

CMat3::CMat3(CVect3 &v0, CVect3 &v1, CVect3 &v2)
{
	e00 = v0.i, e01 = v0.j, e02 = v0.k;
	e10 = v1.i, e11 = v1.j, e12 = v1.k;
	e20 = v2.i, e21 = v2.j, e22 = v2.k;
	return;
}

CMat3::CMat3(CVect3 &vn1, CVect3 &vn2, CVect3 &vb1, CVect3 &vb2)
{
	CVect3 vb1tmp, vb2tmp;
	vb1tmp = (!vn1)/(!vb1)*vb1;
	vb2tmp = (!vn2)/(!vb2)*vb2;
	CVect3 n=vn1*vn2, b=vb1tmp*vb2tmp;
	CMat3 Mn(vn1, n, n*vn1), Mb(vb1tmp, b, b*vb1tmp);
	*this = OptOrtho((Mn^-1)*Mb);	// Cnb
	return;
}

CMat3 CMat3::operator+(CMat3 &mat)
{
	CMat3 mtmp;
	double *p=&mtmp.e00, *p1=&e00, *p2=&mat.e00;
	for(int i=0; i<9; i++)
	{
		*p++ = (*p1++) + (*p2++);
	}
	return mtmp;
}

CMat3& CMat3::operator+=(CMat3 &mat)
{
	double *p1=&e00, *p2=&mat.e00;
	for(int i=0; i<9; i++)
	{
		(*p1++) += *p2++;
	}
	return *this;
}

CMat3 CMat3::operator-(CMat3 &mat)
{
	CMat3 mtmp;
	double *p=&mtmp.e00, *p1=&e00, *p2=&mat.e00;
	for(int i=0; i<9; i++)
	{
		*p++ = (*p1++) - (*p2++);
	}
	return mtmp;
}

CMat3& CMat3::operator-=(CMat3 &mat)
{
	double *p1=&e00, *p2=&mat.e00;
	for(int i=0; i<9; i++)
	{
		(*p1++) -= *p2++;
	}
	return *this;
}

CMat3 CMat3::operator*(CMat3 &mat)
{
	CMat3 mtmp;
	mtmp.e00 = e00*mat.e00 + e01*mat.e10 + e02*mat.e20;
	mtmp.e01 = e00*mat.e01 + e01*mat.e11 + e02*mat.e21;
	mtmp.e02 = e00*mat.e02 + e01*mat.e12 + e02*mat.e22;
	mtmp.e10 = e10*mat.e00 + e11*mat.e10 + e12*mat.e20;
	mtmp.e11 = e10*mat.e01 + e11*mat.e11 + e12*mat.e21;
	mtmp.e12 = e10*mat.e02 + e11*mat.e12 + e12*mat.e22;
	mtmp.e20 = e20*mat.e00 + e21*mat.e10 + e22*mat.e20;
	mtmp.e21 = e20*mat.e01 + e21*mat.e11 + e22*mat.e21;
	mtmp.e22 = e20*mat.e02 + e21*mat.e12 + e22*mat.e22;
	return mtmp;
}

CMat3& CMat3::operator*=(CMat3 &mat)
{
	return (*this=*this*mat);
}
	
CVect3 CMat3::operator*(CVect3 &v)
{
	return CVect3(
		e00*v.i + e01*v.j + e02*v.k,
		e10*v.i + e11*v.j + e12*v.k,
		e20*v.i + e21*v.j + e22*v.k
	);
}

CMat3 CMat3::operator*(double f)
{
	CMat3 mtmp;
	double *p1=&mtmp.e00, *p2=&e00;
	for(int i=0; i<9; i++)
	{
		*p1++ = (*p2++) * f;
	}
	return mtmp;
}

CMat3& CMat3::operator*=(double f)
{
	double *p2=&e00;
	for(int i=0; i<9; i++)
	{
		(*p2++) *= f;
	}
	return *this;
}

CMat3& CMat3::operator++()	// m += I
{ 
	e00 += 1.0, e11 += 1.0, e22 += 1.0;
	return *this;
}

CMat3 operator*(double f, CMat3 &mat)
{
	return mat*f;
}

CMat3 operator-(CMat3 &mat)
{
	CMat3 mtmp;
	double *p1=&mtmp.e00, *p2=&mat.e00;
	for(int i=0; i<9; i++)
	{
		*p1++ = -(*p2++);
	}
	return mtmp;
}

double& CMat3::operator()(int i, int j)
{
//	double *p=&this->e00;
//	return p[i*3+j];
	return (&e00)[i*3+j];
}

CMat3 operator ~(CMat3 &m)
{
/*	CMat3 mtmp;
	double *p1=&mtmp.e00, *p2=&m.e00;
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			*p1++ = p2[j*3+i];
		}
	}
	return mtmp;*/
	return CMat3(m.e00,m.e10,m.e20, m.e01,m.e11,m.e21, m.e02,m.e12,m.e22);
}

double operator!(CMat3 &m)
{
	return (  m.e00*(m.e11*m.e22-m.e12*m.e21)
			 -m.e01*(m.e10*m.e22-m.e12*m.e20)
			 +m.e02*(m.e10*m.e21-m.e11*m.e20) );
}

CMat3 CMat3::operator^(int n)
{
	CMat3 mtmp;
	if(n==0)
	{
		mtmp = CMat3(1.0);
	}
	else if(n<=-1)
	{
		double nm = !(*this);
//		if(!IsZero(nm))
		if(nm>1.0e-45||nm<-1.0e-45)
		{
			mtmp.e00 =  (this->e11*this->e22-this->e12*this->e21)/nm;
			mtmp.e10 = -(this->e10*this->e22-this->e12*this->e20)/nm;
			mtmp.e20 =  (this->e10*this->e21-this->e11*this->e20)/nm;
			mtmp.e01 = -(this->e01*this->e22-this->e02*this->e21)/nm;
			mtmp.e11 =  (this->e00*this->e22-this->e02*this->e20)/nm;
			mtmp.e21 = -(this->e00*this->e21-this->e01*this->e20)/nm;
			mtmp.e02 =  (this->e01*this->e12-this->e02*this->e11)/nm;
			mtmp.e12 = -(this->e00*this->e12-this->e02*this->e10)/nm;
			mtmp.e22 =  (this->e00*this->e11-this->e01*this->e10)/nm;
		}
		else
		{
			mtmp = CMat3(12345678);
		}
		n = -n;
	}
	else
	{
		mtmp = *this;
	}
	CMat3 res = mtmp;
	for(int i=1; i<n; i++)
	{
		res = res*mtmp;
	}
	return res;
}

CMat3 operator/(double f, CMat3 &m)
{
	return (m^-1)*f;
}

CMat3 OptOrtho(CMat3 &m)
{
	CMat3 mtmp = m;
	for(int i=0; i<5; i++)
	{
		mtmp = 1./2*(mtmp+((~mtmp)^-1));
	}
	return mtmp;
}

CMat3 ASkew(CVect3 &v)
{
	return CMat3(0, -v.k, v.j, 
				v.k, 0, -v.i,
				-v.j, v.i, 0);
}

CMat3 Diag(CVect3 &v)
{
	return CMat3(v.i,0,0, 0,v.j,0, 0,0,v.k);
}

CMat3 Cw(CVect3 &a)
{
	double
		si = sin(a.i), ci = cos(a.i),
		sj = sin(a.j), cj = cos(a.j);
	return CMat3(cj,0,-sj*ci, 0,1,si, sj,0,cj*ci);
}

CMat3 Cw_1(CVect3 &a)
{
	double
		si = sin(a.i), ci = cos(a.i),
		sj = sin(a.j), cj = cos(a.j);
	// if(IsZero(cj) error
	return CMat3(cj,0,sj, sj*si/ci,1,-cj*si/ci, -sj/ci,0,cj/ci);
}