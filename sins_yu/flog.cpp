#include "navi.h"

CFlog::CFlog()
{
	f = NULL;
}

CFlog::CFlog(char *fname, char mode)
{
	fmode = mode;
	if(fmode=='t')
		f = fopen(fname, "wt+");
	else
		f = fopen(fname, "wb+");
}

CFlog::~CFlog()
{
	if(f) fclose(f), f=NULL;
}

CFlog& operator<<(CFlog& flog, char c)
{
	if(flog.fmode=='t')
		fprintf(flog.f, "%c", c); 
	else
		NULL;
	return flog;
}

CFlog& operator<<(CFlog& flog, char *str)
{
	fprintf(flog.f, "%s", str);  
	return flog;
}

CFlog& operator<<(CFlog& flog, int i)
{
	fprintf(flog.f, "%d ", i);  
	return flog;
}

CFlog& operator<<(CFlog& flog, unsigned int i)
{
	fprintf(flog.f, "%u ", i);  
	return flog;
}

CFlog& operator<<(CFlog& flog, long l)
{
	fprintf(flog.f, "%ld ", l);  
	return flog;
}

CFlog& operator<<(CFlog& flog, float f)
{
	fprintf(flog.f, "%.8e ", f);  
	return flog;
}

CFlog& operator<<(CFlog& flog, double d)
{
	if(flog.fmode=='t')
		fprintf(flog.f, "%.10e ", d);  
	else
		fwrite(&d, 1, sizeof(double), flog.f);
	return flog;
}

CFlog& operator<<(CFlog& flog, CVect3 &v)
{
	if(flog.fmode=='t')
		fprintf(flog.f, "%.10e %.10e %.10e ", v.i, v.j, v.k);
	else
		fwrite(&v, 3, sizeof(double), flog.f);
	return flog;
}

CFlog& operator<<(CFlog& flog, CQuat &q)
{
	fprintf(flog.f, "%.10e %.10e %.10e %.10e", q.q0, q.q1, q.q2, q.q3);  
	return flog;
}

CFlog& operator<<(CFlog& flog, CMat3 &m)
{
	double *pd = &m.e00;
	for(int i=0; i<3; i++, pd+=3)
		fprintf(flog.f, "%.10e %.10e %.10e\n", pd[0], pd[1], pd[2]);  
	return flog;
}

CFlog& operator<<(CFlog& flog, CMat &m)
{
	double *pd = m.d;
	if(flog.fmode=='t')
	{
		for(int i=0; i<m.row; i++)
		{
			for(int j=0; j<m.clm; j++)
				fprintf(flog.f, "%.10e ", *pd++);
			if(m.row>1)
				fprintf(flog.f, "\n");
		}
	}
	else
	{
		fwrite(pd, m.rc, sizeof(double), flog.f);
	}
	return flog;
}

CFlog& operator<<(CFlog& flog, CSINS &sins)
{
	if(flog.fmode=='t')
		fprintf(flog.f, "%8.6f %8.6f %9.6f %6.4f %6.4f %5.4f %12.10f %13.10f %6.4f ", 
			sins.att.i/NAV_deg,sins.att.j/NAV_deg,sins.att.k/NAV_deg,
			sins.vn.i,sins.vn.j,sins.vn.k, 
			sins.pos.i/NAV_deg,sins.pos.j/NAV_deg,sins.pos.k); 
	else
		flog<<sins.attk<<sins.vnk<<sins.posk;
	return flog;
}

CFlog& operator<<(CFlog& flog, CKalman &kf)
{
	flog<<"Fk=[\n"<<kf.Fk<<"];\n"
		<<"Hk=[\n"<<kf.Hk<<"];\n" 
		<<"Xk=[\n"<<kf.Xk<<"];\n" 
		<<"Pk=[\n"<<kf.Pk<<"];\n" 
		<<"Qk=[\n"<<kf.Qk<<"];\n" 
		<<"Rk=[\n"<<kf.Rk<<"];\n"
		<<"Zk=[\n"<<kf.Zk<<"];\n"
		<<"Pxz=[\n"<<kf.Pxz<<"];\n"
		<<"Pzz=[\n"<<kf.Pzz<<"];\n"
		<<"rk=[\n"<<kf.rk<<"];\n"
		<<"Kk=[\n"<<kf.Kk<<"];\n" ; 
	return flog;
}
