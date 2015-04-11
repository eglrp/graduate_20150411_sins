#include "navi.h"

CKalman::CKalman(int q0, int r0, double forgetting0)
{
	q = q0, r = r0;
	Xk=Xkk_1=fbXk=fbCoef=smXk=smCoef=CMat(q,1);
	Zk=rk=CMat(r,1);
	Ft=Fk=Qt=Qk=Gt=Pk=Iq=Pkk_1=lamda=CMat(q,q);
	Pxz=Kk=CMat(q,r);
	Pzz=Pzz_1=Rk=CMat(r,r);
	Hk=CMat(r,q);
	pd = MatAlloc(&Xk, &fbXk, &fbCoef, &smXk, &smCoef, &Zk, &rk, &Ft, &Fk, &Qt, &Qk, &Gt, &Pxz, &Pzz, &Pzz_1, 
		&Pk, &Hk, &Rk, &Iq, &Kk, &Xkk_1, &Pkk_1, &lamda, NULL);
	Iq = Gt = eye(q);
	forgetting = forgetting0;
	pAPAT = NULL;
}

CKalman::~CKalman()
{
	if(pd)	delete pd, pd = NULL;
}

void CKalman::Discrete(double ts)
{
	Fk = Iq + Ft*ts;
	Qk = Qt*ts;
}

void CKalman::TimeUpdate(void)
{
	if(pAPAT)
	{
		pAPAT(Fk, Pk, Pkk_1, Xk, Xkk_1);
		Pkk_1 = Pkk_1 + Qk;
	}
	else
	{
		Xkk_1 = Fk * Xk;
//		Pkk_1 = Fk*Pk*(~Fk) + Qk;
		Pkk_1 = Fk.SparseMul1(~Fk.SparseMul(Pk)) + Qk;
	}
}

void CKalman::MeasureUpdate(void)
{
	Pxz = Pkk_1*(~Hk);
	Pzz = Hk*Pxz + Rk;
	Pzz_1 = Pzz^-1;
	rk = Zk - Hk*Xkk_1;
	Kk = Pxz*Pzz_1;
	Xk = Xkk_1 + Kk*rk;
	Pk = Pkk_1 - Kk*(~Pxz);	symmetric(Pk);
	if(forgetting>1.0)
		Pk = Pk*forgetting;
}

void CKalman::Update(int obs)
{
	TimeUpdate();
	if(obs)
	{
		MeasureUpdate();
	}
	else
	{
		Xk = Xkk_1;
		Pk = Pkk_1;
	}
}

void CKalman::InitFBCoef(double ts, double coef, ...)
{
	double *pd = fbCoef.d;
	*pd++ = coef;
	va_list vl;
	va_start(vl, coef);
	for(int i=1; i<fbCoef.rc; i++)
		*pd++ = va_arg(vl, double);
	va_end(vl);
	for(i=0,pd=fbCoef.d; i<fbCoef.rc; i++,pd++)
	{
		*pd = *pd<ts ? 0.0 : exp(-ts/(*pd));
	}
}

void CKalman::FeedBack(void)
{
	fbXk = (1.0-fbCoef) & Xk;
	Xk =  fbCoef & Xk;
}

void CKalman::InitSMCoef(double ts, double coef, ...)
{
	double *pd = smCoef.d;
	*pd++ = coef;
	va_list vl;
	va_start(vl, coef);
	for(int i=1; i<smCoef.rc; i++)
		*pd++ = va_arg(vl, double);
	va_end(vl);
	for(i=0,pd=smCoef.d; i<smCoef.rc; i++,pd++)
	{
		*pd = *pd<ts ? 0.0 : exp(-ts/(*pd));
	}
	smXk = 0.0;
}

void CKalman::Smooth(void)
{
	smXk = (smCoef & smXk) + ((1-smCoef) & Xk);
}
