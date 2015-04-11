#include "navi.h"

CEarth::CEarth()
{
}

CEarth::CEarth(double a0, double f0, double wie0)
{
	a = a0;	f = f0; wie = wie0; 
	b = (1-f)*a;
	e = sqrt(a*a-b*b)/a;	e2 = e*e;
	ep = sqrt(a*a-b*b)/b;	ep2 = ep*ep;
}

void CEarth::Update(CVect3 &pos, CVect3 &vn)
{
	sl = sin(pos.LTI), cl = cos(pos.LTI), tl = sl/cl;
	double sq = 1-e2*sl*sl, sq2 = sqrt(sq);
	RMh = a*(1-e2)/sq/sq2+pos.HGT;	f_RMh = 1.0/RMh;
	RNh = a/sq2+pos.HGT;    clRNh = cl*RNh;  f_RNh = 1.0/RNh; f_clRNh = 1.0/clRNh;
	wnie.i = 0,				wnie.j = wie*cl,		wnie.k = wie*sl;
	wnen.i = -vn.NTH*f_RMh,	wnen.j = vn.EST*f_RNh,	wnen.k = wnen.j*tl;//wnen.j = vn.EST*f_RNh ”¶∏√ « vn.EST*f_clRNh ∞…
	wnin = wnie + wnen;
	sl2 = sl*sl, sl4 = sl2*sl2;
	double g = NAV_g0*(1+5.27094e-3*sl2+2.32718e-5*sl4)-3.086e-6*pos.HGT;
	gn = CVect3(0, 0, -g);
}

CVect3 vn2dpos(CEarth &eth, CVect3 &vn, double ts)
{
	return CVect3(vn.NTH*eth.f_RMh*ts, vn.EST*eth.f_clRNh*ts, vn.UP*ts);
}
