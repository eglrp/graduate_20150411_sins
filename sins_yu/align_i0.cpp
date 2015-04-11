#include "navi.h"

#define T_EXC 10

CAlign_i0::CAlign_i0(void)
{
}

CAlign_i0::CAlign_i0(CIMU &imu0, CVect3 &pos00, CEarth &eth0)
{
	imu = imu0; pos0 = pos00;
	eth = eth0; eth.Update(pos0);
	tk = 0;
	t0 = t1 = T_EXC, t2 = 0; 
	vib0 = vi0 = Pib01 = Pib02 = Pi01 = Pi02 = CVect3(0);
	qib0b = CQuat(1.0);
}

int CAlign_i0::Update(CVect3 &wm, CVect3 &vm)
{
	if(imu.Update(wm, vm)!=0)//¼ÆËãÐý×ªÊ¸Á¿¡¢Ô²×¶Îó²î²¹³¥Á¿
		return 0;	
	tk += imu.tss;
	double swiet = sin(tk*eth.wie), cwiet = cos(tk*eth.wie);
	CMat3 Cni0(-swiet,cwiet,0, 
                -eth.sl*cwiet,-eth.sl*swiet,eth.cl, 
                eth.cl*cwiet,eth.cl*swiet,eth.sl);
	vib0 = 0.9999*vib0 - qib0b*imu.dvbm;
	vi0 = 0.9999*vi0 + (~Cni0)*eth.gn*imu.tss;
	Pib02 = 0.999*Pib02 + vib0*imu.tss;
	Pi02 = 0.999*Pi02 + vi0*imu.tss;
	//
	if(++t2>3*t0)
		t0 = t1, Pib01 = tmpPib0, Pi01 = tmpPi0;
	else if(t2>2*t0 && t1==t0)
		t1 = t2, tmpPib0 = Pib02, tmpPi0 = Pi02;
	//
	qib0b *= RV2Q(imu.phim);
	// qnb=qni0*qiib0*qib0b
	if(t2<5*T_EXC)
		qnb = CQuat(1);
	else
		qnb = CQuat(Cni0)*CQuat(CMat3(Pi01, Pi02, Pib01, Pib02))*qib0b;
	
	return 1;
}

int CAlign_i0::IUpdate(CVect3 &wm, CVect3 &vm)
{
	eth.wie=-eth.wie;	//set negative
	int res = Update(-wm, vm);
	eth.wie=-eth.wie;	//restore
	return res;
}


