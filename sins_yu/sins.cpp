#include "navi.h"

CSINS::CSINS(void)
{
}

CSINS::CSINS(CIMU &imu0, CQuat &qnb0, CVect3 &vn0, CVect3 &pos0, CEarth &eth0)
{
	imu = imu0; qnb = qnb0; Cnb = CMat3(qnb); att = CVect3(qnb);
	vn = vn0, pos = pos0;
	attk = att, Cnbk = Cnb, qnbk = qnb, vnk = vn, posk = pos;
	eth = eth0;		eth.Update(pos,vn);
	wb = fb = fn = 0;
}

int CSINS::Update(CVect3 &wm, CVect3 &vm)
{
	imu.Update(wm, vm);
	if(imu.ksample)   // 多子样中间插值，简化版导航
	{
		double kts = imu.ksample*imu.ts;//导航解算周期
		vnk = vn + Cnb*imu.vmm + eth.gn*kts;//速度更新
		posk = pos + vn2dpos(eth, vnk, kts);//位置更新
		qnbk = qnb * RV2Q(imu.wmm);  Cnbk = CMat3(qnbk); attk = CVect3(Cnbk);
	}
	else
	{
		double kts = imu.tss;
		eth.Update(pos, vn);
		CVect3 wbin = (~Cnb)*eth.wnin;
		wb = imu.phim / kts - wbin;	fb = imu.dvbm / kts;	fn = Cnb*fb;
		CVect3 vn1 = vn + (~RV2Q(eth.wnin*(kts/2)))*(Cnb*imu.dvbm) 
			+ (eth.gn-(eth.wnie+eth.wnin)*vn)*kts;
		pos += vn2dpos(eth, vn+vn1, kts/2);	vn = vn1;
		qnb *= RV2Q(wb*kts);	Cnb = CMat3(qnb);	att = CVect3(Cnb);
		attk = att, Cnbk = Cnb; qnbk = qnb, vnk = vn, posk = pos;
	}
	return imu.ksample;
}

int CSINS::IUpdate(CVect3 &wm, CVect3 &vm)
{
	eth.wie=-eth.wie; vn=-vn; //setting reverse
	int k = Update(-wm, vm);
	eth.wie=-eth.wie; vn=-vn; //restore
	return k;
}

CVect3 LeverPos(CSINS &sins, CVect3 &lever, CVect3 *p0)  // 杆臂位置
{
	CVect3 Ln = sins.Cnbk*(lever);
	CVect3 posL = (p0==NULL?sins.posk:*p0) + CVect3(Ln.j/sins.eth.RMh, Ln.i/sins.eth.clRNh, Ln.k);
	return posL;
}

CVect3 LeverVn(CSINS &sins, CVect3 &lever, CVect3 *v0)  // 杆臂速度
{
	CVect3 vnL = (v0==NULL?sins.vnk:*v0) + sins.Cnbk*(sins.wb*lever); 
	return vnL;
}

void CSINS::ErrCoef(CMat3 &S1, CMat3 &S2, 
					CMat3 &M1, CMat3 &M2, CMat3 &M3, CMat3 &M4, CMat3 &M5, CMat3 &M6)
{
	double tl=eth.tl, secl=1.0/eth.cl, secl2=secl*secl, wN=eth.wnie.NTH, 
		wU=eth.wnie.UP, vE=vn.EST, vN=vn.NTH;
	double f_RMh=eth.f_RMh, f_RNh=eth.f_RNh, f_clRNh=eth.f_clRNh, 
		f_RMh2=f_RMh*f_RMh, f_RNh2=f_RNh*f_RNh;
	CMat3 avn=ASkew(vn),
		Mp(0,0,0, -wU,0,0, wN,0,0),
		Mpp(0,0,vN*f_RMh2, 0,0,-vE*f_RNh2, vE*secl2*f_RNh,0,-vE*tl*f_RNh2);
	M1 = CMat3(0,-f_RMh,0, f_RNh,0,0, tl*f_RNh,0,0);
	M2 = Mp+Mpp;
	M3 = avn*M1 - ASkew(eth.wnie+eth.wnin);
	M4 = avn*(Mp+M2);
	M5 = CMat3(0,f_RMh,0, f_clRNh,0,0, 0,0,1);
	M6 = CMat3(0,0,-vN*f_RMh2, vE*tl*f_clRNh,0,-vE*secl*f_RNh2, 0,0,0);
	S1 = -ASkew(eth.wnin);
	S2 = ASkew(fn);
}