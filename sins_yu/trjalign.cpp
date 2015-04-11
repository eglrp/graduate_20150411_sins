#include "navi.h"

CTrjAlign::CTrjAlign(double ts0, double dyaw0, CVect3 &pos0, CEarth &eth0)
{
	firstVng = 0;
	ts = ts0; dyaw = dyaw0; psi = gpspsi = 0.0;
	eth = eth0; eth.Update(pos0); pos = pos0;
	att = vn = 0;	Cnb = I33;
	qnb = Q1; //CQuat(CVect3(10,10,20)*NAV_deg);
	kf = new CKalman(6,3,1.001);
	kf->Fk = eye(6);
	kf->Qk = (Diag(10*NAV_dph,10*NAV_dph,1*NAV_dph, 10*NAV_mg,10*NAV_mg,10*NAV_mg, lfEND)^2)*ts;
	kf->Rk = Diag(.50, .50, 1.00, lfEND)^2;
	kf->Pk = Diag(10*NAV_deg,10*NAV_deg,10*NAV_deg, 10.0,10.0,10.0, lfEND)^2;
	kf->Hk(0,3) = kf->Hk(1,4) = kf->Hk(2,5) = 1.0; 
	kf->InitFBCoef(ts, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0);
}

int CTrjAlign::Update(CVect3 &wm, CVect3 &vm, CVect3 &vng, int vng_valid)
{
	if(vng_valid && (vng.i*vng.i+vng.j*vng.j)>1.0)
	{
		gpspsi = atan2(-vng.i,vng.j)+dyaw;
		if(firstVng==0) 
		{
			qnb = CVect3(att.i, att.j, gpspsi);	//第一次用GPS方位赋值
			vn  = vng;
			firstVng = 1;
		}
	}
	if(firstVng)
	{
		CVect3 dfn = Cnb*vm; 
		kf->Fk(3,1)=-dfn.k; kf->Fk(3,2)=dfn.j; kf->Fk(4,0)=dfn.k; kf->Fk(4,2)=-dfn.i;
		vn  += dfn+eth.gn*ts;
		qnb *= RV2Q(wm-(~Cnb)*eth.wnin*ts);  Cnb = CMat3(qnb); att = CVect3(Cnb);
		if(vng_valid)
		{
			eth.Update(pos, vng);
			*(CVect3*)kf->Zk.d = vn - vng; // moving align observ
		//	*(CVect3*)kf->Zk.d = vn - vng*0;//HYQ modified, static align observ
//			kf->Zk.d[0] = vn.i-vng.i, kf->Zk.d[1] = vn.j-vng.j, kf->Zk.d[2] = vn.k-vng.k;
			kf->Update(1);
			firstVng ++;
		}
		else
		{
			kf->Update(0);
		}
		kf->FeedBack();
		qnb -= *(CVect3*)&kf->fbXk.d[0]; //CVect3(kf->fbXk.d[0], kf->fbXk.d[1], kf->fbXk.d[2]);
		vn  -= *(CVect3*)&kf->fbXk.d[3]; //CVect3(kf->fbXk.d[3], kf->fbXk.d[4], kf->fbXk.d[5]);
	}
	return firstVng;
}

int CTrjAlign::IUpdate(CVect3 &wm, CVect3 &vm, CVect3 &vng, int trjvalid)
{
	eth.wie=-eth.wie;	//set negative
	int res = Update(-wm, vm, -vng, trjvalid);
	eth.wie=-eth.wie;	//restore
	return res;
}


