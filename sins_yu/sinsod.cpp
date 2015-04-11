#include "navi.h"

CSINSOD::CSINSOD(CSINS &sins0, CVect3 &inst0, CVect3 &lever0)
{
	sins = sins0; posDR = LeverPos(sins, lever0);
	lever = lever0; eb = 0; db = 0;
	kf = new CKalman(21,3,1.0);	// 可不需遗忘因子
	pfi=(CVect3*)kf->fbXk.d, pvn=pfi+1, ppos=pfi+2, peb=pfi+3, pdb=pfi+4, pinst=pfi+5, plv=pfi+6;
	psfi=(CVect3*)kf->smXk.d, psvn=psfi+1, pspos=psfi+2, pseb=psfi+3, psdb=psfi+4, psinst=psfi+5, pslv=psfi+6;
	*(CVect3*)&kf->Xk.d[18]=lever;
	pzkv=(CVect3*)&kf->Zk.d[0];
	SVins = SVod = 0; SMod = 0; SMlv = 0; ST = 0; dSm1 = 0;
	inst = inst0;
    prj = ODPrj(inst); prj1=prj;
#ifdef KF_RAPID
	kf->pAPAT = APAT;
#endif
	Init();
}

CVect3 CSINSOD::ODPrj(CVect3 &inst)
{
	double si=sin(inst.i), ci=cos(inst.i), sj=sin(inst.j), cj=cos(inst.j);
	return CVect3(sj*ci, cj*ci, si)*inst.k;
}

CSINSOD::~CSINSOD()
{
	if(kf) delete kf, kf=NULL;
}

int CSINSOD::Update(CVect3 &gyro, CVect3 &acc, double dSm)
{
	CVect3 wm = gyro-eb*sins.imu.ts, vm = acc-db*sins.imu.ts;
	int res=sins.Update(wm, vm);
	posDR +=  vn2dpos(sins.eth, sins.Cnb*prj1*dSm);	//DR
	static runTime=0;
//	cout<<runTime++<<"	"<<
	int obs_valid = 0;
	dSm1 += dSm;
	if(res==0)  //2013. 7.10 AM  HYQ comment
//	if(0)//skip filtering
	{ 
		double ts=sins.imu.tss;
		CMat3 S1, S2, M1, M2, M3, M4, M5, M6, Cnb=sins.Cnb;
		sins.ErrCoef(S1, S2, M1, M2, M3, M4, M5, M6);
		CMat3 Ct = Cnb*ts;
		//  fi			dvn			dpos		eb		db		dins	lv
		kf->Fk.Clear();
		kf->Fk.Setm(
			&(S1*ts),	&(M1*ts),	&(M2*ts),	&-Ct,	&O33,	&O33,	&O33,
			&(S2*ts),	&(M3*ts),	&(M4*ts),	&O33,	&Ct,	&O33,	&O33,
			&O33,		&(M5*ts),	&(M6*ts),	NULL );
		for(int k=0; k<kf->Fk.rc; k+=kf->Fk.clm+1) 
			{	kf->Fk.d[k] += 1.0; kf->Qk.d[k] = kf->Qt.d[k]*ts; }
		SVins += sins.vn*ts;
		CVect3 dSb=prj*dSm1; dSm1=0;
		SVod += Cnb*dSb; SMod += Cnb*ASkew(dSb); SMlv += Cnb*ASkew(sins.wb*ts); ST += ts;
		if(ST>=1.00)
		{
			SMod.e01 = SVod.i, SMod.e11 = SVod.j, SMod.e21 = SVod.k;
			kf->Hk.Setm(&ASkew(-SVod), &(I33*ST), &O33, &O33, &O33, &(-SMod), &(-SMlv));	
			*pzkv = SVins - SVod;
			kf->Update(1);
				CVect3 p=inst; p.i -= psinst->i; p.j -= psinst->j; p.k *= 1-psinst->k;  
				prj1=ODPrj(p);
			if(sins.imu.tt<600)
			{
				posDR = LeverPos(sins, *pslv, &(sins.pos-*pspos));
			}
			obs_valid = 1;
			SVins = SVod = 0; SMod = 0; SMlv = 0; ST = 0;
//			if(sins.imu.tt>543)  { CFlog fg("D:\\ygm\\立得空间\\跑车数据1217\\给严老师测试数据\\Ft.m", 't'); fg<<*kf; exit(0); }
//			if(sins.imu.tt>10) { APATBuild("scr\\kfrapid.cpp", kf->Fk); exit(0); }
		}
		else
		{
			kf->Update(0);
		}
		kf->Smooth();
//		FeedBack();
	}
	
	return obs_valid;
}

int CSINSOD::IUpdate(CVect3 &gyro, CVect3 &acc, double dSm)
{
	sins.eth.wie=-sins.eth.wie; eb=-eb; 
	int res=Update(-gyro, acc, -dSm);
	sins.eth.wie=-sins.eth.wie; eb=-eb;
	return res;
}

void CSINSOD::Init(void)
{
	kf->Qt = Diag(3*NAV_dph,3*NAV_dph,3*NAV_dph, 5.0*NAV_mg,5.0*NAV_mg,5.0*NAV_mg, 
		0.01/NAV_Re,0.01/NAV_Re,0.01, 0.0,0.0,0.0, 0.0,0.0,0.0, 
		0.0,0.0,0.0, 0.001,0.001,0.001, lfEND)^2;
	kf->Rk = Diag(.91, .91, .91, lfEND)^2; 
	kf->Pk = Diag(1*NAV_deg,1*NAV_deg,1*NAV_deg, 1.0,1.0,1.0, 10.0/NAV_Re,10.0/NAV_Re, 10.0, 
		0.1*NAV_dph,0.1*NAV_dph,0.1*NAV_dph, 1000*NAV_ug,1000*NAV_ug,1000*NAV_ug,
		1.0*NAV_deg,0.1,1.0*NAV_deg, 10.,10.,10., lfEND)^2; 	
	kf->InitSMCoef(0.02, 1.0,1.0,1.0, 0.0,0.0,0.0,  0.0,0.0,0.0, 
					1000.0,1000.0,1000.0, 1000.0,1000.0,1000.0,
					10.0,10.0,10.0, 10.0,10.0,10.0);
	InitFBCoef(3);
}

void CSINSOD::InitFBCoef(int idx)
{
	switch(idx)
	{
	case 0:	kf->InitFBCoef(0.02, 10.0,10.0,10.0, 10.0,10.0,100.0, 10.0,10.0,100.0, 
						1000.0,1000.0,1000.0, 1000.0,1000.0,100.0,
						INF,INF,INF, INF,INF,INF); break;
	case 1:	kf->InitFBCoef(0.02, 10.0,10.0,10.0, 10.0,10.0,10.0, 10.0,10.0,10.0,  //反馈姿态、速度、位置
						INF,INF,INF, INF,INF,INF,
						INF,INF,INF, INF,INF,INF); break;
	case 2:	kf->InitFBCoef(0.02, INF,INF,INF, INF,INF,10.0, INF,INF,10.0, 	//反馈天向位置、速度
						INF,INF,INF, INF,INF,INF,
						INF,INF,INF, INF,INF,INF); break;
	case 3:	kf->InitFBCoef(0.02, INF,INF,INF, INF,INF,INF, INF,INF,10.0,	//反馈天向位置
						INF,INF,INF, INF,INF,INF,
						INF,INF,INF, INF,INF,INF); break;
	}
}

void CSINSOD::FeedBack(void)  // 反馈校正
{
	kf->FeedBack();
	sins.qnb-= *pfi; sins.Cnb = CMat3(sins.qnb); sins.att = CVect3(sins.Cnb);
	sins.vn	-= *pvn;
	sins.pos-= *ppos;
	eb		+= *peb;
	db		+= *pdb;		
	lever	-= *plv;
	sins.attk = sins.att, sins.Cnbk = sins.Cnb, sins.qnbk = sins.qnb, 
		sins.vnk = sins.vn, sins.posk = sins.pos;
}

void CSINSOD::Output(CVect3 &lv)  // 输出校正
{
	vnL = LeverVn(sins, lv) - *psvn;
	posL = LeverPos(sins, lv) - *pspos;
}
