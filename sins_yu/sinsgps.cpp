#include "navi.h"

#include "kfrapid.cpp"

CSINSGPS::CSINSGPS(CSINS &sins0, CVect3 &lever0)
{
	lever = lever0; eb = 0; db = 0;
	sins = sins0; //sins.pos=LeverPos(sins,-lever);
	kf = new CKalman(33,6,1.000001);	// 可不需遗忘因子
	pfi=(CVect3*)kf->fbXk.d, pvn=pfi+1, ppos=pfi+2, peb=pfi+3, pdb=pfi+4, 
		pkg1=pfi+5, pkg2=pfi+6, pkg3=pfi+7, pka1=pfi+8, pka23=pfi+9, plv=pfi+10;
	psfi=(CVect3*)kf->smXk.d, psvn=psfi+1, pspos=psfi+2, pseb=psfi+3, psdb=psfi+4, 
		pskg1=psfi+5, pskg2=psfi+6, pskg3=psfi+7, pska1=psfi+8, pska23=psfi+9, pslv=psfi+10;
	pzkv=(CVect3*)&kf->Zk.d[0];
	pzkp=(CVect3*)&kf->Zk.d[3];
#ifdef KF_RAPID
	kf->pAPAT = APAT;
#endif
	Init();
}

CSINSGPS::~CSINSGPS()
{
	if(kf) delete kf, kf=NULL;
}

int CSINSGPS::Update(CVect3 &gyro, CVect3 &acc, CVect3 &vng, CVect3 &gps, int obs)
{
	CVect3 wm = gyro-eb*sins.imu.ts, vm = acc-db*sins.imu.ts;
	int res=sins.Update(wm, vm);
	int obs_valid = obs;
	if(res==0)
	{
		double ts=sins.imu.tss;
		CMat3 S1, S2, M1, M2, M3, M4, M5, M6, Cnb=sins.Cnb;
		sins.ErrCoef(S1, S2, M1, M2, M3, M4, M5, M6);
		CMat3 Ct = Cnb*ts, Cwx=Ct*sins.wb.i, Cwy=Ct*sins.wb.j, Cwz=Ct*sins.wb.k, 
			Cfx=Ct*sins.fb.i, Cfy=Ct*sins.fb.j, Cfz=Ct*sins.fb.k;
		Cfy.e00=Cfz.e02, Cfy.e10=Cfz.e12, Cfy.e20=Cfz.e22;
		//  fi			dvn			dpos		eb		db		dKg(:,1)	dKg(:,2)	dKg(:,3)	dKa(:,1)	dKa(*,2:3)	lv
		kf->Fk.Clear();
		kf->Fk.Setm(
			&(S1*ts),	&(M1*ts),	&(M2*ts),	&-Ct,	&O33,	&Cwx,		&Cwy,		&Cwz,		&O33,		&O33,		&O33,		
			&(S2*ts),	&(M3*ts),	&(M4*ts),	&O33,	&Ct,	&O33,		&O33,		&O33,		&Cfx,		&Cfy,		&O33,
			&O33,		&(M5*ts),	&(M6*ts),	NULL );
		for(int k=0; k<kf->Fk.rc; k+=kf->Fk.clm+1) 
			{	kf->Fk.d[k] += 1.0; kf->Qk.d[k] = kf->Qt.d[k]*ts; }
		if(obs)
		{
//			posL=sins.LeverPos(lever); vnL=sins.LeverVn(lever);
//			kf->Hk.Setm(0,30, -Cnb*ASkew(sins.wb));	*pzkv = LeverVn(sins,lever) - vng;
//			kf->Hk.Setm(3,30, M5*Cnb);	*pzkp = LeverPos(sins,lever) - gps;
//			kf->Update(1);
//			obs_valid = 1;
			sins.vn = 0;

		}
		else
		{
			kf->Update(0);
		}
		kf->Smooth();
		FeedBack();
	}
	return obs_valid;
}

int CSINSGPS::IUpdate(CVect3 &gyro, CVect3 &acc, CVect3 &vng, CVect3 &gps, int obs)
{
	sins.eth.wie=-sins.eth.wie; eb=-eb; //sins.vn=-sins.vn;	//在静止时由正转逆，速度不取反，速度误差会更连续；
	int res=Update(-gyro, acc, -vng, gps, obs);					//但在运动中由正转逆，速度最好取反
	sins.eth.wie=-sins.eth.wie; eb=-eb; //sins.vn=-sins.vn;
	return res;
}

void CSINSGPS::Init(void)
{
	kf->Qt = Diag(3*NAV_dph,3*NAV_dph,3*NAV_dph, 5.0*NAV_mg,5.0*NAV_mg,5.0*NAV_mg, 
		0.01/NAV_Re,0.01/NAV_Re,0.01, 0.0,0.0,0.0, 0.0,0.0,0.0, 
		0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, lfEND)^2;
	kf->Rk = Diag(0.9, 0.9, 0.9, 0.10/NAV_Re, 0.10/NAV_Re, 0.10, lfEND)^2; 
	kf->Pk = Diag(1*NAV_deg,1*NAV_deg,1*NAV_deg, 0.10,0.10,0.10, 10.0/NAV_Re,10.0/NAV_Re, 10.0, 
		0.1*NAV_dph,0.1*NAV_dph,0.1*NAV_dph, 1000*NAV_ug,1000*NAV_ug,1000*NAV_ug,
		100*NAV_ppm,25*NAV_sec,25*NAV_sec, 25*NAV_sec,100*NAV_ppm,25*NAV_sec, 25*NAV_sec,25*NAV_sec,100*NAV_ppm, 
		100*NAV_ppm,25*NAV_sec,25*NAV_sec, 100*NAV_ppm,100*NAV_ppm,25*NAV_sec, 10.,10.,10., lfEND)^2; 	
	kf->Hk.Setm(0, 3, I33); 	kf->Hk.Setm(3, 6, I33);
	kf->InitSMCoef(0.02, 1.0,1.0,1.0, 0.0,0.0,0.0,  0.0,0.0,0.0, 
					1000.0,1000.0,1000.0, 1000.0,1000.0,1000.0,
					1000.0,1000.0,1000.0, 1000.0,1000.0,1000.0, 1000.0,1000.0,1000.0, 
					1000.0,1000.0,1000.0, 1000.0,1000.0,1000.0, 1000.0,1000.0,1000.0);
	InitFBCoef(2);
}

void CSINSGPS::InitFBCoef(int idx)
{
	switch(idx)
	{
	case 0:	kf->InitFBCoef(0.02, 10.0,10.0,10.0, 10.0,10.0,100.0, 10.0,10.0,100.0, 
						1000.0,1000.0,1000.0, 1000.0,1000.0,100.0,
						INF,INF,INF, INF,INF,INF, INF,INF,INF, 
						INF,INF,INF, INF,INF,INF, INF,INF,INF); break;
	case 1:	kf->InitFBCoef(0.02, 10.0,10.0,10.0, 10.0,10.0,10.0, 10.0,10.0,10.0,  //反馈姿态、速度、位置
						INF,INF,INF, INF,INF,INF,
						INF,INF,INF, INF,INF,INF, INF,INF,INF, 
						INF,INF,INF, INF,INF,INF, INF,INF,INF); break;
	case 2:	kf->InitFBCoef(0.02, INF,INF,INF, INF,INF,10.0, INF,INF,10.0, 	//反馈天向位置、速度
						INF,INF,INF, INF,INF,INF,
						INF,INF,INF, INF,INF,INF, INF,INF,INF, 
						INF,INF,INF, INF,INF,INF, INF,INF,INF); break;
	case 3:	kf->InitFBCoef(0.02, INF,INF,INF, INF,INF,INF, INF,INF,10.0,	//反馈天向位置
						INF,INF,INF, INF,INF,INF,
						INF,INF,INF, INF,INF,INF, INF,INF,INF, 
						INF,INF,INF, INF,INF,INF, INF,INF,INF); break;
	}
}

void CSINSGPS::FeedBack(void)  // 反馈校正
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

void CSINSGPS::Output(int isIMU)  // 输出校正
{
//	qnbL = sins.qnbk - *psfi*0;  CnbL = CMat3(qnbL); attL = CVect3(CnbL);
	attL = sins.att;
	if(isIMU)
	{
		vnL  = sins.vnk - *psvn;
		posL = sins.posk - *pspos;
	}
	else
	{
		vnL = LeverVn(sins, lever-*pslv) - *psvn;
		posL = LeverPos(sins, lever-*pslv) - *pspos;
//		vnL  = sins.vnk + CnbL*(sins.wb*(lever-*pslv)) - *psvn;
//		CVect3 Ln = CnbL*(lever-*pslv);
//		posL = sins.posk + CVect3(Ln.j/sins.eth.RMh, Ln.i/sins.eth.clRNh, Ln.k) - *pspos;
	}
}
