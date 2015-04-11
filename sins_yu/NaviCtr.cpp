#include "navi.h"
extern ifstream fRawIMUData;
extern ofstream fNavInfo,fAlignInfo;
extern CVect3 StartPos;
CNaviCtr::CNaviCtr(CVect3 &eb0, CVect3 &db0, CVect3 &lgps0, CVect3 &lod0, CVect3 &inst0)
{
	eb = eb0, db = db0; lgps = lgps0; lod = lod0; inst = inst0;
	step = 0;
	timu = tgps = gpscoef = 0;
	wm = vm = vng = posg = 0;
	pIMU = pOD = pGPS = NULL; 
	plogAlign = plogNavi = NULL;
	pSG = NULL;
	pSOD = NULL;
}

CNaviCtr::~CNaviCtr()
{
}

void CNaviCtr::FileOpen(char *path, char *imu, char *od, char *gps, char *align, char *navres)
{
	char str[256], path1[256];
	strcat(strcpy(path1, path), "\\");
	/* HYQ修改
	if(imu) pIMU = new CBinFile(strcat(strcpy(str, path1), imu), 14);
	if(od)  pOD =new CBinFile(strcat(strcpy(str, path1), od), 1);
	if(gps) pGPS = new CBinFile(strcat(strcpy(str, path1), gps), 15);	*/
	if(align) plogAlign = new CFlog(strcat(strcpy(str, path1), align), 'b');	
	if(navres) plogNavi = new CFlog(strcat(strcpy(str, path1), navres), 'b');	

	if(imu)
	{
		fRawIMUData.open(strcat(strcpy(str, path1), imu));
		fAlignInfo.open(strcat(strcpy(str, path1),align));
		fNavInfo.open(strcat(strcpy(str, path1),navres));
	}

}

int CNaviCtr::StaticAlign(int tStart, int tEnd)
{
	ReadData(tStart*FRQ200);
	posg = StartPos;
	
	CAlign_i0 align(CIMU(TS5,5), posg);
	if(tStart<tEnd)
	{
		for(step=tStart*FRQ200+1; step<=tEnd*FRQ200; step++)
		{
			ReadData(step);
			align.Update(wm, vm);	
			if((step%200==0)&&plogAlign) 
			{
				fAlignInfo<<CVect3(align.qnb).i/NAV_deg<<"	"<<CVect3(align.qnb).j/NAV_deg<<"	"<<CVect3(align.qnb).k/NAV_deg<<endl;
				cout<<CVect3(align.qnb).i/NAV_deg<<"	"<<CVect3(align.qnb).j/NAV_deg<<"	"<<CVect3(align.qnb).k/NAV_deg<<endl;
			}
		}
	}
	else	//reverse align
	{
		for(step=tStart*FRQ200-1; step>=tEnd*FRQ200; step--)
		{
			ReadData(step,1);
			align.IUpdate(wm, vm);	
			if((step%200==0)&&plogAlign) 
				(*plogAlign)<<CVect3(align.qnb)<<timu<<'\n';
		}
	}

	/////////////////
	CVect3 tmpatt=CVect3(align.qnb);

	//CVect3 tmpAtt( (tmpatt.i/NAV_deg)*NAV_deg,(tmpatt.j/NAV_deg-0.0018)*NAV_deg,(-358.973)*NAV_deg);	aqnb = CQuat(tmpAtt);align.qnb = aqnb;/*直接给初始姿态阵赋值  -0.738983   2.202105   359.013275 */
//	CVect3 tmpAtt( (-0.738983)*NAV_deg,(2.202105)*NAV_deg,(-359.01)*NAV_deg);	aqnb = CQuat(tmpAtt);align.qnb = aqnb;/*直接给初始姿态阵赋值  -0.738983   2.202105   359.013275 */

	/////////////////
	 posg = StartPos;

	aqnb = align.qnb; 
	avn = 0;
	CSINS sins(CIMU(TS5,4), aqnb, avn, posg); sins.pos=posg;  // LeverPos(sins,-lgps);
	pSG = new CSINSGPS(sins, lgps); pSOD = new CSINSOD(sins, inst, lod);
	//CVect3 att(align.qnb); printf("%f %f %f\n", att.i/NAV_deg, att.j/NAV_deg, att.k/NAV_deg);//HYQ Modified
	CVect3 att(align.qnb); printf("%f %f %f\n", att.i/NAV_deg, att.j/NAV_deg, att.k/NAV_deg);
	return tEnd;
}

int CNaviCtr::MovingAlign(int tStart, int tEnd)
{
	ReadData(tStart*FRQ200);
	CTrjAlign talign(TS5, 0*NAV_deg, posg);
	if(tStart<tEnd)
	{
		for(step=tStart*FRQ200+1; step<=tEnd*FRQ200; step++)
		{
			ReadData(step);	
			int res=talign.Update(wm, vm, vng, gpsValid&&(gpscoef==1));
			if((step%200)==0&&plogAlign) 
				(*plogAlign)<<talign.att<<talign.vn<<vng<<talign.gpspsi<<timu<<'\n';
			if(res>60) break;
		}
	}
	else
	{
		for(step=tStart*FRQ200-1; step>=tEnd*FRQ200; step--)
		{
			ReadData(step,1);	
			int res=talign.IUpdate(wm, vm, vng, gpsValid&&(gpscoef==1));
			if((step%200)==0&&plogAlign) 
				(*plogAlign)<<talign.att<<talign.vn<<vng<<talign.gpspsi<<timu<<'\n';
			if(res>60) break;
		}
	}
	aqnb = talign.qnb; avn = talign.vn;
	CSINS sins(CIMU(TS5,4), aqnb, avn, posg); sins.pos=LeverPos(sins,-lgps);
	pSG = new CSINSGPS(sins, lgps); 
	CVect3 att(talign.qnb); printf("%f %f %f %f\n", att.i/NAV_deg, att.j/NAV_deg, att.k/NAV_deg, talign.gpspsi/NAV_deg);
	return step/FRQ200;
}

void CNaviCtr::sins_gps(int tStart, int tEnd)
{
	if(tStart<tEnd)
	{
		for(step=tStart*FRQ200+1; step<=tEnd*FRQ200; step++)
		{
			ReadData(step);
			pSG->Update(wm, vm, vng, posg, gpsValid);
			static int m_StopTime=0;
			if(dSm != 0)
				m_StopTime = 0;

			if((step%200)==0&&plogNavi)
			{
				if(dSm == 0) m_StopTime++;
				int temp = m_StopTime;
				if(m_StopTime>3)
				{
				//	CVect3 Zero_vnL(0,0,0);
				//	pSG->vnL = Zero_vnL;
				}
				pSG->Output(0);
				//(*plogNavi)<<pSG->attL<<pSG->vnL<<pSG->posL<<vng<<posg<<pSG->kf->Xk<<timu<<'\n';
				fNavInfo.precision(9);
				fNavInfo<<pSG->attL.i/NAV_deg<<"	"<<pSG->attL.j/NAV_deg<<"	"<<pSG->attL.k/NAV_deg<<"	"<<pSG->vnL.i<<"	"<<pSG->vnL.j<<"	"<<pSG->posL.i/NAV_deg<<"	"<<pSG->posL.j/NAV_deg<<endl;
				//cout<<pSG->attL.i/NAV_deg<<"	"<<pSG->attL.j/NAV_deg<<"	"<<pSG->attL.k/NAV_deg<<"	"<<pSG->vnL.i<<"	"<<pSG->vnL.j<<"	"<<pSG->posL.i/NAV_deg<<"	"<<pSG->posL.j/NAV_deg<<endl;

			}
		}//		step--;
	}
	else	//逆向导航
	{
		for(step=tStart*FRQ200-1; step>=tEnd*FRQ200; step--)
		{
			ReadData(step,1);
			pSG->IUpdate(wm, vm, vng, posg, gpsValid&&(gpscoef<=3));
			if((step%200)==0&&plogNavi)
			{
				pSG->Output(0);
				(*plogNavi)<<pSG->attL<<pSG->vnL<<pSG->posL<<vng<<posg<<pSG->kf->Xk<<timu<<'\n';
			}
		}//		step++;
	}
}

void CNaviCtr::sins_od(int tStart, int tEnd)
{
	if(tStart<tEnd)
	{
		for(step=tStart*FRQ200+1; step<=tEnd*FRQ200; step++)
		{
			ReadData(step);
			pSOD->Update(wm, vm, dSm);
			if((step%200)==0&&plogNavi)
			{
				pSOD->Output(lgps);
				//(*plogNavi)<<pSOD->posDR<<pSOD->vnL<<pSOD->posL<<vng<<posg<<pSOD->kf->Xk<<timu<<'\n';
				fNavInfo.precision(9);
				fNavInfo<<pSG->attL.i/NAV_deg<<"	"<<pSG->attL.j/NAV_deg<<"	"<<pSG->attL.k/NAV_deg<<"	"<<pSG->vnL.i<<"	"<<pSG->vnL.j<<"	"<<pSG->posL.i/NAV_deg<<"	"<<pSG->posL.j/NAV_deg<<endl;

			}
		}//		step--;
	}
	else	//逆向导航
	{
		for(step=tStart*FRQ200-1; step>=tEnd*FRQ200; step--)
		{
			ReadData(step,1);
			pSOD->IUpdate(wm, vm, dSm);
			if((step%200)==0&&plogNavi)
			{
				pSOD->Output(lgps);
				(*plogNavi)<<pSOD->posDR<<pSOD->vnL<<pSOD->posL<<vng<<posg<<pSOD->kf->Xk<<timu<<'\n';
			}
		}//		step++;
	}
}

void CNaviCtr::sins_gps_od(int tStart, int tEnd)
{
	if(tStart<tEnd)
	{
		for(step=tStart*FRQ200+1; step<=tEnd*FRQ200; step++)
		{
			ReadData(step);
			pSG->Update(wm, vm, vng, posg, gpsValid&&(gpscoef<=3));
			pSOD->Update(wm, vm, dSm);
			if((step%200)==0&&plogNavi)
			{
				pSG->Output(0);
				pSOD->Output(lgps);
				(*plogNavi)<<pSG->attL<<pSG->vnL<<pSG->posL<<pSOD->vnL<<pSOD->posL<<vng<<posg<<timu<<'\n';
			}
		}//		step--;
	}
}

int CNaviCtr::ReadData(int i, int reverse)
{
//	if((*pIMU)(i)==NULL) 
//		return 0;
//	else 
	{ 
	//	timu = *(*pIMU)(i); 
		// wm = *(CVect3*)(*pIMU)(i+reverse,1)-eb*TS5; vm = *(CVect3*)(*pIMU)(i+reverse,4)-db*TS5;	// reverse imu: i+1 !!!
		double temp;
		fRawIMUData>>temp>>temp>>wm.i>>wm.j>>wm.k>>vm.i>>vm.j>>vm.k>>dSm>>temp>>temp>>temp>>temp>>temp; 
		temp = wm.j* PI/180 ;// PI/180
		wm.j = wm.i* PI/180;//PI/180;
		wm.i = temp;
		wm.k = -wm.k* PI/180;//;
		temp = vm.j;
		vm.j = vm.i;
		vm.i = temp;
		vm.k = -vm.k;
		wm-=eb*TS5;vm-=db*TS5;
		/*wm.j = -wm.j* PI/180 ;// PI/180
		wm.i = -wm.i* PI/180;//PI/180;
		wm.k = -wm.k* PI/180;//;

		vm.j = -vm.j;
		vm.i = -vm.i;
		vm.k = -vm.k;*/

		//if(pOD) dSm = *(*pOD)(i-10+reverse);
		if(i%FRQ200==0)
		{
			if(dSm == 0)
			{
				gpsValid = 0;//GPS信号频率1Hz，此处屏蔽gps信号
							posg =StartPos;//HYQ 人工赋值gps
							vng = 0;}
			int kk=i/FRQ200; 
			//double *p=(*pGPS)(kk);
			//tgps = *(p+1); 
			//posg = *(CVect3*)(p+2); gpscoef=*(p+5); vng = *(CVect3*)(p+6);



		}
		else
		{
			gpsValid = 0;
		}
	}
	//Display(i,FRQ200,100);
	return 1;
}