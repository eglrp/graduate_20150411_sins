#include "navi.h"

static double conefactors[5][4] = {				//Ô²×¶Îó²î²¹³¥ÏµÊý
	{2./3},										//2
	{9./20, 27./20},							//3
	{54./105, 92./105, 214./105},				//4
	{250./504, 525./504, 650./504, 1375./504}	//5
};

CIMU::CIMU(void)
{
}

CIMU::CIMU(double ts0, int samples0)
{
	ts = ts0; tss = samples0*ts; tt = 0.0;
	samples = samples0;  
	assert(samples<=5);
	pcf = conefactors[samples-2];
	ksample = 0;
}

int CIMU::Update(CVect3 &wmi, CVect3 &vmi)
{
	wm[ksample] = wmi,	vm[ksample] = vmi;
	ksample ++; tt += ts;
	if(ksample==1)	
		{ wmm = wmi, vmm = vmi;	}
	else
		{ wmm += wmi, vmm += vmi; }
	if(ksample==samples)
	{
		CVect3 cm(0), sm(0);
		for(int i=0; i<samples-1; i++)
		{
			cm += pcf[i]*wm[i];
			sm += pcf[i]*vm[i];
		}
		phim = wmm + cm*wm[i];
		dvbm = vmm + 1.0/2*(wmm*vmm) + (cm*vm[i]+sm*wm[i]);
		ksample = 0;
	}
	return ksample;
}
