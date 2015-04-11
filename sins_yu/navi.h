/* navi c++ hearder file navi.h */
/*

*/

/*   */

#ifndef _NAVI_H
#define _NAVI_H

#include <math.h>
#include <stdio.h>
#include <conio.h>
#include <dos.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <fstream.h>

/* type re-define */
#ifndef BOOL
typedef int				BOOL;
#endif

#ifndef uchar
typedef unsigned char   uchar;
typedef unsigned int	uint;
typedef unsigned long   ulong;
#endif

/*constant define*/
#ifndef TRUE
#define TRUE	1
#define FALSE	0
#endif

#ifndef NULL
#define NULL	((void *)0)
#endif

#define PI		3.14159265358979
#define PI_2	(PI/2.0)
#define PI_4	(PI/4.0)
#define _2PI	(2.0*PI)

#define NAV_deg		(PI/180.)
#define NAV_min		(NAV_deg/60.)
#define NAV_sec		(NAV_min/60.)
#define NAV_hur		(3600.)
#define NAV_dph		(NAV_deg/NAV_hur)		/* arcdeg per hour i.e. arcsec per second */
#define NAV_ppm		(1/1000000.)
#define NAV_Re		6378160.0				/* earth coarse radius */
#define NAV_f		(1./298.257223563)
#define NAV_wie		7.2921151467e-5
#define NAV_g0		9.7803267714
#define NAV_mg		(NAV_g0/1000.)
#define NAV_ug		(NAV_g0/1000000.)
#define NAV_mrad	(1/1000.0)
#define NAV_mil		(_2PI/6000.0)			/* or 2*PI/6400.0 */
#define NAV_nm		(1853.0)				/* nautical mile */
#define NAV_km		(1000.0)

#define FRQ400	400							/* sampling frequency in Hz */
#define FRQ200	200
#define FRQ100	100
#define FRQ50	50
#define TS2p5	(1.0/FRQ400)				/* sampling frequency in ms */
#define TS5		(1.0/FRQ200)
#define TS10	(1.0/FRQ100)
#define TS20	(1.0/FRQ50)

#define EPS		2.22044604925031e-016
#define INF		1.0e100
#define lfEND	1234567890.123456

#define EST i
#define NTH j
#define UP  k
#define RGT i
#define FNT j
#define PCH i
#define RLL j
#define YAW k
#define LTI i
#define LGI j
#define HGT k
#define XXX i
#define YYY j
#define ZZZ k

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#define Display(step,frq,n)	{ if((step)%((n)*(frq))==0) printf("%d\n", (int)((step)/(frq))); }

#define MMAX1	33					/* max matrix dimension */
#define MMAX2	(MMAX1*MMAX1)

// class define
class CVect3;	class CMat3;	class CQuat;
class CEarth;	class CIMU;		class CSINS;
class CMat;		class CKalman;	class CUnit;	

// global variables and functions, can not be changed in any way
extern CVect3	O3;
extern CMat3	O33, I33;
extern CQuat	Q1;
extern CEarth	eth84, eth54;

class CVect3 
{
public:
	double i, j, k;

	CVect3(void);
	CVect3(double val);
	CVect3(double xx, double yy, double zz);
	CVect3(CMat3 &Cnb);
	CVect3(CQuat &qnb);

	CVect3 operator+(CVect3 &v);
	CVect3 operator-(CVect3 &v);
	CVect3 operator*(CVect3 &v); // vector cross multiply
	CVect3 operator*(double f);
	friend CVect3 operator*(double f, CVect3 &v);
	CVect3 operator*(CMat3 &m);		// V*M = M^T*V
	CVect3& operator+=(CVect3 &v);
	CVect3& operator-=(CVect3 &v);
	CVect3& operator*=(double f);
//	CVect3 operator=(double f);  // No! for need 'CVect3 operator=(CVect3 &v);' firstly!
	CVect3 operator/(double f);
	CVect3& operator/=(double f);
	friend CVect3 operator-(CVect3 &v);
	double& operator()(int i);

	friend CVect3 operator~(CVect3 &v) { return v; }
	friend double operator!(CVect3 &v);				/* norm */
	friend double VDot(CVect3 v1, CVect3 &v2);
	friend CQuat RV2Q(CVect3 &v);
	friend CMat3 P2Cen(CVect3 &v);
	friend CVect3 PP2Polar(CVect3 &vScr, CVect3 &vObj);
};

class CQuat
{
public:
	double q0, q1, q2, q3;

	CQuat(void);
	CQuat(double val);
	CQuat(double pch, double rll, double yaw);
	CQuat(double qq0, double qq1, double qq2, double qq3);
	CQuat(CVect3 &att);
	CQuat(CMat3 &Cnb);

	CQuat operator+(CVect3 &v);
	CQuat operator-(CVect3 &v);
	CVect3 operator-(CQuat &quat);
	CQuat operator*(CQuat &q);
	CVect3 operator*(CVect3 &v);
	CQuat& operator*=(CQuat &q);
	CQuat& operator+=(CVect3 &v);
	CQuat& operator-=(CVect3 &v);
	CQuat operator*(double f);
	double& operator()(int i);

	friend CQuat operator ~(CQuat &q);		// quaternion conjugate
	friend double operator!(CQuat &q);
	friend CQuat Normlize(CQuat &q);
	friend CVect3 Q2RV(CQuat &q);
};

class CMat3 
{
public:
	double e00, e01, e02, e10, e11, e12, e20, e21, e22;

	CMat3(void);
	CMat3(double val);
	CMat3(double xx, double yy, double zz);
	CMat3(double xx, double xy, double xz, 
		  double yx, double yy, double yz,
		  double zx, double zy, double zz );
	CMat3(CVect3 &att);
	CMat3(CQuat &qnb);
	CMat3(CVect3 &v0, CVect3 &v1, CVect3 &v2);
	CMat3(CVect3 &vn1, CVect3 &vn2, CVect3 &vb1, CVect3 &vb2);

	CMat3 operator+(CMat3 &m);
	CMat3 operator-(CMat3 &m);
	CMat3 operator*(CMat3 &m);
	CVect3 operator*(CVect3 &v);
	CMat3 operator*(double f);
	friend CMat3 operator*(double f, CMat3 &m);
	CMat3& operator+=(CMat3 &m);
	CMat3& operator-=(CMat3 &m);
	CMat3& operator*=(CMat3 &m);
	CMat3& operator*=(double f);
	CMat3& operator++();
	friend CMat3 operator-(CMat3 &m);
	friend CMat3 operator/(double f, CMat3 &m);
	double& operator()(int i, int j);

	friend CMat3 operator ~(CMat3 &m);		// matrix transpose
	CMat3 operator^(int n);					// matrix power
	friend double operator!(CMat3 &m);		/* norm */
	friend CMat3 OptOrtho(CMat3 &m);		/* optimal orthogonality */
	friend CMat3 ASkew(CVect3 &v);
	friend CMat3 Diag(CVect3 &v);
	friend CMat3 Cw(CVect3 &a);
	friend CMat3 Cw_1(CVect3 &a);
};

class CMat
{
public:
	int row, clm, rc;
	double *d;

	CMat(int r=0, int c=0, double *d0=NULL);
//	CMat(CMat &m, double *d0=NULL);  error !

	friend double* MatAlloc(CMat *m, ...);
	double& operator()(int r, int c=0);
	void Clear();
	CMat& Setf(double f, ...);
	CMat& Setv(int r, int c, CVect3 &v);
	CMat& Setm(int r, int c, CMat3 &m);
	CMat& Setm(CMat3 *m, ...);
	CMat& DelRow(int rowN);
	CMat& DelClm(int clmN);
	CMat& DelRC(int rcN);

	CMat operator=(CMat &m);
	CMat operator=(double f);
	CMat operator+(CMat &m);
	CMat operator-(CMat &m);
	CMat operator*(CMat &m);
	CMat operator+(double f);
	CMat operator-(double f);
	CMat operator*(double f);
	CMat& operator+=(CMat &m);
	CMat& operator-=(CMat &m);
	CMat& operator*=(CMat &m);
	CMat& operator*=(double f);
	CMat& operator++();
	friend CMat operator+(double f, CMat &m);
	friend CMat operator-(double f, CMat &m);
	friend CMat operator*(double f, CMat &m);
	friend CMat operator-(CMat &m);

	CMat operator&(CMat &m);		/* element-wise multiply (one by one element) */

	CMat operator^(int n);
	friend CMat operator~(CMat &m);			/* transpose */
	friend CMat Diag(double f, ...);
	friend CMat Diag(CMat &m0);
	friend void symmetric(CMat& m);
	friend CMat eye(int n);
	friend CMat mChol(CMat &m);
	CMat SparseMul(CMat &m0);
	CMat SparseMul1(CMat &m0);
	friend void APATBuild(char *fname, CMat &A);
	friend void APATBuild1(char *fname, CMat &A);
};

class CEarth
{
public:
	double a, b;
	double f, e, e2, ep, ep2;
	double wie;

	double lti, hgt, sl, sl2, sl4, cl, tl, RMh, RNh, clRNh, f_RMh, f_RNh, f_clRNh;;
	CVect3 wnie, wnen, wnin, gn;

	CEarth();
	CEarth(double a0, double f0=NAV_f, double wie0=NAV_wie);

	void Update(CVect3 &pos, CVect3 &vn=CVect3(0));
	friend CVect3 vn2dpos(CEarth &eth, CVect3 &vn, double ts=1.0);
};

class CIMU
{
public:
	double ts, tss, tt, *pcf;
	int samples, ksample;
	CVect3 wm[5], vm[5], wmm, vmm, phim, dvbm;

	CIMU(void);
	CIMU(double ts, int samples0=2);
	int Update(CVect3 &wmi, CVect3 &vmi=CVect3(0));
};

class CSINS
{
public:
	CEarth eth;
	CIMU imu;
	CQuat qnb, qnbk;
	CMat3 Cnb, Cnbk;
	CVect3 	wb, fb, fn, att, vn, pos, attk, vnk, posk;

	CSINS(void);
	CSINS(CIMU &imu0, CQuat &qnb0, CVect3 &vn0, CVect3 &pos0, CEarth &eth0=eth84);

	int Update(CVect3 &wm, CVect3 &vm);
	int IUpdate(CVect3 &wm, CVect3 &vm);
	friend CVect3 LeverPos(CSINS &sins, CVect3 &lever, CVect3 *p0=NULL);
	friend CVect3 LeverVn(CSINS &sins, CVect3 &lever, CVect3 *v0=NULL);
	void ErrCoef(CMat3 &S1, CMat3 &S2, 
				CMat3 &M1, CMat3 &M2, CMat3 &M3, CMat3 &M4, CMat3 &M5, CMat3 &M6);
};

class CDR
{
public:
	double vv, K1, K2, K3, UK1, UK2, UK3, UK4;
	CEarth eth;
	CQuat qnb;
	CVect3 vn, pos, drInst, PJKD, wnc;

	CDR(void);
	CDR(CVect3 &att0, CVect3 &pos0, CVect3 &inst0, CEarth &eth0=eth84);

	void Update(CIMU &imu, double dS);
	void Update0(CIMU &imu, double dS);
	void Update1(CIMU &imu, double dS);
	void IUpdate(CIMU &imu, double dS);
};

class CAlign_i0
{
	double tk;
	long t0, t1, t2;
	CVect3 vib0, vi0, Pib01, Pib02, Pi01, Pi02, tmpPib0, tmpPi0, pos0;
	CQuat qib0b;
	CIMU imu;
	CEarth eth;
public:
	CQuat qnb;
	CAlign_i0(void);
	CAlign_i0(CIMU &imu0, CVect3 &pos0, CEarth &eth0=eth84);

	int Update(CVect3 &wm, CVect3 &vm);
	int IUpdate(CVect3 &wm, CVect3 &vm);
};

class CTrjAlign
{
	int firstVng;
	CEarth eth;
public:
	double ts, dyaw, psi, gpspsi;
	CQuat qnb;
	CMat3 Cnb;
	CVect3 att, vn, pos;
	CKalman *kf;
	CTrjAlign(double ts0, double dyaw0, CVect3 &pos0, CEarth &eth0=eth84);

	int Update(CVect3 &wm, CVect3 &vm, CVect3 &vng, int trjvalid);
	int IUpdate(CVect3 &wm, CVect3 &vm, CVect3 &postrj, int trjvalid);
};

class CFload
{
	FILE *f;
	char *format;
	int clm;
public:
	double *d;

	CFload();
	CFload(char *fname, int clm=0, char separator=' ', char endChar=' ');
	~CFload();

	double operator()(int n);
	CVect3 operator[](int n);

	long Skip(long lines);
	BOOL Load(void);
};

class CBinFile
{
	int row, clm, size;
public:
	double *pd;

	CBinFile(char *fname, int clm, int row=2147483647/sizeof(double));
	~CBinFile();

	double* operator()(int m, int n=0);
};

class CFlog
{
	char fmode;
	FILE *f;
public:
	CFlog();
	CFlog(char *fname, char mode='t');
	~CFlog();

	friend CFlog& operator<<(CFlog& flog, char c);
	friend CFlog& operator<<(CFlog& flog, char *str);
	friend CFlog& operator<<(CFlog& flog, int i);
	friend CFlog& operator<<(CFlog& flog, unsigned int i);
	friend CFlog& operator<<(CFlog& flog, long l);
	friend CFlog& operator<<(CFlog& flog, float f);
	friend CFlog& operator<<(CFlog& flog, double d);
	friend CFlog& operator<<(CFlog& flog, CVect3 &v);
	friend CFlog& operator<<(CFlog& flog, CQuat &q);
	friend CFlog& operator<<(CFlog& flog, CMat3 &m);
	friend CFlog& operator<<(CFlog& flog, CMat &m);
	friend CFlog& operator<<(CFlog& flog, CSINS &sins);
	friend CFlog& operator<<(CFlog& flog, CKalman &kf);
};

class CKalman
{
	double *pd;
public:
	int q, r;
	double chi2, forgetting;
	CMat Xk, fbXk, fbCoef, smXk, smCoef, Zk, rk, Ft, Fk, Qt, Qk, Gt, Pxz, Pzz, Pzz_1, Pk, Hk, Rk, Iq, Kk, Xkk_1, Pkk_1, lamda;
	void (*pAPAT)(CMat &mA, CMat &mPx, CMat &mPy, CMat &vX, CMat &vY);

	CKalman(int q0, int r0, double forgetting0=1.0001);
	~CKalman();

	void Discrete(double ts);
	void TimeUpdate(void);
	void MeasureUpdate(void);
	void Update(int obs=0);
	void InitFBCoef(double ts, double coef, ...);
	void InitSMCoef(double ts, double coef, ...);
	void FeedBack(void);
	void Smooth(void);
};

class CSINSGPS
{
public:
	CMat *psm;
	CVect3 *pfi, *pvn, *ppos, *peb, *pdb, *pkg1, *pkg2, *pkg3, *pka1, *pka23, *plv, *pzkv, *pzkp;
	CVect3 *psfi, *psvn, *pspos, *pseb, *psdb, *pskg1, *pskg2, *pskg3, *pska1, *pska23, *pslv;
	CVect3 attL, vnL, posL, lever, eb, db;
	CQuat qnbL;
	CMat3 CnbL;
	CSINS sins;
	CKalman *kf;
	CSINSGPS(CSINS &sins0, CVect3 &lever0=O3);
	~CSINSGPS();

	void Init(void);
	int Update(CVect3 &gyro, CVect3 &acc, CVect3 &vng, CVect3 &gps, int obs=0);
	int IUpdate(CVect3 &gyro, CVect3 &acc, CVect3 &vng, CVect3 &gps, int obs=0);
	void InitFBCoef(int idx);
	void FeedBack(void);
	void Output(int isIMU=0);
};

class CSINSOD
{
public:
	double ST, dSm1;
	CMat *psm;
	CVect3 *pfi, *pvn, *ppos, *peb, *pdb, *pinst, *plv, *pzkv, *pzkp;
	CVect3 *psfi, *psvn, *pspos, *pseb, *psdb, *psinst, *pslv;
	CVect3 attL, vnL, posL, posDR, inst, lever, eb, db, prj, prj1, SVins, SVod;
	CQuat qnbL;
	CMat3 CnbL, SMod, SMlv;
	CSINS sins;
	CKalman *kf;
	CSINSOD(CSINS &sins0, CVect3 &inst0, CVect3 &lever0=O3);
	~CSINSOD();

	void Init(void);
	int Update(CVect3 &gyro, CVect3 &acc, double dSm);
	int IUpdate(CVect3 &gyro, CVect3 &acc, double dSm);
	CVect3 ODPrj(CVect3 &inst);
	void InitFBCoef(int idx);
	void FeedBack(void);
	void Output(CVect3 &lv=O3);
};

class CUnit
{
	friend double r2d(double rad);
	friend double d2r(double deg);
	friend double r2dm(double rad);
	friend double dm2r(double deg);
	friend double r2dms(double rad);
	friend double dms2r(double deg);
	friend double r2m(double rad);
	friend double m2r(double m);
};

class CTDisp
{
	time_t t0, t;
public:
	CTDisp() { time(&t0); };
	inline long Elasped(void) { return time(&t)-t0; };
	~CTDisp() { time(&t); printf("Time elasped is %lds\n", Elasped()); };
};

class CNaviCtr
{
public:
	int step, gpsValid;
	double timu, tgps, gpscoef, dSm;
	CQuat aqnb;
	CVect3 wm, vm, eb, db, vng, posg, avn, lgps, lod, inst;
	CSINSGPS *pSG;
	CSINSOD *pSOD;
	CBinFile *pIMU, *pOD, *pGPS;
	CFlog *plogAlign, *plogNavi;

	CNaviCtr(CVect3 &eb0, CVect3 &db0, CVect3 &lgps0, CVect3 &lod0, CVect3 &inst0);
	~CNaviCtr();

	void FileOpen(char *path, char *imu, char *od, char *gps, char *align, char *navres);
	int StaticAlign(int tStart, int tEnd);
	int MultiPositionAlign(int tStart, int tEnd);
	int MovingAlign(int tStart, int tEnd);
	void sins_gps(int tStart, int tEnd);
	void sins_od(int tStart, int tEnd);
	void sins_gps_od(int tStart, int tEnd);
	int ReadData(int i, int reverse=0);
};

/* miscellaneous */
BOOL assert(BOOL b);
int signE(double val, double eps);
#define sign(val)			signE(val, EPS)
#define signX(val, x)		signE((val)-(x), EPS);
#define signXE(val, x, eps)	signE((val)-(x), eps);
double range(double val, double minVal, double maxVal);
int nextIndex(int curIdx, int maxIdx);
double asinEx(double x);
double atan2Ex(double y, double x);

#endif