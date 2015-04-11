#include "navi.h"


/* memTempAlloc作为临时中间变量的矩阵内存块分配，使用等号赋值后blocks计数自动归零，仅限于在CMat类中使用 */
int maxBlockUsed = 0, blocks = 0; 

static double* memTempAlloc(void)
{
#define BLOCKS	10
	static double data[MMAX2*BLOCKS], *pd=&data[0];
	pd += MMAX2;
	if(pd>=&data[MMAX2*BLOCKS])
		pd = data;
	blocks ++;
	return pd;
}

/* MatAlloc可用于多个矩阵的内存初始化分配，自动计算所需的内存数量 */
double* MatAlloc(CMat *m, ...)
{
	int len = 0;
	CMat *pm = m;
	// calculate length needed
	va_list vl;
	va_start(vl, m);
	while(1)
	{
		len += pm->rc;
		pm = va_arg(vl, CMat*);
		if(!pm)	break;
	}
	va_end(vl);
	// allocate memory
	double *pd = new double[len]; 
	memset(pd, 0, len*sizeof(double));
	// distribute memory
	len = 0;
	pm = m;
	va_start(vl, m);
	while(1)
	{
		pm->d = &pd[len];
		len += pm->rc;
		pm = va_arg(vl, CMat*);
		if(!pm)	break;
	}
	va_end(vl);
	return pd;
}

CMat::CMat(int r, int c, double *d0)
{
	if(!assert(r<=MMAX1 && c<=MMAX1))
		{ printf("\n\tMatrix dimension error! Try to adjust macro 'MMAX1'.\n\n"); exit(0); }
	row = r, clm = c, rc = row*clm;
	d = d0;
	return;
}

void CMat::Clear()
{
	memset(this->d, 0, this->rc*sizeof(double));
	return;
}

double& CMat::operator()(int r, int c)
{
	return this->d[r*this->clm+c];
}

CMat& CMat::Setf(double f, ...)
{
	double *pd = this->d;
	*pd++ = f;
	va_list vl;
	va_start(vl, f);
	for(int i=1; i<this->rc; i++)
		*pd++ = va_arg(vl, double);
	va_end(vl);
	return *this;
}

CMat& CMat::Setv(int r, int c, CVect3 &v)
{
	assert(r<this->row&&c<this->clm);
	double *p = &d[r*this->clm+c];
					*p = v.i; 
	p += this->clm; *p = v.j; 
	p += this->clm; *p = v.k; 
	return *this;
}

CMat& CMat::Setm(int r, int c, CMat3 &m)
{
	assert(r<this->row&&c<this->clm);
	double *p = &m.e00;
	for(int k=0,kk=r*this->clm+c; k<9; k+=3,kk+=this->clm)
	{
		for(int s=0; s<3; s++)
		{
			this->d[kk+s] = p[k+s];
		}
	}
	return *this;
}

CMat& CMat::Setm(CMat3 *m, ...)
{
	CMat3 *pm = m;
	double *p = &pm->e00;
	va_list vl;
	va_start(vl, m);
	for(int i=0; i<this->row; i+=3)
	{
		for(int j=0; j<this->clm; j+=3)
		{
			if(pm!=&O33)
			{
				for(int k=0,kk=i*this->clm+j; k<9; k+=3,kk+=this->clm)
				{
					for(int s=0; s<3; s++)
					{
						this->d[kk+s] = p[k+s];
					}
				}
			}
			pm = va_arg(vl, CMat3*);
			if(!pm)
				return *this;
			p = &pm->e00;
		}
	}
	va_end(vl);
	return *this;
}

CMat& CMat::DelRow(int rowN)
{
	assert(rowN<this->row && rowN>=0);
	if(rowN != row-1)
	{
		double *p=&this->d[rowN*clm], *pNext=p+clm;
		for(int i=(row-1-rowN)*clm; i>0; i--)
			*p++ = *pNext++;
	}
	this->rc -= clm, this->row--;
	return *this;
}

CMat& CMat::DelClm(int clmN)
{
	assert(clmN<this->clm && clmN>=0);
	double *p=&this->d[clmN], *pNext=p+1;
	for(int i=0; i<this->row; i++)
	{
		for(int j=1; j<this->clm; j++)
			*p++ = *pNext++;
		pNext++;
	}
	this->rc -= row, this->clm--;
	return *this;
}

CMat& CMat::DelRC(int rcN)
{
	DelRow(rcN);
	return DelClm(rcN);
}

CMat CMat::operator=(CMat &m)
{
	row = m.row, clm = m.clm, rc = m.rc;
	if(m.d&&this->d)
		memcpy(this->d, m.d, rc*sizeof(double));
	maxBlockUsed = max(maxBlockUsed,blocks), blocks = 0; assert(maxBlockUsed<BLOCKS);
	return *this;
}

CMat CMat::operator=(double f)
{
	double *p = this->d;
	for(int i=0; i<this->rc; i++) 
		*p++ = f;
	return *this;
}

CMat CMat::operator+(CMat &m)
{
	assert(this->row==m.row&&this->clm==m.clm);
	double *dd = memTempAlloc();
	double *pd=dd, *p=this->d, *pm=m.d;
	for(int i=0; i<m.rc; i++)
		*pd++ = (*p++) + (*pm++);
	return CMat(m.row,m.clm,dd);
}

CMat& CMat::operator+=(CMat &m)
{
	assert(this->row==m.row&&this->clm==m.clm);
	double *p=this->d, *pm=m.d;
	for(int i=0; i<m.rc; i++)
		*p++ += *pm++;
	return *this;
}

CMat CMat::operator-(CMat &m)
{
	assert(this->row==m.row&&this->clm==m.clm);
	double *dd = memTempAlloc();
	double *pd=dd, *p=this->d, *pm=m.d;
	for(int i=0; i<m.rc; i++)
		*pd++ = *p++ - *pm++;
	return CMat(m.row,m.clm,dd);
}

CMat& CMat::operator-=(CMat &m)
{
	assert(this->row==m.row&&this->clm==m.clm);
	double *p=this->d, *pm=m.d;
	for(int i=0; i<m.rc; i++)
		*p++ -= *pm++;
	return *this;
}

CMat CMat::operator*(double f)
{
	double *dd = memTempAlloc();
	double *pd=dd, *p=this->d;
	for(int i=0; i<this->rc; i++)
		*pd++ = (*p++) * f;
	return CMat(this->row,this->clm,dd);
}

CMat& CMat::operator*=(double f)
{
	double *pd=this->d;
	for(int i=0; i<this->rc; i++)
		(*pd++) *= f;
	return *this;
}

CMat operator*(double f, CMat &m)
{
	return m*f;
}

CMat CMat::operator*(CMat &m0)
{
	double *dd = memTempAlloc();
	assert(this->clm==m0.row);
	//dd[m*n]=p1[m*k]*p2[k*n]
	int m=this->row, k=this->clm, n=m0.clm;
	double *p=dd, *p1i=this->d, *p2=m0.d;
	for(int i=0; i<m; i++,p1i+=k)
	{
		for(int j=0; j<n; j++)
		{
			double f=0.0, *p1is=p1i, *p2sj=&p2[j];
			for(int s=0; s<k; s++,p1is++,p2sj+=n)
			{
				f += *p1is * *p2sj;
			}
			*p++ = f;
		}
	}
	return CMat(m,n,dd);
}

CMat CMat::SparseMul(CMat &m0)
{
	double *dd = memTempAlloc();
	assert(this->clm==m0.row);
	//dd[m*n]=p1[m*k]*p2[k*n]
	int m=this->row, k=this->clm, n=m0.clm;
	double *p=dd, *p1i=this->d, *p2=m0.d;
	for(int i=0; i<m; i++,p1i+=k)
	{
		int idx1[MMAX1], idx2[MMAX1];
		double *p1is=p1i;
		for(int s=0, ns=0, j1=0; s<n; s++, ns+=n) // 找出i行不为零的元素，共计j1个，序号存入数组idx1,idx2
		{
			if(*p1is++!=0.0)
			{
				idx1[j1] = s;
				idx2[j1++] = ns; 
			}
		}
		for(int j=0; j<n; j++)
		{
			double f=0.0, *p1is=p1i, *p2sj=&p2[j];
			for(int s=0; s<j1; s++)
			{
				f += p1is[idx1[s]] * p2sj[idx2[s]];
			}
			*p++ = f;
		}
	}
	return CMat(m,n,dd);
}

CMat CMat::SparseMul1(CMat &m0)
{
	double *dd = memTempAlloc();
	assert(this->clm==m0.row);
	//dd[m*n]=p1[m*k]*p2[k*n]
	int m=this->row, k=this->clm, n=m0.clm;
	double *p=dd, *p1i=this->d, *p2=m0.d;
	for(int i=0; i<m; i++,p1i+=k)
	{
		int idx1[MMAX1], idx2[MMAX1];
		double *p1is=p1i;
		for(int s=0, ns=0, j1=0; s<n; s++, ns+=n) // 找出i行不为零的元素，共计j1个，序号存入数组idx1,idx2
		{
			if(*p1is++!=0.0)
			{
				idx1[j1] = s;
				idx2[j1++] = ns; 
			}
		}
		for(int j=0; j<n; j++)
		{
			double f=0.0, *p1is=p1i, *p2sj=&p2[j];
			if(j>=i)
			{
				for(int s=0; s<j1; s++)
				{
					f += p1is[idx1[s]] * p2sj[idx2[s]];
				}
			}
			else
			{
				f = dd[j*n+i];		// 对称
			}
			*p++ = f;
		}
	}
	return CMat(m,n,dd);
}

CMat& CMat::operator*=(CMat &m0)
{
	*this = *this*m0;
	return *this;
}

CMat CMat::operator&(CMat &m0)
{
	double *dd = memTempAlloc();
	CMat m(m0.row, m0.clm, dd);
	double *pd=dd, *pt=this->d, *p0=m0.d;
	for(int i=0; i<m0.rc; i++)
		*pd++ = *pt++ * *p0++;
	return m;
}

CMat operator-(double f, CMat &m)
{
	return (-m)+f;
}

CMat operator+(double f, CMat &m)
{
	return m+f;
}

CMat CMat::operator+(double f)
{
	double *dd = memTempAlloc();
	double *pd=dd, *p=this->d;
	for(int i=0; i<this->rc; i++)
		*pd++ = (*p++) + f;
	return CMat(this->row,this->clm,dd);
}

CMat CMat::operator-(double f)
{
	return *this+(-f);
}

CMat& CMat::operator++()	// m += I
{ 
	double *p=d;
	for(int i=0; i<row; i++, p+=row+1)
		*p += 1.0;
	return *this;
}

static int brinv(double a[], int n);

CMat CMat::operator^(int n)
{
//	static double dd[MMAX], dd1[MMAX];
	double *dd = memTempAlloc(), *dd1 = memTempAlloc();
	CMat mtmp(this->row,this->row,dd), mtmp1(this->row,this->row,dd1);
	if(n==0)
	{
		mtmp = eye(this->row);
	}
	else if(n<=-1)
	{
		mtmp = *this;
		if(!brinv(mtmp.d, mtmp.row))
			mtmp = eye(mtmp.row)*12345678;
		n = -n;
	}
	else
	{
		mtmp = *this;
	}
	mtmp1 = mtmp;
	assert( n<BLOCKS-2 );		// make sure: n < BLOCKS
	for(int i=1; i<n; i++)
		mtmp1 = mtmp1*mtmp;
	return mtmp1;
}

CMat operator-(CMat &m)
{ 
	double *dd = memTempAlloc();
	double *p=dd, *pm=m.d;
	for(int i=0; i<m.rc; i++)
		*p++ = -*pm++;
	return CMat(m.row,m.clm,dd);
}

CMat operator~(CMat &m0)
{
	double *dd = memTempAlloc();
	double *pm=m0.d;
	for(int i=0; i<m0.row; i++)
	{
		for(int j=i; j<m0.rc; j+=m0.row)
		{
			dd[j] = *pm++;
		}
	}
	return CMat(m0.clm,m0.row,dd);
}

CMat Diag(double f, ...)
{
	double *dd = memTempAlloc();
	va_list vl;
	va_start(vl, f);
	for(int i=0; i<MMAX1; i++)
	{
		dd[i] = f;
		f = va_arg(vl, double);
		if(f==lfEND)
			break;
	}
	va_end(vl);
	return Diag(CMat(i+1,1,dd));
}

CMat Diag(CMat &m0)
{
	double *dd = memTempAlloc();
	if(m0.row==1 || m0.clm==1)		// 向量转成对角阵
	{
		CMat m(m0.rc, m0.rc, dd);
		m.Clear();
		for(int i=0; i<m0.rc; i++)
			m.d[i*m0.rc+i] = m0.d[i];
		return m;
	}
	else		// 矩阵提取对角线元素组成向量
	{
		CMat m(m0.row, 1, dd);
		for(int i=0; i<m0.row; i++)
			m.d[i] = m0.d[i*m0.row+i];
		return m;
	}
}

void symmetric(CMat& m)	// m = (m+(~m)) / 2;
{
	int i, j;
	double *pij, *pji;
	for(i=0; i<m.row; i++)
	{
		for(j=i+1,pij=&m.d[i*m.row+j],pji=&m.d[j*m.row+i]; j<m.row; j++,pij++,pji+=m.row)
		{
			*pij = *pji = (*pij+*pji)*0.5;
		}
	}
}

CMat eye(int n)
{
	double *dd = memTempAlloc();
	CMat m(n,n,dd);
	m.Clear();
	for(int i=0; i<m.rc; i+=n+1)
	{
		dd[i] = 1.0;
	}
	return m;
}

int bchol(double a[], int n, double *det);

CMat mChol(CMat &m)
{
	double *dd = memTempAlloc();
	CMat mtmp(m.row,m.clm,dd);
	if(m.d!=dd)
		mtmp = m;
	double det;
	bchol(mtmp.d, mtmp.row, &det);
	return mtmp;
}

static int brinv(double a[], int n)
{ 
	int /**is,*js,*/i,j,k,l,u,v;
	double d,p;

	/*is=malloc(n*sizeof(int));
	js=malloc(n*sizeof(int));*/
	int is[MMAX1], js[MMAX1];
	
	for (k=0; k<=n-1; k++)
	{ 
		d=0.0;
		for (i=k; i<=n-1; i++)
			for (j=k; j<=n-1; j++)
			{ 
				l=i*n+j; p=fabs(a[l]);
				if (p>d) 
				{ 
					d=p; is[k]=i; js[k]=j;
				}
			}
//			if (d+1.0==1.0)
			if (d+EPS==EPS)  // ygm改，适当缩小！
			{ 
				/*free(is); free(js); printf("err**not inv\n");*/
				return(0);
			}
			if (is[k]!=k)
				for (j=0; j<=n-1; j++)
				{ 
					u=k*n+j; v=is[k]*n+j;
					p=a[u]; a[u]=a[v]; a[v]=p;
				}
				if (js[k]!=k)
					for (i=0; i<=n-1; i++)
					{ 
						u=i*n+k; v=i*n+js[k];
						p=a[u]; a[u]=a[v]; a[v]=p;
					}
					l=k*n+k;
					a[l]=1.0/a[l];
					for (j=0; j<=n-1; j++)
						if (j!=k)
						{ 
							u=k*n+j; a[u]=a[u]*a[l];
						}
						for (i=0; i<=n-1; i++)
							if (i!=k)
								for (j=0; j<=n-1; j++)
									if (j!=k)
									{ 
										u=i*n+j;
										a[u]=a[u]-a[i*n+k]*a[k*n+j];
									}
									for (i=0; i<=n-1; i++)
										if (i!=k)
										{ 
											u=i*n+k; a[u]=-a[u]*a[l];
										}
	}
	for (k=n-1; k>=0; k--)
	{ 
		if (js[k]!=k)
			for (j=0; j<=n-1; j++)
			{ 
				u=k*n+j; v=js[k]*n+j;
				p=a[u]; a[u]=a[v]; a[v]=p;
			}
			if (is[k]!=k)
				for (i=0; i<=n-1; i++)
				{ 
					u=i*n+k; v=i*n+is[k];
					p=a[u]; a[u]=a[v]; a[v]=p;
				}
	}
	/*free(is); free(js);*/
	return(1);
}

//return lower triangular matrix
int bchol(double a[], int n, double *det)
{ 
	int i,j,k,u,l;
    double d;
    if((a[0]+1.0==1.0)||(a[0]<0.0))
	{ 
		printf("fail\n"); return(-2);
	}
    a[0]=sqrt(a[0]);
    d=a[0];
    for(i=1; i<=n-1; i++)
	{ 
		u=i*n; a[u]=a[u]/a[0];
	}
    for(j=1; j<=n-1; j++)
	{ 
		l=j*n+j;
		for(k=0; k<=j-1; k++)
		{ 
			u=j*n+k; a[l]=a[l]-a[u]*a[u];
		}
		if((a[l]+1.0==1.0)||(a[l]<0.0))
		{ 
			printf("fail\n"); return(-2);
		}
		a[l]=sqrt(a[l]);
		d=d*a[l];
		for(i=j+1; i<=n-1; i++)
		{ 
			u=i*n+j;
			for (k=0; k<=j-1; k++)
				a[u]=a[u]-a[i*n+k]*a[j*n+k];
			a[u]=a[u]/a[l];
		}
	}
    *det=d*d;
    for (i=0; i<=n-2; i++)
		for (j=i+1; j<=n-1; j++)
			a[i*n+j]=0.0;
	return(2);
}

void APATBuild(char *fname, CMat &A)
{
	int N = A.row, m, n, k, op_add=0, op_mul=0;
	double *pd = A.d;
	FILE *f=fopen(fname, "wt");
    fprintf(f, "#include \"navi.h\"\n\n");
    fprintf(f, "#define KF_RAPID\n\n");
    fprintf(f, "void APAT(CMat &mA, CMat &mPx, CMat &mPy, CMat &vX, CMat &vY)\n{\n");
    fprintf(f, "\tdouble C[%d*%d], *A=mA.d, *P=mPx.d, *B=mPy.d, *X=vX.d, *Y=vY.d;\n\n",N,N);
	// Y = A*X
    for(m=0; m<N; m++)
	{
		fprintf(f, "\tY[%3d]=", m);
		for(n=0; n<N; n++)
		{
			if(pd[m*N+n]==1.0)		{ fprintf(f, "+       X[%3d]", n); op_add++; }
			else if(pd[m*N+n]==-1.0){ fprintf(f, "-       X[%3d]", n); op_add++; }
			else if(pd[m*N+n]!=0.0) { fprintf(f, "+A[%3d]*X[%3d]", m*N+n,n); op_mul++; op_add++; }
		}
		fprintf(f, ";\n");
	}
    fprintf(f, "\n");
	// C = A*Px
    for(m=0; m<N; m++)
	{
        for(n=0; n<N; n++)
		{
            fprintf(f, "\tC[%3d]=", m*N+n);
            for(k=0; k<N; k++)
			{
                if(pd[m*N+k]==1.0)		{ fprintf(f, "+       P[%3d]", k*N+n); op_add++; }
                else if(pd[m*N+k]==-1.0){ fprintf(f, "-       P[%3d]", k*N+n); op_add++; }
                else if(pd[m*N+k]!=0.0) { fprintf(f, "+A[%3d]*P[%3d]", m*N+k,k*N+n); op_mul++; op_add++; }
			}
            fprintf(f, ";\n");
        }
	}
    fprintf(f, "\n");
	// Py = C*A^T
    for(m=0; m<N; m++)       // 先计算上三角
	{
        for(n=m; n<N; n++)
		{
            fprintf(f, "\tB[%3d]=", m*N+n); 
            for(k=0; k<N; k++)
			{
                if(pd[n*N+k]==1.0)		{ fprintf(f, "+C[%3d]       ", m*N+k); op_add++; }
                else if(pd[n*N+k]==-1.0){ fprintf(f, "-C[%3d]       ", m*N+k); op_add++; }
                else if(pd[n*N+k]!=0.0) { fprintf(f, "+C[%3d]*A[%3d]", m*N+k, n*N+k);  op_mul++; op_add++; }
			} 
            fprintf(f, ";\n");
        }
	}
	for(m=1; m<N; m++)    // 再应用对称性获得P的下三角
	{
        for(n=0; n<m; n++)
            fprintf(f, "\tB[%3d]=B[%3d];", m*N+n, n*N+m);   
        fprintf(f, "\n");
    }
    fprintf(f, "\n\treturn;\n}\t//mul=%d,add=%d", op_mul, op_add);
	fclose(f);
}

void APATBuild1(char *fname, CMat &A)
{
	int N=A.row, N2=N*N, m, n, k, i, op_add=0, op_mul=0;
	double *pd = A.d;
	FILE *f=fopen(fname, "wt");
    fprintf(f, "#include \"navi.h\"\n\n");
    fprintf(f, "struct{ double\n\t");
	for(i=0; i<N2; i++) { fprintf(f, "A%d, ", i); if(i%N==N-1) fprintf(f, "\n\t"); };
	for(i=0; i<N2; i++) { fprintf(f, "P%d, ", i); if(i%N==N-1) fprintf(f, "\n\t"); };
	for(i=0; i<N2; i++) { fprintf(f, "C%d, ", i); if(i%N==N-1) fprintf(f, "\n\t"); };
	for(i=0; i<N2; i++) { fprintf(f, "B%d, ", i); if(i%N==N-1) fprintf(f, "\n\t"); };
	for(i=0; i<N; i++)	 fprintf(f, "X%d, ", i);	fprintf(f, "\n\t");
 	for(i=0; i<N-1; i++) fprintf(f, "Y%d, ", i);	fprintf(f, "Y%d;\n} s;\n", i);

    fprintf(f, "void APAT1(CMat &mA, CMat &mPx, CMat &mPy, CMat &vX, CMat &vY)\n{\n");
    fprintf(f, "\tint size1=%d, size2=%d;\n\n", N*sizeof(double), N2*sizeof(double));
    fprintf(f, "\tmemcpy(&s.A0, mA.d, size2);");
    fprintf(f, "\tmemcpy(&s.P0, mPx.d, size2);");
    fprintf(f, "\tmemcpy(&s.X0, vX.d, size1);\n\n");
	// Y = A*X
    for(m=0; m<N; m++)
	{
		fprintf(f, "\ts.Y%-3d=", m);
		for(n=0; n<N; n++)
		{
			if(pd[m*N+n]==1.0)		{ fprintf(f, "+     s.X%-3d", n); op_add++; }
			else if(pd[m*N+n]==-1.0){ fprintf(f, "-     s.X%-3d", n); op_add++; }
			else if(pd[m*N+n]!=0.0) { fprintf(f, "+s.A%-3d*s.X%-3d", m*N+n,n); op_mul++; op_add++; }
		}
		fprintf(f, ";\n");
	}
    fprintf(f, "\n");
	// C = A*Px
    for(m=0; m<N; m++)
	{
        for(n=0; n<N; n++)
		{
            fprintf(f, "\ts.C%-3d=", m*N+n);
            for(k=0; k<N; k++)
			{
                if(pd[m*N+k]==1.0)		{ fprintf(f, "+     s.P%-3d", k*N+n); op_add++; }
                else if(pd[m*N+k]==-1.0){ fprintf(f, "-     s.P%-3d", k*N+n); op_add++; }
                else if(pd[m*N+k]!=0.0) { fprintf(f, "+s.A%-3d*s.P%-3d", m*N+k,k*N+n); op_mul++; op_add++; }
			}
            fprintf(f, ";\n");
        }
	}
    fprintf(f, "\n");
	// Py = C*A^T
    for(m=0; m<N; m++)       // 先计算上三角
	{
        for(n=m; n<N; n++)
		{
            fprintf(f, "\ts.B%-3d=", m*N+n); 
            for(k=0; k<N; k++)
			{
                if(pd[n*N+k]==1.0)		{ fprintf(f, "+s.C%-3d     ", m*N+k); op_add++; }
                else if(pd[n*N+k]==-1.0){ fprintf(f, "-s.C%-3d     ", m*N+k); op_add++; }
                else if(pd[n*N+k]!=0.0) { fprintf(f, "+s.C%-3d*s.A%-3d", m*N+k, n*N+k);  op_mul++; op_add++; }
			} 
            fprintf(f, ";\n");
        }
	}
	for(m=1; m<N; m++)    // 再应用对称性获得P的下三角
	{
        for(n=0; n<m; n++)
            fprintf(f, "\ts.B%-3d=s.B%-3d;", m*N+n, n*N+m);   
        fprintf(f, "\n");
    }
    fprintf(f, "\n\tmemcpy(mPy.d, &s.B0, size2);");
    fprintf(f, "\tmemcpy(vY.d, &s.Y0, size1);\n");
    fprintf(f, "\n\treturn;\n}\t//mul=%d,add=%d", op_mul, op_add);
	fclose(f);
}
