#include "navi.h"

CBinFile::CBinFile(char *fname, int clm, int row)
{
	FILE *f = fopen(fname, "rb");
	fseek(f, 0, SEEK_END);
	int maxRow = ftell(f)/(clm*sizeof(double));
	this->clm = clm;
	this->row = min(row,maxRow);
	int sz = this->row*this->clm*sizeof(double);
	pd = (double*)malloc(sz);
	assert(pd!=NULL);
	fseek(f, 0, SEEK_SET);
	printf("Loading file '%s' ...\n", fname);
	size = fread(pd, 1, sz, f);
	assert(size==sz);
	fclose(f);
}

CBinFile::~CBinFile()
{
	if(pd) delete pd, pd=NULL;
}

double* CBinFile::operator()(int m, int n)
{
	int idx = m*clm+n;
	return (idx<size && idx>=0) ? &pd[idx] : NULL;
}

