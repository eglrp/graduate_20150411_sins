#include "navi.h"

CFload::CFload()
{
}

CFload::CFload(char *fname, int clm, char separator, char endChar)
{
	// In each line of the load file, 'separator' or 'endChar' 
	// must tightly follow each numeric data.
	f = fopen(fname, "rt");
	format = NULL, d = NULL;
	if(f)
	{
		if(clm==0)
		{
			char format_tmp[5]="%lf,"; format_tmp[3] = separator;
			while(1){
				char c;
				fread(&c, 1, 1, f);
				if(c=='\n')	break;
			};
			long pos = ftell(f), pos0=0, pos1;
			rewind(f);
			for(; ; clm++)
			{
				double dtmp;
				fscanf(f, format_tmp, &dtmp);
				pos1 = ftell(f);
				if(pos1==pos0) break;
				if(pos1>=pos) {clm++; break;}
				pos0 = pos1;
			};
			rewind(f);
		}

		this->clm = clm;
		format = new char[5*clm];
		d = new double[clm];
		for(int i=0; i<clm; i++) 
		{
			format[5*i+0] = '%';
			format[5*i+1] = 'l';
			format[5*i+2] = 'f';
			format[5*i+3] = separator;
			format[5*i+4] = '\0';
			d[i] = 0;
		}
		format[5*(i-1)+3] = endChar;
	}
}


CFload::~CFload()
{
	if(f) fclose(f), f=NULL;
	if(format) delete format, format=NULL;
	if(d) delete d, d=NULL;
}

double CFload::operator()(int n)
{
	return d[n];
}

CVect3 CFload::operator[](int n)
{
	return CVect3(d[n],d[n+1],d[n+2]);
}

long CFload::Skip(long lines)
{
	long ln=0;
	if(lines<1) return ln;
	do{
		char c;
		fread(&c, 1, 1, f);		// for read fast, no use 'fscanf'
		if(c=='\n')	ln ++;
		if(feof(f))	break;
	}while(ln<lines);
	return ln;
}

BOOL CFload::Load(void)
{
	for(int i=0; i<clm; i++)
		fscanf(f, &format[5*i], &d[i]);
	return !feof(f);
}
