/* 3.10.97 edited to remove static arrays */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.141592653589793

void isort();
int 	rand_table[98],jindic;
float gfsr4();
void opengfsr(),closegfsr();
void mom();
float qnorm(),pnorm();
int NMAX;
void qcalc();


main()
{
	FILE *input,*out,*dataf;
	int j,ic,n,*index,k,i,ndata,dd1,dd2;
	float *x0,*y0,*x1,*y1;
	float *xmed,*ymed,*ytop,*ybot,*y2;
	float dx0[90000],dy0[90000],*dx1,*dy1,dpval[90000],dnval[90000],dum;
	float mean,sdev,skew,kurt,min,max,g1,d1,l1,xi1,nvar,qval;
	double ppx,ppx2,ppx3,pval,winprop;
	int ifault,itype1,idum,iwin;
	char instr[200],outstr[200],datstr[200];

	printf("data file, output file, simulated density file?:");
	scanf("%s %s %s",datstr,outstr,instr);
	dataf = fopen(datstr,"r");
	input = fopen(instr,"r");
	out = fopen(outstr,"w");
	for(j=0;;++j){
		ic = fscanf(dataf,"%f %f %f %d",&dx0[j],&dy0[j],&dum,&idum);
		if(ic == EOF)break;
	}
	ndata = j;
	dx1 = (float *)malloc(ndata*sizeof(float));
	dy1 = (float *)malloc(ndata*sizeof(float));



	for(j=0;;++j){
		ic = fscanf(input,"%f %f ",&dd1,&dd2);
		if(ic == EOF)break;
	}


	NMAX = n = j;
	fclose(input);

	printf("give p-value (e.g. 0.99)  ");
	scanf("%lf",&pval);
	pval = 1.0 - pval;

	printf("give smoothing proportion (e.g. 0.04)   ");
	scanf("%lf",&winprop);
	iwin = winprop*n;

	index = (int *)malloc(NMAX*sizeof(int));
	x0 = (float *)malloc(NMAX*sizeof(float));
	y0 = (float *)malloc(NMAX*sizeof(float));
	x1 = (float *)malloc(NMAX*sizeof(float));
	y1 = (float *)malloc(NMAX*sizeof(float));

	input = fopen(instr,"r");

	for(j=0;j<n;++j){
		ic = fscanf(input,"%f %f ",&x0[j],&y0[j]);
	}



	isort(n,x0,index);
	for(j=0;j<n;++j){
		x1[j] = x0[index[j]];
		y1[j] = y0[index[j]];
	}
	isort(ndata,dx0,index);
	for(j=0;j<ndata;++j){
		dx1[j] = dx0[index[j]];
		dy1[j] = dy0[index[j]];
	}

	for(j=0;j<ndata;++j) {
		qcalc(dx1[j],dy1[j],x1,y1,pval,n,iwin,&ppx,&ppx2,&ppx3);
		fprintf(out,"%f %f %f %f %f\n",dx1[j],dy1[j],ppx,ppx2,ppx3);
	}

}


void qcalc(het,fst,hvec,fvec,pp,len,win,pval,pval2,pval3)
float het;
float fst;
float hvec[];
float fvec[];
int len;
int win;
double *pval,*pval2,*pval3,pp;

{
	double *dummy,pvala,pvalb;
	int j,i,ic,winstart,winend,ipos,extra,iq;
	double r,pv_l,pv_u,med,f,q_l,q_med,q_u;

	dummy = (double *)malloc((win+5)*sizeof(double));
	for(j=0;j<len;++j){
		if(hvec[j] > het)break;
	}
	if(j==len)j=len-1;
	ipos = j;
	winstart = ipos - win/2;
	winend = ipos+win/2;
	if(winstart < 0){
		extra = -winstart;
		winstart = 0;
		winend = ipos + win/2 + extra;
		if(winend >= len)winend = len-1;
	}
	if(winend >= len){
		extra = winend - len;
		winend = len-1;
		winstart = ipos - win/2 - extra;
		if(win < 0)winstart = 0;
	}

	ic = 0;
	for(j=winstart,i=0;j<=winend;++j,++i){
		dummy[i] = fvec[j];
		++ic;
	}
	dsort('a',ic,dummy);
	pv_l = pp*0.5;
	pv_u = 1-pp*0.5;
	med = 0.5;

	r = 1+(ic-1)*pv_l;
	iq = r;
	f = r-iq;
	if(iq-1 >= ic){printf("iq-1 >= ic\n");exit(1);}
	q_l = (1.0-f)*dummy[iq-1]+f*(iq == ic ? dummy[iq-1] : dummy[iq]);

	r = 1+(ic-1)*med;
	iq = r;
	f = r-iq;
	if(iq-1 >= ic){printf("iq-1 >= ic\n");exit(1);}
	q_med = (1.0-f)*dummy[iq-1]+f*(iq == ic ? dummy[iq-1] : dummy[iq]);

	r = 1+(ic-1)*pv_u;
	iq = r;
	f = r-iq;
	if(iq-1 >= ic){printf("iq-1 >= ic\n");exit(1);}
	q_u = (1.0-f)*dummy[iq-1]+f*(iq == ic ? dummy[iq-1] : dummy[iq]);

	*pval = q_l;
	*pval2 = q_med;
	*pval3 =q_u;


	free(dummy);
}

float qnorm(p)
float p;
{
	static float a0=2.30753,a1=0.27061,b1=0.99229,b2=0.04481;
	if(p>0.5)p=1.0-p;
	p = sqrt(log(1.0/(p*p)));
	return p-(a0+a1*p)/(1+b1*p+b2*p*p);
}

float pnorm(x)
float x;
{
	static double 	p=0.2316419,
					b1=0.319381530,
					b2=-0.356563782,
					b3=1.781477937,
					b4=-1.821255978,
					b5=1.330274429;
	double t,z;

	if(x<-5.0)return 0.0;
	if(x>5.0)return 1.0;
	if(x<0.0) {
		x = -x;
		t = 1.0/(1.0+p*x);
		z = exp(-(x*x)*0.5)/sqrt(2.0*PI);
		return z*(t*b1+t*(b2+t*(b3+t*(b4+t*b5))));
	}
	t = 1.0/(1.0+p*x);
	z = exp(-(x*x)*0.5)/sqrt(2.0*PI);
	return 1.0-z*(t*b1+t*(b2+t*(b3+t*(b4+t*b5))));
}



void mom(x,n,x1,x2,x3,x4,min,max)
int n;
float x[],*x1,*x2,*x3,*x4,*min,*max;
{
	int i;
	float s1,s2,s3,s4,an,an1,dx,dx2,xi,var;

	s1 = x[0];
	s2 = 0.0;
	s3 = 0.0;
	s4 = 0.0;
	*min = s1;
	*max = s1;
	for(i=1;i<n;++i){
		xi = x[i];
		an = i+1;
		an1 = i;
		dx = (xi-s1)/an;
		dx2 = dx*dx;
		s4 -= dx*(4.0*s3-dx*(6.0*s2+an1*(1.0+pow(an1,3.0))*dx2));
		s3 -= dx*(3.0*s2-an*an1*(an-2.0)*dx2);
		s2 += an*an1*dx2;
		s1 += dx;
		if(xi<*min)*min=xi;
		if(xi>*max)*max=xi;
	}
	*x1 = s1;
	var = n>1 ? s2/(n-1) : 0.0;
	*x2 = sqrt(var);
	*x3 = var>0.0 ? 1.0/(n-1)*s3/pow(var,1.5) : 0.0;
	*x4 = var>0.0 ? 1.0/(n-1)*s4/pow(var,2.0) : 0.0;
	return;
}


void isort(n,inar,index)
int n,index[];
float inar[];

{
	int i,j,l,rrb,ir;
	float rra;
	float *ra;
	ra = (float *)malloc(n*sizeof(float));
	for(i=0;i<n;++i){
		ra[i] = inar[i];
		index[i] = i;
	}
	if(n == 1)return;
	l = n/2+1;
	ir = n;
	while(1){
		if(l>1){
			--l;
			rra = ra[l-1];
			rrb = index[l-1];
		}
		else{
			rra = ra[ir-1];
			rrb = index[ir-1];
			ra[ir-1] = ra[0];
			index[ir-1] = index[0];
			--ir;
			if(ir == 0){
				ra[0] = rra;
				index[0] = rrb;
			/*	for(i=1;i<n;++i){
					if(ra[i]<ra[i-1]){
						printf("error in sort\n");
						exit(1);
					}
				} */
				free(ra);
				return;
			}
		}
		i = l;
		j = l+l;
		while(j<=ir){
			if(j < ir){
				if(ra[j-1] < ra[j])++j;
			}

			if(rra < ra[j-1]){
				ra[i-1] = ra[j-1];
				index[i-1] = index[j-1];
				i = j;
				j += j;
			}
			else j = ir+1;
		}
		ra[i-1] = rra;
		index[i-1] = rrb;
	}
}

dsort(dir,n,x)  /* This is adapted from R 0.16 */
char dir;
int n;
double *x;
{
	int i, j, h, asc;
	double xtmp;

	if(dir == 'a' || dir == 'A')asc = 1;
	else asc = 0;
	h = 1;
	do {
		h = 3 * h + 1;
	}
	while (h <= n);

	do {
		h = h / 3;
		for (i = h; i < n; i++) {
			xtmp = x[i];
			j = i;
			if(asc){
				while (x[j - h] > xtmp) {
					x[j] = x[j - h];
					j = j - h;
					if (j < h)
						goto end;
				}
			}
			else{
				while (x[j - h] < xtmp) {
					x[j] = x[j - h];
					j = j - h;
					if (j < h)
						goto end;
				}
			}

		end:	x[j] = xtmp;
		}
	} while (h != 1);
}



int rbit()
{
	static int j = 32;
	static unsigned int orran;
	orran >>= 1;
	++j;
	if(j>32){
		orran = intrand();
		j = 1;
	}
	return orran & 1;
}

int rbit32()
{
	static int j = 32;
	static unsigned int orran;
	orran >>= 5;
	j += 5;
	if(j>32){
		orran = intrand();
		j = 5;
	}
	return orran & 31;
}


int intrand()

{
      int k;
      ++jindic;
      if(jindic>97)jindic=0;
      k = jindic+27;
      if(k>97)k=k-98;
      rand_table[jindic] = rand_table[jindic]^rand_table[k];
      return(rand_table[jindic]);
}

int disrand(l,t)
int l,t;
{
      int k;
      if(t<l){
      	printf("error in disrand\n");
      	exit(1);
      }
      ++jindic;
      if(jindic>97)jindic=0;
      k = jindic+27;
      if(k>97)k=k-98;
      rand_table[jindic] = rand_table[jindic]^rand_table[k];
	return((unsigned)rand_table[jindic]%(t-l+1)+l);
}

float gfsr4()
{
      int k;
      ++jindic;
      if(jindic>97)jindic=0;
      k = jindic+27;
      if(k>97)k=k-98;
      rand_table[jindic] = rand_table[jindic]^rand_table[k];
	return((unsigned)rand_table[jindic]/4294967296.0);
}




float gammln(xx)
float xx;
{
	double 	x,tmp,ser;
	static double cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;

	x=xx-1.0;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.0;
	for (j=0;j<=5;j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}


float expdev()
{

      int k;
      ++jindic;
      if(jindic>97)jindic=0;
      k = jindic+27;
      if(k>97)k=k-98;
      rand_table[jindic] = rand_table[jindic]^rand_table[k];
	return(-log((unsigned)rand_table[jindic]/4294967296.0));
}

void opengfsr()
{

	FILE 	*rt,*fopen();
	int 	j;

	rt = fopen("INTFILE","r");
	for(j=0;j<98;++j)fscanf(rt,"%d",&rand_table[j]);
	fscanf(rt,"%d",&jindic);
	fclose(rt);
}

void closegfsr()
{
	FILE 	*rt,*fopen();
	int 	j;

	rt = fopen("INTFILE","w");
	for(j=0;j<98;++j)fprintf(rt,"%d\n",rand_table[j]);
	fprintf(rt,"%d\n",jindic);
	fclose(rt);
}


