/* 3.10.97 edited to remove static arrays */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.141592653589793
#define WIN 1000
#define BLOCK WIN/2

void isort();
int 	rand_table[98],jindic;
float gfsr4();
void opengfsr(),closegfsr();
void mom();
float qnorm();

int NMAX;

main()
{
	FILE *input,*out;
	int j,ic,n,*index,*index2,k,i;
	float *x0,*y0,*x1,*y1;
	float *xmed,*ymed,*ytop,*ybot,*y2;
	float mean,sdev,skew,kurt,min,max,g1,d1,l1,xi1,nvar,qval;
	float dd1, dd2;
	int ifault,itype1;
	char instr[100],outstr[100];
	
	
	opengfsr();
	printf("infile, outfile ?:");
	scanf("%s %s",instr,outstr);
	input = fopen(instr,"r");
	out = fopen(outstr,"w");
	for(j=0;;++j){
		ic = fscanf(input,"%f %f ",&dd1,&dd2);
		if(ic == EOF)break;
	}

	NMAX = n = j;
	fclose(input);
	
	index = (int *)malloc(NMAX*sizeof(int));
	index2 = (int *)malloc(NMAX*sizeof(int));
	x0 = (float *)malloc(NMAX*sizeof(float));
	y0 = (float *)malloc(NMAX*sizeof(float));
	x1 = (float *)malloc(NMAX*sizeof(float));
	y1 = (float *)malloc(NMAX*sizeof(float));
	xmed = (float *)malloc(NMAX*sizeof(float));
	ymed = (float *)malloc(NMAX*sizeof(float));
	ytop = (float *)malloc(NMAX*sizeof(float));
	ybot = (float *)malloc(NMAX*sizeof(float));
	y2 = (float *)malloc(NMAX*sizeof(float));
	
	
	input = fopen(instr,"r");

	for(j=0;j<n;++j){
		ic = fscanf(input,"%f %f ",&x0[j],&y0[j]);
	}
	
	printf("probability for contour? (e.g. 0.95 or 0.99)   ");
	scanf("%f", &qval);
	qval = qnorm(0.5*(1.0-qval));
	isort(n,x0,index);
	for(j=0;j<n;++j){
		x1[j] = x0[index[j]];
		y1[j] = y0[index[j]];
	}
	for(k=0; k+WIN-1 < n;k+=BLOCK){
		xmed[k] = (x1[WIN/2-1+k]+x1[WIN/2+k])/2.0;
		isort(WIN,&y1[k],index2);
		mom(&y1[k],WIN,&mean,&sdev,&skew,&kurt,&min,&max);
		jnsn(&mean,&sdev,&skew,&kurt,&itype1,&g1,&d1,&l1,&xi1,&ifault); 
/*		printf("%f %f %f %f %d %f %f %f %f %d\n",
		mean,sdev,skew,kurt,itype1,g1,d1,l1,xi1,ifault); */
		if(ifault!=0){
			printf(" jnsn. ifault is %d\n",ifault);
			if(ifault != 3)continue;
		}
		nvar = -qval;
		varj(&ybot[k],&nvar,&itype1,&g1,&d1,&l1,&xi1,&ifault);
/*		printf("%f %f %d %f %f %f %f%d\n",ybot[k],nvar,
		itype1,g1,d1,l1,xi1,ifault); */
		if(ifault!=0){
			printf("varj 1 ifault is %d\n",ifault);
			continue;
		}
		nvar = 0.0;
		varj(&ymed[k],&nvar,&itype1,&g1,&d1,&l1,&xi1,&ifault);
/*		printf("%f %f %d %f %f %f %f%d\n",ymed[k],nvar,
		itype1,g1,d1,l1,xi1,ifault); */
		if(ifault!=0){
			printf("varj 2 ifault is %d\n",ifault);
			continue;
		}
		nvar = qval;
		varj(&ytop[k],&nvar,&itype1,&g1,&d1,&l1,&xi1,&ifault);
/*		printf("%f %f %d %f %f %f %f%d\n",ytop[k],nvar,
		itype1,g1,d1,l1,xi1,ifault); */
		if(ifault!=0){
			printf("varj 3 ifault is %d\n",ifault);
			continue;
		}
		fprintf(out,"%f %f %f %f\n",xmed[k],ybot[k],ymed[k],ytop[k]);
	}
	closegfsr();
}
		
float qnorm(p)
float p;
{
	float a0=2.30753,a1=0.27061,b1=0.99229,b2=0.04481;
	if(p>0.5)p=1.0-p;
	p = sqrt(log(1.0/(p*p)));
	return p-(a0+a1*p)/(1+b1*p+b2*p*p);
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
	if(rt==NULL) {
		printf("I need INTFILE! Where is INTFILE?\n");
		exit(1);
	}
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


