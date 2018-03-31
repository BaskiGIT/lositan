#include <stdlib.h>
#include <stdio.h>
#include <math.h>

	double A_const=0.25,B_const=0.25;
	double Mono;
	double Critlim;


main()
{
	FILE *inp,*out,*pfile,*out_log,*out_ss;
	int ip,init2[2],*pfreq[2],byall,ng=0,*inside;
	int Subs,nloc,**freq_arr,i,j,jj,k,l,*noall,init[1000],ic,*ss,ig,smed,ginitot=0;
	int maxf,isum,too_fixed,*indx,**ssmat;
	float fsum,hsum,h0,h1,h2,fst,*pf,*ph;
	double *fvec,fmed,hcrit,p_trim,pcum,*hvec,hmed;
	char c, str1[10];

	inp = fopen("infile","r");
	out = fopen("data_fst_outfile","w");
	out_log = fopen("datacal_log.txt","w");
	out_ss = fopen("ss_file","w");
	pfile = fopen("pdist.dat","w");
	if(inp == 0){
		printf("no infile\n");
		exit(1);
	}
	fscanf(inp,"%d",&byall); /* 1 is by allele 0 is by population */
	while(!((c=getc(inp)) == '\n' || c == '\f' || c == '\r'));
	fscanf(inp,"%d",&Subs);
	while(!((c=getc(inp)) == '\n' || c == '\f' || c == '\r'));
	fscanf(inp,"%d",&nloc);
	while(!((c=getc(inp)) == '\n' || c == '\f' || c == '\r'));

	noall = (int *)malloc(nloc*sizeof(int));
	freq_arr = (int **)malloc(Subs*sizeof(int *));
	pf = (float *)malloc(Subs*(Subs-1)/2*sizeof(float ));
	ph = (float *)malloc(Subs*(Subs-1)/2*sizeof(float ));
	ss = (int *)malloc(Subs*nloc*sizeof(int));
	ssmat = (int **)malloc(nloc*sizeof(int*));
	for(j=0;j<nloc;++j)ssmat[j] = (int *)malloc(Subs*sizeof(int));
	fvec = (double *)malloc(nloc*sizeof(double));
	hvec = (double *)malloc(nloc*sizeof(double));
	inside = (int *)malloc(nloc*sizeof(int));
	fsum = hsum = 0.0;
	ig = 0;
	for(j=0;j<Subs*(Subs-1)/2;++j) ph[j] = pf[j] = 0.0;
	printf("Give critical frequency for the most common allele (globally)  ");
	scanf("%lf",&Critlim);
	if(Critlim >= 1.0){
		printf("Not a good idea to allow globally monomorphic loci: setting to 0.999\n");
		Critlim =0.999;
	}
	fprintf(out_log,"critical frequency for the most common allele (globally) is %f\n",Critlim);
	printf("give P for trimmed mean (weighted by het_between, 0.5 gives weighted median) ");
	scanf("%lf",&p_trim);
	if(p_trim > 0.5){
		printf("P is > 0.5; setting to 0.5\n");
		p_trim = 0.5;
	}
	fprintf(out_log,"P for trimmed mean (weighted by het_between, 0.5 gives weighted median): %f\n",p_trim);
	printf("give scale for Zhivotovsky estimate  suggest 0.25 and 0.25 ");
	scanf("%lf %lf",&A_const, &B_const);
	fprintf(out_log,"scale for Zhivotovsky estimate: %f %f\n",A_const, B_const);
	printf("infile contains %d loci and %d populations\n",nloc,Subs);
	fprintf(out_log,"infile contains %d loci and %d populations\n",nloc,Subs);
	for(j=0;j<nloc;++j){
		Mono = 0;
		fscanf(inp,"%d",&noall[j]);
		while(!((c=getc(inp)) == '\n' || c == '\f' || c == '\r'));
		for(k=0;k<Subs;++k)freq_arr[k] = (int *)malloc(noall[j]*sizeof(int));
		if(byall){
			for(k=0;k<noall[j];++k){
				for(l=0;l<Subs;++l){
					ic = fscanf(inp,"%d",&freq_arr[l][k]);
					if(ic <= 0 || ic == EOF){
						printf("error reading data\n");
						exit(1);
					}
				}
			}
		}
		else{
			for(k=0;k<Subs;++k){
				for(l=0;l<noall[j];++l){
					ic = fscanf(inp,"%d",&freq_arr[k][l]);
					if(ic == 0 || ic == EOF){
						printf("error reading data\n");
						exit(1);
					}
				}
			}
		}
		ginitot = 0;
		for(k=0;k<Subs;++k){
			init[k] = 0;
			for(i=0;i<noall[j];++i){
				init[k] += freq_arr[k][i];
			}
			ginitot += init[k];
		}


		for(i=0,maxf=0;i<noall[j];++i){
			for(k=0,isum=0;k<Subs;++k){
				isum += freq_arr[k][i];
			}
			if(isum > maxf)maxf = isum;
		}
		if((double)maxf/ginitot > Critlim){
			++too_fixed;
			Mono = 1;
		}
		if(!Mono){
			for(k=0;k<Subs;++k){
				ssmat[ng][k] = init[k];
				ss[ig++] = init[k];
			}
		}

		if(!Mono){
			my_thetacal(freq_arr,noall[j],init,Subs,&h0,&h1,&h2,&fst);
			fprintf(out,"%f %f %f %d\n",h2,fst,h1,j);
			fvec[ng] = fst;
			hvec[ng] = h1;
			++ng;
			fsum += fst*h1;
			hsum += h1;
		}
		ip = 0;
		for(i=0;i<Subs;++i){
			for(k=i+1;k<Subs;++k,++ip){
				pfreq[0] = freq_arr[i];
				pfreq[1] = freq_arr[k];
				init2[0] = init[i];
				init2[1] = init[k];
				if(init2[0] <= 1 || init2[1] <= 1)continue;
				for(jj=0,maxf=0;jj<noall[j];++jj){
					isum = pfreq[0][jj] + pfreq[1][jj];
					if(isum > maxf)maxf = isum;
				}
				if((double)maxf/(init2[0]+init2[1]) <= Critlim){
					my_thetacal(pfreq,noall[j],init2,2,&h0,&h1,&h2,&fst);
					pf[ip] += fst*h1;
					ph[ip] += h1;
				}
			}
		}
		for(k=0;k<Subs;++k)free(freq_arr[k]);
	}

	fprintf(out_ss,"%d\n",ng); /*ng is nloc - too_fixed */
	fprintf(out_ss,"%d\n\n",Subs);
	for(j=0;j<ng;++j){
		for(k=0;k<Subs;++k){
			fprintf(out_ss,"%d ",ssmat[j][k]);
		}
		fprintf(out_ss,"\n");
	}


	printf("(weighted) mean Fst is %f\n",fsum/hsum);
	fprintf(out_log,"(weighted) mean Fst is %f\n",fsum/hsum);
	for(j=0;j<ng;++j)hvec[j] /= hsum;
	if(ng == 0) fmed = -1;
	else{
		indx = (int *)malloc(ng*sizeof(int));
		for(j=0;j<ng;++j)indx[j] = j;
		dsorti('a',ng,fvec,indx);

		pcum = 0;
		for(j=0;j<ng;++j){
			pcum += hvec[indx[j]];
			if(pcum >= p_trim)inside[j] = 1;
			else inside[j] = 0;
		}
		pcum = 0;
		for(j=ng-1;j>=0;--j){
			pcum += hvec[indx[j]];
			if(pcum >= p_trim)inside[j] *= 1;
			else inside[j] *= 0;
		}
		fmed = hmed = 0;
		for(j=0;j<ng;++j){
			if(inside[j]){
				fmed += fvec[indx[j]]*hvec[indx[j]];
				hmed += hvec[indx[j]];
			}
		}
		fmed = fmed/hmed;


	}
	sort(ig,ss);
	if(ig%2 == 0)smed = (ss[ig/2-1] + ss[ig/2])/2.0 + 0.5;
	else smed = ss[ig/2];
	printf("median sample size over all loci and pops is %d\n",smed);
	fprintf(out_log,"median sample size over all loci and pops is %d\n",smed);
	if(ng == 0){
		printf("no observations sufficiently polymorphic\n");
		fprintf(out_log,"no observations sufficiently polymorphic\n");
	}
	else{
		printf("trimmed weighted mean F  %f based on %d observations\n",fmed,ng);
		fprintf(out_log,"trimmed weighted mean F  %f based on %d observations\n",fmed,ng);
	}
	for(i=0,ip=0;i<Subs;++i){
		for(k=i+1;k<Subs;++k,++ip){
			if(ph[ip] == 0)fprintf(pfile,"NA\n");
			else fprintf(pfile,"%f\n",pf[ip]/ph[ip]);
		}
	}
	printf("type in any character and return to close window  ");
	scanf("%s",str1);
}






thetacal(gen,noall,sample_size,no_of_samples,het0,het1,fst)

int *gen[],noall,sample_size[],no_of_samples;
float *het0,*het1,*fst;
{
	int i,j,*psum,ptot;
	double xx,yy,nbar,nc,q2,q3,nbar2;

	psum = (int *)malloc(noall*sizeof(int));

	for(i=0;i<noall;++i)psum[i] = 0;
	nbar = nbar2 = 0;
	for(j=0,xx=0.0;j<no_of_samples;++j) {
		nbar += sample_size[j];
		nbar2 += sample_size[j]*sample_size[j];
		for(i=0;i<noall;++i) {
			psum[i] += gen[j][i];
			if(sample_size[j] > 0)
			xx += ((double)gen[j][i]*gen[j][i])/(double)sample_size[j];
		}
	}
	nc = 1.0/(no_of_samples - 1.0)*(nbar - nbar2/nbar);
	nbar /= no_of_samples;

	for(i=0,yy=0.0,ptot = 0;i<noall;++i) {
		yy += psum[i]*psum[i];
		ptot += psum[i];
	}
	q2 = (xx-no_of_samples)/(no_of_samples*(nbar - 1.0));
	q3 = 1.0/(no_of_samples*(no_of_samples-1.0)*nbar*nc)*
			(yy - nbar*(nc-1.0)/(nbar-1.0)*xx) +
			(nbar-nc)/(nc*(nbar-1.0))*(1.0-1.0/
			(no_of_samples - 1.0)*xx);

	*het0 = 1.0 - q2;
	*het1 = 1.0 - q3;
	if(*het1 < 1.0e-10){
		*fst = -100.0;
		printf("thetacal: fst undefined because of zero heterozygosity. Reduce threshold to define locus as monomorphic");
		exit(1);
	}
	else *fst = 1.0 - (*het0)/(*het1);

	free(psum);

}


my_thetacal(gen,noall,sample_size,no_of_samples,het0,het1,het3,fst)

int *gen[],noall,sample_size[],no_of_samples;
float *het0,*het1,*het3,*fst;
{
	int i,j,k,skip,tot_gen,tot_samp;
	double x2,x0,yy,y1,q2,q3,freq[1000], var[1000],tot_freq,var_tot_freq,x3;
	double a_const=1,b_const=1;

	a_const = A_const;
	b_const = B_const;
	x0 = 0.0;
	skip = 0;
	tot_gen = 0;
	tot_samp = 0;
	for(j=0;j<no_of_samples;++j){
		if(sample_size[j] == 0){++skip;continue;}
		tot_gen += gen[j][0];
		tot_samp += sample_size[j];
		freq[j] = exp(lgamma(gen[j][0]+a_const+0.5)+lgamma(sample_size[j]+a_const+b_const)-
		             lgamma(gen[j][0]+a_const)-lgamma(sample_size[j]+a_const+b_const+0.5));
		             /* from zhivotovsky 10 */
		var[j] =exp(lgamma(gen[j][0]+a_const+1) + lgamma(sample_size[j]+a_const+b_const) -
		             lgamma(gen[j][0]+a_const)-lgamma(sample_size[j]+a_const+b_const+1)) - freq[j]*freq[j];
		             /* from zhivotovsky 11*/
	}
	for(j=0;j<no_of_samples;++j) {
		if(sample_size[j] == 0){continue;}
		x2 = 2*freq[j]*(1-freq[j]) + 2*var[j]; /*from lynch and milligan 4a */
		x0 += 1-x2; /*convert heterozygosity to homozygosity */
	}

	yy = 0.0;
	for(j=0;j<no_of_samples;++j){
		if(sample_size[j] == 0)continue;
		for(k=j+1;k<no_of_samples;++k){
			if(sample_size[k] == 0)continue;

			yy += freq[j]*freq[k] + (1-freq[j])*(1-freq[k]); /* lynch and milligan 9a */
		}
	}



	q2 = x0/(no_of_samples-skip);
	q3 = 2*yy/((no_of_samples-skip)*(no_of_samples-skip-1));


	*het0 = 1.0 - q2;
	*het1 = 1.0 - q3;
	if(*het1 < 1.0e-10){
		*fst = -100.0;
		printf("thetacal: fst undefined because of zero heterozygosity. Reduce threshold to define locus as monomorphic");
		exit(1);
	}
	else *fst = 1.0 - (*het0)/(*het1);


	tot_freq = exp(lgamma(tot_gen+a_const+0.5)+lgamma(tot_samp+a_const+b_const)-
		             lgamma(tot_gen+a_const)-lgamma(tot_samp+a_const+b_const+0.5));
		             /* from zhivotovsky 10 */
	var_tot_freq =exp(lgamma(tot_gen+a_const+1) + lgamma(tot_samp+a_const+b_const) -
		             lgamma(tot_gen+a_const)-lgamma(tot_samp+a_const+b_const+1)) - tot_freq*tot_freq;
		             /* from zhivotovsky 11*/
	x3 = 2*tot_freq*(1-tot_freq) + 2*var_tot_freq;
	*het3 = x3;


}



sort(n,ra)
int n;
int ra[];

{
	int i,j,l,ir;
	int rra;

	if(n == 1)return;
	l = n/2+1;
	ir = n;
	while(1){
		if(l>1){
			--l;
			rra = ra[l-1];
		}
		else{
			rra = ra[ir-1];
			ra[ir-1] = ra[0];
			--ir;
			if(ir == 0){
				ra[0] = rra;
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
				i = j;
				j += j;
			}
			else j = ir+1;
		}
		ra[i-1] = rra;
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

dsorti(dir, n, x, indx)
char dir;
int n;
double *x;
int *indx;
{
	int i, j, h, asc,indtmp;
	double xtmp,*priv;

	priv = (double *)malloc(n*sizeof(double));
	for(j=0;j<n;++j)priv[j] = x[j];
	if(dir == 'a' || dir == 'A')asc = 1;
	else asc = 0;
	for(j=0;j<n;++j)indx[j] = j;
	h = 1;
	do {
		h = 3 * h + 1;
	}
	while (h <= n);

	do {
		h = h / 3;
		for (i = h; i < n; i++) {
			xtmp = priv[i];
			indtmp = indx[i];
			j = i;
			if(asc){
				while (priv[j - h] > xtmp) {
					priv[j] = priv[j - h];
					indx[j] = indx[j-h];
					j = j - h;
					if (j < h)
						goto end;
				}
			}
			else{
				while (priv[j - h] < xtmp) {
					priv[j] = priv[j - h];
					indx[j] = indx[j-h];
					j = j - h;
					if (j < h)
						goto end;
				}
			}

		end:	priv[j] = xtmp;indx[j] = indtmp;
		}
	} while (h != 1);
	free(priv);
}



/* Inserted from R0.16::lgamma.c */

/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Modification by MAB 20.3.97  removed #include "Mathlib.h"
inserted definition of M_PI, copied from Mathlib.h */
#undef M_PI
#define M_PI		3.141592653589793238462643383276

int signgamR1 = 0;

/* log(2*pi)/2 and pi */

static double hl2pi = 0.9189385332046727417803297;
static double xpi = M_PI;

 /* Coefficients from Cheney and Hart */

#define M 6
#define N 8
static double p1[] =
{
	0.83333333333333101837e-1,
	-.277777777735865004e-2,
	0.793650576493454e-3,
	-.5951896861197e-3,
	0.83645878922e-3,
	-.1633436431e-2,
};
static double p2[] =
{
	-.42353689509744089647e5,
	-.20886861789269887364e5,
	-.87627102978521489560e4,
	-.20085274013072791214e4,
	-.43933044406002567613e3,
	-.50108693752970953015e2,
	-.67449507245925289918e1,
	0.0,
};
static double q2[] =
{
	-.42353689509744090010e5,
	-.29803853309256649932e4,
	0.99403074150827709015e4,
	-.15286072737795220248e4,
	-.49902852662143904834e3,
	0.18949823415702801641e3,
	-.23081551524580124562e2,
	0.10000000000000000000e1,
};

static double posarg(double);
static double negarg(double);
static double asform(double);

double lgamma(double arg)
{
	signgamR1 = 1.0;
	if (arg <= 0.0)
		return (negarg(arg));
	if (arg > 8.0)
		return (asform(arg));
	return (log(posarg(arg)));
}

/* Equation 6.1.41 Abramowitz and Stegun */
/* See also ACM algorithm 291 */

static double asform(arg)
double arg;
{
	double log();
	double n, argsq;
	int i;

	argsq = 1. / (arg * arg);
	for (n = 0, i = M - 1; i >= 0; i--) {
		n = n * argsq + p1[i];
	}
	return ((arg - .5) * log(arg) - arg + hl2pi + n / arg);
}

static double negarg(arg)
double arg;
{
	double temp;
	double log(), sin(), posarg();

	arg = -arg;
	temp = sin(xpi * arg);
	if (temp == 0.0)
		/* removed DOMAIN_ERROR  printerr("negarg: temp == 0.0");   */
        return 700.0; /* this is a bad bodge to keep things running */
	if (temp < 0.0)
		temp = -temp;
	else
		signgamR1 = -1;
	return (-log(arg * posarg(arg) * temp / xpi));
}

static double posarg(arg)
double arg;
{
	double n, d, s;
	register i;

	if (arg < 2.0)
		return (posarg(arg + 1.0) / arg);
	if (arg > 3.0)
		return ((arg - 1.0) * posarg(arg - 1.0));

	s = arg - 2.;
	for (n = 0, d = 0, i = N - 1; i >= 0; i--) {
		n = n * s + p2[i];
		d = d * s + q2[i];
	}
	return (n / d);
}

/* end of lgamma insertion */

