#include <stdlib.h>
#include <stdio.h>
#include <math.h>



main()
{
	FILE *inp,*out,*pfile;
	int ip,init2[2],*pfreq[2],byall;
	int Subs,nloc,**freq_arr,i,j,k,l,noall[90000],init[90000],ic,*ss,ig,smed;
	float fsum,hsum,h0,h1,fst,*pf,*ph;
	char c, str1[10];
	
	inp = fopen("infile","r");
	out = fopen("data_fst_outfile","w");
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
	
	freq_arr = (int **)malloc(Subs*sizeof(int *));
	pf = (float *)malloc(Subs*(Subs-1)/2*sizeof(float *));
	ph = (float *)malloc(Subs*(Subs-1)/2*sizeof(float *));
	ss = (int *)malloc(Subs*nloc*sizeof(int));
	fsum = hsum = 0.0;
	ig = 0;
	for(j=0;j<Subs*(Subs-1)/2;++j) ph[j] = pf[j] = 0.0;
	for(j=0;j<nloc;++j){
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
		for(k=0;k<Subs;++k){
			init[k] = 0;
			for(i=0;i<noall[j];++i){
				init[k] += freq_arr[k][i];
			}
			ss[ig++] = init[k];
		}
		my_thetacal(freq_arr,noall[j],init,Subs,&h0,&h1,&fst);
		fprintf(out,"%f %f\n",h1,fst);
		if(fst > -10.0){
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
				my_thetacal(pfreq,noall[j],init2,2,&h0,&h1,&fst);
				if(fst > -10.0){
					pf[ip] += fst*h1;
					ph[ip] += h1;
				}
			}
		}
		for(k=0;k<Subs;++k)free(freq_arr[k]);
	}
	printf("(weighted) mean Fst is %f\n",fsum/hsum);
	sort(Subs*nloc,ss);
	if(ig%2 == 0)smed = (ss[ig/2-1] + ss[ig/2])/2.0 + 0.5;
	else smed = ss[ig/2];
	printf("median sample size over all loci and pops is %d\n",smed);
	for(i=0,ip=0;i<Subs;++i){
		for(k=i+1;k<Subs;++k,++ip){
			fprintf(pfile,"%f\n",pf[ip]/ph[ip]);
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
	if(*het1 < 1.0e-10)*fst = -100.0;
	else *fst = 1.0 - (*het0)/(*het1);
	
	free(psum);
	
}

my_thetacal(gen,noall,sample_size,no_of_samples,het0,het1,fst)
	
int *gen[],noall,sample_size[],no_of_samples;
float *het0,*het1,*fst;
{
	int i,j,k,skip;
	double x2,x0,yy,y1,q2,q3;
	
	
	x0 = 0.0;
	skip = 0;
	for(j=0;j<no_of_samples;++j) {
		if(sample_size[j] == 0){++skip;continue;}
		for(i=0,x2=0.0;i<noall;++i) {
			x2 += gen[j][i]*gen[j][i];
		}
		x0 += (x2-sample_size[j])/(sample_size[j]*(sample_size[j]-1));
	}
	
	yy = 0.0;
	for(j=0;j<no_of_samples;++j){
		if(sample_size[j] == 0)continue;
		for(k=j+1;k<no_of_samples;++k){
			if(sample_size[k] == 0)continue;
			for(i=0,y1=0.0;i<noall;++i){
				y1 += gen[j][i]*gen[k][i];
			}
			yy += y1/(sample_size[j]*sample_size[k]);
		}
	}
	
	
	
	q2 = x0/(no_of_samples-skip);
	q3 = 2*yy/((no_of_samples-skip)*(no_of_samples-skip-1));
	
			
	*het0 = 1.0 - q2;
	*het1 = 1.0 - q3;
	if(*het1 < 1.0e-10)*fst = -100.0;
	else *fst = 1.0 - (*het0)/(*het1);
	
	
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
