#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NOC 100 /* The maximum number of subpopulations possible. Spno - the
				number of subpopulations in the model - must
				be less than or equal to this */
#define PI 3.141592653589793



struct node{
	int sp;
	int osp;
	double time;
	int I;
	int cut;
	int dna;
	struct node *a[2];
	struct node *d[2];
};


typedef struct node *F;

int 	rand_table[98],jindic;

int  Spno;

double  Sd[NOC],Migrate[NOC],Den[NOC][3],Dtop;
int	Ni[NOC],N_n,Ntot;
int	Occ,Occlist[NOC];
float 	Tt;
struct node **List[NOC],**Nlist;
int NEXTMUT,Ms=0;
int Subs,Lmax[NOC];
int Increment1=50,Nmax=50;

void 	opengfsr(),closegfsr(),cpress(),treefill(),killtree();
float gfsr4(),expdev();
int disrand();
double Critlim=0.98;
	double A_const=0.2,B_const=0.2;

main()
{
	double mu,next_mu,mrate,dd,theta;
	float	rr,efst,fsum,wsum;
	float	tcum;
	int	j,jj,j2,i,itno,it,seq,ic,k,ik,ll,initot,init[NOC],ginitot,ginit[NOC],keepmu;
	double rm;
	int i1,j1,kk,**val_arr,**freq_arr,**genotype,noall,in_smp,match_samp,check_subs,lpick;
	int **ss;
	float	vhet,h0,h1,h2,fst;
	int	totall;
	int isum,maxf,mono,gt_two,too_fixed,nloc;
	char dic[100],str1[10];
	FILE *inp,*alls,*out_log,*infl;
	int gt1,gt2,irec,ig;


	opengfsr();
	out_log = fopen("fdist_log.txt","w");
	inp = fopen("Dfdist_params","r");
	fscanf(inp,"%d",&Spno);
	if(Spno > NOC){
		printf("error in parameter file - ");
		printf("number of subpopulation greater than %d\n",NOC);
		exit(1);
	}
	fscanf(inp,"%d",&Subs);
	fscanf(inp,"%f",&efst);
	fscanf(inp,"%d",&in_smp);
	fscanf(inp,"%lf",&theta);
	fscanf(inp,"%d",&itno);
	fscanf(inp,"%lf %lf",&A_const,&B_const);
	fscanf(inp,"%lf",&Critlim);
	fclose(inp);
	printf("number of demes is %d\n",Spno);
	fprintf(out_log,"number of demes is %d\n",Spno);

	printf("number of samples  is %d\n",Subs);
	fprintf(out_log,"number of samples  is %d\n",Subs);

	printf("expected value of W&C estimator is %f\n",efst);
	fprintf(out_log,"expected value of W&C estimator is %f\n",efst);

	efst = 1.0/(1.0+(Spno-1.0)/Spno*(1/efst-1)); /* 15.08.07 now fixed the Subs bug, noted earlier */

	printf("given %d demes, Fst of %f will be used \n",Spno,efst);
	fprintf(out_log,"given %d demes, Fst of %f will be used \n",Spno,efst);

	if(in_smp == 0){
		printf("Sample size is set to 0, therefore taking from file\n");
		fprintf(out_log,"Sample size is set to 0, therefore taking from file\n");
		match_samp = 1;
	}
	else{
		printf("Sample size used is %d \n",in_smp);
		fprintf(out_log,"Sample size uses is %d \n",in_smp);
		match_samp = 0;
	}

	printf("theta (=2*%d*N*mu, where N is deme size) is %f\n",Spno,theta);
	fprintf(out_log,"theta (=2*%d*N*mu, where N is deme size) is %f\n",Spno,theta);

	printf("%d realizations (loci) requested\n",itno);
	fprintf(out_log,"%d realizations (loci) requested\n",itno);

	printf("The parameters (a,b) of the prior for the Zhivotovsky estimate are %f and %f\n",A_const,B_const);
	fprintf(out_log,"The parameters (a,b) of the prior for the Zhivotovsky estimate are %f and %f\n",A_const,B_const);

	printf("The maximum frequency of the most common allele allowed is %f\n",Critlim);
	fprintf(out_log,"The maximum frequency of the most common allele allowed is %f\n",Critlim);


	while(1){
		printf("are these parameters correct? (y/n)  ");
		scanf("%s",dic);
		if(dic[0] == 'y')break;
		else if(dic[0] == 'n')exit(1);
		else printf("que ???\n\n\n");
	}

	alls = fopen("out.dat","w");
	printf("\nold out.dat has now been lost\n");

	if(match_samp){
		printf("sample size matching requested. Reading 'ss_file' \n");
		infl = fopen("ss_file","r");
		if(infl == 0){
			printf("no ss_file\n");
			exit(1);
		}

		fscanf(infl,"%d",&nloc);
		fscanf(infl,"%d",&check_subs);
		if(check_subs != Subs){
			printf("The number of samples in 'ss_file' doesn't match those in parameter file");
			exit(1);
		}


		ss = (int **)malloc(nloc*sizeof(int *));
		for(j=0;j<nloc;++j){
			ss[j] = (int *)malloc(Subs*sizeof(int));
			for(i=0;i<Subs;++i)fscanf(infl,"%d",&ss[j][i]);
		}
	}






	/* this is for dominant markers - in this case the sample size refers to individuals */



	rm = 0.5*(1.0/efst - 1);

	val_arr = (int **)malloc(Subs*sizeof(int *));
	freq_arr = (int **)malloc(Subs*sizeof(int *));
	for(j=0;j<Subs;++j){
		val_arr[j] = (int *)malloc(Nmax*sizeof(int));
		freq_arr[j] = (int *)malloc(Nmax*sizeof(int));
	}
	genotype = (int **)malloc(Subs*sizeof(int *));
	for(j=0;j<Subs;++j){
		genotype[j] = (int *)malloc(2*sizeof(int));
		genotype[j][0]=genotype[j][1]=0;
	}
	fsum = wsum = 0.0;
	keepmu = 0;

	/* NBB for some bizarre reason coalescence time is scaled in units of twice
	the number of subpopulations so divide theta by 2*/
	mu = theta; /*???*/
	mono = 0;
	gt_two = 0;
	too_fixed = 0;
	for(kk=0,i=0;1;++kk){

		initot = 0;
		ginitot = 0;

		if(match_samp)lpick = disrand(0,nloc-1);

/* 14.5.04 have changed sim1 so that it can cope with some of the ss[..][..] being 0. */

		for(j=0;j<Subs;++j){
			if(match_samp){
				init[j] = ss[lpick][j]*2; initot += init[j];
				ginit[j] = ss[lpick][j]; ginitot += ginit[j];
			}
			else{
				init[j] = in_smp*2; initot += init[j];
				ginit[j] = in_smp; ginitot += ginit[j];
			}
		}
		for(j=Subs;j<Spno;++j)init[j]=0;
		for(j=Subs;j<Spno;++j)ginit[j]=0;


		sim1(init,initot,rm,mu,freq_arr,val_arr,&noall);
		if(noall == 1){
			++mono;
			continue;
		}
		if(noall > 2){
			++gt_two;
			continue;
		}
		/* assume we have only 2 alleles*/

		irec = disrand(0,1);
		for(j=0;j<Subs;++j){
			genotype[j][0] = genotype[j][1] = 0;
			while(1){
				if(freq_arr[j][0] == 0 && freq_arr[j][1] == 0)break;
				ig = disrand(1,freq_arr[j][0]+freq_arr[j][1]);
				if(ig <= freq_arr[j][0]){
					--freq_arr[j][0];
					gt1 = 0;
				}
				else{
					--freq_arr[j][1];
					gt1 = 1;
				}

				ig = disrand(1,freq_arr[j][0]+freq_arr[j][1]);
				if(ig <= freq_arr[j][0]){
					--freq_arr[j][0];
					gt2 = 0;
				}
				else{
					--freq_arr[j][1];
					gt2 = 1;
				}
				if(gt1 == irec && gt2 == irec)++genotype[j][0];
				else ++genotype[j][1];
			}
		}


		for(j=0,maxf=0;j<2;++j){
			for(jj=0,isum=0;jj<Subs;++jj){
				isum += genotype[jj][j];
			}
			if(isum > maxf)maxf = isum;
		}
		if((double)maxf/ginitot > Critlim){
			++too_fixed;
			continue;
		}

		my_thetacal(genotype,2,ginit,Subs,&h0,&h1,&h2,&fst);
		fprintf(alls,"%f %f \n",h2,fst);
		fsum += fst*h1;
		wsum += h1;
		++i;
		if(i%10 == 0)fflush(alls);
		if(i==itno)break;
	}
	printf("average Fst is %f\n",fsum/wsum);
	printf("number rejected is %d\n",kk+1-itno);
	printf("         number mono is %d\n",mono);
	printf("         number with > 2 alleles is %d\n",gt_two);
	printf("         number with freq. common allele > %f is %d\n",Critlim,too_fixed);

	fprintf(out_log,"average Fst is %f\n",fsum/wsum);
	fprintf(out_log,"number rejected is %d\n",kk+1-itno);
	fprintf(out_log,"         number mono is %d\n",mono);
	fprintf(out_log,"         number with > 2 alleles is %d\n",gt_two);
	fprintf(out_log,"         number with freq. common allele > %f is %d\n",Critlim,too_fixed);

 	closegfsr();
	printf("type in any character and return to close window");
	scanf("%s",str1);
}

sim1(init,initot,rm,mu,freq_arr,val_arr,noall)
int init[],initot,*freq_arr[],*val_arr[],*noall;
double mu,rm;
{

	int j,nmax,ic,k;
	float rr;

	NEXTMUT = 0;
	*noall = 0;
	for(j=0;j<Spno;++j){
		Sd[j] = Spno/0.5; /* coalescence rate scaled in units of twice the metapopulation size */
		Migrate[j] = Sd[j]*rm;
	}
/*	Occ = Subs;  14.5.04 change this so that allows possibility that some of init[...] are 0 */
	for(j=0;j<Spno;++j)Occlist[j] = 0;
	for(j=0,k=0;j<Subs;++j){
		if(init[j] > 0){
			Occlist[k] = j;
			++k;
		}
	}
	Occ = k;
	/* 15.8.07 bug fix for when there are missing data */
	for(j=0;j<Occ;++j){
		Ni[j] = init[Occlist[j]];
	}
	for(j=Occ;j<Spno;++j)Ni[j] = 0;
	Ntot = initot;
	nmax = 10*Ntot;
	Nlist = (struct node **)malloc(nmax*sizeof(struct node *));
	ic = 0;
	for(k=0;k<Occ;++k){
		Lmax[k] = 2*Ni[k];
		if(Lmax[k] < 10)Lmax[k] = 10;
	}
	for(k=Occ;k<Spno;++k){
		Lmax[k] = 10;
	}
	for(k=0;k<Occ;++k){
		List[k] = (struct node **)malloc(Lmax[k]*sizeof(struct node *));
		for(j=0;j<Ni[k];++j){
			List[k][j] = (struct node *)malloc(sizeof(struct node));
			List[k][j]->d[0]=List[k][j]->d[1]=NULL;
			List[k][j]->a[0]=List[k][j]->a[1] = NULL;
			List[k][j]->time = 0.0;
			List[k][j]->dna = 0;
			List[k][j]->I = 0;
			List[k][j]->sp = Occlist[k];
			List[k][j]->osp = Occlist[k];
			Nlist[ic] = List[k][j];
			++ic;
		}

	}
	for(k=Occ;k<Spno;++k){
		List[k] = (struct node **)malloc(Lmax[k]*sizeof(struct node *));
	}
	N_n = Ntot;
	/* 15.8.07 end of bug fix for when there are missing data */
	Tt = 0.0;
	while(1){
		if(Occ > Spno){
			printf("error Occ > Spno\n");
			exit(1);
		}
		for(k=0;k<Occ;++k){
			if(Ni[k] > Lmax[k]){
				printf("error in Ni/Lmax\n");
				exit(1);
			}
			if(Ni[k] >= Lmax[k]-5){
				Lmax[k] = 2*Ni[k];/*15.8.07 fix possible bug if Ni[k] changes by a large amount */
				if(Ni[k] > Lmax[k]){
					printf("error - Lmax");
					exit(1);
				}
				for(j=0;j<Ni[k];++j){
					if(List[k][j]->sp != Occlist[k]){
						printf("error in sp \n");
						exit(1);
					}
				}
				List[k] = (struct node **)realloc(
					List[k],Lmax[k]*sizeof(struct node *));
				for(j=0;j<Ni[k];++j){
					if(List[k][j]->sp != Occlist[k]){
						printf("error in sp \n");
						exit(1);
					}
				}
			}
		}
		if(N_n >= nmax-1){
			nmax = 2*N_n; /* 15.8.07 fix possible bug if Ni[k] changes by a large amount */
			Nlist = (struct node **)realloc(
				Nlist,nmax*sizeof(struct node *));
		}

		dfill();

		Tt +=  expdev()/Dtop;
		

loopback: rr = gfsr4();
		for(k=0;k<Occ;++k){
			for(j=0;j<2;++j){
				if(rr < Den[k][j]/Dtop)goto loopout;
			}
		}
		goto loopback;
loopout:  if(j==0){
			cnode(k);
			if(Ntot == 1)break;
		}
		else mnode(k);
		


	}

	Nlist[N_n-1]->dna = 0;
	*noall = 0;
	for(j=N_n-1;j>=0;--j){
		treefill(Nlist[j],noall,freq_arr,val_arr,mu);
	}

	for(j=0;j<N_n;++j)killtree(Nlist[j]);
	for(j=0;j<Spno;++j)free(List[j]);
	free(Nlist);

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
	*fst = 1.0 - (*het0)/(*het1);

	free(psum);

}

/*
checkit()
{
	for(j=0;j<N_n;++j){
*/

dfill()
{
	int k;

	Den[0][0] = Sd[Occlist[0]]*Ni[0]*(Ni[0]-1.0)*0.5;
	Den[0][1] = Den[0][0] + Migrate[Occlist[0]]*Ni[0];
	for(k=1;k<Occ;++k){
		Den[k][0] = Den[k-1][1] + Sd[Occlist[k]]*Ni[k]*(Ni[k]-1.0)*0.5;
		Den[k][1] = Den[k][0] + Migrate[Occlist[k]]*Ni[k];
	}
	Dtop = Den[Occ-1][1];
}


cnode(sp)
int sp;
{
	int ind1,ind2,temp,rfs;
	float gfsr4();
	struct node *p1;
	float expdev();
	while(1){
		ind1 = disrand(0,Ni[sp]-1);
		ind2 = disrand(0,Ni[sp]-1);
		if(ind2 != ind1)break;
	}
	if(ind1 > ind2){
		temp = ind1;
		ind1 = ind2;
		ind2 = temp;
	}
	p1 = (F) malloc((unsigned)sizeof(struct node));
	p1->time = Tt;
	p1->d[0] = List[sp][ind1];
	p1->d[1] = List[sp][ind2];
	p1 -> a[0] = p1 -> a[1] = NULL;
	p1->dna = 0;
	p1->I = 0;
	p1->sp = Occlist[sp];
	List[sp][ind1]->a[0] = p1;
	List[sp][ind2]->a[0] = p1;
	List[sp][ind1]->a[1] = List[sp][ind2]->a[1] = NULL;
	List[sp][ind1] = p1;
	List[sp][ind2] = List[sp][Ni[sp]-1];
	--Ni[sp];
	--Ntot;
	Nlist[N_n] = p1;
	++N_n;
	return;
}


mnode(sp)
int sp;
{
	int ind,disrand(),ifs,j,nifs,rn,i,k,click,ii,jj,it;
	int c1,c2,nc1,nc2;
	float gfsr4(),rr;
	struct node *tp,*p;

	ind = disrand(0,Ni[sp]-1);

/* four corners toroidal stepping stone

	c1 = Occlist[sp]/SIDE;
	c2 = Occlist[sp] - c1*SIDE;
	nc1 = disrand(0,1)*2-1;
	nc1 = tau(SIDE,nc1+c1);
	nc2 = disrand(0,1)*2-1;
	nc2 = tau(SIDE,nc2+c2);
	j = nc1*SIDE +nc2; */

/* circular stepping stone

	nc1 = disrand(0,1)*2-1;
	j = tau(SPNO,Occlist[sp]+nc1); */

/* "normal" toroidal stepping stone

	c1 = Occlist[sp]/SIDE;
	c2 = Occlist[sp] - c1*SIDE;
	if(disrand(0,1)){
		nc1 = disrand(0,1)*2-1;
		nc2 = 0;
	}
	else{
		nc1 = 0;
		nc2 = disrand(0,1)*2-1;
	}
	nc1 = tau(SIDE,nc1+c1);
	nc2 = tau(SIDE,nc2+c2);
	j = nc1*SIDE +nc2; */

/* 8 neighbour toroidal stepping stone */

/*	c1 = Occlist[sp]/SIDE;
	c2 = Occlist[sp] - c1*SIDE;
	while(1){
		nc1 = disrand(0,2)-1;
		nc2 = disrand(0,2)-1;
		if(!(nc1 == 0 && nc2 == 0))break;
	}
	nc1 = tau(SIDE,nc1+c1);
	nc2 = tau(SIDE,nc2+c2);
	j = nc1*SIDE +nc2;  */

/* finite island model */

	while(1){
		j = disrand(0,Spno-1);
		if(j != Occlist[sp])break;
	}


/* start of stuff */

	List[sp][ind]->sp = j;
	tp = List[sp][ind];
	List[sp][ind] = List[sp][Ni[sp]-1];
	for(i=0;i<Occ;++i){
		if(Occlist[i] == j)break;
	}
	if(i == Occ){
		if(Ni[sp] == 1){
			i = sp;
			Occlist[i] = j;
			Ni[i] = 0;
			click = 0;
		}
		else{
			++Occ;
			Occlist[i] = j;
			Ni[i] = 0;
			--Ni[sp];
			click = 1;

		}
	}
	else{
		if(Ni[sp] == 1){
			for(k=sp;k<Occ-1;++k){
				Occlist[k] = Occlist[k+1];
				Ni[k] = Ni[k+1];
				if(Ni[k] >= Lmax[k]){
					Lmax[k] = 2*Ni[k]; /* 15.8.07 change to fix bug if Ni[k] changes by a large amount */
					List[k] = (struct node **)realloc(
						List[k],Lmax[k]*sizeof(struct node *));
				}
				for(jj=0;jj<Ni[k];++jj){
					List[k][jj] = List[k+1][jj];
				}
			}
			--Occ;
			if(i>sp)--i;
			click = 2;
		}
		else{
			--Ni[sp];
			click = 3;
		}
	}
	List[i][Ni[i]] = tp;
	++Ni[i];
	return;
}





void treefill(p,noall,freq_arr,val_arr,mu)
struct node *p;
int *noall,*freq_arr[],*val_arr[];
double mu;
{
	struct node *bro;
	int j,mutno,pos,sp,i;
	double time;

	if(p->a[0] == NULL && p->a[1] == NULL){
		return;
	}
	else if(!(p->a[0] == NULL && p->a[1] == NULL)){
		if(p->a[0]->d[0] == p)bro = p->a[0]->d[1];
		else bro = p->a[0]->d[0];
		p->dna = p->a[0]->dna;
	}
	time = p->a[0]->time - p->time;
	mutno = poidev(time*mu);
	for(j=0;j<mutno;++j){
		p->dna = addmut(p->dna);
	}
	if(p->d[0] == NULL && p->d[1] == NULL){
		sp = p->osp;
		for(j=0;j<*noall;++j){
			if(val_arr[sp][j] == p->dna)break;
		}
		if(j<*noall)++freq_arr[sp][j];
		else{
			for(i=0;i<Subs;++i){
				val_arr[i][j] = p->dna;
				freq_arr[i][j] = 0;
			}
			freq_arr[sp][j] = 1;
			++(*noall);
			if(*noall == Nmax){
				Nmax += Increment1;
				for(i=0;i<Subs;++i){
					val_arr[i] = (int *)realloc(val_arr[i],Nmax*sizeof(int));
					freq_arr[i] = (int *)realloc(freq_arr[i],Nmax*sizeof(int));
				}
			}
		}
		return;
	}
	return;
}



void killtree(p)
struct node *p;
{
	free(p);
	return;
}



int addmut(p)
int 	p;
{
	int ic;
	if(Ms){
		ic = disrand(0,1)*2-1;
		return p + ic;
	}
	else return ++NEXTMUT;
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

int poidev(xm)
float xm;
{
	static float sq,alxm,g,oldm=(-1.0);
	float em,t,y;
	float gfsr4(),gammln();

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			em += 1.0;
			t *= gfsr4();
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*gfsr4());
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (gfsr4() > t);
	}
	return (int)(em+0.5);
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



mom(x,n,x1,x2,x3,x4,min,max)
int n;
float x[],*x1,*x2,*x3,*x4,*min,*max;
{
	int i;
	double s1,s2,s3,s4,an,an1,dx,dx2,xi,var,pow(),sqrt();

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
	*x4 = var>0.0 ? 1.0/(n-1)*s4/pow(var,2.0)-3.0 : 0.0;
	return;
}

int tau(lim,val)
int lim,val;
{
      if (val>=0)
            return(val%lim);
      else return(tau(lim,val%lim+lim));
}

void sort(n,ra)
int n;
float ra[];

{
	int i,j,l,ir;
	float rra;

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


