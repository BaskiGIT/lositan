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
int NEXTMUT,Ms,Nmax;
int Smp,Subs,Lmax[NOC];

void 	opengfsr(),closegfsr(),cpress(),treefill(),killtree();
float gfsr4(),expdev();
int disrand();

main()
{
	double mu,next_mu,mrate,dd;
	float	rr,efst,fsum,wsum;
	float	tcum;
	int	j,j2,i,itno,it,seq,ic,k,ik,ll,initot,init[NOC],keepmu;
	double rm;
	int i1,j1,kk,**val_arr,**freq_arr,noall;
	float	vhet,h0,h1,fst;
	int	totall;
	char dic[100],str1[10];
	FILE *inp,*alls;
	
	
	opengfsr();
	inp = fopen("fdist_params2.dat","r");
	fscanf(inp,"%d",&Spno);
	if(Spno > NOC){
		printf("error in parameter file - ");
		printf("number of subpopulation greater than %d\n",NOC);
		exit(1);
	}
	fscanf(inp,"%d",&Subs);
	fscanf(inp,"%f",&efst);
	fscanf(inp,"%d",&Smp);
	fscanf(inp,"%d",&Ms);
	fscanf(inp,"%d",&itno);
	fclose(inp);
	printf("number of demes is %d\n",Spno);
	printf("number of samples  is %d\n",Subs);
	printf("expected (infinite allele, infinite island) Fst is %f\n",efst);
	printf("median sample size is %d - suggest give 50 if >50\n",Smp);
	if(Ms){
		printf("stepwise mutation model assumed\n");
	}
	else printf("infinite alleles mutation model assumed\n");
	printf("%d realizations (loci) requested\n",itno);
	while(1){
		printf("are these parameters correct? (y/n)  ");
		scanf("%s",dic);
		if(dic[0] == 'y')break;
		else if(dic[0] == 'n')exit(1);
		else printf("que ???\n\n\n");
	}
	
	alls = fopen("out.dat","w");
	printf("\nold out.dat has now been lost\n");
	
	rm = 0.5*(1.0/efst - 1);  

	for(j=0;j<Subs;++j){
		init[j] = Smp;
	}
	for(j=Subs;j<Spno;++j)init[j]=0;
	initot = Smp*Subs;
	val_arr = (int **)malloc(Subs*sizeof(int *));
	freq_arr = (int **)malloc(Subs*sizeof(int *));
	for(j=0;j<Subs;++j){
		val_arr[j] = (int *)malloc(Smp*sizeof(int));
		freq_arr[j] = (int *)malloc(Smp*sizeof(int));
	}
	Nmax = Smp;
	fsum = wsum = 0.0;
	keepmu = 0;
	for(kk=0,i=0;1;++kk){
		if(!keepmu){
			rr = gfsr4();
			if(!Ms)mu = rr/(1-rr);
			else{
				dd = 1/(1-rr);
				dd *= dd;
			 	mu = (dd - 1)*0.5;
			}
			if(mu > 100000)mu = 100000;
		}
		sim1(init,initot,rm,mu,freq_arr,val_arr,&noall);
		if(noall > 1){
			my_thetacal(freq_arr,noall,init,Subs,&h0,&h1,&fst);
			fprintf(alls,"%f %f\n",h1,fst);
			fsum += fst*h1;
			wsum += h1;
			++i;
			if(i%10 == 0)fflush(alls);
			if(i==itno)break;
			keepmu = 0;
		}
		else{++keepmu; if(keepmu > 1000)keepmu=0;}
	}
	printf("average Fst is %f\n",fsum/wsum);
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
		Sd[j] = Spno/0.5;
		Migrate[j] = Sd[j]*rm;      
	}
	Occ = Subs;
	for(j=0;j<Spno;++j)Occlist[j] = 0;
	for(j=0;j<Subs;++j){
		Occlist[j] = j;
	}
	for(j=0;j<Spno;++j){
		Ni[j] = init[j];
	}
	Ntot = initot;
	nmax = 10*Ntot;
	Nlist = (struct node **)malloc(nmax*sizeof(struct node *));
	ic = 0;
	for(k=0;k<Spno;++k){
		Lmax[k] = 2*init[k];
		if(Lmax[k] < 10)Lmax[k] = 10;
	}
	for(k=0;k<Spno;++k){
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
	N_n = Ntot;
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
				Lmax[k] = 2*Lmax[k];
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
			nmax = 2*nmax;
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
					Lmax[k] *= 2;
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
				Nmax += Smp;
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
	
my_thetacal(gen,noall,sample_size,no_of_samples,het0,het1,fst)
	
int *gen[],noall,sample_size[],no_of_samples;
float *het0,*het1,*fst;
{
	int i,j,k;
	double x2,x0,yy,y1,q2,q3;
	
	
	x0 = 0.0;
	for(j=0;j<no_of_samples;++j) {
		for(i=0,x2=0.0;i<noall;++i) {
			x2 += gen[j][i]*gen[j][i];
		}
		x0 += (x2-sample_size[j])/(sample_size[j]*(sample_size[j]-1));
	}
	
	yy = 0.0;
	for(j=0;j<no_of_samples;++j){
		for(k=j+1;k<no_of_samples;++k){
			for(i=0,y1=0.0;i<noall;++i){
				y1 += gen[j][i]*gen[k][i];
			}
			yy += y1/(sample_size[j]*sample_size[k]);
		}
	}
	
	
	
	q2 = x0/no_of_samples;
	q3 = 2*yy/(no_of_samples*(no_of_samples-1));
	
			
	*het0 = 1.0 - q2;
	*het1 = 1.0 - q3;
	*fst = 1.0 - (*het0)/(*het1);
	
	
}

