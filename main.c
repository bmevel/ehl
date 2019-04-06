

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

/* FMG Solver of the EHL circular contact */
#define pi 3.1415926535897931

typedef struct
{
double hx;              /*mesh size x*/
double hy;              /*mesh size y*/
int ii;                 /*number of nodes x*/
int jj;                 /*number of nodes y*/
double **p, **f;        /*pressure and right hand side values pressure*/
double **hfi , **hrhs ; /*film thickness and right hand side film thickness*/
double **w ;           /*elastic deformation integrals*/
double **K, **K1 ;       /*kernels*/
double **pconv ;         /*converged pressure for convergence check in FMG */
double **pold ;          /*'old' pressure for use in FAS coarsening*/
double **pjac ;        /*'old' pressure for use in jacobi relaxation*/
double **A, *X, *Y ;     /*for line relaxation*/
double rg ;             /*right hand side of force balance equation*/
double Hm, Hcp, Hc ;     /*minimum and central film thickness for output*/
} Level;



typedef struct
{
int nx0,ny0; /*number of nodes coarsest grid*/
int m1,m2,od; /*number of correction points*/
int maxlev, deep;/*number of grid-levels, grids deep*/
double xa,xb; /*begin,end of computationaI domain*/
double ya,yb; /*begin,end of computational domain*/
double h0,wu; /*global constant and work unit*/
Level *Lk ; /*array of grid levels*/
} Stack;

/********** GLOBAL VARIABLES *****************************/

double	MMoes, LMoes, H_0, rlambda, alphabar;
double	p0r, alpha, eta0, zr;
int	maxl, deepl, starl, order, typecy, ncy, outlev, currfilev;
unsigned long *cputimes,t0;
double	hfact,xi_l,urja,urgs;
int	typecy;


/********ROUTINES FOR DATASTRUCTURE*************************/

double **matrix(int nx, int ny, int shift)
{
int i;

double **m=(double **)calloc(nx+1 ,sizeof(double*));
for (i=0; i<=nx; i++) m[i]=(double *)calloc(ny+1, sizeof(double))+shift ;
return m+shift;
}

void initialize(Stack *U, int nx0, int ny0, int maxl, int deepl, int ord,
               double xa, double xb, double ya, double yb, double h0)
{
/* initialize values in datastructure */


double hx,hy;
Level *L;
int l,ii,jj;

U->xa=xa;
U->xb=xb;
U->ya=ya;
U->yb=yb;
U->maxlev=maxl;
U->deep=deepl;
U->wu=0.0;
U->od = ord-2;
U->h0=h0;

cputimes = (unsigned long *) calloc(maxl+1 , sizeof(unsigned long));

U->Lk=(Level *)calloc(maxl+1,sizeof(Level));

hx=(xb-xa)/nx0;
hy=(yb-ya)/ny0;
ii=nx0;
jj=ny0;
for (l=1;l<=maxl;l++)

{
L=U->Lk+l;
L->hx=hx;
L->hy=hy;
L->ii=ii;
L->jj=jj;
L->p	=matrix(ii+2*U->od, jj+2*U->od,U->od);
L->w	=matrix(ii+2*U->od, jj+2*U->od,U->od);
L->f	=matrix(ii+2*U->od, jj+2*U->od,U->od);
L->pold =matrix(ii+2*U->od, jj+2*U->od,U->od);
L->pjac =matrix(ii+2*U->od, jj+2*U->od,U->od);
L->pconv=matrix(ii+2*U->od, jj+2*U->od,U->od);
L->hfi	=matrix(ii+2*U->od, jj+2*U->od,U->od);
L->hrhs =matrix(ii+2*U->od, jj+2*U->od,U->od);


L->K	=matrix(ii+4*U->od,jj+4*U->od,0);
L->K1   =matrix(ii+4*U->od,jj+4*U->od,0);

L->A=matrix(ii-1,11,0);
L->X=(double *)calloc(ii-1,sizeof (double));
L->Y=(double *)calloc(ii-1,sizeof (double));

printf("\n level: %2d ii=%4d, jj=%4d hx=%f hy=%f",l,ii,jj,hx,hy);

if ((2*(l/2)-l)!=0) {hx*=0.5; ii*=2;} else {hy*=0.5; jj*=2;}

if ((L->p==NULL) || (L->w==NULL) || (L->f==NULL) || (L->K==NULL) || (L->pold==NULL)
                 || (L->pconv==NULL) || (L->hfi==NULL) || (L->hrhs==NULL) || (L->K==NULL)
                 || (L->K1==NULL) || (L->A==NULL) || (L->X==NULL) || (L->Y==NULL))
	{
	printf (" \nprobleme allocating memory\n");
	exit(-1);
	}
}

U->m1=(int)(3+log(1.0*(U->Lk+maxl)->ii)); U->m2=2;
printf("\ncorrection patch in multi-integration: m1=%2d, m2=%2d ",U->m1 ,U->m2);
}


void finalize(Stack *U, int maxlevel)
{
/* free memory at end of program */

Level *L;
int i;
for (i=1 ; i<=U->maxlev; i++)
{
L=U->Lk+i;
free (L->f-U->od);
free (L->w-U->od);
free (L->p-U->od);
free(L->K);
free(L->K1);
free (L->pold-U->od);
free(L->pconv-U->od);
free(L->hfi-U->od);
free(L->hrhs-U->od);
free(L->A);
free(L->X);
free(L->Y);
}

free(cputimes);
free(U->Lk);
}


/****************SPECIAL FUNCTIONS********************/

double reta(double p)
{
/* Barus */

return(exp(-alphabar*p));


/* Roelands */

/*return(exp(-alpha*p0r/zr*(-1.0+pow(1.0+(p/p0r)*(alphabar/alpha),zr)))); */
}


double rho(double p)
{
/* Incompressible */

return(1.0);


/* Compressible */
/* return((5.9e8+1.34*(alphabar/alpba)*p)/(5.9e8+(alphabar/alpha)*p)); */
}

double Lu(double **H, double **P, double rhx, double rhx2, double rhy2, int i,	int j)
/* computes operator in point i,j */
{
double H3, Hx, Qx, Qy,xi_n,xi_s,xi_e,xi_w,r0, r1, r2;

r0=rho(P[i][j]);
H3=r0*H[i][j]*H[i][j]*H[i][j]*reta(P[i][j]);
xi_n=0.5*(rho(P[i  ][j+1])*H[i	][j+1]*H[i  ][j+1]*H[i	][j+1]*reta(P[i  ][j+1])+ H3)*rlambda*rhy2;
xi_s=0.5*(rho(P[i  ][j-1])*H[i	][j-1]*H[i  ][j-1]*H[i	][j-1]*reta(P[i  ][j-1])+ H3)*rlambda*rhy2;
xi_e=0.5*(rho(P[i+1][j  ])*H[i+1][j  ]*H[i+1][j  ]*H[i+1][j  ]*reta(P[i+1][j  ])+ H3)*rlambda*rhx2;
xi_w=0.5*(rho(P[i-1][j	])*H[i-1][j  ]*H[i-1][j  ]*H[i-1][j  ]*reta(P[i-1][j  ])+ H3)*rlambda*rhx2;

Qx=(xi_e*P[i+1][j  ]-(xi_e+xi_w)*P[i][j] + xi_w*P[i-1][j  ]);
Qy=(xi_n*P[i  ][j+1]-(xi_n+xi_s)*P[i][j] + xi_s*P[i  ][j-1]);

r1=rho(P[i-1][j]);
if (i==1) Hx=rhx*(r0*H[i][j]-r1*H[i-1][j]) ;
else
 {
 r2=rho(P[i-2][j]);
 Hx=rhx*(1.5*r0*H[i][j]-2*r1*H[i-1][j]+0.5*r2*H[i-2][j]);
 }
return(Qx+Qy-Hx);
}

/******* SINGLE GRID ROUTINES  *******/

void init_f(Stack *U, int l)
{
	int i, j;
	Level*L;
	double x, y;

	L = U->Lk + l;
	for (i = -1; i <= L->ii; i++)
	{
		x = U->xa+ i*L->hx  ;
		for (j=0;j<=L->jj;j++)
		{
			y= U->ya+ j*L->hy ;
			L->hrhs[i][j] = 0.5*x*x+0.5*y*y;
		}
	}
L->rg = -2.0*pi / 3.0;


}


void init_p(Stack *U, int l)
{
int i,j;
Level *L;
double x,y;
void calchi();

L=U->Lk+l;

for (i=0; i<=L->ii; i++)
 {
 x= U->xa+i*L->hx ;
 for (j=0; j<=L->jj;j++)
 {
  y= U->ya + j*L->hy ;
  L->p[i][j] = 0.0;
  if (x*x+y*y<1.0) L->p[i][j]=sqrt(1.0-x*x-y*y) ; else L->p[i][j]=0.0 ;
  }
 }
calchi(U,l);
}


double resnorm(Stack *U, int l)
{
/* computes the absolute sors of the residual */
/* on level l in points where the pressure>O  */


int i,j;
Level *L;

double herr,rhx2,rhy2,rhx;

L = U->Lk+l;
rhx= 1.0/L->hx;
rhx2= rhx*rhx;
rhy2=1.0/(L->hy*L->hy);
herr = 0.0;

for  (i=1; i<=L->ii-1; i++)
  for (j=1; j<=L->jj-1; j++)
  {
  if (L->p[i][j]>0.0)
  herr+=fabs(L->f[i][j]-Lu(L->hfi,L->p,rhx, rhx2,  rhy2, i, j)) ;
  }
 return (herr/(L->ii-1)/(L->jj-1));

}


void relax(Stack *U,  int l)
/* relaxes the equation on level l */
{
int	i,j;
Level *L;
double **P, **Pj;
double g;

void calchi();
void solve_line();

L=U->Lk+l;
P = L->p;
Pj=L->pjac;

calchi(U,l);

for (i=0; i<=L->ii; i++)
 for (j=0; j<=L->jj; j++)
   Pj[i][j]=P[i][j];


for (j=1; j<L->jj; j++)
solve_line(U,l,j);

for (i=0; i<=L->ii; i++)
 for (j=0; j<=L->jj; j++)
   P[i][j]=Pj[i][j];

g=0.0;
for (i=0; i<=L->ii; i++)
 for (j=0; j<=L->jj; j++)
   g += L->p[i][j];

g = g*L->hx*L->hy+L->rg;


U->wu+=pow(0.25,1.0*((U->maxlev+1)/2-(l+1)/2));
if (l==outlev)
    printf("\n k=%2d, resn=%10.3e, g=%10.3e, h0=%12.8e, wu= %7.3f", (l+1)/2, resnorm(U,l),g,U->h0, U->wu);

}


void relaxh0(Stack *U, int l, int lf)
/* relaxes the force-balance equation */
/* hfact <0.05 gives stable convergence for W cycles */
{
int i,j;
Level *L;
double g;

L=U->Lk+l;
g=0.0;
for (i=0; i<=L->ii; i++)
  for (j=0; j<=L->jj; j++)
  g += L->p[i][j];
g = g*L->hx*L->hy + L->rg;
U->h0 += hfact*g;
}


/***************************INTER GRID ROUTINES**********/


void coarsen_p(Stack*U, int l)
{
/* coarsen the solution from level 1 to level l-2 */
/* and store: coarse grid solution in pold array */

int i, j,iff,jff;
Level *Lc, *L;
double **p ,**hfi;
void calchi();
void init_log();

L  =U->Lk+l;
Lc =U->Lk+l-2;
p = L->p;
hfi=L->hfi;

calchi(U,l);
init_log(U,l-2);

for (i=0; i<=Lc->ii;i++)
	{
	iff=2*i;
	for(j=0; j<=Lc->jj;j++)
	{
	jff=2*j;
	if ((i==0)||(j==0)||(i==Lc->ii)||(j==Lc->jj))
	Lc->p[i][j]=0.0;
	else
	{
	if ((p[iff  ][jff  ]==0) ||
	    (p[iff+1][jff+1]==0) || (p[iff+1][jff-1]==0)||
	    (p[iff-1][jff+1]==0) || (p[iff-1][jff-1]==0)||
	    (p[iff+1][jff  ]==0) || (p[iff-1][jff  ]==0)||
        (p[iff  ][jff+1]==0) || (p[iff  ][jff-1]==0))
	Lc->p[i][j]= p[iff ][jff ];
	else
	Lc->p[i][j]=
	(4.0*p[iff ][jff ]+
	 2.0*(p[iff+1][jff  ]+p[iff-1][jff  ]+p[iff  ][jff+1]+p[iff  ][jff-1])+
	      p[iff+1][jff+1]+p[iff+1][jff-1]+p[iff-1][jff+1]+p[iff-1][jff-1] )/16.0 ;
	      }
	  }
	}

	for (i=0;i<=Lc->ii; i++)
	{
	 iff=2*i;
	 for(j=0;j<=Lc->jj; j++)
	 {
	 jff=2*j ;
 	 Lc->hfi[i][j]=hfi[iff][jff];
	 }
	}

for (i=0;i<=Lc->ii; i++)
 for (j=0;j<=Lc->jj;j++)
	 Lc->pold[i][j]=Lc->p[i][j];
}


void coarsen_f(Stack *U, int l)
{
/**** compute coarse grid right hand side on level l-2: */
/* in coarsening step from level 1	*/

int i,j,iff,jff;
Level *Lc,*L;
double	**f, **p, gf, gc, rhx, rhy, rhx2, rhy2, rhxc, rhyc, rhxc2, rhyc2;
double r0, rn, re, rs, rw, rne, rnw, rse, rsw;
void calchi();


L =U->Lk+l;
Lc=U->Lk+l-2;

p=L->p ;
f=L->f ;

rhx=1.0/L->hx;
rhy=1.0/L->hy;
rhx2=rhx*rhx;
rhy2=rhy*rhy;

rhxc=1.0/Lc->hx;
rhyc=1.0/Lc->hy;
rhxc2=rhxc*rhxc;
rhyc2=rhyc*rhyc;

calchi(U, l-2);

/* right hand side film thickeness equation */

for (i=0; i<=Lc->ii;i++)
 for (j=0;j<=Lc->jj;j++)
  Lc->hrhs[i][j]=-U->h0-2./pi/pi*Lc->w[i][j]+Lc->hfi[i][j];

/* right hand side Reynolds equation */

for (i=1;i<=Lc->ii-1;i++)
{
 iff=2*i;
 for (j=1;j<=Lc->jj-1;j++)
 {
 jff=2*j;
 r0=f[iff][jff]-Lu(L->hfi,L->p, rhx, rhx2, rhy2, iff, jff);
 if (p[iff][jff+1 ]>0)
 rn= f[iff ][jff+1]-Lu(L->hfi,L->p,rhx,rhx2,rhy2,iff,jff+1);
 else  rn= 0.0;
 if (p[iff+1][jff]>0)
 re= f[iff+1][jff ]-Lu(L->hfi,L->p, rhx, rhx2, rhy2, iff+1, jff);
 else re= 0.0 ;
 if (p[iff][jff-1]>0)
 rs=f[iff ][jff-1]-Lu(L->hfi,L->p,rhx,rhx2,rhy2,iff,jff-1);
 else rs=0.0;
 if (p[iff-1][jff]>0)
 rw= f[iff-1][jff]-Lu(L->hfi,L->p,rhx,rhx2,rhy2,iff-1,jff );
 else rw=0.0;
 if (p[iff+1][jff+1]>0)
 rne=f[iff+1][jff+1]-Lu(L->hfi,L->p,rhx,rhx2,rhy2,iff+1,jff+1);
 else rne=0.0;
 if (p[iff-1][jff+1]>0)
  rnw=f[iff-1][jff+1]-Lu(L->hfi,L->p,rhx,rhx2,rhy2,iff-1,jff+1);
 else rnw=0.0;
 if (p[iff+1][jff-1]>0)
  rse=f[iff+1][jff-1]-Lu(L->hfi,L->p,rhx,rhx2,rhy2,iff+1,jff-1);
 else rse=0.0;
 if (p[iff-1][jff-1]>0)
  rsw=f[iff-1][jff-1]-Lu(L->hfi,L->p,rhx,rhx2,rhy2,iff-1, jff-1);
 else rsw=0.0;


 Lc->f[i][j]=(4.0*r0+2.0*(rn+rs+re+rw)+(rne+rse+rsw+rnw))/16.0;
 Lc->f[i][j]+=Lu(Lc->hfi,Lc->p,rhxc,rhxc2,rhyc2,i,j);
 }
}

/* right hand side force balance equation */
gf=0.0;
for (i=0; i<=L->ii; i++)
 for (j=0; j<=L->jj; j++)
  gf += L->p[i][j];


gc = 0.0;
for (i=0; i<=Lc->ii; i++)
  for (j=0; j<=Lc->jj; j++)
  gc += Lc->p[i][j];

Lc->rg=(L->rg+gf*L->hx*L->hy)-gc*Lc->hx*Lc->hy;
}

void refine(Stack *U, int k)
{
/* Interpolation and addition of coarse grid correction f rom grid k72 */
/* to grid k	*/

int i,j,ic,jc,iic,jjc,ii,jj;
Level *Lc,*L;
double **pc,**pco,**p;

void calchi();

L=U->Lk+k;
ii = L->ii; jj=L->jj;
p  = L->p;

Lc= U->Lk+k-2;
iic=Lc->ii; jjc=Lc->jj;
pc = Lc->p;
pco=Lc->pold ;

for (ic=1;ic<=iic; ic++)
 for (jc=1;jc<=jjc; jc++)
 {
 if (p[2*ic ][2*jc ]>0)
		p[2*ic	][2*jc	]+=(pc[ic][jc]-pco[ic][jc]);

 if ((jc<jjc)&&(p[2*ic-1][2*jc ]>0))
		 p[2*ic-1][2*jc	]+=(pc[ic  ][jc]-pco[ic  ][jc ]+
		                    pc[ic-1][jc]-pco[ic-1][jc ])*0.5;

 if ((ic<iic)&&(p[2*ic	][2*jc-1]>0))
		  p[2*ic  ][2*jc-1]+=(pc[ic][jc	 ]-pco[ic][jc  ]+
		                      pc[ic][jc-1]-pco[ic][jc-1])*0.5;

 if  (p[2*ic-1][2*jc-1]> 0)
		  p[2*ic-1][2*jc-1]+=(pc[ic  ][jc  ]-pco[ic  ][jc ]+
		                      pc[ic  ][jc-1]-pco[ic  ][jc-1]+
		                      pc[ic-1][jc  ]-pco[ic-1][jc  ]+
		                      pc[ic-1][jc-1]-pco[ic-1][jc-1])*0.25;
}


for (i=0;i<=ii;i++)
  for (j=0;j<=jj;j++)
      if (p[i][j]<0) p[i][j]=0.;

calchi(U,k);
}


void fmg_interpolate(Stack *U, int k)
{
		/* interpolation of coarse grid k-2 solution to fine grid k */
		/* to serve as first approximation. bi-cubic interpolation	*/

int ic,jc,iic,jjc,i,j,ii,jj;
Level *Lc,*L;
double **pc,**p, **pconv, x,y;


/* set time */

cputimes[k-2]=clock();

L = U->Lk+k;
ii = L->ii; jj = L->jj;
p=L->p ;

Lc =U->Lk+k-2;
iic	= Lc->ii;
jjc = Lc->jj;
pc	= Lc->p	;
pconv=Lc->pconv;

/* store coarse grid solution for inter use in convergence check */
/* and store minimum and central film thickens for output */

jc = jjc/2 ;
ic =1 ;
while ((Lc->p[ic][jc]>Lc->p[ic-1][jc])&&(ic<iic)) ic++ ;
		Lc->Hcp=Lc->hfi[ic][jc] ;

Lc->Hm=1e5; /* arbitrary large value */

for (ic=1; ic<=iic-1; ic++)
{
 x= U->xa +ic*Lc->hx ;
 for (jc=1;jc<=jjc-1;jc++)
 {
  y= U->ya + jc*Lc->hy  ;
  if (Lc->hfi[ic][jc]<Lc->Hm) Lc->Hm = Lc->hfi[ic][jc];
  if ((x==0)&&(y==0)) Lc->Hc = Lc->hfi[ic][jc];
    pconv[ic][jc] =pc[ic][jc] ;
  }
}

/* interpolation*/
/* first inject to points coinciding with coarse grid points */


for (ic=1;ic<=iic-1;ic++)
    for (jc=1;jc<=jjc-1;jc++)
        p[2*ic][2*jc]=pc[ic][jc];

/* interpolate intermediate y direction */

for (i=2;i<=ii-2; i+=2)
{
    p[i][1]=(5.0*p[i][0]+15.0*p[i][2]-5.0*p[i][4]+p[i][6])*0.0625;

    for (j=3;j<=jj-3;j+=2)
        p[i][j]=(-p[i][j-3] +9.0*p[i][j-1] +9.0*p[i][j+1]-p[i][j+3])*0.0625;

    p[i][jj-1]=(5.0*p[i][jj]+15.0*p[i][jj-2]-5.0*p[i][jj-4] +p[i][jj-6])*0.0625;

}
/* interpolate in x direction */

for (j=1;j<=jj-1;j++)
{
    p[1][j]=(5.0*p[0][j]+15.0*p[2][j]-5.0*p[4][j]+p[6][j])*0.0625;

    for (i=3; i<=ii-3 ; i+=2)
        p[i][j]=(-p[i-3][j]+9.0*p[i-1][j]+9.0*p[i+1][j]-p[i+3][j])*0.0625;

    p[ii-1][j]=(5.0*p[ii][j]+15.0*p[ii-2][j]-5.0*p[ii-4][j]+p[ii-6][j])*0.0625;

}
for (i=1 ; i<=ii-1;i++)
 for (j=1;j<=jj-1;j++)
    if (p[i][j]<0) p[i][j]=0.0;
}



double conver(Stack *U, int k)
{
/* convergence check using converged solutions on level k and	*/
/* on next coarser (solution) grid k-2 */

int ic,jc ;
Level *Lc, *L ;
double **pc,**p ;
double err ;


L =U->Lk+k;

if (k==U->maxlev) p=L->p;
  else p=L->pconv;

Lc =U->Lk+k-2;
pc =Lc->pconv;

err= 0.0 ;

for (ic=1;ic<=Lc->ii-1;ic++)
    for (jc=1;jc<=Lc->jj-1;jc++)
        err+=fabs(pc[ic][jc]-p[2*ic][2*jc]) ;

return(err/((Lc->ii-1)*(Lc->jj-1))) ;
}


/********** MULTIGRID DRIVING ROUTINES **********/


void cycle(Stack *U, int k, int nu0, int nu1, int nu2, int gamma)
/* performs coarse grid correction cycle starting on level k */
/* nul pre-relaxations, nu2 postrelaxations, nuO relaxations */
/* on the coarsest grid, cycieindex gamma=l for Vcycle*/
/* gamma=2 for Wcycle*/
{
int i,j;

if (k==1)
 {
 for (i=1;i<=nu0/2;i++) relax(U,k);
  relaxh0(U,k,currfilev) ;
  for (i=1;i<=nu0/2;i++) relax(U,k);
 }
else
 {
 for (i=1;i<=nu1;i++) relax(U,k);
  coarsen_p(U,k) ;
  coarsen_f(U,k);
  for (j=1;j<=gamma;j++) cycle(U,k-2,nu0,nu1,nu2,gamma);
  refine(U,k) ;
  for (i=1;i<=nu2;i++) relax(U,k);
 }
}


void fmg(Stack *U, int k, int ks, int nu0, int nu1, int nu2, int gamma, int ncy)
{
/* performs FMG with k levels and ncy cycles per level */
int i,j ;
void init_log();

if (k==ks)
 {
 if (ks==1)
 {
 for (i=1;i<=nu0/2;i++) relax(U,k);
 relaxh0(U,k,currfilev) ;
 for (i=1;i<=nu0/2;i++) relax(U,k);
 }
else
for (j=1;j<=ncy; j++)
 {
 cycle (U,ks,nu0,nu1,nu2,gamma) ;
 printf("\n") ;
 }
}
else
 if (k>ks)
 {
 fmg(U,k-2,ks,nu0,nu1,nu2,gamma, ncy) ;
 fmg_interpolate(U,k) ;
 outlev=currfilev=k ;
 for (j=1;j<=ncy;j++)
 {
 cycle (U,k,nu0,nu1,nu2,gamma) ;
 printf("\n") ;
 }
 }
}


/********* INPUT ROUTINES * * * * * * * * * * * * * * * * */

void input_loadpar()
/* input of the load conditions for the contact */
{
printf ("\ngive M    ?") ; scanf ("%lf",&MMoes) ;
printf ("\ngive L    ?") ; scanf ("%lf",&LMoes) ;



/* conversion to parameters appearing in equations */

rlambda=1.0/(pi*pow(128.0/3.0/pow(MMoes,4.0),1.0/3.0));
alphabar=LMoes/pi*pow(1.5*MMoes,1.0/3.0);

printf("\nrlambda=%8.5e  alpha*p_h=%8.5e\n",rlambda,alphabar);

/* computation initial guess HO */

H_0=1.67*pow(MMoes,-1.0/9.0)-1.897+0.2*LMoes/50;
if (H_0<-0.99) H_0=-0.99;

printf("First approximation H0=%8.5e",H_0);

/* Optional use of hand - input HO */
/*
printf ("\ngive HO    ?");
scanf ("%lf" ,&H_O) ;
*/


/* parameters Roelands equation */

p0r=1.96e8 ;
alpha=2.2e-8 ;
eta0=40.0e-3 ;
zr= (alpha*p0r) / (log(eta0)+9.67) ;
}


void input_solvepar()
/* input of parameters numerical process */
{
/* account for intermediate grids */
printf("\nhow many levels?") ; scanf("%d",&maxl) ; maxl  = 2*maxl-1 ;
printf("\nstartlevel?")      ; scanf("%d",&starl); starl = 2*starl-1;

printf("\nhow many cycles?") ; scanf("%d",&ncy) ;
printf("\ntype of cycle ?")  ; scanf("%d",&typecy) ;


deepl=2*maxl-2;/* number of coarse grids in multi-integration */
order=6;   /* order of transfer multi-integration*/

xi_l=0.3; /* relaxation switch parameter */
urja=0.2; /* underrelaxation jacobi part */
urgs=0.4; /* underrelaxation Gauss-Seidel part */



/* factor for relaxation of force balance equation */

if (typecy==2)
 hfact=0.05;
else
hfact=0.1;
if ((LMoes>10)||(MMoes>1000)) hfact*=0.25;
}


/********** OUTPUT ROUTINES **********/


void output(Stack *U)
{
/* writes an output file of p and h */
/* output of Hm and Hc to screen */
int i,j ;
Level *L ;
double x,y,dum;
FILE *fp,*fh;

L=U->Lk+U->maxlev ;
fp=fopen("P.dat", "w") ;
fh=fopen("H.dat", "w") ;

L->Hm=1e5; /* arbitrary large value */

/* determine central film thickness */
j=L->jj/2;
i=1;
while ((L->p[i][j] > L->p[i - 1][j]) && (i < L->ii)) i++;
L->Hcp=L->hfi[i][j];


for (i=0;i<=L->ii; i++)
 {
 x=U->xa + i*L->hx ;
 for (j=0;j<=L->jj; j++)
  {
  y=U->ya+ j*L->hy  ;
  if (L->hfi[i][j]<L->Hm) L->Hm=L->hfi[i][j];
  if ((x==0)&&(y==0)) L->Hc=L->hfi[i][j];
  fprintf(fp,"%f %f %f\n",x,y,L->p[i][j]);
  fprintf(fh,"%f %f %e\n",x,y,2.0-L->hfi[i][j]);
  }
  fprintf (fp, "\n") ;
  fprintf (fh, "\n") ;
}
fclose(fp) ;
fclose(fh) ;

/* output of film thickness values */;

dum=sqrt(6*pi*rlambda) ;
printf ("\n\nLevel      Hm       Hm(Moes)     ") ;
printf ("Hc 		Hc(Moes)   Hc(Moes)(p_x=0)");

printf("\n\n");
for (i=starl;i<=U->maxlev;i+=2)
    {
    L=U->Lk+i;
    printf("%d     %8.5e %8.5e    %8.5e   %8.5e    %8.5e\n",
         (i+1)/2, L->Hm,L->Hm*dum,L->Hc,L->Hc*dum,L->Hcp*dum);
    }
printf ("\n\n") ;
}

/* +++++++++++++++++++++++Film Thickness Computation ++++++++++++++++++++*/
/* +++++++++++++++++++++++Multi Level Multi Summation ++++++++++++++++++++*/
double ah(double a, double b)
/*calculates arcsinh(a/b)*/
{
if (a==0.0) return(0.0);
else
if (b==0.0) return(1.0);
else return(log(a/b+sqrt(a*a/b/b+1)));
}


void init_log(Stack *U, int l)
{
/* computes the kernel on level 1 */
int i,j;
Level *L;
double xp,yp,xm,ym;

L =U->Lk+l;

for (i=0; i<=L->ii+4*U->od; i++)
 {
 xp= (i+0.5)*L->hx ; xm=xp-L->hx;
 for (j=0; j<=L->jj+4*U->od; j++)
  {
  yp=(j+0.5)*L->hy ; ym=yp-L->hy ;
  L->K[i][j]=
      fabs(xp)*log(yp/xp + sqrt(1.0+yp*yp/xp/xp))
	- fabs(xp)*log(ym/xp + sqrt(1.0+ym*ym/xp/xp))
	+ fabs(xm)*log(ym/xm + sqrt(1.0+ym*ym/xm/xm))
	- fabs(xm)*log(yp/xm + sqrt(1.0+yp*yp/xm/xm))
	+ fabs(yp)*log(xp/yp + sqrt(1.0+xp*xp/yp/yp))
	- fabs(yp)*log(xm/yp + sqrt(1.0+xm*xm/yp/yp))
	+ fabs(ym)*log(xm/ym + sqrt(1.0+xm*xm/ym/ym))
	- fabs(ym)*log(xp/ym + sqrt(1.0+xp*xp/ym/ym));
	}
  }
}


double kval(Stack *U, int lf, int lc, int i1, int j1)
/* calculates the value of k on the fine grid if in the point il,jl*/
/* and injected to the coarse grid ic */
{
int l,i,j;
double x1,x2,y1,y2,help;
Level *Lf, *Lc;

Lf=U->Lk+lf;
Lc=U->Lk+lc;

i=i1; j=j1;
for (l=lf-1; l>=lc; l--)
  if ((2*(l/2)-l)!=0) i *= 2; else j *= 2;
x1=(i-0.5)*Lf->hx; x2=x1+Lf->hx; y1=(j-0.5)*Lf->hy; y2=y1+Lf->hy;


if ((x1==0.0) || (x2==0.0) || (y1==0.0) || (y2==0.0))
    printf("error in kval lf=%d lc=%d i=%d j=%d",lf,lc,i1,j1);
help=fabs(x2)*ah(y2,x2)+fabs(y2)*ah(x2,y2)-fabs(x2)*ah(y1,x2)-fabs(y2)*ah(x1,y2)
    -fabs(x1)*ah(y2,x1)-fabs(y1)*ah(x2,y1)+fabs(x1)*ah(y1,x1)+fabs(y1)*ah(x1,y1);
 return (help*Lc->hx*Lc->hy/Lf->hx/Lf->hy) ;
}


void fillk(Stack *U, int l)
/* fills the matrix K with integral values*/
/* for matrix multiplication to obtain the integral*/
{
int l1, i,j;
Level *L, *Lc;

for (l1=l; l1>=2; l1--)
 {
  L =U->Lk+l1 ;
  Lc=U->Lk+l1-1 ;
  if ((2*(l1/2)-l1)==0)
   for (j=0; j<=L->jj+4*U->od; j++)
   {
   for (i=0; i<=Lc->ii+2*U->od; i++) Lc->K[i][j]=2.0*L->K[2*i][j];
   for (i=Lc->ii+2*U->od+1;i<=Lc->ii+4*U->od; i++) Lc->K[i][j]=kval(U,l,l1-1,i,j);
   }
   else
   for (i=0; i<=L->ii+4*U->od; i++)
   {
	   for (j=0; j<=Lc->jj+2*U->od; j++) Lc->K[i][j] = 2.0*L->K[i][2*j];
       for (j=Lc->jj+2*U->od+1; j<=Lc->jj+4*U->od; j++) Lc->K[i][j]=kval(U,l,l1-1,i,j);
   }
 }
}


void calcku(Stack *U, int l)
/*calculates the integral on grid I*/
{
int i1,j1,i,j;
double help;
Level *L;

L=U->Lk+l ;
for (i1=-U->od; i1<=L->ii+U->od; i1++)
 for (j1=-U->od; j1<=L->jj+U->od; j1++)
 {
 help=0.0;
 for (i=-U->od; i<=L->ii+U->od; i++)
   for (j=-U->od; j<=L->jj+U->od; j++)
   help+=L->K[abs(i-i1)][abs(j-j1)]*L->p[i][j] ;
 L->w[i1][j1]=help;
 }
}



void sto6k1(Stack *U, int l)
/* stores values of k1 for correct6	*/
/* ki is the value of k minus the 6th */
/* order coarse grid approximation to k*/
{
	int i,j,maxx,maxy;
	Level *L;

	L=U->Lk+l;

	if (4*U->m1<L->ii-1) maxx = 4*U->m1+1; else maxx = L->ii;
	if (4*U->m1<L->jj-1) maxy = 4*U->m1+1; else maxy = L->jj;


	if ((2*(l/2)-l)!= 0)
		for (i=0; i<=maxx + 4*U->od; i++)
		{
			for (j=0; j<=maxy +2*U->od; j++)
			{
				L->K1[i][j] = L->K[i][j]
					- (3 * L->K[i][abs(5-j)] -25 * L->K[i][abs(3-j)]
					+ 150* L->K[i][abs(1-j)] +150* L->K[i][abs(1+j)]
					- 25 * L->K[i][abs(3+j)] +3  * L->K[i][abs(5+j)]) / 256;
			}

			for (j=maxy+2*U->od+1; j<=maxy+4*U->od; j++)
			{
				L->K1[i][j] = L->K[i][j]
					- (3 * L->K[i][abs(5-j)] -  25 * L->K[i][abs(3 - j)]
					+ 150* L->K[i][abs(1-j)] + 150 * kval(U,U->maxlev, l,i,j+1)
					- 25 * kval(U,U->maxlev,l,i,j+3) + 3*kval(U,U->maxlev,l,i,j+5))/256;
			}
		}
	else
	{
		for (j=0; j<=maxy+4*U->od; j++)
		{
			for (i=0; i<=maxx+2*U->od; i++)
			{
				L->K1[i][j] = L->K[i][j] -
					(  3 * L->K[abs(5 - i)][j] - 25 * L->K[abs(3 - i)][j]
					+ 150* L->K[abs(1 - i)][j] + 150* L->K[abs(1 + i)][j]
					- 25 * L->K[abs(3 + i)][j] + 3  * L->K[abs(5 + i)][j])/256;
			}
			for (i=maxx+2*U->od+1; i <= maxx+4*U->od; i++)
			{
				L->K1[i][j] = L->K[i][j] -
					(3 * L->K[abs(5 - i)][j] - 25 * L->K[abs(3-i)][j]
                 + 150 * L->K[abs(1 - i)][j] + 150 * kval(U, U->maxlev, l, i + 1, j)
					-25* kval(U, U->maxlev, l, i+3,j) + 3*kval(U, U->maxlev, l, i+5,j)) / 256;
			}
		}
	}
}
void coarsenp6x(Stack *U, int l)
{
		/* 6th order weighting in x direction */
		/* actually the 6th order approximation is obtained in K */
		/* the stencil is 1/512*(3,0,-25,0,150,256,150,0,-25,0,3) in x direction */


		int i, j;
		double q0, q1, q2, q3, q4, q5, q7, q9;
		Level *L, *Lc;

		L  = U->Lk+l;
		Lc = U->Lk+l-1;

		if ((2*(l/2)-l)!= 0) printf("\nwrong coarsening in coarsenp6x\n");

		for (j=-4; j<= L->jj+4;j++)
		{
			q0 = L->p[-4][j]; q1 = L->p[-3][j]; q2 = L->p[-2][j]; q3 = L->p[-1][j];
			q4 = L->p[ 0][j]; q5 = L->p[ 1][j]; q7 = L->p[ 3][j]; q9 = L->p[ 5][j];

			Lc->p[-4][j] = (                                                      3 * q1) / 512;
			Lc->p[-3][j] = (                                           -25 * q1 + 3 * q3) / 512;
			Lc->p[-2][j] = (                      256 * q0 + 150 * q1 - 25 * q3 + 3 * q5) / 512;
			Lc->p[-1][j] = (           150 * q1 + 256 * q2 + 150 * q3 - 25 * q5 + 3 * q7) / 512;
			Lc->p[ 0][j] = (-25 * q1 + 150 * q3 + 256 * q4 + 150 * q5 - 25 * q7 + 3 * q9) / 512;


			for (i=1;i<=Lc->ii-1; i++)
				Lc->p[i][j] = (3*L->p[2*i-5][j] -25*L->p[2*i-3][j] + 150*L->p[2*i-1][j]
				           + 256*L->p[2*i  ][j]
				             + 3*L->p[2*i+5][j] -25*L->p[2*i+3][j] + 150*L->p[2*i+1][j]) / 512;

			q0 = L->p[L->ii+ 4][j]; q1 = L->p[L->ii+3][j]; q2 = L->p[L->ii+2][j]; q3 = L->p[L->ii+1][j];
			q4 = L->p[L->ii   ][j]; q5 = L->p[L->ii-1][j]; q7 = L->p[L->ii-3][j]; q9 = L->p[L->ii-5][j];

			Lc->p[Lc->ii    ][j] = (-25 * q1 + 150 * q3 + 256 * q4 + 150 * q5 - 25 * q7 + 3 * q9) / 512;
			Lc->p[Lc->ii + 1][j] = (           150 * q1 + 256 * q2 + 150 * q3 - 25 * q5 + 3 * q7) / 512;
			Lc->p[Lc->ii + 2][j] = (                      256 * q0 + 150 * q1 - 25 * q3 + 3 * q5) / 512;
			Lc->p[Lc->ii + 3][j] = (                                           -25 * q1 + 3 * q3) / 512;
			Lc->p[Lc->ii + 4][j] = (                                                      3 * q1) / 512;
		}
}



void coarsenp6y(Stack *U, int l)
/*6th order weighting, in y-direction */
/*actually the 6th order approximation is obtained in K */
/*the stencil is 1/512*93,0,-25,0,150,256,150,0,-25,0,3) in y */
/*direction */

{
int i,j;
double q0,q1,q2,q3,q4,q5,q7,q9;
Level *L, *Lc;


L =U->Lk+l;
Lc=U->Lk+l-1;

if ((2*(l/2)-l)==0) printf("\nwarning:wrong coarsening in coarsenp6y\n") ;

for (i=-4; i<=L->ii+4; i++)
 {
 q0=L->p[i][-4]; q1=L->p[i][-3]; q2=L->p[i][-2]; q3=L->p[i][-1];
 q4=L->p[i][ 0]; q5=L->p[i][ 1]; q7=L->p[i][ 3]; q9=L->p[i][ 5];

 Lc->p[i][-4]=(                                        3*q1)/512;
 Lc->p[i][-3]=(                                -25*q1 +3*q3)/512;
 Lc->p[i][-2]=(                 256*q0 +150*q1 -25*q3 +3*q5)/512;
 Lc->p[i][-1]=(         150*q1 +256*q2 +150*q3 -25*q5 +3*q7)/512;
 Lc->p[i][ 0]=( -25*q1 +150*q3 +256*q4 +150*q5 -25*q7 +3*q9)/512;


  for (j=1; j<=Lc->jj-1;j++)
  {
	Lc->p[i][j]=(	   3*L->p[i][2*j-5]-25*L->p[i][2*j-3]+150*L->p[i][2*j-1]
					+256*L->p[i][2*j  ]
					  +3*L->p[i][2*j+5]-25*L->p[i][2*j+3]+150*L->p[i][2*j+1])/512;
  }

q0=L->p[i][L->jj+4] ; q1=L->p[i][L->jj+3] ; q2=L->p[i][L->jj+2] ; q3=L->p[i][L->jj+1] ;
q4=L->p[i][L->jj ]  ; q5=L->p[i][L->jj-1] ; q7=L->p[i][L->jj-3];  q9=L->p[i][L->jj-5] ;

Lc->p[i][Lc->jj  ]=( -25*q1 +150*q3 +256*q4 +150*q5-25*q7 +3*q9)/512;
Lc->p[i][Lc->jj+1]=(         150*q1 +256*q2 +150*q3-25*q5 +3*q7)/512;
Lc->p[i][Lc->jj+2]=(                 256*q0 +150*q1-25*q3 +3*q5)/512;
Lc->p[i][Lc->jj+3]=(                               -25*q1 +3*q3)/512;
Lc->p[i][Lc->jj+4]=(                                       3*q1)/512;
 }
}


void refine6x(Stack *U, int l)
/* 0(6) correction of values of hfi on the fine grid afterinterpolation */
/* in x direction */
{
int i,j,i1,j1,ibe,ien,jbe,jen;
double help;
Level *L, *Lc ;

L =U->Lk+l;
Lc=U->Lk+l-1;

if ((2*(l/2)-l)!=0) printf("\nwarning:error in refine6x: l=%d\n",l);

for (i=-4;i<=Lc->ii+4; i++) /* fill also the ghost points */
{
ibe=2*i-U->m1; if (ibe<-4)      ibe=-4;      if ((2*(ibe/2)-ibe)==0) ibe++;
ien=2*i+U->m1; if (ien>L->ii+4) ien=L->ii+4; if ((2*(ien/2)-ien)==0) ien--;
for (j=-4;j<=L->jj+4;j++)
 {
 jbe=j-U->m2; if (jbe<-4)      jbe=-4;
 jen=j+U->m2; if (jen>L->jj+4) jen=L->jj+4;
 help=0.0 ;
 for (i1=ibe; i1<=ien;i1 += 2)
  for (j1=jbe; j1<=jen; j1++)
   help += L->p[i1][j1]*L->K1[abs(2*i-i1)][abs(j-j1)];
 Lc->w[i][j] += help;
  }
}


for (i=-2; i<=Lc->ii+2; i++) /* fill also the ghost points */
    for (j=-4; j<=L->jj+4; j++)
       L->w[2*i][j]=Lc->w[i][j];


for (i=-2; i<=Lc->ii+1; i++)
{
ibe=2*i+1-U->m1; if (ibe<-4)      ibe=-4;
ien=2*i+1+U->m1; if (ien>L->ii+4) ien=L->ii+4;
for (j=-4; j<=L->jj+4; j++)
{
jbe=j-U->m2; if (jbe<-4)      jbe=-4;
jen=j+U->m2; if (jen>L->jj+4) jen=L->jj+4;
help=0.0 ;

for (i1=ibe; i1<=ien; i1++)
for (j1=jbe; j1<=jen; j1++)
help += L->p[i1][j1]*L->K1[abs(2*i+1-i1)][abs(j-j1)];
L->w[2*i+1][j]=( 3*Lc->w[i-2][j]-  25*Lc->w[i-1][j]
              +150*Lc->w[i  ][j]+ 150*Lc->w[i+1][j]
              - 25*Lc->w[i+2][j]+   3*Lc->w[i+3][j])/256+help;
  }
 }
}

void refine6y(Stack *U, int l)
/* 0(6) correction of values of WI on the fine grid after interpolation*/
{
int i,j,i1,j1,ibe,ien,jbe,jen;
double help;
Level *L, *Lc;

L =U->Lk+l;
Lc=U->Lk+l-1;

if ((2*(l/2)-l)==0) printf("\nwarning:error in refine6y: %d\n",l);

for (j=-4; j<=Lc->jj+4;j++)
 {
 jbe=2*j-U->m1 ; if (jbe<-4)      jbe=-4     ; if ((2*(jbe/2)-jbe)==0) jbe++;
 jen=2*j+U->m1 ; if (jen>L->jj+4) jen=L->jj+4; if ((2*(jen/2)-jen)==0) jen--;
 for (i=-4; i<=L->ii+4; i++)
 {
 ibe=i-U->m2; if (ibe<-4)      ibe=-4;
 ien=i+U->m2; if (ien>L->ii+4) ien=L->ii+4;
 j1=jbe; help=0.0 ;

for (j1=jbe; j1<=jen; j1 += 2)
    for (i1=ibe; i1<=ien; i1++)
        help += L->p[i1][j1]*L->K1[abs(i-i1)][abs(2*j-j1)] ;
Lc->w[i][j] += help;
 }
}

for (i=-4; i<=L->ii+4; i++)
    for (j=-2; j<=Lc->jj+2; j++)
        L->w[i][2*j]=Lc->w[i][j];

for (j=-2; j<=Lc->jj+1; j++)
{
jbe=2*j+1-U->m1; if (jbe<-4) jbe=-4;
jen=2*j+1+U->m1; if (jen>L->jj+4) jen=L->jj+4;
for (i=-4; i<=L->ii+4; i++)
{
    ibe=i-U->m2; if (ibe<-4) ibe=-4;
    ien=i+U->m2; if (ien>L->ii+4) ien=L->ii+4;
    help=0.0 ;
    for (j1=jbe; j1<=jen; j1++)
        for (i1=ibe; i1<=ien; i1++)
            help+=L->p[i1][j1]*L->K1[abs(i-i1)][abs(2*j+1-j1)] ;
    L->w[i][2*j+1]=
            (  3*Lc->w[i][j-2]-  25*Lc->w[i][j-1]
            +150*Lc->w[i][j  ]+ 150*Lc->w[i][j+1]
            - 25*Lc->w[i][j+2]+   3*Lc->w[i][j+3])/256+help;
  }
 }
}


void calchi(Stack *U, int lev)
/* calculates film thickness on level l */
{
	int i, j, ldeep, l;
	Level *L;

	ldeep = U->deep;
	if (ldeep<0) ldeep=0;
	if (ldeep>lev-1) ldeep=lev-1;
	fillk(U, lev);
	if (U->od == 4) { for (l = lev; l >= lev - ldeep + 1; l--) sto6k1(U, l); }
    for (l=1;l<=ldeep; l+=2)
	{
		if (U->od == 4) { coarsenp6y(U,lev-l+1); coarsenp6x(U,lev-l);}
	}
	calcku(U,lev-ldeep);
	for (l=ldeep;l>=1;l-=2)
	{
		if (U->od == 4) { refine6x(U,lev-l+1); refine6y(U,lev-l+2);}
	}

	L = U->Lk+lev;

		/* compute extra point for 2nd order discretization at first line */

		for (i=-1; i<=L->ii;i++)
			for (j=0; j<=L->jj; j++)
				L->hfi[i][j] = U->h0+2.0/pi/pi*L->w[i][j]+L->hrhs[i][j];

}


/********** ROUTINES FOR LINE RELAXATION **********/


void init_line(Stack *U, int k, int j, int ml, int mr)
/* prepares system to be solved */
{
double H3, Hx, Qx, Qy, xi_n,xi_s,xi_e,xi_w, r0, r1, r2;
double dHx,dHxm,dHxp,dHxmm,dHxpp,dHxm3;
double dK0, dK1, dK2,dK3, dK4, rhx, rhy, rhx2, rhy2, rpi2;
double **P,**Pj, **H, **f, **A, *Y, *X, **K;
Level *L;
int ii,i,m,gs;

L=U->Lk+k;
ii=L->ii;

rhx=1.0/L->hx ;
rhy=1.0/L->hy ;
rhx2 = rhx*rhx ;
rhy2 = rhy*rhy;
rpi2 = 1.0/(pi*pi) ;

P=L->p ;
Pj=L->pjac ;
K=L->K ;
H=L->hfi ;
f=L->f ;
A=L->A ;
Y=L->Y ;
X=L->X ;

for (i=0;i<=L->ii-2;i++)
{
  for (m=0;m<=3*ml+mr ;m++)
     A[i][m] =0.0,
  X[i]=0.0;
  Y[i]=0.0;
}

for (i=1;i<=ii-1;i++)
{
    r0=rho(P[i][j]);
    H3=r0*H[i][j]*H[i][j]*H[i][j]*reta(P[i][j]);

    xi_n=0.5*(rho(P[i  ][j+1]) *H[i ][j+1]*H[i ][j+1]*H[i ][j+1]
            *reta(P[i  ][j+1]) +H3)*rlambda*rhy2 ;
    xi_s=0.5*(rho(P[i  ][j-1]) *H[i ][j-1]*H[i ][j-1]*H[i ][j-1]
            *reta(P[i  ][j-1]) +H3)*rlambda*rhy2 ;
    xi_e=0.5*(rho(P[i+1][j  ]) *H[i+1][j ]*H[i+1][j ]*H[i+1][j ]
            *reta(P[i+1][j  ]) +H3)*rlambda*rhx2 ;
    xi_w=0.5*(rho(P[i-1][j  ]) *H[i-1][j ]*H[i-1][j ]*H[i-1] [j ]
            *reta(P[i-1][j  ]) +H3)*rlambda*rhx2;

gs= ((fabs(xi_n)>xi_l)||(fabs(xi_w) >xi_l )||(fabs (xi_e)>xi_l)||(fabs(xi_s)>xi_l) ) ;

if (gs==1)
{
dK0=K[0][0];
dK1=K[1][0] ;
dK2=K[2][0] ;
dK3=K[3][0] ;
dK4=K[4][0] ;
}
else
{
dK0=K[0][0]-0.25*(2*K[1][0] +2*K[0][1] ) ;
dK1=K[1][0]-0.25*(K[2][0]+K[0][0]+2*K[1][1]);
dK2=K[2][0]-0.25*(K[3][0]+K[1][0]+2*K[2][1]);
dK3=K[3][0]-0.25*(K[4][0]+K[2][0]+2*K[3][1]);
dK4=K[4][0]-0.25*(K[5][0]+K[3][0]+2*K[4][1]);
}
r1=rho(P[i-1][j]) ;
if (i==1) Hx=rhx*(r0*H[i][j]-r1*H[i-1][j]);
else
{
r2=rho(P[i-2][j]);
Hx=rhx*(1.5*r0*H[i][j]-2*r1*H[i-1][j]+0.5*r2*H[i-2][j]);
}


dHx=dHxm=dHxmm=dHxm3=0.0;
dHxp = 0.0;
dHxpp  = 0.0;

if (i==1)
{
dHxm3=0.0*rhx*(r0*dK3-r1*dK2);
dHxmm=0.0*rhx*(r0*dK2-r1*dK1);
dHxm= 0.0*rhx*(r0*dK1-r1*dK0);
dHx=  rhx*(r0*dK0-r1*dK1);
dHxp= rhx*(r0*dK1-r1*dK2);
dHxpp=rhx*(r0*dK2-r1*dK3);
}
else
{
if (i>3) dHxm3=rhx*(1.5*r0*dK3-2.0*r1*dK2+0.5*r2*dK1);
if (i>2) dHxmm=rhx*(1.5*r0*dK2-2.0*r1*dK1+0.5*r2*dK0);
if (i>1) dHxm= rhx*(1.5*r0*dK1-2.0*r1*dK0+0.5*r2*dK1);
dHx= rhx*(1.5*r0*dK0-2.0*r1*dK1+0.5*r2*dK2) ;
if (i<ii-1) dHxp = rhx*(1.5*r0*dK1-2.0*r1*dK2+0.5*r2*dK3);
if (i<ii-2) dHxpp= rhx*(1.5*r0*dK2-2.0*r1*dK3+0.5*r2*dK4) ;
}

if (gs==1)
{
 Qx=(xi_e*Pj[i+1][j  ] -(xi_e+xi_w)*Pj[i][j] + xi_w*Pj[i-1][j ]) ;
 Qy=(xi_n*Pj[i  ][j+1] -(xi_n+xi_s)*Pj[i][j] + xi_s*Pj[i  ][j-1]);


Y[i-1]=f[i][j] -Qx-Qy + Hx;
A[i-1][ml] =-(xi_w+xi_e)-(xi_s+xi_n)-2*rpi2*dHx ;

if ((Pj[i-2][j]>0)&&(Pj[i-1][j]>0)&&(Pj[i+1][j]>0)&&(Pj[i+2][j]>0))
{
if (i>1)    A[i-1][ml-1] = xi_w-2*rpi2*dHxm;
if (i>2)    A[i-1][ml-2] = -2*rpi2*dHxmm;
if (i>3)    A[i-1][ml-3] = -2*rpi2*dHxm3;
if (i<ii-1) A[i-1][ml+1] =  xi_e-2*rpi2*dHxp;
if (i<ii-2) A[i-1][ml+2] = -2*rpi2*dHxpp;
 }
}
    else
    {
    Qx=(xi_e*P[i+1][j  ]  -(xi_e+xi_w)*P[i][j] + xi_w*P[i-1][j ]);
    Qy=(xi_n*P[i  ][j+1]  -(xi_n+xi_s)*P[i][j] + xi_s*P[i  ][j-1]);


    Y[i-1]=f[i][j] -Qx -Qy + Hx;
    A[i-1][ml] =-1.25*(xi_w+xi_e+xi_s+xi_w)-2*rpi2*dHx ;
    if((Pj[i-2][j]>0)&&(Pj[i-1][j]>0)&&(Pj[i+1][j]>0)&&(Pj[i+2][j]>0))
    {
    if (i>1)    A[i-1][ml-1]= xi_w+(xi_n+xi_e+xi_e+xi_w)*0.25-2*rpi2*dHxm;
    if (i>2)    A[i-1][ml-2]=-0.25*xi_w-2*rpi2*dHxmm;
    if (i>3)    A[i-1][ml-3]=-2*rpi2*dHxm3 ;
    if (i<ii-1) A[i-1][ml+1]= xi_e+(xi_n+xi_e+xi_s+xi_w)*0.25-2*rpi2*dHxp;
    if (i<ii-2) A[i-1][ml+2]=-0.25*xi_e-2*rpi2*dHxpp;
    }
   }
}
}

int pivot(double **A, int ri, int ii, int ml)
{
int i,j ,I,end;
double piv ;

piv=fabs(A[ri][ml]) ;
I=ri ;

end=ml; if ((ri+ml)>ii)  end=ii-ri;
for (i=1;i<=end; i++)
{
j=ml-i ;
  if (fabs(A[ri+i][j]) > piv)
  I=ri+i ;
}
return (I) ;
}



void swap(double **A, double *y, int ri, int rI, int ml, int m)
{
int j,cJ ;
double temp ;

for (j=ml;j<=m;j++)
   {
   cJ=j-(rI-ri) ;
   if (cJ<0) continue;
   temp=A[ri][j]; A[ri][j]=A[rI][cJ]; A[rI][cJ]=temp;
   }

temp=y[ri]; y[ri]=y[rI]; y[rI]=temp;
}



void bcksb(double **A,double *Y,double *X, int mm, int ml, int N)
{
int i,j;
double sum;

X[N]=Y[N]/(A[N][ml]) ;
for (i=N-1;i>=0;i--)
{
sum=0.0 ;
for (j=1 ; j<=mm; j++){
 if(i+j<=N) sum+=A[i][ml+j]*X[i+j];
 }
 X[i]=(Y[i]-sum)/(A[i][ml]) ;
 }
}

void solve_system(Stack *U, int l, int m, int ml, int mr)
{
int N,end,i,j,k,rI,J;
double dum;
Level *L;
double **A, *X, *Y;

L=U->Lk+l;
A=L->A;
X=L->X;
Y=L->Y;
N=L->ii-2;

for (i=0;i<=N;i++)
 {
    rI = pivot(A, i, N, ml);
    if (i!=rI) swap(A,Y,i, rI, ml,m);

end = i+ml ; if ( end>N ) end = N ;
for (k=i+1; k<=end; k++)
 {
 J=ml-(k-i) ;
 if (J<0) continue ;
 dum = -(A[k][J]/A[i][ml]);
 for(j=J+1;j<=J+(mr+ml);j++)
     A[k][j] += dum*A[i][j+(ml-J)];
 A[k][J] = 0.0 ;
 Y[k]+=dum*Y[i] ;
 }
}
bcksb(A,Y,X,m,ml,N);
}


void return_line(Stack *U, int k , int j)
{
double **H, H3,xi_n,xi_s,xi_e,xi_w,r0;
int ii,jj,i,gs;
double rhx,rhy,rhx2,rhy2,*X,del0,p0;
Level *L ;

L=U->Lk+k;
ii=L->ii;
jj=L->jj;

rhx=1.0/L->hx;
rhy=1.0/L->hy;
rhx2=rhx*rhx;
rhy2=rhy*rhy;
X=L->X;
H=L->hfi;



for (i=1; i<=ii-1;i++)
   {
   r0=rho(L->p[i][j]);
   H3=r0*H[i][j]*H[i][j]*H[i][j] *reta(L->p[i][j]);

   xi_n=0.5*(rho(L->p[i  ][j+1]) *H[i   ][j+1]*H[i   ][j+1]*H[i   ][j+1]
          * reta(L->p[i  ][j+1]) +H3) * rlambda*rhy2 ;
   xi_s=0.5*(rho(L->p[i  ][j-1]) *H[i   ][j-1]*H[i   ][j-1]*H[i   ][j-1]
          * reta(L->p[i  ][j-1]) +H3) * rlambda*rhy2 ;
   xi_e=0.5*(rho(L->p[i+1][j  ]) *H[i+1][j   ]*H[i+1][j   ]*H[i+1][j   ]
          * reta(L->p[i+1][j  ]) +H3) * rlambda*rhx2 ;
   xi_w=0.5*(rho(L->p[i-1][j  ]) *H[i-1][j   ]*H[i-1][j   ]*H[i-1][j   ]
          * reta(L->p[i-1][j  ]) +H3) * rlambda*rhx2 ;

gs=((fabs(xi_n)>xi_l)||(fabs(xi_w)>xi_l)||(fabs(xi_e)>xi_l)||(fabs(xi_s)>xi_l));

if(gs==1)
{
p0=L->pjac[i][j] ;
L->pjac[i][j] +=urgs*X[i-1] ;
if (L->pjac[i][j]<0.0) L->pjac[i][j]=0.0;
del0=L->pjac[i][j]-p0;
}
else
{
p0=L->pjac[i][j] ;
L->pjac[i][j]+=urja*X[i-1] ;
if (L->pjac[i][j]<0) L->pjac[i][j]=0.0;
del0=L->pjac[i][j]-p0;
if ((L->pjac[i][j]>0)&&(L->pjac[i-1][j]>0)&&(L->pjac[i+1][j]>0)
                     &&(L->pjac[i][j-1]>0)&&(L->pjac[i][j+1]>0))
{
if (i>1)
{
L->pjac[i-1][j] -=0.25*del0;
if (L->pjac[i-1][j]<0) L->pjac[i-1][j]=0.0;
}
if (i<ii-1)
{
L->pjac[i+1][j] -=0.25*del0;
if (L->pjac[i+1][j]<0) L->pjac[i+1][j]=0.0;
}
if (j>1)
{
L->pjac [i][j-1]-=0.25*del0;
if (L->pjac[i][j-1]<0) L->pjac[i][j-1]=0.0;
}
if (j<jj-1)
{
L->pjac [i][j+1]-=0.25*del0;
if (L->pjac[i][j+1]<0) L->pjac[i][j+1]=0.0;
}
}
}
}
}

void solve_line (Stack *U, int k, int j)
{
int ml,m,mr ;
ml=3;
mr=2 ;
m=2*ml+mr ;
init_line(U, k, j, ml,mr);
solve_system(U ,k ,m,ml ,mr) ;
return_line (U,k, j) ;
}
/********** MAIN PROGRAM **********/

void main()
{

int l ;
Stack U;

input_loadpar();
input_solvepar();

initialize(&U,16,16,maxl,deepl,order,-4.5,1.5,-3.0,3.0,H_0);

for (l=1; l<=maxl; l+=2) { init_log(&U,l); init_f(&U,l) ; }
outlev=currfilev=starl ;
t0=clock();
init_p (&U, starl) ;


fmg (&U, maxl, starl, 40,2, 1, typecy, ncy) ;
cputimes[maxl] =clock() ;

output(&U) ;


if (maxl>1){
 for (l=starl+2;l<=maxl;l+=2) printf(" aen(%2d %2d)=%8.5e \n",(l+1)/2,(l-1)/2,conver (&U,l)) ;
}

printf ("\n\n") ;
for (l=starl;l<=maxl;l+=2) printf("cpu times level %d = %8.5e \n",(l+1)/2,
                                  (cputimes [l]-t0)/1000000.0) ;

finalize (&U, maxl) ;
}

































































