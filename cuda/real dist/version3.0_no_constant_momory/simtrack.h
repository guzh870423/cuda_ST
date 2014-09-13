#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cdf.h>

using namespace std;
//================================================================
//
//  Global  constants 
//
//===============================================================
double PI      = 3.14159265;
double light_speed =  299792458;
double proton_mass =  938.272046;     // in units of MeV
double electron_mass =  0.510998928;  // in units of MeV
double Fdrift1 = 1.e0 / (2.e0 * (2.e0 - pow(2.e0, 1.0/3.0 )));
double Fdrift2 = 0.5e0 - Fdrift1;
double Fkick1  = 2.e0 * Fdrift1;
double Fkick2  = 1.e0 - 2e0 * Fkick1;
double QLslice = 0.02;
double BLslice = 0.02;
enum ps_index { x_ = 0, px_ = 1, y_ = 2, py_ = 3,  z_ = 4, delta_ = 5 };  //  z > 0 earlier than synchronous particle
double seed=1232346;

//================================================================
//
//      Global variables and flags
//
//===============================================================

//-----global beam and tracking parameters
class  GlobalVariables
{
 public:
  GlobalVariables()
    {
      energy=250000.;          // particle energy of one nuclei in MeV
      gamma =266;              // beambeam needs it
      beta = 0.999955981502;   // particle volecity dived by light speed = sqrt(1-1.0/gamma/gamma) 
      brho = 831.763013;       // SOLEN needs it  brho= e/p_0
      harm = 360;             // for accelerating RF cavities
      A= 1;                   // Atom index  197 
      Q= 1;                   // Particle charge   79 
      bbscale =1 ;             // force scaling for beambeam, compared to p-p collision, 
                               // for example: bbscale =1.0 for proton-proton,  bbscale =79*79/197=31.68 for Au-Au,.. 
      turn = 0 ;               // for tracking turn control ( used by ACKICK, ACDIPOLE, ACMULT )
      step_deltap=0.0001;      // used in numeric chrom calculations
      twiss_6d = 0;            // for 5D or 6D Twiss calculation  , 0 will treat RF as drift in Twiss calculation      
    }
  int    turn, A, Q, harm, twiss_6d ;             
  double energy, gamma, beta, brho, bbscale, step_deltap;
};
GlobalVariables GP; 

void Print_GlobalVariables()
{
  cout<<"-------------------------------"<<endl;
  cout<<"GP.energy  :  "<<GP.energy<<endl;
  cout<<"GP.gamma   :  "<<GP.gamma<<endl;
  cout<<"GP.beta    :  "<<GP.beta<<endl;
  cout<<"GP.brho    :  "<<GP.brho<<endl;
  cout<<"GP.harm    :  "<<GP.harm<<endl;
  cout<<"GP.A       :  "<<GP.A<<endl;
  cout<<"GP.Q       :  "<<GP.Q<<endl;
  cout<<"GP.bbscale :  "<<GP.bbscale<<endl;
  cout<<"GP.turn    :  "<<GP.turn<<endl; 
  cout<<"GP.step_deltap :  "<<GP.step_deltap<<endl;
  cout<<"GP.twiss_6d:  "<<GP.twiss_6d<<endl;
  cout<<"-------------------------------"<<endl;
}

//---containers for elements for Track_Fast()

int     nelement;
int     Ndrift, Nsbend, Nquad,  Nsext,     Nmult;
int     Nrf,    Nbb,    Nelens, Ncooling,  Nacmult, Nmatrix;
int     Nkicker, Nhkicker, Nvkicker, Ndiffuse, Ngmult, Nlrbb;

vector <int>     type;
vector <double>  nelement_dx, nelement_dy, nelement_dt, nelement_APx, nelement_APy;

vector <double>  drift_l;

vector <int>     sbend_nint;
vector <double>  sbend_l, sbend_angle, sbend_e1, sbend_e2;

vector <int>     quad_nint;
vector <double>  quad_l,  quad_k1l;

vector <int>     sext_nint;
vector <double>  sext_l,  sext_k2l;

vector <int>     mult_nint, mult_norder;
vector <double>  mult_l, mult_knl,  mult_knsl;

vector <double>  bb_np,  bb_emitx,  bb_emity, bb_sigmal, bb_betax, bb_alfax, bb_betay, bb_alfay;
vector <int>     bb_treat, bb_nslice;

vector <int>     elens_nslice;
vector <double>  elens_l, elens_ne, elens_betae, elens_sigmax, elens_sigmay;

vector <double>  rf_l, rf_vrf,  rf_frf,  rf_phi0;

vector <double>  cooling_l, cooling_alpha;

vector <int>     acmult_norder, acmult_tturns;
vector <double>  acmult_l, acmult_kl, acmult_phi0; 

vector <double>  matrix_l, matrix_m66, matrix_xco_in, matrix_xco_out;

vector <double>  kick_l, kicker_hkick, kicker_vkick;

vector <double>  hkick_l, hkicker_hkick;

vector <double>  vkick_l, vkicker_vkick;

vector <double>  diffuse_l, diffuse_diff_x, diffuse_diff_y, diffuse_diff_delta;

vector <int>     gmult_nint, gmult_norder;
vector <double>  gmult_l,gmult_angle, gmult_knl,  gmult_knsl, gmult_e1, gmult_e2;

vector <double>  lrbb_np, lrbb_sepx, lrbb_sepy, lrbb_sigmax, lrbb_sigmay, lrbb_kickx0, lrbb_kicky0;

//===========================================
//
//          math functions
//
//===========================================
long int  fac(int n)
{
  if (n<=0) {
    cout<<"n should be non-negative."<<endl;
    exit(1);
  }
  if(n==1) return 1;
  if(n==2) return 2;
  if(n==3) return 6;
  if(n==4) return 24;
  if(n==5) return 120;
  if(n==6) return 720;
  if(n==7) return 5040;
  if(n==8) return 40320;
  if(n==9) return 362880;
  if(n==10) return 3628800;  
  if(n >10) {
    cout<<" fac(n), n>10. exit."<<endl;
    exit(1);
  }   
  return 0;
}

//---random number generator
double rnd(double & r)
{
  int m;
  double s,u,v,p;
  s=65536.0; u=2053.0; v=13849.0;
  m=int(r/s);  r=r-m*s;  r=u*r+v; 
  m=(int)(r/s); r=r-m*s; p=r/s;
  return(p);
}

double guass(double u,double g, double & r)
{ int i,m;
 double s,w,v,t;
 s=65536.0; w=2053.0; v=13849.0;
 t=0.0;
 for (i=1; i<=12; i++)
   { r=r*w+v; m=(int)(r/s);
   r=r-m*s; t=t+(r)/s;
   }
 t=u+g*(t-6.0);
 return(t);
}

//---polynominal fitting
void pfit(double x[],double y[], int n, double a[], int m)
{ 
  int m1,i,j,l,ii,k,im,ix[21];
  double h[21],ha,hh,y1,y2,h1,h2,d,hm;
  for (i=0; i<=m; i++) a[i]=0.0;
  if (m>=n) m=n-1;
  if (m>=20) m=19;
  m1=m+1;
  ha=0.0;
  ix[0]=0; ix[m]=n-1;
  l=(n-1)/m; j=l;
  for (i=1; i<=m-1; i++)
    { ix[i]=j; j=j+l;}
  while (1==1)
    { hh=1.0;
    for (i=0; i<=m; i++)
      { a[i]=y[ix[i]]; h[i]=-hh; hh=-hh;}
    for (j=1; j<=m; j++)
      { ii=m1; y2=a[ii-1]; h2=h[ii-1];
      for (i=j; i<=m; i++)
	{ d=x[ix[ii-1]]-x[ix[m1-i-1]];
	y1=a[m-i+j-1];
	h1=h[m-i+j-1];
	a[ii-1]=(y2-y1)/d;
	h[ii-1]=(h2-h1)/d;
	ii=m-i+j; y2=y1; h2=h1;
	}
      }
    hh=-a[m]/h[m];
    for (i=0; i<=m; i++)
      a[i]=a[i]+h[i]*hh;
    for (j=1; j<=m-1; j++)
      { ii=m-j; d=x[ix[ii-1]];
      y2=a[ii-1];
      for (k=m1-j; k<=m; k++)
	{ y1=a[k-1]; a[ii-1]=y2-d*y1;
	y2=y1; ii=k;
	}
      }
    hm=fabs(hh);
    if (hm<=ha) { a[m]=-hm; return;}
    a[m]=hm; ha=hm; im=ix[0]; h1=hh;
    j=0;
    for (i=0; i<=n-1; i++)
      { if (i==ix[j])
	{ if (j<m) j=j+1;}
      else
	{ h2=a[m-1];
	for (k=m-2; k>=0; k--)
	  h2=h2*x[i]+a[k];
	h2=h2-y[i];
	if (fabs(h2)>hm)
	  { hm=fabs(h2); h1=h2; im=i;}
	}
      }
    if (im==ix[0]) return;
    i=0;l=1;
    while (l==1)
      { l=0;
      if (im>=ix[i])
	{ i=i+1;
	if (i<=m) l=1;
	}
      }
    if (i>m) i=m;
    if (i==(i/2)*2) h2=-hh;
    else h2=hh;
    if (h1*h2>=0.0) ix[i]=im;
    else
      { if (im<ix[0])
	{ for (j=m-1; j>=0; j--)
	  ix[j+1]=ix[j];
	ix[0]=im;
	}
      else
	{ if (im>ix[m])
	  { for (j=1; j<=m; j++)
	    ix[j-1]=ix[j];
	  ix[m]=im;
	  }
	else ix[i-1]=im;
	}
      }
    }
}

//---matrix kits
void mat_mult(double A[],  double B[], double C[], int m, int n, int k)
{
  int i, j, l;

  for(i=0; i<m;i++)
    for(j=0;j<k;j++) C[i*m+j]=0.;

  for(i=0; i<m;i++)
    for(j=0;j<k;j++)
      for (l=0;l<n;l++) C[i*m+j]=C[i*m+j] + A[i*m+l]*B[l*n+j];
}

double mat_det(double a[],int n)
{
  int i,j,k,is=0,js=0,l,u,v;
  double f,det,q,d;
  f=1.0; det=1.0;
  for (k=0; k<=n-2; k++)
    { q=0.0;
    for (i=k; i<=n-1; i++)
      for (j=k; j<=n-1; j++)
	{ l=i*n+j; d=fabs(a[l]);
	if (d>q) { q=d; is=i; js=j;}
	}
    if (q+1.0==1.0)
      { det=0.0; return(det);}
    if (is!=k)
      { f=-f;
      for (j=k; j<=n-1; j++)
	{ u=k*n+j; v=is*n+j;
	d=a[u]; a[u]=a[v]; a[v]=d;
	}
      }
    if (js!=k)
      { f=-f;
      for (i=k; i<=n-1; i++)
	{ u=i*n+js; v=i*n+k;
	d=a[u]; a[u]=a[v]; a[v]=d;
	}
      }
    l=k*n+k;
    det=det*a[l];
    for (i=k+1; i<=n-1; i++)
      { d=a[i*n+k]/a[l];
      for (j=k+1; j<=n-1; j++)
	{ u=i*n+j;
	a[u]=a[u]-d*a[k*n+j];
	}
      }
    }
  det=f*det*a[n*n-1];
  return(det);
}

int mat_inv(double a[],int n)
{
  int *is,*js,i,j,k,l,u,v;
  double d,p;
  is=new int[n];
  js=new int[n];
  for (k=0; k<=n-1; k++)
    { d=0.0;
    for (i=k; i<=n-1; i++)
      for (j=k; j<=n-1; j++)
	{ l=i*n+j; p=abs(a[l]);
	if (p>d) { d=p; is[k]=i; js[k]=j;}
	}
    if (d+1.0==1.0)
      { free(is); free(js); printf("err**not inv\n"); exit(1);
      return(0);
      }
    if (is[k]!=k)
      for (j=0; j<=n-1; j++)
	{ u=k*n+j; v=is[k]*n+j;
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
    if (js[k]!=k)
      for (i=0; i<=n-1; i++)
	{ u=i*n+k; v=i*n+js[k];
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
    l=k*n+k;
    a[l]=1.0/a[l];
    for (j=0; j<=n-1; j++)
      if (j!=k)
	{ u=k*n+j; a[u]=a[u]*a[l];}
    for (i=0; i<=n-1; i++)
      if (i!=k)
	for (j=0; j<=n-1; j++)
	  if (j!=k)
	    { u=i*n+j;
	    a[u]=a[u]-a[i*n+k]*a[k*n+j];
	    }
    for (i=0; i<=n-1; i++)
      if (i!=k)
	{ u=i*n+k; a[u]=-a[u]*a[l];}
    }
  for (k=n-1; k>=0; k--)
    { if (js[k]!=k)
      for (j=0; j<=n-1; j++)
	{ u=k*n+j; v=js[k]*n+j;
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
    if (is[k]!=k)
      for (i=0; i<=n-1; i++)
	{ u=i*n+k; v=i*n+is[k];
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
    }
  free(is); free(js);
  return(1);
}

//----eigen tune solver
void mat_change_hessenberg(double a[], int n)
{ 
  int i,j,k,u,v;
  double d,t;
  for (k=1; k<=n-2; k++)
    { d=0.0;
    for (j=k; j<=n-1; j++)
      { u=j*n+k-1; t=a[u];
      if (fabs(t)>fabs(d))
	{ d=t; i=j;}
      }
    if (fabs(d)+1.0!=1.0)
      { if (i!=k)
	{ for (j=k-1; j<=n-1; j++)
	  { u=i*n+j; v=k*n+j;
	  t=a[u]; a[u]=a[v]; a[v]=t;
	  }
	for (j=0; j<=n-1; j++)
	  { u=j*n+i; v=j*n+k;
	  t=a[u]; a[u]=a[v]; a[v]=t;
	  }
	}
      for (i=k+1; i<=n-1; i++)
	{ u=i*n+k-1; t=a[u]/d; a[u]=0.0;
	for (j=k; j<=n-1; j++)
	  { v=i*n+j;
	  a[v]=a[v]-t*a[k*n+j];
	  }
	for (j=0; j<=n-1; j++)
	  { v=j*n+k;
	  a[v]=a[v]+t*a[j*n+i];
	  }
	}
      }
    }
  return;
}

int mat_root_hessenberg(double a[],int n,double u[],double v[],double eps,int jt)
{ 
  int m,it,i,j,k,l,ii,jj,kk,ll;
  double b,c,w,g,xy,p,q,r,x,s,e,f,z,y;
  it=0; m=n;
  while (m!=0)
    { l=m-1;
    while ((l>0)&&(fabs(a[l*n+l-1])>eps*
		   (fabs(a[(l-1)*n+l-1])+fabs(a[l*n+l])))) l=l-1;
    ii=(m-1)*n+m-1; jj=(m-1)*n+m-2;
    kk=(m-2)*n+m-1; ll=(m-2)*n+m-2;
    if (l==m-1)
      { u[m-1]=a[(m-1)*n+m-1]; v[m-1]=0.0;
      m=m-1; it=0;
      }
    else if (l==m-2)
      { b=-(a[ii]+a[ll]);
      c=a[ii]*a[ll]-a[jj]*a[kk];
      w=b*b-4.0*c;
      y=sqrt(fabs(w));
      if (w>0.0)
	{ xy=1.0;
	if (b<0.0) xy=-1.0;
	u[m-1]=(-b-xy*y)/2.0;
	u[m-2]=c/u[m-1];
	v[m-1]=0.0; v[m-2]=0.0;
	}
      else
	{ u[m-1]=-b/2.0; u[m-2]=u[m-1];
	v[m-1]=y/2.0; v[m-2]=-v[m-1];
	}
      m=m-2; it=0;
      }
    else
      { if (it>=jt)
	{ printf("fail\n");
	return(-1);
	}
      it=it+1;
      for (j=l+2; j<=m-1; j++)
	a[j*n+j-2]=0.0;
      for (j=l+3; j<=m-1; j++)
	a[j*n+j-3]=0.0;
      for (k=l; k<=m-2; k++)
	{ if (k!=l)
	  { p=a[k*n+k-1]; q=a[(k+1)*n+k-1];
	  r=0.0;
	  if (k!=m-2) r=a[(k+2)*n+k-1];
	  }
	else
	  { x=a[ii]+a[ll];
	  y=a[ll]*a[ii]-a[kk]*a[jj];
	  ii=l*n+l; jj=l*n+l+1;
	  kk=(l+1)*n+l; ll=(l+1)*n+l+1;
	  p=a[ii]*(a[ii]-x)+a[jj]*a[kk]+y;
	  q=a[kk]*(a[ii]+a[ll]-x);
	  r=a[kk]*a[(l+2)*n+l+1];
	  }
	if ((fabs(p)+fabs(q)+fabs(r))!=0.0)
	  { xy=1.0;
	  if (p<0.0) xy=-1.0;
	  s=xy*sqrt(p*p+q*q+r*r);
	  if (k!=l) a[k*n+k-1]=-s;
	  e=-q/s; f=-r/s; x=-p/s;
	  y=-x-f*r/(p+s);
	  g=e*r/(p+s);
	  z=-x-e*q/(p+s);
	  for (j=k; j<=m-1; j++)
	    { ii=k*n+j; jj=(k+1)*n+j;
	    p=x*a[ii]+e*a[jj];
	    q=e*a[ii]+y*a[jj];
	    r=f*a[ii]+g*a[jj];
	    if (k!=m-2)
	      { kk=(k+2)*n+j;
	      p=p+f*a[kk];
	      q=q+g*a[kk];
	      r=r+z*a[kk]; a[kk]=r;
	      }
	    a[jj]=q; a[ii]=p;
	    }
	  j=k+3;
	  if (j>=m-1) j=m-1;
	  for (i=l; i<=j; i++)
	    { ii=i*n+k; jj=i*n+k+1;
	    p=x*a[ii]+e*a[jj];
	    q=e*a[ii]+y*a[jj];
	    r=f*a[ii]+g*a[jj];
	    if (k!=m-2)
	      { kk=i*n+k+2;
	      p=p+f*a[kk];
	      q=q+g*a[kk];
	      r=r+z*a[kk]; a[kk]=r;
	      }
	    a[jj]=q; a[ii]=p;
	    }
	  }
	}
      }
    }
  return(1);
}

void LinearEquations(double a1, double b1, double c1, double a2, double b2, double c2, double & x, double & y  )
// --the equations to be solved
//     a1*x +b1*y= c1
//     a2*x +b2*y= c2
{
  if(a1*b2-a2*b1 != 0. ){
    x=(c1*b2-c2*b1) / ( a1*b2-a2*b1);
    y=(c1*a2-c2*a1) / ( a2*b1-a1*b2);
  }
}

//---eigen solver from GSL
void  EigenSolver(double Matrix[4][4] , double wr[4], double  wi[4], double vr[4][4], double vi[4][4]) 
{
  int i, j;
  double data[16];

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      data[i*4+j]= Matrix[i][j];
  
  gsl_matrix_view m  = gsl_matrix_view_array (data, 4, 4);
  gsl_vector_complex *eval = gsl_vector_complex_alloc (4);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (4, 4);
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (4);
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
  gsl_eigen_nonsymmv_free (w);
  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  
  for (i = 0; i < 4; i++)
    {
      gsl_complex eval_i = gsl_vector_complex_get (eval, i);
      gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
      wr[i]=  	GSL_REAL(eval_i);
      wi[i]=  	GSL_IMAG(eval_i);
      
      for (j = 0; j < 4; j++)
	{
	  gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);	
	  vr[i][j]=GSL_REAL(z);
	  vi[i][j]=GSL_IMAG(z);
	}
    }
  
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
}

void  EigenSolver_6D(double Matrix[6][6] , double wr[6], double  wi[6], double vr[6][6], double vi[6][6]) 
{
  int i, j;
  double data[36];

  for(i=0;i<6;i++)
    for(j=0;j<6;j++)
      data[i*6+j]= Matrix[i][j];
  
  gsl_matrix_view m  = gsl_matrix_view_array (data, 6, 6);
  gsl_vector_complex *eval = gsl_vector_complex_alloc (6);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (6, 6);
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (6);
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
  gsl_eigen_nonsymmv_free (w);
  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  
  for (i = 0; i < 6; i++)
    {
      gsl_complex eval_i = gsl_vector_complex_get (eval, i);
      gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
      wr[i]=  	GSL_REAL(eval_i);
      wi[i]=  	GSL_IMAG(eval_i);
      
      for (j = 0; j < 6; j++)
	{
	  gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);	
	  vr[i][j]=GSL_REAL(z);
	  vi[i][j]=GSL_IMAG(z);
	}
    }
  
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
}
//---SVD stuffs
void brmul(double a[], double b[] , int m, int n, int k, double c[] )
{ 
  int i,j,l,u;
  for (i=0; i<=m-1; i++)
    for (j=0; j<=k-1; j++)
      { u=i*k+j; c[u]=0.0;
        for (l=0; l<=n-1; l++)
          c[u]=c[u]+a[i*n+l]*b[l*k+j];
      }
  return;
}

static void ppp(double *a,double *e,double *s,double *v,int m,int n)
{ 
  int i,j,p,q;
  double d;
  if (m>=n) i=n;
 else i=m;
 for (j=1; j<=i-1; j++)
   { a[(j-1)*n+j-1]=s[j-1];
   a[(j-1)*n+j]=e[j-1];
      }
    a[(i-1)*n+i-1]=s[i-1];
    if (m<n) a[(i-1)*n+i]=e[i-1];
    for (i=1; i<=n-1; i++)
    for (j=i+1; j<=n; j++)
      { p=(i-1)*n+j-1; q=(j-1)*n+i-1;
        d=v[p]; v[p]=v[q]; v[q]=d;
      }
    return;
  }

static void sss(double fg[2],double cs[2])
{ double r,d;
 if ((fabs(fg[0])+fabs(fg[1]))==0.0)
      { cs[0]=1.0; cs[1]=0.0; d=0.0;}
    else 
      { d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
        if (fabs(fg[0])>fabs(fg[1]))
          { d=fabs(d);
            if (fg[0]<0.0) d=-d;
          }
        if (fabs(fg[1])>=fabs(fg[0]))
          { d=fabs(d);
            if (fg[1]<0.0) d=-d;
          }
        cs[0]=fg[0]/d; cs[1]=fg[1]/d;
      }
    r=1.0;
    if (fabs(fg[0])>fabs(fg[1])) r=cs[1];
    else
      if (cs[0]!=0.0) r=1.0/cs[0];
    fg[0]=d; fg[1]=r;
    return;
  }

int bmuav(double a[],int m,int n,double *u,double *v,double eps,int ka)
{ int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
 double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh,fg[2],cs[2];
 double *s,*e,*w;
 
 s=new  double[ka];
 e=new  double[ka];
 w=new  double[ka];
 
 it=60; k=n;
 if (m-1<n) k=m-1;
 l=m;
 if (n-2<m) l=n-2;
 if (l<0) l=0;
 ll=k;
 if (l>k) ll=l;
 if (ll>=1)
   { for (kk=1; kk<=ll; kk++)
     { if (kk<=k)
       { d=0.0;
       for (i=kk; i<=m; i++)
	 { ix=(i-1)*n+kk-1; d=d+a[ix]*a[ix];}
       s[kk-1]=sqrt(d);
       if (s[kk-1]!=0.0)
	 { ix=(kk-1)*n+kk-1;
	 if (a[ix]!=0.0)
	   { s[kk-1]=fabs(s[kk-1]);
	   if (a[ix]<0.0) s[kk-1]=-s[kk-1];
	   }
	 for (i=kk; i<=m; i++)
	   { iy=(i-1)*n+kk-1;
	   a[iy]=a[iy]/s[kk-1];
	   }
	 a[ix]=1.0+a[ix];
	 }
       s[kk-1]=-s[kk-1];
       }
     if (n>=kk+1)
       { for (j=kk+1; j<=n; j++)
	 { if ((kk<=k)&&(s[kk-1]!=0.0))
	   { d=0.0;
	   for (i=kk; i<=m; i++)
	     { ix=(i-1)*n+kk-1;
	     iy=(i-1)*n+j-1;
	     d=d+a[ix]*a[iy];
	     }
	   d=-d/a[(kk-1)*n+kk-1];
	   for (i=kk; i<=m; i++)
	     { ix=(i-1)*n+j-1;
	     iy=(i-1)*n+kk-1;
	     a[ix]=a[ix]+d*a[iy];
	     }
	   }
                    e[j-1]=a[(kk-1)*n+j-1];
	 }
              }
            if (kk<=k)
              { for (i=kk; i<=m; i++)
                  { ix=(i-1)*m+kk-1; iy=(i-1)*n+kk-1;
                    u[ix]=a[iy];
                  }
              }
            if (kk<=l)
              { d=0.0;
                for (i=kk+1; i<=n; i++)
                  d=d+e[i-1]*e[i-1];
                e[kk-1]=sqrt(d);
                if (e[kk-1]!=0.0)
                  { if (e[kk]!=0.0)
                      { e[kk-1]=fabs(e[kk-1]);
                        if (e[kk]<0.0) e[kk-1]=-e[kk-1];
                      }
                    for (i=kk+1; i<=n; i++)
                      e[i-1]=e[i-1]/e[kk-1];
                    e[kk]=1.0+e[kk];
                  }
                e[kk-1]=-e[kk-1];
                if ((kk+1<=m)&&(e[kk-1]!=0.0))
                  { for (i=kk+1; i<=m; i++) w[i-1]=0.0;
                    for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                        w[i-1]=w[i-1]+e[j-1]*a[(i-1)*n+j-1];
                    for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                        { ix=(i-1)*n+j-1;
                          a[ix]=a[ix]-w[i-1]*e[j-1]/e[kk];
                        }
                  }
                for (i=kk+1; i<=n; i++)
                  v[(i-1)*n+kk-1]=e[i-1];
              }
          }
      }
    mm=n;
    if (m+1<n) mm=m+1;
    if (k<n) s[k]=a[k*n+k];
    if (m<mm) s[mm-1]=0.0;
    if (l+1<mm) e[l]=a[l*n+mm-1];
    e[mm-1]=0.0;
    nn=m;
    if (m>n) nn=n;
    if (nn>=k+1)
      { for (j=k+1; j<=nn; j++)
          { for (i=1; i<=m; i++)
              u[(i-1)*m+j-1]=0.0;
            u[(j-1)*m+j-1]=1.0;
          }
      }
    if (k>=1)
      { for (ll=1; ll<=k; ll++)
          { kk=k-ll+1; iz=(kk-1)*m+kk-1;
            if (s[kk-1]!=0.0)
              { if (nn>=kk+1)
                  for (j=kk+1; j<=nn; j++)
                    { d=0.0;
                      for (i=kk; i<=m; i++)
                        { ix=(i-1)*m+kk-1;
                          iy=(i-1)*m+j-1;
                          d=d+u[ix]*u[iy]/u[iz];
                        }
                      d=-d;
                      for (i=kk; i<=m; i++)
                        { ix=(i-1)*m+j-1;
                          iy=(i-1)*m+kk-1;
                          u[ix]=u[ix]+d*u[iy];
                        }
                    }
                  for (i=kk; i<=m; i++)
                    { ix=(i-1)*m+kk-1; u[ix]=-u[ix];}
                  u[iz]=1.0+u[iz];
                  if (kk-1>=1)
                    for (i=1; i<=kk-1; i++)
                      u[(i-1)*m+kk-1]=0.0;
              }
            else
              { for (i=1; i<=m; i++)
                  u[(i-1)*m+kk-1]=0.0;
                u[(kk-1)*m+kk-1]=1.0;
              }
          }
      }
    for (ll=1; ll<=n; ll++)
      { kk=n-ll+1; iz=kk*n+kk-1;
        if ((kk<=l)&&(e[kk-1]!=0.0))
          { for (j=kk+1; j<=n; j++)
              { d=0.0;
                for (i=kk+1; i<=n; i++)
                  { ix=(i-1)*n+kk-1; iy=(i-1)*n+j-1;
                    d=d+v[ix]*v[iy]/v[iz];
                  }
                d=-d;
                for (i=kk+1; i<=n; i++)
                  { ix=(i-1)*n+j-1; iy=(i-1)*n+kk-1;
                    v[ix]=v[ix]+d*v[iy];
                  }
              }
          }
        for (i=1; i<=n; i++)
          v[(i-1)*n+kk-1]=0.0;
        v[iz-n]=1.0;
      }
    for (i=1; i<=m; i++)
    for (j=1; j<=n; j++)
      a[(i-1)*n+j-1]=0.0;
    m1=mm; it=60;
    while (1==1)
      { if (mm==0)
          { ppp(a,e,s,v,m,n);
            delete s; delete e; delete w; return(1);
          }
        if (it==0)
          { ppp(a,e,s,v,m,n);
            delete s; delete e;  delete w; return(-1);
          }
        kk=mm-1;
	while ((kk!=0)&&(fabs(e[kk-1])!=0.0))
          { d=fabs(s[kk-1])+fabs(s[kk]);
            dd=fabs(e[kk-1]);
            if (dd>eps*d) kk=kk-1;
            else e[kk-1]=0.0;
          }
        if (kk==mm-1)
          { kk=kk+1;
            if (s[kk-1]<0.0)
              { s[kk-1]=-s[kk-1];
                for (i=1; i<=n; i++)
                  { ix=(i-1)*n+kk-1; v[ix]=-v[ix];}
              }
            while ((kk!=m1)&&(s[kk-1]<s[kk]))
              { d=s[kk-1]; s[kk-1]=s[kk]; s[kk]=d;
                if (kk<n)
                  for (i=1; i<=n; i++)
                    { ix=(i-1)*n+kk-1; iy=(i-1)*n+kk;
                      d=v[ix]; v[ix]=v[iy]; v[iy]=d;
                    }
                if (kk<m)
                  for (i=1; i<=m; i++)
                    { ix=(i-1)*m+kk-1; iy=(i-1)*m+kk;
                      d=u[ix]; u[ix]=u[iy]; u[iy]=d;
                    }
                kk=kk+1;
              }
            it=60;
            mm=mm-1;
          }
        else
          { ks=mm;
            while ((ks>kk)&&(fabs(s[ks-1])!=0.0))
              { d=0.0;
                if (ks!=mm) d=d+fabs(e[ks-1]);
                if (ks!=kk+1) d=d+fabs(e[ks-2]);
                dd=fabs(s[ks-1]);
                if (dd>eps*d) ks=ks-1;
                else s[ks-1]=0.0;
              }
            if (ks==kk)
              { kk=kk+1;
                d=fabs(s[mm-1]);
                t=fabs(s[mm-2]);
                if (t>d) d=t;
                t=fabs(e[mm-2]);
                if (t>d) d=t;
                t=fabs(s[kk-1]);
                if (t>d) d=t;
                t=fabs(e[kk-1]);
                if (t>d) d=t;
                sm=s[mm-1]/d; sm1=s[mm-2]/d;
                em1=e[mm-2]/d;
                sk=s[kk-1]/d; ek=e[kk-1]/d;
                b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
                c=sm*em1; c=c*c; shh=0.0;
                if ((b!=0.0)||(c!=0.0))
                  { shh=sqrt(b*b+c);
                    if (b<0.0) shh=-shh;
                    shh=c/(b+shh);
                  }
                fg[0]=(sk+sm)*(sk-sm)-shh;
                fg[1]=sk*ek;
                for (i=kk; i<=mm-1; i++)
                  { sss(fg,cs);
                    if (i!=kk) e[i-2]=fg[0];
                    fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
                    e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
                    fg[1]=cs[1]*s[i];
                    s[i]=cs[0]*s[i];
                    if ((cs[0]!=1.0)||(cs[1]!=0.0))
                      for (j=1; j<=n; j++)
                        { ix=(j-1)*n+i-1;
                          iy=(j-1)*n+i;
                          d=cs[0]*v[ix]+cs[1]*v[iy];
                          v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                          v[ix]=d;
                        }
                    sss(fg,cs);
                    s[i-1]=fg[0];
                    fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
                    s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
                    fg[1]=cs[1]*e[i];
                    e[i]=cs[0]*e[i];
                    if (i<m)
                      if ((cs[0]!=1.0)||(cs[1]!=0.0))
                        for (j=1; j<=m; j++)
                          { ix=(j-1)*m+i-1;
                            iy=(j-1)*m+i;
                            d=cs[0]*u[ix]+cs[1]*u[iy];
                            u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                            u[ix]=d;
                          }
                  }
                e[mm-2]=fg[0];
                it=it-1;
              }
            else
              { if (ks==mm)
                  { kk=kk+1;
                    fg[1]=e[mm-2]; e[mm-2]=0.0;
                    for (ll=kk; ll<=mm-1; ll++)
                      { i=mm+kk-ll-1;
                        fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        if (i!=kk)
                          { fg[1]=-cs[1]*e[i-2];
                            e[i-2]=cs[0]*e[i-2];
                          }
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=n; j++)
                            { ix=(j-1)*n+i-1;
                              iy=(j-1)*n+mm-1;
                              d=cs[0]*v[ix]+cs[1]*v[iy];
                              v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                              v[ix]=d;
                            }
                      }
                  }
                else
                  { kk=ks+1;
                    fg[1]=e[kk-2];
                    e[kk-2]=0.0;
                    for (i=kk; i<=mm; i++)
                      { fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        fg[1]=-cs[1]*e[i-1];
                        e[i-1]=cs[0]*e[i-1];
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=m; j++)
                            { ix=(j-1)*m+i-1;
                              iy=(j-1)*m+kk-2;
                              d=cs[0]*u[ix]+cs[1]*u[iy];
                              u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                              u[ix]=d;
                            }
                      }
                  }
              }
          }
      }
  }

//------FFT, tune finders
void fft(int m, double*x, double*y)
{
   short int dir=1;
   int n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   n = 1;
   for (i=0;i<m;i++) 
      n *= 2;

   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }
}

void  Sin2FFT(int n, double*xr, double*xi, int & nfft, double*famp, double*fphi)
{
  int i, m;
  double temp1,temp2;

  nfft= 1;
  m=0;
  do {
    m=m+1;
    nfft *= 2;} while (nfft <= n);  
  m=m-1;
  nfft=nfft/2;
  
  for (i=0;i<=nfft-1;i++){
    famp[i]=xr[i]*2.0*pow( sin(PI*(i+1)/nfft),2.); 
    fphi[i]=xi[i]*2.0*pow( sin(PI*(i+1)/nfft),2.);
  }
  
  fft(m,famp,fphi);
  
  for (i=0;i<nfft-1;i++){
    temp1=famp[i];
    temp2=fphi[i];
    famp[i]=sqrt(temp1*temp1+temp2*temp2);
    fphi[i]=atan2(temp2,temp1);
  }
}

void  FindMaxPeak(int n, double*xr, double*xi, int & nfft, double*famp, double*fphi, double& peakf, double& peakamp, double& peakphi, int flag)
{
  int i, k;
  double temp1,temp2;
  double peak1;
  double a,b,c;

  peak1=0.0;
  k=0;
  for (i=600;i<=750;i++){
    if(famp[i] > peak1){
      k=i;
      peak1=famp[i];
    }  
  }
  
  a=famp[k];
  if( famp[k+1] >  famp[k-1] ) {
    b=famp[k+1];
  }
  else{
    b=famp[k-1];
  }
  c=cos(2.0*PI/nfft);
  
  temp1=pow( c*(a+b),2.0)-2*a*b*(2*c*c-c-1);
  temp1=( -(a+b*c)*(a-b)+b*sqrt(temp1) )/( a*a+b*b+2.0*a*b*c );
  peakf=1.0*k/nfft+(1.0/2.0/PI)*asin(temp1 *sin(2.0*PI/nfft));
  //if(flag==1) peakf=1.-peakf;
  
  peakamp=0.;
  peakphi=0.;
  for (i=0;i<n;i++){
    peakamp=peakamp+xr[i]*cos(-2.0*PI*i*peakf) -xi[i]*sin(-2.0*PI*i*peakf);
    peakphi=peakphi+xr[i]*sin(-2.0*PI*i*peakf) + xi[i]*cos(-2.0*PI*i*peakf);
  }
  peakamp=peakamp/n;
  peakphi=peakphi/n;
  temp1=sqrt(peakamp*peakamp + peakphi*peakphi)*2.0;
  temp2=atan2(peakphi, peakamp);
  peakamp=temp1;
  peakphi=temp2;
}

void FineTuneFinder(int nfft, double *x, double & fpeak)
{
  int i;
  double xr[nfft], xi[nfft], famp[nfft], fphi[nfft];
  double peakfr0,peakamp, peakphi;
  
  for(i=0;i<nfft;i++){
    xr[i]=x[i];
    xi[i]=0.;
  }
  Sin2FFT(nfft, xr, xi, nfft, famp, fphi);
  FindMaxPeak(1024, xr, xi, nfft, famp, fphi, peakfr0, peakamp, peakphi,1);
  fpeak=peakfr0;  
}

//==============================================================
//
//         linear  tpsa: x[7]  and  linear map x[6][7]
//
//===============================================================
class tps
{
 public:
  tps() 
    {
      int i;
      for (i=0;i<7;i++) sample[i]=0.;  
    }
  tps(double d )
    {
      int i;
      for (i=0;i<7;i++) sample[i]=0;
      sample[0]=d;  
    }
  tps(double d[7] )
    {
      int i;
      for (i=0;i<7;i++) sample[i]=d[i];  
    }
  double &operator[](int i)
    { 
      return sample[i]; 
    }
  friend ostream &operator<<(ostream &stream, tps x)
    {
      int i;
      for(i=0;i<7;i++) stream<<setw(12)<<setprecision(9)<<scientific<<x[i]<<"  ";
      return stream;
    }
  friend tps operator+(tps x, tps y)
    {
      int i;
      tps z;

      for (i = 0; i < 7; i++)  z[i] = x[i] + y[i];
      return z;
    }
  friend tps operator-(tps x, tps y)
    {
      int i;
      tps z;
      
      for (i = 0; i < 7; i++)  z[i] = x[i] - y[i];
      return z;
    }
  friend tps operator*(tps x, tps y)
    {
      int i;
      tps z;
      
      z[0] = x[0] * y[0];
      for (i = 1; i < 7; i++)
	z[i] = x[0] * y[i] + x[i] * y[0];
      return z;
    }
  friend tps DAinv(tps x)
    {
      int i;
      double a, temp;
      tps z;      

      z[0] = 1.0 / x[0];
      temp = x[0];
      a = -1.0 / (temp * temp);
      for (i = 1; i < 7; i++) z[i] = a * x[i];
      return z;
    }

  friend tps operator/(tps x, tps y)
    {
      return x*DAinv(y);
    }
  friend tps sqr(tps x)
    {
      return x*x;
    }
  friend tps sqrt(tps x)
    {
      int i;
      double a;
      tps  z;
      
      a = sqrt(x[0]);
      z[0] = a;
      a = 0.5 / a;
      for (i = 1; i < 7; i++) z[i] = a * x[i];
      return z;
    }
  friend tps sin(tps x)
    {
      int i;
      double a;
      tps z;

      z[0] = sin(x[0]); 
      a =    cos(x[0]);
      for (i = 1; i < 7; i++) z[i] = a * x[i];
      return z;
    }
  friend tps cos(tps x)
    {
      int i;
      double a;
      tps z;
      
      z[0] = cos(x[0]); 
      a =   -sin(x[0]);
      for (i = 1; i < 7; i++) z[i] = a * x[i];
      return z;
    }
  friend tps tan(tps x)
    {
      return sin(x)/cos(x);
    }
  friend tps atan(tps x)
    {
      int  i;
      double a;
      tps z;
      
      a = x[0];
      z[0] = atan(a);
      a = 1 / (1 + a * a);
      for (i = 1; i < 7; i++)  z[i] = a * x[i];
      return z;
    }
  friend tps sinh(tps x)
    {
      int i;
      double a;
      tps z;
      
      z[0] = sinh(x[0]); 
      a = cosh(x[0]);
      for (i = 1; i < 7; i++) z[i] = a * x[i];
      return z;
    }
  friend tps cosh(tps x)
    {
      int i;
      double a;
      tps z;
      
      z[0] = cosh(x[0]); 
      a = sinh(x[0]);
      for (i = 1; i <7; i++) z[i] = a * x[i];
      return z;
    }
  friend tps exp(tps x)
    {
      int i;
      double a;
      tps z;
      
      a = exp(x[0]);
      z[0] = a;
      for (i = 1; i <7; i++) z[i] = a * x[i];
      return z;
    }
  friend  tps ln(tps x)
    {
      int i;
      tps z;
      
      z[0] = log(x[0]);
      for (i = 1; i < 7; i++) z[i] = x[i] / x[0];
      return z;
    } 
 private:
  double  sample[7];
};

class linmap
{
 public:
  linmap()
    {
      int i,j;
      for(i=0;i<6;i++)
	for(j=0;j<7;j++) map0[i][j]=0.;
    }
  void identity()
    {
      int i,j;
      for(i=0;i<6;i++)
	for(j=0;j<7;j++) map0[i][j]=0.;
      map0[0][1]=1;
      map0[1][2]=1;
      map0[2][3]=1;
      map0[3][4]=1;
      map0[4][5]=1;
      map0[5][6]=1;
    }
  linmap(double x[6])
    {
      int i,j;
      for(i=0;i<6;i++)
	for(j=0;j<7;j++) map0[i][j]=0.;
      for(i=0;i<6;i++) map0[i][0]=x[i];
    }
  tps &operator[](int i) 
    { 
      return map0[i]; 
    }
  
  friend linmap operator+(linmap x, linmap y)
  {
    int i,j;
    linmap z;
    for(i=0;i<6;i++)
      for(j=0;j<7;j++) z[i][j]=x[i][j]+y[i][j];
    return z;
  }
  void print()
  {
    int i;
    for(i=0;i<6;i++)  cout<<map0[i]<<endl;  
  }
 private:
  tps map0[6];
};

void Getmat(linmap map0, double x[36])
{
  int i,j;
  for (i=0;i<6;i++)
    for(j=0;j<6;j++) x[i*6+j]=map0[i][j+1];
}

void Getpos(linmap map0, double x[6])
{
  int i;
  for (i=0;i<6;i++)  x[i]=map0[i][0];
}

//=========================================
//
//            Magnete Alignment
//
//=========================================
template<typename T>
void  GtoL (T x[6], double DX, double DY, double DT)
{

  int i;
  T xtemp[6];
  
  if ( abs(DX) +  abs(DY)  > 1.0e-10 ) {
    x[0] = x[0] - DX;   x[2] = x[2] - DY;
  }
  if (  abs(DT)  > 1.0e-10 ) {
    double cosT=cos(DT), sinT=sin(DT);
    for (i=0; i<6;i++) xtemp[i]=x[i];
    x[0] = cosT * xtemp[0] + sinT * xtemp[2];
    x[1] = cosT * xtemp[1] + sinT * xtemp[3];
    x[2] = cosT * xtemp[2] - sinT * xtemp[0];
    x[3] = cosT * xtemp[3] - sinT * xtemp[1];
  }
}

template<typename T>
void LtoG(T x[6], double DX, double DY, double DT)
{
  int i;
  T xtemp[6];

  if (  abs(DT)  > 1.0e-10 ) {
    double cosT=cos(DT), sinT=sin(DT);
    for(i=0; i<6;i++) xtemp[i]=x[i];
    x[0] = cosT * xtemp[0] - sinT * xtemp[2];
    x[1] = cosT * xtemp[1] - sinT * xtemp[3];
    x[2] = sinT * xtemp[0] + cosT * xtemp[2];
    x[3] = sinT * xtemp[1] + cosT * xtemp[3];
  }
  if ( abs(DX) +  abs(DY)  > 1.0e-10 ) {
    x[0] = x[0] + DX  ;  x[2] = x[2]+ DY;
  }
}

//=================================================================
//
//    4th order Symplectic Integrator  (S.I.)
//
//=================================================================
template <class T>
void DRIFT_Pass(T x[6], double L)
{
  T u;
  if( L !=0. ) {
    u=L/(1+x[delta_]);
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
  }
}

template <class T>
void bend_kick_pass(T x[6], double L, double href  )
{
  x[px_]=x[px_]+(href*x[delta_]-href*href*x[x_])*L;
  x[z_]=x[z_]-href*x[x_]*L;
}

template <class T>
void quad_kick_pass(T x[6], double k1l, double k1sl)
{
  int i, Norder;
  T  Xn,  Yn, Xn0, Yn0;
  T  By, Bx;

  By=0.;  Bx=0.;  Xn=1.;  Yn=0.;
  Norder=1;
  for(i=1;i<Norder+1;i++){
    Xn0=Xn;
    Yn0=Yn;
    Xn=Xn0*x[x_]-Yn0*x[y_];
    Yn=Xn0*x[y_]+Yn0*x[x_];
  }
  By=By+(k1l*Xn-k1sl*Yn);
  Bx=Bx+(k1l*Yn+k1sl*Xn);
  x[px_]=x[px_]-By;
  x[py_]=x[py_]+Bx;
}

template <class T>
void sext_kick_pass(T x[6], double k2l, double k2sl)
{
  int i, Norder;
  T  Xn,  Yn, Xn0, Yn0;
  T  By, Bx;

  By=0.;  Bx=0.;  Xn=1.;  Yn=0.;
  Norder=2;
  for(i=1;i<Norder+1;i++){
    Xn0=Xn;
    Yn0=Yn;
    Xn=Xn0*x[x_]-Yn0*x[y_];
    Yn=Xn0*x[y_]+Yn0*x[x_];
  }
  By=By+(k2l*Xn-k2sl*Yn)/2;
  Bx=Bx+(k2l*Yn+k2sl*Xn)/2;
  
  x[px_]=x[px_]-By;
  x[py_]=x[py_]+Bx;
}

template <class T>
void oct_kick_pass(T x[6], double k3l, double k3sl)
{
  int i, Norder;
  T  Xn,  Yn, Xn0, Yn0;
  T  By, Bx;

  By=0.;   Bx=0.;  Xn=1.;  Yn=0.;
  Norder=3;
  for(i=1;i<Norder+1;i++){
    Xn0=Xn;
    Yn0=Yn;
    Xn=Xn0*x[x_]-Yn0*x[y_];
    Yn=Xn0*x[y_]+Yn0*x[x_];
  }
  By=By+(k3l*Xn-k3sl*Yn)/6;
  Bx=Bx+(k3l*Yn+k3sl*Xn)/6;
  
  x[px_]=x[px_]-By;
  x[py_]=x[py_]+Bx;
}

template <class T>
void mult_kick_pass(T x[6], int Norder, double KNL[11], double KNSL[11])
{
  int i;
  int fac=1;
  T   Xn,  Yn, Xn0, Yn0;
  T   By, Bx;

  By=KNL[0];
  Bx=KNSL[0];
  Xn=1.;
  Yn=0.;

  for(i=1;i<Norder+1;i++){
    Xn0=Xn;
    Yn0=Yn;
    Xn=Xn0*x[x_]-Yn0*x[y_];
    Yn=Xn0*x[y_]+Yn0*x[x_];
    fac=fac*i;
    if ( KNL[i] != 0. || KNSL[i] !=0. ) {
      By=By+(KNL[i]*Xn-KNSL[i]*Yn)/fac;
      Bx=Bx+(KNL[i]*Yn+KNSL[i]*Xn)/fac;
    }
  }
  x[px_]=x[px_]-By;
  x[py_]=x[py_]+Bx;
}

template <class T>
void bend_mult_kick_pass(T x[6], double L, double href, int Norder, double KNL[11], double KNSL[11])
{
  int i;
  int fac=1;
  T   Xn,  Yn, Xn0, Yn0;
  T   By, Bx;

  By=KNL[0];
  Bx=KNSL[0];
  Xn=1.;
  Yn=0.;

  for(i=1;i<Norder+1;i++){
    Xn0=Xn;
    Yn0=Yn;
    Xn=Xn0*x[x_]-Yn0*x[y_];
    Yn=Xn0*x[y_]+Yn0*x[x_];
    fac=fac*i;
    if ( KNL[i] != 0. || KNSL[i] !=0. ) {
      By=By+(KNL[i]*Xn-KNSL[i]*Yn)/fac;
      Bx=Bx+(KNL[i]*Yn+KNSL[i]*Xn)/fac;
    }
  }
  x[px_]=x[px_]-By  + (href*x[delta_]-href*href*x[x_])*L;;
  x[py_]=x[py_]+Bx;
  x[z_]=x[z_]-href*x[x_]*L;
}

template <class T> 
void SBEND_Pass(T x[6], double L, int Nint, double Angle, double E1, double E2) 
{
  int i;
  double href=Angle/L;
  double Lint=L/Nint;
  T  u;
  
  x[1] = x[1]+ tan(E1)*x[0]*href;   
  x[3] = x[3]- tan(E1)*x[2]*href;

  for(i=0;i<Nint;i++){
    //DRIFT_Pass(x,Fdrift1*Lint);
    L=Fdrift1*Lint;
    u=L/(1+x[delta_]);
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
    
    //bend_kick_pass(x, Fkick1*Lint, href);
    L=Fkick1*Lint;
    x[px_]=x[px_]+(href*x[delta_]-href*href*x[x_])*L;
    x[z_]=x[z_]-href*x[x_]*L;
    
    //DRIFT_Pass(x,Fdrift2*Lint);
    L=Fdrift2*Lint;
    u=L/(1+x[delta_]);
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
    
    //bend_kick_pass(x, Fkick2*Lint, href);
    L=Fkick2*Lint;
    x[px_]=x[px_]+(href*x[delta_]-href*href*x[x_])*L;
    x[z_]=x[z_]-href*x[x_]*L;
    
    //DRIFT_Pass(x,Fdrift2*Lint);
    L=Fdrift2*Lint;
    u=L/(1+x[delta_]);
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
    
    //bend_kick_pass(x, Fkick1*Lint, href);
    L=Fkick1*Lint;
    x[px_]=x[px_]+(href*x[delta_]-href*href*x[x_])*L;
    x[z_]=x[z_]-href*x[x_]*L;
    
    //DRIFT_Pass(x,Fdrift1*Lint); 
    L=Fdrift1*Lint;
    u=L/(1+x[delta_]);
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
  }

  x[1] = x[1]+ tan(E2)*x[0]*href;   
  x[3] = x[3]- tan(E2)*x[2]*href;
}

template <class T> 
void QUAD_Pass(T x[6], double L, int Nint, double k1l, double k1sl)
{
  int i;
  double Lint=L/Nint;
  
  if(L==0.) 
    {
      quad_kick_pass(x, k1l, k1sl);
    }
  else 
    {
      double k1l_kick1,k1sl_kick1;
      double k1l_kick2,k1sl_kick2;
      k1l_kick1 =Fkick1*k1l/Nint;
      k1sl_kick1=Fkick1*k1sl/Nint;
      k1l_kick2 =Fkick2*k1l/Nint;
      k1sl_kick2=Fkick2*k1sl/Nint;

      int j, Norder=1;
      T  Xn,  Yn, Xn0, Yn0;
      T  By, Bx;
      T  u;     
      for(i=0;i<Nint;i++){
	
	//DRIFT_Pass(x,Fdrift1*Lint);
	L=Fdrift1*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
	
	//quad_kick_pass(x, k1l_kick1, k1sl_kick1);
        k1l=k1l_kick1;
        k1sl=k1sl_kick1;
	By=0.;  Bx=0.;  Xn=1.;  Yn=0.;
	for(j=1;j<Norder+1;j++){
	  Xn0=Xn;
	  Yn0=Yn;
	  Xn=Xn0*x[x_]-Yn0*x[y_];
	  Yn=Xn0*x[y_]+Yn0*x[x_];
	}
	By=By+(k1l*Xn-k1sl*Yn);
	Bx=Bx+(k1l*Yn+k1sl*Xn);
	x[px_]=x[px_]-By;
	x[py_]=x[py_]+Bx;
	
	//DRIFT_Pass(x,Fdrift2*Lint);
	L=Fdrift2*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
	
	//quad_kick_pass(x, k1l_kick2, k1sl_kick2);
        k1l=k1l_kick2;
        k1sl=k1sl_kick2;
	By=0.;  Bx=0.;  Xn=1.;  Yn=0.;
	for(j=1;j<Norder+1;j++){
	  Xn0=Xn;
	  Yn0=Yn;
	  Xn=Xn0*x[x_]-Yn0*x[y_];
	  Yn=Xn0*x[y_]+Yn0*x[x_];
	}
	By=By+(k1l*Xn-k1sl*Yn);
	Bx=Bx+(k1l*Yn+k1sl*Xn);
	x[px_]=x[px_]-By;
	x[py_]=x[py_]+Bx;
	
	//DRIFT_Pass(x,Fdrift2*Lint);
	L=Fdrift2*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
	
	//quad_kick_pass(x, k1l_kick1, k1sl_kick1);
        k1l=k1l_kick1;
        k1sl=k1sl_kick1;
	By=0.;  Bx=0.;  Xn=1.;  Yn=0.;
	for(j=1;j<Norder+1;j++){
	  Xn0=Xn;
	  Yn0=Yn;
	  Xn=Xn0*x[x_]-Yn0*x[y_];
	  Yn=Xn0*x[y_]+Yn0*x[x_];
	}
	By=By+(k1l*Xn-k1sl*Yn);
	Bx=Bx+(k1l*Yn+k1sl*Xn);
	x[px_]=x[px_]-By;
	x[py_]=x[py_]+Bx;
	
	//DRIFT_Pass(x,Fdrift1*Lint);
	L=Fdrift1*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
      }
    }
} 

template <class T> 
void SEXT_Pass(T x[6], double L, int Nint, double k2l, double k2sl)
{
  int i;
  double Lint=L/Nint;

  if(L==0.) 
    {
      sext_kick_pass(x, k2l, k2sl);
    }
  else 
    {
      double k2l_kick1,k2sl_kick1;
      double k2l_kick2,k2sl_kick2;
      k2l_kick1 =Fkick1*k2l/Nint;
      k2sl_kick1=Fkick1*k2sl/Nint;
      k2l_kick2 =Fkick2*k2l/Nint;
      k2sl_kick2=Fkick2*k2sl/Nint;
      
      int j, Norder=2;
      T  Xn,  Yn, Xn0, Yn0;
      T  By, Bx;
      T  u;  
      for(i=0;i<Nint;i++){
	
	//DRIFT_Pass(x,Fdrift1*Lint);
	L=Fdrift1*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
	
	//sext_kick_pass(x, k2l_kick1, k2sl_kick1);
	k2l= k2l_kick1; k2sl= k2sl_kick1;
	By=0.;  Bx=0.;  Xn=1.;  Yn=0.;
	for(j=1;j<Norder+1;j++){
	  Xn0=Xn;
	  Yn0=Yn;
	  Xn=Xn0*x[x_]-Yn0*x[y_];
	  Yn=Xn0*x[y_]+Yn0*x[x_];
	}
	By=By+(k2l*Xn-k2sl*Yn)/2;
	Bx=Bx+(k2l*Yn+k2sl*Xn)/2;
	x[px_]=x[px_]-By;   
	x[py_]=x[py_]+Bx;
	
	//DRIFT_Pass(x,Fdrift2*Lint);
	L=Fdrift2*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
	
	//sext_kick_pass(x, k2l_kick2, k2sl_kick2);
	k2l= k2l_kick2;   k2sl= k2sl_kick2;
	By=0.;  Bx=0.;  Xn=1.;  Yn=0.;
	for(j=1;j<Norder+1;j++){
	  Xn0=Xn;
	  Yn0=Yn;
	  Xn=Xn0*x[x_]-Yn0*x[y_];
	  Yn=Xn0*x[y_]+Yn0*x[x_];
	}
	By=By+(k2l*Xn-k2sl*Yn)/2;
	Bx=Bx+(k2l*Yn+k2sl*Xn)/2;
	x[px_]=x[px_]-By;   
	x[py_]=x[py_]+Bx;
	
	//DRIFT_Pass(x,Fdrift2*Lint);
	L=Fdrift2*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
	
	//sext_kick_pass(x, k2l_kick1, k2sl_kick1);
	k2l= k2l_kick1; k2sl= k2sl_kick1;
	By=0.;  Bx=0.;  Xn=1.;  Yn=0.;
	for(j=1;j<Norder+1;j++){
	  Xn0=Xn;
	  Yn0=Yn;
	  Xn=Xn0*x[x_]-Yn0*x[y_];
	  Yn=Xn0*x[y_]+Yn0*x[x_];
	}
	By=By+(k2l*Xn-k2sl*Yn)/2;
	Bx=Bx+(k2l*Yn+k2sl*Xn)/2;
	x[px_]=x[px_]-By;   
	x[py_]=x[py_]+Bx;
	
	//DRIFT_Pass(x,Fdrift1*Lint);
	L=Fdrift1*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
      }
    }
} 

template <class T> 
void OCT_Pass(T x[6], double L, int Nint, double k3l, double k3sl)
{
  int i;
  double Lint=L/Nint;
  
  if(L==0.) 
    {
      oct_kick_pass(x, k3l, k3sl);
    }
  else 
    {
      double k3l_kick1,k3sl_kick1;
      double k3l_kick2,k3sl_kick2;
      k3l_kick1 =Fkick1*k3l/Nint;
      k3sl_kick1=Fkick1*k3sl/Nint;
      k3l_kick2 =Fkick2*k3l/Nint;
      k3sl_kick2=Fkick2*k3sl/Nint;
      for(i=0;i<Nint;i++){
	DRIFT_Pass(x,Fdrift1*Lint);
	oct_kick_pass(x, k3l_kick1, k3sl_kick1);
	DRIFT_Pass(x,Fdrift2*Lint);
	oct_kick_pass(x, k3l_kick2, k3sl_kick2);
	DRIFT_Pass(x,Fdrift2*Lint);
	oct_kick_pass(x, k3l_kick1, k3sl_kick1);
	DRIFT_Pass(x,Fdrift1*Lint);
      }
    }
} 

template <class T> 
void MULT_Pass(T x[6], double L, int Nint, int Norder, double KNL[11], double KNSL[11]) 
{
  int i;
  double Lint=L/Nint;
  
  if(L==0.) 
    {
      mult_kick_pass(x, Norder, KNL, KNSL);
    }
  else 
    {
      double knl_kick1[11],knsl_kick1[11];
      double knl_kick2[11],knsl_kick2[11];
      for(i=0;i<11;i++) {
	knl_kick1[i] =Fkick1*KNL[i]/Nint;
	knsl_kick1[i]=Fkick1*KNSL[i]/Nint;
      }
      for(i=0;i<11;i++) {
	knl_kick2[i] =Fkick2*KNL[i]/Nint;
	knsl_kick2[i]=Fkick2*KNSL[i]/Nint;
      }
      
      for(i=0;i<Nint;i++){
	DRIFT_Pass(x,Fdrift1*Lint);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift2*Lint);
	mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
	DRIFT_Pass(x,Fdrift2*Lint);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift1*Lint);
      }
    }
} 

template <class T> 
void GMULT_Pass(T x[6], double L, int Nint,  double Angle, int Norder, double KNL[11], double KNSL[11], double E1, double E2) 
{
   if(L==0.) 
    {
      mult_kick_pass(x, Norder, KNL, KNSL);
    }
   else 
     {
       int i;
       double href=Angle/L;
       double Lint=L/Nint;
       double knl_kick1[11],knsl_kick1[11];
       double knl_kick2[11],knsl_kick2[11];

       for(i=0;i<11;i++) {
	 knl_kick1[i] =Fkick1*KNL[i]/Nint;
	 knsl_kick1[i]=Fkick1*KNSL[i]/Nint;
       }
       for(i=0;i<11;i++) {
	 knl_kick2[i] =Fkick2*KNL[i]/Nint;
	 knsl_kick2[i]=Fkick2*KNSL[i]/Nint;
       }
       
       x[1] = x[1]+ tan(E1)*x[0]*href;   
       x[3] = x[3]- tan(E1)*x[2]*href;       
       for(i=0;i<Nint;i++){
	 DRIFT_Pass(x,Fdrift1*Lint);
	 bend_mult_kick_pass(x, Fdrift1*Lint, href, Norder,knl_kick1, knsl_kick1);
	 DRIFT_Pass(x,Fdrift2*Lint);
	 bend_mult_kick_pass(x, Fdrift2*Lint, href, Norder,knl_kick2, knsl_kick2);
	 DRIFT_Pass(x,Fdrift2*Lint);
	 bend_mult_kick_pass(x, Fdrift1*Lint, href,  Norder,knl_kick1, knsl_kick1);
	 DRIFT_Pass(x,Fdrift1*Lint);
       }
       x[1] = x[1]+ tan(E2)*x[0]*href;   
       x[3] = x[3]- tan(E2)*x[2]*href; 
     }
} 

//===========================================
//
//       Element classes
//
//============================================

class Element
{
 public:
  Element(string name)
    {
      int i;
      NAME=name;
      DX=0.;      // offset of magnet center w.r.t x-y frame, DX>0 means magnet moved to positive x
      DY=0.;      // offset of magnet center w.r.t x-y frame, DY>0 means magnet moved to positive y
      DT=0.;      // tilt angle of magnet w.r.t x-y frame, roll from e_x-->e_y, DT > 0 ;  Normal Q rolled -PI/4 to get a skewQ, k1s=k1
      for(i=0;i<6;i++)   X[i]=0.;
      for(i=0;i<36;i++)  T[i]=0.;
      for(i=0;i<36;i++)  M[i]=0.;
      for(i=0;i<36;i++)  A[i]=0.;
      Beta1=0.;
      Beta2=0.;
      Beta3=0.;
      Alfa1=0.;
      Alfa2=0.;
      Alfa3=0.;
      r=0;
      c11=0.;
      c12=0.;
      c21=0.;
      c22=0.;
      Etax=0.;
      Etay=0.;
      Etaxp=0.;
      Etayp=0.;
      Mu1=0.;
      Mu2=0.;
      Mu3=0.;
      APx=1.;
      APy=1.;
    }
  virtual void    SetP(const char *name, double value)=0;
  virtual double  GetP(const char *name)=0;
  virtual void    Pass(double x[6])=0;
  virtual void    DAPass(tps x[6])=0; 
  string  NAME, TYPE;
  double  S, L, DX, DY, DT;
  double  X[6], T[36], M[36], A[36];
  double  Beta1, Alfa1, Beta2, Alfa2,  Beta3, Alfa3, Mu1, Mu2, Mu3;
  double  r, c11, c12, c21, c22;
  double  Etax, Etay, Etaxp,  Etayp;
  double  APx, APy;
};

//---------------DRIFT-----------------------------------
class DRIFT: public Element
{
 public:
  DRIFT(string name, double l): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("DRIFT");
	  L=l;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void  SetP(const char *name, double value) 
    {
      cout<<"No parameter to be set for DRIFT."<<endl;
      exit(0);
    }
  double GetP(const char *name) 
    {
      cout<<"No parameter to be returned for DRIFT."<<endl;
      exit(0);
    }
  void   Pass(double x[6]){
    DRIFT_Pass(x, L);
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x, L); 
  }
};

//--------------------------SBEND------------------------------------
class SBEND: public Element
{
 public:
  SBEND(string name, double l, double angle, double e1, double e2): Element(name)
    {
      if (l >0 )
	{ 
	  TYPE=string("SBEND");
	  L=l;
	  ANGLE=angle;
          E1=e1;
          E2=e2;
	  Nint=int(L/BLslice)+1;
	}
      else
	{
	  cout<<"Error: SBEND length can not be zero."<<endl;
          exit(1); 
	}
    }
  void SetP(const char *name, double value)
    {
      if (strcmp( name, "ANGLE" ) == 0) 
	{
	  ANGLE=value;
	}
      else if (strcmp( name, "E1" ) == 0) 
	{
	  E1=value;
	}
      else if (strcmp( name, "E2" ) == 0) 
	{
	  E2=value;
	}
      else if (strcmp( name, "Nint" ) == 0) 
	{
	  Nint=int(value);
	}
      else 
	{
	  cout<<"SBEND does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 
  double GetP(const char *name)
    {
      if (strcmp( name, "ANGLE" ) == 0) 
	{
	  return ANGLE;
	}
      else if (strcmp( name, "E1" ) == 0) 
	{
	  return E1;
	}
      else if (strcmp( name, "E2" ) == 0) 
	{
	  return E2;
	}
      else if (strcmp( name, "Nint" ) == 0) 
	{
	  return Nint;
	}
      else 
	{
	  cout<<"SBEND does not have parameter of "<<name<<endl; 
          exit(0); 
	} 
    }  
  void   Pass(double x[6]){
    GtoL(x, DX,DY,DT);  SBEND_Pass(x,L,Nint,ANGLE,E1,E2);  LtoG(x,DX,DY,DT); 
  }
  void   DAPass(tps x[6]){
    GtoL(x,DX,DY,DT);   SBEND_Pass(x,L,Nint,ANGLE,E1,E2);  LtoG(x,DX,DY,DT);  
  }
 private:
  double  ANGLE, E1,E2;
  int Nint;
};

//---------------------QUAD-----------------------------------------
class QUAD: public Element
{
 public:
  QUAD(string name, double l, double k1l ): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("QUAD");
	  L=l;
	  K1L=k1l;
          Nint=int(L/QLslice)+1;
          Norder=1;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void SetP(const char *name, double value)
    { 
      if (strcmp( name, "K1L" ) == 0) 
	{
	  K1L=value;
	}
      else if (strcmp( name, "Nint" ) == 0) 
	{
	  Nint=int(value);
	}
      else 
	{
	  cout<<"QUAD does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 
  double GetP(const char *name)
    {
      if (strcmp( name, "K1L" ) == 0) 
	{
	  return K1L;
	}
      else if (strcmp( name, "Nint" ) == 0) 
	{
	  return Nint;
	}
      else 
	{
	  cout<<"QUAD does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }  
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);    QUAD_Pass(x,L,Nint,K1L, 0.);    LtoG(x,DX,DY,DT);
  }
  void   DAPass(tps x[6]){
    GtoL(x,DX,DY,DT);    QUAD_Pass(x,L,Nint,K1L, 0.);    LtoG(x,DX,DY,DT);
  }
 private:
  double  K1L;
  int Nint, Norder;
};

//--------------------------------SKEWQ------------------------------
class SKEWQ: public Element
{
 public:
  SKEWQ(string name, double l, double k1sl): Element(name)
    { 
      if (l >=0 ) { 
	TYPE=string("SKEWQ");  
	L=l; 
	K1SL= k1sl; 
	Nint=int(L/QLslice)+1;
	Norder=1;}
      else { 
	cout<<"Error: check signs of L and S"<<endl; 
	exit(1); 
      }
    }
  void SetP(const char *name, double value)     
    {
      if (strcmp( name, "K1SL" ) == 0) 	
	{  
	  K1SL=value;	
	}
      else if (strcmp( name, "Nint" ) == 0) 
	{
	  Nint=int(value);
	}
      else  { 
	cout<<"SKEWQ does not have  parameter of "<<name<<endl; 
	exit(0);
      } 
    } 
  double GetP(const char *name)
    {
      if (strcmp( name, "K1SL" ) == 0) 	
	{
	  return K1SL;	
	}
      else if (strcmp( name, "Nint" ) == 0) 
	{
	  return Nint;
	}
      else  
	{
	  cout<<"SKEWQ does not have  parameter of "<<name<<endl; 
	  exit(0); 
	}
    }  
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);    QUAD_Pass(x,L,Nint, 0., K1SL);    LtoG(x,DX,DY,DT);
  }
  void   DAPass(tps x[6]){
    GtoL(x,DX,DY,DT);    QUAD_Pass(x,L,Nint, 0., K1SL);    LtoG(x,DX,DY,DT);
  }
 private:
  double  K1SL;
  int Nint, Norder;
};

//----------------------------SEXT-------------------------------------------
class SEXT: public Element
{
 public:
  SEXT(string name, double l, double k2l): Element(name)
    { 
      if (l >=0)
	{
	  TYPE=string("SEXT");
	  L=l;
          K2L=k2l;
	  Nint=int(L/QLslice)+1;
	  Norder=2;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void SetP(const char *name, double value)
    { 
      if (strcmp( name, "K2L" ) == 0) 
	{
	  K2L=value;
	}
      else if (strcmp( name, "Nint" ) == 0) 
	{
	  Nint=int(value);
	}
      else 
	{
	  cout<<"SEXT does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 
  double GetP(const char *name)
    {
      if (strcmp( name, "K2L" ) == 0) 
	{
	  return K2L;
	}
      else if (strcmp( name, "Nint" ) == 0) 
	{
	  return Nint;
	}
      else 
	{
	  cout<<"SEXT does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }  
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);    SEXT_Pass(x,L,Nint,K2L, 0.);      LtoG(x, DX,DY,DT);
  }
  void   DAPass(tps x[6]){
    GtoL(x,DX,DY,DT);    SEXT_Pass(x,L,Nint,K2L, 0.);      LtoG(x, DX,DY,DT);
  }
 private:
  double  K2L;
  int Nint, Norder;
};

//-------------------------OCT---------------------------------------------
class OCT: public Element
{
 public:
  OCT(string name, double l, double k3l): Element(name)
    { 
      if (l >=0  )
	{
	  TYPE=string("OCT");
	  L=l;
          K3L=k3l;
	  Nint=int(L/QLslice)+1;
	  Norder=3;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void SetP(const char *name, double value)
    { 
     if (strcmp( name, "K3L" ) == 0) 
	{
	  K3L=value;
	}
      else if (strcmp( name, "Nint" ) == 0) 
	{
	  Nint=int(value);
	}
      else 
	{
	  cout<<"OCT does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 
  double GetP(const char *name)
    {
      if (strcmp( name, "K3L" ) == 0) 
	{
	  return K3L;
	}
      else if (strcmp( name, "Nint" ) == 0) 
	{
	  return Nint;
	}
      else 
	{
	  cout<<"OCT does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }  
  void   Pass(double x[6]){
      GtoL(x,DX,DY,DT);    OCT_Pass(x,L,Nint,K3L, 0.);      LtoG(x, DX,DY,DT);
  }
  void   DAPass(tps x[6]){
      GtoL(x,DX,DY,DT);    OCT_Pass(x,L,Nint,K3L, 0.);      LtoG(x, DX,DY,DT);
  }
 private:
  double  K3L;
  int Nint,Norder;
};

//-------------------------------MULT--------------------------------------
class MULT: public Element
{
 public:
  MULT(string name, double l, double knl[11], double knsl[11]): Element(name)
    { 
      int i;
      if (l >=0  )
	{
	  TYPE=string("MULT");
	  L=l;
          for(i=0;i<11;i++) {
	    KNL[i]=knl[i];  KNSL[i]=knsl[i]; }
	  Nint=int(L/QLslice)+1;
	  Norder=1;
	  for(i=0;i<10;i++) {
	    if( KNL[10-i] != 0. ||  KNSL[10-i] != 0. ){
	      Norder=10-i;
	      break;
	    }
	  }
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
   void SetP(const char *name, double value)
     { 
       if (     strcmp( name, "K0L" ) == 0) 	{  KNL[0]=value; }
       else if (strcmp( name, "K1L" ) == 0) 	{  KNL[1]=value; }      
       else if (strcmp( name, "K2L" ) == 0) 	{  KNL[2]=value; }
       else if (strcmp( name, "K3L" ) == 0) 	{  KNL[3]=value; }
       else if (strcmp( name, "K4L" ) == 0) 	{  KNL[4]=value; }      
       else if (strcmp( name, "K5L" ) == 0) 	{  KNL[5]=value; }
       else if (strcmp( name, "K6L" ) == 0) 	{  KNL[6]=value; }
       else if (strcmp( name, "K7L" ) == 0) 	{  KNL[7]=value; }      
       else if (strcmp( name, "K8L" ) == 0) 	{  KNL[8]=value; }
       else if (strcmp( name, "K9L" ) == 0) 	{  KNL[9]=value; }      
       else if (strcmp( name, "K10L" ) == 0)    {  KNL[10]=value; }
       else if (strcmp( name, "K0SL" ) == 0)    {  KNSL[0]=value; }
       else if (strcmp( name, "K1SL" ) == 0)    {  KNSL[1]=value; }      
       else if (strcmp( name, "K2SL" ) == 0)    {  KNSL[2]=value; }
       else if (strcmp( name, "K3SL" ) == 0)    {  KNSL[3]=value; }
       else if (strcmp( name, "K4SL" ) == 0)    {  KNSL[4]=value; }      
       else if (strcmp( name, "K5SL" ) == 0)    {  KNSL[5]=value; }
       else if (strcmp( name, "K6SL" ) == 0)    {  KNSL[6]=value; }
       else if (strcmp( name, "K7SL" ) == 0)    {  KNSL[7]=value; }      
       else if (strcmp( name, "K8SL" ) == 0)    {  KNSL[8]=value; }
       else if (strcmp( name, "K9SL" ) == 0)    {  KNSL[9]=value; }      
       else if (strcmp( name, "K10SL" ) == 0)   {  KNSL[10]=value;}
       else if (strcmp( name, "Nint" ) == 0)    {  Nint=int(value);}
       else 
	 {
	   cout<<"MULT does not have parameter of  "<<name<<endl; 
	   exit(0); 
	 } 
     }
  double GetP(const char *name)
    {
      if (     strcmp( name, "K0L" ) == 0)  { return   KNL[0]; }
      else if (strcmp( name, "K1L" ) == 0)  { return   KNL[1]; }      
      else if (strcmp( name, "K2L" ) == 0)  { return   KNL[2]; }
      else if (strcmp( name, "K3L" ) == 0)  { return   KNL[3]; }
      else if (strcmp( name, "K4L" ) == 0)  { return   KNL[4]; }      
      else if (strcmp( name, "K5L" ) == 0)  { return   KNL[5]; }
      else if (strcmp( name, "K6L" ) == 0)  { return   KNL[6]; }
      else if (strcmp( name, "K7L" ) == 0)  { return   KNL[7]; }      
      else if (strcmp( name, "K8L" ) == 0)  { return   KNL[8]; }
      else if (strcmp( name, "K9L" ) == 0)  { return   KNL[9]; }      
      else if (strcmp( name, "K10L" ) == 0) { return   KNL[10]; }
      else if (strcmp( name, "K0SL" ) == 0) { return   KNSL[0]; }
      else if (strcmp( name, "K1SL" ) == 0) { return   KNSL[1]; }      
      else if (strcmp( name, "K2SL" ) == 0) { return   KNSL[2]; }
      else if (strcmp( name, "K3SL" ) == 0) { return   KNSL[3]; }
      else if (strcmp( name, "K4SL" ) == 0) { return   KNSL[4]; }      
      else if (strcmp( name, "K5SL" ) == 0) { return   KNSL[5]; }
      else if (strcmp( name, "K6SL" ) == 0) { return   KNSL[6]; }
      else if (strcmp( name, "K7SL" ) == 0) { return   KNSL[7]; }      
      else if (strcmp( name, "K8SL" ) == 0) { return   KNSL[8]; }
      else if (strcmp( name, "K9SL" ) == 0) { return   KNSL[9]; }      
      else if (strcmp( name, "K10SL" ) == 0){ return   KNSL[10];}
      else if (strcmp( name, "Norder" ) == 0){ return  Norder;}
      else if (strcmp( name, "Nint" ) == 0)  { return  Nint;}
      else 
	{
	  cout<<"MULT does not have parameter of  "<<name<<endl; 
	  exit(0); 
	}   
    }
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);   MULT_Pass(x,L,Nint, Norder,KNL,KNSL);    LtoG(x, DX,DY,DT);
  }
  void   DAPass(tps x[6]){
    GtoL(x,DX,DY,DT);   MULT_Pass(x,L,Nint, Norder,KNL,KNSL);    LtoG(x, DX,DY,DT);
  }
 private:
  double  KNL[11], KNSL[11];
  int Nint, Norder;
};

//-------------------------------GMULT--------------------------------------
class GMULT: public Element
{
 public:
 GMULT(string name, double l, double angle, double e1, double e2, double knl[11], double knsl[11]): Element(name)
    { 
      int i;
      if (l >=0  )
	{
	  TYPE=string("GMULT");
	  L=l;
	  ANGLE=angle;
          E1=e1;
          E2=e2;
	  Nint=int(L/QLslice)+1;
          for(i=0;i<11;i++) {
	    KNL[i]=knl[i];  KNSL[i]=knsl[i]; }
	  Norder=1;
	  for(i=0;i<10;i++) {
	    if( KNL[10-i] != 0. ||  KNSL[10-i] != 0. ){
	      Norder=10-i;
	      break;
	    }
	  }
	}
      else
	{
	  cout<<"Error: GMULT length can not be negative."<<endl;
          exit(1); 
	}
    }
   void SetP(const char *name, double value)
     { 
       if (strcmp( name, "ANGLE" ) == 0) 
	 {
	   ANGLE=value;
	 }
       else if (strcmp( name, "E1" ) == 0) 
	 {
	   E1=value;
	 }
       else if (strcmp( name, "E2" ) == 0) 
	 {
	   E2=value;
	 }
       else if (strcmp( name, "Nint" ) == 0) 
	 {
	   Nint=int(value);
	 }
       else if (strcmp( name, "K0L" ) == 0) 	{  KNL[0]=value; }
       else if (strcmp( name, "K1L" ) == 0) 	{  KNL[1]=value; }      
       else if (strcmp( name, "K2L" ) == 0) 	{  KNL[2]=value; }
       else if (strcmp( name, "K3L" ) == 0) 	{  KNL[3]=value; }
       else if (strcmp( name, "K4L" ) == 0) 	{  KNL[4]=value; }      
       else if (strcmp( name, "K5L" ) == 0) 	{  KNL[5]=value; }
       else if (strcmp( name, "K6L" ) == 0) 	{  KNL[6]=value; }
       else if (strcmp( name, "K7L" ) == 0) 	{  KNL[7]=value; }      
       else if (strcmp( name, "K8L" ) == 0) 	{  KNL[8]=value; }
       else if (strcmp( name, "K9L" ) == 0) 	{  KNL[9]=value; }      
       else if (strcmp( name, "K10L" ) == 0)    {  KNL[10]=value; }
       else if (strcmp( name, "K0SL" ) == 0)    {  KNSL[0]=value; }
       else if (strcmp( name, "K1SL" ) == 0)    {  KNSL[1]=value; }      
       else if (strcmp( name, "K2SL" ) == 0)    {  KNSL[2]=value; }
       else if (strcmp( name, "K3SL" ) == 0)    {  KNSL[3]=value; }
       else if (strcmp( name, "K4SL" ) == 0)    {  KNSL[4]=value; }      
       else if (strcmp( name, "K5SL" ) == 0)    {  KNSL[5]=value; }
       else if (strcmp( name, "K6SL" ) == 0)    {  KNSL[6]=value; }
       else if (strcmp( name, "K7SL" ) == 0)    {  KNSL[7]=value; }      
       else if (strcmp( name, "K8SL" ) == 0)    {  KNSL[8]=value; }
       else if (strcmp( name, "K9SL" ) == 0)    {  KNSL[9]=value; }      
       else if (strcmp( name, "K10SL" ) == 0)   {  KNSL[10]=value;}
       else if (strcmp( name, "Nint" ) == 0)    {  Nint=int(value);}
       else 
	 {
	   cout<<"MULT does not have parameter of  "<<name<<endl; 
	   exit(0); 
	 } 
     }
   double GetP(const char *name)
   {
     if (strcmp( name, "ANGLE" ) == 0) 
       {
	 return ANGLE;
       }
     else if (strcmp( name, "E1" ) == 0) 
       {
	 return E1;
       }
     else if (strcmp( name, "E2" ) == 0) 
       {
	 return E2;
       }
     else if (strcmp( name, "Nint" ) == 0) 
       {
	 return Nint;
       }
     else if (strcmp( name, "K0L" ) == 0)  { return   KNL[0]; }
     else if (strcmp( name, "K1L" ) == 0)  { return   KNL[1]; }      
     else if (strcmp( name, "K2L" ) == 0)  { return   KNL[2]; }
     else if (strcmp( name, "K3L" ) == 0)  { return   KNL[3]; }
     else if (strcmp( name, "K4L" ) == 0)  { return   KNL[4]; }      
     else if (strcmp( name, "K5L" ) == 0)  { return   KNL[5]; }
     else if (strcmp( name, "K6L" ) == 0)  { return   KNL[6]; }
     else if (strcmp( name, "K7L" ) == 0)  { return   KNL[7]; }      
     else if (strcmp( name, "K8L" ) == 0)  { return   KNL[8]; }
     else if (strcmp( name, "K9L" ) == 0)  { return   KNL[9]; }      
     else if (strcmp( name, "K10L" ) == 0) { return   KNL[10]; }
     else if (strcmp( name, "K0SL" ) == 0) { return   KNSL[0]; }
     else if (strcmp( name, "K1SL" ) == 0) { return   KNSL[1]; }      
     else if (strcmp( name, "K2SL" ) == 0) { return   KNSL[2]; }
     else if (strcmp( name, "K3SL" ) == 0) { return   KNSL[3]; }
     else if (strcmp( name, "K4SL" ) == 0) { return   KNSL[4]; }      
     else if (strcmp( name, "K5SL" ) == 0) { return   KNSL[5]; }
     else if (strcmp( name, "K6SL" ) == 0) { return   KNSL[6]; }
     else if (strcmp( name, "K7SL" ) == 0) { return   KNSL[7]; }      
     else if (strcmp( name, "K8SL" ) == 0) { return   KNSL[8]; }
     else if (strcmp( name, "K9SL" ) == 0) { return   KNSL[9]; }      
     else if (strcmp( name, "K10SL" ) == 0){ return   KNSL[10];}
     else if (strcmp( name, "Norder" ) == 0){ return  Norder;}
     else if (strcmp( name, "Nint" ) == 0)  { return  Nint;}
     else 
       {
	 cout<<"MULT does not have parameter of  "<<name<<endl; 
	 exit(0); 
       }   
   }
   void   Pass(double x[6]){
     GtoL(x,DX,DY,DT);   GMULT_Pass(x,L,Nint,ANGLE,Norder,KNL,KNSL,E1,E2);    LtoG(x, DX,DY,DT);
   }
   void   DAPass(tps x[6]){
     GtoL(x,DX,DY,DT);   GMULT_Pass(x,L,Nint,ANGLE,Norder,KNL,KNSL,E1,E2);    LtoG(x, DX,DY,DT);
   }
 private:
   double  ANGLE, E1,E2, KNL[11], KNSL[11];
   int Nint, Norder;
};

//---------------------------KICKER-----------------------------------------
template <class T> 
void KICK_Pass(T x[6], double L, double HKICK, double VKICK) 
{
    DRIFT_Pass(x,L/2.0); 
    x[1]=x[1]+HKICK;     
    x[3]=x[3]+VKICK;
    DRIFT_Pass(x,L/2.0);
}

class KICKER: public Element
{
 public:
  KICKER(string name, double l, double hkick, double vkick): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("KICKER");
	  L=l;
          HKICK=hkick;
          VKICK=vkick;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void SetP(const char *name, double value)
    { 
      if (strcmp( name, "HKICK" ) == 0) 
	{
	  HKICK=value;
	}
      else if (strcmp( name, "VKICK" ) == 0) 
	{
	  VKICK=value;
	}
      else 
	{
	  cout<<"KICKRE does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 
  double GetP(const char *name)
    {
      if (strcmp( name, "HKICK" ) == 0) 
	{
	  return HKICK;
	}
      else if (strcmp( name, "VKICK" ) == 0) 
	{
	  return VKICK;
	}
      else 
	{
	  cout<<"KICKRE does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }  
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);   KICK_Pass(x,L,HKICK, VKICK);   LtoG(x, DX,DY,DT);
  }
  void   DAPass(tps x[6]){
    GtoL(x,DX,DY,DT);  KICK_Pass(x,L,HKICK, VKICK);    LtoG(x, DX,DY,DT);
 }
 private:
  double  HKICK, VKICK;
};

//----------------------------HKICKER--------------------------------------
class HKICKER: public Element
{
 public:
  HKICKER(string name, double l, double hkick): Element(name)
    { 
      if (l >=0  )
	{
	  TYPE=string("HKICKER");
	  L=l;
          HKICK=hkick;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void SetP(const char *name, double value)
    { 
      if (strcmp( name, "HKICK" ) == 0) 
	{
	  HKICK=value;
	}
      else 
	{
	  cout<<"HKICKER does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 
  double GetP(const char *name)
    {
      if (strcmp( name, "HKICK" ) == 0) 
	{
	  return HKICK;
	}
      else 
	{
	  cout<<"HKICKER does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }  
  void   Pass(double x[6]){ 
    GtoL(x,DX,DY,DT);  KICK_Pass(x,L, HKICK, 0.);  LtoG(x, DX,DY,DT);
  }
  void   DAPass(tps x[6]){
    GtoL(x,DX,DY,DT);  KICK_Pass(x,L, HKICK, 0.);  LtoG(x, DX,DY,DT);
  }
 private:
  double  HKICK;   // kick in x'
};

//------------------------------VKICKER------------------------------------------
class VKICKER: public Element
{
 public:
  VKICKER(string name, double l, double vkick): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("VKICKER");
	  L=l;
          VKICK=vkick;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void SetP(const char *name, double value)
    { 
      if (strcmp( name, "VKICK" ) == 0) 
	{
	  VKICK=value;
	}
      else 
	{
	  cout<<"VKICKER does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 
  double GetP(const char *name)
    {
      if (strcmp( name, "VKICK" ) == 0) 
	{
	  return VKICK;
	}
      else 
	{
	  cout<<"KICKER does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }  
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);   KICK_Pass(x,L, 0., VKICK);  LtoG(x, DX,DY,DT);
  }
  void   DAPass(tps x[6]){
    GtoL(x,DX,DY,DT);   KICK_Pass(x,L, 0., VKICK);  LtoG(x, DX,DY,DT);
  }
 private:
  double  VKICK;  // kick in y'
};

//---------------------------HACKICK ( amplitude constant ac dipole kicking)------------------------------------
template <class T> 
void ACKICK_Pass(T x[6], double L, double HKICK, double VKICK, double TTURNS, double PHI0) 
{
    DRIFT_Pass(x,L/2.0); 
    x[1]=x[1]+HKICK *sin( 2.0 * PI * GP.turn / TTURNS + PHI0);
    x[3]=x[3]+VKICK *sin( 2.0 * PI * GP.turn / TTURNS + PHI0);
    DRIFT_Pass(x,L/2.0);
}

class HACKICK: public Element
{
 public:
 HACKICK(string name, double l, double hkick, int tturns, double phi0): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("HACKICK");
	  L=l;
          HKICK=hkick;
	  TTURNS=tturns;
	  PHI0=phi0;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void  SetP(const char *name, double value) 
    {
      if (strcmp( name, "HKICK" ) == 0) 
	{
	  HKICK=value;
	}
      else if (strcmp( name, "TTURNS" ) == 0) 
	{
	  TTURNS=int(value);
	}
      else if (strcmp( name, "PHI0" ) == 0) 
	{
	  PHI0=value;
	}
      else 
	{
	  cout<<"HACKICK does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }
  double GetP(const char *name) 
  {
    if (strcmp( name, "HKICK" ) == 0) 
      {
	return HKICK;
      }
    else if (strcmp( name, "TTURNS" ) == 0) 
      {
	return TTURNS;
      }
    else if (strcmp( name, "PHI0" ) == 0) 
      {
	return PHI0;
      }
    else 
      {
	cout<<"HACKICK  does not have  parameter of "<<name<<endl; 
	exit(0); 
      } 
  }
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);  ACKICK_Pass(X,L,HKICK,0.,TTURNS,PHI0); LtoG(x,DX,DY,DT); 
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x, L);
  }
 private:
  int    TTURNS;              //  T period in unit of turns
  double HKICK, PHI0;         //  HKICK is kick amplitude
};

//---------------------------VACKICK ( amplitude constant ac dipole kicking)------------------------------------
class VACKICK: public Element
{
 public:
 VACKICK(string name, double l, double vkick, int tturns, double phi0): Element(name)
  { 
    if (l >=0 )
      {
	TYPE=string("VACKICK");
	L=l;
	VKICK=vkick;
	TTURNS=int(tturns);
	PHI0=phi0;
      }
    else
      {
	cout<<"Error: check signs of L and S"<<endl;
	exit(1); 
      }
  }
  void  SetP(const char *name, double value) 
  {
    if (strcmp( name, "VKICK" ) == 0) 
      {
	VKICK=value;
      }
    else if (strcmp( name, "TTURNS" ) == 0) 
      {
	TTURNS=int(value);
      }
    else if (strcmp( name, "PHI0" ) == 0) 
      {
	PHI0=value;
      }
    else 
      {
	cout<<"VACKICK does not have  parameter of "<<name<<endl; 
	exit(0); 
      } 
  }
  double GetP(const char *name) 
  {
    if (strcmp( name, "VKICK" ) == 0) 
      {
	return VKICK;
      }
    else if (strcmp( name, "TTURNS" ) == 0) 
      {
	return TTURNS;
      }
    else if (strcmp( name, "PHI0" ) == 0) 
      {
	return PHI0;
      }
    else 
      {
	cout<<"VACKICK  does not have  parameter of "<<name<<endl; 
	exit(0); 
      } 
  }
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);  ACKICK_Pass(X,L,0.,VKICK,TTURNS,PHI0);  LtoG(x,DX,DY,DT); 
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x, L);
  }
 private:
  int    TTURNS;            // T period in unit of turns 
  double VKICK, PHI0;       
};

//------------------------------HACDIP ( slowly ramp up the kick amplitude)-----------------------------------------
template <class T> 
void ACDIP_Pass(T x[6], double L, double HKICKMAX, double VKICKMAX, double NUD, double TURNS, double TURNE, double PHID) 
{
  if( GP.turn < TURNS ) {
    DRIFT_Pass(x, L);
  }
  else if ( GP.turn >= TURNS  and GP.turn <= TURNE  ){
    DRIFT_Pass(x, L/2.);
    x[1] += ((GP.turn-TURNS)*1.0 * HKICKMAX /(TURNE-TURNS )) * sin( 2.0* 3.14159265*NUD*(GP.turn-TURNS)+PHID);
    x[3] += ((GP.turn-TURNS)*1.0 * VKICKMAX /(TURNE-TURNS )) * sin( 2.0* 3.14159265*NUD*(GP.turn-TURNS)+PHID);
    DRIFT_Pass(x, L/2.);
  }
  else {
    DRIFT_Pass(x, L/2.);
    x[1] +=  HKICKMAX* sin( 2.0* 3.14159265 * NUD *(GP.turn-TURNS) +PHID );
    x[3] +=  VKICKMAX* sin( 2.0* 3.14159265 * NUD *(GP.turn-TURNS) +PHID );
    DRIFT_Pass(x, L/2.);
  }
}

class HACDIP: public Element
{
 public:
  HACDIP(string name, double l, double hkickmax, double nud, double phid, int turns, int turne ): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("HACDIP");
	  L=l;
	  HKICKMAX=hkickmax;
          NUD=nud;
	  PHID=phid;
	  TURNS=turns;
          TURNE=turne;
	}
      else
	{
	  cout<<"Error: check signs of L. "<<endl;
          exit(1); 
	}
    }
  void SetP(const char *name, double value)
    { 
      if (strcmp( name, "HKICKMAX" ) == 0) 
	{
	  HKICKMAX=value;
	}
      else if (strcmp( name, "NUD" ) == 0) 
	{
	  NUD=value;
	}
      else if (strcmp( name, "PHID" ) == 0) 
	{
	  PHID=value;
	}
      else if (strcmp( name, "TURNS" ) == 0) 
	{
	  TURNS=int(value);
	}
      else if (strcmp( name, "TURNE" ) == 0) 
	{
	  TURNE=int(value);
	}

      else 
	{
	  cout<<"HACDIP does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 
  double GetP(const char *name)
    { 
      if (strcmp( name, "HKICKMAX" ) == 0) 
	{
	  return HKICKMAX;
	}
      else if (strcmp( name, "TURNS" ) == 0) 
	{
	  return TURNS;
	}
      else if (strcmp( name, "TURNE" ) == 0) 
	{
	  return TURNE;
	}
      else if (strcmp( name, "NUD" ) == 0) 
	{
	  return NUD;
	}
      else if (strcmp( name, "PHID" ) == 0) 
	{
	  return PHID;
	}
      else 
	{
	  cout<<"HACDIP does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 
  void   Pass(double x[6]){
     GtoL(x,DX,DY,DT); ACDIP_Pass(x, L, HKICKMAX, 0., NUD, TURNS, TURNE, PHID);  LtoG(x,DX,DY,DT); 
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x,L);
  }
 private:
  int     TURNS, TURNE;          // kick strength ramping between  TURNS and TURNE turns
  double  HKICKMAX, NUD, PHID;   // HKICKMAX, maximum kick in x'
};

//------------------------------VACDIP------------------------------------------
class VACDIP: public Element
{
 public:
  VACDIP(string name, double l, double vkickmax, double nud, double phid, int turns, int turne ): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("VACDIP");
	  L=l;
	  VKICKMAX=vkickmax;
          NUD=nud;
	  PHID=phid;
	  TURNS=turns;
          TURNE=turne;
	  
	}
      else
	{
	  cout<<"Error: check signs of L. "<<endl;
          exit(1); 
	}
    }
    void SetP(const char *name, double value)
    { 
      if (strcmp( name, "VKICKMAX" ) == 0) 
	{
	  VKICKMAX=value;
	}
      else if (strcmp( name, "NUD" ) == 0) 
	{
	  NUD=value;
	}
      else if (strcmp( name, "PHID" ) == 0) 
	{
	  PHID=value;
	}
      else if (strcmp( name, "TURNS" ) == 0) 
	{
	  TURNS=int(value);
	}
      else if (strcmp( name, "TURNE" ) == 0) 
	{
	  TURNE=int(value);
	}
      else 
	{
	  cout<<"VACDIP does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 
    double GetP(const char *name)
    { 
      if (strcmp( name, "VKICKMAX" ) == 0) 
	{
	  return VKICKMAX;
	}
      else if (strcmp( name, "TURNS" ) == 0) 
	{
	  return TURNS;
	}
      else if (strcmp( name, "TURNE" ) == 0) 
	{
	  return TURNE;
	}
      else if (strcmp( name, "NUD" ) == 0) 
	{
	  return NUD;
	}
      else if (strcmp( name, "PHID" ) == 0) 
	{
	  return PHID;
	}
      else 
	{
	  cout<<"VACDIP does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 
    void   Pass(double x[6]){
      GtoL(x,DX,DY,DT);  ACDIP_Pass(x, L, 0.,VKICKMAX,NUD, TURNS, TURNE, PHID);   LtoG(x,DX,DY,DT);
    }
    void   DAPass(tps x[6]){
      DRIFT_Pass(x,L);
    }
 private:
    int     TURNS, TURNE;
    double  VKICKMAX, NUD, PHID;  // HKICKMAX, maximum kick in y'
};

//----------------------------ACMULT ( single order, order >= quadrupole )-----------------
void  ACMULT_Pass(double x[6], double L, int Norder, double KL, int TTURNS, double PHI0)
{
  int     i;
  int     fac=1;
  double  Xn,  Yn, Xn0, Yn0;
  double  KNL, KNSL;
  double  By, Bx;

  DRIFT_Pass(x, L/2.0);

  if(Norder < 0 ){
     KNL=0.;
     KNSL=KL;
  }
  else{
     KNL=KL;
     KNSL=0.;
  }
  
  By=0.;    
  Bx=0.;
  Xn=1.;
  Yn=0.;

  for(i=1;i<Norder+1;i++){
    Xn0=Xn;
    Yn0=Yn;
    Xn=Xn0*x[x_]-Yn0*x[y_];
    Yn=Xn0*x[y_]+Yn0*x[x_];
    fac=fac*i;
  }
  By=By+(KNL*Xn-KNSL*Yn)/fac;
  Bx=Bx+(KNL*Yn+KNSL*Xn)/fac;
  
  x[px_]=x[px_]- By * sin( 2.0 * PI * GP.turn / TTURNS + PHI0);
  x[py_]=x[py_]+ Bx * sin( 2.0 * PI * GP.turn / TTURNS + PHI0);

  DRIFT_Pass(x, L/2.0);
}

class ACMULT: public Element
{
 public:
 ACMULT(string name, double l, int norder, double kl, int tturns, double phi0 ): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("ACMULT");
	  L=l;
          Norder=norder;
          KL=kl;
          TTURNS=tturns;
          PHI0=phi0;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void  SetP(const char *name, double value) 
    {
      if (strcmp( name, "KL" ) == 0) 
	{
	  KL=value;
	}
      else if (strcmp( name, "Norder" ) == 0) 
	{
          if ( value == 0) {
	    cout<<" Use ACKICK for Norder = 0."<<endl;
	    exit(0);
	  }
	  else if ( abs(value) > 10 ){
            cout<<" The maximum order of ACMULT is 10."<<endl;
            exit(1);
	  }
	  else{
	    Norder=int(value);
	  }
	}
      else if (strcmp( name, "TTURNS" ) == 0) 
	{
	  TTURNS=int(value);
	}
      else if (strcmp( name, "PHI0" ) == 0) 
	{
	  PHI0=value;
	}
      else 
	{
	  cout<<"ACMULT does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }
  double GetP(const char *name) 
  {
    if (strcmp( name, "KL" ) == 0) 
      {
	return KL;
      }
    else if (strcmp( name, "Norder" ) == 0) 
      {
	return Norder;
      }
    else if (strcmp( name, "TTURNS" ) == 0) 
      {
	return TTURNS;
      }
    else if (strcmp( name, "PHI0" ) == 0) 
      {
	return PHI0;
      }
    else 
      {
	cout<<"ACMULT  does not have  parameter of "<<name<<endl; 
	exit(0); 
      } 
  }
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);  ACMULT_Pass(x, L, Norder, KL, TTURNS, PHI0);   LtoG(x,DX,DY,DT); 
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x, L);
  }
 private:
  int    Norder, TTURNS;            // TTURNs period in numbers of turns, 
  double KL, PHI0;                  //  Norder >0, KL->KNL, Norder <0, KL->KNSL 
};

//---------------------------COOLING--------------------------------
template <class T> 
void COOLING_Pass(T x[6], double L, double ALPHA)
{
  int i;
  DRIFT_Pass(x, L/2.);   
  for(i=0;i<4;i++) x[i]= (1.- ALPHA )* x[i];
  DRIFT_Pass(x, L/2.);
}

class COOLING: public Element
{
 public:
  COOLING(string name, double l, double alpha ): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("COOLING");
	  L=l;
	  ALPHA=alpha;
	}
      else
	{
	  cout<<"Error: check signs of L and S of COOLING element."<<endl;
          exit(1);
	}
    }
    void  SetP(const char *name, double value) 
    {
      if (strcmp( name, "ALPHA" ) == 0) 
	{
	  ALPHA=value;
	}
      else 
	{
	  cout<<"COOLING does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }
  double GetP(const char *name) 
  {
    if (strcmp( name, "ALPHA" ) == 0) 
      {
	return ALPHA;
      }
    else 
      {
	cout<<"COOLING does not have  parameter of "<<name<<endl; 
	exit(0); 
      } 
  }
  void   Pass(double x[6]){
    int i;
    GtoL(x,DX,DY,DT);    COOLING_Pass(x, L, ALPHA);    LtoG(x,DX,DY,DT); 
  }
  void   DAPass(tps x[6]){ 
    DRIFT_Pass(x, L);
  }
 private:
  double ALPHA;  
};

//---------------------------DIFFUSE--------------------------------
template <class T> 
void DIFFUSE_Pass(T x[6], double L,  double DIFF_X, double  DIFF_Y, double DIFF_DELTA)
{
  DRIFT_Pass(x, L/2.);   
  x[1]= x[1]+ DIFF_X * rnd(seed);
  x[3]= x[3]+ DIFF_Y * rnd(seed);
  x[5]= x[5]+ DIFF_DELTA *rnd(seed);
  DRIFT_Pass(x, L/2.);
}

class DIFFUSE: public Element
{
 public:
  DIFFUSE(string name, double l, double diff_x, double diff_y, double diff_delta ): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("DIFFUSE");
	  L=l;
	  DIFF_X = diff_x;  DIFF_Y = diff_y; DIFF_DELTA = diff_delta;
	}
      else
	{
	  cout<<"Error: check signs of L and S of DIFFUSE element."<<endl;
          exit(1);
	}
    }
    void  SetP(const char *name, double value) 
    {
      if (strcmp( name, "DIFF_X" ) == 0) 
	{
	 DIFF_X = value;
	}
      else if (strcmp( name, "DIFF_Y" ) == 0) 
	{
	 DIFF_Y = value;
	}
      else if (strcmp( name, "DIFF_delta" ) == 0) 
	{
	 DIFF_DELTA = value;
	}
      else 
	{
	  cout<<"DIFFUSE does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }
  double GetP(const char *name) 
  {
    if (strcmp( name, "DIFF_X" ) == 0) 
      {
	return DIFF_X;
      }
    else if (strcmp( name, "DIFF_Y" ) == 0) 
      {
	return DIFF_Y;
      }
    else if (strcmp( name, "DIFF_DELTA" ) == 0) 
      {
	return DIFF_DELTA;
      }
    else 
      {
	cout<<"DIFFUSE does not have  parameter of "<<name<<endl; 
	exit(0); 
      } 
  }
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);   DIFFUSE_Pass(x, L, DIFF_X, DIFF_Y, DIFF_DELTA);   LtoG(x,DX,DY,DT); 
  }
  void   DAPass(tps x[6]){ 
    DRIFT_Pass(x, L);
  }
 private:
  double DIFF_X, DIFF_Y, DIFF_DELTA;  
};

//---------------------------------BPM-------------------------
class BPM: public Element
{
 public:
 BPM(string name, double l): Element(name)
  { 
    if (l >=0 )
      {
	TYPE=string("BPM");
	L=l;
      }
    else
      {
	cout<<"Error: check signs of L and S"<<endl;
	exit(1); 
      }
  }
  void  SetP(const char *name, double value) 
  {
    cout<<"No parameter to be set for BPM."<<endl;
    exit(1);
  }
  double GetP(const char *name) 
  {
    cout<<"No parameter to be returned for BPM."<<endl;
    exit(1);
  }
  void   Pass(double x[6]){
    DRIFT_Pass(x, L);
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x, L);
  }
};

//----------------------------HBPM-----------------------------------
class HBPM: public Element
{
 public:
  HBPM(string name, double l): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("HBPM");
	  L=l;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void  SetP(const char *name, double value) 
    {
      cout<<"No parameter to be set for HBPM."<<endl;
      exit(1);
    }
  double GetP(const char *name) 
    {
      cout<<"No parameter to be returned for HBPM."<<endl;
      exit(1);
    }
  void   Pass(double x[6]){
    DRIFT_Pass(x, L);
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x, L);
  }
};

//---------------------------------VBPM------------------------------
class VBPM: public Element
{
 public:
  VBPM(string name, double l): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("VBPM");
	  L=l;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void  SetP(const char *name, double value) 
    {
      cout<<"No parameter to be set for VBPM."<<endl;
      exit(1);
    }
  double GetP(const char *name) 
    {
      cout<<"No parameter to be returned for VBPM."<<endl;
      exit(0);
    }
  void   Pass(double x[6]){
    DRIFT_Pass(x, L);
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x, L);
  }
};

//-----------------------------MARKER--------------------------
class MARKER: public Element
{
 public:
  MARKER(string name, double l): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("MARKER");
	  L=l;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void  SetP(const char *name, double value) 
    {
      cout<<"No parameter to be set for MARKER."<<endl;
      exit(1);
    }
  double GetP(const char *name) 
    {
      cout<<"No parameter to be returned for MARKER."<<endl;
      exit(1);
    }
  void   Pass(double x[6]){
    DRIFT_Pass(x, L);   
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x, L);   
  }
};

//---------------------------RCOLL--------------------------------
class RCOLL: public Element
{
 public:
  RCOLL(string name, double l): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("RCOLL");
	  L=l;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void  SetP(const char *name, double value) 
  {
    cout<<"No parameter to be set for MARKER."<<endl;
    exit(1);
  }
  double GetP(const char *name) 
  {
    cout<<"No parameter to be returned for MARKER."<<endl;
    exit(1);
  }
  void   Pass(double x[6]){
    DRIFT_Pass(x, L);   
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x, L);   
  }
};

//----------------------------ECOLL-------------------------------
class ECOLL: public Element
{
 public:
  ECOLL(string name, double l): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("ECOLL");
	  L=l;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void  SetP(const char *name, double value) 
  {
    cout<<"No parameter to be set for MARKER."<<endl;
    exit(1);
  }
  double GetP(const char *name) 
  {
    cout<<"No parameter to be returned for MARKER."<<endl;
    exit(1);
  }
  void   Pass(double x[6]){
    DRIFT_Pass(x, L);  
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x, L);  
  }
};

//----------------------INSTR------------------------------------
class INSTR: public Element
{
 public:
  INSTR(string name, double l): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("INSTR");
	  L=l;
	}
      else
	{
	  cout<<"Error: check signs of L and S"<<endl;
          exit(1); 
	}
    }
  void   SetP(const char *name, double value) 
    {
      cout<<"No parameter is set for INSTR."<<endl;
      exit(0);
    }
  double GetP(const char *name) 
    {
      cout<<"No parameter is returned for INSTR."<<endl;
      exit(1);
    }
  void   Pass(double x[6]){
    DRIFT_Pass(x, L); 
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x, L);  
  }
};

//----------------------------SOLEN-----------------------------------
template <class T> 
void SOLEN_Pass(T x[6], double L, double KS, double M[36])
{
    int i, j, k;
    T xtemp[6];

    if(KS == 0 ) {
      DRIFT_Pass(x, L);}
    else{
      x[1]=x[1]/(1+x[5]);  // x'=px / (1+delta)
      x[3]=x[3]/(1+x[5]);
      for(i=0;i<6;i++) xtemp[i]=x[i];

      for(k=0;k<6;k++) {
	x[k]=0.;
	for(j=0;j<6;j++) x[k]= x[k]+ M[k*6+j]*xtemp[j];
      }
      x[1]=x[1]*(1+x[5]); // px=x' * (1+delta)
      x[3]=x[3]*(1+x[5]);
    }
}

class SOLEN: public Element
{
 public:
  SOLEN(string name, double l, double ks): Element(name)
    { 
      if ( l >0 && ks !=0.)
	{
	  TYPE=string("SOLEN");
	  L=l;
          KS=ks;
	  
	  int i, j;	  
	  double g=KS/2, theta=KS*L/2, costheta= cos( theta ), sintheta=  sin( theta ) ;
	  
          if( KS != 0. ) {
	    for(i=0;i<6;i++)
	      for(j=0;j<6;j++) M[i*6+j]= 0.0;
	    
	    M[0*6+0] = costheta * costheta;
	    M[0*6+1] = sintheta * costheta / g; 
	    M[0*6+2] = sintheta * costheta;
	    M[0*6+3] = sintheta *  sintheta / g;
	    
	    M[1*6+0] = -g * sintheta * costheta ;
	    M[1*6+1] =      costheta * costheta ;
	    M[1*6+2] = -g*  sintheta * sintheta ;
	    M[1*6+3] =      sintheta *  costheta;
	    
	    M[2*6+0] = -sintheta * costheta ;
	    M[2*6+1] = -sintheta * sintheta / g ;
	    M[2*6+2] =  costheta * costheta ;
	    M[2*6+3] =  sintheta * costheta / g;
	    
	    M[3*6+0] = g *  sintheta * sintheta ;
	    M[3*6+1] =     -sintheta * costheta ;
	    M[3*6+2] = -g * sintheta * costheta ;
	    M[3*6+3] =      costheta * costheta;
	    
	    M[4*6+4] = 1.0;
	    M[5*6+5] = 1.0;
	  }
	}
      else if (l >0 && ks ==0.)
	{ 
	  TYPE=string("SOLEN");
	  L=l;
	  KS=0;
	}
      else
	{
	  cout<<"Solenoid: L must be positive. "<<endl;
          exit(1); 
	}
    }
  void   SetP(const char *name, double value) 
    {
     if (strcmp( name, "KS" ) == 0) 
	{
	  KS=value;
	  int i, j;	  
	  double g=KS/2, theta=KS*L/2, costheta= cos( theta ), sintheta=  sin( theta ) ;
	  
          if( KS  != 0. ) {
	    for(i=0;i<6;i++)
	      for(j=0;j<6;j++) M[i*6+j]= 0.0;
	    
	    M[0*6+0] = costheta * costheta;
	    M[0*6+1] = sintheta * costheta / g; 
	    M[0*6+2] = sintheta * costheta;
	    M[0*6+3] = sintheta *  sintheta / g;
	    
	    M[1*6+0] = -g * sintheta * costheta ;
	    M[1*6+1] =      costheta * costheta ;
	    M[1*6+2] = -g*  sintheta * sintheta ;
	    M[1*6+3] =      sintheta *  costheta;
	    
	    M[2*6+0] = -sintheta * costheta ;
	    M[2*6+1] = -sintheta * sintheta / g ;
	    M[2*6+2] =  costheta * costheta ;
	    M[2*6+3] =  sintheta * costheta / g;
	    
	    M[3*6+0] = g *  sintheta * sintheta ;
	    M[3*6+1] =     -sintheta * costheta ;
	    M[3*6+2] = -g * sintheta * costheta ;
	    M[3*6+3] =      costheta * costheta;
	    
	    M[4*6+4] = 1.0;
	    M[5*6+5] = 1.0;
	  }
	  else
	    {
	      KS=0.;
	    }
	}
     else 
       {
	 cout<<"SOLEN does not have  parameter of "<<name<<endl; 
	 exit(0); 
       } 
    }
  double GetP(const char *name) 
    {
     if (strcmp( name, "KS" ) == 0) 
	{
	  return KS;
	}
      else 
	{
	  cout<<"SOLEN does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);   SOLEN_Pass(x, L, KS, M);    LtoG(x,DX,DY,DT); 
  }
  void   DAPass(tps x[6]){
    GtoL(x,DX,DY,DT);   SOLEN_Pass(x, L, KS, M);    LtoG(x,DX,DY,DT); 
  }
 private:
  double  KS, M[36];
};

//-----------------------MATRIX-----------------------------------
template<typename T>
void MATRIX_Pass( T x[6], double XCO_IN[6], double XCO_OUT[6], double M66[36])
{
  int i,j;
  T x1[6];

  for(i=0;i<6;i++) x[i]=x[i]-XCO_IN[i];

  for(i=0;i<6;i++) {
    x1[i]=0.0;
    for(j=0;j<6;j++)  x1[i] = x1[i] + M66[i*6+j] * x[j];
  }   
  for(i=0;i<6;i++)    x[i]=x1[i];
  
  for(i=0;i<6;i++) x[i]=x[i]+XCO_OUT[i];
}

class MATRIX: public Element
{
 public:
  MATRIX(string name, double l, double m66[36], double xco_in[6], double xco_out[6]): Element(name)
    {
      int i; 
      TYPE=string("MATRIX");
      L=l;
      for(i=0;i<36;i++) M66[i]=m66[i];
      for(i=0;i<6;i++)  XCO_IN[i]=xco_in[i];
      for(i=0;i<6;i++)  XCO_OUT[i]=xco_out[i];
    }
    void   SetP(const char *name, double value) 
    {
      if (     strcmp( name, "M11" ) == 0) 	{  M66[0]=value; }
      else if (strcmp( name, "M12" ) == 0) 	{  M66[1]=value; }      
      else if (strcmp( name, "M13" ) == 0) 	{  M66[2]=value; }
      else if (strcmp( name, "M14" ) == 0) 	{  M66[3]=value; }
      else if (strcmp( name, "M15" ) == 0) 	{  M66[4]=value; }      
      else if (strcmp( name, "M16" ) == 0) 	{  M66[5]=value; }
      else if (strcmp( name, "M21" ) == 0) 	{  M66[6]=value; }   
      else if (strcmp( name, "M22" ) == 0) 	{  M66[7]=value; }      
      else if (strcmp( name, "M23" ) == 0) 	{  M66[8]=value; }
      else if (strcmp( name, "M24" ) == 0) 	{  M66[9]=value; }
      else if (strcmp( name, "M25" ) == 0) 	{  M66[10]=value; }      
      else if (strcmp( name, "M26" ) == 0) 	{  M66[11]=value; }
      else if (strcmp( name, "M31" ) == 0) 	{  M66[12]=value; }   
      else if (strcmp( name, "M32" ) == 0) 	{  M66[13]=value; }      
      else if (strcmp( name, "M33" ) == 0) 	{  M66[14]=value; }
      else if (strcmp( name, "M34" ) == 0) 	{  M66[15]=value; }
      else if (strcmp( name, "M35" ) == 0) 	{  M66[16]=value; }      
      else if (strcmp( name, "M36" ) == 0) 	{  M66[17]=value; }
      else if (strcmp( name, "M41" ) == 0) 	{  M66[18]=value; }   
      else if (strcmp( name, "M42" ) == 0) 	{  M66[19]=value; }      
      else if (strcmp( name, "M43" ) == 0) 	{  M66[20]=value; }
      else if (strcmp( name, "M44" ) == 0) 	{  M66[21]=value; }
      else if (strcmp( name, "M45" ) == 0) 	{  M66[22]=value; }      
      else if (strcmp( name, "M46" ) == 0) 	{  M66[23]=value; }
      else if (strcmp( name, "M51" ) == 0) 	{  M66[24]=value; }   
      else if (strcmp( name, "M52" ) == 0) 	{  M66[25]=value; }      
      else if (strcmp( name, "M53" ) == 0) 	{  M66[26]=value; }
      else if (strcmp( name, "M54" ) == 0) 	{  M66[27]=value; }
      else if (strcmp( name, "M55" ) == 0) 	{  M66[28]=value; }      
      else if (strcmp( name, "M56" ) == 0) 	{  M66[29]=value; }
      else if (strcmp( name, "M61" ) == 0) 	{  M66[30]=value; }   
      else if (strcmp( name, "M62" ) == 0) 	{  M66[31]=value; }      
      else if (strcmp( name, "M63" ) == 0) 	{  M66[32]=value; }
      else if (strcmp( name, "M64" ) == 0) 	{  M66[33]=value; }
      else if (strcmp( name, "M65" ) == 0) 	{  M66[34]=value; }      
      else if (strcmp( name, "M66" ) == 0) 	{  M66[35]=value; }
      else if (strcmp( name, "XCO_IN_X" ) == 0)     {  XCO_IN[0]=value; }
      else if (strcmp( name, "XCO_IN_PX" ) == 0)    {  XCO_IN[1]=value; }
      else if (strcmp( name, "XCO_IN_Y" ) == 0)     {  XCO_IN[2]=value; }
      else if (strcmp( name, "XCO_IN_PY" ) == 0)    {  XCO_IN[3]=value; }
      else if (strcmp( name, "XCO_IN_Z" ) == 0)     {  XCO_IN[4]=value; }
      else if (strcmp( name, "XCO_IN_DELTA" ) == 0) {  XCO_IN[5]=value; }
      else if (strcmp( name, "XCO_OUT_X" ) == 0)    {  XCO_OUT[0]=value; }
      else if (strcmp( name, "XCO_OUT_PX" ) == 0)   {  XCO_OUT[1]=value; }
      else if (strcmp( name, "XCO_OUT_Y" ) == 0)    {  XCO_OUT[2]=value; }
      else if (strcmp( name, "XCO_OUT_PY" ) == 0)   {  XCO_OUT[3]=value; }
      else if (strcmp( name, "XCO_OUT_Z" ) == 0)    {  XCO_OUT[4]=value; }
      else if (strcmp( name, "XCO_OUT_DELTA" ) == 0){  XCO_OUT[5]=value; }
      else 
	{
	  cout<<"Matrix does not have parameter of  "<<name<<endl; 
	  exit(0); 
	} 
    }
    double GetP(const char *name) 
    {
      if (     strcmp( name, "M11" ) == 0) 	{  return M66[0]; }
      else if (strcmp( name, "M12" ) == 0) 	{  return M66[1]; }      
      else if (strcmp( name, "M13" ) == 0) 	{  return M66[2]; }
      else if (strcmp( name, "M14" ) == 0) 	{  return M66[3]; }
      else if (strcmp( name, "M15" ) == 0) 	{  return M66[4]; }      
      else if (strcmp( name, "M16" ) == 0) 	{  return M66[5]; }
      else if (strcmp( name, "M21" ) == 0) 	{  return M66[6]; }   
      else if (strcmp( name, "M22" ) == 0) 	{  return M66[7]; }      
      else if (strcmp( name, "M23" ) == 0) 	{  return M66[8]; }
      else if (strcmp( name, "M24" ) == 0) 	{  return M66[9]; }
      else if (strcmp( name, "M25" ) == 0) 	{  return M66[10]; }      
      else if (strcmp( name, "M26" ) == 0) 	{  return M66[11]; }
      else if (strcmp( name, "M31" ) == 0) 	{  return M66[12]; }   
      else if (strcmp( name, "M32" ) == 0) 	{  return M66[13]; }      
      else if (strcmp( name, "M33" ) == 0) 	{  return M66[14]; }
      else if (strcmp( name, "M34" ) == 0) 	{  return M66[15]; }
      else if (strcmp( name, "M35" ) == 0) 	{  return M66[16]; }      
      else if (strcmp( name, "M36" ) == 0) 	{  return M66[17]; }
      else if (strcmp( name, "M41" ) == 0) 	{  return M66[18]; }   
      else if (strcmp( name, "M42" ) == 0) 	{  return M66[19]; }      
      else if (strcmp( name, "M43" ) == 0) 	{  return M66[20]; }
      else if (strcmp( name, "M44" ) == 0) 	{  return M66[21]; }
      else if (strcmp( name, "M45" ) == 0) 	{  return M66[22]; }      
      else if (strcmp( name, "M46" ) == 0) 	{  return M66[23]; }
      else if (strcmp( name, "M51" ) == 0) 	{  return M66[24]; }   
      else if (strcmp( name, "M52" ) == 0) 	{  return M66[25]; }      
      else if (strcmp( name, "M53" ) == 0) 	{  return M66[26]; }
      else if (strcmp( name, "M54" ) == 0) 	{  return M66[27]; }
      else if (strcmp( name, "M55" ) == 0) 	{  return M66[28]; }      
      else if (strcmp( name, "M56" ) == 0) 	{  return M66[29]; }
      else if (strcmp( name, "M61" ) == 0) 	{  return M66[30]; }   
      else if (strcmp( name, "M62" ) == 0) 	{  return M66[31]; }      
      else if (strcmp( name, "M63" ) == 0) 	{  return M66[32]; }
      else if (strcmp( name, "M64" ) == 0) 	{  return M66[33]; }
      else if (strcmp( name, "M65" ) == 0) 	{  return M66[34]; }      
      else if (strcmp( name, "M66" ) == 0) 	{  return M66[35]; }
      else if (strcmp( name, "XCO_IN_X" ) == 0)     {  return XCO_IN[0]; }
      else if (strcmp( name, "XCO_IN_PX" ) == 0)    {  return XCO_IN[1]; }
      else if (strcmp( name, "XCO_IN_Y" ) == 0)     {  return XCO_IN[2]; }
      else if (strcmp( name, "XCO_IN_PY" ) == 0)    {  return XCO_IN[3]; }
      else if (strcmp( name, "XCO_IN_Z" ) == 0)     {  return XCO_IN[4]; }
      else if (strcmp( name, "XCO_IN_DELTA" ) == 0) {  return XCO_IN[5]; }
      else if (strcmp( name, "XCO_OUT_X" ) == 0)    {  return XCO_OUT[0]; }
      else if (strcmp( name, "XCO_OUT_PX" ) == 0)   {  return XCO_OUT[1]; }
      else if (strcmp( name, "XCO_OUT_Y" ) == 0)    {  return XCO_OUT[2]; }
      else if (strcmp( name, "XCO_OUT_PY" ) == 0)   {  return XCO_OUT[3]; }
      else if (strcmp( name, "XCO_OUT_Z" ) == 0)    {  return XCO_OUT[4]; }
      else if (strcmp( name, "XCO_OUT_DELTA" ) == 0){  return XCO_OUT[5]; }
      else 
	{
	  cout<<"Matrix does not have parameter of  "<<name<<endl; 
	  exit(0); 
	} 
    }
    void   Pass(double x[6]){
      GtoL(x,DX,DY,DT); MATRIX_Pass(x,XCO_IN, XCO_OUT, M66);  LtoG(x,DX,DY,DT); 
    }
    void   DAPass(tps x[6]){
      GtoL(x,DX,DY,DT); MATRIX_Pass(x, XCO_IN, XCO_OUT, M66); LtoG(x,DX,DY,DT);
    }
 private:
    double  M66[36];
    double  XCO_IN[6], XCO_OUT[6];
};

//--------------------ELSEP----------------------------------
class ELSEP: public Element
{
 public:
 ELSEP(string name, double l): Element(name)
  { 
    if (l >=0 )
      {
	TYPE=string("ELSEP");
	L=l;
      }
    else
      {
	cout<<"Error: check signs of L and S"<<endl;
	exit(1); 
      }
  }
  void   SetP(const char *name, double value) 
  {
    cout<<"No parameter to be set for ELSEP."<<endl;
    exit(1);
  }
  double GetP(const char *name) 
  {
    cout<<"No parameter to be returned for ELSEP."<<endl;
    exit(1);
  }
  void   Pass(double x[6]){
    DRIFT_Pass(x, L);
  }
  void   DAPass(tps x[6]){
    DRIFT_Pass(x, L); 
  }
};

//--------------------OFFSET (coordinate system change)-------------------------------
template<typename T>
void OFFSET_Pass( T x[6], double DX, double DPX,double DY, double DPY, double TILT)
{
  int i;
  double cosT=cos(-TILT), sinT=sin(-TILT);
  T xtemp[6];
  x[0]=x[0]-DX;
  x[1]=x[1]-DPX;
  x[2]=x[2]-DY;
  x[3]=x[3]-DPY;
  for (i=0; i<6;i++) xtemp[i]=x[i];
  x[0] = cosT * xtemp[0] + sinT * xtemp[2];
  x[1] = cosT * xtemp[1] + sinT * xtemp[3];
  x[2] = cosT * xtemp[2] - sinT * xtemp[0];
  x[3] = cosT * xtemp[3] - sinT * xtemp[1];
}

class OFFSET: public Element
{
 public:
  OFFSET(string name, double l, double dx, double dpx, double dy, double dpy, double tilt): Element(name)
    { 
      if (l ==0 )
	{
	  TYPE=string("OFFSET");
	  L=l;
          DX=  dx;
          DPX= dpx;
          DY=  dy;
          DPY= dpy;
          TILT=tilt;
	}
      else
	{
	  cout<<"Error: L for OFFSET must be zero."<<endl;
          exit(1); 
	}
    }
  void   SetP(const char *name, double value) 
    {
      if (strcmp( name, "DX" ) == 0) 
	{
	  DX=value;
	}
      else if (strcmp( name, "DPX" ) == 0) 
	{
	  DPX=value;
	}
      if (strcmp( name, "DY" ) == 0) 
	{
	  DY=value;
	}
      else if (strcmp( name, "DPY" ) == 0) 
	{
	  DPY=value;
	}
      else if (strcmp( name, "TILT" ) == 0) 
	{
	  TILT=value;
	}
      else 
	{
	  cout<<"OFFSET does not have  parameter of "<<name<<endl; 
          exit(0); 
	}        
    }
  double GetP(const char *name) 
    {
      if (strcmp( name, "DX" ) == 0) 
	{
	  return DX;
	}
      else if (strcmp( name, "DPX" ) == 0) 
	{
	  return DPX;
	}
      if (strcmp( name, "DY" ) == 0) 
	{
	  return DY;
	}
      else if (strcmp( name, "DPY" ) == 0) 
	{
	  return DPY;
	}
      else if (strcmp( name, "TILT" ) == 0) 
	{
	  return TILT;
	}
      else 
	{
	  cout<<"OFFSET does not have  parameter of "<<name<<endl; 
          exit(0); 
	}        
    }
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);  OFFSET_Pass(x, DX,DPX,DY,DPY,TILT);  LtoG(x,DX,DY,DT);
  }
  void   DAPass(tps x[6]){
    GtoL(x,DX,DY,DT);  OFFSET_Pass(x, DX,DPX,DY,DPY,TILT);  LtoG(x,DX,DY,DT);
  }
 private:
  double DX, DPX, DY, DPY, TILT;
};

//----------------------------RFCAV--------------------------------
template <class T> 
void RFCAV_Pass(T x[6], double L, double VRF, double FRF, double PHASE0)
{
  DRIFT_Pass(x, L/2.);
  x[5] =x[5] + (VRF*GP.Q/GP.A/GP.energy)*sin(2.0*PI*FRF*x[z_]/3.0e8);
  DRIFT_Pass(x, L/2.);
}

class RFCAV: public Element
{
 public:
  RFCAV(string name, double l, double vrf, double frf, double phase0): Element(name)
    { 
      TYPE=string("RFCAV");
      VRF= vrf; 
      FRF= frf;
      PHASE0=phase0;
      L=l;
    }
  void   SetP(const char *name, double value) 
    {
     if (strcmp( name, "VRF" ) == 0) 
	{
	  VRF=value;
	}
     else if (strcmp( name, "FRF" ) == 0) 
	{
	  FRF=value;
	}
     else if (strcmp( name, "PHASE0" ) == 0) 
	{
	  PHASE0=value;
	}
     else 
       {
	  cout<<"RFCAV does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }
  double GetP(const char *name) 
    {
     if (strcmp( name, "VRF" ) == 0) 
	{
	  return VRF;
	}
     else if (strcmp( name, "FRF" ) == 0) 
	{
	  return FRF;
	}
     else if (strcmp( name, "PHASE0" ) == 0) 
	{
	  return PHASE0;
	}
      else 
	{
	  cout<<"RFCAV does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT);  RFCAV_Pass(x, L, VRF,FRF,PHASE0); LtoG(x,DX,DY,DT); 
  }
  void   DAPass(tps x[6]){
    if(GP.twiss_6d == 1 ) {
      GtoL(x,DX,DY,DT);  RFCAV_Pass(x, L, VRF,FRF,PHASE0); LtoG(x,DX,DY,DT); 
    }
    else{
      DRIFT_Pass(x, L);
    }
  }
 private:
  double VRF, FRF, PHASE0;
};

//----------------------BEAMBEAM------------------------------
void cerrf( double xx, double yy, double & wx, double & wy )
{
  int n,nc,nu ;
  double x,y,q,h,xl,xh,yh,tx,ty,tn,sx,sy,saux;
  double rx[33],ry[33];
  double cc=1.12837916709551, zero=0, one=1, two=2, half=.5, 
         xlim=5.33, ylim=4.29, fac1=3.2, fac2=23, fac3=21;

  x = abs(xx);
  y = abs(yy);
  if (y < ylim  and  x  < xlim) 
    {
      q  = (one - y / ylim) * sqrt(one - (x/xlim)*(x/xlim) );
      h  = one / (fac1 * q);
      nc = 7 + int(fac2*q);
      xl = pow(h,(1 - nc));
      xh = y + half/h;
      yh = x;
      nu = 10 + int(fac3*q);
      rx[nu+1] = zero;
      ry[nu+1] = zero;
      for ( n=nu; n>=1; n--) {
	tx = xh + n * rx[n+1];
	ty = yh - n * ry[n+1];
	tn = tx*tx + ty*ty;
	rx[n] = half * tx / tn;
	ry[n] = half * ty / tn;
      }
      sx = zero;
      sy = zero;
      for( n=nc;n>=1; n--) {
	saux = sx + xl;
	sx = rx[n] * saux - ry[n] * sy;
	sy = rx[n] * sy + ry[n] * saux;
	xl = h * xl;
      }
      wx = cc * sx;
      wy = cc * sy;
    }
  else 
    {
      xh = y;
      yh = x;
      rx[1] = zero;
      ry[1] = zero;
      for ( n = 9; n>= 1; n--) {
	tx = xh + n * rx[1];
	ty = yh - n * ry[1];
	tn = tx*tx + ty*ty;
	rx[1] = half * tx / tn;
	ry[1] = half * ty / tn;
      }
      wx = cc * rx[1];
      wy = cc * ry[1];
    }
  if(yy < zero) 
    {
      wx =   two * exp(y*y-x*x) * cos(two*x*y) - wx;
      wy = - two * exp(y*y-x*x) * sin(two*x*y) - wy;
      if(xx >  zero) 	wy = -wy;
    }
  else
    {
      if(xx <  zero) wy = -wy;
    }	
}

void BB4D(double x, double y, double gamma, double N, double sigmax, double sigmay, double & Dpx, double & Dpy)
{
  double rp=1.534698e-18;
  double x1, y1, x2, y2, signx, signy;

  signx = ((x>=0.0)?(1):(-1));
  signy = ((y>=0.0)?(1):(-1));

  if (x==0. && y==0. )
    {
      Dpx=0.;
      Dpy=0.;
    }
  else
    {
      if ( abs(sigmax -sigmay)/sigmax < 1.0e-6 )
	{  
	  double r2= x*x+ y*y;
	  double temp1= 2*N*rp/gamma/r2;
	  double temp2= 1-exp(-r2/2/sigmax/sigmax);
	  Dpx=temp1 * x * temp2;
	  Dpy=temp1 * y * temp2;
	}
      else if ( sigmax > sigmay) 
	{
	  double temp1=(2.*N*rp/gamma)*sqrt(PI/2./(sigmax*sigmax-sigmay*sigmay));
	  double temp2=exp(-x*x/2./sigmax/sigmax - y*y/2./sigmay/sigmay);
	  double temp3=sqrt(2.*(sigmax*sigmax-sigmay*sigmay));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs(x)/temp3;
	  y1 = fabs(y)/temp3;
	  x2 = fabs(x)*sigmay/sigmax/temp3;
	  y2 = fabs(y)*sigmax/sigmay/temp3;
	  cerrf( x1, y1, w1_real, w1_imag);
	  cerrf( x2, y2, w2_real, w2_imag);
	  Dpy=temp1 *( w1_real-temp2* w2_real)*signy;
	  Dpx=temp1 *( w1_imag-temp2* w2_imag)*signx;  
	}
      else if ( sigmax < sigmay) 
	{
	  double temp1=(2*N*rp/gamma)*sqrt(PI/2/(sigmay*sigmay-sigmax*sigmax));      
	  double temp2=exp(-x*x/2/sigmax/sigmax - y*y/2/sigmay/sigmay);
	  double temp3=sqrt(2*(sigmay*sigmay-sigmax*sigmax));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs(x)/temp3;
	  y1 = fabs(y)/temp3;
	  x2 = fabs(x)*sigmay/sigmax/temp3;
	  y2 = fabs(y)*sigmax/sigmay/temp3;
	  cerrf( y1, x1, w1_real, w1_imag);
	  cerrf( y2, x2, w2_real, w2_imag);
	  Dpy=temp1 *( w1_imag-temp2* w2_imag)*signy; 
	  Dpx=temp1 *( w1_real-temp2* w2_real)*signx;
	}
    }
}

void BB6D(double x[6], double gamma, double Np, double sigma_l, int N_slice, 
                       double emitx_rms,  double betax_star, double alfx_star,
	               double emity_rms,  double betay_star, double alfy_star)
 // emitx_rms, emity_rms:  un-normalized rms emittance,  sigma=SQRT[ emitx_rms * betax ]  
{
  int i,j; 
  double rp= 1.534698e-18;
  double x0[6];
  double Nsigma=4, Np_slice[N_slice], z_star[N_slice];
  double S, gx_star, gy_star, betax, alfx, betay, alfy;
  double sigmax, sigmay, dsigmax2ds, dsigmay2ds, dUdsigmax2, dUdsigmay2; 
  double X, Y,  Dpx, Dpy, temp1, temp2, temp3, temp4;
  
  //----center of each slice of strong bunch, each slice has not same particle population.
  
  for(i=0;i<N_slice;i++)
    z_star[i]= -1.0*Nsigma + (Nsigma*2.0/N_slice)*(2*i+1)/2;
  
  if(N_slice == 11 ) {
    Np_slice[0]=Np*  0.000500905  ;
    Np_slice[1]=Np*  0.0049242    ;
    Np_slice[2]=Np*  0.0290614    ;
    Np_slice[3]=Np*  0.103138     ;
    Np_slice[4]=Np*  0.220408     ;
    Np_slice[5]=Np*  0.28387      ;
    Np_slice[6]=Np*  0.220408     ;
    Np_slice[7]=Np*  0.103138     ;
    Np_slice[8]=Np*  0.0290614    ;
    Np_slice[9]=Np*  0.0049242    ;
    Np_slice[10]=Np* 0.000500905  ;  
  }
  else{
    for(i=0;i<N_slice;i++)
      Np_slice[i]=Np* (gsl_cdf_ugaussian_P( z_star[i]+ Nsigma*2.0/N_slice/2 )- gsl_cdf_ugaussian_P( z_star[i]- Nsigma*2.0/N_slice/2  ) );
  }  
  
  for(i=0; i<N_slice;i++ ) z_star[i]=z_star[i]*sigma_l;
  
  //for(i=0;i<N_slice;i++) cout << z_star[i]/sigma_l<<"   "<< Np_slice[i]/Np<<endl;

  //-----calculate the changes in x[6]
  for(i=0; i<N_slice;i++) {

    for (j=0;j<6;j++) x0[j]=x[j];
    S=(x0[4]-z_star[i])/2.0;

    gx_star= (1.0 + alfx_star * alfx_star )/ betax_star;
    gy_star= (1.0 + alfy_star * alfy_star )/ betay_star;
    betax= ( 1.0/gx_star + gx_star *(S-alfx_star/gx_star) *(S-alfx_star/gx_star) );
    alfx =(alfx_star -gx_star * S );
    betay= ( 1.0/gy_star + gy_star *(S-alfy_star/gy_star) *(S-alfy_star/gy_star) );
    alfy =(alfy_star -gy_star * S );

    sigmax  =  sqrt( emitx_rms * betax ) ;
    sigmay  =  sqrt( emity_rms * betay ) ;
    dsigmax2ds =-2.0*alfx* ( emitx_rms );
    dsigmay2ds =-2.0*alfy* ( emity_rms );
    
    //----calculate the kicks for each slice
    X=x0[0] + x0[1]*S;
    Y=x0[2] + x0[3]*S;
    BB4D(X, Y, gamma, Np_slice[i], sigmax, sigmay, Dpx, Dpy);
    temp1 = 1.0/2./(sigmax*sigmax-sigmay*sigmay) ;
    temp2 = X *  Dpx + Y *  Dpy;
    temp3 = 2.* Np_slice[i] *rp /gamma ;
    temp4 = exp ( - X * X / 2. / sigmax /sigmax -  Y * Y / 2. / sigmay /sigmay ) ;
    
    x[0]=x0[0] - S * Dpx * GP.bbscale ;
    x[1]=x0[1] + Dpx  * GP.bbscale;   
    x[2]=x0[2] - S * Dpy * GP.bbscale ;
    x[3]=x0[3] + Dpy * GP.bbscale;   
    x[4]=x0[4];
    
    if (sigmax == sigmay ) 
      {
      	x[5]=x0[5]+ 0.5 * Dpx * GP.bbscale * ( x0[1] + 0.5* Dpx * GP.bbscale ) 
                  + 0.5*  Dpy * GP.bbscale * ( x0[3] + 0.5* Dpy * GP.bbscale ) 
                  + (1.0/sigmax/sigmax ) * dsigmax2ds  * ( temp3 * GP.bbscale/ 4. ) * temp4;
      }
    else if (sigmax > sigmay ) 
      {
	dUdsigmax2 =  temp1 * GP.bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	dUdsigmay2 = -temp1 * GP.bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	x[5]= x0[5] + 0.5 * Dpx * GP.bbscale * ( x0[1] + 0.5* Dpx * GP.bbscale) 
                    + 0.5 * Dpy *GP.bbscale  * ( x0[3] + 0.5* Dpy *  GP.bbscale) 
                    - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
    else
      {
	dUdsigmax2 = -temp1 *  GP.bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	dUdsigmay2 =  temp1 *  GP.bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	x[5]= x0[5] + 0.5 * Dpx *GP.bbscale * ( x0[1] + 0.5* Dpx* GP.bbscale) 
                    + 0.5 * Dpy *GP.bbscale *  ( x0[3] + 0.5* Dpy * GP.bbscale) 
                    - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
  }
}

template <class T> 
void BEAMBEAM_Pass(T x[6], int TREATMENT, double NP, double SIGMAL, int NSLICE, double EMITX,  double BETAX, double ALFAX, double EMITY, double BETAY, double ALFAY)
{
  if ( NP != 0. ) {
      if(int(TREATMENT) == 6){
	BB6D(x, GP.gamma, NP, SIGMAL, NSLICE, EMITX,  BETAX, ALFAX, EMITY, BETAY, ALFAY); 
      }
      else{
	double Dpx, Dpy;
	BB4D(x[0], x[2], GP.gamma, NP, sqrt( EMITX * BETAX), sqrt( EMITY * BETAY ), Dpx,  Dpy);
	x[1]=x[1]+Dpx * GP.bbscale ;
	x[3]=x[3]+Dpy * GP.bbscale;
      }
  } 
}

class BEAMBEAM: public Element
{
 public:
 BEAMBEAM(string name, int treatment, double np, double sigmal, int nslice, double emitx, double emity, 
	  double betax, double alfax, double betay, double alfay): Element(name)
 //     emitx, emity:  un-normalized rms emittance,  sigma=SQRT[ emitx * betax ]  
  { 
    TREATMENT=treatment;
    NP=np;
    SIGMAL=sigmal;
    NSLICE=nslice;
    EMITX=emitx;
    EMITY=emity;
    BETAX=betax;
    ALFAX=alfax;
    BETAY=betay;
    ALFAY=alfay;
    TYPE=string("BEAMBEAM");
    L=0.;
  }
  void   SetP(const char *name, double value) 
  {
    if (strcmp( name, "TREATMENT" ) == 0) 
      { 
	int temp;
	temp= int(value);
	if(temp == 4 || temp == 6 ) {
	  TREATMENT=temp;}
	else{
	  cout<<"Error: beambeam treatment only can be 4 or 6"<<endl;
	  exit(1);
	}
      }
    else if (strcmp( name, "NP" ) == 0) 
      {
	NP=value;
      }
    else if (strcmp( name, "SIGMAL" ) == 0) 
      {
	SIGMAL=value;
      }
    else if (strcmp( name, "NSLICE" ) == 0) 
      {
	NSLICE=int(value);
      }
    else if (strcmp( name, "EMITX" ) == 0) 
      {
	EMITX=value;
      }
    else if (strcmp( name, "EMITY" ) == 0) 
      {
	EMITY=value;
      }
    else if (strcmp( name, "BETAX" ) == 0) 
      {
	BETAX=value;
      }
    else if (strcmp( name, "ALFAX" ) == 0) 
      {
	ALFAX=value;
      }
    else if (strcmp( name, "BETAY" ) == 0) 
      {
	BETAY=value;
      }
     else if (strcmp( name, "ALFAY" ) == 0) 
       {
	 ALFAY=value;
       }
     else 
       {
	 cout<<"BEAMBEAM does not have  parameter of "<<name<<endl; 
	 exit(0); 
       }    
  }
  double GetP(const char *name) 
  {
    if (strcmp( name, "TREATMENT" ) == 0) 
      {
	return TREATMENT;
      }
    else if (strcmp( name, "NP" ) == 0) 
      {
	return NP;
      }
    else if (strcmp( name, "SIGMAL" ) == 0) 
      {
	return SIGMAL;
      }
    else if (strcmp( name, "NSLICE" ) == 0) 
      {
	return NSLICE;
      }
    else if (strcmp( name, "EMITX" ) == 0) 
      {
	return EMITX;
      }
    else if (strcmp( name, "EMITY" ) == 0) 
      {
	return EMITY;
      }
    else if (strcmp( name, "BETAX" ) == 0) 
      {
	return BETAX;
      }
    else if (strcmp( name, "ALFAX" ) == 0) 
      {
	return ALFAX;
      }
    else if (strcmp( name, "BETAY" ) == 0) 
      {
	return BETAY;
      }
    else if (strcmp( name, "ALFAY" ) == 0) 
      {
	return ALFAY;
      }
    else 
      {
	cout<<"BEAMBEAM does not have  parameter of "<<name<<endl; 
	exit(0); 
      }    
  }
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT); 
    BEAMBEAM_Pass(x, TREATMENT,NP,SIGMAL,NSLICE,EMITX,BETAX,ALFAX,EMITY,BETAY,ALFAY);
    LtoG(x,DX,DY,DT);  
  }
  void   DAPass(tps x[6]){
    double sigmax, sigmay;
    double ksi_x, ksi_y;
    double rp= 1.534698e-18;
    sigmax=  sqrt( EMITX * BETAX );
    sigmay=  sqrt( EMITY * BETAY );
    ksi_x=   2.0*NP*rp/sigmax/GP.gamma/(sigmax+sigmay);  
    ksi_y=   2.0*NP*rp/sigmay/GP.gamma/(sigmax+sigmay);  
    GtoL(x,DX,DY,DT); 
    x[1]=x[1]+ksi_x * x[0] * GP.bbscale ; 
    x[3]=x[3]+ksi_y * x[2] * GP.bbscale ;
    LtoG(x,DX,DY,DT); 
  }
 private:
  int NSLICE, TREATMENT;
  double  NP, EMITX, EMITY, SIGMAL, BETAX, ALFAX, BETAY, ALFAY;
};

//---------------------LRBB ( long-range BEAMBEAM )------------
template <class T> 
void LRBB_Pass(T x[6], double gamma, double NP,double SEPX, double SEPY,double SIGMAX, double SIGMAY, double KICKX0, double KICKY0)
//  here we assume the particle's coordinates from the strong beam center is ( SEPX+x[1], SEPY+x[3] )
{
  double dpx1, dpy1;
  BB4D(SEPX+x[1], SEPY+x[3], GP.gamma, NP, SIGMAX, SIGMAY, dpx1,  dpy1);
  x[1]=x[1]+(dpx1 - KICKX0);
  x[3]=x[3]+(dpy1 - KICKY0);
}

class LRBB: public Element
{
 public:
 LRBB(string name, double np, double sepx, double sepy, double sigmax, double sigmay): Element(name)
 //     sepx, sepy are the offset of the weak beam center w.r.t to the center strong beam  
  { 
    NP=np;
    SEPX=sepx;
    SEPY=sepy;
    SIGMAX=sigmax;
    SIGMAY=sigmay;
    TYPE=string("LRBB");
    L=0.;
    BB4D(SEPX, SEPY, GP.gamma, NP, SIGMAX, SIGMAY, KICKX0,  KICKY0);
  }
  void   SetP(const char *name, double value) 
  {
    if (strcmp( name, "NP" ) == 0) 
      {
	NP=value;
      }
    else if (strcmp( name, "SEPX" ) == 0) 
      {
	SEPX=value;
      }
    else if (strcmp( name, "SEPY" ) == 0) 
      {
	SEPY=value;
      }
    else if (strcmp( name, "SIGMAX" ) == 0) 
      {
	SIGMAX=value;
      }
    else if (strcmp( name, "SIGMAY" ) == 0) 
      {
	SIGMAY=value;
      }
    else 
      {
	cout<<"LRBB does not have parameter of "<<name<<endl; 
	exit(0); 
      } 
    BB4D(SEPX, SEPY, GP.gamma, NP, SIGMAX, SIGMAY, KICKX0,  KICKY0);
  }
  double GetP(const char *name) 
  {
    if (strcmp( name, "NP" ) == 0) 
      {
	return NP;
      }
    else if (strcmp( name, "SEPX" ) == 0) 
      {
	return SEPX;
      }
    else if (strcmp( name, "SEPY" ) == 0) 
      {
	return SEPY;
      }
    else if (strcmp( name, "SIGMAX" ) == 0) 
      {
	return SIGMAX;
      }
    else if (strcmp( name, "SIGMAY" ) == 0) 
      {
	return SIGMAY;
      }
    else if (strcmp( name, "KICKX0" ) == 0) 
      {
	return KICKX0;
      }
    else if (strcmp( name, "KICKY0" ) == 0) 
      {
	return KICKY0;
      }
    else 
      {
	cout<<"LRBB does not have  parameter of "<<name<<endl; 
	exit(0); 
      }    
  }
  void   Pass(double x[6]){
    GtoL(x,DX,DY,DT); LRBB_Pass(x, GP.gamma, NP,SEPX,SEPY,SIGMAX, SIGMAY, KICKX0, KICKY0);  LtoG(x,DX,DY,DT);  
  }
  void   DAPass(tps x[6]){
    //no effect to any parameters.
  }
 private:
  double  NP, SEPX, SEPY, SIGMAX, SIGMAY;
  double  KICKX0, KICKY0;
};

//----------------------E-LENS------------------------------
template <class T> 
void ELENS_Pass(T x[6], double Ne, double Le, double beta_e, int N_slice, 
		double sigmax, double sigmay)
{
  int i;
  double Le_slice=Le*1.0/ N_slice; 
  double Ne_slice=Ne*1.0/N_slice;
  double Dpx, Dpy;
  
  for(i=0;i<N_slice;i++){
    DRIFT_Pass(x, Le_slice/2.);
    BB4D(x[0], x[2], GP.gamma, Ne_slice, sigmax, sigmay, Dpx, Dpy);
    x[1]=x[1]- Dpx * ( 1.+ beta_e);
    x[3]=x[3]- Dpy * ( 1.+ beta_e);
    DRIFT_Pass(x,Le_slice/2.);
  }
}

void elens_pass_round_Gaussian_topoff(double x[6], double gamma, 
	   double Ne, double Le, double beta_e, int N_slice, 
	   double sigmax, double sigmay) 
// I assume sigmax=sigmay, top off [-sigmax*a, sigmax*a ], Ne from perfect Gaussian
{
  int i;
  double Ne_slice= Ne*1.0/ N_slice; 
  double Le_slice= Le*1.0/ N_slice; 
  double a=0.4;                
  double scale;
  double r, rp=1.534698e-18;

  for(i=0;i<N_slice;i++){
    DRIFT_Pass(x, Le_slice/2.);
    r=  sqrt(x[0]*x[0]+x[2]*x[2]);
    if(  r <=  a * sigmax  ) 
      {
	scale = -2 * Ne_slice *(  exp(-a*a/2) * (r/sigmax) * (r/sigmax) /2  ) * rp  / gamma ;
	x[1]=x[1] + scale * x[0]/ r/ r ;
	x[3]=x[3] + scale * x[2]/ r/ r;
      }
    else
      {
	scale = -2 * Ne_slice * ( 1 -exp(-r*r/2/sigmax/sigmax) -  (1 -exp(-a*a/2) ) + exp(-a*a/2) * a * a /2  ) *  rp  / gamma ;
	x[1]=x[1] +  scale * x[0]/ r/r ;
	x[3]=x[3] +  scale * x[2]/ r/r ;
      }
    DRIFT_Pass(x,Le_slice/2.);
  }
}

void elens_pass_round_Gaussian_truncated(double x[6], double gamma, 
	   double Ne, double Le, double beta_e, int N_slice, 
	   double sigmax, double sigmay) 
// I assume sigmax=sigmay, Gaussian tail cut off from Nc*sigmax
{
  int i;
  double Le_slice=Le*1.0/ N_slice; 
  double Ne_slice=Ne*1.0/N_slice;
  double Dpx, Dpy;
  double Nc = 100, Nsigma;   
  
  for(i=0;i<N_slice;i++){
    DRIFT_Pass(x, Le_slice/2.);
    Nsigma= sqrt(x[0]*x[0]+ x[2]*x[2]) /  sigmax ;
    if(  Nsigma <=  Nc  ) 
      {
	BB4D(x[0], x[2], gamma, Ne_slice, sigmax, sigmay, Dpx, Dpy);
      }
    else
      {
	BB4D(x[0]* Nc/Nsigma, x[2]* Nc/Nsigma, gamma, Ne_slice, sigmax, sigmay, Dpx, Dpy);      
	Dpx =  Dpx * Nc / Nsigma;
	Dpy =  Dpy * Nc / Nsigma;
      }
    x[1]=x[1]- Dpx * ( 1.+ beta_e);
    x[3]=x[3]- Dpy * ( 1.+ beta_e);
    DRIFT_Pass(x,Le_slice/2.);
  }
}

void elens_pass_round_uniform(double x[6], double gamma, 
	   double Ne, double Le, double beta_e, int N_slice, 
	   double sigmax, double sigmay) 
// I assume round uniform distribution in r<=a
{
  int i;
  double Le_slice=Le*1.0/ N_slice; 
  double Ne_slice=Ne*1.0/N_slice;
  double a = 0.31e-3 * 5 ;
  double scale, r2;
  double rp=1.534698e-18;
  
  for(i=0;i<N_slice;i++){
    DRIFT_Pass(x, Le_slice/2.);
    scale=-2*Ne_slice * rp * ( 1 + beta_e ) / GP.gamma;
    r2=(x[0]*x[0]+x[2]*x[2]);
    if(  sqrt(r2) < a  ) 
      {
	x[1]=x[1] + scale * x[0]/ a/a ;
	x[3]=x[3] + scale * x[2]/ a/a ;
      }
    else
      {
	x[1]=x[1] + scale * x[0]/ r2 ;
	x[3]=x[3] + scale * x[2]/ r2 ;
      }
    DRIFT_Pass(x,Le_slice/2.);
  }
}

class ELENS: public Element
{
 public:
  ELENS(string name, double l, double ne, int nslice, double betae, double sigmax, double sigmay): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("ELENS");
	  L=l;
	  NE=ne;
	  NSLICE=nslice;
	  BETAE=betae;
	  SIGMAX=sigmax;
	  SIGMAY=sigmay;
	}
      else
	{
	  cout<<"Error: check signs of L. "<<endl;
          exit(1); 
	}
    }
    void   SetP(const char *name, double value) 
    {
      if (strcmp( name, "NE" ) == 0) 
	{
	  NE=value;
	}
      else if (strcmp( name, "NSLICE" ) == 0) 
	{
	  NSLICE=int(value);
	}
      else if (strcmp( name, "BETAE" ) == 0) 
	{
	  BETAE=value;
	}
      else if (strcmp( name, "SIGMAX" ) == 0) 
	{
	  SIGMAX=value;
	}
      else if (strcmp( name, "SIGMAY" ) == 0) 
	{
	  SIGMAY=value;
	}
      else 
	{
	  cout<<"ELENS does not have  parameter of "<<name<<endl; 
          exit(0); 
	}    
    }
    double GetP(const char *name) 
    {
      if (strcmp( name, "NE" ) == 0) 
	{
	  return NE;
	}
     else if (strcmp( name, "NSLICE" ) == 0) 
	{
	  return NSLICE;
	}
     else if (strcmp( name, "BETAE" ) == 0) 
	{
	  return BETAE;
	}
     else if (strcmp( name, "SIGMAX" ) == 0) 
	{
	  return SIGMAX;
	}
     else if (strcmp( name, "SIGMAY" ) == 0) 
	{
	  return SIGMAY;
	}
      else 
	{
	  cout<<"ELENS does not have  parameter of "<<name<<endl; 
          exit(0); 
	}    
    }
    void   Pass(double x[6]){
	GtoL(x,DX,DY,DT);  ELENS_Pass(x, NE, L, BETAE, NSLICE, SIGMAX, SIGMAY); LtoG(x,DX,DY,DT); 
    }
    void   DAPass(tps x[6]){
      int i;
      double ksi_x, ksi_y;
      double rp= 1.534698e-18;
      ksi_x= -2.0*(1.0*NE/NSLICE)*rp*(1+BETAE)/SIGMAX/GP.gamma/(SIGMAX+SIGMAY);  
      ksi_y= -2.0*(1.0*NE/NSLICE)*rp*(1+BETAE)/SIGMAY/GP.gamma/(SIGMAX+SIGMAY); 
      GtoL(x,DX,DY,DT); 
      for(i=0;i<NSLICE;i++){
	DRIFT_Pass(x, L*1.0/NSLICE/2.0);
	x[1]=x[1]+ksi_x * x[0] ;
	x[3]=x[3]+ksi_y * x[2] ;
	DRIFT_Pass(x, L*1.0/NSLICE/2.0);
      }
      LtoG(x,DX,DY,DT); 
    }
 private:
    double NE, BETAE, SIGMAX, SIGMAY;
    int    NSLICE;
};


//=========================================
//
//            Line  class
//
//=========================================
class Line
{
 public:
  Line()
    {
      Ncell=0;
      Length=0;
      Tune1=0.;
      Tune2=0.;
      Tune3=0.;
      Chromx1=0.;
      Chromy1=0.;
      Chromx2=0.;
      Chromy2=0.;
      Chromx3=0.;
      Chromy3=0.;
    }
  void Update()
  {
    int i;
    double spointer=0.;
    for(i=0;i<Cell.size();i++){
      spointer=Cell[i]->L+spointer;
      Cell[i]->S=spointer;
    }
    Ncell=Cell.size();  
    Length=spointer;
    frev0 = ( light_speed * GP.beta ) / Length ;
    frf   =  frev0  * GP.harm;
  }
  void Append(Element * x)
  {
    Cell.push_back( x );
    Update();
  }
  void Delete(int i)
  {
    vector<Element *>::iterator start;
    start=Cell.begin()+i;
    Cell.erase(start);
    Update();
  }
  void Insert( int i, Element * temp)
  {
    vector<Element *>::iterator start;
    start=Cell.begin()+i;
    Cell.insert(start,1,temp);
    Update();
  }  
  void Empty()
  {
    vector<Element *>::iterator start;
    vector<Element *>::iterator end;      
    start=Cell.begin();
    end=Cell.end();      
    Cell.erase(start,end);
    Update();
  }
  vector  <Element *> Cell;
  double Length;            //---length of reference orbit
  double frev0;             //---decided by Length
  double frf ;              //---rf frequency
  double Vrf_tot;           //---total rf voltage
  double Orbit_Length ;     //---length of closed orbit length.
  long   Ncell;             //---number of elements 
  double Tune1, Tune2, Tune3;  //----tunes
  double Chromx1, Chromy1, Chromx2, Chromy2, Chromx3, Chromy3;  //   chromaticities
  double Alfa0, Alfa1, Alfa2, Gammat, Slip;  //  transition parameters
  double Qs;                //  longitudinal tune
  double Bucket_length;     //  in unit of ns 
  double Bucket_height;     //  (dp/p0)_max
  double Bucket_area;       //  in ( phi_rf, dE/(hw_rev) ) phase space per nucleon
  double Bunch_length;      //  in unit of m, 6*sigma_l, in unit of ns
  double Bunch_area;        //  in ( phi_rf, dE/(hw_rev) ) phase space per nucleon
  double Bunch_height;      //  (dp/p0)_max for the given bunch area
};

void Line_Rewind(Line & linename,  int k )
{
  int i;
  Element * new_element;
  Line  temp_line;
  
  for(i=k;i<linename.Ncell;i++) {
    new_element = linename.Cell[i];
    temp_line.Append( new_element );
  }
  for(i=0; i<k;i++) {
    new_element = linename.Cell[i];
    temp_line.Append( new_element );
  }
  linename=temp_line;
}

void Line_Invert(Line & linename)
{
  int i;
  Element * new_element;
  Line  temp_line;
  
  for(i=0;i<linename.Ncell;i++) {
    new_element = linename.Cell[linename.Ncell-1 - i];
    temp_line.Append( new_element );
  }
  linename=temp_line;
}

void Line_Repeat(Line & linename1, Line & linename2, int n)
{
  int i,j;
  Element * new_element;
  
  for(j=0;j<n;j++){
    for(i=0;i<linename2.Ncell;i++) {
      new_element = linename2.Cell[i];
      linename1.Append( new_element );
    }
  }
}

void Line_Connect(Line & linename1, Line & linename2, Line & linename3)
{
  int i;
  Element * new_element;

  for(i=0;i<linename2.Ncell;i++) {
    new_element = linename2.Cell[i];
    linename1.Append( new_element );
  }
  for(i=0;i<linename3.Ncell;i++) {
    new_element = linename3.Cell[i];
    linename1.Append( new_element );
  }
}

int Get_Index(Line & linename, const char * name, int k)
{
  int i;
  int count=0;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->NAME==string(name)){
      count++;
      if(count == k ){
	return i;
      }
    }
  }
  if(count == 0) {
    cout<<"Element not found."<<endl;
    return 0;
  }
  return 0;
}

void Split_Quad(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length,k1l;
  Element * new_element;
  char  index[125];
  
  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}

  if(linename.Cell[i]->TYPE==string("QUAD")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    k1l=  linename.Cell[i]->GetP("K1L");
    linename.Delete(i);
    for(j=0;j<m;j++) {
      sprintf(index,"_%d",j+1); 
      new_element= new QUAD(old_name,length/m, k1l/m);
      linename.Insert(i,new_element);
    }
  }
  else{
    cout<<" The i-th element is not QUAD." <<endl;
    exit(1);
  }
}

void Split_Sext(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length,k2l;
  Element * new_element;
  char  index[125];
  string new_name;
  
  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}
  
  if(linename.Cell[i]->TYPE==string("SEXT")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    k2l=  linename.Cell[i]->GetP("K2L");
    linename.Delete(i);    
    for(j=0;j<m;j++) {
      sprintf(index,"_%d",j+1); 
      new_element= new SEXT(old_name,length/m, k2l/m);
      linename.Insert(i,new_element);
    }
  }
  else{
    cout<<" The i-th element is not SEXT." <<endl;
    exit(1);
  }
}

void Split_Quad_Sext(Line & linename)
{
  int i;
  
  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->TYPE==string("SEXT") )
      Split_Sext(linename,i,int(linename.Cell[i]->L /0.05)+1  );
  }
  
  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->TYPE==string("QUAD") )
      Split_Quad(linename,i,int(linename.Cell[i]->L /0.05)+1  );
  }
}

void Split_Sbend(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length,angle, e1, e2;
  Element * new_element;
  char index[125];

  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}
  
  if(linename.Cell[i]->TYPE==string("SBEND")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    angle=  linename.Cell[i]->GetP("ANGLE");
    e1   =  linename.Cell[i]->GetP("E1");
    e2   =  linename.Cell[i]->GetP("E2");
    linename.Delete(i);    
    for(j=0;j<m;j++) {
      sprintf(index,"_%d",j+1); 
      new_element= new SBEND(old_name,length/m, angle/m, e1,e2);
      linename.Insert(i,new_element);
    }
  }
  else{
    cout<<" The i-th element is not SBEND ."<<endl;
    exit(1);
  }
}

void Split_Drift(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length;
  Element * new_element;
  char index[125];

  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}
  
  if(linename.Cell[i]->TYPE==string("DRIFT")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    linename.Delete(i);    
    for(j=0;j<m;j++) {
      sprintf(index,"_%d",j+1); 
      new_element= new DRIFT(old_name,length/m);
      linename.Insert(i,new_element);
    }
  }
  else{
    cout<<" The i-th element is not DRIFT."<<endl;
    exit(1);
  }
}

void Split_Mult(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length;
  double knl[11],knsl[11];
  Element * new_element;
  char index[125];

  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}
  
  if(linename.Cell[i]->TYPE==string("DRIFT")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    
    for(j=0;j<11;j++){
      char name1[125], name2[125];
      sprintf(name1, "K%dL",j);
      sprintf(name2, "K%dSL",j);
      knl[j]= linename.Cell[i]->GetP(name1)/m;
      knsl[j]=linename.Cell[i]->GetP(name2)/m;
    }

    linename.Delete(i);    
    for(j=0;j<m;j++) {
      sprintf(index,"_%d",j+1); 
      new_element= new MULT(old_name,length/m,knl, knsl);
      linename.Insert(i,new_element);
    }
  }
  else{
    cout<<" The i-th element is not MULT."<<endl;
    exit(1);
  }
}

//----function: concate adjacient DRIFT, doesn't change Twiss
void Concat_Drift(Line & linename)
{
  int i;
  double length;
  
  i=1;
  while(i < linename.Ncell) {
    if( linename.Cell[i-1]->TYPE == string("DRIFT")  and  linename.Cell[i]->TYPE == string("DRIFT") ){
      length=linename.Cell[i]->L;
      linename.Delete(i);  
      linename.Cell[i-1]->L=linename.Cell[i-1]->L + length;
      linename.Update();
    }
    else{
      i++;
    }
  }
}

//----get rid of zero strength element, doesn't change Twiss
void Clean_Up(Line & linename)
{ 
  int i;  
  double length;
  Element *temp_element;

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("SEXT") and linename.Cell[i]->GetP("K2L") ==0  ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);
      temp_element= new DRIFT("TEMPD",length);
      linename.Insert(i,temp_element);
    }
  }

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("OCT") and linename.Cell[i]->GetP("K3L") ==0  ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);
      temp_element= new DRIFT("TEMPD",length);
      linename.Insert(i,temp_element);
    }
  }

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("SOLEN") and linename.Cell[i]->GetP("KS") ==0  ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);
      temp_element= new DRIFT("TEMPD",length);
      linename.Insert(i,temp_element);
    }
  }

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("RFCAV") and linename.Cell[i]->GetP("VRF") ==0  ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);
      temp_element= new DRIFT("TEMPD",length);
      linename.Insert(i,temp_element);
    }
  }
  
  for(i=0;i<linename.Ncell;i++)
    if( linename.Cell[i]->TYPE==string("KICKER") and linename.Cell[i]->GetP("HKICK") ==0.  and  linename.Cell[i]->GetP("VKICK") ==0.  ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);  
      if (length != 0. ) {
	temp_element=new DRIFT("TEMPD", length);
	linename.Insert(i, temp_element);
      }
    }
  
  for(i=0;i<linename.Ncell;i++)
    if( linename.Cell[i]->TYPE==string("HKICKER") and linename.Cell[i]->GetP("HKICK") ==0. ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);  
      if (length != 0. ) {
	temp_element=new DRIFT("TEMPD", length);
	linename.Insert(i, temp_element);
      }
    }
  
  for(i=0;i<linename.Ncell;i++)
    if( linename.Cell[i]->TYPE==string("VKICKER") and linename.Cell[i]->GetP("VKICK") ==0. ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);  
      if (length != 0. ) {
	temp_element=new DRIFT("TEMPD", length);
	linename.Insert(i, temp_element);
      }
    }

  for(i=0;i<linename.Ncell;i++)
    if( linename.Cell[i]->TYPE==string("MARKER")    or linename.Cell[i]->TYPE==string("BPM") 
	or linename.Cell[i]->TYPE==string("HBPM")   or linename.Cell[i]->TYPE==string("VBPM") 
	or linename.Cell[i]->TYPE==string("HBPM")   or linename.Cell[i]->TYPE==string("VBPM") ){
      length= linename.Cell[i]->L;
      linename.Delete(i);  
      temp_element=new DRIFT("TEMPD", length);
      linename.Insert(i, temp_element);
    }
  
  for(i=0;i<linename.Ncell;i++)
    if( linename.Cell[i]->TYPE==string("MULT") ) {
      length= linename.Cell[i]->L;
      int m;
      double kl=0.;
      double  knl[11], knsl[11];
      for(m=0;m<11;m++){
        char name1[125], name2[125];
        sprintf(name1, "K%dL",m);
        sprintf(name2, "K%dSL",m);
        knl[m] =linename.Cell[i]->GetP(name1);
        knsl[m]=linename.Cell[i]->GetP(name2);
      }
      for(m=0;m<11;m++)  kl=kl+abs(knl[m])+abs(knsl[m]);
      
      if(kl==0.){
	linename.Delete(i);  
	temp_element=new DRIFT("TEMPD",length);
	linename.Insert(i, temp_element);
      }
    }
  Concat_Drift(linename);
}

//----prepare for fast tarcking: will change Twiss
void Make_Thin(Line & linename) 
{
  int i;
  double length;
  Element * temp_element;

  //---convert Track_Fast() not-acceptable elements into drift
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE !=string("DRIFT")  &   linename.Cell[i]->TYPE !=string("SBEND")   &
	linename.Cell[i]->TYPE !=string("QUAD")   &   linename.Cell[i]->TYPE !=string("SEXT") &
	linename.Cell[i]->TYPE !=string("MULT")   &   linename.Cell[i]->TYPE !=string("BEAMBEAM")  &
        linename.Cell[i]->TYPE !=string("ELENS")  &   linename.Cell[i]->TYPE !=string("RFCAV")  &
	linename.Cell[i]->TYPE !=string("COOLING") &  linename.Cell[i]->TYPE !=string("ACMULT") &
	linename.Cell[i]->TYPE !=string("MATRIX")  &  linename.Cell[i]->TYPE !=string("KICKER") & 
        linename.Cell[i]->TYPE !=string("HKICKER") &  linename.Cell[i]->TYPE !=string("VKICKER") & 
        linename.Cell[i]->TYPE !=string("DIFFUSE") &  linename.Cell[i]->TYPE !=string("GMULT")  &
        linename.Cell[i]->TYPE !=string("LRBB") )  {
      length=linename.Cell[i]->L;
      temp_element= new  DRIFT("TEMPD",length);
      linename.Delete(i);
      linename.Insert(i,temp_element);
    }
  }

  //-----further make thin of some elements  among  Track_Fast() acceptable elements
  i=0;
  while(i < linename.Ncell-1 ) {
    if( linename.Cell[i]->TYPE ==string("SEXT")    || linename.Cell[i]->TYPE ==string("MULT")   ||
	linename.Cell[i]->TYPE ==string("HKICKER") || linename.Cell[i]->TYPE ==string("VKICKER")||
	linename.Cell[i]->TYPE ==string("RFCAV")   || linename.Cell[i]->TYPE ==string("ACMULT") ||
	linename.Cell[i]->TYPE ==string("COOLING") || linename.Cell[i]->TYPE ==string("DIFFUSE")||
        linename.Cell[i]->TYPE ==string("MATRIX")  || linename.Cell[i]->TYPE ==string("KICKER")  ) 
      {
	length=linename.Cell[i]->L;
	if(length != 0. ) {
	  temp_element=new DRIFT("TEMPD", length/2.0);
	  linename.Insert(i, temp_element);
	  linename.Cell[i+1]->L =0.;
	  temp_element=new DRIFT("TEMPD", length/2.0);
	  linename.Insert(i+2, temp_element);
	}
      }
    i++;  
  }

  //---even reduce the integration steps for SBEND  and QUAD to minimum 3
  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("SBEND") ) {
      if( linename.Cell[i]->GetP("Nint") > 3)  linename.Cell[i]->SetP("Nint",3); 
    }
  }
  
  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("QUAD") ) {
      if( linename.Cell[i]->GetP("Nint") > 3 )  linename.Cell[i]->SetP("Nint",3); 
    }
  }

  //---concat all created drifts in above processes
  Concat_Drift(linename);
}

double Get_KL(Line & linename, const char * name, const char *  kl)
{
  int i;
  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->NAME==string(name) ) {
      return linename.Cell[i]->GetP(kl);
    }
  }  
  return 0.;
}

void Set_KL(Line & linename, const char * name, const char *  kl, double strength)
{
  int i;
  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->NAME==string(name) )  linename.Cell[i]->SetP(kl, strength);
  }  
}

void Set_dKL(Line & linename, const char * name, const char * kl, double dstrength)
{
  int i;
  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->NAME==string(name) ) {
      linename.Cell[i]->SetP(kl, linename.Cell[i]->GetP(kl) + dstrength);
    }
  } 
}

//=====================================
//
//     interfaces to MADX
//
//=====================================
void  Read_MADXLattice(const char * filename, Line & linename)
{
  int i;
  fstream f1;
  string line;
  string temp_name;
  string temp_type;
  double temp_s;
  double temp_l;
  double temp_angle;
  double temp_e1;
  double temp_e2;
  double temp_tilt;
  double knl[11], knsl[11];
  double temp_sk;
  Element * temp_element;

  f1.open(filename,ios::in);
  if(!f1)
    {
      cout<<"error in opening the file: "<<filename<<endl;
      exit(0);
    }
  for(i=1;i<=47;i++) getline(f1,line,'\n');
  
  double s1=0;
  while( getline(f1,line,'\n')){
    istringstream ss(line);
    ss>>temp_name>>temp_type>>temp_s>>temp_l>>temp_angle>>temp_e1>>temp_e2>>temp_tilt
      >>knl[0]>>knsl[0]
      >>knl[1]>>knsl[1]
      >>knl[2]>>knsl[2]
      >>knl[3]>>knsl[3]
      >>knl[4]>>knsl[4]
      >>knl[5]>>knsl[5]
      >>knl[6]>>knsl[6]
      >>knl[7]>>knsl[7]
      >>knl[8]>>knsl[8]
      >>knl[9]>>knsl[9]
      >>knl[10]>>knsl[10]
      >>temp_sk;
    temp_name=string(temp_name,1,temp_name.size()-2 );     
    temp_type=string(temp_type,1,temp_type.size()-2 );
    
    s1=s1+temp_l;
    if(temp_type==string("DRIFT")){
      temp_element=new DRIFT(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("SBEND"))    {  // MADX also supports sbends with K1 and K2
      if(knl[1] == 0 and  knl[2] == 0  ) {
	temp_element= new SBEND(temp_name,temp_l,temp_angle, temp_e1, temp_e2);
	linename.Cell.push_back(temp_element); }
      else {
	temp_element=new GMULT(temp_name,temp_l, temp_angle, temp_e1, temp_e2, knl, knsl);
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("RBEND"))    {
      if(knl[1] == 0 ) {
	temp_element= new SBEND(temp_name,temp_l,temp_angle, temp_e1, temp_e2);
	linename.Cell.push_back(temp_element); }
      else {
	temp_element=new GMULT(temp_name,temp_l, temp_angle, temp_e1, temp_e2, knl, knsl);
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("QUADRUPOLE")){  // MAD8 quads only have K1, MADX also has K1S
      if ( knsl[1] != 0 ) {
	cout<<" Alert: QUAD with k1s detected. "<<endl;
      }
      temp_element= new QUAD(temp_name,temp_l, knl[1]);
      linename.Cell.push_back(temp_element); 
    }
    else if (temp_type==string("SEXTUPOLE")){   // MAD8 sexts only have K2, MADX also has K2S
      if ( knsl[2] != 0 ) {
	cout<<" Alert: SEXT with kss detected. "<<endl;
      }
      temp_element= new SEXT(temp_name,temp_l, knl[2]);
      linename.Cell.push_back(temp_element); 
    }
    else if (temp_type==string("OCTUPOLE")){    // MAD8 octs only have K3, MADX also has K2S
      if ( knsl[3] != 0 ) {
	cout<<" Alert: OCT with kss detected. "<<endl;
      }
      temp_element= new OCT(temp_name,temp_l, knl[3]);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("HKICKER")){
      temp_element= new HKICKER(temp_name,temp_l, -knl[0] );
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("VKICKER")){
      temp_element= new VKICKER(temp_name,temp_l, knsl[0]);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("KICKER")){
      temp_element= new KICKER(temp_name,temp_l, -knl[0], knsl[0]);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("MULTIPOLE")){  // MADX only support zero-length multipoles
      int m;
      double kl=0.;
      for(m=0;m<11;m++)  kl=kl+abs(knl[m])+abs(knsl[m]);
      if(kl==0.){
	temp_element=new DRIFT(temp_name,temp_l);
	linename.Cell.push_back(temp_element);	
      }
      else{
	temp_element=new MULT(temp_name,temp_l, knl, knsl);
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("HMONITOR")){
      temp_element= new HBPM(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("VMONITOR")){
      temp_element= new VBPM(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("MONITOR")){
      temp_element= new BPM(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("MARKER")){
      temp_element= new MARKER(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("SOLENOID")){  
      if(temp_sk !=0.0  and temp_l !=0. ){    // meaning  input for solenoid 
	temp_element= new SOLEN(temp_name,temp_l, temp_sk);
	linename.Cell.push_back(temp_element);
      }
      else{
	temp_element=new DRIFT(temp_name,temp_l);
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("RFCAVITY")){
      temp_element= new RFCAV(temp_name,temp_l, 0.0, 78250.42279*GP.harm, 0.);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("BEAMBEAM")){
      if(temp_l !=0.){ 
	cout<<"Beam-Beam element should have zero length!"<<endl;
        exit(1);  
      }
      temp_element= new BEAMBEAM(temp_name,6, 0., 0.4545, 11, 2.5e-06, 2.5e-06, 0.53, 0.0, 0.53, 0.0);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("ELENS")){
      temp_element= new ELENS(temp_name, temp_l, 0.0, 5, 0.14, .31e-03, .31e-03);
      linename.Cell.push_back(temp_element);
    }
    else
      {
	temp_element=new DRIFT(temp_name,temp_l);
	linename.Cell.push_back(temp_element);
      }
  }
  linename.Ncell=linename.Cell.size();
  linename.Update();
  f1.close();
}

void  Print_MADXLattice(const char * filename, Line & linename)
{
  int i, j, flag;
  fstream f2;
  char str[125];
  vector <string> namelist;
  
  f2.open(filename,ios::out);
  if(!f2)
    {
      cout<<"error in opening the file: "<<filename<<endl;
      exit(0);
    }
  
  for(i=0;i<linename.Ncell; i++){
    flag=0;
    for(j=0;j<namelist.size();j++) {
      if(linename.Cell[i]->NAME  == namelist [j] ) {
	flag=1;break;
      }
    }
    if(flag==0) {
      namelist.push_back(linename.Cell[i]->NAME);
      if(linename.Cell[i]->TYPE==string("DRIFT") ){
      	f2<<linename.Cell[i]->NAME<<" : DRIFT, L = "<<setprecision(10)<<linename.Cell[i]->L<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("QUAD") )  {
	f2<<linename.Cell[i]->NAME<<" : QUADRUPOLE, L = "<<setprecision(10)<<linename.Cell[i]->L<<", K1 =  "<<setprecision(10)<<linename.Cell[i]->GetP("K1L")/ linename.Cell[i]->L<<" ;"<<endl;}
      else if(linename.Cell[i]->TYPE==string("SEXT") ) {
	f2<<linename.Cell[i]->NAME<<" : SEXTUPOLE,  L = "<<setprecision(10)<<linename.Cell[i]->L<<", K2 =  "<<setprecision(10)<<linename.Cell[i]->GetP("K2L")/ linename.Cell[i]->L<<" ;"<<endl; }
      else if(linename.Cell[i]->TYPE==string("SBEND") ) { 
	f2<<linename.Cell[i]->NAME<<" : SBEND,      L = "<<setprecision(10)<<linename.Cell[i]->L<<", ANGLE =  "<<setprecision(10)<<linename.Cell[i]->GetP("ANGLE")<<" ;"<<endl; }
      else if(linename.Cell[i]->TYPE==string("RFCAV") ) {  
	f2<<linename.Cell[i]->NAME<<" : RFCAVITY,   L = "<<setprecision(10)<<linename.Cell[i]->L<<" ;"<<endl; }
      else if(linename.Cell[i]->TYPE==string("MULT") ) {
	f2<<linename.Cell[i]->NAME<<" : MULTIPOLE,  "<<endl;
	f2<<" knl:={ "<<linename.Cell[i]->GetP("K0L")<<","<<linename.Cell[i]->GetP("K1L")<<","<<linename.Cell[i]->GetP("K2L")<<","<<linename.Cell[i]->GetP("K3L")<<","<<linename.Cell[i]->GetP("K4L")<<","<<linename.Cell[i]->GetP("K5L")<<",";  
	f2<<linename.Cell[i]->GetP("K6L")<<","<<linename.Cell[i]->GetP("K7L")<<","<<linename.Cell[i]->GetP("K8L")<<","<<linename.Cell[i]->GetP("K9L")<<","<<linename.Cell[i]->GetP("K10L")<<"},"<<endl;
	f2<<" ksl:={ "<<linename.Cell[i]->GetP("K0SL")<<","<<linename.Cell[i]->GetP("K1SL")<<","<<linename.Cell[i]->GetP("K2SL")<<","<<linename.Cell[i]->GetP("K3SL")<<","<<linename.Cell[i]->GetP("K4SL")<<","<<linename.Cell[i]->GetP("K5SL")<<",";  
	f2<<linename.Cell[i]->GetP("K6SL")<<","<<linename.Cell[i]->GetP("K7SL")<<","<<linename.Cell[i]->GetP("K8SL")<<","<<linename.Cell[i]->GetP("K9SL")<<","<<linename.Cell[i]->GetP("K10SL")<<"};"<<endl; 
        if (linename.Cell[i]->L != 0. ) cout<<"Warning: "<< linename.Cell[i]->NAME<<"  "<<linename.Cell[i]->TYPE<<"   "<<linename.Cell[i]->L<<"  transferred to DRIFT ."<<endl;   }
      else {
	f2<<linename.Cell[i]->NAME<<" : DRIFT, L = "<<setprecision(10)<<linename.Cell[i]->L<<";"<<endl;
        cout<<"Warning: "<< linename.Cell[i]->NAME<<"  "<<linename.Cell[i]->TYPE<<"   "<<setprecision(10)<<linename.Cell[i]->L<<"  transferred to DRIFT ."<<endl;
      }
    }
  }
  
  int istart=0;
  int inumber10=0; 
  int nline=linename.Ncell / 1000;

  for(i=0;i<nline+1;i++) {
    sprintf(str,"LIN%d",i);
    f2<<str<<" : LINE=("<<endl;
    inumber10=0;

    do{
      f2<<linename.Cell[istart]->NAME<<",";
      istart++;
      inumber10++;
      if (inumber10 ==10){ f2<< endl; inumber10=0;}
    } while ( istart < 1000*(i+1)  &&  istart < linename.Ncell-1  );
    
    if(i < nline ){
      f2<<linename.Cell[istart]->NAME<<" ); "<<endl;
      istart++;
    }
    else{
      f2<<linename.Cell[linename.Ncell-1]->NAME<<" ); "<<endl;
    }
  }

  f2<<"rhic:  LINE = ( ";
  for(i=0;i<nline;i++) {
    sprintf(str,"LIN%d",i);
    f2<<str<<",";
  }
  sprintf(str,"LIN%d",i);
  f2<<str<<" );"<<endl;

  f2<<"beam, mass:=0.93827, charge:=1, gamma:=268.2, exn:=20.0e-06, eyn:=20.0e-06, sige:=0.001;"<<endl;
  f2<<"use, period=rhic;"<<endl;
  f2<<"select, flag=twiss,  column=NAME, KEYWORD,S,BETX, BETY, DX, DY, MUX, MUY;"<<endl;
  f2<<"twiss,table=twiss,file=twiss.table;"<<endl;
  f2<<"stop;"<<endl;
  f2.close();
}

//================================================================
//
//    PHYSICS  calculation
//
//================================================================

//---5D orbit and Twiss 
void Cal_Orbit_Num(Line & linename, double deltap)
{
  int i,j,k,iter=0;
  double x0[6],x01[6], x1[6], dx[4];
  double d=1.0e-09;
  double mat[4][4];
  double chi;
  double codeps=1e-12;
  int Max_iter=20;

  for (i=0;i<6;i++) x0[i]=0.;
  do{
    
    iter++;
    x0[4]=0;
    x0[5]=deltap;
    
    for(j=0;j<6;j++) x01[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x01);
    
    
    for(k=0;k<4;k++){
      for(j=0;j<6;j++) x1[j]=x0[j];
      x1[k]=x1[k]+d;
      for(j=0;j<linename.Ncell;j++) linename.Cell[j]->Pass(x1);
      for(j=0;j<4;j++) mat[j][k]=(x1[j]-x01[j])/d;
    }
    
    for(j=0;j<6;j++) x1[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x1);
    for(i=0;i<4;i++) {
      dx[i]=x0[i]-x1[i];
    }
    chi=0;
    for(i=0;i<4;i++) chi+=dx[i]*dx[i];
    chi=sqrt(chi/4.);
    
    if( chi > codeps ){
      for(i=0;i<4;i++) mat[i][i]=mat[i][i]-1.0000001;
      mat_inv(&mat[0][0],4); 
      for(i=0;i<4;i++){
	for(j=0;j<4;j++) x0[i]+=mat[i][j]*dx[j];  
      } 
    } 
    
  }while (chi>codeps && iter <Max_iter );
  
  if(iter == Max_iter-1 ) 
    {
      cout<<"Failed to find COD."<<endl;
      exit(1);
    }
  else
    {
      x0[4]=0.000;
      x0[5]=deltap;
      for(j=0;j<linename.Ncell; j++) {
	linename.Cell[j]->Pass(x0);
	for(i=0;i<6;i++) linename.Cell[j]->X[i]=x0[i];
      } 
    }
}

void Cal_Orbit_DA(Line & linename, double deltap)
{
  int i,j,iter=0;
  double x0[6], x1[6], dx[6];
  tps tps1[6];
  linmap m1;
  double mat[6][6];
  double chi;
  double codeps=1e-12;
  int Max_iter=20;
  
  for (i=0;i<6;i++) x0[i]=0.;

  do{
  
    iter++;

    x0[4]=0;
    x0[5]=deltap;
    m1.identity();
    m1=m1+x0;
    
    for(i=0;i<6;i++) tps1[i]=m1[i];
    for(j=0;j<linename.Ncell; j++) linename.Cell[j]->DAPass(tps1);
    for(i=0;i<6;i++) m1[i]=tps1[i];
    Getmat(m1, &mat[0][0]);
    Getpos(m1, x1);
    
    for(i=0;i<6;i++) {
      dx[i]=x0[i]-x1[i];
    }

    chi=0;
    for(i=0;i<4;i++) chi+=dx[i]*dx[i];
    chi=sqrt(chi/4.);
 
    if( chi > codeps ){
      dx[4]=0.;
      dx[5]=0.;
      for(i=0;i<6;i++) mat[i][i]=mat[i][i]-1.0000001;
      
      mat_inv(&mat[0][0],6); 
      for(i=0;i<6;i++){
	for(j=0;j<6;j++) x0[i]+=mat[i][j]*dx[j];  
      } 
    } 
    
  }while (chi>codeps && iter <Max_iter );
  
  if(iter == Max_iter-1 ) 
    {
      cout<<"Failed to find COD."<<endl;
      exit(1);
    }
  else
    {
      x0[4]=0.000;
      x0[5]=deltap;
      for(j=0;j<linename.Ncell; j++) {
	linename.Cell[j]->Pass(x0);
	for(i=0;i<6;i++) linename.Cell[j]->X[i]=x0[i];
      } 
    }
}

void Cal_OneTurnMap(Line & linename, double deltap)
{
  int i,j;
  double x[6];
  tps tps1[6];
  linmap m1;
  double m66[36], b66[36];
  int flag;
  double u[6], v[6];

  for(i=0;i<6;i++) x[i]=linename.Cell[linename.Ncell-1]->X[i];
  x[4]=0.000;
  x[5]=deltap; 
  m1.identity();
  m1=m1+x;
  
  for(i=0;i<6;i++) tps1[i]=m1[i];
  for(j=0;j<linename.Ncell; j++) linename.Cell[j]->DAPass(tps1);
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, m66); 

  for(i=0;i<36;i++) linename.Cell[linename.Ncell-1]->M[i]=m66[i];
  
  if( false ) {
    cout<<"One turn map M:"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<m66[i*6+j]<<"  ";
      cout<<endl;
    }
    for(i=0;i<36;i++) b66[i]=m66[i];
    cout<<"Det of M = "<<mat_det(b66,6)<<endl;
  }
  
  mat_change_hessenberg(m66, 6);
  flag=mat_root_hessenberg(m66,6,u,v,1.0e-10,60);
  if(flag < 0){
    cout<<"Unstable motion. "<<endl;
    exit(1);
  }
}

void Cal_ElementMap(Line & linename, double deltap)
{
  int i,j,k, istart;
  double x[6];
  tps tps1[6];
  linmap m1;
  double t66[36];

  for(k=0; k<linename.Ncell; k++) 
    {
      istart=k-1;
      if (istart <0) istart=linename.Ncell-1;

      for(i=0;i<6;i++) x[i]=linename.Cell[istart]->X[i];
      x[4]=0.000;
      x[5]=deltap; 
      m1.identity();
      m1=m1+x;
      
      for(i=0;i<6;i++) tps1[i]=m1[i];
      linename.Cell[k]->DAPass(tps1);
      for(i=0;i<6;i++) m1[i]=tps1[i];
      Getmat(m1, t66);
      
      for(j=0;j<36;j++) linename.Cell[k]->T[j]=t66[j]; 
    }
}

void Cal_SectionMap(Line & linename, int i1, int i2, double deltap, double t66[36] )
{
  int i;
  double x[6];
  tps tps1[6];
  linmap m1;
  
  if (i1== 0 ) { 
    for(i=0;i<6;i++) x[i] = linename.Cell[linename.Ncell-1]->X[i]; }
  else {
    for(i=0;i<6;i++) x[i] = linename.Cell[i1-1]->X[i]; }
  x[4]=0.000;
  x[5]=deltap; 
  m1.identity();
  m1=m1+x;
  
  for(i=0;i<6;i++) tps1[i]=m1[i];
  for(i=i1;i<=i2;i++)  linename.Cell[i]->DAPass(tps1);
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, t66);
}

//------original 4-D approach
void Cal_A0(Line & linename, double deltap)
{
  int i,j;
  double m66[6][6], m44[4][4];

  double wr[4], wi[4], vr[4][4], vi[4][4], temp_wr[4], temp_wi[4], temp_vr[4][4], temp_vi[4][4];  
  double temp1, temp2;
  int flag;
  double theta1, theta2, tempx, tempy, x1,y1;
  double a44[4][4], b44[4][4];
  double det, scale;
  
  //---solve eigen problem to construct A[4*4]
  for(i=0;i<6;i++)
   for(j=0;j<6;j++) m66[i][j]=linename.Cell[linename.Ncell-1]->M[i*6+j];

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      m44[i][j]=m66[i][j];
  EigenSolver(m44, wr, wi, vr, vi);
  if(false){
    cout<<"Eigen values and eigen vectors: "<<endl;
    for(i=0;i<4;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<4;j++){
        cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl; 
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

  //---link v and v* together
  temp1=wr[0];
  temp2=wi[0];
  for(i=0; i<4;i++)     {
    if (wr[i]==temp1 && wi[i]== -temp2 ) { 
      temp_wr[0]= wr[0];
      temp_wi[0]= wi[0];
      for(j=0;j<4;j++) temp_vr[0][j]= vr[0][j];
      for(j=0;j<4;j++) temp_vi[0][j]= vi[0][j];
      temp_wr[1]= wr[i];
      temp_wi[1]= wi[i];       
      for(j=0;j<4;j++) temp_vr[1][j]= vr[1][j];
      for(j=0;j<4;j++) temp_vi[1][j]= vi[1][j];
    }
  }
  
  flag=0;
  for(i=0; i<4;i++)     {
    if (wr[i] != temp1 ) { 
      if(flag ==0 ){
	temp_wr[2]= wr[i];
	temp_wi[2]= wi[i];
	for(j=0;j<4;j++) temp_vr[2][j]= vr[i][j];
	for(j=0;j<4;j++) temp_vi[2][j]= vi[i][j];
	flag=1;
      }       
      if(flag==1){
	temp_wr[3]= wr[i];
	temp_wi[3]= wi[i];       
	for(j=0;j<4;j++) temp_vr[3][j]= vr[i][j];
	for(j=0;j<4;j++) temp_vi[3][j]= vi[i][j];
      }
    }
  }

  //---link eigen vetcors to  H/V planes 
  temp1=abs( temp_vr[0][0] ) + abs( temp_vi[0][0] ) + abs( temp_vr[0][1] ) + abs( temp_vi[0][1] );
  temp2=abs( temp_vr[2][0] ) + abs( temp_vi[2][0] ) + abs( temp_vr[2][1] ) + abs( temp_vi[2][1] );
  
  if (temp1 > temp2 ) {
    for(i=0;i<4;i++){
      wr[i]= temp_wr[i];
      wi[i]= temp_wi[i];       
      for(j=0;j<4;j++) vr[i][j]= temp_vr[i][j];
      for(j=0;j<4;j++) vi[i][j]= temp_vi[i][j];        
    }
  }
  else
    {
      for(i=0;i<4;i++){
	wr[3-i]= temp_wr[i];
	wi[3-i]= temp_wi[i];       
	for(j=0;j<4;j++) vr[3-i][j]= temp_vr[i][j];
	for(j=0;j<4;j++) vi[3-i][j]= temp_vi[i][j];        
      }
    }

  //----rotating to make A12=0, A34=0
  theta1=atan2(vi[0][0], vr[0][0]);
  x1= cos(-PI/2-theta1);
  y1=sin(-PI/2-theta1);
  for(i=0;i<4;i++){
    tempx=vr[0][i]*x1-vi[0][i]*y1;
    tempy=vi[0][i]*x1+vr[0][i]*y1;
    vr[0][i]=tempx;
    vi[0][i]=tempy; }
  x1= cos(PI/2+theta1);
  y1=sin(PI/2+theta1);
  for(i=0;i<4;i++){
    tempx=vr[1][i]*x1-vi[1][i]*y1;
    tempy=vi[1][i]*x1+vr[1][i]*y1;
    vr[1][i]=tempx;
    vi[1][i]=tempy;}

  theta2=atan2(vi[2][2], vr[2][2]);
  x1= cos(-PI/2-theta2);
  y1=sin(-PI/2-theta2);
  for(i=0;i<4;i++){
    tempx=vr[2][i]*x1-vi[2][i]*y1;
    tempy=vi[2][i]*x1+vr[2][i]*y1;
    vr[2][i]=tempx;
    vi[2][i]=tempy; }
   x1  = cos(PI/2+theta2);
   y1=sin(PI/2+theta2);
  for(i=0;i<4;i++){
    tempx=vr[3][i]*x1-vi[3][i]*y1;
    tempy=vi[3][i]*x1+vr[3][i]*y1;
    vr[3][i]=tempx;
    vi[3][i]=tempy;}

  //----produce A
  for(i=0;i<4;i++) a44[i][0]=-(vi[0][i]-vi[1][i]);
  for(i=0;i<4;i++) a44[i][1]=vr[0][i]+vr[1][i];
  for(i=0;i<4;i++) a44[i][2]=-(vi[2][i]-vi[3][i]);
  for(i=0;i<4;i++) a44[i][3]=vr[2][i]+vr[3][i];

  if(false) {
    cout<<"A matrix :"<<endl;
    for (i=0;i<4;i++) {
      for(j=0;j<4;j++) cout<<setw(12)<<a44[i][j]<<"  ";
      cout<<endl;
    }
  }

  if(a44[1][1] < 0. ) {
    //cout<<"... Eigenmode I: reverted sign. "<<endl;
    for(i=0;i<4;i++) a44[i][1]=- a44[i][1];
  }
  if(a44[3][3] < 0. ) {
    //cout<<"... Eigenmode II: reverted sign. "<<endl;
    for(i=0;i<4;i++) a44[i][3]=- a44[i][3];
  }
  
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) b44[i][j]=a44[i][j];
  }
  det= mat_det( &b44[0][0], 4);

  scale=sqrt(sqrt(abs(det)));
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) a44[i][j]=a44[i][j]/scale;
  }
  
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) 
      linename.Cell[linename.Ncell-1]->A[i*6+j]=a44[i][j];
}

//-----adapted from 6-D general approach
void Cal_A(Line & linename, double deltap)
{
  
  int i,j;
  double m44[4][4];

  double wr[4], wi[4], vr[4][4], vi[4][4], temp_wr, temp_wi, temp_vr[4], temp_vi[4];  
  double temp1, temp2;
  double theta1, theta2, tempx, tempy, x1,y1;
  double a44[4][4], b44[4][4];
  double det, scale;
  
  //---solve eigen problem to construct A[4*4]
  for(i=0;i<4;i++)
   for(j=0;j<4;j++) m44[i][j]=linename.Cell[linename.Ncell-1]->M[i*6+j];
  EigenSolver(m44, wr, wi, vr, vi);

  //---link eigen vetcors to  H/V planes 
  temp1=abs(vr[0][0] ) + abs( vi[0][0] ) + abs( vr[0][1] ) + abs( vi[0][1] );
  temp2=abs(vr[2][0] ) + abs( vi[2][0] ) + abs( vr[2][1] ) + abs( vi[2][1] );
  if (temp2 > temp1 ) {
      temp_wr    =  wr[0] ;
      temp_wi    =  wi[0] ;
      wr[0]      =  wr[2] ;
      wi[0]      =  wi[2] ;
      wr[2]      =  temp_wr;
      wi[2]      =  temp_wi;
      wr[1]      =  wr[0];
      wi[1]      = -wi[0]; 
      wr[3]      =  wr[2];
      wi[3]      = -wi[2]; 
    for(i=0;i<4;i++){
      temp_vr[i]    = vr[0][i];
      temp_vi[i]    = vi[0][i];
      vr[0][i]      = vr[2][i];
      vi[0][i]      = vi[2][i];
      vr[2][i]      = temp_vr[i];
      vi[2][i]      = temp_vi[i];
      vr[1][i]      = vr[0][i];
      vi[1][i]      =-vi[0][i]; 
      vr[3][i]      = vr[2][i];
      vi[3][i]      =-vi[2][i]; 
    }
  }

  //----rotating to make A12=0, A34=0
  theta1=atan2(vi[0][0], vr[0][0]);
  x1  = cos(-PI/2-theta1);
  y1  =sin(-PI/2-theta1);
  for(i=0;i<4;i++){
    tempx=vr[0][i]*x1-vi[0][i]*y1;
    tempy=vi[0][i]*x1+vr[0][i]*y1;
    vr[0][i]=tempx;
    vi[0][i]=tempy; }
  x1  = cos(PI/2+theta1);
  y1 =sin(PI/2+theta1);
  for(i=0;i<4;i++){
    tempx=vr[1][i]*x1-vi[1][i]*y1;
    tempy=vi[1][i]*x1+vr[1][i]*y1;
    vr[1][i]=tempx;
    vi[1][i]=tempy;}

  theta2=atan2(vi[2][2], vr[2][2]);
  x1  = cos(-PI/2-theta2);
  y1 =sin(-PI/2-theta2);
  for(i=0;i<4;i++){
    tempx=vr[2][i]*x1-vi[2][i]*y1;
    tempy=vi[2][i]*x1+vr[2][i]*y1;
    vr[2][i]=tempx;
    vi[2][i]=tempy; }
  x1  = cos(PI/2+theta2);
  y1=sin(PI/2+theta2);
  for(i=0;i<4;i++){
    tempx=vr[3][i]*x1-vi[3][i]*y1;
    tempy=vi[3][i]*x1+vr[3][i]*y1;
    vr[3][i]=tempx;
    vi[3][i]=tempy;}

  //----produce A
  for(i=0;i<4;i++) a44[i][0]=-(vi[0][i]-vi[1][i]);
  for(i=0;i<4;i++) a44[i][1]=vr[0][i]+vr[1][i];
  for(i=0;i<4;i++) a44[i][2]=-(vi[2][i]-vi[3][i]);
  for(i=0;i<4;i++) a44[i][3]=vr[2][i]+vr[3][i];

  if(true){
    if(a44[1][1] < 0. ) {
      for(i=0;i<4;i++) a44[i][1]=- a44[i][1];
      //cout<<"Mode I : negative beta reverted."<<endl; 
    }
    if(a44[3][3] < 0. ) {
      for(i=0;i<4;i++) a44[i][3]=- a44[i][3];
      //cout<<"Mode II : negative beta reverted."<<endl; 
    }
  }
  
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) b44[i][j]=a44[i][j];
  }
  det= mat_det( &b44[0][0], 4);

  scale=sqrt(sqrt(abs(det)));
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) a44[i][j]=a44[i][j]/scale;
  }
  
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) 
      linename.Cell[linename.Ncell-1]->A[i*6+j]=a44[i][j];

 if(false) {
    cout<<"A matrix :"<<endl;
    for (i=0;i<4;i++) {
      for(j=0;j<4;j++) cout<<setw(12)<<a44[i][j]<<"  ";
      cout<<endl;
    }
    
    for(i=0;i<4;i++) {
      for(j=0;j<4;j++) b44[i][j]=a44[i][j];
    }
    det= mat_det( &b44[0][0], 4);
    cout<<" Det of A  =  "<< det <<endl;
  }

}

void Trace_A(Line & linename, double deltap)
{
  int i,j,k;
  double T[4][4], A[4][4], G[4][4];
  double dphi1, dphi2, scale1, scale2;
  double mu1, mu2;
 
  for(i=0;i<4;i++)
    for(j=0;j<4;j++) A[i][j]= linename.Cell[linename.Ncell-1]->A[i*6+j];

  mu1=0.;
  mu2=0.;
  for(k=0;k<linename.Ncell;k++){
    
    for (i=0;i<4;i++){
      for(j=0;j<4;j++) 	T[i][j]=linename.Cell[k]->T[i*6+j];
    }
    
    mat_mult( &T[0][0], &A[0][0], &G[0][0], 4,4,4 );
    
    dphi1=atan2(G[0][1], G[0][0]);
    dphi2=atan2(G[2][3], G[2][2]);
    scale1=sqrt( G[0][0]*G[0][0] + G[0][1]*G[0][1]);
    scale2=sqrt( G[2][2]*G[2][2] + G[2][3]*G[2][3]);

    A[0][0]= ( G[0][0]*G[0][0]+G[0][1]*G[0][1] ) / scale1;
    A[0][1]= (-G[0][0]*G[0][1]+G[0][1]*G[0][0] ) / scale1;
    A[1][0]= ( G[1][0]*G[0][0]+G[1][1]*G[0][1] ) / scale1;
    A[1][1]= (-G[1][0]*G[0][1]+G[1][1]*G[0][0] ) / scale1;
                                                                          
    A[2][0]= ( G[2][0]*G[0][0]+G[2][1]*G[0][1] ) / scale1;
    A[2][1]= (-G[2][0]*G[0][1]+G[2][1]*G[0][0] ) / scale1;
    A[3][0]= ( G[3][0]*G[0][0]+G[3][1]*G[0][1] ) / scale1;
    A[3][1]= (-G[3][0]*G[0][1]+G[3][1]*G[0][0] ) / scale1; 
                                                                          
    A[0][2]= ( G[0][2]*G[2][2]+G[0][3]*G[2][3] ) / scale2;
    A[0][3]= (-G[0][2]*G[2][3]+G[0][3]*G[2][2] ) / scale2;
    A[1][2]= ( G[1][2]*G[2][2]+G[1][3]*G[2][3] ) / scale2;
    A[1][3]= (-G[1][2]*G[2][3]+G[1][3]*G[2][2] ) / scale2;
                                                                          
    A[2][2]= ( G[2][2]*G[2][2]+G[2][3]*G[2][3] ) / scale2;
    A[2][3]= (-G[2][2]*G[2][3]+G[2][3]*G[2][2] ) / scale2;
    A[3][2]= ( G[3][2]*G[2][2]+G[3][3]*G[2][3] ) / scale2;
    A[3][3]= (-G[3][2]*G[2][3]+G[3][3]*G[2][2] ) / scale2; 

    for(i=0;i<4;i++){
      for(j=0;j<4;j++) linename.Cell[k]->A[i*6+j]=A[i][j];
    }
    
    mu1=mu1+dphi1/ 2./PI;
    mu2=mu2+dphi2/ 2./PI;
    linename.Cell[k]->Mu1= mu1;
    linename.Cell[k]->Mu2= mu2;
  }

  linename.Tune1=  mu1;
  linename.Tune2=  mu2;
}

void Cal_Twiss(Line & linename, double deltap)
{
  int i,j,k;
  double temp;
  double A12[4], A22[4], C[4];

  Cal_Orbit_Num(linename, deltap);
  Cal_OneTurnMap(linename, deltap);
  Cal_A(linename, deltap);
  Cal_ElementMap(linename, deltap);
  Trace_A(linename, deltap);
  
  for(k=0; k<linename.Ncell;k++){
    linename.Cell[k]->Beta1= linename.Cell[k]->A[0*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Alfa1=-linename.Cell[k]->A[1*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Beta2= linename.Cell[k]->A[2*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Alfa2=-linename.Cell[k]->A[3*6+2]/linename.Cell[k]->A[3*6+3];
    
    temp=linename.Cell[k]->A[0*6+0]*linename.Cell[k]->A[1*6+1];
    temp=temp*linename.Cell[k]->A[2*6+2]*linename.Cell[k]->A[3*6+3];
    temp=sqrt(sqrt(temp));
    linename.Cell[k]->r=temp;
    
    for(i=0;i<2;i++){
      for(j=0;j<2;j++) {
	A12[i*2+j]= linename.Cell[k]->A[6*i+(j+2)];
	A22[i*2+j]= linename.Cell[k]->A[6*(i+2)+(j+2)];
      }    
    }
    
    mat_inv(A22,2);
    mat_mult(A12, A22,C,2,2,2);
    linename.Cell[k]->c11=temp*C[0];
    linename.Cell[k]->c12=temp*C[1];
    linename.Cell[k]->c21=temp*C[2];
    linename.Cell[k]->c22=temp*C[3];
  }
}

//------trace the Twiss parameters through the ring to the end
void Trace_Twiss(Line & linename, double deltap, double x[4], double Beta1, double Beta2, double Alfa1, double Alfa2, double c11, double c12, double c21, double c22)
{
  double r;
  double x0[6];
  int i,j,k;
  double temp;
  double A12[4], A22[4], C[4];

  for(i=0;i<4;i++)  x0[i]=x[i];
  x0[4]=0.000;
  x0[5]=deltap;
  for(j=0;j<linename.Ncell; j++) {
    linename.Cell[j]->Pass(x0);
    for(i=0;i<6;i++) linename.Cell[j]->X[i]=x0[i];
  } 
  
  r= sqrt(1 - ( c11*c22-c12*c21));
  linename.Cell[linename.Ncell-1]->A[0*6+0]= r * sqrt(Beta1);
  linename.Cell[linename.Ncell-1]->A[0*6+1]= 0; 
  linename.Cell[linename.Ncell-1]->A[0*6+2]= c11*sqrt(Beta2)- c12*Alfa2/sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[0*6+3]= c12 / sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[1*6+0]= -Alfa1*r / sqrt(Beta1);
  linename.Cell[linename.Ncell-1]->A[1*6+1]= r/sqrt(Beta1); 
  linename.Cell[linename.Ncell-1]->A[1*6+2]= c21*sqrt(Beta2)-c22*Alfa2/sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[1*6+3]= c22 / sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[2*6+0]=  -c12*Alfa1 / sqrt(Beta1) - c22*sqrt(Beta1);
  linename.Cell[linename.Ncell-1]->A[2*6+1]=  c12/sqrt(Beta1); 
  linename.Cell[linename.Ncell-1]->A[2*6+2]= r*sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[2*6+3]= 0.;
  linename.Cell[linename.Ncell-1]->A[3*6+0]= c11*Alfa1 / sqrt(Beta1) + c21*sqrt(Beta1);
  linename.Cell[linename.Ncell-1]->A[3*6+1]= -c11/sqrt(Beta1); 
  linename.Cell[linename.Ncell-1]->A[3*6+2]= -Alfa2*r/sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[3*6+3]= r/sqrt(Beta2);
  
  Cal_ElementMap(linename, deltap);
  Trace_A(linename, deltap);

  for(k=0; k<linename.Ncell;k++){
    linename.Cell[k]->Beta1= linename.Cell[k]->A[0*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Alfa1=-linename.Cell[k]->A[1*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Beta2= linename.Cell[k]->A[2*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Alfa2=-linename.Cell[k]->A[3*6+2]/linename.Cell[k]->A[3*6+3];
    
    temp=linename.Cell[k]->A[0*6+0]*linename.Cell[k]->A[1*6+1];
    temp=temp*linename.Cell[k]->A[2*6+2]*linename.Cell[k]->A[3*6+3];
    temp=sqrt(sqrt(temp));
    linename.Cell[k]->r=temp;
    
    for(i=0;i<2;i++){
      for(j=0;j<2;j++) {
	A12[i*2+j]= linename.Cell[k]->A[6*i+(j+2)];
	A22[i*2+j]= linename.Cell[k]->A[6*(i+2)+(j+2)];
      }    
    }
    
    mat_inv(A22,2);
    mat_mult(A12, A22,C,2,2,2);
    linename.Cell[k]->c11=temp*C[0];
    linename.Cell[k]->c12=temp*C[1];
    linename.Cell[k]->c21=temp*C[2];
    linename.Cell[k]->c22=temp*C[3];
  }
}

//-----trace the orbit through the ring to the end 
void Trace_Orbit(Line & linename, double x[6], int i1)
{
  int i,j;
  for(i=i1;i<linename.Ncell;i++){
    linename.Cell[i]->Pass(x);
    for(j=0;j<6;j++) linename.Cell[i]->X[j]=x[j];
  }
}

//---6D orbit and Twiss
void Cal_Orbit_Num_6D(Line & linename)
{
  int i,j,k,iter=0;
  double x0[6],x01[6], x1[6], dx[6];
  double d=1.0e-09;
  double mat[6][6];
  double chi;
  double codeps=1e-12;
  int Max_iter=20;
  
  for (i=0;i<6;i++) x0[i]=0.;

  do{
    iter++;
    
    for(j=0;j<6;j++) x01[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x01);
    
    for(k=0;k<6;k++){
      for(j=0;j<6;j++) x1[j]=x0[j];
      x1[k]=x1[k]+d;
      for(j=0;j<linename.Ncell;j++) linename.Cell[j]->Pass(x1);
      for(j=0;j<6;j++) mat[j][k]=(x1[j]-x01[j])/d;
    }
    
    for(j=0;j<6;j++) x1[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x1);
    for(i=0;i<6;i++) {
      dx[i]=x0[i]-x1[i];
    }
    chi=0;
    for(i=0;i<6;i++) chi+=dx[i]*dx[i];
    chi=sqrt(chi/6.);
    
    if( chi > codeps ){
      for(i=0;i<6;i++) mat[i][i]=mat[i][i]-1.0000001;
      mat_inv(&mat[0][0],6); 
      for(i=0;i<6;i++){
	for(j=0;j<6;j++) x0[i]+=mat[i][j]*dx[j];  
      } 
    } 
    
  }while (chi>codeps && iter <Max_iter );
  
  if(iter == Max_iter-1 ) 
    {
      cout<<"Failed to find COD."<<endl;
      exit(1);
    }
  else
    {
      for(j=0;j<linename.Ncell; j++) {
	linename.Cell[j]->Pass(x0);
	for(i=0;i<6;i++) linename.Cell[j]->X[i]=x0[i];
      } 
    }
}

void Cal_OneTurnMap_6D(Line & linename)
// there is not such clsoed orbits for different constant deltaps like in 5-d simulation
{
  int i,j;
  double x[6];
  tps tps1[6];
  linmap m1;
  double m66[36], b66[36];
  int flag;
  double u[6], v[6];

  for(i=0;i<6;i++) x[i]=linename.Cell[linename.Ncell-1]->X[i];
  m1.identity();
  m1=m1+x;
  
  for(i=0;i<6;i++) tps1[i]=m1[i];
  for(j=0;j<linename.Ncell; j++) linename.Cell[j]->DAPass(tps1);
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, m66); 

  for(i=0;i<36;i++) linename.Cell[linename.Ncell-1]->M[i]=m66[i];
  
  if( false ) {
    cout<<"One turn map M:"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<m66[i*6+j]<<"  ";
      cout<<endl;
    }
    for(i=0;i<36;i++) b66[i]=m66[i];
    cout<<"Det of M = "<<mat_det(b66,6)<<endl;
  }

  mat_change_hessenberg(m66, 6);
  flag=mat_root_hessenberg(m66,6,u,v,1.0e-10,60);
  if(flag < 0){
    cout<<"Unstable motion. "<<endl;
    exit(1);
  }
}

void Cal_ElementMap_6D(Line & linename)
{
  int i,j,k, istart;
  double x[6];
  tps tps1[6];
  linmap m1;
  double t66[36];

  for(k=0; k<linename.Ncell; k++) 
    {
      istart=k-1;
      if (istart <0) istart=linename.Ncell-1;

      for(i=0;i<6;i++) x[i]=linename.Cell[istart]->X[i];
      m1.identity();
      m1=m1+x;
      
      for(i=0;i<6;i++) tps1[i]=m1[i];
      linename.Cell[k]->DAPass(tps1);
      for(i=0;i<6;i++) m1[i]=tps1[i];
      Getmat(m1, t66);
      
      for(j=0;j<36;j++) linename.Cell[k]->T[j]=t66[j]; 
    }
}

void Cal_A_6D(Line & linename)
{
  
  int i,j;
  double m66[6][6];
  double wr[6], wi[6], vr[6][6], vi[6][6], temp_wr, temp_wi, temp_vr[6], temp_vi[6];  
  double temp1, temp2, temp3;
  double theta1, theta2, theta3, tempx, tempy, x1,y1;
  double a66[6][6], b66[6][6], c44[4][4];
  double det, scale;
  
  //---solve eigen problem to construct A[6*6]
  for(i=0;i<6;i++)
   for(j=0;j<6;j++) m66[i][j]=linename.Cell[linename.Ncell-1]->M[i*6+j];

  EigenSolver_6D(m66, wr, wi, vr, vi);
  if(false){
    cout<<"Eigen values and eigen vectors: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++){
        cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl; 
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

 //---link eigen vetcors S planes 
  temp1=abs(atan2(wi[0], wr[0]))/2/PI;
  temp2=abs(atan2(wi[2], wr[2]))/2/PI; 
  temp3=abs(atan2(wi[4], wr[4]))/2/PI;
  cout<<"Tunes: "<<temp1<<"  "<<temp2<<"   "<<temp3<<endl;
  if( temp1 < temp2 and temp1 < temp3) {
      temp_wr    =  wr[4] ;
      temp_wi    =  wi[4] ;
      wr[4]      =  wr[0] ;
      wi[4]      =  wi[0] ;
      wr[0]      =  temp_wr;
      wi[0]      =  temp_wi;
      wr[1]      =  wr[0];
      wi[1]      = -wi[0]; 
      wr[5]      =  wr[4];
      wi[5]      = -wi[4]; 
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[4][i];
      temp_vi[i]    = vi[4][i];
      vr[4][i]      = vr[0][i];
      vi[4][i]      = vi[0][i];
      vr[0][i]      = temp_vr[i];
      vi[0][i]      = temp_vi[i];
      vr[1][i]      = vr[0][i];
      vi[1][i]      =-vi[0][i]; 
      vr[5][i]      = vr[4][i];
      vi[5][i]      =-vi[4][i]; 
    }
  }
  
  if( temp2 < temp1 and temp2 < temp3) {
      temp_wr    =  wr[4] ;
      temp_wi    =  wi[4] ;
      wr[4]      =  wr[2] ;
      wi[4]      =  wi[2] ;
      wr[2]      =  temp_wr;
      wi[2]      =  temp_wi;
      wr[3]      = wr[2];
      wi[3]      =-wi[2]; 
      wr[5]      = wr[4];
      wi[5]      =-wi[4]; 
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[4][i];
      temp_vi[i]    = vi[4][i];
      vr[4][i]      = vr[2][i];
      vi[4][i]      = vi[2][i];
      vr[2][i]      = temp_vr[i];
      vi[2][i]      = temp_vi[i];
      vr[3][i]      = vr[2][i];
      vi[3][i]      =-vi[2][i]; 
      vr[5][i]      = vr[4][i];
      vi[5][i]      =-vi[4][i]; 
    }
  }

  if(false){
    cout<<"Eigen values and eigen vectors after s mode fixed: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++) cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
    }
  }

 //---link eigen vetcors H/V planes 
  temp1=abs(vr[0][0] ) + abs( vi[0][0] ) + abs( vr[0][1] ) + abs( vi[0][1] );
  temp2=abs(vr[2][0] ) + abs( vi[2][0] ) + abs( vr[2][1] ) + abs( vi[2][1] );
  if (temp2 > temp1 ) {
      temp_wr    =  wr[0] ;
      temp_wi    =  wi[0] ;
      wr[0]      =  wr[2] ;
      wi[0]      =  wi[2] ;
      wr[2]      =  temp_wr;
      wi[2]      =  temp_wi;
      wr[1]      =  wr[0];
      wi[1]      = -wi[0]; 
      wr[3]      =  wr[2];
      wi[3]      = -wi[2]; 
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[0][i];
      temp_vi[i]    = vi[0][i];
      vr[0][i]      = vr[2][i];
      vi[0][i]      = vi[2][i];
      vr[2][i]      = temp_vr[i];
      vi[2][i]      = temp_vi[i];
      vr[1][i]      = vr[0][i];
      vi[1][i]      =-vi[0][i]; 
      vr[3][i]      = vr[2][i];
      vi[3][i]      =-vi[2][i]; 
    }
  }

  if(false){
    cout<<"Eigen values and eigen vectors after H/V/S modes fixed: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++){
	cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl;
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

  //----if here we swap eigen vetcor 4 and 5 
  if(false){
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[4][i];
      temp_vi[i]    = vi[4][i];
      vr[4][i]      = vr[5][i];
      vi[4][i]      = vi[5][i];
      vr[5][i]      = temp_vr[i];
      vi[5][i]      = temp_vi[i];
    }
  }

  //----rotating to make A12=0, A34=0, A56=0
  theta1=atan2(vi[0][0], vr[0][0]);
  x1  =cos(-PI/2-theta1);
  y1  =sin(-PI/2-theta1);
  for(i=0;i<6;i++){
    tempx=vr[0][i]*x1-vi[0][i]*y1;
    tempy=vi[0][i]*x1+vr[0][i]*y1;
    vr[0][i]=tempx;
    vi[0][i]=tempy; }
  x1 =cos(PI/2+theta1);
  y1 =sin(PI/2+theta1);
  for(i=0;i<6;i++){
    tempx=vr[1][i]*x1-vi[1][i]*y1;
    tempy=vi[1][i]*x1+vr[1][i]*y1;
    vr[1][i]=tempx;
    vi[1][i]=tempy;}

  theta2=atan2(vi[2][2], vr[2][2]);
  x1 =cos(-PI/2-theta2);
  y1 =sin(-PI/2-theta2);
  for(i=0;i<6;i++){
    tempx=vr[2][i]*x1-vi[2][i]*y1;
    tempy=vi[2][i]*x1+vr[2][i]*y1;
    vr[2][i]=tempx;
    vi[2][i]=tempy; }
  x1=cos(PI/2+theta2);
  y1=sin(PI/2+theta2);
  for(i=0;i<6;i++){
    tempx=vr[3][i]*x1-vi[3][i]*y1;
    tempy=vi[3][i]*x1+vr[3][i]*y1;
    vr[3][i]=tempx;
    vi[3][i]=tempy;}

  theta3=atan2(vi[4][4], vr[4][4]);
  x1 =cos(-PI/2-theta3);
  y1 =sin(-PI/2-theta3);
  for(i=0;i<6;i++){
    tempx=vr[4][i]*x1-vi[4][i]*y1;
    tempy=vi[4][i]*x1+vr[4][i]*y1;
    vr[4][i]=tempx;
    vi[4][i]=tempy; }
  x1=cos(PI/2+theta3);
  y1=sin(PI/2+theta3);
  for(i=0;i<6;i++){
    tempx=vr[5][i]*x1-vi[5][i]*y1;
    tempy=vi[5][i]*x1+vr[5][i]*y1;
    vr[5][i]=tempx;
    vi[5][i]=tempy;}

  if(false){
    cout<<"Eigen values and eigen vectors after phase fixed: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++){
	cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl;
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

  //----produce A
  for(i=0;i<6;i++) a66[i][0]=-(vi[0][i]-vi[1][i]);
  for(i=0;i<6;i++) a66[i][1]=  vr[0][i]+vr[1][i];
  for(i=0;i<6;i++) a66[i][2]=-(vi[2][i]-vi[3][i]);
  for(i=0;i<6;i++) a66[i][3]=  vr[2][i]+vr[3][i];
  for(i=0;i<6;i++) a66[i][4]=-(vi[4][i]-vi[5][i]);
  for(i=0;i<6;i++) a66[i][5]=  vr[4][i]+vr[5][i];

  if(true){
    if(a66[1][1] < 0. ) {
      for(i=0;i<6;i++) a66[i][1]=- a66[i][1];
      //cout<<"Mode I : negative beta reverted."<<endl; 
   }
    if(a66[3][3] < 0. ) {
      for(i=0;i<6;i++) a66[i][3]=- a66[i][3];
      //cout<<"Mode II : negative beta reverted."<<endl; 
    }
    if(a66[5][5] < 0. ) {
      for(i=0;i<6;i++) a66[i][5]=- a66[i][5];
      //cout<<"Mode III: negative beta reverted."<<endl; 
    }
  }
  
  if(false) {
    cout<<"A matrix :"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<a66[i][j]<<"  ";
      cout<<endl;
    }
    
    for(i=0;i<6;i++) {
      for(j=0;j<6;j++) b66[i][j]=a66[i][j];
    }
    det= mat_det( &b66[0][0], 6);
    cout<<" Det of A66  =  "<< det <<endl;
  }
  
  //---normalizing 
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)  c44[i][j]=a66[i][j];
  det= mat_det( &c44[0][0], 4);
  scale=sqrt(sqrt(abs(det)));
  for(i=0;i<6;i++) {
    for(j=0;j<6;j++) a66[i][j]=a66[i][j]/scale;
  }
  
  for(i=0;i<6;i++) 
    for(j=0;j<6;j++) 
      linename.Cell[linename.Ncell-1]->A[i*6+j]=a66[i][j];
  
 if(false) {
    cout<<"A matrix :"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<a66[i][j]<<"  ";
      cout<<endl;
    }
    
    for(i=0;i<6;i++) {
      for(j=0;j<6;j++) b66[i][j]=a66[i][j];
    }
    det= mat_det( &b66[0][0], 6);
    cout<<" Det of A  =  "<< det <<endl;
  }

}

void Trace_A_6D(Line & linename)
{
  int i,j,k;
  double T[6][6], A[6][6], G[6][6], RI[6][6];
  double dphi1, dphi2, dphi3;
  double mu1, mu2, mu3;
 
  for(i=0;i<6;i++)
    for(j=0;j<6;j++) A[i][j]= linename.Cell[linename.Ncell-1]->A[i*6+j];

  mu1=0.;
  mu2=0.;
  for(k=0;k<linename.Ncell;k++){
    for (i=0;i<6;i++){
      for(j=0;j<6;j++) 	T[i][j]=linename.Cell[k]->T[i*6+j];
    }
  
    mat_mult( &T[0][0], &A[0][0], &G[0][0], 6,6,6 );
    dphi1=atan2(G[0][1], G[0][0]);
    dphi2=atan2(G[2][3], G[2][2]);
    dphi3=atan2(G[4][5], G[4][4]);

    for(i=0;i<6;i++)
      for(j=0;j<6;j++) RI[i][j]=0;

    RI[0][0]  =  cos(dphi1);   
    RI[0][1]  = -sin(dphi1);  
    RI[1][0]  =  sin(dphi1);  
    RI[1][1]  =  cos(dphi1);

    RI[2][2]  =  cos(dphi2);   
    RI[2][3]  = -sin(dphi2);  
    RI[3][2]  =  sin(dphi2);  
    RI[3][3]  =  cos(dphi2);

    RI[4][4]  =  cos(dphi3);   
    RI[4][5]  = -sin(dphi3);  
    RI[5][4]  =  sin(dphi3);  
    RI[5][5]  =  cos(dphi3);

   mat_mult( &G[0][0], &RI[0][0], &A[0][0], 6,6,6 );
   for(i=0;i<6;i++){
     for(j=0;j<6;j++) linename.Cell[k]->A[i*6+j]=A[i][j];
   }
    
    mu1=mu1+dphi1/ 2./PI;
    mu2=mu2+dphi2/ 2./PI;
    mu3=mu3+dphi3/ 2./PI;
    linename.Cell[k]->Mu1= mu1;
    linename.Cell[k]->Mu2= mu2;
    linename.Cell[k]->Mu3= mu3;
  }

  linename.Tune1=  mu1;
  linename.Tune2=  mu2;
  linename.Tune3=  mu3;
}

void Cal_Twiss_6D(Line & linename)
{
  int i,j,k;
  double temp;
  double A12[4], A22[4], C[4];

  Cal_Orbit_Num_6D(linename);
  Cal_OneTurnMap_6D(linename);
  Cal_A_6D(linename);
  Cal_ElementMap_6D(linename);
  Trace_A_6D(linename);
  
  for(k=0; k<linename.Ncell;k++){
    linename.Cell[k]->Beta1= linename.Cell[k]->A[0*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Alfa1=-linename.Cell[k]->A[1*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Beta2= linename.Cell[k]->A[2*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Alfa2=-linename.Cell[k]->A[3*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Beta3= linename.Cell[k]->A[4*6+4]/linename.Cell[k]->A[5*6+5];
    linename.Cell[k]->Alfa3=-linename.Cell[k]->A[5*6+4]/linename.Cell[k]->A[5*6+5];
    
    temp=linename.Cell[k]->A[0*6+0]*linename.Cell[k]->A[1*6+1];
    temp=temp*linename.Cell[k]->A[2*6+2]*linename.Cell[k]->A[3*6+3];
    temp=sqrt(sqrt(temp));
    linename.Cell[k]->r=temp;
    
    for(i=0;i<2;i++){
      for(j=0;j<2;j++) {
	A12[i*2+j]= linename.Cell[k]->A[6*i+(j+2)];
	A22[i*2+j]= linename.Cell[k]->A[6*(i+2)+(j+2)];
      }    
    }
    
    mat_inv(A22,2);
    //C[1]=-A22[1];
    //C[2]=-A22[2];
    //C[3]= A22[0];
    //det = (A22[0]*A22[3]-A22[1]*A22[2]);  
    //for(i=0;i<3;i++) A22[i]= C[0] / det ;
    //if(det==0){
    //  cout<<" Can't invert A22 in Cal_Twiss(). "<<endl;
    //  exit(1);
    //}

    mat_mult(A12, A22,C,2,2,2);
    linename.Cell[k]->c11=temp*C[0];
    linename.Cell[k]->c12=temp*C[1];
    linename.Cell[k]->c21=temp*C[2];
    linename.Cell[k]->c22=temp*C[3];
  }
}

void Cal_Chrom( Line & linename)
{
  double deltap;
  double qx0, qy0, qxp, qyp, qxm, qym;

  deltap=0.0003;
  Cal_Twiss(linename, deltap);
  qxp=linename.Tune1;
  qyp=linename.Tune2;
  
  deltap=-0.0003;
  Cal_Twiss(linename, deltap);
  qxm=linename.Tune1;
  qym=linename.Tune2;

  deltap=0.00;
  Cal_Twiss(linename,deltap);
  qx0=linename.Tune1;
  qy0=linename.Tune2;
  
  linename.Chromx1=(qxp-qxm)/2./0.0003;
  linename.Chromy1=(qyp-qym)/2./0.0003;
  linename.Chromx2=(qxp+qxm-2*qx0)/2./0.0003/0.0003;
  linename.Chromy2=(qyp+qym-2*qy0)/2./0.0003/0.0003;
}

void Cal_Dispersion(Line & linename)
{
  int i,j, k;
  double m66[6][6];
  double r44[4][4];
  double x[6], xtemp[6];  
  double T66[6][6];
   
  for(i=0;i<6;i++) 
    for(j=0;j<6;j++) m66[i][j]= linename.Cell[linename.Ncell-1]->M[i*6+j];
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) r44[i][j]= linename.Cell[linename.Ncell-1]->M[i*6+j];
  for(i=0;i<4;i++) r44[i][i]  = r44[i][i]-1.0;
  mat_inv(&r44[0][0],4); 
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) r44[i][j]= -1.0 * r44[i][j];
  linename.Cell[linename.Ncell-1]->Etax  = r44[0][0] * m66[0][5] +  r44[0][1] * m66[1][5] +   r44[0][2] * m66[2][5] +  r44[0][3] * m66[3][5];
  linename.Cell[linename.Ncell-1]->Etaxp = r44[1][0] * m66[0][5] +  r44[1][1] * m66[1][5] +   r44[1][2] * m66[2][5] +  r44[1][3] * m66[3][5];
  linename.Cell[linename.Ncell-1]->Etay  = r44[2][0] * m66[0][5] +  r44[2][1] * m66[1][5] +   r44[2][2] * m66[2][5] +  r44[2][3] * m66[3][5];
  linename.Cell[linename.Ncell-1]->Etayp = r44[3][0] * m66[0][5] +  r44[3][1] * m66[1][5] +   r44[3][2] * m66[2][5] +  r44[3][3] * m66[3][5];

  x[0]= linename.Cell[linename.Ncell-1]->Etax  ;
  x[1]= linename.Cell[linename.Ncell-1]->Etaxp ;
  x[2]= linename.Cell[linename.Ncell-1]->Etay  ;
  x[3]= linename.Cell[linename.Ncell-1]->Etayp ;
  for(k=0;k<linename.Ncell; k++){
    x[4]= 0.;  x[5]= 1;
    for(i=0;i<6;i++) 
      for(j=0;j<6;j++) T66[i][j]= linename.Cell[k]->T[i*6+j];
    xtemp[0]=  x[0] * T66[0][0] +  x[1] * T66[0][1] +   x[2] * T66[0][2] +  x[3] * T66[0][3]  +   x[5] * T66[0][5]   ;
    xtemp[1]=  x[0] * T66[1][0] +  x[1] * T66[1][1] +   x[2] * T66[1][2] +  x[3] * T66[1][3]  +   x[5] * T66[1][5]   ;
    xtemp[2]=  x[0] * T66[2][0] +  x[1] * T66[2][1] +   x[2] * T66[2][2] +  x[3] * T66[2][3]  +   x[5] * T66[2][5]   ;
    xtemp[3]=  x[0] * T66[3][0] +  x[1] * T66[3][1] +   x[2] * T66[3][2] +  x[3] * T66[3][3]  +   x[5] * T66[3][5]   ;
    linename.Cell[k]->Etax  =  xtemp[0];
    linename.Cell[k]->Etaxp =  xtemp[1];
    linename.Cell[k]->Etay  =  xtemp[2];
    linename.Cell[k]->Etayp =  xtemp[3];
    for(i=0;i<4;i++) x[i]= xtemp[i];
  }
}

//-----------------------------------------
//   fitting tunes and linear chroms
//-----------------------------------------
void Fit_Tune(Line & linename, double q1, double q2, const char * qf_name, const char * qd_name)
{
  double tunex0, tuney0, tunex1=q1, tuney1=q2, dtunex, dtuney;
  double qf_k1l_0, qd_k1l_0;
  double dk1l_qf, dtunex_qf,  dtuney_qf, dk1l_qd, dtunex_qd,  dtuney_qd;
  double scale_qf, scale_qd;

  Cal_Twiss(linename,0.0);
  tunex0=linename.Tune1;
  tuney0=linename.Tune2;

  while( (tunex1-tunex0)*(tunex1-tunex0) + (tuney1-tuney0)*(tuney1-tuney0)  > 1.0e-10 ) {
    qf_k1l_0= Get_KL(linename,qf_name, "K1L");
    qd_k1l_0= Get_KL(linename, qd_name,"K1L"); 
    
    dk1l_qf=  qf_k1l_0 * 0.001;
    Set_dKL(linename,qf_name, "K1L", dk1l_qf);
    Cal_Twiss(linename,0.0);
    dtunex_qf=linename.Tune1 - tunex0;
    dtuney_qf=linename.Tune2 - tuney0;
    Set_dKL(linename,qf_name, "K1L",-dk1l_qf);
    
    dk1l_qd=  qd_k1l_0 * 0.001;
    Set_dKL(linename,qd_name, "K1L", dk1l_qd);
    Cal_Twiss(linename,0.0);
    dtunex_qd=linename.Tune1 - tunex0;
    dtuney_qd=linename.Tune2 - tuney0;
    Set_dKL(linename,qd_name, "K1L",-dk1l_qd);
    
    dtunex=tunex1- tunex0;
    dtuney=tuney1- tuney0;
    
    LinearEquations(dtunex_qf, dtunex_qd,dtunex, dtuney_qf, dtuney_qd, dtuney, scale_qf, scale_qd);
    Set_dKL(linename,qf_name, "K1L", dk1l_qf * scale_qf);
    Set_dKL(linename,qd_name, "K1L", dk1l_qd * scale_qd);
    
    Cal_Twiss(linename,0.0);
    tunex0=linename.Tune1;
    tuney0=linename.Tune2;
  }
}

void Fit_Tune_RHICelens(Line & linename, double q1, double q2)
{
  double tunex0, tuney0, tunex1=q1, tuney1=q2, dtunex, dtuney;
  double qf_k1l_0, qd_k1l_0;
  double dk1l_qf, dtunex_qf,  dtuney_qf, dk1l_qd, dtunex_qd,  dtuney_qd;
  double scale_qf, scale_qd;

  Cal_Twiss(linename,0.0);
  tunex0=linename.Tune1;
  tuney0=linename.Tune2;

  while( (tunex1-tunex0)*(tunex1-tunex0) + (tuney1-tuney0)*(tuney1-tuney0)  > 1.0e-10 ) {
    qf_k1l_0= Get_KL(linename,"QF", "K1L");
    qd_k1l_0= Get_KL(linename, "QD","K1L"); 
    
    dk1l_qf=  qf_k1l_0 * 0.001;
    Set_dKL(linename,"QF", "K1L", dk1l_qf);
    Set_dKL(linename,"QFSHFT", "K1L", dk1l_qf);
    Set_dKL(linename,"QFSHFT2", "K1L", dk1l_qf/2);

    Cal_Twiss(linename,0.0);
    dtunex_qf=linename.Tune1 - tunex0;
    dtuney_qf=linename.Tune2 - tuney0;
    Set_dKL(linename,"QF", "K1L",-dk1l_qf);
    Set_dKL(linename,"QFSHFT", "K1L", -dk1l_qf);
    Set_dKL(linename,"QFSHFT2", "K1L", -dk1l_qf/2);
    
    dk1l_qd=  qd_k1l_0 * 0.001;
    Set_dKL(linename,"QD", "K1L", dk1l_qd);
    Set_dKL(linename,"QDSHFT", "K1L", dk1l_qd);
    Set_dKL(linename,"QDSHFT2", "K1L", dk1l_qd/2);
    Cal_Twiss(linename,0.0);
    dtunex_qd=linename.Tune1 - tunex0;
    dtuney_qd=linename.Tune2 - tuney0;
    Set_dKL(linename,"QD", "K1L",-dk1l_qd);
    Set_dKL(linename,"QDSHFT", "K1L", -dk1l_qd);
    Set_dKL(linename,"QDSHFT2", "K1L", -dk1l_qd/2);
   
    dtunex=tunex1- tunex0;
    dtuney=tuney1- tuney0;
    
    LinearEquations(dtunex_qf, dtunex_qd,dtunex, dtuney_qf, dtuney_qd, dtuney, scale_qf, scale_qd);
    Set_dKL(linename,"QF", "K1L", dk1l_qf * scale_qf);
    Set_dKL(linename,"QFSHFT", "K1L",  dk1l_qf * scale_qf);
    Set_dKL(linename,"QFSHFT2", "K1L",  dk1l_qf * scale_qf /2 );

    Set_dKL(linename,"QD", "K1L", dk1l_qd * scale_qd);
    Set_dKL(linename,"QDSHFT", "K1L", dk1l_qd * scale_qd);
    Set_dKL(linename,"QDSHFT2", "K1L", dk1l_qd * scale_qd/2);
  
    Cal_Twiss(linename,0.0);
    tunex0=linename.Tune1;
    tuney0=linename.Tune2;
  }
}

void Fit_Chrom(Line & linename, double chrom1x_want, double chrom1y_want, const char * sf_name, const char * sd_name )
{
  double chrom1x0, chrom1y0,  dchrom1x, dchrom1y;
  double dk2l_sf, dk2l_sd, dchrom1x_sf,  dchrom1y_sf, dchrom1x_sd, dchrom1y_sd;
  double scale_sf, scale_sd;

  Cal_Chrom(linename) ; 
  chrom1x0= linename.Chromx1;
  chrom1y0= linename.Chromy1;

  while( (chrom1x_want-chrom1x0)*(chrom1x_want-chrom1x0) + (chrom1y_want-chrom1y0)*(chrom1y_want-chrom1y0)  > 0.0001 ) {
    //sf_k2l_0= Get_KL(linename,sf_name,"K2L");
    //sd_k2l_0= Get_KL(linename,sd_name,"K2L"); 
    
    dk2l_sf= 0.3 * 0.005;
    Set_dKL(linename,sf_name, "K2L", dk2l_sf);
    Cal_Chrom(linename);
    dchrom1x_sf=linename.Chromx1 - chrom1x0;
    dchrom1y_sf=linename.Chromy1 - chrom1y0;
    Set_dKL(linename,sf_name, "K2L",-dk2l_sf);
    
    dk2l_sd= -0.5 * 0.005;
    Set_dKL(linename,sd_name, "K2L", dk2l_sd);
    Cal_Chrom(linename);
    dchrom1x_sd=linename.Chromx1 - chrom1x0;
    dchrom1y_sd=linename.Chromy1 - chrom1y0;
    Set_dKL(linename,sd_name, "K2L",-dk2l_sd);

    dchrom1x=chrom1x_want- chrom1x0;
    dchrom1y=chrom1y_want- chrom1y0;
    
    LinearEquations(dchrom1x_sf, dchrom1x_sd, dchrom1x, dchrom1y_sf, dchrom1y_sd, dchrom1y, scale_sf, scale_sd);
    Set_dKL(linename,sf_name, "K2L", dk2l_sf * scale_sf);
    Set_dKL(linename,sd_name, "K2L", dk2l_sd * scale_sd);
    
    Cal_Chrom(linename);
    chrom1x0=linename.Chromx1;
    chrom1y0=linename.Chromy1;
  }
}

void Fit_Chrom_RHIC8fam(Line & linename, double chrom1x_want, double chrom1y_want )
{
  int i;
  double chrom1x0, chrom1y0,  dchrom1x, dchrom1y;
  double dk2l_sf, dk2l_sd, dchrom1x_sf,  dchrom1y_sf, dchrom1x_sd, dchrom1y_sd;
  double scale_sf, scale_sd;

  Cal_Chrom(linename);
  chrom1x0=linename.Chromx1;
  chrom1y0=linename.Chromy1;

  while( (chrom1x_want-chrom1x0)*(chrom1x_want-chrom1x0) + (chrom1y_want-chrom1y0)*(chrom1y_want-chrom1y0)  > 0.0001 ) {
    
    dk2l_sf =  0.3 * 0.005;
    for(i=0;i<linename.Ncell;i++) 
      if( linename.Cell[i]->NAME.find("SF") != string::npos and linename.Cell[i]->TYPE ==string("SEXT") )
	linename.Cell[i]->SetP("K2L",  linename.Cell[i]->GetP("K2L") +  dk2l_sf );
    Cal_Chrom(linename);
    dchrom1x_sf=linename.Chromx1 - chrom1x0;
    dchrom1y_sf=linename.Chromy1 - chrom1y0;  
    for(i=0;i<linename.Ncell;i++) 
      if(  linename.Cell[i]->NAME.find("SF") != string::npos and linename.Cell[i]->TYPE ==string("SEXT") ) 
	linename.Cell[i]->SetP( "K2L",linename.Cell[i]->GetP("K2L") -  dk2l_sf );

    dk2l_sd =  -0.5 * 0.005;
    for(i=0;i<linename.Ncell;i++) 
      if(  linename.Cell[i]->NAME.find("SD") != string::npos and linename.Cell[i]->TYPE ==string("SEXT")  ) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") +  dk2l_sd );
    Cal_Chrom(linename);
    dchrom1x_sd=linename.Chromx1 - chrom1x0;
    dchrom1y_sd=linename.Chromy1 - chrom1y0;
    for(i=0;i<linename.Ncell;i++) 
       if(  linename.Cell[i]->NAME.find("SD") != string::npos and linename.Cell[i]->TYPE ==string("SEXT") ) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") -  dk2l_sd );

    dchrom1x=chrom1x_want- chrom1x0;
    dchrom1y=chrom1y_want- chrom1y0;
    LinearEquations(dchrom1x_sf, dchrom1x_sd, dchrom1x, dchrom1y_sf, dchrom1y_sd, dchrom1y, scale_sf, scale_sd);
    for(i=0;i<linename.Ncell;i++) 
      if( linename.Cell[i]->NAME.find("SF") != string::npos and linename.Cell[i]->TYPE ==string("SEXT") ) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") + dk2l_sf * scale_sf );
    for(i=0;i<linename.Ncell;i++) 
      if( linename.Cell[i]->NAME.find("SD") != string::npos and linename.Cell[i]->TYPE ==string("SEXT") ) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") +dk2l_sd * scale_sd );
    
    Cal_Chrom(linename);
    chrom1x0=linename.Chromx1;
    chrom1y0=linename.Chromy1;
  }
}

//----------------------------------------
//    Chromatic Calculation 
//----------------------------------------
void  chrom_fit(double qx[],double qy[],double & chromx1,double & chromy1,double & chromx2,double & chromy2,double & chromx3,double & chromy3 )
{
  int i;
  double temp;
  double xa[21],ya[21];
  double coeff[8];

  temp=qx[10];
  for (i=0; i<21;i++){
    xa[i]=(i-10)*GP.step_deltap;
    ya[i]=qx[i]-temp;
  }
  pfit(xa, ya, 21, coeff, 7);
  chromx1=coeff[1];
  chromx2=coeff[2];
  chromx3=coeff[3];  
  
  temp=qy[10];
  for (i=0; i<21;i++){
    xa[i]=(i-10)*GP.step_deltap;
    ya[i]=qy[i]-temp;
  }
  pfit(xa, ya, 21, coeff,7);
  chromy1=coeff[1];
  chromy2=coeff[2];
  chromy3=coeff[3]; 
}

void Cal_Chrom_Num( Line & linename)
{
  int i;
  double qx[21],qy[21];
  double deltap;
  double chromx1,chromy1,chromx2,chromy2,chromx3,chromy3;  

  for(i=0;i<21;i++){
    deltap=GP.step_deltap*(i-10);
    Cal_Twiss(linename, deltap);
    qx[i]=linename.Tune1;
    qy[i]=linename.Tune2;
  }

  chrom_fit(qx,qy,chromx1,chromy1,chromx2,chromy2,chromx3,chromy3);
  
  linename.Chromx1 =  chromx1;
  linename.Chromy1 =  chromy1;
  linename.Chromx2 =  chromx2;
  linename.Chromy2 =  chromy2;
  linename.Chromx3 =  chromx3;
  linename.Chromy3 =  chromy3;

  Cal_Twiss(linename, 0.0);
}

void Correct_Chrom_Manual( Line & linename)
{
  int i;
  char fam1[16], fam2[16];
  double step;
  
  for(i=0;i<8;i++) {
    cout<<"Input two Sextupole families and change step:"<<endl;
    cin>>fam1>>fam2>>step;
    Set_dKL(linename,fam1,"K2L",step);
    Set_dKL(linename,fam2,"K2L",-step);
    Fit_Chrom_RHIC8fam(linename, 1.0, 1.0);
    Cal_Chrom(linename);  //  Cal_Chrom_Num(linename);
    cout<<linename.Chromx1<<"  "<<linename.Chromy1<<endl;
    cout<<linename.Chromx2<<"  "<<linename.Chromy2<<endl;
    //cout<<linename.Chromx3<<"  "<<linename.Chromy3<<endl;
  }
}

void Cal_Tune_vs_Deltap(Line & linename, const char *filename)
{
  int i;
  double qx[21],qy[21];
  double deltap;
  fstream fout;

  for(i=0;i<21;i++){
    deltap=GP.step_deltap*(i-10);
    Cal_Twiss(linename, deltap);
    qx[i]=linename.Tune1;
    qy[i]=linename.Tune2;
  }

  fout.open(filename, ios::out);
  for (i = 0; i < 21; i++)
    fout << setw(10) <<GP.step_deltap*(i-10)
  	 << scientific << setw(15) << qx[i]
  	 << scientific << setw(15) << qy[i]<<endl;
   fout.close();
}

void Plot_Tune_vs_Deltap(Line & linename, const char* filename)
{
  char command[256];
  fstream fout;

  fout.open("temp222.p", ios::out);
  fout<<"set term post color enhanced 20 "<<endl;
  fout<<"set output 'tune_vs_delta.ps' "<<endl;
  fout<<"set xlabel  'dp/p_0 [10^{-3}] '  " <<endl;
  fout<<"set ylabel 'Q_x' " <<endl;
  fout<<"set y2label 'Q_y' "<<endl;
  fout<<"set ytics nomirror"<<endl;
  fout<<"set y2tics"<<endl;
  
  sprintf(command, "plot '%s' u ($1*1000):2 tit 'Q_x' w l lt 1  lw 2,\\", filename );
  fout<<command<<endl;
  sprintf(command, "     '%s' u ($1*1000):3 axes x1y2 tit 'Q_y' w l lt 3  lw 2", filename );
  fout<<command<<endl;
  fout<<"exit"<<endl;
  fout.close();
  system("gnuplot temp222.p");
  system("rm temp222.p");
}

void Cal_Beta_Star_vs_Deltap(Line & linename, const char *filename)
{
  int i;
  double betx[21],bety[21];
  double deltap;
  fstream fout;

  for(i=0;i<21;i++){
    deltap=GP.step_deltap*(i-10);
    Cal_Twiss(linename, deltap);
    betx[i]=linename.Cell[linename.Ncell-1]->Beta1;
    bety[i]=linename.Cell[linename.Ncell-1]->Beta2;
  }

  fout.open(filename, ios::out);
  for (i = 0; i < 21; i++)
    fout << setw(10) <<GP.step_deltap*(i-10)
  	 << scientific << setw(15) << betx[i]
  	 << scientific << setw(15) << bety[i]<<endl;
   fout.close();
}

void Plot_Beta_Star_vs_Deltap(Line & linename, const char* filename)
{
  char command[256];
  fstream fout;

  fout.open("temp222.p", ios::out);
  fout<<"set term post color enhanced 20 "<<endl;
  fout<<"set output 'beta_star_vs_delta.ps' "<<endl;
  fout<<"set xlabel  'dp/p_0 [10^{-3}]'  " <<endl;
  fout<<"set ylabel '{/Symbol \142}*_x' " <<endl;
  fout<<"set y2label '{/Symbol \142}*_y' "<<endl;
  fout<<"set ytics nomirror"<<endl;
  fout<<"set y2tics"<<endl;
  
  sprintf(command, "plot '%s' u ($1*1000):2 tit '{/Symbol \142}*_x' w l lt 1  lw 2,\\", filename );
  fout<<command<<endl;
  sprintf(command, "     '%s' u ($1*1000):3 axes x1y2 tit '{/Symbol \142}*_y' w l lt 3  lw 2", filename );
  fout<<command<<endl;
  fout<<"exit"<<endl;
  fout.close();
  system("gnuplot temp222.p");
  system("rm temp222.p");
 }

void Cal_Beta_vs_Deltap(Line & linename, const char *filename)
{
  int i,j;
  double delta[21];
  double betax[21][linename.Ncell], betay[21][linename.Ncell];
  double DBXDD[linename.Ncell], DBYDD[linename.Ncell], DBXDD2[linename.Ncell], DBYDD2[linename.Ncell], DBXDD3[linename.Ncell], DBYDD3[linename.Ncell];
  double input1[21],input2[21];
  double term11, term12,term21, term22, term31, term32; 
  fstream fout;
  
  for(i=0;i<21;i++) delta[i]=0.00002*(i-10);
  
  for(i=0;i<21;i++){
    Cal_Twiss(linename, delta[i]);
    for(j=0;j<linename.Ncell;j++){
      betax[i][j]=linename.Cell[j]->Beta1;
      betay[i][j]=linename.Cell[j]->Beta2;
    }
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  betax[i][j];
      input2[i]=  betay[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DBXDD[j]=term11;   DBXDD2[j]=term21*2;  DBXDD3[j]=term31*6;
    DBYDD[j]=term12;   DBYDD2[j]=term22*2;  DBYDD3[j]=term32*6;
  }
  
  fout.open(filename, ios::out);
  for(i=0;i<linename.Ncell;i++) 
    fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->S
        <<setw(20)<<DBXDD[i]<<setw(20)<<DBXDD2[i]<<setw(20)<<DBXDD3[i]<<setw(20)<<DBYDD[i]<<setw(20)<<DBYDD2[i]<<setw(20)<<DBYDD3[i]<<endl;
  fout.close();
}

void Cal_Dispersion_vs_Deltap(Line & linename, const char *filename)
{
  int i,j;
  double delta[21];
  double xco[21][linename.Ncell], yco[21][linename.Ncell];
  double DX[linename.Ncell], DX2[linename.Ncell], DX3[linename.Ncell];
  double input1[21],input2[21];
  double term11, term12,term21, term22, term31, term32; 
  fstream fout;
 
  for(i=0;i<21;i++) delta[i]=GP.step_deltap*(i-10);
  
  for(i=0;i<21;i++){
    Cal_Twiss(linename, delta[i]);
    for(j=0;j<linename.Ncell;j++){
      xco[i][j]=linename.Cell[j]->X[0];
    }
  }
  
  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  xco[i][j];
      input2[i]=  yco[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DX[j]=term11;  DX2[j]=term21*2;  DX3[j]=term31*6;
  }

  fout.open(filename, ios::out);
  for(i=0;i<linename.Ncell;i++) 
    fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->S<<setw(20)<< DX[i]<<setw(20)<<DX2[i]<<setw(20)<<DX3[i]<<endl;
  fout.close();
}

void Cal_Half_Integer_RDT(Line & linename,  const char* filename)
//  horizontal: h20001, vertical: h00021
{
  int i,j;
  fstream fout;
  double h_real, h_imag, v_real, v_imag;
  double k1l, k2l, betx, bety, Dx,dphix, dphiy; 

  Cal_Twiss(linename,0.0);
  Cal_Dispersion(linename);
  fout.open(filename, ios::out);

  for(i=0;i<linename.Ncell;i++){
    h_real=0; h_imag=0;
    v_real=0; v_imag=0;
    
    for(j=0;j<linename.Ncell;j++){
      if( linename.Cell[j]->TYPE == string("QUAD") ){
	k1l=linename.Cell[j]->GetP("K1L");
	betx=linename.Cell[j]->Beta1;
	bety=linename.Cell[j]->Beta2;
	Dx=linename.Cell[j]->Etax;
	dphix= abs(linename.Cell[j]->Mu1 - linename.Cell[i]->Mu1);
	dphiy= abs(linename.Cell[j]->Mu2 - linename.Cell[i]->Mu2);        
	h_real = h_real  -  k1l*betx*cos(  2*dphix * 2 * PI ) ; 
	h_imag = h_imag  -  k1l*betx*sin(  2*dphix * 2 * PI );	
	v_real = v_real  +  k1l*bety*cos(  2*dphiy * 2 * PI ) ; 
	v_imag = v_imag  +  k1l*bety*sin(  2*dphiy * 2 * PI ) ; 
      }
      else if(linename.Cell[j]->TYPE == string("SEXT") ){
	k2l=linename.Cell[j]->GetP("K2L");
	betx=linename.Cell[j]->Beta1;
	bety=linename.Cell[j]->Beta2;
	Dx=linename.Cell[j]->Etax;
	dphix= abs(linename.Cell[j]->Mu1 - linename.Cell[i]->Mu1);
	dphiy= abs(linename.Cell[j]->Mu2 - linename.Cell[i]->Mu2);  
	h_real = h_real +  k2l*betx*Dx*cos(  2*dphix * 2 * PI ) ;  
	h_imag = h_imag +  k2l*betx*Dx*sin(  2*dphix * 2 * PI );	
	v_real = v_real -  k2l*bety*Dx*cos(  2*dphiy * 2 * PI ) ;  
	v_imag = v_imag -  k2l*bety*Dx*sin(  2*dphiy * 2 * PI ) ;  
      }
      else{}
    }
    fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->S<<setw(20)<<linename.Cell[i]->Beta1<<setw(20)<<linename.Cell[i]->Beta2<<setw(20)<<sqrt(h_real*h_real+h_imag*h_imag)<<"  "<<sqrt(v_real*v_real+v_imag*v_imag)<<endl;
  }

  fout.close();
}

void Cal_Q2_Source(Line & linename, const char* filename )
{
  int i,j;
  double delta[21];
  double xco[21][linename.Ncell], yco[21][linename.Ncell];
  double betax[21][linename.Ncell], betay[21][linename.Ncell];
  double DX[linename.Ncell], DX2[linename.Ncell], DX3[linename.Ncell];
  double DBXDD[linename.Ncell], DBYDD[linename.Ncell], DBXDD2[linename.Ncell], DBYDD2[linename.Ncell], DBXDD3[linename.Ncell], DBYDD3[linename.Ncell];
  double input1[21],input2[21];
  double term11, term12,term21, term22, term31, term32; 
  fstream fout;  

  fout.open(filename, ios::out);

  for(i=0;i<21;i++) delta[i]=0.00002*(i-10);
  for(i=0;i<21;i++){
    Cal_Twiss(linename, delta[i]);
    for(j=0;j<linename.Ncell;j++){
      xco[i][j]=linename.Cell[j]->X[0];
      betax[i][j]=linename.Cell[j]->Beta1;
      betay[i][j]=linename.Cell[j]->Beta2;
    }
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  xco[i][j];
      input2[i]=  yco[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DX[j]=term11;  DX2[j]=term21*2;  DX3[j]=term31*6;
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  betax[i][j];
      input2[i]=  betay[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DBXDD[j]=term11;   DBXDD2[j]=term21*2;  DBXDD3[j]=term31*6;
    DBYDD[j]=term12;   DBYDD2[j]=term22*2;  DBYDD3[j]=term32*6;
  }

  double chromx1,   chromy1, chromx2,  chromy2,  chromx3, chromy3;
  Cal_Twiss(linename,0.);

  for(i=1;i<linename.Ncell;i++){
    
    chromx1=0;  chromy1=0 ; chromx2=0 ;  chromy2=0;  chromx3=0; chromy3=0;
    
    if(linename.Cell[i]->TYPE ==string("SEXT") ) {
      
      chromx1 +=  1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/PI;
      chromy1 +=- 1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/PI;
      
      chromx2 += -1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/PI/2;
      chromy2 += +1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/PI/2;
      
      chromx2 +=  linename.Cell[i]->GetP("K2L") * DX[i]  * DBXDD[i]/4/PI/2;
      chromy2 += -linename.Cell[i]->GetP("K2L") * DX[i]  * DBYDD[i]/4/PI/2;
      
      chromx2 +=  linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/4/PI/2;
      chromy2 += -linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/4/PI/2;
      
      chromx3 += +1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/8/PI/3;
      chromy3 += -1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/8/PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K2L") * DX[i]  * DBXDD[i]/8/PI/3;
      chromy3 += +linename.Cell[i]->GetP("K2L") * DX[i]  * DBYDD[i]/8/PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/8/PI/3;
      chromy3 += +linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/8/PI/3; 
      
      chromx3 +=  linename.Cell[i]->GetP("K2L") *  DX[i] * DBXDD2[i]/8/PI/3;
      chromy3 += -linename.Cell[i]->GetP("K2L") *  DX[i] * DBYDD2[i]/8/PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K2L") *  DX[i] * DBXDD[i]/8/PI/3;
      chromy3 += +linename.Cell[i]->GetP("K2L") *  DX[i] * DBYDD[i]/8/PI/3;       
      
      
      chromx3 += +linename.Cell[i]->GetP("K2L") *  DX2[i] * DBXDD[i]/8/PI/3;
      chromy3 += -linename.Cell[i]->GetP("K2L") *  DX2[i] * DBYDD[i]/8/PI/3; 
      
      chromx3 += +linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX3[i]/8/PI/3;
      chromy3 += -linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX3[i]/8/PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/8/PI/3;
      chromy3 += +linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/8/PI/3;
      
      chromx3 += +linename.Cell[i]->GetP("K2L") *  DBXDD[i] * DX2[i]/8/PI/3;
      chromy3 += -linename.Cell[i]->GetP("K2L") *  DBYDD[i] * DX2[i]/8/PI/3;        
      
      fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->TYPE<<setw(20)<<linename.Cell[i]->S
          <<setw(20)<<chromx1<<setw(20)<<chromy1<<setw(20)<<chromx2<<setw(20)<<chromy2<<setw(20)<<chromx3<<setw(20)<<chromy3<<endl;
    }
    
    if(linename.Cell[i]->TYPE ==string("QUAD") ) {
      chromx1 += - 1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI;
      chromy1 += + 1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI;
      
      chromx2 += +1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI/2;
      chromy2 += -1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI/2;
      
      chromx2 += -linename.Cell[i]->GetP("K1L") * DBXDD[i]/4/PI/2;
      chromy2 += +linename.Cell[i]->GetP("K1L") * DBYDD[i]/4/PI/2;
      
      chromx3 += -1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/8/PI/3;
      chromy3 += +1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/8/PI/3;
      
      chromx3 += +linename.Cell[i]->GetP("K1L") * DBXDD[i]/8/PI/3;
      chromy3 += -linename.Cell[i]->GetP("K1L") * DBYDD[i]/8/PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K1L") * DBXDD2[i]/8/PI/3;
      chromy3 += +linename.Cell[i]->GetP("K1L") * DBYDD2[i]/8/PI/3;
      
      chromx3 += +linename.Cell[i]->GetP("K1L") * DBXDD[i]/8/PI/3;
      chromy3 += -linename.Cell[i]->GetP("K1L") * DBYDD[i]/8/PI/3;
      
      fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->TYPE<<setw(20)<<linename.Cell[i]->S
          <<setw(20)<<chromx1<<setw(20)<<chromy1<<setw(20)<<chromx2<<setw(20)<<chromy2<<setw(20)<<chromx3<<setw(20)<<chromy3<<endl;
    }
  }
  fout.close();
} 

//--------------------------------------------------
//  RDTs and detunings
//--------------------------------------------------
void Cal_Coupling_Coefficient( Line & linename )
// coupling coefficient from whole ring
{
  int i;
  double creal=0.0, cimag=0.0;
  double angle, sin_angle, cos_angle;
  
  for(i=0; i< linename.Ncell;i++)
    if( linename.Cell[i]->TYPE == string("SKEWQ") ){
      creal= creal + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) *
        linename.Cell[i]->GetP("K1SL") / 2. / PI  ;
      cimag= cimag + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) *
        linename.Cell[i]->GetP("K1SL") / 2. / PI  ;
    }

   for(i=0; i< linename.Ncell;i++)
    if( linename.Cell[i]->TYPE == string("QUAD")  and linename.Cell[i]->DT !=0 ){
      creal= creal + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) *
        linename.Cell[i]->GetP("K1L") * sin(-1.0 * linename.Cell[i]->DT*2) / 2. / PI  ;
      cimag= cimag + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) *
        linename.Cell[i]->GetP("K1L") * sin(-1.0* linename.Cell[i]->DT*2) / 2. / PI  ;
    }
 
  angle=  atan( abs(cimag/creal) );
  sin_angle=cimag/ sqrt(cimag*cimag+ creal*creal) ;
  cos_angle=creal/ sqrt(cimag*cimag+ creal*creal);
  
  if( sin_angle  >= 0. and  cos_angle >=0 ) angle=angle;
  if( sin_angle >= 0. and cos_angle <=0 )  angle=3.14159265 - angle;
  if( sin_angle < 0.  and cos_angle <=0)   angle=3.14159265 + angle;
  if( sin_angle < 0.  and cos_angle >=0)   angle=2*3.14159265-angle;

  cout<<" C_real  = "<<creal<<endl;  
  cout<<" C_imag  = "<<cimag<<endl;
  cout<<" C_amp   =  "<<sqrt(cimag*cimag+ creal*creal)<<endl;
  cout<<" C_phase =  "<< angle <<endl;
}

void Cal_Coupling_Coefficient_Single( Line & linename, const char * filename)
// single element's coupling contribution
{
  int i;
  double creal=0.0, cimag=0.0;
  double angle, sin_angle, cos_angle;
  fstream  f1;

  f1.open(filename, ios::out);
      f1<<setw(15)<<" NAME "<<setw(15)<<" TYPE "
	<<setw(15)<<" creal "<<setw(15)<<" cimag " <<setw(15)<<" C_amp "<<setw(15)<<" C_phase "<<endl;

  for(i=0; i< linename.Ncell;i++) {

    if( linename.Cell[i]->TYPE == string("SKEWQ") ){
      creal= sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) / 2. / PI  ;
      cimag= sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) / 2. / PI  ;
      
      angle=  atan( abs(cimag/creal) );
      sin_angle=cimag/ sqrt(cimag*cimag+ creal*creal) ;
      cos_angle=creal/ sqrt(cimag*cimag+ creal*creal);
      if( sin_angle  >= 0. and  cos_angle >=0 ) angle=angle;
      if( sin_angle >= 0. and cos_angle <=0 )  angle=3.14159265 - angle;
      if( sin_angle < 0.  and cos_angle <=0)   angle=3.14159265 + angle;
      if( sin_angle < 0.  and cos_angle >=0)   angle=2*3.14159265-angle;
      f1<<setw(15)<<linename.Cell[i]->NAME<<setw(15)<<setw(15)<<linename.Cell[i]->TYPE
	  <<setw(15)<<creal<<setw(15)<<cimag<<setw(15)<<sqrt(cimag*cimag+ creal*creal)<<setw(15)<<angle<<endl;
    }
    
    if( linename.Cell[i]->TYPE == string("QUAD")  and linename.Cell[i]->DT !=0 ){
      creal= sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) / 2. / PI  ;
      cimag= sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) / 2. / PI  ;
      
      angle=  atan( abs(cimag/creal) );
      sin_angle=cimag/ sqrt(cimag*cimag+ creal*creal) ;
      cos_angle=creal/ sqrt(cimag*cimag+ creal*creal);
      if( sin_angle  >= 0. and  cos_angle >=0 ) angle=angle;
      if( sin_angle >= 0. and cos_angle <=0 )  angle=3.14159265 - angle;
      if( sin_angle < 0.  and cos_angle <=0)   angle=3.14159265 + angle;
      if( sin_angle < 0.  and cos_angle >=0)   angle=2*3.14159265-angle;
      f1<<setw(15)<<linename.Cell[i]->NAME<<setw(15)<<setw(15)<<linename.Cell[i]->TYPE
	<<setw(15)<<creal<<setw(15)<<cimag<<setw(15)<<sqrt(cimag*cimag+ creal*creal)<<setw(15)<<angle<<endl;
    }

  }
   f1.close();
}

void Cal_Sext_RDTs( Line & linename )
// sextupole linear RDT at the staring point
{
  int i;
  double h2100_c, h2100_s;
  double h3000_c, h3000_s;
  double h1011_c, h1011_s;  
  double h1002_c, h1002_s;
  double h1020_c, h1020_s;
  double amp, phi, k2l, betx, bety, mux, muy;

  h2100_c=0 ; h2100_s =0 ;
  h3000_c=0 ; h3000_s =0 ;
  h1011_c=0 ; h1011_s =0 ;  
  h1002_c=0 ; h1002_s =0 ;
  h1020_c=0 ; h1020_s =0 ;

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("SEXT") ||  linename.Cell[i]->TYPE==string("MULT") ){
      k2l=linename.Cell[i]->GetP("K2L");
      betx=linename.Cell[i]->Beta1;
      bety=linename.Cell[i]->Beta2;
      mux =linename.Cell[i]->Mu1 * 2.0 * PI;
      muy =linename.Cell[i]->Mu2 * 2.0 * PI;
      
      amp = k2l * pow(betx,1.5);      phi = mux;
      h2100_c += amp * cos(phi);      h2100_s += amp * sin(phi);

      amp = k2l * pow(betx,1.5);      phi = 3* mux;
      h3000_c += amp * cos(phi);      h3000_s += amp * sin(phi);

      amp = k2l * pow(betx,0.5) * bety;      phi = mux;
      h1011_c += amp * cos(phi);      h1011_s += amp * sin(phi);

      amp = k2l * pow(betx,0.5) * bety;    phi = mux - 2* muy;
      h1002_c += amp * cos(phi);      h1002_s += amp * sin(phi);

      amp = k2l * pow(betx,0.5) * bety;    phi = mux + 2* muy;
      h1020_c += amp * cos(phi);      h1020_s += amp * sin(phi);
     }
  }

  h2100_c=   h2100_c * (-1.0/8);  h2100_s=   h2100_s * (-1.0/8);
  h3000_c=   h3000_c * (-1.0/24); h3000_s=   h3000_s * (-1.0/24);
  h1011_c=   h1011_c * ( 1.0/ 4); h1011_s=   h1011_s * ( 1.0/ 4);
  h1002_c=   h1002_c * ( 1.0/ 8); h1002_s=   h1002_s * ( 1.0/ 8);
  h1020_c=   h1020_c * ( 1.0/ 8); h1020_s=   h1020_s * ( 1.0/ 8);

  cout<<">>>1st order sextupole geometric resonance driving terms: "<<endl;
  cout<<" h2100 :  "<<h2100_c <<" +i  "<<h2100_s<<" , "<<sqrt(h2100_c*h2100_c + h2100_s*h2100_s)<<" / "<<atan2(h2100_s,h2100_c)<<endl;
  cout<<" h3000 :  "<<h3000_c <<" +i  "<<h3000_s<<" , "<<sqrt(h3000_c*h3000_c + h3000_s*h3000_s)<<" / "<<atan2(h3000_s,h3000_c)<<endl;
  cout<<" h1011 :  "<<h1011_c <<" +i  "<<h1011_s<<" , "<<sqrt(h1011_c*h1011_c + h1011_s*h1011_s)<<" / "<<atan2(h1011_s,h1011_c)<<endl;
  cout<<" h1002 :  "<<h1002_c <<" +i  "<<h1002_s<<" , "<<sqrt(h1002_c*h1002_c + h1002_s*h1002_s)<<" / "<<atan2(h1002_s,h1002_c)<<endl;
  cout<<" h1020 :  "<<h1020_c <<" +i  "<<h1020_s<<" , "<<sqrt(h1020_c*h1020_c + h1020_s*h1020_s)<<" / "<<atan2(h1020_s,h1020_c)<<endl;
}

void Cal_Detuning_Sext( Line & linename )
{
  int i,j;
  double axx=0, axy=0, ayy=0;
  double k2l_1, betx_1, bety_1, mux_1, muy_1;
  double k2l_2, betx_2, bety_2, mux_2, muy_2;
  double dmux, dmuy;
  double Qx0=linename.Tune1, Qy0=linename.Tune2;

  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->TYPE==string("SEXT") ||  linename.Cell[i]->TYPE==string("MULT") ) {
      k2l_1=linename.Cell[i]->GetP("K2L");
      betx_1=linename.Cell[i]->Beta1;
      bety_1=linename.Cell[i]->Beta2;
      mux_1 =linename.Cell[i]->Mu1 * 2.0 * PI;
      muy_1 =linename.Cell[i]->Mu2 * 2.0 * PI;
      
      for(j=0;j<linename.Ncell;j++) {
	if(linename.Cell[j]->TYPE==string("SEXT") ||  linename.Cell[j]->TYPE==string("MULT") ) {
	  k2l_2=linename.Cell[j]->GetP("K2L");
	  betx_2=linename.Cell[j]->Beta1;
	  bety_2=linename.Cell[j]->Beta2;
	  mux_2 =linename.Cell[j]->Mu1 * 2.0 * PI;
	  muy_2 =linename.Cell[j]->Mu2 * 2.0 * PI;
	  
	  dmux=abs(   mux_1 -   mux_2 );   dmuy=abs(   muy_1 -   muy_2 );
	  axx += k2l_1 *k2l_2 * pow(betx_1, 1.5) * pow(betx_2, 1.5) 
	    *(      cos(3 *(PI*Qx0-dmux))/sin(3*PI*Qx0) 
		    + 3* cos(PI*Qx0 -dmux )/sin(PI*Qx0) );
	  axy += k2l_1 *k2l_2 *sqrt(betx_1*betx_2) * bety_1  
	    *( -bety_2*cos(dmux+2*dmuy-PI*(Qx0+2*Qy0))/sin(PI*(Qx0+2*Qy0))
	       +bety_2*cos(dmux-2*dmuy-PI*(Qx0-2*Qy0))/sin(PI*(Qx0-2*Qy0))  
	       +2*betx_2*cos(dmux-PI*Qx0)/sin(PI*Qx0) );
	  ayy += k2l_1 *k2l_2 *sqrt(betx_1 *betx_2) * bety_1 * bety_2 
	    *(  cos(dmux+2*dmuy-PI*(Qx0+2*Qy0)) / sin(PI *(Qx0+2*Qy0)) 
		+ cos(dmux-2*dmuy-PI*(Qx0-2*Qy0))/sin(PI*(Qx0-2*Qy0)) 
		+4 * cos(dmux-PI*Qx0) /sin(PI*Qx0)  );
	}
      }
    }
  }
  axx=axx*(-1.0/64/PI); ayy=ayy*(-1.0/64/PI); axy= axy*(1.0/32/PI);
  cout<<">>>Amplitude dependent tune shifts from sextupoles:" <<endl;
  cout<<"   2PI * dQx =  axx * (2Jx) - axy *( 2Jy) "<<endl;
  cout<<"   2PI * dQy = -axy * (2Jx) + ayy *( 2Jy) "<<endl;
  cout<<"  "<<endl;
  cout<<"  a_xx = "<<axx<<endl;
  cout<<"  a_xy = "<<axy<<endl;
  cout<<"  a_yy = "<<ayy<<endl;
}

void Cal_Detuning_Oct( Line & linename )
{
  int i;
  double axx=0, axy=0, ayy=0;
  double k3l, betx, bety;

  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->TYPE==string("OCT") ||  linename.Cell[i]->TYPE==string("MULT") ) {
      k3l=linename.Cell[i]->GetP("K3L");
      betx=linename.Cell[i]->Beta1;
      bety=linename.Cell[i]->Beta2;
      axx +=  k3l*betx*betx;
      axy -=  2*k3l*betx*bety;
      ayy +=  k3l*bety*bety;
    }
  }
  axx=axx*(1.0/32/PI); ayy=ayy*(1.0/32/PI); axy= axy*(1.0/32/PI);
  cout<<">>>Amplitude dependent tune shifts from octupoles:" <<endl;
  cout<<"   dQx =  axx * (2Jx) - axy *( 2Jy) "<<endl;
  cout<<"   dQy = -axy * (2Jx) + ayy *( 2Jy) "<<endl;
  cout<<"  "<<endl;
  cout<<"  a_xx = "<<axx<<endl;
  cout<<"  a_xy = "<<axy<<endl;
  cout<<"  a_yy = "<<ayy<<endl;
}

//-----------------------------------------
//   path length and gamma-t 
//-----------------------------------------
void Cal_Gammat(Line & linename)
{

  int i,j;
  double x[6];
  double deltap;
  double input1[21],input2[21];
  double coeff[8];
  double alfa0, alfa1, alfa2;
  fstream fout;

  for(i=0;i<21;i++){
    deltap=GP.step_deltap*(i-10);
    Cal_Orbit_Num(linename, deltap);
    x[0]=linename.Cell[linename.Ncell-1]->X[0];
    x[1]=linename.Cell[linename.Ncell-1]->X[1];
    x[2]=linename.Cell[linename.Ncell-1]->X[2];
    x[3]=linename.Cell[linename.Ncell-1]->X[3];
    x[4]=0.000;
    x[5]=deltap;
    for(j=0;j<linename.Ncell; j++) linename.Cell[j]->Pass(x);
    input1[i] = deltap;
    input2[i] = -x[4] / linename.Length  ;
  } 
  
  fout.open("DLL0_vers_deltap.dat", ios::out);
  for (i = 0; i < 21; i++)
    fout << scientific << setw(15) << input1[i]<< scientific << setw(15) << input2[i]<<endl;
  fout.close();

  pfit(input1, input2, 21, coeff, 7);
  alfa0=coeff[1];
  alfa1=coeff[2];
  alfa2=coeff[3]; 
  
  linename.Alfa0=alfa0  ;
  linename.Alfa1=alfa1  ;
  linename.Alfa2=alfa2  ;
  linename.Gammat =  1.0/ sqrt(alfa0);
  linename.Slip=(linename.Alfa0 - 1.0/GP.gamma/GP.gamma);
}

void Cal_Orbit_Length(Line & linename, double deltap)
{
  int j;
  double x[6];
  
  Cal_Orbit_Num(linename, deltap);
  x[0]=linename.Cell[linename.Ncell-1]->X[0];
  x[1]=linename.Cell[linename.Ncell-1]->X[1];
  x[2]=linename.Cell[linename.Ncell-1]->X[2];
  x[3]=linename.Cell[linename.Ncell-1]->X[3];
  x[4]=0.000;
  x[5]=deltap;
  for(j=0;j<linename.Ncell; j++) linename.Cell[j]->Pass(x);
  linename.Orbit_Length =  -x[4] * 3.0e8 +  linename.Length;
}

//-------------------------------
//  longitudinal calculations
//--------------------------------

void Cal_Qs(Line & linename)
{
  int i;
  double Vrf_tot=0,  phi_s =0;
  
  for(i=0;i<linename.Ncell;i++)
    if(linename.Cell[i]->TYPE==string("RFCAV") ){
      Vrf_tot +=linename.Cell[i]->GetP("VRF");
    }
  linename.Vrf_tot= Vrf_tot;
  
  Cal_Gammat(linename);
  linename.Qs= sqrt( GP.harm * GP.Q * Vrf_tot * abs(linename.Slip * cos(phi_s) )/ 2/PI/GP.beta/GP.beta/GP.energy/GP.A ) ;
}

void Cal_Bucket_Area(Line & linename)
{
  int i;
  double Vrf_tot=0;
  
  for(i=0;i<linename.Ncell;i++)
    if(linename.Cell[i]->TYPE==string("RFCAV") ){
      Vrf_tot +=linename.Cell[i]->GetP("VRF");
    }

  linename.Bucket_length= linename.Length * 1.0e9 / GP.harm / ( GP.beta * light_speed) ;
  linename.Bucket_height= 2*sqrt(GP.Q*Vrf_tot/2/PI/GP.beta/GP.beta/GP.energy/GP.A/GP.harm/abs(linename.Slip)); 
  linename.Bucket_area=16*1.0e6*sqrt( GP.beta*GP.beta* GP.A *GP.energy * GP.Q * Vrf_tot /2/PI/(2*PI*linename.frev0)/(2*PI*linename.frev0)/GP.harm/GP.harm/GP.harm/linename.Slip ) / GP.A ;
}

double RF_F_function(double phi_s, double phi_right, double phi_left)
{
  int i;
  double dphi,phi;
  double sum=0;

  dphi=(phi_right - phi_left) / 1000;
  for(i=0;i<1000;i++){
    phi= phi_left + dphi * (i+1);
    sum += sqrt( abs(cos(phi_left)-cos(phi) + (phi_left-phi) * sin(phi_s)  ) ) * dphi;
  }
  return sum*sqrt(2.0)/8;
}

void Cal_Bunch_Area(Line & linename, double full_length)
// full length is +/-3sigma_l, that is, 6 sigma_l, in units of ns
{
  double phi_s=0, delta_phi, phi_right, phi_left;

  linename.Bunch_length =  full_length ; 
  delta_phi =  full_length * 2 *PI / linename.Bucket_length ;
  phi_right =   delta_phi /2  ;
  phi_left =   -delta_phi /2  ;
  linename.Bunch_area  =   RF_F_function( 0, phi_right, phi_left) * linename.Bucket_area ;
  linename.Bunch_height =  sqrt(abs(cos(phi_right)  - cos(phi_s) + ( phi_right - phi_s ) * sin(phi_s) ) )
    *  linename.Bucket_area * GP.harm * ( 2*PI*linename.frev0)/ 8/ sqrt(2.0) / (GP.energy*1.0e6)*GP.beta*GP.beta;
}

void Cal_Bunch_Height(Line & linename, double bunch_area)
{
  int i;
  double guess,  bunch_area_0,  bunch_area_1, scale;
  
  Cal_Qs(linename);
  Cal_Bucket_Area(linename);

  i=0;
  guess = linename.Bucket_length*0.5; 
  do{
    Cal_Bunch_Area( linename, guess);
    bunch_area_0  = linename.Bunch_area ;   
    Cal_Bunch_Area( linename, guess+0.05);
    bunch_area_1  = linename.Bunch_area ; 
    scale = 0.05 / ( bunch_area_1 -  bunch_area_0 );
    
    guess = ( bunch_area -  bunch_area_0 ) * scale + guess; 
    Cal_Bunch_Area( linename, guess);
    i++;
  } while (i < 30 && abs(bunch_area - linename.Bunch_area ) > 0.1 );  

}

void Print_Longitudinal_Summary( Line & linename)
{
  cout<<"========================================="<<endl;
  cout<<"Beam energy  = "<<GP.energy<<"  MeV "<<endl;
  cout<<"gamma        = "<<GP.gamma<<endl;
  cout<<"beta         = "<<GP.beta<<endl;
  cout<<"circumference= "<<linename.Length<<"  m   "<<endl;
  cout<<"Revolution frequency = "<< linename.frev0 <<"  Hz  "<<endl;
  cout<<"particle's A : "<<GP.A<<endl;
  cout<<"particle's Q  : "<<GP.Q<<endl;
  cout<<"harmnic number = "<<GP.harm<<endl;
  cout<<"RF freq = "<< linename.frf <<"  Hz "<<endl;
  cout<<"RF total voltage = "<< linename.Vrf_tot <<"  MV "<<endl; 
  cout<<"GammaT  = "<< linename.Gammat<<endl;
  cout<<"Phase slip factor = " << linename.Slip<<endl;
  cout<<"Qs= "<<linename.Qs<<endl;  
  cout<<"bucket length = "<<linename.Bucket_length<<"  ns   "<<endl; 
  cout<<"bucket height  (dp/p0_max) = "<<linename.Bucket_height<<endl; 
  cout<<"bucket area (un-normalized) = "<<linename.Bucket_area << "  eV.s/n  "<<endl;  
  cout<<"bunch length = "<<linename.Bunch_length<<"  ns  "<<endl; 
  cout<<"bunch area (un-normalized) = "<<linename.Bunch_area<<"  eV.s/n   "<<endl; 
  cout<<"bunch height (dp/p0_max) = "<<linename.Bunch_height<<endl; 
}

//------------------------------------------------
//        Print and Plot
//------------------------------------------------
void Print_Twiss(Line & linename, const char* filename)
{
 int i;
 fstream fout;
 
 fout.open(filename, ios::out);
 fout <<setw(15) <<"NAME"<<setw(15) <<"TYPE"<<setw(15) <<"S"
      <<setw(15) <<"X   "<<setw(15)  <<"PX    "<<setw(15)<<"Y    "<<setw(15)<<"PY"
      <<setw(15) <<"BETA1"<<setw(15) <<"ALFA1"<<setw(15) <<"MU1"
      <<setw(15) <<"BETA2"<<setw(15) <<"ALFA2"<<setw(15) <<"MU2"
      <<setw(15) <<"ETAX"<<setw(15)  <<"ETAXP"<<setw(15) <<"ETAY"<<setw(15)  <<"ETAYP"<<endl;
 for(i=0;i<linename.Ncell;i++)
   fout <<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->S
	<<setw(15) <<linename.Cell[i]->X[0]<<setw(15) <<linename.Cell[i]->X[1]<<setw(15) <<linename.Cell[i]->X[2]<<setw(15) <<linename.Cell[i]->X[3]
	<<setw(15) <<linename.Cell[i]->Beta1<<setw(15) <<linename.Cell[i]->Alfa1<<setw(15) <<linename.Cell[i]->Mu1
	<<setw(15) <<linename.Cell[i]->Beta2<<setw(15) <<linename.Cell[i]->Alfa2<<setw(15) <<linename.Cell[i]->Mu2
	<<setw(15) <<linename.Cell[i]->Etax<<setw(15) <<linename.Cell[i]->Etaxp<<setw(15) <<linename.Cell[i]->Etay<<setw(15) <<linename.Cell[i]->Etayp<<endl;
 fout.close();
}

void Print_Twiss_6D(Line & linename, const char* filename)
{
 int i;
 fstream fout;
 
 fout.open(filename, ios::out);
 fout <<setw(15) <<"NAME"<<setw(15) <<"TYPE"<<setw(15) <<"S"
      <<setw(15) <<"X   "<<setw(15)  <<"Y    "<<setw(15)<<"Z    "<<setw(15)<<"Delta"
      <<setw(15) <<"BETA1"<<setw(15) <<"ALFA1"<<setw(15) <<"MU1"
      <<setw(15) <<"BETA2"<<setw(15) <<"ALFA2"<<setw(15) <<"MU2"
      <<setw(15) <<"BETA3"<<setw(15) <<"ALFA3"<<setw(15) <<"MU3"
      <<setw(15) <<"ETAX"<<setw(15)  <<"ETAXP"<<setw(15) <<"ETAY"<<setw(15)  <<"ETAYP"<<endl;
 for(i=0;i<linename.Ncell;i++)
   fout <<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->S
	<<setw(15) <<linename.Cell[i]->X[0]<<setw(15) <<linename.Cell[i]->X[2]<<setw(15) <<linename.Cell[i]->X[4]<<setw(15) <<linename.Cell[i]->X[5]
	<<setw(15) <<linename.Cell[i]->Beta1<<setw(15) <<linename.Cell[i]->Alfa1<<setw(15) <<linename.Cell[i]->Mu1
	<<setw(15) <<linename.Cell[i]->Beta2<<setw(15) <<linename.Cell[i]->Alfa2<<setw(15) <<linename.Cell[i]->Mu2
	<<setw(15) <<linename.Cell[i]->Beta3<<setw(15) <<linename.Cell[i]->Alfa3<<setw(15) <<linename.Cell[i]->Mu3
	<<setw(15) <<linename.Cell[i]->Etax<<setw(15) <<linename.Cell[i]->Etaxp<<setw(15) <<linename.Cell[i]->Etay<<setw(15) <<linename.Cell[i]->Etayp<<endl;
 fout.close();
}

void Print_A_Matrix(Line & linename, const char* filename)
{
  int i,j;
 fstream fout;
 
 fout.open(filename, ios::out);
 for(i=0;i<linename.Ncell;i++) {
   fout<<">>>>"<<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->L<<setw(15) <<linename.Cell[i]->S<<endl;
   for(j=0;j<6;j++)
     fout<<setw(15) <<linename.Cell[i]->A[j*6+0]<<setw(15) <<linename.Cell[i]->A[j*6+1]<<setw(15) <<linename.Cell[i]->A[j*6+2]<<setw(15)
	 <<setw(15) <<linename.Cell[i]->A[j*6+3]<<setw(15) <<linename.Cell[i]->A[j*6+4]<<setw(15) <<linename.Cell[i]->A[j*6+5]<<endl;
 }
 fout.close();
}

void Print_Optics_Summary(Line & linename)
{
  cout<<"-----------------------------------------------------------"<<endl;
  cout<<setw(12)<<" Optics Summary:"<<endl;
  cout<<setw(12)<<" Length:" <<setw(20)<<setprecision(15)<<linename.Length<<endl; 
  cout<<setw(12)<<" Tunes :" <<setw(20)<<setprecision(15)<<linename.Tune1<<setw(20)<<setprecision(15)<<linename.Tune2<<endl;
  cout<<setw(12)<<" Chrom1:" <<setw(20)<<setprecision(15)<<linename.Chromx1<<setw(20)<<setprecision(15)<<linename.Chromy1<<endl;
  cout<<setw(12)<<" Chrom2:" <<setw(20)<<setprecision(15)<<linename.Chromx2<<setw(20)<<setprecision(15)<<linename.Chromy2<<endl;
  cout<<setw(12)<<" Chrom3:" <<setw(20)<<setprecision(15)<<linename.Chromx3<<setw(20)<<setprecision(15)<<linename.Chromy3<<endl;
  cout<<setw(12)<<" Beta* :" <<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Beta1<<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Beta2<<endl;
  cout<<setw(12)<<" Alfa* :" <<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Alfa1<<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Alfa2<<endl;
  cout<<setw(12)<<" Etax* :" <<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Etax<<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Etaxp<<endl;
  cout<<setw(12)<<" Etay* :" <<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Etay<<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Etayp<<endl;
  int  i;
  double betax_max=0, betay_max=0;
  for(i=0;i<linename.Ncell;i++){
    if ( linename.Cell[i]->Beta1 > betax_max ) betax_max = linename.Cell[i]->Beta1 ;
    if ( linename.Cell[i]->Beta2 > betay_max ) betay_max = linename.Cell[i]->Beta2 ;
  }
  cout<<setw(12)<<" Beta_max:" <<setw(20)<<setprecision(15)<<betax_max<<setw(20)<<setprecision(15)<<betay_max<<endl; 
}

void Plot_Twiss(Line & linename, double s1, double s2)
{
 int i;
 fstream fout;
 fout.open("temp_twiss", ios::out);
 for(i=0;i<linename.Ncell;i++)
   if(linename.Cell[i]->S < s2 and linename.Cell[i]->S > s1 )
     fout<<linename.Cell[i]->TYPE<<"  "<<linename.Cell[i]->L<<"  "<<linename.Cell[i]->S<<"    "<<linename.Cell[i]->Beta1<<"    "<<linename.Cell[i]->Beta2<<"  "<<linename.Cell[i]->Etax<<endl;
 fout.close();

 fout.open("temp1.p", ios::out);
 fout<<"set term post color enhanced 20 "<<endl;
 fout<<"set output 'twiss.ps' "<<endl;
 fout<<"set xlabel 's [m]'  " <<endl;
 fout<<"set ylabel '{/Symbol \142}_x , {/Symbol \142}_y [m]  " <<endl;
 fout<<"set y2label 'Dx [m]' "<<endl;
 fout<<"set ytics nomirror"<<endl;
 fout<<"set y2tics"<<endl;
 fout<<"plot 'temp_twiss' u 3:4 tit '{/Symbol \142}_x' w l lt 1  lw 2,\\"<<endl;
 fout<<"    'temp_twiss' u 3:5 tit '{/Symbol \142}_y' w l lt 3 lw 2,\\"<<endl;
 fout<<"    'temp_twiss' u 3:6 axes x1y2 tit 'D_x'  w l lt 2 lw 2"<<endl;
 fout<<"exit"<<endl;
 fout.close();

 system("gnuplot temp1.p");
 system("rm temp1.p");
 system("rm temp_twiss");
}

void Plot_Orbit(Line & linename, double s1, double s2)
{
 int i;
 fstream fout;
 fout.open("temp_orbit", ios::out);
 for(i=0;i<linename.Ncell;i++)
   if(linename.Cell[i]->S < s2 and linename.Cell[i]->S > s1 )
     fout<<linename.Cell[i]->TYPE<<"  "<<linename.Cell[i]->L<<"  "<<linename.Cell[i]->S<<"    "<<linename.Cell[i]->X[0]<<"    "<<linename.Cell[i]->X[1]<<endl;
 fout.close();

 fout.open("temp1.p", ios::out);
 fout<<"set term post color enhanced 20 "<<endl;
 fout<<"set output 'orbit.ps' "<<endl;
 fout<<"set xlabel 's [m]'  " <<endl;
 fout<<"set ylabel 'x_{co}, y_{co}  [m]  " <<endl;
 fout<<"plot 'temp_orbit' u 3:4 tit 'x_{co}' w l lw 2,\\"<<endl;
 fout<<"     'temp_orbit' u 3:5 tit 'y_{co}'  w l lt 2 lw 2"<<endl;
 fout<<"exit"<<endl;
 fout.close();

 system("gnuplot temp1.p");
 system("rm temp1.p");
 system("rm temp_orbit");
}

//---------------------------------------
//   Artifical phase rotator matrix
//--------------------------------------
void Add_Phaser(Line & linename, int loc, const char * name,  double mux, double muy)
{
  int i;
  Element * temp_element;

  double bx, ax, gx,  by, ay, gy, dx, dxp;
  double R11, R12, R21,R22, R33, R34,R43,R44, A,B,C;
  double m[6][6];
  double xco_in[6], xco_out[6];

  for(i=0;i<6;i++) xco_in[i]=0.;
  for(i=0;i<6;i++) xco_out[i]=0.;
  
  bx= linename.Cell[loc-1]->Beta1;
  ax= linename.Cell[loc-1]->Alfa1;
  gx= (1+ax*ax)/bx; 
  
  by= linename.Cell[loc-1]->Beta2;
  ay= linename.Cell[loc-1]->Alfa2;
  gy= (1+ay*ay)/by; 
  
  dx = linename.Cell[loc-1]->Etax;
  dxp= linename.Cell[loc-1]->Etaxp; 
  
  R11=cos(mux) + ax*sin(mux);
  R12=bx*sin(mux);
  R21=-gx*sin(mux);
  R22=cos(mux)-ax*sin(mux);
  
  R33=cos(muy) + ay*sin(muy);
  R34=by*sin(muy);
  R43=-gy*sin(muy);
  R44=cos(muy)-ay*sin(muy);
  
  A=R21*dx-R11*dxp+dxp;
  B=-R12*dxp+R22*dx-dx;
  C=0;
  
  m[0][0]=R11;
  m[0][1]=R12;
  m[0][2]=0;
  m[0][3]=0;
  m[0][4]=0;
  m[0][5]=dx-R11*dx-R12*dxp;
  
  m[1][0]=R21;
  m[1][1]=R22;
  m[1][2]=0;
  m[1][3]=0;
  m[1][4]=0;
  m[1][5]=dxp-R21*dx-R22*dxp;
  
  m[2][0]=0;
  m[2][1]=0;
  m[2][2]=R33;
  m[2][3]=R34;
  m[2][4]=0;
  m[2][5]=0;
  
  m[3][0]=0;
  m[3][1]=0;
  m[3][2]=R43;
  m[3][3]=R44;
  m[3][4]=0;
  m[3][5]=0;
  
  m[4][0]=A;
  m[4][1]=B;
  m[4][2]=0;
  m[4][3]=0;
  m[4][4]=1;
  m[4][5]=C-A*dx-B*dxp;
  
  m[5][0]=0;
  m[5][1]=0;
  m[5][2]=0;
  m[5][3]=0;
  m[5][4]=0;
  m[5][5]=1;
  
  temp_element= new MATRIX(name, 0.0, &m[0][0], xco_in, xco_out);
  linename.Insert(loc, temp_element);
}

//----------------------------------------
//  ORBIT correction
//----------------------------------------
void Correct_Orbit_SVD(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane)
{
  if( n > m ) {
    cout<<" Numer of BPM should be larger than the number of correctors."<<endl;
    exit(0);
  }
  
  int     i,j, k, p;
  double  A[m][n], U[m][m], VT[n][n], AI[n][m];
  double  eps, cut_scale=1000;
  vector<double> temp_vector; 
  vector< vector<double> > V1, U1;
  double read[m], dkick[n];
  
  if( plane ==0 ) {
      for(i=0;i<m;i++) 
	for(j=0;j<n;j++)  A[i][j]=sqrt( linename.Cell[bpm_index[i]]->Beta1 * linename.Cell[kicker_index[j]]->Beta1 ) 
	  *cos( abs(  linename.Cell[bpm_index[i]]->Mu1 - linename.Cell[kicker_index[j]]->Mu1 )*2.0*PI - PI* linename.Tune1 )
	  /2.0/sin( PI * linename.Tune1);
  }
  else{
    for(i=0;i<m;i++) 
      for(j=0;j<n;j++)  A[i][j]=sqrt( linename.Cell[bpm_index[i]]->Beta2 * linename.Cell[kicker_index[j]]->Beta2 ) 
	*cos( abs(  linename.Cell[bpm_index[i]]->Mu2*2*PI    - linename.Cell[kicker_index[j]]->Mu2*2*PI  ) - PI* linename.Tune2 )
	/2.0/sin( PI * linename.Tune2);
  }
  
  eps=1.0e-10;
  i=bmuav(&A[0][0],m,n,&U[0][0], &VT[0][0],eps,m+1);
  if( i<= 0) { cout<<" SVD failed"<<endl;    exit(0);  }
  
  for(i=0;i<n;i++) {
    if( A[i][i] < A[0][0]/cut_scale ) {
      p=i;
      break;
    }
  }
  for(i=0;i<m;i++){
    for(j=0;j<p;j++) temp_vector.push_back(U[i][j]);
    U1.push_back(temp_vector);
    temp_vector.clear();
  }
  for(i=0;i<n;i++){
    for(j=0;j<p;j++) temp_vector.push_back(VT[j][i]);
    V1.push_back(temp_vector);
    temp_vector.clear();
  } 
  for(i=0;i<n;i++) 
    for(j=0;j<m;j++) {
      AI[i][j] = 0.0;
      for(k=0;k<p;k++) AI[i][j]=AI[i][j] + V1[i][k]/ A[k][k] * U1[j][k];
    }

  if (plane ==0 ) {
    for(i=0;i<m;i++) read[i]=linename.Cell[ bpm_index[i]]->X[0]; 
    for(i=0;i<n;i++) {
      dkick[i]=0;
      for(j=0;j<m;j++) dkick[i]= dkick[i] - AI[i][j]*read[j]; 
      linename.Cell[ kicker_index[i] ]->SetP( "HKICK", linename.Cell[ kicker_index[i] ]->GetP("HKICK")+ dkick[i] ); 
    }
  }
  else{
    for(i=0;i<m;i++) read[i]=linename.Cell[ bpm_index[i]]->X[2]; 
    for(i=0;i<n;i++) {
      dkick[i]=0;
      for(j=0;j<m;j++) dkick[i]= dkick[i] - AI[i][j]*read[j]; 
      linename.Cell[ kicker_index[i] ]->SetP( "VKICK", linename.Cell[ kicker_index[i] ]->GetP("VKICK")+ dkick[i]); 
    }
  }
}

void Local_Three_Bump(Line linename, int plane,  const char *corr1,   const char *corr2,   const char *corr3, double kick1)
{
  int loc1, loc2, loc3;
  double beta1, beta2, beta3, scale2, scale3, phi_21,phi_31, phi_32,  kick2, kick3;

  loc1=Get_Index(linename,corr1, 1);
  loc2=Get_Index(linename,corr2, 1);
  loc3=Get_Index(linename,corr3, 1);
  if(plane == 1) {
    beta1 =  linename.Cell[loc1]->Beta1;
    beta2 =  linename.Cell[loc2]->Beta1;
    beta3 =  linename.Cell[loc3]->Beta1;
    phi_21 = ( linename.Cell[loc2]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * PI ; 
    phi_31 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * PI ; 
    phi_32 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc2]->Mu1 )  * 2.0 * PI ; 
    scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
    scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
    kick2=scale2*kick1; 
    kick3=scale3*kick1;
    Set_KL(linename,corr1,"HKICK", kick1);
    Set_KL(linename,corr2,"HKICK", kick2);   
    Set_KL(linename,corr3,"HKICK", kick3);
 }
  else if ( plane == 0 ){
    beta1 =  linename.Cell[loc1]->Beta2;
    beta2 =  linename.Cell[loc2]->Beta2;
    beta3 =  linename.Cell[loc3]->Beta2;
    phi_21 = ( linename.Cell[loc2]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * PI ; 
    phi_31 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * PI ; 
    phi_32 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc2]->Mu2 )  * 2.0 * PI ; 
    kick1 = 0.2e-03;
    scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
    scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
    kick2=scale2*kick1; 
    kick3=scale3*kick1;
    Set_KL(linename,corr1,"VKICK", kick1);
    Set_KL(linename,corr2,"VKICK", kick2);   
    Set_KL(linename,corr3,"VKICK", kick3);
 }
  else{
  }
  if(false){
    Cal_Twiss(linename, 0);
    Print_Twiss(linename,"./twiss");
    exit(0);
  }
}

double RMS_Leakage_Orbit( Line linename, int plane, int i1, int i2 )
{
  int i, count;
  double mean, sum, rms;
  
  sum=0.;
  count=0;
  for(i=0;i<linename.Ncell;i++){
    if( (i-i1)*(i-i2) > 0 )  {
      count++;
      sum += linename.Cell[i]->X[plane*2];
    }
  }
  mean=sum/count;
  
  sum=0.;
  count=0;
  for(i=0;i<linename.Ncell;i++) {
    if( (i-i1)*(i-i2) > 0 ){
      count++;
      sum += (linename.Cell[i]->X[plane*2]-mean) *  (linename.Cell[i]->X[plane*2]-mean) ;
    }
  }
  rms=sqrt ( sum / count );
  
  return rms; 
}

void Correct_Orbit_SlidingBump1(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane)
{
  int i,j;
  double  reading , exceed=0.5e-03;
  int     loc1, loc2, loc3, loc_bpm;
  double  s1, s2, s3, phi_31, phi_21, phi_32, dphi_bpm, beta1,  beta2,  beta3, kick1, kick2, kick3, scale2, scale3;
  int     dir;
  double distance;
  
  for(i=0;i<kicker_index.size();i++ ){

    if(i== kicker_index.size()-1 ){
      loc1=kicker_index[kicker_index.size()-1];
      loc2=kicker_index[0];
      loc3=kicker_index[1];
    }
    else if (i== kicker_index.size()-2 ){
      loc1=kicker_index[kicker_index.size()-2];
      loc2=kicker_index[kicker_index.size()-1];
      loc3=kicker_index[0];
    }
    else{
      loc1=kicker_index[i];
      loc2=kicker_index[i+1];
      loc3=kicker_index[i+2];
    }
    s1 = linename.Cell[loc1]->S;
    s2 = linename.Cell[loc2]->S;
    s3 = linename.Cell[loc3]->S;

    if( s2 > s1 and s2 > s3 and s1 > s3  ) {
      loc_bpm=  loc2 ;
      dir=0;
      if ( plane ==0 ) {
	dphi_bpm=( linename.Tune1 - linename.Cell[loc_bpm]->Mu1 +  linename.Cell[loc3]->Mu1 )   * 2 * PI;
      }
      else{
	dphi_bpm=( linename.Tune2 - linename.Cell[loc_bpm]->Mu2 +  linename.Cell[loc3]->Mu2 )   * 2 * PI;
      }
    }
    else if ( s1> s2 and s1 > s3 and s3 > s2) {
      loc_bpm= loc2 ;
      dir=1;     
      if ( plane ==0 ) {
	dphi_bpm=( linename.Cell[loc_bpm]->Mu1 + linename.Tune1 -  linename.Cell[loc1]->Mu1 )   * 2 * PI;
      }
      else{
	dphi_bpm=( linename.Cell[loc_bpm]->Mu2 + linename.Tune2 -  linename.Cell[loc1]->Mu2 )   * 2 * PI;
      }
    }
    else if ( s3 > s2 and s2 > s1 ) {
      distance=100;
      loc_bpm=0;
      for(j=0;j<bpm_index.size();j++){
	if( abs( linename.Cell[ bpm_index[j] ]->S - linename.Cell[loc2]->S )  < distance ){
	  distance=  abs( linename.Cell[ bpm_index[j] ]->S - linename.Cell[loc2]->S );
	  loc_bpm   =   bpm_index[j];
	}
      }
      if(  linename.Cell[loc_bpm ]->S > linename.Cell[loc2]->S ) {
	dir=0; 
	if (plane == 0 ) {dphi_bpm= linename.Cell[loc3]->Mu1*2*PI - linename.Cell[ loc_bpm ]->Mu1 * 2.0 * PI;}
	else{ dphi_bpm= linename.Cell[loc3]->Mu2*2*PI - linename.Cell[ loc_bpm ]->Mu2 * 2.0 * PI;}
      }
      else{  
	dir=1; 
	if(plane==0 ){ dphi_bpm=  linename.Cell[ loc_bpm ]->Mu1 * 2.0 * PI - linename.Cell[loc1]->Mu1*2*PI; }
        else{ dphi_bpm=  linename.Cell[ loc_bpm ]->Mu2 * 2.0 * PI - linename.Cell[loc1]->Mu2*2*PI; }
      }
    }
    else{ cout<<"something wrong in sliding bump."<<endl; exit(1); }

    if ( plane == 0 ) {
      reading= linename.Cell[ loc_bpm ]->X[0];
      if(abs(reading) > exceed ) { 
	beta1= linename.Cell[loc1]->Beta1 ; 
	beta2= linename.Cell[loc2]->Beta1 ;  
	beta3= linename.Cell[loc3]->Beta1 ; 
	phi_21 = ( linename.Cell[loc2]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * PI ; 
	phi_31 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * PI ; 
	phi_32 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc2]->Mu1 )  * 2.0 * PI ; 
        if ( phi_21 < 0 ) phi_21 = phi_21 + linename.Tune1*2*PI;
        if ( phi_31 < 0 ) phi_31 = phi_31 + linename.Tune1*2*PI;
        if ( phi_32 < 0 ) phi_32 = phi_32 + linename.Tune1*2*PI;
	if( dir ==1 ) {
	  kick1=-reading
	    / sqrt(beta1* linename.Cell[ loc_bpm ]->Beta1 ) / sin( dphi_bpm  );
	  scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	  scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	  kick2=scale2*kick1;  kick3=scale3 * kick1;
	}
	else{
	  kick3=-reading 
	    / sqrt(beta3* linename.Cell[ loc_bpm ]->Beta1 ) / sin( dphi_bpm ); 
	  scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	  scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	  kick1=kick3/scale3; kick2=scale2*kick1;
	}
	linename.Cell[loc1]->SetP("HKICK", linename.Cell[loc1]->GetP("HKICK")+ kick1);
	linename.Cell[loc2]->SetP("HKICK", linename.Cell[loc2]->GetP("HKICK")+ kick2);
	linename.Cell[loc3]->SetP("HKICK", linename.Cell[loc3]->GetP("HKICK")+ kick3);
      }
    }
    else{
      reading= linename.Cell[ loc_bpm ]->X[2];
      if(abs(reading) > exceed ) { 
	beta1= linename.Cell[loc1]->Beta2 ; 
	beta2= linename.Cell[loc2]->Beta2 ;  
	beta3= linename.Cell[loc3]->Beta2 ; 
	phi_21 = ( linename.Cell[loc2]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * PI ; 
	phi_31 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * PI ; 
	phi_32 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc2]->Mu2 )  * 2.0 * PI ; 
        if ( phi_21 < 0 ) phi_21 = phi_21 + linename.Tune2*2*PI;
        if ( phi_31 < 0 ) phi_31 = phi_31 + linename.Tune2*2*PI;
        if ( phi_32 < 0 ) phi_32 = phi_32 + linename.Tune2*2*PI;
	if( dir ==1 ) {
	  kick1=-reading
	    / sqrt(beta1* linename.Cell[ loc_bpm ]->Beta2 ) / sin( dphi_bpm  );
	  scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	  scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	  kick2=scale2*kick1;  kick3=scale3 * kick1;
	}
	else{
	  kick3=-reading 
	    / sqrt(beta3* linename.Cell[ loc_bpm ]->Beta2 ) / sin( dphi_bpm ); 
	  scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	  scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	  kick1=kick3/scale3; kick2=scale2*kick1;
	}
	linename.Cell[loc1]->SetP("VKICK", linename.Cell[loc1]->GetP("VKICK")+ kick1);
	linename.Cell[loc2]->SetP("VKICK", linename.Cell[loc2]->GetP("VKICK")+ kick2);
	linename.Cell[loc3]->SetP("VKICK", linename.Cell[loc3]->GetP("VKICK")+ kick3);
      }
    }
  }
}

void Correct_Orbit_SlidingBump2(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane)
{
  int i;
  double  reading , exceed=0.5e-03;
  int     loc1, loc2, loc3, loc_bpm;
  double  phi_31, phi_21, phi_32, dphi_bpm, beta1,  beta2,  beta3, kick1, kick2, kick3, scale2, scale3;
  
  for(i=0;i<kicker_index.size();i++ ){

    if(i== kicker_index.size()-1 ){
      loc1=kicker_index[kicker_index.size()-1];
      loc2=kicker_index[0];
      loc3=kicker_index[1];
    }
    else if (i== kicker_index.size()-2 ){
      loc1=kicker_index[kicker_index.size()-2];
      loc2=kicker_index[kicker_index.size()-1];
      loc3=kicker_index[0];
    }
    else{
      loc1=kicker_index[i];
      loc2=kicker_index[i+1];
      loc3=kicker_index[i+2];
    }
    
    loc_bpm= loc2;

    if ( plane == 0 ) {
      reading= linename.Cell[ loc_bpm ]->X[0];
      if(abs(reading) > exceed ) { 
	beta1= linename.Cell[loc1]->Beta1 ; 
	beta2= linename.Cell[loc2]->Beta1 ;  
	beta3= linename.Cell[loc3]->Beta1 ; 
	phi_21 = ( linename.Cell[loc2]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * PI ; 
	phi_31 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * PI ; 
	phi_32 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc2]->Mu1 )  * 2.0 * PI ; 
	dphi_bpm=( linename.Cell[loc_bpm]->Mu1 -  linename.Cell[loc1]->Mu1 )   * 2 * PI;
        if ( phi_21 < 0 ) phi_21 = phi_21 + linename.Tune1*2*PI;
        if ( phi_31 < 0 ) phi_31 = phi_31 + linename.Tune1*2*PI;
        if ( phi_32 < 0 ) phi_32 = phi_32 + linename.Tune1*2*PI;
	if ( dphi_bpm < 0 ) dphi_bpm = dphi_bpm + linename.Tune1*2*PI;
	kick1=-reading
	  / sqrt(beta1* linename.Cell[ loc_bpm ]->Beta1 ) / sin( dphi_bpm  );
	scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	kick2=scale2*kick1;  kick3=scale3 * kick1;
	linename.Cell[loc1]->SetP("HKICK", linename.Cell[loc1]->GetP("HKICK")+ kick1);
	linename.Cell[loc2]->SetP("HKICK", linename.Cell[loc2]->GetP("HKICK")+ kick2);
	linename.Cell[loc3]->SetP("HKICK", linename.Cell[loc3]->GetP("HKICK")+ kick3);
      }
    }
    else{
      reading= linename.Cell[ loc_bpm ]->X[2];
      if(abs(reading) > exceed ) { 
	beta1= linename.Cell[loc1]->Beta2 ; 
	beta2= linename.Cell[loc2]->Beta2 ;  
	beta3= linename.Cell[loc3]->Beta2 ; 
	phi_21 = ( linename.Cell[loc2]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * PI ; 
	phi_31 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * PI ; 
	phi_32 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc2]->Mu2 )  * 2.0 * PI ; 
	dphi_bpm=( linename.Cell[loc_bpm]->Mu2  -  linename.Cell[loc1]->Mu2 )   * 2 * PI;
        if ( phi_21 < 0 ) phi_21 = phi_21 + linename.Tune2*2*PI;
        if ( phi_31 < 0 ) phi_31 = phi_31 + linename.Tune2*2*PI;
        if ( phi_32 < 0 ) phi_32 = phi_32 + linename.Tune2*2*PI;
	if ( dphi_bpm < 0 ) dphi_bpm = dphi_bpm + linename.Tune2*2*PI;
	kick1=-reading
	  / sqrt(beta1* linename.Cell[ loc_bpm ]->Beta2 ) / sin( dphi_bpm  );
	scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	kick2=scale2*kick1;  kick3=scale3 * kick1;
	linename.Cell[loc1]->SetP("VKICK", linename.Cell[loc1]->GetP("VKICK")+ kick1);
	linename.Cell[loc2]->SetP("VKICK", linename.Cell[loc2]->GetP("VKICK")+ kick2);
	linename.Cell[loc3]->SetP("VKICK", linename.Cell[loc3]->GetP("VKICK")+ kick3);
      }
    }
  }
}

void Orbit_Status( Line linename, vector<int> bpm_index, int plane, double &orbit_mean, double & orbit_rms )
{
  int i;
  double sum;
  sum=0.;
  for(i=0;i<bpm_index.size();i++) sum += linename.Cell[ bpm_index[i] ]->X[plane*2];
  orbit_mean=sum/bpm_index.size();

  sum=0.;
  for(i=0;i<bpm_index.size();i++) sum += (linename.Cell[ bpm_index[i] ]->X[plane*2]-orbit_mean) *  (linename.Cell[ bpm_index[i] ]->X[plane*2]-orbit_mean) ;
  orbit_rms=sqrt ( sum / bpm_index.size() );
}

//==============================================================
//
//                Tracking 
//
//===============================================================

void Track(Line & linename, double x[6], int nturn, int & stable, int & lost_turn, int & lost_post)
{
  int j;
  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;    lost_turn= 0;    lost_post = 0;    return; }
  
  //----now we do tracking
  for(GP.turn=0; GP.turn < nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy  ) { 
      	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      }
      //each element to output
      // cout<<j<<"   ";
      //for(k=0;k<6;k++) cout<<x[k]<<" ";
      //cout<<endl;
    } 
    // each turn to output
    // for(k=0;k<6;k++) cout<<x[k]<<" ";
    // cout<<endl;
  }
}

void Track_tbt(Line & linename, double x[6], int nturn, double  x_tbt[], int & stable, int & lost_turn, int & lost_post)
{
  int j;

  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;     lost_turn= 0;     lost_post = 0;     return;   }

  //----now we track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    for(j=0;j<6;j++) x_tbt[GP.turn*6+j]=x[j];
    
    //for(j=0;j<6;j++) cout<<x[j]<<"  ";
    //cout<<endl;
    
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy  ) { 
	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      } } }
}

void Track_tbt(Line & linename, double x[6], int nturn, double  x_tbt[], int &bpm_index, int & stable, int & lost_turn, int & lost_post)
{
  int j,k;

  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;     lost_turn= 0;     lost_post = 0;     return;   }

  //----now we track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);
      if(j == bpm_index ){
	for(k=0;k<6;k++) x_tbt[GP.turn*6+k]=x[k];
	//for(k=0;k<6;k++) cout<<x[k]<<"  ";
	//cout<<endl;
      }
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy  ) { 
	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      } } }
}

void Cal_Tune_Num(Line & linename, double deltap0)
{
  int k;
  double x[6];

  int Nturn=1024;
  int stable = 1, lost_turn=0, lost_post=0;
  double  x_tbt[Nturn*6], xtemp[1024];
  double tunex1, tuney1;
  
  Cal_Orbit_Num(linename, deltap0);

  x[0]=linename.Cell[linename.Ncell-1]->X[0] ;
  x[1]=linename.Cell[linename.Ncell-1]->X[1] ;
  x[2]=linename.Cell[linename.Ncell-1]->X[2] ;
  x[3]=linename.Cell[linename.Ncell-1]->X[3] ;
  x[4]=0. ;
  x[5]= deltap0;
  
  if(x[0] < 1.0e-9 ) x[0]=x[0]+1.0e-9;
  if(x[2] < 1.0e-9 ) x[2]=x[2]+1.0e-9; 

  stable = 1; lost_turn=0; lost_post=0;
  Track_tbt(linename, x, Nturn, x_tbt, stable, lost_turn, lost_post);

  if(stable==1) {
    for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+0];
    FineTuneFinder(1024, xtemp, tunex1);
    for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+2];
    FineTuneFinder(1024, xtemp, tuney1);
  }
  linename.Tune1= tunex1;
  linename.Tune2= tuney1;
}

//----------------------------------------------------------------
//Tracking without  linename.Cell[i]->Pass() to speed up tracking
//----------------------------------------------------------------
void Prepare_Track_Fast(Line & linename)
//  this function only accepts the allowed elements and change others to drift
//  this function accepts thick lenses: SBEND, QUAD, SEXT, MULT, ELENS, GMULT
{
  int i,j,k;
  //-----------do some statistics
  nelement=linename.Ncell;
  Ndrift=0;
  Nsbend=0;
  Nquad=0;
  Nsext=0;
  Nmult=0;
  Nrf=0;
  Nbb=0;
  Nelens=0;
  Ncooling=0;
  Nacmult =0;
  Nmatrix =0;
  Nkicker =0;
  Nhkicker=0;
  Nvkicker=0;
  Ndiffuse=0;
  Ngmult  =0;
  Nlrbb   =0;

  for(i=0; i<linename.Ncell; i++){
    if(linename.Cell[i]->TYPE==string("DRIFT") ) {
      Ndrift++;
    }
    else if(linename.Cell[i]->TYPE==string("SBEND") ) {
      Nsbend++;
    }
    else if(linename.Cell[i]->TYPE==string("QUAD") ){
      Nquad++;
    }
    else if(linename.Cell[i]->TYPE==string("SEXT") ){
      Nsext++;
    }
    else if(linename.Cell[i]->TYPE==string("MULT") ){
      Nmult++;
    }
    else if(linename.Cell[i]->TYPE==string("BEAMBEAM") ){
      Nbb++;
    }
    else if(linename.Cell[i]->TYPE==string("ELENS") ){
      Nelens++;
    }
    else if(linename.Cell[i]->TYPE==string("RFCAV") ){
      Nrf++;
    }
    else if(linename.Cell[i]->TYPE==string("COOLING") ){
      Ncooling++;
    }
    else if(linename.Cell[i]->TYPE==string("ACMULT") ){
      Nacmult++;
    }
    else if(linename.Cell[i]->TYPE==string("MATRIX") ){
      Nmatrix++;
    }
    else if(linename.Cell[i]->TYPE==string("KICKER") ){
      Nkicker++;
    }
    else if(linename.Cell[i]->TYPE==string("HKICKER") ){
      Nhkicker++;
    }
    else if(linename.Cell[i]->TYPE==string("VKICKER") ){
      Nvkicker++;
    }
    else if(linename.Cell[i]->TYPE==string("DIFFUSE") ){
      Ndiffuse++;
    }
    else if(linename.Cell[i]->TYPE==string("GMULT") ){
      Ngmult++;
    }
    else if(linename.Cell[i]->TYPE==string("LRBB") ){
      Nlrbb++;
    }
    else{ 
      Ndrift++;
    }
  }
  
  if( true ) {
    cout<<"---------------------------------"<<endl;
    cout<<"Statistics before Track_Fast(): "<<endl;
    cout<<"There are "<<Ndrift<< " DRIFT"<<endl;
    cout<<"There are "<<Nsbend<< " SEBDN"<<endl;
    cout<<"There are "<<Nquad<< "  QUAD"<<endl;
    cout<<"There are "<<Nsext<< "  SEXT"<<endl;
    cout<<"There are "<<Nmult<< "  MULT"<<endl;
    cout<<"There are "<<Nbb<< "  BEAMBEAM"<<endl;
    cout<<"There are "<<Nelens<< "  ELENS"<<endl;
    cout<<"There are "<<Nrf<< "  RFCAV"<<endl;
    cout<<"There are "<<Ncooling<< "  COOLING"<<endl;
    cout<<"There are "<<Nacmult<< "  ACMULT"<<endl;
    cout<<"There are "<<Nmatrix<< "  MATRIX"<<endl;
    cout<<"There are "<<Nkicker<< "  KICKER"<<endl;
    cout<<"There are "<<Nhkicker<< "  HKICKER"<<endl;
    cout<<"There are "<<Nvkicker<< "  VKICKER"<<endl;
    cout<<"There are "<<Ndiffuse<< "  DIFFUSE"<<endl;
    cout<<"There are "<<Ngmult<< "  GMULT"<<endl;
    cout<<"There are "<<Nlrbb<< "  LRBB"<<endl;
    cout<<"Totally "<<linename.Ncell<<"  elements"<<endl;
  }
 
  //--- store strengths
  for(i=0; i<nelement; i++){
    
    nelement_dx.push_back(linename.Cell[i]->DX );
    nelement_dy.push_back(linename.Cell[i]->DY );
    nelement_dt.push_back(linename.Cell[i]->DT );
    nelement_APx.push_back(linename.Cell[i]->APx );
    nelement_APy.push_back(linename.Cell[i]->APy ); 

    if(linename.Cell[i]->TYPE==string("DRIFT") ) {
      type.push_back( 0 );
      drift_l.push_back(linename.Cell[i]->L);
    }
    else if(linename.Cell[i]->TYPE==string("SBEND") ) {
      type.push_back( 1  );
      sbend_l.push_back( linename.Cell[i]->L );
      sbend_angle.push_back( linename.Cell[i]->GetP("ANGLE") );
      sbend_nint.push_back( int(linename.Cell[i]->GetP("Nint")) );
      sbend_e1.push_back( linename.Cell[i]->GetP("E1") );
      sbend_e2.push_back( linename.Cell[i]->GetP("E2") );
    }
    else if(linename.Cell[i]->TYPE==string("QUAD") ){
      type.push_back( 2 );
      quad_l.push_back( linename.Cell[i]->L );
      quad_k1l.push_back( linename.Cell[i]->GetP("K1L") );
      quad_nint.push_back( int(linename.Cell[i]->GetP("Nint")) );
    }
    else if(linename.Cell[i]->TYPE==string("SEXT") ){
      type.push_back( 3 );
      sext_l.push_back(linename.Cell[i]->L);
      sext_nint.push_back( int(linename.Cell[i]->GetP("Nint")));
      sext_k2l.push_back(linename.Cell[i]->GetP("K2L"));
    }
    else if(linename.Cell[i]->TYPE==string("MULT") ){
      type.push_back( 4 );
      mult_l.push_back(linename.Cell[i]->L);
      mult_nint.push_back(int(linename.Cell[i]->GetP("Nint")));
      mult_norder.push_back( int( linename.Cell[i]->GetP("Norder") )  );
      for(j=0;j<11;j++){
        char name1[125], name2[125];
        sprintf(name1, "K%dL",j);
        sprintf(name2, "K%dSL",j);
        mult_knl.push_back(linename.Cell[i]->GetP(name1));
        mult_knsl.push_back(linename.Cell[i]->GetP(name2));
      }
    }
    else if(linename.Cell[i]->TYPE==string("BEAMBEAM") ){
      type.push_back( 5 );
      bb_np.push_back( linename.Cell[i]->GetP("NP") );
      bb_treat.push_back( int( linename.Cell[i]->GetP("TREATMENT") ) );
      bb_nslice.push_back( int( linename.Cell[i]->GetP("NSLICE") ) );
      bb_emitx.push_back( linename.Cell[i]->GetP("EMITX") );
      bb_emity.push_back( linename.Cell[i]->GetP("EMITY") );
      bb_sigmal.push_back( linename.Cell[i]->GetP("SIGMAL") );
      bb_betax.push_back( linename.Cell[i]->GetP("BETAX") );
      bb_alfax.push_back( linename.Cell[i]->GetP("ALFAX") );
      bb_betay.push_back( linename.Cell[i]->GetP("BETAY") );
      bb_alfay.push_back( linename.Cell[i]->GetP("ALFAY") );
    }
    else if(linename.Cell[i]->TYPE==string("ELENS") ){
      type.push_back( 6  );
      elens_ne.push_back( linename.Cell[i]->GetP("NE") );
      elens_l.push_back( linename.Cell[i]->L );     
      elens_nslice.push_back( int( linename.Cell[i]->GetP("NSLICE") ) );
      elens_sigmax.push_back( linename.Cell[i]->GetP("SIGMAX") );
      elens_sigmay.push_back( linename.Cell[i]->GetP("SIGMAY") );
      elens_betae.push_back( linename.Cell[i]->GetP("BETAE") );
    }
    else if(linename.Cell[i]->TYPE==string("RFCAV") ){
      type.push_back( 7 );
      rf_l.push_back( linename.Cell[i]->L);
      rf_vrf.push_back( linename.Cell[i]->GetP("VRF") );
      rf_frf.push_back( linename.Cell[i]->GetP("FRF") );
      rf_phi0.push_back( linename.Cell[i]->GetP("PHASE0") );
    }
    else if(linename.Cell[i]->TYPE==string("COOLING") ){
      type.push_back( 8 );
      cooling_l.push_back( linename.Cell[i]->L );
      cooling_alpha.push_back( linename.Cell[i]->GetP("ALPHA") );
    }
    else if(linename.Cell[i]->TYPE==string("ACMULT") ){
      type.push_back( 9 );
      acmult_l.push_back( linename.Cell[i]->L );     
      acmult_norder.push_back( int(linename.Cell[i]->GetP("Norder") ));
      acmult_tturns.push_back( int(linename.Cell[i]->GetP("Tturns") ));
      acmult_kl.push_back( linename.Cell[i]->GetP("KL") );
      acmult_phi0.push_back( linename.Cell[i]->GetP("PHI0") );
    }
   else if(linename.Cell[i]->TYPE==string("MATRIX") ){
      type.push_back( 10 );
      matrix_l.push_back( linename.Cell[i]->L); 
      for(j=0;j<6;j++)
	for(k=0;k<6;k++)
	  {
	    char name1[125];
	    sprintf(name1, "M%d%d",j+1,k+1);
	    matrix_m66.push_back( linename.Cell[i]->GetP(name1) ); 
	  }
      
      matrix_xco_in.push_back(  linename.Cell[i]->GetP("XCO_IN_X") );
      matrix_xco_in.push_back(  linename.Cell[i]->GetP("XCO_IN_PX") );
      matrix_xco_in.push_back(  linename.Cell[i]->GetP("XCO_IN_Y") );
      matrix_xco_in.push_back(  linename.Cell[i]->GetP("XCO_IN_PY") );
      matrix_xco_in.push_back(  linename.Cell[i]->GetP("XCO_IN_Z") );
      matrix_xco_in.push_back(  linename.Cell[i]->GetP("XCO_IN_DELTA") );

      matrix_xco_out.push_back(  linename.Cell[i]->GetP("XCO_OUT_X") );
      matrix_xco_out.push_back(  linename.Cell[i]->GetP("XCO_OUT_PX") );
      matrix_xco_out.push_back(  linename.Cell[i]->GetP("XCO_OUT_Y") );
      matrix_xco_out.push_back(  linename.Cell[i]->GetP("XCO_OUT_PY") );
      matrix_xco_out.push_back(  linename.Cell[i]->GetP("XCO_OUT_Z") );
      matrix_xco_out.push_back(  linename.Cell[i]->GetP("XCO_OUT_DELTA") );
    }
    else if(linename.Cell[i]->TYPE==string("KICKER") ){
      type.push_back( 11 );
      kicker_hkick.push_back( linename.Cell[i]->GetP("HKICK") );
      kicker_vkick.push_back( linename.Cell[i]->GetP("VKICK") );
    }
    else if(linename.Cell[i]->TYPE==string("HKICKER") ){
      type.push_back( 12 );
      hkicker_hkick.push_back( linename.Cell[i]->GetP("HKICK") );
    }
    else if(linename.Cell[i]->TYPE==string("VKICKER") ){
      type.push_back( 13 );
      vkicker_vkick.push_back( linename.Cell[i]->GetP("VKICK") );
    }
    else if(linename.Cell[i]->TYPE==string("DIFFUSE") ){
      type.push_back( 14 );
      diffuse_l.push_back( linename.Cell[i]->L);
      diffuse_diff_x.push_back( linename.Cell[i]->GetP("DIFF_X") );
      diffuse_diff_y.push_back( linename.Cell[i]->GetP("DIFF_Y") );
      diffuse_diff_delta.push_back( linename.Cell[i]->GetP("DIFF_DELTA") );
    }
    else if(linename.Cell[i]->TYPE==string("GMULT") ){
      type.push_back( 15 );
      gmult_l.push_back(linename.Cell[i]->L); 
      gmult_nint.push_back( int(linename.Cell[i]->GetP("Nint")) ); 
      gmult_norder.push_back( int(linename.Cell[i]->GetP("Norder")) ); 
      gmult_l.push_back( linename.Cell[i]->L );
      gmult_angle.push_back( linename.Cell[i]->GetP("ANGLE") );
      gmult_e1.push_back( linename.Cell[i]->GetP("E1") );
      gmult_e2.push_back( linename.Cell[i]->GetP("E2") );
      for(j=0;j<11;j++){
        char name1[125], name2[125];
        sprintf(name1, "K%dL",j);
        sprintf(name2, "K%dSL",j);
        gmult_knl.push_back(linename.Cell[i]->GetP(name1));
        gmult_knsl.push_back(linename.Cell[i]->GetP(name2));
      }
    }
    else if(linename.Cell[i]->TYPE==string("LRBB") ){
      type.push_back( 16 );
      lrbb_np.push_back( linename.Cell[i]->GetP("NP") );
      lrbb_sepx.push_back( linename.Cell[i]->GetP("SEPX") );
      lrbb_sepy.push_back( linename.Cell[i]->GetP("SEPY") );
      lrbb_sigmax.push_back( linename.Cell[i]->GetP("SIGMAX") );
      lrbb_sigmay.push_back( linename.Cell[i]->GetP("SIGMAY") );
      lrbb_kickx0.push_back( linename.Cell[i]->GetP("KICKX0") );
      lrbb_kicky0.push_back( linename.Cell[i]->GetP("KICKY0") );
    }
    else{ 
      //-----I hope no more Warning here. Should be prevented by Make_Thin()
      cout<<" Warning: "<<linename.Cell[i]->NAME<<" :  "<<linename.Cell[i]->TYPE<<"  at  "<<linename.Cell[i]->S<<"  converted to DRIFT."<<endl; 
      type.push_back( 0 );
      double temp=linename.Cell[i]->L;
      drift_l.push_back( temp  );
    }
  }
}

void Track_Fast(Line & linename, double x[6], int nturn, int & stable, int & lost_turn, int & lost_post, double & sum_x2, double & sum_y2, double & sum_z2)
{
  int i,j,k;
  double  Dpx, Dpy;
  double  knl[11], knsl[11];
  double  m66[36], xco_in[6], xco_out[6];
  double  x1[6];

  sum_x2=0.;  sum_y2=0.;  sum_z2=0.;

  //-----quick check before track
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;      lost_turn= 0;       lost_post = 0;      return ;   }
  
  //---now we do track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    
    Ndrift=0;
    Nsbend=0;
    Nquad=0;
    Nsext=0;
    Nmult=0;
    Nrf=0;
    Nbb=0;
    Nelens=0;
    Ncooling=0;
    Nacmult=0;
    Nmatrix=0;   
    Nkicker =0;
    Nhkicker=0;
    Nvkicker=0;
    Ndiffuse=0;
    Ngmult=0;
    Nlrbb=0;
    
    for(i=0; i<nelement; i++){
      
      //--can be used for each element or some type of element
      //GtoL(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 
      
      if(type[i]==0 ) {
	DRIFT_Pass(x, drift_l[Ndrift]);
	Ndrift++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==1 ) {
	SBEND_Pass(x,sbend_l[Nsbend],sbend_nint[Nsbend], sbend_angle[Nsbend],sbend_e1[Nsbend],sbend_e2[Nsbend]); 
	Nsbend++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==2 ){
	QUAD_Pass(x,  quad_l[Nquad], quad_nint[Nquad], quad_k1l[Nquad], 0.); 
	Nquad++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==3 ){
	SEXT_Pass(x, sext_l[Nsext], sext_nint[Nsext], sext_k2l[Nsext], 0.);
	Nsext++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==4 ){
	for(j=0;j<11;j++) {
	  knl[j]=mult_knl[Nmult*11+j]; 
	  knsl[j]=mult_knsl[Nmult*11+j];
	}
	MULT_Pass(x, mult_l[Nmult], mult_nint[Nmult], mult_norder[Nmult], knl, knsl);
	Nmult++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==5 ){
	BEAMBEAM_Pass(x,bb_treat[Nbb], bb_np[Nbb], bb_sigmal[Nbb], bb_nslice[Nbb], bb_emitx[Nbb], 
		      bb_betax[Nbb], bb_alfax[Nbb], bb_emity[Nbb], bb_betay[Nbb], bb_alfay[Nbb]);
	Nbb++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==6 ){
	GtoL(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 
        ELENS_Pass(x, elens_ne[Nelens], elens_l[Nelens], elens_betae[Nelens],
		   elens_nslice[Nelens], elens_sigmax[Nelens], elens_sigmay[Nelens]);
	LtoG(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 
	Nelens++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if (type[i]==7 ){
        //RFCAV_Pass(x, rf_l[Nrf], rf_vrf[Nrf], rf_frf[Nrf], 0.);
	x[5] =x[5] + (rf_vrf[Nrf]*1.0/GP.energy)*sin(2.0*PI*rf_frf[Nrf]*x[z_]/3.0e8);
	Nrf++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if (type[i]==8 ){
        //COOLING_Pass(x, cooling_L[Ncooling], cooling_alpha[Ncooling]);
	for(j=0;j<4;j++) x[j]= (1.- cooling_alpha[Ncooling] )* x[j];
	Ncooling++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if (type[i]==9 ){
	ACMULT_Pass(x,  acmult_l[Nacmult], acmult_norder[Nacmult], acmult_kl[Nacmult], acmult_tturns[Nacmult], acmult_phi0[Nacmult]);
	Nacmult++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==10 ){
	for(j=0;j<36;j++) {
	  m66[j]=matrix_m66[Nmatrix*36+j]; 
	}
	for(j=0;j<6;j++) xco_in[j] = matrix_xco_in[Nmatrix*6+j];
	for(j=0;j<6;j++) xco_out[j]= matrix_xco_out[Nmatrix*6+j];
	
	for(j=0;j<6;j++) x[j]=x[j]-xco_in[j];
	for(j=0;j<6;j++)    x1[j]=x[j];
	for(j=0;j<6;j++) {
	  x[j]=0.0;
	  for(k=0;k<6;k++)  x[j] = x[j] + m66[j*6+k] * x1[k];
	}  
	for(j=0;j<6;j++) x[j]=x[j]+xco_out[j];
	
	Nmatrix++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==11 ){
	x[1]=x[1]+kicker_hkick[Nkicker];   
	x[3]=x[3]+kicker_vkick[Nkicker];
	Nkicker++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==12 ){
	x[1]=x[1]+hkicker_hkick[Nhkicker];   
	Nhkicker++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==13 ){
	x[1]=x[1]+vkicker_vkick[Nvkicker];   
	Nvkicker++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==14 ){
	x[1]= x[1]+ diffuse_diff_x[Ndiffuse] * rnd(seed); 
	x[3]= x[1]+ diffuse_diff_y[Ndiffuse] * rnd(seed); 
	x[5]= x[1]+ diffuse_diff_delta[Ndiffuse] * rnd(seed); 
	Ndiffuse++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==15 ){
	for(j=0;j<11;j++) {
	  knl[j] =gmult_knl[Ngmult*11+j]; 
	  knsl[j]=gmult_knsl[Ngmult*11+j];
	}
	GMULT_Pass(x, gmult_l[Ngmult], gmult_nint[Ngmult], gmult_angle[Ngmult], gmult_norder[Ngmult],knl,knsl,
                      gmult_e1[Ngmult],gmult_e2[Ngmult]);
	Ngmult++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==16 ){
        LRBB_Pass(x, GP.gamma, lrbb_np[Nlrbb],lrbb_sepx[Nlrbb],lrbb_sepy[Nlrbb],lrbb_sigmax[Nlrbb],lrbb_sigmay[Nlrbb],lrbb_kickx0[Nlrbb],lrbb_kicky0[Nlrbb] );
	Nlrbb++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else{ }

     //LtoG(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 

       //print out each element 
      // cout<<i<<"   ";
      // for(j=0;j<6;j++) cout<<x[j]<<" ";
      // cout<<endl;
    }
    sum_x2 += x[0]*x[0];
    sum_y2 += x[2]*x[2]; 
    sum_z2 += x[4]*x[4]; 

    // print out turn-by-turn data
    // for(j=0;j<6;j++) cout<<x[j]<<" ";
    // cout<<endl;
  }
}

//--------------------------------------------------
//      Nonlinear tools: FOOTPRINT, FMA
//--------------------------------------------------
void Track_tbt_FMA( Line & linename, double deltap0, double sigmax0, double sigmay0 ) // I prefer RF off for this 
{
  int i,j,k;
  int Nturn=2048;
  double sigma_step=0.05;
  int sigma0=0, sigma1=6,  nsigma= int((sigma1-sigma0 ) *1.0 / sigma_step  )+1  ;  

  double x[6];
  double  x_tbt[Nturn*6];
  int stable = 1, lost_turn=0, lost_post=0;
  
  double xtemp[1024];
  double tunex1, tunex2, tuney1, tuney2;

  fstream  f1;
  f1.open("./FMA-output.dat", ios::out);
  f1.close();
  
  Cal_Orbit_Num(linename, deltap0);

  for(i=0; i<nsigma;  i++) {
    for(j=0; j<nsigma; j++ ) {
      if(  sqrt( 1.0*i*i +1.0*j*j) <= (6./sigma_step) ) {
        
        x[0] = i*sigma_step*sigmax0;
        x[1] = 0.;
        x[2] = j*sigma_step*sigmay0;
        x[3] = 0.;
        x[4] = 0.;
        x[5] = deltap0;
	x[0]= x[0] + linename.Cell[linename.Ncell-1]->X[0]  ;
	x[1]= x[1] + linename.Cell[linename.Ncell-1]->X[1]  ;
	x[2]= x[2] + linename.Cell[linename.Ncell-1]->X[2]  ;
	x[3]= x[3] + linename.Cell[linename.Ncell-1]->X[3]  ; 
        
        if(x[0] < 1.0e-09 ) x[0]=x[0]+1.0e-09;
        if(x[2] < 1.0e-09 ) x[2]=x[2]+1.0e-09; 

        stable = 1; lost_turn=0; lost_post=0;
        Track_tbt(linename, x, Nturn, x_tbt, stable, lost_turn, lost_post);
        if(stable==1) {
          for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+0];
          FineTuneFinder(1024,xtemp,tunex1);
          for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+2];
          FineTuneFinder(1024,xtemp,tuney1);
          for(k=0;k<1024;k++) xtemp[k]=x_tbt[(k+1023)*6+0];
          FineTuneFinder(1024,xtemp,tunex2);
          for(k=0;k<1024;k++) xtemp[k]=x_tbt[(k+1023)*6+2];
          FineTuneFinder(1024,xtemp,tuney2);
          f1.open("./FMA-output.dat",  ios::out | ios::app);
          f1<<i*sigma_step<<" "<<j*sigma_step<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tunex1<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tunex2<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tuney1<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tuney2<<"  "
            <<setw(25)<<setprecision(18)<<scientific<<tunex1-tunex2<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tuney1-tuney2<<endl;
          f1.close();
        }
       } 
    } 
  }
}

void Track_tbt_tune_footprint( Line & linename, double deltap0, double emitx, double emity)  // I prefer RF off for this 
{
  double betax0, alfax0, betay0, alfay0;

  int i,j,k;
  double nsigma, nsigmax, nsigmay;
  double xn, pxn, yn, pyn, x[6];

  int Nturn=2048;
  int stable = 1, lost_turn=0, lost_post=0;
  double  x_tbt[Nturn*6];
  double xtemp[1024];
  double tunex1, tunex2, tuney1, tuney2;

  fstream  f1;
  f1.open("./footprint-output.dat", ios::out);
  f1.close();

  Cal_Twiss(linename, deltap0);
  betax0= linename.Cell[linename.Ncell-1]->Beta1; 
  betay0= linename.Cell[linename.Ncell-1]->Beta2;
  alfax0= linename.Cell[linename.Ncell-1]->Alfa1; 
  alfay0= linename.Cell[linename.Ncell-1]->Alfa2;

  Cal_Orbit_Num(linename,deltap0);

  for( i=0;i<5;i++)
    for(j=0;j<8;j++){

      nsigma=  0.+j*(7-0.)/7 ;
      nsigmax= nsigma * cos((i+1)*15 *PI/180);
      nsigmay= nsigma * sin((i+1)*15 *PI/180);

      xn=  sqrt(emitx) * nsigmax ; 
      pxn= 0.;
      yn=  sqrt(emity) * nsigmay ; 
      pyn= 0;
      
      x[0]= sqrt(betax0) * xn ;                                          
      x[1]=-alfax0 * xn / sqrt( betax0 ) + 1.0* pxn /  sqrt( betax0 );   
      x[2]= sqrt(betay0) * yn ;                                          
      x[3]=-alfay0 * yn / sqrt( betay0 ) + 1.0* pyn /  sqrt( betay0 );   
      x[4]= 0. ;
      x[5]=  deltap0;
      x[0]= x[0] + linename.Cell[linename.Ncell-1]->X[0]  ;
      x[1]= x[1] + linename.Cell[linename.Ncell-1]->X[1]  ;
      x[2]= x[2] + linename.Cell[linename.Ncell-1]->X[2]  ;
      x[3]= x[3] + linename.Cell[linename.Ncell-1]->X[3]  ; 

      if(x[0] == 0. ) x[0]=x[0]+1.0e-09;
      if(x[2] == 0. ) x[2]=x[2]+1.0e-09; 
      
      stable = 1; lost_turn=0; lost_post=0;
      Track_tbt(linename, x, Nturn, x_tbt, stable, lost_turn, lost_post);
      if(stable==1) {
	for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+0];
	FineTuneFinder(1024,xtemp,tunex1);
	for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+2];
	FineTuneFinder(1024,xtemp,tuney1);
	for(k=0;k<1024;k++) xtemp[k]=x_tbt[(k+1023)*6+0];
	FineTuneFinder(1024,xtemp,tunex2);
	for(k=0;k<1024;k++) xtemp[k]=x_tbt[(k+1023)*6+2];
	FineTuneFinder(1024,xtemp,tuney2);
	f1.open("./footprint-output.dat",  ios::out | ios::app);
	f1<<(i+1)*15<<" "<<j<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tunex1<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tunex2<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tuney1<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tuney2<<"  "
	  <<setw(25)<<setprecision(18)<<scientific<<tunex1-tunex2<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tuney1-tuney2<<endl;
	f1.close();
      }
    }
}

void Track_DA( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0 ) // Track with ->Pass()
{
  int i;
  double angle, sigma_s, sigma_e, sigma, sigma1, sigma2, sigma_max=100.0;  
  int    stable, stable1=1, stable2=1,  lost_turn=0, lost_post=0; 
  FILE *f2;
  f2=fopen("./DA-output.dat","w");
  fclose(f2);
  double x[6];

  for(i=0; i<10;i++){
    
    angle=( 90/11.)*(i+1)* PI/180;
    sigma_s=0.;
    sigma_e=sigma_max;
    
    do {
      sigma = (sigma_s + sigma_e )/2.0;
      sigma1= sigma;
      sigma2= sigma + 0.05;

      x[0]= sigma1 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma1 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable1=1; lost_turn=0; lost_post=0;
      Track(linename, x, nturn, stable1, lost_turn, lost_post);
      
      x[0]= sigma2 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma2 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable2=1; lost_turn=0; lost_post=0;
      Track(linename, x, nturn, stable2, lost_turn, lost_post);
      
      if ( stable1* stable2 ==1 ) {
	stable =1 ; }
      else{
	stable=0;
      }
      if( stable ==0 ) 	sigma_e= sigma;
      if( stable ==1 ) 	sigma_s= sigma;
    } while ( abs( sigma_e- sigma_s) > 0.1  );

    f2=fopen("DA-output.dat" ,"a");
    fprintf(f2,"%f %f \n", ( 90/11.)*(i+1),  sigma);
    fclose(f2);
  }
} 

void Track_Fast_DA( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0 ) // Track without class
{
  int i;
  double angle, sigma_s, sigma_e, sigma, sigma1, sigma2, sigma_max=100.0;  
  int    stable, stable1=1, stable2=1,  lost_turn=0, lost_post=0; 
  double sum_x2, sum_y2, sum_z2;
  FILE *f2;
  f2=fopen("./DA-output.dat","w");
  fclose(f2);
  double x[6];

  for(i=0; i<10;i++){

    angle=( 90/11.)*(i+1)* PI/180;
    sigma_s=0.;
    sigma_e=sigma_max;
    
    do {
      sigma = (sigma_s + sigma_e )/2.0;
      sigma1= sigma;
      sigma2= sigma + 0.05;

      x[0]= sigma1 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma1 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable1=1; lost_turn=0; lost_post=0;
      Track_Fast(linename, x, nturn, stable1, lost_turn, lost_post, sum_x2, sum_y2, sum_z2);
      
      x[0]= sigma2 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma2 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable2=1; lost_turn=0; lost_post=0;
      Track_Fast(linename, x, nturn, stable2, lost_turn, lost_post, sum_x2, sum_y2, sum_z2);
      
      if ( stable1* stable2 ==1 ) {
	stable =1 ; }
      else{
	stable=0;
      }
      if( stable ==0 ) 	sigma_e= sigma;
      if( stable ==1 ) 	sigma_s= sigma;
    } while ( abs( sigma_e- sigma_s) > 0.1  );

    f2=fopen("DA-output.dat" ,"a");
    fprintf(f2,"%f %f \n", ( 90/11.)*(i+1),  sigma);
    fclose(f2);
  }
} 
