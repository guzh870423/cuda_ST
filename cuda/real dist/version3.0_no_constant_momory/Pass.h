#include "constants.h"
//1.0/(2*(2-pow(2.0,1.0/3)))
#define _CFdrift1 0.675603596 
//0.5-CFdrift1
#define _CFdrift2 -0.175603596 
//2*CFdrift1
#define _CFkick1 1.351207192 
//1.e0 - 2e0 * CFkick1
#define _CFkick2 -1.702414384 



using namespace std;

enum ps_index { x_ = 0, px_ = 1, y_ = 2, py_ = 3,  z_ = 4, delta_ = 5 };  //  z > 0 earlier than synchronous particle


class COORD
{
  public:
	REAL x[6];
};

class Element
{
  public:

  	REAL   L,DX, DY;




};

//---------------DRIFT-----------------------------------
class DRIFT: public Element
{

 public:

	__host__ void init(REAL l=0)
	{
		L=l;

	}


};

//--------------------------SBEND------------------------------------
class SBEND: public Element
{
    public:
	REAL ANGLE,E1,E2;
	int Nint;
	
	 void init(REAL l=0, REAL angle=0, REAL e1=0, REAL e2=0, int N_kicks=1)
	{
		L=l;
		ANGLE=angle;
        		E1=e1;
         	E2=e2;
	 	Nint=N_kicks;

	}

};


//---------------------QUAD-----------------------------------------
class QUAD: public Element
{
 public:
    REAL K1L,K1SL;
     int Nint;
    void init(REAL l=0, REAL k1=0, REAL tilt=0, int N_kicks=1 )
    	{
		L=l;
		K1L=k1*l*cos(2*tilt);
		K1SL=k1*l*sin(2*tilt);
	 	Nint=N_kicks;

	}
};

//---------------------SEXT-----------------------------------------
class SEXT: public Element
{
 public:
    REAL K2L,K2SL;
     int Nint;
    void init(REAL l=0, REAL k2=0, REAL tilt=0, int N_kicks=1 )
    	{
		L=l;
		K2L=k2*l*cos(3*tilt);
		K2SL=k2*l*sin(3*tilt);
	 	Nint=N_kicks;

	}
};

//-------------------------------CMB_SBEND--------------------------------------
class CMB_SBEND: public Element
{
 public:
   REAL ANGLE,E1,E2,KNL[11];
   int Nint,Norder;
    void init( REAL l=0, REAL angle=0, REAL e1=0, REAL e2=0, int N_kicks=1, REAL kn[11]=NULL)
    { 
      int i;

	  L=l;
	  ANGLE=angle;
          E1=e1;
          E2=e2;
	  Nint=N_kicks;
          for(i=0;i<11;i++) {
	    KNL[i]=kn[i]*l;  }
	  Norder=1;
	  for(i=0;i<10;i++) {
	    if( KNL[10-i] != 0. ){
	      Norder=10-i;
	      break;
	    }
	  }
	}
};

//----------------------------RFCAV--------------------------------
class RFCAV: public Element
{
 public:
 void init( REAL l, REAL vrf, REAL frf, REAL phase0)
    { 

      VRF= vrf; 
      FRF= frf;
      PHASE0=phase0;
      L=l;
    }


  REAL VRF, FRF, PHASE0;
};

template <class T> 
__device__ void RFCAV_Pass(T x[6], RFCAV *CA)
{
	T u,L;
//  DRIFT_Pass(x, L/2.);
    L=CA->L/2;
    u=L/(1+x[delta_]);
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
//RF Kick
  x[5] =x[5] + (CA->VRF/E0)*sin(2.0*M_PI*CA->FRF*x[z_]/v_c+CA->PHASE0);
//  DRIFT_Pass(x, L/2.);
    L=CA->L/2;
    u=L/(1+x[delta_]);
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
}




//=================================================================
//
//    4th order Symplectic Integrator  (S.I.)
//
//=================================================================

//---------------------DRIFT_PASS------------------------
template <class T>
__device__ void DRIFT_Pass(T x[6], DRIFT *D)
{
  T u;
  if( D->L !=0. ) {
    u=D->L/(1+x[delta_]);
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
  }
}

//--------------------------SBEND_PASS---------------------
template <class T> 
__device__ void SBEND_Pass(T x[6], SBEND *B) 
{
  int i;
  T href=B->ANGLE/B->L;
  T Lint=B->L/B->Nint;
  T  u,L;
  
  x[1] = x[1]+ tan(B->E1)*x[0]*href;   
  x[3] = x[3]- tan(B->E1)*x[2]*href;

  for(i=0;i<B->Nint;i++){
    //DRIFT_Pass(x,Fdrift1*Lint);
    L=_CFdrift1*Lint;
    u=L/(1+x[delta_]);
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
    
    //bend_kick_pass(x, Fkick1*Lint, href);
    L=_CFkick1*Lint;
    x[px_]=x[px_]+(href*x[delta_]-href*href*x[x_])*L;
    x[z_]=x[z_]-href*x[x_]*L;
    
    //DRIFT_Pass(x,Fdrift2*Lint);
    L=_CFdrift2*Lint;
    u=L/(1+x[delta_]);
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
    
    //bend_kick_pass(x, Fkick2*Lint, href);
    L=_CFkick2*Lint;
    x[px_]=x[px_]+(href*x[delta_]-href*href*x[x_])*L;
    x[z_]=x[z_]-href*x[x_]*L;
    
    //DRIFT_Pass(x,Fdrift2*Lint);
    L=_CFdrift2*Lint;
    u=L/(1+x[delta_]);
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
    
    //bend_kick_pass(x, Fkick1*Lint, href);
    L=_CFkick1*Lint;
    x[px_]=x[px_]+(href*x[delta_]-href*href*x[x_])*L;
    x[z_]=x[z_]-href*x[x_]*L;
    
    //DRIFT_Pass(x,Fdrift1*Lint); 
    L=_CFdrift1*Lint;
    u=L/(1+x[delta_]);
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
  }

  x[1] = x[1]+ tan(B->E2)*x[0]*href;   
  x[3] = x[3]- tan(B->E2)*x[2]*href;
}

template <class T> 
__device__ void QUAD_Pass(T x[6], QUAD *Q)
{
  int i;
  T Lint=Q->L/Q->Nint;
  

      T k1l_kick1,k1sl_kick1;
      T k1l_kick2,k1sl_kick2;
      k1l_kick1 =_CFkick1*Q->K1L/Q->Nint;
      k1sl_kick1=_CFkick1*Q->K1SL/Q->Nint;
      k1l_kick2 =_CFkick2*Q->K1L/Q->Nint;
      k1sl_kick2=_CFkick2*Q->K1SL/Q->Nint;

      int j, Norder=1;
      T  Xn,  Yn, Xn0, Yn0;
      T  By, Bx;
      T  u,k1l,k1sl,L;     
      for(i=0;i<Q->Nint;i++){
	
	//DRIFT_Pass(x,Fdrift1*Lint);
	L=_CFdrift1*Lint;
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
	L=_CFdrift2*Lint;
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
	L=_CFdrift2*Lint;
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
	L=_CFdrift1*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
      }
 
} 

template <class T> 
void __device__ SEXT_Pass(T x[6], SEXT *S)
{
  int i;
  T Lint=S->L/S->Nint;

      T k2l_kick1,k2sl_kick1;
      T k2l_kick2,k2sl_kick2;
      k2l_kick1 =_CFkick1*S->K2L/S->Nint;
      k2sl_kick1=_CFkick1*S->K2SL/S->Nint;
      k2l_kick2 =_CFkick2*S->K2L/S->Nint;
      k2sl_kick2=_CFkick2*S->K2SL/S->Nint;
      
      int j, Norder=2;
      T  Xn,  Yn, Xn0, Yn0;
      T  By, Bx;
      T  u,L,k2l,k2sl;  
      for(i=0;i<S->Nint;i++){
	
	//DRIFT_Pass(x,Fdrift1*Lint);
	L=_CFdrift1*Lint;
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
	L=_CFdrift2*Lint;
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
	L=_CFdrift2*Lint;
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
	L=_CFdrift1*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);
      }

} 

template <class T>
void __device__ bend_mult_kick_pass(T x[6], T L, T href, int Norder, T KNL[11])
{
  int i;
  int fac=1;
  T   Xn,  Yn, Xn0, Yn0;
  T   By, Bx;

  By=KNL[0];
  Bx=0;
  Xn=1.;
  Yn=0.;

  for(i=1;i<Norder+1;i++){
    Xn0=Xn;
    Yn0=Yn;
    Xn=Xn0*x[x_]-Yn0*x[y_];
    Yn=Xn0*x[y_]+Yn0*x[x_];
    fac=fac*i;
    if ( KNL[i] != 0. ) {
      By=By+(KNL[i]*Xn-0*Yn)/fac;
      Bx=Bx+(KNL[i]*Yn+0*Xn)/fac;
    }
  }
  x[px_]=x[px_]-By  + (href*x[delta_]-href*href*x[x_])*L;;
  x[py_]=x[py_]+Bx;
  x[z_]=x[z_]-href*x[x_]*L;
}

template <class T> 
void __device__ CMB_SBEND_Pass(T x[6], CMB_SBEND *C, curandState *s,T *dqueue) 
{

       int i;
       T href=C->ANGLE/C->L,Es1=x[0]*C->E1 / (1+x[0]*href),Es2=0;
       T Lint=C->L/C->Nint;
       T knl_kick1[11];
       T knl_kick2[11];
	T L,u;

       for(i=0;i<11;i++) {
	 knl_kick1[i] =_CFkick1*C->KNL[i]/C->Nint;

       }
       for(i=0;i<11;i++) {
	 knl_kick2[i] =_CFkick2*C->KNL[i]/C->Nint;
       }
       
//edge angle
       x[1] = x[1]+ tan(C->E1)*x[0]*href;   
       x[3] = x[3]- tan(C->E1)*x[2]*href;
	if(Es1<0)
	{
		radiation(x, href, -Es1, s, dqueue);	
	}
       
      for(i=0;i<C->Nint;i++){
	Es2=x[0]*C->E2 / (1+x[0]*href);

//	 DRIFT_Pass(x,Fdrift1*Lint);
	L=_CFdrift1*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);

	 bend_mult_kick_pass(x, (REAL)_CFkick1*Lint, href, C->Norder,knl_kick1);
	if( Es1>i * Lint ||  Es2>Lint* (0.6666666666f + C->Nint -1 - i) )
	{
		if(Es1>i * Lint && Es1<Lint * (i + 0.333333333f) )
		{
			radiation(x, href, Lint*(i + 0.333333333f)-Es1, s, dqueue);
		}

		if(Es2>Lint* (0.6666666666f + C->Nint -1 - i ) && Es2<Lint * (1+C->Nint -1 - i))
		{
			radiation(x, href, Lint* (1+C->Nint -1 - i)-Es2, s, dqueue);
		}
	}
	else
	{
		radiation(x, href, Lint/3, s, dqueue);
	}
//	 DRIFT_Pass(x,Fdrift2*Lint);
	L=_CFdrift2*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);

	 bend_mult_kick_pass(x, (REAL)_CFkick2*Lint, href, C->Norder,knl_kick2);
	if( Es1> Lint * (i+0.33333333f) ||  Es2>Lint* (0.33333333f + C->Nint -1 - i ) )
	{
		if(Es1> Lint * (i+0.33333333f)&& Es1<Lint * (i + 0.666666666f) )
		{
			radiation(x, href, Lint*(i + 0.666666666f)-Es1, s, dqueue);
		}

		if(Es2>Lint* (0.33333333f + C->Nint -1 - i) && Es2<Lint * (0.666666666f+C->Nint -1 - i) )
		{
			radiation(x, href, Lint* (0.666666666f+C->Nint -1 - i)-Es2, s, dqueue);
		}
	}
	else
	{
		radiation(x, href, Lint/3, s, dqueue);
	}

//	 DRIFT_Pass(x,Fdrift2*Lint);
	L=_CFdrift2*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);

	 bend_mult_kick_pass(x, (REAL)_CFkick1*Lint, href,  C->Norder,knl_kick1);
	if( Es1> Lint * (i+0.666666666f) ||  Es2>Lint* ( C->Nint -1 - i ) )
	{
		if(Es1> Lint * (i+0.666666666f)&& Es1<Lint * (i + 1) )
		{
			radiation(x, href, Lint*(i + 1)-Es1, s, dqueue);
		}

		if(Es2>Lint* (C->Nint -1 - i) && Es2<Lint * (0.33333333f+C->Nint -1 - i) )
		{
			radiation(x, href, Lint* (0.33333333f+C->Nint -1 - i)-Es2, s, dqueue);
		}
	}
	else
	{
		radiation(x, href, Lint/3, s, dqueue);
	}

//	 DRIFT_Pass(x,Fdrift1*Lint);
	L=_CFdrift1*Lint;
	u=L/(1+x[delta_]);
	x[x_]=x[x_]+x[px_]*u;
	x[y_]=x[y_]+x[py_]*u;
	x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[delta_]);

       }

//Edge angle 2
       x[1] = x[1]+ tan(C->E2)*x[0]*href;   
       x[3] = x[3]- tan(C->E2)*x[2]*href;
	if(Es2<0)
	{
		radiation(x, href, -Es2, s, dqueue);	
	} 

} 

template <class T> 
void __device__ radiation(T x[6], T href, T dsISR, curandState *s,T *dqueue)
{
	T ddel, A_ph,n_ph;
	int j=0,rand_num;
	T  f = (1+x[0]*href)/sqrt((1+x[5])*(1+x[5])-x[1]*x[1]-x[3]*x[3]);
	T  xp = x[1]*f;
	T  yp = x[3]*f;
	T  dsFactor = sqrt((1+x[0]*href)*(1+x[0]*href) + xp*xp + yp*yp);
		A_ph=1.44337567297 * _alpha * _gamma * dsISR * href *dsFactor  ; //1.4433756729740644112728719512549=5/2/sqrt(3)


		n_ph = curand_poisson(s, A_ph);

		for(j=0,ddel=0; j < n_ph; j++)
		{
			rand_num= floor(curand_uniform(s)*_pool);

			ddel=dqueue[rand_num] * (1+x[5]) * (1+x[5]) * href;
			x[1]=(1+x[5]-ddel) * x[1] /(1+ x[5]);
			x[3]=(1+x[5]-ddel) * x[3] /(1+ x[5]);
			x[5]-=ddel;

		}



}


