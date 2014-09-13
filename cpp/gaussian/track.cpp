#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/piecewise_linear_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include "Pass.h"
#include "constants.h"
#include  "lattice_original.h"




using namespace std;




  	const gsl_rng_type * T;
  	gsl_rng * r;


void radiation( double *x, double rho);
void emittance(COORD *y, double *PEx, double *PEy, double *PEdelta)
{

	double avg_x=0,avg_xp=0,avg_y=0,avg_yp=0,avg_delta=0,sig_xx=0,sig_xpxp=0,sig_xxp=0,sig_yy=0,sig_ypyp=0,sig_yyp=0,sig_delta=0;
	for(int i=0;i<_Npart;i++)
	{
		avg_x+=y[i].x[x_]/_Npart;
		avg_xp+=y[i].x[px_]/(1+y[i].x[delta_])/_Npart;
		avg_y+=y[i].x[y_]/_Npart;
		avg_yp+=y[i].x[py_]/(1+y[i].x[delta_])/_Npart;
		avg_delta+=y[i].x[delta_]/_Npart;
	}

	for(int i=0;i<_Npart;i++)
	{
		sig_xx+=(y[i].x[x_]-avg_x)*(y[i].x[x_]-avg_x)/_Npart;
		sig_xpxp+=(y[i].x[px_]/(1+y[i].x[delta_])-avg_xp)*(y[i].x[px_]/(1+y[i].x[delta_])-avg_xp)/_Npart;
		sig_xxp+=(y[i].x[x_]-avg_x)*(y[i].x[px_]/(1+y[i].x[delta_])-avg_xp)/_Npart;
		sig_yy+=(y[i].x[y_]-avg_y)*(y[i].x[y_]-avg_y)/_Npart;
		sig_ypyp+=(y[i].x[py_]/(1+y[i].x[delta_])-avg_yp)*(y[i].x[py_]/(1+y[i].x[delta_])-avg_yp)/_Npart;
		sig_yyp+=(y[i].x[y_]-avg_y)*(y[i].x[py_]/(1+y[i].x[delta_])-avg_yp)/_Npart;
		sig_delta+=(y[i].x[delta_]-avg_delta)*(y[i].x[delta_]-avg_delta)/_Npart;
	}
	
	*PEx=sqrt(sig_xx*sig_xpxp-sig_xxp*sig_xxp);
	*PEy=sqrt(sig_yy*sig_ypyp-sig_yyp*sig_yyp);
	*PEdelta=sqrt(sig_delta);
}

int main(int argc, char** argv)
{





//initialization

//**********************************************************


  	gsl_rng_env_setup();

  	T = gsl_rng_default;
 	r = gsl_rng_alloc (T);
//**********************************************************




  	COORD *part;

	part=(COORD*)malloc(_Npart*sizeof(COORD));

	
	
	double phi_x,phi_y,Jx,Jy,Ex,Ey,Sdelta;
	int i,n;

	for(i=0;i<_Npart;i++)
	{
	     
	     do {Jx=gsl_ran_exponential(r, 2*E_x);}
		while(Jx>E_x*6);
	     do {Jy=gsl_ran_exponential(r, 2*E_y);}
		while(Jy>E_y*6);
		phi_x=gsl_ran_flat(r,0,2*M_PI);
		phi_y=gsl_ran_flat(r,0,2*M_PI);

		part[i].x[x_]=sqrt(Jx*Beta_x)*cos(phi_x);
		part[i].x[px_]=sqrt(Jx/Beta_x)*sin(phi_x);
		part[i].x[y_]=sqrt(Jy*Beta_y)*cos(phi_y);
		part[i].x[py_]=sqrt(Jy/Beta_y)*sin(phi_y);
		part[i].x[z_]=0;
		part[i].x[delta_]=0.00;			
	
	}
//	part[0].x[0]=0.000;part[0].x[1]=0.000;part[0].x[2]=0.000;part[0].x[3]=0.000;part[0].x[5]=0.01;




       ofstream outfile("abc.txt");
        outfile.close();

	for(n=0;n<_Nturn1;n++)
	{


		Track(part);

	emittance(part,&Ex,&Ey,&Sdelta);

       ofstream outfile("abc.txt",ios::app);
     outfile<<n<<"    "<<Ex<<"    "<<Ey<<"    "<<Sdelta<<endl;

 //     outfile<<part[0].x[0]<<"  "<<part[0].x[1]<<"  "<<part[0].x[2]<<"  "<<part[0].x[3]<<"  "<<part[0].x[4]<<"  "<<part[0].x[5]<<endl;
 //	for(int k=0;k<_Npart;k++) {  if(abs(part[k].x[0])>10||abs(part[k].x[1])>10||part[k].x[2]>10||part[k].x[3]>10||part[k].x[5]>2) {cout<<n<<"  "<<part[k].x[0]<<"  "<<part[k].x[1]<<"  "<<part[k].x[2]<<"  "<<part[k].x[3]<<"  "<<part[k].x[4]<<"  "<<part[k].x[5]<<endl;}  }
        outfile.close();

	}


	cout<<part[0].x[0]<<"  "<<part[0].x[1]<<"  "<<part[0].x[2]<<"  "<<part[0].x[3]<<"  "<<part[0].x[4]<<"  "<<part[0].x[5]<<endl;
 

	free(part);

	gsl_rng_free (r);


//delta.close();
}
void radiation(double x[6], double rho, double L, double dsISR)
{
	double ddel, rand_num,sigma;
	int j=0;
	double  f = (1+x[0]/rho)/sqrt((1+x[5])*(1+x[5])-x[1]*x[1]-x[3]*x[3]);
	double  xp = x[1]*f;
	double  yp = x[3]*f;
	double  dsFactor = sqrt((1+x[0]/rho)*(1+x[0]/rho) + xp*xp + yp*yp);

//average relative energy loss
	ddel=_beta*_beta*_beta * C_gamma /2/3.1415926 *E0*E0*E0 * (1+x[5])*(1+x[5]) /rho/rho * L * dsFactor;

//gaussian rms
	sigma = sqrt( 0.315865513*C_gamma*hbar*v_c*_gamma*_gamma*_gamma *E0*E0/(rho*rho*rho)/_beta*dsISR*dsFactor) * (1+x[5])*(1+x[5]); 
//	do {rand_num=gsl_ran_gaussian(r,sigma);} while(abs(rand_num)>9*sigma);
	rand_num=gsl_ran_gaussian(r,sigma);
	ddel+=rand_num;
	x[1]=(1+x[5]-ddel) * x[1] /(1+ x[5]);
	x[3]=(1+x[5]-ddel) * x[3] /(1+ x[5]);
	x[5]-=ddel;

		


} 

