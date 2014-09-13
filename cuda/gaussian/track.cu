#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <curand.h>
#include <curand_kernel.h>
#include "Pass.h"
#include "constants.h"
#include "lattice_PAR.h"



using namespace std;



void emittance(COORD *y, REAL *PEx, REAL *PEy, REAL *PEdelta)
{

	REAL avg_x=0,avg_xp=0,avg_y=0,avg_yp=0,avg_delta=0,sig_xx=0,sig_xpxp=0,sig_xxp=0,sig_yy=0,sig_ypyp=0,sig_yyp=0,sig_delta=0;
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
  	const gsl_rng_type * T;
  	gsl_rng * r;

  	gsl_rng_env_setup();

  	T = gsl_rng_default;
 	r = gsl_rng_alloc (T);

  	COORD *part, *dpart;
	int size = _Npart * sizeof(COORD);
	part=(COORD*)malloc(size);
	cudaMalloc(&dpart,size);
	
	
	REAL phi_x,phi_y,Jx,Jy,Ex,Ey,Sdelta;
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
//	part[0].x[0]=0.000;part[0].x[1]=0.000;part[0].x[2]=0.000;part[0].x[3]=0.000;part[0].x[5]=0.00;


	cudaMemcpy(dpart,part,size,cudaMemcpyHostToDevice);

       ofstream outfile("abc.txt");
        outfile.close();
	for(n=0;n<_Nturn1;n++)
	{
		Track<<<_BlockNum,_ThreadNum>>>(dpart,n);

		cudaMemcpy(part,dpart,size,cudaMemcpyDeviceToHost);
	emittance(part,&Ex,&Ey,&Sdelta);

       ofstream outfile("abc.txt",ios::app);
     outfile<<n<<"    "<<Ex<<"    "<<Ey<<"    "<<Sdelta<<endl;
 //     outfile<<part[0].x[0]<<"  "<<part[0].x[1]<<"  "<<part[0].x[2]<<"  "<<part[0].x[3]<<"  "<<part[0].x[4]<<"  "<<part[0].x[5]<<endl;
//	for(int k=0;k<_Npart;k++) {  if(abs(part[k].x[0])>10||abs(part[k].x[1])>10||part[k].x[2]>10||part[k].x[3]>10||part[k].x[5]>2) {cout<<n<<"  "<<part[k].x[0]<<"  "<<part[k].x[1]<<"  "<<part[k].x[2]<<"  "<<part[k].x[3]<<"  "<<part[k].x[4]<<"  "<<part[k].x[5]<<endl;}  }
        outfile.close();

	}


	cout<<part[0].x[0]<<"  "<<part[0].x[1]<<"  "<<part[0].x[2]<<"  "<<part[0].x[3]<<"  "<<part[0].x[4]<<"  "<<part[0].x[5]<<endl;
 

	free(part);
	cudaFree(dpart);

	gsl_rng_free (r);

}

